# AGENTS.md

This file provides guidance to AI coding agents when working with code in this repository.

## Build and Test Commands

```bash
# Build the project
cargo build

# Run tests
cargo test

# Run a single test
cargo test <test_name>

# Run examples (visualization requires Rerun viewer running)
cargo run --example draw_faces
cargo run --example draw_many
cargo run --example draw_shapes
cargo run --example floor_plan
cargo run --example ray_2_boxes
cargo run --example ray_teapot
cargo run --example bench_teapot
cargo run --example sim_acoustics
cargo run --example sim_lighting
cargo run --example sim_lighting_heatmap
cargo run --example sim_energy
cargo run --example bestest_600_energy
cargo run --example bestest_energy_suite
cargo run --example bras_cr2
cargo run --example pipeline_solar_energy
cargo run --example pipeline_weather_shaded_energy

# Check code without building
cargo check

# Format code
cargo fmt

# Lint
cargo clippy
```

## Architecture

This is a 3D building modeling and simulation library with a strict hierarchical composition model:

```
Building → Zone → Solid → Wall → Polygon → Mesh
```

Each entity has:
- `name: String` (path separator `/` not allowed)
- `uid: UID` (UUID v4 wrapper)
- `parent: Option<UID>` (reference to parent)

### Core Geometry Types (src/geom/)

- **Point**: 3D point (x, y, z as f64) with arithmetic operations against Vector
- **Vector**: 3D displacement (dx, dy, dz) with cross/dot products, normalization
- **Mesh**: vertices (`Vec<Point>`) + optional faces (`Vec<TriangleIndex>`)
- **Polygon**: planar multi-sided face, auto-triangulates via ear-clipping
- **Wall**: container of polygons (`HashMap<String, Polygon>`)
- **Solid**: 3D volume containing walls; constructable via `from_box()` or `from_floor_plan(FloorPlan)`
- **Zone**: container of solids, groups related volumes
- **Building**: top-level container of zones
- **PlaneBasis** (`projection.rs`): orthonormal basis for projecting 3D points onto a 2D plane and back

### Key Traits

- **HasMesh**: provides `copy_mesh()` for polymorphic mesh access (used by drawing)
- **HasName**: polymorphic name access
- **SortByName**: in-place stable sort by name
- **IsClose**: epsilon-based f64/Point/Vector comparison (uses `EPS = 1e-13`)

### Transformations

`rotate()` and `translate()` methods cascade down the hierarchy. Rotation uses Rodrigues' formula (src/geom/rotation.rs).

### Simulation (src/sim/)

The simulation layer provides ray-based and thermal simulation capabilities, sharing a common engine, material system, and composable framework.

#### Simulation Framework (`sim/framework/`)

Domain-agnostic runtime for composing simulation modules:
- **Bus**: typed message bus for inter-module communication (any `'static` type can be published/consumed)
- **SimContext**: shared read-only context holding `&Building` and `&SurfaceIndex`
- **SimModule**: trait for composable simulation modules (`init()` + `step()` lifecycle)
- **Pipeline**: sequences `SimModule`s with shared `Bus` and `SimContext`

#### Surface Index (`sim/index/`)

- **SurfaceIndex**: domain-agnostic lookup table for polygon surfaces keyed by `UID`; provides path-based and UID-based access to `SurfaceRef` metadata (path, polygon UID, zone UID, area)
- **SurfaceRef**: reference information for a polygon surface in the building hierarchy

#### Surface Semantics (`sim/surfaces.rs`)

Cross-domain surface classification using the polygon-facing graph:
- **SurfaceKind**: `Exterior`, `SameZoneInterface`, `InterZoneInterface`, `UnknownInterface`
- **SurfaceInterface**: tracks same-zone and cross-zone relationships
- **SurfaceSemantics**: classifies all polygon surfaces; used by energy, lighting, and acoustics domains consistently

#### Coupling Payloads (`sim/coupling.rs`)

Stable contracts for inter-module communication via Bus (radiometric SI units):
- **WeatherHourIndex**: step counter for weather datasets
- **OutdoorAirTemperatureC**: outdoor air temperature [C]
- **InternalGainsWPerZone** / **InternalGainsWTotal**: heat gains [W] keyed by zone UID
- **ShortwaveAbsorbedWPerPolygon**: absorbed solar on surfaces [W]
- **ShortwaveTransmittedWPerZone**: transmitted shortwave into zones [W]

#### Materials (`sim/materials.rs`)

Unified material system with domain-specific properties:
- **AcousticMaterial**: frequency-dependent absorption/scattering (6 octave bands: 125-4000 Hz)
- **OpticalMaterial**: diffuse/specular reflectance and transmittance (RGB)
- **ThermalMaterial**: U-value, thermal capacity, construction layers
- **Material**: composite type holding optional acoustic, optical, and thermal properties
- **MaterialLibrary**: named material registry with path-pattern assignment to building surfaces; includes presets (concrete, glass, gypsum, carpet, wood, metal)

#### Simulation Engine (`sim/engine/`)

Shared infrastructure for ray-based simulations:
- **FlatScene**: flattened polygon scene from building hierarchy with voxel-grid spatial acceleration, transparent surface detection, and ray-surface intersection queries
- **RayBatch** / **RayState**: batch of rays with position, velocity, energy
- **VoxelGrid**: spatial hash grid for fast polygon lookup
- Pluggable models via traits:
  - **AbsorptionModel**: `ScalarAbsorption`, `FrequencyDependentAbsorption`, `AirAbsorption`
  - **ReflectionModel**: `Specular`, `Diffuse`, `Hybrid`
  - **PropagationModel**: `FixedTimeStep`, `EventDriven`

#### Acoustic Ray Tracing (`sim/rays/`)

- **SimulationConfig**: ray count, speed, source/absorber positions, scalar or frequency-dependent mode, material library integration, `store_ray_history` flag (set `false` to avoid ~800 MB RAM for large runs)
- **Simulation**: builds `FlatScene`, runs time-stepped ray tracing with configurable absorption/reflection/propagation
- **SimulationResult**: per-step ray positions/energies (optional, controlled by `store_ray_history`), absorber hits (scalar and per-band)

#### Acoustics (`sim/acoustics/`)

- **Receiver**: spherical receiver accumulating energy-time-frequency histograms
- **ImpulseResponse**: extracted from receiver data (per-band and broadband)
- **Metrics**: Schroeder backward integration for reverberation time (RT20/RT30/EDT), clarity (C50/C80), definition (D50), STI; `RoomAcousticReport` aggregates all metrics
- **Source directivity**: `Omnidirectional`, `Cardioid` patterns via `SourceDirectivity` trait
- **Auralization** (`auralization/`): WAV export of impulse responses (`write_ir_wav`), convolution with dry audio (`auralize`, `auralize_per_band`), octave-band filtering, resampling utilities

#### Lighting (`sim/lighting/`)

- **LightingConfig**: point/directional lights, ray count, bounce limit, material library
- **LightingSimulation**: forward ray tracing with RGB illuminance accumulation
- **LightingModule**: `SimModule` wrapper for ray-tracing in composed pipelines
- **Light sources**: `PointLight`, `DirectionalLight` via `LightSource` trait
- **SensorGrid**: auto-generated sensor points on polygon surfaces
- **Sky models**: `CIEOvercast`, `CIEClearSky` via `SkyModel` trait
- **Solar**: `SolarPosition` calculation from latitude/longitude/time
- **Backward ray tracing**: `BackwardTracer` for sensor-to-sky illuminance
- **LightingResult**: per-polygon RGB illuminance map (UID-keyed)
- **Shortwave coupling** (`shortwave.rs`): modules converting lighting/solar results to thermal coupling payloads:
  - `LightingToShortwaveModule`: converts `LightingResult` to shortwave payloads
  - `SolarShortwaveStepModule` / `SolarShortwaveModule`: deterministic SHGC-based solar (step-based or single-shot)
  - `SolarShortwaveShadedStepModule`: adds direct-sun occlusion via FlatScene ray casting
  - `SolarEpwModule` / `SolarEpwBusModule`: EPW-driven solar with bus integration
  - `SolarEpwShadedModule` / `SolarEpwShadedBusModule`: EPW-driven solar with shadowing

#### Energy (`sim/energy/`)

- **ThermalConfig**: wall constructions, U-values, indoor/outdoor temperatures, infiltration, gains
- **WallConstruction**: layered wall constructions with R-value/U-value calculation
- **WeatherData**: EPW (EnergyPlus Weather) file parser with hourly records
- **InternalGainsProfile**: hourly/daily schedules for occupancy, equipment, lighting
- **HvacIdealLoads**: simplified HVAC controller with heating/cooling setpoints
- **SolarBridge**: coupling between lighting and energy simulations (legacy)
- **run_annual_simulation()**: 8760-hour annual simulation producing heating/cooling demand, peak loads, monthly breakdown
- **ThermalResult** / **AnnualResult**: simulation output types
- **ThermalBoundaries** (`boundary.rs`): re-exports `SurfaceSemantics` as thermal-specific type alias
- **EnergyModule** (`module.rs`): `SimModule` wrapper for step-based multi-zone thermal simulation
  - `EnergyModelKind`: `AirOnly` or `EnvelopeRc2R1C`
  - `EnergyModuleConfig`: combines ThermalConfig, HVAC, timestep, model selection
- **WeatherModule** (`weather_module.rs`): publishes weather data to Bus per timestep
- **InternalGainsModule** (`gains_module.rs`): publishes time-varying gains from profile to Bus
- **MultiZoneRecorderData** (`recorder.rs`): buffer for annual simulation results
- **Thermal Network** (`network/`): graph-based thermal representation
  - `ThermalNetwork`: inter-zone and exterior conductances (W/K)
  - `MultiZoneAirModel`: zone-air-only thermal model with Gaussian elimination solver
  - `MultiZoneEnvelopeRcModel`: 2R1C envelope model (thermal mass in wall assemblies)
  - `solve_dense()`: dense linear system solver

### Visualization (src/draw/)

Uses Rerun (localhost:9876). `RerunConfig` controls session name, entity prefix, colors, sizes, and simulation drawing parameters. Drawing functions accept any `T: HasMesh + HasName`:
- `draw_faces()`: renders triangulated mesh
- `draw_edges()`: draws triangle edges
- `draw_points()`: renders vertices
- `draw_simulation()`: animated ray tracing visualization
- `draw_receivers()`, `draw_impulse_response()` (`draw/acoustics.rs`): acoustic receiver and IR visualization
- `draw_illuminance_heatmap()` (`draw/lighting.rs`): illuminance heatmap on building surfaces
- `draw_heat_loss_heatmap()`, `draw_annual_profile()` (`draw/thermal.rs`): thermal simulation visualization

### File I/O (src/io/)

- **B3D** (`b3d.rs`): Native JSON format, preserves full hierarchy and UIDs
- **STL** (`stl.rs`): Triangulated mesh format (ASCII/Binary), loses hierarchy
- **BIM** (`bim.rs`): dotbim JSON format for BIM interoperability
- **AC3D** (`ac3d.rs`): AC3D `.ac` format reader; materials become Walls, surfaces become Polygons; supports Y-up and Z-up coordinate systems

### Module Layout

```
src/
├── lib.rs              # Public API
├── uid.rs              # UUID wrapper
├── name.rs             # Naming traits
├── vecutils.rs         # Array utilities (min, max, roll, flip)
├── world.rs            # Placeholder for global registry
├── geom/               # Core geometry
│   ├── point.rs        # + point/check.rs, point/convert.rs
│   ├── vector.rs       # + vector/check.rs
│   ├── polygon.rs      # Polygon struct
│   │   ├── containment.rs  # is_point_inside()
│   │   ├── relations.rs    # are_polygons_facing(), touching, coplanar
│   │   ├── slice.rs        # slice_polygon()
│   │   └── boolean.rs      # polygon_intersection(), difference
│   ├── wall.rs
│   ├── solid.rs        # includes FloorPlan
│   │   ├── containment.rs  # is_point_inside() for solids
│   │   └── adjacency.rs    # is_solid_adjacent_to(), get_shared_polygons()
│   ├── zone.rs         # Zone container
│   ├── building/
│   │   ├── mod.rs      # Building struct + path-based access
│   │   └── graph.rs    # get_graph(), stitch_solids()
│   ├── mesh/
│   │   ├── mod.rs      # Mesh struct + HasMesh trait
│   │   ├── quality.rs  # analyze_triangle(), analyze_mesh()
│   │   └── tetrahedralize.rs  # tetrahedralize_centroid()
│   ├── projection.rs   # PlaneBasis for 3D ↔ 2D projection
│   ├── segment.rs      # Line segment operations, intersections
│   ├── triangles.rs    # Ear-clipping triangulation
│   ├── rotation.rs     # Rodrigues rotation, rotate_points_to_plane()
│   ├── ray.rs          # Ray struct, ray-polygon intersection
│   ├── visibility.rs   # are_points_visible(), visibility_matrix()
│   ├── tetrahedron.rs  # Tetrahedron volume/centroid
│   └── bboxes.rs       # Bounding box operations
├── sim/                # Simulation
│   ├── materials.rs    # Material types, MaterialLibrary, presets
│   ├── coupling.rs     # Cross-domain coupling payloads (Bus contracts)
│   ├── surfaces.rs     # SurfaceSemantics, SurfaceKind, SurfaceInterface
│   ├── framework/      # Domain-agnostic simulation runtime
│   │   ├── mod.rs      # Re-exports Bus, SimContext, SimModule, Pipeline
│   │   ├── bus.rs      # Typed message bus
│   │   ├── context.rs  # SimContext (shared building + index)
│   │   ├── module.rs   # SimModule trait
│   │   └── pipeline.rs # Pipeline (module sequencer)
│   ├── index/
│   │   └── mod.rs      # SurfaceIndex, SurfaceRef
│   ├── engine/
│   │   ├── mod.rs      # FlatScene, RayBatch, RayState
│   │   ├── absorption.rs   # AbsorptionModel trait + impls
│   │   ├── reflection.rs   # ReflectionModel trait + impls
│   │   ├── propagation.rs  # PropagationModel trait + impls
│   │   └── voxel_grid.rs   # Spatial hash grid
│   ├── rays/
│   │   ├── mod.rs      # Re-exports AcousticMode, SimulationConfig, Simulation, SimulationResult
│   │   ├── config.rs   # SimulationConfig, AcousticMode
│   │   └── simulation.rs   # Simulation, SimulationResult
│   ├── acoustics/
│   │   ├── source.rs   # SourceDirectivity, Omnidirectional, Cardioid
│   │   ├── receiver.rs # Receiver (spherical energy-time collector)
│   │   ├── impulse_response.rs  # ImpulseResponse extraction
│   │   ├── metrics.rs  # RT, EDT, C50, C80, D50, STI, RoomAcousticReport
│   │   └── auralization/  # Audio output from simulation results
│   │       ├── mod.rs     # write_ir_wav(), auralize(), auralize_per_band()
│   │       ├── convolve.rs    # Fast convolution
│   │       ├── filters.rs     # Octave-band filters (125-4000 Hz)
│   │       └── wav.rs         # WAV file I/O
│   ├── lighting/
│   │   ├── config.rs   # LightingConfig
│   │   ├── simulation.rs   # LightingSimulation (forward ray tracing)
│   │   ├── module.rs   # LightingModule (SimModule wrapper)
│   │   ├── shortwave.rs    # Solar/shortwave coupling modules
│   │   ├── sources.rs  # PointLight, DirectionalLight, LightSource trait
│   │   ├── sensor.rs   # SensorGrid
│   │   ├── sky.rs      # CIEOvercast, CIEClearSky
│   │   ├── solar.rs    # SolarPosition
│   │   ├── backward.rs # BackwardTracer
│   │   └── result.rs   # LightingResult (UID-keyed)
│   └── energy/
│       ├── config.rs   # ThermalConfig
│       ├── simulation.rs   # run_annual_simulation(), AnnualResult
│       ├── module.rs   # EnergyModule, EnergyModuleConfig (SimModule wrapper)
│       ├── boundary.rs # Re-exports SurfaceSemantics as ThermalBoundaries
│       ├── construction.rs # WallConstruction, layer presets
│       ├── weather.rs  # WeatherData, EPW parser
│       ├── weather_module.rs # WeatherModule (publishes weather to Bus)
│       ├── schedule.rs # InternalGainsProfile
│       ├── gains_module.rs  # InternalGainsModule (publishes gains to Bus)
│       ├── hvac.rs     # HvacIdealLoads
│       ├── solar_bridge.rs # Lighting-energy coupling (legacy)
│       ├── zone.rs     # Zone-level heat balance
│       ├── result.rs   # ThermalResult
│       ├── recorder.rs # MultiZoneRecorderData
│       └── network/    # Thermal network graph
│           ├── mod.rs  # ThermalNetwork, InterZoneConductance
│           ├── multizone.rs        # MultiZoneAirModel, MultiZoneStepResult
│           ├── multizone_envelope.rs # MultiZoneEnvelopeRcModel (2R1C)
│           └── solve.rs            # solve_dense() linear solver
├── io/
│   ├── mod.rs          # I/O module exports
│   ├── b3d.rs          # Native JSON format
│   ├── stl.rs          # STL mesh format
│   ├── bim.rs          # dotbim BIM format
│   └── ac3d.rs         # AC3D format reader
└── draw/
    ├── config.rs       # RerunConfig struct
    ├── rerun.rs        # Rerun visualization (faces, edges, points, simulation)
    ├── acoustics.rs    # Receiver and impulse response drawing
    ├── lighting.rs     # Illuminance heatmap
    └── thermal.rs      # Heat loss heatmap, annual profile
```

### Key Algorithms

- **Triangulation**: Ear-clipping algorithm in `triangles.rs`
- **Point-in-polygon**: Ray casting in `polygon/containment.rs`
- **Point-in-solid**: Ray casting counting polygon crossings in `solid/containment.rs`
- **Polygon intersection**: Sutherland-Hodgman clipping in `polygon/boolean.rs`
- **Visibility**: Ray-polygon intersection tests in `visibility.rs`
- **Tetrahedralization**: Centroid-based decomposition in `mesh/tetrahedralize.rs`
- **Acoustic ray tracing**: Time-stepped ray propagation with absorption/reflection in `sim/rays/`
- **Reverberation time**: Schroeder backward integration in `sim/acoustics/metrics.rs`
- **Forward lighting**: RGB ray tracing with bounce accumulation in `sim/lighting/`
- **Energy simulation**: Hourly steady-state heat balance in `sim/energy/`
- **Multi-zone thermal network**: Gaussian elimination solver in `sim/energy/network/solve.rs`
- **Solar shading**: Direct-sun occlusion via FlatScene ray casting in `sim/lighting/shortwave.rs`

### Conventions

- Containers use `HashMap<String, ChildType>` for O(1) access by name
- Methods returning collections (e.g., `walls()`, `polygons()`) return sorted vectors
- Point validation functions in `src/geom/point/check.rs`: coplanarity, collinearity, segment containment
- Tests are inline within modules (`#[cfg(test)] mod tests`)
- Path-based access uses `/` separator: `"zone/solid/wall/polygon"`
- Material assignment uses substring path matching via `MaterialLibrary`
- Simulation modules communicate via `Bus` with UID-keyed payloads in SI units
- Surface semantics (exterior/interior/inter-zone) come from `SurfaceSemantics`, not stored on geometry types

### Path-Based Access

Building provides path-based access to nested entities:
```rust
building.get_zone("zone_name")
building.get_solid("zone/solid")
building.get_wall("zone/solid/wall")
building.get_polygon("zone/solid/wall/polygon")
```
