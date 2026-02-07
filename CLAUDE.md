# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

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

The simulation layer provides ray-based and thermal simulation capabilities, sharing a common engine and material system.

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

- **SimulationConfig**: ray count, speed, source/absorber positions, scalar or frequency-dependent mode, material library integration
- **Simulation**: builds `FlatScene`, runs time-stepped ray tracing with configurable absorption/reflection/propagation
- **SimulationResult**: per-step ray positions, energies, absorber hits (scalar and per-band)

#### Acoustics (`sim/acoustics/`)

- **Receiver**: spherical receiver accumulating energy-time-frequency histograms
- **ImpulseResponse**: extracted from receiver data (per-band and broadband)
- **Metrics**: Schroeder backward integration for reverberation time (RT20/RT30/EDT), clarity (C50/C80), definition (D50)
- **Source directivity**: `Omnidirectional`, `Cardioid` patterns via `SourceDirectivity` trait

#### Lighting (`sim/lighting/`)

- **LightingConfig**: point/directional lights, ray count, bounce limit, material library
- **LightingSimulation**: forward ray tracing with RGB illuminance accumulation
- **Light sources**: `PointLight`, `DirectionalLight` via `LightSource` trait
- **SensorGrid**: auto-generated sensor points on polygon surfaces
- **Sky models**: `CIEOvercast`, `CIEClearSky` via `SkyModel` trait
- **Solar**: `SolarPosition` calculation from latitude/longitude/time
- **Backward ray tracing**: `BackwardTracer` for sensor-to-sky illuminance
- **LightingResult**: per-polygon RGB illuminance map

#### Energy (`sim/energy/`)

- **ThermalConfig**: wall constructions, U-values, indoor/outdoor temperatures, infiltration, gains
- **WallConstruction**: layered wall constructions with R-value/U-value calculation
- **WeatherData**: EPW (EnergyPlus Weather) file parser with hourly records
- **InternalGainsProfile**: hourly/daily schedules for occupancy, equipment, lighting
- **HvacIdealLoads** / **LumpedThermalModel**: simplified HVAC and thermal models
- **SolarBridge**: coupling between lighting and energy simulations
- **run_annual_simulation()**: 8760-hour annual simulation producing heating/cooling demand, peak loads, monthly breakdown
- **ThermalResult** / **AnnualResult**: simulation output types

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
│   ├── distance.rs     # Point-to-line, point-to-polygon distances
│   ├── triangles.rs    # Ear-clipping triangulation
│   ├── rotation.rs     # Rodrigues rotation, rotate_points_to_plane()
│   ├── ray.rs          # Ray struct, ray-polygon intersection
│   ├── visibility.rs   # are_points_visible(), visibility_matrix()
│   ├── tetrahedron.rs  # Tetrahedron volume/centroid
│   └── bboxes.rs       # Bounding box operations
├── sim/                # Simulation
│   ├── materials.rs    # Material types, MaterialLibrary, presets
│   ├── engine/
│   │   ├── mod.rs      # FlatScene, RayBatch, RayState
│   │   ├── absorption.rs   # AbsorptionModel trait + impls
│   │   ├── reflection.rs   # ReflectionModel trait + impls
│   │   ├── propagation.rs  # PropagationModel trait + impls
│   │   ├── find_transparent.rs  # Transparent surface detection
│   │   └── voxel_grid.rs   # Spatial hash grid
│   ├── rays/
│   │   ├── config.rs   # SimulationConfig, AcousticMode
│   │   └── simulation.rs   # Simulation, SimulationResult
│   ├── acoustics/
│   │   ├── source.rs   # SourceDirectivity, Omnidirectional, Cardioid
│   │   ├── receiver.rs # Receiver (spherical energy-time collector)
│   │   ├── impulse_response.rs  # ImpulseResponse extraction
│   │   └── metrics.rs  # RT, EDT, C50, C80, D50
│   ├── lighting/
│   │   ├── config.rs   # LightingConfig
│   │   ├── simulation.rs   # LightingSimulation (forward ray tracing)
│   │   ├── sources.rs  # PointLight, DirectionalLight, LightSource trait
│   │   ├── sensor.rs   # SensorGrid
│   │   ├── sky.rs      # CIEOvercast, CIEClearSky
│   │   ├── solar.rs    # SolarPosition
│   │   ├── backward.rs # BackwardTracer
│   │   └── result.rs   # LightingResult
│   └── energy/
│       ├── config.rs   # ThermalConfig
│       ├── simulation.rs   # run_annual_simulation(), AnnualResult
│       ├── construction.rs # WallConstruction, layer presets
│       ├── weather.rs  # WeatherData, EPW parser
│       ├── schedule.rs # InternalGainsProfile
│       ├── hvac.rs     # HvacIdealLoads, LumpedThermalModel
│       ├── solar_bridge.rs # Lighting-energy coupling
│       ├── zone.rs     # Zone-level heat balance
│       └── result.rs   # ThermalResult
├── io/
│   ├── mod.rs          # I/O module exports
│   ├── b3d.rs          # Native JSON format
│   ├── stl.rs          # STL mesh format
│   └── bim.rs          # dotbim BIM format
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

### Conventions

- Containers use `HashMap<String, ChildType>` for O(1) access by name
- Methods returning collections (e.g., `walls()`, `polygons()`) return sorted vectors
- Point validation functions in `src/geom/point/check.rs`: coplanarity, collinearity, segment containment
- Tests are inline within modules (`#[cfg(test)] mod tests`)
- Path-based access uses `/` separator: `"zone/solid/wall/polygon"`
- Material assignment uses substring path matching via `MaterialLibrary`

### Path-Based Access

Building provides path-based access to nested entities:
```rust
building.get_zone("zone_name")
building.get_solid("zone/solid")
building.get_wall("zone/solid/wall")
building.get_polygon("zone/solid/wall/polygon")
```
