# building3d-rs

3D building modeling and simulation library in Rust. A rewrite of [building3d](https://github.com/krzysztofarendt/building3d) (Python).

This crate provides a strict hierarchical composition model, geometry algorithms for building-like solids, and simulation capabilities for acoustics, lighting, and energy.

## Hierarchy

`Building → Zone → Solid → Wall → Polygon → Mesh`

Each entity has:
- `name: String` (must not contain `/`, which is used as a path separator)
- `uid: UID`
- `parent: Option<UID>`

## Features

- Core geometry: `Point`, `Vector`, `Mesh`, `PlaneBasis`
- Polygon operations: triangulation (ear-clipping), area/centroid, plane coefficients, containment, relations, slicing, basic boolean ops
- Solids: `from_box()`, `from_floor_plan(FloorPlan)`, volume, point-in-solid (ray casting), adjacency detection
- Building analysis: path-based access (`zone/solid/wall/polygon`), adjacency graph, stitching report
- Simulation:
  - **Materials**: unified `MaterialLibrary` with acoustic, optical, and thermal properties; built-in presets (concrete, glass, gypsum, carpet, wood, metal)
  - **Engine**: shared ray tracing infrastructure with voxel-grid acceleration, pluggable absorption/reflection/propagation models, optional ray history storage (`store_ray_history`)
  - **Acoustic ray tracing**: scalar or frequency-dependent (6 octave bands), configurable source/absorber positions (`building3d::sim::rays`)
  - **Acoustics**: spherical receivers, impulse response extraction, room acoustic metrics (RT20/RT30/EDT, C50/C80, D50, STI), source directivity patterns (`building3d::sim::acoustics`)
  - **Auralization**: WAV export of impulse responses, convolution with dry audio, per-band frequency-dependent auralization (`building3d::sim::acoustics::auralization`)
  - **Lighting**: forward ray tracing with RGB illuminance, point/directional lights, CIE sky models, solar position, sensor grids, backward tracing (`building3d::sim::lighting`)
  - **Energy**: steady-state heat balance, layered wall constructions, EPW weather data, internal gains schedules, HVAC models, annual hourly simulation (`building3d::sim::energy`)
- File I/O:
  - **B3D**: native JSON format preserving hierarchy + UIDs (`building3d::io::b3d`)
  - **STL**: triangulated mesh import/export (`building3d::io::stl`) *(hierarchy is not preserved)*
  - **BIM**: dotbim import/export (`building3d::io::bim`) *(hierarchy is not preserved)*
- Visualization: Rerun-based drawing helpers (`building3d::draw`) with customizable `RerunConfig`, including acoustic, lighting, and thermal heatmap visualizations

## Quick Start

Add as a dependency (path or git, depending on how you consume it), then:

```rust
use building3d::{Building, Solid};

fn main() -> anyhow::Result<()> {
    let solid = Solid::from_box(2.0, 3.0, 4.0, None, "box")?;
    let building = Building::from_solids("building", vec![solid])?;

    println!("Volume: {}", building.volume());
    println!("Zones: {}", building.zones().len());
    Ok(())
}
```

## Examples

Run examples (visualization requires a Rerun viewer running at `localhost:9876`):

```bash
cargo run --example draw_faces
cargo run --example draw_many
cargo run --example draw_shapes
cargo run --example floor_plan
cargo run --example ray_2_boxes
cargo run --example ray_teapot
cargo run --example bench_teapot
cargo run --example sim_acoustics
cargo run --example sim_lighting
cargo run --example sim_energy
```

## I/O

### B3D (native JSON)

```rust
use building3d::{Building, Solid};
use building3d::io::b3d::{read_b3d, write_b3d};
use std::path::Path;

let solid = Solid::from_box(1.0, 1.0, 1.0, None, "box")?;
let building = Building::from_solids("my_building", vec![solid])?;
write_b3d(Path::new("model.b3d"), &building)?;
let loaded = read_b3d(Path::new("model.b3d"))?;
assert!((loaded.volume() - building.volume()).abs() < 1e-10);
```

### STL / BIM

STL and dotbim formats do **not** preserve the full hierarchy. Importers build a simplified structure.

## Development

```bash
cargo fmt
cargo test
cargo clippy --all-targets -- -D warnings
```

## Status

This is an actively evolving library. Expect API changes.
