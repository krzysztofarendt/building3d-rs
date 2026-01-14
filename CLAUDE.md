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
cargo run --example draw_shapes
cargo run --example floor_plan

# Check code without building
cargo check

# Format code
cargo fmt

# Lint
cargo clippy
```

## Architecture

This is a 3D building modeling library with a strict hierarchical composition model:

```
Building → Solid → Wall → Polygon → Mesh
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
- **Building**: top-level container of solids

### Key Traits

- **HasMesh**: provides `copy_mesh()` for polymorphic mesh access (used by drawing)
- **HasName**: polymorphic name access
- **SortByName**: in-place stable sort by name
- **IsClose**: epsilon-based f64/Point/Vector comparison (uses `EPS = 1e-13`)

### Transformations

`rotate()` and `translate()` methods cascade down the hierarchy. Rotation uses Rodrigues' formula (src/geom/rotation.rs).

### Visualization (src/draw/)

Uses Rerun (localhost:9876). Drawing functions accept any `T: HasMesh + HasName`:
- `draw_faces()`: renders triangulated mesh
- `draw_edges()`: draws triangle edges
- `draw_points()`: renders vertices

### Module Layout

```
src/
├── lib.rs              # Public API
├── uid.rs              # UUID wrapper
├── name.rs             # Naming traits
├── vecutils.rs         # Array utilities (min, max, roll, flip)
├── geom/               # Core geometry
│   ├── point.rs        # + point/check.rs, point/convert.rs
│   ├── vector.rs       # + vector/check.rs
│   ├── polygon.rs
│   ├── wall.rs
│   ├── solid.rs        # includes FloorPlan
│   ├── building.rs
│   ├── mesh.rs
│   ├── triangles.rs    # Ear-clipping triangulation
│   ├── rotation.rs     # Rodrigues rotation
│   └── bboxes.rs
└── draw/
    └── rerun.rs
```

### Conventions

- Containers use `HashMap<String, ChildType>` for O(1) access by name
- Methods returning collections (e.g., `walls()`, `polygons()`) return sorted vectors
- Point validation functions in `src/geom/point/check.rs`: coplanarity, collinearity, segment containment
- Tests are inline within modules (`#[cfg(test)] mod tests`)
