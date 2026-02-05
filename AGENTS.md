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
│   ├── segment.rs      # Line segment operations, intersections
│   ├── distance.rs     # Point-to-line, point-to-polygon distances
│   ├── triangles.rs    # Ear-clipping triangulation
│   ├── rotation.rs     # Rodrigues rotation, rotate_points_to_plane()
│   ├── ray.rs          # Ray struct, ray-polygon intersection
│   ├── visibility.rs   # are_points_visible(), visibility_matrix()
│   ├── tetrahedron.rs  # Tetrahedron volume/centroid
│   └── bboxes.rs       # Bounding box operations
├── io/
│   ├── mod.rs          # I/O module exports
│   ├── b3d.rs          # Native JSON format
│   ├── stl.rs          # STL mesh format
│   └── bim.rs          # dotbim BIM format
└── draw/
    └── rerun.rs        # Rerun visualization
```

### Key Algorithms

- **Triangulation**: Ear-clipping algorithm in `triangles.rs`
- **Point-in-polygon**: Ray casting in `polygon/containment.rs`
- **Point-in-solid**: Ray casting counting polygon crossings in `solid/containment.rs`
- **Polygon intersection**: Sutherland-Hodgman clipping in `polygon/boolean.rs`
- **Visibility**: Ray-polygon intersection tests in `visibility.rs`
- **Tetrahedralization**: Centroid-based decomposition in `mesh/tetrahedralize.rs`

### Conventions

- Containers use `HashMap<String, ChildType>` for O(1) access by name
- Methods returning collections (e.g., `walls()`, `polygons()`) return sorted vectors
- Point validation functions in `src/geom/point/check.rs`: coplanarity, collinearity, segment containment
- Tests are inline within modules (`#[cfg(test)] mod tests`)
- Path-based access uses `/` separator: `"zone/solid/wall/polygon"`

### Path-Based Access

Building provides path-based access to nested entities:
```rust
building.get_zone("zone_name")
building.get_solid("zone/solid")
building.get_wall("zone/solid/wall")
building.get_polygon("zone/solid/wall/polygon")
```
