# Building3D Rust Rewrite Plan

This document outlines the staged implementation plan to complete the Rust rewrite of building3d.

## Current State

The Rust codebase has the core geometry hierarchy implemented:
- Point, Vector, Mesh, Polygon, Wall, Solid, Building
- Basic transformations (rotate, translate)
- Triangulation (ear-clipping)
- Visualization (Rerun)

## Missing Features (from Python)

### Critical Missing
1. Zone level in hierarchy
2. Point-in-polygon test
3. Point-in-solid test (ray casting)
4. Polygon area/centroid (public API)
5. File I/O (B3D, STL)

### Medium Priority
6. Polygon relationships (facing, crossing, touching)
7. Line segment operations
8. Distance calculations
9. Solid adjacency
10. Building graph analysis

### Lower Priority
11. Polygon slicing
12. Visibility/ray tracing
13. Path-based access (World)
14. .BIM format I/O

---

## Stage 1: Polygon Properties and Point-in-Polygon

**Goal**: Expose polygon geometric properties and implement point containment test.

### Tasks

1.1. **Add `area()` method to Polygon**
   - Calculate using cross product of edges
   - Formula: `0.5 * |sum(cross(v[i], v[i+1]))|`
   - File: `src/geom/polygon.rs`

1.2. **Add `centroid()` method to Polygon**
   - Weighted average of triangle centroids by area
   - File: `src/geom/polygon.rs`

1.3. **Add `edges()` method to Polygon**
   - Return `Vec<(Point, Point)>` of consecutive vertex pairs
   - File: `src/geom/polygon.rs`

1.4. **Add `plane_coefficients()` method to Polygon**
   - Return `(a, b, c, d)` for plane equation `ax + by + cz + d = 0`
   - Normal vector gives (a, b, c), solve for d using any point
   - File: `src/geom/polygon.rs`

1.5. **Implement `is_point_inside()` for Polygon**
   - Ray casting or point-in-triangle accumulation
   - Handle boundary cases with `boundary_in: bool` parameter
   - File: `src/geom/polygon/containment.rs` (new)

### Tests
```rust
#[test] fn test_polygon_area_square()
#[test] fn test_polygon_area_triangle()
#[test] fn test_polygon_centroid()
#[test] fn test_polygon_edges()
#[test] fn test_plane_coefficients()
#[test] fn test_point_inside_polygon()
#[test] fn test_point_outside_polygon()
#[test] fn test_point_on_polygon_boundary()
```

---

## Stage 2: Line Segment Operations

**Goal**: Implement line segment intersection and related operations.

### Tasks

2.1. **Create `src/geom/segment.rs` module**
   - Define `Segment` struct (two Points) or use tuple `(Point, Point)`

2.2. **Implement `line_segment_intersection()`**
   - 3D line-line intersection (closest points if skew)
   - Return `Option<Point>` or `Option<(Point, Point)>` for skew lines
   - File: `src/geom/segment.rs`

2.3. **Implement `are_segments_parallel()`**
   - Check if direction vectors are parallel (cross product near zero)
   - File: `src/geom/segment.rs`

2.4. **Implement `is_line_segment_crossing_polygon()`**
   - Check plane intersection, then point-in-polygon test
   - File: `src/geom/segment.rs` or `src/geom/polygon/intersect.rs`

### Tests
```rust
#[test] fn test_segment_intersection_2d()
#[test] fn test_segment_intersection_3d()
#[test] fn test_parallel_segments()
#[test] fn test_skew_segments()
#[test] fn test_segment_crosses_polygon()
#[test] fn test_segment_misses_polygon()
```

---

## Stage 3: Distance Calculations

**Goal**: Implement distance functions for geometric queries.

### Tasks

3.1. **Implement `distance_point_to_point()`**
   - Simple Euclidean distance (may already exist as Vector::length)
   - File: `src/geom/point.rs`

3.2. **Implement `distance_point_to_line()`**
   - Perpendicular distance from point to infinite line
   - File: `src/geom/distance.rs` (new)

3.3. **Implement `distance_point_to_edge()`**
   - Distance to line segment (clamped to endpoints)
   - File: `src/geom/distance.rs`

3.4. **Implement `closest_point_on_line()`**
   - Project point onto line, return closest point
   - File: `src/geom/distance.rs`

3.5. **Implement `distance_point_to_polygon()`**
   - If inside: 0, else: min distance to edges
   - File: `src/geom/distance.rs`

### Tests
```rust
#[test] fn test_point_to_point_distance()
#[test] fn test_point_to_line_distance()
#[test] fn test_point_to_edge_distance()
#[test] fn test_closest_point_on_line()
#[test] fn test_distance_to_polygon_inside()
#[test] fn test_distance_to_polygon_outside()
```

---

## Stage 4: Polygon Relationships

**Goal**: Determine spatial relationships between polygons.

### Tasks

4.1. **Implement `are_polygons_facing()`**
   - Same plane, opposite normals, overlapping areas
   - File: `src/geom/polygon/relations.rs` (new)

4.2. **Implement `are_polygons_coplanar()`**
   - Check if all points satisfy same plane equation
   - File: `src/geom/polygon/relations.rs`

4.3. **Implement `are_polygons_crossing()`**
   - Check for edge intersections (excluding exact overlap)
   - File: `src/geom/polygon/relations.rs`

4.4. **Implement `are_polygons_touching()`**
   - Coplanar + share edge segments without crossing
   - File: `src/geom/polygon/relations.rs`

4.5. **Implement `polygon_overlap_area()`**
   - Calculate intersection area of two coplanar polygons
   - File: `src/geom/polygon/relations.rs`

### Tests
```rust
#[test] fn test_polygons_facing()
#[test] fn test_polygons_not_facing()
#[test] fn test_polygons_coplanar()
#[test] fn test_polygons_crossing()
#[test] fn test_polygons_touching()
#[test] fn test_polygons_separate()
```

---

## Stage 5: Point-in-Solid (Ray Casting)

**Goal**: Implement 3D point containment test for Solid.

### Tasks

5.1. **Implement `is_point_inside()` for Solid**
   - Ray casting algorithm: cast ray from point, count polygon crossings
   - Odd count = inside
   - Handle edge cases (ray through edge/vertex)
   - File: `src/geom/solid/containment.rs` (new)

5.2. **Implement `is_point_at_boundary()` for Solid**
   - Check if point lies on any polygon surface
   - File: `src/geom/solid/containment.rs`

5.3. **Implement `are_bboxes_overlapping()`**
   - Quick rejection test before expensive containment checks
   - File: `src/geom/bboxes.rs`

### Tests
```rust
#[test] fn test_point_inside_box_solid()
#[test] fn test_point_outside_solid()
#[test] fn test_point_on_solid_face()
#[test] fn test_point_on_solid_edge()
#[test] fn test_bboxes_overlapping()
#[test] fn test_bboxes_separate()
```

---

## Stage 6: Zone Container

**Goal**: Add Zone level to hierarchy (Building → Zone → Solid).

### Tasks

6.1. **Create `src/geom/zone.rs`**
   - Struct with `name`, `uid`, `parent`, `solids: HashMap<String, Solid>`
   - Implement `HasName`, `HasMesh`

6.2. **Implement Zone methods**
   - `new()`, `add_solid()`, `solids()`, `walls()`, `polygons()`
   - `rotate()`, `translate()` (cascade to children)
   - `volume()` (sum of solid volumes)
   - `bbox()`

6.3. **Update Building to use Zone**
   - Change `solids: HashMap<String, Solid>` to `zones: HashMap<String, Zone>`
   - Add methods: `zones()`, `add_zone()`
   - Update `walls()`, `polygons()` to traverse through zones
   - Maintain backward compatibility or provide migration path

6.4. **Update visualization**
   - Ensure `HasMesh` works correctly with Zone

### Tests
```rust
#[test] fn test_zone_creation()
#[test] fn test_zone_volume()
#[test] fn test_building_with_zones()
#[test] fn test_zone_transform_cascade()
#[test] fn test_building_volume_with_zones()
```

---

## Stage 7: Solid and Zone Adjacency

**Goal**: Detect when solids share faces.

### Tasks

7.1. **Implement `is_adjacent_to_solid()` for Solid**
   - Check if any polygon faces another solid's polygon
   - Use `are_polygons_facing()` from Stage 4
   - File: `src/geom/solid/adjacency.rs` (new)

7.2. **Implement `has_correct_interface()` for Solid**
   - Validate that interface polygons match exactly
   - File: `src/geom/solid/adjacency.rs`

7.3. **Implement `get_shared_polygons()` for Solid**
   - Return list of polygon pairs that form interfaces
   - File: `src/geom/solid/adjacency.rs`

### Tests
```rust
#[test] fn test_adjacent_solids()
#[test] fn test_non_adjacent_solids()
#[test] fn test_correct_interface()
#[test] fn test_incorrect_interface()
#[test] fn test_shared_polygons()
```

---

## Stage 8: Building Graph

**Goal**: Build and query adjacency graph for building analysis.

### Tasks

8.1. **Implement `get_graph()` for Building**
   - Parameters: `level` (polygon/wall/solid/zone), `facing`, `overlapping`, `touching`
   - Return `HashMap<String, Vec<String>>` mapping paths to adjacent paths
   - File: `src/geom/building/graph.rs` (new)

8.2. **Implement path-based access**
   - `get(path: &str)` method to retrieve any object by path
   - Path format: "zone/solid/wall/polygon"
   - File: `src/geom/building.rs` or `src/geom/paths.rs`

8.3. **Implement `stitch_solids()`**
   - Find and record adjacent solid interfaces
   - File: `src/geom/building/graph.rs`

### Tests
```rust
#[test] fn test_building_graph_polygon_level()
#[test] fn test_building_graph_solid_level()
#[test] fn test_path_access()
#[test] fn test_stitch_solids()
```

---

## Stage 9: File I/O - B3D Format

**Goal**: Implement native JSON serialization.

### Tasks

9.1. **Add serde dependencies**
   - `serde`, `serde_json` in Cargo.toml

9.2. **Derive Serialize/Deserialize for core types**
   - Point, Vector, UID
   - Polygon, Wall, Solid, Zone, Building
   - Handle Mesh specially (vertices + faces)

9.3. **Create `src/io/mod.rs` and `src/io/b3d.rs`**
   - `write_b3d(path: &Path, building: &Building) -> Result<()>`
   - `read_b3d(path: &Path) -> Result<Building>`

9.4. **Maintain hierarchy in JSON**
   - Nested structure matching Python implementation
   - Store UIDs for reference integrity

### Tests
```rust
#[test] fn test_write_b3d()
#[test] fn test_read_b3d()
#[test] fn test_b3d_roundtrip()
#[test] fn test_b3d_preserves_uids()
```

---

## Stage 10: File I/O - STL Format

**Goal**: Import/export STL mesh files.

### Tasks

10.1. **Create `src/io/stl.rs`**
   - Use existing crate like `stl_io` or implement manually

10.2. **Implement `write_stl()`**
   - Export all triangles with normals
   - ASCII or binary format option

10.3. **Implement `read_stl()`**
   - Parse triangles, group into walls/solids
   - Note: STL loses hierarchy information

### Tests
```rust
#[test] fn test_write_stl_ascii()
#[test] fn test_write_stl_binary()
#[test] fn test_read_stl()
#[test] fn test_stl_roundtrip()
```

---

## Stage 11: Polygon Slicing

**Goal**: Cut polygons with a line.

### Tasks

11.1. **Implement `slice_polygon()`**
   - Input: polygon, slicing points (line on polygon plane)
   - Output: two new polygons
   - Handle edge cases (slice through vertex)
   - File: `src/geom/polygon/slice.rs` (new)

### Tests
```rust
#[test] fn test_slice_square_horizontally()
#[test] fn test_slice_square_diagonally()
#[test] fn test_slice_complex_polygon()
#[test] fn test_slice_through_vertex()
```

---

## Stage 12: Visibility and Ray Tracing

**Goal**: Implement visibility analysis.

### Tasks

12.1. **Implement `are_points_visible()`**
   - Check if line segment between points crosses any blocking geometry
   - File: `src/geom/visibility.rs` (new)

12.2. **Implement `visibility_matrix()`**
   - N×N matrix of pairwise visibility
   - File: `src/geom/visibility.rs`

12.3. **Implement basic ray casting infrastructure**
   - Ray struct, intersection tests
   - File: `src/geom/ray.rs` (new)

### Tests
```rust
#[test] fn test_points_visible()
#[test] fn test_points_blocked()
#[test] fn test_visibility_matrix()
#[test] fn test_ray_polygon_intersection()
```

---

## Stage 13: rotate_points_to_plane and Cleanup

**Goal**: Complete unimplemented functions and cleanup.

### Tasks

13.1. **Implement `rotate_points_to_plane()`**
   - Rotate points so they lie on a target plane
   - File: `src/geom/rotation.rs`

13.2. **Implement World struct (if needed)**
   - Global registry for objects by UID
   - Path-based lookup
   - File: `src/world.rs`

13.3. **Code cleanup**
   - Remove TODO comments for completed features
   - Add documentation
   - Ensure all public API is consistent

### Tests
```rust
#[test] fn test_rotate_points_to_xy_plane()
#[test] fn test_rotate_points_to_arbitrary_plane()
#[test] fn test_world_registry()
```

---

## Stage 14 (Future): Advanced Features

These features from Python's `_to_be_ported` directory are deferred:

- **Tetrahedralization** (3D volumetric meshing)
- **Advanced triangulation** (Delaunay with constraints)
- **Mesh quality analysis**
- **Boolean operations** (polygon holes, imprinting)
- **.BIM format** I/O

---

## Summary Timeline

| Stage | Focus | Estimated Complexity |
|-------|-------|---------------------|
| 1 | Polygon properties + point-in-polygon | Medium |
| 2 | Line segment operations | Medium |
| 3 | Distance calculations | Easy |
| 4 | Polygon relationships | Medium |
| 5 | Point-in-solid (ray casting) | Medium |
| 6 | Zone container | Easy |
| 7 | Solid adjacency | Medium |
| 8 | Building graph | Medium |
| 9 | B3D file I/O | Easy |
| 10 | STL file I/O | Easy |
| 11 | Polygon slicing | Hard |
| 12 | Visibility/ray tracing | Medium |
| 13 | Cleanup and finish | Easy |

---

## Running Tests

After each stage:
```bash
cargo test
cargo clippy
cargo fmt --check
```

For specific stage tests:
```bash
cargo test polygon_area
cargo test point_inside
# etc.
```
