# ISSUES.md

Review of the Rust rewrite against REWRITE_PLAN.md and the Python source.

---

## Critical Bugs

### 1. Vector / Vector performs multiplication

**File:** `src/geom/vector.rs:224-226`

`Div for Vector` uses `*` instead of `/`:

```rust
dx: self.dx * other.dx,  // should be /
dy: self.dy * other.dy,
dz: self.dz * other.dz,
```

Any code path using `Vector / Vector` silently returns wrong results.
The `f64 / Vector` impl (line 252) is correct.

**Fix:** Change `*` to `/` on lines 224-226.

---

### 2. Bounding box containment check is wrong

**File:** `src/geom/bboxes.rs:7-8`

`is_point_inside_bbox` uses `&&` within each branch:

```rust
if (ptest.x < pmin.x && ptest.y < pmin.y && ptest.z < pmin.z)
    || (ptest.x > pmax.x && ptest.y > pmax.y && ptest.z > pmax.z)
```

A point at `(pmin.x - 1, pmax.y, pmax.z)` is outside the box but passes
this check. Each coordinate must be checked independently:

```rust
if ptest.x < pmin.x || ptest.x > pmax.x
    || ptest.y < pmin.y || ptest.y > pmax.y
    || ptest.z < pmin.z || ptest.z > pmax.z
```

This affects every containment test that uses bbox pre-filtering.

**Fix:** Replace `&&` with `||` so any single out-of-range coordinate rejects.

---

### 3. Vector / f64 rejects negative divisors

**File:** `src/geom/vector.rs:234`

```rust
if other < EPS {  // negative values pass, but so does -1e-14
```

Should be `other.abs() < EPS`. Currently `Vector / (-0.5)` returns `None`.

**Fix:** Use `other.abs() < EPS`.

---

### 4. `from_box()` doesn't set wall parent UIDs

**File:** `src/geom/solid.rs:237-251`

`Solid::new()` (line 62-63) sets `w.parent = Some(uid.clone())` for each wall.
`Solid::from_box()` builds the HashMap directly and never sets parent UIDs.
All walls created via `from_box()` have `parent = None`, breaking the hierarchy.

**Fix:** Loop over walls and set parent before inserting into HashMap,
or delegate to `Solid::new()`.

---

### 5. `volume()` missing `abs()`

**File:** `src/geom/solid.rs:142`

```rust
total_volume  // Python returns abs(total_volume)
```

For solids with inverted winding order, this returns a negative volume.

**Fix:** Return `total_volume.abs()`.

---

## Panic-Prone Library Code

### 6. `copy_mesh()` unwraps `Option<faces>` at every hierarchy level

**Files:**
- `src/geom/wall.rs:36`
- `src/geom/solid.rs:43, 124`
- `src/geom/zone.rs:44`
- `src/geom/building/mod.rs:46`

All do `.faces.unwrap()` on `Mesh.faces: Option<Vec<TriangleIndex>>`.
If any polygon has `faces = None`, the library panics.

**Fix:** Handle `None` with `unwrap_or_default()` or propagate as `Result`.

---

### 7. `bounding_box()` panics on empty or NaN input

**File:** `src/geom/bboxes.rs:50-55`

```rust
let xmin = *x.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
```

Two panic sources: `partial_cmp()` returns `None` for NaN, outer `unwrap()`
panics on empty iterator. Same issue in `vecutils::min/max` (lines 2-8).

**Fix:** Return `Option<(Point, Point)>` or validate input.

---

### 8. `are_points_close()` panics on empty slice

**File:** `src/geom/point/check.rs:78`

```rust
let p0 = pts[0];  // panics if pts is empty
```

Same pattern in `are_vectors_close()` (`src/geom/vector/check.rs:6`).

**Fix:** Return `true` for empty (vacuously true) or guard with `is_empty()`.

---

## Algorithmic Issues

### 9. Triangulation retry doesn't reset failure counter

**File:** `src/geom/triangles.rs:41`

When ear-clipping fails and retries with flipped winding, `num_fail` carries
over from the previous attempt. The existing TODO comment confirms this.
Python resets it. This causes premature failure on valid non-convex polygons.

**Fix:** Reset `num_fail` to 0 on retry, or pass it as 0 in the recursive call.

---

### 10. Sutherland-Hodgman unprojection assumes axis-aligned planes

**File:** `src/geom/polygon/boolean.rs:168-187`

When projecting 3D polygons to 2D for clipping, the unprojection uses a
constant coordinate from `subject[0]` (e.g., `z_val = subject[0].z`).
For tilted polygons, intersection points will be placed at the wrong
z-coordinate instead of on the actual plane.

**Fix:** After unprojection, reproject each result point onto the plane
equation of the original polygon.

---

### 11. Polygon slicing fails when edges are collinear with slice line

**File:** `src/geom/polygon/slice.rs:84-98`

When polygon edges are collinear with the slice line, both edge endpoints
are added as intersections. Two consecutive collinear edges produce 3-4
intersections, causing the `intersections.len() != 2` guard (line 104)
to reject the slice.

**Fix:** After collecting intersections, deduplicate and keep only the
two extreme points along the slice direction.

---

### 12. STL reader duplicates all vertices

**File:** `src/io/stl.rs:232-234` (ASCII), `285-286` (binary)

Every triangle gets 3 fresh vertices. A cube mesh produces 36 vertices
instead of 8. Adjacent triangles don't share vertices, so topology is lost.
`stl_to_solid()` then creates one Polygon per triangle, producing a Solid
with no meaningful wall structure.

**Fix:** Deduplicate vertices by position (within epsilon) during read,
and group connected triangles into walls during `stl_to_solid()`.

---

### 13. `copy_mesh()` never deduplicates vertices

**Files:** Wall, Solid, Zone, Building `copy_mesh()` implementations

Each polygon's vertices are concatenated with offset indices. A 6-face cube
Solid produces 24 mesh vertices instead of 8. This wastes memory and breaks
algorithms that expect shared topology.

Additionally, `poly.copy_mesh()` is called **twice** per polygon in each
implementation (once for vertices, once for faces).

**Fix:** Store the mesh from one `copy_mesh()` call and reuse. Optionally
deduplicate vertices across polygons.

---

## Missing Features vs Python

### 14. `stitch_solids()` only detects adjacency, doesn't modify geometry

**Rust:** `src/geom/building/graph.rs` returns `StitchInfo` structs.
**Python:** Recursively slices adjacent polygons so they share vertices/edges,
actually modifying the building geometry. The Rust version is informational
only.

---

### 15. Missing polygon methods

- `contains_polygon(other, margin)` - checks if one polygon is entirely
  inside another. Required for proper stitching.
- `get_some_interior_point()` - returns a point guaranteed inside the polygon.
- `is_point_inside_projection()` / `is_point_inside_ortho_projection()` -
  projects a 3D point onto the polygon's plane before containment check.

---

### 16. Simplified polygon slicing

**Python:** Full pipeline with point classification (VERTEX/EDGE/INTERIOR),
redundant point removal, and edge intersection detection across multiple
helper modules.

**Rust:** Basic two-point slice only. Missing classification, edge cases,
and multi-edge handling.

---

### 17. Missing `distance_point_to_edge()`

Plan Stage 3.3 specifies this (distance to line segment, clamped to
endpoints). Python has it. Rust has `distance_point_to_line()` (infinite
line) but not the segment-clamped variant.

---

### 18. Path access stops at polygon level

**Python:** Supports 6 levels including individual point access via index.
**Rust:** `src/geom/building/mod.rs:145-180` only goes 4 levels deep
(zone/solid/wall/polygon).

---

## Design Issues

### 19. No name validation

CLAUDE.md states `/` is forbidden in entity names (path separator).
No constructor (`Polygon::new`, `Wall::new`, `Solid::new`, `Zone::new`,
`Building::new`) validates this. Users can silently break path-based access.

**Fix:** Validate in constructors and return `Result` (or at least assert).

---

### 20. B3D deserialization doesn't validate hierarchy integrity

**File:** `src/io/b3d.rs:55-63`

After `serde_json::from_reader`, parent UIDs are not checked for:
- Existence of referenced parents
- Uniqueness of UIDs
- Consistency of the hierarchy

A corrupted B3D file creates a building with broken references.

**Fix:** Add a `validate()` method and call it after deserialization.

---

### 21. BIM round-trip destroys zone structure

**File:** `src/io/bim.rs:255-272, 317-323`

Export flattens all zones into a flat list of solids. Import puts everything
into a single zone named `"imported"`. Original zone names and grouping
are lost.

---

### 22. `stitch_solids` applies one interface check to all shared polygons

**File:** `src/geom/solid/adjacency.rs:454`

`has_correct_interface()` returns one boolean for the entire solid pair.
The loop applies that same result to every shared polygon between those
solids, even if only one interface is incorrect.

**Fix:** Check interface correctness per polygon pair, not per solid pair.

---

## Summary

| #  | Severity     | Category        | Location                      |
|----|-------------|-----------------|-------------------------------|
| 1  | Critical    | Wrong operator  | `vector.rs:224`               |
| 2  | Critical    | Wrong logic     | `bboxes.rs:7`                 |
| 3  | Critical    | Missing check   | `vector.rs:234`               |
| 4  | Critical    | Missing parent  | `solid.rs:237`                |
| 5  | Critical    | Missing abs()   | `solid.rs:142`                |
| 6  | High        | Panic           | wall/solid/zone/building      |
| 7  | High        | Panic           | `bboxes.rs:50`                |
| 8  | High        | Panic           | `point/check.rs:78`           |
| 9  | High        | Algorithm       | `triangles.rs:41`             |
| 10 | High        | Precision       | `polygon/boolean.rs:168`      |
| 11 | High        | Edge case       | `polygon/slice.rs:84`         |
| 12 | High        | Data loss       | `stl.rs:232`                  |
| 13 | Medium      | Waste/topology  | `copy_mesh()` everywhere      |
| 14 | Medium      | Missing feature | `graph.rs` stitch             |
| 15 | Medium      | Missing feature | polygon methods               |
| 16 | Medium      | Missing feature | polygon slicing               |
| 17 | Medium      | Missing feature | distance to edge              |
| 18 | Low         | Missing feature | path access depth             |
| 19 | Low         | Validation      | name validation               |
| 20 | Low         | Validation      | B3D deserialization           |
| 21 | Low         | Data loss       | BIM round-trip                |
| 22 | Low         | Logic           | stitch interface check        |
