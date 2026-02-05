# ISSUES.md

Review of the Rust rewrite against REWRITE_PLAN.md and the Python source.

---

## Critical Bugs (ALL FIXED)

### 1. ~~Vector / Vector performs multiplication~~ FIXED (commit 3aaf2e6)

Removed the `Div<Vector> for Vector` impl entirely. Component-wise division
of vectors is not a meaningful geometric operation.

### 2. ~~Bounding box containment check is wrong~~ FIXED (commit fe4054f)

Rewritten with correct `||` logic and empty-input guard.

### 3. ~~Vector / f64 rejects negative divisors~~ FIXED (commit 3aaf2e6)

Removed the EPS guard entirely. `Div<f64> for Vector` now does plain division.

### 4. ~~`from_box()` doesn't set wall parent UIDs~~ FIXED (commit 3aaf2e6)

`from_box()` delegates to `Solid::new()` which sets parent UIDs.

### 5. ~~`volume()` missing `abs()`~~ FIXED (commit fe4054f)

Returns `total_volume.abs()`.

---

## Panic-Prone Library Code (ALL FIXED)

### 6. ~~`copy_mesh()` unwraps `Option<faces>`~~ FIXED (commit fe4054f)

All `copy_mesh()` implementations now use `if let Some(faces)` pattern.

### 7. ~~`bounding_box()` panics on empty or NaN input~~ FIXED (commit fe4054f)

Returns `(zero, zero)` for empty input, uses `total_cmp` for NaN safety.

### 8. ~~`are_points_close()` panics on empty slice~~ FIXED (commit fe4054f)

Guards for empty input with early return `true`.

---

## Algorithmic Issues (MOSTLY FIXED)

### 9. ~~Triangulation retry doesn't reset failure counter~~ NOT A BUG

The recursive call creates a new stack frame where `num_fail` is initialized
to 0. The counter does not actually carry over between attempts.

### 10. ~~Sutherland-Hodgman unprojection assumes axis-aligned planes~~ FIXED

Rewritten to use a proper `PlaneBasis` with orthonormal vectors in the polygon
plane. Projection and unprojection work correctly for any plane orientation.

### 11. ~~Polygon slicing fails when edges are collinear with slice line~~ FIXED

Added deduplication and extreme-point reduction for collinear intersections
(keeps only the two most extreme points along the slice direction).

---

## Open Issues

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

### 19. ~~No name validation~~ FIXED (commit 3aaf2e6)

All constructors and `add_*` methods now call `geom::validate_name()` and
return `Result`. Names containing `/` are rejected.

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

| #  | Severity     | Category        | Status  |
|----|-------------|-----------------|---------|
| 1  | Critical    | Wrong operator  | FIXED   |
| 2  | Critical    | Wrong logic     | FIXED   |
| 3  | Critical    | Missing check   | FIXED   |
| 4  | Critical    | Missing parent  | FIXED   |
| 5  | Critical    | Missing abs()   | FIXED   |
| 6  | High        | Panic           | FIXED   |
| 7  | High        | Panic           | FIXED   |
| 8  | High        | Panic           | FIXED   |
| 9  | High        | Algorithm       | Not a bug |
| 10 | High        | Precision       | FIXED   |
| 11 | High        | Edge case       | FIXED   |
| 12 | High        | Data loss       | Open    |
| 13 | Medium      | Waste/topology  | Open    |
| 14 | Medium      | Missing feature | Open    |
| 15 | Medium      | Missing feature | Open    |
| 16 | Medium      | Missing feature | Open    |
| 17 | Medium      | Missing feature | Open    |
| 18 | Low         | Missing feature | Open    |
| 19 | Low         | Validation      | FIXED   |
| 20 | Low         | Validation      | Open    |
| 21 | Low         | Data loss       | Open    |
| 22 | Low         | Logic           | Open    |
