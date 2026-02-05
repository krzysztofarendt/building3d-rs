# Next Steps

Prioritized improvements and next steps for building3d-rs.

## Critical: Fix Known Bugs

There are 17 documented issues in `ISSUES_CLAUDE.md`, several of which cause silent incorrect results:

1. **Vector/Vector division** uses `*` instead of `/` (`vector.rs:224-226`)
2. **Bounding box containment** has wrong boolean logic (`&&` should be `||`)
3. **Vector/f64 division** incorrectly rejects negative divisors
4. **`from_box()` doesn't set wall parent UIDs** — breaks hierarchy invariant
5. **`volume()` missing `abs()`** — returns negative for inverted winding

These should be the top priority since they produce wrong results silently.

## High Priority

### 1. Remove panic-prone code in library functions

~20 `unwrap()` calls in non-test code. Functions like `bounding_box()`, `are_points_close()`, and `copy_mesh()` will panic on empty inputs. A library should return `Result` or `Option` instead.

### 2. Add CI/CD

No `.github/workflows/` exists. A basic pipeline running `cargo test`, `cargo clippy`, and `cargo fmt --check` on PRs would catch regressions early.

### 3. Audit/remove unused dependency

`three-d = "0.18"` appears in `Cargo.toml` but doesn't seem to be used in the code. Unused deps increase compile time significantly.

## Medium Priority

### 4. Vertex deduplication in mesh hierarchy

`copy_mesh()` duplicates vertices at every level — a cube ends up with 24 vertices instead of 8. For large buildings this is a significant memory waste. A deduplication pass (or shared vertex pool) would help.

### 5. Improve STL import

Currently creates one `Polygon` per triangle with no topology reconstruction. Importing a real-world STL produces degenerate structures. Grouping coplanar adjacent triangles into polygons would make round-tripping useful.

### 6. Fix Sutherland-Hodgman unprojection

Polygon boolean operations (`polygon/boolean.rs`) assume axis-aligned planes during unprojection. Tilted polygons get wrong z-coordinates. Needs reprojection onto the actual plane equation.

### 7. Documentation

- No module-level docs on `lib.rs`
- Public API items (~226) are moderately documented (~120 doc comments)
- Advanced features (tetrahedralization, visibility, adjacency) lack usage examples
- Running `cargo doc` with `#![warn(missing_docs)]` would highlight gaps

## Lower Priority / Feature Ideas

### 8. Implement the `world` module

Currently a placeholder. A global registry of buildings would enable multi-building scenes and cross-building queries.

### 9. More examples

Current 4 examples all require the Rerun viewer. Adding CLI-friendly examples for I/O round-tripping, containment queries, and adjacency analysis would lower the onboarding barrier.

### 10. Benchmarks

No `benches/` directory. Adding benchmarks for hot paths (triangulation, containment, mesh operations) would help track performance as the library grows.

### 11. Stitch solids properly

`stitch_solids()` currently only detects adjacency — it doesn't actually slice/merge shared polygon edges like the Python version does. This limits its usefulness for thermal simulation workflows.

### 12. Property-based testing

With 194 unit tests passing, the next level would be property-based tests (e.g., with `proptest`) for geometry invariants: "translating then inverse-translating returns the original point", "triangulated mesh area equals polygon area", etc.

## Summary

| Priority | Area | Effort |
|----------|------|--------|
| Critical | Fix 5 silent-wrong-result bugs | Small |
| High | Replace `unwrap()` with `Result` | Medium |
| High | Add CI pipeline | Small |
| High | Remove unused `three-d` dep | Trivial |
| Medium | Vertex deduplication | Medium |
| Medium | STL topology reconstruction | Large |
| Medium | Fix polygon boolean unprojection | Medium |
| Medium | Documentation pass | Medium |
| Low | World module, benchmarks, more examples | Ongoing |
