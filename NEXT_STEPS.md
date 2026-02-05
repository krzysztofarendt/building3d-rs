# Next Steps

Prioritized improvements and next steps for building3d-rs.

All critical bugs (issues #1-#5) and most high-severity issues (#6-#11, #19)
from ISSUES_CLAUDE.md have been fixed. High-priority items (CI, unused deps,
O(1) path access) are also done. The items below are what remains.

## Medium Priority

### 4. Vertex deduplication in mesh hierarchy

`copy_mesh()` duplicates vertices at every level — a cube ends up with 24 vertices instead of 8. For large buildings this is a significant memory waste. A deduplication pass (or shared vertex pool) would help.

### 5. Improve STL import

Currently creates one `Polygon` per triangle with no topology reconstruction. Importing a real-world STL produces degenerate structures. Grouping coplanar adjacent triangles into polygons would make round-tripping useful.

### 6. Documentation

- No module-level docs on `lib.rs`
- Public API items (~226) are moderately documented (~120 doc comments)
- Advanced features (tetrahedralization, visibility, adjacency) lack usage examples
- Running `cargo doc` with `#![warn(missing_docs)]` would highlight gaps

### 7. B3D deserialization validation

After `serde_json::from_reader`, parent UIDs are not checked for existence,
uniqueness, or hierarchy consistency. A corrupted B3D file creates a building
with broken references. Add a `validate()` method called after deserialization.

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

### 13. BIM round-trip fidelity

Export flattens zones into a flat list of solids; import puts everything into
a single zone named `"imported"`. Zone structure is lost on round-trip.

## Summary

| Priority | Area | Status |
|----------|------|--------|
| ~~High~~ | ~~Add CI pipeline~~ | Done |
| ~~High~~ | ~~Remove unused `three-d` dep~~ | Done |
| ~~High~~ | ~~O(1) path-based access~~ | Done |
| Medium | Vertex deduplication | Open |
| Medium | STL topology reconstruction | Open |
| Medium | Documentation pass | Open |
| Medium | B3D validation | Open |
| Low | World module, benchmarks, more examples | Open |
