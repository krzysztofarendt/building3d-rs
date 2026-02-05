# Review: Rust rewrite of `building3d` (branch `dev`)

This document records issues found during a code review of the current branch implementation against `REWRITE_PLAN.md`.

Status:
- `cargo test`: PASS (194 unit tests + doc-tests)
- `cargo fmt --check`: FAIL (formatting diffs in `src/geom/polygon/boolean.rs`)
- `cargo clippy --all-targets -- -D warnings`: FAIL (several clippy lints; see “Tooling”)

---

## 1) High-Risk Correctness & API-Safety Issues (Fix First)

### 1.1 Name invariants are not enforced consistently

The repo’s own architecture says `/` is a path separator and must not appear in any `name`. In practice, only `Polygon::new` validates names (`geom::validate_name`), while `Building/Zone/Solid/Wall` do not validate, allowing invalid names that break path-based access and can serialize into invalid B3D.

Where:
- `src/geom.rs`: `validate_name` exists but is private to `geom`.
- `src/geom/polygon.rs`: `Polygon::new` calls `geom::validate_name(name)?;`
- `src/geom/building/mod.rs`, `src/geom/zone.rs`, `src/geom/solid.rs`, `src/geom/wall.rs`: constructors accept arbitrary `name`.

Fix suggestions:
- Apply `validate_name` uniformly:
  - `Building::new`, `Zone::new`, `Solid::new`, `Wall::new`, plus `add_zone/add_solid/add_wall/add_polygon`.
- Consider moving name validation into a shared helper (e.g. `crate::name` or public `geom::validate_name`) and/or a `Name(String)` newtype that enforces invariants at construction.

### 1.2 Public APIs can panic due to unvalidated normals / rotation axes

There are `panic!`/`assert!` paths reachable through public APIs because normal vectors and rotation axes are not normalized/validated consistently.

Where:
- `src/geom/triangles.rs`: `is_corner_convex` asserts the polygon normal length is ~1.
  - `triangulate` uses caller-provided `vn` (from `Polygon::new`) without normalizing it.
- `src/geom/polygon.rs`: `Polygon::new` explicitly says “If [normal] is provided, its validity isn’t checked.”
- `src/geom/rotation.rs`: `rotation_matrix` panics if axis is not unit; `rotate_points_around_vector` calls it without normalizing.

Impact:
- A user can pass `Some(Vector::new(0.0, 0.0, 2.0))` as a polygon normal and trigger panics during triangulation.
- A user can call `rotate(angle, &Vector::new(0.0, 0.0, 2.0))` and panic.

Fix suggestions:
- Make panics unreachable from library public methods:
  - Normalize normals and axes inside implementation (preferred), and/or
  - Change `rotation_matrix` to return `Result<Array2<f64>>` and propagate errors.
- In `Polygon::new`, if `normal: Some(v)` is passed:
  - validate non-zero and normalize, or reject with `Err`.

### 1.3 `Solid::from_floor_plan` uses `assert_eq!` on user input

`Solid::from_floor_plan(fp: FloorPlan) -> Result<Self>` has input checks that panic rather than return `Err`.

Where:
- `src/geom/solid.rs`: `assert_eq!` checks on `wall_names` length and uniqueness.

Fix suggestion:
- Replace `assert_eq!` with `anyhow::bail!` (or `Err(anyhow!(...))`) with context.

### 1.4 Point-in-solid ray casting uses a ray length derived only from `bbox_max`

`cast_ray` uses `bbox_max` (and the point under test) to compute the “diagonal”, which can be near zero when the point is near `bbox_max` and the ray points toward `bbox_min`. That can produce a ray that does not reliably leave the solid, causing incorrect parity counts.

Where:
- `src/geom/solid/containment.rs`: `cast_ray(..., bbox_max: &Point)` and `diagonal` calculation uses `bbox_max - ptest`.

Fix suggestions:
- Compute ray length based on full bounding box diagonal (`bbox_max - bbox_min`), or pick an endpoint guaranteed outside by using both corners.
- Consider using a single well-chosen direction + robust tie-breaking (or jitter) rather than a fixed set plus majority vote.

---

## 2) Robustness / Numerical Stability

### 2.1 Absolute epsilon (`EPS = 1e-13`) is too rigid for scale-varying models

Many geometric predicates use a single absolute epsilon, which is fragile for:
- large coordinate magnitudes (meters vs kilometers, etc.),
- very small features.

Where:
- `src/geom.rs`: `const EPS: f64 = 1e-13;`
- Plane checks, segment-plane classification, containment, polygon relations, etc.

Fix suggestions:
- Introduce a relative tolerance scheme (e.g. `abs(a-b) <= abs_eps + rel_eps * max(|a|,|b|)`).
- For plane distances, scale epsilon by normal magnitude and/or bbox size.

### 2.2 Deserialization does not restore invariants / may enable later panics

The library relies on invariants like “faces are present” in many code paths, yet `read_b3d` can deserialize objects with `Mesh.faces = None` or invalid geometry. Subsequent calls may `unwrap()` and panic.

Where:
- `src/io/b3d.rs`: `serde_json::from_reader` directly into `Building`.
- Several `copy_mesh` implementations `unwrap()` on `faces`.

Fix suggestions:
- Add `validate()` methods for `Building/Zone/Solid/Wall/Polygon/Mesh` and call them after load.
- Consider custom `Deserialize` implementations or “builder” types that validate on construction.

---

## 3) API Design & Performance

### 3.1 Constructors silently overwrite duplicate child names

Constructors convert `Vec<T>` into `HashMap<String, T>` using `collect()`. If there are duplicates, earlier entries are overwritten silently.

Where:
- `Building::new` (`zones`), `Zone::new` (`solids`), `Solid::new` (`walls`), `Wall::new` (`polygons`).

Fix suggestions:
- Change constructors to `Result<Self>` and explicitly error on duplicate names.
- Alternatively, accept `HashMap` as input and provide `try_from(Vec<...>)`.

### 3.2 Path-based access methods are avoidably O(n) and allocate

`Building::get_solid/get_wall/get_polygon` split the path then use `zone.solids()` which allocates and sorts `Vec<&Solid>`, then linear search by name.

Where:
- `src/geom/building/mod.rs`: `get_solid/get_wall/get_polygon`.

Fix suggestions:
- Add direct getters that hit internal `HashMap`s:
  - `Zone::get_solid(&self, name: &str) -> Option<&Solid>`
  - `Solid::get_wall(&self, name: &str) -> Option<&Wall>`
  - `Wall::get_polygon(&self, name: &str) -> Option<&Polygon>`
- Then implement path getters using direct lookup without sorting/allocations.

### 3.3 `copy_mesh()` clones aggressively

Many operations build a fully cloned combined mesh for bbox/containment/graph operations. This is convenient but can be expensive for large models.

Fix suggestions:
- Add “mesh views” / iterators to traverse vertices/triangles without allocation.
- Make `bbox()` compute via iterating all vertices without creating a concatenated mesh.

---

## 4) Tooling / Maintenance

### 4.1 Formatting: `cargo fmt --check` fails

Where:
- `src/geom/polygon/boolean.rs`: at least two blocks around the Sutherland–Hodgman implementation (Rustfmt diff indicates indentation/bracing).

Fix:
- Run `cargo fmt` and commit the formatting changes.

### 4.2 Clippy: `cargo clippy --all-targets -- -D warnings` fails

Current failures include:
- `src/geom/polygon/boolean.rs`: `clippy::type_complexity` (boxed closure tuple); suggest factoring into `type` aliases or a small struct/enum.
- `src/geom/polygon/boolean.rs`: `clippy::needless_range_loop` (index-only loops); switch to iterators.
- `src/geom/building/graph.rs` tests: `clippy::for_kv_map`; iterate `graph.values()` instead of `for (_k, v) in &graph`.

Fix:
- Address the lints and add a CI/pre-push recommendation in README (optional).

---

## 5) Smaller Notes / Nice-to-Haves

- `src/geom/segment.rs`: coplanar segment–polygon intersection is explicitly TODO; this affects `segment_crosses_polygon` and thus point-in-solid if a ray lies in a face plane (edge cases). Consider implementing coplanar handling or ensuring rays avoid coplanarity (jitter direction).
- Consider exposing `EPS`/tolerance configuration publicly (or per-operation parameters) if this is intended as a general geometry library.
- `src/world.rs` is a placeholder; ok, but ensure public API/docs don’t imply it exists yet.

