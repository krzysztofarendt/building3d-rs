# Plan: Add RerunConfig for Customizable Drawing Settings

## Context
All Rerun visualization settings (colors, sizes, labels, animation thresholds) are hardcoded in `src/draw/rerun.rs`. This adds a `RerunConfig` struct following the `SimulationConfig` pattern so users can customize everything programmatically.

## Files to Create/Modify

### 1. Create `src/draw/config.rs` (new file)
- Define `pub type Rgba = (f32, f32, f32, f32);`
- Define `RerunConfig` struct with all pub fields:
  - **Labels**: `session_name: String`, `entity_prefix: String`
  - **General drawing**: `face_color`, `edge_color`, `edge_radius`, `point_color`, `point_radius`
  - **Simulation**: `sim_building_color`, `sim_absorber_color`, `sim_source_color`, `sim_source_radius`, `sim_ray_color_high`, `sim_ray_color_low`, `sim_ray_radius`, `sim_ray_energy_threshold`
- `new()` returns defaults matching current hardcoded values
- `impl Default` delegates to `new()`
- Inline tests for defaults and custom values

### 2. Update `src/draw.rs`
- Add `pub mod config;` alongside existing `pub mod rerun;`

### 3. Update `src/lib.rs`
- Add `pub use draw::config::RerunConfig;` to prelude

### 4. Update `src/draw/rerun.rs`
- Import `RerunConfig` and `Rgba` from `super::config`
- Remove `const SESSION_NAME`
- Add `lerp_color(low, high, t) -> Rgba` helper for ray color interpolation
- Update function signatures (add `config: &RerunConfig` as last param):
  - `start_session(config)` — use `config.session_name`
  - `draw_faces(session, model, rgba, config)` — use `config.entity_prefix`
  - `draw_edges(session, model, radius, rgba, config)` — use `config.entity_prefix`
  - `draw_points(session, model, radius, rgba, config)` — use `config.entity_prefix`
  - `draw_simulation(session, result, building, config)` — use all `sim_*` fields, `lerp_color` for ray energy mapping
- Add inline tests for `lerp_color`

### 5. Update all 6 examples
Each gets: `use building3d::RerunConfig;`, `let config = RerunConfig::new();`, pass `&config` to draw calls.
- `draw_faces.rs` — 2 call sites (start_session, draw_faces)
- `draw_many.rs` — 5 call sites (start_session, draw_faces, 2x draw_edges, draw_points)
- `draw_shapes.rs` — 3 call sites (start_session, draw_edges, draw_faces)
- `floor_plan.rs` — 4 call sites (start_session, draw_faces, draw_edges, draw_points)
- `ray_2_boxes.rs` — 2 call sites (start_session, draw_simulation)
- `ray_teapot.rs` — 2 call sites (start_session, draw_simulation)

## Design Decisions
- **Per-call colors/radii kept** — examples pass different colors per model, so config doesn't replace those params, it adds session-wide settings (entity prefix, session name) and simulation defaults
- **Two ray colors** (`sim_ray_color_high/low`) with lerp replaces hardcoded `(e, 0.0, 0.0, 0.8)`
- **No absorber radius in RerunConfig** — already in `SimulationConfig.absorber_radius`
- **No new dependencies** — pure Rust struct, no serde

## Verification
```bash
cargo test           # config tests + lerp_color tests pass
cargo clippy         # no warnings
cargo fmt --check    # formatted
cargo build --examples  # all examples compile
```
