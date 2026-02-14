# Next Step: Iterative Surface Heat Balance

## Goal

Replace the current sequential surface → air solve with an **iterative
surface heat balance** that couples all interior surfaces via view-factor
radiation within each substep. Convert the floor from a special-cased
internal mass slab into a normal FVM wall. This enables:

1. **Convection-only air gain** (`q_to_air = h_conv * (T_surf - T_air)`),
   which is physically correct — radiation stays between surfaces.
2. **Wall mass participation** in thermal storage via radiation from floor
   to wall concrete and back.
3. **Non-lagged MRT**: surface temperatures and MRT are self-consistent
   within each substep (not lagged by one substep).
4. **Unified surface model**: all surfaces (walls, ceiling, floor) are
   geometric FVM walls — no special-cased mass slabs for BESTEST.

Expected impact: Case 900 heating from +30% → ~15%, cooling from +15% → ~8%.

---

## Current Architecture (sequential)

Each substep currently does:

```
1. Read T_air from previous substep
2. Read surface temps from previous substep → compute MRT (LAGGED)
3. For each FVM wall:  step solver with BC(h_total, t_eff)
4. For each mass slab: step solver with BC(h_total, t_eff)
5. Sum q_to_air = h_total * (T_surf - t_eff) * A for all surfaces
6. Step zone air model: T_air_new = f(T_air_old, gains, q_to_air, HVAC)
```

Problems:
- MRT uses **previous substep** temperatures → energy conservation only
  holds if `h_total * (T_surf - t_eff)` is used for air gain (radiative
  portion nets to zero only with simultaneous temperatures).
- Cannot separate conv/rad air gains because the rad portion does NOT
  net to zero with lagged temperatures.
- All solar goes to floor mass slab; walls never receive radiative heat
  from floor.
- Floor is a non-geometric internal mass slab, disconnected from the
  actual floor polygon. View factors must match it by `cos_tilt` heuristic.

## New Architecture (iterative)

Each substep will do:

```
1. Read T_air from previous substep (initial guess)
2. OUTER LOOP (max 3-5 iterations, convergence < 0.1°C):
   a. Collect all interior surface temps → compute per-surface MRT
   b. For each FVM surface (walls, ceiling, floor):
        step solver with BC(h_conv + h_rad, t_eff_i)
   c. Compute convective-only air gain:
        q_conv = Σ h_conv_i * (T_surf_i - T_air) * A_i
   d. Step zone air model with q_conv → get new T_air
   e. Check convergence: max|T_surf_new - T_surf_old| < tol
   f. If not converged, RESTORE solver states, update T_air, repeat
3. Accept final surface temps, T_air, HVAC power
```

Key differences:
- **All surfaces are FVM walls** — floor is a regular geometric FVM wall
  with ground-temperature exterior BC, not a separate mass slab.
- **Surfaces and air iterate together** until self-consistent.
- **MRT is current-step** (not lagged) — surfaces exchange radiation
  within the iteration loop.
- **Air gain is convective-only** — radiation stays between surfaces
  and nets to zero by reciprocity (now guaranteed because temperatures
  are simultaneous).
- FVM solver states are saved/restored cheaply between iterations.

---

## Implementation Plan

### Step 0: Convert floor from internal mass slab to FVM wall

Currently the floor polygon exists in the building geometry but is
skipped by `collect_fvm_exterior_walls` because `is_ground_coupled`
returns `true`. It is then re-created as a non-geometric
`InternalMassSurface` with `OneSidedAdiabatic` BC. This is redundant
and disconnects the floor from its actual polygon.

#### 0a. Allow ground-coupled surfaces into FVM collection

**File: `src/sim/energy/simulation.rs`, `collect_fvm_exterior_walls()`**

Remove the early `continue` for ground-coupled surfaces (line 202-204).
Instead, mark the `FvmExteriorWall` as ground-coupled and use a
ground-temperature exterior BC in the step function:

```rust
// Before (skips ground-coupled):
if is_ground_coupled { continue; }

// After (includes ground-coupled):
let is_ground_coupled = is_ground_coupled_exterior_surface(config, &s.path, &poly.vn);
// Don't skip — let it become an FVM wall with ground BC.
```

The step function already handles ground-coupled walls correctly
(lines 579-583): `Convective { h: h_out, t_fluid: T_ground }`.

#### 0b. Remove floor from internal mass config in BESTEST

**File: `examples/bestest_energy_suite/main.rs`**

Remove the `InternalMassSurface` for `floor_mass`. The floor polygon's
own `WallConstruction` (R-25 insulation + concrete/timber) provides the
FVM mesh directly.

```rust
// Remove this block:
// if let Some(floor) = cfg.constructions.get("floor").cloned() {
//     cfg.internal_mass_surfaces.push(InternalMassSurface { ... });
// }
```

#### 0c. Ensure floor construction is applied to floor polygon

The floor polygon path is `zone_one/space/floor/floor`. The construction
is already registered as `cfg.constructions.insert("floor", ...)` with
substring matching, so it will be resolved by `resolve_construction()`
when the FVM collector encounters the floor polygon. No change needed.

#### 0d. Assign solar to the floor FVM wall

With the floor as a regular FVM wall, transmitted solar can be directed
to it via the existing `distribute_transmitted_solar_to_fvm_walls`
mechanism. The beam solar distribution should target the floor polygon
specifically (downward-facing interior surface).

**Modify solar distribution logic** to identify the floor FVM wall(s)
by normal direction (`vn.dz <= -0.5`) and route beam solar there:

```rust
// In the solar distribution block:
// Beam solar → floor FVM wall(s) only (interior face, vn.dz <= -0.5)
// Diffuse solar → area-proportional across all interior FVM surfaces
```

This replaces the current `internal_mass_sources_w_by_zone_uid` path
with direct flux to the floor FVM wall's interior BC.

#### 0e. Update view factor surface enumeration

With the floor as a `SurfaceHandle::Polygon(floor_uid)` instead of
`SurfaceHandle::InternalMass { index }`, the view factor computation
naturally includes it — no more `cos_tilt` matching heuristic needed.

Remove the `match_internal_mass_to_polygons` logic from
`compute_building_view_factors` for the floor case. The floor polygon
is already in the zone's polygon list.

#### 0f. Handle remaining internal mass surfaces

`InternalMassSurface` remains available for truly non-geometric mass
(e.g., furniture, partition walls not in the 3D model). But for BESTEST,
all mass is geometric — no `internal_mass_surfaces` needed.

The `step_internal_mass_surfaces_fill_gains_by_zone_uid` function stays
unchanged for cases that still use it. The iteration loop simply has an
empty mass list for BESTEST.

### Step 1: Add `save_state` / `restore_state` to FVM surfaces

`FvmWallSolver` already implements `Clone`. But cloning the entire
`FvmExteriorWall` vector each iteration is wasteful (copies mesh data).

Add a lightweight state snapshot:

**File: `src/sim/heat_transfer/solver.rs`**

```rust
/// Lightweight snapshot of solver state (temperatures only).
#[derive(Clone)]
pub struct FvmSolverSnapshot {
    temperatures: Vec<f64>,
}

impl FvmWallSolver {
    pub fn save_state(&self) -> FvmSolverSnapshot {
        FvmSolverSnapshot {
            temperatures: self.temperatures.clone(),
        }
    }

    pub fn restore_state(&mut self, snapshot: &FvmSolverSnapshot) {
        self.temperatures.copy_from_slice(&snapshot.temperatures);
    }
}
```

### Step 2: Config additions

**File: `src/sim/energy/config.rs`**

```rust
/// Enable iterative surface heat balance (simultaneous solve of surface
/// temperatures, MRT, and zone air temperature within each substep).
pub use_iterative_surface_balance: bool,  // default: false

/// Maximum iterations per substep for the iterative surface balance.
pub surface_balance_max_iterations: usize,  // default: 4

/// Convergence tolerance for surface temperatures [°C].
pub surface_balance_tolerance_c: f64,  // default: 0.1
```

### Step 3: Restructure the substep loop

**File: `src/sim/energy/simulation.rs`**

This is the core change. Inside the `for _substep in 0..substeps_per_hour`
loop (lines 1756-2259), wrap the FVM wall step + mass step + air step in
an iteration loop.

#### 3a. Before the iteration loop: save solver states

```rust
// Save FVM wall solver states for potential restoration.
let wall_snapshots: Vec<FvmSolverSnapshot> = fvm_walls
    .iter()
    .map(|w| w.solver.save_state())
    .collect();
let mass_snapshots: Vec<FvmSolverSnapshot> = internal_mass_surfaces
    .iter()
    .map(|m| m.solver.save_state())
    .collect();
```

#### 3b. The iteration loop

```rust
let max_iter = if config.use_iterative_surface_balance {
    config.surface_balance_max_iterations
} else {
    1  // no iteration, behaves like current code
};

let mut t_air_iter = t_air_start_c;

for iter in 0..max_iter {
    // Restore solver states if this is a retry (iter > 0).
    if iter > 0 {
        for (w, snap) in fvm_walls.iter_mut().zip(&wall_snapshots) {
            w.solver.restore_state(snap);
        }
        for (m, snap) in internal_mass_surfaces.iter_mut().zip(&mass_snapshots) {
            m.solver.restore_state(snap);
        }
    }

    // 1. Compute per-surface MRT from CURRENT surface temps.
    let (vf_mrt_map, vf_h_r) = compute_vf_mrt(...);

    // 2. Step all FVM surfaces (walls, ceiling, floor) with current
    //    T_air and MRT.
    step_fvm_exterior_walls_fill_gains_by_zone_uid(
        &mut fvm_walls, ..., |_| t_air_iter, ...,
        vf_mrt_map.as_ref(), vf_h_r,
    );

    // 3. Step any remaining internal mass (empty for BESTEST).
    step_internal_mass_surfaces_fill_gains_by_zone_uid(
        &mut internal_mass_surfaces, ..., |_| t_air_iter, ...,
        vf_mrt_map.as_ref(), vf_h_r,
    );

    // 4. Compute CONVECTIVE-ONLY air gains from all surfaces.
    //    (step functions fill gains_out with h_conv-only terms)

    // 5. Step zone air model with q_conv → get new T_air.
    let total_gains = gains_air_w + solar_total_air_w + q_ground
                    + q_fvm_conv + q_mass_conv;
    // ... (free-float → HVAC check → setpoint re-step, same as now)

    // 6. Check convergence.
    let max_dt_surf = max_surface_temp_change(
        &fvm_walls, &internal_mass_surfaces, &prev_surface_temps,
    );
    if max_dt_surf < config.surface_balance_tolerance_c {
        break;  // Converged
    }

    // Save current surface temps for convergence check.
    update_prev_surface_temps(...);

    // Update T_air for next iteration.
    t_air_iter = new_t_air;
}
```

#### 3c. Change air gain computation

Currently `q_to_air` uses `h_total * (T_surf - t_eff)` with VF, which
conflates convective and radiative transport. Change to:

**When iterative balance is active:**
```rust
let q_conv_w_per_m2 = h_conv * (t_surf - t_air);
```

This is correct because within the converged iteration, the radiative
terms between surfaces are self-consistent and net to zero:
`Σ A_i * h_rad * (T_i - T_mrt_i) = 0` (by reciprocity, with
simultaneous temperatures).

**When iterative balance is NOT active (backward compat):**
Keep existing behavior: `h_total * (T_surf - t_eff)` with lagged MRT.

### Step 4: Modify step functions for convective-only air gain

**File: `src/sim/energy/simulation.rs`**

The `step_fvm_exterior_walls_fill_gains_by_zone_uid` function currently
computes `q_to_air` and fills `gains_out`. Add a flag parameter:

```rust
convective_only_air_gain: bool
```

When `true`, air gain is `h_conv * (T_surf - T_air)`.
When `false`, existing behavior: `h_total * (T_surf - t_eff)`.

```rust
// In the gains computation block (lines 602-615):
let q_w_per_m2 = if convective_only_air_gain {
    h_in_conv * (t_surf - t_air)
} else if per_surface_mrt.is_some() {
    h_in_total * (t_surf - t_eff)
} else {
    h_in_conv * (t_surf - t_air)
};
```

Same change in `step_internal_mass_surfaces_fill_gains_by_zone_uid`.

### Step 5: Solar distribution to floor FVM wall

With the floor as a regular FVM wall, solar distribution changes:

**Beam solar → floor FVM wall interior face.**
Identify floor walls by `normal.dz <= -0.5` (downward-facing = floor
interior faces upward into the zone). Apply beam solar as
`ConvectiveWithFluxToDomain` BC on the interior face.

**Diffuse solar → area-proportional across all interior FVM surfaces.**
Same as current `distribute_transmitted_solar_to_fvm_walls` but now
naturally includes the floor polygon alongside walls and ceiling.

**Internal gains → radiant fraction to all interior surfaces.**
Same area-proportional split, now including floor.

The key change is in the solar distribution block (lines 1662-1744):
replace `interior_sources_mass_w_by_zone_uid` (which targeted mass
slabs) with routing beam solar to floor FVM walls. The
`interior_sources_walls_w_by_zone_uid` map already distributes to FVM
wall interior faces — just ensure floor walls are included and beam
solar is directed specifically to them.

```rust
// Beam solar: identify floor FVM walls and route 100% there.
let floor_area: f64 = fvm_walls.iter()
    .filter(|w| w.normal.dz <= -0.5)
    .map(|w| w.area_m2)
    .sum();

// Diffuse solar + internal gains: area-proportional across ALL
// interior FVM surfaces (walls + ceiling + floor).
let all_interior_area: f64 = fvm_walls.iter()
    .map(|w| w.area_m2)
    .sum();
```

### Step 6: HVAC integration within iteration

The free-float → setpoint two-pass strategy stays the same, but now
operates inside the iteration loop:

```
Iteration 1: free-float T_air → check setpoints → compute HVAC
Iteration 2: use HVAC T_air as new guess → re-step surfaces → re-check
Iteration 3: usually converged (HVAC clamps T_air, surfaces adjust)
```

The key insight: with HVAC active, T_air is clamped to the setpoint,
so surface temperatures converge quickly (1-2 extra iterations).

### Step 7: Enable in BESTEST config

**Files: `examples/bestest_energy_suite/main.rs`, `tests/bestest_energy_suite.rs`**

```rust
// Re-enable view factors with iterative balance.
cfg.use_view_factor_radiation = true;
cfg.use_iterative_surface_balance = true;
cfg.surface_balance_max_iterations = 4;
cfg.surface_balance_tolerance_c = 0.1;

// Solar to all FVM surfaces (including floor).
cfg.distribute_transmitted_solar_to_fvm_walls = true;
cfg.use_beam_solar_distribution = true;

// Remove internal mass slab for floor (it's now an FVM wall).
cfg.internal_mass_surfaces.clear();
```

For Case 900, the floor construction (R-25 insulation + concrete slab)
is already registered and will be picked up by the FVM collector via
the existing `resolve_construction("floor")` path matching.

### Step 8: Update `module.rs` (pipeline path)

**File: `src/sim/energy/module.rs`**

Apply the same changes:
- Floor as FVM wall (ground-coupled with FVM)
- Iteration loop in `step()`
- Convective-only air gain when iterative balance is active

### Step 9: UA computation update

**File: `src/sim/energy/simulation.rs`**

Currently the floor polygon is NOT in `fvm_skip_polygons` (it was
skipped before FVM collection). With the floor now becoming an FVM wall,
its polygon UID must be added to `fvm_skip_polygons` so that the
steady-state UA computation excludes it (same as other FVM walls).

This happens automatically — `collect_fvm_exterior_walls` already adds
every collected polygon's UID to the skip set (line 250). Since the
floor polygon is no longer skipped during collection, it will be added
to the skip set and excluded from UA. No code change needed.

---

## Convergence Analysis

Typical BESTEST zone: 8 surfaces (4 walls + floor + ceiling + 2 windows).
Interior h_rad ≈ 5 W/m²K, h_conv ≈ 2.5 W/m²K.

Surface temperature response to 1°C change in T_air:
- FVM wall interior cell: ΔT_surf ≈ h / (h + K_face) ≈ 0.3°C
- The coupling is weak (thermal mass dominates), so iteration converges
  fast: typically 2-3 iterations for 0.1°C tolerance.

Expected overhead: ~3x per substep (3 iterations × 1 FVM solve each).
With 6 substeps/hour × 8760 hours = 52,560 substeps, FVM solve is O(N)
Thomas algorithm with N ≈ 5-10 cells. Total: negligible.

---

## Energy Conservation Verification

After implementing, verify:

```
Σ_i A_i * h_rad * (T_i - T_mrt_i) = 0
```

This should hold to machine precision with converged temperatures.
Add a debug assertion in the iteration loop:

```rust
#[cfg(debug_assertions)]
{
    let rad_imbalance: f64 = zone_surfaces.iter()
        .map(|s| s.area * h_rad * (s.t_surf - s.t_mrt))
        .sum();
    debug_assert!(
        rad_imbalance.abs() < 1e-6,
        "Radiative imbalance: {rad_imbalance} W"
    );
}
```

---

## Files Modified

| File | Change |
|------|--------|
| `src/sim/heat_transfer/solver.rs` | Add `FvmSolverSnapshot`, `save_state()`, `restore_state()` |
| `src/sim/energy/config.rs` | Add 3 iterative balance config fields |
| `src/sim/energy/simulation.rs` | Floor as FVM wall, iteration loop, conv-only air gain, solar to floor FVM |
| `src/sim/energy/module.rs` | Same iteration pattern for pipeline path |
| `src/sim/energy/view_factors.rs` | Remove internal mass polygon matching (floor is now geometric) |
| `examples/bestest_energy_suite/main.rs` | Remove floor mass slab, re-enable VF, enable iterative balance |
| `tests/bestest_energy_suite.rs` | Same config changes + tighten tolerances |
| `SIMULATION.md` | Update results and document iterative balance |

---

## Test Plan

1. **Backward compatibility**: `use_iterative_surface_balance = false`
   with floor as internal mass must produce identical results to current
   code (run before removing the mass slab config).
2. **Floor FVM wall correctness**: floor surface temperature should match
   the old internal mass slab temperature within 0.5°C (same construction,
   same BCs, ground BC ≈ adiabatic given R-25).
3. **Energy conservation**: radiative imbalance < 1e-6 W per zone per
   substep with converged iteration.
4. **Convergence**: typical substeps converge in ≤ 3 iterations at 0.1°C.
5. **BESTEST 600 regression**: should remain within 5% heating, 10% cooling.
6. **BESTEST 900 improvement**: heating should decrease from +30% toward 15%.
7. **BESTEST no-solar cases**: should remain nearly identical (radiation
   coupling has minimal effect when all surfaces are near T_air).

---

## Why This Will Work

The previous VF attempt failed because:
- MRT was lagged → `h_conv * (T_surf - T_air)` for air gain lost energy
- Using `h_total * (T_surf - t_eff)` preserved energy but with wrong
  convective/radiative split, making the model behave as if all heat
  transfer (including radiation) went directly to air

The iterative approach fixes both:
- **MRT is current-step** → radiation between surfaces is self-consistent
- **Air gain is convective-only** → correct physics, and radiation nets
  to zero by reciprocity with simultaneous temperatures
- **All surfaces are FVM walls** → floor is a real geometric polygon with
  its own FVM solver, participating directly in the VF enclosure
- **Wall mass participates** → floor-to-wall radiation redistributes
  stored solar, providing the multi-rate time constant spectrum that
  EnergyPlus has

The floor at 25°C radiates to the 20°C wall interior surface. The wall
concrete (Case 900: 0.10m, k=0.51) absorbs this heat and releases it
slowly through the foam insulation (outward, τ > 10h) and back to the
room via convection (inward, τ ≈ 2.8h). Meanwhile, the floor concrete
(0.08m, k=1.13) releases via convection with τ ≈ 4.4h. Together, these
provide the spectrum of release rates that covers evening through
pre-dawn — exactly what a single mass slab cannot do.
