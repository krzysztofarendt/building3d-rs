# Energy Simulation Physics

This document describes the physical models, assumptions, and known gaps in the
building3d energy simulation engine, with a systematic comparison to EnergyPlus.

---

## 1. Current Physical Model

### 1.1 Wall Conduction

Two conduction models are available, selected per surface:

**Steady-state U*A** (for surfaces without layered constructions or with U-value overrides):
```
Q = U * A * (T_indoor - T_outdoor)     [W]
```

U-values are computed per ISO 6946:
```
R_total = R_se + sum(t_i / lambda_i) + R_si
U = 1 / R_total
```

**1D Finite Volume Method** (for eligible opaque exterior surfaces with `WallConstruction`):

Each surface gets its own FVM solver instance. The heat equation is discretized
over cells with Backward Euler (unconditionally stable at any dt):
```
(rho * c_p * V_i / dt) * (T_i^{n+1} - T_i^n) = sum_faces k_f * A_f / d_f * (T_j^{n+1} - T_i^{n+1}) + Q_i
```

Boundary conditions are convective on both sides:
- Interior: `q = h_in * (T_surface - T_air)`
- Exterior: `q = h_out * (T_surface - T_outdoor)`, with optional absorbed solar flux

The tridiagonal system is solved with the Thomas algorithm (O(N)).

FVM eligibility: exterior, has resolved `WallConstruction`, not glazing, not U-value
override, not ground-coupled.

### 1.2 Zone Thermal Models

**AirOnly (1R1C)**: Single zone air node with lumped capacity.
```
C = V * thermal_capacity_j_per_m3_k     [J/K]
T_zone(t+dt) = (T_zone(t) + dt/C * (Q_gains + Q_hvac + K*T_out)) / (1 + dt*K/C)
```

**EnvelopeRc2R1C**: Two-node model (air + envelope mass) per zone.
Exterior conductance is split into two equal resistances with an envelope mass node
between them. Heat sources are split between nodes (internal gains to air, exterior
absorbed solar to envelope).

**Multi-zone coupling**: Inter-zone partitions become conductances `K_ij = U_eq * A`.
Backward Euler yields a linear system solved per timestep (Gaussian elimination).

### 1.3 Internal Mass Slabs

Non-geometric 1D FVM slabs (e.g., floor slab) coupled to zone air via convective BCs.
Can receive transmitted solar and radiant internal gains when
`use_surface_aware_solar_distribution` is enabled. This provides thermal lag for
solar gains absorbed by interior surfaces.

### 1.4 Solar Gains

Transmitted solar through glazing:
```
Q_solar = (DNI * cos(theta_inc) + DHI * sky_view + ground_reflected) * A * SHGC
```

Exterior opaque sol-air absorption:
```
Q_opaque = (U / h_out) * alpha * I_incident * A
```

Solar position: Spencer (1971) with equation-of-time and longitude correction.
EPW hour-ending convention handled via mid-hour timestamp (`hour - 0.5`).

Angular transmittance modifier (two options, selected via config):
- **ASHRAE IAM:** `IAM = 1 - a * (1/cos(theta) - 1)`, clamped to [0, 1].
- **Polynomial:** `SHGC(theta) = SHGC_0 * sum(c[i] * cos^i(theta))` for i=0..5,
  with preset coefficients derived from Fresnel + absorption fit for single-pane
  clear glass (n=1.526, BESTEST Glass Type 1). The polynomial is more physically
  accurate at high incidence angles (70-80°) but gives nearly identical annual
  results to ASHRAE IAM a=0.1 for single-pane glass.

### 1.5 HVAC

Ideal loads with heating/cooling setpoints. Implicit formulation accounts for
concurrent envelope losses:
```
Q_hvac = C*(T_set - T_zone)/dt + (UA + K_inf)*(T_set - T_out) - Q_gains
```

### 1.6 Infiltration

Constant ACH model:
```
Q_inf = rho * c_p * V * ACH / 3600 * (T_indoor - T_outdoor)
```

Fixed air properties: rho = 1.2 kg/m^3, c_p = 1005 J/(kg*K).

### 1.7 Surface Classification

Polygons are classified at simulation time (not stored on geometry) using a
facing-graph overlay keyed by polygon UID:
- **Exterior**: faces nothing (contributes to envelope UA)
- **SameZoneInterface**: internal partition within a zone (excluded from envelope)
- **InterZoneInterface**: partition between zones (coupling conductance)

### 1.8 Key Assumptions

- **1D heat conduction** only (no lateral flow within wall layers)
- **No per-surface interior temperature nodes** (except FVM wall faces)
- **Simplified interior longwave radiation exchange** (area-weighted MRT from FVM
  wall faces, internal mass slabs, and steady-state window surfaces)
- **Simplified exterior LW radiation** (sky/ground exchange enabled on FVM walls;
  uses air temperature for outgoing radiation, not surface temperature)
- **Constant h_out** from ISO 6946 R_se (wind-dependent model available but disabled
  by default)
- **Constant SHGC** (no spectral dependence; polynomial angular model available)
- **Constant infiltration** (no wind/stack pressure dependence)
- **Constant ground temperature** boundary
- **No latent loads** (sensible only)
- **No thermal bridging** (clear-field U-values only)
- **Ideal HVAC** (instantaneous response, perfect setpoint control)

---

## 2. Comparison to EnergyPlus

### 2.1 Side-by-Side

| Aspect | EnergyPlus | building3d |
|--------|-----------|------------|
| **Exterior surface balance** | Full separate-term balance: absorbed solar + LW radiation (sky/ground/air/surrounding surfaces) + wind-dependent convection + conduction. Each term computed independently. No sol-air lumping. | FVM walls: absorbed solar + LW sky/ground exchange as BC fluxes. Non-FVM: sol-air `(U/h_out) * (alpha*I - eps*delta_R)`. Constant h_out from ISO 6946 R_se (wind model available). |
| **Interior surface balance** | Full grey-interchange LW model (ScriptF/CarrollMRT) between all zone surfaces. Each surface has its own temperature node. | Simplified area-weighted MRT from FVM wall faces, internal mass slabs, and steady-state window surfaces. |
| **Wall conduction** | Conduction Transfer Functions (CTF) from state-space pre-computation, or Conduction Finite Difference (CondFD). Both capture exact distributed thermal mass. | 1D FVM (similar to CondFD) for eligible exterior opaque surfaces. Windows and inter-zone partitions use steady U*A. |
| **Solar distribution** | Beam solar distributed to surfaces proportional to `Area * Absorptance`. Full ray-tracing option available. | Fraction-based: `transmitted_solar_to_air_fraction` + remainder to mass node or internal mass slabs. Not geometry-aware. |
| **Zone air balance** | `C_z * dT/dt = convective_gains + sum(h_i*A_i*(T_si - T_z)) + infiltration + HVAC`. 3rd-order backward difference. Predictor-corrector with HVAC. | Backward Euler on lumped or 2-node model. HVAC applied as ideal loads. |
| **Window thermal** | ISO 15099 layer-by-layer spectral calculation. Angular transmittance/reflectance per layer. | Bulk U-value override + constant SHGC. Polynomial angular transmittance (Fresnel-based) or ASHRAE IAM. |
| **Ground coupling** | Kiva 2D finite difference solver with detailed soil domain. | Constant ground temperature boundary. |
| **Infiltration** | Design flow rate, effective leakage area (Sherman-Grimsrud), flow coefficient (AIM-2), or full airflow network. | Constant ACH. |

### 2.2 Architectural Differences

EnergyPlus solves **coupled** surface heat balances: at each timestep, exterior
surface temperature, interior surface temperature, and zone air temperature are
solved simultaneously via iteration. Each surface has explicit inside and outside
temperature nodes. This means:

1. Solar absorbed on an exterior opaque surface raises the outside surface temperature,
   which drives conduction inward through the wall (with thermal lag from CTF/CondFD),
   which raises the inside surface temperature, which transfers heat to zone air via
   convection and to other surfaces via LW radiation.

2. Transmitted solar hitting the floor is absorbed at the floor's inside surface
   temperature node, stored in the floor's thermal mass, and released slowly via
   convection + LW radiation over subsequent hours.

building3d approximates these paths: FVM walls provide the conduction path (1) for
eligible surfaces, but the interior coupling (step 2) uses fraction-based splits
rather than explicit surface temperature nodes.

---

## 3. BESTEST Validation Status

We validate against ASHRAE Standard 140 / BESTEST cases 600 (lightweight) and 900
(heavyweight), comparing to OpenStudio/EnergyPlus reference values using a Boston
Logan TMY3 EPW file.

### 3.1 Current Discrepancies

| Case | Metric | Annual Error vs E+ | Primary Cause |
|------|--------|-------------------|---------------|
| 600 | Heating | -29.3% | Remaining gap from simplified h_out + no interior LW exchange |
| 600 | Cooling | +18.3% | Solar over-prediction (simplified distribution, no interior LW) |
| 900 | Heating | +20.7% | LW sky loss + thermal mass dynamics over-predicting heating |
| 900 | Cooling | +9.2% | Improved by routing solar to internal mass slabs |

### 3.2 Monthly Pattern Analysis

**Case 600 heating (-29.3%):** building3d under-predicts heating. Exterior LW sky
exchange is now included on FVM walls, which improved this from -34.3% to -29.3%.
Remaining gap likely from constant h_out and missing interior LW radiation exchange.

**Case 600 cooling (+18.3%):** Over-prediction in summer from solar gains being too
high. The angular SHGC polynomial (Phase 3) had negligible effect because the
Fresnel+absorption curve for single-pane glass closely matches the ASHRAE IAM.
Remaining gap is primarily from simplified solar distribution and missing interior
LW radiation exchange.

**Case 900 heating (+20.7%):** Over-predicts heating. The LW sky correction increased
heating demand beyond the reference. May indicate the LW term is too strong without
a coupled surface temperature solve, or that the thermal mass seasonal buffering
needs tuning. Slightly improved from +21.7% by the polynomial angular model.

**Case 900 cooling (+9.2%):** Significantly improved from +22.2% by routing
transmitted solar to internal mass slabs (FVM floor) rather than directly to zone
air. The thermal lag from the FVM slab delays solar release.

### 3.3 Compensating Errors

Several BESTEST totals appear "close" for the wrong reasons. For example, a too-low
window U-value can partially cancel missing opaque exterior solar absorption. Fixes
must be validated on **monthly shapes + peaks**, not just annual totals.

---

## 4. Gap Analysis and Fix Plan

### Gap 1: Exterior Surface Heat Balance -- DONE

**Status:** Implemented. LW sky/ground radiation exchange added to FVM wall exterior
BCs. Wind-dependent h_out infrastructure available but kept disabled by default.

**Effect:** Case 600 heating improved from -34.3% to -29.2%. Case 600 cooling
improved from +23.2% to +18.4%. Remaining gap likely from constant h_out and
missing interior LW radiation exchange.

### Gap 2: Interior Solar Distribution -- DONE

**Status:** Implemented. Transmitted solar routed to internal mass slabs (floor)
as surface heat flux via `ConvectiveWithFluxToDomain` BCs.

**Effect:** Case 900 cooling improved from +22.2% to +10.0%. The FVM slab provides
realistic thermal lag for absorbed solar.

### Gap 3: Interior Longwave Radiation Exchange (MEDIUM)

**Problem:** No proper radiation exchange between interior surfaces. EnergyPlus uses
a full grey-interchange model (ScriptF) that redistributes heat from warm surfaces
(sun-heated floor) to cooler surfaces (walls, ceiling) and then to air.

**Symptom:** Affects the rate at which stored solar heat is released from mass.

**Fix:**
1. Short term: improve the existing `use_interior_radiative_exchange` model by
   including all interior surface temperatures (not just FVM wall faces).
2. Medium term: implement a radiation star network where each surface exchanges LW
   with a zone mean radiant node weighted by area*emissivity.
3. Long term: full view-factor calculation (can reuse ray infrastructure).

### Gap 4: Angular SHGC -- DONE

**Status:** Implemented. 5th-order polynomial angular transmittance model
`SHGC(theta) = SHGC_0 * sum(c[i]*cos^i(theta))` with preset coefficients derived
from Fresnel reflectance + absorption fit for single-pane clear glass (n=1.526).
Applied consistently in all three SHGC computation paths (unshaded, shaded, legacy).
Falls back to ASHRAE IAM when polynomial coefficients are not set.

**Effect:** Negligible for single-pane clear glass (BESTEST Glass Type 1). The
Fresnel+absorption transmittance curve closely matches ASHRAE IAM a=0.1 at most
angles. Annual changes were <0.2% for Case 600 and ~1% for Case 900. The angular
model is not a significant contributor to the remaining cooling over-prediction;
that gap comes from other sources (solar distribution, interior LW exchange).

**Future value:** The polynomial infrastructure will have more impact with multi-pane
low-e glazing coefficients, where angular drop-off is steeper than single-pane.

### Gap 5: Ground Coupling (LOW for BESTEST)

**Problem:** Constant ground temperature. Real ground temperature varies seasonally
and depends on building footprint geometry.

**Impact:** BESTEST floors are heavily insulated (R-25), so this is a small term.

**Fix:** Implement Kasuda correlation for monthly ground temperature profile, or a
simplified 2D finite difference model.

### Gap 6: Zone Air Capacitance (LOW)

**Problem:** When FVM walls + internal mass slabs handle structural mass, the
remaining lumped capacity should represent only air + furniture. The capacity
deduction logic may not be fully correct.

**Fix:** Verify the deduction logic. Consider using `C_z = rho_air * Cp_air * V`
(pure air) when FVM walls + internal mass slabs handle all structural mass.

### 4.1 Recommended Implementation Sequence

| Phase | Gap | Target | Status |
|-------|-----|--------|--------|
| **1** | Gap 2: Interior solar distribution | Fix Case 900 cooling | DONE |
| **2** | Gap 1: Exterior surface heat balance | Fix Case 600 heating | DONE |
| **3** | Gap 4: Angular SHGC | Improve seasonal accuracy | DONE (negligible effect) |
| **4** | Gap 3: Interior LW exchange (windows in MRT) | Improve Case 600 heating | DONE |
| **5** | Gap 6: Zone capacity audit | Prevent double-counting | |
| **6** | Gap 5: Ground coupling | Marginal improvement | |

After Phase 4 (window surfaces in MRT), Case 600 is within ~17% heating / ~1%
cooling and Case 900 within ~54% heating / ~11% cooling. Phase 4 added
steady-state interior surface temperatures for non-FVM exterior surfaces
(windows) to the area-weighted MRT. This correctly increases heating demand by
accounting for radiative exchange between warm walls and cold window surfaces.
Case 600 heating improved significantly (-29% → -17%). Case 900 heating
worsened (+21% → +54%) because the heavyweight model already over-predicted
heating; the root cause is likely elsewhere (thermal mass dynamics, h_in
calibration, or envelope capacity accounting).

### 4.2 Validation Strategy

After each phase, rerun:
- `cargo test` (CI validation with `tests/bestest_energy_suite.rs`)
- `examples/bestest_energy_suite` (monthly loads + reference comparison)
- Validate **monthly shapes + peaks**, not just annual totals

Target tolerances (vs EnergyPlus single-run reference):
- Case 600: heating +/- 10%, cooling +/- 10%
- Case 900: heating +/- 15%, cooling +/- 15%

These are tighter than the ASHRAE 140 inter-program spread but realistic for a
single-reference comparison with the same EPW file.
