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

### 1.2 Interior Convection

Two models available, selected via `InteriorConvectionModel`:

**Fixed** (default 3.0 W/m²K): constant convection coefficient.

**TARP** (Walton 1983, buoyancy-driven): temperature- and tilt-dependent:
- Vertical (|cos_tilt| < 0.707): `h = 1.31 * |dT|^(1/3)`
- Horizontal unstable (dT * cos_tilt < 0): `h = 9.482 * |dT|^(1/3) / (7.238 - |cos_tilt|)`
- Horizontal stable: `h = 1.810 * |dT|^(1/3) / (1.382 + |cos_tilt|)`
- Physical minimum: H_MIN = 0.1 W/m²K

When TARP is active, radiative exchange is bypassed (h_conv = h_total = TARP,
t_eff = T_air). The simplified area-weighted MRT model cannot properly return
energy to air; a proper view-factor model would be needed to split conv/rad.

### 1.3 Exterior Convection

Two models available, selected via `ExteriorConvectionModel`:

**Fixed**: constant from ISO 6946 R_se (~25 W/m²K).

**DOE-2 simplified**: combines natural + wind-forced convection:
- Natural: same TARP correlations as interior
- Forced: `h_f = 3.26 + 3.89 * V_wind` (rough surface)
- Combined: `h_ext = sqrt(h_n² + h_f²)`

Wind speed read from EPW weather data via Bus coupling.

### 1.5 Zone Thermal Models

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

### 1.6 Internal Mass Slabs

Non-geometric 1D FVM slabs (e.g., floor slab) coupled to zone air via convective BCs.
Can receive transmitted solar and radiant internal gains when
`use_surface_aware_solar_distribution` is enabled. This provides thermal lag for
solar gains absorbed by interior surfaces.

### 1.7 Solar Gains

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

### 1.8 HVAC

Ideal loads with heating/cooling setpoints. Implicit formulation accounts for
concurrent envelope losses:
```
Q_hvac = C*(T_set - T_zone)/dt + (UA + K_inf)*(T_set - T_out) - Q_gains
```

### 1.9 Infiltration

Constant ACH model:
```
Q_inf = rho * c_p * V * ACH / 3600 * (T_indoor - T_outdoor)
```

Fixed air properties: rho = 1.2 kg/m^3, c_p = 1005 J/(kg*K).

### 1.10 Surface Classification

Polygons are classified at simulation time (not stored on geometry) using a
facing-graph overlay keyed by polygon UID:
- **Exterior**: faces nothing (contributes to envelope UA)
- **SameZoneInterface**: internal partition within a zone (excluded from envelope)
- **InterZoneInterface**: partition between zones (coupling conductance)

### 1.11 Key Assumptions

- **1D heat conduction** only (no lateral flow within wall layers)
- **No per-surface interior temperature nodes** (except FVM wall faces)
- **Simplified interior longwave radiation exchange** (area-weighted MRT from FVM
  wall faces, internal mass slabs, and steady-state window surfaces; uses proper
  half-cell-interpolated surface temperatures, not cell centroids)
- **Simplified exterior LW radiation** (sky/ground exchange enabled on FVM walls;
  uses air temperature for outgoing radiation, not surface temperature)
- **Dynamic or constant h_out**: DOE-2 exterior convection available (natural +
  wind-forced via sqrt combination); ISO 6946 R_se fixed value as fallback
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
| **Exterior surface balance** | Full separate-term balance: absorbed solar + LW radiation (sky/ground/air/surrounding surfaces) + wind-dependent convection + conduction. Each term computed independently. No sol-air lumping. | FVM walls: absorbed solar + LW sky/ground exchange as BC fluxes + DOE-2 wind-dependent convection. Non-FVM: sol-air `(U/h_out) * (alpha*I - eps*delta_R)`. |
| **Interior surface balance** | Full grey-interchange LW model (ScriptF/CarrollMRT) between all zone surfaces. Each surface has its own temperature node. | TARP convection (dT/tilt-dependent) or simplified area-weighted MRT from FVM wall faces, internal mass slabs, and steady-state window surfaces. When TARP is active, radiative exchange is bypassed. |
| **Wall conduction** | Conduction Transfer Functions (CTF) from state-space pre-computation, or Conduction Finite Difference (CondFD). Both capture exact distributed thermal mass. | 1D FVM (similar to CondFD) for eligible exterior opaque surfaces. Windows and inter-zone partitions use steady U*A. |
| **Solar distribution** | Beam solar distributed to surfaces proportional to `Area * Absorptance`. Full ray-tracing option available. | Fraction-based: `transmitted_solar_to_air_fraction` (default 0.4) to air, remainder to internal mass slabs (floor). Not geometry-aware. |
| **Zone air balance** | `C_z * dT/dt = convective_gains + sum(h_i*A_i*(T_si - T_z)) + infiltration + HVAC`. 3rd-order backward difference. Predictor-corrector with HVAC. | Backward Euler on lumped or 2-node model. HVAC applied as ideal loads. |
| **Convection** | TARP, DOE-2, MoWiTT, adaptive; windward/leeward; surface roughness | TARP interior + DOE-2 exterior (simplified, no windward/leeward) |
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
| 600 | Heating | -0.7% | Nearly matched. Small residual from simplified solar distribution and no coupled surface solve. |
| 600 | Cooling | -7.2% | Minor under-prediction. Missing geometry-aware solar distribution and no full view-factor radiation exchange. |
| 900 | Heating | +42.7% | Heavyweight thermal mass dynamics: over-predicts heating losses through massive walls. No coupled interior surface temperature iteration. |
| 900 | Cooling | +26.8% | Over-predicts cooling. Thermal mass stores/releases heat differently than E+'s coupled CTF solve with per-surface radiation exchange. |

### 3.2 Monthly Pattern Analysis

**Case 600 heating (-0.7%):** Nearly perfect annual match. TARP dynamic convection
and corrected insulation thickness (0.066m) bring the model close to EnergyPlus.
Monthly shape shows slight over-prediction in winter (Jan +31%) compensated by
under-prediction in shoulder months (Apr -15%).

**Case 600 cooling (-7.2%):** Good overall match. Summer months (Jul, Aug) are
close; shoulder months (Jan-Mar, Oct-Dec) under-predict cooling due to
simplified solar distribution not tracking beam direction through windows.

**Case 900 heating (+42.7%):** Significantly improved from +106% after fixing TARP
h_min_iso clamp and insulation thickness. Remaining gap is in winter months where
heavy concrete walls should store more daytime solar heat to offset nighttime
losses. The sequential (non-iterative) thermal solve and TARP radiative bypass
likely contribute.

**Case 900 cooling (+26.8%):** Over-predicts cooling, especially in summer (Jun-Aug).
TARP dynamic convection increases coupling between warm surfaces and air, but
without a proper view-factor radiation model, the split between convective and
radiative pathways differs from EnergyPlus's coupled solve.

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

**Status:** Partially addressed. Area-weighted MRT model exists but is bypassed when
TARP is active because the simplified MRT loses energy through incomplete feedback.
TARP treats all interior coupling as convection-only (h_conv = h_total = TARP).

**Symptom:** Affects the rate at which stored solar heat is released from mass.
Main contributor to remaining Case 900 deviation.

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
| **5** | Gap 6: Zone capacity audit | Prevent double-counting | DONE |
| **6** | MRT surface temperature fix | Fix centroid-vs-surface bug | DONE |
| **7** | Dynamic convection (TARP interior + DOE-2 exterior) | Fix Case 600 cooling, Case 900 | DONE |
| **8** | TARP h_min_iso fix + insulation thickness correction | Fix TARP clamp bug + wall R-value | DONE |
| **9** | Gap 5: Ground coupling | Marginal improvement | |
| **10** | Gap 3: Full view-factor radiation exchange | Fix Case 900 deviation | |

After correcting the window U-value (1.8 → 2.8 W/m²K) and updating glass
properties to ASHRAE 140-2020, Case 600 was within ~7% heating / ~7% cooling.
The MRT centroid fix (Phase 6) then moved Case 600 to +12.5% heating / -33%
cooling as it removed a compensating error in the radiation exchange.

Phase 7 (dynamic convection) implemented TARP interior and DOE-2 exterior
convection models, plus `transmitted_solar_to_air_fraction=0.4`. This improved
Case 600 cooling from -53% to -25% and Case 900 cooling from -44% to -5%.

Phase 8 (bug fixes) corrected three issues:
- TARP was clamped to h_min_iso = 7.69 W/m²K (ISO combined coefficient),
  defeating dynamic convection entirely. Removed the clamp.
- Case 600 fiberglass thickness was 0.0776m (RSI 1.94) instead of correct
  0.066m (RSI 1.65) from BESTEST-GSR reference. Fixed.
- Linearized radiation through MRT lost energy. Changed TARP mode to bypass
  radiative exchange (h_conv = h_total = TARP, t_eff = T_air).

After Phase 8: Case 600 heating -0.7%, cooling -7.2% (nearly matched).
Case 900 heating +42.7%, cooling +26.8% (improved but still significant).

The dominant remaining gap for Case 900 is the absence of a proper view-factor
interior radiation exchange model. TARP bypasses radiation entirely, so heat
stored in massive walls can only return to air via convection, missing the
radiative pathway that EnergyPlus models with its ScriptF grey-interchange model.

### 4.2 Validation Strategy

After each phase, rerun:
- `cargo test` (CI validation with `tests/bestest_energy_suite.rs`)
- `examples/bestest_energy_suite` (monthly loads + reference comparison)
- Validate **monthly shapes + peaks**, not just annual totals

Target tolerances (vs EnergyPlus single-run reference):
- Case 600: heating +/- 10%, cooling +/- 10% -- **MET** (H: -0.7%, C: -7.2%)
- Case 900: heating +/- 15%, cooling +/- 15% -- not yet met (H: +42.7%, C: +26.8%)

These are tighter than the ASHRAE 140 inter-program spread but realistic for a
single-reference comparison with the same EPW file.
