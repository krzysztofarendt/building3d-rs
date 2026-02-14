# Energy Simulation: Physics, Methods & Validation

This document describes the physical models, numerical methods, EnergyPlus comparison,
and validation status of the building3d energy simulation engine.

Reference: ASHRAE Standard 140 / BESTEST Cases 600 (lightweight) and 900 (heavyweight),
validated against OpenStudio/EnergyPlus using Boston Logan TMY3 EPW.

---

## 1. Current Status & Gaps

### 1.1 Current Results

Reference annual values (OpenStudio/EnergyPlus):
- Case 600: heating = 4324.8 kWh, cooling = 6044.1 kWh
- Case 900: heating = 1661.2 kWh, cooling = 2498.2 kWh

| Case | Heating | Cooling | Status |
|------|---------|---------|--------|
| 600 | **+11.4%** (4816 vs 4325 kWh) | **-0.5%** (6013 vs 6044 kWh) | Good |
| 900 | **+58.6%** (2635 vs 1661 kWh) | **+31.8%** (3293 vs 2498 kWh) | Major gap |

Target tolerances: +/- 15% (tighter than ASHRAE 140 inter-program spread).
Case 600 meets target. Case 900 does not.

### 1.2 Monthly Error Patterns

The monthly pattern reveals a **seasonal phase error** diagnostic of incorrect
solar distribution:

**Case 900 heating (annual +58.6%):**

| Month | Simulated | Reference | Error | Pattern |
|-------|-----------|-----------|-------|---------|
| Jan | 784 kWh | 255 kWh | +207% | Winter: massive over-prediction |
| Mar | 364 kWh | 125 kWh | +190% | Winter: massive over-prediction |
| Apr | 132 kWh | 268 kWh | -51% | Shoulder: under-prediction |
| Oct | 36 kWh | 81 kWh | -56% | Shoulder: under-prediction |
| Dec | 499 kWh | 308 kWh | +62% | Winter: over-prediction |

In winter, the model loses too much heat (solar leaking through exterior walls
instead of being stored in mass). In shoulder months, mass should release stored
winter solar to offset heating, but was never sufficiently charged.

**Case 900 cooling (annual +31.8%):**

| Month | Simulated | Reference | Error | Pattern |
|-------|-----------|-----------|-------|---------|
| Jan | 12 kWh | 53 kWh | -77% | Cold months: under-prediction |
| Jul | 754 kWh | 400 kWh | +89% | Summer peak: over-prediction |
| Aug | 705 kWh | 513 kWh | +37% | Summer: over-prediction |
| Jun | 409 kWh | 404 kWh | +1% | Nearly perfect |
| Sep | 531 kWh | 493 kWh | +8% | Good |

Summer peaks are over-predicted (too much solar reaching zone air without mass
buffering). Shoulder months (Jun, Sep) are nearly perfect.

### 1.3 Error Budget

| Error Source | 900 Heating (+59pp) | 900 Cooling (+32pp) |
|-------------|---------------------|---------------------|
| Remaining wall solar leakage (diffuse) | ~15pp | ~-10pp (helps cooling) |
| Solar timing / surface selection | ~20pp | ~20pp |
| Interior film coefficient (h_in slightly low) | ~5pp | ~5pp |
| Remaining / other | ~19pp | ~17pp |

**Beam solar distribution** (beam to floor, diffuse area-proportional) is now
active. This reduced wall leakage vs the prior uniform area split, improving
heating by 5pp. However, the reduced leakage also worsened cooling by 11pp
because wall leakage was a beneficial heat-dumping path in summer.

**Trade-off:** Wall leakage simultaneously helps cooling and hurts heating.
Without geometry-aware per-surface distribution, this tension cannot be fully
resolved. The remaining gap requires distributing solar to the RIGHT surface
at the RIGHT time (floor in winter, north wall mass in summer).

**Thermal mass model is sound:** no-solar cases show only 1.7% difference
between Case 600 (7508 kWh) and Case 900 (7378 kWh).

**Envelope UA is correct:** 80.96 W/K matches hand-calculated 80.9 W/K.

### 1.4 Prioritized Next Steps

| Priority | Improvement | Expected Impact | Effort |
|----------|------------|-----------------|--------|
| 1 | **Sun-angle solar projection** | Break heating/cooling trade-off | Medium |
| | Hourly profile angle → beam to floor (winter) or north wall (summer). | | |
| | North wall has mass inboard of insulation = minimal leakage. | | |
| 2 | **Exclude window-side wall from diffuse** | 900H: -5pp | Low |
| | Solar through south windows can't hit south wall interior. Remove it | | |
| | from diffuse distribution denominator. | | |
| 3 | **Higher-order time integration** | ~2-3pp on all cases | Medium |
| | 3rd-order BDF or predictor-corrector for zone air temperature. | | |

---

## 2. Physical Models

### 2.1 Wall Conduction

Two conduction models, selected per surface:

**Steady-state U*A** (glazing, ground-coupled, U-value overrides):
```
Q = U * A * (T_indoor - T_outdoor)
R_total = R_se + sum(t_i / k_i) + R_si       [ISO 6946]
```

**1D Finite Volume Method** (eligible opaque exterior surfaces with `WallConstruction`):

Each surface gets its own FVM solver. Backward Euler, unconditionally stable:
```
(rho * Cp * V_i / dt) * (T_i^{n+1} - T_i^n) = sum K_face * (T_j - T_i) + Q_i
```

Face conductance at material interfaces uses series resistance:
```
K_face = A / (half_dx_L / k_L + half_dx_R / k_R)
```

Boundary convective BC accounts for half-cell resistance:
```
K_eff = 1 / (1/(h*A) + 1/K_face)
```

Tridiagonal system solved with Thomas algorithm (O(N)). Max cell size 0.05m.

FVM eligibility: exterior, has `WallConstruction`, not glazing, not U-value
override, not ground-coupled.

### 2.2 Interior Convection

Two models, selected via `InteriorConvectionModel`:

**Fixed** (default 3.0 W/m²K): constant coefficient.

**TARP** (Walton 1983, buoyancy-driven, temperature- and tilt-dependent):
- Vertical (|cos_tilt| < 0.707): `h = 1.31 * |dT|^(1/3)`
- Horizontal unstable (dT * cos_tilt < 0): `h = 9.482 * |dT|^(1/3) / (7.238 - |cos_tilt|)`
- Horizontal stable: `h = 1.810 * |dT|^(1/3) / (1.382 + |cos_tilt|)`
- Minimum: H_MIN = 0.1 W/m²K

TARP returns **convection-only** values (typically 1.5-4.0 W/m²K). When view
factors are disabled, radiative exchange is bypassed (h_total = h_conv, t_eff = T_air).
When view factors are enabled, radiation is handled separately (see 2.4).

### 2.3 Exterior Convection

Two models, selected via `ExteriorConvectionModel`:

**Fixed**: ISO 6946 R_se (~25 W/m²K).

**DOE-2 simplified**: natural + wind-forced convection:
- Natural: TARP correlations (same as interior)
- Forced: `h_f = 3.26 + 3.89 * V_wind` (rough surface)
- Combined: `h_ext = sqrt(h_n^2 + h_f^2)`

No windward/leeward distinction or surface roughness categories.

### 2.4 Interior Longwave Radiation (View Factors)

Optional, enabled via `use_view_factor_radiation`. Provides proper radiative
exchange between all interior zone surfaces.

**View factor computation:**
- Monte Carlo ray casting: N cosine-weighted rays per surface (default 10,000)
- Per-zone FlatScene with voxel-accelerated intersection
- Building polygons have outward normals; VF code negates them for interior rays
- Reciprocity enforced: `A_i * F_ij = A_j * F_ji`; rows normalized to sum ~1.0

**Per-timestep radiation exchange:**
```
h_rad  = eps * 4 * sigma * T_mean^3              # uniform (~5.1 W/m²K at 20C)
T_mrt_i = sum_j(F_ij * T_j)                      # per-surface MRT
h_total = h_conv(TARP) + h_rad
t_eff   = (h_conv * T_air + h_rad * T_mrt_i) / h_total
```

**Energy conservation:** Guaranteed by reciprocity:
`sum(A_i * h_rad * (T_i - T_mrt_i)) = 0`

Combined h_in ~7 W/m²K (TARP ~2.5 + radiation ~4.6), matching ISO 6946
R_si = 0.13 m²K/W = 7.69 W/m²K.

**Air gain:** Must use `h_total * (T_surf - t_eff)`, not `h_conv * (T_surf - T_air)`.
The radiative portion nets to zero by reciprocity, but omitting it from the zone
balance loses energy.

### 2.5 Solar Gains & Distribution

**Transmitted solar through glazing:**
```
Q_direct  = DNI * cos(theta_inc) * A * SHGC * IAM
Q_diffuse = DHI * sky_view * A * SHGC
Q_ground  = GHI * rho_ground * ground_view * A * SHGC
```

**Angular transmittance (two options):**
- Polynomial: `SHGC(theta) = SHGC_0 * sum(c[i] * cos^i(theta))` for i=0..5
  (Fresnel + absorption fit for single-pane clear glass, n=1.526)
- ASHRAE IAM: `IAM = 1 - a * (1/cos(theta) - 1)`, clamped to [0, 1]

**Exterior opaque sol-air absorption:**
```
Q_opaque = (U / h_out) * alpha * I_incident * A
```

Applied as `ConvectiveWithFlux` BC on FVM walls (fraction alpha entering domain).

**Interior solar distribution** (`transmitted_solar_to_air_fraction = 0.0`):

When `use_beam_solar_distribution = true` (recommended):
- **Beam (direct)** component → 100% to internal mass slabs (floor)
- **Diffuse** component → area-proportional split: floor (~30%) + FVM walls (~70%)
- `distribute_transmitted_solar_to_fvm_walls = true`

This approximates EnergyPlus beam tracking: beam through south windows mostly
hits the floor. Diffuse remains area-proportional. Reduces wall solar leakage
compared to uniform area split, improving heating accuracy.

When `use_beam_solar_distribution = false` (legacy):
- All transmitted solar split area-proportionally between floor and FVM walls.

**Solar position:** Spencer (1971) with equation-of-time and longitude correction.
EPW hour-ending convention handled via mid-hour timestamp.

### 2.6 Zone Thermal Models

**AirOnly (1R1C):** Single zone air node with lumped capacity.
```
T(t+dt) = [T(t) + (dt/C) * (Q_gains + Q_hvac + K*T_out)] / (1 + dt*K/C)
```

**EnvelopeRc2R1C:** Two-node model (air + envelope mass) per zone. Exterior
conductance split into two equal resistances with envelope mass node between them.

**Multi-zone coupling:** Inter-zone partitions become conductances `K_ij = U * A`.
Backward Euler yields a linear system solved per timestep (Gaussian elimination).

### 2.7 Internal Mass Slabs

Non-geometric 1D FVM slabs (e.g., floor slab, concrete block mass) coupled to
zone air via convective BCs. Receive transmitted solar and radiant internal gains
via `ConvectiveWithFluxToDomain`. Provides thermal lag for absorbed solar.

When view factors are enabled, mass slabs are matched to building polygons by
zone and tilt (floor → polygon with vn.dz <= -0.5) and inherit their view factor
row in the F_ij matrix.

### 2.8 HVAC, Infiltration, Ground

**HVAC:** Ideal loads with heating/cooling setpoints. Implicit formulation:
```
Q_hvac = C*(T_set - T_zone)/dt + (UA + K_inf)*(T_set - T_out) - Q_gains
```

**Infiltration:** Constant ACH model:
```
Q_inf = rho * Cp * V * ACH / 3600 * (T_indoor - T_outdoor)
```
Fixed air properties: rho = 1.2 kg/m³, Cp = 1005 J/(kg*K). No wind or stack
dependence. BESTEST specifies 0.5 ACH.

**Ground coupling:** Constant ground temperature boundary (default 10C). Applied
as correction to outdoor temperature for ground-coupled surfaces. BESTEST floors
are heavily insulated (R-25), so this is a negligible term (UA_ground ~1.89 W/K).

### 2.9 Key Assumptions

- **1D heat conduction** only (no lateral flow within wall layers)
- **Constant material properties** (no temperature-dependent k or Cp)
- **Simplified solar distribution** (beam to floor, diffuse area-proportional; not geometry-traced)
- **Simplified exterior LW radiation** (uses air temperature, not surface temperature)
- **Constant SHGC** (no spectral dependence; polynomial angular model available)
- **Constant infiltration** (no wind/stack pressure dependence)
- **Constant ground temperature** boundary
- **No latent loads** (sensible only)
- **No thermal bridging** (clear-field U-values only)
- **Ideal HVAC** (instantaneous response, perfect setpoint control)
- **1st-order time integration** (Backward Euler, vs EnergyPlus 3rd-order BDF)

---

## 3. Comparison to EnergyPlus

### 3.1 Core Methods

| Aspect | EnergyPlus | building3d | Gap |
|--------|-----------|------------|-----|
| Wall conduction | CTF (default) or CondFD (Crank-Nicolson) | 1D FVM, Backward Euler, Thomas algorithm | Comparable |
| Interior convection | TARP, Ceiling Diffuser, adaptive | TARP natural convection | Minor |
| Exterior convection | DOE-2, MoWiTT; windward/leeward | DOE-2 simplified; no wind direction | Moderate |
| Interior LW radiation | View-factor matrix, ScriptF grey-interchange | Monte Carlo view factors, uniform h_rad | Comparable |
| Exterior LW radiation | Sky temp model, cloud cover | EPW horizontal IR, linearized h_rad | OK |
| Time integration | 3rd-order BDF, predictor-corrector | Backward Euler (1st order) | Minor |
| Air capacitance | `V * rho * Cp / dt` with multiplier | Configurable J/(m³*K) | OK |

### 3.2 Solar & Glazing

| Aspect | EnergyPlus | building3d | Gap |
|--------|-----------|------------|-----|
| Solar distribution | Full beam tracking to specific surfaces | Beam to floor, diffuse area-proportional | **Moderate** |
| Window thermal | ISO 15099 multi-layer, spectral | Bulk U-value + SHGC | Large |
| Window angular | Per-layer spectral polynomial | Single SHGC with IAM polynomial | Moderate |
| Shading | Full 3D shadow with self-shading | FlatScene ray casting (optional) | Comparable |
| Ground reflectance | View-factor based | `GHI * rho * ground_view` | OK |

EnergyPlus `FullInteriorAndExterior` traces beam solar through windows to specific
interior surfaces based on geometry. Low winter sun hits the floor; high summer
sun hits the back wall. This is the single most important difference for Case 900.

### 3.3 Infiltration, Ground, HVAC

| Aspect | EnergyPlus | building3d | Gap |
|--------|-----------|------------|-----|
| Infiltration | ELA, crack method, wind/stack | Constant ACH | Moderate |
| Ground coupling | Kiva 2D transient FD | Constant temperature | Moderate (low for BESTEST) |
| HVAC | Full component-based | Ideal loads | By design |
| Moisture | Full moisture balance | Not modeled | By design |

### 3.4 Architectural Difference

EnergyPlus solves **coupled** surface heat balances: exterior surface temperature,
interior surface temperature, and zone air temperature are solved simultaneously
via iteration. Each surface has explicit inside and outside temperature nodes.

Solar absorbed on an exterior opaque surface raises outside surface temperature,
drives conduction inward (with thermal lag from CTF/CondFD), raises inside surface
temperature, transfers to zone air via convection and to other surfaces via LW
radiation.

building3d approximates this path: FVM walls provide the conduction path for
eligible surfaces. Interior coupling uses view-factor radiation + TARP convection.
The main gap is solar distribution: which surfaces receive the solar flux.

---

## 4. Validation History

BESTEST Cases 600 & 900 vs OpenStudio/EnergyPlus reference.

### 4.1 Results by Commit

| Commit | Description | 600 Heat | 600 Cool | 900 Heat | 900 Cool |
|--------|-------------|----------|----------|----------|----------|
| 536a19f | Merge heat into fvm | -3.8% | -53.2% | +87.6% | -43.7% |
| bf4ae52 | TARP + DOE-2 + solar-to-air=0.4 | +13.1% | -25.0% | +105.8% | -5.3% |
| e9c8993 | Fix TARP clamp + insulation + radiative bypass | -0.7% | -7.2% | +42.7% | +26.8% |
| 5d9645f | View-factor interior radiation | +9.9% | +4.1% | +63.7% | +50.3% |
| 6c777c3 | Solar-to-air=0.0 (all solar to mass) | +9.2% | +3.6% | +53.5% | +42.0% |
| cd243bf | Distribute solar to all interior surfaces | +13.1% | -5.7% | +63.5% | +20.6% |
| PENDING | Beam solar to floor, diffuse area-proportional | +11.4% | -0.5% | +58.6% | +31.8% |

### 4.2 Ablation: Dynamic Convection (on top of 536a19f)

| Change | 600 Heat | 600 Cool | 900 Heat | 900 Cool |
|--------|----------|----------|----------|----------|
| Baseline (Fixed h_in=3.0, Fixed h_out) | -3.8% | -53.2% | +87.6% | -43.7% |
| + TARP interior convection | +12.5% | -33.0% | +101.9% | -24.5% |
| + DOE-2 exterior convection | +13.0% | -32.9% | +102.3% | -24.8% |
| + transmitted_solar_to_air_fraction=0.4 | +13.1% | -25.0% | +105.8% | -5.3% |

TARP is the biggest single lever: +20pp improvement on cooling, but worsens
heating by ~15pp. DOE-2 exterior is negligible (~0.3pp). Solar-to-air=0.4
improves cooling further (900C from -24.5% to -5.3%).

### 4.3 Ablation: View Factors + Solar Distribution (on top of e9c8993)

| Change | 600 Heat | 600 Cool | 900 Heat | 900 Cool |
|--------|----------|----------|----------|----------|
| Baseline (TARP bypass, no VF) | -0.7% | -7.2% | +42.7% | +26.8% |
| + View factors (VF on) | +9.9% | +4.1% | +63.7% | +50.3% |
| + Solar-to-air=0.0 | +9.2% | +3.6% | +53.5% | +42.0% |
| + Solar to FVM walls (area-proportional) | +13.1% | -5.7% | +63.5% | +20.6% |

View factors removed "phantom insulation" from TARP-only (R_si 0.40 → 0.14 m²K/W),
exposing other errors. Solar-to-air=0.0 matches EnergyPlus FullInteriorAndExterior.
Wall solar dilutes floor flux 3.3x, greatly improving 900C but worsening 900H
via 7% wall leakage.

### 4.4 Ablation: Solar Distribution Alternatives (on top of 6c777c3)

| Change | 600 Heat | 600 Cool | 900 Heat | 900 Cool |
|--------|----------|----------|----------|----------|
| Floor-only solar (baseline) | +9.2% | +3.6% | +53.5% | +42.0% |
| + Wall solar (ConvectiveWithFluxToDomain) | +13.1% | -5.7% | +63.5% | +20.6% |
| + Wall solar to air (redirect to zone air) | +11.5% | +5.5% | +82.2% | +67.0% |
| + ConvectiveWithFlux (alpha=85% into domain) | +13.4% | -7.8% | +69.5% | +18.1% |
| + Beam to floor, diffuse area-proportional | +11.4% | -0.5% | +58.6% | +31.8% |

Beam solar distribution: beam (direct) component goes 100% to floor mass, diffuse
component is area-proportional between floor and walls. Reduces wall leakage in
winter (heating improves), but also reduces beneficial leakage in summer (cooling
worsens). Case 600 nearly matched. Net total absolute error improves from 2482 to
2291 kWh (-8%).

---

## 5. Lessons & Rejected Approaches

### 5.1 Bug: TARP Clamped to ISO Combined Coefficient

`interior_convection_h()` applied `h.max(h_min_iso)` where h_min_iso = 1/R_si = 7.69.
ISO R_si = 0.13 is a **combined** (convective + radiative) coefficient; TARP returns
convection-only values (1.5-4.0 W/m²K). The clamp made TARP return 7.69 always.
Fix: removed h_min_iso clamp, keeping only H_MIN = 0.1 as physical minimum.

### 5.2 Bug: Case 600 Fiberglass Over-Insulated 16%

Fiberglass thickness was 0.0776m (RSI 1.94 * k=0.04) instead of correct 0.066m
(RSI 1.65 from BESTEST-GSR reference). This artificially reduced envelope losses.

### 5.3 Rejected: Linearized Radiation Through Area-Weighted MRT

Adding h_rad = 4*eps*sigma*T^3 as separate radiative component routed heat through
area-weighted MRT, but the simplified model cannot properly return energy to air.
Results worsened: 600H +21.6%, 900H +134.6%. Fix: TARP bypasses radiation when
view factors are disabled. Proper view-factor model needed to split conv/rad.

### 5.4 Rejected: Wall Solar to Air

`fvm_wall_solar_to_air = true` redirects 70% of solar directly to zone air.
Bypasses ALL thermal mass buffering. Both heating and cooling dramatically
worsened (900H +82.2%, 900C +67.0%).

### 5.5 Rejected: ConvectiveWithFlux for Interior Wall Solar

Uses sol-air BC (alpha ~85% into domain) instead of 100% domain injection.
Marginally less leakage, but heating worse (+69.5% vs +63.5%). The 15%
surface-retained fraction eventually absorbs into the wall anyway.
ConvectiveWithFluxToDomain remains correct for interior absorbed shortwave.

### 5.6 Rejected: Solar to FVM Walls Without View Factors

Before view factors, depositing solar on FVM wall interiors leaked heat outward
without the compensating benefit of proper interior radiation exchange. Heating
worsened for both cases. Only viable when view factors provide correct h_in.

### 5.7 Key Insight: Phantom Insulation from TARP-Only

TARP convection-only gives h_in ~2.5 W/m²K → R_si = 0.40 m²K/W (vs correct
0.13). This "phantom insulation" accidentally compensated for other model
deficiencies. Adding view factors (h_in ~7 W/m²K, R_si ~0.14) removed this
compensation, revealing the solar distribution as the dominant remaining error.

### 5.8 Key Insight: Energy Conservation Requires h_total in Air Gain

When view-factor radiation is active, the air gain must be computed as
`h_total * (T_surf - t_eff)`, not `h_conv * (T_surf - T_air)`. The radiative
portion nets to zero across all surfaces by reciprocity, but omitting it from
the individual surface-to-air gain calculation breaks energy conservation and
loses heat from the zone balance.

### 5.9 Key Insight: Building Polygon Normal Direction

Building polygons have **outward-facing** normals. View-factor ray casting must
negate them to shoot rays into the zone interior. Missing this causes rays to
escape the zone and hit nothing.
