# Building3D vs EnergyPlus: Physics Comparison

## 1. Zone Air Heat Balance

| Aspect | EnergyPlus | Building3D | Gap |
|--------|-----------|------------|-----|
| **Time integration** | 3rd-order backward difference (default), Euler, or exact analytical | Backward Euler (1st order) | Minor - Euler is stable but less accurate per timestep |
| **Predictor-Corrector** | Yes - predicts load, then corrects zone temp | No - single solve per timestep | Moderate - E+ converges better with large timesteps |
| **Air capacitance** | `V * rho * Cp / dt` with adjustable multiplier (`ZoneVolCapMultpSens`) | `50,000 J/(m³·K)` default or configurable | OK - both adjustable, but B3D default is high |
| **Surface coupling** | Explicit `SumHA + SumHATsurf` per surface with per-surface `h_c` | Via FVM boundary or lumped `U*A` | See convection below |

### EnergyPlus Details

Zone air temperature solved using predictor-corrector with key coefficients:

```
TempDepCoef = SumHA + SumMCp
TempIndCoef = SumIntGain + SumHATsurf - SumHATref + SumMCpT + SysDepZoneLoadsLagged
```

Where:
- `SumHA` = Sum of convection coefficient x area for all surfaces
- `SumHATsurf` = Sum of h_c x A x T_surface (sensible heat from surfaces)
- `SumIntGain` = Sum of convective internal gains (people, lights, equipment)
- `SumMCp` = Sum of mass flow rate x specific heat (infiltration/ventilation)
- `SumMCpT` = Sum of mass flow rate x Cp x inlet air temperature

3rd-order backward difference (default):

```
TempHistoryTerm = AirPowerCap * (3*ZTM[0] - 1.5*ZTM[1] + (1/3)*ZTM[2])
tempDepLoad = (11/6)*AirPowerCap + TempDepCoef
tempIndLoad = TempHistoryTerm + TempIndCoef
T_zone = tempIndLoad / tempDepLoad
```

Air power capacity: `AirPowerCap = V * ZoneVolCapMultpSens * rho * Cp / dt`

### Building3D Details

Backward Euler single-node:

```
T(t+dt) = [T(t) + (dt/C)*(Q_gains + Q_hvac + K*T_out)] / (1 + dt*K/C)
```

Two-node envelope model (2R1C):

```
[1 + dt(K_env+K_am)/Ca    -dt*K_am/Ca  ] [T_a(t+dt)]   [...]
[-dt*K_am/Cm              1 + dt*K_am/Cm] [T_m(t+dt)] = [...]
```

---

## 2. Wall Conduction

| Aspect | EnergyPlus | Building3D | Gap |
|--------|-----------|------------|-----|
| **Method** | CTF (default, fast) or Finite Difference (Crank-Nicolson / Implicit) | FVM with Implicit Euler + Thomas algorithm | Comparable - both solve 1D transient conduction |
| **Discretization** | Variable node spacing, temperature-dependent properties | Fixed max cell size 0.05m, constant properties | Minor - E+ handles variable conductivity |
| **Material properties** | Temperature-dependent k, Cp, rho (in FiniteDiff mode) | Constant per layer | Gap - no temperature-dependent properties in B3D |
| **Phase change** | Supported (enthalpy-temperature function) | Not supported | Minor for most cases |

### EnergyPlus Details

**CTF (Conduction Transfer Functions):** Pre-computed response factors.

```
Q = X0*T_in + X1*T_in_prev + ... + Y0*T_out + Y1*T_out_prev + ...
```

**Finite Difference:** Crank-Nicolson (2nd order) or Fully Implicit (1st order). Variable node spacing with temperature-dependent properties. Stability monitored via Fourier number.

### Building3D Details

**FVM (Finite Volume Method):** Implicit backward Euler with Thomas algorithm for tridiagonal system.

Series-resistance face conductance at material interfaces:

```
K_face = A / (half_dx_L / k_L + half_dx_R / k_R)
```

Boundary convective BC with series resistance:

```
K_eff = 1 / (1/(h*A) + 1/K_face)
```

---

## 3. Convection Coefficients

| Aspect | EnergyPlus | Building3D | Gap |
|--------|-----------|------------|-----|
| **Interior** | Multiple models: TARP (buoyancy), Ceiling Diffuser, Beausoleil-Morrison (mixed), Khalifa, adaptive selection | TARP natural convection (buoyancy-driven, tilt-dependent) or fixed value | Minor - TARP covers most cases; missing ceiling diffuser and mixed convection |
| **Exterior** | DOE-2, TARP, MoWiTT, adaptive; wind speed + direction dependent; windward/leeward distinction | DOE-2 simplified (natural + wind-forced, sqrt combination) or fixed | Moderate - no windward/leeward distinction, no surface roughness |
| **Natural convection** | `h = 1.31 * abs(dT)^(1/3)` (vertical), tilt-dependent (horizontal) | Same TARP/Walton correlations | OK |

### EnergyPlus Interior Convection Details

**TARP/Walton natural convection (buoyancy-driven):**

- Vertical walls (ASHRAE): `h = 1.31 * |dT|^(1/3)`
- Horizontal unstable (heated below): `h = 9.482 * |dT|^(1/3) / (7.238 - |cos(tilt)|)`
- Horizontal stable: `h = 1.810 * |dT|^(1/3) / (1.382 + |cos(tilt)|)`

**Ceiling diffuser (Fisher-Pedersen):** `h = f(ACH)` for different surface orientations.

**Mixed convection (Beausoleil-Morrison):** Combines natural and forced effects for aided/opposing flows with different correlations for walls, floors, ceilings.

**Adaptive algorithm:** Maps surface conditions (orientation, temperature difference, air movement) to appropriate correlation.

### EnergyPlus Exterior Convection Details

**Natural:** Same correlations as interior, applied to exterior surface.

**Forced (wind-driven):**
- DOE-2 model with windward/leeward distinction
- MoWiTT model combining buoyancy + wind
- Surface roughness factors

**Combined:** `h_ext = h_natural + h_forced` or weighted combination.

### Building3D Details

**Interior:** TARP/Walton natural convection model (default) or fixed value (default 3.0 W/m²K).
TARP correlations are the same as EnergyPlus:
- Vertical: `h = 1.31 * |dT|^(1/3)`
- Horizontal unstable: `h = 9.482 * |dT|^(1/3) / (7.238 - |cos(tilt)|)`
- Horizontal stable: `h = 1.810 * |dT|^(1/3) / (1.382 + |cos(tilt)|)`
- Floor minimum: H_MIN = 0.1 W/m²K

When TARP is active, radiative exchange is bypassed (h_conv = h_total = TARP, t_eff = T_air) because the simplified area-weighted MRT model loses energy through incomplete feedback. A proper view-factor model would be needed to split convection and radiation.

**Exterior:** DOE-2 simplified model (default) or fixed (~25 W/m²K from ISO R_se).
- Natural: same TARP correlations as interior
- Forced: `h_f = a + b * V_wind` (a=3.26, b=3.89 for rough surfaces)
- Combined: `h_ext = sqrt(h_n² + h_f²)`

No windward/leeward distinction or surface roughness categories.

---

## 4. Longwave Radiation

| Aspect | EnergyPlus | Building3D | Gap |
|--------|-----------|------------|-----|
| **Interior LW** | Full view-factor matrix, linearized radiation exchange between all surfaces | Area-weighted MRT approximation (optional) | Moderate - B3D approximation is reasonable for simple geometries |
| **Exterior LW** | Sky temperature model, cloud cover, atmospheric radiation | Uses EPW horizontal IR directly, linearized `h_rad = 4*eps*sigma*T^3` | OK - using EPW IR is valid |
| **Emissivity** | Per-surface, per-layer | Configurable per material | OK |

### EnergyPlus Details

**Interior radiation:** Pre-computed view-factor matrix between all zone surfaces. Linearized radiative exchange:

```
Q_rad = h_r * A * (T_surface - T_MRT)
h_r = 4 * sigma * eps * T_mean^3
T_MRT = (Sum F_ij * eps_j * T_j) / (Sum F_ij * eps_j)
```

**Exterior radiation:** Sky temperature modeling with cloud cover effects. Ground temperature-based ground radiation.

### Building3D Details

**Interior radiation (optional):** Area-weighted MRT:

```
T_mrt = Sum(A_i * T_i) / Sum(A_i)
h_conv = h_total * (1 - f_rad)
h_rad = h_total * f_rad
T_eff = (h_conv * T_air + h_rad * T_mrt) / h_total
```

**Exterior radiation:** Uses EPW `horizontal_infrared_radiation` directly with isotropic sky/ground view:

```
sky_view = 0.5 * (1 + n_z)
ground_view = 1 - sky_view
L_in = sky_view * L_sky + ground_view * L_ground
q_lw = eps * (L_in - sigma * T_air^4)
h_rad = 4 * eps * sigma * T_air^3
```

---

## 5. Solar Radiation

| Aspect | EnergyPlus | Building3D | Gap |
|--------|-----------|------------|-----|
| **Shading** | Full 3D shadow calculation with self-shading, neighboring buildings | Optional FlatScene ray casting | Comparable when ray casting enabled |
| **Window angular properties** | Spectral data, polynomial fits, BSDF for complex glazing | Simple SHGC with optional IAM polynomial or ASHRAE formula | Moderate - B3D lacks spectral detail |
| **Solar distribution** | Full tracking: direct beam hits specific interior surfaces, shortwave reflected/absorbed per surface | Fraction-based: configurable `transmitted_solar_to_air_fraction` to air, remainder to internal mass slabs (floor) | Moderate - B3D distributes to floor mass but not geometry-aware |
| **Opaque absorption** | Computed per surface with sol-air temperature concept | Applied as ConvectiveWithFlux BC on FVM walls | OK - similar physics, different implementation |
| **Ground reflectance** | Configurable, view-factor based | `GHI * rho_ground * ground_view` | OK |

### EnergyPlus Details

**Solar distribution:** Full tracking of beam solar through windows onto interior surfaces. Direct beam hits specific surfaces based on geometry. Shortwave reflected and absorbed per surface with multiple bounces.

**Window optical model:** Multi-layer spectral calculation at each incidence angle. Polynomial angular fits (POLYF). BSDF for complex fenestration. Separate visible vs infrared handling.

**Shading:** Full 3D shadow calculation considering self-shading and neighboring buildings.

### Building3D Details

**Glazing transmittance:**

```
Q_direct = DNI * max(0, sun_dir . normal) * area * SHGC * IAM
Q_diffuse = DHI * sky_view * area * SHGC
Q_ground_reflect = GHI * rho_ground * ground_view * area * SHGC
```

**IAM options:**
- Polynomial: `SHGC(theta) = SHGC_0 * Sum(c[i] * cos^i(theta))` (5th order)
- ASHRAE: `IAM = 1 - a * (1/cos(theta) - 1)` (a ~ 0.1)

**Opaque absorption:** Applied as ConvectiveWithFlux BC on FVM walls with fraction entering wall:

```
alpha_wall = K_face / (K_face + h_out*A)
```

**Solar distribution:** Transmitted solar split by `transmitted_solar_to_air_fraction` (default 0.4): fraction goes directly to zone air, remainder distributed to internal mass slabs (floor) as surface heat flux via `ConvectiveWithFluxToDomain` BCs. Provides thermal lag for absorbed solar but is not geometry-aware (does not track beam direction through windows).

---

## 6. Infiltration/Ventilation

| Aspect | EnergyPlus | Building3D | Gap |
|--------|-----------|------------|-----|
| **Models** | ELA (effective leakage area), crack method, flow coefficient, scheduled | Simple ACH model only | Moderate - ACH is adequate for many cases |
| **Wind/stack effects** | `Q = ELA * sqrt(C_s*dT + C_w*V^2)` | Not modeled | Moderate |
| **Density basis** | Outdoor, standard, or indoor | Fixed 1.2 kg/m³ | Minor |

### EnergyPlus Details

**Effective Leakage Area (ELA):**

```
Q = ELA * sqrt(C_stack * dT + C_wind * V_wind^2)
```

**Flow coefficient method:** Pressure-dependent power law. **Crack width method** (DOE-2 style). **Scheduled ventilation:** Constant or wind/temperature controlled.

### Building3D Details

Simple ACH model:

```
Q_infiltration = rho * cp * V * ACH / 3600 * dT
= 1.2 * 1005 * V * ACH / 3600 * dT
~ 0.335 * V * ACH * dT
```

---

## 7. HVAC / Ideal Loads

| Aspect | EnergyPlus | Building3D | Gap |
|--------|-----------|------------|-----|
| **Ideal loads** | Full predictor-corrector with system feedback | Implicit formula with envelope coupling | Similar approach |
| **Capacity limits** | Max heating/cooling power | Supported | OK |
| **COP** | Supported | Supported | OK |
| **Humidity control** | Latent loads, dehumidification | Not modeled | Moderate for humid climates |

### EnergyPlus Details

Predictor-corrector: predicts system load in prediction step, corrects zone temperature in correction step. Supports thermostat scheduling, dual setpoint, humidity control.

### Building3D Details

Implicit formulation:

```
Q_hvac = C*(T_set - T_zone)/dt + (UA + Inf_cond)*(T_set - T_out) - Q_gains
Q_delivered = min(Q_required, max_capacity)
E_electric = Q_delivered / COP
```

Setpoint logic: heat if T_free < T_heat, cool if T_free > T_cool.

---

## 8. Window/Glazing Model

| Aspect | EnergyPlus | Building3D | Gap |
|--------|-----------|------------|-----|
| **Thermal model** | Multi-layer with gap gas convection (Nusselt correlations), per-layer absorption | Simple SHGC-based | **LARGE** - E+ models each glass layer and gas gap |
| **Optical model** | Spectral, angular (polynomial), multi-bounce between layers | Single SHGC value with optional angular modifier | **LARGE** |
| **Frame/divider** | Modeled separately with thermal bridge | Not modeled | Moderate |

### EnergyPlus Details

**Multi-layer window heat balance:** Glass layers + air/gas gap layers. Gap convection modeled via Nusselt number correlations (Rayleigh number determines conduction vs convection regime). Per-layer absorption and heat balance.

**Optical model:** Spectral data at multiple incidence angles. Multi-bounce between layers. Separate visible vs infrared. BSDF for complex glazing (screens, louvers).

### Building3D Details

Single SHGC value applied to entire window assembly. Optional angular modifier (IAM). No multi-layer thermal or optical modeling.

---

## 9. Ground Coupling

| Aspect | EnergyPlus | Building3D | Gap |
|--------|-----------|------------|-----|
| **Model** | Kiva (2D transient FD), F-factor, or ground temp from weather | Fixed temperature BC | Moderate - simple approach is common |

### EnergyPlus Details

**Kiva Foundation Model:** 2D transient finite difference conduction in ground, coupled to surface boundary conditions. Accounts for slab-on-grade and basement effects.

**F-factor Method:** Simplified U-value and F-value approach without transient ground effects.

### Building3D Details

Simple fixed-temperature boundary condition for ground-coupled surfaces. No deep-soil model or 2D ground conduction.

---

## 10. Additional Differences

| Aspect | EnergyPlus | Building3D |
|--------|-----------|------------|
| **Moisture/latent loads** | Full moisture balance, dehumidification | Not modeled |
| **Airflow network** | Multi-zone pressure-driven airflow | Not modeled |
| **Daylighting** | Detailed daylight factor calculation | Separate lighting module (ray tracing) |
| **HVAC systems** | Full component-based HVAC (coils, fans, chillers, boilers) | Ideal loads only |
| **Schedules** | Arbitrary schedules for all parameters | 24h/168h/8760h repeating |

---

## Priority Recommendations to Align with EnergyPlus

### IMPLEMENTED

#### 1. Dynamic Interior Convection Coefficients -- DONE

TARP/Walton natural convection model implemented in `src/sim/energy/convection.rs`.
Opt-in via `InteriorConvectionModel::Tarp`. Applied to FVM walls, internal mass slabs,
and steady-state surfaces. TARP bypasses radiative exchange (simplified MRT loses energy).

#### 2. Dynamic Exterior Convection Coefficients -- DONE

DOE-2 simplified model implemented in `src/sim/energy/convection.rs`.
Opt-in via `ExteriorConvectionModel::Doe2`. Combines natural (TARP) + forced
(wind-driven) via `sqrt(h_n² + h_f²)`. Wind speed read from EPW via Bus.

#### 3. Interior Solar Distribution -- DONE

Transmitted solar split by `transmitted_solar_to_air_fraction` (default 0.0, recommended 0.4).
Remainder distributed to internal mass slabs as surface heat flux. Provides thermal lag
for absorbed solar via FVM slab conduction.

### MEDIUM PRIORITY

#### 4. Higher-Order Time Integration

Implement 3rd-order backward difference for zone air temperature:

```
TempHistoryTerm = AirPowerCap * (3*ZTM[0] - 1.5*ZTM[1] + (1/3)*ZTM[2])
T_zone = (TempHistoryTerm + TempIndCoef) / ((11/6)*AirPowerCap + TempDepCoef)
```

This improves accuracy at hourly timesteps without needing substeps.

#### 5. Predictor-Corrector for Zone Air

Two-step approach per timestep:
1. **Predict:** Estimate system load needed to reach setpoint
2. **Correct:** Apply actual system output and compute resulting zone temperature

Better convergence especially with large timesteps.

#### 6. Multi-Layer Window Thermal Model

Model U-value from glass + gap + glass with gas fill properties:

- Glass conduction: `k_glass / thickness`
- Gap convection: Nusselt correlations based on Rayleigh number
- Surface radiation in gap: linearized exchange between glass surfaces

At minimum, compute window U-value from layer properties rather than using a single prescribed value.

#### 7. Wind/Stack Infiltration

Replace pure ACH with physics-based model:

```
Q = ELA * sqrt(C_stack * |T_in - T_out| + C_wind * V_wind^2)
```

Where:
- `ELA` = effective leakage area (m²)
- `C_stack` = stack coefficient (depends on building height)
- `C_wind` = wind coefficient (depends on shielding)

### LOWER PRIORITY

#### 8. Temperature-Dependent Material Properties

Allow k and Cp to vary with temperature in FVM solver. Most important for PCM (phase change materials) but also affects accuracy for large temperature swings.

#### 9. Full View-Factor Radiation Exchange

Replace area-weighted MRT with computed view factors between zone surfaces. Current MRT approximation is reasonable for convex rooms but breaks down for L-shaped or complex geometries.

#### 10. Latent Loads / Moisture

Model moisture balance for humidity-sensitive applications. Includes dehumidification loads in HVAC sizing.
