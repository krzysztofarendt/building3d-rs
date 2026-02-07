# Physics Models in building3d-rs

This document describes the physics equations, assumptions, and algorithms underlying the
acoustics, lighting, and energy simulation modules.

---

## Table of Contents

1. [Shared Infrastructure](#1-shared-infrastructure)
2. [Acoustics](#2-acoustics)
3. [Lighting](#3-lighting)
4. [Energy](#4-energy)
5. [Potential Improvements](#5-potential-improvements)

---

## 1. Shared Infrastructure

### 1.1 Ray-Triangle Intersection (Moller-Trumbore Algorithm)

All ray-based simulations use the Moller-Trumbore algorithm (`src/geom/ray.rs`) for
ray-triangle intersection. Given a ray with origin **O** and direction **D**, and a triangle
with vertices **P0**, **P1**, **P2**:

```
E1 = P1 - P0
E2 = P2 - P0
H  = D x E2
a  = E1 . H

f = 1 / a
S = O - P0
u = f * (S . H)
Q = S x E1
v = f * (D . Q)
t = f * (E2 . Q)
```

The intersection is valid when:

- `|a| > 1e-10` (ray not parallel to triangle)
- `-tol <= u <= 1 + tol`
- `-tol <= v <= 1 + tol` and `u + v <= 1 + tol`
- `t > 1e-6` (intersection in front of ray origin)

A barycentric tolerance `tol = 1e-3` is used to prevent rays from slipping through cracks
between adjacent triangles in a mesh.

### 1.2 Voxel Grid Spatial Acceleration

A uniform hash grid (`src/sim/engine/voxel_grid.rs`) accelerates ray-surface queries:

- Cell size defaults to 0.1 m
- Grid extends 1 cell beyond the scene bounding box
- Each cell stores indices of polygons whose bounding boxes overlap it
- Queries return polygon candidates from 3x3x3 = 27 neighboring cells

### 1.3 Flat Scene

The `FlatScene` (`src/sim/engine/mod.rs`) flattens the building hierarchy into a list of
polygons with:

- Precomputed vertex normals
- Material assignments resolved from the `MaterialLibrary`
- Transparent surface detection (internal interfaces between solids in the same zone)
- Voxel grid for spatial queries

---

## 2. Acoustics

The acoustic simulation uses geometric ray tracing to estimate room impulse responses and
derive standard room acoustic metrics.

### 2.1 Ray Initialization

Rays are emitted from a point source with:

- **Initial energy**: 1.0 per ray (scalar) or `[1.0; 6]` per octave band
- **Direction**: uniformly distributed on the unit sphere via rejection sampling
- **Speed**: 343.0 m/s (speed of sound at 20 C in air)

Optional source directivity patterns scale initial energy by a gain factor:

| Pattern | Gain(theta) | Range |
|---------|-------------|-------|
| Omnidirectional | 1.0 | constant |
| Cardioid | 0.5 * (1 + cos(theta)) | [0, 1] |
| Sub-cardioid | 0.75 + 0.25 * cos(theta) | [0.5, 1] |

where theta is the angle between the ray direction and the source forward axis.

### 2.2 Propagation

Two propagation models are available:

**Fixed time step** (default):
```
P(t + dt) = P(t) + V * dt
```
where dt = 2.5e-5 s (default). Reflections are searched within `1.5 * speed * dt` of the
current position, providing a buffer for time-stepping artifacts.

**Event-driven**: Same position update, but sets reflection search distance to infinity,
forcing exact intersection tests for all rays.

### 2.3 Absorption

#### Surface absorption

At each surface interaction, energy is reduced:

```
E_remaining = E * (1 - alpha)
```

In frequency-dependent mode, this is applied per octave band:

```
E_remaining[b] = E[b] * (1 - alpha[b])    for b in {125, 250, 500, 1k, 2k, 4k} Hz
```

Absorption coefficients for preset materials:

| Material | 125 Hz | 250 Hz | 500 Hz | 1 kHz | 2 kHz | 4 kHz |
|----------|--------|--------|--------|-------|-------|-------|
| Concrete | 0.01 | 0.01 | 0.02 | 0.02 | 0.02 | 0.03 |
| Glass | 0.18 | 0.06 | 0.04 | 0.03 | 0.02 | 0.02 |
| Gypsum | 0.29 | 0.10 | 0.05 | 0.04 | 0.07 | 0.09 |
| Carpet | 0.02 | 0.06 | 0.14 | 0.37 | 0.60 | 0.65 |
| Wood | 0.15 | 0.11 | 0.10 | 0.07 | 0.06 | 0.07 |
| Metal | 0.04 | 0.04 | 0.03 | 0.03 | 0.03 | 0.03 |

#### Air absorption (ISO 9613-1)

Optional atmospheric attenuation at standard conditions (20 C, 50% RH):

| Band (Hz) | 125 | 250 | 500 | 1000 | 2000 | 4000 |
|-----------|------|------|------|------|------|------|
| Attenuation (dB/m) | 0.0003 | 0.0011 | 0.0027 | 0.0067 | 0.024 | 0.084 |

Applied as a linear multiplier over traveled distance d:

```
E[b] *= 10^(-attenuation[b] * d / 10)
```

### 2.4 Reflection

Three reflection models:

**Specular** (mirror-like):
```
V_reflected = V_incident - 2 * (V_incident . N) * N
```

**Diffuse** (Lambertian): Random direction uniformly sampled in the hemisphere opposite
to the incident side, generated via rejection sampling.

**Hybrid**: Blends specular and diffuse using the scattering coefficient s:
```
if random() < s:
    use diffuse reflection
else:
    use specular reflection
```

In frequency-dependent mode, scattering is averaged across bands:
```
s_mean = sum(scattering[b]) / 6
```

### 2.5 Receiver Model

A spherical receiver collects ray energy into time-frequency bins:

- **Radius**: configurable (default 0.1 m)
- **Time resolution**: configurable bin width (seconds)
- **Histogram**: `Vec<[f64; 6]>` for per-band energy, `Vec<f64>` for broadband

A ray contributes energy to the receiver when:
```
|P_ray - P_receiver|^2 <= radius^2
```

Energy is accumulated into the time bin:
```
bin = floor(t / time_resolution)
```

### 2.6 Impulse Response

The impulse response is extracted from the receiver histogram. For audio output, the
energy histogram is resampled to a target sample rate by distributing energy density
uniformly within each bin:

```
energy_density = E / time_resolution
sample_value = energy_density * overlap_duration
```

### 2.7 Schroeder Backward Integration

The energy decay curve is computed by backward integration:

```
EDC(t) = 10 * log10( sum(E[t:]) / sum(E[:]) )
```

The result is a monotonically decreasing curve in dB, normalized to 0 dB at t=0.

### 2.8 Room Acoustic Metrics

All metrics are derived from the Schroeder decay curve or the energy histogram.

**Reverberation time (RT)**:

Linear regression on the decay curve between two dB thresholds, extrapolated to -60 dB:

```
slope = (n * sum(x*y) - sum(x) * sum(y)) / (n * sum(x^2) - sum(x)^2)
RT = -60 / slope   [seconds]
```

| Metric | Start (dB) | End (dB) |
|--------|-----------|----------|
| T20 | -5 | -25 |
| T30 | -5 | -35 |
| EDT | 0 | -10 |

RT60 prefers T30; falls back to T20 if the decay range is insufficient.

**Clarity (C80)**:

Ratio of early to late energy at the 80 ms boundary:

```
C80 = 10 * log10( E[0..80ms] / E[80ms..] )   [dB]
```

**Definition (D50)**:

Fraction of energy arriving within 50 ms:

```
D50 = E[0..50ms] / E_total
```

**Speech Transmission Index (STI)**:

Simplified single-modulation-frequency estimate:

1. Per-band modulation transfer function (f_mod = 1 Hz):
   ```
   m[b] = 1 / sqrt(1 + (2*pi*f_mod * RT60[b] / 13.8)^2)
   ```

2. Apparent signal-to-noise ratio:
   ```
   SNR[b] = 10 * log10(m[b] / (1 - m[b]))     clipped to [-15, 15] dB
   ```

3. Weighted average with male-speech band weights:
   ```
   weights = [0.129, 0.143, 0.114, 0.114, 0.186, 0.171]
   mean_SNR = sum(weight[b] * SNR[b]) / sum(weight)
   ```

4. Final index:
   ```
   STI = (mean_SNR + 15) / 30     clipped to [0, 1]
   ```

### 2.9 Key Assumptions (Acoustics)

- **Energy-based model**: rays carry scalar energy (or per-band intensity), no phase or
  wave interference. Valid in the geometric acoustics regime (room dimensions >> wavelength).
- **No diffraction**: rays either hit or miss surfaces; no edge diffraction modeling.
- **Lambertian diffuse scattering**: uniform hemisphere distribution, no angle-dependent
  BRDF.
- **Uniform material properties per polygon**: no spatial variation within a surface.
- **Omnidirectional receiver**: no listener directivity.
- **Fixed atmospheric conditions**: air absorption assumes 20 C, 50% RH.
- **Simplified STI**: uses single modulation frequency (1 Hz) instead of the full
  14-frequency-per-band IEC 60268-16 method.

---

## 3. Lighting

The lighting simulation uses forward and backward ray tracing to compute illuminance
on building surfaces, with physically-based sky and solar models.

### 3.1 Light Sources

**Point light** (omnidirectional):

```
I(direction) = [I_R, I_G, I_B]     (constant in all directions)
total_flux = (I_R + I_G + I_B) * 4*pi / 3
```

Constructor from total lumens: `I_per_channel = lumens / (4*pi)`.

**Area light** (Lambertian emitter):

```
I(direction) = I_base * |N . direction|
total_flux = pi * I_per_area * width * height
```

where N is the surface normal. The pi factor comes from integrating the cosine distribution
over the hemisphere.

**Directional light** (parallel beam, e.g. sunlight):

```
I(direction) = [E_R, E_G, E_B]     (constant irradiance)
total_flux = infinity
```

Rays originate from random points on the bounding box face perpendicular to the light
direction.

### 3.2 Forward Ray Tracing

For each light source, `num_rays` rays are traced through the scene:

1. **Initialize** energy per ray:
   - Point light: `E[c] = intensity[c] / num_rays`
   - Directional: `E[c] = irradiance[c] / num_rays`

2. **Trace** each ray up to `max_bounces` (default 5):
   - Find closest surface intersection
   - Accumulate incident flux at the hit polygon
   - Apply surface reflectance: `E'[c] = E[c] * reflectance[c]`
   - Russian roulette termination:
     ```
     p_survive = max(reflectance[R], reflectance[G], reflectance[B])
     if random() > p_survive: terminate ray
     else: E[c] = E[c] / p_survive    (energy compensation)
     ```
   - Reflect diffusely (Lambertian) for next bounce
   - Terminate if total energy < 1e-6

3. **Compute illuminance**: `illuminance[surface] = flux[surface] / area[surface]` (lux)

### 3.3 Backward Ray Tracing

Sensor-centric calculation for daylight analysis:

For each sensor on a surface with normal **N**:

1. Cast `num_rays` random directions into the hemisphere above the sensor
2. For unobstructed rays (no surface hit), evaluate sky/light intensity:
   ```
   contribution[c] += intensity[c] * |direction . N|
   ```
3. Scale by hemisphere solid angle:
   ```
   illuminance[c] = contribution[c] * 2*pi / num_rays
   ```

The `|direction . N|` factor is Lambert's cosine law for diffuse surface reception.

### 3.4 Solar Position (Spencer 1971)

Solar position is calculated from latitude, day of year, and hour (solar time).

**Day angle**:
```
gamma = 2*pi * (day_of_year - 1) / 365
```

**Declination** (Spencer's Fourier approximation):
```
delta = 0.006918 - 0.399912*cos(gamma) + 0.070257*sin(gamma)
        - 0.006758*cos(2*gamma) + 0.000907*sin(2*gamma)
        - 0.002697*cos(3*gamma) + 0.00148*sin(3*gamma)
```

**Hour angle**:
```
h = (hour - 12) * 15 deg
```

**Altitude**:
```
sin(alt) = sin(lat)*sin(delta) + cos(lat)*cos(delta)*cos(h)
```

**Azimuth** (from north, clockwise):
```
cos(azi) = (sin(delta)*cos(lat) - cos(delta)*sin(lat)*cos(h)) / cos(alt)
if h > 0: azi = 360 - azi
```

**Direction vector** (North = +Y, East = +X, Up = +Z):
```
x = cos(alt) * sin(azi)
y = cos(alt) * cos(azi)
z = sin(alt)
```

### 3.5 Sky Luminance Models

Three sky models compute luminance L (cd/m^2) as a function of direction:

**CIE Standard Overcast Sky**:
```
L(direction) = L_z * (1 + 2*sin(altitude)) / 3
```
where L_z is the zenith luminance. For directions below the horizon, L = 0.

Zenith luminance from diffuse horizontal illuminance:
```
L_z = E_h * 9 / (7*pi)
```

**CIE Clear Sky**:
```
L = L_z * f(gamma) * g(altitude) / (f(gamma_z) * g(pi/2))
```

Indicatrix (circumsolar brightening):
```
f(gamma) = 0.91 + 10*exp(-3*gamma) + 0.45*cos^2(gamma)
```

Gradation (altitude dependence):
```
g(alt) = 0.274 * (0.45 - sin(alt) + sin^2(alt))
```

where gamma is the angular distance from the sun.

**Perez All-Weather Model**:
```
L = L_z * f(theta, gamma) / f(0, theta_z)
```

with:
```
f(theta, gamma) = (1 + a*exp(b / max(cos(theta), 0.01)))
                * (1 + c*exp(-d*gamma) + e*cos^2(gamma))
```

Parameterized by 5 coefficients `[a, b, c, d, e]` that vary with sky conditions.

### 3.6 Sensor Grid

Sensors are placed on polygon surfaces using a rectangular grid in local 2D coordinates:

1. Build orthonormal basis on the polygon plane (u_axis from first edge, v_axis = normal x u_axis)
2. Project vertices to 2D, find bounding box
3. Place sensors at `(u_min + spacing/2 + i*spacing, v_min + spacing/2 + j*spacing)`
4. Reject sensors outside the polygon boundary

**Daylight factor**:
```
DF = E_interior / E_exterior
```

### 3.7 Key Assumptions (Lighting)

- **Diffuse-only reflections**: no specular or glossy materials in the tracing loop.
- **RGB color model**: three broadband channels, no spectral wavelength dependence.
- **Isotropic surface properties**: reflectance independent of incident angle (no Fresnel).
- **No transparency**: surfaces are fully opaque (except transparent surface detection
  for internal interfaces).
- **Russian roulette is unbiased**: energy compensation preserves expected value.
- **Uniform hemisphere sampling**: no importance sampling toward light sources (increases
  variance for small/distant lights).
- **Solar time only**: longitude is not used for equation-of-time correction.
- **No atmospheric scattering**: sky models provide direct luminance, no volumetric effects.

---

## 4. Energy

The energy simulation calculates annual heating and cooling demand using steady-state or
simplified transient thermal models.

### 4.1 Wall Construction and U-Value

Thermal resistance is calculated per ISO 6946:

```
R_total = R_se + sum(t_i / lambda_i) + R_si
U = 1 / R_total     [W/(m^2*K)]
```

where:
- `t_i` = layer thickness (m)
- `lambda_i` = thermal conductivity (W/(m*K))
- `R_se` = external surface resistance (m^2*K/W)
- `R_si` = internal surface resistance (m^2*K/W)

**Surface resistances (ISO 6946)**:

| Surface type | R_se | R_si |
|-------------|------|------|
| Wall | 0.04 | 0.13 |
| Floor (downward flow) | 0.04 | 0.17 |
| Roof (upward flow) | 0.04 | 0.10 |

**Thermal capacity per unit area**:
```
C = sum(rho_i * c_p_i * t_i)     [J/(m^2*K)]
```

### 4.2 Steady-State Heat Balance

For each zone at each hour:

**Transmission loss**:
```
Q_transmission = sum(U_i * A_i) * (T_indoor - T_outdoor)     [W]
```

**Infiltration loss**:
```
Q_infiltration = rho * c_p * V * ACH / 3600 * (T_indoor - T_outdoor)     [W]
```

where:
- `rho` = 1.2 kg/m^3 (air density at ~20 C)
- `c_p` = 1005 J/(kg*K) (specific heat of air at constant pressure)
- `V` = zone volume (m^3)
- `ACH` = air changes per hour

**Net demand**:
```
Q_net = Q_transmission + Q_infiltration - Q_internal_gains - Q_solar_gains

Heating demand = max(0, Q_net)
Cooling demand = max(0, -Q_net)
```

### 4.3 Transient Thermal Model (1R1C Lumped Model)

A single-node resistor-capacitor model captures thermal mass effects:

**Total conductance**:
```
UA_total = sum(U_i * A_i)     [W/K]
Infiltration_cond = rho * c_p * V * ACH / 3600     [W/K]
```

**Thermal capacity**:
```
C = V * 50000     [J/K]
```

The factor 50 kJ/(m^3*K) represents a medium-weight construction estimate.

**Temperature evolution** (explicit Euler):
```
Q_losses = (UA_total + Infiltration_cond) * (T_zone - T_outdoor)
Q_net = Q_gains + Q_hvac - Q_losses
T_zone(t + dt) = T_zone(t) + (dt / C) * Q_net
```

with dt = 3600 s (1-hour time step).

**Zero-capacity fallback** (instantaneous steady-state):
```
T_zone = T_outdoor + (Q_gains + Q_hvac) / (UA_total + Infiltration_cond)
```

### 4.4 HVAC Ideal Loads

The HVAC model maintains zone temperature between heating and cooling setpoints:

| Parameter | Default |
|-----------|---------|
| Heating setpoint | 20 C |
| Cooling setpoint | 26 C |
| Heating COP | 1.0 (resistive) |
| Cooling COP | 3.0 (heat pump) |

**Required thermal power**:
```
Q_required = C * (T_setpoint - T_free_floating) / dt
```

**Delivered power** (clamped to capacity):
```
Q_delivered = min(|Q_required|, Q_max_capacity)
```

**Electrical power**:
```
P_electric = Q_delivered / COP
```

A deadband exists between the heating and cooling setpoints where no HVAC operates.

### 4.5 Solar Gains (Lighting-Energy Coupling)

The `SolarBridge` converts lighting simulation results to thermal gains:

```
Q_solar = (Phi_R + Phi_G + Phi_B) / eta * SHGC     [W]
```

where:
- `Phi_R, Phi_G, Phi_B` = incident luminous flux on glazing (lm)
- `eta` = luminous efficacy, default 120 lm/W (solar spectrum)
- `SHGC` = solar heat gain coefficient, default 0.6 (double glazing)

Glazing surfaces are identified by pattern matching on path names ("window", "glazing",
"glass").

### 4.6 Weather Data

EPW (EnergyPlus Weather) files provide hourly data for 8760 hours:

- Dry-bulb temperature (C)
- Relative humidity (%)
- Global horizontal radiation (Wh/m^2)
- Direct normal radiation (Wh/m^2)
- Diffuse horizontal radiation (Wh/m^2)
- Wind speed (m/s) and direction (degrees)

A synthetic weather generator is available using:

**Annual temperature profile**:
```
T(day) = T_mean + amplitude * cos(2*pi * (day - 200) / 365)
```

**Daily variation**:
```
T(hour) = T_annual + 3 * cos(2*pi * (hour - 14) / 24)
```

**Solar radiation** (simplified parabolic model):
```
solar_hour = (hour - 12) / 6
factor = max(0, 1 - solar_hour^2)     for hours 7-19
GHR = 800 * factor     [W/m^2]
Direct = 0.6 * GHR
Diffuse = 0.4 * GHR
```

### 4.7 Internal Gains Schedules

Gains are computed from occupancy, equipment, and lighting schedules:

```
Q_gains(t) = q_person * n_occupants * f_occ(t)
           + P_equipment * f_equip(t)
           + P_lighting * f_light(t)
```

**Typical office values** (per 100 m^2):
- Heat per person: 120 W (seated, light work)
- Occupant density: 1 per 10 m^2
- Equipment power density: 15 W/m^2
- Lighting power density: 10 W/m^2

### 4.8 Key Assumptions (Energy)

- **Steady-state or simple transient**: no detailed finite-element thermal modeling or
  multi-zone air flow.
- **Single-node thermal mass**: entire zone lumped into one capacitance; no temperature
  gradients within the zone.
- **Fixed air properties**: density 1.2 kg/m^3, c_p = 1005 J/(kg*K); no humidity or
  pressure dependence.
- **Ideal HVAC**: instantaneous response, perfect setpoint control (within capacity limits).
- **No latent loads**: only sensible heating/cooling; no humidity control.
- **Uniform infiltration**: constant ACH regardless of wind pressure or stack effect.
- **Solar gains via luminous efficacy conversion**: simplified relationship between
  lighting flux and thermal power.
- **No thermal bridging**: U-values are for clear-field construction only.
- **1D heat conduction**: no lateral heat flow within wall layers.

---

## 5. Potential Improvements

### 5.1 Acoustics

- **Edge diffraction**: Add UTD (Uniform Theory of Diffraction) or Biot-Tolstoy-Medwin
  diffraction around edges, important for low frequencies and barriers.
- **Full STI calculation**: Use all 14 modulation frequencies per band per IEC 60268-16.
- **Importance sampling**: Weight initial ray directions toward receivers to reduce variance.
- **Angle-dependent absorption**: Model absorption as a function of incidence angle
  (impedance-based model).
- **Variable atmospheric conditions**: Allow temperature and humidity inputs instead of
  fixed ISO 9613-1 values.
- **Phase-aware modeling**: Implement image-source method for early reflections to capture
  interference effects.
- **Listener directivity**: HRTF-weighted receiver for binaural metrics.

### 5.2 Lighting

- **Specular and glossy reflections**: Add a Cook-Torrance or Phong BRDF for materials
  like polished metal and glass.
- **Importance sampling**: Direct light sampling (next-event estimation) to reduce variance
  for small or distant sources.
- **Spectral rendering**: Replace RGB with spectral wavelength sampling for accurate color
  rendering and metameric analysis.
- **Transparency and refraction**: Model glass transmission with Fresnel equations and
  Snell's law for windows.
- **Equation of time**: Use longitude for proper solar time correction (currently uses
  solar time directly).
- **Participating media**: Volumetric scattering for fog, smoke, or hazy conditions.

### 5.3 Energy

- **Multi-zone air flow**: Implement pressure network or CFD coupling for inter-zone
  air movement.
- **Detailed thermal mass**: Multi-node RC networks (e.g. ISO 13790 5R1C) or finite
  difference through wall layers.
- **Latent loads**: Humidity balance for HVAC sizing and comfort analysis.
- **Thermal bridging**: Linear and point thermal bridge corrections (psi and chi values).
- **Ground coupling**: Ground temperature model for slab-on-grade and basement walls.
- **Variable infiltration**: Pressure-dependent air leakage model (e.g. power law).
- **Detailed HVAC**: Part-load curves, variable COP, multi-stage systems, heat recovery.
- **Radiant exchange**: View-factor-based long-wave radiation between interior surfaces.
- **Moisture transport**: Hygrothermal modeling for condensation risk assessment.
