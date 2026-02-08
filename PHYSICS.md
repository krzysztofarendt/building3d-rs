# Physics Models in building3d-rs

This document describes the physics equations, assumptions, and algorithms underlying the
acoustics, lighting, and energy simulation modules.

---

## Table of Contents

1. [Shared Infrastructure](#1-shared-infrastructure)
2. [Acoustics](#2-acoustics)
3. [Lighting](#3-lighting)
4. [Energy](#4-energy)
5. [Known Issues and Fix Plan](#5-known-issues-and-fix-plan)
6. [Potential Improvements](#6-potential-improvements)

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

> **Fixed**: Tolerance reduced from `1e-3` to `1e-6`. See [5.1](#51-shared-infrastructure).

### 1.2 Voxel Grid Spatial Acceleration

A uniform hash grid (`src/sim/engine/voxel_grid.rs`) accelerates ray-surface queries:

- Cell size defaults to 0.1 m
- Grid extends 1 cell beyond the scene bounding box
- Each cell stores indices of polygons whose bounding boxes overlap it
- Two query modes:
  - `find_nearby(pos)`: returns candidates from 3x3x3 = 27 neighboring cells (fast, local)
  - `find_along_ray(origin, direction, max_distance)`: Amanatides-Woo 3D-DDA ray marching
    that steps through all grid cells along the ray path, collecting polygon candidates
    from each visited cell. Used by `find_target_surface()` so that rays originating far
    from geometry can still find intersections.

### 1.3 Flat Scene

The `FlatScene` (`src/sim/engine/mod.rs`) flattens the building hierarchy into a list of
polygons with:

- Precomputed vertex normals
- Material assignments resolved from the `MaterialLibrary`
- Transparent surface detection (internal interfaces between solids in the same zone)
- Voxel grid for spatial queries

### 1.4 Composable Simulation Pipeline (Bus + Modules)

Multi-physics workflows (lighting ↔ thermal, acoustics-only runs, etc.) are intended to be
composed from small modules rather than hard-wired into a monolithic “simulation app”.

The shared runtime lives in `src/sim/framework/`:
- `SimContext`: immutable inputs shared by modules (the `Building` plus `SurfaceIndex`)
- `Bus`: typed message/value store connecting modules (keyed by concrete Rust type)
- `Pipeline`: ordered list of `SimModule`s with `init()` and `step()`

Cross-module “contracts” (payload types carried on the `Bus`) live in `src/sim/coupling.rs`.
Examples:
- `OutdoorAirTemperatureC` (weather boundary for step-based thermal)
- `ShortwaveTransmittedWPerZone` (solar shortwave gains per zone)
- `InternalGainsWPerZone` / `InternalGainsWTotal`

This enables step-based pipelines such as:
1) a weather/solar producer publishes `OutdoorAirTemperatureC` and `ShortwaveTransmittedWPerZone`
2) the thermal module consumes those inputs each step (see `sim::energy::module::EnergyModule`)

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

**Diffuse** (Lambertian): Cosine-weighted random direction in the hemisphere opposite to
the incident side, generated via Malley's method (uniform disk projected onto hemisphere).

> **Fixed**: Now uses Malley's method for cosine-weighted Lambertian sampling.
> See [5.2](#52-acoustics).

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

> **Fixed**: Normalization methods added to `Receiver` that convert raw energy to
> energy density (W/m^2), independent of receiver size and ray count.
> See [5.2](#52-acoustics).

### 2.6 Impulse Response

The impulse response is extracted from the receiver histogram. For audio output, the
energy histogram is resampled to a target sample rate by distributing energy density
uniformly within each bin:

```
energy_density = E / time_resolution
sample_value = energy_density * overlap_duration
```

> **Note**: This produces an *energy* time series (proportional to pressure^2), not
> a signed pressure impulse response. It cannot be used for convolution reverb. For metrics
> (RT, C80, D50) this is acceptable. The documentation now clearly states this.
> See [5.2](#52-acoustics).

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
- **Cosine-weighted (Lambertian) diffuse scattering**: implemented via Malley's method.
- **Uniform material properties per polygon**: no spatial variation within a surface.
- **Omnidirectional receiver**: no listener directivity.
- **Fixed atmospheric conditions**: air absorption assumes 20 C, 50% RH.
- **Simplified STI-like index**: uses single modulation frequency (1 Hz) instead of the
  full 14-frequency-per-band IEC 60268-16 method. Uses 6 bands (125-4k Hz); many STI
  variants include 7 bands up to 8 kHz. Should not be labeled "STI" without qualification.

---

## 3. Lighting

The lighting simulation uses forward and backward ray tracing to compute illuminance
on building surfaces, with physically-based sky and solar models.

### 3.1 Light Sources

**Point light** (omnidirectional):

```
I(direction) = [I_R, I_G, I_B]     (constant in all directions, W/sr per channel)
total_flux = (I_R + I_G + I_B) * 4*pi     [W]
```

Constructor from total radiant flux: `I_per_channel = flux / (3 * 4*pi)`, so that
`total_flux = 3 * I_per_channel * 4*pi = flux`.

**Area light** (Lambertian emitter):

```
I(direction) = I_base * |N . direction|     [W/(sr*m^2) per channel]
total_flux = (I_R + I_G + I_B) * pi * area     [W]
```

where N is the surface normal. The pi factor comes from integrating the cosine distribution
over the hemisphere.

**Directional light** (parallel beam, e.g. sunlight):

```
I(direction) = [E_R, E_G, E_B]     (constant irradiance, W/m^2 per channel)
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

3. **Compute irradiance**: `irradiance[surface] = flux[surface] / area[surface]` (W/m^2)

The 1/r^2 falloff for point lights emerges naturally from the Monte Carlo sampling:
rays are uniformly distributed on the sphere, so a surface at distance r subtends solid
angle `A*cos(theta)/r^2`, and the number of rays hitting it is proportional to that.
This has been verified with a dedicated test.

RGB channels use radiometric units (W/sr per channel for sources, W/m^2 per channel for
irradiance) and are treated as independent spectral bands.

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

The sky model returns radiance (W/(sr*m^2)). The backward tracer variable is named
`radiance` accordingly. The Monte Carlo estimator
`E = (2*pi/N) * sum(L * cos(theta))` is correct for uniform hemisphere sampling
(pdf = 1/(2*pi)).

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

> **Fixed**: Now uses `atan2(sin_azi, cos_azi)` for numerically robust azimuth calculation.
> See [5.3](#53-lighting).

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

### 3.8 Next Steps (Lighting Engine Roadmap)

This project’s long-term goal is to be a **composable simulation core** (acoustics, heat,
light) with **no hard dependency on a specific GUI**. Visualization is handled via Rerun,
and the user-facing interface is expected to be **natural language through an AI agent**
communicating with a small, deterministic MCP tool API.

For the lighting engine specifically, the next steps should be organized around two
principles:

1. **Correctness first, then performance**: establish a validation harness against a
   trusted reference (Radiance-like workflows) before accelerating or adding features.
2. **Composable building blocks**: keep the lighting engine assembled from small pieces
   (scene, materials/BSDF, sky/sun, sampler, integrator, outputs) so an AI agent can
   “compose” a simulation for a specific use-case.

#### 3.8.1 Decide the “first product” (deliverables)

Avoid trying to land everything at once. Pick one of these as the first-class output and
design the integrator around it:

- **Daylight metrics mode (recommended first)**: sensor-based irradiance/illuminance on
  workplanes + annual metrics (DF, DA/sDA/UDI/ASE) driven by sky + sun time series.
- **View/image mode**: luminance/HDR rendering from a camera, for visual evaluation and
  glare analysis.

Both can share almost all infrastructure, but the sampling strategy, outputs, and “what is
considered correct” differ enough that one should lead.

#### 3.8.2 Create a reference-backed validation harness (Radiance as oracle)

Before changing integrators, create a repeatable way to compare against Radiance-style
results. The goal is not a byte-identical match; it is to detect regressions and reduce
unknown bias.

Step-by-step:

1. **Canonical scenes** (small set, hand-auditable):
   - empty box + point light (already partially validated via inverse-square test)
   - box + window + overcast sky (daylight factor sanity)
   - box + window + sun (hard shadows and direct component)
   - simple “L-room” (multiple bounces and occlusion)
   - glass pane variants (clear vs diffuse) once transmissive materials exist
2. **Scene export** from `FlatScene`:
   - geometry: triangulated polygons (consistent winding and normals)
   - materials: map `MaterialLibrary` optical properties to a minimal Radiance material set
   - sky/sun: export CIE sky as an environment source; export sun as a directional source
3. **Golden outputs**:
   - store reference sensor values (and/or images) in `validation/` with configuration
     metadata (ray counts, bounces, random seed, sky parameters)
4. **Tolerance policy**:
   - define acceptable error bands per deliverable (e.g., relative error on mean sensor
     illuminance, percentile error across sensors, and “shadow edge” checks for sun cases)
5. **Regression tests**:
   - add non-flaky tests that run fast locally (small ray counts) and a slower “validation”
     suite for CI/manual runs (large ray counts).

This harness becomes the backbone of future refactors (sampling, BVH, glazing, MIS).

#### 3.8.3 Make the lighting pipeline explicitly modular (composability)

The code already has `LightingSimulation` (forward) and `BackwardTracer`. To move toward a
Radiance-like model while remaining composable, restructure conceptually into:

- **Scene**: intersection, normals, area, and material lookup (currently `FlatScene`)
- **Emitter model**: point/directional/area lights + sky + sun (treat sky/sun as emitters)
- **Surface model (BSDF)**: mapping from `OpticalMaterial` → reflect/transmit sampling
- **Sampler**: deterministic random/low-discrepancy sequences (reproducible per run)
- **Integrator**: a small set of algorithms (direct-only, path tracing, bidirectional later)
- **Measurements**: sensors and/or camera film; post-processing to metrics

Key constraint for agentic use: each component should have a clear config struct that can
be serialized, diffed, and re-run deterministically.

#### 3.8.4 Implement a Radiance-style backward integrator (incremental milestones)

Radiance’s strength in building daylighting comes from **backward sampling from sensors**
(and images), plus aggressive variance reduction. Implement this as a staged sequence:

1. **Direct-only backward estimator (sky + sun)**:
   - keep the current hemisphere sampling, but add explicit evaluation of the chosen sky
     and sun models
   - add a “separate sun” option (sun as delta light) so the estimator handles sharp shadows
2. **Next-event estimation (NEE)**:
   - for each bounce, sample direct contribution from explicit emitters:
     - sun direction (delta) if visible
     - one or more sky patches / environment samples with proper pdf
   - keep the existing path continuation for indirect light
3. **Multiple bounces with energy conservation**:
   - ensure diffuse bounce sampling is cosine-weighted (Lambertian BSDF)
   - keep Russian roulette, but tie termination probability to path throughput
4. **Multiple importance sampling (MIS)**:
   - combine BSDF sampling and light sampling (especially important for sun/bright sky)
   - implement power heuristic weights
5. **Specular/glossy reflection and transmissive glazing**:
   - specular reflection: deterministic reflection direction
   - refraction through glazing: direction change + transmittance/absorption model
   - “thin glass” approximation as a first step (no thickness), then optionally thickness
6. **BSDF hooks**:
   - keep a trait boundary so a future measured/window-system BSDF can plug in without
     changing the integrator’s core logic.

At each milestone, add/extend one canonical validation scene so correctness stays anchored.

#### 3.8.5 Sky, sun, and climate integration (annual workflows)

For building simulation, the biggest value is often **annual daylight availability** and
coupling to thermal loads. A practical progression:

1. **Keep CIE skies for unit tests** (stable and parameter-light).
2. **Add a weather-driven sky/sun pipeline**:
   - reuse EPW parsing from `src/sim/energy/weather.rs` to get per-hour DNI/DHI/GHI
   - compute sun position per timestamp (already present in `src/sim/lighting/solar.rs`)
   - generate an environment distribution (sky patches) driven by DHI and sun from DNI
3. **Introduce luminous efficacy boundaries**:
   - clearly separate *radiometric* engine units from *photometric* reporting (lux)
   - convert only at outputs, with documented assumptions and user-overridable constants
4. **Thermal coupling**:
   - use lighting results as a potential source for solar gains, but prefer a single
     authoritative pipeline (avoid “two ways to compute the same gains” drifting apart)
   - if a lighting→thermal bridge remains, document it as an approximation with limits.

#### 3.8.6 Outputs, metrics, and post-processing (what the engine returns)

To keep the core composable, treat “metrics” as a layer on top of a small set of primary
outputs:

- **Primary outputs**:
  - per-sensor irradiance/illuminance (direct + indirect separated if possible)
  - per-surface irradiance maps (for visualization and debugging)
  - optional path statistics: bounce counts, visibility ratios, variance estimates
- **Derived daylight metrics** (computed from time series + thresholds):
  - daylight factor (single sky condition)
  - DA/sDA, UDI, ASE (annual)
  - glare metrics (later; typically needs view/camera mode and luminance distributions)

Prefer returning both:
1) machine-friendly data (tables), and
2) a small set of Rerun artifacts (heatmaps, falsecolor, debug rays) for inspection.

#### 3.8.6a Interface guidelines for thermal coupling (keep worktrees aligned)

To keep lighting results usable as inputs to thermal simulation (without hardcoding any
thermal metadata into geometry), adopt the following conventions:

- **Stable identifiers**: any per-surface output intended for other modules should be keyed
  by polygon `UID` (with optional `zone/solid/wall/polygon` path strings for reporting).
- **Payload types (code)**: `sim::coupling::ShortwaveAbsorbedWPerPolygon` and
  `sim::coupling::ShortwaveTransmittedWPerZone` define the default cross-module contracts.
- **Producers**: choose exactly one shortwave producer in a composed pipeline:
  - deterministic EPW-driven producer: `sim::lighting::shortwave::SolarShortwaveModule`, or
  - ray-based producer: `sim::lighting::shortwave::LightingToShortwaveModule` fed by a lighting run.
- **Step-based pipelines**: for hour-by-hour composition, use `sim::lighting::shortwave::SolarShortwaveStepModule`
  to publish `OutdoorAirTemperatureC` and `ShortwaveTransmittedWPerZone` each step, then consume them in
  `sim::energy::module::EnergyModule`.
- **Bus inputs**: step-based thermal simulations consume weather and gains from the `Bus`:
  - `sim::coupling::OutdoorAirTemperatureC`
  - `sim::coupling::InternalGainsWPerZone` (preferred) or `sim::coupling::InternalGainsWTotal` (fallback)
  - `sim::coupling::ShortwaveTransmittedWPerZone`
- **Separation of concerns**: `sim::lighting::module::LightingModule` publishes `LightingResult`
  only; shortwave coupling payloads are produced explicitly by the chosen producer module.
- **Units**: keep the integrator in radiometric units (W, W/m², W/sr) and convert to
  photometric units (lux, cd/m²) only at output/reporting boundaries.
- **Single source of truth for shortwave gains**: in a composed simulation pipeline, do not
  compute solar gains in two different places. Either:
  - a lighting/solar module produces absorbed/transmitted shortwave power, or
  - the thermal module uses EPW + SHGC as a fallback approximation.

#### 3.8.7 Performance roadmap (after correctness)

Once validated, address performance in a way that preserves determinism:

1. **Better acceleration**:
   - keep `VoxelGrid` for “broadphase,” but add a BVH for triangle intersections (narrowphase)
   - consider per-solid/per-zone BVHs to align with the hierarchy and allow incremental updates
2. **Sampling efficiency**:
   - importance sample the sky (bright patches) and the sun
   - stratified / low-discrepancy sequences per sensor pixel (reduces noise at same ray count)
3. **Parallelism**:
   - parallelize over sensors (and/or paths) with stable chunking to keep reproducibility
4. **Caching** (optional, but important for agentic iteration):
   - cache `FlatScene` and derived structures (BVH, patch distributions) keyed by config hash

#### 3.8.8 MCP + agentic “simulation composition” interface

The lighting engine should be controllable via a small MCP tool surface that an AI agent can
call safely and deterministically.

Recommended design steps:

1. **Define a “case spec”** (serializable struct) containing:
   - building/model reference (or procedural instructions to generate one)
   - material assignments (via `MaterialLibrary` patterns and/or explicit overrides)
   - sky/sun/weather configuration (including timestamp series for annual mode)
   - sensor definition (workplane height, grids, polygon attachments)
   - integrator selection (direct-only vs path tracing), ray counts, bounces, seeds, budgets
2. **Expose idempotent MCP tools**:
   - create/update/delete model entities and materials
   - run lighting simulation as an asynchronous job with explicit resource limits
   - fetch results (tables + references to generated Rerun logs/artifacts)
3. **Make every run explainable**:
   - return the full resolved configuration (after defaults) and a minimal “diff” from
     the previous run so an agent can justify changes in natural language
4. **Keep the core GUI-free**:
   - never return “UI actions” from the engine; return data + optional Rerun recordings only.

This mirrors the composability goal: the agent composes the *case spec*, the core runs it,
and Rerun is just one output consumer.

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

Where the sum is taken over **exterior** envelope surfaces only (see 4.2.1).

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

#### 4.2.1 Envelope classification and multi-zone coupling (steady-state)

Polygons are classified using a facing-graph overlay (keyed by polygon `UID`) into:
- **exterior** (faces “nothing”),
- **same-zone interface** (internal partitions between solids in the same zone; excluded from envelope),
- **inter-zone interface** (partition between two zones; becomes a coupling term).

This classification is imposed *at simulation time* and is not stored on geometry types.

For two zones `i` and `j` with an inter-zone partition conductance `K_ij` (W/K):
```
Q_ij = K_ij * (T_i - T_j)   [W]
```

Zone-level steady-state balance with coupling:
```
0 = -K_out,i*(T_i - T_out) - sum_j K_ij*(T_i - T_j) + Q_gains,i + Q_hvac,i
```

The conductance `K_ij` is computed from geometric overlap area and an explicit policy for
combining the two assigned U-values (one per facing polygon). Two common policies:

- **Mean** (default): `U_eq = (U1 + U2)/2`
- **Series**: `U_eq = 1 / (1/U1 + 1/U2)`

Then:
```
K_ij = U_eq * A_overlap
```

In code this policy is controlled by `ThermalConfig::interzone_u_value_policy`
(`sim::energy::config::InterZoneUValuePolicy`).

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

> **Note**: This is a tuning parameter with large influence on dynamics. Typical values
> range from ~30 kJ/(m^3*K) (lightweight) to ~80 kJ/(m^3*K) (heavyweight). Should be
> documented as such and ideally derived from actual construction layers.

**Temperature evolution** (backward Euler / implicit):
```
K = UA_total + Infiltration_cond
T_zone(t + dt) = (T_zone(t) + dt/C * (Q_gains + Q_hvac + K*T_outdoor)) / (1 + dt*K/C)
```

with dt = 3600 s (1-hour time step).

**Zero-capacity fallback** (instantaneous steady-state):
```
T_zone = T_outdoor + (Q_gains + Q_hvac) / (UA_total + Infiltration_cond)
```

#### 4.3.1 Multi-zone transient air-node model (network solve)

For `N` zones, represent each zone air temperature as a state `T_i`. Let:
- `C_i` be zone thermal capacity (J/K), typically `C_i = V_i * C_vol`,
- `K_out,i = UA_i + K_inf,i` (W/K),
- `K_ij` be inter-zone conductance (W/K),
- `Q_gains,i` be total gains in zone `i` (W),
- `Q_hvac,i` be HVAC heat input to zone `i` (W, positive heating, negative cooling).

Backward Euler discretization yields a linear system each timestep:
```
C_i/dt * (T_i^{n+1} - T_i^n)
  = -K_out,i * (T_i^{n+1} - T_out)
    - sum_j K_ij * (T_i^{n+1} - T_j^{n+1})
    + Q_gains,i + Q_hvac,i
```

In “ideal loads” mode, some zones may be clamped to heating/cooling setpoints, turning those
temperatures into fixed boundary conditions for the solve.

Implementation notes (code):
- Multi-zone model: `MultiZoneAirModel` (`src/sim/energy/network/multizone.rs`)
- Annual runners: `run_multizone_transient_simulation()` and `run_multizone_steady_simulation()`
  (`src/sim/energy/simulation.rs`)
- Per-zone default solar: `compute_solar_gains_per_zone()` (`src/sim/energy/solar_bridge.rs`)

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

> **Fixed**: `calculate_with_losses()` uses the implicit formulation that accounts for
> concurrent envelope losses. The transient simulation now uses this method.
> See [5.4](#54-energy).

**Delivered power** (clamped to capacity):
```
Q_delivered = min(|Q_required|, Q_max_capacity)
```

**Electrical power**:
```
P_electric = Q_delivered / COP
```

A deadband exists between the heating and cooling setpoints where no HVAC operates.

### 4.5 Solar Gains

Solar gains are computed directly from EPW weather radiation data combined with sun
geometry and window properties (`SolarGainConfig`):

```
Q_solar_window = (DNI * cos(theta_inc) + DHI * sky_view) * A_window * SHGC     [W]
```

where:
- `DNI` = direct normal irradiance from weather record (W/m^2)
- `DHI` = diffuse horizontal irradiance from weather record (W/m^2)
- `cos(theta_inc) = max(0, sun_direction . polygon_normal)` (only sun-facing surfaces)
- `sky_view = 0.5 * (1 + max(0, normal_z))` (isotropic sky view factor)
- `A_window` = glazing polygon area (m^2)
- `SHGC` = solar heat gain coefficient, default 0.6 (double glazing)

Solar position is computed per hour using the Spencer algorithm (see Section 3.4).
When the sun is below the horizon, only the diffuse component contributes.

Glazing surfaces are identified by pattern matching on path names ("window", "glazing",
"glass") or by explicit per-surface SHGC entries in `SolarGainConfig`.

For multi-zone simulations, the same calculation can be aggregated **per zone** by summing
glazing gains within each zone:
```
Q_solar,zone = sum_{glazing in zone} Q_solar_window
```

A legacy `SolarBridgeConfig` and `lighting_to_solar_gains()` function remain available
for converting lighting simulation results to thermal gains via luminous efficacy, but
the recommended approach is the direct weather-based calculation.

> **Note**: Window detection by name pattern is fragile. Should use an explicit material
> flag or construction type property.

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
  pressure-based multi-zone air flow. Multi-zone heat coupling through partitions is supported.
- **Single-node thermal mass**: entire zone lumped into one capacitance; no temperature
  gradients within the zone.
- **Fixed air properties**: density 1.2 kg/m^3, c_p = 1005 J/(kg*K); no humidity or
  pressure dependence.
- **Ideal HVAC**: instantaneous response, perfect setpoint control (within capacity limits).
- **No latent loads**: only sensible heating/cooling; no humidity control.
- **Uniform infiltration**: constant ACH regardless of wind pressure or stack effect.
- **Isotropic sky model for diffuse solar**: diffuse component uses a simple tilt factor
  rather than an anisotropic sky model (e.g. Perez).
- **No thermal bridging**: U-values are for clear-field construction only.
- **1D heat conduction**: no lateral heat flow within wall layers.

### 4.9 Next Steps (Thermal Simulation Roadmap)

This section is an **interface document**: it describes the intended model boundaries and
data contracts so that different worktrees (e.g. lighting-focused vs energy-focused) can
advance independently without drifting.

The guiding principle is the same as for lighting (Section 3.8): **composable building
blocks** with **deterministic configs**, where domain-specific semantics are imposed as
late-stage *overlays* keyed by stable IDs (UIDs), not stored directly on geometry types.

#### 4.9.1 First milestone: correct envelope accounting + multi-zone readiness

1. **Boundary classification overlay**
   - Classify each polygon as:
     - exterior (to weather boundary),
     - interface to same zone (ignore for transmission),
     - inter-zone partition (coupled to another zone air node),
     - (later) ground / adiabatic / fixed-temperature boundary.
   - This must remain an overlay keyed by polygon `UID` (not a `Polygon` field).

2. **Stop double-counting internal partitions**
   - Exterior transmission should sum only `U*A` for exterior polygons.
   - Inter-zone partitions should not appear in “to-outdoor” transmission; they become
     coupling conductances between zone air nodes.

#### 4.9.2 Thermal network core (EnergyPlus-like heat balance kernel)

Move from “ad-hoc formulas” to a reusable **thermal network representation**:

- **Nodes** (initially): one air temperature node per `Zone` (a deliberate simplification).
- **Components**:
  - exterior conductances `UA_zone = sum(U*A)` (W/K),
  - inter-zone conductances `K_zone_a_zone_b` (W/K),
  - infiltration/ventilation conductance (W/K),
  - heat sources: internal + solar + HVAC (W).

Inter-zone conductances require an explicit policy for combining two facing U-values. This
must be configurable so different modeling interpretations can be used without changing
geometry or solver code (see 4.2.1).

This network should be able to run in:
- **steady-state** mode (instantaneous loads), and
- **transient** mode (implicit timestep update; stable at 1-hour steps).

#### 4.9.3 Pluggable wall/partition models (ISO 13790 → FD → 3D)

Keep the network API stable while swapping internal component models:

1. **Steady U-value** (current baseline): partitions and envelope are pure conductances.
2. **Low-order RC networks** (recommended next): 2R1C/3R2C per construction for dynamic
   surface temperatures and thermal mass effects.
3. **1D finite-difference through layers** (optional): higher fidelity transient conduction.
4. **3D conduction** (future): couple to tetrahedral meshes / FEM; should still expose
   the same boundary heat-flow interface to the zone network.

#### 4.9.4 Comfort outputs and radiant exchange (incremental)

Once interior surface temperatures exist (via RC/FD models), add:
- **MRT / operative temperature** per zone (start with area-weighted approximation),
- **long-wave radiant exchange** (view-factor-based, later; can reuse ray infrastructure).

#### 4.9.5 HVAC and airflow as replaceable boundary components (Modelica-style boundary)

Treat HVAC and airflow as components connected to the zone air node:
- keep **ideal loads** as the default actuator,
- add optional explicit systems (fan-coil, heat pump with part-load curves, ERV/HRV),
- add optional **airflow coupling** (pressure network) without changing geometry.

#### 4.9.6 Cross-domain coupling contract (lighting ↔ thermal)

To prevent “two ways to compute the same gains” drifting apart:

- In any composed simulation, **exactly one module** should be responsible for producing
  *shortwave solar absorption*.
- Downstream thermal modules should consume a common payload keyed by polygon `UID`
  (with optional path strings only for reporting/debugging).

Suggested payload (conceptual):
- `ShortwaveAbsorbedWPerPolygon { polygon_uid -> watts }`
- `ShortwaveTransmittedWPerZone { zone_uid -> watts }` (optional simplification)

In code, these contracts live in `sim::coupling` as `ShortwaveAbsorbedWPerPolygon` and
`ShortwaveTransmittedWPerZone`.

For multi-zone thermal models, `ShortwaveTransmittedWPerZone` is the most direct input: it
maps cleanly onto the per-zone gains vector used by the zone-air solver. Per-polygon absorbed
shortwave is optional (useful for future surface-temperature / radiant models).

In step-based composed simulations, `sim::energy::module::EnergyModule` consumes these payloads
along with:
- `sim::coupling::OutdoorAirTemperatureC`
- `sim::coupling::InternalGainsWPerZone` (or `sim::coupling::InternalGainsWTotal`)

Thermal should support a fallback path (EPW + SHGC) when no lighting/solar module is present,
but the composed pipeline should designate a single authoritative producer.

#### 4.9.7 Determinism requirements (for agentic iteration)

All thermal simulation configs should be:
- serializable and diffable,
- reproducible given a seed and a fixed case spec,
- explicit about units (SI internally).

Where randomness is introduced (e.g. Monte Carlo view factors), it must be seeded and the
results should be cacheable by a config hash.

---

## 5. Known Issues and Fix Plan

Issues are grouped by severity: **bugs** (produce incorrect results and must be fixed),
**biases** (produce results that depend on numerical parameters in non-physical ways), and
**labeling** (code works but documentation or naming is misleading).

### 5.1 Shared Infrastructure

#### ~~BUG: Voxel grid misses distant intersections~~ (FIXED)

**Fixed**: `find_along_ray()` implements Amanatides-Woo 3D-DDA ray marching that steps
through grid cells along the ray direction. `find_target_surface()` now uses this instead
of `find_nearby()`, so rays originating far from geometry find intersections correctly.

#### ~~BIAS: Barycentric tolerance too large~~ (FIXED)

**Fixed**: Reduced `BARY_TOLERANCE` from `1e-3` to `1e-6` in `src/sim/engine/mod.rs`.
This prevents energy leaks from double-counted hits and wrong normals near triangle edges.

### 5.2 Acoustics

#### ~~BIAS: Diffuse reflection is uniform, not Lambertian~~ (FIXED)

**Fixed**: Replaced uniform hemisphere rejection sampling with Malley's method in the
`Diffuse` reflection model (`src/sim/engine/reflection.rs`). The new implementation
samples a uniform disk and projects onto the hemisphere, producing a cosine-weighted
distribution (`pdf = cos(theta) / pi`). Verified by a test checking that the mean
`cos(theta)` equals 2/3.

#### ~~BIAS: Receiver energy depends on radius and timestep~~ (FIXED)

**Fixed**: Added `normalization_factor(num_rays)`, `normalized_scalar_histogram(num_rays)`,
and `normalized_histogram(num_rays)` methods to `Receiver` (`src/sim/acoustics/receiver.rs`).
The normalization factor `4*pi / (pi * r^2 * N_rays)` converts raw collected energy to
energy density (W/m^2) independent of receiver size and ray count. Raw recording is
unchanged; normalization is applied as a post-processing step.

#### ~~LABELING: Impulse response is energy, not pressure~~ (FIXED)

**Fixed**: Updated `ImpulseResponse` documentation in `src/sim/acoustics/impulse_response.rs`
to clearly state it contains energy values (proportional to pressure^2), not signed pressure.
Doc comments on the struct and its `to_time_series` / `band_to_time_series` methods now
note the energy nature and that it is suitable for metrics (RT, C80, D50) but not for
convolution reverb.

#### ~~LABELING: STI should be called "STI-like index"~~ (FIXED)

**Fixed**: Renamed `sti()` to `sti_approximate()` in `src/sim/acoustics/metrics.rs`.
The struct field `RoomAcousticReport::sti` is now `sti_approximate`. The docstring
documents the limitations (1 modulation frequency, no noise/masking, 6 bands).

### 5.3 Lighting

#### ~~BUG: Point light total flux formula is wrong~~ (FIXED)

**Fixed**: Removed the erroneous `/3` from `total_flux()`. The formula is now
`total_flux = (I_R + I_G + I_B) * 4*pi`. The `white()` constructor was updated to
`I_per_channel = flux / (3 * 4*pi)` so that `total_flux()` recovers the original flux.
`AreaLight::total_flux()` similarly fixed.

#### ~~BUG: RGB photometric/radiometric unit confusion~~ (FIXED)

**Fixed**: All lighting module documentation now uses radiometric units:
- Light source intensity: W/sr per channel
- Surface irradiance: W/m^2 per channel
- Incident flux: W per channel
- Backward tracer variable renamed from `intensity` to `radiance`

#### ~~BUG: Missing inverse-square verification for point lights~~ (FIXED)

**Fixed**: A unit test (`test_inverse_square_falloff`) verifies that 1/r^2 falloff
emerges naturally from the Monte Carlo solid angle sampling. A point light placed
off-center in a box confirms that the near surface receives significantly more flux
than the far surface.

#### ~~BIAS: Solar azimuth numerically fragile~~ (FIXED)

**Fixed**: Replaced `acos`-based azimuth calculation with `atan2(sin_azi, cos_azi)` in
`SolarPosition::calculate()` (`src/sim/lighting/solar.rs`). This preserves quadrant
information and is numerically robust for all latitude/declination combinations.

#### ~~LABELING: Backward tracer variable naming~~ (FIXED)

**Fixed**: Renamed `intensity` to `radiance` in the backward tracer code, consistent
with the radiometric unit convention adopted across the lighting module.

### 5.4 Energy

#### ~~BUG: Solar gains computed from lighting lumens instead of solar radiation~~ (FIXED)

**Fixed**: Solar gains are now computed directly from EPW weather data (DNI + DHI) using
sun position, window orientation (`cos(theta_incidence)`), isotropic sky view factor, and
SHGC. The `solar_gain_factor: f64` parameter was replaced with
`solar_config: Option<&SolarGainConfig>` in both `run_annual_simulation()` and
`run_transient_simulation()`. The legacy `SolarBridgeConfig`/`lighting_to_solar_gains()`
functions remain for backward compatibility.

#### ~~BIAS: HVAC load calculation ignores concurrent losses~~ (FIXED)

**Fixed**: Added `HvacIdealLoads::calculate_with_losses()` method in `src/sim/energy/hvac.rs`
that uses the implicit formula:
`Q_hvac = C*(T_set - T_zone)/dt + (UA + Inf_cond)*(T_set - T_out) - Q_gains`.
Updated `run_transient_simulation()` in `src/sim/energy/simulation.rs` to use this method.
The original `calculate_for_timestep()` is preserved for standalone use cases.

#### ~~LABELING: Window detection by name pattern~~ (FIXED)

**Fixed**: Added `is_glazing: bool` field to `Material` (`src/sim/materials.rs`) with a
`with_glazing()` builder method. The glass preset is marked as glazing. Added
`compute_solar_gains_with_materials()` in `src/sim/energy/solar_bridge.rs` that accepts
an optional `MaterialLibrary` and checks `is_glazing` before falling back to name pattern
matching. The original `compute_solar_gains()` delegates with `None` for backward
compatibility.

### 5.5 Prioritized Fix Order

Issues are ordered by impact on result correctness:

| Priority | Issue | Severity | Status |
|----------|-------|----------|--------|
| 1 | Solar gains from weather data, not lumens | BUG | **FIXED** |
| 2 | Point light total flux formula | BUG | **FIXED** |
| 3 | RGB unit convention (radiometric) | BUG | **FIXED** |
| 4 | Verify inverse-square for point lights | BUG | **FIXED** |
| 5 | Voxel grid ray marching (3D-DDA) | BUG | **FIXED** |
| 6 | Cosine-weighted diffuse reflection | BIAS | **FIXED** |
| 7 | Receiver normalization | BIAS | **FIXED** |
| 8 | HVAC concurrent-loss formula | BIAS | **FIXED** |
| 9 | Solar azimuth atan2 | BIAS | **FIXED** |
| 10 | Barycentric tolerance reduction | BIAS | **FIXED** |
| 11 | Rename STI to STI-like | LABELING | **FIXED** |
| 12 | Rename IR to energy IR | LABELING | **FIXED** |
| 13 | Backward tracer variable names | LABELING | **FIXED** |
| 14 | Window detection by material flag | LABELING | **FIXED** |
| 15 | Document 1R1C capacity as tuning param | LABELING | **FIXED** |

### 5.6 Acceptable Simplifications

The following are deliberate simplifications that do not need fixing, but should remain
clearly documented:

- Energy-only acoustics (no phase/interference) -- standard for geometric room acoustics
- No edge diffraction -- acceptable at mid/high frequencies
- Diffuse-only lighting reflections -- acceptable for matte interiors
- 1R1C zone thermal model -- standard simplified method (ISO 13790)
- Steady-state option for quick estimates
- Fixed atmospheric conditions for air absorption
- No latent loads in energy model

---

## 6. Potential Improvements

### 6.1 Acoustics

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

### 6.2 Lighting

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

### 6.3 Energy

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
