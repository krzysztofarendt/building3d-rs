# BRAS Benchmark - Summary for building3d-rs Validation

## Overview

The BRAS (Benchmark for Room Acoustical Simulation) is a CC BY-SA 4.0 licensed dataset
from RWTH Aachen / TU Berlin containing 11 acoustic scenes with measured impulse responses,
geometry, and material data for validating room acoustic simulation software.

**Source**: https://depositonce.tu-berlin.de/items/version/64

## Complex Rooms Relevant for Validation

| Scene | Name | Volume | Materials | Feasibility |
|-------|------|--------|-----------|-------------|
| CR2 (Scene 09) | Small seminar room, RWTH Aachen | 145 m³ | 5 | **Best candidate** |
| CR3 (Scene 10) | Chamber music hall, Konzerthaus Berlin | 2,350 m³ | 7 | Good, complex geometry |
| CR4 (Scene 11) | Large auditorium, TU Berlin | 8,650 m³ | 8 | Good, complex geometry |

## CR2: Small Seminar Room (Primary Target)

### Room Properties
- **Volume**: 145 m³ (S = 202.53 m², V = 146.1 m³ from 3D model)
- **Dimensions**: ~8.4 m × 6.7 m × 3 m (roughly rectangular)
- **Temperature**: 19.5°C
- **Humidity**: 41.7%
- **Sampling rate**: 44100 Hz

### Source Positions (from RAVEN .rpf file, x,y,z in meters)
The RAVEN coordinate system has Y=up:
- **LS1-mid**: (0.931, 0.723, 2.547) → in BRAS global coords
- **LS2-mid**: (0.119, 0.723, -2.880) → in BRAS global coords

Note: Source heights include mid-freq tweeter at 0.723m (in RAVEN Y-up coords).
The Treble docs report: LS1=(0.93, -2.55, 0.72), LS2=(0.19, 2.88, 0.72) in a different
coordinate convention (likely z-up). The RAVEN file uses Y-up, so positions are:
- LS1: x=0.931, height=0.723, z=2.547
- LS2: x=0.119, height=0.723, z=-2.880

### Receiver Positions (from RAVEN .rpf file, x,y,z)
All receivers at height 1.230 m (Y-up):
- **MP1**: (-0.993, 1.230, -1.426)
- **MP2**: (0.439, 1.230, 0.147)
- **MP3**: (1.361, 1.230, 0.603)
- **MP4**: (-1.110, 1.230, 0.256)
- **MP5**: (-0.998, 1.230, 1.409)

### Materials (5 materials, fitted estimates)
All values are 31 third-octave bands from 20 Hz to 20 kHz.
Format: line 1 = frequencies, line 2 = absorption, line 3 = scattering.

| Material | Surface Area | Characteristic Depth |
|----------|-------------|---------------------|
| plaster | 34.93 m² | 0.006 m |
| concrete | 56.94 m² | 0.005 m |
| windows | 9.75 m² | 0.002 m |
| floor | 49.29 m² | 0.003 m |
| ceiling | 51.62 m² | 0.060 m |

#### Absorption at octave bands (fitted estimates):
| Material | 125 Hz | 250 Hz | 500 Hz | 1000 Hz | 2000 Hz | 4000 Hz |
|----------|--------|--------|--------|---------|---------|---------|
| plaster  | 0.033  | 0.050  | 0.039  | 0.044   | 0.048   | 0.036   |
| concrete | 0.085  | 0.075  | 0.056  | 0.059   | 0.059   | 0.044   |
| windows  | 0.175  | 0.073  | 0.049  | 0.057   | 0.133   | 0.055   |
| floor    | 0.071  | 0.091  | 0.070  | 0.065   | 0.062   | 0.043   |
| ceiling  | 0.083  | 0.104  | 0.048  | 0.049   | 0.047   | 0.062   |

#### Scattering at octave bands (same for all materials except ceiling):
Scattering is frequency-dependent, estimated from characteristic depth via:
s(f) = 0.5 * sqrt(d_char / (c/f)), clamped to [0.05, 0.99]

### Geometry
- Available as AC3D format (ASCII text) in `informed_sim/RavenModels/scene9.ac`
- 5 OBJECT poly groups, one per material:
  - Object 0: mat_scene09_concrete (mat 0) - 194 vertices, 147 surfaces
  - Object 1: mat_scene09_windows (mat 1) - 12 vertices, 3 surfaces (3 window panes)
  - Object 2: mat_scene09_ceiling (mat 2) - 102 vertices, 134 surfaces
  - Object 3: mat_scene09_plaster (mat 3) - 25 vertices, 17 surfaces
  - Object 4: mat_scene09_floor (mat 4) - 31 vertices, 29 surfaces
- Materials indexed 0-4: concrete, windows, ceiling, plaster, floor

### Measured Room Acoustic Parameters (averaged across all source-receiver pairs)
Data in `additional_data/4 Additional data/Room acoustic parameters/data/`
24 third-octave bands from 80 Hz to 16 kHz.

**CR2 at standard octave bands:**

| Band | EDT (s) | T20 (s) | C80 (dB) | D50 (%) |
|------|---------|---------|----------|---------|
| 125 Hz | 1.501 | 1.453 | 0.551 | 29.5 |
| 250 Hz | 1.306 | 1.345 | 0.902 | 39.4 |
| 500 Hz | 2.018 | 2.018 | -1.862 | 27.8 |
| 1000 Hz | 1.919 | 1.934 | -0.606 | 32.1 |
| 2000 Hz | 1.758 | 1.715 | -0.483 | 35.6 |
| 4000 Hz | 1.649 | 1.604 | 0.164 | 36.6 |

### RAVEN Simulation Settings (for reference)
- Image source order: 2
- Number of ray tracing particles: 500,000
- Filter length: 2800 ms
- Energy loss threshold: 75 dB
- Detection sphere radius: 0.8 m
- Time resolution: 2 ms
- Third-octave band resolution
- Air absorption: enabled

### WARNING
The level of the RIRs was manually corrected by 6 dB to meet the expected direct sound level.

## CR3: Chamber Music Hall (Konzerthaus Berlin)

### Room Properties
- Volume: 2,350 m³ (excluding attic)
- Temperature: 22.4°C, Humidity: 40.9%
- 7 materials: ceiling, floor, plaster, seating, stagePanels, structuredPlaster, windows
- Note: Considerable volume above ceiling (attic), coupled via connections

### Measured Parameters at octave bands:
| Band | EDT (s) | T20 (s) | C80 (dB) | D50 (%) |
|------|---------|---------|----------|---------|
| 125 Hz | 1.694 | 1.563 | -0.955 | 27.7 |
| 250 Hz | 1.564 | 1.497 | 0.566 | 40.3 |
| 500 Hz | 1.457 | 1.258 | 1.241 | 42.3 |
| 1000 Hz | 1.316 | 1.342 | 1.348 | 39.9 |
| 2000 Hz | 1.273 | 1.320 | 1.398 | 41.2 |
| 4000 Hz | 1.087 | 1.049 | 2.724 | 49.0 |

## CR4: Large Auditorium (TU Berlin)

### Room Properties
- Volume: 8,650 m³
- Temperature: 20.9°C, Humidity: 37.5%
- 8 materials: brickwall, concrete, linoleum, parquet, seating, whitePanels, windows, woodPanels

### Measured Parameters at octave bands:
| Band | EDT (s) | T20 (s) | C80 (dB) | D50 (%) |
|------|---------|---------|----------|---------|
| 125 Hz | 2.231 | 2.263 | -1.407 | 31.8 |
| 250 Hz | 2.599 | 2.378 | -1.730 | 27.3 |
| 500 Hz | 2.108 | 2.044 | -0.341 | 34.4 |
| 1000 Hz | 2.256 | 2.145 | -2.151 | 27.3 |
| 2000 Hz | 1.903 | 1.887 | -0.350 | 36.3 |
| 4000 Hz | 1.502 | 1.405 | 1.680 | 45.4 |

## Data Files Available Locally

```
validation/bras/
├── Documentation.pdf                    # Full BRAS documentation (28 pages)
├── DocumentationInformedSim.pdf         # RWTH informed simulation docs
├── SUMMARY.md                           # This file
├── additional_data/
│   └── 4 Additional data/
│       └── Room acoustic parameters/
│           ├── calculateRoomAcousticParameters.m
│           ├── data/                    # EDT, T20, C80, D50 for CR2-CR4
│           └── plots/                   # Pre-rendered comparison plots
├── surface_descriptions/
│   └── 3 Surface descriptions/
│       ├── AllMaterials.pdf
│       ├── _csv/
│       │   ├── fitted_estimates/        # Absorption/scattering for CR1-CR4
│       │   └── initial_estimates/       # Initial estimates
│       ├── _descr/                      # Material description text files
│       ├── _img/                        # Material photos
│       └── _plots/                      # Material property plots
└── informed_sim/
    ├── RavenModels/
    │   ├── scene9.ac                    # CR2 geometry (AC3D, ASCII)
    │   ├── scene9_RR.ac                 # CR2 round-robin variant
    │   └── scene11.ac                   # CR4 geometry (AC3D, ASCII)
    ├── RavenInput/
    │   ├── Scene09_Informed.rpf         # CR2 simulation config with positions
    │   └── Scene11_Informed.rpf         # CR4 simulation config
    └── RavenDatabase/
        └── MaterialDatabase/
            └── BRASv2_fitted_rev2019v3_oldNames/  # .mat material files
```

## AC3D File Format Reference

The `.ac` files are ASCII text with this structure:
```
AC3Db
MATERIAL "name" rgb r g b amb r g b emis r g b spec r g b shi N trans T
...
OBJECT world
kids N
OBJECT poly
name "polygon_object"
numvert N
x y z        # vertex coordinates (Y-up in RAVEN)
...
numsurf N
SURF 0x10
mat M        # material index (0-based, references MATERIAL declarations)
refs N       # number of vertices in this face
idx u v      # vertex index + UV coords (0 0 for no texture)
...
kids 0
```

## Coordinate System Notes

- **RAVEN/AC3D**: Y-up coordinate system
- **BRAS global**: x/y horizontal plane, z up (φ=0 → +x, φ=90 → +y, ϑ=90 → +z)
- **building3d-rs**: Typically z-up
- Conversion: RAVEN (x, y_up, z) → building3d (x, z, -y) or similar remapping needed

## Next Steps for Validation

1. **Parse AC3D geometry** → Extract vertices and faces per material
2. **Map to building3d hierarchy**: Create Building → Zone → Solid → Wall → Polygon
3. **Assign materials**: Map the 5 CR2 materials to AcousticMaterial with per-band absorption/scattering
4. **Set source/receiver positions**: From RAVEN .rpf coordinates (transform to building3d coords)
5. **Run simulation**: Match RAVEN settings (500k rays, detection sphere, etc.)
6. **Compare**: EDT, T20, C80, D50 at octave bands against measured values
