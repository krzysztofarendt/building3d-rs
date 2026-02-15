# BESTEST ASHRAE 140 Energy Suite

## What is BESTEST?

BESTEST (Building Energy Simulation Test) is a standardized diagnostic methodology
developed at NREL by Judkoff and Neymark (1995), later codified as
[ANSI/ASHRAE Standard 140](https://www.ashrae.org/technical-resources/bookstore/standard-140-2020)
("Standard Method of Test for the Evaluation of Building Energy Analysis Computer
Programs"). The current edition is Standard 140-2020.

It uses a three-tier validation framework:

1. **Analytical verification** -- test cases with exact mathematical solutions
   (e.g., steady-state conduction through a known wall).
2. **Empirical validation** -- comparison to measured data from instrumented test
   cells (e.g., ETNA test cells in France, LBNL's FLEXLAB).
3. **Comparative testing** -- the dominant method: multiple reference programs
   (EnergyPlus, TRNSYS, ESP-r, DOE-2, BLAST, etc.) run the same cases, and
   their result spread defines an acceptance range.

The acceptance range reflects the consensus of state-of-the-art programs, not
ground truth. If all programs share a systematic error, BESTEST won't catch it.
BESTEST is best understood as a diagnostic and quality-assurance tool -- it catches
bugs and outlier behavior, but does not guarantee a match to physical reality.

## References

- ANSI/ASHRAE Standard 140-2020, *Standard Method of Test for the Evaluation
  of Building Energy Analysis Computer Programs*
  ([bookstore](https://www.ashrae.org/technical-resources/bookstore/standard-140-2020),
  [resource files](https://data.ashrae.org/standard140/))
- Judkoff, R. and Neymark, J. (1995), *International Energy Agency Building
  Energy Simulation Test (BESTEST) and Diagnostic Method*, NREL/TP-472-6231
  ([PDF](https://www.nrel.gov/docs/legosti/old/6231.pdf))
- Neymark, J. et al. (2020), *Update of ASHRAE Standard 140 Section 5.2 and
  Related Sections*, ANL-20/26
  ([PDF](https://publications.anl.gov/anlpubs/2020/05/158451.pdf))
- NREL BESTEST-GSR automated test framework
  ([GitHub](https://github.com/NREL/BESTEST-GSR))
- DOE ASHRAE 140 Maintenance and Development
  ([website](https://www.energy.gov/eere/buildings/ashrae-standard-140-maintenance-and-development))

## This suite

Runs 17 BESTEST Section 5.2 cases covering the building thermal fabric tests
(cases 600--960). Reference annual heating/cooling loads and peak loads are
compared against OpenStudio/EnergyPlus results from the
[NREL/BESTEST-GSR](https://github.com/NREL/BESTEST-GSR) workflow, using the
Boston Logan TMY3 EPW weather file.

## Common specifications

All cases share the following base geometry and parameters
(per ASHRAE 140-2020, Section 5.2.1):

| Parameter | Value |
|-----------|-------|
| Floor plan | 8 m x 6 m (48 m² floor area) |
| Height | 2.7 m |
| Zone volume | 129.6 m³ |
| Total window area | 12 m² (two windows, each 3 m x 2 m) |
| Window sill height | 0.2 m above floor |
| Glazing | Double-pane clear, SHGC ~ 0.76 |
| Internal gains | 200 W continuous (60% radiant, 40% convective) |
| Infiltration | 0.5 ACH (constant, no wind/temperature dependence) |
| Heating setpoint | 20 C (unless noted) |
| Cooling setpoint | 27 C (unless noted) |
| HVAC system | 100% convective, ideal (infinite capacity) |
| Ground temperature | 10 C (constant) |
| Site | Boston Logan TMY3 (42.36 N, 71.01 W) |

## Case descriptions

### 600-series (lightweight construction)

Exterior walls: wood siding (9 mm) + fiberglass insulation (66 mm) +
plasterboard (12 mm). Roof: wood deck + fiberglass + plasterboard.
Floor: timber flooring over R-25 insulation.

| Case | Description | Key modification from Case 600 |
|------|-------------|-------------------------------|
| **600** | Base case, lightweight | South-facing windows (2 x 3m x 2m), all defaults |
| **610** | South overhang | 1.0 m deep horizontal overhang at roof level above each south window; no lateral extension beyond window edges |
| **620** | East/west windows | Windows moved from south wall to east and west walls (one 3m x 2m window on each); same total 12 m² glazing area |
| **630** | East/west windows + overhang + fins | Case 620 windows with 1.0 m deep overhang and 1.0 m deep vertical side fins on each window |
| **640** | Thermostat setback | Heating setpoint = 10 C from 23:00--07:00, 20 C from 07:00--23:00; cooling setpoint unchanged (27 C) |
| **650** | Night ventilation cooling | No mechanical heating; cooling at 27 C from 07:00--18:00 only; night ventilation: 1703.16 m³/h (~13.14 ACH) from 18:00--07:00 in addition to base 0.5 ACH |
| **600FF** | Free-float | No mechanical heating or cooling; zone temperature floats freely |
| **650FF** | Free-float + night ventilation | No HVAC; same night ventilation schedule as Case 650 |

### 900-series (heavyweight construction)

Exterior walls: wood siding (9 mm) + foam insulation (61.5 mm) + concrete
block (100 mm). Floor: concrete slab (80 mm) over R-25 insulation.
Roof: same as 600-series (lightweight). Overall wall R-values are similar to
the 600-series but thermal mass is significantly higher.

| Case | Description | Key modification from Case 900 |
|------|-------------|-------------------------------|
| **900** | Base case, heavyweight | Same geometry and windows as Case 600; heavyweight wall and floor constructions |
| **910** | South overhang | Same overhang as Case 610 applied to heavyweight building |
| **920** | East/west windows | Same window relocation as Case 620 with heavyweight constructions |
| **930** | East/west windows + overhang + fins | Same shading as Case 630 with heavyweight constructions |
| **940** | Thermostat setback | Same setback schedule as Case 640 with heavyweight constructions |
| **950** | Night ventilation cooling | Same ventilation schedule as Case 650 with heavyweight constructions |
| **900FF** | Free-float | No HVAC; heavyweight building floats freely |
| **950FF** | Free-float + night ventilation | No HVAC; same night ventilation as Case 950 |

### Case 960 (two-zone sunspace)

| Case | Description |
|------|-------------|
| **960** | Sunlit zone (south, 8 m x 2 m) with heavyweight construction and 12 m² south windows, coupled via partition wall to a back zone (6 m x 6 m, lightweight, no windows, unconditioned). Only the south zone is conditioned. Tests inter-zone heat transfer. **Not yet implemented** -- requires multi-zone geometry and separate HVAC zones. |

## What each case tests

The cases are organized as a diagnostic sequence. Each case modifies one
parameter from a base case, isolating a specific physical mechanism:

| Mechanism tested | Cases |
|-----------------|-------|
| Basic conduction, solar gain, internal loads | 600, 900 |
| Effect of thermal mass (light vs heavy) | 600 vs 900 |
| Overhang shading of beam solar | 610, 910 |
| Window orientation (east/west vs south) | 620, 920 |
| Combined shading (overhang + fins) | 630, 930 |
| Thermostat setback control | 640, 940 |
| Night ventilation (passive cooling) | 650, 950 |
| Free-float temperature prediction | 600FF, 900FF |
| Free-float + ventilation | 650FF, 950FF |
| Multi-zone coupling | 960 |

## Data download (EPW)

```bash
bash examples/bestest_energy_suite/download_data.sh
```

## Run

```bash
cargo run --example bestest_energy_suite --release
```

### Outputs

| File | Content |
|------|---------|
| `results.csv` | Per-case monthly and annual heating/cooling loads with E+ reference comparison |
| `diagnostics.csv` | UA breakdown, solar gains, peak load timing per case |
| `diagnostics_monthly.csv` | Monthly solar transmitted, opaque sol-air, ground correction |

### Plotting

```bash
python3 examples/bestest_energy_suite/plot_results.py
```

Generates comparison figures in the same directory:

| Figure | Content |
|--------|---------|
| `bestest_600_series.png` | 600-series annual heating & cooling vs E+ reference |
| `bestest_900_series.png` | 900-series annual heating & cooling vs E+ reference |
| `bestest_peaks.png` | Peak heating & cooling loads (2x2 grid) |
| `bestest_freefloat.png` | Free-float min/max zone temperatures |
| `bestest_error_summary.png` | Percentage error overview (all HVAC cases) |
| `bestest_NNN_monthly.png` | Monthly heating & cooling profile for case NNN |

## Environment variables

| Variable | Default | Description |
|----------|---------|-------------|
| `BESTEST_600_EPW` | `data/USA_MA_Boston-Logan...epw` | Path to EPW weather file |
| `BESTEST_WARMUP_DAYS` | `7` | Simulation warmup period |
| `BESTEST_SUBSTEPS_PER_HOUR` | `6` | Timesteps per hour (10-min intervals at default) |
| `BESTEST_VF` | `0` | Enable view-factor interior radiation |
| `BESTEST_GLOBAL_VF` | `0` | Enable view-factor radiation in global FVM solve |
| `BESTEST_SOLAR_ALPHA` | `0.6` | Interior solar absorptance (0--1) |

## Implementation notes

- Wall conduction uses a finite volume method (FVM) with multi-layer constructions.
- Exterior convection uses the DOE-2 model (wind-dependent combined coefficient).
- Interior convection uses the TARP model (buoyancy-driven, tilt-dependent).
- Solar gains use a polynomial angular SHGC model for clear glass, with
  ground-reflected diffuse and exterior opaque sol-air absorption.
- Overhang and fin shading is computed analytically from sun position and
  device geometry (see `src/sim/energy/shading.rs`).
- HVAC schedule overrides (setback, night ventilation, free-float) are
  implemented as hourly vectors passed to `TransientSimulationOptions`.
- The 900-series uses the same roof construction as 600 (lightweight),
  matching the ASHRAE 140 specification.
- Case 960 (two-zone sunspace) is listed but skipped pending multi-zone
  geometry support.
