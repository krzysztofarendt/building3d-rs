# BESTEST Energy Suite (600 + 900)

## What is BESTEST?

BESTEST (Building Energy Simulation Test) is a standardized diagnostic methodology
developed at NREL by Judkoff and Neymark (1995), later codified as
ANSI/ASHRAE Standard 140. It uses a three-tier validation framework:

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

## This suite

Runs a small BESTEST-inspired validation suite for the energy simulation:

- **600 - Base Case** (light-mass)
- **600 - Base Case (no solar)** (internal sanity check; no external reference)
- **900 - High-Mass Base Case** (approximated via increased zone thermal capacity)

Reference monthly/annual heating and cooling loads for 600 and 900 are taken from the
NREL/BESTEST-GSR `workflow_results.csv` (OpenStudio/EnergyPlus), using the same Boston Logan
TMY3 EPW.

## Data download (EPW)

```bash
bash examples/bestest_energy_suite/download_data.sh
```

## Run

```bash
cargo run --example bestest_energy_suite --release
```

Outputs:
- `examples/bestest_energy_suite/results.csv`

## Notes

- The 900 case is **not a perfect mapping** to EnergyPlus BESTEST 900 (which uses a detailed
  high-mass construction and internal mass). In `building3d`, this suite approximates “high
  mass” by scaling the zone thermal capacity (`ThermalConfig::thermal_capacity_j_per_m3_k`).
- Tune the high-mass capacity scale with:
  - `BESTEST_900_CAPACITY_SCALE` (default: `8.0`)

