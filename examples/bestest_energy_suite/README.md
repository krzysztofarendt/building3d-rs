# BESTEST Energy Suite (600 + 900)

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

