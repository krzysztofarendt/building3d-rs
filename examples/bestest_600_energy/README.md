# BESTEST 600 (Energy) — Reference-Backed Validation Example

This example implements a **BESTEST Case 600**-inspired single-zone building and compares
`building3d` energy results against **OpenStudio/EnergyPlus reference outputs** published by
NREL’s BESTEST-GSR workflow.

This is **not** a perfect apples-to-apples comparison (our thermal model is intentionally
simplified), but it provides a reproducible “anchor” case that can catch regressions and
quantify drift as the energy model evolves.

## What it runs

- Geometry: 8 m × 6 m × 2.7 m single zone with two south windows (per `case_en_600.idf`).
- HVAC: ideal loads with dual setpoints (20°C heating, 27°C cooling).
- Infiltration: 0.5 ACH (constant).
- Internal gains: 200 W (constant).
- Solar: simple SHGC model applied to polygons whose path contains `"window"`.
- Weather: TMY3 EPW used by the BESTEST-GSR Case 600 workflow (Boston Logan).

## Data download (EPW)

The EPW is not committed to this repo. Download it once:

```bash
bash examples/bestest_600_energy/download_data.sh
```

## Run

```bash
cargo run --example bestest_600_energy --release
```

This writes:
- `examples/bestest_600_energy/results.csv`

Optional plot:
```bash
python3 examples/bestest_600_energy/generate_figures.py
```

This writes:
- `examples/bestest_600_energy/results.png`

## Reference values

Reference monthly/annual heating and cooling loads are taken from the
`results/workflow_results.csv` line for **“600 - Base Case”** in NREL/BESTEST-GSR
(same commit pinned by `download_data.sh`).

The benchmark stores only derived numeric results (no ASHRAE spreadsheets).
