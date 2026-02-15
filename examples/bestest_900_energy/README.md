# BESTEST 900 Feature Sweep

This example runs a curated feature sweep for BESTEST Case 900 (heavy-mass) against EnergyPlus reference values.

## Run

Weather resolution order:

1. `BESTEST_900_EPW`
2. `BESTEST_600_EPW` (backward-compatible fallback)
3. `examples/bestest_energy_suite/data/USA_MA_Boston-Logan.Intl.AP.725090_TMY3.epw`

To place the EPW at the default location:

```bash
bash examples/bestest_energy_suite/download_data.sh
```

Run:

```bash
cargo run --example bestest_900_energy --release
```

## Outputs

- `examples/bestest_900_energy/results.csv` (baseline monthly/annual errors)
- `examples/bestest_900_energy/feature_sweep.csv` (all sweep cases)
- `examples/bestest_900_energy/feature_sweep_barplot.svg` (annual heating/cooling error bars)

## Environment Controls

- `BESTEST_900_SWEEP_CASES`: comma-separated case IDs to run
- `BESTEST_SWEEP_CASES`: fallback filter (for compatibility)
- `BESTEST_900_VF_RAYS`: view-factor rays for `suite_view_factor_radiation`
- `BESTEST_VF_RAYS`: fallback VF rays (for compatibility)
- `BESTEST_WARMUP_DAYS`: default warmup days (default: `7`)
- `BESTEST_SUBSTEPS_PER_HOUR`: default substeps per hour (default: `6`)

Example:

```bash
BESTEST_900_SWEEP_CASES=suite_baseline,suite_no_surface_aware cargo run --example bestest_900_energy --release
```

## Cases

All cases are one-at-a-time variants around the suite baseline.

| Case ID | Label | What changes |
|---|---|---|
| `suite_baseline` | Suite baseline | Baseline matching `bestest_energy_suite` Case 900 settings (heavy constructions). |
| `suite_fixed_convection` | Suite + fixed convection | Replaces TARP/DOE2 convection with fixed coefficients. |
| `suite_no_interior_radiation` | Suite + no interior radiation | Disables interior radiative exchange (and view-factor radiation). |
| `suite_view_factor_radiation` | Suite + view-factor LW | Enables per-surface view-factor longwave radiation. |
| `suite_no_beam_split` | Suite + no beam split | Disables beam/diffuse split for interior transmitted solar handling. |
| `suite_no_surface_aware` | Suite + no surface-aware solar | Disables surface-aware transmitted-solar distribution path. |
| `suite_transmit_to_air_30` | Suite + 30% transmitted solar to air | Sets `transmitted_solar_to_air_fraction = 0.3`. |
| `suite_distribute_to_fvm_walls` | Suite + distribute to FVM walls | Includes FVM wall area for interior source distribution. |
| `suite_wall_share_to_air` | Suite + wall share rerouted to air | When distributing to FVM walls, routes wall share to air node. |
| `suite_alpha_1p0` | Suite + alpha = 1.0 | Sets interior solar absorptance to 1.0 (no bounce pool). |
| `suite_no_warmup` | Suite + no warmup | Solver control: warmup `0` days. |
| `suite_substeps_1` | Suite + 1 substep | Solver control: `1` substep/hour. |
