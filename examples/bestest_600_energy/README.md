# BESTEST 600 Feature Sweep

This example runs a curated feature sweep for BESTEST Case 600 against EnergyPlus reference values.

## Run

Make sure the EPW is available. The example defaults to:

- `examples/bestest_energy_suite/data/USA_MA_Boston-Logan.Intl.AP.725090_TMY3.epw`

You can either:

- Download to the default location: `bash examples/bestest_energy_suite/download_data.sh`
- Or set an explicit weather file: `BESTEST_600_EPW=/path/to/file.epw`

Then run:

```bash
cargo run --example bestest_600_energy --release
```

## Outputs

- `examples/bestest_600_energy/results.csv` (baseline monthly/annual errors)
- `examples/bestest_600_energy/feature_sweep.csv` (all sweep cases)
- `examples/bestest_600_energy/feature_sweep_barplot.svg` (annual heating/cooling error bars)

## Environment Controls

- `BESTEST_SWEEP_CASES`: comma-separated case IDs to run (default: all)
- `BESTEST_VF_RAYS`: view-factor rays for `suite_view_factor_radiation` (default: `10000`)
- `BESTEST_WARMUP_DAYS`: default warmup days (default: `7`)
- `BESTEST_SUBSTEPS_PER_HOUR`: default substeps per hour (default: `6`)

Example:

```bash
BESTEST_SWEEP_CASES=suite_baseline,suite_no_surface_aware cargo run --example bestest_600_energy --release
```

## Cases

All cases are one-at-a-time variants around the suite baseline unless noted.

| Case ID | Label | What changes |
|---|---|---|
| `legacy_fixed_floor_mass` | Legacy fixed + floor mass | Legacy-style setup: fixed interior/exterior convection, no beam split, add explicit floor internal mass. |
| `suite_baseline` | Suite baseline | Baseline matching `bestest_energy_suite` Case 600 settings. |
| `suite_plus_floor_mass` | Suite + explicit floor mass | Adds one-sided floor internal mass slab. |
| `suite_fixed_convection` | Suite + fixed convection | Replaces TARP/DOE2 convection with fixed coefficients. |
| `suite_no_interior_radiation` | Suite + no interior radiation | Disables interior radiative exchange (and view-factor radiation). |
| `suite_view_factor_radiation` | Suite + view-factor LW | Enables per-surface view-factor longwave radiation. |
| `suite_no_beam_split` | Suite + no beam split | Disables beam/diffuse split for interior transmitted solar handling. |
| `suite_no_surface_aware` | Suite + no surface-aware solar | Disables surface-aware transmitted-solar distribution path. |
| `suite_transmit_to_air_30` | Suite + 30% transmitted solar to air | Sets `transmitted_solar_to_air_fraction = 0.3`. |
| `suite_distribute_to_fvm_walls` | Suite + distribute to FVM walls | Includes FVM wall area for interior source distribution. |
| `suite_wall_share_to_air` | Suite + wall share rerouted to air | When distributing to FVM walls, routes wall share to air node. |
| `suite_alpha_1p0` | Suite + alpha = 1.0 | Sets interior solar absorptance to 1.0 (no bounce pool). |
| `suite_no_warmup_substeps1` | Suite + no warmup + 1 substep | Coarse solver controls: warmup `0` days and `1` substep/hour. |
