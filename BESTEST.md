# BESTEST Notes (building3d)

## Scope

This repository currently tracks a BESTEST-inspired regression subset:
- Case 600 (light mass)
- Case 900 (heavy mass)
- Diagnostic variants: `600_no_solar`, `900_no_solar`

Weather source for regression runs:
- Boston Logan TMY3 EPW (same EPW used by the local OpenStudio/EnergyPlus reference extraction)

## Runtime Modes in building3d

Only two transient runtime modes are supported:

1. `Global solve + VF=0` (default)
2. `Global solve + VF=1`

`Global solve` means the simulator solves wall cells, surface nodes, and air node together in one coupled system each substep.

Legacy simplified transient runtime code is removed. Historical behavior is documented in `LEGACY_SIMULATION.md`.

## Current Annual Results Snapshot

Reference annual values (OpenStudio/E+):
- 600: heating `4324.8`, cooling `6044.1` kWh
- 900: heating `1661.2`, cooling `2498.2` kWh

building3d current modes:

| Mode | 600 Heating | 600 Cooling | 900 Heating | 900 Cooling |
|---|---:|---:|---:|---:|
| Global + VF=0 (default) | -10.8% | -15.1% | +15.7% | +3.8% |
| Global + VF=1 | -6.0% | -21.0% | +28.6% | -10.7% |

Default recommendation:
- Use `VF=0` for best overall annual agreement across the 600/900 pair.

## Why Differences vs EnergyPlus Remain

Main gap is window/surface thermal fidelity:
- EnergyPlus: dynamic layered glazing and iterative inside/outside surface heat balances.
- building3d: simplified steady interior window node + `U*A` path to outdoor.

This difference affects indoor longwave redistribution and annual heating/cooling partition.

## Running the Example

Download EPW:

```bash
bash examples/bestest_energy_suite/download_data.sh
```

Run (default mode, VF off):

```bash
cargo run --example bestest_energy_suite --release
```

Run with VF on:

```bash
BESTEST_GLOBAL_VF=1 cargo run --example bestest_energy_suite --release
```

Outputs:
- `examples/bestest_energy_suite/results.csv`
- `examples/bestest_energy_suite/diagnostics.csv`
- `examples/bestest_energy_suite/diagnostics_monthly.csv`

## Test Intent

`tests/bestest_energy_suite.rs` is a regression stability guard for current model behavior under the two supported global modes, not an ASHRAE 140 conformance certification harness.

## References

- Judkoff, R.; Neymark, J. (1995), NREL/TP-472-6231
- ANSI/ASHRAE Standard 140
- NREL BESTEST-GSR workflow reference outputs
