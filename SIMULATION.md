# Energy Simulation: Current Runtime Model

This document describes the active thermal simulation runtime in `building3d`.

## Runtime Policy

Only two transient runtime modes are supported:

1. `global + VF=0` (default)
2. `global + VF=1`

Where:
- `global` = global simultaneous FVM solve (wall cells + interior surface nodes + zone air node in one linear system)
- `VF` = interior longwave view-factor radiation coupling

Default behavior:
- global solve is always used in `run_transient_simulation*`
- `use_view_factor_radiation = false`

## Legacy Paths

Legacy simplified transient paths have been removed from runtime code:
- lumped envelope-only wall capacitance path (1R1C wall treatment)
- two-node / three-node lumped envelope mass variants
- non-global sequential wall/air coupling and iterative fallback path

Historical notes are preserved in `LEGACY_SIMULATION.md`.

## Active Numerical Model

### 1. Conduction

Opaque envelope surfaces with resolved `WallConstruction` are solved with 1D implicit FVM:

```
(C_i / dt) * (T_i^{n+1} - T_i^n) = sum_f K_f * (T_nb^{n+1} - T_i^{n+1}) + S_i
```

- Layer interfaces use series-resistance face conductance.
- Time integration is backward Euler.
- Each eligible polygon has its own wall solver mesh.

### 2. Global Coupled Solve

At each substep, a single global dense system is assembled:
- wall cell unknowns
- interior surface nodes
- steady exterior surface nodes (e.g. glazing interior nodes)
- zone air node(s)

The solve includes:
- wall conduction
- interior convection (TARP/Fixed)
- infiltration and glazing UA terms on air node
- solar/internal gains source terms
- optional interior longwave radiation matrix (VF mode)

### 3. Interior Longwave Radiation

`VF=0`:
- no surface-to-surface longwave enclosure matrix is applied.

`VF=1`:
- Monte Carlo geometric view factors are computed per zone.
- ScriptF exchange factors are derived from geometric F and emissivity.
- radiation coupling is embedded directly in the global matrix.

### 4. HVAC

Ideal loads HVAC clamps air nodes to heating/cooling setpoints by constrained re-solve and reports required heating/cooling power.

## BESTEST (Boston TMY3) Current Reference Snapshot

Reference annual values (OpenStudio/EnergyPlus):
- Case 600: heating `4324.8 kWh`, cooling `6044.1 kWh`
- Case 900: heating `1661.2 kWh`, cooling `2498.2 kWh`

Active runtime modes:

| Mode | 600 Heat | 600 Cool | 900 Heat | 900 Cool |
|---|---:|---:|---:|---:|
| Global, VF=0 (default) | -10.8% | -15.1% | +15.7% | +3.8% |
| Global, VF=1 | -6.0% | -21.0% | +28.6% | -10.7% |

Interpretation:
- `VF=0` is the default because it gives the best overall annual agreement across 600/900.
- `VF=1` changes interior longwave redistribution and shifts heating/cooling tradeoffs.

## Known Model Gap vs EnergyPlus

Largest remaining source of divergence is window/surface physics fidelity:
- EnergyPlus solves dynamic inside/outside window surface heat balances with layered glazing/gap physics.
- `building3d` uses simplified steady window interior nodes coupled by indoor film and outdoor `U*A`.

This difference materially affects annual heating/cooling partitioning in BESTEST.

## User-Facing Controls

For BESTEST example runs:
- Global solve is always ON.
- Toggle only VF:
  - default: `BESTEST_GLOBAL_VF=0`
  - optional: `BESTEST_GLOBAL_VF=1`

For API users:
- `use_view_factor_radiation` remains the active mode selector.
