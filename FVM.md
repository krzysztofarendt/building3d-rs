# FVM Heat Transfer Runtime (Current)

This document describes the active FVM heat-transfer runtime used by energy simulations.

## Runtime Policy

Transient building energy simulation uses:
- FVM walls enabled
- global simultaneous solve enabled
- optional VF interior longwave coupling

Supported runtime modes:
1. `global + VF=0` (default)
2. `global + VF=1`

Legacy simplified transient runtime code paths were removed. Historical notes are in `LEGACY_SIMULATION.md`.

## Implemented FVM Stack

`src/sim/heat_transfer/`:
- `mesh.rs`: `FvmCell`, `FvmFace`, `FvmMesh` core mesh types
- `mesh_1d.rs`: layered wall -> 1D FVM mesh
- `solver.rs`: implicit tridiagonal wall solver (Thomas)
- `mesh_3d.rs`, `solver_sparse.rs`: generic APIs (not primary BESTEST runtime path)
- `boundary.rs`: boundary condition types

## Governing Discretization

For each cell:

```
(C_i / dt) * (T_i^{n+1} - T_i^n) = sum_f K_f * (T_nb^{n+1} - T_i^{n+1}) + S_i
```

where:
- `C_i = rho * cp * V`
- `K_f` uses series resistance at interfaces
- backward Euler in time

## Boundary Conditions Used in Energy Runtime

- `Convective`
- `ConvectiveWithFlux`
- `ConvectiveWithFluxToDomain`
- `Neumann` (for adiabatic internal-mass side where applicable)

## Global Coupled Assembly

In `src/sim/energy/global_solve.rs`, one system is assembled per substep with:
- wall cells
- interior surface nodes
- explicit steady exterior surface nodes (e.g. windows)
- air node(s)

Optional VF mode adds interior longwave matrix terms (`h_rad * ScriptF`) directly into the same global linear system.

## Surface Eligibility

FVM exterior wall eligibility (energy runtime):
- exterior surface
- resolvable `WallConstruction`
- not glazing-like
- no explicit per-surface U override
- not ground-coupled skip case

Non-FVM exterior surfaces (notably windows) are represented as steady interior nodes in the global solve.

## Configuration Notes

Active transient mode selector:
- `use_view_factor_radiation` (`false` default, `true` optional)

## Validation Snapshot (BESTEST)

Current annual deviations vs local OpenStudio/E+ reference:

| Mode | 600 Heat | 600 Cool | 900 Heat | 900 Cool |
|---|---:|---:|---:|---:|
| Global + VF=0 | -10.8% | -15.1% | +15.7% | +3.8% |
| Global + VF=1 | -6.0% | -21.0% | +28.6% | -10.7% |

Default mode is `VF=0`.
