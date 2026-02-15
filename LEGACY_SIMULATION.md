# Legacy Simulation History

This file preserves documentation of transient runtime paths that were removed from active code.

## Removal Summary

Legacy transient runtime paths were removed and replaced by a single global FVM runtime with optional view-factor radiation:
- `global + VF=0` (default)
- `global + VF=1`

## Removed Legacy Transient Paths

1. Lumped envelope-only wall capacitance path (1R1C wall treatment)
- Building envelope losses and thermal mass were represented by a coarse lumped formulation.

2. Two-node lumped envelope mass variants
- Air + mass 2R2C lumped model.
- Envelope-to-mass 2R2C variant.

3. Three-node lumped envelope variant
- Air + interior-surface + envelope coarse model.

4. Non-global sequential wall/air coupling path
- Per-wall sequential stepping with separate air update.
- Iterative surface-balance fallback loop.

## Removed Legacy APIs / Config Knobs

The following transient legacy controls were removed from `ThermalConfig`:
- `two_node_mass_fraction`
- `solar_gains_to_mass_fraction`
- `two_node_envelope_to_mass`
- `three_node_envelope_mass_fraction`
- `use_iterative_surface_balance`
- `surface_balance_max_iterations`
- `surface_balance_tolerance_c`
- `use_global_fvm_solve`

## Where to Find Old Behavior

Historical implementations remain available in git history before the legacy-removal change.

Recommended lookup strategy:
- inspect `src/sim/energy/simulation.rs` history for transient runtime branching,
- inspect `src/sim/energy/hvac.rs` history for removed 2R2C/3R3C model implementations.
