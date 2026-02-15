# Next Step: Improve Case 900 BESTEST Results

## Current State

**Baseline results** (VF OFF, TARP interior, 1R1C model):
- 600: heating -1.4%, cooling -7.7% (excellent)
- 900: heating +30.0%, cooling +15.4% (structural)
- Total absolute deviation: 54.5pp

## What Was Tried (and failed)

The iterative surface heat balance plan was **fully implemented** but does not
improve BESTEST results with the 1R1C model. See `memory/iterative_balance_findings.md`
for detailed experimental results (7+ configurations tested).

**Root cause**: Adding interior radiation coupling (h_rad ≈ 4.6 W/m²K) to the
1R1C model increases effective building U-value by ~14% and reduces floor mass
time constant from 22.5h to 7.4h. Both effects worsen Case 900. The 1R1C model
was implicitly calibrated with TARP-only h ≈ 2.5 W/m²K.

All infrastructure remains in place but disabled by default:
- `FvmSolverSnapshot` save/restore
- `use_iterative_surface_balance` config flag
- Iteration loop in 1R1C substep
- Ground-coupled floor-as-FVM guard
- `convective_only_air_gain` and `floor_beam_sources` step function params

## Option A: Per-Surface Heat Balance Model

Replace 1R1C with a **surface-by-surface simultaneous solve** (like EnergyPlus CTF):

Each surface i has a heat balance equation:
```
q_cond_i + h_conv_i*(T_air - T_i) + h_rad*sum_j(F_ij*(T_j - T_i)) + q_solar_i = 0
```

All N surface equations + 1 zone air equation solved simultaneously via direct
matrix solve (N+1 × N+1 system). No iteration needed.

**Pros**: Physically correct, naturally handles radiation coupling, matches E+.
**Cons**: Major refactor of the thermal model — new module, not a patch on 1R1C.

## Option B: Empirical Correction for Case 900

Adjust thermal mass parameters or add a second mass node to better represent
the distributed mass in heavyweight construction:
- Two-mass model: floor slab + wall mass with different time constants
- Calibrate against Case 900 reference data
- Quick implementation, limited generality

## Option C: Accept Current Results

54.5pp total deviation is reasonable for a simplified model. Case 600 is excellent.
Case 900's +30% heating is a known limitation of single-mass models without
interior radiation coupling.

Focus development effort on other simulation domains instead.

---

## Recommendation

Option A is the correct long-term path but is a significant effort (~2-3 sessions).
Option C is reasonable if BESTEST validation is not the primary goal.
