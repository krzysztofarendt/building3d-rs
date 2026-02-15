  Bugs / wrong assumptions (high confidence)

  1. BESTEST regression harness drift.
     tests/bestest_energy_suite.rs:478 and tests/bestest_energy_suite.rs:529 describe an older baseline, but current config path no
     longer matches it, and with EPW present both reference tests fail.
  2. Global VF path excludes windows from the radiation matrix.
     VF builder includes all zone polygons (src/sim/energy/view_factors.rs:527), but global assembly only keeps fvm_walls +
     internal_mass_surfaces (src/sim/energy/simulation.rs:2499), so window handles are dropped.
  3. Sequential iterative VF coupling is not fully coupled.
     Inside the iteration loop, MRT is recomputed but still uses t_air_start_c (src/sim/energy/simulation.rs:2685) and wall
     callbacks still use t_rad_guess_c (src/sim/energy/simulation.rs:2713, src/sim/energy/simulation.rs:2771), so the “iterative”
     solve is partially stale.
  4. Inconsistent solar absorptance defaults across codepaths.
     ThermalConfig::new() default is 0.6 (src/sim/energy/config.rs:325), while BESTEST example defaults to 1.0 unless env override
     (examples/bestest_energy_suite/main.rs:414).

  Model limitations vs EnergyPlus (not simple bugs)

  1. No dynamic window surface-node heat balance.
     Glazing is a UA sink on air node (src/sim/energy/global_solve.rs:86, src/sim/energy/global_solve.rs:590), unlike E+ inside/
     outside window surface state solution.
  2. Interior shortwave distribution is heuristic.
     Beam/diffuse split and absorptance redistribution are simplified rules, not E+ full interior solar distribution.
  3. Heavyweight mass representation is still coarse relative to E+ per-surface mass dynamics.
  4. ScriptF implementation is approximate, not full radiosity solve (src/sim/energy/view_factors.rs:452).

  Concrete patch plan (ordered)

  1. Patch A (quick): fix harness determinism and honesty.
     Make BESTEST tests pin all critical toggles explicitly (alpha, VF, solver mode) and update assertions/comments to current
     intended baseline.
  2. Patch B (medium): include steady-state exterior surfaces (windows) in global radiation enclosure.
     Add radiative nodes for ss_exterior_surfaces to global VF matrix assembly and couple each node to air/outdoor with its
     existing U/h_in approximation.
  3. Patch C (quick): fix iterative coupling bug.
     In iteration loop, use current iterate air temp for MRT and stop passing stale t_rad_guess_c in VF mode.
  4. Patch D (larger): add explicit window interior surface nodes (and optionally exterior nodes) to thermal state equations.
     This is the first real step toward EnergyPlus-like behavior.
  5. Patch E (larger): replace heuristic interior solar deposition with geometric per-surface distribution (ray/visibility based),
     then compare again.

