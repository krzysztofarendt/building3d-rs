# FVM Heat Transfer: Current Implementation

This document describes the implemented Finite Volume Method (FVM) heat-transfer
stack in `building3d` as it exists in code today.

## Scope

Implemented:
- 1D transient conduction through layered wall constructions.
- 3D transient conduction on tetrahedral meshes.
- Implicit solvers for both dense 1D tridiagonal and generic sparse topologies.
- Integration of 1D per-surface FVM walls into energy simulations.
- Integration of optional internal mass slabs (also 1D FVM) into energy simulations.

Not currently wired into building energy runtime:
- 3D FVM is available as mesh+solver APIs and example code, but energy simulation
  integration uses 1D per-surface wall solvers.

## Governing Discretization

Heat equation:

```
rho * cp * dT/dt = div(k * grad(T)) + Q
```

Cell-centered implicit FVM update:

```
(C_i / dt) * (T_i^{n+1} - T_i^n) = sum_f K_f * (T_nb^{n+1} - T_i^{n+1}) + S_i
```

where:
- `C_i = rho_i * cp_i * V_i` [J/K]
- `K_f` is face conductance [W/K]
- `S_i` is per-cell source power [W]

After assembly, this gives a linear system solved each timestep.

## Module Layout

```
src/sim/heat_transfer/
├── boundary.rs       # BoundaryCondition enum
├── mesh.rs           # FvmCell, FvmFace, FvmMesh, BOUNDARY sentinel
├── mesh_1d.rs        # build_1d_mesh() from WallConstruction
├── mesh_3d.rs        # build_3d_mesh(), build_3d_mesh_uniform()
├── solver.rs         # FvmWallSolver (1D tridiagonal + Thomas)
└── solver_sparse.rs  # FvmSparseSolver (generic sparse + PCG)
```

Public exports are re-exported from `src/sim/heat_transfer/mod.rs`.

## Core Data Model

`FvmMesh` is dimension-agnostic:
- `FvmCell { volume, conductivity, density, specific_heat }`
- `FvmFace { cell_left, cell_right, area, distance, conductance }`
- boundary side is marked by `BOUNDARY = usize::MAX`

Face conductance is precomputed and stored in each face.

## Boundary Conditions

`BoundaryCondition` currently supports:
- `Dirichlet { temperature }`
- `Neumann { heat_flux }` where heat flux is positive into domain
- `Convective { h, t_fluid }`
- `ConvectiveWithFlux { h, t_fluid, heat_flux }`
- `ConvectiveWithFluxToDomain { h, t_fluid, heat_flux }`

`ConvectiveWithFlux` behavior:
- surface source is split between convection and conduction using
  `alpha = K_face / (K_face + h*A)`
- only `alpha * heat_flux` enters the domain

`ConvectiveWithFluxToDomain` behavior:
- full imposed flux enters the domain (no split)

## 1D Wall Meshing

`build_1d_mesh(construction, wall_area)` in `mesh_1d.rs`:
- layers are ordered exterior -> interior
- each layer is subdivided with max cell thickness `0.05 m`
- each layer gets at least one cell
- interior interface conductance uses series resistance:

```
K = A / (half_dx_L / k_L + half_dx_R / k_R)
```

- boundary faces use half-cell distance:
  `K = k * A / half_dx`

## 1D Solver (`FvmWallSolver`)

`solver.rs` implements:
- backward-Euler implicit time stepping
- tridiagonal matrix assembly
- Thomas algorithm solve (O(N))

API highlights:
- `step(dt, bc_exterior, bc_interior, sources)`
- `temperatures()`
- `total_capacity_j_per_k()`
- `interior_heat_flux(...)`, `exterior_heat_flux(...)`
- `interior_surface_temperature(...)`, `exterior_surface_temperature(...)`
- `boundary_conduction_fraction(h, is_interior)`

Notes:
- `interior_surface_temp()` / `exterior_surface_temp()` return boundary-adjacent
  cell centroid temperatures.
- `interior_surface_temperature(...)` / `exterior_surface_temperature(...)`
  reconstruct actual boundary surface temperatures using half-cell resistance.

## 3D Meshing (`mesh_3d.rs`)

`build_3d_mesh_uniform()` and `build_3d_mesh()` convert a tetrahedral mesh to `FvmMesh`:
- one FVM cell per tetrahedron
- shared tetra faces become interior faces
- unmatched tetra faces become boundary faces

Interior face conductance:

```
K = A / (d_L / k_L + d_R / k_R)
```

Boundary face conductance:

```
K = k * A / d
```

with `d` equal to centroid-to-face distance of the adjacent tetrahedron.

The builder validates:
- property count matches tetra count
- positive tetra volumes / face areas / distances
- non-manifold face sharing is rejected

## Generic Sparse Solver (`FvmSparseSolver`)

`solver_sparse.rs` provides a topology-agnostic implicit solver for `FvmMesh`:
- sparse assembly (`diag + adjacency off-diagonals`)
- PCG solve with Jacobi preconditioning

Default config (`SparseSolverConfig`):
- `max_iterations = 2000`
- `rel_tolerance = 1e-9`
- `abs_tolerance = 1e-12`

Boundary handling:
- BCs are passed as `(face_index, BoundaryCondition)` pairs
- omitted boundary faces are adiabatic
- applying BC on interior faces is rejected

## Energy Integration (Implemented)

1D FVM is integrated in both:
- `src/sim/energy/module.rs` (step-based `EnergyModule`)
- `src/sim/energy/simulation.rs` (`run_transient_simulation*`, multizone transient)

Key behavior:
- `ThermalConfig.use_fvm_walls` defaults to `true`
- eligible FVM exterior surfaces are:
  exterior + non-glazing + no explicit U-value override + resolvable construction + not ground-coupled
- each eligible polygon gets its own `FvmWallSolver`
- exterior absorbed shortwave/sol-air terms are applied via `ConvectiveWithFlux`
- interior surface sources (transmitted solar/radiant gains where configured) use
  `ConvectiveWithFluxToDomain`
- wall-to-zone contribution is convective flux at interior surface
- FVM polygons are excluded from steady UA network via
  `ThermalNetwork::build_with_ignored_exterior_polygons(...)` to avoid double counting
- explicit FVM/internal-mass capacities are subtracted from lumped zone capacity where needed

Related configuration flags:
- `use_fvm_walls`
- `internal_mass_surfaces`
- `distribute_transmitted_solar_to_fvm_walls`
- `use_surface_aware_solar_distribution`
- `use_interior_radiative_exchange`

## Tests and Validation Coverage

Implemented test coverage includes:
- `mesh_1d.rs`: layer subdivision, thickness conservation, conductance checks
- `solver.rs`: Thomas solve, steady-state Dirichlet/convective/multilayer, transient
  step response, periodic forcing, Neumann sign convention, convective-with-flux variants
- `mesh_3d.rs`: tetra-to-FVM topology, heterogeneous conductance, manifold checks
- `solver_sparse.rs`: generic BC behavior, flux/surface temperature reconstruction,
  3D mesh-builder compatibility, PCG-driven solves
- energy integration tests in `module.rs` and `simulation.rs`, including FVM-specific behavior
- BESTEST suite paths with FVM enabled by default (`BESTEST_USE_FVM_WALLS`)

## Examples

- `cargo run --example sim_fvm_wall`
- `cargo run --example sim_fvm_wall_viz`
- `cargo run --example sim_fvm_3d`
