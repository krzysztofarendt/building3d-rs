# FVM Heat Transfer: Implementation Plan

## Goal

Add a Finite Volume Method (FVM) heat transfer solver to building3d, starting
with 1D conduction through wall layers and designed so the same method extends
to 3D volumetric conduction on tetrahedral meshes.

## Why FVM

- **Same method in 1D and 3D.** The algorithm (cells, faces, conductances,
  implicit solve) is identical — only the mesh changes.
- **Exact energy conservation** by construction: every watt leaving one cell
  enters its neighbor.
- **Natural fit for heat transfer.** FVM is derived from the integral form of
  the conservation law, not from function approximation (FEM) or pointwise
  derivatives (FDM).
- **Simple implementation.** The diffusion equation on any mesh reduces to
  assembling a sparse symmetric matrix and solving a linear system.

## The heat equation

```
ρ c_p ∂T/∂t = ∇·(k ∇T) + Q
```

FVM discretizes this over each cell i with volume V_i:

```
(ρ c_p V_i / Δt) (T_i^{n+1} - T_i^n) = Σ_{faces f} k_f A_f / d_f (T_j^{n+1} - T_i^{n+1}) + Q_i
```

where:
- `k_f` = face conductivity (harmonic mean of neighbors)
- `A_f` = face area
- `d_f` = distance between cell centroids across face f
- `Q_i` = volumetric source term in cell i (W)

Rearranging gives a sparse linear system: `M T^{n+1} = b`.

## Module structure

```
src/sim/heat_transfer/
├── mod.rs          # Public API, trait definitions
├── mesh.rs         # FvmMesh: cells, faces, connectivity (dimension-agnostic)
├── mesh_1d.rs      # Build 1D mesh from WallConstruction layers
├── solver.rs       # FVM assembly + implicit time step (shared by 1D and 3D)
├── boundary.rs     # Boundary condition types (Dirichlet, Neumann, convective)
└── (future) mesh_3d.rs  # Build 3D mesh from TetrahedralMesh
```

## Data structures

### FvmMesh (dimension-agnostic)

```rust
pub struct FvmCell {
    pub volume: f64,            // m^3 (1D: A_wall * dx)
    pub conductivity: f64,      // W/(m·K)
    pub density: f64,           // kg/m^3
    pub specific_heat: f64,     // J/(kg·K)
}

pub struct FvmFace {
    pub cell_left: usize,       // index into cells (or BOUNDARY sentinel)
    pub cell_right: usize,      // index into cells (or BOUNDARY sentinel)
    pub area: f64,              // m^2 (1D: A_wall, the polygon area)
    pub distance: f64,          // m (centroid-to-centroid distance across face)
    pub conductance: f64,       // W/K = k_f * A / d (precomputed)
}

pub struct FvmMesh {
    pub cells: Vec<FvmCell>,
    pub faces: Vec<FvmFace>,
}
```

### Boundary conditions

```rust
pub enum BoundaryCondition {
    /// Fixed temperature (T = T_prescribed)
    Dirichlet { temperature: f64 },
    /// Fixed heat flux (q = q_prescribed, W/m^2)
    Neumann { heat_flux: f64 },
    /// Convective (q = h * (T_surface - T_fluid))
    Convective { h: f64, t_fluid: f64 },
}
```

Convective BCs are the primary mode — exterior surfaces see outdoor air with
h_out (~10-25 W/(m^2·K)), interior surfaces see zone air with h_in (~3-8).

### Solver interface

```rust
pub trait HeatTransferSolver {
    /// Advance one time step. Returns new cell temperatures.
    fn step(
        &mut self,
        dt: f64,
        bc_exterior: &BoundaryCondition,
        bc_interior: &BoundaryCondition,
        sources: &[f64],           // per-cell volumetric source (W)
    ) -> &[f64];

    /// Interior surface temperature (last cell).
    fn interior_surface_temp(&self) -> f64;

    /// Exterior surface temperature (first cell).
    fn exterior_surface_temp(&self) -> f64;

    /// Heat flux into zone through interior surface (W/m^2).
    fn interior_heat_flux(&self) -> f64;
}
```

This is the interface the zone thermal model calls. It doesn't know if the
solver is 1D or 3D.

## Step-by-step implementation

### Step 1: FvmMesh and 1D mesh builder ✅

Build a 1D mesh from a `WallConstruction`:

```
Exterior BC ─┤ cell_0 │ cell_1 │ ... │ cell_N ├─ Interior BC
              layer 0    layer 1        layer M
```

Each `Layer` is subdivided into one or more cells (at least one per layer;
thin layers get one cell, thick layers can be split for accuracy). The face
conductance between cells in different layers uses the **series-resistance
formula**:

```
K = A / (half_dx_L/k_L + half_dx_R/k_R)
```

(reduces to harmonic mean when half-distances are equal).

Input: `WallConstruction` (already has `Vec<Layer>` with thickness, k, ρ, c_p)
and wall polygon area.

Output: `FvmMesh` with cells and faces.

Implemented in `mesh.rs` and `mesh_1d.rs`. Tests:
- `test_single_layer_mesh`: 0.20m concrete → 4 cells, 5 faces
- `test_three_layer_mesh`: layer subdivision and total thickness
- `test_thin_layer_gets_one_cell`: 0.001m layer → 1 cell
- `test_face_conductance_same_material`: K ≈ 28.0 for uniform material

### Step 2: Implicit FVM solver ✅

Assemble and solve the linear system using Backward Euler:

```
(C/Δt + K) T^{n+1} = (C/Δt) T^n + Q + BC contributions
```

In 1D this is a **tridiagonal system** — solved with Thomas algorithm (O(N)).

Implemented in `solver.rs`. Tests:
- `test_thomas_algorithm`: 3×3 system
- `test_steady_state_dirichlet`: linear profile and q = k·ΔT/L

### Step 3: Boundary condition handling ✅

Four BC types implemented in `boundary.rs`:

- **Dirichlet**: modify the matrix row for the boundary cell to fix temperature.
- **Neumann**: add the prescribed flux to the RHS of the boundary cell.
- **Convective**: series resistance `1/(h·A) + 1/K_face` for half-cell accuracy.
- **ConvectiveWithFlux**: convective + imposed surface source (e.g. absorbed
  solar), with flux split between conduction and convection via
  `boundary_conduction_fraction()`.

Tests:
- `test_steady_state_convective`: q = ΔT/R_total within 1%
- `test_convective_multilayer_interface_temps`: surface temps with half-cell offset
- `test_neumann_flux_reporting_sign`: sign convention verified

### Step 4: Validation against analytical solutions ✅

1. **Steady-state multi-layer wall** ✅: 3 layers with different k, convective
   BCs. Test: `test_steady_state_multilayer`.

2. **Transient step response** ✅: semi-infinite solid with sudden surface
   temperature change. Compared to `T(x,t) = T_s·erfc(x/(2√(αt)))`.
   Test: `test_transient_step_response`.

3. **Periodic forcing** ✅: 24h sinusoidal Dirichlet BC on a 1m concrete wall.
   Verified against `T(x,t) = T_mean + A₀·exp(-x/δ)·sin(ωt - x/δ)` where
   `δ = √(2α/ω)` is the penetration depth.
   Test: `test_periodic_forcing` — checks amplitude decay (<5% error),
   phase lag (<500s absolute), and mean temperature at 5 cell depths.

### Step 5: Integration with zone thermal model ✅

Connected FVM wall solver to the energy simulation in both `module.rs`
(SimModule pipeline) and `simulation.rs` (annual simulation):

1. `ThermalConfig.use_fvm_walls: bool` (default: **true**). Eligible surfaces:
   exterior, resolved `WallConstruction`, not glazing, not U-value override,
   not ground-coupled.
2. Each timestep: convective BCs from zone air / outdoor temps, absorbed solar
   via `ConvectiveWithFlux`, interior heat flux fed to zone energy balance.
3. `ThermalNetwork::build_with_ignored_exterior_polygons()` excludes FVM
   surfaces from steady-state UA to avoid double-counting.
4. FVM wall capacity deducted from zone lumped capacity.

### Step 6: BESTEST validation ✅

BESTEST 600 and 900 run with FVM walls enabled by default. The
`bestest_energy_suite` example supports `BESTEST_USE_FVM_WALLS` env var.

### Step 7 (future): 3D mesh builder

Build a 3D `FvmMesh` from a `TetrahedralMesh` (already exists in the codebase
at `src/geom/mesh/tetrahedralize.rs`):

- Cells = tetrahedra (volume from `tetrahedron_volume()`)
- Faces = shared triangle faces between adjacent tetrahedra
- Face area = triangle area
- Distance = centroid-to-centroid distance
- Conductance = `k_f * A_f / d_f` (harmonic mean at material interfaces)

Build adjacency from the `TetrahedronIndex` list: two tetrahedra are neighbors
if they share exactly 3 vertex indices.

The solver (`solver.rs`) is reused without changes — only the matrix is no
longer tridiagonal, so swap the Thomas algorithm for a sparse CG solver.

## Dependencies

- No external crates needed for 1D (Thomas algorithm is trivial)
- For 3D sparse solve: either implement a simple CG, or add `nalgebra-sparse`
  (evaluate later)

## Key design decisions

1. **FvmMesh is dimension-agnostic.** The solver sees cells and faces, never
   coordinates. This is what makes 1D → 3D seamless.

2. **One solver per surface.** Each exterior polygon gets its own 1D solver
   instance. This is embarrassingly parallel and avoids coupling complexity
   (lateral conduction between surfaces is negligible for most buildings).

3. **Implicit time stepping only.** Backward Euler is unconditionally stable
   at any Δt. No CFL constraint, compatible with the hourly timestep of the
   energy simulation. Crank-Nicolson can be added later as an option for
   second-order accuracy.

4. **Series-resistance formula for face conductivity.** At interfaces between
   layers with different k, `K = A / (half_dx_L/k_L + half_dx_R/k_R)` ensures
   correct steady-state heat flux. Reduces to harmonic mean when cell sizes
   are equal. This is standard for FVM on heterogeneous media.

5. **Concrete solver type.** `FvmWallSolver` is used directly (no trait
   abstraction). The zone model calls `step()` and reads
   `interior_heat_flux()`. A trait can be added later if alternative solver
   implementations are needed.
