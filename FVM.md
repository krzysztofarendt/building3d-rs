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

### Step 1: FvmMesh and 1D mesh builder

Build a 1D mesh from a `WallConstruction`:

```
Exterior BC ─┤ cell_0 │ cell_1 │ ... │ cell_N ├─ Interior BC
              layer 0    layer 1        layer M
```

Each `Layer` is subdivided into one or more cells (at least one per layer;
thin layers get one cell, thick layers can be split for accuracy). The face
conductance between cells in different layers uses the **harmonic mean**:

```
k_face = 2 * k_L * k_R / (k_L + k_R)
```

Input: `WallConstruction` (already has `Vec<Layer>` with thickness, k, ρ, c_p)
and wall polygon area.

Output: `FvmMesh` with cells and faces.

Deliverable: unit test — build a 3-layer wall mesh, verify cell count, face
count, and total thickness.

### Step 2: Implicit FVM solver

Assemble and solve the linear system using Backward Euler:

```
(C/Δt + K) T^{n+1} = (C/Δt) T^n + Q + BC contributions
```

where:
- `C` = diagonal capacity matrix: `C_ii = ρ_i c_p_i V_i`
- `K` = stiffness matrix from face conductances (symmetric, tridiagonal in 1D)
- `Q` = source vector
- BC contributions modify the relevant rows

In 1D this is a **tridiagonal system** — solve with Thomas algorithm (O(N)).
In 3D it becomes a sparse symmetric positive-definite system — solve with
conjugate gradient (future).

Use the Thomas algorithm for 1D first. When extending to 3D, swap in a sparse
solver (e.g. `nalgebra` sparse, or a simple CG implementation).

Deliverable: unit test — single-layer wall, steady-state analytical solution.
Apply T_ext = 0, T_int = 20 (Dirichlet), verify linear temperature profile
and correct heat flux (q = k * ΔT / L).

### Step 3: Boundary condition handling

Implement the three BC types:

- **Dirichlet**: modify the matrix row for the boundary cell to fix temperature.
- **Neumann**: add the prescribed flux to the RHS of the boundary cell.
- **Convective**: add `h * A` to the diagonal and `h * A * T_fluid` to the RHS.
  This is the most common mode in building simulation.

Deliverable: unit test — convective BCs on both sides, verify that steady-state
matches the analytical solution: `q = (T_out - T_in) / (1/h_out + R_wall + 1/h_in)`.

### Step 4: Validation against analytical solutions

1. **Steady-state multi-layer wall**: 3 layers with different k, convective BCs.
   Verify interface temperatures and heat flux against hand calculation.

2. **Transient step response**: semi-infinite solid with sudden surface
   temperature change. Compare to error-function analytical solution
   `T(x,t) = T_s * erfc(x / (2√(αt)))` at selected times.

3. **Periodic forcing**: sinusoidal outdoor temperature on a concrete wall.
   Verify amplitude decay and phase lag against analytical solution for
   periodic conduction through a slab.

### Step 5: Integration with zone thermal model

Connect the FVM wall solver to the existing energy simulation:

1. For each exterior polygon with a `WallConstruction`, instantiate an
   `FvmWallSolver` (1D mesh from layers, area from polygon).
2. Each simulation timestep:
   - Set exterior BC: `Convective { h: h_out, t_fluid: T_outdoor }` (or sol-air)
   - Set interior BC: `Convective { h: h_in, t_fluid: T_zone_air }`
   - Call `solver.step(dt, bc_ext, bc_int, &sources)`
   - Read back `interior_heat_flux()` → contributes to zone energy balance
   - Read back `interior_surface_temp()` → available for MRT/comfort (future)
3. The per-surface heat flux replaces the current `U * A * ΔT` steady-state
   calculation, giving proper thermal lag and surface temperatures.

This step requires a new `ThermalConfig` option to enable FVM walls (e.g.
`use_fvm_walls: bool`). The existing U-value path remains as fallback for
surfaces without layer definitions.

### Step 6: BESTEST validation

Run Cases 600 and 900 with FVM walls enabled. Compare:
- Annual heating/cooling totals (should stay within tolerance)
- **Monthly shape** (should improve, especially for Case 900 where thermal
  mass dynamics matter)
- Interior surface temperatures (new diagnostic output)
- Peak load timing (should improve with proper thermal lag)

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

4. **Harmonic mean for face conductivity.** At interfaces between layers with
   different k, the harmonic mean ensures correct steady-state heat flux.
   This is standard for FVM on heterogeneous media.

5. **The trait boundary (`HeatTransferSolver`) hides the solver.** The zone
   model calls `step()` and reads `interior_heat_flux()`. It never knows
   whether the solver is 1D FVM, 3D FVM, or a future FEM implementation.
