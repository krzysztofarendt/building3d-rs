# 3D FVM Heat Transfer Example

Transient heat conduction in a solid aluminum block with an internal heat
source, solved with the 3D finite-volume method (FVM) on a tetrahedral mesh.

## Case Description

A unit cube is heated from within and cooled convectively on all six faces.
The simulation starts from a uniform temperature and tracks how the central
hot spot develops and approaches thermal equilibrium.

## Geometry

| Parameter | Value |
|-----------|-------|
| Shape | Rectangular block |
| Dimensions | 1.0 x 1.0 x 1.0 m |
| Origin | (0, 0, 0) |

## Material (Aluminum)

| Property | Value | Unit |
|----------|-------|------|
| Thermal conductivity (k) | 200 | W/(m-K) |
| Density (rho) | 2700 | kg/m^3 |
| Specific heat (cp) | 900 | J/(kg-K) |
| Thermal diffusivity (alpha) | 8.23e-5 | m^2/s |

## Grid

The mesh is generated through a three-stage pipeline:

1. **Surface mesh**: The box solid provides a triangulated surface mesh (12
   triangles on 6 faces).
2. **Point enrichment**: Extra points are added to improve mesh quality:
   - Surface points via barycentric subdivision of each triangle face
     (spacing = 0.2 m).
   - Interior points on a regular 3D grid (spacing = 0.2 m), filtered to
     include only points inside the solid.
   - Near-duplicate points are removed (tolerance = 0.02 m).
3. **Delaunay tetrahedralization**: The Bowyer-Watson algorithm produces a
   conforming tetrahedral mesh (~1300-1400 tetrahedra, ~450 vertices).
   Tetrahedra outside the surface boundary are automatically culled via
   ray casting.

Each tetrahedron becomes one FVM cell. Shared triangular faces between
tetrahedra become interior faces; unshared faces on the surface become
boundary faces (~770 boundary faces).

Interior face conductance uses the series-resistance formula:
`K = A / (d_L/k_L + d_R/k_R)`, where `d` is the centroid-to-face distance
and `k` is the cell conductivity (uniform here, so this reduces to the
harmonic mean).

## Boundary Conditions

All boundary faces receive the same convective boundary condition:

| Parameter | Value |
|-----------|-------|
| Convective coefficient (h) | 10 W/(m^2-K) |
| Fluid temperature (T_fluid) | 20 C |

This represents natural convection in still air. The effective surface
resistance is `1/h = 0.1 m^2-K/W`.

## Heat Source

A volumetric heat source is applied to cells whose centroid lies within
0.15 m of the box center (0.5, 0.5, 0.5). The total source power (5000 W)
is distributed equally among these cells (~6 cells, ~833 W each).

## Solver

- **Method**: Implicit Euler time integration with preconditioned conjugate
  gradient (PCG) linear solver (Jacobi preconditioner).
- **Time step**: dt = 1.0 s
- **Duration**: 600 steps (10 minutes simulated time)
- **Initial temperature**: 20 C (uniform, equal to fluid temperature)

The implicit scheme is unconditionally stable, allowing the relatively large
time step without numerical oscillation.

## Visualization

The results are streamed to a Rerun viewer (`localhost:9876`) with:

- **Cross-section slice**: Only tetrahedra whose centroid y-coordinate is in
  [0.4, 0.6] are rendered, giving a see-through view of the interior
  temperature field. Each cell is colored uniformly using a blue (cold) to
  white to red (hot) color ramp.
- **Box wireframe**: The original solid edges provide spatial reference.
- **Temperature time series**: Min, mean, and max cell temperatures are
  logged as scalar values on each frame.

Frames are logged every 10 time steps (61 frames total).

## Running

```bash
# Start the Rerun viewer (if not already running)
rerun

# Run the example
cargo run --example sim_fvm_3d
```

## Expected Results

The simulation shows heat diffusing outward from the center:

- **T_max** rises quickly at first, then asymptotically approaches ~38 C as
  the source power balances the convective losses.
- **T_mean** climbs slowly as energy accumulates in the block.
- **T_min** stays near 20 C (the fluid temperature) at the surfaces.

At steady state the temperature profile is governed by `q_source = h * A *
(T_surface - T_fluid)`, giving a bulk energy balance of
`5000 W = 10 * 6 * (T_avg_surface - 20)`, or `T_avg_surface ~ 103 C` for
the surface average. The actual peak is lower because the source is
concentrated and the block's high conductivity redistributes heat
efficiently.
