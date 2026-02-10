# FVM Wall Heat Transfer -- Rerun Visualization

Visualizes the transient 1D FVM heat transfer simulation in Rerun, showing a
temperature wave propagating through a multi-layer wall over 24 hours.

## What it shows

- **Animated wall cross-section**: Each FVM cell drawn as a colored quad,
  blue (cold) to red (hot), updated every 15 minutes of simulation time.
- **Layer boundaries**: Vertical lines at material interfaces.
- **Temperature chart**: Outdoor air, exterior surface, and interior surface
  temperatures over the 24-hour output period.
- **Heat flux chart**: Interior surface heat flux (W/m^2) over time.

## Wall construction

Same 4-layer insulated exterior wall as `sim_fvm_wall`:

| Layer        | Thickness (mm) | k (W/m*K) |
|--------------|---------------:|----------:|
| plaster_ext  |             20 |      0.87 |
| insulation   |            100 |      0.04 |
| concrete     |            150 |      1.40 |
| plaster_int  |             15 |      0.87 |

## FVM discretization

Each layer is subdivided into cells no thicker than 50 mm. This produces 7
cells across the 285 mm wall:

| Layer        | Thickness (mm) | Cells | Cell size (mm) |
|--------------|---------------:|------:|---------------:|
| plaster_ext  |             20 |     1 |             20 |
| insulation   |            100 |     2 |             50 |
| concrete     |            150 |     3 |             50 |
| plaster_int  |             15 |     1 |             15 |

Between each pair of adjacent cells there is a **face** carrying a precomputed
conductance. At material interfaces the conductance uses the series-resistance
formula `K = A / (half_dx_L/k_L + half_dx_R/k_R)` so that the jump in
conductivity is handled correctly. The two boundary faces (exterior and
interior) connect the outermost cells to the convective boundary conditions.

Total: 7 cells, 6 interior faces, 2 boundary faces.

## Simulation parameters

| Parameter        | Value              |
|------------------|--------------------|
| Time step        | 60 s               |
| Warmup           | 7 days (no output) |
| Output period    | 24 h               |
| T_outdoor        | sinusoidal, mean -5 C, amplitude 5 C |
| T_indoor         | 20 C (constant)    |
| h_ext / h_int    | 25 / 7.7 W/(m^2*K)|
| Log interval     | 15 minutes (96 frames) |
| Color range      | -10 to 22 C        |

## Run

Start the Rerun viewer first, then run the example:

```bash
rerun &
cargo run --example sim_fvm_wall_viz
```
