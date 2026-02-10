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
