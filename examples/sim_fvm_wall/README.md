# FVM Wall Heat Transfer -- Analytical Verification

This example demonstrates the implicit Finite Volume Method (FVM) solver for 1D
heat conduction through a multi-layer wall and verifies it against closed-form
analytical solutions.

## Wall Construction

A four-layer insulated exterior wall (exterior to interior):

| Layer        | Thickness (mm) | Conductivity k (W/m*K) | Density (kg/m^3) | Specific heat (J/kg*K) |
|--------------|---------------:|----------------------:|------------------:|-----------------------:|
| plaster_ext  |             20 |                  0.87 |              1800 |                    840 |
| insulation   |            100 |                  0.04 |                30 |                   1030 |
| concrete     |            150 |                  1.40 |              2300 |                    880 |
| plaster_int  |             15 |                  0.87 |              1800 |                    840 |

Total thickness: 285 mm. Wall area: 12 m^2.

## Boundary Conditions

Convective (Robin) boundary conditions on both surfaces:

| Boundary | Coefficient h (W/m^2*K) |
|----------|------------------------:|
| Exterior |                    25.0 |
| Interior |                     7.7 |

## Part 1: Steady-State Verification

**Conditions**: T_out = -10 C, T_in = 20 C (constant).

### Analytical Solution -- Resistance Network

At steady state the heat flux through the wall is uniform and determined by the
total thermal resistance of the series circuit:

```
R_total = 1/h_ext + SUM(dx_i / k_i) + 1/h_int
```

where dx_i and k_i are the thickness and conductivity of layer i.

The steady-state heat flux per unit area:

```
q = (T_out - T_in) / R_total    [W/m^2]
```

The temperature at any point inside the wall follows from the linear temperature
profile within each homogeneous layer. Walking inward from the exterior air:

```
T_ext_surface = T_out - q / h_ext
T(x)          = T_layer_start - q * x / k_layer
```

where x is the distance from the start of the current layer.

### FVM Approach

Each material layer is subdivided into cells of at most 50 mm. The solver uses
a fully implicit (backward Euler) time discretization. To reach steady state,
a very large time step (dt = 10^8 s) is used for 30 iterations, which drives
the transient capacity term to zero and yields the stationary solution.

### Acceptance Criteria

- Heat flux error < 1%
- Maximum cell-centroid temperature error < 0.15 C

## Part 2: Transient 24-Hour Simulation

**Conditions**: T_in = 20 C (constant), T_out varies sinusoidally.

### Outdoor Temperature Forcing

```
T_out(t) = -5 + 5 * sin(pi * (t - 6) / 12)    [C]
```

Mean = -5 C, amplitude = 5 C, minimum at hour 6 (00 h shifted), maximum at
hour 12.

### Simulation Parameters

| Parameter        | Value    |
|------------------|----------|
| Time step (dt)   | 60 s     |
| Warmup period    | 7 days   |
| Output period    | 24 h     |
| Initial temp     | 15 C     |

The 7-day warmup lets transients from the uniform initial condition die out so
that the final 24-hour output represents a quasi-periodic steady state.

### Analytical Checks

**Mean flux**: Over one full period of a sinusoidal forcing, the mean heat flux
must equal the steady-state flux at the mean outdoor temperature:

```
q_mean = (T_out_mean - T_in) / R_total
```

**Decrement factor**: Thermal mass attenuates the amplitude of the flux
oscillation at the interior surface relative to the forcing amplitude:

```
df = q_amplitude_interior / q_amplitude_forcing
```

where q_amplitude_forcing = T_amplitude / R_total. A value df < 1 confirms that
the wall's thermal mass dampens the signal.

**Phase lag**: The peak interior heat flux occurs later than the peak outdoor
temperature due to thermal inertia. A positive phase lag (in hours) confirms
that the solver correctly models heat storage in the wall.

### Acceptance Criteria

- Mean flux error < 2%
- Decrement factor < 1 (thermal mass attenuates)
- Phase lag > 0 hours (thermal inertia delays the response)

## Run

```bash
cargo run --example sim_fvm_wall
```
