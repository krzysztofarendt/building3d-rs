use anyhow::Result;
use building3d::sim::energy::construction::WallConstruction;
use building3d::sim::heat_transfer::boundary::BoundaryCondition;
use building3d::sim::heat_transfer::mesh_1d::build_1d_mesh;
use building3d::sim::heat_transfer::solver::FvmWallSolver;
use building3d::sim::materials::Layer;

/// Demonstrate the FVM heat transfer solver and verify against analytical solutions.
///
/// Part 1: Steady-state multi-layer wall — compare FVM temperature profile and
///         heat flux against the resistance-network analytical solution.
/// Part 2: Transient 24-hour simulation with sinusoidal outdoor temperature.
fn main() -> Result<()> {
    // --- Wall construction (exterior to interior) ---
    let construction = WallConstruction::new(
        "insulated_exterior_wall",
        vec![
            Layer {
                name: "plaster_ext".into(),
                thickness: 0.02,
                conductivity: 0.87,
                density: 1800.0,
                specific_heat: 840.0,
            },
            Layer {
                name: "insulation".into(),
                thickness: 0.10,
                conductivity: 0.04,
                density: 30.0,
                specific_heat: 1030.0,
            },
            Layer {
                name: "concrete".into(),
                thickness: 0.15,
                conductivity: 1.4,
                density: 2300.0,
                specific_heat: 880.0,
            },
            Layer {
                name: "plaster_int".into(),
                thickness: 0.015,
                conductivity: 0.87,
                density: 1800.0,
                specific_heat: 840.0,
            },
        ],
    );

    let wall_area = 12.0; // m^2
    let h_ext = 25.0; // W/(m2*K)
    let h_int = 7.7; // W/(m2*K)

    // Resistance chain: ext film | plaster | insulation | concrete | plaster | int film
    let r_layers: Vec<f64> = construction
        .layers
        .iter()
        .map(|l| l.thickness / l.conductivity)
        .collect();
    let r_wall: f64 = r_layers.iter().sum();
    let r_total = 1.0 / h_ext + r_wall + 1.0 / h_int;

    let total_thickness: f64 = construction.layers.iter().map(|l| l.thickness).sum();

    println!("FVM Wall Heat Transfer — Analytical Verification");
    println!("{:=<60}", "");
    println!();
    println!("Wall: {}", construction.name);
    println!(
        "  Layers: {}",
        construction
            .layers
            .iter()
            .map(|l| format!("{} ({:.0}mm)", l.name, l.thickness * 1000.0))
            .collect::<Vec<_>>()
            .join(" | ")
    );
    println!("  Total thickness: {:.0} mm", total_thickness * 1000.0);
    println!("  R-wall: {r_wall:.4} m2*K/W");
    println!("  R-total (with films): {r_total:.4} m2*K/W");
    println!("  U-value: {:.3} W/(m2*K)", 1.0 / r_total);
    println!();

    // =====================================================================
    // PART 1: Steady-state verification
    // =====================================================================
    println!("PART 1: Steady-State Verification");
    println!("{:-<60}", "");
    println!();

    let t_out = -10.0;
    let t_in = 20.0;

    // Analytical solution: heat flux and interface temperatures
    let q_analytical = (t_out - t_in) / r_total; // W/m^2 (negative: heat leaves zone)

    // Temperature at each interface, walking from exterior air inward:
    //   T_ext_surface = T_out - q * (1/h_ext)    [q is negative, so surface > T_out]
    //   T_after_plaster = T_ext_surface - q * R_plaster
    //   ... etc.
    let mut t_interfaces = Vec::new();
    let mut t = t_out;
    t -= q_analytical / h_ext; // exterior surface
    t_interfaces.push(("ext surface", t));
    for (i, layer) in construction.layers.iter().enumerate() {
        t -= q_analytical * r_layers[i];
        t_interfaces.push((&layer.name, t));
    }
    // t should now equal t_in + q/h_int ≈ t_in - |q|/h_int

    println!(
        "  Boundary conditions: T_out = {t_out:.1} C, T_in = {t_in:.1} C"
    );
    println!(
        "  h_ext = {h_ext:.0} W/(m2*K), h_int = {h_int:.1} W/(m2*K)"
    );
    println!();
    println!("  Analytical heat flux: {q_analytical:.4} W/m2");
    println!();
    println!("  Analytical interface temperatures:");
    println!("    {:>20}  {:>10}", "Interface", "T [C]");
    for (name, temp) in &t_interfaces {
        println!("    {:>20}  {:>10.4}", name, temp);
    }
    println!();

    // Run FVM to steady state
    let mesh = build_1d_mesh(&construction, wall_area);
    let n_cells = mesh.cells.len();
    let mut solver = FvmWallSolver::new(mesh, 10.0);
    let sources = vec![0.0; n_cells];

    let bc_ext = BoundaryCondition::Convective {
        h: h_ext,
        t_fluid: t_out,
    };
    let bc_int = BoundaryCondition::Convective {
        h: h_int,
        t_fluid: t_in,
    };

    // Large dt drives capacity term to zero -> steady state
    let dt_ss = 1e8;
    for _ in 0..30 {
        solver.step(dt_ss, &bc_ext, &bc_int, &sources);
    }

    let q_fvm = solver.interior_heat_flux(&bc_int);
    let q_err = ((q_fvm - q_analytical) / q_analytical * 100.0).abs();

    println!("  FVM heat flux:        {q_fvm:.4} W/m2");
    println!("  Error:                {q_err:.2}%");
    println!();

    // Compare temperature profile
    // Build analytical cell-centroid temperatures by interpolating within layers
    println!(
        "  Temperature profile comparison (cell centroids):"
    );
    println!(
        "    {:>4}  {:>12}  {:>12}  {:>8}",
        "Cell", "FVM [C]", "Exact [C]", "Err [C]"
    );
    println!("    {:-<44}", "");

    let temps = solver.temperatures();
    let analytical_temps = analytical_cell_temperatures(
        &construction,
        q_analytical,
        t_out,
        h_ext,
        wall_area,
    );

    let mut max_temp_err = 0.0_f64;
    for i in 0..n_cells {
        let err = (temps[i] - analytical_temps[i]).abs();
        max_temp_err = max_temp_err.max(err);
        println!(
            "    {:>4}  {:>12.4}  {:>12.4}  {:>8.4}",
            i, temps[i], analytical_temps[i], err
        );
    }
    println!("    {:-<44}", "");
    println!("  Max temperature error: {max_temp_err:.4} C");

    let pass_ss = q_err < 1.0 && max_temp_err < 0.15;
    println!();
    if pass_ss {
        println!("  PASS: steady-state flux error < 1%, temperature error < 0.15 C");
    } else {
        println!(
            "  FAIL: flux error = {q_err:.2}%, max temp error = {max_temp_err:.4} C"
        );
    }
    println!();

    // =====================================================================
    // PART 2: Transient 24-hour simulation
    // =====================================================================
    println!("PART 2: 24-Hour Transient Simulation");
    println!("{:-<60}", "");
    println!();

    // Re-create solver with fresh initial conditions
    let mesh = build_1d_mesh(&construction, wall_area);
    let n_cells = mesh.cells.len();
    let mut solver = FvmWallSolver::new(mesh, 15.0);
    let sources = vec![0.0; n_cells];

    let dt = 60.0; // 1-minute time steps
    let warmup_days = 7; // let transient die out
    let total_hours = warmup_days * 24 + 24; // warmup + 1 day of output
    let steps_per_hour = (3600.0 / dt) as usize;
    let output_start = warmup_days * 24;

    println!(
        "  Running {warmup_days} days warmup + 24h output (dt={dt:.0}s)"
    );
    println!("  Indoor air:  {t_in:.1} C (constant)");
    println!(
        "  Outdoor air: sinusoidal, mean -5 C, amplitude 5 C"
    );
    println!();

    // Store results for the last 24 hours
    let mut hourly_results: Vec<(f64, f64, f64, f64)> = Vec::new();

    for hour in 0..total_hours {
        let t_outdoor = outdoor_temperature(hour as f64);

        let bc_ext = BoundaryCondition::Convective {
            h: h_ext,
            t_fluid: t_outdoor,
        };
        let bc_int = BoundaryCondition::Convective {
            h: h_int,
            t_fluid: t_in,
        };

        for _ in 0..steps_per_hour {
            solver.step(dt, &bc_ext, &bc_int, &sources);
        }

        if hour >= output_start {
            let q_int = solver.interior_heat_flux(&bc_int);
            hourly_results.push((
                t_outdoor,
                solver.exterior_surface_temp(),
                solver.interior_surface_temp(),
                q_int,
            ));
        }
    }

    println!(
        "  {:>4}  {:>7}  {:>7}  {:>7}  {:>9}  {:>9}",
        "Hour", "T_out", "T_ext", "T_int", "q_int", "Q_wall"
    );
    println!(
        "  {:>4}  {:>7}  {:>7}  {:>7}  {:>9}  {:>9}",
        "", "[C]", "[C]", "[C]", "[W/m2]", "[W]"
    );
    println!("  {:-<55}", "");

    for (h, (t_outdoor, t_ext, t_int_s, q_int)) in
        hourly_results.iter().enumerate()
    {
        println!(
            "  {:>4}  {:>7.1}  {:>7.2}  {:>7.2}  {:>9.2}  {:>9.1}",
            h,
            t_outdoor,
            t_ext,
            t_int_s,
            q_int,
            q_int * wall_area,
        );
    }
    println!("  {:-<55}", "");
    println!();

    // Verify: mean heat flux over the quasi-periodic day should match
    // the steady-state flux at the mean outdoor temperature
    let t_out_mean = -5.0;
    let q_mean_analytical = (t_out_mean - t_in) / r_total;
    let q_mean_fvm: f64 =
        hourly_results.iter().map(|(_, _, _, q)| q).sum::<f64>()
            / hourly_results.len() as f64;
    let q_mean_err =
        ((q_mean_fvm - q_mean_analytical) / q_mean_analytical * 100.0).abs();

    println!("  Periodic verification:");
    println!(
        "    Mean outdoor temp:       {t_out_mean:.1} C"
    );
    println!(
        "    Analytical mean flux:    {q_mean_analytical:.4} W/m2"
    );
    println!(
        "    FVM mean flux (day 2):   {q_mean_fvm:.4} W/m2"
    );
    println!("    Error:                   {q_mean_err:.2}%");
    println!();

    // Interior flux amplitude: peak-to-peak / 2
    let q_max = hourly_results
        .iter()
        .map(|(_, _, _, q)| *q)
        .fold(f64::NEG_INFINITY, f64::max);
    let q_min = hourly_results
        .iter()
        .map(|(_, _, _, q)| *q)
        .fold(f64::INFINITY, f64::min);
    let q_amplitude_fvm = (q_max - q_min) / 2.0;

    // Outdoor forcing amplitude in equivalent flux
    let forcing_amplitude = 5.0; // C
    let q_amplitude_no_lag = forcing_amplitude / r_total;

    // The wall attenuates the amplitude. For a multi-layer wall, the
    // decrement factor df = q_amplitude_actual / q_amplitude_forcing.
    let decrement_factor = q_amplitude_fvm / q_amplitude_no_lag;

    println!("    Forcing flux amplitude:  {q_amplitude_no_lag:.4} W/m2");
    println!(
        "    Interior flux amplitude: {q_amplitude_fvm:.4} W/m2"
    );
    println!("    Decrement factor:        {decrement_factor:.4}");
    println!(
        "    (df < 1 confirms thermal mass attenuates the signal)"
    );
    println!();

    // Find phase lag: hour of max outdoor temp vs hour of max interior flux
    let hour_max_tout = 12; // sin peaks at hour 12 (6 + 6)
    let hour_max_q = hourly_results
        .iter()
        .enumerate()
        .max_by(|a, b| a.1 .3.partial_cmp(&b.1 .3).unwrap())
        .map(|(i, _)| i)
        .unwrap();
    // q is negative (heat loss), so "max" flux = least negative = warmest hour
    let phase_lag = if hour_max_q >= hour_max_tout {
        hour_max_q - hour_max_tout
    } else {
        hour_max_q + 24 - hour_max_tout
    };

    println!(
        "    Peak outdoor temp at:    hour {hour_max_tout}"
    );
    println!(
        "    Peak interior flux at:   hour {hour_max_q}"
    );
    println!(
        "    Phase lag:               {phase_lag} hours"
    );
    println!(
        "    (positive lag confirms thermal inertia delays the response)"
    );
    println!();

    let pass_periodic = q_mean_err < 2.0 && decrement_factor < 1.0 && phase_lag > 0;
    if pass_periodic {
        println!(
            "  PASS: mean flux error < 2%, df < 1, phase lag > 0h"
        );
    } else {
        println!("  FAIL: check results above");
    }

    Ok(())
}

/// Sinusoidal outdoor temperature with mean -5 C, amplitude 5 C.
/// Minimum at hour 6, maximum at hour 12.
fn outdoor_temperature(hour: f64) -> f64 {
    let mean = -5.0;
    let amplitude = 5.0;
    mean + amplitude * ((hour - 6.0) * std::f64::consts::PI / 12.0).sin()
}

/// Compute analytical cell-centroid temperatures for steady-state conduction
/// through a multi-layer wall.
///
/// Uses the resistance-chain: at steady state, temperature varies linearly
/// within each layer. Each FVM cell centroid sits at the midpoint of its slab.
fn analytical_cell_temperatures(
    construction: &WallConstruction,
    q: f64,     // steady-state heat flux [W/m2] (negative = heat loss)
    t_out: f64, // outdoor air temperature
    h_ext: f64, // exterior convective coefficient
    wall_area: f64,
) -> Vec<f64> {
    let max_cell_thickness = 0.05;
    let mut temps = Vec::new();

    // Temperature at exterior wall surface
    let mut t_surface = t_out - q / h_ext;

    for layer in &construction.layers {
        let n = (layer.thickness / max_cell_thickness).ceil().max(1.0) as usize;
        let dx = layer.thickness / n as f64;
        let k = layer.conductivity;

        for j in 0..n {
            // Distance from layer start to cell centroid
            let x_centroid = (j as f64 + 0.5) * dx;
            // T(x) = T_surface_of_layer - q * x / k
            // (q is negative, so temperature increases going inward)
            let t_centroid = t_surface - q * x_centroid / k;
            temps.push(t_centroid);
        }

        // Update surface temperature for next layer
        t_surface -= q * layer.thickness / k;
    }

    // Verify: remaining temperature drop through interior film should reach t_in
    // t_surface here is the interior wall surface
    // t_in = t_surface - q / h_int (should be ~20 C)
    let _r_check = wall_area; // just to suppress unused warning

    temps
}
