use anyhow::Result;
use rerun as rr;

use building3d::draw::config::RerunConfig;
use building3d::draw::rerun::start_session;
use building3d::sim::energy::construction::WallConstruction;
use building3d::sim::heat_transfer::boundary::BoundaryCondition;
use building3d::sim::heat_transfer::mesh_1d::build_1d_mesh;
use building3d::sim::heat_transfer::solver::FvmWallSolver;
use building3d::sim::materials::Layer;

/// Visualize FVM wall heat transfer in Rerun.
///
/// Shows an animated wall cross-section colored by temperature, plus
/// temperature and heat flux time-series charts over a 24-hour period.
///
/// Requires a running Rerun viewer (`rerun`).
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
    let t_in = 20.0; // C

    // --- Build FVM mesh and compute cell geometry ---
    let mesh = build_1d_mesh(&construction, wall_area);
    let n_cells = mesh.cells.len();

    // Compute cell x-positions (start, end) for drawing quads
    let cell_dx: Vec<f64> = mesh.cells.iter().map(|c| c.volume / wall_area).collect();
    let mut cell_x_start = Vec::with_capacity(n_cells);
    let mut x = 0.0;
    for dx in &cell_dx {
        cell_x_start.push(x);
        x += dx;
    }

    // Compute layer boundary x-positions for vertical lines
    let mut layer_boundaries = Vec::new();
    let mut lx = 0.0;
    for layer in &construction.layers {
        lx += layer.thickness;
        layer_boundaries.push(lx);
    }
    // Remove last boundary (it's the wall end, not an interior interface)
    layer_boundaries.pop();

    // --- Start Rerun session ---
    let mut config = RerunConfig::new();
    config.session_name = "FVM Wall Heat Transfer".to_string();
    let session = start_session(&config)?;

    // --- Log static elements ---

    // Layer boundary lines (vertical lines at material interfaces)
    let quad_h: f64 = 0.5; // height of cell quads
    if !layer_boundaries.is_empty() {
        let lines: Vec<Vec<rr::Vec3D>> = layer_boundaries
            .iter()
            .map(|&bx| {
                vec![
                    rr::Vec3D([bx as f32, 0.0, -0.02]),
                    rr::Vec3D([bx as f32, 0.0, quad_h as f32 + 0.02]),
                ]
            })
            .collect();
        session.log_static(
            "wall/boundaries",
            &rr::LineStrips3D::new(lines)
                .with_colors(vec![
                    rr::Color::from_rgb(180, 180, 180);
                    layer_boundaries.len()
                ])
                .with_radii(vec![rr::Radius::new_ui_points(1.0); layer_boundaries.len()]),
        )?;
    }

    // Configure time-series appearance (static, logged once)
    session.log_static(
        "temperature/outdoor",
        &rr::SeriesLines::new()
            .with_colors([[80, 80, 255]])
            .with_names(["T outdoor (C)"])
            .with_widths([2.0]),
    )?;
    session.log_static(
        "temperature/ext_surface",
        &rr::SeriesLines::new()
            .with_colors([[255, 165, 0]])
            .with_names(["T ext surface (C)"])
            .with_widths([1.5]),
    )?;
    session.log_static(
        "temperature/int_surface",
        &rr::SeriesLines::new()
            .with_colors([[255, 80, 80]])
            .with_names(["T int surface (C)"])
            .with_widths([1.5]),
    )?;
    session.log_static(
        "heat_flux/interior",
        &rr::SeriesLines::new()
            .with_colors([[200, 50, 50]])
            .with_names(["q interior (W/m2)"])
            .with_widths([2.0]),
    )?;

    // --- Run simulation ---
    let mut solver = FvmWallSolver::new(mesh, 15.0);
    let sources = vec![0.0; n_cells];
    let dt = 60.0; // 1-minute time steps
    let warmup_days = 7;
    let steps_per_hour = (3600.0 / dt) as usize;

    // Temperature color range
    let t_min = -10.0_f64;
    let t_max = 22.0_f64;

    // Warmup: run silently (no Rerun logging)
    println!("Running {warmup_days}-day warmup...");
    for hour in 0..(warmup_days * 24) {
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
    }

    // Output day: log every 15 minutes (96 frames over 24h)
    println!("Logging 24-hour output to Rerun...");
    let log_interval_minutes = 15;
    let steps_per_log = log_interval_minutes * 60 / dt as usize;
    let output_start_hour = warmup_days * 24;
    let mut minute = 0_i64;
    let mut step_in_output = 0_usize;

    for hour_offset in 0..24 {
        let hour = (output_start_hour + hour_offset) as f64;
        let t_outdoor = outdoor_temperature(hour);
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
            step_in_output += 1;

            if step_in_output.is_multiple_of(steps_per_log) {
                session.set_time_sequence("minute", minute);

                // Log wall cross-section: one Mesh3D per cell
                let temps = solver.temperatures();
                for i in 0..n_cells {
                    let xs = cell_x_start[i] as f32;
                    let xe = (cell_x_start[i] + cell_dx[i]) as f32;
                    let h = quad_h as f32;

                    let vertices: Vec<rr::Vec3D> = vec![
                        rr::Vec3D([xs, 0.0, 0.0]),
                        rr::Vec3D([xe, 0.0, 0.0]),
                        rr::Vec3D([xe, 0.0, h]),
                        rr::Vec3D([xs, 0.0, h]),
                    ];
                    let indices: Vec<rr::TriangleIndices> = vec![
                        rr::TriangleIndices(rr::datatypes::UVec3D([0, 1, 2])),
                        rr::TriangleIndices(rr::datatypes::UVec3D([0, 2, 3])),
                    ];

                    let t_norm = ((temps[i] - t_min) / (t_max - t_min)).clamp(0.0, 1.0) as f32;
                    let (r, g, b) = heat_loss_color(t_norm);

                    session.log(
                        format!("wall/cell_{i}"),
                        &rr::Mesh3D::new(vertices)
                            .with_triangle_indices(indices)
                            .with_albedo_factor(rr::Rgba32::from_linear_unmultiplied_rgba_f32(
                                r, g, b, 1.0,
                            )),
                    )?;
                }

                // Log temperature time series
                session.log("temperature/outdoor", &rr::Scalars::single(t_outdoor))?;
                session.log(
                    "temperature/ext_surface",
                    &rr::Scalars::single(solver.exterior_surface_temp()),
                )?;
                session.log(
                    "temperature/int_surface",
                    &rr::Scalars::single(solver.interior_surface_temp()),
                )?;

                // Log heat flux time series
                let q_int = solver.interior_heat_flux(&bc_int);
                session.log("heat_flux/interior", &rr::Scalars::single(q_int))?;

                minute += log_interval_minutes as i64;
            }
        }
    }

    println!("Done. View results in the Rerun viewer.");
    Ok(())
}

/// Sinusoidal outdoor temperature with mean -5 C, amplitude 5 C.
fn outdoor_temperature(hour: f64) -> f64 {
    let mean = -5.0;
    let amplitude = 5.0;
    mean + amplitude * ((hour - 6.0) * std::f64::consts::PI / 12.0).sin()
}

/// Maps a normalized temperature value (0=cold, 1=hot) to a blue-white-red ramp.
fn heat_loss_color(t: f32) -> (f32, f32, f32) {
    if t < 0.5 {
        let s = t * 2.0;
        (s, s, 1.0) // blue -> white
    } else {
        let s = (t - 0.5) * 2.0;
        (1.0, 1.0 - s, 1.0 - s) // white -> red
    }
}
