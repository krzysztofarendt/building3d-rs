use anyhow::Result;
use building3d::draw::acoustics::{draw_impulse_response, draw_receivers};
use building3d::draw::rerun::{draw_simulation, start_session};
use building3d::sim::acoustics::impulse_response::ImpulseResponse;
use building3d::sim::acoustics::metrics::RoomAcousticReport;
use building3d::sim::acoustics::receiver::Receiver;
use building3d::sim::materials::{
    AcousticMaterial, Material, MaterialLibrary, OCTAVE_BAND_FREQUENCIES,
};
use building3d::sim::rays::{AcousticMode, Simulation, SimulationConfig};
use building3d::{Building, FloorPlan, Point, RerunConfig, Solid, Zone};

/// Build an L-shaped room (single solid) in one zone.
///
/// ```text
///     (0,8)-----(4,8)
///       |         |
///       |         |
///       |         |
///     (0,4)-----(4,4)-----(6,4)
///       |                   |
///       |                   |
///       |                   |
///     (0,0)---------------(6,0)
/// ```
fn build_l_shaped() -> Result<Building> {
    let room = Solid::from_floor_plan(FloorPlan {
        plan: vec![
            (0.0, 0.0),
            (6.0, 0.0),
            (6.0, 4.0),
            (4.0, 4.0),
            (4.0, 8.0),
            (0.0, 8.0),
        ],
        height: 3.0,
        name: "room".to_string(),
        ..Default::default()
    })?;

    let zone = Zone::new("zone", vec![room])?;
    Building::new("L-building", vec![zone])
}

fn main() -> Result<()> {
    let building = build_l_shaped()?;

    // Set up material library with presets
    let mut lib = MaterialLibrary::with_presets();

    // Add a custom absorptive material for the ceiling
    lib.add(
        Material::new("acoustic_ceiling").with_acoustic(AcousticMaterial::new(
            "acoustic_ceiling",
            [0.15, 0.20, 0.30, 0.45, 0.55, 0.60],
            [0.30, 0.30, 0.30, 0.30, 0.30, 0.30],
        )),
    );

    // Assign materials to surfaces
    lib.assign("wall", "concrete"); // All walls get concrete
    lib.assign("floor", "carpet"); // Floors get carpet
    lib.assign("ceiling", "acoustic_ceiling"); // Ceilings get absorptive panels

    // Configure acoustic simulation
    let mut config = SimulationConfig::new();
    // Keep all surfaces; this example uses a single (non-convex) solid with no stitched interfaces.
    config.search_transparent = false;
    config.acoustic_mode = AcousticMode::FrequencyDependent;
    config.material_library = Some(lib);
    config.enable_air_absorption = true;

    // Source in lower part of the L-shape
    config.source = Point::new(3.0, 2.0, 1.5);
    // Absorber (receiver location) in upper part of the L-shape
    config.absorbers = vec![Point::new(2.0, 6.0, 1.5)];
    config.absorber_radius = 0.5;

    config.num_rays = 5000;
    config.num_steps = 2000;
    config.store_ray_history = false; // saves ~800 MB RAM

    println!("Running frequency-dependent acoustic ray tracing...");
    println!(
        "  Source: ({:.1}, {:.1}, {:.1})",
        config.source.x, config.source.y, config.source.z
    );
    println!(
        "  Receiver: ({:.1}, {:.1}, {:.1})",
        config.absorbers[0].x, config.absorbers[0].y, config.absorbers[0].z
    );
    println!("  Rays: {}", config.num_rays);
    println!("  Steps: {}", config.num_steps);
    println!();

    // Run simulation
    let sim = Simulation::new(&building, config)?;
    let result = sim.run();

    // Print absorber hit summary
    let total_hits: f64 = result.hits.iter().map(|h| h[0]).sum();
    println!(
        "Simulation complete: {} steps, {} rays",
        result.config.num_steps, result.config.num_rays
    );
    println!("Total energy absorbed: {:.4}", total_hits);
    println!();

    // Create receiver at the absorber location
    let max_time = result.config.num_steps as f64 * result.config.time_step;
    let mut receiver = Receiver::new(
        result.config.absorbers[0],
        result.config.absorber_radius,
        result.config.time_step,
        max_time,
    );

    // Feed hits from simulation result into receiver
    if let Some(ref band_hits) = result.band_hits {
        for (step, step_band_hits) in band_hits.iter().enumerate() {
            let time = step as f64 * result.config.time_step;
            // Feed per-absorber band hits (absorber index 0)
            if !step_band_hits.is_empty() {
                let energy_sum: f64 = step_band_hits[0].iter().sum();
                if energy_sum > 0.0 {
                    receiver.record_band_hit(time, &step_band_hits[0]);
                }
            }
        }
    }

    println!("Receiver total energy: {:.6}", receiver.total_energy());
    println!();

    // Extract impulse response
    let ir = ImpulseResponse::from_receiver(&receiver);
    println!(
        "Impulse response: {} bins, {:.3} ms resolution",
        ir.len(),
        ir.time_resolution * 1000.0
    );
    println!("Total IR energy: {:.6}", ir.total_energy());
    println!();

    // Compute room acoustic metrics
    let report = RoomAcousticReport::from_ir(&ir);

    println!("Room Acoustic Metrics:");
    println!("{:-<60}", "");
    println!(
        "{:>8} {:>10} {:>10} {:>10} {:>10}",
        "Freq", "RT60 (s)", "EDT (s)", "C80 (dB)", "D50"
    );
    println!("{:-<60}", "");

    for (band, &freq) in OCTAVE_BAND_FREQUENCIES.iter().enumerate() {
        let rt60_str = report.rt60[band]
            .map(|v| format!("{:.3}", v))
            .unwrap_or_else(|| "  ---".to_string());
        let edt_str = report.edt[band]
            .map(|v| format!("{:.3}", v))
            .unwrap_or_else(|| "  ---".to_string());
        let c80_str = report.c80[band]
            .map(|v| format!("{:.1}", v))
            .unwrap_or_else(|| "  ---".to_string());
        let d50_str = report.d50[band]
            .map(|v| format!("{:.3}", v))
            .unwrap_or_else(|| "  ---".to_string());

        println!(
            "{:>7.0} Hz {:>10} {:>10} {:>10} {:>10}",
            freq, rt60_str, edt_str, c80_str, d50_str
        );
    }
    println!("{:-<60}", "");

    if let Some(sti_val) = report.sti {
        println!("STI: {:.3}", sti_val);
    }
    println!();

    // Visualize with Rerun
    let draw_config = RerunConfig::new();
    let session = start_session(&draw_config)?;
    if !result.positions.is_empty() {
        draw_simulation(&session, &result, &building, &draw_config)?;
    }
    draw_receivers(&session, &[receiver])?;
    draw_impulse_response(&session, &ir, "receiver")?;

    println!("Visualization sent to Rerun (localhost:9876)");

    Ok(())
}
