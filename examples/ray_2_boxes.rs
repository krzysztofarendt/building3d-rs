use anyhow::Result;
use building3d::RerunConfig;
use building3d::draw::rerun::{draw_simulation, start_session};
use building3d::sim::rays::{Simulation, SimulationConfig};
use building3d::{Building, Point, Solid, Zone};

fn main() -> Result<()> {
    // Create building: two adjacent boxes in the same zone
    let s0 = Solid::from_box(1.0, 1.0, 1.0, None, "s0")?;
    let s1 = Solid::from_box(1.0, 1.0, 1.0, Some((1.0, 0.0, 0.0)), "s1")?;
    let zone = Zone::new("z", vec![s0, s1])?;
    let building = Building::new("b", vec![zone])?;

    // Configure simulation
    let mut config = SimulationConfig::new();
    config.num_steps = 200;
    config.num_rays = 100;
    config.source = Point::new(0.3, 0.3, 0.3);
    config.absorbers = vec![Point::new(0.6, 0.6, 0.6), Point::new(0.1, 0.1, 0.6)];
    config.store_ray_history = true; // needed for draw_simulation

    // Run simulation
    let sim = Simulation::new(&building, config)?;
    let result = sim.run();

    // Print summary
    let total_hits: Vec<f64> = (0..result.config.absorbers.len())
        .map(|ai| result.hits.iter().map(|h| h[ai]).sum())
        .collect();
    println!(
        "Simulation complete: {} steps, {} rays",
        result.config.num_steps,
        result.positions[0].len()
    );
    for (i, hit) in total_hits.iter().enumerate() {
        println!("Absorber {}: total energy absorbed = {:.4}", i, hit);
    }

    // Visualize with Rerun
    let draw_config = RerunConfig::new();
    let session = start_session(&draw_config)?;
    draw_simulation(&session, &result, &building, &draw_config)?;

    Ok(())
}
