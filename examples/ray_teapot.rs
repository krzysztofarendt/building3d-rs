use std::path::Path;
use std::time::Instant;

use anyhow::Result;
use building3d::RerunConfig;
use building3d::draw::rerun::{draw_simulation, start_session};
use building3d::io::stl::{read_stl, stl_to_solid};
use building3d::sim::rays::{Simulation, SimulationConfig};
use building3d::{Building, Point, Zone};

fn main() -> Result<()> {
    println!("Ray simulation in the Utah teapot model.");

    // Load teapot from STL
    let mesh = read_stl(Path::new("resources/stl/utah_teapot.stl"))?;
    println!(
        "Loaded teapot: {} vertices, {} triangles",
        mesh.vertices.len(),
        mesh.faces.as_ref().map_or(0, |f| f.len()),
    );

    let solid = stl_to_solid(&mesh, "teapot")?;
    let zone = Zone::new("z", vec![solid])?;
    let building = Building::new("b", vec![zone])?;

    // Configure simulation
    let mut config = SimulationConfig::new();
    config.num_steps = 700;
    config.num_rays = 5000;
    config.source = Point::new(0.0, 0.0, 8.0);
    config.absorbers = vec![Point::new(0.0, 0.0, 4.0)];
    config.default_absorption = 0.1;
    config.voxel_size = 0.25;
    config.search_transparent = false; // Too slow for large meshes
    config.store_ray_history = true; // needed for draw_simulation

    // Run simulation
    let t0 = Instant::now();
    let sim = Simulation::new(&building, config)?;
    let result = sim.run();
    let elapsed = t0.elapsed();
    println!("Simulation finished in {:.2}s", elapsed.as_secs_f64());

    // Print summary
    let total_hits: Vec<f64> = (0..result.config.absorbers.len())
        .map(|ai| result.hits.iter().map(|h| h[ai]).sum())
        .collect();
    println!(
        "{} steps, {} rays",
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
