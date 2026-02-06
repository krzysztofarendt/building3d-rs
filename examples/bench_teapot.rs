use std::path::Path;
use std::time::Instant;

use anyhow::Result;
use building3d::io::stl::{read_stl, stl_to_solid};
use building3d::sim::rays::{Simulation, SimulationConfig};
use building3d::{Building, Point, Zone};

fn main() -> Result<()> {
    println!("Benchmark: Ray simulation in the Utah teapot model.");

    // Load teapot from STL
    let t_load = Instant::now();
    let mesh = read_stl(Path::new("resources/stl/utah_teapot.stl"))?;
    println!(
        "Loaded teapot: {} vertices, {} triangles ({:.2}s)",
        mesh.vertices.len(),
        mesh.faces.as_ref().map_or(0, |f| f.len()),
        t_load.elapsed().as_secs_f64(),
    );

    // Build model
    let t_setup = Instant::now();
    let solid = stl_to_solid(&mesh, "teapot")?;
    let zone = Zone::new("z", vec![solid]);
    let building = Building::new("b", vec![zone]);

    // Configure simulation
    let mut config = SimulationConfig::new();
    config.num_steps = 700;
    config.num_rays = 5000;
    config.source = Point::new(0.0, 0.0, 8.0);
    config.absorbers = vec![Point::new(0.0, 0.0, 4.0)];
    config.default_absorption = 0.1;
    config.voxel_size = 0.25;
    config.search_transparent = false;

    let sim = Simulation::new(&building, config)?;
    let setup_time = t_setup.elapsed().as_secs_f64();
    println!("Setup time: {:.2}s", setup_time);

    // Run simulation
    println!("Running simulation: 700 steps, 5000 rays...");
    let t_sim = Instant::now();
    let result = sim.run();
    let sim_time = t_sim.elapsed().as_secs_f64();
    println!("Simulation time: {:.2}s", sim_time);

    // Summary
    let total_hits: Vec<f64> = (0..result.config.absorbers.len())
        .map(|ai| result.hits.iter().map(|h| h[ai]).sum())
        .collect();
    for (i, hit) in total_hits.iter().enumerate() {
        println!("Absorber {}: total energy absorbed = {:.4}", i, hit);
    }

    println!("\nTotal time (setup + sim): {:.2}s", setup_time + sim_time);
    println!("  Setup: {:.2}s", setup_time);
    println!("  Simulation: {:.2}s", sim_time);

    Ok(())
}
