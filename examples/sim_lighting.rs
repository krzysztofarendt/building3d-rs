use anyhow::Result;
use building3d::draw::lighting::draw_illuminance_heatmap;
use building3d::draw::rerun::start_session;
use building3d::sim::lighting::config::LightingConfig;
use building3d::sim::lighting::simulation::LightingSimulation;
use building3d::sim::lighting::sources::PointLight;
use building3d::sim::materials::MaterialLibrary;
use building3d::{Building, FloorPlan, Point, RerunConfig, Solid, Zone};

/// Build an L-shaped room (single solid) in one zone.
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

    // Set up material library with optical properties
    let mut lib = MaterialLibrary::with_presets();

    // Assign materials: concrete walls, gypsum ceiling, carpet floor
    lib.assign("wall", "concrete");
    lib.assign("floor", "carpet");
    lib.assign("ceiling", "gypsum");

    // Configure lighting simulation
    let mut config = LightingConfig::new();
    config.material_library = Some(lib);
    config.num_rays = 50000;
    config.max_bounces = 5;

    // Point light in lower part of the L-shape (ceiling-mounted)
    config
        .point_lights
        .push(PointLight::white(Point::new(3.0, 2.0, 2.8), 3000.0));
    // Point light in upper part of the L-shape (ceiling-mounted)
    config
        .point_lights
        .push(PointLight::white(Point::new(2.0, 6.0, 2.8), 2000.0));

    println!("Running forward lighting simulation...");
    println!("  Lights: {}", config.point_lights.len());
    println!("  Rays per light: {}", config.num_rays);
    println!("  Max bounces: {}", config.max_bounces);
    println!();

    // Run simulation
    let sim = LightingSimulation::new(&building, config)?;
    let result = sim.run();

    // Print illuminance summary
    let avg = result.average_illuminance();
    let avg_lux = (avg[0] + avg[1] + avg[2]) / 3.0;
    println!("Illuminance Summary:");
    println!("{:-<50}", "");
    println!(
        "  Average: {:.1} lux (R={:.1}, G={:.1}, B={:.1})",
        avg_lux, avg[0], avg[1], avg[2]
    );

    let mut min_lux = f64::MAX;
    let mut max_lux = 0.0_f64;
    let mut min_path = String::new();
    let mut max_path = String::new();

    for (path, ill) in &result.illuminance {
        let lux = (ill[0] + ill[1] + ill[2]) / 3.0;
        if lux < min_lux {
            min_lux = lux;
            min_path = path.clone();
        }
        if lux > max_lux {
            max_lux = lux;
            max_path = path.clone();
        }
    }

    println!("  Min: {:.1} lux ({})", min_lux, min_path);
    println!("  Max: {:.1} lux ({})", max_lux, max_path);
    println!("  Surfaces hit: {}", result.hit_count.len());
    println!("{:-<50}", "");
    println!();

    // Print per-surface illuminance
    println!("Per-surface illuminance:");
    let mut surfaces: Vec<_> = result.illuminance.iter().collect();
    surfaces.sort_by(|a, b| a.0.cmp(b.0));
    for (path, ill) in &surfaces {
        let lux = (ill[0] + ill[1] + ill[2]) / 3.0;
        let hits = result.hit_count.get(*path).copied().unwrap_or(0);
        println!("  {:50} {:>8.1} lux  ({} hits)", path, lux, hits);
    }
    println!();

    // Visualize with Rerun
    let draw_config = RerunConfig::new();
    let session = start_session(&draw_config)?;
    draw_illuminance_heatmap(&session, &result, &building)?;

    println!("Visualization sent to Rerun (localhost:9876)");

    Ok(())
}
