use anyhow::Result;
use building3d::draw::lighting::{draw_illuminance_heatmap, draw_sensor_grid};
use building3d::draw::rerun::start_session;
use building3d::sim::lighting::config::LightingConfig;
use building3d::sim::lighting::simulation::LightingSimulation;
use building3d::sim::lighting::sources::DirectionalLight;
use building3d::sim::materials::MaterialLibrary;
use building3d::{Building, Point, Polygon, RerunConfig, Solid, Vector, Wall, Zone};

/// Build a 6x4x3m room with a window opening in the south wall.
fn build_room_with_window() -> Result<Building> {
    // Room dimensions
    let lx = 6.0;
    let ly = 4.0;
    let lz = 3.0;

    // Floor corners (z=0)
    let p0 = Point::new(0.0, 0.0, 0.0);
    let p1 = Point::new(lx, 0.0, 0.0);
    let p2 = Point::new(lx, ly, 0.0);
    let p3 = Point::new(0.0, ly, 0.0);

    // Ceiling corners (z=lz)
    let p4 = Point::new(0.0, 0.0, lz);
    let p5 = Point::new(lx, 0.0, lz);
    let p6 = Point::new(lx, ly, lz);
    let p7 = Point::new(0.0, ly, lz);

    // Floor (normal pointing down into the room interior when viewed from outside,
    // but the polygon normal auto-calculation will handle orientation)
    let floor_poly = Polygon::new("floor", vec![p0, p3, p2, p1], None)?;
    let floor_wall = Wall::new("floor", vec![floor_poly])?;

    // Ceiling
    let ceil_poly = Polygon::new("ceiling", vec![p4, p5, p6, p7], None)?;
    let ceil_wall = Wall::new("ceiling", vec![ceil_poly])?;

    // South wall (y=0) — has a 2m x 1.5m glass window pane centered at x=3, z=1.5
    // Window: x from 2.0 to 4.0, z from 0.75 to 2.25
    let south_outer = vec![p0, p1, p5, p4];
    let south_normal = Vector::new(0.0, -1.0, 0.0);
    let window_hole = vec![
        Point::new(2.0, 0.0, 0.75),
        Point::new(4.0, 0.0, 0.75),
        Point::new(4.0, 0.0, 2.25),
        Point::new(2.0, 0.0, 2.25),
    ];
    let south_poly = Polygon::with_holes(
        "south",
        south_outer,
        vec![window_hole.clone()],
        Some(south_normal),
    )?;
    // Glass pane filling the window opening
    let glass_poly = Polygon::new("glass_pane", window_hole, Some(south_normal))?;
    let south_wall = Wall::new("wall_south", vec![south_poly, glass_poly])?;

    // North wall (y=ly)
    let north_poly = Polygon::new("north", vec![p2, p3, p7, p6], None)?;
    let north_wall = Wall::new("wall_north", vec![north_poly])?;

    // East wall (x=lx)
    let east_poly = Polygon::new("east", vec![p1, p2, p6, p5], None)?;
    let east_wall = Wall::new("wall_east", vec![east_poly])?;

    // West wall (x=0)
    let west_poly = Polygon::new("west", vec![p3, p0, p4, p7], None)?;
    let west_wall = Wall::new("wall_west", vec![west_poly])?;

    let room = Solid::new(
        "room",
        vec![
            floor_wall, ceil_wall, south_wall, north_wall, east_wall, west_wall,
        ],
    )?;
    let zone = Zone::new("zone", vec![room])?;
    Building::new("window-room", vec![zone])
}

fn main() -> Result<()> {
    let building = build_room_with_window()?;

    // Material setup
    let mut lib = MaterialLibrary::with_presets();
    lib.assign("wall", "concrete");
    lib.assign("floor", "carpet");
    lib.assign("ceiling", "gypsum");
    lib.assign("glass_pane", "glass");
    lib.assign("wall_north", "metal");

    // Lighting config with directional light (sunlight from south, angled down)
    let mut config = LightingConfig::new();
    config.material_library = Some(lib);
    config.num_rays = 200_000;
    config.max_bounces = 5;

    // Sunlight entering from the south through the window
    // Direction: slightly into the room (+y) and downward (-z)
    config.directional_lights.push(DirectionalLight::new(
        Vector::new(0.0, 1.0, -0.5),
        [5000.0, 5000.0, 4500.0],
    ));

    // Sensor grid on floor with 0.2m spacing
    config.sensor_spacing = Some(0.2);
    config.sensor_patterns = vec!["floor".to_string()];

    println!("Running lighting simulation with sensor grid...");
    println!("  Directional lights: {}", config.directional_lights.len());
    println!("  Rays: {}", config.num_rays);
    println!("  Max bounces: {}", config.max_bounces);
    println!("  Sensor spacing: 0.2 m");
    println!();

    let sim = LightingSimulation::new(&building, config)?;
    let result = sim.run();

    // Print summary
    let avg = result.average_illuminance();
    let avg_lux = (avg[0] + avg[1] + avg[2]) / 3.0;
    println!("Illuminance Summary:");
    println!("{:-<50}", "");
    println!("  Average (per polygon): {:.1} lux", avg_lux);
    println!("  Surfaces hit: {}", result.hit_count.len());
    println!("  Sensor grids: {}", result.sensor_grids.len());

    for grid in &result.sensor_grids {
        let avg_s = grid.average_illuminance();
        let avg_s_lux = (avg_s[0] + avg_s[1] + avg_s[2]) / 3.0;
        println!(
            "  Grid '{}': {} sensors, avg {:.1} lux",
            grid.polygon_path,
            grid.sensors.len(),
            avg_s_lux
        );
    }
    println!("{:-<50}", "");
    println!();

    // Visualize
    let draw_config = RerunConfig::new();
    let session = start_session(&draw_config)?;

    // Draw per-polygon heatmap (walls/ceiling)
    draw_illuminance_heatmap(&session, &result, &building)?;

    // Draw per-point sensor grids on floor
    // Find max lux across all sensor grids for consistent color scale
    let max_sensor_lux = result
        .sensor_grids
        .iter()
        .flat_map(|g| g.sensors.iter())
        .map(|s| (s.illuminance[0] + s.illuminance[1] + s.illuminance[2]) / 3.0)
        .fold(0.0_f64, f64::max);

    for grid in &result.sensor_grids {
        draw_sensor_grid(&session, grid, max_sensor_lux)?;
    }

    println!("Visualization sent to Rerun (localhost:9876)");

    Ok(())
}
