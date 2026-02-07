use anyhow::Result;
use building3d::draw::rerun::start_session;
use building3d::draw::thermal::{draw_annual_timeline, draw_heat_loss_heatmap};
use building3d::sim::energy::config::ThermalConfig;
use building3d::sim::energy::construction::{concrete_wall, double_glazing, insulated_wall};
use building3d::sim::energy::hvac::HvacIdealLoads;
use building3d::sim::energy::schedule::InternalGainsProfile;
use building3d::sim::energy::simulation::run_transient_simulation;
use building3d::sim::energy::weather::WeatherData;
use building3d::sim::energy::zone::calculate_heat_balance;
use building3d::{Building, FloorPlan, RerunConfig, Solid, Zone};

/// Build an L-shaped building with two rooms in one zone.
fn build_l_shaped() -> Result<Building> {
    let room1 = Solid::from_floor_plan(FloorPlan {
        plan: vec![(0.0, 0.0), (6.0, 0.0), (6.0, 4.0), (0.0, 4.0)],
        height: 3.0,
        name: "room1".to_string(),
        wall_names: Some(vec![
            "south".to_string(),
            "east".to_string(),
            "north".to_string(),
            "west".to_string(),
        ]),
        floor_name: None,
        ceiling_name: None,
    })?;

    let room2 = Solid::from_floor_plan(FloorPlan {
        plan: vec![(0.0, 4.0), (4.0, 4.0), (4.0, 8.0), (0.0, 8.0)],
        height: 3.0,
        name: "room2".to_string(),
        wall_names: Some(vec![
            "south".to_string(),
            "east".to_string(),
            "north".to_string(),
            "west".to_string(),
        ]),
        floor_name: None,
        ceiling_name: None,
    })?;

    let zone = Zone::new("zone", vec![room1, room2])?;
    Building::new("L-building", vec![zone])
}

fn main() -> Result<()> {
    let building = build_l_shaped()?;

    // Set up thermal configuration with wall constructions
    let mut config = ThermalConfig::new();
    config.indoor_temperature = 20.0;
    config.infiltration_ach = 0.5;

    // Assign constructions to surfaces
    config
        .constructions
        .insert("wall".to_string(), insulated_wall());
    config
        .constructions
        .insert("floor".to_string(), concrete_wall());
    config
        .constructions
        .insert("ceiling".to_string(), insulated_wall());
    // South wall of room1 gets double glazing (window wall)
    config
        .constructions
        .insert("room1/south".to_string(), double_glazing());

    println!("U-values:");
    println!(
        "  Insulated wall: {:.3} W/(m2*K)",
        insulated_wall().u_value()
    );
    println!(
        "  Concrete floor: {:.3} W/(m2*K)",
        concrete_wall().u_value()
    );
    println!(
        "  Double glazing: {:.3} W/(m2*K)",
        double_glazing().u_value()
    );
    println!();

    // Create synthetic weather data (temperate climate, e.g., Central Europe)
    let weather = WeatherData::synthetic("Central Europe", 52.0, 13.0, 10.0, 15.0);
    println!(
        "Weather: {} ({} hours)",
        weather.location,
        weather.num_hours()
    );
    println!("  Mean temperature: {:.1} C", weather.mean_temperature());
    println!(
        "  Design heating temp: {:.1} C",
        weather.design_heating_temperature()
    );
    println!();

    // Create internal gains profile for office use
    // Total floor area: Room1 (6*4=24) + Room2 (4*4=16) = 40 m2
    let floor_area = 40.0;
    let gains = InternalGainsProfile::office(floor_area);

    // HVAC system
    let hvac = HvacIdealLoads::new();
    println!(
        "HVAC setpoints: heating={:.0} C, cooling={:.0} C",
        hvac.heating_setpoint, hvac.cooling_setpoint
    );
    println!();

    // Run transient annual simulation
    let solar_gain_factor = 0.3; // 30% of solar radiation enters as gains
    println!("Running transient annual simulation...");
    let annual = run_transient_simulation(
        &building,
        &config,
        &weather,
        &hvac,
        Some(&gains),
        solar_gain_factor,
    );

    // Print annual results
    println!();
    println!("Annual Energy Summary:");
    println!("{:-<50}", "");
    println!(
        "  Heating demand: {:>10.1} kWh/a",
        annual.annual_heating_kwh
    );
    println!(
        "  Cooling demand: {:>10.1} kWh/a",
        annual.annual_cooling_kwh
    );
    println!("  Peak heating:   {:>10.1} W", annual.peak_heating);
    println!("  Peak cooling:   {:>10.1} W", annual.peak_cooling);
    println!(
        "  Heating EUI:    {:>10.1} kWh/(m2*a)",
        annual.annual_heating_kwh / floor_area
    );
    println!(
        "  Cooling EUI:    {:>10.1} kWh/(m2*a)",
        annual.annual_cooling_kwh / floor_area
    );
    println!("{:-<50}", "");
    println!();

    // Monthly breakdown
    let months = [
        "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
    ];
    println!("Monthly Breakdown:");
    println!("{:>5} {:>12} {:>12}", "Month", "Heating kWh", "Cooling kWh");
    println!("{:-<35}", "");
    for (i, month) in months.iter().enumerate() {
        println!(
            "{:>5} {:>12.1} {:>12.1}",
            month, annual.monthly_heating_kwh[i], annual.monthly_cooling_kwh[i]
        );
    }
    println!("{:-<35}", "");
    println!();

    // Run single-hour heat balance for winter design conditions
    let mut winter_config = config.clone();
    winter_config.outdoor_temperature = weather.design_heating_temperature();
    winter_config.internal_gains = 0.0;
    winter_config.solar_gains = 0.0;

    let heat_balance = calculate_heat_balance(&building, &winter_config);
    println!(
        "Winter Design Heat Balance (outdoor={:.1} C):",
        winter_config.outdoor_temperature
    );
    println!("{:-<50}", "");
    println!(
        "  Transmission loss: {:>10.1} W",
        heat_balance.transmission_loss
    );
    println!(
        "  Infiltration loss: {:>10.1} W",
        heat_balance.infiltration_loss
    );
    println!(
        "  Heating demand:    {:>10.1} W",
        heat_balance.heating_demand
    );
    println!("{:-<50}", "");
    println!();

    // Print per-surface heat loss (top 10)
    let mut surface_losses: Vec<_> = heat_balance.surface_heat_loss.iter().collect();
    surface_losses.sort_by(|a, b| b.1.abs().partial_cmp(&a.1.abs()).unwrap());
    println!("Top surface heat losses (winter design):");
    for (path, loss) in surface_losses.iter().take(10) {
        println!("  {:50} {:>8.1} W", path, loss);
    }
    println!();

    // Visualize with Rerun
    let draw_config = RerunConfig::new();
    let session = start_session(&draw_config)?;
    draw_heat_loss_heatmap(&session, &heat_balance, &building)?;
    draw_annual_timeline(&session, &annual)?;

    println!("Visualization sent to Rerun (localhost:9876)");

    Ok(())
}
