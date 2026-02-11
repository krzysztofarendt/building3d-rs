use building3d::sim::energy::config::ThermalConfig;
use building3d::sim::energy::construction::WallConstruction;
use building3d::sim::energy::hvac::HvacIdealLoads;
use building3d::sim::energy::simulation::{
    TransientSimulationOptions, run_transient_simulation, run_transient_simulation_with_options,
};
use building3d::sim::energy::solar_bridge::SolarGainConfig;
use building3d::sim::energy::weather::WeatherData;
use building3d::sim::materials::Layer;
use building3d::{Building, Point, Polygon, Solid, Wall, Zone};
use std::path::PathBuf;

fn bestest_600_constructions() -> (
    WallConstruction,
    WallConstruction,
    WallConstruction,
    WallConstruction,
) {
    let lt_wall = WallConstruction::new(
        "LTWALL",
        vec![
            Layer {
                name: "WOOD SIDING-1".to_string(),
                thickness: 0.009,
                conductivity: 0.14,
                density: 530.0,
                specific_heat: 900.0,
            },
            Layer {
                name: "FIBERGLASS QUILT-1".to_string(),
                // Target RSI ≈ 1.94 m²K/W (R-11 in IP units).
                thickness: 1.94 * 0.04,
                conductivity: 0.04,
                density: 12.0,
                specific_heat: 840.0,
            },
            Layer {
                name: "PLASTERBOARD-1".to_string(),
                thickness: 0.012,
                conductivity: 0.16,
                density: 950.0,
                specific_heat: 840.0,
            },
        ],
    );

    let lt_roof = WallConstruction::roof(
        "LTROOF",
        vec![
            Layer {
                name: "ROOF DECK".to_string(),
                thickness: 0.019,
                conductivity: 0.14,
                density: 530.0,
                specific_heat: 900.0,
            },
            Layer {
                name: "FIBERGLASS QUILT-2".to_string(),
                // Target RSI ≈ 3.35 m²K/W (R-19 in IP units).
                thickness: 3.35 * 0.04,
                conductivity: 0.04,
                density: 12.0,
                specific_heat: 840.0,
            },
            Layer {
                name: "PLASTERBOARD-2".to_string(),
                thickness: 0.01,
                conductivity: 0.16,
                density: 950.0,
                specific_heat: 840.0,
            },
        ],
    );

    // Represent Material:NoMass "R-25 INSULATION" as a pure resistance layer.
    // In the BESTEST IDF this value is already in m²*K/W.
    let r25 = 25.075_f64;
    let lt_floor = WallConstruction::floor(
        "LTFLOOR",
        vec![
            Layer {
                name: "R-25 INSULATION (no mass)".to_string(),
                thickness: 1.0,
                conductivity: 1.0 / r25,
                density: 0.0,
                specific_heat: 0.0,
            },
            Layer {
                name: "TIMBER FLOORING".to_string(),
                thickness: 0.025,
                conductivity: 0.14,
                density: 650.0,
                specific_heat: 1200.0,
            },
        ],
    );

    let window = WallConstruction::new(
        "Double Pane Window (approx)",
        vec![
            Layer {
                name: "Glass".to_string(),
                thickness: 0.003175,
                conductivity: 1.06,
                density: 2500.0,
                specific_heat: 750.0,
            },
            Layer {
                name: "Air gap (approx)".to_string(),
                thickness: 0.013,
                conductivity: 0.026,
                density: 0.0,
                specific_heat: 0.0,
            },
            Layer {
                name: "Glass".to_string(),
                thickness: 0.003175,
                conductivity: 1.06,
                density: 2500.0,
                specific_heat: 750.0,
            },
        ],
    );

    (lt_wall, lt_roof, lt_floor, window)
}

fn bestest_900_constructions() -> (
    WallConstruction,
    WallConstruction,
    WallConstruction,
    WallConstruction,
) {
    // From BESTEST-GSR shared `bestest_resources.osm` (construction set "BESTEST HW").
    //
    // NOTE: In that resource model, the heavyweight case uses:
    // - HWWALL (WOOD SIDING-1 + FOAM INSULATION + CONCRETE BLOCK)
    // - HWFLOOR (R-25 INSULATION + CONCRETE SLAB)
    // - LTROOF (same as the light-mass case)
    // - Double Pane Window (same as the light-mass case)
    let (_lt_wall, lt_roof, _lt_floor, window) = bestest_600_constructions();

    let hw_wall = WallConstruction::new(
        "HWWALL",
        vec![
            Layer {
                name: "WOOD SIDING-1".to_string(),
                thickness: 0.009,
                conductivity: 0.14,
                density: 530.0,
                specific_heat: 900.0,
            },
            Layer {
                name: "FOAM INSULATION".to_string(),
                thickness: 0.0615,
                conductivity: 0.04,
                density: 10.0,
                specific_heat: 1400.0,
            },
            Layer {
                name: "CONCRETE BLOCK".to_string(),
                thickness: 0.1,
                conductivity: 0.51,
                density: 1400.0,
                specific_heat: 1000.0,
            },
        ],
    );

    // HWFLOOR: R-25 INSULATION (no mass) + CONCRETE SLAB
    // In BESTEST-GSR resources, R-25 is 25.175 m²*K/W.
    let r25 = 25.175_f64;
    let hw_floor = WallConstruction::floor(
        "HWFLOOR",
        vec![
            Layer {
                name: "R-25 INSULATION (no mass)".to_string(),
                thickness: 1.0,
                conductivity: 1.0 / r25,
                density: 0.0,
                specific_heat: 0.0,
            },
            Layer {
                name: "CONCRETE SLAB".to_string(),
                thickness: 0.08,
                conductivity: 1.13,
                density: 1400.0,
                specific_heat: 1000.0,
            },
        ],
    );

    (hw_wall, lt_roof, hw_floor, window)
}

fn poly_rect_y0(name: &str, x0: f64, x1: f64, z0: f64, z1: f64) -> Polygon {
    let y = 0.0;
    Polygon::new(
        name,
        vec![
            Point::new(x0, y, z0),
            Point::new(x1, y, z0),
            Point::new(x1, y, z1),
            Point::new(x0, y, z1),
        ],
        Some(building3d::Vector::new(0.0, -1.0, 0.0)),
    )
    .unwrap()
}

fn build_bestest_600_geometry() -> Building {
    let x = 8.0;
    let y = 6.0;
    let z = 2.7;

    let p0 = Point::new(0.0, 0.0, 0.0);
    let p1 = Point::new(x, 0.0, 0.0);
    let p2 = Point::new(x, y, 0.0);
    let p3 = Point::new(0.0, y, 0.0);
    let p4 = Point::new(0.0, 0.0, z);
    let p5 = Point::new(x, 0.0, z);
    let p6 = Point::new(x, y, z);
    let p7 = Point::new(0.0, y, z);

    let poly_floor = Polygon::new("floor", vec![p0, p3, p2, p1], None).unwrap();
    let poly_roof = Polygon::new("ceiling", vec![p4, p5, p6, p7], None).unwrap();

    let mut south_polys: Vec<Polygon> = Vec::new();
    south_polys.push(poly_rect_y0("opaque_left", 0.0, 0.5, 0.0, z));
    south_polys.push(poly_rect_y0("opaque_mid", 3.5, 4.5, 0.0, z));
    south_polys.push(poly_rect_y0("opaque_right", 7.5, 8.0, 0.0, z));

    south_polys.push(poly_rect_y0("opaque_below_w1", 0.5, 3.5, 0.0, 0.2));
    south_polys.push(poly_rect_y0("opaque_above_w1", 0.5, 3.5, 2.2, z));
    south_polys.push(poly_rect_y0("window_1", 0.5, 3.5, 0.2, 2.2));

    south_polys.push(poly_rect_y0("opaque_below_w2", 4.5, 7.5, 0.0, 0.2));
    south_polys.push(poly_rect_y0("opaque_above_w2", 4.5, 7.5, 2.2, z));
    south_polys.push(poly_rect_y0("window_2", 4.5, 7.5, 0.2, 2.2));

    let wall_floor = Wall::new("floor", vec![poly_floor]).unwrap();
    let wall_0 = Wall::new("wall_0", south_polys).unwrap();
    let wall_1 = Wall::new(
        "wall_1",
        vec![Polygon::new("poly_1", vec![p1, p2, p6, p5], None).unwrap()],
    )
    .unwrap();
    let wall_2 = Wall::new(
        "wall_2",
        vec![Polygon::new("poly_2", vec![p3, p7, p6, p2], None).unwrap()],
    )
    .unwrap();
    let wall_3 = Wall::new(
        "wall_3",
        vec![Polygon::new("poly_3", vec![p0, p4, p7, p3], None).unwrap()],
    )
    .unwrap();
    let wall_roof = Wall::new("ceiling", vec![poly_roof]).unwrap();

    let solid = Solid::new(
        "space",
        vec![wall_floor, wall_0, wall_1, wall_2, wall_3, wall_roof],
    )
    .unwrap();
    let zone = Zone::new("zone_one", vec![solid]).unwrap();
    Building::new("bestest_case_600", vec![zone]).unwrap()
}

fn estimate_zone_capacity_j_per_m3_k(building: &Building, cfg: &ThermalConfig) -> f64 {
    let volume_m3: f64 = building.zones().iter().map(|z| z.volume()).sum();
    if volume_m3 <= 0.0 {
        return cfg.thermal_capacity_j_per_m3_k;
    }

    let mut c_total_j_per_k = 0.0;
    for zone in building.zones() {
        for solid in zone.solids() {
            for wall in solid.walls() {
                for polygon in wall.polygons() {
                    let path = format!(
                        "{}/{}/{}/{}",
                        zone.name, solid.name, wall.name, polygon.name
                    );
                    let c_j_per_m2_k =
                        cfg.resolve_envelope_capacity_j_per_m2_k(Some(&polygon.uid), &path, 0.0);
                    c_total_j_per_k += c_j_per_m2_k * polygon.area();
                }
            }
        }
    }

    (c_total_j_per_k / volume_m3).max(0.0)
}

fn make_cfg_600(building: &Building) -> ThermalConfig {
    let (lt_wall, lt_roof, lt_floor, window) = bestest_600_constructions();
    let mut cfg = ThermalConfig::new();
    cfg.default_u_value = 0.0;
    cfg.infiltration_ach = 0.5;
    cfg.internal_gains = 200.0;
    cfg.indoor_temperature = 20.0;

    cfg.use_fvm_walls = true;

    cfg.constructions.insert("window".to_string(), window);
    cfg.constructions.insert("ceiling".to_string(), lt_roof);
    cfg.constructions.insert("floor".to_string(), lt_floor);
    cfg.constructions.insert("wall".to_string(), lt_wall);

    // Treat windows as a whole-window U-value (the layered "air gap as pure conduction"
    // approximation is far too insulating).
    cfg.u_value_overrides_by_path_pattern
        .insert("window".to_string(), 1.8);

    // BESTEST case floors are ground-coupled in the reference model.
    cfg.ground_temperature_c = Some(10.0);

    // Surface-aware policy for transmitted solar + radiant internal gains.
    cfg.use_surface_aware_solar_distribution = true;
    cfg.transmitted_solar_to_air_fraction = 0.0;
    cfg.internal_gains_to_mass_fraction = 0.6; // from BESTEST-GSR "OtherEquipment" radiant fraction

    cfg.thermal_capacity_j_per_m3_k = estimate_zone_capacity_j_per_m3_k(building, &cfg);
    cfg
}

fn make_cfg_900(building: &Building) -> ThermalConfig {
    let (wall, roof, floor, window) = bestest_900_constructions();
    let mut cfg = make_cfg_600(building);

    cfg.constructions.insert("window".to_string(), window);
    cfg.constructions.insert("ceiling".to_string(), roof);
    cfg.constructions.insert("floor".to_string(), floor);
    cfg.constructions.insert("wall".to_string(), wall);

    cfg.thermal_capacity_j_per_m3_k = estimate_zone_capacity_j_per_m3_k(building, &cfg);
    cfg
}

fn solar_cfg() -> SolarGainConfig {
    let mut solar = SolarGainConfig::new();
    solar.glazing_patterns = vec!["window".to_string()];
    solar.default_shgc = 0.86156;
    solar.include_ground_reflection = true;
    solar.ground_reflectance = 0.2;
    solar.include_incidence_angle_modifier = true;
    solar.incidence_angle_modifier_a = 0.1;
    solar.include_exterior_opaque_absorption = true;
    solar.default_opaque_absorptance = 0.6;
    // Phase 2.2 refinements (optional): exterior longwave exchange + wind-based h_out.
    solar.include_exterior_longwave_exchange = false;
    solar.use_wind_speed_for_h_out = false;
    solar
}

fn assert_rel_close(name: &str, got: f64, expected: f64, rel_tol: f64) {
    let err = (got - expected).abs();
    let scale = expected.abs().max(1.0);
    assert!(
        err <= rel_tol * scale,
        "{name}: got {got}, expected {expected}, rel err {:.3}%",
        100.0 * err / scale
    );
}

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn find_bestest_epw() -> Option<PathBuf> {
    let candidates = [
        repo_root()
            .join("examples/bestest_energy_suite/data/USA_MA_Boston-Logan.Intl.AP.725090_TMY3.epw"),
        repo_root()
            .join("examples/bestest_600_energy/data/USA_MA_Boston-Logan.Intl.AP.725090_TMY3.epw"),
    ];
    candidates.into_iter().find(|p| p.exists())
}

#[test]
fn test_synthetic_weather_outputs_are_finite() {
    let building = build_bestest_600_geometry();
    let cfg = make_cfg_600(&building);
    let hvac = HvacIdealLoads::with_setpoints(20.0, 27.0);
    let solar = solar_cfg();
    let weather = WeatherData::synthetic("bestest_synth", 42.0, 0.0, 10.0, 15.0);

    let annual = run_transient_simulation(&building, &cfg, &weather, &hvac, None, Some(&solar));

    assert!(annual.annual_heating_kwh.is_finite());
    assert!(annual.annual_cooling_kwh.is_finite());
    assert!(annual.annual_heating_kwh >= 0.0);
    assert!(annual.annual_cooling_kwh >= 0.0);
}

#[test]
fn test_no_solar_increases_heating_decreases_cooling() {
    let building = build_bestest_600_geometry();
    let cfg = make_cfg_600(&building);
    let hvac = HvacIdealLoads::with_setpoints(20.0, 27.0);
    let solar = solar_cfg();
    let weather = WeatherData::synthetic("bestest_synth", 42.0, 0.0, 10.0, 15.0);

    let with_solar = run_transient_simulation(&building, &cfg, &weather, &hvac, None, Some(&solar));
    let no_solar = run_transient_simulation(&building, &cfg, &weather, &hvac, None, None);

    assert!(
        no_solar.annual_heating_kwh > with_solar.annual_heating_kwh,
        "Expected heating to increase without solar gains"
    );
    assert!(
        no_solar.annual_cooling_kwh < with_solar.annual_cooling_kwh,
        "Expected cooling to decrease without solar gains"
    );
}

#[test]
fn test_bestest_600_epw_reference_within_tolerance_if_present() {
    let Some(epw_path) = find_bestest_epw() else {
        eprintln!("Skipping EPW validation (no EPW found). Run download_data.sh in examples.");
        return;
    };

    // Reference annual totals from BESTEST-GSR (OpenStudio/EnergyPlus), in kWh.
    let ref_heating_kwh = 4324.76;
    let ref_cooling_kwh = 6044.07;

    let building = build_bestest_600_geometry();
    let cfg = make_cfg_600(&building);
    let hvac = HvacIdealLoads::with_setpoints(20.0, 27.0);
    let solar = solar_cfg();

    let epw_content = std::fs::read_to_string(&epw_path).unwrap();
    let weather = WeatherData::from_epw(&epw_content).unwrap();

    let options = TransientSimulationOptions {
        warmup_hours: 7 * 24,
    };
    let annual = run_transient_simulation_with_options(
        &building,
        &cfg,
        &weather,
        &hvac,
        None,
        Some(&solar),
        &options,
    );

    // Wide tolerances: this is a simplified model (1R1C + simple solar gains),
    // but we still want to catch large regressions.
    assert_rel_close(
        "epw_annual_heating_kwh",
        annual.annual_heating_kwh,
        ref_heating_kwh,
        0.20,
    );
    assert_rel_close(
        "epw_annual_cooling_kwh",
        annual.annual_cooling_kwh,
        ref_cooling_kwh,
        0.20,
    );
}

#[test]
fn test_bestest_900_epw_reference_within_tolerance_if_present() {
    let Some(epw_path) = find_bestest_epw() else {
        eprintln!("Skipping EPW validation (no EPW found). Run download_data.sh in examples.");
        return;
    };

    // Reference annual totals from BESTEST-GSR (OpenStudio/EnergyPlus), in kWh.
    let ref_heating_kwh = 1661.17;
    let ref_cooling_kwh = 2498.16;

    let building = build_bestest_600_geometry();
    let cfg = make_cfg_900(&building);

    let hvac = HvacIdealLoads::with_setpoints(20.0, 27.0);
    let solar = solar_cfg();

    let epw_content = std::fs::read_to_string(&epw_path).unwrap();
    let weather = WeatherData::from_epw(&epw_content).unwrap();

    let options = TransientSimulationOptions {
        warmup_hours: 7 * 24,
    };
    let annual = run_transient_simulation_with_options(
        &building,
        &cfg,
        &weather,
        &hvac,
        None,
        Some(&solar),
        &options,
    );

    assert_rel_close(
        "epw_900_annual_heating_kwh",
        annual.annual_heating_kwh,
        ref_heating_kwh,
        0.30,
    );
    assert_rel_close(
        "epw_900_annual_cooling_kwh",
        annual.annual_cooling_kwh,
        ref_cooling_kwh,
        0.40,
    );
}
