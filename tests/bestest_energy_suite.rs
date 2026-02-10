use building3d::sim::energy::config::ThermalConfig;
use building3d::sim::energy::construction::WallConstruction;
use building3d::sim::energy::hvac::HvacIdealLoads;
use building3d::sim::energy::simulation::run_transient_simulation;
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

fn make_cfg(_building: &Building) -> ThermalConfig {
    let (lt_wall, lt_roof, lt_floor, window) = bestest_600_constructions();
    let mut cfg = ThermalConfig::new();
    cfg.default_u_value = 0.0;
    cfg.infiltration_ach = 0.5;
    cfg.internal_gains = 200.0;
    cfg.indoor_temperature = 20.0;

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

    // Keep a fixed deterministic thermal capacity for regression (J/(m3*K)).
    cfg.thermal_capacity_j_per_m3_k = 22_000.0;
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
    let cfg = make_cfg(&building);
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
    let cfg = make_cfg(&building);
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
    let cfg = make_cfg(&building);
    let hvac = HvacIdealLoads::with_setpoints(20.0, 27.0);
    let solar = solar_cfg();

    let epw_content = std::fs::read_to_string(&epw_path).unwrap();
    let weather = WeatherData::from_epw(&epw_content).unwrap();

    let annual = run_transient_simulation(&building, &cfg, &weather, &hvac, None, Some(&solar));

    // Wide tolerances: this is a simplified model (1R1C + simple solar gains),
    // but we still want to catch large regressions.
    assert_rel_close(
        "epw_annual_heating_kwh",
        annual.annual_heating_kwh,
        ref_heating_kwh,
        0.10,
    );
    assert_rel_close(
        "epw_annual_cooling_kwh",
        annual.annual_cooling_kwh,
        ref_cooling_kwh,
        0.15,
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
    let mut cfg = make_cfg(&building);
    cfg.thermal_capacity_j_per_m3_k *= 8.0;
    cfg.two_node_mass_fraction = 0.95;
    cfg.interior_heat_transfer_coeff_w_per_m2_k = 8.0;
    cfg.solar_gains_to_mass_fraction = 0.9;
    cfg.internal_gains_to_mass_fraction = 0.0;
    cfg.two_node_envelope_to_mass = true;

    let hvac = HvacIdealLoads::with_setpoints(20.0, 27.0);
    let solar = solar_cfg();

    let epw_content = std::fs::read_to_string(&epw_path).unwrap();
    let weather = WeatherData::from_epw(&epw_content).unwrap();

    let annual = run_transient_simulation(&building, &cfg, &weather, &hvac, None, Some(&solar));

    assert_rel_close(
        "epw_900_annual_heating_kwh",
        annual.annual_heating_kwh,
        ref_heating_kwh,
        0.25,
    );
    assert_rel_close(
        "epw_900_annual_cooling_kwh",
        annual.annual_cooling_kwh,
        ref_cooling_kwh,
        0.25,
    );
}
