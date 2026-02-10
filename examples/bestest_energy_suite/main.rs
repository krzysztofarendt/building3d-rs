use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};

use building3d::sim::energy::config::ThermalConfig;
use building3d::sim::energy::construction::WallConstruction;
use building3d::sim::energy::hvac::HvacIdealLoads;
use building3d::sim::energy::simulation::{
    AnnualResult, TransientSimulationOptions, run_transient_simulation_with_options,
};
use building3d::sim::energy::solar_bridge::{
    SolarGainConfig, SolarHourParams, compute_solar_gains_with_materials,
};
use building3d::sim::energy::weather::WeatherData;
use building3d::sim::materials::Layer;
use building3d::{Building, Point, Polygon, Solid, Wall, Zone};

// Reference outputs from NREL/BESTEST-GSR (`results/workflow_results.csv`)
// (OpenStudio/EnergyPlus) for the Boston Logan TMY3 EPW.
//
// Values below are derived numeric results (monthly ideal loads, kWh).
const REF_600_MONTHLY_HEATING_KWH: [f64; 12] = [
    712.08, 682.83, 472.18, 509.87, 136.75, 10.08, 12.02, 6.58, 73.33, 347.66, 625.41, 735.97,
];
const REF_600_MONTHLY_COOLING_KWH: [f64; 12] = [
    523.18, 388.87, 485.77, 241.13, 340.91, 542.19, 537.21, 650.81, 735.08, 676.95, 438.99, 482.98,
];

const REF_900_MONTHLY_HEATING_KWH: [f64; 12] = [
    255.36, 293.66, 125.29, 268.29, 33.28, 0.0, 0.0, 0.0, 0.11, 80.96, 296.63, 307.59,
];
const REF_900_MONTHLY_COOLING_KWH: [f64; 12] = [
    53.41, 12.48, 65.79, 18.63, 114.91, 404.25, 399.74, 512.94, 493.44, 293.16, 68.29, 61.12,
];

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn default_epw_path() -> PathBuf {
    repo_root()
        .join("examples/bestest_energy_suite/data")
        .join("USA_MA_Boston-Logan.Intl.AP.725090_TMY3.epw")
}

fn sum(monthly: &[f64; 12]) -> f64 {
    monthly.iter().sum()
}

fn day_of_year(month: u8, day: u8) -> u16 {
    const DAYS_BEFORE_MONTH: [u16; 12] = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];
    let m = (month as usize).saturating_sub(1).min(11);
    DAYS_BEFORE_MONTH[m] + day as u16
}

fn bestest_600_constructions() -> (
    WallConstruction,
    WallConstruction,
    WallConstruction,
    WallConstruction,
) {
    // From `case_en_600.idf` (BESTEST-GSR workflow resource).

    // LTWALL: WOOD SIDING-1 + FIBERGLASS QUILT-1 + PLASTERBOARD-1
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

    // LTROOF: ROOF DECK + FIBERGLASS QUILT-2 + PLASTERBOARD-2
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

    // LTFLOOR: R-25 INSULATION (no mass) + TIMBER FLOORING
    //
    // Represent Material:NoMass "R-25 INSULATION" as a pure resistance layer:
    //   R = thickness / conductivity  => choose thickness=1, conductivity=1/R.
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

    // Double Pane Window: approximate layers.
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

fn poly_rect_y0(name: &str, x0: f64, x1: f64, z0: f64, z1: f64) -> Result<Polygon> {
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
}

fn build_case_600_geometry() -> Result<Building> {
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

    let poly_floor = Polygon::new("floor", vec![p0, p3, p2, p1], None)?;
    let poly_roof = Polygon::new("ceiling", vec![p4, p5, p6, p7], None)?;

    let mut south_polys: Vec<Polygon> = Vec::new();
    south_polys.push(poly_rect_y0("opaque_left", 0.0, 0.5, 0.0, z)?);
    south_polys.push(poly_rect_y0("opaque_mid", 3.5, 4.5, 0.0, z)?);
    south_polys.push(poly_rect_y0("opaque_right", 7.5, 8.0, 0.0, z)?);

    south_polys.push(poly_rect_y0("opaque_below_w1", 0.5, 3.5, 0.0, 0.2)?);
    south_polys.push(poly_rect_y0("opaque_above_w1", 0.5, 3.5, 2.2, z)?);
    south_polys.push(poly_rect_y0("window_1", 0.5, 3.5, 0.2, 2.2)?);

    south_polys.push(poly_rect_y0("opaque_below_w2", 4.5, 7.5, 0.0, 0.2)?);
    south_polys.push(poly_rect_y0("opaque_above_w2", 4.5, 7.5, 2.2, z)?);
    south_polys.push(poly_rect_y0("window_2", 4.5, 7.5, 0.2, 2.2)?);

    let wall_floor = Wall::new("floor", vec![poly_floor])?;
    let wall_0 = Wall::new("wall_0", south_polys)?;
    let wall_1 = Wall::new(
        "wall_1",
        vec![Polygon::new("poly_1", vec![p1, p2, p6, p5], None)?],
    )?;
    let wall_2 = Wall::new(
        "wall_2",
        vec![Polygon::new("poly_2", vec![p3, p7, p6, p2], None)?],
    )?;
    let wall_3 = Wall::new(
        "wall_3",
        vec![Polygon::new("poly_3", vec![p0, p4, p7, p3], None)?],
    )?;
    let wall_roof = Wall::new("ceiling", vec![poly_roof])?;

    let solid = Solid::new(
        "space",
        vec![wall_floor, wall_0, wall_1, wall_2, wall_3, wall_roof],
    )?;
    let zone = Zone::new("zone_one", vec![solid])?;
    Building::new("bestest_case_600", vec![zone])
}

fn config_for_case_600(building: &Building) -> ThermalConfig {
    let (lt_wall, lt_roof, lt_floor, window) = bestest_600_constructions();

    let mut cfg = ThermalConfig::new();
    cfg.default_u_value = 0.0; // fail closed
    cfg.infiltration_ach = 0.5;
    cfg.internal_gains = 200.0;
    cfg.indoor_temperature = 20.0;

    cfg.constructions.insert("window".to_string(), window);
    cfg.constructions.insert("ceiling".to_string(), lt_roof);
    cfg.constructions.insert("floor".to_string(), lt_floor);
    cfg.constructions.insert("wall".to_string(), lt_wall);

    // BESTEST window conduction should be treated as a whole-window U-value.
    // Using a layered construction with a still-air gap modeled as pure conduction
    // is far too insulating (missing convection+radiation in the gap).
    cfg.u_value_overrides_by_path_pattern
        .insert("window".to_string(), 1.8);

    // BESTEST case floors are ground-coupled in the reference model.
    cfg.ground_temperature_c = Some(10.0);

    // Surface-aware policy for transmitted solar + radiant internal gains.
    //
    // With FVM envelope walls, routing all transmitted shortwave into the zone air node
    // over-predicts peaks and fails to use the envelope thermal mass. Instead, deposit
    // most of it on interior surfaces (floor-first), and similarly deposit the radiant
    // fraction of internal gains to surfaces.
    cfg.use_surface_aware_solar_distribution = true;
    cfg.transmitted_solar_to_air_fraction = 0.0;
    cfg.internal_gains_to_mass_fraction = 0.6; // from BESTEST-GSR "OtherEquipment" radiant fraction

    let _ = building;
    cfg
}

fn solar_config_for_case_600() -> SolarGainConfig {
    let mut solar = SolarGainConfig::new();
    solar.glazing_patterns = vec!["window".to_string()];
    // From `Glass Type 1` solar transmittance in the IDF (approximate SHGC).
    solar.default_shgc = 0.86156;
    // Small physics refinements (still simplified vs EnergyPlus):
    // - ground-reflected shortwave using a constant albedo
    // - a simple incidence-angle modifier to reduce gains at grazing angles
    solar.include_ground_reflection = true;
    solar.ground_reflectance = 0.2;
    solar.include_incidence_angle_modifier = true;
    solar.incidence_angle_modifier_a = 0.1;
    solar.include_exterior_opaque_absorption = true;
    solar.default_opaque_absorptance = 0.6;
    // Phase 2.2 refinements: exterior longwave exchange + wind-based h_out.
    solar.include_exterior_longwave_exchange = false;
    solar.use_wind_speed_for_h_out = false;
    solar
}

#[derive(Debug, Clone, Copy)]
enum RefKind {
    Ref600,
    Ref900,
    None,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CaseMassKind {
    Light,
    Heavy,
}

#[derive(Debug, Clone)]
struct CaseSpec {
    name: &'static str,
    ref_kind: RefKind,
    mass_kind: CaseMassKind,
    solar: bool,
}

fn reference_monthly(metric: &str, kind: RefKind) -> Option<[f64; 12]> {
    match (kind, metric) {
        (RefKind::Ref600, "heating") => Some(REF_600_MONTHLY_HEATING_KWH),
        (RefKind::Ref600, "cooling") => Some(REF_600_MONTHLY_COOLING_KWH),
        (RefKind::Ref900, "heating") => Some(REF_900_MONTHLY_HEATING_KWH),
        (RefKind::Ref900, "cooling") => Some(REF_900_MONTHLY_COOLING_KWH),
        _ => None,
    }
}

fn write_suite_csv(path: &Path, cases: &[(CaseSpec, AnnualResult)]) -> Result<()> {
    let mut f = fs::File::create(path).context("create results.csv")?;
    writeln!(
        f,
        "case,month,metric,simulated_kwh,reference_kwh,error_kwh,error_pct"
    )?;

    for (spec, annual) in cases {
        for metric in ["heating", "cooling"] {
            let sim_monthly: &[f64; 12] = match metric {
                "heating" => &annual.monthly_heating_kwh,
                "cooling" => &annual.monthly_cooling_kwh,
                _ => unreachable!(),
            };
            let reference = reference_monthly(metric, spec.ref_kind);

            for i in 0..12 {
                let sim = sim_monthly[i];
                let (ref_v, err, pct) = if let Some(reference) = reference {
                    let ref_v = reference[i];
                    let err = sim - ref_v;
                    let pct = if ref_v != 0.0 {
                        100.0 * err / ref_v
                    } else {
                        0.0
                    };
                    (ref_v, err, pct)
                } else {
                    (f64::NAN, f64::NAN, f64::NAN)
                };

                writeln!(
                    f,
                    "{},{},{},{:.3},{:.3},{:.3},{:.2}",
                    spec.name,
                    i + 1,
                    metric,
                    sim,
                    ref_v,
                    err,
                    pct
                )?;
            }

            if let Some(reference) = reference {
                let sim_a = sim_monthly.iter().sum::<f64>();
                let ref_a = reference.iter().sum::<f64>();
                let err = sim_a - ref_a;
                let pct = if ref_a != 0.0 {
                    100.0 * err / ref_a
                } else {
                    0.0
                };
                writeln!(
                    f,
                    "{},{},{},{:.3},{:.3},{:.3},{:.2}",
                    spec.name, "annual", metric, sim_a, ref_a, err, pct
                )?;
            } else {
                let sim_a = sim_monthly.iter().sum::<f64>();
                writeln!(
                    f,
                    "{},{},{},{:.3},{},{},{}",
                    spec.name, "annual", metric, sim_a, "", "", ""
                )?;
            }
        }
    }

    Ok(())
}

#[derive(Debug, Clone, Copy)]
struct UaBreakdown {
    ua_total: f64,
    ua_glazing: f64,
    ua_opaque: f64,
    ua_ground: f64,
}

fn ua_breakdown(
    building: &Building,
    cfg: &ThermalConfig,
    solar: &SolarGainConfig,
) -> Result<UaBreakdown> {
    use building3d::sim::energy::boundary::ThermalBoundaries;
    use building3d::sim::index::SurfaceIndex;

    let index = SurfaceIndex::new(building);
    let boundaries = ThermalBoundaries::classify(building, &index);

    let mut ua_total = 0.0;
    let mut ua_glazing = 0.0;
    let mut ua_opaque = 0.0;
    let mut ua_ground = 0.0;

    for s in &index.surfaces {
        if !boundaries.is_exterior(&s.polygon_uid) {
            continue;
        }
        let u = cfg.resolve_u_value_for_surface(&s.polygon_uid, &s.path);
        let ua = u * s.area_m2;
        ua_total += ua;
        if solar
            .resolve_shgc(&s.path, cfg.material_library.as_ref())
            .is_some()
        {
            ua_glazing += ua;
        } else {
            ua_opaque += ua;
        }

        if cfg.ground_temperature_c.is_some()
            && cfg
                .ground_surface_patterns
                .iter()
                .any(|p| s.path.contains(p.as_str()))
        {
            let Some(poly) = building.get_polygon(&s.path) else {
                continue;
            };
            if poly.vn.dz <= -0.5 {
                ua_ground += ua;
            }
        }
    }

    Ok(UaBreakdown {
        ua_total,
        ua_glazing,
        ua_opaque,
        ua_ground,
    })
}

fn print_ua_breakdown(
    building: &Building,
    cfg: &ThermalConfig,
    solar: &SolarGainConfig,
) -> Result<()> {
    let ua = ua_breakdown(building, cfg, solar)?;
    println!(
        "UA breakdown: total={:.2} W/K (glazing={:.2}, opaque={:.2}, ground={:.2})",
        ua.ua_total, ua.ua_glazing, ua.ua_opaque, ua.ua_ground
    );
    Ok(())
}

#[derive(Debug, Clone)]
struct CaseDiagnostics {
    ua: UaBreakdown,
    annual_solar_transmitted_kwh: f64,
    annual_solar_opaque_sol_air_kwh: f64,
    annual_ground_correction_kwh: f64,
    monthly_solar_transmitted_kwh: [f64; 12],
    monthly_solar_opaque_sol_air_kwh: [f64; 12],
    monthly_ground_correction_kwh: [f64; 12],
    peak_heating_hour_idx: usize,
    peak_cooling_hour_idx: usize,
}

fn compute_opaque_sol_air_gain_total_w(
    building: &Building,
    cfg: &ThermalConfig,
    solar: &SolarGainConfig,
    index: &building3d::sim::index::SurfaceIndex,
    boundaries: &building3d::sim::energy::boundary::ThermalBoundaries,
    params: &SolarHourParams,
) -> f64 {
    use building3d::sim::lighting::solar::SolarPosition;

    let solar_pos = SolarPosition::calculate_from_local_time(
        params.latitude,
        params.longitude,
        params.timezone,
        params.day_of_year,
        params.local_time_hours,
    );
    let sun_above = solar_pos.is_above_horizon();
    let sun_dir = solar_pos.to_direction();

    let h_out = solar.exterior_heat_transfer_coeff_w_per_m2_k.max(1e-9);

    let mut total = 0.0;
    for surface in &index.surfaces {
        if !boundaries.is_exterior(&surface.polygon_uid) {
            continue;
        }
        if solar
            .resolve_shgc(&surface.path, cfg.material_library.as_ref())
            .is_some()
        {
            continue;
        }

        let area = surface.area_m2;
        if area <= 0.0 {
            continue;
        }
        let Some(poly) = building.get_polygon(&surface.path) else {
            continue;
        };

        let normal = poly.vn;
        let sky_view = 0.5 * (1.0 + normal.dz.max(0.0));

        let mut incident = params.diffuse_horizontal_irradiance.max(0.0) * sky_view;
        if sun_above && params.direct_normal_irradiance > 0.0 {
            let cos_incidence = sun_dir.dot(&normal).max(0.0);
            incident += params.direct_normal_irradiance.max(0.0) * cos_incidence;
        }
        if incident <= 0.0 {
            continue;
        }

        let absorptance = solar.default_opaque_absorptance.clamp(0.0, 1.0);
        if absorptance <= 0.0 {
            continue;
        }

        let u = cfg.resolve_u_value_for_surface(&surface.polygon_uid, &surface.path);
        if !(u.is_finite() && u > 0.0) {
            continue;
        }

        total += (u / h_out) * incident * absorptance * area;
    }
    total
}

fn compute_case_diagnostics(
    building: &Building,
    cfg: &ThermalConfig,
    solar_cfg: &SolarGainConfig,
    weather: &WeatherData,
    annual: &AnnualResult,
) -> Result<CaseDiagnostics> {
    use building3d::sim::energy::boundary::ThermalBoundaries;
    use building3d::sim::index::SurfaceIndex;

    let index = SurfaceIndex::new(building);
    let boundaries = ThermalBoundaries::classify(building, &index);
    let ua = ua_breakdown(building, cfg, solar_cfg)?;

    let mut solar_transmitted_wh = 0.0;
    let mut solar_opaque_wh = 0.0;
    let mut ground_correction_wh = 0.0;
    let mut monthly_solar_transmitted_wh = [0.0; 12];
    let mut monthly_solar_opaque_wh = [0.0; 12];
    let mut monthly_ground_correction_wh = [0.0; 12];

    for record in &weather.records {
        let month_idx = (record.month as usize).saturating_sub(1).min(11);
        let params = SolarHourParams {
            outdoor_air_temperature_c: record.dry_bulb_temperature,
            global_horizontal_irradiance: record.global_horizontal_radiation,
            direct_normal_irradiance: record.direct_normal_radiation,
            diffuse_horizontal_irradiance: record.diffuse_horizontal_radiation,
            horizontal_infrared_radiation: record.horizontal_infrared_radiation,
            wind_speed: record.wind_speed,
            day_of_year: day_of_year(record.month, record.day),
            local_time_hours: record.hour as f64 - 0.5,
            latitude: weather.latitude,
            longitude: weather.longitude,
            timezone: weather.timezone,
        };

        let transmitted_w = compute_solar_gains_with_materials(
            building,
            &params,
            solar_cfg,
            cfg.material_library.as_ref(),
        );
        solar_transmitted_wh += transmitted_w;
        monthly_solar_transmitted_wh[month_idx] += transmitted_w;

        if solar_cfg.include_exterior_opaque_absorption {
            let opaque_w = compute_opaque_sol_air_gain_total_w(
                building,
                cfg,
                solar_cfg,
                &index,
                &boundaries,
                &params,
            );
            solar_opaque_wh += opaque_w;
            monthly_solar_opaque_wh[month_idx] += opaque_w;
        }

        if let Some(tg) = cfg.ground_temperature_c {
            let tout = record.dry_bulb_temperature;
            let q_ground = ua.ua_ground * (tg - tout);
            ground_correction_wh += q_ground;
            monthly_ground_correction_wh[month_idx] += q_ground;
        }
    }

    let (peak_heating_hour_idx, _) = annual
        .hourly_heating
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.total_cmp(b))
        .unwrap_or((0, &0.0));
    let (peak_cooling_hour_idx, _) = annual
        .hourly_cooling
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.total_cmp(b))
        .unwrap_or((0, &0.0));

    Ok(CaseDiagnostics {
        ua,
        annual_solar_transmitted_kwh: solar_transmitted_wh / 1000.0,
        annual_solar_opaque_sol_air_kwh: solar_opaque_wh / 1000.0,
        annual_ground_correction_kwh: ground_correction_wh / 1000.0,
        monthly_solar_transmitted_kwh: monthly_solar_transmitted_wh.map(|v| v / 1000.0),
        monthly_solar_opaque_sol_air_kwh: monthly_solar_opaque_wh.map(|v| v / 1000.0),
        monthly_ground_correction_kwh: monthly_ground_correction_wh.map(|v| v / 1000.0),
        peak_heating_hour_idx,
        peak_cooling_hour_idx,
    })
}

fn weather_timestamp(weather: &WeatherData, hour_idx: usize) -> String {
    match weather.records.get(hour_idx) {
        Some(r) => format!("{:02}-{:02} {:02}:00", r.month, r.day, r.hour),
        None => "n/a".to_string(),
    }
}

fn write_diagnostics_csv(
    path: &Path,
    weather: &WeatherData,
    cases: &[(CaseSpec, AnnualResult, CaseDiagnostics)],
) -> Result<()> {
    let mut f = fs::File::create(path).context("create diagnostics.csv")?;
    writeln!(
        f,
        "case,ua_total_w_per_k,ua_glazing_w_per_k,ua_opaque_w_per_k,ua_ground_w_per_k,annual_solar_transmitted_kwh,annual_solar_opaque_sol_air_kwh,annual_ground_correction_kwh,peak_heating_w,peak_heating_time,peak_cooling_w,peak_cooling_time"
    )?;

    for (spec, annual, diag) in cases {
        writeln!(
            f,
            "{},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{},{:.3},{}",
            spec.name,
            diag.ua.ua_total,
            diag.ua.ua_glazing,
            diag.ua.ua_opaque,
            diag.ua.ua_ground,
            diag.annual_solar_transmitted_kwh,
            diag.annual_solar_opaque_sol_air_kwh,
            diag.annual_ground_correction_kwh,
            annual.peak_heating,
            weather_timestamp(weather, diag.peak_heating_hour_idx),
            annual.peak_cooling,
            weather_timestamp(weather, diag.peak_cooling_hour_idx),
        )?;
    }

    Ok(())
}

fn write_diagnostics_monthly_csv(path: &Path, cases: &[(CaseSpec, CaseDiagnostics)]) -> Result<()> {
    let mut f = fs::File::create(path).context("create diagnostics_monthly.csv")?;
    writeln!(f, "case,month,metric,value_kwh")?;

    for (spec, diag) in cases {
        for (metric, arr) in [
            ("solar_transmitted", &diag.monthly_solar_transmitted_kwh),
            (
                "solar_opaque_sol_air",
                &diag.monthly_solar_opaque_sol_air_kwh,
            ),
            ("ground_correction", &diag.monthly_ground_correction_kwh),
        ] {
            for (i, v) in arr.iter().enumerate() {
                writeln!(f, "{},{},{},{:.3}", spec.name, i + 1, metric, v)?;
            }
        }
    }

    Ok(())
}

fn main() -> Result<()> {
    let use_fvm_walls: bool = std::env::var("BESTEST_USE_FVM_WALLS")
        .ok()
        .as_deref()
        .map(|s| s == "1" || s.eq_ignore_ascii_case("true"))
        .unwrap_or(true);

    let epw_path = std::env::var("BESTEST_600_EPW")
        .map(PathBuf::from)
        .unwrap_or_else(|_| default_epw_path());
    let epw_content = fs::read_to_string(&epw_path)
        .with_context(|| format!("read EPW at {} (run download_data.sh?)", epw_path.display()))?;
    let weather = WeatherData::from_epw(&epw_content).context("parse EPW")?;

    let building = build_case_600_geometry()?;
    let mut base_cfg = config_for_case_600(&building);
    base_cfg.use_fvm_walls = use_fvm_walls;
    let solar_cfg = solar_config_for_case_600();
    let hvac = HvacIdealLoads::with_setpoints(20.0, 27.0);

    let enable_two_node_600: bool = std::env::var("BESTEST_600_ENABLE_TWO_NODE")
        .ok()
        .as_deref()
        .map(|s| s == "1" || s.eq_ignore_ascii_case("true"))
        .unwrap_or(false);
    let two_node_mass_fraction_600: f64 = std::env::var("BESTEST_600_TWO_NODE_MASS_FRACTION")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(0.6);
    let solar_to_mass_600: f64 = std::env::var("BESTEST_600_SOLAR_TO_MASS_FRACTION")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(0.6);
    let interior_h_600: f64 = std::env::var("BESTEST_600_INTERIOR_H_W_PER_M2_K")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(3.0);

    let suite = vec![
        CaseSpec {
            name: "600",
            ref_kind: RefKind::Ref600,
            mass_kind: CaseMassKind::Light,
            solar: true,
        },
        CaseSpec {
            name: "600_no_solar",
            ref_kind: RefKind::None,
            mass_kind: CaseMassKind::Light,
            solar: false,
        },
        CaseSpec {
            name: "900",
            ref_kind: RefKind::Ref900,
            mass_kind: CaseMassKind::Heavy,
            solar: true,
        },
        CaseSpec {
            name: "900_no_solar",
            ref_kind: RefKind::None,
            mass_kind: CaseMassKind::Heavy,
            solar: false,
        },
    ];

    println!("BESTEST energy suite (building3d vs OpenStudio/E+ reference)");
    println!(
        "Weather: {} ({} hours)",
        weather.location,
        weather.num_hours()
    );
    println!(
        "Reference annual: 600 h={:.1} c={:.1} | 900 h={:.1} c={:.1} (kWh)",
        sum(&REF_600_MONTHLY_HEATING_KWH),
        sum(&REF_600_MONTHLY_COOLING_KWH),
        sum(&REF_900_MONTHLY_HEATING_KWH),
        sum(&REF_900_MONTHLY_COOLING_KWH),
    );
    println!("FVM walls: {}", if use_fvm_walls { "ON" } else { "off" });
    if enable_two_node_600 {
        println!(
            "Light-mass 2R2C enabled (600): mass_fraction={:.2}, solar→mass={:.2}, interior_h={:.2} W/(m²·K)",
            two_node_mass_fraction_600, solar_to_mass_600, interior_h_600
        );
    }
    println!();

    let warmup_days: usize = std::env::var("BESTEST_WARMUP_DAYS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(7);
    let options = TransientSimulationOptions {
        warmup_hours: warmup_days.saturating_mul(24),
    };

    println!("Warmup: {warmup_days} days");
    print_ua_breakdown(&building, &base_cfg, &solar_cfg)?;
    println!();

    let mut outputs: Vec<(CaseSpec, AnnualResult)> = Vec::new();
    let mut diag_rows: Vec<(CaseSpec, AnnualResult, CaseDiagnostics)> = Vec::new();
    for spec in suite {
        let mut cfg = base_cfg.clone();
        if spec.mass_kind == CaseMassKind::Heavy {
            let (wall, roof, floor, window) = bestest_900_constructions();
            cfg.constructions.insert("window".to_string(), window);
            cfg.constructions.insert("ceiling".to_string(), roof);
            cfg.constructions.insert("floor".to_string(), floor);
            cfg.constructions.insert("wall".to_string(), wall);
        }
        if enable_two_node_600 && spec.name.starts_with("600") {
            cfg.two_node_mass_fraction = two_node_mass_fraction_600;
            cfg.interior_heat_transfer_coeff_w_per_m2_k = interior_h_600;
            cfg.solar_gains_to_mass_fraction = solar_to_mass_600;
            cfg.internal_gains_to_mass_fraction = 0.0;
        }

        cfg.thermal_capacity_j_per_m3_k = estimate_zone_capacity_j_per_m3_k(&building, &cfg);

        let solar = if spec.solar { Some(&solar_cfg) } else { None };
        let annual = run_transient_simulation_with_options(
            &building, &cfg, &weather, &hvac, None, solar, &options,
        );

        if spec.solar {
            let diag = compute_case_diagnostics(&building, &cfg, &solar_cfg, &weather, &annual)?;
            diag_rows.push((spec.clone(), annual.clone(), diag));
        }

        let sim_h = annual.annual_heating_kwh;
        let sim_c = annual.annual_cooling_kwh;

        print!(
            "Case {:>12}: heating={:8.1} kWh, cooling={:8.1} kWh",
            spec.name, sim_h, sim_c
        );
        if let Some(ref_h) = reference_monthly("heating", spec.ref_kind) {
            let ref_h = sum(&ref_h);
            let ref_c = sum(&reference_monthly("cooling", spec.ref_kind).unwrap());
            print!(
                " | vs ref: heating {:+6.1}%, cooling {:+6.1}%",
                100.0 * (sim_h - ref_h) / ref_h,
                100.0 * (sim_c - ref_c) / ref_c
            );
        }
        println!();

        outputs.push((spec, annual));
    }

    let out_path = repo_root().join("examples/bestest_energy_suite/results.csv");
    write_suite_csv(&out_path, &outputs)?;
    println!();
    println!("Wrote {}", out_path.display());
    let diag_path = repo_root().join("examples/bestest_energy_suite/diagnostics.csv");
    write_diagnostics_csv(&diag_path, &weather, &diag_rows)?;
    println!("Wrote {}", diag_path.display());
    let diag_monthly_path =
        repo_root().join("examples/bestest_energy_suite/diagnostics_monthly.csv");
    let diag_monthly_rows: Vec<(CaseSpec, CaseDiagnostics)> = diag_rows
        .iter()
        .map(|(spec, _annual, diag)| (spec.clone(), diag.clone()))
        .collect();
    write_diagnostics_monthly_csv(&diag_monthly_path, &diag_monthly_rows)?;
    println!("Wrote {}", diag_monthly_path.display());
    Ok(())
}
