use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};

use building3d::sim::energy::config::ThermalConfig;
use building3d::sim::energy::construction::WallConstruction;
use building3d::sim::energy::hvac::HvacIdealLoads;
use building3d::sim::energy::simulation::{AnnualResult, run_transient_simulation};
use building3d::sim::energy::solar_bridge::SolarGainConfig;
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
                thickness: 0.066,
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
                thickness: 0.1118,
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
        None,
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

    cfg.thermal_capacity_j_per_m3_k = estimate_zone_capacity_j_per_m3_k(building, &cfg);
    cfg
}

fn solar_config_for_case_600() -> SolarGainConfig {
    let mut solar = SolarGainConfig::new();
    solar.glazing_patterns = vec!["window".to_string()];
    // From `Glass Type 1` solar transmittance in the IDF (approximate SHGC).
    solar.default_shgc = 0.86156;
    solar
}

#[derive(Debug, Clone, Copy)]
enum RefKind {
    Ref600,
    Ref900,
    None,
}

#[derive(Debug, Clone)]
struct CaseSpec {
    name: &'static str,
    ref_kind: RefKind,
    high_mass_capacity_scale: f64,
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

fn main() -> Result<()> {
    let epw_path = std::env::var("BESTEST_600_EPW")
        .map(PathBuf::from)
        .unwrap_or_else(|_| default_epw_path());
    let epw_content = fs::read_to_string(&epw_path)
        .with_context(|| format!("read EPW at {} (run download_data.sh?)", epw_path.display()))?;
    let weather = WeatherData::from_epw(&epw_content).context("parse EPW")?;

    let building = build_case_600_geometry()?;
    let base_cfg = config_for_case_600(&building);
    let solar_cfg = solar_config_for_case_600();
    let hvac = HvacIdealLoads::with_setpoints(20.0, 27.0);

    let cap_scale_900: f64 = std::env::var("BESTEST_900_CAPACITY_SCALE")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(8.0);

    let suite = vec![
        CaseSpec {
            name: "600",
            ref_kind: RefKind::Ref600,
            high_mass_capacity_scale: 1.0,
            solar: true,
        },
        CaseSpec {
            name: "600_no_solar",
            ref_kind: RefKind::None,
            high_mass_capacity_scale: 1.0,
            solar: false,
        },
        CaseSpec {
            name: "900",
            ref_kind: RefKind::Ref900,
            high_mass_capacity_scale: cap_scale_900,
            solar: true,
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
    println!("High-mass capacity scale (900): {cap_scale_900}");
    println!();

    let mut outputs: Vec<(CaseSpec, AnnualResult)> = Vec::new();
    for spec in suite {
        let mut cfg = base_cfg.clone();
        cfg.thermal_capacity_j_per_m3_k *= spec.high_mass_capacity_scale.max(0.0);

        let solar = if spec.solar { Some(&solar_cfg) } else { None };
        let annual = run_transient_simulation(&building, &cfg, &weather, &hvac, None, solar);

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
    Ok(())
}
