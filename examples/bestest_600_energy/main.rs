use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use building3d::sim::energy::config::{InternalMassBoundary, InternalMassSurface, ThermalConfig};
use building3d::sim::energy::construction::WallConstruction;
use building3d::sim::energy::hvac::HvacIdealLoads;
use building3d::sim::energy::simulation::{
    TransientSimulationOptions, run_transient_simulation_with_options,
};
use building3d::sim::energy::solar_bridge::SolarGainConfig;
use building3d::sim::energy::weather::WeatherData;
use building3d::sim::materials::Layer;
use building3d::{Building, Point, Polygon, Solid, Vector, Wall, Zone};

// Reference outputs from NREL/BESTEST-GSR (`results/workflow_results.csv`)
// line for "600 - Base Case" (OpenStudio/EnergyPlus).
//
// These are *derived* numeric results (monthly ideal loads, kWh).
const REF_MONTHLY_HEATING_KWH: [f64; 12] = [
    712.08, 682.83, 472.18, 509.87, 136.75, 10.08, 12.02, 6.58, 73.33, 347.66, 625.41, 735.97,
];
const REF_MONTHLY_COOLING_KWH: [f64; 12] = [
    523.18, 388.87, 485.77, 241.13, 340.91, 542.19, 537.21, 650.81, 735.08, 676.95, 438.99, 482.98,
];

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn default_epw_path() -> PathBuf {
    repo_root()
        .join("examples/bestest_600_energy/data")
        .join("USA_MA_Boston-Logan.Intl.AP.725090_TMY3.epw")
}

fn annual_kwh(monthly: &[f64; 12]) -> f64 {
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
                // BESTEST-GSR: thickness = 0.066 m, k = 0.04 → RSI = 1.65 m²K/W.
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

    // Double Pane Window: ASHRAE 140-2020 glass properties.
    // The layered construction is kept for capacity estimation; actual heat
    // transfer uses the U-value override below (ISO 15099 convection+radiation
    // in the gap cannot be represented by a simple conductivity layer).
    let window = WallConstruction::new(
        "Double Pane Window",
        vec![
            Layer {
                name: "Glass".to_string(),
                thickness: 0.003048,
                conductivity: 1.0,
                density: 2500.0,
                specific_heat: 750.0,
            },
            Layer {
                name: "Air gap (approx)".to_string(),
                thickness: 0.012,
                conductivity: 0.026,
                density: 0.0,
                specific_heat: 0.0,
            },
            Layer {
                name: "Glass".to_string(),
                thickness: 0.003048,
                conductivity: 1.0,
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
        Some(Vector::new(0.0, -1.0, 0.0)),
    )
}

/// Builds a BESTEST 600-like single-zone model.
///
/// Matches the zone box (8 x 6 x 2.7) and two south windows from the
/// BESTEST-GSR `case_en_600.idf` resource.
fn build_bestest_600() -> Result<Building> {
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

    // Floor (normal down)
    let poly_floor = Polygon::new("floor", vec![p0, p3, p2, p1], None)?;
    // Roof (normal up)
    let poly_roof = Polygon::new("ceiling", vec![p4, p5, p6, p7], None)?;

    // South wall (y=0) is split into opaque panels + window polygons to avoid overlap.
    //
    // Windows (from IDF):
    // - Window 1: x=0.5..3.5, z=0.2..2.2
    // - Window 2: x=4.5..7.5, z=0.2..2.2
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
    Building::new("bestest_600", vec![zone])
}

fn write_results_csv(
    path: &Path,
    annual: &building3d::sim::energy::simulation::AnnualResult,
) -> Result<()> {
    let mut f = fs::File::create(path).context("create results.csv")?;
    writeln!(
        f,
        "month,metric,simulated_kwh,reference_kwh,error_kwh,error_pct"
    )?;

    for (i, (sim, reference)) in annual
        .monthly_heating_kwh
        .iter()
        .zip(REF_MONTHLY_HEATING_KWH.iter())
        .enumerate()
    {
        let err = sim - reference;
        let pct = if *reference != 0.0 {
            100.0 * err / reference
        } else {
            0.0
        };
        writeln!(
            f,
            "{},{},{:.3},{:.3},{:.3},{:.2}",
            i + 1,
            "heating",
            sim,
            reference,
            err,
            pct
        )?;
    }

    for (i, (sim, reference)) in annual
        .monthly_cooling_kwh
        .iter()
        .zip(REF_MONTHLY_COOLING_KWH.iter())
        .enumerate()
    {
        let err = sim - reference;
        let pct = if *reference != 0.0 {
            100.0 * err / reference
        } else {
            0.0
        };
        writeln!(
            f,
            "{},{},{:.3},{:.3},{:.3},{:.2}",
            i + 1,
            "cooling",
            sim,
            reference,
            err,
            pct
        )?;
    }

    let ref_annual_h = annual_kwh(&REF_MONTHLY_HEATING_KWH);
    let ref_annual_c = annual_kwh(&REF_MONTHLY_COOLING_KWH);

    let err_h = annual.annual_heating_kwh - ref_annual_h;
    let err_c = annual.annual_cooling_kwh - ref_annual_c;
    let pct_h = 100.0 * err_h / ref_annual_h;
    let pct_c = 100.0 * err_c / ref_annual_c;

    writeln!(
        f,
        "annual,heating,{:.3},{:.3},{:.3},{:.2}",
        annual.annual_heating_kwh, ref_annual_h, err_h, pct_h
    )?;
    writeln!(
        f,
        "annual,cooling,{:.3},{:.3},{:.3},{:.2}",
        annual.annual_cooling_kwh, ref_annual_c, err_c, pct_c
    )?;

    Ok(())
}

fn main() -> Result<()> {
    let building = build_bestest_600()?;

    let epw_path = std::env::var("BESTEST_600_EPW")
        .map(PathBuf::from)
        .unwrap_or_else(|_| default_epw_path());
    let epw_content = fs::read_to_string(&epw_path)
        .with_context(|| format!("read EPW at {} (run download_data.sh?)", epw_path.display()))?;
    let weather = WeatherData::from_epw(&epw_content).context("parse EPW")?;

    let (lt_wall, lt_roof, lt_floor, window) = bestest_600_constructions();

    let mut cfg = ThermalConfig::new();
    cfg.default_u_value = 0.0; // fail closed: expect all surfaces to be matched by constructions
    cfg.infiltration_ach = 0.5;
    cfg.internal_gains = 200.0;
    cfg.indoor_temperature = 20.0; // initial state only

    // Assign constructions by deterministic substring match (longest wins).
    cfg.constructions.insert("window".to_string(), window);
    cfg.constructions.insert("ceiling".to_string(), lt_roof);
    cfg.constructions.insert("floor".to_string(), lt_floor);
    cfg.constructions.insert("wall".to_string(), lt_wall);

    // BESTEST window: whole-window U-value from ISO 15099 gap calculation.
    // Double-pane clear glass (3mm, eps=0.84) with 12mm air gap: center-of-glass
    // U ranges from ~2.8 (annual average) to ~3.0 (NFRC winter peak). Since
    // building3d uses a fixed U year-round, 2.8 best represents the annual average.
    // The layered construction using still-air conductivity gives only ~1.5,
    // missing convection + radiation in the gap.
    cfg.u_value_overrides_by_path_pattern
        .insert("window".to_string(), 2.8);

    // BESTEST case floors are ground-coupled in the reference model.
    cfg.ground_temperature_c = Some(10.0);

    // Surface-aware policy for transmitted solar + radiant internal gains.
    cfg.use_surface_aware_solar_distribution = true;
    cfg.distribute_transmitted_solar_to_fvm_walls = false;
    cfg.transmitted_solar_to_air_fraction = 0.0;
    cfg.internal_gains_to_mass_fraction = 0.6;

    // Interior radiative exchange (convective vs radiative split) for FVM walls/mass.
    cfg.interior_heat_transfer_coeff_w_per_m2_k = 3.0;
    cfg.use_interior_radiative_exchange = true;
    cfg.interior_radiation_fraction = 0.6;

    // Model the floor as an internal mass slab (one-sided, insulated/adiabatic underside).
    let floor_area_m2 = building
        .get_polygon("zone_one/space/floor/floor")
        .map(|p| p.area())
        .unwrap_or(48.0);
    if let Some(floor) = cfg.constructions.get("floor").cloned() {
        cfg.internal_mass_surfaces.push(InternalMassSurface {
            name: "floor_mass".to_string(),
            zone_path_pattern: "zone_one".to_string(),
            area_m2: floor_area_m2,
            construction: floor,
            boundary: InternalMassBoundary::OneSidedAdiabatic,
            cos_tilt: -1.0,
        });
    }

    // Derive a zone capacity estimate from envelope layer capacities.
    cfg.thermal_capacity_j_per_m3_k = estimate_zone_capacity_j_per_m3_k(&building, &cfg);

    let hvac = HvacIdealLoads::with_setpoints(20.0, 27.0);

    let mut solar = SolarGainConfig::new();
    solar.glazing_patterns = vec!["window".to_string()];
    // BESTEST Glass Type 1 solar transmittance (single-pane, ASHRAE 140-2020).
    // Used as an approximation of the double-pane assembly SHGC (~0.76 at
    // normal incidence) because the single-pane angular polynomial already
    // accounts for most of the angular transmittance reduction.
    solar.default_shgc = 0.86156;
    solar.include_ground_reflection = true;
    solar.ground_reflectance = 0.2;
    solar.include_incidence_angle_modifier = true;
    solar.incidence_angle_modifier_a = 0.1;
    solar.include_exterior_opaque_absorption = true;
    solar.default_opaque_absorptance = 0.6;
    solar.include_exterior_longwave_exchange = true;
    solar.use_wind_speed_for_h_out = false;
    solar.angular_shgc_coefficients = Some(SolarGainConfig::single_pane_clear_coefficients());

    let warmup_days: usize = 7;
    let substeps_per_hour: usize = 6;
    let options = TransientSimulationOptions {
        warmup_hours: warmup_days * 24,
        substeps_per_hour,
    };

    println!("BESTEST 600 energy benchmark (building3d vs OpenStudio/E+ reference)");
    println!(
        "Weather: {} ({} hours)",
        weather.location,
        weather.num_hours()
    );
    println!(
        "Estimated zone capacity: {:.0} kJ/(m3·K)",
        cfg.thermal_capacity_j_per_m3_k / 1000.0
    );
    println!("Warmup: {warmup_days} days, substeps/hour: {substeps_per_hour}");
    println!(
        "Reference annual: heating={:.1} kWh, cooling={:.1} kWh",
        annual_kwh(&REF_MONTHLY_HEATING_KWH),
        annual_kwh(&REF_MONTHLY_COOLING_KWH)
    );
    println!();

    let annual = run_transient_simulation_with_options(
        &building, &cfg, &weather, &hvac, None, Some(&solar), &options,
    );

    let ref_annual_h = annual_kwh(&REF_MONTHLY_HEATING_KWH);
    let ref_annual_c = annual_kwh(&REF_MONTHLY_COOLING_KWH);

    println!(
        "Simulated annual: heating={:.1} kWh ({:+.1}%), cooling={:.1} kWh ({:+.1}%)",
        annual.annual_heating_kwh,
        100.0 * (annual.annual_heating_kwh - ref_annual_h) / ref_annual_h,
        annual.annual_cooling_kwh,
        100.0 * (annual.annual_cooling_kwh - ref_annual_c) / ref_annual_c
    );

    let out_path = repo_root().join("examples/bestest_600_energy/results.csv");
    write_results_csv(&out_path, &annual)?;
    println!("Wrote {}", out_path.display());

    Ok(())
}
