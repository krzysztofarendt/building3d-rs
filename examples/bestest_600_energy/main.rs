use std::collections::HashSet;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, anyhow};
use building3d::sim::energy::config::{InternalMassBoundary, InternalMassSurface, ThermalConfig};
use building3d::sim::energy::construction::WallConstruction;
use building3d::sim::energy::convection::{ExteriorConvectionModel, InteriorConvectionModel};
use building3d::sim::energy::hvac::HvacIdealLoads;
use building3d::sim::energy::simulation::{
    AnnualResult, TransientSimulationOptions, run_transient_simulation_with_options,
};
use building3d::sim::energy::solar_bridge::SolarGainConfig;
use building3d::sim::energy::weather::WeatherData;
use building3d::sim::materials::Layer;
use building3d::{Building, Point, Polygon, Solid, Vector, Wall, Zone};

// Reference outputs from NREL/BESTEST-GSR (`results/workflow_results.csv`)
// line for "600 - Base Case" (OpenStudio/EnergyPlus).
//
// These are derived numeric results (monthly ideal loads, kWh).
const REF_MONTHLY_HEATING_KWH: [f64; 12] = [
    712.08, 682.83, 472.18, 509.87, 136.75, 10.08, 12.02, 6.58, 73.33, 347.66, 625.41, 735.97,
];
const REF_MONTHLY_COOLING_KWH: [f64; 12] = [
    523.18, 388.87, 485.77, 241.13, 340.91, 542.19, 537.21, 650.81, 735.08, 676.95, 438.99, 482.98,
];

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum SweepKind {
    LegacyFixedWithFloorMass,
    SuiteBaseline,
    SuitePlusFloorMass,
    SuiteFixedConvection,
    SuiteNoInteriorRadiation,
    SuiteViewFactorRadiation,
    SuiteNoBeamSplit,
    SuiteNoSurfaceAwareSolar,
    SuiteTransmittedToAir30Pct,
    SuiteDistributeToFvmWalls,
    SuiteDistributeToFvmWallsAirReroute,
    SuiteInteriorSolarAbsorptance100Pct,
    SuiteNoWarmupSubsteps1,
}

#[derive(Clone, Copy)]
struct SweepCase {
    id: &'static str,
    label: &'static str,
    description: &'static str,
    kind: SweepKind,
}

#[derive(Clone)]
struct SweepSummary {
    warmup_days: usize,
    substeps_per_hour: usize,
    use_view_factor_radiation: bool,
    interior_convection_model: InteriorConvectionModel,
    exterior_convection_model: ExteriorConvectionModel,
    use_interior_radiative_exchange: bool,
    use_surface_aware_solar_distribution: bool,
    use_beam_solar_distribution: bool,
    distribute_transmitted_solar_to_fvm_walls: bool,
    fvm_wall_solar_to_air: bool,
    transmitted_solar_to_air_fraction: f64,
    internal_gains_to_mass_fraction: f64,
    interior_solar_absorptance: f64,
    floor_internal_mass: bool,
}

#[derive(Clone)]
struct SweepResult {
    case: SweepCase,
    summary: SweepSummary,
    annual: AnnualResult,
    annual_heating_err_kwh: f64,
    annual_heating_err_pct: f64,
    annual_cooling_err_kwh: f64,
    annual_cooling_err_pct: f64,
}

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn default_epw_path() -> PathBuf {
    // Keep this aligned with the suite weather dataset location.
    repo_root()
        .join("examples/bestest_energy_suite/data")
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
                // BESTEST-GSR: thickness = 0.066 m, k = 0.04 -> RSI = 1.65 m2K/W.
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
                // Target RSI ~= 3.35 m2K/W (R-19 in IP units).
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

    // LTFLOOR: R-25 INSULATION (no mass) + TIMBER FLOORING.
    // In the BESTEST IDF this R-value is already in m2*K/W.
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
    // The layered construction is kept for capacity estimation; actual heat transfer
    // uses the U-value override below.
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
    Building::new("bestest_600", vec![zone])
}

fn suite_case_600_config(building: &Building) -> ThermalConfig {
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

    cfg.u_value_overrides_by_path_pattern
        .insert("window".to_string(), 2.8);

    cfg.ground_temperature_c = Some(10.0);

    cfg.use_surface_aware_solar_distribution = true;
    cfg.distribute_transmitted_solar_to_fvm_walls = false;
    cfg.use_beam_solar_distribution = true;
    cfg.fvm_wall_solar_to_air = false;
    cfg.transmitted_solar_to_air_fraction = 0.0;
    cfg.internal_gains_to_mass_fraction = 0.6;

    cfg.interior_convection_model = InteriorConvectionModel::Tarp;
    cfg.exterior_convection_model = ExteriorConvectionModel::Doe2;

    cfg.interior_heat_transfer_coeff_w_per_m2_k = 3.0;
    cfg.use_interior_radiative_exchange = true;
    cfg.interior_radiation_fraction = 0.6;

    cfg.use_view_factor_radiation = false;
    cfg.view_factor_rays_per_surface = 10_000;
    cfg.interior_emissivity = 0.9;

    cfg.interior_solar_absorptance = 0.6;

    cfg.thermal_capacity_j_per_m3_k = estimate_zone_capacity_j_per_m3_k(building, &cfg);
    cfg
}

fn case_600_solar_config() -> SolarGainConfig {
    let mut solar = SolarGainConfig::new();
    solar.glazing_patterns = vec!["window".to_string()];
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
    solar
}

fn base_options() -> TransientSimulationOptions {
    let warmup_days: usize = std::env::var("BESTEST_WARMUP_DAYS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(7);
    let substeps_per_hour: usize = std::env::var("BESTEST_SUBSTEPS_PER_HOUR")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(6)
        .max(1);
    TransientSimulationOptions {
        warmup_hours: warmup_days.saturating_mul(24),
        substeps_per_hour,
        ..Default::default()
    }
}

fn add_floor_internal_mass(building: &Building, cfg: &mut ThermalConfig) {
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
}

fn sweep_cases() -> Vec<SweepCase> {
    vec![
        SweepCase {
            id: "legacy_fixed_floor_mass",
            label: "Legacy fixed + floor mass",
            description: "bestest_600_energy style before suite-style upgrades",
            kind: SweepKind::LegacyFixedWithFloorMass,
        },
        SweepCase {
            id: "suite_baseline",
            label: "Suite baseline",
            description: "bestest_energy_suite Case 600 modeling style",
            kind: SweepKind::SuiteBaseline,
        },
        SweepCase {
            id: "suite_plus_floor_mass",
            label: "Suite + explicit floor mass",
            description: "add one-sided floor internal mass slab",
            kind: SweepKind::SuitePlusFloorMass,
        },
        SweepCase {
            id: "suite_fixed_convection",
            label: "Suite + fixed convection",
            description: "replace TARP/DOE2 with fixed interior/exterior films",
            kind: SweepKind::SuiteFixedConvection,
        },
        SweepCase {
            id: "suite_no_interior_radiation",
            label: "Suite + no interior radiation",
            description: "disable interior radiative split",
            kind: SweepKind::SuiteNoInteriorRadiation,
        },
        SweepCase {
            id: "suite_view_factor_radiation",
            label: "Suite + view-factor LW",
            description: "enable per-surface MRT view-factor radiation",
            kind: SweepKind::SuiteViewFactorRadiation,
        },
        SweepCase {
            id: "suite_no_beam_split",
            label: "Suite + no beam split",
            description: "disable beam-vs-diffuse interior solar handling",
            kind: SweepKind::SuiteNoBeamSplit,
        },
        SweepCase {
            id: "suite_no_surface_aware",
            label: "Suite + no surface-aware solar",
            description: "route transmitted solar by non-surface-aware path",
            kind: SweepKind::SuiteNoSurfaceAwareSolar,
        },
        SweepCase {
            id: "suite_transmit_to_air_30",
            label: "Suite + 30% transmitted solar to air",
            description: "shift part of transmitted solar directly to air node",
            kind: SweepKind::SuiteTransmittedToAir30Pct,
        },
        SweepCase {
            id: "suite_distribute_to_fvm_walls",
            label: "Suite + distribute to FVM walls",
            description: "include FVM wall area for transmitted/radiant source splitting",
            kind: SweepKind::SuiteDistributeToFvmWalls,
        },
        SweepCase {
            id: "suite_wall_share_to_air",
            label: "Suite + wall share rerouted to air",
            description: "when distributing to walls, redirect wall share into air gains",
            kind: SweepKind::SuiteDistributeToFvmWallsAirReroute,
        },
        SweepCase {
            id: "suite_alpha_1p0",
            label: "Suite + alpha = 1.0",
            description: "100% first-hit interior solar absorptance (no bounce pool)",
            kind: SweepKind::SuiteInteriorSolarAbsorptance100Pct,
        },
        SweepCase {
            id: "suite_no_warmup_substeps1",
            label: "Suite + no warmup + 1 substep",
            description: "coarse solver controls: 0 warmup days and 1 substep per hour",
            kind: SweepKind::SuiteNoWarmupSubsteps1,
        },
    ]
}

fn configure_case(
    building: &Building,
    case: SweepCase,
    cfg: &mut ThermalConfig,
    _solar: &mut SolarGainConfig,
    options: &mut TransientSimulationOptions,
) {
    match case.kind {
        SweepKind::LegacyFixedWithFloorMass => {
            cfg.interior_convection_model = InteriorConvectionModel::Fixed(3.0);
            cfg.exterior_convection_model = ExteriorConvectionModel::Fixed;
            cfg.use_beam_solar_distribution = false;
            add_floor_internal_mass(building, cfg);
        }
        SweepKind::SuiteBaseline => {}
        SweepKind::SuitePlusFloorMass => {
            add_floor_internal_mass(building, cfg);
        }
        SweepKind::SuiteFixedConvection => {
            cfg.interior_convection_model = InteriorConvectionModel::Fixed(3.0);
            cfg.exterior_convection_model = ExteriorConvectionModel::Fixed;
        }
        SweepKind::SuiteNoInteriorRadiation => {
            cfg.use_interior_radiative_exchange = false;
            cfg.use_view_factor_radiation = false;
        }
        SweepKind::SuiteViewFactorRadiation => {
            cfg.use_view_factor_radiation = true;
            cfg.view_factor_rays_per_surface = std::env::var("BESTEST_VF_RAYS")
                .ok()
                .and_then(|s| s.parse::<usize>().ok())
                .unwrap_or(10_000)
                .max(200);
        }
        SweepKind::SuiteNoBeamSplit => {
            cfg.use_beam_solar_distribution = false;
        }
        SweepKind::SuiteNoSurfaceAwareSolar => {
            cfg.use_surface_aware_solar_distribution = false;
        }
        SweepKind::SuiteTransmittedToAir30Pct => {
            cfg.transmitted_solar_to_air_fraction = 0.3;
        }
        SweepKind::SuiteDistributeToFvmWalls => {
            cfg.distribute_transmitted_solar_to_fvm_walls = true;
            cfg.fvm_wall_solar_to_air = false;
        }
        SweepKind::SuiteDistributeToFvmWallsAirReroute => {
            cfg.distribute_transmitted_solar_to_fvm_walls = true;
            cfg.fvm_wall_solar_to_air = true;
        }
        SweepKind::SuiteInteriorSolarAbsorptance100Pct => {
            cfg.interior_solar_absorptance = 1.0;
        }
        SweepKind::SuiteNoWarmupSubsteps1 => {
            options.warmup_hours = 0;
            options.substeps_per_hour = 1;
        }
    }
}

fn case_filter(all_cases: &[SweepCase]) -> Result<Vec<SweepCase>> {
    let Some(filter_raw) = std::env::var("BESTEST_SWEEP_CASES").ok() else {
        return Ok(all_cases.to_vec());
    };

    let wanted: HashSet<String> = filter_raw
        .split(',')
        .map(str::trim)
        .filter(|s| !s.is_empty())
        .map(str::to_string)
        .collect();

    if wanted.is_empty() {
        return Ok(all_cases.to_vec());
    }

    let mut selected = Vec::new();
    for c in all_cases {
        if wanted.contains(c.id) {
            selected.push(*c);
        }
    }

    if selected.is_empty() {
        return Err(anyhow!(
            "BESTEST_SWEEP_CASES matched no cases. Available ids: {}",
            all_cases
                .iter()
                .map(|c| c.id)
                .collect::<Vec<_>>()
                .join(", ")
        ));
    }

    Ok(selected)
}

fn summarize(cfg: &ThermalConfig, options: &TransientSimulationOptions) -> SweepSummary {
    SweepSummary {
        warmup_days: options.warmup_hours / 24,
        substeps_per_hour: options.substeps_per_hour,
        use_view_factor_radiation: cfg.use_view_factor_radiation,
        interior_convection_model: cfg.interior_convection_model,
        exterior_convection_model: cfg.exterior_convection_model,
        use_interior_radiative_exchange: cfg.use_interior_radiative_exchange,
        use_surface_aware_solar_distribution: cfg.use_surface_aware_solar_distribution,
        use_beam_solar_distribution: cfg.use_beam_solar_distribution,
        distribute_transmitted_solar_to_fvm_walls: cfg.distribute_transmitted_solar_to_fvm_walls,
        fvm_wall_solar_to_air: cfg.fvm_wall_solar_to_air,
        transmitted_solar_to_air_fraction: cfg.transmitted_solar_to_air_fraction,
        internal_gains_to_mass_fraction: cfg.internal_gains_to_mass_fraction,
        interior_solar_absorptance: cfg.interior_solar_absorptance,
        floor_internal_mass: !cfg.internal_mass_surfaces.is_empty(),
    }
}

fn interior_model_name(m: InteriorConvectionModel) -> &'static str {
    match m {
        InteriorConvectionModel::Fixed(_) => "fixed",
        InteriorConvectionModel::Tarp => "tarp",
    }
}

fn exterior_model_name(m: ExteriorConvectionModel) -> &'static str {
    match m {
        ExteriorConvectionModel::Fixed => "fixed",
        ExteriorConvectionModel::Doe2 => "doe2",
    }
}

fn write_monthly_results_csv(path: &Path, annual: &AnnualResult) -> Result<()> {
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

fn write_feature_sweep_csv(path: &Path, rows: &[SweepResult], suite_ref_idx: usize) -> Result<()> {
    let mut f = fs::File::create(path).context("create feature_sweep.csv")?;
    writeln!(
        f,
        "case_id,case_label,description,warmup_days,substeps_per_hour,use_view_factor_radiation,interior_convection_model,exterior_convection_model,use_interior_radiative_exchange,use_surface_aware_solar_distribution,use_beam_solar_distribution,distribute_transmitted_solar_to_fvm_walls,fvm_wall_solar_to_air,transmitted_solar_to_air_fraction,internal_gains_to_mass_fraction,interior_solar_absorptance,floor_internal_mass,annual_heating_kwh,annual_cooling_kwh,err_heating_kwh,err_heating_pct,err_cooling_kwh,err_cooling_pct,delta_heating_vs_suite_kwh,delta_cooling_vs_suite_kwh"
    )?;

    let suite_h = rows[suite_ref_idx].annual.annual_heating_kwh;
    let suite_c = rows[suite_ref_idx].annual.annual_cooling_kwh;

    for r in rows {
        let delta_h = r.annual.annual_heating_kwh - suite_h;
        let delta_c = r.annual.annual_cooling_kwh - suite_c;
        writeln!(
            f,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{:.6},{:.6},{:.6},{},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3},{:.3}",
            r.case.id,
            r.case.label,
            r.case.description,
            r.summary.warmup_days,
            r.summary.substeps_per_hour,
            r.summary.use_view_factor_radiation,
            interior_model_name(r.summary.interior_convection_model),
            exterior_model_name(r.summary.exterior_convection_model),
            r.summary.use_interior_radiative_exchange,
            r.summary.use_surface_aware_solar_distribution,
            r.summary.use_beam_solar_distribution,
            r.summary.distribute_transmitted_solar_to_fvm_walls,
            r.summary.fvm_wall_solar_to_air,
            r.summary.transmitted_solar_to_air_fraction,
            r.summary.internal_gains_to_mass_fraction,
            r.summary.interior_solar_absorptance,
            r.summary.floor_internal_mass,
            r.annual.annual_heating_kwh,
            r.annual.annual_cooling_kwh,
            r.annual_heating_err_kwh,
            r.annual_heating_err_pct,
            r.annual_cooling_err_kwh,
            r.annual_cooling_err_pct,
            delta_h,
            delta_c,
        )?;
    }

    Ok(())
}

fn xml_escape(raw: &str) -> String {
    raw.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}

fn write_feature_sweep_barplot_svg(path: &Path, rows: &[SweepResult]) -> Result<()> {
    let width: i32 = 1800;
    let left_margin: i32 = 470;
    let right_margin: i32 = 120;
    let top_margin: i32 = 130;
    let bottom_margin: i32 = 80;
    let row_h: i32 = 34;
    let bar_h: i32 = 10;

    let chart_h = (rows.len() as i32) * row_h;
    let height = top_margin + chart_h + bottom_margin;
    let plot_w = width - left_margin - right_margin;

    let max_abs_err = rows
        .iter()
        .map(|r| {
            r.annual_heating_err_pct
                .abs()
                .max(r.annual_cooling_err_pct.abs())
        })
        .fold(5.0_f64, f64::max);
    let max_abs_err = (max_abs_err / 5.0).ceil() * 5.0;
    let x_min = -max_abs_err;
    let x_max = max_abs_err;

    let x_of = |v: f64| -> f64 {
        let t = if (x_max - x_min).abs() > 0.0 {
            (v - x_min) / (x_max - x_min)
        } else {
            0.5
        };
        left_margin as f64 + t.clamp(0.0, 1.0) * plot_w as f64
    };

    let mut f = fs::File::create(path).context("create feature_sweep_barplot.svg")?;

    writeln!(
        f,
        "<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>"
    )?;
    writeln!(
        f,
        "<rect x='0' y='0' width='{width}' height='{height}' fill='#fbfbfd'/>"
    )?;

    writeln!(
        f,
        "<text x='{}' y='44' font-family='monospace' font-size='30' fill='#111827'>BESTEST 600 Feature Sweep</text>",
        left_margin
    )?;
    writeln!(
        f,
        "<text x='{}' y='76' font-family='monospace' font-size='16' fill='#374151'>Bars show annual error vs EnergyPlus reference (percent). Red = heating, blue = cooling.</text>",
        left_margin
    )?;

    // Grid and tick labels.
    let tick_count = 8;
    for i in 0..=tick_count {
        let v = x_min + (x_max - x_min) * (i as f64 / tick_count as f64);
        let x = x_of(v);
        let stroke = if v.abs() < 1e-9 { "#111827" } else { "#d1d5db" };
        let stroke_w = if v.abs() < 1e-9 { 2 } else { 1 };
        writeln!(
            f,
            "<line x1='{:.2}' y1='{}' x2='{:.2}' y2='{}' stroke='{}' stroke-width='{}'/>",
            x,
            top_margin - 6,
            x,
            top_margin + chart_h + 6,
            stroke,
            stroke_w
        )?;
        writeln!(
            f,
            "<text x='{:.2}' y='{}' text-anchor='middle' font-family='monospace' font-size='13' fill='#4b5563'>{:+.1}%</text>",
            x,
            top_margin + chart_h + 28,
            v
        )?;
    }

    // Legend.
    let legend_x = width - 330;
    writeln!(
        f,
        "<rect x='{}' y='34' width='14' height='14' fill='#dc2626'/>",
        legend_x
    )?;
    writeln!(
        f,
        "<text x='{}' y='46' font-family='monospace' font-size='14' fill='#111827'>Heating error %</text>",
        legend_x + 24
    )?;
    writeln!(
        f,
        "<rect x='{}' y='58' width='14' height='14' fill='#2563eb'/>",
        legend_x
    )?;
    writeln!(
        f,
        "<text x='{}' y='70' font-family='monospace' font-size='14' fill='#111827'>Cooling error %</text>",
        legend_x + 24
    )?;

    for (i, row) in rows.iter().enumerate() {
        let y0 = top_margin + (i as i32) * row_h;
        let yc = y0 + row_h / 2;

        writeln!(
            f,
            "<text x='{}' y='{}' text-anchor='end' font-family='monospace' font-size='13' fill='#111827'>{}</text>",
            left_margin - 10,
            yc + 5,
            xml_escape(row.case.label)
        )?;

        let x_zero = x_of(0.0);

        // Heating bar (upper)
        let xh = x_of(row.annual_heating_err_pct);
        let hx = x_zero.min(xh);
        let hw = (xh - x_zero).abs();
        writeln!(
            f,
            "<rect x='{:.2}' y='{}' width='{:.2}' height='{}' fill='#dc2626' opacity='0.85'/>",
            hx,
            yc - bar_h - 2,
            hw.max(1.0),
            bar_h
        )?;

        // Cooling bar (lower)
        let xc = x_of(row.annual_cooling_err_pct);
        let cx = x_zero.min(xc);
        let cw = (xc - x_zero).abs();
        writeln!(
            f,
            "<rect x='{:.2}' y='{}' width='{:.2}' height='{}' fill='#2563eb' opacity='0.85'/>",
            cx,
            yc + 2,
            cw.max(1.0),
            bar_h
        )?;

        writeln!(
            f,
            "<text x='{:.2}' y='{}' font-family='monospace' font-size='11' fill='#7f1d1d'>{:+.1}%</text>",
            if xh >= x_zero { xh + 6.0 } else { xh - 48.0 },
            yc - 2,
            row.annual_heating_err_pct
        )?;
        writeln!(
            f,
            "<text x='{:.2}' y='{}' font-family='monospace' font-size='11' fill='#1e3a8a'>{:+.1}%</text>",
            if xc >= x_zero { xc + 6.0 } else { xc - 48.0 },
            yc + bar_h + 4,
            row.annual_cooling_err_pct
        )?;
    }

    writeln!(f, "</svg>")?;
    Ok(())
}

fn run_sweep_case(
    building: &Building,
    weather: &WeatherData,
    hvac: &HvacIdealLoads,
    case: SweepCase,
    ref_annual_h: f64,
    ref_annual_c: f64,
) -> SweepResult {
    let mut cfg = suite_case_600_config(building);
    let mut solar = case_600_solar_config();
    let mut options = base_options();

    configure_case(building, case, &mut cfg, &mut solar, &mut options);

    let annual = run_transient_simulation_with_options(
        building,
        &cfg,
        weather,
        hvac,
        None,
        Some(&solar),
        &options,
    );

    let heating_err_kwh = annual.annual_heating_kwh - ref_annual_h;
    let cooling_err_kwh = annual.annual_cooling_kwh - ref_annual_c;

    SweepResult {
        case,
        summary: summarize(&cfg, &options),
        annual,
        annual_heating_err_kwh: heating_err_kwh,
        annual_heating_err_pct: 100.0 * heating_err_kwh / ref_annual_h,
        annual_cooling_err_kwh: cooling_err_kwh,
        annual_cooling_err_pct: 100.0 * cooling_err_kwh / ref_annual_c,
    }
}

fn main() -> Result<()> {
    let building = build_bestest_600()?;

    let epw_path = std::env::var("BESTEST_600_EPW")
        .map(PathBuf::from)
        .unwrap_or_else(|_| default_epw_path());
    let epw_content = fs::read_to_string(&epw_path)
        .with_context(|| format!("read EPW at {} (run download_data.sh?)", epw_path.display()))?;
    let weather = WeatherData::from_epw(&epw_content).context("parse EPW")?;

    let hvac = HvacIdealLoads::with_setpoints(20.0, 27.0);

    let ref_annual_h = annual_kwh(&REF_MONTHLY_HEATING_KWH);
    let ref_annual_c = annual_kwh(&REF_MONTHLY_COOLING_KWH);

    let all_cases = sweep_cases();
    let cases = case_filter(&all_cases)?;

    println!("BESTEST 600 feature sweep (building3d vs OpenStudio/E+ reference)");
    println!(
        "Weather: {} ({} hours)",
        weather.location,
        weather.num_hours()
    );
    println!(
        "Reference annual: heating={:.1} kWh, cooling={:.1} kWh",
        ref_annual_h, ref_annual_c
    );
    println!("Sweep cases: {}", cases.len());
    println!();

    let mut rows = Vec::with_capacity(cases.len());
    for (i, case) in cases.iter().enumerate() {
        println!("[{}/{}] {}", i + 1, cases.len(), case.label);
        let row = run_sweep_case(
            &building,
            &weather,
            &hvac,
            *case,
            ref_annual_h,
            ref_annual_c,
        );
        println!(
            "  annual: heating={:.1} kWh ({:+.1}%), cooling={:.1} kWh ({:+.1}%)",
            row.annual.annual_heating_kwh,
            row.annual_heating_err_pct,
            row.annual.annual_cooling_kwh,
            row.annual_cooling_err_pct,
        );
        rows.push(row);
    }

    let suite_ref_idx = rows
        .iter()
        .position(|r| r.case.kind == SweepKind::SuiteBaseline)
        .context("suite baseline case missing from sweep")?;

    // Keep legacy monthly CSV behavior for quick continuity checks, but now from suite baseline.
    let baseline_monthly_path = repo_root().join("examples/bestest_600_energy/results.csv");
    write_monthly_results_csv(&baseline_monthly_path, &rows[suite_ref_idx].annual)?;

    let sweep_csv_path = repo_root().join("examples/bestest_600_energy/feature_sweep.csv");
    write_feature_sweep_csv(&sweep_csv_path, &rows, suite_ref_idx)?;

    let sweep_plot_path = repo_root().join("examples/bestest_600_energy/feature_sweep_barplot.svg");
    write_feature_sweep_barplot_svg(&sweep_plot_path, &rows)?;

    println!();
    println!("Wrote {}", baseline_monthly_path.display());
    println!("Wrote {}", sweep_csv_path.display());
    println!("Wrote {}", sweep_plot_path.display());

    Ok(())
}
