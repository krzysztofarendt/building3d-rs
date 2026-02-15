use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};

use building3d::sim::energy::config::ThermalConfig;
use building3d::sim::energy::construction::WallConstruction;
use building3d::sim::energy::convection::{ExteriorConvectionModel, InteriorConvectionModel};
use building3d::sim::energy::hvac::HvacIdealLoads;
use building3d::sim::energy::shading::{FinGeometry, OverhangGeometry};
use building3d::sim::energy::simulation::{
    AnnualResult, TransientSimulationOptions, run_transient_simulation_with_options,
};
use building3d::sim::energy::solar_bridge::{
    SolarGainConfig, SolarHourParams, WindowShading, compute_solar_gains_with_materials,
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
                // BESTEST-GSR: thickness = 0.066 m, k = 0.04 → RSI = 1.65 m²K/W.
                // Previously used 0.0776 m (RSI 1.94 from R-11 IP), which
                // over-insulated Case 600 walls by ~16% vs the reference.
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

fn config_for_case_600(_building: &Building) -> ThermalConfig {
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
    //
    // With FVM envelope walls, routing all transmitted shortwave into the zone air node
    // over-predicts peaks and fails to use the envelope thermal mass. Instead, deposit
    // most of it on interior surfaces (floor-first), and similarly deposit the radiant
    // fraction of internal gains to surfaces.
    cfg.use_surface_aware_solar_distribution = true;
    // EnergyPlus FullInteriorAndExterior distributes ALL transmitted solar to interior
    // surfaces (0% direct-to-air). Include FVM wall areas in the distribution to dilute
    // floor flux and spread solar across all interior surfaces including walls.
    cfg.distribute_transmitted_solar_to_fvm_walls = false;
    cfg.use_beam_solar_distribution = true;
    cfg.fvm_wall_solar_to_air = false;
    cfg.transmitted_solar_to_air_fraction = 0.0;
    cfg.internal_gains_to_mass_fraction = 0.6; // from BESTEST-GSR "OtherEquipment" radiant fraction

    // Dynamic interior convection (TARP/Walton): h depends on dT and surface tilt.
    cfg.interior_convection_model = InteriorConvectionModel::Tarp;
    // Dynamic exterior convection (DOE-2): combined natural + wind-forced.
    cfg.exterior_convection_model = ExteriorConvectionModel::Doe2;

    // Use a representative combined interior coefficient for explicit internal mass slabs.
    cfg.interior_heat_transfer_coeff_w_per_m2_k = 3.0;
    // Approximate interior longwave exchange (convective vs radiative split) for FVM walls/mass.
    cfg.use_interior_radiative_exchange = true;
    cfg.interior_radiation_fraction = 0.6;

    // View-factor interior longwave radiation exchange (per-surface MRT with uniform h_rad).
    cfg.use_view_factor_radiation = std::env::var("BESTEST_VF")
        .ok()
        .as_deref()
        .map(|s| s == "1" || s.eq_ignore_ascii_case("true"))
        .unwrap_or(false);
    cfg.view_factor_rays_per_surface = 10_000;
    cfg.interior_emissivity = 0.9;

    // Solar absorption/reflection: α=1.0 disables multi-bounce (100% absorbed on first hit).
    // Set to 0.6 for EnergyPlus-like model with 60% absorption, 40% reflected/redistributed.
    cfg.interior_solar_absorptance = std::env::var("BESTEST_SOLAR_ALPHA")
        .ok()
        .and_then(|s| s.parse::<f64>().ok())
        .unwrap_or(0.6);

    cfg.distribute_transmitted_solar_to_fvm_walls = false;

    // Floor is modeled as a layered FVM wall with ground-coupled exterior BC,
    // auto-collected from building geometry in collect_fvm_exterior_walls().
    cfg
}

fn solar_config_for_case_600() -> SolarGainConfig {
    let mut solar = SolarGainConfig::new();
    solar.glazing_patterns = vec!["window".to_string()];
    // BESTEST Glass Type 1 solar transmittance (single-pane, ASHRAE 140-2020).
    // Used as an approximation of the double-pane assembly SHGC (~0.76 at
    // normal incidence) because the single-pane angular polynomial already
    // accounts for most of the angular transmittance reduction.
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
    // Phase 2.2: exterior longwave exchange with sky/ground.
    solar.include_exterior_longwave_exchange = true;
    solar.use_wind_speed_for_h_out = false;
    // Phase 3: polynomial angular SHGC for single-pane clear glass.
    solar.angular_shgc_coefficients = Some(SolarGainConfig::single_pane_clear_coefficients());
    solar
}

// ── E+ annual reference totals (kWh and W) from BESTEST-GSR workflow_results.csv ──
#[allow(dead_code)]
struct AnnualRef {
    heating_kwh: f64,
    cooling_kwh: f64,
    peak_heating_w: f64,
    peak_cooling_w: f64,
}

struct FreefloatRef {
    min_temp_c: f64,
    max_temp_c: f64,
}

fn annual_ref(case: &str) -> Option<AnnualRef> {
    match case {
        "600" => Some(AnnualRef {
            heating_kwh: 4325.0,
            cooling_kwh: 6042.0,
            peak_heating_w: 3204.0,
            peak_cooling_w: 6352.0,
        }),
        "610" => Some(AnnualRef {
            heating_kwh: 4375.0,
            cooling_kwh: 4344.0,
            peak_heating_w: 3192.0,
            peak_cooling_w: 6137.0,
        }),
        "620" => Some(AnnualRef {
            heating_kwh: 4483.0,
            cooling_kwh: 4069.0,
            peak_heating_w: 3229.0,
            peak_cooling_w: 4797.0,
        }),
        "630" => Some(AnnualRef {
            heating_kwh: 4781.0,
            cooling_kwh: 2842.0,
            peak_heating_w: 3207.0,
            peak_cooling_w: 4212.0,
        }),
        "640" => Some(AnnualRef {
            heating_kwh: 2658.0,
            cooling_kwh: 5778.0,
            peak_heating_w: 4547.0,
            peak_cooling_w: 6299.0,
        }),
        "650" => Some(AnnualRef {
            heating_kwh: 0.0,
            cooling_kwh: 4839.0,
            peak_heating_w: 0.0,
            peak_cooling_w: 6141.0,
        }),
        "900" => Some(AnnualRef {
            heating_kwh: 1661.0,
            cooling_kwh: 2492.0,
            peak_heating_w: 2688.0,
            peak_cooling_w: 3042.0,
        }),
        "910" => Some(AnnualRef {
            heating_kwh: 1953.0,
            cooling_kwh: 1386.0,
            peak_heating_w: 2699.0,
            peak_cooling_w: 2224.0,
        }),
        "920" => Some(AnnualRef {
            heating_kwh: 3331.0,
            cooling_kwh: 2733.0,
            peak_heating_w: 2770.0,
            peak_cooling_w: 3261.0,
        }),
        "930" => Some(AnnualRef {
            heating_kwh: 3989.0,
            cooling_kwh: 1922.0,
            peak_heating_w: 2785.0,
            peak_cooling_w: 2782.0,
        }),
        "940" => Some(AnnualRef {
            heating_kwh: 1064.0,
            cooling_kwh: 2428.0,
            peak_heating_w: 3142.0,
            peak_cooling_w: 3041.0,
        }),
        "950" => Some(AnnualRef {
            heating_kwh: 0.0,
            cooling_kwh: 694.0,
            peak_heating_w: 0.0,
            peak_cooling_w: 2370.0,
        }),
        "960" => Some(AnnualRef {
            heating_kwh: 2703.0,
            cooling_kwh: 903.0,
            peak_heating_w: 2263.0,
            peak_cooling_w: 1480.0,
        }),
        _ => None,
    }
}

fn freefloat_ref(case: &str) -> Option<FreefloatRef> {
    match case {
        "600FF" => Some(FreefloatRef {
            min_temp_c: -12.56,
            max_temp_c: 63.92,
        }),
        "650FF" => Some(FreefloatRef {
            min_temp_c: -17.07,
            max_temp_c: 62.58,
        }),
        "900FF" => Some(FreefloatRef {
            min_temp_c: 1.25,
            max_temp_c: 44.30,
        }),
        "950FF" => Some(FreefloatRef {
            min_temp_c: -12.81,
            max_temp_c: 36.69,
        }),
        _ => None,
    }
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum GeometryKind {
    South,
    EastWest,
    Sunspace,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ShadingKind {
    None,
    Overhang,
    OverhangFins,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum HvacMode {
    Normal,
    Setback,
    NightVent,
    FreeFloat,
    FreeFloatNightVent,
}

#[derive(Debug, Clone)]
struct CaseSpec {
    name: &'static str,
    ref_kind: RefKind,
    mass_kind: CaseMassKind,
    solar: bool,
    geometry_kind: GeometryKind,
    shading: ShadingKind,
    hvac_mode: HvacMode,
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

// ── Case 620: East/west windows (no south windows) ──

fn poly_rect_east(name: &str, y0: f64, y1: f64, z0: f64, z1: f64, x: f64) -> Result<Polygon> {
    Polygon::new(
        name,
        vec![
            Point::new(x, y0, z0),
            Point::new(x, y1, z0),
            Point::new(x, y1, z1),
            Point::new(x, y0, z1),
        ],
        Some(building3d::Vector::new(1.0, 0.0, 0.0)),
    )
}

fn poly_rect_west(name: &str, y0: f64, y1: f64, z0: f64, z1: f64, x: f64) -> Result<Polygon> {
    Polygon::new(
        name,
        vec![
            Point::new(x, y1, z0),
            Point::new(x, y0, z0),
            Point::new(x, y0, z1),
            Point::new(x, y1, z1),
        ],
        Some(building3d::Vector::new(-1.0, 0.0, 0.0)),
    )
}

fn build_case_620_geometry() -> Result<Building> {
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

    let poly_floor = Polygon::new("floor", vec![p0, p3, p2, p1], Option::None)?;
    let poly_roof = Polygon::new("ceiling", vec![p4, p5, p6, p7], Option::None)?;

    // South wall: fully opaque (no windows)
    let wall_0 = Wall::new("wall_0", vec![poly_rect_y0("poly_0", 0.0, x, 0.0, z)?])?;

    // East wall (x=8): two windows, each 3m wide x 2m high, centered vertically
    // Total glazing = 2 * 3 * 2 = 12 m² (same as case 600)
    let mut east_polys = Vec::new();
    east_polys.push(poly_rect_east("opaque_bottom", 0.0, y, 0.0, 0.2, x)?);
    east_polys.push(poly_rect_east("opaque_top", 0.0, y, 2.2, z, x)?);
    east_polys.push(poly_rect_east("opaque_left", 0.0, 1.5, 0.2, 2.2, x)?);
    east_polys.push(poly_rect_east("window_1", 1.5, 4.5, 0.2, 2.2, x)?);
    east_polys.push(poly_rect_east("opaque_right", 4.5, y, 0.2, 2.2, x)?);
    let wall_1 = Wall::new("wall_1", east_polys)?;

    // North wall: opaque
    let wall_2 = Wall::new(
        "wall_2",
        vec![Polygon::new("poly_2", vec![p3, p7, p6, p2], Option::None)?],
    )?;

    // West wall (x=0): same window layout as east
    let mut west_polys = Vec::new();
    west_polys.push(poly_rect_west("opaque_bottom", 0.0, y, 0.0, 0.2, 0.0)?);
    west_polys.push(poly_rect_west("opaque_top", 0.0, y, 2.2, z, 0.0)?);
    west_polys.push(poly_rect_west("opaque_left", 0.0, 1.5, 0.2, 2.2, 0.0)?);
    west_polys.push(poly_rect_west("window_2", 1.5, 4.5, 0.2, 2.2, 0.0)?);
    west_polys.push(poly_rect_west("opaque_right", 4.5, y, 0.2, 2.2, 0.0)?);
    let wall_3 = Wall::new("wall_3", west_polys)?;

    let wall_floor = Wall::new("floor", vec![poly_floor])?;
    let wall_roof = Wall::new("ceiling", vec![poly_roof])?;

    let solid = Solid::new(
        "space",
        vec![wall_floor, wall_0, wall_1, wall_2, wall_3, wall_roof],
    )?;
    let zone = Zone::new("zone_one", vec![solid])?;
    Building::new("bestest_case_620", vec![zone])
}

// ── Case 960: Two-zone sunspace geometry ──

/// Helper to create a rectangle polygon at a given y position facing north (+y).
fn poly_rect_y_north(name: &str, x0: f64, x1: f64, z0: f64, z1: f64, y: f64) -> Result<Polygon> {
    Polygon::new(
        name,
        vec![
            Point::new(x1, y, z0),
            Point::new(x0, y, z0),
            Point::new(x0, y, z1),
            Point::new(x1, y, z1),
        ],
        Some(building3d::Vector::new(0.0, 1.0, 0.0)),
    )
}

/// Helper to create a rectangle polygon at a given y position facing south (-y).
fn poly_rect_y_south(name: &str, x0: f64, x1: f64, z0: f64, z1: f64, y: f64) -> Result<Polygon> {
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

fn build_case_960_geometry() -> Result<Building> {
    let x = 8.0;
    let z = 2.7;
    let y_sun = 2.0; // sunspace depth (south)
    let y_back = 6.0; // back zone depth (north)

    // ── Sunspace (south zone, UNCONDITIONED): 8m x 2m x 2.7m ──
    // y = 0..2, south wall at y=0, common wall at y=2
    let sun_floor = Polygon::new(
        "floor",
        vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(0.0, y_sun, 0.0),
            Point::new(x, y_sun, 0.0),
            Point::new(x, 0.0, 0.0),
        ],
        None,
    )?;
    let sun_ceiling = Polygon::new(
        "ceiling",
        vec![
            Point::new(0.0, 0.0, z),
            Point::new(x, 0.0, z),
            Point::new(x, y_sun, z),
            Point::new(0.0, y_sun, z),
        ],
        None,
    )?;

    // South wall (y=0): two 3m x 2m windows (same as Case 600/900 layout)
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

    // Common wall (y=2, facing north from sunspace side)
    let sun_common = poly_rect_y_north("common_wall", 0.0, x, 0.0, z, y_sun)?;

    // East wall (x=8)
    let sun_east = Polygon::new(
        "poly_east",
        vec![
            Point::new(x, 0.0, 0.0),
            Point::new(x, y_sun, 0.0),
            Point::new(x, y_sun, z),
            Point::new(x, 0.0, z),
        ],
        Some(building3d::Vector::new(1.0, 0.0, 0.0)),
    )?;
    // West wall (x=0)
    let sun_west = Polygon::new(
        "poly_west",
        vec![
            Point::new(0.0, y_sun, 0.0),
            Point::new(0.0, 0.0, 0.0),
            Point::new(0.0, 0.0, z),
            Point::new(0.0, y_sun, z),
        ],
        Some(building3d::Vector::new(-1.0, 0.0, 0.0)),
    )?;

    let sun_solid = Solid::new(
        "sunspace",
        vec![
            Wall::new("floor", vec![sun_floor])?,
            Wall::new("wall_0", south_polys)?,
            Wall::new("wall_1", vec![sun_east])?,
            Wall::new("wall_2", vec![sun_common])?,
            Wall::new("wall_3", vec![sun_west])?,
            Wall::new("ceiling", vec![sun_ceiling])?,
        ],
    )?;

    // ── Back zone (north, CONDITIONED): 8m x 6m x 2.7m ──
    // y = 2..8, common wall at y=2, north wall at y=8
    let y_north = y_sun + y_back; // 8.0
    let back_floor = Polygon::new(
        "floor",
        vec![
            Point::new(0.0, y_sun, 0.0),
            Point::new(0.0, y_north, 0.0),
            Point::new(x, y_north, 0.0),
            Point::new(x, y_sun, 0.0),
        ],
        None,
    )?;
    let back_ceiling = Polygon::new(
        "ceiling",
        vec![
            Point::new(0.0, y_sun, z),
            Point::new(x, y_sun, z),
            Point::new(x, y_north, z),
            Point::new(0.0, y_north, z),
        ],
        None,
    )?;

    // Common wall (y=2, facing south from back zone side)
    let back_common = poly_rect_y_south("common_wall", 0.0, x, 0.0, z, y_sun)?;

    // East wall (x=8)
    let back_east = Polygon::new(
        "poly_east",
        vec![
            Point::new(x, y_sun, 0.0),
            Point::new(x, y_north, 0.0),
            Point::new(x, y_north, z),
            Point::new(x, y_sun, z),
        ],
        Some(building3d::Vector::new(1.0, 0.0, 0.0)),
    )?;
    // North wall (y=8)
    let back_north = Polygon::new(
        "poly_north",
        vec![
            Point::new(x, y_north, 0.0),
            Point::new(0.0, y_north, 0.0),
            Point::new(0.0, y_north, z),
            Point::new(x, y_north, z),
        ],
        Some(building3d::Vector::new(0.0, 1.0, 0.0)),
    )?;
    // West wall (x=0)
    let back_west = Polygon::new(
        "poly_west",
        vec![
            Point::new(0.0, y_north, 0.0),
            Point::new(0.0, y_sun, 0.0),
            Point::new(0.0, y_sun, z),
            Point::new(0.0, y_north, z),
        ],
        Some(building3d::Vector::new(-1.0, 0.0, 0.0)),
    )?;

    let back_solid = Solid::new(
        "back",
        vec![
            Wall::new("floor", vec![back_floor])?,
            Wall::new("wall_0", vec![back_common])?,
            Wall::new("wall_1", vec![back_east])?,
            Wall::new("wall_2", vec![back_north])?,
            Wall::new("wall_3", vec![back_west])?,
            Wall::new("ceiling", vec![back_ceiling])?,
        ],
    )?;

    // Zone names sorted: z_back (idx 0), z_sun (idx 1)
    let z_sun = Zone::new("z_sun", vec![sun_solid])?;
    let z_back = Zone::new("z_back", vec![back_solid])?;
    Building::new("bestest_960", vec![z_sun, z_back])
}

// ── Overhang/fin shading configurations ──

/// BESTEST case 610/910: south overhang 1m deep, no gap.
fn south_overhang_shading(window_height: f64, window_width: f64) -> WindowShading {
    WindowShading {
        overhang: Some(OverhangGeometry {
            depth: 1.0,
            gap: 0.0,
            window_height,
            window_width,
            extension_left: 0.0,
            extension_right: 0.0,
        }),
        fin_left: Option::None,
        fin_right: Option::None,
        surface_azimuth_deg: 180.0, // south-facing
    }
}

/// BESTEST case 630/930: east/west windows with overhang (1m) and side fins (1m).
fn east_overhang_fins_shading(window_height: f64, window_width: f64) -> WindowShading {
    WindowShading {
        overhang: Some(OverhangGeometry {
            depth: 1.0,
            gap: 0.0,
            window_height,
            window_width,
            extension_left: 0.0,
            extension_right: 0.0,
        }),
        fin_left: Some(FinGeometry {
            depth: 1.0,
            height: window_height,
        }),
        fin_right: Some(FinGeometry {
            depth: 1.0,
            height: window_height,
        }),
        surface_azimuth_deg: 90.0, // east-facing
    }
}

fn west_overhang_fins_shading(window_height: f64, window_width: f64) -> WindowShading {
    WindowShading {
        overhang: Some(OverhangGeometry {
            depth: 1.0,
            gap: 0.0,
            window_height,
            window_width,
            extension_left: 0.0,
            extension_right: 0.0,
        }),
        fin_left: Some(FinGeometry {
            depth: 1.0,
            height: window_height,
        }),
        fin_right: Some(FinGeometry {
            depth: 1.0,
            height: window_height,
        }),
        surface_azimuth_deg: 270.0, // west-facing
    }
}

// ── HVAC schedules ──

/// Case 640/940: thermostat setback schedule.
/// Heating setpoint = 10C during 23:00-07:00, 20C otherwise.
fn make_setback_heating_schedule() -> Vec<f64> {
    (0..8760)
        .map(|h| {
            let hour_of_day = h % 24;
            if hour_of_day >= 23 || hour_of_day < 7 {
                10.0
            } else {
                20.0
            }
        })
        .collect()
}

/// Case 650/950: night ventilation schedule.
/// Night ventilation of 1703.16 m³/h from 18:00-07:00, plus base 0.5 ACH.
fn make_night_vent_infiltration_schedule(volume: f64) -> Vec<f64> {
    let night_ach = 1703.16 / volume; // ~13.14 ACH for BESTEST zone
    (0..8760)
        .map(|h| {
            let hour_of_day = h % 24;
            if hour_of_day >= 18 || hour_of_day < 7 {
                0.5 + night_ach
            } else {
                0.5
            }
        })
        .collect()
}

/// Case 650/950: cooling setpoint schedule (daytime only).
/// Cooling at 27C from 07:00-18:00, disabled (999C) otherwise.
fn make_daytime_cooling_schedule() -> Vec<f64> {
    (0..8760)
        .map(|h| {
            let hour_of_day = h % 24;
            if hour_of_day >= 7 && hour_of_day < 18 {
                27.0
            } else {
                999.0
            }
        })
        .collect()
}

/// Free-float cases (no HVAC): heating disabled, cooling disabled.
fn make_free_float_heating_schedule() -> Vec<f64> {
    vec![-999.0; 8760]
}

fn make_free_float_cooling_schedule() -> Vec<f64> {
    vec![999.0; 8760]
}

fn write_suite_csv(path: &Path, cases: &[(CaseSpec, AnnualResult)]) -> Result<()> {
    let mut f = fs::File::create(path).context("create results.csv")?;
    writeln!(
        f,
        "case,month,metric,simulated_kwh,reference_kwh,error_kwh,error_pct"
    )?;

    for (spec, annual) in cases {
        let is_ff = matches!(
            spec.hvac_mode,
            HvacMode::FreeFloat | HvacMode::FreeFloatNightVent
        );

        if is_ff {
            // Free-float cases: write min/max zone temperature instead of loads.
            let ff_ref = freefloat_ref(spec.name);
            let (ref_min, err_min, pct_min) = if let Some(r) = &ff_ref {
                let e = annual.min_zone_temp_c - r.min_temp_c;
                let p = if r.min_temp_c.abs() > 0.01 {
                    100.0 * e / r.min_temp_c.abs()
                } else {
                    0.0
                };
                (
                    format!("{:.3}", r.min_temp_c),
                    format!("{:.3}", e),
                    format!("{:.2}", p),
                )
            } else {
                ("".to_string(), "".to_string(), "".to_string())
            };
            writeln!(
                f,
                "{},annual,min_zone_temp_c,{:.3},{},{},{}",
                spec.name, annual.min_zone_temp_c, ref_min, err_min, pct_min
            )?;

            let (ref_max, err_max, pct_max) = if let Some(r) = &ff_ref {
                let e = annual.max_zone_temp_c - r.max_temp_c;
                let p = if r.max_temp_c.abs() > 0.01 {
                    100.0 * e / r.max_temp_c.abs()
                } else {
                    0.0
                };
                (
                    format!("{:.3}", r.max_temp_c),
                    format!("{:.3}", e),
                    format!("{:.2}", p),
                )
            } else {
                ("".to_string(), "".to_string(), "".to_string())
            };
            writeln!(
                f,
                "{},annual,max_zone_temp_c,{:.3},{},{},{}",
                spec.name, annual.max_zone_temp_c, ref_max, err_max, pct_max
            )?;
        } else {
            // HVAC cases: write monthly and annual heating/cooling loads.
            let a_ref = annual_ref(spec.name);
            let monthly_ref = |metric: &str, kind: RefKind| -> Option<[f64; 12]> {
                reference_monthly(metric, kind)
            };

            for metric in ["heating", "cooling"] {
                let sim_monthly: &[f64; 12] = match metric {
                    "heating" => &annual.monthly_heating_kwh,
                    "cooling" => &annual.monthly_cooling_kwh,
                    _ => unreachable!(),
                };
                let reference = monthly_ref(metric, spec.ref_kind);

                // Monthly rows (with reference if available).
                for i in 0..12 {
                    let sim = sim_monthly[i];
                    if let Some(ref_arr) = reference {
                        let ref_v = ref_arr[i];
                        let err = sim - ref_v;
                        let pct = if ref_v != 0.0 {
                            100.0 * err / ref_v
                        } else {
                            0.0
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
                    } else {
                        writeln!(f, "{},{},{},{:.3},,,", spec.name, i + 1, metric, sim,)?;
                    }
                }

                // Annual row (use annual_ref for all cases).
                let sim_a: f64 = sim_monthly.iter().sum();
                let ref_annual = a_ref.as_ref().map(|r| match metric {
                    "heating" => r.heating_kwh,
                    "cooling" => r.cooling_kwh,
                    _ => unreachable!(),
                });
                if let Some(ref_a) = ref_annual {
                    let err = sim_a - ref_a;
                    let pct = if ref_a != 0.0 {
                        100.0 * err / ref_a
                    } else {
                        0.0
                    };
                    writeln!(
                        f,
                        "{},annual,{},{:.3},{:.3},{:.3},{:.2}",
                        spec.name, metric, sim_a, ref_a, err, pct
                    )?;
                } else {
                    writeln!(f, "{},annual,{},{:.3},,,", spec.name, metric, sim_a,)?;
                }
            }

            // Peak loads row (with reference if available).
            if let Some(r) = &a_ref {
                let pk_h_err = annual.peak_heating - r.peak_heating_w;
                let pk_h_pct = if r.peak_heating_w > 0.0 {
                    100.0 * pk_h_err / r.peak_heating_w
                } else {
                    0.0
                };
                writeln!(
                    f,
                    "{},annual,peak_heating_w,{:.1},{:.1},{:.1},{:.2}",
                    spec.name, annual.peak_heating, r.peak_heating_w, pk_h_err, pk_h_pct
                )?;
                let pk_c_err = annual.peak_cooling - r.peak_cooling_w;
                let pk_c_pct = if r.peak_cooling_w > 0.0 {
                    100.0 * pk_c_err / r.peak_cooling_w
                } else {
                    0.0
                };
                writeln!(
                    f,
                    "{},annual,peak_cooling_w,{:.1},{:.1},{:.1},{:.2}",
                    spec.name, annual.peak_cooling, r.peak_cooling_w, pk_c_err, pk_c_pct
                )?;
            } else {
                writeln!(
                    f,
                    "{},annual,peak_heating_w,{:.1},,,",
                    spec.name, annual.peak_heating
                )?;
                writeln!(
                    f,
                    "{},annual,peak_cooling_w,{:.1},,,",
                    spec.name, annual.peak_cooling
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
        )
        .total();
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

fn make_case_spec(
    name: &'static str,
    mass: CaseMassKind,
    geom: GeometryKind,
    shading: ShadingKind,
    hvac_mode: HvacMode,
) -> CaseSpec {
    let ref_kind = match name {
        "600" => RefKind::Ref600,
        "900" => RefKind::Ref900,
        _ => RefKind::None,
    };
    CaseSpec {
        name,
        ref_kind,
        mass_kind: mass,
        solar: true,
        geometry_kind: geom,
        shading,
        hvac_mode,
    }
}

fn main() -> Result<()> {
    use CaseMassKind::*;
    use GeometryKind::*;
    use HvacMode::*;
    use ShadingKind as SK;

    let epw_path = std::env::var("BESTEST_600_EPW")
        .map(PathBuf::from)
        .unwrap_or_else(|_| default_epw_path());
    let epw_content = fs::read_to_string(&epw_path)
        .with_context(|| format!("read EPW at {} (run download_data.sh?)", epw_path.display()))?;
    let weather = WeatherData::from_epw(&epw_content).context("parse EPW")?;

    // Build geometry variants.
    let building_south = build_case_600_geometry()?;
    let building_ew = build_case_620_geometry()?;
    let building_960 = build_case_960_geometry()?;

    let base_cfg_south = {
        let mut c = config_for_case_600(&building_south);
        c.use_fvm_walls = true;
        c
    };
    let base_cfg_ew = {
        let mut c = config_for_case_600(&building_ew);
        c.use_fvm_walls = true;
        c
    };
    let solar_cfg = solar_config_for_case_600();
    let hvac = HvacIdealLoads::with_setpoints(20.0, 27.0);

    let suite: Vec<CaseSpec> = vec![
        // ── 600 series (lightweight) ──
        make_case_spec("600", Light, South, SK::None, Normal),
        make_case_spec("610", Light, South, SK::Overhang, Normal),
        make_case_spec("620", Light, EastWest, SK::None, Normal),
        make_case_spec("630", Light, EastWest, SK::OverhangFins, Normal),
        make_case_spec("640", Light, South, SK::None, Setback),
        make_case_spec("650", Light, South, SK::None, NightVent),
        make_case_spec("600FF", Light, South, SK::None, FreeFloat),
        make_case_spec("650FF", Light, South, SK::None, FreeFloatNightVent),
        // ── 900 series (heavyweight) ──
        make_case_spec("900", Heavy, South, SK::None, Normal),
        make_case_spec("910", Heavy, South, SK::Overhang, Normal),
        make_case_spec("920", Heavy, EastWest, SK::None, Normal),
        make_case_spec("930", Heavy, EastWest, SK::OverhangFins, Normal),
        make_case_spec("940", Heavy, South, SK::None, Setback),
        make_case_spec("950", Heavy, South, SK::None, NightVent),
        make_case_spec("900FF", Heavy, South, SK::None, FreeFloat),
        make_case_spec("950FF", Heavy, South, SK::None, FreeFloatNightVent),
        // ── Case 960: two-zone sunspace ──
        make_case_spec("960", Heavy, Sunspace, SK::None, Normal),
    ];

    println!("BESTEST ASHRAE 140 energy suite (building3d vs OpenStudio/E+ reference)");
    println!(
        "Weather: {} ({} hours)",
        weather.location,
        weather.num_hours()
    );
    println!(
        "Cases: {} total ({} HVAC + {} free-float)",
        suite.len(),
        suite
            .iter()
            .filter(|s| !matches!(s.hvac_mode, FreeFloat | FreeFloatNightVent))
            .count(),
        suite
            .iter()
            .filter(|s| matches!(s.hvac_mode, FreeFloat | FreeFloatNightVent))
            .count(),
    );
    println!("FVM walls: ON");
    println!("Global FVM solve: ON");
    println!();

    let warmup_days: usize = std::env::var("BESTEST_WARMUP_DAYS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(7);
    let substeps_per_hour: usize = std::env::var("BESTEST_SUBSTEPS_PER_HOUR")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(6)
        .max(1);

    println!("Warmup: {warmup_days} days");
    println!("Substeps per hour: {substeps_per_hour}");
    print_ua_breakdown(&building_south, &base_cfg_south, &solar_cfg)?;
    println!();

    let mut outputs: Vec<(CaseSpec, AnnualResult)> = Vec::new();
    let mut diag_rows: Vec<(CaseSpec, AnnualResult, CaseDiagnostics)> = Vec::new();

    // Case 960 config: sunspace uses heavyweight walls/floor (900-series)
    // for the sunspace zone, lightweight walls for the back zone, and a
    // concrete block partition for the common wall.
    let base_cfg_960 = {
        let mut c = config_for_case_600(&building_960);
        c.use_fvm_walls = true;

        // Back zone (z_back): lightweight constructions (Case 600)
        // These match by substring on paths containing "z_back"
        let (lt_wall, lt_roof, lt_floor, window) = bestest_600_constructions();
        c.constructions
            .insert("z_back/back/wall".to_string(), lt_wall);
        c.constructions
            .insert("z_back/back/floor".to_string(), lt_floor.clone());

        // Sunspace (z_sun): heavyweight constructions (Case 900)
        let (hw_wall, _hw_roof, hw_floor, _hw_window) = bestest_900_constructions();
        c.constructions
            .insert("z_sun/sunspace/wall_0".to_string(), hw_wall.clone());
        c.constructions
            .insert("z_sun/sunspace/wall_1".to_string(), hw_wall.clone());
        c.constructions
            .insert("z_sun/sunspace/wall_3".to_string(), hw_wall);
        c.constructions
            .insert("z_sun/sunspace/floor".to_string(), hw_floor);

        // Both zones: lightweight roof
        c.constructions.insert("ceiling".to_string(), lt_roof);
        // Windows (only on sunspace south wall)
        c.constructions.insert("window".to_string(), window);

        // Common wall (partition): concrete block
        // Thickness: 0.200 m, k = 0.510 W/(mK), c = 1000 J/(kgK), rho = 1400 kg/m³
        let partition = WallConstruction::new(
            "COMMON_WALL",
            vec![Layer {
                name: "CONCRETE BLOCK".to_string(),
                thickness: 0.200,
                conductivity: 0.510,
                density: 1400.0,
                specific_heat: 1000.0,
            }],
        );
        c.constructions.insert("common_wall".to_string(), partition);

        // Ground-coupled floors
        c.ground_temperature_c = Some(10.0);

        c
    };

    for spec in &suite {
        // Select geometry.
        let building = match spec.geometry_kind {
            South => &building_south,
            EastWest => &building_ew,
            Sunspace => &building_960,
        };
        let base_cfg = match spec.geometry_kind {
            South => &base_cfg_south,
            EastWest => &base_cfg_ew,
            Sunspace => &base_cfg_960,
        };

        // Apply heavyweight constructions (skip for Sunspace; base_cfg_960 is pre-configured).
        let mut cfg = base_cfg.clone();
        if spec.mass_kind == Heavy && spec.geometry_kind != Sunspace {
            let (wall, roof, floor, window) = bestest_900_constructions();
            cfg.constructions.insert("window".to_string(), window);
            cfg.constructions.insert("ceiling".to_string(), roof);
            cfg.constructions.insert("floor".to_string(), floor);
            cfg.constructions.insert("wall".to_string(), wall);
        }
        cfg.use_view_factor_radiation = std::env::var("BESTEST_GLOBAL_VF")
            .ok()
            .as_deref()
            .map(|s| s == "1" || s.eq_ignore_ascii_case("true"))
            .unwrap_or(false);
        cfg.view_factor_rays_per_surface = 10_000;
        cfg.interior_emissivity = 0.9;
        cfg.thermal_capacity_j_per_m3_k = estimate_zone_capacity_j_per_m3_k(building, &cfg);

        // Apply shading to solar config.
        let mut case_solar_cfg = solar_cfg.clone();
        match spec.shading {
            SK::Overhang => {
                // South-facing window overhang (cases 610/910).
                let sh = south_overhang_shading(2.0, 3.0);
                case_solar_cfg
                    .window_shading
                    .insert("window".to_string(), sh);
            }
            SK::OverhangFins => {
                // East/west window overhang + fins (cases 630/930).
                let sh_east = east_overhang_fins_shading(2.0, 3.0);
                let sh_west = west_overhang_fins_shading(2.0, 3.0);
                case_solar_cfg
                    .window_shading
                    .insert("window_1".to_string(), sh_east);
                case_solar_cfg
                    .window_shading
                    .insert("window_2".to_string(), sh_west);
            }
            SK::None => {}
        }

        // Build HVAC schedules.
        let volume: f64 = building.zones().iter().map(|z| z.volume()).sum();
        let mut hourly_heating_setpoint: Option<Vec<f64>> = Option::None;
        let mut hourly_cooling_setpoint: Option<Vec<f64>> = Option::None;
        let mut hourly_infiltration_ach: Option<Vec<f64>> = Option::None;

        match spec.hvac_mode {
            Normal => {}
            Setback => {
                hourly_heating_setpoint = Some(make_setback_heating_schedule());
            }
            NightVent => {
                // No heating, daytime cooling only, night ventilation.
                hourly_heating_setpoint = Some(make_free_float_heating_schedule());
                hourly_cooling_setpoint = Some(make_daytime_cooling_schedule());
                hourly_infiltration_ach = Some(make_night_vent_infiltration_schedule(volume));
            }
            FreeFloat => {
                hourly_heating_setpoint = Some(make_free_float_heating_schedule());
                hourly_cooling_setpoint = Some(make_free_float_cooling_schedule());
            }
            FreeFloatNightVent => {
                hourly_heating_setpoint = Some(make_free_float_heating_schedule());
                hourly_cooling_setpoint = Some(make_free_float_cooling_schedule());
                hourly_infiltration_ach = Some(make_night_vent_infiltration_schedule(volume));
            }
        }

        // Per-zone HVAC for Case 960: z_back (idx 0) conditioned, z_sun (idx 1) free-float.
        let per_zone_hvac = if spec.geometry_kind == Sunspace {
            Some(vec![
                HvacIdealLoads::with_setpoints(20.0, 27.0), // z_back: conditioned
                HvacIdealLoads::with_setpoints(-999.0, 999.0), // z_sun: free-float
            ])
        } else {
            None
        };

        let options = TransientSimulationOptions {
            warmup_hours: warmup_days.saturating_mul(24),
            substeps_per_hour,
            hourly_heating_setpoint,
            hourly_cooling_setpoint,
            hourly_infiltration_ach,
            per_zone_hvac,
        };

        let solar = if spec.solar {
            Some(&case_solar_cfg)
        } else {
            Option::None
        };
        let annual = run_transient_simulation_with_options(
            building,
            &cfg,
            &weather,
            &hvac,
            Option::None,
            solar,
            &options,
        );

        let is_freefloat = matches!(spec.hvac_mode, FreeFloat | FreeFloatNightVent);

        if is_freefloat {
            // Free-float output: min/max zone temperature.
            let sim_min = annual.min_zone_temp_c;
            let sim_max = annual.max_zone_temp_c;
            print!(
                "Case {:>12}: Tmin={:7.2} C, Tmax={:7.2} C",
                spec.name, sim_min, sim_max
            );
            if let Some(ff_ref) = freefloat_ref(spec.name) {
                print!(
                    " | vs ref: Tmin {:+6.2}C, Tmax {:+6.2}C",
                    sim_min - ff_ref.min_temp_c,
                    sim_max - ff_ref.max_temp_c
                );
            }
            println!();
        } else {
            // HVAC output: heating/cooling.
            let sim_h = annual.annual_heating_kwh;
            let sim_c = annual.annual_cooling_kwh;

            print!(
                "Case {:>12}: heating={:8.1} kWh, cooling={:8.1} kWh",
                spec.name, sim_h, sim_c
            );
            if let Some(a_ref) = annual_ref(spec.name) {
                let h_pct = if a_ref.heating_kwh > 0.0 {
                    100.0 * (sim_h - a_ref.heating_kwh) / a_ref.heating_kwh
                } else {
                    0.0
                };
                let c_pct = if a_ref.cooling_kwh > 0.0 {
                    100.0 * (sim_c - a_ref.cooling_kwh) / a_ref.cooling_kwh
                } else {
                    0.0
                };
                print!(
                    " | vs ref: heating {:+6.1}%, cooling {:+6.1}%",
                    h_pct, c_pct
                );
            }
            println!();

            // For Case 960, also show sunspace free-float temperature range.
            if spec.geometry_kind == Sunspace && annual.num_zones >= 2 {
                let sun_min = annual.per_zone_min_temp_c.get(1).copied().unwrap_or(0.0);
                let sun_max = annual.per_zone_max_temp_c.get(1).copied().unwrap_or(0.0);
                println!(
                    "             sunspace: Tmin={:7.2} C, Tmax={:7.2} C",
                    sun_min, sun_max
                );
            }

            if spec.solar {
                if let Ok(diag) =
                    compute_case_diagnostics(building, &cfg, &case_solar_cfg, &weather, &annual)
                {
                    diag_rows.push((spec.clone(), annual.clone(), diag));
                }
            }
        }

        outputs.push((spec.clone(), annual));
    }

    // ── Summary table ──
    println!();
    println!("=== BESTEST ASHRAE 140 Summary ===");
    println!();
    println!("─── HVAC Cases ───");
    println!(
        "{:>8} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}",
        "Case",
        "Htg kWh",
        "Clg kWh",
        "Pk Htg W",
        "Pk Clg W",
        "Htg err%",
        "Clg err%",
        "PkH err%",
        "PkC err%"
    );
    let mut total_abs_dev = 0.0_f64;
    let mut num_dev_terms = 0_u32;

    for (spec, annual) in &outputs {
        let is_ff = matches!(spec.hvac_mode, FreeFloat | FreeFloatNightVent);
        if is_ff {
            continue;
        }
        if let Some(a_ref) = annual_ref(spec.name) {
            let h_pct = if a_ref.heating_kwh > 0.0 {
                let v = 100.0 * (annual.annual_heating_kwh - a_ref.heating_kwh) / a_ref.heating_kwh;
                total_abs_dev += v.abs();
                num_dev_terms += 1;
                format!("{:+.1}", v)
            } else {
                "-".to_string()
            };
            let c_pct = if a_ref.cooling_kwh > 0.0 {
                let v = 100.0 * (annual.annual_cooling_kwh - a_ref.cooling_kwh) / a_ref.cooling_kwh;
                total_abs_dev += v.abs();
                num_dev_terms += 1;
                format!("{:+.1}", v)
            } else {
                "-".to_string()
            };
            let pk_h_pct = if a_ref.peak_heating_w > 0.0 {
                let v = 100.0 * (annual.peak_heating - a_ref.peak_heating_w) / a_ref.peak_heating_w;
                format!("{:+.1}", v)
            } else {
                "-".to_string()
            };
            let pk_c_pct = if a_ref.peak_cooling_w > 0.0 {
                let v = 100.0 * (annual.peak_cooling - a_ref.peak_cooling_w) / a_ref.peak_cooling_w;
                format!("{:+.1}", v)
            } else {
                "-".to_string()
            };
            println!(
                "{:>8} {:>10.1} {:>10.1} {:>10.0} {:>10.0} {:>10} {:>10} {:>10} {:>10}",
                spec.name,
                annual.annual_heating_kwh,
                annual.annual_cooling_kwh,
                annual.peak_heating,
                annual.peak_cooling,
                h_pct,
                c_pct,
                pk_h_pct,
                pk_c_pct,
            );
        } else {
            println!(
                "{:>8} {:>10.1} {:>10.1} {:>10.0} {:>10.0} {:>10} {:>10} {:>10} {:>10}",
                spec.name,
                annual.annual_heating_kwh,
                annual.annual_cooling_kwh,
                annual.peak_heating,
                annual.peak_cooling,
                "-",
                "-",
                "-",
                "-",
            );
        }
    }

    println!();
    println!("─── Free-Float Cases ───");
    println!(
        "{:>8} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}",
        "Case", "Tmin C", "Tmax C", "ref Tmin", "ref Tmax", "dTmin C", "dTmax C"
    );
    for (spec, annual) in &outputs {
        let is_ff = matches!(spec.hvac_mode, FreeFloat | FreeFloatNightVent);
        if !is_ff {
            continue;
        }
        if let Some(ff_ref) = freefloat_ref(spec.name) {
            println!(
                "{:>8} {:>10.2} {:>10.2} {:>10.2} {:>10.2} {:>10} {:>10}",
                spec.name,
                annual.min_zone_temp_c,
                annual.max_zone_temp_c,
                ff_ref.min_temp_c,
                ff_ref.max_temp_c,
                format!("{:+.2}", annual.min_zone_temp_c - ff_ref.min_temp_c),
                format!("{:+.2}", annual.max_zone_temp_c - ff_ref.max_temp_c),
            );
        } else {
            println!(
                "{:>8} {:>10.2} {:>10.2} {:>10} {:>10} {:>10} {:>10}",
                spec.name, annual.min_zone_temp_c, annual.max_zone_temp_c, "-", "-", "-", "-",
            );
        }
    }

    println!();
    println!(
        "Total absolute deviation (HVAC annual loads): {:.1}pp across {} terms",
        total_abs_dev, num_dev_terms
    );

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
