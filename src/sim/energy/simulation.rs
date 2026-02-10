use crate::Building;

use crate::sim::index::SurfaceIndex;

use super::boundary::ThermalBoundaries;
use super::config::ThermalConfig;
use super::hvac::{HvacIdealLoads, LumpedThermalModel, TwoNodeThermalModel};
use super::network::{MultiZoneAirModel, ThermalNetwork};
use super::schedule::InternalGainsProfile;
use super::solar_bridge::{
    SolarGainConfig, SolarHourParams, compute_solar_gains_per_zone_with_materials,
    compute_solar_gains_with_materials,
};
use super::weather::WeatherData;
use super::zone::calculate_heat_balance_with_boundaries;
use crate::UID;
use crate::sim::lighting::solar::SolarPosition;

#[cfg(test)]
use super::zone::calculate_heat_balance;

/// Result of an annual energy simulation.
#[derive(Debug, Clone)]
pub struct AnnualResult {
    /// Hourly heating demand in W for each hour.
    pub hourly_heating: Vec<f64>,
    /// Hourly cooling demand in W for each hour.
    pub hourly_cooling: Vec<f64>,
    /// Total annual heating energy in kWh.
    pub annual_heating_kwh: f64,
    /// Total annual cooling energy in kWh.
    pub annual_cooling_kwh: f64,
    /// Peak heating demand in W.
    pub peak_heating: f64,
    /// Peak cooling demand in W.
    pub peak_cooling: f64,
    /// Monthly heating energy in kWh (12 values).
    pub monthly_heating_kwh: [f64; 12],
    /// Monthly cooling energy in kWh (12 values).
    pub monthly_cooling_kwh: [f64; 12],
}

/// Multi-zone transient simulation output (zone air node per `Zone`).
#[derive(Debug, Clone)]
pub struct MultiZoneAnnualResult {
    pub zone_uids: Vec<UID>,
    pub zone_names: Vec<String>,
    /// Temperatures per zone per hour, indexed as `[zone][hour]`.
    pub hourly_zone_temperatures_c: Vec<Vec<f64>>,
    /// Thermal HVAC heating power per zone per hour, indexed as `[zone][hour]`.
    pub hourly_zone_heating_w: Vec<Vec<f64>>,
    /// Thermal HVAC cooling power per zone per hour, indexed as `[zone][hour]`.
    pub hourly_zone_cooling_w: Vec<Vec<f64>>,
    /// Aggregated building-level totals (sums over zones).
    pub annual: AnnualResult,
}

/// Computes day of year from month and day.
fn day_of_year(month: u8, day: u8) -> u16 {
    const DAYS_BEFORE_MONTH: [u16; 12] = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];
    let m = (month as usize).saturating_sub(1).min(11);
    DAYS_BEFORE_MONTH[m] + day as u16
}

fn opaque_absorptance_for_path(
    path: &str,
    material_library: Option<&crate::sim::materials::MaterialLibrary>,
    default_absorptance: f64,
) -> f64 {
    let Some(lib) = material_library else {
        return default_absorptance;
    };
    let Some(mat) = lib.lookup(path) else {
        return default_absorptance;
    };
    let Some(opt) = mat.optical.as_ref() else {
        return default_absorptance;
    };

    let mut a = 0.0;
    for c in 0..3 {
        let absorb =
            (1.0 - opt.diffuse_reflectance[c] - opt.specular_reflectance[c] - opt.transmittance[c])
                .clamp(0.0, 1.0);
        a += absorb;
    }
    (a / 3.0).clamp(0.0, 1.0)
}

fn compute_exterior_opaque_sol_air_gain_total_w(
    building: &Building,
    base_config: &ThermalConfig,
    index: &SurfaceIndex,
    boundaries: &ThermalBoundaries,
    params: &SolarHourParams,
    solar_config: &SolarGainConfig,
) -> f64 {
    let solar_pos = SolarPosition::calculate(
        params.latitude,
        params.longitude,
        params.day_of_year,
        params.hour,
    );
    let sun_above = solar_pos.is_above_horizon();
    let sun_dir = solar_pos.to_direction();

    let h_out = solar_config
        .exterior_heat_transfer_coeff_w_per_m2_k
        .max(1e-9);

    let mut total = 0.0;
    for surface in &index.surfaces {
        if !boundaries.is_exterior(&surface.polygon_uid) {
            continue;
        }
        if solar_config
            .resolve_shgc_with_materials(&surface.path, base_config.material_library.as_ref())
            .is_some()
        {
            continue; // glazing handled via transmitted-to-zone path
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

        let absorptance = opaque_absorptance_for_path(
            &surface.path,
            base_config.material_library.as_ref(),
            solar_config.default_opaque_absorptance,
        );
        if absorptance <= 0.0 {
            continue;
        }

        let u = base_config.resolve_u_value_for_surface(&surface.polygon_uid, &surface.path);
        if !(u.is_finite() && u > 0.0) {
            continue;
        }

        // Sol-air style coupling: absorbed shortwave heats the exterior surface, but only
        // a fraction conducts inward instantaneously. Approximate that fraction as U/h_out.
        total += (u / h_out) * incident * absorptance * area;
    }
    total
}

fn compute_exterior_opaque_sol_air_gains_by_zone_w(
    building: &Building,
    base_config: &ThermalConfig,
    index: &SurfaceIndex,
    boundaries: &ThermalBoundaries,
    params: &SolarHourParams,
    solar_config: &SolarGainConfig,
) -> std::collections::HashMap<UID, f64> {
    let solar_pos = SolarPosition::calculate(
        params.latitude,
        params.longitude,
        params.day_of_year,
        params.hour,
    );
    let sun_above = solar_pos.is_above_horizon();
    let sun_dir = solar_pos.to_direction();

    let h_out = solar_config
        .exterior_heat_transfer_coeff_w_per_m2_k
        .max(1e-9);

    let mut gains: std::collections::HashMap<UID, f64> = std::collections::HashMap::new();
    for surface in &index.surfaces {
        if !boundaries.is_exterior(&surface.polygon_uid) {
            continue;
        }
        if solar_config
            .resolve_shgc_with_materials(&surface.path, base_config.material_library.as_ref())
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

        let absorptance = opaque_absorptance_for_path(
            &surface.path,
            base_config.material_library.as_ref(),
            solar_config.default_opaque_absorptance,
        );
        if absorptance <= 0.0 {
            continue;
        }

        let u = base_config.resolve_u_value_for_surface(&surface.polygon_uid, &surface.path);
        if !(u.is_finite() && u > 0.0) {
            continue;
        }

        let q = (u / h_out) * incident * absorptance * area;
        if q != 0.0 {
            *gains.entry(surface.zone_uid.clone()).or_insert(0.0) += q;
        }
    }
    gains
}

fn looks_like_glazing(
    path: &str,
    base_config: &ThermalConfig,
    solar_config: Option<&SolarGainConfig>,
) -> bool {
    if let Some(sc) = solar_config
        && sc
            .resolve_shgc_with_materials(path, base_config.material_library.as_ref())
            .is_some()
    {
        return true;
    }
    if let Some(lib) = base_config.material_library.as_ref()
        && let Some(mat) = lib.lookup(path)
        && mat.is_glazing
    {
        return true;
    }
    let p = path.to_ascii_lowercase();
    p.contains("window") || p.contains("glazing") || p.contains("glass")
}

fn estimate_opaque_envelope_area_m2(
    index: &SurfaceIndex,
    boundaries: &ThermalBoundaries,
    base_config: &ThermalConfig,
    solar_config: Option<&SolarGainConfig>,
) -> f64 {
    let mut a = 0.0;
    for s in &index.surfaces {
        if !boundaries.is_exterior(&s.polygon_uid) {
            continue;
        }
        if looks_like_glazing(&s.path, base_config, solar_config) {
            continue;
        }
        a += s.area_m2.max(0.0);
    }
    a.max(0.0)
}

/// Runs an annual hourly energy simulation.
///
/// For each hour:
///   1. Set outdoor temperature from weather data
///   2. Calculate internal gains from schedule
///   3. Compute solar gains from DNI/DHI + sun geometry + window SHGC
///   4. Run steady-state heat balance
///   5. Accumulate heating/cooling demand
///
/// If `solar_config` is `None`, solar gains are zero for all hours.
pub fn run_annual_simulation(
    building: &Building,
    base_config: &ThermalConfig,
    weather: &WeatherData,
    gains_profile: Option<&InternalGainsProfile>,
    solar_config: Option<&SolarGainConfig>,
) -> AnnualResult {
    let index = SurfaceIndex::new(building);
    let boundaries = ThermalBoundaries::classify(building, &index);

    let num_hours = weather.num_hours();
    let mut hourly_heating = Vec::with_capacity(num_hours);
    let mut hourly_cooling = Vec::with_capacity(num_hours);

    let mut annual_heating = 0.0;
    let mut annual_cooling = 0.0;
    let mut peak_heating = 0.0_f64;
    let mut peak_cooling = 0.0_f64;
    let mut monthly_heating = [0.0; 12];
    let mut monthly_cooling = [0.0; 12];

    for (hour_idx, record) in weather.records.iter().enumerate() {
        let mut config = base_config.clone();
        config.outdoor_temperature = record.dry_bulb_temperature;

        // Internal gains
        if let Some(profile) = gains_profile {
            config.internal_gains = profile.gains_at(hour_idx);
        }

        // Solar gains from weather data + sun geometry + glazing properties (+ optional opaque sol-air coupling)
        config.solar_gains = match solar_config {
            Some(sc) => {
                let params = SolarHourParams {
                    direct_normal_irradiance: record.direct_normal_radiation,
                    diffuse_horizontal_irradiance: record.diffuse_horizontal_radiation,
                    day_of_year: day_of_year(record.month, record.day),
                    hour: record.hour as f64,
                    latitude: weather.latitude,
                    longitude: weather.longitude,
                };
                let mut solar = compute_solar_gains_with_materials(
                    building,
                    &params,
                    sc,
                    config.material_library.as_ref(),
                );
                if sc.include_exterior_opaque_absorption {
                    solar += compute_exterior_opaque_sol_air_gain_total_w(
                        building,
                        &config,
                        &index,
                        &boundaries,
                        &params,
                        sc,
                    );
                }
                solar
            }
            None => 0.0,
        };

        let result = calculate_heat_balance_with_boundaries(building, &config, &boundaries);

        hourly_heating.push(result.heating_demand);
        hourly_cooling.push(result.cooling_demand);

        annual_heating += result.heating_demand;
        annual_cooling += result.cooling_demand;
        peak_heating = peak_heating.max(result.heating_demand);
        peak_cooling = peak_cooling.max(result.cooling_demand);

        let month_idx = (record.month as usize).saturating_sub(1).min(11);
        monthly_heating[month_idx] += result.heating_demand;
        monthly_cooling[month_idx] += result.cooling_demand;
    }

    // Convert W*h to kWh
    let to_kwh = 1.0 / 1000.0;

    AnnualResult {
        hourly_heating,
        hourly_cooling,
        annual_heating_kwh: annual_heating * to_kwh,
        annual_cooling_kwh: annual_cooling * to_kwh,
        peak_heating,
        peak_cooling,
        monthly_heating_kwh: monthly_heating.map(|v| v * to_kwh),
        monthly_cooling_kwh: monthly_cooling.map(|v| v * to_kwh),
    }
}

/// Runs a transient annual simulation with thermal mass and HVAC.
///
/// Unlike `run_annual_simulation` which is steady-state, this models
/// the zone temperature state between time steps using a 1R1C model.
///
/// If `solar_config` is `None`, solar gains are zero for all hours.
pub fn run_transient_simulation(
    building: &Building,
    base_config: &ThermalConfig,
    weather: &WeatherData,
    hvac: &HvacIdealLoads,
    gains_profile: Option<&InternalGainsProfile>,
    solar_config: Option<&SolarGainConfig>,
) -> AnnualResult {
    let index = SurfaceIndex::new(building);
    let boundaries = ThermalBoundaries::classify(building, &index);

    let num_hours = weather.num_hours();
    let mut hourly_heating = Vec::with_capacity(num_hours);
    let mut hourly_cooling = Vec::with_capacity(num_hours);

    // Compute building-level UA and thermal capacity
    let steady = calculate_heat_balance_with_boundaries(building, base_config, &boundaries);
    let dt = base_config.indoor_temperature - base_config.outdoor_temperature;
    let ua_total = if dt.abs() > 1e-10 {
        steady.transmission_loss / dt
    } else {
        // Fall back to calculating from default U-value and areas
        let mut ua = 0.0;
        for zone in building.zones() {
            for solid in zone.solids() {
                for wall in solid.walls() {
                    for polygon in wall.polygons() {
                        if !boundaries.is_exterior(&polygon.uid) {
                            continue;
                        }
                        let path = format!(
                            "{}/{}/{}/{}",
                            zone.name, solid.name, wall.name, polygon.name
                        );
                        let u = base_config.resolve_u_value_for_surface(&polygon.uid, &path);
                        ua += u * polygon.area();
                    }
                }
            }
        }
        ua
    };

    let infiltration_cond = if dt.abs() > 1e-10 {
        steady.infiltration_loss / dt
    } else {
        // rho * cp * V * ACH / 3600
        let volume: f64 = building.zones().iter().map(|z| z.volume()).sum();
        1.2 * 1005.0 * volume * base_config.infiltration_ach / 3600.0
    };

    // Estimate thermal capacity from building volume.
    // The factor 50 kJ/(m^3*K) is a tuning parameter for medium-weight construction.
    // Typical range: ~30 kJ/(m^3*K) (lightweight) to ~80 kJ/(m^3*K) (heavyweight).
    // For more accurate results, derive from actual construction layer properties.
    let volume: f64 = building.zones().iter().map(|z| z.volume()).sum();
    let thermal_capacity = volume * base_config.thermal_capacity_j_per_m3_k; // J/K

    let k_env = ua_total + infiltration_cond;
    let dt_s = 3600.0;

    let use_two_node = base_config.two_node_mass_fraction > 0.0
        && base_config.interior_heat_transfer_coeff_w_per_m2_k > 0.0
        && thermal_capacity > 0.0
        && k_env.is_finite()
        && k_env >= 0.0;

    let opaque_area_m2 = if use_two_node {
        estimate_opaque_envelope_area_m2(&index, &boundaries, base_config, solar_config)
    } else {
        0.0
    };

    let (mut model_1r1c, mut model_2r2c) = if use_two_node && opaque_area_m2 > 0.0 {
        let f_mass = base_config.two_node_mass_fraction.clamp(0.0, 1.0);
        let c_total = thermal_capacity;
        let air_capacity_min = 1.2 * 1005.0 * volume; // rho*cp*V
        let mut c_air = (1.0 - f_mass) * c_total;
        c_air = c_air.max(air_capacity_min).min(c_total);
        let c_mass = (c_total - c_air).max(0.0);

        let k_am = base_config.interior_heat_transfer_coeff_w_per_m2_k.max(0.0) * opaque_area_m2;
        (
            None,
            Some(TwoNodeThermalModel::new(
                base_config.indoor_temperature,
                base_config.indoor_temperature,
                k_env,
                k_am,
                c_air,
                c_mass.max(1.0),
            )),
        )
    } else {
        (
            Some(LumpedThermalModel::new(
                base_config.indoor_temperature,
                ua_total,
                infiltration_cond,
                thermal_capacity,
            )),
            None,
        )
    };

    let mut annual_heating = 0.0;
    let mut annual_cooling = 0.0;
    let mut peak_heating = 0.0_f64;
    let mut peak_cooling = 0.0_f64;
    let mut monthly_heating = [0.0; 12];
    let mut monthly_cooling = [0.0; 12];

    for (hour_idx, record) in weather.records.iter().enumerate() {
        let gains = gains_profile.map(|p| p.gains_at(hour_idx)).unwrap_or(0.0);
        let (solar_transmitted, solar_opaque_sol_air) = match solar_config {
            Some(sc) => {
                let params = SolarHourParams {
                    direct_normal_irradiance: record.direct_normal_radiation,
                    diffuse_horizontal_irradiance: record.diffuse_horizontal_radiation,
                    day_of_year: day_of_year(record.month, record.day),
                    hour: record.hour as f64,
                    latitude: weather.latitude,
                    longitude: weather.longitude,
                };
                let transmitted = compute_solar_gains_with_materials(
                    building,
                    &params,
                    sc,
                    base_config.material_library.as_ref(),
                );
                let opaque = if sc.include_exterior_opaque_absorption {
                    compute_exterior_opaque_sol_air_gain_total_w(
                        building,
                        base_config,
                        &index,
                        &boundaries,
                        &params,
                        sc,
                    )
                } else {
                    0.0
                };
                (transmitted, opaque)
            }
            None => (0.0, 0.0),
        };
        let solar_total = solar_transmitted + solar_opaque_sol_air;

        let (heating_power, cooling_power) = if let Some(model) = model_2r2c.as_mut() {
            let f_solar_mass = base_config.solar_gains_to_mass_fraction.clamp(0.0, 1.0);
            let f_internal_mass = base_config.internal_gains_to_mass_fraction.clamp(0.0, 1.0);

            let gains_mass = gains * f_internal_mass + solar_transmitted * f_solar_mass;
            let gains_air = gains * (1.0 - f_internal_mass)
                + solar_opaque_sol_air
                + solar_transmitted * (1.0 - f_solar_mass);

            let t_air = model.air_temperature_c;
            let setpoint = hvac.active_setpoint(t_air);
            let need_hvac = (setpoint - t_air).abs() > 1e-10;
            let q_hvac = if need_hvac {
                hvac.required_hvac_power_two_node(
                    setpoint,
                    t_air,
                    model.mass_temperature_c,
                    record.dry_bulb_temperature,
                    model.envelope_conductance_w_per_k,
                    model.air_mass_conductance_w_per_k,
                    gains_air,
                    gains_mass,
                    model.air_capacity_j_per_k,
                    model.mass_capacity_j_per_k,
                    dt_s,
                )
            } else {
                0.0
            };

            model.step(
                record.dry_bulb_temperature,
                gains_air,
                gains_mass,
                q_hvac,
                dt_s,
            );

            (q_hvac.max(0.0), (-q_hvac).max(0.0))
        } else {
            let model = model_1r1c.as_mut().unwrap();
            let total_gains = gains + solar_total;

            // Calculate HVAC power using implicit formula that accounts for
            // concurrent envelope losses during the timestep.
            let total_conductance = model.ua_total + model.infiltration_conductance;
            let (heating_power, cooling_power) = hvac.calculate_with_losses(
                model.zone_temperature,
                record.dry_bulb_temperature,
                total_conductance,
                total_gains,
                model.thermal_capacity,
                dt_s,
            );

            let hvac_net = heating_power - cooling_power;
            model.step(record.dry_bulb_temperature, total_gains, hvac_net, dt_s);
            (heating_power, cooling_power)
        };

        hourly_heating.push(heating_power);
        hourly_cooling.push(cooling_power);

        annual_heating += heating_power;
        annual_cooling += cooling_power;
        peak_heating = peak_heating.max(heating_power);
        peak_cooling = peak_cooling.max(cooling_power);

        let month_idx = (record.month as usize).saturating_sub(1).min(11);
        monthly_heating[month_idx] += heating_power;
        monthly_cooling[month_idx] += cooling_power;
    }

    let to_kwh = 1.0 / 1000.0;

    AnnualResult {
        hourly_heating,
        hourly_cooling,
        annual_heating_kwh: annual_heating * to_kwh,
        annual_cooling_kwh: annual_cooling * to_kwh,
        peak_heating,
        peak_cooling,
        monthly_heating_kwh: monthly_heating.map(|v| v * to_kwh),
        monthly_cooling_kwh: monthly_cooling.map(|v| v * to_kwh),
    }
}

/// Runs a multi-zone transient annual simulation with zone air nodes and ideal HVAC.
///
/// Differences vs [`run_transient_simulation`]:
/// - one air temperature state per `Zone` (coupled through inter-zone partitions),
/// - exterior transmission excludes internal interfaces by construction,
/// - HVAC is applied per-zone using ideal setpoints (unlimited capacity).
///
/// Gains are distributed across zones proportional to zone volume (placeholder policy
/// until per-zone schedules and solar distributions are provided by upstream modules).
pub fn run_multizone_transient_simulation(
    building: &Building,
    base_config: &ThermalConfig,
    weather: &WeatherData,
    hvac: &HvacIdealLoads,
    gains_profile: Option<&InternalGainsProfile>,
    solar_config: Option<&SolarGainConfig>,
) -> anyhow::Result<MultiZoneAnnualResult> {
    let index = SurfaceIndex::new(building);
    let boundaries = ThermalBoundaries::classify(building, &index);
    let network = ThermalNetwork::build(building, base_config, &index, &boundaries);

    let mut model = MultiZoneAirModel::new(
        building,
        &network,
        base_config.infiltration_ach,
        base_config.thermal_capacity_j_per_m3_k,
        base_config.indoor_temperature,
    );

    let zones = building.zones();
    let zone_volumes_m3: Vec<f64> = zones.iter().map(|z| z.volume()).collect();
    let total_volume: f64 = zone_volumes_m3.iter().sum();

    let num_hours = weather.num_hours();
    let n_zones = zones.len();

    let mut hourly_zone_temperatures_c = vec![Vec::with_capacity(num_hours); n_zones];
    let mut hourly_zone_heating_w = vec![Vec::with_capacity(num_hours); n_zones];
    let mut hourly_zone_cooling_w = vec![Vec::with_capacity(num_hours); n_zones];

    let mut hourly_heating = Vec::with_capacity(num_hours);
    let mut hourly_cooling = Vec::with_capacity(num_hours);

    let mut annual_heating = 0.0;
    let mut annual_cooling = 0.0;
    let mut peak_heating = 0.0_f64;
    let mut peak_cooling = 0.0_f64;
    let mut monthly_heating = [0.0; 12];
    let mut monthly_cooling = [0.0; 12];

    for (hour_idx, record) in weather.records.iter().enumerate() {
        let gains_internal = gains_profile.map(|p| p.gains_at(hour_idx)).unwrap_or(0.0);
        let solar_by_zone = match solar_config {
            Some(sc) => {
                let params = SolarHourParams {
                    direct_normal_irradiance: record.direct_normal_radiation,
                    diffuse_horizontal_irradiance: record.diffuse_horizontal_radiation,
                    day_of_year: day_of_year(record.month, record.day),
                    hour: record.hour as f64,
                    latitude: weather.latitude,
                    longitude: weather.longitude,
                };
                let mut solar = compute_solar_gains_per_zone_with_materials(
                    building,
                    &params,
                    sc,
                    base_config.material_library.as_ref(),
                );
                if sc.include_exterior_opaque_absorption {
                    let opaque = compute_exterior_opaque_sol_air_gains_by_zone_w(
                        building,
                        base_config,
                        &index,
                        &boundaries,
                        &params,
                        sc,
                    );
                    for (z, q) in opaque {
                        *solar.entry(z).or_insert(0.0) += q;
                    }
                }
                solar
            }
            None => std::collections::HashMap::new(),
        };
        let mut gains_by_zone = vec![0.0; n_zones];
        if total_volume > 1e-14 {
            for i in 0..n_zones {
                let internal_i = gains_internal * (zone_volumes_m3[i] / total_volume);
                let solar_i = solar_by_zone
                    .get(&model.zone_uids()[i])
                    .cloned()
                    .unwrap_or(0.0);
                gains_by_zone[i] = internal_i + solar_i;
            }
        }

        let step = model.step(record.dry_bulb_temperature, &gains_by_zone, hvac, 3600.0)?;

        let hour_heating: f64 = step.zone_heating_w.iter().sum();
        let hour_cooling: f64 = step.zone_cooling_w.iter().sum();

        hourly_heating.push(hour_heating);
        hourly_cooling.push(hour_cooling);

        annual_heating += hour_heating;
        annual_cooling += hour_cooling;
        peak_heating = peak_heating.max(hour_heating);
        peak_cooling = peak_cooling.max(hour_cooling);

        let month_idx = (record.month as usize).saturating_sub(1).min(11);
        monthly_heating[month_idx] += hour_heating;
        monthly_cooling[month_idx] += hour_cooling;

        for z in 0..n_zones {
            hourly_zone_temperatures_c[z].push(step.zone_temperatures_c[z]);
            hourly_zone_heating_w[z].push(step.zone_heating_w[z]);
            hourly_zone_cooling_w[z].push(step.zone_cooling_w[z]);
        }
    }

    let to_kwh = 1.0 / 1000.0;
    let annual = AnnualResult {
        hourly_heating,
        hourly_cooling,
        annual_heating_kwh: annual_heating * to_kwh,
        annual_cooling_kwh: annual_cooling * to_kwh,
        peak_heating,
        peak_cooling,
        monthly_heating_kwh: monthly_heating.map(|v| v * to_kwh),
        monthly_cooling_kwh: monthly_cooling.map(|v| v * to_kwh),
    };

    Ok(MultiZoneAnnualResult {
        zone_uids: model.zone_uids().to_vec(),
        zone_names: model.zone_names().to_vec(),
        hourly_zone_temperatures_c,
        hourly_zone_heating_w,
        hourly_zone_cooling_w,
        annual,
    })
}

/// Runs a multi-zone *steady-state* annual simulation with zone air nodes and ideal HVAC.
///
/// This uses the same multi-zone coupling model as [`run_multizone_transient_simulation`],
/// but with zero thermal capacity (instantaneous steady-state each hour). It is useful for
/// quick annual load estimates and for debugging inter-zone coupling without dynamics.
pub fn run_multizone_steady_simulation(
    building: &Building,
    base_config: &ThermalConfig,
    weather: &WeatherData,
    hvac: &HvacIdealLoads,
    gains_profile: Option<&InternalGainsProfile>,
    solar_config: Option<&SolarGainConfig>,
) -> anyhow::Result<MultiZoneAnnualResult> {
    let index = SurfaceIndex::new(building);
    let boundaries = ThermalBoundaries::classify(building, &index);
    let network = ThermalNetwork::build(building, base_config, &index, &boundaries);

    // Zero-capacity zones: steady-state per hour.
    let mut model = MultiZoneAirModel::new(
        building,
        &network,
        base_config.infiltration_ach,
        0.0,
        base_config.indoor_temperature,
    );

    let zones = building.zones();
    let zone_volumes_m3: Vec<f64> = zones.iter().map(|z| z.volume()).collect();
    let total_volume: f64 = zone_volumes_m3.iter().sum();

    let num_hours = weather.num_hours();
    let n_zones = zones.len();

    let mut hourly_zone_temperatures_c = vec![Vec::with_capacity(num_hours); n_zones];
    let mut hourly_zone_heating_w = vec![Vec::with_capacity(num_hours); n_zones];
    let mut hourly_zone_cooling_w = vec![Vec::with_capacity(num_hours); n_zones];

    let mut hourly_heating = Vec::with_capacity(num_hours);
    let mut hourly_cooling = Vec::with_capacity(num_hours);

    let mut annual_heating = 0.0;
    let mut annual_cooling = 0.0;
    let mut peak_heating = 0.0_f64;
    let mut peak_cooling = 0.0_f64;
    let mut monthly_heating = [0.0; 12];
    let mut monthly_cooling = [0.0; 12];

    for (hour_idx, record) in weather.records.iter().enumerate() {
        let gains_internal = gains_profile.map(|p| p.gains_at(hour_idx)).unwrap_or(0.0);
        let solar_by_zone = match solar_config {
            Some(sc) => {
                let params = SolarHourParams {
                    direct_normal_irradiance: record.direct_normal_radiation,
                    diffuse_horizontal_irradiance: record.diffuse_horizontal_radiation,
                    day_of_year: day_of_year(record.month, record.day),
                    hour: record.hour as f64,
                    latitude: weather.latitude,
                    longitude: weather.longitude,
                };
                compute_solar_gains_per_zone_with_materials(
                    building,
                    &params,
                    sc,
                    base_config.material_library.as_ref(),
                )
            }
            None => std::collections::HashMap::new(),
        };

        let mut gains_by_zone = vec![0.0; n_zones];
        if total_volume > 1e-14 {
            for i in 0..n_zones {
                let internal_i = gains_internal * (zone_volumes_m3[i] / total_volume);
                let solar_i = solar_by_zone
                    .get(&model.zone_uids()[i])
                    .cloned()
                    .unwrap_or(0.0);
                gains_by_zone[i] = internal_i + solar_i;
            }
        }

        let step = model.step(record.dry_bulb_temperature, &gains_by_zone, hvac, 3600.0)?;

        let hour_heating: f64 = step.zone_heating_w.iter().sum();
        let hour_cooling: f64 = step.zone_cooling_w.iter().sum();

        hourly_heating.push(hour_heating);
        hourly_cooling.push(hour_cooling);

        annual_heating += hour_heating;
        annual_cooling += hour_cooling;
        peak_heating = peak_heating.max(hour_heating);
        peak_cooling = peak_cooling.max(hour_cooling);

        let month_idx = (record.month as usize).saturating_sub(1).min(11);
        monthly_heating[month_idx] += hour_heating;
        monthly_cooling[month_idx] += hour_cooling;

        for z in 0..n_zones {
            hourly_zone_temperatures_c[z].push(step.zone_temperatures_c[z]);
            hourly_zone_heating_w[z].push(step.zone_heating_w[z]);
            hourly_zone_cooling_w[z].push(step.zone_cooling_w[z]);
        }
    }

    let to_kwh = 1.0 / 1000.0;
    let annual = AnnualResult {
        hourly_heating,
        hourly_cooling,
        annual_heating_kwh: annual_heating * to_kwh,
        annual_cooling_kwh: annual_cooling * to_kwh,
        peak_heating,
        peak_cooling,
        monthly_heating_kwh: monthly_heating.map(|v| v * to_kwh),
        monthly_cooling_kwh: monthly_cooling.map(|v| v * to_kwh),
    };

    Ok(MultiZoneAnnualResult {
        zone_uids: model.zone_uids().to_vec(),
        zone_names: model.zone_names().to_vec(),
        hourly_zone_temperatures_c,
        hourly_zone_heating_w,
        hourly_zone_cooling_w,
        annual,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::energy::weather::HourlyRecord;
    use crate::{Solid, Zone};

    #[test]
    fn test_annual_simulation_basic() {
        let s = Solid::from_box(5.0, 5.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let config = ThermalConfig::new();
        let weather = WeatherData::synthetic("Test", 52.0, 13.0, 10.0, 12.0);

        let result = run_annual_simulation(&building, &config, &weather, None, None);

        assert_eq!(result.hourly_heating.len(), 8760);
        assert_eq!(result.hourly_cooling.len(), 8760);
        assert!(
            result.annual_heating_kwh > 0.0,
            "Should need heating in cold climate"
        );
        assert!(result.peak_heating > 0.0, "Should have peak heating load");
    }

    #[test]
    fn test_annual_with_gains() {
        let s = Solid::from_box(10.0, 10.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let config = ThermalConfig::new();
        let weather = WeatherData::synthetic("Test", 52.0, 13.0, 10.0, 12.0);
        let gains = InternalGainsProfile::office(100.0);

        let result_no_gains = run_annual_simulation(&building, &config, &weather, None, None);
        let result_gains = run_annual_simulation(&building, &config, &weather, Some(&gains), None);

        // Internal gains should reduce heating demand
        assert!(
            result_gains.annual_heating_kwh < result_no_gains.annual_heating_kwh,
            "Internal gains should reduce heating: {} vs {}",
            result_gains.annual_heating_kwh,
            result_no_gains.annual_heating_kwh
        );
    }

    #[test]
    fn test_monthly_sums() {
        let s = Solid::from_box(3.0, 3.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let config = ThermalConfig::new();
        let weather = WeatherData::synthetic("Test", 52.0, 13.0, 10.0, 12.0);

        let result = run_annual_simulation(&building, &config, &weather, None, None);

        let monthly_sum: f64 = result.monthly_heating_kwh.iter().sum();
        assert!(
            (monthly_sum - result.annual_heating_kwh).abs() < 1.0,
            "Monthly heating sum should equal annual: {monthly_sum} vs {}",
            result.annual_heating_kwh
        );
    }

    #[test]
    fn test_day_of_year_basic() {
        assert_eq!(day_of_year(1, 1), 1);
        assert_eq!(day_of_year(3, 1), 60); // non-leap-year convention (Jan31+Feb28+1)
        assert_eq!(day_of_year(12, 31), 365);
    }

    #[test]
    fn test_annual_simulation_with_solar_config_small_weather() {
        let s = Solid::from_box(5.0, 5.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let config = ThermalConfig::new();
        let weather = WeatherData {
            location: "Mini".to_string(),
            latitude: 0.0,
            longitude: 0.0,
            timezone: 0.0,
            elevation: 0.0,
            records: vec![
                HourlyRecord {
                    month: 3,
                    day: 20,
                    hour: 12,
                    dry_bulb_temperature: 0.0,
                    relative_humidity: 50.0,
                    global_horizontal_radiation: 0.0,
                    direct_normal_radiation: 800.0,
                    diffuse_horizontal_radiation: 100.0,
                    wind_speed: 0.0,
                    wind_direction: 0.0,
                },
                HourlyRecord {
                    month: 3,
                    day: 20,
                    hour: 1,
                    dry_bulb_temperature: 0.0,
                    relative_humidity: 50.0,
                    global_horizontal_radiation: 0.0,
                    direct_normal_radiation: 0.0,
                    diffuse_horizontal_radiation: 50.0,
                    wind_speed: 0.0,
                    wind_direction: 0.0,
                },
            ],
        };

        let mut solar = SolarGainConfig::new();
        solar.default_shgc = 1.0;
        solar.glazing_patterns = vec!["floor".to_string()];

        let no_solar = run_annual_simulation(&building, &config, &weather, None, None);
        let with_solar = run_annual_simulation(&building, &config, &weather, None, Some(&solar));

        assert_eq!(no_solar.hourly_heating.len(), 2);
        assert_eq!(with_solar.hourly_heating.len(), 2);
        assert!(
            with_solar.hourly_heating[0] <= no_solar.hourly_heating[0],
            "Solar gains should not increase heating demand"
        );
    }

    #[test]
    fn test_multizone_steady_simulation_hits_solar_per_zone_branch() {
        let s = Solid::from_box(5.0, 5.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config = ThermalConfig::new();
        config.default_u_value = 1.0;
        config.infiltration_ach = 0.0;
        config.indoor_temperature = 20.0;
        config.outdoor_temperature = 20.0;

        let weather = WeatherData {
            location: "Mini".to_string(),
            latitude: 0.0,
            longitude: 0.0,
            timezone: 0.0,
            elevation: 0.0,
            records: vec![HourlyRecord {
                month: 6,
                day: 1,
                hour: 12,
                dry_bulb_temperature: 20.0,
                relative_humidity: 50.0,
                global_horizontal_radiation: 0.0,
                direct_normal_radiation: 0.0,
                diffuse_horizontal_radiation: 200.0,
                wind_speed: 0.0,
                wind_direction: 0.0,
            }],
        };

        let hvac = HvacIdealLoads::new();

        let no_solar =
            run_multizone_steady_simulation(&building, &config, &weather, &hvac, None, None)
                .unwrap();

        let mut solar = SolarGainConfig::new();
        solar.default_shgc = 1.0;
        solar.glazing_patterns = vec!["floor".to_string()];

        let with_solar = run_multizone_steady_simulation(
            &building,
            &config,
            &weather,
            &hvac,
            None,
            Some(&solar),
        )
        .unwrap();

        assert_eq!(no_solar.annual.hourly_heating.len(), 1);
        assert_eq!(with_solar.annual.hourly_heating.len(), 1);
        assert_eq!(no_solar.annual.hourly_cooling.len(), 1);
        assert_eq!(with_solar.annual.hourly_cooling.len(), 1);
        assert!(
            with_solar.annual.hourly_cooling[0] > no_solar.annual.hourly_cooling[0],
            "Solar gains should increase cooling demand in steady-state"
        );
    }

    #[test]
    fn test_transient_simulation() {
        let s = Solid::from_box(5.0, 5.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let config = ThermalConfig::new();
        let weather = WeatherData::synthetic("Test", 52.0, 13.0, 10.0, 12.0);
        let hvac = HvacIdealLoads::new();

        let result = run_transient_simulation(&building, &config, &weather, &hvac, None, None);

        assert_eq!(result.hourly_heating.len(), 8760);
        assert!(
            result.annual_heating_kwh > 0.0,
            "Should need heating in cold climate"
        );
        assert!(result.peak_heating > 0.0, "Should have peak heating load");

        // Transient should have different results from steady-state due to thermal mass
        let steady = run_annual_simulation(&building, &config, &weather, None, None);
        // They won't be identical, but both should indicate heating is needed
        assert!(steady.annual_heating_kwh > 0.0);
    }

    #[test]
    fn test_transient_simulation_dt_zero_fallback_ua_and_solar() {
        let s = Solid::from_box(2.0, 2.0, 2.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config = ThermalConfig::new();
        config.indoor_temperature = 20.0;
        config.outdoor_temperature = 20.0; // dt=0 -> triggers UA/infiltration fallback branches

        let weather = WeatherData {
            location: "Mini".to_string(),
            latitude: 0.0,
            longitude: 0.0,
            timezone: 0.0,
            elevation: 0.0,
            records: vec![
                HourlyRecord {
                    month: 6,
                    day: 1,
                    hour: 12,
                    dry_bulb_temperature: 20.0,
                    relative_humidity: 50.0,
                    global_horizontal_radiation: 0.0,
                    direct_normal_radiation: 500.0,
                    diffuse_horizontal_radiation: 50.0,
                    wind_speed: 0.0,
                    wind_direction: 0.0,
                },
                HourlyRecord {
                    month: 6,
                    day: 1,
                    hour: 1,
                    dry_bulb_temperature: 20.0,
                    relative_humidity: 50.0,
                    global_horizontal_radiation: 0.0,
                    direct_normal_radiation: 0.0,
                    diffuse_horizontal_radiation: 10.0,
                    wind_speed: 0.0,
                    wind_direction: 0.0,
                },
            ],
        };

        let hvac = HvacIdealLoads::new();
        let gains = InternalGainsProfile::office(100.0);
        let mut solar = SolarGainConfig::new();
        solar.default_shgc = 1.0;
        solar.glazing_patterns = vec!["floor".to_string()];

        let result = run_transient_simulation(
            &building,
            &config,
            &weather,
            &hvac,
            Some(&gains),
            Some(&solar),
        );
        assert_eq!(result.hourly_heating.len(), 2);
        assert_eq!(result.hourly_cooling.len(), 2);
    }

    #[test]
    fn test_multizone_transient_simulation_basic() {
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "s0").unwrap();
        let s1 = Solid::from_box(2.0, 2.0, 2.0, Some((2.0, 0.0, 0.0)), "s1").unwrap();
        let z0 = Zone::new("z0", vec![s0]).unwrap();
        let z1 = Zone::new("z1", vec![s1]).unwrap();
        let building = Building::new("b", vec![z0, z1]).unwrap();

        let config = ThermalConfig::new();
        let weather = WeatherData::synthetic("Test", 52.0, 13.0, 10.0, 12.0);
        let hvac = HvacIdealLoads::new();

        let result =
            run_multizone_transient_simulation(&building, &config, &weather, &hvac, None, None)
                .unwrap();

        assert_eq!(result.zone_names.len(), 2);
        assert_eq!(result.hourly_zone_temperatures_c.len(), 2);
        assert_eq!(result.hourly_zone_temperatures_c[0].len(), 8760);
        assert_eq!(result.annual.hourly_heating.len(), 8760);
    }

    #[test]
    fn test_multizone_transient_with_solar_distribution() {
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "s0").unwrap();
        let s1 = Solid::from_box(2.0, 2.0, 2.0, Some((2.0, 0.0, 0.0)), "s1").unwrap();
        let z0 = Zone::new("z0", vec![s0]).unwrap();
        let z1 = Zone::new("z1", vec![s1]).unwrap();
        let building = Building::new("b", vec![z0, z1]).unwrap();

        let config = ThermalConfig::new();
        let weather = WeatherData {
            location: "Mini".to_string(),
            latitude: 0.0,
            longitude: 0.0,
            timezone: 0.0,
            elevation: 0.0,
            records: vec![
                HourlyRecord {
                    month: 6,
                    day: 1,
                    hour: 12,
                    dry_bulb_temperature: 0.0,
                    relative_humidity: 50.0,
                    global_horizontal_radiation: 0.0,
                    direct_normal_radiation: 800.0,
                    diffuse_horizontal_radiation: 100.0,
                    wind_speed: 0.0,
                    wind_direction: 0.0,
                },
                HourlyRecord {
                    month: 6,
                    day: 1,
                    hour: 1,
                    dry_bulb_temperature: 0.0,
                    relative_humidity: 50.0,
                    global_horizontal_radiation: 0.0,
                    direct_normal_radiation: 0.0,
                    diffuse_horizontal_radiation: 20.0,
                    wind_speed: 0.0,
                    wind_direction: 0.0,
                },
            ],
        };
        let hvac = HvacIdealLoads::new();
        let gains = InternalGainsProfile::office(100.0);
        let mut solar = SolarGainConfig::new();
        solar.default_shgc = 1.0;
        solar.glazing_patterns = vec!["floor".to_string()];

        let result = run_multizone_transient_simulation(
            &building,
            &config,
            &weather,
            &hvac,
            Some(&gains),
            Some(&solar),
        )
        .unwrap();

        assert_eq!(result.zone_names.len(), 2);
        assert_eq!(result.hourly_zone_temperatures_c.len(), 2);
        assert_eq!(result.hourly_zone_temperatures_c[0].len(), 2);
        assert_eq!(result.annual.hourly_heating.len(), 2);
        assert_eq!(result.annual.hourly_cooling.len(), 2);
    }

    #[test]
    fn test_multizone_steady_simulation_basic() {
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "s0").unwrap();
        let s1 = Solid::from_box(2.0, 2.0, 2.0, Some((2.0, 0.0, 0.0)), "s1").unwrap();
        let z0 = Zone::new("z0", vec![s0]).unwrap();
        let z1 = Zone::new("z1", vec![s1]).unwrap();
        let building = Building::new("b", vec![z0, z1]).unwrap();

        let config = ThermalConfig::new();
        let weather = WeatherData::synthetic("Test", 52.0, 13.0, 10.0, 12.0);
        let hvac = HvacIdealLoads::new();

        let result =
            run_multizone_steady_simulation(&building, &config, &weather, &hvac, None, None)
                .unwrap();

        assert_eq!(result.zone_names.len(), 2);
        assert_eq!(result.hourly_zone_temperatures_c.len(), 2);
        assert_eq!(result.hourly_zone_temperatures_c[0].len(), 8760);
        assert_eq!(result.annual.hourly_heating.len(), 8760);
    }

    // ── Physics verification tests ──────────────────────────────────────

    #[test]
    fn test_annual_energy_conservation() {
        // With a 1-hour timestep, summing hourly demands (W) over all hours is numerically
        // equivalent to the annual energy in Wh, i.e. annual_kwh * 1000.
        let s = Solid::from_box(4.0, 4.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let config = ThermalConfig::new();
        let weather = WeatherData::synthetic("Test", 50.0, 10.0, 8.0, 10.0);

        let result = run_annual_simulation(&building, &config, &weather, None, None);

        let sum_heating_wh: f64 = result.hourly_heating.iter().sum();
        let sum_cooling_wh: f64 = result.hourly_cooling.iter().sum();

        assert!(
            (sum_heating_wh / 1000.0 - result.annual_heating_kwh).abs() < 0.01,
            "Hourly heating sum={:.1} Wh vs annual={:.1} kWh",
            sum_heating_wh,
            result.annual_heating_kwh
        );
        assert!(
            (sum_cooling_wh / 1000.0 - result.annual_cooling_kwh).abs() < 0.01,
            "Hourly cooling sum={:.1} Wh vs annual={:.1} kWh",
            sum_cooling_wh,
            result.annual_cooling_kwh
        );
    }

    #[test]
    fn test_peak_is_max_of_hourly() {
        let s = Solid::from_box(5.0, 5.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let config = ThermalConfig::new();
        let weather = WeatherData::synthetic("Test", 52.0, 13.0, 10.0, 15.0);

        let result = run_annual_simulation(&building, &config, &weather, None, None);

        let max_heating = result
            .hourly_heating
            .iter()
            .cloned()
            .fold(0.0_f64, f64::max);
        let max_cooling = result
            .hourly_cooling
            .iter()
            .cloned()
            .fold(0.0_f64, f64::max);

        assert!(
            (result.peak_heating - max_heating).abs() < 1e-10,
            "peak_heating={} vs max(hourly)={}",
            result.peak_heating,
            max_heating
        );
        assert!(
            (result.peak_cooling - max_cooling).abs() < 1e-10,
            "peak_cooling={} vs max(hourly)={}",
            result.peak_cooling,
            max_cooling
        );
    }

    #[test]
    fn test_annual_hour_matches_steady_state_heat_balance() {
        let s = Solid::from_box(5.0, 5.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config = ThermalConfig::new();
        config.indoor_temperature = 20.0;

        // Synthetic weather with 0 annual amplitude (still has daily variation)
        let weather = WeatherData::synthetic("Const", 50.0, 10.0, 5.0, 0.0);

        let result = run_annual_simulation(&building, &config, &weather, None, None);

        // Annual simulation runs a steady-state heat balance every hour. Spot-check one hour.
        let record = &weather.records[0];
        let mut hour_config = config.clone();
        hour_config.outdoor_temperature = record.dry_bulb_temperature;
        hour_config.internal_gains = 0.0;
        hour_config.solar_gains = 0.0;
        let hour_result = calculate_heat_balance(&building, &hour_config);

        assert!(
            (result.hourly_heating[0] - hour_result.heating_demand).abs() < 1e-12,
            "Hour 0 heating mismatch: annual={} vs heat_balance={}",
            result.hourly_heating[0],
            hour_result.heating_demand
        );
        assert!(
            (result.hourly_cooling[0] - hour_result.cooling_demand).abs() < 1e-12,
            "Hour 0 cooling mismatch: annual={} vs heat_balance={}",
            result.hourly_cooling[0],
            hour_result.cooling_demand
        );
    }

    #[test]
    fn test_synthetic_weather_amplitude_zero_repeats_daily_cycle() {
        // Synthetic weather always has a ±3°C daily swing; with `temp_amplitude=0.0`,
        // the annual component is removed, so the day-to-day pattern repeats.
        let s = Solid::from_box(5.0, 5.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config = ThermalConfig::new();
        config.indoor_temperature = 20.0;

        let weather = WeatherData::synthetic("Const", 50.0, 10.0, 5.0, 0.0);
        let result = run_annual_simulation(&building, &config, &weather, None, None);

        // Compare day 1 vs day 2 (hours 0..24 vs 24..48).
        for h in 0..24 {
            let day1 = result.hourly_heating[h];
            let day2 = result.hourly_heating[h + 24];
            assert!(
                (day1 - day2).abs() < 1e-6,
                "Hour-of-day {h}: day1={day1}, day2={day2}"
            );
        }

        // All heating values should be non-negative
        for (i, &val) in result.hourly_heating.iter().enumerate() {
            assert!(
                val >= 0.0,
                "Hour {i}: heating should be non-negative, got {val}"
            );
        }
    }

    #[test]
    fn test_transient_energy_conservation() {
        // For the transient simulation, hourly sums should also match annual kWh.
        let s = Solid::from_box(5.0, 5.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let config = ThermalConfig::new();
        let weather = WeatherData::synthetic("Test", 52.0, 13.0, 10.0, 12.0);
        let hvac = HvacIdealLoads::new();

        let result = run_transient_simulation(&building, &config, &weather, &hvac, None, None);

        let sum_heating: f64 = result.hourly_heating.iter().sum();
        let sum_cooling: f64 = result.hourly_cooling.iter().sum();

        assert!(
            (sum_heating / 1000.0 - result.annual_heating_kwh).abs() < 0.01,
            "Transient hourly heating sum mismatch"
        );
        assert!(
            (sum_cooling / 1000.0 - result.annual_cooling_kwh).abs() < 0.01,
            "Transient hourly cooling sum mismatch"
        );
    }

    #[test]
    fn test_more_insulation_less_heating() {
        // Better insulated building should need less heating.
        let s = Solid::from_box(5.0, 5.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let weather = WeatherData::synthetic("Cold", 55.0, 10.0, 5.0, 12.0);

        let mut config_poor = ThermalConfig::new();
        config_poor.default_u_value = 3.0;

        let mut config_good = ThermalConfig::new();
        config_good.default_u_value = 0.5;

        let result_poor = run_annual_simulation(&building, &config_poor, &weather, None, None);
        let result_good = run_annual_simulation(&building, &config_good, &weather, None, None);

        assert!(
            result_good.annual_heating_kwh < result_poor.annual_heating_kwh,
            "Better insulation should reduce heating: good={}, poor={}",
            result_good.annual_heating_kwh,
            result_poor.annual_heating_kwh
        );
        // With U=0.5 vs U=3.0, the ratio of transmission losses is 6:1.
        // Infiltration is unchanged, so the total ratio is less than 6,
        // but the well-insulated building should use significantly less.
        assert!(
            result_good.annual_heating_kwh < result_poor.annual_heating_kwh * 0.5,
            "U=0.5 should use well under half the energy of U=3.0"
        );
    }
}
