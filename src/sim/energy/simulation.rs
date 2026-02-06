use crate::Building;

use super::config::ThermalConfig;
use super::hvac::{HvacIdealLoads, LumpedThermalModel};
use super::schedule::InternalGainsProfile;
use super::weather::WeatherData;
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

/// Runs an annual hourly energy simulation.
///
/// For each hour:
///   1. Set outdoor temperature from weather data
///   2. Calculate internal gains from schedule
///   3. Estimate solar gains from weather radiation
///   4. Run steady-state heat balance
///   5. Accumulate heating/cooling demand
pub fn run_annual_simulation(
    building: &Building,
    base_config: &ThermalConfig,
    weather: &WeatherData,
    gains_profile: Option<&InternalGainsProfile>,
    solar_gain_factor: f64,
) -> AnnualResult {
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

        // Solar gains (simplified: fraction of global horizontal radiation * window area)
        config.solar_gains = record.global_horizontal_radiation * solar_gain_factor;

        let result = calculate_heat_balance(building, &config);

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
pub fn run_transient_simulation(
    building: &Building,
    base_config: &ThermalConfig,
    weather: &WeatherData,
    hvac: &HvacIdealLoads,
    gains_profile: Option<&InternalGainsProfile>,
    solar_gain_factor: f64,
) -> AnnualResult {
    let num_hours = weather.num_hours();
    let mut hourly_heating = Vec::with_capacity(num_hours);
    let mut hourly_cooling = Vec::with_capacity(num_hours);

    // Compute building-level UA and thermal capacity
    let steady = calculate_heat_balance(building, base_config);
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
                        let path = format!(
                            "{}/{}/{}/{}",
                            zone.name, solid.name, wall.name, polygon.name
                        );
                        let u = base_config.resolve_u_value(&path);
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

    // Estimate thermal capacity from building volume (rough: 50 kJ/(m^3*K) for medium weight)
    let volume: f64 = building.zones().iter().map(|z| z.volume()).sum();
    let thermal_capacity = volume * 50000.0; // J/K

    let mut model = LumpedThermalModel::new(
        base_config.indoor_temperature,
        ua_total,
        infiltration_cond,
        thermal_capacity,
    );

    let mut annual_heating = 0.0;
    let mut annual_cooling = 0.0;
    let mut peak_heating = 0.0_f64;
    let mut peak_cooling = 0.0_f64;
    let mut monthly_heating = [0.0; 12];
    let mut monthly_cooling = [0.0; 12];

    for (hour_idx, record) in weather.records.iter().enumerate() {
        let gains = gains_profile.map(|p| p.gains_at(hour_idx)).unwrap_or(0.0);
        let solar = record.global_horizontal_radiation * solar_gain_factor;
        let total_gains = gains + solar;

        // Determine HVAC mode from free-floating temperature
        let setpoint = hvac.active_setpoint(model.zone_temperature);

        // Calculate required HVAC power to reach setpoint
        let q_loss = (model.ua_total + model.infiltration_conductance)
            * (model.zone_temperature - record.dry_bulb_temperature);
        let q_needed = q_loss - total_gains;

        let (heating_power, cooling_power) = if model.zone_temperature < hvac.heating_setpoint {
            (q_needed.max(0.0), 0.0)
        } else if model.zone_temperature > hvac.cooling_setpoint {
            (0.0, (-q_needed).max(0.0))
        } else {
            (0.0, 0.0)
        };

        let hvac_net = heating_power - cooling_power;

        // Advance the thermal model
        model.step(record.dry_bulb_temperature, total_gains, hvac_net, 3600.0);

        // Clamp temperature to setpoints if HVAC is active
        if heating_power > 0.0 && model.zone_temperature < setpoint {
            model.zone_temperature = setpoint;
        }
        if cooling_power > 0.0 && model.zone_temperature > setpoint {
            model.zone_temperature = setpoint;
        }

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Solid, Zone};

    #[test]
    fn test_annual_simulation_basic() {
        let s = Solid::from_box(5.0, 5.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let config = ThermalConfig::new();
        let weather = WeatherData::synthetic("Test", 52.0, 13.0, 10.0, 12.0);

        let result = run_annual_simulation(&building, &config, &weather, None, 0.0);

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

        let result_no_gains = run_annual_simulation(&building, &config, &weather, None, 0.0);
        let result_gains = run_annual_simulation(&building, &config, &weather, Some(&gains), 0.0);

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

        let result = run_annual_simulation(&building, &config, &weather, None, 0.0);

        let monthly_sum: f64 = result.monthly_heating_kwh.iter().sum();
        assert!(
            (monthly_sum - result.annual_heating_kwh).abs() < 1.0,
            "Monthly heating sum should equal annual: {monthly_sum} vs {}",
            result.annual_heating_kwh
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

        let result = run_transient_simulation(&building, &config, &weather, &hvac, None, 0.0);

        assert_eq!(result.hourly_heating.len(), 8760);
        assert!(
            result.annual_heating_kwh > 0.0,
            "Should need heating in cold climate"
        );
        assert!(result.peak_heating > 0.0, "Should have peak heating load");

        // Transient should have different results from steady-state due to thermal mass
        let steady = run_annual_simulation(&building, &config, &weather, None, 0.0);
        // They won't be identical, but both should indicate heating is needed
        assert!(steady.annual_heating_kwh > 0.0);
    }
}
