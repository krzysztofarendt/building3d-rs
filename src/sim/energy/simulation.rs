use crate::Building;

use super::config::ThermalConfig;
use super::hvac::{HvacIdealLoads, LumpedThermalModel};
use super::schedule::InternalGainsProfile;
use super::solar_bridge::{SolarGainConfig, SolarHourParams, compute_solar_gains};
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

/// Computes day of year from month and day.
fn day_of_year(month: u8, day: u8) -> u16 {
    const DAYS_BEFORE_MONTH: [u16; 12] = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];
    let m = (month as usize).saturating_sub(1).min(11);
    DAYS_BEFORE_MONTH[m] + day as u16
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

        // Solar gains from weather data + sun geometry + glazing properties
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
                compute_solar_gains(building, &params, sc)
            }
            None => 0.0,
        };

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

    // Estimate thermal capacity from building volume.
    // The factor 50 kJ/(m^3*K) is a tuning parameter for medium-weight construction.
    // Typical range: ~30 kJ/(m^3*K) (lightweight) to ~80 kJ/(m^3*K) (heavyweight).
    // For more accurate results, derive from actual construction layer properties.
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
        let solar = match solar_config {
            Some(sc) => {
                let params = SolarHourParams {
                    direct_normal_irradiance: record.direct_normal_radiation,
                    diffuse_horizontal_irradiance: record.diffuse_horizontal_radiation,
                    day_of_year: day_of_year(record.month, record.day),
                    hour: record.hour as f64,
                    latitude: weather.latitude,
                    longitude: weather.longitude,
                };
                compute_solar_gains(building, &params, sc)
            }
            None => 0.0,
        };
        let total_gains = gains + solar;

        // Calculate HVAC power using implicit formula that accounts for
        // concurrent envelope losses during the timestep.
        let total_conductance = model.ua_total + model.infiltration_conductance;
        let (heating_power, cooling_power) = hvac.calculate_with_losses(
            model.zone_temperature,
            record.dry_bulb_temperature,
            total_conductance,
            total_gains,
            model.thermal_capacity,
            3600.0,
        );

        let hvac_net = heating_power - cooling_power;

        // Advance the thermal model
        model.step(record.dry_bulb_temperature, total_gains, hvac_net, 3600.0);

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
