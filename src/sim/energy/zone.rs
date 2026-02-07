use crate::Building;

use super::config::ThermalConfig;
use super::result::{ThermalResult, ZoneResult};

/// Specific heat capacity of air in J/(kg*K).
const AIR_SPECIFIC_HEAT: f64 = 1005.0;
/// Density of air in kg/m^3 (at ~20°C).
const AIR_DENSITY: f64 = 1.2;

/// Performs a steady-state thermal zone heat balance.
///
/// For each zone in the building:
///   Q_transmission = sum(U_i * A_i * dT) for all exterior surfaces
///   Q_infiltration = rho * cp * V * ACH / 3600 * dT
///   Q_gains = internal + solar
///   Q_net = Q_transmission + Q_infiltration - Q_gains
pub fn calculate_heat_balance(building: &Building, config: &ThermalConfig) -> ThermalResult {
    let mut result = ThermalResult::new();
    let dt = config.indoor_temperature - config.outdoor_temperature;

    for zone in building.zones() {
        let zone_name = zone.name.clone();
        let volume = zone.volume();
        let mut zone_transmission = 0.0;
        let mut zone_area = 0.0;

        for solid in zone.solids() {
            for wall in solid.walls() {
                for polygon in wall.polygons() {
                    let path = format!(
                        "{}/{}/{}/{}",
                        zone.name, solid.name, wall.name, polygon.name
                    );
                    let u_value = config.resolve_u_value(&path);
                    let area = polygon.area();
                    let q_surface = u_value * area * dt;

                    result.surface_heat_loss.insert(path, q_surface);
                    zone_transmission += q_surface;
                    zone_area += area;
                }
            }
        }

        // Infiltration loss: rho * cp * V * ACH / 3600 * dT
        let zone_infiltration =
            AIR_DENSITY * AIR_SPECIFIC_HEAT * volume * config.infiltration_ach / 3600.0 * dt;

        result.transmission_loss += zone_transmission;
        result.infiltration_loss += zone_infiltration;

        let net = zone_transmission + zone_infiltration;
        result.zone_results.insert(
            zone_name.clone(),
            ZoneResult {
                zone_name,
                volume,
                envelope_area: zone_area,
                transmission_loss: zone_transmission,
                infiltration_loss: zone_infiltration,
                net_demand: net,
            },
        );
    }

    result.total_gains = config.internal_gains + config.solar_gains;

    let net = result.transmission_loss + result.infiltration_loss - result.total_gains;
    if net > 0.0 {
        result.heating_demand = net;
        result.cooling_demand = 0.0;
    } else {
        result.heating_demand = 0.0;
        result.cooling_demand = -net;
    }

    result
}

/// Computes the mean U-value of the building envelope.
///
/// U_mean = sum(U_i * A_i) / sum(A_i)
pub fn mean_u_value(building: &Building, config: &ThermalConfig) -> f64 {
    let mut sum_ua = 0.0;
    let mut sum_a = 0.0;

    for zone in building.zones() {
        for solid in zone.solids() {
            for wall in solid.walls() {
                for polygon in wall.polygons() {
                    let path = format!(
                        "{}/{}/{}/{}",
                        zone.name, solid.name, wall.name, polygon.name
                    );
                    let u = config.resolve_u_value(&path);
                    let a = polygon.area();
                    sum_ua += u * a;
                    sum_a += a;
                }
            }
        }
    }

    if sum_a > 0.0 { sum_ua / sum_a } else { 0.0 }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Solid, Zone};

    #[test]
    fn test_heat_balance_simple_box() {
        // 5x5x3m room, default U=2.0, dT=20°C
        let s = Solid::from_box(5.0, 5.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let config = ThermalConfig::new(); // dT = 20 - 0 = 20°C
        let result = calculate_heat_balance(&building, &config);

        // Total area: 2*(5*5 + 5*3 + 5*3) = 2*(25+15+15) = 110 m^2
        // Q_transmission = 2.0 * 110 * 20 = 4400 W
        let expected_transmission = 2.0 * 110.0 * 20.0;
        assert!(
            (result.transmission_loss - expected_transmission).abs() < 1.0,
            "Transmission loss should be ~{expected_transmission}, got {}",
            result.transmission_loss
        );

        // Q_infiltration = 1.2 * 1005 * 75 * 0.5 / 3600 * 20 ≈ 251.25 W
        let expected_infiltration = AIR_DENSITY * AIR_SPECIFIC_HEAT * 75.0 * 0.5 / 3600.0 * 20.0;
        assert!(
            (result.infiltration_loss - expected_infiltration).abs() < 1.0,
            "Infiltration loss should be ~{expected_infiltration}, got {}",
            result.infiltration_loss
        );

        assert!(
            result.heating_demand > 0.0,
            "Should need heating when outdoor < indoor"
        );
        assert!(
            (result.cooling_demand - 0.0).abs() < 1e-10,
            "Should not need cooling"
        );
    }

    #[test]
    fn test_heat_balance_no_dt() {
        let s = Solid::from_box(3.0, 3.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config = ThermalConfig::new();
        config.outdoor_temperature = 20.0; // same as indoor
        let result = calculate_heat_balance(&building, &config);

        assert!(
            result.transmission_loss.abs() < 1e-10,
            "No temperature difference, no loss"
        );
        assert!(
            result.infiltration_loss.abs() < 1e-10,
            "No temperature difference, no infiltration loss"
        );
    }

    #[test]
    fn test_cooling_demand() {
        let s = Solid::from_box(3.0, 3.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config = ThermalConfig::new();
        config.outdoor_temperature = 35.0; // hot outside
        config.solar_gains = 5000.0;
        let result = calculate_heat_balance(&building, &config);

        // When outdoor > indoor, transmission is negative (heat flows in)
        // Plus large solar gains -> cooling demand
        assert!(
            result.cooling_demand > 0.0,
            "Should need cooling with high outdoor temp and solar gains"
        );
    }

    #[test]
    fn test_mean_u_value() {
        let s = Solid::from_box(3.0, 3.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let config = ThermalConfig::new();
        let u_mean = mean_u_value(&building, &config);
        assert!(
            (u_mean - config.default_u_value).abs() < 1e-10,
            "All surfaces use default, mean should equal default"
        );
    }

    // ── Physics verification tests ──────────────────────────────────────

    #[test]
    fn test_doubling_insulation_halves_transmission() {
        let s = Solid::from_box(5.0, 5.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config_a = ThermalConfig::new();
        config_a.default_u_value = 2.0;
        let result_a = calculate_heat_balance(&building, &config_a);

        let mut config_b = ThermalConfig::new();
        config_b.default_u_value = 1.0; // half the U-value
        let result_b = calculate_heat_balance(&building, &config_b);

        let ratio = result_a.transmission_loss / result_b.transmission_loss;
        assert!(
            (ratio - 2.0).abs() < 1e-10,
            "Halving U-value should halve transmission loss, ratio={ratio}"
        );
    }

    #[test]
    fn test_doubling_volume_doubles_infiltration() {
        // Room A: 5x5x3 = 75 m³
        let s_a = Solid::from_box(5.0, 5.0, 3.0, None, "room").unwrap();
        let z_a = Zone::new("z", vec![s_a]).unwrap();
        let b_a = Building::new("b", vec![z_a]).unwrap();

        // Room B: 10x5x3 = 150 m³ (double the volume)
        let s_b = Solid::from_box(10.0, 5.0, 3.0, None, "room").unwrap();
        let z_b = Zone::new("z", vec![s_b]).unwrap();
        let b_b = Building::new("b", vec![z_b]).unwrap();

        let config = ThermalConfig::new();
        let result_a = calculate_heat_balance(&b_a, &config);
        let result_b = calculate_heat_balance(&b_b, &config);

        let ratio = result_b.infiltration_loss / result_a.infiltration_loss;
        assert!(
            (ratio - 2.0).abs() < 1e-6,
            "Double volume should double infiltration loss, ratio={ratio}"
        );
    }

    #[test]
    fn test_heating_cooling_symmetry() {
        // With no gains, heating demand at dT=+20 should equal cooling demand at dT=-20.
        let s = Solid::from_box(4.0, 4.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config_heat = ThermalConfig::new();
        config_heat.indoor_temperature = 20.0;
        config_heat.outdoor_temperature = 0.0; // dT = +20
        config_heat.internal_gains = 0.0;
        config_heat.solar_gains = 0.0;

        let mut config_cool = ThermalConfig::new();
        config_cool.indoor_temperature = 20.0;
        config_cool.outdoor_temperature = 40.0; // dT = -20
        config_cool.internal_gains = 0.0;
        config_cool.solar_gains = 0.0;

        let result_heat = calculate_heat_balance(&building, &config_heat);
        let result_cool = calculate_heat_balance(&building, &config_cool);

        assert!(
            (result_heat.heating_demand - result_cool.cooling_demand).abs() < 1e-10,
            "Symmetric dT should give equal demand: heating={}, cooling={}",
            result_heat.heating_demand,
            result_cool.cooling_demand
        );
    }

    #[test]
    fn test_energy_conservation_steady_state() {
        // In steady-state: Q_heating = Q_transmission + Q_infiltration - Q_gains
        let s = Solid::from_box(6.0, 4.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config = ThermalConfig::new();
        config.indoor_temperature = 21.0;
        config.outdoor_temperature = -5.0;
        config.internal_gains = 500.0;
        config.solar_gains = 200.0;

        let result = calculate_heat_balance(&building, &config);

        let net_losses = result.transmission_loss + result.infiltration_loss;
        let net = net_losses - result.total_gains;
        // net > 0 → heating, net < 0 → cooling
        assert!(
            (result.heating_demand - result.cooling_demand - net).abs() < 1e-10,
            "Energy balance: heating={}, cooling={}, net_losses={}, gains={}",
            result.heating_demand,
            result.cooling_demand,
            net_losses,
            result.total_gains
        );
    }

    #[test]
    fn test_transmission_proportional_to_dt() {
        // Transmission loss should be exactly proportional to (T_in - T_out).
        let s = Solid::from_box(3.0, 3.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config_10 = ThermalConfig::new();
        config_10.indoor_temperature = 20.0;
        config_10.outdoor_temperature = 10.0; // dT = 10

        let mut config_30 = ThermalConfig::new();
        config_30.indoor_temperature = 20.0;
        config_30.outdoor_temperature = -10.0; // dT = 30

        let r10 = calculate_heat_balance(&building, &config_10);
        let r30 = calculate_heat_balance(&building, &config_30);

        let ratio = r30.transmission_loss / r10.transmission_loss;
        assert!(
            (ratio - 3.0).abs() < 1e-10,
            "Transmission should scale linearly with dT, ratio={ratio}"
        );
    }
}
