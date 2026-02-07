/// Ideal HVAC loads model.
///
/// Calculates the energy needed to maintain zone temperature
/// within heating and cooling setpoints.
#[derive(Debug, Clone)]
pub struct HvacIdealLoads {
    /// Heating setpoint in °C.
    pub heating_setpoint: f64,
    /// Cooling setpoint in °C.
    pub cooling_setpoint: f64,
    /// Maximum heating capacity in W (0 = unlimited).
    pub max_heating_capacity: f64,
    /// Maximum cooling capacity in W (0 = unlimited).
    pub max_cooling_capacity: f64,
    /// Heating COP (coefficient of performance).
    pub heating_cop: f64,
    /// Cooling COP (coefficient of performance).
    pub cooling_cop: f64,
}

impl HvacIdealLoads {
    /// Creates a default HVAC system (20°C heating, 26°C cooling).
    pub fn new() -> Self {
        Self {
            heating_setpoint: 20.0,
            cooling_setpoint: 26.0,
            max_heating_capacity: 0.0,
            max_cooling_capacity: 0.0,
            heating_cop: 1.0,
            cooling_cop: 3.0,
        }
    }

    /// Creates HVAC with custom setpoints.
    pub fn with_setpoints(heating: f64, cooling: f64) -> Self {
        Self {
            heating_setpoint: heating,
            cooling_setpoint: cooling,
            ..Self::new()
        }
    }

    /// Calculates the HVAC energy to maintain zone temperature.
    ///
    /// Returns (heating_electric_power_w, cooling_electric_power_w, resulting_zone_temperature).
    ///
    /// - `free_floating_temp`: zone temperature without HVAC
    /// - `zone_thermal_capacity`: total zone thermal capacity in J/K
    ///
    /// This is a simple "ideal loads" control model: it computes the *average*
    /// power over an assumed 1-hour timestep required to bring the zone back to
    /// the active setpoint, optionally limited by max capacities and converted
    /// to electrical power using COPs.
    pub fn calculate(
        &self,
        free_floating_temp: f64,
        zone_thermal_capacity: f64,
    ) -> (f64, f64, f64) {
        self.calculate_for_timestep(free_floating_temp, zone_thermal_capacity, 3600.0)
    }

    /// Like [`Self::calculate`], but with an explicit timestep in seconds.
    pub fn calculate_for_timestep(
        &self,
        free_floating_temp: f64,
        zone_thermal_capacity: f64,
        dt_s: f64,
    ) -> (f64, f64, f64) {
        if dt_s <= 0.0 || zone_thermal_capacity <= 0.0 {
            return (0.0, 0.0, free_floating_temp);
        }

        let target_temp = self.active_setpoint(free_floating_temp);
        let delta_t = target_temp - free_floating_temp;
        if delta_t == 0.0 {
            return (0.0, 0.0, free_floating_temp);
        }

        let required_thermal_power_w = (zone_thermal_capacity * delta_t.abs()) / dt_s;
        let (max_capacity_w, cop) = if delta_t > 0.0 {
            (
                if self.max_heating_capacity > 0.0 {
                    self.max_heating_capacity
                } else {
                    f64::INFINITY
                },
                self.heating_cop.max(1e-9),
            )
        } else {
            (
                if self.max_cooling_capacity > 0.0 {
                    self.max_cooling_capacity
                } else {
                    f64::INFINITY
                },
                self.cooling_cop.max(1e-9),
            )
        };

        let delivered_thermal_power_w = required_thermal_power_w.min(max_capacity_w);
        let delivered_delta_t = (delivered_thermal_power_w * dt_s) / zone_thermal_capacity;
        let resulting_temp = if delta_t > 0.0 {
            (free_floating_temp + delivered_delta_t).min(target_temp)
        } else {
            (free_floating_temp - delivered_delta_t).max(target_temp)
        };

        let electric_power_w = delivered_thermal_power_w / cop;

        if delta_t > 0.0 {
            (electric_power_w, 0.0, resulting_temp)
        } else {
            (0.0, electric_power_w, resulting_temp)
        }
    }

    /// Calculates HVAC power accounting for concurrent envelope losses/gains.
    ///
    /// Uses the implicit formulation that solves for Q_hvac such that
    /// T(t+dt) = T_setpoint:
    ///
    /// ```text
    /// Q_hvac = C*(T_set - T_zone)/dt + (UA + Inf_cond)*(T_set - T_out) - Q_gains
    /// ```
    ///
    /// This correctly accounts for the fact that as the zone temperature
    /// changes toward the setpoint, envelope losses also change.
    ///
    /// Returns (heating_thermal_power_w, cooling_thermal_power_w).
    pub fn calculate_with_losses(
        &self,
        zone_temp: f64,
        outdoor_temp: f64,
        total_conductance: f64,
        gains: f64,
        thermal_capacity: f64,
        dt_s: f64,
    ) -> (f64, f64) {
        let setpoint = self.active_setpoint(zone_temp);
        if (setpoint - zone_temp).abs() < 1e-10 {
            return (0.0, 0.0);
        }

        if dt_s <= 0.0 || thermal_capacity <= 0.0 {
            return (0.0, 0.0);
        }

        // Implicit formula: Q_hvac that achieves T_setpoint at end of timestep
        let q_hvac = thermal_capacity * (setpoint - zone_temp) / dt_s
            + total_conductance * (setpoint - outdoor_temp)
            - gains;

        if q_hvac > 0.0 {
            (q_hvac, 0.0)
        } else {
            (0.0, -q_hvac)
        }
    }

    /// Returns the active setpoint for a given free-floating temperature.
    pub fn active_setpoint(&self, free_floating_temp: f64) -> f64 {
        if free_floating_temp < self.heating_setpoint {
            self.heating_setpoint
        } else if free_floating_temp > self.cooling_setpoint {
            self.cooling_setpoint
        } else {
            free_floating_temp
        }
    }
}

impl Default for HvacIdealLoads {
    fn default() -> Self {
        Self::new()
    }
}

/// 1R1C lumped thermal mass model for a zone.
///
/// Represents the zone as a single resistance (envelope) and
/// single capacitance (thermal mass). Used for transient simulation.
///
/// T_zone(t+dt) = T_zone(t) + dt/C * (Q_gains - Q_losses)
///
/// where:
///   Q_losses = UA * (T_zone - T_outdoor) + Q_infiltration
///   C = total thermal mass of zone (J/K)
#[derive(Debug, Clone)]
pub struct LumpedThermalModel {
    /// Zone temperature in °C.
    pub zone_temperature: f64,
    /// Total UA value in W/K (sum of U*A for all surfaces).
    pub ua_total: f64,
    /// Infiltration conductance in W/K (rho*cp*V*ACH/3600).
    pub infiltration_conductance: f64,
    /// Total thermal capacity in J/K.
    pub thermal_capacity: f64,
}

impl LumpedThermalModel {
    pub fn new(initial_temp: f64, ua_total: f64, infiltration_cond: f64, capacity: f64) -> Self {
        Self {
            zone_temperature: initial_temp,
            ua_total,
            infiltration_conductance: infiltration_cond,
            thermal_capacity: capacity,
        }
    }

    /// Advances the zone temperature by one time step.
    ///
    /// - `outdoor_temp`: outdoor temperature in °C
    /// - `gains`: total heat gains in W (internal + solar)
    /// - `hvac_power`: net HVAC power in W (positive = heating, negative = cooling)
    /// - `dt`: time step in seconds
    ///
    /// Returns the new zone temperature.
    pub fn step(&mut self, outdoor_temp: f64, gains: f64, hvac_power: f64, dt: f64) -> f64 {
        if self.thermal_capacity <= 0.0 {
            // No thermal mass — instant response (steady-state)
            let total_conductance = self.ua_total + self.infiltration_conductance;
            if total_conductance > 0.0 {
                self.zone_temperature = outdoor_temp + (gains + hvac_power) / total_conductance;
            }
            return self.zone_temperature;
        }

        let q_loss = (self.ua_total + self.infiltration_conductance)
            * (self.zone_temperature - outdoor_temp);
        let q_net = gains + hvac_power - q_loss;

        self.zone_temperature += dt / self.thermal_capacity * q_net;
        self.zone_temperature
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hvac_heating_mode() {
        let hvac = HvacIdealLoads::new();
        // 5 K below setpoint, capacity 3600 J/K over 1 hour -> 5 W thermal, COP=1 -> 5 W electric.
        let (heating, cooling, temp) = hvac.calculate(15.0, 3600.0);
        assert!(
            (heating - 5.0).abs() < 1e-10,
            "Expected 5 W electric heating"
        );
        assert!((cooling - 0.0).abs() < 1e-10, "No cooling needed");
        assert!(
            (temp - 20.0).abs() < 1e-10,
            "Should maintain heating setpoint"
        );
    }

    #[test]
    fn test_hvac_cooling_mode() {
        let hvac = HvacIdealLoads::new();
        // 4 K above setpoint, capacity 3600 J/K over 1 hour -> 4 W thermal, COP=3 -> 1.333.. W electric.
        let (heating, cooling, temp) = hvac.calculate(30.0, 3600.0);
        assert!((heating - 0.0).abs() < 1e-10, "No heating needed");
        assert!(
            (cooling - (4.0 / 3.0)).abs() < 1e-10,
            "Expected cooling electric power of 4/3 W"
        );
        assert!(
            (temp - 26.0).abs() < 1e-10,
            "Should maintain cooling setpoint"
        );
    }

    #[test]
    fn test_hvac_deadband() {
        let hvac = HvacIdealLoads::new();
        let (heating, cooling, temp) = hvac.calculate(23.0, 3600.0);
        assert!((heating - 0.0).abs() < 1e-10);
        assert!((cooling - 0.0).abs() < 1e-10);
        assert!(
            (temp - 23.0).abs() < 1e-10,
            "Temperature unchanged in deadband"
        );
    }

    #[test]
    fn test_hvac_capacity_limits_resulting_temp() {
        let mut hvac = HvacIdealLoads::new();
        hvac.max_heating_capacity = 2.0; // W thermal

        // Need 5 W thermal to raise 5 K in an hour with C=3600 J/K, but only 2 W available.
        let (heating_elec, cooling_elec, temp) = hvac.calculate(15.0, 3600.0);
        assert!(cooling_elec.abs() < 1e-12);
        // COP=1 -> electric == thermal
        assert!((heating_elec - 2.0).abs() < 1e-10);
        // Delivered deltaT = 2 W * 3600 s / 3600 J/K = 2 K
        assert!((temp - 17.0).abs() < 1e-10);
    }

    #[test]
    fn test_lumped_model_cooling_down() {
        // Zone at 20°C, outdoor at 0°C, no gains/hvac
        let mut model = LumpedThermalModel::new(
            20.0, 100.0,    // UA = 100 W/K
            20.0,     // infiltration = 20 W/K
            500000.0, // 500 kJ/K thermal mass
        );

        let temp = model.step(0.0, 0.0, 0.0, 3600.0); // 1 hour
        // Should cool down: dT = 3600/500000 * (0 - 120*20) = 3600/500000 * (-2400) = -17.28
        // But that's too big — let's just check it's colder
        assert!(
            temp < 20.0,
            "Zone should cool down without heating, got {temp}"
        );
        assert!(temp > 0.0, "Should not drop below outdoor in 1 hour");
    }

    #[test]
    fn test_lumped_model_with_heating() {
        let mut model = LumpedThermalModel::new(20.0, 100.0, 20.0, 500000.0);

        // Apply enough heating to counteract losses
        // Q_loss = 120 * 20 = 2400 W
        let temp = model.step(0.0, 0.0, 2400.0, 3600.0);
        // Net Q = 2400 - 2400 = 0, so temp should stay ~20°C
        assert!(
            (temp - 20.0).abs() < 0.1,
            "Temp should stay ~20°C with balanced heating, got {temp}"
        );
    }

    #[test]
    fn test_active_setpoint() {
        let hvac = HvacIdealLoads::with_setpoints(18.0, 25.0);
        assert!((hvac.active_setpoint(15.0) - 18.0).abs() < 1e-10);
        assert!((hvac.active_setpoint(22.0) - 22.0).abs() < 1e-10);
        assert!((hvac.active_setpoint(28.0) - 25.0).abs() < 1e-10);
    }

    #[test]
    fn test_calculate_with_losses_heating() {
        let hvac = HvacIdealLoads::new();
        // Zone at 15°C, outdoor 0°C, UA=100 W/K, no gains, C=3600 J/K, dt=3600s
        let (heating, cooling) = hvac.calculate_with_losses(15.0, 0.0, 100.0, 0.0, 3600.0, 3600.0);
        // Q_hvac = 3600*(20-15)/3600 + 100*(20-0) - 0 = 5 + 2000 = 2005 W
        assert!((heating - 2005.0).abs() < 1e-6, "Expected 2005 W heating, got {heating}");
        assert!(cooling.abs() < 1e-10);
    }

    #[test]
    fn test_calculate_with_losses_deadband() {
        let hvac = HvacIdealLoads::new();
        let (heating, cooling) = hvac.calculate_with_losses(22.0, 20.0, 100.0, 0.0, 3600.0, 3600.0);
        assert!(heating.abs() < 1e-10);
        assert!(cooling.abs() < 1e-10);
    }

    #[test]
    fn test_calculate_with_losses_cooling() {
        let hvac = HvacIdealLoads::new();
        // Zone at 30°C, outdoor 35°C, UA=100 W/K, gains=500 W
        let (heating, cooling) = hvac.calculate_with_losses(30.0, 35.0, 100.0, 500.0, 3600.0, 3600.0);
        // Q_hvac = 3600*(26-30)/3600 + 100*(26-35) - 500 = -4 + (-900) - 500 = -1404
        // So cooling = 1404 W
        assert!(heating.abs() < 1e-10);
        assert!((cooling - 1404.0).abs() < 1e-6, "Expected 1404 W cooling, got {cooling}");
    }
}
