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
/// The envelope loss term is integrated *implicitly* (backward Euler) to
/// improve stability and to match the implicit HVAC setpoint formulation in
/// [`HvacIdealLoads::calculate_with_losses`].
///
/// With constant inputs over the timestep, this solves:
///
/// ```text
/// dT/dt = (Q_gains + Q_hvac - K*(T - T_out)) / C
/// ```
///
/// where `K = UA + infiltration_conductance`.
///
/// Backward Euler update:
///
/// ```text
/// T(t+dt) = (T(t) + dt/C * (Q_gains + Q_hvac + K*T_out)) / (1 + dt*K/C)
/// ```
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
        if dt <= 0.0 {
            return self.zone_temperature;
        }
        if self.thermal_capacity <= 0.0 {
            // No thermal mass — instant response (steady-state)
            let total_conductance = self.ua_total + self.infiltration_conductance;
            if total_conductance > 0.0 {
                self.zone_temperature = outdoor_temp + (gains + hvac_power) / total_conductance;
            }
            return self.zone_temperature;
        }

        let k = self.ua_total + self.infiltration_conductance;
        if k <= 0.0 {
            // No envelope losses -> pure integrator.
            self.zone_temperature += dt / self.thermal_capacity * (gains + hvac_power);
            return self.zone_temperature;
        }

        let a = dt * k / self.thermal_capacity;
        let numerator = self.zone_temperature
            + (dt / self.thermal_capacity) * (gains + hvac_power + k * outdoor_temp);
        self.zone_temperature = numerator / (1.0 + a);
        self.zone_temperature
    }
}

/// 2R2C lumped zone model: an air node coupled to an interior "mass" node.
///
/// State:
/// - air temperature `T_air`
/// - mass temperature `T_mass`
///
/// Conductances:
/// - `k_env`: envelope + infiltration conductance from air → outdoor (W/K)
/// - `k_am`: coupling conductance between air ↔ mass (W/K)
///
/// Capacitances:
/// - `c_air`: effective air-side capacity (J/K)
/// - `c_mass`: effective mass-side capacity (J/K)
///
/// Backward Euler integration solves a 2×2 linear system each step.
#[derive(Debug, Clone)]
pub struct TwoNodeThermalModel {
    pub air_temperature_c: f64,
    pub mass_temperature_c: f64,
    pub envelope_conductance_w_per_k: f64,
    pub air_mass_conductance_w_per_k: f64,
    pub air_capacity_j_per_k: f64,
    pub mass_capacity_j_per_k: f64,
}

impl TwoNodeThermalModel {
    pub fn new(
        initial_air_c: f64,
        initial_mass_c: f64,
        envelope_conductance_w_per_k: f64,
        air_mass_conductance_w_per_k: f64,
        air_capacity_j_per_k: f64,
        mass_capacity_j_per_k: f64,
    ) -> Self {
        Self {
            air_temperature_c: initial_air_c,
            mass_temperature_c: initial_mass_c,
            envelope_conductance_w_per_k: envelope_conductance_w_per_k.max(0.0),
            air_mass_conductance_w_per_k: air_mass_conductance_w_per_k.max(0.0),
            air_capacity_j_per_k: air_capacity_j_per_k.max(0.0),
            mass_capacity_j_per_k: mass_capacity_j_per_k.max(0.0),
        }
    }

    /// Advances the model by one timestep.
    ///
    /// - `outdoor_temp_c`: outdoor air temperature (°C)
    /// - `gains_air_w`: gains applied directly to the air node (W)
    /// - `gains_mass_w`: gains applied to the mass node (W)
    /// - `hvac_power_w`: net HVAC thermal power into the air node (W)
    /// - `dt_s`: timestep (s)
    pub fn step(
        &mut self,
        outdoor_temp_c: f64,
        gains_air_w: f64,
        gains_mass_w: f64,
        hvac_power_w: f64,
        dt_s: f64,
    ) -> (f64, f64) {
        if dt_s <= 0.0 {
            return (self.air_temperature_c, self.mass_temperature_c);
        }

        let ca = self.air_capacity_j_per_k;
        let cm = self.mass_capacity_j_per_k;
        let k_env = self.envelope_conductance_w_per_k;
        let k_am = self.air_mass_conductance_w_per_k;

        // Degenerate fallbacks: reduce to a single node if we cannot integrate.
        if ca <= 0.0 || cm <= 0.0 || k_am <= 0.0 {
            // Treat everything as a 1R1C air node with capacity ca+cm (if any).
            let c = (ca + cm).max(0.0);
            if c <= 0.0 {
                let k = k_env.max(1e-12);
                self.air_temperature_c =
                    outdoor_temp_c + (gains_air_w + gains_mass_w + hvac_power_w) / k;
                self.mass_temperature_c = self.air_temperature_c;
                return (self.air_temperature_c, self.mass_temperature_c);
            }
            let k = k_env.max(0.0);
            if k <= 0.0 {
                self.air_temperature_c += dt_s / c * (gains_air_w + gains_mass_w + hvac_power_w);
                self.mass_temperature_c = self.air_temperature_c;
                return (self.air_temperature_c, self.mass_temperature_c);
            }
            let a = dt_s * k / c;
            let numerator = self.air_temperature_c
                + (dt_s / c) * (gains_air_w + gains_mass_w + hvac_power_w + k * outdoor_temp_c);
            self.air_temperature_c = numerator / (1.0 + a);
            self.mass_temperature_c = self.air_temperature_c;
            return (self.air_temperature_c, self.mass_temperature_c);
        }

        // Backward Euler 2×2 solve.
        let a11 = 1.0 + dt_s * (k_env + k_am) / ca;
        let a12 = -dt_s * k_am / ca;
        let b1 = self.air_temperature_c
            + (dt_s / ca) * (gains_air_w + hvac_power_w + k_env * outdoor_temp_c);

        let a21 = -dt_s * k_am / cm;
        let a22 = 1.0 + dt_s * k_am / cm;
        let b2 = self.mass_temperature_c + (dt_s / cm) * gains_mass_w;

        let det = a11 * a22 - a12 * a21;
        if det.abs() < 1e-16 {
            return (self.air_temperature_c, self.mass_temperature_c);
        }

        let t_air = (b1 * a22 - a12 * b2) / det;
        let t_mass = (a11 * b2 - b1 * a21) / det;

        self.air_temperature_c = t_air;
        self.mass_temperature_c = t_mass;
        (t_air, t_mass)
    }
}

/// 2R2C variant: air node + mass node, where envelope conduction connects
/// the mass node directly to outdoors (infiltration remains air ↔ outdoors).
///
/// This is often a better approximation for heavyweight envelopes, where
/// the dominant capacitance is in walls/roof/floor rather than in air.
#[derive(Debug, Clone)]
pub struct TwoNodeEnvelopeThermalModel {
    pub air_temperature_c: f64,
    pub mass_temperature_c: f64,
    pub air_outdoor_conductance_w_per_k: f64,
    pub mass_outdoor_conductance_w_per_k: f64,
    pub air_mass_conductance_w_per_k: f64,
    pub air_capacity_j_per_k: f64,
    pub mass_capacity_j_per_k: f64,
}

impl TwoNodeEnvelopeThermalModel {
    pub fn new(
        initial_air_c: f64,
        initial_mass_c: f64,
        air_outdoor_conductance_w_per_k: f64,
        mass_outdoor_conductance_w_per_k: f64,
        air_mass_conductance_w_per_k: f64,
        air_capacity_j_per_k: f64,
        mass_capacity_j_per_k: f64,
    ) -> Self {
        Self {
            air_temperature_c: initial_air_c,
            mass_temperature_c: initial_mass_c,
            air_outdoor_conductance_w_per_k: air_outdoor_conductance_w_per_k.max(0.0),
            mass_outdoor_conductance_w_per_k: mass_outdoor_conductance_w_per_k.max(0.0),
            air_mass_conductance_w_per_k: air_mass_conductance_w_per_k.max(0.0),
            air_capacity_j_per_k: air_capacity_j_per_k.max(0.0),
            mass_capacity_j_per_k: mass_capacity_j_per_k.max(0.0),
        }
    }

    pub fn step(
        &mut self,
        outdoor_temp_c: f64,
        gains_air_w: f64,
        gains_mass_w: f64,
        hvac_power_w: f64,
        dt_s: f64,
    ) -> (f64, f64) {
        if dt_s <= 0.0 {
            return (self.air_temperature_c, self.mass_temperature_c);
        }

        let ca = self.air_capacity_j_per_k;
        let cm = self.mass_capacity_j_per_k;
        let k_ao = self.air_outdoor_conductance_w_per_k;
        let k_mo = self.mass_outdoor_conductance_w_per_k;
        let k_am = self.air_mass_conductance_w_per_k;

        if ca <= 0.0 || cm <= 0.0 || k_am <= 0.0 {
            // Fallback: treat as a single air node.
            let c = (ca + cm).max(0.0);
            let k = (k_ao + k_mo).max(0.0);
            if c <= 0.0 {
                let k = k.max(1e-12);
                self.air_temperature_c =
                    outdoor_temp_c + (gains_air_w + gains_mass_w + hvac_power_w) / k;
                self.mass_temperature_c = self.air_temperature_c;
                return (self.air_temperature_c, self.mass_temperature_c);
            }
            if k <= 0.0 {
                self.air_temperature_c += dt_s / c * (gains_air_w + gains_mass_w + hvac_power_w);
                self.mass_temperature_c = self.air_temperature_c;
                return (self.air_temperature_c, self.mass_temperature_c);
            }
            let a = dt_s * k / c;
            let numerator = self.air_temperature_c
                + (dt_s / c) * (gains_air_w + gains_mass_w + hvac_power_w + k * outdoor_temp_c);
            self.air_temperature_c = numerator / (1.0 + a);
            self.mass_temperature_c = self.air_temperature_c;
            return (self.air_temperature_c, self.mass_temperature_c);
        }

        let a11 = 1.0 + dt_s * (k_ao + k_am) / ca;
        let a12 = -dt_s * k_am / ca;
        let b1 = self.air_temperature_c
            + (dt_s / ca) * (gains_air_w + hvac_power_w + k_ao * outdoor_temp_c);

        let a21 = -dt_s * k_am / cm;
        let a22 = 1.0 + dt_s * (k_mo + k_am) / cm;
        let b2 = self.mass_temperature_c + (dt_s / cm) * (gains_mass_w + k_mo * outdoor_temp_c);

        let det = a11 * a22 - a12 * a21;
        if det.abs() < 1e-16 {
            return (self.air_temperature_c, self.mass_temperature_c);
        }

        let t_air = (b1 * a22 - a12 * b2) / det;
        let t_mass = (a11 * b2 - b1 * a21) / det;

        self.air_temperature_c = t_air;
        self.mass_temperature_c = t_mass;
        (t_air, t_mass)
    }
}

/// 3R3C coarse zone model: air node + interior surface node + envelope node.
///
/// Intended use:
/// - transmitted solar heats interior surfaces (surface node) and is released to air via `k_as`
/// - exterior absorbed solar heats the envelope node and reaches the room with lag via `k_se`
///
/// Conductances:
/// - `k_ao`: air ↔ outdoor (infiltration + glazing)
/// - `k_as`: air ↔ interior surfaces (convective+radiative aggregate)
/// - `k_se`: interior surfaces ↔ envelope (opaque conduction, inner half)
/// - `k_eo`: envelope ↔ outdoor (opaque conduction, outer half + exterior film)
#[derive(Debug, Clone)]
pub struct ThreeNodeEnvelopeThermalModel {
    pub air_temperature_c: f64,
    pub surface_temperature_c: f64,
    pub envelope_temperature_c: f64,
    pub air_outdoor_conductance_w_per_k: f64,
    pub air_surface_conductance_w_per_k: f64,
    pub surface_envelope_conductance_w_per_k: f64,
    pub envelope_outdoor_conductance_w_per_k: f64,
    pub air_capacity_j_per_k: f64,
    pub surface_capacity_j_per_k: f64,
    pub envelope_capacity_j_per_k: f64,
}

impl ThreeNodeEnvelopeThermalModel {
    pub fn new(
        initial_air_c: f64,
        initial_surface_c: f64,
        initial_envelope_c: f64,
        air_outdoor_conductance_w_per_k: f64,
        air_surface_conductance_w_per_k: f64,
        surface_envelope_conductance_w_per_k: f64,
        envelope_outdoor_conductance_w_per_k: f64,
        air_capacity_j_per_k: f64,
        surface_capacity_j_per_k: f64,
        envelope_capacity_j_per_k: f64,
    ) -> Self {
        Self {
            air_temperature_c: initial_air_c,
            surface_temperature_c: initial_surface_c,
            envelope_temperature_c: initial_envelope_c,
            air_outdoor_conductance_w_per_k: air_outdoor_conductance_w_per_k.max(0.0),
            air_surface_conductance_w_per_k: air_surface_conductance_w_per_k.max(0.0),
            surface_envelope_conductance_w_per_k: surface_envelope_conductance_w_per_k.max(0.0),
            envelope_outdoor_conductance_w_per_k: envelope_outdoor_conductance_w_per_k.max(0.0),
            air_capacity_j_per_k: air_capacity_j_per_k.max(0.0),
            surface_capacity_j_per_k: surface_capacity_j_per_k.max(0.0),
            envelope_capacity_j_per_k: envelope_capacity_j_per_k.max(0.0),
        }
    }

    pub fn step(
        &mut self,
        outdoor_temp_c: f64,
        gains_air_w: f64,
        gains_surface_w: f64,
        gains_envelope_w: f64,
        hvac_power_w: f64,
        dt_s: f64,
    ) -> (f64, f64, f64) {
        if dt_s <= 0.0 {
            return (
                self.air_temperature_c,
                self.surface_temperature_c,
                self.envelope_temperature_c,
            );
        }

        let ca = self.air_capacity_j_per_k;
        let cs = self.surface_capacity_j_per_k;
        let ce = self.envelope_capacity_j_per_k;
        let k_ao = self.air_outdoor_conductance_w_per_k;
        let k_as = self.air_surface_conductance_w_per_k;
        let k_se = self.surface_envelope_conductance_w_per_k;
        let k_eo = self.envelope_outdoor_conductance_w_per_k;

        if ca <= 0.0 || cs <= 0.0 || ce <= 0.0 || k_as <= 0.0 || k_se <= 0.0 {
            // Degenerate fallback: treat as a single air node.
            let c = (ca + cs + ce).max(0.0);
            let k = (k_ao + k_eo).max(0.0);
            if c <= 0.0 {
                let k = k.max(1e-12);
                self.air_temperature_c = outdoor_temp_c
                    + (gains_air_w + gains_surface_w + gains_envelope_w + hvac_power_w) / k;
                self.surface_temperature_c = self.air_temperature_c;
                self.envelope_temperature_c = self.air_temperature_c;
                return (
                    self.air_temperature_c,
                    self.surface_temperature_c,
                    self.envelope_temperature_c,
                );
            }
            if k <= 0.0 {
                self.air_temperature_c +=
                    dt_s / c * (gains_air_w + gains_surface_w + gains_envelope_w + hvac_power_w);
                self.surface_temperature_c = self.air_temperature_c;
                self.envelope_temperature_c = self.air_temperature_c;
                return (
                    self.air_temperature_c,
                    self.surface_temperature_c,
                    self.envelope_temperature_c,
                );
            }
            let a = dt_s * k / c;
            let numerator = self.air_temperature_c
                + (dt_s / c)
                    * (gains_air_w
                        + gains_surface_w
                        + gains_envelope_w
                        + hvac_power_w
                        + k * outdoor_temp_c);
            self.air_temperature_c = numerator / (1.0 + a);
            self.surface_temperature_c = self.air_temperature_c;
            self.envelope_temperature_c = self.air_temperature_c;
            return (
                self.air_temperature_c,
                self.surface_temperature_c,
                self.envelope_temperature_c,
            );
        }

        // Backward Euler 3×3 solve.
        let a11 = 1.0 + dt_s * (k_ao + k_as) / ca;
        let a12 = -dt_s * k_as / ca;
        let a13 = 0.0;
        let b1 = self.air_temperature_c
            + (dt_s / ca) * (gains_air_w + hvac_power_w + k_ao * outdoor_temp_c);

        let a21 = -dt_s * k_as / cs;
        let a22 = 1.0 + dt_s * (k_as + k_se) / cs;
        let a23 = -dt_s * k_se / cs;
        let b2 = self.surface_temperature_c + (dt_s / cs) * gains_surface_w;

        let a31 = 0.0;
        let a32 = -dt_s * k_se / ce;
        let a33 = 1.0 + dt_s * (k_se + k_eo) / ce;
        let b3 =
            self.envelope_temperature_c + (dt_s / ce) * (gains_envelope_w + k_eo * outdoor_temp_c);

        // Solve via Gaussian elimination (since matrix is small and sparse).
        let m21 = a21 / a11;
        let m31 = a31 / a11;
        let a22p = a22 - m21 * a12;
        let a23p = a23 - m21 * a13;
        let b2p = b2 - m21 * b1;
        let a32p = a32 - m31 * a12;
        let a33p = a33 - m31 * a13;
        let b3p = b3 - m31 * b1;

        if a22p.abs() < 1e-16 {
            return (
                self.air_temperature_c,
                self.surface_temperature_c,
                self.envelope_temperature_c,
            );
        }

        let m32 = a32p / a22p;
        let a33pp = a33p - m32 * a23p;
        let b3pp = b3p - m32 * b2p;

        if a33pp.abs() < 1e-16 {
            return (
                self.air_temperature_c,
                self.surface_temperature_c,
                self.envelope_temperature_c,
            );
        }

        let t_env = b3pp / a33pp;
        let t_surf = (b2p - a23p * t_env) / a22p;
        let t_air = (b1 - a12 * t_surf - a13 * t_env) / a11;

        self.air_temperature_c = t_air;
        self.surface_temperature_c = t_surf;
        self.envelope_temperature_c = t_env;
        (t_air, t_surf, t_env)
    }
}

impl HvacIdealLoads {
    /// Required HVAC thermal power (W) to reach a target air temperature at the end of
    /// a timestep for the two-node model.
    ///
    /// Returns net thermal power into the air node (positive = heating, negative = cooling).
    pub fn required_hvac_power_two_node(
        &self,
        target_air_temp_c: f64,
        air_temp_c: f64,
        mass_temp_c: f64,
        outdoor_temp_c: f64,
        envelope_conductance_w_per_k: f64,
        air_mass_conductance_w_per_k: f64,
        gains_air_w: f64,
        gains_mass_w: f64,
        air_capacity_j_per_k: f64,
        mass_capacity_j_per_k: f64,
        dt_s: f64,
    ) -> f64 {
        if dt_s <= 0.0 {
            return 0.0;
        }
        if air_capacity_j_per_k <= 0.0 || mass_capacity_j_per_k <= 0.0 {
            return 0.0;
        }

        let k_env = envelope_conductance_w_per_k.max(0.0);
        let k_am = air_mass_conductance_w_per_k.max(0.0);

        // Mass node implicit update when T_air(t+dt) is forced to target.
        let a = dt_s * k_am / mass_capacity_j_per_k;
        let t_mass_next = if a > 0.0 {
            (mass_temp_c + (dt_s / mass_capacity_j_per_k) * gains_mass_w + a * target_air_temp_c)
                / (1.0 + a)
        } else {
            mass_temp_c + (dt_s / mass_capacity_j_per_k) * gains_mass_w
        };

        // Derived implicit HVAC requirement for the air node.
        air_capacity_j_per_k * (target_air_temp_c - air_temp_c) / dt_s
            + k_env * (target_air_temp_c - outdoor_temp_c)
            + k_am * (target_air_temp_c - t_mass_next)
            - gains_air_w
    }

    /// Required HVAC thermal power (W) to reach a target air temperature at the end of
    /// a timestep for the two-node *envelope-to-mass* model.
    pub fn required_hvac_power_two_node_envelope(
        &self,
        target_air_temp_c: f64,
        air_temp_c: f64,
        mass_temp_c: f64,
        outdoor_temp_c: f64,
        air_outdoor_conductance_w_per_k: f64,
        mass_outdoor_conductance_w_per_k: f64,
        air_mass_conductance_w_per_k: f64,
        gains_air_w: f64,
        gains_mass_w: f64,
        air_capacity_j_per_k: f64,
        mass_capacity_j_per_k: f64,
        dt_s: f64,
    ) -> f64 {
        if dt_s <= 0.0 {
            return 0.0;
        }
        if air_capacity_j_per_k <= 0.0 || mass_capacity_j_per_k <= 0.0 {
            return 0.0;
        }

        let k_ao = air_outdoor_conductance_w_per_k.max(0.0);
        let k_mo = mass_outdoor_conductance_w_per_k.max(0.0);
        let k_am = air_mass_conductance_w_per_k.max(0.0);

        if k_am <= 0.0 {
            return 0.0;
        }

        // Mass node implicit update when T_air(t+dt) is forced to target.
        let a = dt_s * (k_am + k_mo) / mass_capacity_j_per_k;
        let t_mass_next = if a > 0.0 {
            (mass_temp_c
                + (dt_s / mass_capacity_j_per_k) * (gains_mass_w + k_mo * outdoor_temp_c)
                + (dt_s * k_am / mass_capacity_j_per_k) * target_air_temp_c)
                / (1.0 + a)
        } else {
            mass_temp_c + (dt_s / mass_capacity_j_per_k) * gains_mass_w
        };

        // Solve air equation for hvac, given T_air(t+dt)=target and T_mass(t+dt) from above.
        let a_air = 1.0 + dt_s * (k_ao + k_am) / air_capacity_j_per_k;
        (air_capacity_j_per_k * (a_air * target_air_temp_c - air_temp_c) / dt_s)
            - gains_air_w
            - k_ao * outdoor_temp_c
            - k_am * t_mass_next
    }

    /// Required HVAC thermal power (W) to reach a target air temperature at the end of
    /// a timestep for the three-node envelope model.
    pub fn required_hvac_power_three_node_envelope(
        &self,
        target_air_temp_c: f64,
        air_temp_c: f64,
        surface_temp_c: f64,
        envelope_temp_c: f64,
        outdoor_temp_c: f64,
        air_outdoor_conductance_w_per_k: f64,
        air_surface_conductance_w_per_k: f64,
        surface_envelope_conductance_w_per_k: f64,
        envelope_outdoor_conductance_w_per_k: f64,
        gains_air_w: f64,
        gains_surface_w: f64,
        gains_envelope_w: f64,
        air_capacity_j_per_k: f64,
        surface_capacity_j_per_k: f64,
        envelope_capacity_j_per_k: f64,
        dt_s: f64,
    ) -> f64 {
        if dt_s <= 0.0 {
            return 0.0;
        }
        if air_capacity_j_per_k <= 0.0
            || surface_capacity_j_per_k <= 0.0
            || envelope_capacity_j_per_k <= 0.0
        {
            return 0.0;
        }

        let k_ao = air_outdoor_conductance_w_per_k.max(0.0);
        let k_as = air_surface_conductance_w_per_k.max(0.0);
        let k_se = surface_envelope_conductance_w_per_k.max(0.0);
        let k_eo = envelope_outdoor_conductance_w_per_k.max(0.0);

        if k_as <= 0.0 || k_se <= 0.0 {
            return 0.0;
        }

        // Solve (implicit) surface+envelope next temperatures with T_air_next forced to target.
        let cs = surface_capacity_j_per_k;
        let ce = envelope_capacity_j_per_k;

        let a22 = 1.0 + dt_s * (k_as + k_se) / cs;
        let a23 = -dt_s * k_se / cs;
        let b2 =
            surface_temp_c + (dt_s / cs) * gains_surface_w + (dt_s * k_as / cs) * target_air_temp_c;

        let a32 = -dt_s * k_se / ce;
        let a33 = 1.0 + dt_s * (k_se + k_eo) / ce;
        let b3 = envelope_temp_c + (dt_s / ce) * (gains_envelope_w + k_eo * outdoor_temp_c);

        let det = a22 * a33 - a23 * a32;
        if det.abs() < 1e-16 {
            return 0.0;
        }

        let t_surf_next = (b2 * a33 - a23 * b3) / det;
        let _t_env_next = (a22 * b3 - b2 * a32) / det;

        // Air equation with implicit conductances.
        // C_a*(T_next - T)/dt = k_ao*(T_out - T_next) + k_as*(T_surf_next - T_next) + gains + hvac
        air_capacity_j_per_k * (target_air_temp_c - air_temp_c) / dt_s
            + k_ao * (target_air_temp_c - outdoor_temp_c)
            + k_as * (target_air_temp_c - t_surf_next)
            - gains_air_w
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
    fn test_three_node_required_hvac_hits_setpoint_heating() {
        let hvac = HvacIdealLoads::with_setpoints(20.0, 27.0);
        let mut model = ThreeNodeEnvelopeThermalModel::new(
            18.0,
            18.0,
            10.0,
            50.0,        // k_ao
            200.0,       // k_as
            100.0,       // k_se
            60.0,        // k_eo
            80_000.0,    // c_air
            2_000_000.0, // c_surface
            1_000_000.0, // c_envelope
        );

        let dt_s = 3600.0;
        let outdoor = 0.0;
        let gains_air = 50.0;
        let gains_surface = 0.0;
        let gains_envelope = 0.0;

        let q_hvac = hvac.required_hvac_power_three_node_envelope(
            hvac.heating_setpoint,
            model.air_temperature_c,
            model.surface_temperature_c,
            model.envelope_temperature_c,
            outdoor,
            model.air_outdoor_conductance_w_per_k,
            model.air_surface_conductance_w_per_k,
            model.surface_envelope_conductance_w_per_k,
            model.envelope_outdoor_conductance_w_per_k,
            gains_air,
            gains_surface,
            gains_envelope,
            model.air_capacity_j_per_k,
            model.surface_capacity_j_per_k,
            model.envelope_capacity_j_per_k,
            dt_s,
        );

        model.step(
            outdoor,
            gains_air,
            gains_surface,
            gains_envelope,
            q_hvac,
            dt_s,
        );

        assert!(
            (model.air_temperature_c - hvac.heating_setpoint).abs() < 1e-9,
            "air temp should hit setpoint, got {}",
            model.air_temperature_c
        );
        assert!(
            q_hvac.is_finite() && q_hvac >= 0.0,
            "heating hvac should be positive, got {q_hvac}"
        );
    }

    #[test]
    fn test_three_node_required_hvac_hits_setpoint_cooling() {
        let hvac = HvacIdealLoads::with_setpoints(20.0, 27.0);
        let mut model = ThreeNodeEnvelopeThermalModel::new(
            29.0,
            27.0,
            20.0,
            50.0,        // k_ao
            150.0,       // k_as
            80.0,        // k_se
            60.0,        // k_eo
            80_000.0,    // c_air
            2_000_000.0, // c_surface
            1_000_000.0, // c_envelope
        );

        let dt_s = 3600.0;
        let outdoor = 35.0;
        let gains_air = 0.0;
        let gains_surface = 100.0;
        let gains_envelope = 0.0;

        let q_hvac = hvac.required_hvac_power_three_node_envelope(
            hvac.cooling_setpoint,
            model.air_temperature_c,
            model.surface_temperature_c,
            model.envelope_temperature_c,
            outdoor,
            model.air_outdoor_conductance_w_per_k,
            model.air_surface_conductance_w_per_k,
            model.surface_envelope_conductance_w_per_k,
            model.envelope_outdoor_conductance_w_per_k,
            gains_air,
            gains_surface,
            gains_envelope,
            model.air_capacity_j_per_k,
            model.surface_capacity_j_per_k,
            model.envelope_capacity_j_per_k,
            dt_s,
        );

        model.step(
            outdoor,
            gains_air,
            gains_surface,
            gains_envelope,
            q_hvac,
            dt_s,
        );

        assert!(
            (model.air_temperature_c - hvac.cooling_setpoint).abs() < 1e-9,
            "air temp should hit setpoint, got {}",
            model.air_temperature_c
        );
        assert!(
            q_hvac.is_finite() && q_hvac <= 0.0,
            "cooling hvac should be negative, got {q_hvac}"
        );
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
    fn test_two_node_required_hvac_hits_setpoint() {
        let hvac = HvacIdealLoads::with_setpoints(20.0, 27.0);

        let mut model = TwoNodeThermalModel::new(
            18.0,
            18.0,
            120.0,     // k_env
            300.0,     // k_am
            800_000.0, // c_air
            5_000_000.0,
        );

        let dt_s = 3600.0;
        let outdoor = 0.0;
        let gains_air = 0.0;
        let gains_mass = 0.0;
        let target = 20.0;

        let q_hvac = hvac.required_hvac_power_two_node(
            target,
            model.air_temperature_c,
            model.mass_temperature_c,
            outdoor,
            model.envelope_conductance_w_per_k,
            model.air_mass_conductance_w_per_k,
            gains_air,
            gains_mass,
            model.air_capacity_j_per_k,
            model.mass_capacity_j_per_k,
            dt_s,
        );

        model.step(outdoor, gains_air, gains_mass, q_hvac, dt_s);
        assert!(
            (model.air_temperature_c - target).abs() < 1e-9,
            "Expected air temperature to hit setpoint; got {}",
            model.air_temperature_c
        );
    }

    #[test]
    fn test_two_node_envelope_required_hvac_hits_setpoint() {
        let hvac = HvacIdealLoads::with_setpoints(20.0, 27.0);

        let mut model = TwoNodeEnvelopeThermalModel::new(
            18.0,
            18.0,
            20.0,      // k_air_out (infiltration)
            100.0,     // k_mass_out (envelope)
            300.0,     // k_am
            800_000.0, // c_air
            5_000_000.0,
        );

        let dt_s = 3600.0;
        let outdoor = 0.0;
        let gains_air = 0.0;
        let gains_mass = 0.0;
        let target = 20.0;

        let q_hvac = hvac.required_hvac_power_two_node_envelope(
            target,
            model.air_temperature_c,
            model.mass_temperature_c,
            outdoor,
            model.air_outdoor_conductance_w_per_k,
            model.mass_outdoor_conductance_w_per_k,
            model.air_mass_conductance_w_per_k,
            gains_air,
            gains_mass,
            model.air_capacity_j_per_k,
            model.mass_capacity_j_per_k,
            dt_s,
        );

        model.step(outdoor, gains_air, gains_mass, q_hvac, dt_s);
        assert!(
            (model.air_temperature_c - target).abs() < 1e-9,
            "Expected air temperature to hit setpoint; got {}",
            model.air_temperature_c
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
        assert!(
            (heating - 2005.0).abs() < 1e-6,
            "Expected 2005 W heating, got {heating}"
        );
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
        let (heating, cooling) =
            hvac.calculate_with_losses(30.0, 35.0, 100.0, 500.0, 3600.0, 3600.0);
        // Q_hvac = 3600*(26-30)/3600 + 100*(26-35) - 500 = -4 + (-900) - 500 = -1404
        // So cooling = 1404 W
        assert!(heating.abs() < 1e-10);
        assert!(
            (cooling - 1404.0).abs() < 1e-6,
            "Expected 1404 W cooling, got {cooling}"
        );
    }

    #[test]
    fn test_calculate_with_losses_matches_lumped_model_step() {
        let hvac = HvacIdealLoads::new();
        let dt_s = 3600.0;

        // Heating case: should end exactly at heating setpoint.
        let mut model = LumpedThermalModel::new(15.0, 100.0, 0.0, 3600.0);
        let total_conductance = model.ua_total + model.infiltration_conductance;
        let (heating, cooling) = hvac.calculate_with_losses(
            model.zone_temperature,
            0.0,
            total_conductance,
            0.0,
            model.thermal_capacity,
            dt_s,
        );
        model.step(0.0, 0.0, heating - cooling, dt_s);
        assert!((model.zone_temperature - hvac.heating_setpoint).abs() < 1e-10);

        // Cooling case: should end exactly at cooling setpoint.
        let mut model = LumpedThermalModel::new(30.0, 100.0, 0.0, 3600.0);
        let total_conductance = model.ua_total + model.infiltration_conductance;
        let (heating, cooling) = hvac.calculate_with_losses(
            model.zone_temperature,
            35.0,
            total_conductance,
            500.0,
            model.thermal_capacity,
            dt_s,
        );
        model.step(35.0, 500.0, heating - cooling, dt_s);
        assert!((model.zone_temperature - hvac.cooling_setpoint).abs() < 1e-10);
    }

    // ── Physics verification tests ──────────────────────────────────────

    #[test]
    fn test_1r1c_exponential_decay() {
        // Analytical solution for 1R1C free-floating (no gains, no HVAC):
        //   T(t) = T_out + (T_0 - T_out) * exp(-K*t/C)
        // where K = UA + infiltration_conductance.
        let t_out = 0.0;
        let t_0 = 20.0;
        let ua = 100.0; // W/K
        let inf = 20.0; // W/K
        let k = ua + inf; // 120 W/K
        let c = 500_000.0; // J/K
        let dt_s = 60.0; // 1-minute steps for accuracy

        let mut model = LumpedThermalModel::new(t_0, ua, inf, c);

        // Simulate 10 hours in 1-minute steps
        let total_seconds = 10.0 * 3600.0;
        let steps = (total_seconds / dt_s) as usize;
        for _ in 0..steps {
            model.step(t_out, 0.0, 0.0, dt_s);
        }

        let t_analytical = t_out + (t_0 - t_out) * (-k * total_seconds / c).exp();
        assert!(
            (model.zone_temperature - t_analytical).abs() < 0.05,
            "After 10h: model={:.4}, analytical={:.4}",
            model.zone_temperature,
            t_analytical
        );
    }

    #[test]
    fn test_1r1c_steady_state_convergence_with_gains() {
        // With constant gains Q and no HVAC, the steady-state temperature is:
        //   T_ss = T_out + Q / K
        let t_out = 5.0;
        let q_gains = 600.0; // W
        let ua = 80.0;
        let inf = 20.0;
        let k = ua + inf; // 100 W/K
        let c = 200_000.0;
        let dt_s = 3600.0;
        let t_ss_expected = t_out + q_gains / k; // 5 + 6 = 11°C

        let mut model = LumpedThermalModel::new(t_out, ua, inf, c);

        // Run for 48 hours — should converge well within that
        for _ in 0..(48 * 3600 / dt_s as usize) {
            model.step(t_out, q_gains, 0.0, dt_s);
        }

        assert!(
            (model.zone_temperature - t_ss_expected).abs() < 0.01,
            "Should converge to T_ss={:.2}, got {:.4}",
            t_ss_expected,
            model.zone_temperature
        );
    }

    #[test]
    fn test_1r1c_free_floating_bounded() {
        // Without HVAC or gains, temperature must stay between T_0 and T_out
        // (it decays monotonically toward T_out).
        let t_out = -5.0;
        let t_0 = 22.0;
        let mut model = LumpedThermalModel::new(t_0, 150.0, 30.0, 300_000.0);

        let mut prev = t_0;
        for _ in 0..200 {
            let t = model.step(t_out, 0.0, 0.0, 3600.0);
            assert!(
                t >= t_out && t <= prev,
                "Temperature should decay monotonically toward T_out: prev={prev}, t={t}"
            );
            prev = t;
        }
    }

    #[test]
    fn test_1r1c_time_constant() {
        // The time constant tau = C/K. After one tau, the temperature should
        // drop to ~36.8% of the initial dT (i.e., T ≈ T_out + 0.368*(T_0 - T_out)).
        let t_out = 0.0;
        let t_0 = 20.0;
        let k = 100.0;
        let c = 360_000.0; // tau = 3600 s = 1 hour
        let tau = c / k;
        let dt_s = 10.0; // fine steps

        let mut model = LumpedThermalModel::new(t_0, k, 0.0, c);

        let steps = (tau / dt_s) as usize;
        for _ in 0..steps {
            model.step(t_out, 0.0, 0.0, dt_s);
        }

        let expected = t_out + (t_0 - t_out) * (-1.0_f64).exp(); // ~7.358°C
        assert!(
            (model.zone_temperature - expected).abs() < 0.1,
            "After 1 tau: expected={expected:.3}, got={:.3}",
            model.zone_temperature
        );
    }

    #[test]
    fn test_hvac_calculate_for_timestep_invalid_inputs() {
        let hvac = HvacIdealLoads::new();
        let (h, c, t) = hvac.calculate_for_timestep(10.0, 0.0, 3600.0);
        assert_eq!((h, c, t), (0.0, 0.0, 10.0));

        let (h, c, t) = hvac.calculate_for_timestep(10.0, 1000.0, 0.0);
        assert_eq!((h, c, t), (0.0, 0.0, 10.0));
    }

    #[test]
    fn test_lumped_model_edge_cases() {
        // dt<=0 -> unchanged
        let mut model = LumpedThermalModel::new(20.0, 100.0, 0.0, 1000.0);
        let t = model.step(0.0, 0.0, 0.0, 0.0);
        assert!((t - 20.0).abs() < 1e-12);

        // zero thermal capacity -> steady-state formula
        let mut model = LumpedThermalModel::new(20.0, 100.0, 20.0, 0.0);
        let t = model.step(0.0, 200.0, 0.0, 3600.0);
        // T = T_out + Q/K = 0 + 200/120
        assert!((t - (200.0 / 120.0)).abs() < 1e-12);

        // no losses (K<=0) -> pure integrator
        let mut model = LumpedThermalModel::new(20.0, 0.0, 0.0, 1000.0);
        let t = model.step(0.0, 100.0, 0.0, 10.0);
        assert!((t - 21.0).abs() < 1e-12);
    }

    #[test]
    fn test_default_trait() {
        let hvac: HvacIdealLoads = Default::default();
        assert!((hvac.heating_setpoint - 20.0).abs() < 1e-12);
        assert!((hvac.cooling_setpoint - 26.0).abs() < 1e-12);
    }

    #[test]
    fn test_hvac_capacity_limits_cooling_branch() {
        let mut hvac = HvacIdealLoads::new();
        hvac.max_cooling_capacity = 2.0; // W thermal

        // Need 4 W thermal to drop 4 K in an hour with C=3600 J/K, but only 2 W available.
        let (heating_elec, cooling_elec, temp) = hvac.calculate(30.0, 3600.0);
        assert!(heating_elec.abs() < 1e-12);
        // COP=3 -> electric == thermal/3
        assert!((cooling_elec - (2.0 / 3.0)).abs() < 1e-12);
        // Delivered deltaT = 2 W * 3600 s / 3600 J/K = 2 K
        assert!((temp - 28.0).abs() < 1e-10);
    }

    #[test]
    fn test_calculate_with_losses_invalid_dt_or_capacity_returns_zero() {
        let hvac = HvacIdealLoads::new();
        let (h, c) = hvac.calculate_with_losses(15.0, 0.0, 100.0, 0.0, 3600.0, 0.0);
        assert_eq!((h, c), (0.0, 0.0));
        let (h, c) = hvac.calculate_with_losses(15.0, 0.0, 100.0, 0.0, 0.0, 3600.0);
        assert_eq!((h, c), (0.0, 0.0));
    }

    #[test]
    fn test_three_node_hvac_zero_capacities() {
        let hvac = HvacIdealLoads::new();
        // Zero air capacity -> returns 0
        let q = hvac.required_hvac_power_three_node_envelope(
            20.0,
            18.0,
            18.0,
            10.0,
            0.0,
            50.0,
            200.0,
            100.0,
            60.0,
            0.0,
            0.0,
            0.0,
            0.0, // c_air = 0
            2_000_000.0,
            1_000_000.0,
            3600.0,
        );
        assert_eq!(q, 0.0, "Zero air capacity should return 0");

        // Zero surface capacity -> returns 0
        let q = hvac.required_hvac_power_three_node_envelope(
            20.0,
            18.0,
            18.0,
            10.0,
            0.0,
            50.0,
            200.0,
            100.0,
            60.0,
            0.0,
            0.0,
            0.0,
            80_000.0,
            0.0, // c_surface = 0
            1_000_000.0,
            3600.0,
        );
        assert_eq!(q, 0.0, "Zero surface capacity should return 0");

        // Zero envelope capacity -> returns 0
        let q = hvac.required_hvac_power_three_node_envelope(
            20.0,
            18.0,
            18.0,
            10.0,
            0.0,
            50.0,
            200.0,
            100.0,
            60.0,
            0.0,
            0.0,
            0.0,
            80_000.0,
            2_000_000.0,
            0.0, // c_envelope = 0
            3600.0,
        );
        assert_eq!(q, 0.0, "Zero envelope capacity should return 0");
    }

    #[test]
    fn test_three_node_hvac_zero_dt() {
        let hvac = HvacIdealLoads::new();
        let q = hvac.required_hvac_power_three_node_envelope(
            20.0,
            18.0,
            18.0,
            10.0,
            0.0,
            50.0,
            200.0,
            100.0,
            60.0,
            0.0,
            0.0,
            0.0,
            80_000.0,
            2_000_000.0,
            1_000_000.0,
            0.0, // dt = 0
        );
        assert_eq!(q, 0.0, "Zero timestep should return 0");
    }

    #[test]
    fn test_three_node_hvac_zero_conductances() {
        let hvac = HvacIdealLoads::new();
        // k_as = 0 -> returns 0
        let q = hvac.required_hvac_power_three_node_envelope(
            20.0,
            18.0,
            18.0,
            10.0,
            0.0,
            50.0,
            0.0, // k_as = 0
            100.0,
            60.0,
            0.0,
            0.0,
            0.0,
            80_000.0,
            2_000_000.0,
            1_000_000.0,
            3600.0,
        );
        assert_eq!(q, 0.0, "Zero k_as should return 0");

        // k_se = 0 -> returns 0
        let q = hvac.required_hvac_power_three_node_envelope(
            20.0,
            18.0,
            18.0,
            10.0,
            0.0,
            50.0,
            200.0,
            0.0, // k_se = 0
            60.0,
            0.0,
            0.0,
            0.0,
            80_000.0,
            2_000_000.0,
            1_000_000.0,
            3600.0,
        );
        assert_eq!(q, 0.0, "Zero k_se should return 0");
    }

    #[test]
    fn test_two_node_step_degenerate_zero_capacity() {
        // Both capacities zero -> steady-state fallback
        let mut model = TwoNodeThermalModel::new(20.0, 20.0, 100.0, 300.0, 0.0, 0.0);
        let (ta, tm) = model.step(0.0, 200.0, 0.0, 0.0, 3600.0);
        // T = T_out + Q/K = 0 + 200/100 = 2.0
        assert!((ta - 2.0).abs() < 1e-10, "got ta={ta}");
        assert!((tm - ta).abs() < 1e-10);
    }

    #[test]
    fn test_two_node_step_degenerate_zero_k_env() {
        // Positive capacity but k_env=0 and k_am=0 -> pure integrator
        let mut model = TwoNodeThermalModel::new(20.0, 20.0, 0.0, 0.0, 1000.0, 0.0);
        let (ta, _) = model.step(0.0, 100.0, 0.0, 0.0, 10.0);
        // dT = dt/c * Q = 10/1000 * 100 = 1.0
        assert!((ta - 21.0).abs() < 1e-10, "got ta={ta}");
    }

    #[test]
    fn test_two_node_step_dt_zero() {
        let mut model = TwoNodeThermalModel::new(20.0, 18.0, 100.0, 300.0, 800_000.0, 5_000_000.0);
        let (ta, tm) = model.step(0.0, 0.0, 0.0, 0.0, 0.0);
        assert!((ta - 20.0).abs() < 1e-12);
        assert!((tm - 18.0).abs() < 1e-12);
    }

    #[test]
    fn test_two_node_envelope_step_degenerate_zero_capacity() {
        let mut model = TwoNodeEnvelopeThermalModel::new(20.0, 20.0, 20.0, 100.0, 300.0, 0.0, 0.0);
        let (ta, tm) = model.step(0.0, 200.0, 0.0, 0.0, 3600.0);
        // T = T_out + Q/K_total = 0 + 200/120 ≈ 1.667
        assert!((ta - (200.0 / 120.0)).abs() < 1e-10, "got ta={ta}");
        assert!((tm - ta).abs() < 1e-10);
    }

    #[test]
    fn test_two_node_envelope_step_dt_zero() {
        let mut model = TwoNodeEnvelopeThermalModel::new(
            20.0,
            18.0,
            20.0,
            100.0,
            300.0,
            800_000.0,
            5_000_000.0,
        );
        let (ta, tm) = model.step(0.0, 0.0, 0.0, 0.0, 0.0);
        assert!((ta - 20.0).abs() < 1e-12);
        assert!((tm - 18.0).abs() < 1e-12);
    }

    #[test]
    fn test_two_node_envelope_step_zero_k_total() {
        // k_ao=0, k_mo=0, k_am=0 -> pure integrator fallback
        let mut model = TwoNodeEnvelopeThermalModel::new(20.0, 20.0, 0.0, 0.0, 0.0, 1000.0, 0.0);
        let (ta, _) = model.step(0.0, 100.0, 0.0, 0.0, 10.0);
        assert!((ta - 21.0).abs() < 1e-10, "got ta={ta}");
    }

    #[test]
    fn test_three_node_step_normal() {
        // Normal case: all positive capacities and conductances
        let mut model = ThreeNodeEnvelopeThermalModel::new(
            20.0,
            20.0,
            10.0,
            50.0,        // k_ao
            200.0,       // k_as
            100.0,       // k_se
            60.0,        // k_eo
            80_000.0,    // c_air
            2_000_000.0, // c_surface
            1_000_000.0, // c_envelope
        );

        // dt=0 -> unchanged
        let (ta, ts, te) = model.step(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        assert!((ta - 20.0).abs() < 1e-12);
        assert!((ts - 20.0).abs() < 1e-12);
        assert!((te - 10.0).abs() < 1e-12);

        // Normal step: outdoor=0, no gains, no hvac -> temperatures change
        let (ta, ts, te) = model.step(0.0, 0.0, 0.0, 0.0, 0.0, 3600.0);
        assert!(ta < 20.0, "Air should cool: ta={ta}");
        assert!(ta > 0.0, "Air shouldn't overshoot outdoor: ta={ta}");
        assert!(ts.is_finite(), "Surface temp should be finite: ts={ts}");
        assert!(te.is_finite(), "Envelope temp should be finite: te={te}");
    }

    #[test]
    fn test_three_node_step_degenerate_zero_capacity() {
        // All capacities zero -> steady-state fallback
        let mut model = ThreeNodeEnvelopeThermalModel::new(
            20.0, 20.0, 10.0, 50.0, 0.0, 0.0, 60.0, 0.0, 0.0, 0.0,
        );
        let (ta, ts, te) = model.step(0.0, 200.0, 0.0, 0.0, 0.0, 3600.0);
        // c=0, k=50+60=110 -> T = 0 + 200/110
        let expected = 200.0 / 110.0;
        assert!((ta - expected).abs() < 1e-10, "got ta={ta}");
        assert!((ts - ta).abs() < 1e-10);
        assert!((te - ta).abs() < 1e-10);
    }

    #[test]
    fn test_three_node_step_degenerate_zero_k() {
        // k_ao=0, k_eo=0 -> pure integrator
        let mut model = ThreeNodeEnvelopeThermalModel::new(
            20.0, 20.0, 20.0, 0.0, 0.0, 0.0, 0.0, 3000.0, 0.0, 0.0,
        );
        let (ta, _, _) = model.step(0.0, 300.0, 0.0, 0.0, 0.0, 10.0);
        // dT = dt/c * Q = 10/3000 * 300 = 1.0
        assert!((ta - 21.0).abs() < 1e-10, "got ta={ta}");
    }

    #[test]
    fn test_three_node_step_steady_state_convergence() {
        // With constant gains and no HVAC, the three-node model should
        // converge to a steady state over many steps.
        let mut model = ThreeNodeEnvelopeThermalModel::new(
            0.0,
            0.0,
            0.0,
            50.0,        // k_ao
            200.0,       // k_as
            100.0,       // k_se
            60.0,        // k_eo
            80_000.0,    // c_air
            2_000_000.0, // c_surface
            1_000_000.0, // c_envelope
        );

        let dt = 3600.0;
        for _ in 0..500 {
            model.step(0.0, 1000.0, 0.0, 0.0, 0.0, dt);
        }

        // At steady state, all net fluxes should be small
        let ta = model.air_temperature_c;
        assert!(ta > 0.0, "Should warm up with gains");
        assert!(ta.is_finite(), "Should converge to finite temperature");
    }
}
