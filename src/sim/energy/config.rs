use std::collections::HashMap;

use super::construction::WallConstruction;

/// Policy for computing inter-zone partition conductance from two assigned U-values.
///
/// When two adjacent solids belong to different zones, the interface is represented
/// by two facing polygons. It is ambiguous whether these represent:
/// - the same physical construction (duplicated surfaces), or
/// - two constructions in series (each side contributes its own resistance).
///
/// This policy makes the behavior explicit and configurable.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InterZoneUValuePolicy {
    /// Use the mean of the two U-values (default; preserves historical behavior).
    Mean,
    /// Treat U-values as two resistances in series: `U_eq = 1 / (1/U1 + 1/U2)`.
    Series,
}

/// Configuration for thermal simulation.
#[derive(Debug, Clone)]
pub struct ThermalConfig {
    /// Wall constructions assigned to polygon paths (pattern matching).
    /// Key is a path pattern, value is a WallConstruction.
    pub constructions: HashMap<String, WallConstruction>,
    /// Default U-value for surfaces without assigned construction (W/(m^2*K)).
    pub default_u_value: f64,
    /// Outdoor temperature in °C (for steady-state calculation).
    pub outdoor_temperature: f64,
    /// Indoor setpoint temperature in °C.
    pub indoor_temperature: f64,
    /// Infiltration rate in air changes per hour (ACH).
    pub infiltration_ach: f64,
    /// Internal heat gains in W (from people, equipment, lighting).
    pub internal_gains: f64,
    /// Solar heat gains in W (simplified, can be overridden by lighting sim).
    pub solar_gains: f64,
    /// Thermal capacity per zone air volume in J/(m^3*K).
    ///
    /// Used by transient zone-air models to estimate total zone thermal capacity:
    /// `C_zone = V_zone * thermal_capacity_j_per_m3_k`.
    pub thermal_capacity_j_per_m3_k: f64,
    /// Policy for computing inter-zone partition conductance from two assigned U-values.
    pub interzone_u_value_policy: InterZoneUValuePolicy,
}

impl ThermalConfig {
    pub fn new() -> Self {
        Self {
            constructions: HashMap::new(),
            default_u_value: 2.0,
            outdoor_temperature: 0.0,
            indoor_temperature: 20.0,
            infiltration_ach: 0.5,
            internal_gains: 0.0,
            solar_gains: 0.0,
            thermal_capacity_j_per_m3_k: 50_000.0,
            interzone_u_value_policy: InterZoneUValuePolicy::Mean,
        }
    }

    /// Resolves U-value for a polygon path.
    pub fn resolve_u_value(&self, path: &str) -> f64 {
        // Check for exact match first
        if let Some(construction) = self.constructions.get(path) {
            return construction.u_value();
        }
        // Check for partial match (path starts with pattern)
        for (pattern, construction) in &self.constructions {
            if path.starts_with(pattern) || path.ends_with(pattern) || path.contains(pattern) {
                return construction.u_value();
            }
        }
        self.default_u_value
    }

    /// Computes an equivalent inter-zone conductance (W/K) for a partition.
    pub fn interzone_conductance_w_per_k(&self, u1: f64, u2: f64, area_m2: f64) -> f64 {
        if area_m2 <= 0.0 {
            return 0.0;
        }
        let u_eq = match self.interzone_u_value_policy {
            InterZoneUValuePolicy::Mean => {
                if u1.is_finite() && u2.is_finite() {
                    0.5 * (u1 + u2)
                } else if u1.is_finite() {
                    u1
                } else if u2.is_finite() {
                    u2
                } else {
                    return 0.0;
                }
            }
            InterZoneUValuePolicy::Series => {
                if !(u1.is_finite() && u2.is_finite() && u1 > 0.0 && u2 > 0.0) {
                    return 0.0;
                }
                1.0 / (1.0 / u1 + 1.0 / u2)
            }
        };

        (u_eq * area_m2).max(0.0)
    }
}

impl Default for ThermalConfig {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::energy::construction::insulated_wall;

    #[test]
    fn test_config_defaults() {
        let config = ThermalConfig::new();
        assert!((config.default_u_value - 2.0).abs() < 1e-10);
        assert!((config.indoor_temperature - 20.0).abs() < 1e-10);
        assert!((config.thermal_capacity_j_per_m3_k - 50_000.0).abs() < 1e-10);
        assert_eq!(config.interzone_u_value_policy, InterZoneUValuePolicy::Mean);
    }

    #[test]
    fn test_resolve_u_value() {
        let mut config = ThermalConfig::new();
        config
            .constructions
            .insert("wall".to_string(), insulated_wall());

        // Exact match
        let u = config.resolve_u_value("wall");
        assert!(u < 0.5, "Should use insulated wall U-value");

        // Pattern match
        let u = config.resolve_u_value("zone/solid/wall/poly");
        assert!(u < 0.5, "Should match pattern containing 'wall'");

        // No match - default
        let u = config.resolve_u_value("zone/solid/roof/poly");
        assert!((u - 2.0).abs() < 1e-10, "Should use default U-value");
    }

    #[test]
    fn test_interzone_conductance_policies() {
        let mut config = ThermalConfig::new();

        // Mean: (2 + 2)/2 * 1 = 2
        config.interzone_u_value_policy = InterZoneUValuePolicy::Mean;
        let k_mean = config.interzone_conductance_w_per_k(2.0, 2.0, 1.0);
        assert!((k_mean - 2.0).abs() < 1e-12);

        // Series: 1 / (1/2 + 1/2) * 1 = 1
        config.interzone_u_value_policy = InterZoneUValuePolicy::Series;
        let k_series = config.interzone_conductance_w_per_k(2.0, 2.0, 1.0);
        assert!((k_series - 1.0).abs() < 1e-12);
    }
}
