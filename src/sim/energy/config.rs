use std::collections::HashMap;

use super::construction::WallConstruction;

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
}
