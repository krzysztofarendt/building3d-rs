use std::collections::HashMap;

/// Result of a thermal simulation.
#[derive(Debug, Clone)]
pub struct ThermalResult {
    /// Heat loss through each surface in W.
    /// Key is polygon path.
    pub surface_heat_loss: HashMap<String, f64>,
    /// Total transmission heat loss in W (through surfaces).
    pub transmission_loss: f64,
    /// Infiltration heat loss in W.
    pub infiltration_loss: f64,
    /// Total heat gains in W (internal + solar).
    pub total_gains: f64,
    /// Net heating demand in W (positive = needs heating).
    pub heating_demand: f64,
    /// Net cooling demand in W (positive = needs cooling).
    pub cooling_demand: f64,
    /// Per-zone results.
    pub zone_results: HashMap<String, ZoneResult>,
}

/// Per-zone thermal result.
#[derive(Debug, Clone)]
pub struct ZoneResult {
    pub zone_name: String,
    /// Zone volume in m^3.
    pub volume: f64,
    /// Total envelope area in m^2.
    pub envelope_area: f64,
    /// Transmission loss for this zone in W.
    pub transmission_loss: f64,
    /// Infiltration loss for this zone in W.
    pub infiltration_loss: f64,
    /// Net demand in W.
    pub net_demand: f64,
}

impl ThermalResult {
    pub fn new() -> Self {
        Self {
            surface_heat_loss: HashMap::new(),
            transmission_loss: 0.0,
            infiltration_loss: 0.0,
            total_gains: 0.0,
            heating_demand: 0.0,
            cooling_demand: 0.0,
            zone_results: HashMap::new(),
        }
    }
}

impl Default for ThermalResult {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_result_creation() {
        let result = ThermalResult::new();
        assert!(result.surface_heat_loss.is_empty());
        assert!((result.transmission_loss - 0.0).abs() < 1e-10);
    }
}
