use std::collections::HashMap;

use crate::sim::lighting::result::LightingResult;

// Converts lighting simulation results to thermal solar heat gains.
//
// Luminous flux (lm) → Radiometric flux (W) using luminous efficacy.
// For solar radiation, typical luminous efficacy ≈ 120 lm/W.
//
// Solar heat gain for a surface = transmitted_flux / luminous_efficacy * SHGC
// where SHGC = Solar Heat Gain Coefficient (fraction of solar energy entering the zone).

/// Default luminous efficacy for solar radiation (lm/W).
const DEFAULT_LUMINOUS_EFFICACY: f64 = 120.0;

/// Default Solar Heat Gain Coefficient for glazing.
const DEFAULT_SHGC: f64 = 0.6;

/// Configuration for solar bridge conversion.
#[derive(Debug, Clone)]
pub struct SolarBridgeConfig {
    /// Luminous efficacy in lm/W.
    pub luminous_efficacy: f64,
    /// Solar heat gain coefficient per polygon path.
    /// Only surfaces with entries are considered as glazing.
    pub shgc: HashMap<String, f64>,
    /// Default SHGC for surfaces matched by pattern.
    pub default_shgc: f64,
    /// Patterns identifying glazing surfaces.
    pub glazing_patterns: Vec<String>,
}

impl SolarBridgeConfig {
    pub fn new() -> Self {
        Self {
            luminous_efficacy: DEFAULT_LUMINOUS_EFFICACY,
            shgc: HashMap::new(),
            default_shgc: DEFAULT_SHGC,
            glazing_patterns: vec![
                "window".to_string(),
                "glazing".to_string(),
                "glass".to_string(),
            ],
        }
    }

    /// Returns the SHGC for a polygon path, or None if it's not glazing.
    fn resolve_shgc(&self, path: &str) -> Option<f64> {
        // Exact match
        if let Some(&shgc) = self.shgc.get(path) {
            return Some(shgc);
        }
        // Pattern match
        for pattern in &self.glazing_patterns {
            if path.contains(pattern.as_str()) {
                return Some(self.default_shgc);
            }
        }
        None
    }
}

impl Default for SolarBridgeConfig {
    fn default() -> Self {
        Self::new()
    }
}

/// Converts lighting simulation illuminance into total solar heat gain (W).
///
/// Only considers surfaces identified as glazing (by pattern or explicit SHGC).
/// For each glazing surface:
///   Q_solar = incident_flux_total / luminous_efficacy * SHGC
pub fn lighting_to_solar_gains(
    lighting_result: &LightingResult,
    config: &SolarBridgeConfig,
) -> f64 {
    if config.luminous_efficacy <= 0.0 {
        return 0.0;
    }

    let mut total_gains = 0.0;

    for (path, flux) in &lighting_result.incident_flux {
        if let Some(shgc) = config.resolve_shgc(path) {
            // Total luminous flux hitting this surface
            let total_flux = flux[0] + flux[1] + flux[2];
            // Convert lm to W and apply SHGC
            let q = total_flux / config.luminous_efficacy * shgc;
            total_gains += q;
        }
    }

    total_gains
}

/// Returns per-surface solar gains in W.
pub fn lighting_to_solar_gains_per_surface(
    lighting_result: &LightingResult,
    config: &SolarBridgeConfig,
) -> HashMap<String, f64> {
    if config.luminous_efficacy <= 0.0 {
        return HashMap::new();
    }

    let mut gains = HashMap::new();

    for (path, flux) in &lighting_result.incident_flux {
        if let Some(shgc) = config.resolve_shgc(path) {
            let total_flux = flux[0] + flux[1] + flux[2];
            let q = total_flux / config.luminous_efficacy * shgc;
            gains.insert(path.clone(), q);
        }
    }

    gains
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_basic_solar_gains() {
        let mut result = LightingResult::new();
        // Simulate light hitting a window surface
        result.record_hit("zone/room/window/glass", [1200.0, 1200.0, 1200.0]);

        let config = SolarBridgeConfig::new();
        let gains = lighting_to_solar_gains(&result, &config);

        // Total flux = 3600 lm
        // Q = 3600 / 120 * 0.6 = 18.0 W
        assert!(
            (gains - 18.0).abs() < 0.1,
            "Expected ~18 W solar gain, got {gains}"
        );
    }

    #[test]
    fn test_non_glazing_surface_ignored() {
        let mut result = LightingResult::new();
        result.record_hit("zone/room/wall/concrete", [1000.0, 1000.0, 1000.0]);

        let config = SolarBridgeConfig::new();
        let gains = lighting_to_solar_gains(&result, &config);

        assert!(
            gains.abs() < 1e-10,
            "Non-glazing surfaces should not contribute solar gains"
        );
    }

    #[test]
    fn test_custom_shgc() {
        let mut result = LightingResult::new();
        result.record_hit("zone/room/special/pane", [600.0, 600.0, 600.0]);

        let mut config = SolarBridgeConfig::new();
        config
            .shgc
            .insert("zone/room/special/pane".to_string(), 0.3);

        let gains = lighting_to_solar_gains(&result, &config);
        // Total flux = 1800 lm
        // Q = 1800 / 120 * 0.3 = 4.5 W
        assert!(
            (gains - 4.5).abs() < 0.1,
            "Expected ~4.5 W with SHGC=0.3, got {gains}"
        );
    }

    #[test]
    fn test_per_surface_gains() {
        let mut result = LightingResult::new();
        result.record_hit("zone/room/window/w1", [120.0, 120.0, 120.0]);
        result.record_hit("zone/room/window/w2", [240.0, 240.0, 240.0]);

        let config = SolarBridgeConfig::new();
        let per_surface = lighting_to_solar_gains_per_surface(&result, &config);

        assert_eq!(per_surface.len(), 2);
        // w1: 360/120*0.6 = 1.8
        // w2: 720/120*0.6 = 3.6
        assert!((per_surface["zone/room/window/w1"] - 1.8).abs() < 0.1);
        assert!((per_surface["zone/room/window/w2"] - 3.6).abs() < 0.1);
    }

    #[test]
    fn test_zero_efficacy_guard() {
        let mut result = LightingResult::new();
        result.record_hit("zone/room/window/w1", [120.0, 120.0, 120.0]);

        let mut config = SolarBridgeConfig::new();
        config.luminous_efficacy = 0.0;

        let gains = lighting_to_solar_gains(&result, &config);
        assert!(gains == 0.0);

        let per_surface = lighting_to_solar_gains_per_surface(&result, &config);
        assert!(per_surface.is_empty());
    }
}
