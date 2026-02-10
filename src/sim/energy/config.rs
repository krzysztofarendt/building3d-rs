use std::collections::HashMap;

use super::construction::WallConstruction;
use crate::UID;
use crate::sim::materials::MaterialLibrary;

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
    /// Optional shared material library used for thermal U-value/capacity resolution.
    ///
    /// When set, and when a surface is assigned a material that contains
    /// `ThermalMaterial`, that material becomes the default source of truth for
    /// U-values/capacities (unless overridden by exact per-surface overrides).
    pub material_library: Option<MaterialLibrary>,
    /// Exact U-value overrides keyed by polygon `UID` (W/(m^2*K)).
    pub u_value_overrides_by_polygon_uid: HashMap<UID, f64>,
    /// U-value overrides keyed by polygon path or path substring pattern.
    ///
    /// Resolution is deterministic and mirrors `constructions`:
    /// - exact path match wins
    /// - otherwise, the longest substring match wins (ties: lexicographic key order)
    ///
    /// This is useful for cases where a surface should be specified directly by a
    /// manufacturer-style U-value (e.g., glazing) rather than by a layered construction.
    pub u_value_overrides_by_path_pattern: HashMap<String, f64>,
    /// Exact envelope capacity overrides keyed by polygon `UID` (J/(m^2*K)).
    pub envelope_capacity_overrides_j_per_m2_k_by_polygon_uid: HashMap<UID, f64>,
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
    /// Optional ground temperature boundary for ground-coupled surfaces (°C).
    ///
    /// When set, surfaces matched by [`Self::ground_surface_patterns`] contribute an
    /// extra conductance term to a "ground" boundary rather than outdoor air. The
    /// transient solvers approximate this by adding a per-hour gain:
    /// `Q_ground = UA_ground * (T_ground - T_outdoor)`.
    pub ground_temperature_c: Option<f64>,
    /// Path substring patterns identifying ground-coupled exterior surfaces.
    ///
    /// Default: `["floor"]`.
    pub ground_surface_patterns: Vec<String>,
    /// Optional two-node (air + mass) thermal model enable knob.
    ///
    /// - `0.0` (default): use the historical 1R1C zone model.
    /// - `> 0.0`: enable a 2R2C model by splitting the total thermal capacity
    ///   into an air node and a mass node.
    ///
    /// Interpretation: fraction of total capacity assigned to the **mass** node.
    pub two_node_mass_fraction: f64,
    /// Interior heat transfer coefficient used to couple air ↔ mass (W/(m²·K)).
    ///
    /// This is a coarse aggregate that lumps convection+radiation exchange between
    /// zone air and interior surfaces.
    pub interior_heat_transfer_coeff_w_per_m2_k: f64,
    /// Fraction (0..1) of **solar** gains applied to the mass node when the two-node
    /// model is enabled. The remainder is applied to the air node.
    pub solar_gains_to_mass_fraction: f64,
    /// Fraction (0..1) of **internal** gains applied to the mass node when the two-node
    /// model is enabled. The remainder is applied to the air node.
    pub internal_gains_to_mass_fraction: f64,
    /// If true, attach the envelope conductance to the **mass** node (instead of the air node)
    /// when the 2R2C model is enabled.
    ///
    /// Rationale: in heavyweight buildings (e.g. BESTEST 900), most of the thermal mass is
    /// in the envelope and should be directly coupled to outdoors. The legacy 2R2C model
    /// couples the air node directly to outdoors, which can distort lag/peak behavior.
    ///
    /// When enabled:
    /// - `UA` (conduction) is connected mass ↔ outdoors
    /// - infiltration remains air ↔ outdoors
    pub two_node_envelope_to_mass: bool,
    /// Policy for computing inter-zone partition conductance from two assigned U-values.
    pub interzone_u_value_policy: InterZoneUValuePolicy,
}

impl ThermalConfig {
    pub fn new() -> Self {
        Self {
            constructions: HashMap::new(),
            material_library: None,
            u_value_overrides_by_polygon_uid: HashMap::new(),
            u_value_overrides_by_path_pattern: HashMap::new(),
            envelope_capacity_overrides_j_per_m2_k_by_polygon_uid: HashMap::new(),
            default_u_value: 2.0,
            outdoor_temperature: 0.0,
            indoor_temperature: 20.0,
            infiltration_ach: 0.5,
            internal_gains: 0.0,
            solar_gains: 0.0,
            thermal_capacity_j_per_m3_k: 50_000.0,
            ground_temperature_c: None,
            ground_surface_patterns: vec!["floor".to_string()],
            two_node_mass_fraction: 0.0,
            interior_heat_transfer_coeff_w_per_m2_k: 3.0,
            solar_gains_to_mass_fraction: 0.0,
            internal_gains_to_mass_fraction: 0.0,
            two_node_envelope_to_mass: false,
            interzone_u_value_policy: InterZoneUValuePolicy::Mean,
        }
    }

    /// Resolves U-value for a polygon path.
    pub fn resolve_u_value(&self, path: &str) -> f64 {
        self.resolve_u_value_with_uid(None, path)
    }

    /// Resolves U-value for a specific polygon surface, enabling UID-keyed overrides.
    pub fn resolve_u_value_for_surface(&self, polygon_uid: &UID, path: &str) -> f64 {
        self.resolve_u_value_with_uid(Some(polygon_uid), path)
    }

    fn resolve_u_value_with_uid(&self, polygon_uid: Option<&UID>, path: &str) -> f64 {
        // 1) Exact override by UID
        if let Some(uid) = polygon_uid
            && let Some(&u) = self.u_value_overrides_by_polygon_uid.get(uid)
        {
            return u;
        }

        // 2) U-value override by path/pattern (manufacturer-style override)
        if let Some(&u) = self.u_value_overrides_by_path_pattern.get(path) {
            return u;
        }
        if let Some(u) = self.best_matching_u_value_override(path) {
            return u;
        }

        // 3) Exact override by path via explicit construction entry
        if let Some(construction) = self.constructions.get(path) {
            return construction.u_value();
        }

        // 4) MaterialLibrary assignment (if provided and contains ThermalMaterial)
        if let Some(lib) = self.material_library.as_ref()
            && let Some(mat) = lib.lookup(path)
            && let Some(t) = mat.thermal.as_ref()
        {
            return t.u_value;
        }

        // 5) Deterministic best-match construction pattern (fallback)
        if let Some(construction) = self.best_matching_construction(path) {
            return construction.u_value();
        }

        // 6) Default
        self.default_u_value
    }

    /// Resolves envelope thermal capacity per unit area (J/(m²·K)) for a surface.
    ///
    /// Precedence:
    /// 1) exact override by polygon UID
    /// 2) exact path match in `constructions`
    /// 3) material library ThermalMaterial capacity
    /// 4) best-matching construction pattern
    /// 5) provided fallback default
    pub fn resolve_envelope_capacity_j_per_m2_k(
        &self,
        polygon_uid: Option<&UID>,
        path: &str,
        default_capacity_j_per_m2_k: f64,
    ) -> f64 {
        if let Some(uid) = polygon_uid
            && let Some(&c) = self
                .envelope_capacity_overrides_j_per_m2_k_by_polygon_uid
                .get(uid)
        {
            return c.max(0.0);
        }

        if let Some(construction) = self.constructions.get(path) {
            return construction.thermal_capacity().max(0.0);
        }

        if let Some(lib) = self.material_library.as_ref()
            && let Some(mat) = lib.lookup(path)
            && let Some(t) = mat.thermal.as_ref()
        {
            return t.thermal_capacity.max(0.0);
        }

        if let Some(construction) = self.best_matching_construction(path) {
            return construction.thermal_capacity().max(0.0);
        }

        default_capacity_j_per_m2_k.max(0.0)
    }

    fn best_matching_construction(&self, path: &str) -> Option<&WallConstruction> {
        let mut best: Option<(&str, &WallConstruction)> = None;
        for (pattern, construction) in &self.constructions {
            if pattern == path {
                continue;
            }
            if !path.contains(pattern.as_str()) {
                continue;
            }
            match best {
                None => best = Some((pattern.as_str(), construction)),
                Some((best_pat, _)) => {
                    let better = pattern.len() > best_pat.len()
                        || (pattern.len() == best_pat.len() && pattern.as_str() < best_pat);
                    if better {
                        best = Some((pattern.as_str(), construction));
                    }
                }
            }
        }
        best.map(|(_, c)| c)
    }

    fn best_matching_u_value_override(&self, path: &str) -> Option<f64> {
        let mut best: Option<(&str, f64)> = None;
        for (pattern, &u) in &self.u_value_overrides_by_path_pattern {
            if pattern == path {
                continue;
            }
            if !path.contains(pattern.as_str()) {
                continue;
            }
            match best {
                None => best = Some((pattern.as_str(), u)),
                Some((best_pat, _)) => {
                    let better = pattern.len() > best_pat.len()
                        || (pattern.len() == best_pat.len() && pattern.as_str() < best_pat);
                    if better {
                        best = Some((pattern.as_str(), u));
                    }
                }
            }
        }
        best.map(|(_, u)| u)
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
        assert!(config.ground_temperature_c.is_none());
        assert_eq!(config.ground_surface_patterns, vec!["floor".to_string()]);
        assert!((config.two_node_mass_fraction - 0.0).abs() < 1e-12);
        assert!((config.interior_heat_transfer_coeff_w_per_m2_k - 3.0).abs() < 1e-12);
        assert!((config.solar_gains_to_mass_fraction - 0.0).abs() < 1e-12);
        assert!((config.internal_gains_to_mass_fraction - 0.0).abs() < 1e-12);
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
    fn test_resolve_u_value_override_by_path_pattern() {
        let mut cfg = ThermalConfig::new();
        cfg.constructions
            .insert("window".to_string(), insulated_wall());

        cfg.u_value_overrides_by_path_pattern
            .insert("window".to_string(), 3.0);

        let u = cfg.resolve_u_value("zone/solid/wall/window_1");
        assert!(
            (u - 3.0).abs() < 1e-12,
            "Expected U override to win, got {u}"
        );
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

    #[test]
    fn test_interzone_conductance_edge_cases() {
        let mut config = ThermalConfig::new();

        // Non-positive area -> 0.
        assert_eq!(config.interzone_conductance_w_per_k(2.0, 2.0, 0.0), 0.0);
        assert_eq!(config.interzone_conductance_w_per_k(2.0, 2.0, -1.0), 0.0);

        // Mean policy: prefer finite U-values.
        config.interzone_u_value_policy = InterZoneUValuePolicy::Mean;
        let k = config.interzone_conductance_w_per_k(f64::NAN, 3.0, 2.0);
        assert!((k - 6.0).abs() < 1e-12);
        let k = config.interzone_conductance_w_per_k(4.0, f64::INFINITY, 2.0);
        assert!((k - 8.0).abs() < 1e-12);
        assert_eq!(
            config.interzone_conductance_w_per_k(f64::NAN, f64::NAN, 2.0),
            0.0
        );

        // Series policy: invalid/non-positive U-values -> 0.
        config.interzone_u_value_policy = InterZoneUValuePolicy::Series;
        assert_eq!(config.interzone_conductance_w_per_k(-1.0, 2.0, 1.0), 0.0);
        assert_eq!(config.interzone_conductance_w_per_k(0.0, 2.0, 1.0), 0.0);
        assert_eq!(
            config.interzone_conductance_w_per_k(2.0, f64::NAN, 1.0),
            0.0
        );
    }

    #[test]
    fn test_default_trait() {
        let cfg: ThermalConfig = Default::default();
        assert!((cfg.default_u_value - 2.0).abs() < 1e-12);
    }
}
