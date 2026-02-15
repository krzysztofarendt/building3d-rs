use std::collections::HashMap;

use super::construction::WallConstruction;
use super::convection::{ExteriorConvectionModel, InteriorConvectionModel};
use crate::UID;
use crate::sim::materials::MaterialLibrary;

/// Boundary type for an internal mass surface.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InternalMassBoundary {
    /// Both sides of the mass exchange heat with the zone air (e.g., an interior partition).
    TwoSided,
    /// Only one side exchanges heat with the zone air; the other side is adiabatic.
    ///
    /// This is a simple way to represent a massive floor slab sitting on high-R insulation.
    OneSidedAdiabatic,
}

/// A non-geometric internal thermal mass surface assigned to one or more zones.
///
/// These surfaces are modeled as 1D FVM slabs with convective coupling to zone air.
/// They can receive a portion of transmitted solar and radiant internal gains when
/// `use_surface_aware_solar_distribution` is enabled.
#[derive(Debug, Clone)]
pub struct InternalMassSurface {
    /// Human-readable name for diagnostics.
    pub name: String,
    /// Zone name/path substring pattern (matched against `Zone::name`).
    pub zone_path_pattern: String,
    /// Exposed face area (m²) for the slab.
    pub area_m2: f64,
    /// Layer stack (outside → inside). The *zone-exposed* face is treated as the "inside".
    pub construction: WallConstruction,
    /// Boundary type (one-sided vs two-sided).
    pub boundary: InternalMassBoundary,
    /// Cosine of the surface tilt from horizontal (polygon normal dz equivalent).
    ///
    /// - `0.0` (default) → vertical partition
    /// - `-1.0` → floor (interior side faces up)
    /// - `+1.0` → ceiling (interior side faces down)
    ///
    /// Used by dynamic convection models (TARP) to select the appropriate
    /// natural convection correlation.
    pub cos_tilt: f64,
}

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
    /// If true, approximate interior longwave exchange using a simple "radiant node"
    /// (MRT-style) per zone.
    ///
    /// This affects only the FVM-based interior boundary conditions by splitting the
    /// interior film coefficient into:
    /// - a convective portion to zone air, and
    /// - a radiative portion to a zone radiant temperature estimate.
    ///
    /// The radiant temperature is estimated from the current interior surface
    /// temperatures (area-weighted), and does not add or remove energy from the zone
    /// (it only redistributes among surfaces).
    pub use_interior_radiative_exchange: bool,
    /// Fraction (0..1) of the interior film coefficient assigned to the radiative path
    /// when [`Self::use_interior_radiative_exchange`] is enabled.
    pub interior_radiation_fraction: f64,
    /// Fraction (0..1) of **solar** gains applied to the mass node when the two-node
    /// model is enabled. The remainder is applied to the air node.
    pub solar_gains_to_mass_fraction: f64,
    /// If true, use a surface-aware policy for distributing **transmitted** solar gains
    /// in the two-node model.
    ///
    /// When enabled, transmitted shortwave is treated primarily as an **interior surface**
    /// heat input (mass node) rather than as an instantaneous air gain, which improves
    /// lag/peak behavior in heavyweight cases.
    pub use_surface_aware_solar_distribution: bool,
    /// Fraction (0..1) of transmitted solar gains applied directly to the **air** node
    /// when [`Self::use_surface_aware_solar_distribution`] is enabled.
    ///
    /// The remainder is applied to the mass node (interior surfaces).
    pub transmitted_solar_to_air_fraction: f64,
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
    /// Optional 3-node (air + interior surface + envelope) model enable knob.
    ///
    /// When `> 0.0` and when the 2R2C model is otherwise enabled, the transient solver
    /// uses a coarse 3-node model:
    /// - air node
    /// - interior surface/mass node (receives most transmitted solar)
    /// - envelope node (receives exterior absorbed solar / ground correction)
    ///
    /// Interpretation: fraction (0..1) of the total **mass-side** capacity assigned to the
    /// envelope node. The remainder is assigned to the interior surface node.
    pub three_node_envelope_mass_fraction: f64,
    /// If true, replace steady `U*A*ΔT` exterior conduction for eligible opaque
    /// exterior surfaces with per-surface 1D FVM wall solvers.
    ///
    /// Eligibility (current policy):
    /// - exterior surface,
    /// - has a resolved [`WallConstruction`] (layer stack),
    /// - does **not** have an explicit U-value override (UID or path-pattern),
    /// - is not treated as glazing (material `is_glazing` or name heuristics),
    /// - is not ground-coupled (when `ground_temperature_c` is set).
    ///
    /// This is intended for transient simulations and step-based pipelines.
    pub use_fvm_walls: bool,
    /// Optional internal thermal mass surfaces (non-geometric).
    ///
    /// These are modeled as 1D FVM slabs coupled to zone air, and can receive a
    /// portion of transmitted solar and radiant internal gains when
    /// `use_surface_aware_solar_distribution` is enabled.
    pub internal_mass_surfaces: Vec<InternalMassSurface>,
    /// Policy for computing inter-zone partition conductance from two assigned U-values.
    pub interzone_u_value_policy: InterZoneUValuePolicy,
    /// If true, distribute transmitted solar to FVM wall interior faces (area-proportional
    /// with internal mass slabs). If false (default), transmitted solar goes only to
    /// internal mass slabs, which avoids an unrealistic heat-loss path through wall
    /// insulation to outdoors.
    pub distribute_transmitted_solar_to_fvm_walls: bool,
    /// If true, use beam vs diffuse solar split for interior distribution.
    ///
    /// Beam (direct) solar is sent 100% to internal mass slabs (floor), while
    /// diffuse solar is distributed area-proportionally across all interior
    /// surfaces (mass + FVM walls). This approximates the geometric reality
    /// that beam solar through windows mostly hits the floor, not walls.
    pub use_beam_solar_distribution: bool,
    /// If true AND `distribute_transmitted_solar_to_fvm_walls` is true, the wall
    /// portion of interior surface sources (transmitted solar + radiant gains) is
    /// redirected to zone air instead of being injected into FVM wall domains.
    ///
    /// This includes wall area in the solar distribution denominator (diluting the
    /// flux on mass surfaces), while avoiding the heat-loss path through wall
    /// insulation to outdoors. The wall mass still participates in thermal buffering
    /// via its normal convective coupling to zone air.
    pub fvm_wall_solar_to_air: bool,
    /// Interior surface convection model.
    ///
    /// - `Fixed(3.0)` (default): legacy constant coefficient.
    /// - `Tarp`: temperature- and tilt-dependent natural convection (Walton/TARP).
    pub interior_convection_model: InteriorConvectionModel,
    /// Exterior surface convection model.
    ///
    /// - `Fixed` (default): coefficient from ISO R_se.
    /// - `Doe2`: DOE-2 simplified combined natural + wind-forced convection.
    pub exterior_convection_model: ExteriorConvectionModel,
    /// If true, compute per-surface view factors and use them for interior
    /// longwave radiation exchange instead of the simplified area-weighted MRT.
    ///
    /// This replaces the `use_interior_radiative_exchange` f_rad split with a
    /// proper energy-conserving model: uniform `h_rad` + per-surface MRT from
    /// Monte Carlo view factors.
    pub use_view_factor_radiation: bool,
    /// Number of cosine-weighted rays per surface for Monte Carlo view factor
    /// computation. More rays → more accurate F_ij but slower init.
    pub view_factor_rays_per_surface: usize,
    /// Interior surface emissivity for linearized radiative exchange.
    ///
    /// Used to compute `h_rad = 4 * eps * sigma * T_mean^3`.
    /// Default: 0.9 (typical for building interior surfaces).
    pub interior_emissivity: f64,
    /// Enable iterative surface heat balance (simultaneous solve of surface
    /// temperatures, MRT, and zone air temperature within each substep).
    ///
    /// When enabled, the substep loop wraps FVM wall steps + air model in an
    /// outer iteration that converges surface temperatures and MRT together.
    /// This eliminates the one-substep lag in MRT and allows convective-only
    /// air gain (radiation stays between surfaces, netting to zero by reciprocity).
    pub use_iterative_surface_balance: bool,
    /// Maximum iterations per substep for the iterative surface balance.
    pub surface_balance_max_iterations: usize,
    /// Convergence tolerance for surface temperatures [°C].
    pub surface_balance_tolerance_c: f64,
    /// If true, assemble all wall FVM cells, surface nodes, and air nodes into a
    /// single global matrix and solve simultaneously each substep.
    ///
    /// This enables radiation coupling to be embedded directly in the matrix
    /// (no iteration needed) and produces self-consistent surface temperatures,
    /// air temperatures, and radiative exchange within each timestep.
    ///
    /// Requires `use_fvm_walls = true`. When enabled, the sequential
    /// per-wall Thomas solve + separate air model is replaced by a dense
    /// global solve.
    pub use_global_fvm_solve: bool,
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
            use_interior_radiative_exchange: false,
            interior_radiation_fraction: 0.6,
            solar_gains_to_mass_fraction: 0.0,
            use_surface_aware_solar_distribution: false,
            transmitted_solar_to_air_fraction: 0.0,
            internal_gains_to_mass_fraction: 0.0,
            two_node_envelope_to_mass: false,
            three_node_envelope_mass_fraction: 0.0,
            use_fvm_walls: true,
            internal_mass_surfaces: vec![],
            interzone_u_value_policy: InterZoneUValuePolicy::Mean,
            distribute_transmitted_solar_to_fvm_walls: false,
            use_beam_solar_distribution: false,
            fvm_wall_solar_to_air: false,
            interior_convection_model: InteriorConvectionModel::default(),
            exterior_convection_model: ExteriorConvectionModel::default(),
            use_view_factor_radiation: false,
            view_factor_rays_per_surface: 10_000,
            interior_emissivity: 0.9,
            use_iterative_surface_balance: false,
            surface_balance_max_iterations: 4,
            surface_balance_tolerance_c: 0.1,
            use_global_fvm_solve: false,
        }
    }

    /// Resolves a wall construction for a path (if any).
    ///
    /// Resolution is deterministic and mirrors U-value resolution for constructions:
    /// - exact path match wins
    /// - otherwise, the longest substring match wins (ties: lexicographic key order)
    pub fn resolve_construction(&self, path: &str) -> Option<&WallConstruction> {
        if let Some(c) = self.constructions.get(path) {
            return Some(c);
        }
        self.best_matching_construction(path)
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

    /// Returns true if this surface has an explicit U-value override (by UID or path/pattern).
    ///
    /// This is used by higher-fidelity wall models to decide whether a surface should
    /// follow an explicit manufacturer-style U-value rather than a layered construction.
    pub fn has_u_value_override_for_surface(&self, polygon_uid: &UID, path: &str) -> bool {
        if self
            .u_value_overrides_by_polygon_uid
            .contains_key(polygon_uid)
        {
            return true;
        }
        if self.u_value_overrides_by_path_pattern.contains_key(path) {
            return true;
        }
        self.best_matching_u_value_override(path).is_some()
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
        assert!(config.use_fvm_walls);
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

    #[test]
    fn test_resolve_envelope_capacity_precedence() {
        let uid = crate::UID::new();
        let mut config = ThermalConfig::new();
        let fallback = 50_000.0;

        // Default fallback
        let c = config.resolve_envelope_capacity_j_per_m2_k(
            Some(&uid),
            "zone/solid/wall/poly",
            fallback,
        );
        assert!((c - fallback).abs() < 1e-12);

        // Construction match wins over default
        let wall = insulated_wall();
        // insulated_wall() has layers -> thermal_capacity() returns their sum
        let expected_cap = wall.thermal_capacity();
        config.constructions.insert("wall".to_string(), wall);
        let c = config.resolve_envelope_capacity_j_per_m2_k(None, "zone/solid/wall/poly", fallback);
        assert!(
            (c - expected_cap.max(0.0)).abs() < 1e-6,
            "Construction should win: got {c}, expected {expected_cap}"
        );

        // UID override wins over construction
        config
            .envelope_capacity_overrides_j_per_m2_k_by_polygon_uid
            .insert(uid.clone(), 99_999.0);
        let c = config.resolve_envelope_capacity_j_per_m2_k(
            Some(&uid),
            "zone/solid/wall/poly",
            fallback,
        );
        assert!(
            (c - 99_999.0).abs() < 1e-12,
            "UID override should win: got {c}"
        );
    }
}
