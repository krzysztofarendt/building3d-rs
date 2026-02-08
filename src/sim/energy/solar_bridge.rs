use std::collections::HashMap;

use crate::Building;
use crate::sim::index::SurfaceIndex;
use crate::sim::lighting::result::LightingResult;
use crate::sim::lighting::solar::SolarPosition;
use crate::sim::materials::MaterialLibrary;

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
    surface_index: &SurfaceIndex,
    config: &SolarBridgeConfig,
) -> f64 {
    if config.luminous_efficacy <= 0.0 {
        return 0.0;
    }

    let mut total_gains = 0.0;

    for (polygon_uid, flux) in &lighting_result.incident_flux {
        let Some(path) = surface_index.path_by_polygon_uid(polygon_uid) else {
            continue;
        };
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
    surface_index: &SurfaceIndex,
    config: &SolarBridgeConfig,
) -> HashMap<String, f64> {
    if config.luminous_efficacy <= 0.0 {
        return HashMap::new();
    }

    let mut gains = HashMap::new();

    for (polygon_uid, flux) in &lighting_result.incident_flux {
        let Some(path) = surface_index.path_by_polygon_uid(polygon_uid) else {
            continue;
        };
        if let Some(shgc) = config.resolve_shgc(path) {
            let total_flux = flux[0] + flux[1] + flux[2];
            let q = total_flux / config.luminous_efficacy * shgc;
            gains.insert(path.to_string(), q);
        }
    }

    gains
}

/// Configuration for physics-based solar gain calculation using EPW weather data.
///
/// Instead of converting lighting lumens to watts, this computes solar gains
/// directly from DNI/DHI radiation data combined with sun geometry and window
/// orientation.
#[derive(Debug, Clone)]
pub struct SolarGainConfig {
    /// Solar heat gain coefficient per polygon path.
    /// Only surfaces with entries are considered as glazing.
    pub shgc: HashMap<String, f64>,
    /// Default SHGC for surfaces matched by pattern.
    pub default_shgc: f64,
    /// Patterns identifying glazing surfaces (substring match on polygon path).
    pub glazing_patterns: Vec<String>,
}

impl SolarGainConfig {
    pub fn new() -> Self {
        Self {
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
    ///
    /// Resolution order:
    /// 1. Exact match in `shgc` map
    /// 2. Material library `is_glazing` flag (if `material_library` provided)
    /// 3. Substring pattern match (fallback)
    fn resolve_shgc_with_materials(
        &self,
        path: &str,
        material_library: Option<&MaterialLibrary>,
    ) -> Option<f64> {
        // 1. Exact match in per-surface SHGC map
        if let Some(&shgc) = self.shgc.get(path) {
            return Some(shgc);
        }
        // 2. Check material library is_glazing flag
        if let Some(lib) = material_library
            && let Some(mat) = lib.lookup(path)
            && mat.is_glazing
        {
            return Some(self.default_shgc);
        }
        // 3. Fallback: substring pattern match
        for pattern in &self.glazing_patterns {
            if path.contains(pattern.as_str()) {
                return Some(self.default_shgc);
            }
        }
        None
    }
}

impl Default for SolarGainConfig {
    fn default() -> Self {
        Self::new()
    }
}

/// Parameters describing the solar conditions for a single hour.
pub struct SolarHourParams {
    /// Direct normal irradiance (W/m^2).
    pub direct_normal_irradiance: f64,
    /// Diffuse horizontal irradiance (W/m^2).
    pub diffuse_horizontal_irradiance: f64,
    /// Day of year (1-365).
    pub day_of_year: u16,
    /// Solar hour (0-24).
    pub hour: f64,
    /// Site latitude in degrees.
    pub latitude: f64,
    /// Site longitude in degrees.
    pub longitude: f64,
}

/// Computes solar heat gains (W) from EPW weather data and sun geometry.
///
/// For each glazing polygon in the building:
///   Q = DNI * cos(incidence) * area * SHGC + DHI * sky_view * area * SHGC
///
/// where:
///   - DNI = direct normal irradiance from weather record (W/m^2)
///   - DHI = diffuse horizontal irradiance from weather record (W/m^2)
///   - cos(incidence) = max(0, sun_direction . polygon_normal)
///   - sky_view = 0.5 * (1 + max(0, normal_z)) isotropic sky view factor
///
/// Glazing surfaces are identified by (in order of priority):
/// 1. Explicit per-surface SHGC entries in the config
/// 2. `is_glazing` flag on the material (if `material_library` is provided)
/// 3. Substring pattern matching on polygon path names (fallback)
pub fn compute_solar_gains(
    building: &Building,
    params: &SolarHourParams,
    config: &SolarGainConfig,
) -> f64 {
    compute_solar_gains_with_materials(building, params, config, None)
}

/// Like [`compute_solar_gains`], but also accepts a material library for
/// `is_glazing`-based surface identification.
pub fn compute_solar_gains_with_materials(
    building: &Building,
    params: &SolarHourParams,
    config: &SolarGainConfig,
    material_library: Option<&MaterialLibrary>,
) -> f64 {
    let solar_pos = SolarPosition::calculate(
        params.latitude,
        params.longitude,
        params.day_of_year,
        params.hour,
    );
    if !solar_pos.is_above_horizon() {
        // Sun below horizon: only diffuse contribution
        return compute_diffuse_only(
            building,
            params.diffuse_horizontal_irradiance,
            config,
            material_library,
        );
    }

    let sun_dir = solar_pos.to_direction();
    let mut total_gains = 0.0;

    for zone in building.zones() {
        for solid in zone.solids() {
            for wall in solid.walls() {
                for polygon in wall.polygons() {
                    let path = format!(
                        "{}/{}/{}/{}",
                        zone.name, solid.name, wall.name, polygon.name
                    );
                    if let Some(shgc) = config.resolve_shgc_with_materials(&path, material_library)
                    {
                        let normal = polygon.vn;
                        let area = polygon.area();

                        // Direct component: DNI * cos(incidence_angle) * area * SHGC
                        // cos(incidence) = sun_direction . surface_normal
                        // Only positive values (sun facing the surface)
                        let cos_incidence = sun_dir.dot(&normal).max(0.0);
                        let q_direct =
                            params.direct_normal_irradiance * cos_incidence * area * shgc;

                        // Diffuse component: DHI * sky_view_factor * area * SHGC
                        // Isotropic sky model: view factor = 0.5 for vertical, 1.0 for horizontal up
                        let sky_view = 0.5 * (1.0 + normal.dz.max(0.0));
                        let q_diffuse =
                            params.diffuse_horizontal_irradiance * sky_view * area * shgc;

                        total_gains += q_direct + q_diffuse;
                    }
                }
            }
        }
    }

    total_gains
}

/// Computes diffuse-only solar gains when the sun is below the horizon.
fn compute_diffuse_only(
    building: &Building,
    diffuse_horizontal_irradiance: f64,
    config: &SolarGainConfig,
    material_library: Option<&MaterialLibrary>,
) -> f64 {
    let mut total = 0.0;
    for zone in building.zones() {
        for solid in zone.solids() {
            for wall in solid.walls() {
                for polygon in wall.polygons() {
                    let path = format!(
                        "{}/{}/{}/{}",
                        zone.name, solid.name, wall.name, polygon.name
                    );
                    if let Some(shgc) = config.resolve_shgc_with_materials(&path, material_library)
                    {
                        let area = polygon.area();
                        let normal = polygon.vn;
                        let sky_view = 0.5 * (1.0 + normal.dz.max(0.0));
                        total += diffuse_horizontal_irradiance * sky_view * area * shgc;
                    }
                }
            }
        }
    }
    total
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Point, Polygon, Solid, Wall, Zone};

    fn make_building_with_surfaces(wall_and_polygon_names: Vec<(&str, &str)>) -> Building {
        let mut walls = Vec::new();
        for (wall_name, polygon_name) in wall_and_polygon_names {
            let poly = Polygon::new(
                polygon_name,
                vec![
                    Point::new(0.0, 0.0, 0.0),
                    Point::new(1.0, 0.0, 0.0),
                    Point::new(1.0, 1.0, 0.0),
                    Point::new(0.0, 1.0, 0.0),
                ],
                None,
            )
            .unwrap();
            walls.push(Wall::new(wall_name, vec![poly]).unwrap());
        }
        let solid = Solid::new("room", walls).unwrap();
        let zone = Zone::new("zone", vec![solid]).unwrap();
        Building::new("b", vec![zone]).unwrap()
    }

    #[test]
    fn test_basic_solar_gains() {
        let building = make_building_with_surfaces(vec![("window", "glass")]);
        let index = SurfaceIndex::new(&building);
        let polygon_uid = index
            .polygon_uid_by_path("zone/room/window/glass")
            .unwrap()
            .clone();

        let mut result = LightingResult::new();
        result
            .incident_flux
            .insert(polygon_uid, [1200.0, 1200.0, 1200.0]);

        let config = SolarBridgeConfig::new();
        let gains = lighting_to_solar_gains(&result, &index, &config);

        // Total flux = 3600 lm
        // Q = 3600 / 120 * 0.6 = 18.0 W
        assert!(
            (gains - 18.0).abs() < 0.1,
            "Expected ~18 W solar gain, got {gains}"
        );
    }

    #[test]
    fn test_non_glazing_surface_ignored() {
        let building = make_building_with_surfaces(vec![("wall", "concrete")]);
        let index = SurfaceIndex::new(&building);
        let polygon_uid = index
            .polygon_uid_by_path("zone/room/wall/concrete")
            .unwrap()
            .clone();

        let mut result = LightingResult::new();
        result
            .incident_flux
            .insert(polygon_uid, [1000.0, 1000.0, 1000.0]);

        let config = SolarBridgeConfig::new();
        let gains = lighting_to_solar_gains(&result, &index, &config);

        assert!(
            gains.abs() < 1e-10,
            "Non-glazing surfaces should not contribute solar gains"
        );
    }

    #[test]
    fn test_custom_shgc() {
        let building = make_building_with_surfaces(vec![("special", "pane")]);
        let index = SurfaceIndex::new(&building);
        let polygon_uid = index
            .polygon_uid_by_path("zone/room/special/pane")
            .unwrap()
            .clone();

        let mut result = LightingResult::new();
        result
            .incident_flux
            .insert(polygon_uid, [600.0, 600.0, 600.0]);

        let mut config = SolarBridgeConfig::new();
        config
            .shgc
            .insert("zone/room/special/pane".to_string(), 0.3);

        let gains = lighting_to_solar_gains(&result, &index, &config);
        // Total flux = 1800 lm
        // Q = 1800 / 120 * 0.3 = 4.5 W
        assert!(
            (gains - 4.5).abs() < 0.1,
            "Expected ~4.5 W with SHGC=0.3, got {gains}"
        );
    }

    #[test]
    fn test_per_surface_gains() {
        let building = make_building_with_surfaces(vec![("window", "w1"), ("window2", "w2")]);
        let index = SurfaceIndex::new(&building);
        let uid1 = index
            .polygon_uid_by_path("zone/room/window/w1")
            .unwrap()
            .clone();
        let uid2 = index
            .polygon_uid_by_path("zone/room/window2/w2")
            .unwrap()
            .clone();

        let mut result = LightingResult::new();
        result.incident_flux.insert(uid1, [120.0, 120.0, 120.0]);
        result.incident_flux.insert(uid2, [240.0, 240.0, 240.0]);

        let config = SolarBridgeConfig::new();
        let per_surface = lighting_to_solar_gains_per_surface(&result, &index, &config);

        assert_eq!(per_surface.len(), 2);
        // w1: 360/120*0.6 = 1.8
        // w2: 720/120*0.6 = 3.6
        assert!((per_surface["zone/room/window/w1"] - 1.8).abs() < 0.1);
        assert!((per_surface["zone/room/window2/w2"] - 3.6).abs() < 0.1);
    }

    #[test]
    fn test_zero_efficacy_guard() {
        let building = make_building_with_surfaces(vec![("window", "w1")]);
        let index = SurfaceIndex::new(&building);
        let polygon_uid = index
            .polygon_uid_by_path("zone/room/window/w1")
            .unwrap()
            .clone();

        let mut result = LightingResult::new();
        result
            .incident_flux
            .insert(polygon_uid, [120.0, 120.0, 120.0]);

        let mut config = SolarBridgeConfig::new();
        config.luminous_efficacy = 0.0;

        let gains = lighting_to_solar_gains(&result, &index, &config);
        assert!(gains == 0.0);

        let per_surface = lighting_to_solar_gains_per_surface(&result, &index, &config);
        assert!(per_surface.is_empty());
    }

    #[test]
    fn test_compute_solar_gains_south_facing_window() {
        use crate::{Solid, Zone};

        // Create a box with a "window" named polygon
        // We need a building with a glazing surface
        let s = Solid::from_box(5.0, 5.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let config = SolarGainConfig::new();

        // At solar noon in summer (day 172), latitude 45N
        let params = SolarHourParams {
            direct_normal_irradiance: 500.0,
            diffuse_horizontal_irradiance: 200.0,
            day_of_year: 172,
            hour: 12.0,
            latitude: 45.0,
            longitude: 0.0,
        };
        let gains = compute_solar_gains(&building, &params, &config);

        // No surfaces have "window" in their name in a basic box,
        // so gains should be zero
        assert!(
            gains.abs() < 1e-10,
            "No glazing surfaces in basic box, gains should be 0, got {gains}"
        );
    }

    #[test]
    fn test_compute_solar_gains_zero_radiation() {
        use crate::{Solid, Zone};

        let s = Solid::from_box(5.0, 5.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let config = SolarGainConfig::new();

        let params = SolarHourParams {
            direct_normal_irradiance: 0.0,
            diffuse_horizontal_irradiance: 0.0,
            day_of_year: 172,
            hour: 12.0,
            latitude: 45.0,
            longitude: 0.0,
        };
        let gains = compute_solar_gains(&building, &params, &config);

        assert!(gains.abs() < 1e-10, "Zero radiation should give zero gains");
    }
}
