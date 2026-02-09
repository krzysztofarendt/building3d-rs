use std::collections::HashMap;

use crate::Building;
use crate::sim::lighting::solar::SolarPosition;
use crate::sim::materials::MaterialLibrary;

/// Default Solar Heat Gain Coefficient for glazing.
const DEFAULT_SHGC: f64 = 0.6;

/// Configuration for physics-based solar gain calculation using EPW weather data.
///
/// Computes solar gains directly from DNI/DHI radiation data combined with sun
/// geometry and window orientation.
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
    pub(crate) fn resolve_shgc_with_materials(
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
#[derive(Debug, Clone, Copy)]
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

/// Computes solar gains per zone UID (W) from EPW radiation data and sun geometry.
///
/// This is the zone-resolved counterpart of [`compute_solar_gains`]. It is useful for
/// multi-zone thermal simulations that need a reasonable default distribution of solar gains.
pub fn compute_solar_gains_per_zone(
    building: &Building,
    params: &SolarHourParams,
    config: &SolarGainConfig,
) -> HashMap<crate::UID, f64> {
    compute_solar_gains_per_zone_with_materials(building, params, config, None)
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

/// Like [`compute_solar_gains_per_zone`], but also accepts a material library for
/// `is_glazing`-based surface identification.
pub fn compute_solar_gains_per_zone_with_materials(
    building: &Building,
    params: &SolarHourParams,
    config: &SolarGainConfig,
    material_library: Option<&MaterialLibrary>,
) -> HashMap<crate::UID, f64> {
    let solar_pos = SolarPosition::calculate(
        params.latitude,
        params.longitude,
        params.day_of_year,
        params.hour,
    );
    if !solar_pos.is_above_horizon() {
        return compute_diffuse_only_per_zone(
            building,
            params.diffuse_horizontal_irradiance,
            config,
            material_library,
        );
    }

    let sun_dir = solar_pos.to_direction();
    let mut gains: HashMap<crate::UID, f64> = HashMap::new();

    for zone in building.zones() {
        let mut zone_gains = 0.0;
        for solid in zone.solids() {
            for wall in solid.walls() {
                for polygon in wall.polygons() {
                    let path = format!(
                        "{}/{}/{}/{}",
                        zone.name, solid.name, wall.name, polygon.name
                    );
                    let Some(shgc) = config.resolve_shgc_with_materials(&path, material_library)
                    else {
                        continue;
                    };

                    let area = polygon.area();
                    if area <= 0.0 {
                        continue;
                    }

                    let normal = polygon.vn;
                    let cos_incidence = sun_dir.dot(&normal).max(0.0);
                    let q_direct = params.direct_normal_irradiance * cos_incidence * area * shgc;

                    let sky_view = 0.5 * (1.0 + normal.dz.max(0.0));
                    let q_diffuse = params.diffuse_horizontal_irradiance * sky_view * area * shgc;

                    zone_gains += q_direct + q_diffuse;
                }
            }
        }
        if zone_gains != 0.0 {
            gains.insert(zone.uid.clone(), zone_gains);
        }
    }

    gains
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

fn compute_diffuse_only_per_zone(
    building: &Building,
    diffuse_horizontal_irradiance: f64,
    config: &SolarGainConfig,
    material_library: Option<&MaterialLibrary>,
) -> HashMap<crate::UID, f64> {
    let mut gains: HashMap<crate::UID, f64> = HashMap::new();
    for zone in building.zones() {
        let mut zone_total = 0.0;
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
                        zone_total += diffuse_horizontal_irradiance * sky_view * area * shgc;
                    }
                }
            }
        }
        if zone_total != 0.0 {
            gains.insert(zone.uid.clone(), zone_total);
        }
    }
    gains
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::materials::{Material, MaterialLibrary};
    use crate::{Solid, Zone};

    #[test]
    fn test_compute_solar_gains_south_facing_window() {
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

    #[test]
    fn test_compute_solar_gains_per_zone_sums_to_total() {
        let s0 = Solid::from_box(1.0, 1.0, 1.0, None, "s0").unwrap();
        let s1 = Solid::from_box(1.0, 1.0, 1.0, Some((2.0, 0.0, 0.0)), "s1").unwrap();
        let z0 = Zone::new("z0", vec![s0]).unwrap();
        let z1 = Zone::new("z1", vec![s1]).unwrap();
        let building = Building::new("b", vec![z0, z1]).unwrap();

        let params = SolarHourParams {
            direct_normal_irradiance: 800.0,
            diffuse_horizontal_irradiance: 100.0,
            day_of_year: 81, // ~equinox
            hour: 9.0,       // morning sun (east-ish)
            latitude: 0.0,
            longitude: 0.0,
        };

        let mut cfg = SolarGainConfig::new();
        cfg.shgc.insert("z0/s0/wall_1/poly_1".to_string(), 1.0);
        cfg.shgc.insert("z1/s1/wall_1/poly_1".to_string(), 1.0);

        let total = compute_solar_gains(&building, &params, &cfg);
        let by_zone = compute_solar_gains_per_zone(&building, &params, &cfg);
        let sum: f64 = by_zone.values().sum();

        assert!(
            (sum - total).abs() < 1e-9,
            "Per-zone sum should match total: sum={sum}, total={total}"
        );
        assert_eq!(by_zone.len(), 2);
    }

    #[test]
    fn test_default_traits() {
        let cfg: SolarGainConfig = Default::default();
        assert!(cfg.default_shgc > 0.0);
    }

    #[test]
    fn test_resolve_shgc_with_materials_uses_is_glazing_flag() {
        let mut lib = MaterialLibrary::new();
        lib.add(Material::new("glass").with_glazing());
        lib.assign("/", "glass");

        let mut cfg = SolarGainConfig::new();
        cfg.default_shgc = 0.75;

        let shgc = cfg.resolve_shgc_with_materials("zone/room/wall/window_0", Some(&lib));
        assert_eq!(shgc, Some(0.75));
    }
}
