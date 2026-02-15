use std::collections::HashMap;

use crate::sim::lighting::solar::SolarPosition;
use crate::sim::materials::MaterialLibrary;
use crate::{Building, Point, UID, Vector};

use super::shading::{FinGeometry, OverhangGeometry};

/// Default Solar Heat Gain Coefficient for glazing.
const DEFAULT_SHGC: f64 = 0.6;
const DEFAULT_OPAQUE_ABSORPTANCE: f64 = 0.7;
const DEFAULT_EXTERIOR_HEAT_TRANSFER_COEFF_W_PER_M2_K: f64 = 17.0;
const DEFAULT_GROUND_REFLECTANCE: f64 = 0.2;
const DEFAULT_INC_ANGLE_MODIFIER_A: f64 = 0.1;

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
    /// If true, also model absorbed shortwave on exterior *opaque* surfaces as an
    /// additional thermal gain coupled through the envelope (sol-air approximation).
    pub include_exterior_opaque_absorption: bool,
    /// Default solar absorptance for opaque surfaces (0..1).
    pub default_opaque_absorptance: f64,
    /// Exterior combined heat transfer coefficient `h_out` (W/(m²·K)) used for the
    /// sol-air coupling term. Typical range: ~15–25 W/(m²·K).
    pub exterior_heat_transfer_coeff_w_per_m2_k: f64,
    /// If true, add a simple ground-reflected component to incident irradiance.
    ///
    /// Approximation: `I_ground = GHI * rho_g * ground_view`, where
    /// `ground_view = 0.5 * (1 - normal_z)`.
    pub include_ground_reflection: bool,
    /// Ground shortwave reflectance (albedo), 0..1.
    pub ground_reflectance: f64,
    /// If true, apply a simple incidence-angle modifier (IAM) to glazing gains.
    ///
    /// Approximation (ASHRAE-style):
    /// `IAM = 1 - a * (1/cos(theta) - 1)` clamped to [0, 1].
    pub include_incidence_angle_modifier: bool,
    /// IAM shape coefficient `a` (typical: ~0.05–0.2).
    pub incidence_angle_modifier_a: f64,
    /// If true, include exterior longwave exchange to sky/ground in the
    /// sol-air coupling term for opaque surfaces.
    pub include_exterior_longwave_exchange: bool,
    /// Effective longwave emissivity of exterior opaque surfaces (0..1).
    pub exterior_opaque_emissivity: f64,
    /// Effective longwave emissivity of the ground (0..1).
    pub ground_emissivity: f64,
    /// If true, compute `h_out` from wind speed (instead of a constant).
    pub use_wind_speed_for_h_out: bool,
    /// Wind-model base term for `h_out` in W/(m²·K).
    pub h_out_base_w_per_m2_k: f64,
    /// Wind-model slope for `h_out` in W/(m²·K) per (m/s).
    pub h_out_wind_coeff_w_per_m2_k_per_m_s: f64,
    /// Optional tilt scaling applied to `h_out`: `h_out *= 1 + tilt_scale*|n_z|`.
    pub h_out_tilt_scale: f64,
    /// Optional 5th-order polynomial coefficients for angular SHGC.
    /// `SHGC(theta) = SHGC_0 * sum(c[i] * cos^i(theta))` for i=0..5.
    /// When `Some`, overrides the ASHRAE IAM modifier.
    /// Coefficients should sum to 1.0 (normalization at theta=0).
    pub angular_shgc_coefficients: Option<[f64; 6]>,
    /// Per-window-pattern shading devices. Key is a substring pattern matched
    /// against the polygon path. The shading factor multiplies only the beam
    /// component; diffuse is unaffected.
    pub window_shading: HashMap<String, WindowShading>,
}

/// Shading devices (overhang and/or side fins) associated with a window.
#[derive(Debug, Clone)]
pub struct WindowShading {
    /// Optional overhang above the window.
    pub overhang: Option<OverhangGeometry>,
    /// Optional left side fin.
    pub fin_left: Option<FinGeometry>,
    /// Optional right side fin.
    pub fin_right: Option<FinGeometry>,
    /// Outward normal azimuth of the wall containing this window [degrees from north, clockwise].
    pub surface_azimuth_deg: f64,
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
            include_exterior_opaque_absorption: false,
            default_opaque_absorptance: DEFAULT_OPAQUE_ABSORPTANCE,
            exterior_heat_transfer_coeff_w_per_m2_k:
                DEFAULT_EXTERIOR_HEAT_TRANSFER_COEFF_W_PER_M2_K,
            include_ground_reflection: false,
            ground_reflectance: DEFAULT_GROUND_REFLECTANCE,
            include_incidence_angle_modifier: false,
            incidence_angle_modifier_a: DEFAULT_INC_ANGLE_MODIFIER_A,
            include_exterior_longwave_exchange: false,
            exterior_opaque_emissivity: 0.9,
            ground_emissivity: 0.95,
            use_wind_speed_for_h_out: false,
            h_out_base_w_per_m2_k: 5.0,
            h_out_wind_coeff_w_per_m2_k_per_m_s: 4.0,
            h_out_tilt_scale: 0.0,
            angular_shgc_coefficients: None,
            window_shading: HashMap::new(),
        }
    }

    /// Polynomial coefficients for single-pane clear glass (n ~ 1.526).
    /// Derived from Fresnel transmittance + absorption fit for 3mm clear float glass
    /// (BESTEST Glass Type 1, T_sol(0) = 0.834).
    pub fn single_pane_clear_coefficients() -> [f64; 6] {
        [0.0, 3.2695, -2.4987, -3.1251, 5.7123, -2.3581]
    }

    /// Returns the SHGC for a polygon path, or None if it's not glazing.
    ///
    /// Resolution order:
    /// 1. Exact match in `shgc` map
    /// 2. Material library `is_glazing` flag (if `material_library` provided)
    /// 3. Substring pattern match (fallback)
    pub fn resolve_shgc(
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

/// Computes the sunlit fraction for a window given its shading devices and solar position.
fn compute_sunlit_fraction(shading: &WindowShading, solar_pos: &SolarPosition) -> f64 {
    use super::shading::{overhang_and_fins_sunlit_fraction, overhang_sunlit_fraction};

    match (&shading.overhang, &shading.fin_left, &shading.fin_right) {
        (Some(oh), Some(fl), Some(fr)) => overhang_and_fins_sunlit_fraction(
            oh,
            fl,
            fr,
            solar_pos.altitude,
            solar_pos.azimuth,
            shading.surface_azimuth_deg,
        ),
        (Some(oh), _, _) => overhang_sunlit_fraction(
            oh,
            solar_pos.altitude,
            solar_pos.azimuth,
            shading.surface_azimuth_deg,
        ),
        _ => 1.0,
    }
}

/// Angular transmittance modifier for glazing.
///
/// When polynomial coefficients are provided, evaluates:
///   T(theta)/T(0) = sum(c[i] * cos^i(theta))
///
/// Otherwise falls back to ASHRAE IAM:
///   IAM = 1 - a * (1/cos(theta) - 1)
pub fn angular_transmittance_modifier(cos_incidence: f64, config: &SolarGainConfig) -> f64 {
    if !(cos_incidence.is_finite() && cos_incidence > 0.0) {
        return 0.0;
    }
    if let Some(c) = &config.angular_shgc_coefficients {
        let x = cos_incidence;
        let val = c[0]
            + c[1] * x
            + c[2] * x * x
            + c[3] * x * x * x
            + c[4] * x * x * x * x
            + c[5] * x * x * x * x * x;
        val.clamp(0.0, 1.0)
    } else if config.include_incidence_angle_modifier {
        let a = config.incidence_angle_modifier_a.max(0.0);
        (1.0 - a * (1.0 / cos_incidence - 1.0)).clamp(0.0, 1.0)
    } else {
        1.0
    }
}

/// Transmitted solar split into beam (direct) and diffuse components.
#[derive(Debug, Clone, Copy, Default)]
pub struct TransmittedSolarSplit {
    /// Direct (beam) component through glazing (W).
    pub beam_w: f64,
    /// Diffuse (sky + ground-reflected) component through glazing (W).
    pub diffuse_w: f64,
}

impl TransmittedSolarSplit {
    /// Total transmitted solar (W).
    pub fn total(&self) -> f64 {
        self.beam_w + self.diffuse_w
    }
}

/// Per-glazing transmitted solar contribution for one timestep.
#[derive(Debug, Clone)]
pub struct GlazingTransmission {
    /// Zone UID containing the glazing polygon.
    pub zone_uid: UID,
    /// Glazing polygon UID.
    pub polygon_uid: UID,
    /// Polygon outward normal.
    pub outward_normal: Vector,
    /// Polygon centroid.
    pub centroid: Point,
    /// Polygon area [m²].
    pub area_m2: f64,
    /// Transmitted direct (beam) power [W].
    pub beam_w: f64,
    /// Transmitted diffuse power [W].
    pub diffuse_w: f64,
}

/// Parameters describing the solar conditions for a single hour.
#[derive(Debug, Clone, Copy)]
pub struct SolarHourParams {
    /// Outdoor air temperature (°C).
    pub outdoor_air_temperature_c: f64,
    /// Global horizontal irradiance (W/m^2).
    pub global_horizontal_irradiance: f64,
    /// Direct normal irradiance (W/m^2).
    pub direct_normal_irradiance: f64,
    /// Diffuse horizontal irradiance (W/m^2).
    pub diffuse_horizontal_irradiance: f64,
    /// Horizontal infrared radiation intensity from sky (W/m^2).
    pub horizontal_infrared_radiation: f64,
    /// Wind speed (m/s).
    pub wind_speed: f64,
    /// Day of year (1-365).
    pub day_of_year: u16,
    /// Solar hour (0-24).
    ///
    /// Interpretation: local standard time (clock time) in hours.
    /// For EPW "hour-ending" records, pass `hour - 0.5` to use the mid-hour timestamp.
    pub local_time_hours: f64,
    /// Site latitude in degrees.
    pub latitude: f64,
    /// Site longitude in degrees.
    pub longitude: f64,
    /// Timezone (hours from UTC), as provided by EPW.
    pub timezone: f64,
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
    compute_solar_gains_with_materials(building, params, config, None).total()
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
) -> TransmittedSolarSplit {
    let mut total_beam = 0.0_f64;
    let mut total_diffuse = 0.0_f64;
    for t in
        compute_glazing_transmissions_with_materials(building, params, config, material_library)
    {
        total_beam += t.beam_w;
        total_diffuse += t.diffuse_w;
    }
    TransmittedSolarSplit {
        beam_w: total_beam,
        diffuse_w: total_diffuse,
    }
}

/// Computes per-glazing transmitted solar contributions.
pub fn compute_glazing_transmissions_with_materials(
    building: &Building,
    params: &SolarHourParams,
    config: &SolarGainConfig,
    material_library: Option<&MaterialLibrary>,
) -> Vec<GlazingTransmission> {
    let solar_pos = SolarPosition::calculate_from_local_time(
        params.latitude,
        params.longitude,
        params.timezone,
        params.day_of_year,
        params.local_time_hours,
    );
    let sun_above = solar_pos.is_above_horizon();
    let sun_dir = solar_pos.to_direction();
    let mut transmissions = Vec::new();

    let rho_g = config.ground_reflectance.clamp(0.0, 1.0);
    let use_ground = config.include_ground_reflection && rho_g > 0.0;
    let iam_diffuse = angular_transmittance_modifier(0.5, config); // cos(60°)

    for zone in building.zones() {
        for solid in zone.solids() {
            for wall in solid.walls() {
                for polygon in wall.polygons() {
                    let path = format!(
                        "{}/{}/{}/{}",
                        zone.name, solid.name, wall.name, polygon.name
                    );
                    let Some(shgc) = config.resolve_shgc(&path, material_library) else {
                        continue;
                    };

                    let area = polygon.area();
                    if area <= 0.0 {
                        continue;
                    }
                    let normal = polygon.vn;

                    // Diffuse sky component: isotropic view factor.
                    let sky_view = 0.5 * (1.0 + normal.dz.max(0.0));
                    let mut incident_diffuse =
                        params.diffuse_horizontal_irradiance.max(0.0) * sky_view * iam_diffuse;

                    // Ground reflected component: simple view factor approximation.
                    if use_ground {
                        let ground_view = 0.5 * (1.0 - normal.dz.clamp(-1.0, 1.0));
                        incident_diffuse += params.global_horizontal_irradiance.max(0.0)
                            * rho_g
                            * ground_view
                            * iam_diffuse;
                    }

                    // Direct beam component.
                    let mut incident_direct = 0.0;
                    if sun_above && params.direct_normal_irradiance > 0.0 {
                        let cos_incidence = sun_dir.dot(&normal).max(0.0);
                        if cos_incidence > 0.0 {
                            let mut i = params.direct_normal_irradiance.max(0.0) * cos_incidence;
                            i *= angular_transmittance_modifier(cos_incidence, config);
                            incident_direct = i;
                        }
                    }

                    let mut beam_w = incident_direct * area * shgc;
                    let diffuse_w = incident_diffuse * area * shgc;

                    // Apply overhang/fin shading to beam component only.
                    if beam_w > 0.0 && !config.window_shading.is_empty() {
                        for (pattern, shading) in &config.window_shading {
                            if path.contains(pattern.as_str()) {
                                let sunlit = compute_sunlit_fraction(shading, &solar_pos);
                                beam_w *= sunlit;
                                break;
                            }
                        }
                    }

                    if beam_w <= 0.0 && diffuse_w <= 0.0 {
                        continue;
                    }

                    transmissions.push(GlazingTransmission {
                        zone_uid: zone.uid.clone(),
                        polygon_uid: polygon.uid.clone(),
                        outward_normal: normal,
                        centroid: polygon.centroid(),
                        area_m2: area,
                        beam_w,
                        diffuse_w,
                    });
                }
            }
        }
    }

    transmissions
}

/// Like [`compute_solar_gains_per_zone`], but also accepts a material library for
/// `is_glazing`-based surface identification.
pub fn compute_solar_gains_per_zone_with_materials(
    building: &Building,
    params: &SolarHourParams,
    config: &SolarGainConfig,
    material_library: Option<&MaterialLibrary>,
) -> HashMap<crate::UID, f64> {
    let mut gains: HashMap<crate::UID, f64> = HashMap::new();
    for t in
        compute_glazing_transmissions_with_materials(building, params, config, material_library)
    {
        *gains.entry(t.zone_uid).or_insert(0.0) += t.beam_w + t.diffuse_w;
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
            outdoor_air_temperature_c: 20.0,
            global_horizontal_irradiance: 700.0,
            direct_normal_irradiance: 500.0,
            diffuse_horizontal_irradiance: 200.0,
            horizontal_infrared_radiation: 300.0,
            wind_speed: 3.0,
            day_of_year: 172,
            local_time_hours: 12.0,
            latitude: 45.0,
            longitude: 0.0,
            timezone: 0.0,
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
            outdoor_air_temperature_c: 20.0,
            global_horizontal_irradiance: 0.0,
            direct_normal_irradiance: 0.0,
            diffuse_horizontal_irradiance: 0.0,
            horizontal_infrared_radiation: 300.0,
            wind_speed: 3.0,
            day_of_year: 172,
            local_time_hours: 12.0,
            latitude: 45.0,
            longitude: 0.0,
            timezone: 0.0,
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
            outdoor_air_temperature_c: 20.0,
            global_horizontal_irradiance: 900.0,
            direct_normal_irradiance: 800.0,
            diffuse_horizontal_irradiance: 100.0,
            horizontal_infrared_radiation: 300.0,
            wind_speed: 3.0,
            day_of_year: 81,       // ~equinox
            local_time_hours: 9.0, // morning sun (east-ish)
            latitude: 0.0,
            longitude: 0.0,
            timezone: 0.0,
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

        let shgc = cfg.resolve_shgc("zone/room/wall/window_0", Some(&lib));
        assert_eq!(shgc, Some(0.75));
    }
}
