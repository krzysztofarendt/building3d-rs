//! Dynamic convection coefficient models aligned with EnergyPlus.
//!
//! Two opt-in models are provided:
//!
//! - **TARP** (interior): temperature- and tilt-dependent natural convection
//!   correlations from Walton (1983) as implemented in EnergyPlus.
//!
//! - **DOE-2** (exterior): combined natural + forced (wind-driven) convection
//!   using the DOE-2 simplified model.
//!
//! Default behaviour is fixed coefficients (`Fixed` variants) which preserves
//! backward compatibility with existing simulations.

/// Minimum interior convection coefficient for floor surfaces [W/(m²·K)].
///
/// EnergyPlus uses 1.0 for consolidated interior algorithm. We use a
/// conservative 0.1 to avoid numerical instabilities with nearly-adiabatic
/// slabs while still allowing very low natural convection on stable floors.
const H_MIN_FLOOR: f64 = 0.1;

// ----- Interior convection model -----

/// Interior surface convection model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum InteriorConvectionModel {
    /// Fixed coefficient [W/(m²·K)]. This is the legacy default (3.0).
    Fixed(f64),
    /// TARP (Walton) natural convection correlations.
    ///
    /// The coefficient depends on the surface-to-air temperature difference
    /// and the surface tilt (via `cos_tilt`).
    Tarp,
}

impl Default for InteriorConvectionModel {
    fn default() -> Self {
        Self::Fixed(3.0)
    }
}

/// Exterior surface convection model.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExteriorConvectionModel {
    /// Fixed coefficient derived from ISO R_se (legacy default).
    Fixed,
    /// DOE-2 simplified model: combined natural + wind-forced convection.
    Doe2,
}

impl Default for ExteriorConvectionModel {
    fn default() -> Self {
        Self::Fixed
    }
}

/// TARP natural convection coefficient [W/(m²·K)].
///
/// Correlations from Walton (1983) as used in EnergyPlus:
///
/// - Vertical (|cos_tilt| < 0.707):
///   `h = 1.31 * |dT|^(1/3)`
///
/// - Tilted/horizontal unstable (heated floor or cooled ceiling,
///   i.e. `dT * cos_tilt < 0`):
///   `h = 9.482 * |dT|^(1/3) / (7.238 - |cos_tilt|)`
///
/// - Tilted/horizontal stable (cooled floor or heated ceiling):
///   `h = 1.810 * |dT|^(1/3) / (1.382 + |cos_tilt|)`
///
/// # Arguments
/// - `dt`: surface temperature minus air temperature [K or °C difference]
/// - `cos_tilt`: cosine of the surface tilt from horizontal (`polygon.vn.dz`).
///   `cos_tilt = -1` → floor facing down (interior side up),
///   `cos_tilt = +1` → ceiling facing up (interior side down for ext. roof).
pub fn tarp_natural_h(dt: f64, cos_tilt: f64) -> f64 {
    let abs_dt = dt.abs();
    if abs_dt < 1e-15 {
        return H_MIN_FLOOR;
    }
    let dt_third = abs_dt.powf(1.0 / 3.0);
    let abs_cos = cos_tilt.abs();

    if abs_cos < 0.707 {
        // Near-vertical surface.
        1.31 * dt_third
    } else if dt * cos_tilt < 0.0 {
        // Unstable: buoyancy-driven enhanced convection
        // (heated surface facing up, or cooled surface facing down).
        9.482 * dt_third / (7.238 - abs_cos)
    } else {
        // Stable: suppressed convection.
        1.810 * dt_third / (1.382 + abs_cos)
    }
}

/// Compute interior convection coefficient for a surface.
///
/// When model is `Fixed`, returns the fixed value. When model is `Tarp`,
/// uses the TARP correlation but clamps to at least `h_min_iso` (the
/// ISO 6946 film coefficient from the construction's R_si).
///
/// # Arguments
/// - `model`: selected interior convection model
/// - `dt`: `T_surface - T_air` [K]
/// - `cos_tilt`: surface tilt cosine (polygon normal dz)
/// - `h_min_iso`: minimum coefficient from ISO R_si [W/(m²·K)]
pub fn interior_convection_h(
    model: &InteriorConvectionModel,
    dt: f64,
    cos_tilt: f64,
    h_min_iso: f64,
) -> f64 {
    match model {
        InteriorConvectionModel::Fixed(h) => *h,
        InteriorConvectionModel::Tarp => {
            let h = tarp_natural_h(dt, cos_tilt);
            h.max(h_min_iso).max(H_MIN_FLOOR)
        }
    }
}

// ----- Exterior convection model -----

/// DOE-2 wind-driven forced convection coefficients for rough surfaces.
const DOE2_A_ROUGH: f64 = 3.26;
const DOE2_B_ROUGH: f64 = 3.89;

/// DOE-2 simplified exterior convection coefficient [W/(m²·K)].
///
/// Combined natural + forced (wind):
/// - Natural component: same TARP correlation
/// - Forced component: `a + b * V_wind` (rough surface coefficients)
/// - Combined: `sqrt(h_n² + h_f²)`
///
/// # Arguments
/// - `dt`: `T_surface - T_outdoor_air` [K]
/// - `cos_tilt`: surface tilt cosine (polygon normal dz)
/// - `wind_speed`: outdoor wind speed [m/s]
pub fn doe2_exterior_h(dt: f64, cos_tilt: f64, wind_speed: f64) -> f64 {
    let h_natural = tarp_natural_h(dt, cos_tilt);
    let v = wind_speed.max(0.0);
    let h_forced = DOE2_A_ROUGH + DOE2_B_ROUGH * v;
    (h_natural * h_natural + h_forced * h_forced).sqrt()
}

/// Compute exterior convection coefficient for a surface.
///
/// When model is `Fixed`, returns the provided `h_fixed` (from ISO R_se).
/// When model is `Doe2`, uses the DOE-2 combined correlation.
///
/// # Arguments
/// - `model`: selected exterior convection model
/// - `h_fixed`: fixed coefficient [W/(m²·K)] (used only for `Fixed` model)
/// - `dt`: `T_surface - T_outdoor_air` [K]
/// - `cos_tilt`: surface tilt cosine (polygon normal dz)
/// - `wind_speed`: outdoor wind speed [m/s]
pub fn exterior_convection_h(
    model: &ExteriorConvectionModel,
    h_fixed: f64,
    dt: f64,
    cos_tilt: f64,
    wind_speed: f64,
) -> f64 {
    match model {
        ExteriorConvectionModel::Fixed => h_fixed,
        ExteriorConvectionModel::Doe2 => doe2_exterior_h(dt, cos_tilt, wind_speed),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tarp_vertical_surface() {
        // Vertical wall: cos_tilt ≈ 0 → h = 1.31 * |dT|^(1/3)
        let dt = 8.0_f64;
        let h = tarp_natural_h(dt, 0.0);
        let expected = 1.31 * dt.powf(1.0 / 3.0);
        assert!(
            (h - expected).abs() < 1e-10,
            "Vertical: got {h}, expected {expected}"
        );
    }

    #[test]
    fn test_tarp_horizontal_unstable() {
        // Heated floor facing down (interior side up): cos_tilt = -1.0, dt > 0 → dt*cos < 0 → unstable
        let dt = 5.0;
        let cos_tilt = -1.0;
        let h = tarp_natural_h(dt, cos_tilt);
        let expected = 9.482 * dt.powf(1.0 / 3.0) / (7.238 - 1.0);
        assert!(
            (h - expected).abs() < 1e-10,
            "Unstable: got {h}, expected {expected}"
        );
    }

    #[test]
    fn test_tarp_horizontal_stable() {
        // Cooled floor (interior side up): cos_tilt = -1.0, dt < 0 → dt*cos > 0 → stable
        let dt = -5.0;
        let cos_tilt = -1.0;
        let h = tarp_natural_h(dt, cos_tilt);
        let expected = 1.810 * 5.0_f64.powf(1.0 / 3.0) / (1.382 + 1.0);
        assert!(
            (h - expected).abs() < 1e-10,
            "Stable: got {h}, expected {expected}"
        );
    }

    #[test]
    fn test_tarp_zero_dt() {
        let h = tarp_natural_h(0.0, 0.0);
        assert!(
            (h - H_MIN_FLOOR).abs() < 1e-12,
            "Zero dT should return H_MIN_FLOOR"
        );
    }

    #[test]
    fn test_interior_convection_fixed() {
        let h = interior_convection_h(&InteriorConvectionModel::Fixed(7.69), 5.0, 0.0, 1.0);
        assert!((h - 7.69).abs() < 1e-12);
    }

    #[test]
    fn test_interior_convection_tarp_clamps_to_min() {
        // Very small dT → TARP gives tiny h, but should be clamped to h_min_iso.
        let h_min = 5.0;
        let h = interior_convection_h(&InteriorConvectionModel::Tarp, 0.001, 0.0, h_min);
        assert!(
            h >= h_min - 1e-12,
            "Should be at least h_min_iso={h_min}, got {h}"
        );
    }

    #[test]
    fn test_doe2_zero_wind() {
        // No wind → forced = a, combined = sqrt(natural² + a²)
        let dt = 10.0;
        let cos_tilt = 0.0;
        let h = doe2_exterior_h(dt, cos_tilt, 0.0);
        let h_n = tarp_natural_h(dt, cos_tilt);
        let h_f = DOE2_A_ROUGH;
        let expected = (h_n * h_n + h_f * h_f).sqrt();
        assert!(
            (h - expected).abs() < 1e-10,
            "DOE-2 zero wind: got {h}, expected {expected}"
        );
    }

    #[test]
    fn test_doe2_with_wind() {
        let dt = 5.0;
        let cos_tilt = -0.5;
        let wind = 4.0;
        let h = doe2_exterior_h(dt, cos_tilt, wind);
        let h_n = tarp_natural_h(dt, cos_tilt);
        let h_f = DOE2_A_ROUGH + DOE2_B_ROUGH * wind;
        let expected = (h_n * h_n + h_f * h_f).sqrt();
        assert!(
            (h - expected).abs() < 1e-10,
            "DOE-2 with wind: got {h}, expected {expected}"
        );
        // With 4 m/s wind, forced component should dominate.
        assert!(h > 15.0, "Expected significant h with 4 m/s wind, got {h}");
    }

    #[test]
    fn test_exterior_convection_fixed() {
        let h = exterior_convection_h(&ExteriorConvectionModel::Fixed, 25.0, 10.0, 0.0, 5.0);
        assert!((h - 25.0).abs() < 1e-12, "Fixed model ignores dT/wind");
    }

    #[test]
    fn test_exterior_convection_doe2() {
        let h = exterior_convection_h(&ExteriorConvectionModel::Doe2, 25.0, 10.0, 0.0, 5.0);
        let expected = doe2_exterior_h(10.0, 0.0, 5.0);
        assert!(
            (h - expected).abs() < 1e-12,
            "Doe2 model should use doe2_exterior_h"
        );
    }

    #[test]
    fn test_tarp_slightly_tilted() {
        // Surface tilted 50° from horizontal → cos_tilt ≈ 0.643 (< 0.707 → vertical branch)
        let cos_tilt = 0.643;
        let dt = 3.0;
        let h = tarp_natural_h(dt, cos_tilt);
        let expected = 1.31 * dt.powf(1.0 / 3.0);
        assert!(
            (h - expected).abs() < 1e-10,
            "Slightly tilted should use vertical: got {h}, expected {expected}"
        );
    }

    #[test]
    fn test_tarp_ceiling_heated_stable() {
        // Heated ceiling (interior side faces down): cos_tilt = 1.0, dt > 0 → dt*cos > 0 → stable
        let dt = 3.0;
        let cos_tilt = 1.0;
        let h = tarp_natural_h(dt, cos_tilt);
        let expected = 1.810 * dt.powf(1.0 / 3.0) / (1.382 + 1.0);
        assert!(
            (h - expected).abs() < 1e-10,
            "Heated ceiling stable: got {h}, expected {expected}"
        );
    }

    #[test]
    fn test_defaults() {
        assert_eq!(
            InteriorConvectionModel::default(),
            InteriorConvectionModel::Fixed(3.0)
        );
        assert_eq!(
            ExteriorConvectionModel::default(),
            ExteriorConvectionModel::Fixed
        );
    }

    #[test]
    fn test_doe2_negative_wind_clamped() {
        // Negative wind speed should be clamped to zero.
        let h_neg = doe2_exterior_h(5.0, 0.0, -3.0);
        let h_zero = doe2_exterior_h(5.0, 0.0, 0.0);
        assert!(
            (h_neg - h_zero).abs() < 1e-12,
            "Negative wind should equal zero wind"
        );
    }
}
