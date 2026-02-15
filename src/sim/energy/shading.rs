/// Analytical overhang and fin shading models for windows.
///
/// Computes the fraction of a window area that is sunlit given the geometry
/// of shading devices (overhangs, side fins) and the current sun position.
///
/// Conventions:
/// - Solar azimuth: degrees from north, clockwise (0=N, 90=E, 180=S, 270=W).
/// - Solar altitude: degrees above horizon (0=horizon, 90=zenith).
/// - Surface azimuth: outward normal azimuth in the same convention.

/// Geometry of a horizontal overhang above a window.
#[derive(Debug, Clone, Copy)]
pub struct OverhangGeometry {
    /// Projection depth from the wall plane [m].
    pub depth: f64,
    /// Vertical gap between the top of the window and the overhang [m].
    pub gap: f64,
    /// Window height [m].
    pub window_height: f64,
    /// Window width [m].
    pub window_width: f64,
    /// Overhang lateral extension beyond the left edge of the window [m].
    pub extension_left: f64,
    /// Overhang lateral extension beyond the right edge of the window [m].
    pub extension_right: f64,
}

/// Geometry of a vertical side fin adjacent to a window.
#[derive(Debug, Clone, Copy)]
pub struct FinGeometry {
    /// Projection depth from the wall plane [m].
    pub depth: f64,
    /// Vertical extent of the fin below the overhang (or top of window) [m].
    pub height: f64,
}

/// Returns the fraction of the window that is sunlit (0.0 = fully shaded, 1.0 = fully sunlit).
///
/// Uses a simple geometric projection of the overhang shadow onto the window plane.
/// Only the beam (direct) component is affected; diffuse sky radiation is not modified.
///
/// # Arguments
/// * `overhang` - Overhang geometry
/// * `solar_altitude_deg` - Solar altitude angle in degrees
/// * `solar_azimuth_deg` - Solar azimuth angle in degrees from north, clockwise
/// * `surface_azimuth_deg` - Outward normal azimuth of the wall in degrees from north, clockwise
pub fn overhang_sunlit_fraction(
    overhang: &OverhangGeometry,
    solar_altitude_deg: f64,
    solar_azimuth_deg: f64,
    surface_azimuth_deg: f64,
) -> f64 {
    if overhang.depth <= 0.0 || overhang.window_height <= 0.0 || overhang.window_width <= 0.0 {
        return 1.0;
    }
    if solar_altitude_deg <= 0.0 {
        return 1.0; // sun below horizon
    }

    let alt_rad = solar_altitude_deg.to_radians();
    let gamma_rad = (solar_azimuth_deg - surface_azimuth_deg).to_radians();
    let cos_gamma = gamma_rad.cos();

    // If the sun is behind the surface, the window is fully sunlit (no shadow from overhang).
    if cos_gamma <= 0.0 {
        return 1.0;
    }

    // Vertical shadow depth: how far below the overhang the shadow extends on the wall.
    let tan_alt = alt_rad.tan();
    if tan_alt <= 0.0 {
        return 1.0;
    }
    let vertical_shadow = overhang.depth * cos_gamma / tan_alt;

    // The shadow starts at the overhang position, which is `gap` above the window top.
    // So the shadow on the window = vertical_shadow - gap.
    let shadow_on_window = (vertical_shadow - overhang.gap).clamp(0.0, overhang.window_height);
    if shadow_on_window <= 0.0 {
        return 1.0;
    }

    // Horizontal shadow offset from the overhang edge.
    let sin_gamma = gamma_rad.sin();
    let horizontal_shadow_offset = overhang.depth * sin_gamma / cos_gamma; // = depth * tan(gamma)

    // The overhang extends from -extension_left to (window_width + extension_right) relative
    // to the left edge of the window. The shadow is shifted horizontally.
    let overhang_total_width =
        overhang.window_width + overhang.extension_left + overhang.extension_right;

    // Shadow strip on the wall plane (relative to window left edge = 0):
    //   The overhang left edge is at -extension_left on the window.
    //   Shadow shifts by horizontal_shadow_offset.
    let shadow_left = -overhang.extension_left + horizontal_shadow_offset;
    let shadow_right = shadow_left + overhang_total_width;

    // Overlap with window [0, window_width]
    let overlap_left = shadow_left.max(0.0);
    let overlap_right = shadow_right.min(overhang.window_width);
    let overlap_width = (overlap_right - overlap_left).max(0.0);

    let shaded_area = shadow_on_window * overlap_width;
    let window_area = overhang.window_height * overhang.window_width;

    (1.0 - shaded_area / window_area).clamp(0.0, 1.0)
}

/// Returns the fraction of the window that is sunlit considering both an overhang and two side fins.
///
/// The overhang and fin shadows are computed independently and combined conservatively
/// (union of shaded regions approximated as the sum minus overlap).
pub fn overhang_and_fins_sunlit_fraction(
    overhang: &OverhangGeometry,
    fin_left: &FinGeometry,
    fin_right: &FinGeometry,
    solar_altitude_deg: f64,
    solar_azimuth_deg: f64,
    surface_azimuth_deg: f64,
) -> f64 {
    if solar_altitude_deg <= 0.0 {
        return 1.0;
    }

    let window_area = overhang.window_height * overhang.window_width;
    if window_area <= 0.0 {
        return 1.0;
    }

    let alt_rad = solar_altitude_deg.to_radians();
    let gamma_rad = (solar_azimuth_deg - surface_azimuth_deg).to_radians();
    let cos_gamma = gamma_rad.cos();
    let sin_gamma = gamma_rad.sin();

    // If the sun is behind the surface, fully sunlit.
    if cos_gamma <= 0.0 {
        return 1.0;
    }

    // ── Overhang shadow ──
    let overhang_shaded_fraction = 1.0
        - overhang_sunlit_fraction(
            overhang,
            solar_altitude_deg,
            solar_azimuth_deg,
            surface_azimuth_deg,
        );

    // ── Fin shadows ──
    // Fins project a horizontal shadow across the window.
    // Left fin is at x=0 (left edge of window), right fin at x=window_width.
    let tan_gamma = sin_gamma / cos_gamma; // positive = sun from left, negative = sun from right
    let _tan_alt = alt_rad.tan().max(1e-10);

    // Left fin shadow (fin is on the left edge of the window, projecting rightward when sun is from left)
    let mut left_fin_shaded = 0.0;
    if fin_left.depth > 0.0 && fin_left.height > 0.0 {
        // Horizontal projection of the fin shadow on the wall
        let h_shadow = fin_left.depth * tan_gamma.abs();
        // The fin shadow projects rightward when tan_gamma > 0 (sun from left side)
        if tan_gamma > 0.0 && h_shadow > 0.0 {
            let shadow_width = h_shadow.min(overhang.window_width);
            // Vertical extent: the fin shadow height on the wall
            let v_shadow = (fin_left.height).min(overhang.window_height);
            left_fin_shaded = shadow_width * v_shadow / window_area;
        }
    }

    // Right fin shadow (fin is on the right edge of the window, projecting leftward when sun is from right)
    let mut right_fin_shaded = 0.0;
    if fin_right.depth > 0.0 && fin_right.height > 0.0 {
        let h_shadow = fin_right.depth * tan_gamma.abs();
        if tan_gamma < 0.0 && h_shadow > 0.0 {
            let shadow_width = h_shadow.min(overhang.window_width);
            let v_shadow = (fin_right.height).min(overhang.window_height);
            right_fin_shaded = shadow_width * v_shadow / window_area;
        }
    }

    // Combine: approximate union as min(1, sum) since overlaps are typically small.
    let total_shaded = (overhang_shaded_fraction + left_fin_shaded + right_fin_shaded).min(1.0);
    (1.0 - total_shaded).max(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_no_overhang_fully_sunlit() {
        let oh = OverhangGeometry {
            depth: 0.0,
            gap: 0.0,
            window_height: 2.0,
            window_width: 3.0,
            extension_left: 0.0,
            extension_right: 0.0,
        };
        let f = overhang_sunlit_fraction(&oh, 45.0, 180.0, 180.0);
        assert!(
            (f - 1.0).abs() < 1e-10,
            "No overhang should be fully sunlit"
        );
    }

    #[test]
    fn test_sun_behind_surface() {
        let oh = OverhangGeometry {
            depth: 1.0,
            gap: 0.0,
            window_height: 2.0,
            window_width: 3.0,
            extension_left: 0.0,
            extension_right: 0.0,
        };
        // Sun at azimuth 0 (north), surface facing south (180).
        let f = overhang_sunlit_fraction(&oh, 45.0, 0.0, 180.0);
        assert!(
            (f - 1.0).abs() < 1e-10,
            "Sun behind surface should be fully sunlit"
        );
    }

    #[test]
    fn test_overhang_partial_shading() {
        // 1m deep overhang, no gap, 2m tall window, sun at 45° altitude directly facing.
        // Vertical shadow = 1.0 * cos(0) / tan(45°) = 1.0m on 2m window = 50% shaded
        let oh = OverhangGeometry {
            depth: 1.0,
            gap: 0.0,
            window_height: 2.0,
            window_width: 3.0,
            extension_left: 10.0, // wide enough to cover full width
            extension_right: 10.0,
        };
        let f = overhang_sunlit_fraction(&oh, 45.0, 180.0, 180.0);
        assert!((f - 0.5).abs() < 0.01, "Expected ~50% sunlit, got {f}");
    }

    #[test]
    fn test_overhang_with_gap() {
        // 1m deep overhang, 0.5m gap, 2m window, sun at 45° facing.
        // Vertical shadow = 1.0m, but gap = 0.5m, so shadow on window = 0.5m.
        // Fraction shaded = 0.5/2.0 = 25%
        let oh = OverhangGeometry {
            depth: 1.0,
            gap: 0.5,
            window_height: 2.0,
            window_width: 3.0,
            extension_left: 10.0,
            extension_right: 10.0,
        };
        let f = overhang_sunlit_fraction(&oh, 45.0, 180.0, 180.0);
        assert!((f - 0.75).abs() < 0.01, "Expected ~75% sunlit, got {f}");
    }

    #[test]
    fn test_sun_below_horizon() {
        let oh = OverhangGeometry {
            depth: 1.0,
            gap: 0.0,
            window_height: 2.0,
            window_width: 3.0,
            extension_left: 0.0,
            extension_right: 0.0,
        };
        let f = overhang_sunlit_fraction(&oh, -5.0, 180.0, 180.0);
        assert!(
            (f - 1.0).abs() < 1e-10,
            "Sun below horizon should be fully sunlit"
        );
    }

    #[test]
    fn test_high_sun_fully_shaded() {
        // Very deep overhang, low sun altitude = deep shadow that covers entire window.
        let oh = OverhangGeometry {
            depth: 5.0,
            gap: 0.0,
            window_height: 2.0,
            window_width: 3.0,
            extension_left: 10.0,
            extension_right: 10.0,
        };
        // tan(15°) ≈ 0.268, shadow = 5.0/0.268 ≈ 18.7m >> 2m window
        let f = overhang_sunlit_fraction(&oh, 15.0, 180.0, 180.0);
        assert!(
            f < 0.01,
            "Deep overhang at low sun should fully shade, got {f}"
        );
    }

    #[test]
    fn test_overhang_and_fins_no_devices() {
        let oh = OverhangGeometry {
            depth: 0.0,
            gap: 0.0,
            window_height: 2.0,
            window_width: 3.0,
            extension_left: 0.0,
            extension_right: 0.0,
        };
        let fl = FinGeometry {
            depth: 0.0,
            height: 0.0,
        };
        let fr = FinGeometry {
            depth: 0.0,
            height: 0.0,
        };
        let f = overhang_and_fins_sunlit_fraction(&oh, &fl, &fr, 45.0, 180.0, 180.0);
        assert!((f - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_fins_add_shading() {
        let oh = OverhangGeometry {
            depth: 1.0,
            gap: 0.0,
            window_height: 2.0,
            window_width: 3.0,
            extension_left: 10.0,
            extension_right: 10.0,
        };
        let fl = FinGeometry {
            depth: 1.0,
            height: 2.0,
        };
        let fr = FinGeometry {
            depth: 1.0,
            height: 2.0,
        };

        // Overhang only (sun directly facing)
        let f_oh = overhang_sunlit_fraction(&oh, 45.0, 180.0, 180.0);

        // With fins, sun slightly from the left (azimuth 170 vs surface 180)
        let f_full = overhang_and_fins_sunlit_fraction(&oh, &fl, &fr, 45.0, 170.0, 180.0);

        // Fins should add some shading when sun is from the side
        assert!(f_full <= f_oh + 0.01, "Fins should add shading or be equal");
    }
}
