use crate::Vector;

/// Solar position (azimuth and elevation angles).
#[derive(Debug, Clone, Copy)]
pub struct SolarPosition {
    /// Solar altitude angle in degrees (0 = horizon, 90 = zenith).
    pub altitude: f64,
    /// Solar azimuth angle in degrees from north, clockwise (0=N, 90=E, 180=S, 270=W).
    pub azimuth: f64,
}

impl SolarPosition {
    /// Calculates the solar position using the Spencer algorithm.
    ///
    /// - `latitude`: in degrees (positive north)
    /// - `longitude`: in degrees (positive east)
    /// - `day_of_year`: 1-365
    /// - `hour`: solar time in hours (0-24)
    pub fn calculate(latitude: f64, longitude: f64, day_of_year: u16, hour: f64) -> Self {
        let lat = latitude.to_radians();
        let _lon = longitude; // Longitude only affects solar time correction

        // Day angle (Spencer)
        let gamma = 2.0 * std::f64::consts::PI * (day_of_year as f64 - 1.0) / 365.0;

        // Solar declination (Spencer approximation)
        let declination = 0.006918 - 0.399912 * gamma.cos() + 0.070257 * gamma.sin()
            - 0.006758 * (2.0 * gamma).cos()
            + 0.000907 * (2.0 * gamma).sin()
            - 0.002697 * (3.0 * gamma).cos()
            + 0.00148 * (3.0 * gamma).sin();

        // Hour angle (15 degrees per hour from solar noon)
        let hour_angle = (hour - 12.0) * 15.0_f64.to_radians();

        // Solar altitude
        let sin_alt =
            lat.sin() * declination.sin() + lat.cos() * declination.cos() * hour_angle.cos();
        let altitude = sin_alt.asin().to_degrees();

        // Solar azimuth
        let cos_azimuth = (declination.sin() * lat.cos()
            - declination.cos() * lat.sin() * hour_angle.cos())
            / altitude.to_radians().cos().max(1e-10);

        let mut azimuth = cos_azimuth.clamp(-1.0, 1.0).acos().to_degrees();
        if hour_angle > 0.0 {
            azimuth = 360.0 - azimuth;
        }

        Self { altitude, azimuth }
    }

    /// Returns true if the sun is above the horizon.
    pub fn is_above_horizon(&self) -> bool {
        self.altitude > 0.0
    }

    /// Converts solar position to a direction vector (pointing toward the sun).
    pub fn to_direction(&self) -> Vector {
        let alt = self.altitude.to_radians();
        let azi = self.azimuth.to_radians();

        // Convention: azimuth from north clockwise
        // North = +Y, East = +X
        let x = alt.cos() * azi.sin();
        let y = alt.cos() * azi.cos();
        let z = alt.sin();

        Vector::new(x, y, z)
    }

    /// Returns the direction vector from scene toward the sun (for use as light direction).
    /// This is the opposite of `to_direction()` — it gives the direction light travels.
    pub fn light_direction(&self) -> Vector {
        let d = self.to_direction();
        Vector::new(-d.dx, -d.dy, -d.dz)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_solar_noon_equator_equinox() {
        // At solar noon on the equinox, sun should be directly overhead at equator
        let pos = SolarPosition::calculate(0.0, 0.0, 80, 12.0); // March equinox ~ day 80
        assert!(
            pos.altitude > 80.0,
            "Sun should be near zenith at equator on equinox noon"
        );
        assert!(pos.is_above_horizon());
    }

    #[test]
    fn test_solar_midnight() {
        // At midnight, sun should be below horizon
        let _pos = SolarPosition::calculate(45.0, 0.0, 172, 0.0);
        // At high latitudes in summer, sun might still be up at midnight
        // Use winter instead
        let pos_winter = SolarPosition::calculate(45.0, 0.0, 355, 0.0);
        assert!(
            !pos_winter.is_above_horizon(),
            "Sun should be below horizon at midnight in winter"
        );
    }

    #[test]
    fn test_direction_vector() {
        let pos = SolarPosition {
            altitude: 90.0,
            azimuth: 0.0,
        };
        let dir = pos.to_direction();
        // Sun at zenith should point straight up
        assert!((dir.dz - 1.0).abs() < 1e-6);
        assert!(dir.dx.abs() < 1e-6);
        assert!((dir.dy - 0.0).abs() < 0.1); // cos(90) is approximately 0
    }

    #[test]
    fn test_light_direction_opposite() {
        let pos = SolarPosition {
            altitude: 45.0,
            azimuth: 180.0,
        };
        let sun_dir = pos.to_direction();
        let light_dir = pos.light_direction();
        // Should be opposite
        assert!((sun_dir.dx + light_dir.dx).abs() < 1e-10);
        assert!((sun_dir.dy + light_dir.dy).abs() < 1e-10);
        assert!((sun_dir.dz + light_dir.dz).abs() < 1e-10);
    }
}
