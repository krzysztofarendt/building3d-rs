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

        // Solar azimuth (atan2-based for numerical robustness)
        let cos_alt = altitude.to_radians().cos().max(1e-10);
        let sin_azi = -declination.cos() * hour_angle.sin() / cos_alt;
        let cos_azi = (declination.sin() * lat.cos()
            - declination.cos() * lat.sin() * hour_angle.cos())
            / cos_alt;

        let mut azimuth = sin_azi.atan2(cos_azi).to_degrees();
        if azimuth < 0.0 {
            azimuth += 360.0;
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

    /// Calculates solar position from local standard time (clock time) using an
    /// equation-of-time + longitude/timezone-meridian correction.
    ///
    /// - `longitude`: degrees (positive east)
    /// - `timezone`: hours from UTC (EPW header value)
    /// - `local_time_hours`: local standard time in hours (0-24)
    ///
    /// Note: EPW records are typically "hour-ending"; callers often want to pass
    /// `hour - 0.5` to represent the mid-hour timestamp.
    pub fn calculate_from_local_time(
        latitude: f64,
        longitude: f64,
        timezone: f64,
        day_of_year: u16,
        local_time_hours: f64,
    ) -> Self {
        // Equation of time (Spencer), minutes.
        let gamma = 2.0 * std::f64::consts::PI * (day_of_year as f64 - 1.0) / 365.0;
        let eot_min = 229.18
            * (0.000075 + 0.001868 * gamma.cos()
                - 0.032077 * gamma.sin()
                - 0.014615 * (2.0 * gamma).cos()
                - 0.040849 * (2.0 * gamma).sin());

        // Local standard time meridian, degrees.
        let lstm = 15.0 * timezone;

        // Time correction factor, minutes.
        let tc_min = 4.0 * (longitude - lstm) + eot_min;

        let solar_time_hours = local_time_hours + tc_min / 60.0;
        Self::calculate(latitude, longitude, day_of_year, solar_time_hours)
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
    fn test_solar_azimuth_morning_vs_afternoon() {
        // At 45°N on equinox, morning sun should be east (azi < 180),
        // afternoon sun should be west (azi > 180).
        let morning = SolarPosition::calculate(45.0, 0.0, 80, 8.0);
        let afternoon = SolarPosition::calculate(45.0, 0.0, 80, 16.0);
        assert!(
            morning.azimuth < 180.0,
            "Morning azimuth should be east (<180), got {}",
            morning.azimuth
        );
        assert!(
            afternoon.azimuth > 180.0,
            "Afternoon azimuth should be west (>180), got {}",
            afternoon.azimuth
        );
        // Azimuths should be roughly symmetric around south (180)
        let diff = (morning.azimuth + afternoon.azimuth) / 2.0;
        assert!(
            (diff - 180.0).abs() < 5.0,
            "Mean azimuth should be ~180 (south), got {diff}"
        );
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
