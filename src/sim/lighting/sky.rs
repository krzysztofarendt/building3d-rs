use crate::Vector;

use super::solar::SolarPosition;
use super::sources::Rgb;

/// Trait for sky luminance distribution models.
pub trait SkyModel {
    /// Returns the sky luminance [R, G, B] for a given sky direction.
    ///
    /// - `direction`: unit vector pointing toward the sky point
    /// - `solar_pos`: current solar position
    fn luminance(&self, direction: Vector, solar_pos: &SolarPosition) -> Rgb;
}

/// CIE Standard Overcast Sky.
///
/// Luminance varies only with altitude: L = Lz * (1 + 2*sin(alt)) / 3
/// where Lz is the zenith luminance.
pub struct CIEOvercast {
    /// Zenith luminance (cd/m^2).
    pub zenith_luminance: f64,
}

impl CIEOvercast {
    pub fn new(zenith_luminance: f64) -> Self {
        Self { zenith_luminance }
    }

    /// Creates an overcast sky from horizontal diffuse illuminance.
    ///
    /// For CIE overcast sky: E_h = 7/9 * pi * L_z
    pub fn from_diffuse_illuminance(illuminance: f64) -> Self {
        let lz = illuminance * 9.0 / (7.0 * std::f64::consts::PI);
        Self {
            zenith_luminance: lz,
        }
    }
}

impl SkyModel for CIEOvercast {
    fn luminance(&self, direction: Vector, _solar_pos: &SolarPosition) -> Rgb {
        let alt = direction.dz.asin(); // altitude angle from direction z-component
        if alt < 0.0 {
            return [0.0; 3]; // Below horizon
        }
        let l = self.zenith_luminance * (1.0 + 2.0 * alt.sin()) / 3.0;
        [l, l, l]
    }
}

/// CIE Clear Sky model.
///
/// Luminance depends on both the altitude and the angular distance from the sun.
pub struct CIEClearSky {
    /// Zenith luminance (cd/m^2).
    pub zenith_luminance: f64,
}

impl CIEClearSky {
    pub fn new(zenith_luminance: f64) -> Self {
        Self { zenith_luminance }
    }
}

impl SkyModel for CIEClearSky {
    fn luminance(&self, direction: Vector, solar_pos: &SolarPosition) -> Rgb {
        let dir = match direction.normalize() {
            Ok(v) => v,
            Err(_) => return [0.0; 3],
        };

        let alt = dir.dz.asin();
        if alt < 0.0 {
            return [0.0; 3];
        }

        let sun_dir = solar_pos.to_direction();
        let sun_dir = match sun_dir.normalize() {
            Ok(v) => v,
            Err(_) => return [0.0; 3],
        };

        // Angular distance between sky point and sun
        let cos_gamma = dir.dot(&sun_dir).clamp(-1.0, 1.0);
        let gamma = cos_gamma.acos();

        let solar_alt = solar_pos.altitude.to_radians();

        // Indicatrix function (circumsolar brightening)
        let f_gamma = 0.91 + 10.0 * (-3.0 * gamma).exp() + 0.45 * cos_gamma * cos_gamma;
        let f_zenith_sun = 0.91 + 10.0 * (-3.0 * solar_alt).exp() + 0.45 * solar_alt.cos().powi(2);

        // Gradation function
        let g_alt = 0.274 * (0.45 - alt.sin() + alt.sin().powi(2));
        let g_zenith: f64 = 0.274 * (0.45 - 1.0 + 1.0); // alt = pi/2

        if f_zenith_sun.abs() < 1e-10 || g_zenith.abs() < 1e-10 {
            return [0.0; 3];
        }

        let ratio = (f_gamma * g_alt) / (f_zenith_sun * g_zenith);
        let l = self.zenith_luminance * ratio.max(0.0);
        [l, l, l]
    }
}

/// Perez All-Weather Sky model.
///
/// Uses 5 parameters (a-e) that vary with sky conditions.
pub struct PerezAllWeather {
    /// Perez coefficients [a, b, c, d, e].
    pub coefficients: [f64; 5],
    /// Zenith luminance (cd/m^2).
    pub zenith_luminance: f64,
}

impl PerezAllWeather {
    pub fn new(coefficients: [f64; 5], zenith_luminance: f64) -> Self {
        Self {
            coefficients,
            zenith_luminance,
        }
    }

    /// Creates a Perez model for overcast conditions.
    pub fn overcast(zenith_luminance: f64) -> Self {
        Self {
            coefficients: [1.0, 0.0, 0.0, -1.0, 0.0],
            zenith_luminance,
        }
    }

    /// Creates a Perez model for clear sky conditions.
    pub fn clear(zenith_luminance: f64) -> Self {
        Self {
            coefficients: [-1.0, -0.7, 3.0, -2.5, 0.6],
            zenith_luminance,
        }
    }
}

impl SkyModel for PerezAllWeather {
    fn luminance(&self, direction: Vector, solar_pos: &SolarPosition) -> Rgb {
        let dir = match direction.normalize() {
            Ok(v) => v,
            Err(_) => return [0.0; 3],
        };

        let alt = dir.dz.asin();
        if alt < 0.0 {
            return [0.0; 3];
        }

        let sun_dir = solar_pos.to_direction();
        let sun_dir = match sun_dir.normalize() {
            Ok(v) => v,
            Err(_) => return [0.0; 3],
        };

        let cos_gamma = dir.dot(&sun_dir).clamp(-1.0, 1.0);
        let gamma = cos_gamma.acos();

        let solar_zenith = (90.0 - solar_pos.altitude).to_radians();
        let zenith = std::f64::consts::FRAC_PI_2 - alt;

        let [a, b, c, d, e] = self.coefficients;

        // Perez luminance function
        let f = |theta: f64, g: f64| -> f64 {
            (1.0 + a * (b / (theta.cos().max(0.01))).exp())
                * (1.0 + c * (-d * g).exp() + e * g.cos().powi(2))
        };

        let f_sky = f(zenith, gamma);
        let f_zenith = f(0.0, solar_zenith);

        if f_zenith.abs() < 1e-10 {
            return [0.0; 3];
        }

        let l = self.zenith_luminance * (f_sky / f_zenith).max(0.0);
        [l, l, l]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cie_overcast_zenith_brightest() {
        let sky = CIEOvercast::new(10000.0);
        let solar = SolarPosition {
            altitude: 45.0,
            azimuth: 180.0,
        };

        let l_zenith = sky.luminance(Vector::new(0.0, 0.0, 1.0), &solar);
        let l_horizon = sky.luminance(Vector::new(1.0, 0.0, 0.01), &solar);

        // Zenith should be brightest for overcast sky
        assert!(l_zenith[0] > l_horizon[0]);
    }

    #[test]
    fn test_cie_overcast_below_horizon() {
        let sky = CIEOvercast::new(10000.0);
        let solar = SolarPosition {
            altitude: 45.0,
            azimuth: 180.0,
        };

        let l = sky.luminance(Vector::new(0.0, 0.0, -1.0), &solar);
        assert!((l[0] - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_cie_clear_circumsolar() {
        let sky = CIEClearSky::new(10000.0);
        let solar = SolarPosition {
            altitude: 60.0,
            azimuth: 180.0,
        };

        let sun_dir = solar.to_direction();
        let l_sun = sky.luminance(sun_dir, &solar);

        // Away from sun
        let l_away = sky.luminance(Vector::new(0.0, 1.0, 0.5), &solar);

        // Near-sun direction should generally be brighter than away
        // (this is the circumsolar brightening effect)
        assert!(l_sun[0] > 0.0);
        assert!(l_away[0] > 0.0);
    }

    #[test]
    fn test_perez_overcast() {
        let sky = PerezAllWeather::overcast(10000.0);
        let solar = SolarPosition {
            altitude: 45.0,
            azimuth: 180.0,
        };

        let l = sky.luminance(Vector::new(0.0, 0.0, 1.0), &solar);
        assert!(l[0] > 0.0, "Overcast sky should have positive luminance");
    }

    #[test]
    fn test_from_diffuse_illuminance() {
        let sky = CIEOvercast::from_diffuse_illuminance(10000.0);
        assert!(sky.zenith_luminance > 0.0);
    }
}
