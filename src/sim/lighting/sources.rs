use crate::{Point, Vector};

/// RGB color/intensity triplet.
pub type Rgb = [f64; 3];

/// Trait for all light sources.
pub trait LightSource {
    /// Returns the position of the light source (or representative point).
    fn position(&self) -> Point;

    /// Returns the intensity [R, G, B] (in lumens or watts) at a given direction from the source.
    fn intensity(&self, direction: Vector) -> Rgb;

    /// Total luminous flux (lumens) of the source.
    fn total_flux(&self) -> f64;
}

/// An omnidirectional point light source.
pub struct PointLight {
    pub position: Point,
    /// RGB intensity (lumens per steradian in each channel).
    pub intensity: Rgb,
}

impl PointLight {
    pub fn new(position: Point, intensity: Rgb) -> Self {
        Self {
            position,
            intensity,
        }
    }

    /// Creates a white point light with given total luminous flux (lumens).
    pub fn white(position: Point, lumens: f64) -> Self {
        // Distribute equally across RGB, over full sphere (4*pi steradians)
        let i = lumens / (4.0 * std::f64::consts::PI);
        Self {
            position,
            intensity: [i, i, i],
        }
    }
}

impl LightSource for PointLight {
    fn position(&self) -> Point {
        self.position
    }

    fn intensity(&self, _direction: Vector) -> Rgb {
        self.intensity
    }

    fn total_flux(&self) -> f64 {
        (self.intensity[0] + self.intensity[1] + self.intensity[2]) * 4.0 * std::f64::consts::PI
            / 3.0
    }
}

/// A rectangular area light source.
pub struct AreaLight {
    /// Center position.
    pub position: Point,
    /// Normal direction (emission direction).
    pub normal: Vector,
    /// Width of the area light.
    pub width: f64,
    /// Height of the area light.
    pub height: f64,
    /// RGB intensity per unit area.
    pub intensity: Rgb,
}

impl AreaLight {
    pub fn new(position: Point, normal: Vector, width: f64, height: f64, intensity: Rgb) -> Self {
        Self {
            position,
            normal,
            width,
            height,
            intensity,
        }
    }
}

impl LightSource for AreaLight {
    fn position(&self) -> Point {
        self.position
    }

    fn intensity(&self, direction: Vector) -> Rgb {
        // Lambertian emission: intensity * cos(theta)
        let n = match self.normal.normalize() {
            Ok(v) => v,
            Err(_) => return [0.0; 3],
        };
        let d = match direction.normalize() {
            Ok(v) => v,
            Err(_) => return [0.0; 3],
        };
        let cos_theta = n.dot(&d).abs();
        [
            self.intensity[0] * cos_theta,
            self.intensity[1] * cos_theta,
            self.intensity[2] * cos_theta,
        ]
    }

    fn total_flux(&self) -> f64 {
        let area = self.width * self.height;
        // Lambertian emitter: flux = pi * intensity * area
        (self.intensity[0] + self.intensity[1] + self.intensity[2]) / 3.0
            * std::f64::consts::PI
            * area
    }
}

/// A directional (parallel) light source, such as sunlight.
pub struct DirectionalLight {
    /// Direction the light travels (from source toward scene).
    pub direction: Vector,
    /// RGB irradiance (W/m^2 or lux in each channel).
    pub irradiance: Rgb,
}

impl DirectionalLight {
    pub fn new(direction: Vector, irradiance: Rgb) -> Self {
        Self {
            direction,
            irradiance,
        }
    }
}

impl LightSource for DirectionalLight {
    fn position(&self) -> Point {
        // Directional lights don't have a position; return origin
        Point::new(0.0, 0.0, 0.0)
    }

    fn intensity(&self, _direction: Vector) -> Rgb {
        self.irradiance
    }

    fn total_flux(&self) -> f64 {
        // Infinite for truly directional lights
        f64::INFINITY
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_point_light_white() {
        let light = PointLight::white(Point::new(0.0, 0.0, 0.0), 1000.0);
        let i = light.intensity(Vector::new(1.0, 0.0, 0.0));
        assert!(i[0] > 0.0);
        assert!((i[0] - i[1]).abs() < 1e-10);
        assert!((i[1] - i[2]).abs() < 1e-10);
    }

    #[test]
    fn test_area_light_lambertian() {
        let light = AreaLight::new(
            Point::new(0.0, 0.0, 0.0),
            Vector::new(0.0, 0.0, -1.0),
            1.0,
            1.0,
            [100.0, 100.0, 100.0],
        );
        // Normal direction should give full intensity
        let i_normal = light.intensity(Vector::new(0.0, 0.0, -1.0));
        // 90 degree angle should give zero
        let i_perp = light.intensity(Vector::new(1.0, 0.0, 0.0));
        assert!(i_normal[0] > i_perp[0]);
        assert!(i_perp[0] < 1e-10);
    }

    #[test]
    fn test_directional_light() {
        let light = DirectionalLight::new(Vector::new(0.0, 0.0, -1.0), [500.0, 500.0, 500.0]);
        let i = light.intensity(Vector::new(1.0, 0.0, 0.0));
        assert!((i[0] - 500.0).abs() < 1e-10);
        assert!(light.total_flux().is_infinite());
    }
}
