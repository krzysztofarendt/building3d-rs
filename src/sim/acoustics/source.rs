use crate::Vector;

/// Directivity pattern for sound sources.
pub trait SourceDirectivity {
    /// Returns the energy gain factor for a given emission direction.
    ///
    /// The direction is relative to the source's forward axis.
    /// Returns a value >= 0, where 1.0 is the reference level.
    fn gain(&self, direction: Vector, forward: Vector) -> f64;
}

/// Omnidirectional source: equal energy in all directions.
pub struct Omnidirectional;

impl SourceDirectivity for Omnidirectional {
    fn gain(&self, _direction: Vector, _forward: Vector) -> f64 {
        1.0
    }
}

/// Cardioid directivity pattern: gain = 0.5 * (1 + cos(theta)).
///
/// Maximum gain (1.0) in the forward direction, zero gain behind.
pub struct Cardioid;

impl SourceDirectivity for Cardioid {
    fn gain(&self, direction: Vector, forward: Vector) -> f64 {
        let dir_norm = match direction.normalize() {
            Ok(v) => v,
            Err(_) => return 0.0,
        };
        let fwd_norm = match forward.normalize() {
            Ok(v) => v,
            Err(_) => return 1.0,
        };
        let cos_theta = dir_norm.dot(&fwd_norm);
        0.5 * (1.0 + cos_theta)
    }
}

/// Subcardioid directivity: gain = 0.75 + 0.25 * cos(theta).
///
/// More omnidirectional than cardioid but still directional.
pub struct SubCardioid;

impl SourceDirectivity for SubCardioid {
    fn gain(&self, direction: Vector, forward: Vector) -> f64 {
        let dir_norm = match direction.normalize() {
            Ok(v) => v,
            Err(_) => return 0.0,
        };
        let fwd_norm = match forward.normalize() {
            Ok(v) => v,
            Err(_) => return 1.0,
        };
        let cos_theta = dir_norm.dot(&fwd_norm);
        0.75 + 0.25 * cos_theta
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_omnidirectional() {
        let omni = Omnidirectional;
        let fwd = Vector::new(1.0, 0.0, 0.0);
        assert!((omni.gain(Vector::new(1.0, 0.0, 0.0), fwd) - 1.0).abs() < 1e-10);
        assert!((omni.gain(Vector::new(-1.0, 0.0, 0.0), fwd) - 1.0).abs() < 1e-10);
        assert!((omni.gain(Vector::new(0.0, 1.0, 0.0), fwd) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_cardioid_forward() {
        let card = Cardioid;
        let fwd = Vector::new(1.0, 0.0, 0.0);
        // Forward direction: cos(0) = 1, gain = 0.5*(1+1) = 1.0
        let g = card.gain(Vector::new(1.0, 0.0, 0.0), fwd);
        assert!((g - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_cardioid_backward() {
        let card = Cardioid;
        let fwd = Vector::new(1.0, 0.0, 0.0);
        // Backward direction: cos(pi) = -1, gain = 0.5*(1-1) = 0.0
        let g = card.gain(Vector::new(-1.0, 0.0, 0.0), fwd);
        assert!(g.abs() < 1e-10);
    }

    #[test]
    fn test_cardioid_perpendicular() {
        let card = Cardioid;
        let fwd = Vector::new(1.0, 0.0, 0.0);
        // Perpendicular: cos(pi/2) = 0, gain = 0.5*(1+0) = 0.5
        let g = card.gain(Vector::new(0.0, 1.0, 0.0), fwd);
        assert!((g - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_subcardioid() {
        let sub = SubCardioid;
        let fwd = Vector::new(1.0, 0.0, 0.0);
        // Forward: 0.75 + 0.25 = 1.0
        assert!((sub.gain(fwd, fwd) - 1.0).abs() < 1e-10);
        // Backward: 0.75 - 0.25 = 0.5
        let g = sub.gain(Vector::new(-1.0, 0.0, 0.0), fwd);
        assert!((g - 0.5).abs() < 1e-10);
    }
}
