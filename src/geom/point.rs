use crate::Vector;
use crate::geom::EPS;
use std::fmt;
use std::ops::Add;

pub mod check;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Point {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Returns true if both points are very close to each other.
    pub fn is_close(&self, other: &Self) -> bool {
        (self.x - other.x).abs() < EPS
            && (self.y - other.y).abs() < EPS
            && (self.z - other.z).abs() < EPS
    }

    /// Multiplies all coordinates by a scalar and returns a copy.
    pub fn scale(&self, scale: f64) -> Self {
        Self {
            x: self.x * scale,
            y: self.y * scale,
            z: self.z * scale,
        }
    }
}

impl fmt::Display for Point {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let prec = f.precision().unwrap_or(2); // Default 2 decimals
        write!(
            f,
            "Point({:.prec$}, {:.prec$}, {:.prec$})",
            self.x,
            self.y,
            self.z,
            prec = prec
        )
    }
}

// Implement +
// (Sub is NOT implemented)
impl Add<Vector> for Point {
    type Output = Point;
    fn add(self, other: Vector) -> Self {
        Self {
            x: self.x + other.dx,
            y: self.y + other.dy,
            z: self.z + other.dz,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_is_close() {
        let pa = Point::new(5., 5., 5.);
        let pb = Point::new(5.00000000000001, 5., 5.);
        let pc = Point::new(5.0001, 5., 5.);
        assert!(pa.is_close(&pb));
        assert!(!pa.is_close(&pc));
    }

    #[test]
    fn test_scale() {
        let p1 = Point::new(1., 2., 3.);
        let p2 = p1.scale(10.);
        assert!(p2.is_close(&Point::new(10., 20., 30.)));
    }
}
