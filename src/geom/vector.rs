use crate::Point;
use crate::geom::EPS;
use crate::geom::IsClose;
use anyhow::{Result, anyhow};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;
use std::fmt;
use std::ops::{Add, Div, Mul, Sub};

pub mod check;

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Vector {
    pub dx: f64,
    pub dy: f64,
    pub dz: f64,
}

impl Vector {
    pub fn new(dx: f64, dy: f64, dz: f64) -> Self {
        Self { dx, dy, dz }
    }

    pub fn from_a_point(pt: Point) -> Self {
        Self::from_points(Point::new(0., 0., 0.), pt)
    }

    pub fn from_points(beg: Point, end: Point) -> Self {
        Self {
            dx: end.x - beg.x,
            dy: end.y - beg.y,
            dz: end.z - beg.z,
        }
    }

    /// Cross product between 2 vectors.
    pub fn cross(&self, other: &Self) -> Self {
        Self {
            dx: self.dy * other.dz - self.dz * other.dy,
            dy: self.dz * other.dx - self.dx * other.dz,
            dz: self.dx * other.dy - self.dy * other.dx,
        }
    }

    /// Dot product between 2 vectors.
    pub fn dot(&self, other: &Self) -> f64 {
        self.dx * other.dx + self.dy * other.dy + self.dz * other.dz
    }

    /// Returns the length of the vector.
    pub fn length(&self) -> f64 {
        (self.dx.powi(2) + self.dy.powi(2) + self.dz.powi(2)).sqrt()
    }

    pub fn is_close(&self, other: &Self) -> bool {
        self.dx.is_close(other.dx) && self.dy.is_close(other.dy) && self.dz.is_close(other.dz)
    }

    /// Normalizes the vector (divides by its length) and returns a copy.
    pub fn normalize(&self) -> Result<Self> {
        let len = self.length();
        if len < EPS {
            Err(anyhow!("Cannot normalize a zero-length vector"))
        } else {
            Ok(Self {
                dx: self.dx / len,
                dy: self.dy / len,
                dz: self.dz / len,
            })
        }
    }

    /// Calculates vector normal to the surface defined with 3 points.
    ///
    /// If the normal does not exist, returns Err.
    /// The normal does not exist if the points are collinear.
    pub fn normal(p0: Point, p1: Point, p2: Point) -> Result<Self> {
        let v01 = Self::from_points(p0, p1);
        let v02 = Self::from_points(p0, p2);
        let vn = v01.cross(&v02);
        vn.normalize()
    }

    /// Calculates angle in radians between two vectors.
    ///
    /// The output range is within 0 and PI.
    /// Returns None if any of the vectors is zero-length.
    pub fn angle(&self, other: &Self) -> Option<f64> {
        let self_opt_norm = self.normalize();
        let other_opt_norm = other.normalize();
        let self_norm: Vector;
        let other_norm: Vector;
        if let (Ok(v1), Ok(v2)) = (self_opt_norm, other_opt_norm) {
            self_norm = v1;
            other_norm = v2;
        } else {
            return None;
        }
        let dot = self_norm.dot(&other_norm).clamp(-1.0, 1.0);
        let rad = dot.acos();
        assert!((0.0..=PI).contains(&rad));
        Some(rad)
    }
}

impl fmt::Display for Vector {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let prec = f.precision().unwrap_or(2); // Default 2 decimals
        write!(
            f,
            "Vector({:.prec$}, {:.prec$}, {:.prec$})",
            self.dx,
            self.dy,
            self.dz,
            prec = prec
        )
    }
}

// Implement +
impl Add for Vector {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self {
            dx: self.dx + other.dx,
            dy: self.dy + other.dy,
            dz: self.dz + other.dz,
        }
    }
}
impl Add<f64> for Vector {
    type Output = Self;
    fn add(self, other: f64) -> Self {
        Self {
            dx: self.dx + other,
            dy: self.dy + other,
            dz: self.dz + other,
        }
    }
}
impl Add<Vector> for f64 {
    type Output = Vector;
    fn add(self, other: Vector) -> Vector {
        Vector {
            dx: self + other.dx,
            dy: self + other.dy,
            dz: self + other.dz,
        }
    }
}

// Implement +
impl Sub for Vector {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        Self {
            dx: self.dx - other.dx,
            dy: self.dy - other.dy,
            dz: self.dz - other.dz,
        }
    }
}
impl Sub<f64> for Vector {
    type Output = Self;
    fn sub(self, other: f64) -> Self {
        Self {
            dx: self.dx - other,
            dy: self.dy - other,
            dz: self.dz - other,
        }
    }
}
impl Sub<Vector> for f64 {
    type Output = Vector;
    fn sub(self, other: Vector) -> Vector {
        Vector {
            dx: self - other.dx,
            dy: self - other.dy,
            dz: self - other.dz,
        }
    }
}

// Implement *
impl Mul<f64> for Vector {
    type Output = Self;
    fn mul(self, other: f64) -> Self {
        Self {
            dx: self.dx * other,
            dy: self.dy * other,
            dz: self.dz * other,
        }
    }
}
impl Mul<Vector> for f64 {
    type Output = Vector;
    fn mul(self, other: Vector) -> Vector {
        Vector {
            dx: self * other.dx,
            dy: self * other.dy,
            dz: self * other.dz,
        }
    }
}

// Implement /
impl Div<f64> for Vector {
    type Output = Self;
    fn div(self, other: f64) -> Self {
        Self {
            dx: self.dx / other,
            dy: self.dy / other,
            dz: self.dz / other,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_points() {
        let p0 = Point::new(1., 1., 1.);
        let p1 = Point::new(0., 0., 0.);
        let va = Vector::from_points(p0, p1);
        let vb = Vector::from_points(p1, p0);
        assert_eq!(va, vb * -1.);
    }

    #[test]
    fn test_cross() {
        let vx = Vector::new(1., 0., 0.);
        let vy = Vector::new(0., 1., 0.);
        let v_cross = vx.cross(&vy);
        assert_eq!(v_cross, Vector::new(0., 0., 1.));
        let len = v_cross.length();
        assert_eq!(len, 1.);
    }

    #[test]
    fn test_normalize() {
        // Non-zero-length vector
        let v = Vector::new(9., 0., 0.);
        let vnorm = v.normalize();
        assert!(vnorm.is_ok());
        assert_eq!(vnorm.unwrap(), Vector::new(1., 0., 0.));
        // Zero-length vector
        let v = Vector::new(0., 0., 0.);
        assert!(v.normalize().is_err());
    }

    #[test]
    fn test_normal() {
        let p0 = Point::new(1., 0., 0.);
        let p1 = Point::new(1., 1., 0.);
        let p2 = Point::new(0., 1., 0.);
        let vn = Vector::normal(p0, p1, p2);
        assert!(vn.is_ok());
        let expected = Vector::new(0., 0., 1.);
        assert!(vn.unwrap().is_close(&expected));
    }

    #[test]
    fn test_ops() {
        let v1 = Vector::new(1., 1., 1.);
        let v2 = Vector::new(0., 0., 0.);
        let _ = v1 + v2;
        let _ = v1 - v2;
        let _ = v1 * 2.;
        let _ = v1 / 2.;
        let _ = 2. * v1;
        let result = v1 / 2.;
        assert!(result.is_close(&Vector::new(0.5, 0.5, 0.5)));
    }

    #[test]
    fn test_angle() {
        let vx = Vector::new(1., 0., 0.); // 3 orthogonal vectors
        let vy = Vector::new(0., 1., 0.); // 3 orthogonal vectors
        let vz = Vector::new(0., 0., 1.); // 3 orthogonal vectors
        let vxy = Vector::new(1., 1., 0.); // 45 degrees w.r.t. vx
        let deg90 = PI / 2.0;
        let deg45 = PI / 4.0;
        let angle = vx.angle(&vy).unwrap();
        assert!(angle.is_close(deg90));
        let angle = vx.angle(&vz).unwrap();
        assert!(angle.is_close(deg90));
        let angle = vy.angle(&vz).unwrap();
        assert!(angle.is_close(deg90));
        let angle = vx.angle(&vxy).unwrap();
        assert!(angle.is_close(deg45));
    }
}
