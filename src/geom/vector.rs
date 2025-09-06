use crate::Point;
use crate::geom::EPS;
use std::fmt;
use std::ops::{Add, Div, Mul, Sub};

pub mod check;

#[derive(Debug, Clone, Copy, PartialEq)]
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
    pub fn cross(self, other: Self) -> Self {
        Self {
            dx: self.dy * other.dz - self.dz * other.dy,
            dy: self.dz * other.dx - self.dx * other.dz,
            dz: self.dx * other.dy - self.dy * other.dx,
        }
    }

    /// Dot product between 2 vectors.
    pub fn dot(self, other: Self) -> f64 {
        self.dx * other.dx + self.dy * other.dy + self.dz * other.dz
    }

    /// Returns the length of the vector.
    pub fn length(&self) -> f64 {
        (self.dx.powi(2) + self.dy.powi(2) + self.dz.powi(2)).sqrt()
    }

    pub fn is_close(&self, other: &Self) -> bool {
        (self.dx - other.dx).abs() < EPS
            && (self.dy - other.dy).abs() < EPS
            && (self.dz - other.dz).abs() < EPS
    }

    /// Normalizes the vector (divides by its length) and returns a copy.
    pub fn normalize(&self) -> Option<Self> {
        let len = self.length();
        if len < EPS {
            None
        } else {
            Some(Self {
                dx: self.dx / len,
                dy: self.dy / len,
                dz: self.dz / len,
            })
        }
    }

    /// Calculates vector normal to the surface defined with 3 points.
    ///
    /// If the normal does not exist, returns None.
    /// The normal does not exist if the points are collinear.
    pub fn normal(pt0: Point, pt1: Point, pt2: Point) -> Option<Self> {
        let v01 = Self::from_points(pt0, pt1);
        let v02 = Self::from_points(pt0, pt2);
        let vn = v01.cross(v02);
        vn.normalize()
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

// Implement /
impl Div<f64> for Vector {
    type Output = Option<Self>;
    fn div(self, other: f64) -> Option<Self> {
        if other < EPS {
            None
        } else {
            Some(Self {
                dx: self.dx / other,
                dy: self.dy / other,
                dz: self.dz / other,
            })
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
        let v_cross = vx.cross(vy);
        assert_eq!(v_cross, Vector::new(0., 0., 1.));
        let len = v_cross.length();
        assert_eq!(len, 1.);
    }

    #[test]
    fn test_normalize() {
        // Non-zero-length vector
        let v = Vector::new(9., 0., 0.);
        let vnorm = v.normalize();
        assert!(vnorm.is_some());
        assert_eq!(vnorm.unwrap(), Vector::new(1., 0., 0.));
        // Zero-length vector
        let v = Vector::new(0., 0., 0.);
        assert!(v.normalize().is_none());
    }

    #[test]
    fn test_normal() {
        let p0 = Point::new(1., 0., 0.);
        let p1 = Point::new(1., 1., 0.);
        let p2 = Point::new(0., 1., 0.);
        let vn = Vector::normal(p0, p1, p2);
        assert!(vn.is_some());
        let expected = Vector::new(0., 0., 1.);
        assert!(vn.unwrap().is_close(&expected));
    }
}
