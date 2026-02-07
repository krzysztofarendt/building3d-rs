use crate::Point;
use crate::Vector;

/// Orthonormal basis for projecting 3D points onto a 2D plane and back.
#[derive(Debug, Clone, Copy)]
pub struct PlaneBasis {
    pub origin: Point,
    pub u: Vector,
    pub v: Vector,
}

impl PlaneBasis {
    /// Creates a `PlaneBasis` from the first three vertices of a polygon.
    pub fn from_polygon(poly: &[Point]) -> Option<Self> {
        if poly.len() < 3 {
            return None;
        }

        let v1 = poly[1] - poly[0];
        let v2 = poly[2] - poly[0];
        let n = v1.cross(&v2).normalize().ok()?;

        Self::build(poly[0], n)
    }

    /// Creates a `PlaneBasis` from an origin point and a normal vector.
    pub fn from_normal(origin: Point, normal: Vector) -> Option<Self> {
        let n = normal.normalize().ok()?;
        Self::build(origin, n)
    }

    fn build(origin: Point, n: Vector) -> Option<Self> {
        let helper = if n.dz.abs() < 0.9 {
            Vector::new(0.0, 0.0, 1.0)
        } else {
            Vector::new(0.0, 1.0, 0.0)
        };

        let u = helper.cross(&n).normalize().ok()?;
        let v = n.cross(&u).normalize().ok()?;

        Some(Self { origin, u, v })
    }

    /// Projects a 3D point onto the 2D plane, returning (u, v) coordinates.
    pub fn project(&self, p: Point) -> (f64, f64) {
        let r = p - self.origin;
        (r.dot(&self.u), r.dot(&self.v))
    }

    /// Unprojects 2D (u, v) coordinates back to a 3D point on the plane.
    pub fn unproject(&self, x: f64, y: f64) -> Point {
        self.origin + self.u * x + self.v * y
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_polygon_roundtrip() {
        let pts = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
        ];
        let basis = PlaneBasis::from_polygon(&pts).unwrap();
        for &p in &pts {
            let (u, v) = basis.project(p);
            let back = basis.unproject(u, v);
            assert!(p.is_close(&back));
        }
    }

    #[test]
    fn test_from_normal_roundtrip() {
        let origin = Point::new(1.0, 2.0, 3.0);
        let normal = Vector::new(0.0, 0.0, 1.0);
        let basis = PlaneBasis::from_normal(origin, normal).unwrap();
        let p = Point::new(2.0, 3.0, 3.0);
        let (u, v) = basis.project(p);
        let back = basis.unproject(u, v);
        assert!(p.is_close(&back));
    }
}
