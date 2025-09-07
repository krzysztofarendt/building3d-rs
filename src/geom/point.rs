use crate::Vector;
use crate::geom::{EPS, IsClose};
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
        self.x.is_close(&other.x) && self.y.is_close(&other.y) && self.z.is_close(&other.z)
    }

    /// Multiplies all coordinates by a scalar and returns a copy.
    pub fn scale(&self, scale: f64) -> Self {
        Self {
            x: self.x * scale,
            y: self.y * scale,
            z: self.z * scale,
        }
    }

    // Creates a new point along the edge pt1->pt2 with some relative distance from pt1.
    pub fn new_between_2_points(pt1: Self, pt2: Self, rel_d: f64) -> Self {
        let alpha_v = Vector::new(rel_d, rel_d, rel_d);
        let v_pt1 = Vector::from_a_point(pt1);
        let v_pt2 = Vector::from_a_point(pt2);
        let new_v = v_pt1 * (1. - alpha_v) + v_pt2 * alpha_v;
        Self::new(new_v.dx, new_v.dy, new_v.dz)
    }

    /// Creates `num` new points along the edge spanning from `pt1` to `pt2`.
    ///
    /// Points are uniformly distributed based on `num`.
    /// Returns a vector of all points: `[pt1, new_1, new_2, ..., new_num, pt2]`
    pub fn many_new_between_2_points(pt1: Self, pt2: Self, num: usize) -> Vec<Point> {
        let mut vecs: Vec<Point> = vec![pt1];
        let vp1 = Vector::from_a_point(pt1);
        let vp2 = Vector::from_a_point(pt2);

        for i in 1..num + 1 {
            let alpha = (i as f64) / (num + 1) as f64;
            let alpha_v = Vector::new(alpha, alpha, alpha);
            let new_vec = vp1 * (1. - alpha_v) + vp2 * alpha_v;
            let new_pt = Point::new(new_vec.dx, new_vec.dy, new_vec.dz);
            vecs.push(new_pt);
        }
        vecs.push(pt2);

        vecs
    }

    /// Checks if the point lies on the line segment defined by points `pt1` and `pt2`.
    pub fn is_on_segment(self, pt1: Self, pt2: Self) -> bool {
        if self.is_close(&pt1) || self.is_close(&pt2) {
            return true;
        }
        if !check::are_points_collinear(&[self, pt1, pt2]) {
            return false;
        }
        // Check if the point lies within the segment bounds
        let v_self = Vector::from_a_point(self);
        let v_pt1 = Vector::from_a_point(pt1);
        let v_seg = Vector::from_points(pt1, pt2);
        let dot_prod = (v_self - v_pt1).dot(&v_seg);
        let sq_len_seg = v_seg.length().powi(2);

        if dot_prod < -EPS || dot_prod > (sq_len_seg + EPS) {
            return false;
        }

        true
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

    #[test]
    fn test_new_between_2_points() {
        let p0 = Point::new(0., 0., 0.);
        let p1 = Point::new(1., 1., 1.);
        let ptest = Point::new_between_2_points(p0, p1, 0.5);
        assert!(ptest.is_close(&Point::new(0.5, 0.5, 0.5)));
        let ptest = Point::new_between_2_points(p0, p1, 0.0);
        assert!(ptest.is_close(&p0));
        let ptest = Point::new_between_2_points(p0, p1, 1.0);
        assert!(ptest.is_close(&p1));
        let ptest = Point::new_between_2_points(p0, p0, 0.5);
        assert!(ptest.is_close(&p0));
        let ptest = Point::new_between_2_points(p0, p0, 0.0);
        assert!(ptest.is_close(&p0));
        let ptest = Point::new_between_2_points(p0, p0, 1.0);
        assert!(ptest.is_close(&p0));
    }

    #[test]
    fn test_many_new_between_2_points() {
        let p_beg = Point::new(0., 0., 0.);
        let p_end = Point::new(1., 1., 1.);
        let num: usize = 3;
        let pt_all = Point::many_new_between_2_points(p_beg, p_end, num);
        let expected = vec![
            p_beg,
            Point::new(0.25, 0.25, 0.25),
            Point::new(0.5, 0.5, 0.5),
            Point::new(0.75, 0.75, 0.75),
            p_end,
        ];
        // Check if output points are same as expected
        for (pa, pb) in pt_all.iter().zip(expected) {
            assert!(pa.is_close(&pb));
        }
    }

    #[test]
    fn test_is_on_segment() {
        let p_beg = Point::new(0., 0., 0.);
        let p_end = Point::new(1., 1., 1.);
        // True
        let p_test = Point::new(0.5, 0.5, 0.5);
        assert!(p_test.is_on_segment(p_beg, p_end));
        let p_test = Point::new(0.0, 0.0, 0.0);
        assert!(p_test.is_on_segment(p_beg, p_end));
        let p_test = Point::new(1.0, 1.0, 1.0);
        assert!(p_test.is_on_segment(p_beg, p_end));
        // False
        let p_test = Point::new(0.6, 0.5, 0.5);
        assert!(!p_test.is_on_segment(p_beg, p_end));
        let p_test = Point::new(-0.1, -0.1, -0.1);
        assert!(!p_test.is_on_segment(p_beg, p_end));
        let p_test = Point::new(1.1, 1.1, 1.1);
        assert!(!p_test.is_on_segment(p_beg, p_end));
    }
}
