use super::*;
use crate::Vector;
use crate::geom::vector::check::are_vectors_close;

pub fn are_points_coplanar(pts: Vec<Point>) -> bool {
    if pts.len() <= 3 {
        return true;
    }
    // TODO: Finish. It will require a dot product of 2 vectors.
    true
}

/// Checks if (multiple) points are collinear
pub fn are_points_collinear(pts: &[Point]) -> bool {
    if pts.len() <= 2 {
        return true; // 1 or 2 points are always collinear
    }
    // Calculate unit vectors using the first point as origin
    let mut unit_vectors: Vec<Vector> = Vec::new();
    for p in pts.iter().skip(1) {
        let v = Vector::from_points(pts[0], *p)
            .normalize()
            .unwrap_or(Vector::new(0., 0., 0.));
        unit_vectors.push(v);
    }
    // Unit vectors must be equal if points are collinear
    are_vectors_close(&unit_vectors)
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_are_points_collinear() {
        // Collinear
        let p0 = Point::new(0., 0., 0.);
        let p1 = Point::new(1., 0., 0.);
        let p2 = Point::new(2., 0., 0.);
        assert!(are_points_collinear(&[p0]));
        assert!(are_points_collinear(&[p0, p1]));
        assert!(are_points_collinear(&[p0, p1, p2]));
        // Not collinear
        let p3 = Point::new(1., 1., 0.);
        assert!(!are_points_collinear(&[p0, p1, p2, p3]));
    }
}
