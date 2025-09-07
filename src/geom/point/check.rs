use super::*;
use crate::Vector;
use crate::geom::EPS;

/// Checks if all points are (almost) equal
pub fn are_points_close(pts: &[Point]) -> bool {
    let mut all_close = true;
    let p0 = pts[0];
    for p in pts.iter().skip(1) {
        if !p.is_close(&p0) {
            all_close = false;
            break;
        }
    }
    all_close
}

/// Checks if (multiple) points are coplanar.
pub fn are_points_coplanar(pts: &[Point]) -> bool {
    if pts.len() <= 3 {
        return true;
    }

    // Calculate normal vector based on 3 points from pts,
    // however some points may be collinear, so we need to search
    // for a triplet of non-collinear points.
    let triplet = first_3_noncollinear(pts);
    if triplet.is_none() {
        // All points must be collinear
        // Collinear points are also coplanar
        return true;
    }
    let (i, j, k) = triplet.unwrap();

    // I can unwrap, because I already know the points are not collinear
    // so the normal must exist
    let vn = Vector::normal(pts[i], pts[j], pts[k]).unwrap();

    // Plane equation:
    // ax + by + cz + d = 0
    // (a, b, c) are taken from the normal vector
    // d is calculated by substituting one of the points
    let r: usize = 0; // Reference point index
    let pt_r = pts[r];
    let d = -vn.dot(&Vector::from_a_point(pt_r));

    // Check if all points lay on the same plane
    for pt in pts.iter() {
        let coplanar: bool = (d + vn.dot(&Vector::from_a_point(*pt))).abs() < EPS;
        if !coplanar {
            return false;
        }
    }
    true
}

/// Returns indices of the first 3 noncollinear points from the slice.
///
/// If all points are collinear, returns None.
/// Complexity: O(n^2)
fn first_3_noncollinear(pts: &[Point]) -> Option<(usize, usize, usize)> {
    // It works by first finding a pair of different points (i, j)
    // and then searching for the third noncollinear point (k).
    let n = pts.len();
    for i in 0..n {
        for j in i + 1..n {
            if pts[i].is_close(&pts[j]) {
                continue;
            }
            for k in j + 1..n {
                if !are_points_collinear(&[pts[i], pts[j], pts[k]]) {
                    return Some((i, j, k));
                }
            }
        }
    }
    None
}

/// Checks if (multiple) points are collinear.
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
    // Unit vectors must be equal or exactly negative if points are collinear
    let mut are_collinear = true;
    let v0 = &unit_vectors[0];
    let v0_rev = &(-1. * unit_vectors[0]);
    for uv in unit_vectors.iter().skip(1) {
        if !(uv.is_close(v0) || uv.is_close(v0_rev)) {
            are_collinear = false;
            break;
        }
    }

    are_collinear
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_are_points_close() {
        // Close
        let p0 = Point::new(0., 0., 0.);
        let p1 = Point::new(EPS * 0.5, 0., 0.);
        let p2 = Point::new(EPS * 0.5, EPS * 0.5, EPS * 0.5);
        assert!(are_points_close(&[p0, p1, p2]));
        // Not close
        let p0 = Point::new(0., 0., 0.);
        let p1 = Point::new(1., 0., 0.);
        let p2 = Point::new(2., 0., 0.);
        assert!(!are_points_close(&[p0, p1, p2]));
    }

    #[test]
    fn test_are_points_collinear() {
        // Collinear
        let p0 = Point::new(0., 0., 0.);
        let p1 = Point::new(1., 0., 0.);
        let p2 = Point::new(2., 0., 0.);
        let p3 = Point::new(3., 0., 0.);
        let p4 = Point::new(4., 0., 0.);
        let p5 = Point::new(5., 0., 0.);
        assert!(are_points_collinear(&[p0]));
        assert!(are_points_collinear(&[p0, p1]));
        assert!(are_points_collinear(&[p0, p1, p2]));
        assert!(are_points_collinear(&[p0, p1, p2, p3]));
        assert!(are_points_collinear(&[p0, p1, p2, p3, p4]));
        assert!(are_points_collinear(&[p0, p1, p2, p3, p4, p5]));
        assert!(are_points_collinear(&[p5, p1, p2, p4, p3, p0]));  // Shuffled order
        assert!(first_3_noncollinear(&[p0, p1, p2, p3, p4, p5]).is_none());
        // Not collinear
        let p3 = Point::new(1., 1., 0.);
        assert!(!are_points_collinear(&[p0, p1, p2, p3]));
    }

    #[test]
    fn test_are_points_coplanar() {
        // Coplanar (z = 0)
        let p0 = Point::new(0., 0., 0.);
        let p1 = Point::new(1., 0., 0.);
        let p2 = Point::new(0., 1., 0.);
        let p3 = Point::new(2., 3., 0.);
        assert!(are_points_coplanar(&[p0, p1]));
        assert!(are_points_coplanar(&[p0, p1, p2]));
        assert!(are_points_coplanar(&[p0, p1, p2, p3]));

        // Non-coplanar: p4 is off the plane z = 0
        let p4 = Point::new(3., 1., 1.);
        assert!(!are_points_coplanar(&[p0, p1, p2, p3, p4]));
    }
}
