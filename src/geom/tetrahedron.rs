use crate::Point;
use crate::Vector;

/// Calculates volume using the Cayley-Menger formula.
pub fn tetrahedron_volume(pt0: Point, pt1: Point, pt2: Point, pt3: Point) -> f64 {
    // Edge lengths
    let a = (pt1 - pt0).length();
    let b = (pt2 - pt0).length();
    let c = (pt3 - pt0).length();
    let d = (pt2 - pt1).length();
    let e = (pt3 - pt2).length();
    let f = (pt3 - pt1).length();

    let a2 = a.powi(2);
    let b2 = b.powi(2);
    let c2 = c.powi(2);
    let d2 = d.powi(2);
    let e2 = e.powi(2);
    let f2 = f.powi(2);

    let x = b2 + c2 - e2;
    let y = a2 + c2 - f2;
    let z = a2 + b2 - d2;

    let x2 = x.powi(2);
    let y2 = y.powi(2);
    let z2 = z.powi(2);

    let nominator = 4. * a2 * b2 * c2 - a2 * x2 - b2 * y2 - c2 * z2 + x * y * z;

    // Use max(0.0) to handle negative values from floating-point errors
    // for degenerate tetrahedrons (e.g., all 4 points coplanar)
    nominator.max(0.0).sqrt() / 12.
}

/// Returns tetrahedron centroid (i.e. average of each vertices)
pub fn tetrahedron_centroid(pt0: Point, pt1: Point, pt2: Point, pt3: Point) -> Point {
    let x = (pt0.x + pt1.x + pt2.x + pt3.x) / 4.;
    let y = (pt0.y + pt1.y + pt2.y + pt3.y) / 4.;
    let z = (pt0.z + pt1.z + pt2.z + pt3.z) / 4.;
    Point::new(x, y, z)
}

/// Circumsphere center and radius-squared of a tetrahedron.
///
/// Returns `None` for degenerate (coplanar/coincident) vertices.
pub fn circumsphere(p0: Point, p1: Point, p2: Point, p3: Point) -> Option<(Point, f64)> {
    let a = Vector::from_points(p0, p1);
    let b = Vector::from_points(p0, p2);
    let c = Vector::from_points(p0, p3);

    let a_sq = a.dot(&a);
    let b_sq = b.dot(&b);
    let c_sq = c.dot(&c);

    let ab = a.dot(&b);
    let ac = a.dot(&c);
    let bc = b.dot(&c);

    // Solve via Cramer's rule for the 3x3 system:
    // 2 * [[a·a, a·b, a·c], [a·b, b·b, b·c], [a·c, b·c, c·c]] * [α, β, γ] = [|a|², |b|², |c|²]
    let det =
        a_sq * (b_sq * c_sq - bc * bc) - ab * (ab * c_sq - bc * ac) + ac * (ab * bc - b_sq * ac);

    // `det` scales with length^6, so use a scale-aware tolerance to avoid
    // classifying valid small tetrahedra as degenerate.
    let scale = a_sq
        .max(b_sq)
        .max(c_sq)
        .max(ab.abs())
        .max(ac.abs())
        .max(bc.abs());
    if scale <= f64::EPSILON {
        return None;
    }
    let det_tol = 1e-12 * scale.powi(3);
    if det.abs() <= det_tol {
        return None;
    }

    let inv_det = 0.5 / det;

    let alpha = inv_det
        * (a_sq * (b_sq * c_sq - bc * bc)
            + b_sq * (ac * bc - ab * c_sq)
            + c_sq * (ab * bc - ac * b_sq));
    let beta = inv_det
        * (a_sq * (bc * ac - ab * c_sq)
            + b_sq * (a_sq * c_sq - ac * ac)
            + c_sq * (ab * ac - a_sq * bc));
    let gamma = inv_det
        * (a_sq * (ab * bc - b_sq * ac)
            + b_sq * (ab * ac - a_sq * bc)
            + c_sq * (a_sq * b_sq - ab * ab));

    let center = Point::new(
        p0.x + alpha * a.dx + beta * b.dx + gamma * c.dx,
        p0.y + alpha * a.dy + beta * b.dy + gamma * c.dy,
        p0.z + alpha * a.dz + beta * b.dz + gamma * c.dz,
    );

    let dx = center.x - p0.x;
    let dy = center.y - p0.y;
    let dz = center.z - p0.z;
    let radius_sq = dx * dx + dy * dy + dz * dz;

    Some((center, radius_sq))
}

/// True if `test` lies strictly inside the circumsphere of (p0, p1, p2, p3).
pub fn is_point_in_circumsphere(p0: Point, p1: Point, p2: Point, p3: Point, test: Point) -> bool {
    if let Some((center, radius_sq)) = circumsphere(p0, p1, p2, p3) {
        let dx = center.x - test.x;
        let dy = center.y - test.y;
        let dz = center.z - test.z;
        let dist_sq = dx * dx + dy * dy + dz * dz;
        dist_sq < radius_sq - 1e-10
    } else {
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_circumsphere_right_corner() {
        let p0 = Point::new(0.0, 0.0, 0.0);
        let p1 = Point::new(1.0, 0.0, 0.0);
        let p2 = Point::new(0.0, 1.0, 0.0);
        let p3 = Point::new(0.0, 0.0, 1.0);

        let (center, r_sq) = circumsphere(p0, p1, p2, p3).unwrap();
        assert!((center.x - 0.5).abs() < 1e-10);
        assert!((center.y - 0.5).abs() < 1e-10);
        assert!((center.z - 0.5).abs() < 1e-10);
        // radius = sqrt(3)/2, radius² = 3/4
        assert!((r_sq - 0.75).abs() < 1e-10);
    }

    #[test]
    fn test_circumsphere_degenerate() {
        // Coplanar points -> None
        let p0 = Point::new(0.0, 0.0, 0.0);
        let p1 = Point::new(1.0, 0.0, 0.0);
        let p2 = Point::new(0.0, 1.0, 0.0);
        let p3 = Point::new(1.0, 1.0, 0.0);

        assert!(circumsphere(p0, p1, p2, p3).is_none());
    }

    #[test]
    fn test_circumsphere_small_scale_non_degenerate() {
        let s = 0.01;
        let p0 = Point::new(0.0, 0.0, 0.0);
        let p1 = Point::new(s, 0.0, 0.0);
        let p2 = Point::new(0.0, s, 0.0);
        let p3 = Point::new(0.0, 0.0, s);

        let (center, r_sq) = circumsphere(p0, p1, p2, p3).expect("small tetra should be valid");
        assert!((center.x - 0.5 * s).abs() < 1e-12);
        assert!((center.y - 0.5 * s).abs() < 1e-12);
        assert!((center.z - 0.5 * s).abs() < 1e-12);
        assert!((r_sq - 0.75 * s * s).abs() < 1e-12);
    }

    #[test]
    fn test_in_circumsphere() {
        let p0 = Point::new(0.0, 0.0, 0.0);
        let p1 = Point::new(1.0, 0.0, 0.0);
        let p2 = Point::new(0.0, 1.0, 0.0);
        let p3 = Point::new(0.0, 0.0, 1.0);

        // Center is inside
        let center = Point::new(0.5, 0.5, 0.5);
        assert!(is_point_in_circumsphere(p0, p1, p2, p3, center));

        // Point far away is outside
        let far = Point::new(10.0, 10.0, 10.0);
        assert!(!is_point_in_circumsphere(p0, p1, p2, p3, far));

        // Vertex is on the sphere (not strictly inside)
        assert!(!is_point_in_circumsphere(p0, p1, p2, p3, p0));
    }
}
