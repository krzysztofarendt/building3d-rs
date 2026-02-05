use crate::Point;

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
