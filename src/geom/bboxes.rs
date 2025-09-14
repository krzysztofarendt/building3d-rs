use crate::geom::point::Point;

/// Checks whether a point is inside the bounding box holding all points `pts`.
pub fn is_point_inside_bbox(ptest: Point, pts: &[Point]) -> bool {
    let (pmin, pmax) = bounding_box(pts);
    if (ptest.x < pmin.x && ptest.y < pmin.y && ptest.z < pmin.z)
        || (ptest.x > pmax.x && ptest.y > pmax.y && ptest.z > pmax.z)
    {
        // Point outside bbox
        false
    } else {
        // Point inside bbox
        true
    }
}

pub fn bounding_box(pts: &[Point]) -> (Point, Point) {
    let x: Vec<f64> = pts.iter().map(|v| v.x).collect();
    let y: Vec<f64> = pts.iter().map(|v| v.y).collect();
    let z: Vec<f64> = pts.iter().map(|v| v.z).collect();
    let xmin = *x.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
    let xmax = *x.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
    let ymin = *y.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
    let ymax = *y.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
    let zmin = *z.iter().min_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();
    let zmax = *z.iter().max_by(|a, b| a.partial_cmp(b).unwrap()).unwrap();

    let pmin = Point::new(xmin, ymin, zmin);
    let pmax = Point::new(xmax, ymax, zmax);

    (pmin, pmax)
}
