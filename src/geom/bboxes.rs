use crate::geom::EPS;
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

/// Checks whether a point is strictly inside a bounding box (not on boundary).
pub fn is_point_strictly_inside_bbox(ptest: Point, pmin: Point, pmax: Point) -> bool {
    ptest.x > pmin.x + EPS
        && ptest.x < pmax.x - EPS
        && ptest.y > pmin.y + EPS
        && ptest.y < pmax.y - EPS
        && ptest.z > pmin.z + EPS
        && ptest.z < pmax.z - EPS
}

/// Checks whether two bounding boxes overlap.
///
/// Takes min and max corners of each bbox.
/// Returns true if boxes overlap (including touching).
pub fn are_bboxes_overlapping(min1: Point, max1: Point, min2: Point, max2: Point) -> bool {
    // Boxes don't overlap if separated along any axis
    if max1.x < min2.x - EPS || min1.x > max2.x + EPS {
        return false;
    }
    if max1.y < min2.y - EPS || min1.y > max2.y + EPS {
        return false;
    }
    if max1.z < min2.z - EPS || min1.z > max2.z + EPS {
        return false;
    }
    true
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
