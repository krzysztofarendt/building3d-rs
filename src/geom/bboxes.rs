use crate::geom::EPS;
use crate::geom::point::Point;

/// Checks whether a point is inside the bounding box holding all points `pts`.
pub fn is_point_inside_bbox(ptest: Point, pts: &[Point]) -> bool {
    if pts.is_empty() {
        return false;
    }
    let (pmin, pmax) = bounding_box(pts);
    // A point is outside if ANY coordinate is out of range.
    !(ptest.x < pmin.x
        || ptest.x > pmax.x
        || ptest.y < pmin.y
        || ptest.y > pmax.y
        || ptest.z < pmin.z
        || ptest.z > pmax.z)
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
    if pts.is_empty() {
        let zero = Point::new(0.0, 0.0, 0.0);
        return (zero, zero);
    }

    let mut xmin = pts[0].x;
    let mut xmax = pts[0].x;
    let mut ymin = pts[0].y;
    let mut ymax = pts[0].y;
    let mut zmin = pts[0].z;
    let mut zmax = pts[0].z;

    for p in pts.iter().skip(1) {
        if xmin.is_nan() || p.x < xmin {
            xmin = p.x;
        }
        if xmax.is_nan() || p.x > xmax {
            xmax = p.x;
        }
        if ymin.is_nan() || p.y < ymin {
            ymin = p.y;
        }
        if ymax.is_nan() || p.y > ymax {
            ymax = p.y;
        }
        if zmin.is_nan() || p.z < zmin {
            zmin = p.z;
        }
        if zmax.is_nan() || p.z > zmax {
            zmax = p.z;
        }
    }

    (Point::new(xmin, ymin, zmin), Point::new(xmax, ymax, zmax))
}
