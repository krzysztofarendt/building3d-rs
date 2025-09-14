use crate::Point;
use crate::geom::EPS;
use crate::geom::IsClose;
use crate::geom::bboxes::is_point_inside_bbox;
use crate::geom::point::check::are_points_collinear;
use crate::geom::point::check::is_point_on_same_side;
use crate::geom::vector::Vector;
use crate::vecutils::{max, min};
use anyhow::{Result, anyhow};

/// Type for holding vertex indices for a triangle.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TriangleIndex(pub usize, pub usize, pub usize);

/// Triangulates the polygon defined by points `pts` and normal `vn`.
pub fn triangulate(
    mut pts: Vec<Point>,
    vn: Vector,
    num_try: usize,
) -> Result<(Vec<Point>, Vec<TriangleIndex>)> {
    if num_try >= 2 {
        return Err(anyhow!("Ear-clipping algorithm failed."));
    }
    if vn.length().is_close(0.) {
        return Err(anyhow!("Normal vector cannot have zero length"));
    }

    let mut vertices: Vec<(usize, &Point)> = Vec::new();
    for (i, p) in pts.iter().enumerate() {
        vertices.push((i, p));
    }
    let mut triangles: Vec<TriangleIndex> = Vec::new();
    let mut pos: usize = 0;
    let mut num_fail: usize = 0;

    while vertices.len() > 2 {
        if num_fail > pts.len() {
            // Try with flipped points
            pts = pts.iter().rev().cloned().collect();
            return triangulate(pts, vn, num_try + 1); // TODO: shouldn't num_try be resetted here?
        }

        // If last vertex, start from the beginning
        if pos > vertices.len() - 1 {
            pos = 0;
        }

        let prev_pos = if pos > 0 { pos - 1 } else { vertices.len() - 1 };
        let next_pos = if pos < vertices.len() - 1 { pos + 1 } else { 0 };

        let (prev_id, prev_pt) = vertices[prev_pos];
        let (curr_id, curr_pt) = vertices[pos];
        let (next_id, next_pt) = vertices[next_pos];

        let convex_corner = is_corner_convex(prev_pt, curr_pt, next_pt, &vn);

        if convex_corner {
            // Check if no other point is within this triangle
            // Needed for non-convex polygons
            let mut any_point_inside = false;
            for (test_id, _) in vertices.iter() {
                if ![prev_id, curr_id, next_id].contains(test_id) {
                    let point_inside = is_point_inside_triangle(
                        pts[*test_id],
                        pts[prev_id],
                        pts[curr_id],
                        pts[next_id],
                    );
                    if point_inside {
                        any_point_inside = true;
                        break;
                    }
                }
            }
            if !any_point_inside {
                // Add triangle
                triangles.push(TriangleIndex(prev_id, curr_id, next_id));
                // Remove pos from index
                vertices.remove(pos); // TODO: Consider linked list
                continue;
            } else {
                // There is some point inside this triangle
                // so it is not an ear
                num_fail += 1;
            }
        } else {
            // Non-convex corner
            num_fail += 1;
        }
        pos += 1;
    }

    Ok((pts, triangles))
}

/// Checks if the angle between p2->p1 and p2->p3 is less than 180 degrees
///
/// It is done by comparing the polygon normal vector with the cross
/// product p2->p3 x p2->p1. The points p1, p2, p3 should be ordered
/// counter-clockwise with respect to the surface front side.
///
/// # Panics
/// It panics if the length of the normal vector vn is not 1.
pub fn is_corner_convex(p1: &Point, p2: &Point, p3: &Point, vn: &Vector) -> bool {
    assert!((vn.length() - 1.0).abs() < EPS);

    let v1: Vector = *p2 - *p1;
    let v2: Vector = *p3 - *p2;
    let v1v2_n = v1.cross(&v2).normalize();
    if v1v2_n.is_err() {
        return false; // Collinear points p1, p2, p3
    }
    let v1v2_n = v1v2_n.unwrap();

    // If v1v2_n is equal to the vn, the corner must be convex
    v1v2_n.is_close(vn)
}

/// Tests if point `ptest` is inside the triangle `(p1, p2, p3)`.
///
/// Using the "same side technique" described at:
/// https://blackpawn.com/texts/pointinpoly/
/// This function does not test if the point is coplanar with the triangle.
pub fn is_point_inside_triangle(ptest: Point, p1: Point, p2: Point, p3: Point) -> bool {
    if !is_point_inside_bbox(ptest, &[p1, p2, p3]) {
        return false;
    }
    // Test if the point is at any of the three vertices
    if ptest.is_close(&p1) || ptest.is_close(&p2) || ptest.is_close(&p3) {
        return true;
    }
    // Test if it's at any of the edges
    for (pa, pb) in [(p1, p2), (p2, p3), (p3, p1)].iter() {
        let pts = vec![*pa, *pb, ptest];
        if are_points_collinear(&pts) {
            // ptest is collinear, but is it on the edge or outside the triangle?
            return !(ptest.x > max(&[pa.x, pb.x]) + EPS
                || ptest.y > max(&[pa.y, pb.y]) + EPS
                || ptest.z > max(&[pa.z, pb.z]) + EPS
                || ptest.x < min(&[pa.x, pb.x]) - EPS
                || ptest.y < min(&[pa.y, pb.y]) - EPS
                || ptest.z < min(&[pa.z, pb.z]) - EPS);
        }
    }

    // Test if ptest is inside
    let side1 = is_point_on_same_side(p1, p2, ptest, p3).unwrap_or(false);
    let side2 = is_point_on_same_side(p2, p3, ptest, p1).unwrap_or(false);
    let side3 = is_point_on_same_side(p3, p1, ptest, p2).unwrap_or(false);

    side1 && side2 && side3
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::vecutils::roll;

    #[test]
    fn test_triangulate_square() -> Result<()> {
        let p0 = Point::new(0., 0., 0.);
        let p1 = Point::new(1., 0., 0.);
        let p2 = Point::new(1., 1., 0.);
        let p3 = Point::new(0., 1., 0.);
        let pts = vec![p0, p1, p2, p3];
        let vn = Vector::new(0., 0., 1.);
        let (pts, tri) = triangulate(pts, vn, 0)?;
        assert!(tri.len() == 2);
        assert!(pts.len() == 4);

        Ok(())
    }

    #[test]
    fn test_triangulate_l_shape() -> Result<()> {
        let p0 = Point::new(0., 0., 0.);
        let p1 = Point::new(1., 0., 0.);
        let p2 = Point::new(1., 1., 0.);
        let p3 = Point::new(2., 1., 0.);
        let p4 = Point::new(2., 2., 0.);
        let p5 = Point::new(0., 2., 0.);
        let mut pts = vec![p0, p1, p2, p3, p4, p5];
        let vn = Vector::new(0., 0., 1.);
        let mut tri: Vec<TriangleIndex>;

        // Test at different starting points
        let expected_num_triangles = 4;
        for i in 0..pts.len() {
            if i > 0 {
                roll(&mut pts, 1);
            }
            (pts, tri) = triangulate(pts, vn, 0)?;
            assert!(tri.len() == expected_num_triangles);
            for ix in tri.iter() {
                let tri_vn = Vector::normal(pts[ix.0], pts[ix.1], pts[ix.2])?;
                assert!(tri_vn.is_close(&vn));
            }
        }
        Ok(())
    }

    #[test]
    fn test_triangulate_u_shape() -> Result<()> {
        let p0 = Point::new(0., 0., 0.);
        let p1 = Point::new(1., 0., 0.);
        let p2 = Point::new(1., 1., 0.);
        let p3 = Point::new(2., 1., 0.);
        let p4 = Point::new(2., 0., 0.);
        let p5 = Point::new(3., 0., 0.);
        let p6 = Point::new(3., 2., 0.);
        let p7 = Point::new(0., 2., 0.);
        let mut pts = vec![p0, p1, p2, p3, p4, p5, p6, p7];
        let vn = Vector::new(0., 0., 1.);
        let mut tri: Vec<TriangleIndex>;

        // Test at different starting points
        let expected_num_triangles = 6;
        for i in 0..pts.len() {
            if i > 0 {
                roll(&mut pts, 1);
            }
            (pts, tri) = triangulate(pts, vn, 0)?;
            assert!(tri.len() == expected_num_triangles);
            for ix in tri.iter() {
                let tri_vn = Vector::normal(pts[ix.0], pts[ix.1], pts[ix.2])?;
                assert!(tri_vn.is_close(&vn));
            }
        }
        Ok(())
    }

    #[test]
    fn test_triangulate_c_shape() -> Result<()> {
        let p0 = Point::new(0.75, 0.75, 1.0);
        let p1 = Point::new(0.75, 0.25, 1.0);
        let p2 = Point::new(0.25, 0.25, 1.0);
        let p3 = Point::new(0.0, 0.0, 1.0);
        let p4 = Point::new(1.0, 0.0, 1.0);
        let p5 = Point::new(1.0, 1.0, 1.0);
        let p6 = Point::new(0.0, 1.0, 1.0);
        let p7 = Point::new(0.25, 0.75, 1.0);
        let mut pts = vec![p0, p1, p2, p3, p4, p5, p6, p7];
        let num_pts = pts.len();
        let tri: Vec<TriangleIndex>;
        let vn = Vector::new(0., 0., 1.);
        (pts, tri) = triangulate(pts, vn, 0)?;
        let expected_num_triangles = 6;
        assert!(tri.len() == expected_num_triangles);
        assert!(pts.len() == num_pts);
        Ok(())
    }

    #[test]
    fn test_triangulate_one_more() -> Result<()> {
        // This example was failing long time ago.
        let p0 = Point::new(0.5, 0.5, 0.0);
        let p1 = Point::new(0.5, 1.0, 0.0);
        let p2 = Point::new(1.0, 1.0, 0.0);
        let p3 = Point::new(1.0, 0.0, 0.0);
        let p4 = Point::new(0.0, 0.0, 0.0);
        let mut pts = vec![p0, p1, p2, p3, p4];
        let vn = Vector::new(0., 0., 1.);
        let mut tri: Vec<TriangleIndex>;
        let expected_num_triangles = 3;
        for _ in 0..pts.len() {
            roll(&mut pts, 1);
            (pts, tri) = triangulate(pts, vn, 0)?;
            assert!(tri.len() == expected_num_triangles);
        }
        Ok(())
    }

    #[test]
    fn test_is_point_inside_triangle() {
        let p1 = Point::new(1., 0., 0.);
        let p2 = Point::new(0., 0., 0.);
        let p3 = Point::new(0., 1., 0.);

        let ptest = Point::new(0.1, 0.1, 0.0); // inside
        assert!(is_point_inside_triangle(ptest, p1, p2, p3));

        let ptest = Point::new(0.0, 0.0, 0.0); // inside (at the corner)
        assert!(is_point_inside_triangle(ptest, p1, p2, p3));

        let ptest = Point::new(1.0, 0.0, 0.0); // inside (at the corner)
        assert!(is_point_inside_triangle(ptest, p1, p2, p3));

        let ptest = Point::new(0.5, 0.5, 0.0); // inside (at the edge)
        assert!(is_point_inside_triangle(ptest, p1, p2, p3));

        let ptest = Point::new(0.51, 0.51, 0.0); // outside
        assert!(!is_point_inside_triangle(ptest, p1, p2, p3));
    }
}
