use crate::Point;
use crate::Vector;
use crate::geom::EPS;
use crate::geom::triangles::{TriangleIndex, is_point_inside_triangle};

/// Checks if a point lies inside a polygon.
///
/// The polygon is defined by its vertices `pts`, triangulation `tri`, and normal vector `vn`.
/// If `boundary_in` is true, points on the boundary (edges or vertices) are considered inside.
/// `outer_boundary` and `holes` are used for boundary checks instead of `pts` when non-empty,
/// to avoid treating bridge edges (from hole merging) as real boundaries.
///
/// This function first projects the test point onto the polygon's plane, then checks
/// if the projection lies within any of the polygon's triangles.
pub fn is_point_inside_polygon(
    ptest: Point,
    pts: &[Point],
    tri: &[TriangleIndex],
    vn: &Vector,
    boundary_in: bool,
    outer_boundary: &[Point],
    holes: &[Vec<Point>],
) -> bool {
    // First check if the point lies on the polygon's plane
    if !is_point_on_plane(ptest, pts, vn) {
        return false;
    }

    // Check boundary against the real boundaries (outer + holes), not the merged mesh
    let boundary_pts = if !outer_boundary.is_empty() {
        outer_boundary
    } else {
        pts
    };

    // Check if point is on a hole boundary
    for hole in holes {
        if is_point_on_boundary(ptest, hole) {
            return boundary_in;
        }
    }

    // Check if point is on the outer boundary
    if is_point_on_boundary(ptest, boundary_pts) {
        return boundary_in;
    }

    // Check if point is inside any of the triangles
    for triangle in tri {
        let p1 = pts[triangle.0];
        let p2 = pts[triangle.1];
        let p3 = pts[triangle.2];

        if is_point_inside_triangle(ptest, p1, p2, p3) {
            return true;
        }
    }

    false
}

/// Checks if a point lies on the plane defined by the polygon.
///
/// Uses the plane equation: ax + by + cz + d = 0
/// where (a, b, c) is the normal vector and d is computed from any point on the plane.
fn is_point_on_plane(ptest: Point, pts: &[Point], vn: &Vector) -> bool {
    if pts.is_empty() {
        return false;
    }

    // Compute d from first point: d = -(ax + by + cz)
    let p0 = pts[0];
    let d = -(vn.dx * p0.x + vn.dy * p0.y + vn.dz * p0.z);

    // Check if test point satisfies the plane equation
    let distance = vn.dx * ptest.x + vn.dy * ptest.y + vn.dz * ptest.z + d;

    distance.abs() < EPS
}

/// Checks if a point lies on the boundary of the polygon (vertices or edges).
fn is_point_on_boundary(ptest: Point, pts: &[Point]) -> bool {
    let n = pts.len();
    if n < 2 {
        return false;
    }

    // Check if point is at any vertex
    for pt in pts {
        if ptest.is_close(pt) {
            return true;
        }
    }

    // Check if point is on any edge
    for i in 0..n {
        let p1 = pts[i];
        let p2 = pts[(i + 1) % n];

        if ptest.is_on_segment(p1, p2) {
            return true;
        }
    }

    false
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::triangles::triangulate;

    fn make_square() -> (Vec<Point>, Vec<TriangleIndex>, Vector) {
        let pts = vec![
            Point::new(0., 0., 0.),
            Point::new(1., 0., 0.),
            Point::new(1., 1., 0.),
            Point::new(0., 1., 0.),
        ];
        let vn = Vector::new(0., 0., 1.);
        let (pts, tri) = triangulate(pts, vn, 0).unwrap();
        (pts, tri, vn)
    }

    fn make_triangle() -> (Vec<Point>, Vec<TriangleIndex>, Vector) {
        let pts = vec![
            Point::new(0., 0., 0.),
            Point::new(1., 0., 0.),
            Point::new(0.5, 1., 0.),
        ];
        let vn = Vector::new(0., 0., 1.);
        let (pts, tri) = triangulate(pts, vn, 0).unwrap();
        (pts, tri, vn)
    }

    #[test]
    fn test_point_inside_square() {
        let (pts, tri, vn) = make_square();
        let no_holes: &[Vec<Point>] = &[];

        // Center point - should be inside
        let ptest = Point::new(0.5, 0.5, 0.0);
        assert!(is_point_inside_polygon(
            ptest, &pts, &tri, &vn, true, &pts, no_holes
        ));
        assert!(is_point_inside_polygon(
            ptest, &pts, &tri, &vn, false, &pts, no_holes
        ));
    }

    #[test]
    fn test_point_outside_square() {
        let (pts, tri, vn) = make_square();
        let no_holes: &[Vec<Point>] = &[];

        // Point outside - should be false
        let ptest = Point::new(1.5, 0.5, 0.0);
        assert!(!is_point_inside_polygon(
            ptest, &pts, &tri, &vn, true, &pts, no_holes
        ));
        assert!(!is_point_inside_polygon(
            ptest, &pts, &tri, &vn, false, &pts, no_holes
        ));

        // Point above plane - should be false
        let ptest = Point::new(0.5, 0.5, 1.0);
        assert!(!is_point_inside_polygon(
            ptest, &pts, &tri, &vn, true, &pts, no_holes
        ));
        assert!(!is_point_inside_polygon(
            ptest, &pts, &tri, &vn, false, &pts, no_holes
        ));
    }

    #[test]
    fn test_point_on_vertex() {
        let (pts, tri, vn) = make_square();
        let no_holes: &[Vec<Point>] = &[];

        // Point at vertex
        let ptest = Point::new(0., 0., 0.);
        assert!(is_point_inside_polygon(
            ptest, &pts, &tri, &vn, true, &pts, no_holes
        ));
        assert!(!is_point_inside_polygon(
            ptest, &pts, &tri, &vn, false, &pts, no_holes
        ));
    }

    #[test]
    fn test_point_on_edge() {
        let (pts, tri, vn) = make_square();
        let no_holes: &[Vec<Point>] = &[];

        // Point on edge (midpoint of bottom edge)
        let ptest = Point::new(0.5, 0., 0.);
        assert!(is_point_inside_polygon(
            ptest, &pts, &tri, &vn, true, &pts, no_holes
        ));
        assert!(!is_point_inside_polygon(
            ptest, &pts, &tri, &vn, false, &pts, no_holes
        ));
    }

    #[test]
    fn test_point_inside_triangle() {
        let (pts, tri, vn) = make_triangle();
        let no_holes: &[Vec<Point>] = &[];

        // Centroid - should be inside
        let ptest = Point::new(0.5, 0.33, 0.0);
        assert!(is_point_inside_polygon(
            ptest, &pts, &tri, &vn, true, &pts, no_holes
        ));
    }

    #[test]
    fn test_point_outside_triangle() {
        let (pts, tri, vn) = make_triangle();
        let no_holes: &[Vec<Point>] = &[];

        // Outside
        let ptest = Point::new(1.0, 1.0, 0.0);
        assert!(!is_point_inside_polygon(
            ptest, &pts, &tri, &vn, true, &pts, no_holes
        ));
    }

    #[test]
    fn test_point_not_on_plane() {
        let (pts, tri, vn) = make_square();
        let no_holes: &[Vec<Point>] = &[];

        // Point not on the plane
        let ptest = Point::new(0.5, 0.5, 0.1);
        assert!(!is_point_inside_polygon(
            ptest, &pts, &tri, &vn, true, &pts, no_holes
        ));
    }

    #[test]
    fn test_l_shaped_polygon() {
        // L-shaped polygon
        let pts = vec![
            Point::new(0., 0., 0.),
            Point::new(1., 0., 0.),
            Point::new(1., 1., 0.),
            Point::new(2., 1., 0.),
            Point::new(2., 2., 0.),
            Point::new(0., 2., 0.),
        ];
        let vn = Vector::new(0., 0., 1.);
        let (pts, tri) = triangulate(pts, vn, 0).unwrap();
        let no_holes: &[Vec<Point>] = &[];

        // Inside the L
        let ptest = Point::new(0.5, 0.5, 0.0);
        assert!(is_point_inside_polygon(
            ptest, &pts, &tri, &vn, true, &pts, no_holes
        ));

        let ptest = Point::new(0.5, 1.5, 0.0);
        assert!(is_point_inside_polygon(
            ptest, &pts, &tri, &vn, true, &pts, no_holes
        ));

        let ptest = Point::new(1.5, 1.5, 0.0);
        assert!(is_point_inside_polygon(
            ptest, &pts, &tri, &vn, true, &pts, no_holes
        ));

        // Outside the L (in the cutout)
        let ptest = Point::new(1.5, 0.5, 0.0);
        assert!(!is_point_inside_polygon(
            ptest, &pts, &tri, &vn, true, &pts, no_holes
        ));

        // Outside the L (completely outside)
        let ptest = Point::new(3.0, 1.0, 0.0);
        assert!(!is_point_inside_polygon(
            ptest, &pts, &tri, &vn, true, &pts, no_holes
        ));
    }
}
