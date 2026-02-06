use crate::Point;
use crate::geom::EPS;
use crate::geom::bboxes::is_point_inside_bbox;
use crate::geom::point::check::are_points_collinear;
use crate::geom::point::check::is_point_on_same_side;
use crate::geom::projection::PlaneBasis;
use crate::geom::vector::Vector;
use crate::vecutils::{max, min};
use anyhow::{Result, anyhow};
use serde::{Deserialize, Serialize};

/// Type for holding vertex indices for a triangle.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
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
    let vn = vn
        .normalize()
        .map_err(|_| anyhow!("Normal vector cannot have zero length"))?;

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
            return triangulate(pts, vn, num_try + 1);
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
                    // Also skip points that coincide with any triangle vertex
                    // (can happen with bridge edges in hole-merged polygons)
                    let test_pt = pts[*test_id];
                    if test_pt.is_close(&pts[prev_id])
                        || test_pt.is_close(&pts[curr_id])
                        || test_pt.is_close(&pts[next_id])
                    {
                        continue;
                    }
                    let point_inside =
                        is_point_inside_triangle(test_pt, pts[prev_id], pts[curr_id], pts[next_id]);
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
    let vn = match vn.normalize() {
        Ok(v) => v,
        Err(_) => return false,
    };

    let v1: Vector = *p2 - *p1;
    let v2: Vector = *p3 - *p2;
    let v1v2_n = v1.cross(&v2).normalize();
    if v1v2_n.is_err() {
        return false; // Collinear points p1, p2, p3
    }
    let v1v2_n = v1v2_n.unwrap();

    // If v1v2_n is equal to the vn, the corner must be convex
    v1v2_n.is_close(&vn)
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

/// Triangulates a polygon with holes by merging holes into the outer boundary
/// using bridge edges, then delegating to the standard ear-clipping `triangulate()`.
///
/// The outer boundary should be wound counter-clockwise w.r.t. the normal.
/// Holes should be wound clockwise (opposite to outer).
pub fn triangulate_with_holes(
    outer: Vec<Point>,
    holes: Vec<Vec<Point>>,
    vn: Vector,
    num_try: usize,
) -> Result<(Vec<Point>, Vec<TriangleIndex>)> {
    if holes.is_empty() {
        return triangulate(outer, vn, num_try);
    }

    let basis = PlaneBasis::from_normal(outer[0], vn)
        .ok_or_else(|| anyhow!("Cannot create plane basis for hole merging"))?;

    // Project outer and holes to 2D
    let outer_2d: Vec<(f64, f64)> = outer.iter().map(|p| basis.project(*p)).collect();
    let holes_2d: Vec<Vec<(f64, f64)>> = holes
        .iter()
        .map(|h| h.iter().map(|p| basis.project(*p)).collect())
        .collect();

    // Merge holes into outer boundary in 2D
    let merged_2d = merge_holes_into_boundary(outer_2d, holes_2d);

    // Unproject back to 3D
    let merged_3d: Vec<Point> = merged_2d
        .iter()
        .map(|&(u, v)| basis.unproject(u, v))
        .collect();

    triangulate(merged_3d, vn, num_try)
}

/// 2D cross product (z-component of the 3D cross product of 2D vectors).
fn cross_2d(ox: f64, oy: f64, ax: f64, ay: f64, bx: f64, by: f64) -> f64 {
    (ax - ox) * (by - oy) - (ay - oy) * (bx - ox)
}

/// Finds the intersection x-coordinate of a horizontal ray from `(mx, my)`
/// going rightward (+x) with the segment `(ax, ay)-(bx, by)`.
/// Returns `Some(x)` if the ray hits the segment, `None` otherwise.
fn ray_intersect_segment_x(mx: f64, my: f64, ax: f64, ay: f64, bx: f64, by: f64) -> Option<f64> {
    // The segment must straddle the ray's y coordinate
    let (min_y, max_y) = if ay < by { (ay, by) } else { (by, ay) };
    if my < min_y - EPS || my > max_y + EPS {
        return None;
    }
    // Avoid division by zero for horizontal segments
    let dy = by - ay;
    if dy.abs() < EPS {
        return None;
    }
    let t = (my - ay) / dy;
    if !(-EPS..=1.0 + EPS).contains(&t) {
        return None;
    }
    let ix = ax + t * (bx - ax);
    // Must be to the right of (or at) the ray origin
    if ix >= mx - EPS { Some(ix) } else { None }
}

/// Checks if a 2D point `(px, py)` is inside triangle `(ax, ay), (bx, by), (cx, cy)`
/// using the sign-of-cross-product method.
#[allow(clippy::too_many_arguments)]
fn is_point_in_triangle_2d(
    px: f64,
    py: f64,
    ax: f64,
    ay: f64,
    bx: f64,
    by: f64,
    cx: f64,
    cy: f64,
) -> bool {
    let d1 = cross_2d(px, py, ax, ay, bx, by);
    let d2 = cross_2d(px, py, bx, by, cx, cy);
    let d3 = cross_2d(px, py, cx, cy, ax, ay);

    let has_neg = (d1 < 0.0) || (d2 < 0.0) || (d3 < 0.0);
    let has_pos = (d1 > 0.0) || (d2 > 0.0) || (d3 > 0.0);

    !(has_neg && has_pos)
}

/// Merges all holes into the outer boundary using the bridge-vertex algorithm.
///
/// 1. Sort holes by their rightmost vertex (max x), descending.
/// 2. For each hole, find the rightmost vertex M, cast a ray rightward,
///    find the nearest visible vertex on the outer boundary, and splice
///    the hole in via bridge edges.
fn merge_holes_into_boundary(
    mut outer: Vec<(f64, f64)>,
    holes: Vec<Vec<(f64, f64)>>,
) -> Vec<(f64, f64)> {
    if holes.is_empty() {
        return outer;
    }

    // For each hole, find index of rightmost vertex and the x value
    let mut hole_info: Vec<(usize, usize, f64)> = holes
        .iter()
        .enumerate()
        .map(|(hi, h)| {
            let (mi, mx) = h
                .iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.0.partial_cmp(&b.0).unwrap())
                .unwrap();
            (hi, mi, mx.0)
        })
        .collect();

    // Sort by rightmost x, descending (process rightmost holes first)
    hole_info.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap());

    for (hole_idx, m_idx, _) in hole_info {
        let hole = &holes[hole_idx];
        let m = hole[m_idx];

        // Cast ray from M rightward; find nearest edge intersection on outer boundary
        let mut best_ix: Option<f64> = None;
        let mut best_edge_idx: Option<usize> = None;

        let n = outer.len();
        for i in 0..n {
            let a = outer[i];
            let b = outer[(i + 1) % n];
            if let Some(ix) = ray_intersect_segment_x(m.0, m.1, a.0, a.1, b.0, b.1)
                && (best_ix.is_none() || ix < best_ix.unwrap())
            {
                best_ix = Some(ix);
                best_edge_idx = Some(i);
            }
        }

        let Some(ix) = best_ix else {
            continue; // hole is outside outer boundary (shouldn't happen for valid input)
        };
        let edge_idx = best_edge_idx.unwrap();

        // The intersection point I is at (ix, m.1)
        // Candidate bridge vertex P is the endpoint of the intersected edge
        // with the larger x. If both have the same x, pick the one closer to I.
        let a = outer[edge_idx];
        let b = outer[(edge_idx + 1) % n];
        let p_idx = if (a.0 - b.0).abs() < EPS {
            // Same x — pick the one closer to I in y
            if (a.1 - m.1).abs() < (b.1 - m.1).abs() {
                edge_idx
            } else {
                (edge_idx + 1) % n
            }
        } else if a.0 > b.0 {
            edge_idx
        } else {
            (edge_idx + 1) % n
        };
        let mut bridge_idx = p_idx;
        let mut bridge_pt = outer[bridge_idx];

        // Check if any outer vertex inside the triangle M-I-P is closer
        // (i.e. would block the bridge). If so, pick the one with smallest
        // angle to the horizontal ray as the bridge vertex.
        let i_pt = (ix, m.1);
        let n = outer.len();
        for (idx, &v) in outer.iter().enumerate().take(n) {
            if is_point_in_triangle_2d(v.0, v.1, m.0, m.1, i_pt.0, i_pt.1, bridge_pt.0, bridge_pt.1)
            {
                // Pick the vertex with the smallest angle from M to horizontal
                let angle_v = (v.1 - m.1).abs().atan2(v.0 - m.0);
                let angle_best = (bridge_pt.1 - m.1).abs().atan2(bridge_pt.0 - m.0);
                if angle_v < angle_best {
                    bridge_idx = idx;
                    bridge_pt = v;
                }
            }
        }

        // Splice hole into outer boundary at bridge_idx
        // New boundary: outer[0..=bridge_idx] + hole[m_idx..] + hole[..=m_idx] + outer[bridge_idx..]
        let hole_len = hole.len();
        let mut new_outer = Vec::with_capacity(outer.len() + hole_len + 2);

        // Part of outer up to and including bridge vertex
        new_outer.extend(&outer[..=bridge_idx]);

        // Hole starting from M, wrapping around
        for j in 0..=hole_len {
            new_outer.push(hole[(m_idx + j) % hole_len]);
        }

        // Bridge back: duplicate bridge vertex
        new_outer.push(outer[bridge_idx]);

        // Rest of outer after bridge vertex
        if bridge_idx + 1 < outer.len() {
            new_outer.extend(&outer[bridge_idx + 1..]);
        }

        outer = new_outer;
    }

    outer
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

    #[test]
    fn test_triangulate_square_with_square_hole() -> Result<()> {
        // 4x4 outer square (CCW)
        let outer = vec![
            Point::new(0., 0., 0.),
            Point::new(4., 0., 0.),
            Point::new(4., 4., 0.),
            Point::new(0., 4., 0.),
        ];
        // 1x1 hole centered at (2, 2) (CW winding)
        let hole = vec![
            Point::new(1.5, 1.5, 0.),
            Point::new(1.5, 2.5, 0.),
            Point::new(2.5, 2.5, 0.),
            Point::new(2.5, 1.5, 0.),
        ];
        let vn = Vector::new(0., 0., 1.);
        let (pts, tri) = triangulate_with_holes(outer, vec![hole], vn, 0)?;

        // After merging, we should have triangles covering the area outside the hole
        assert!(!tri.is_empty());

        // No triangle should have its centroid inside the hole
        for t in &tri {
            let cx = (pts[t.0].x + pts[t.1].x + pts[t.2].x) / 3.0;
            let cy = (pts[t.0].y + pts[t.1].y + pts[t.2].y) / 3.0;
            let in_hole = cx > 1.5 + EPS && cx < 2.5 - EPS && cy > 1.5 + EPS && cy < 2.5 - EPS;
            assert!(
                !in_hole,
                "Triangle centroid ({}, {}) is inside the hole",
                cx, cy
            );
        }

        Ok(())
    }

    #[test]
    fn test_triangulate_with_multiple_holes() -> Result<()> {
        // 10x10 outer square
        let outer = vec![
            Point::new(0., 0., 0.),
            Point::new(10., 0., 0.),
            Point::new(10., 10., 0.),
            Point::new(0., 10., 0.),
        ];
        // Hole 1 at (1,1)-(3,3) (CW)
        let hole1 = vec![
            Point::new(1., 1., 0.),
            Point::new(1., 3., 0.),
            Point::new(3., 3., 0.),
            Point::new(3., 1., 0.),
        ];
        // Hole 2 at (5,5)-(7,7) (CW)
        let hole2 = vec![
            Point::new(5., 5., 0.),
            Point::new(5., 7., 0.),
            Point::new(7., 7., 0.),
            Point::new(7., 5., 0.),
        ];
        let vn = Vector::new(0., 0., 1.);
        let (_pts, tri) = triangulate_with_holes(outer, vec![hole1, hole2], vn, 0)?;

        assert!(!tri.is_empty());
        Ok(())
    }

    #[test]
    fn test_merge_holes_into_boundary_vertex_count() {
        // 4x4 outer square
        let outer = vec![(0., 0.), (4., 0.), (4., 4.), (0., 4.)];
        // 2x2 hole
        let hole = vec![(1., 1.), (1., 3.), (3., 3.), (3., 1.)];
        let merged = merge_holes_into_boundary(outer.clone(), vec![hole.clone()]);

        // Merged boundary: outer(4) + hole(4) + 1 (return to hole start) + 1 (bridge back) = 10
        // = outer vertices + hole vertices + 2 bridge duplicates
        assert_eq!(merged.len(), outer.len() + hole.len() + 2);
    }
}
