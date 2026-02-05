//! Polygon relationship detection.
//!
//! This module provides functions for determining spatial relationships
//! between polygons: facing, coplanar, crossing, and touching.

use crate::geom::EPS;
use crate::geom::polygon::Polygon;
use crate::geom::segment::SegmentIntersection;
use crate::geom::segment::segment_intersection;

/// Checks if two polygons are facing each other.
///
/// Two polygons are facing if:
/// 1. They are coplanar (lie in the same plane)
/// 2. Their normal vectors point in opposite directions
/// 3. They have overlapping area (share some common region)
///
/// This is typically used to identify interfaces between adjacent solids.
pub fn are_polygons_facing(poly1: &Polygon, poly2: &Polygon) -> bool {
    // Check if normals are opposite
    let dot = poly1.vn.dot(&poly2.vn);
    if dot > -1.0 + EPS {
        return false; // Normals not opposite
    }

    // Check if coplanar
    if !are_polygons_coplanar(poly1, poly2) {
        return false;
    }

    // Check if they have overlapping area
    // For facing polygons, at least one vertex of each should be inside the other
    // or their edges should intersect
    has_overlap(poly1, poly2)
}

/// Checks if two polygons lie in the same plane.
///
/// Two polygons are coplanar if all vertices of both polygons satisfy
/// the same plane equation (within tolerance).
pub fn are_polygons_coplanar(poly1: &Polygon, poly2: &Polygon) -> bool {
    let (a, b, c, d) = poly1.plane_coefficients();

    // Check if all vertices of poly2 satisfy poly1's plane equation
    for pt in poly2.vertices() {
        let dist = a * pt.x + b * pt.y + c * pt.z + d;
        if dist.abs() > EPS {
            return false;
        }
    }

    true
}

/// Checks if two polygons cross (intersect) each other.
///
/// Two polygons cross if they are not coplanar and their edges
/// intersect with each other's surface.
pub fn are_polygons_crossing(poly1: &Polygon, poly2: &Polygon) -> bool {
    // If coplanar, they can't cross in the 3D sense
    // (they might overlap, but that's different)
    if are_polygons_coplanar(poly1, poly2) {
        return false;
    }

    // Check if any edge of poly1 crosses poly2's surface
    for (p1, p2) in poly1.edges() {
        if crate::geom::segment::segment_crosses_polygon(p1, p2, poly2).is_some() {
            return true;
        }
    }

    // Check if any edge of poly2 crosses poly1's surface
    for (p1, p2) in poly2.edges() {
        if crate::geom::segment::segment_crosses_polygon(p1, p2, poly1).is_some() {
            return true;
        }
    }

    false
}

/// Checks if two polygons are touching.
///
/// Two polygons are touching if they are coplanar and share one or more
/// edges or edge segments, but do not overlap in their interior regions.
pub fn are_polygons_touching(poly1: &Polygon, poly2: &Polygon) -> bool {
    // Must be coplanar to touch
    if !are_polygons_coplanar(poly1, poly2) {
        return false;
    }

    // Check for interior overlap by testing if centroids are inside
    let centroid1 = poly1.centroid();
    let centroid2 = poly2.centroid();

    if poly2.is_point_inside(centroid1, false) {
        return false;
    }
    if poly1.is_point_inside(centroid2, false) {
        return false;
    }

    // Check that no vertex of one polygon is strictly inside the other
    for pt in poly1.vertices() {
        if poly2.is_point_inside(*pt, false) {
            return false;
        }
    }

    for pt in poly2.vertices() {
        if poly1.is_point_inside(*pt, false) {
            return false;
        }
    }

    // Check for edge crossings that indicate overlap
    // If edges cross at interior points (not at their endpoints), it's overlap
    let edges1 = poly1.edges();
    let edges2 = poly2.edges();

    let mut has_shared_edge = false;

    for (a1, a2) in &edges1 {
        for (b1, b2) in &edges2 {
            let result = segment_intersection(*a1, *a2, *b1, *b2);
            match result {
                SegmentIntersection::Point(pt) => {
                    // If intersection is at endpoints of both segments, it's touching
                    // If intersection is interior to either segment, it might be overlap
                    let at_a_endpoint = pt.is_close(a1) || pt.is_close(a2);
                    let at_b_endpoint = pt.is_close(b1) || pt.is_close(b2);

                    if at_a_endpoint && at_b_endpoint {
                        has_shared_edge = true;
                    } else {
                        // Interior crossing - this is overlap, not touching
                        return false;
                    }
                }
                SegmentIntersection::Collinear(start, end) => {
                    // Collinear overlap - check if it's the full edge or partial
                    // If edges are exactly the same, it's touching
                    // If one edge partially overlaps another, it's overlap
                    let same_edge = (start.is_close(a1) && end.is_close(a2))
                        || (start.is_close(a2) && end.is_close(a1))
                        || (start.is_close(b1) && end.is_close(b2))
                        || (start.is_close(b2) && end.is_close(b1));

                    if same_edge {
                        has_shared_edge = true;
                    } else {
                        // Partial collinear overlap - this is overlap, not touching
                        return false;
                    }
                }
                _ => {}
            }
        }
    }

    has_shared_edge
}

/// Checks if two coplanar polygons have overlapping area.
fn has_overlap(poly1: &Polygon, poly2: &Polygon) -> bool {
    // Check if any vertex of poly1 is inside poly2
    for pt in poly1.vertices() {
        if poly2.is_point_inside(*pt, true) {
            return true;
        }
    }

    // Check if any vertex of poly2 is inside poly1
    for pt in poly2.vertices() {
        if poly1.is_point_inside(*pt, true) {
            return true;
        }
    }

    // Check if any edges intersect (for cases where vertices don't overlap
    // but edges cross)
    let edges1 = poly1.edges();
    let edges2 = poly2.edges();

    for (a1, a2) in &edges1 {
        for (b1, b2) in &edges2 {
            let result = segment_intersection(*a1, *a2, *b1, *b2);
            match result {
                SegmentIntersection::Point(pt) => {
                    // Make sure intersection is not just at endpoints
                    if !pt.is_close(a1) && !pt.is_close(a2) && !pt.is_close(b1) && !pt.is_close(b2)
                    {
                        return true;
                    }
                }
                SegmentIntersection::Collinear(_, _) => {
                    return true;
                }
                _ => {}
            }
        }
    }

    false
}

/// Calculates the approximate overlap area between two coplanar polygons.
///
/// This is a simplified implementation that returns a boolean-like value:
/// - 0.0 if no overlap
/// - positive value if overlap exists
///
/// For exact overlap area calculation, a more sophisticated polygon
/// clipping algorithm would be needed (e.g., Sutherland-Hodgman).
pub fn polygon_overlap_area(poly1: &Polygon, poly2: &Polygon) -> f64 {
    if !are_polygons_coplanar(poly1, poly2) {
        return 0.0;
    }

    if !has_overlap(poly1, poly2) {
        return 0.0;
    }

    // Simplified: return the smaller polygon's area as upper bound
    // A full implementation would compute the actual intersection polygon
    poly1.area().min(poly2.area())
}

/// Checks if two polygons are exactly equal (same vertices, possibly rotated).
pub fn are_polygons_equal(poly1: &Polygon, poly2: &Polygon) -> bool {
    poly1 == poly2
}

/// Calculates the minimum distance between two polygons.
///
/// Returns 0.0 if polygons touch or overlap.
pub fn polygon_distance(poly1: &Polygon, poly2: &Polygon) -> f64 {
    // Quick check: if they overlap or touch, distance is 0
    if are_polygons_coplanar(poly1, poly2) && has_overlap(poly1, poly2) {
        return 0.0;
    }

    let mut min_dist = f64::MAX;

    // Check distance from each vertex of poly1 to poly2's edges
    for v1 in poly1.vertices() {
        for (e1, e2) in poly2.edges() {
            let dist = crate::geom::segment::distance_point_to_segment(*v1, e1, e2);
            min_dist = min_dist.min(dist);
        }
    }

    // Check distance from each vertex of poly2 to poly1's edges
    for v2 in poly2.vertices() {
        for (e1, e2) in poly1.edges() {
            let dist = crate::geom::segment::distance_point_to_segment(*v2, e1, e2);
            min_dist = min_dist.min(dist);
        }
    }

    // Check edge-to-edge distances
    for (a1, a2) in poly1.edges() {
        for (b1, b2) in poly2.edges() {
            let result = segment_intersection(a1, a2, b1, b2);
            match result {
                SegmentIntersection::Point(_) | SegmentIntersection::Collinear(_, _) => {
                    return 0.0;
                }
                SegmentIntersection::Skew(p1, p2) => {
                    let dist = (p1 - p2).length();
                    min_dist = min_dist.min(dist);
                }
                _ => {}
            }
        }
    }

    min_dist
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::IsClose;
    use crate::{Point, Vector};
    use anyhow::Result;

    fn make_xy_square(x: f64, y: f64, size: f64, z: f64, name: &str) -> Result<Polygon> {
        let pts = vec![
            Point::new(x, y, z),
            Point::new(x + size, y, z),
            Point::new(x + size, y + size, z),
            Point::new(x, y + size, z),
        ];
        Polygon::new(name, pts, None)
    }

    #[test]
    fn test_coplanar_same_plane() -> Result<()> {
        let poly1 = make_xy_square(0.0, 0.0, 1.0, 0.0, "p1")?;
        let poly2 = make_xy_square(2.0, 0.0, 1.0, 0.0, "p2")?;

        assert!(are_polygons_coplanar(&poly1, &poly2));
        Ok(())
    }

    #[test]
    fn test_coplanar_different_planes() -> Result<()> {
        let poly1 = make_xy_square(0.0, 0.0, 1.0, 0.0, "p1")?;
        let poly2 = make_xy_square(0.0, 0.0, 1.0, 1.0, "p2")?; // Different z

        assert!(!are_polygons_coplanar(&poly1, &poly2));
        Ok(())
    }

    #[test]
    fn test_facing_polygons() -> Result<()> {
        // Two squares at z=0, one with normal up, one with normal down
        let pts1 = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
        ];
        let poly1 = Polygon::new("up", pts1, Some(Vector::new(0.0, 0.0, 1.0)))?;

        // Same square but reversed winding (normal down)
        let pts2 = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
        ];
        let poly2 = Polygon::new("down", pts2, Some(Vector::new(0.0, 0.0, -1.0)))?;

        assert!(are_polygons_facing(&poly1, &poly2));
        Ok(())
    }

    #[test]
    fn test_not_facing_same_normal() -> Result<()> {
        // Two squares with same normal direction
        let poly1 = make_xy_square(0.0, 0.0, 1.0, 0.0, "p1")?;
        let poly2 = make_xy_square(0.5, 0.5, 1.0, 0.0, "p2")?;

        assert!(!are_polygons_facing(&poly1, &poly2));
        Ok(())
    }

    #[test]
    fn test_not_facing_no_overlap() -> Result<()> {
        // Two facing squares but no overlap
        let pts1 = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
        ];
        let poly1 = Polygon::new("up", pts1, Some(Vector::new(0.0, 0.0, 1.0)))?;

        let pts2 = vec![
            Point::new(5.0, 5.0, 0.0),
            Point::new(5.0, 6.0, 0.0),
            Point::new(6.0, 6.0, 0.0),
            Point::new(6.0, 5.0, 0.0),
        ];
        let poly2 = Polygon::new("down", pts2, Some(Vector::new(0.0, 0.0, -1.0)))?;

        assert!(!are_polygons_facing(&poly1, &poly2));
        Ok(())
    }

    #[test]
    fn test_crossing_polygons() -> Result<()> {
        // Horizontal square at z=0.5
        let pts1 = vec![
            Point::new(0.0, 0.0, 0.5),
            Point::new(1.0, 0.0, 0.5),
            Point::new(1.0, 1.0, 0.5),
            Point::new(0.0, 1.0, 0.5),
        ];
        let poly1 = Polygon::new("horizontal", pts1, None)?;

        // Vertical square crossing through
        let pts2 = vec![
            Point::new(0.5, 0.5, 0.0),
            Point::new(0.5, 0.5, 1.0),
            Point::new(0.5, -0.5, 1.0),
            Point::new(0.5, -0.5, 0.0),
        ];
        let poly2 = Polygon::new("vertical", pts2, None)?;

        assert!(are_polygons_crossing(&poly1, &poly2));
        Ok(())
    }

    #[test]
    fn test_not_crossing_parallel() -> Result<()> {
        let poly1 = make_xy_square(0.0, 0.0, 1.0, 0.0, "p1")?;
        let poly2 = make_xy_square(0.0, 0.0, 1.0, 1.0, "p2")?;

        assert!(!are_polygons_crossing(&poly1, &poly2));
        Ok(())
    }

    #[test]
    fn test_touching_shared_edge() -> Result<()> {
        // Two adjacent squares sharing an edge
        let poly1 = make_xy_square(0.0, 0.0, 1.0, 0.0, "p1")?;
        let poly2 = make_xy_square(1.0, 0.0, 1.0, 0.0, "p2")?;

        assert!(are_polygons_touching(&poly1, &poly2));
        Ok(())
    }

    #[test]
    fn test_not_touching_separated() -> Result<()> {
        // Two squares with a gap between them
        let poly1 = make_xy_square(0.0, 0.0, 1.0, 0.0, "p1")?;
        let poly2 = make_xy_square(2.0, 0.0, 1.0, 0.0, "p2")?;

        assert!(!are_polygons_touching(&poly1, &poly2));
        Ok(())
    }

    #[test]
    fn test_not_touching_overlapping() -> Result<()> {
        // Two overlapping squares (not just touching)
        let poly1 = make_xy_square(0.0, 0.0, 1.0, 0.0, "p1")?;
        let poly2 = make_xy_square(0.5, 0.0, 1.0, 0.0, "p2")?;

        // They overlap, so they're not just "touching"
        assert!(!are_polygons_touching(&poly1, &poly2));
        Ok(())
    }

    #[test]
    fn test_polygon_distance_touching() -> Result<()> {
        let poly1 = make_xy_square(0.0, 0.0, 1.0, 0.0, "p1")?;
        let poly2 = make_xy_square(1.0, 0.0, 1.0, 0.0, "p2")?;

        let dist = polygon_distance(&poly1, &poly2);
        assert!(dist.is_close(0.0));
        Ok(())
    }

    #[test]
    fn test_polygon_distance_separated() -> Result<()> {
        let poly1 = make_xy_square(0.0, 0.0, 1.0, 0.0, "p1")?;
        let poly2 = make_xy_square(2.0, 0.0, 1.0, 0.0, "p2")?;

        let dist = polygon_distance(&poly1, &poly2);
        assert!(dist.is_close(1.0)); // Gap of 1 unit
        Ok(())
    }

    #[test]
    fn test_overlap_area_overlapping() -> Result<()> {
        let poly1 = make_xy_square(0.0, 0.0, 2.0, 0.0, "p1")?;
        let poly2 = make_xy_square(1.0, 1.0, 2.0, 0.0, "p2")?;

        let area = polygon_overlap_area(&poly1, &poly2);
        assert!(area > 0.0);
        Ok(())
    }

    #[test]
    fn test_overlap_area_no_overlap() -> Result<()> {
        let poly1 = make_xy_square(0.0, 0.0, 1.0, 0.0, "p1")?;
        let poly2 = make_xy_square(5.0, 5.0, 1.0, 0.0, "p2")?;

        let area = polygon_overlap_area(&poly1, &poly2);
        assert!(area.is_close(0.0));
        Ok(())
    }
}
