//! Line segment operations for 3D geometry.
//!
//! This module provides functions for working with line segments in 3D space,
//! including intersection tests, parallelism checks, and polygon crossing detection.

use crate::geom::EPS;
use crate::geom::polygon::Polygon;
use crate::{Point, Vector};

/// Result of a line segment intersection test.
#[derive(Debug, Clone, PartialEq)]
pub enum SegmentIntersection {
    /// Segments intersect at a single point
    Point(Point),
    /// Segments are parallel and do not intersect
    Parallel,
    /// Segments are collinear and overlap (returns the overlap segment)
    Collinear(Point, Point),
    /// Segments are skew (non-parallel, non-intersecting in 3D)
    /// Returns the closest points on each segment
    Skew(Point, Point),
    /// No intersection (segments don't reach each other)
    None,
}

/// Checks if two line segments are parallel.
///
/// Two segments are parallel if their direction vectors are parallel
/// (cross product is zero or near-zero).
pub fn are_segments_parallel(p1: Point, p2: Point, p3: Point, p4: Point) -> bool {
    let d1 = p2 - p1; // Direction of first segment
    let d2 = p4 - p3; // Direction of second segment

    // Segments are parallel if cross product is zero
    let cross = d1.cross(&d2);
    cross.length() < EPS
}

/// Finds the intersection of two line segments in 3D space.
///
/// Returns the type of intersection:
/// - `Point`: segments intersect at a single point
/// - `Parallel`: segments are parallel and don't intersect
/// - `Collinear`: segments are collinear and may overlap
/// - `Skew`: segments are skew lines (closest points returned)
/// - `None`: segments don't intersect
pub fn segment_intersection(p1: Point, p2: Point, p3: Point, p4: Point) -> SegmentIntersection {
    let d1 = p2 - p1; // Direction of segment 1
    let d2 = p4 - p3; // Direction of segment 2
    let r = p3 - p1; // Vector from p1 to p3

    let d1_len = d1.length();
    let d2_len = d2.length();

    // Handle degenerate segments (points)
    if d1_len < EPS && d2_len < EPS {
        // Both segments are points
        if p1.is_close(&p3) {
            return SegmentIntersection::Point(p1);
        }
        return SegmentIntersection::None;
    }

    if d1_len < EPS {
        // First segment is a point, check if it's on second segment
        if is_point_on_segment(p1, p3, p4) {
            return SegmentIntersection::Point(p1);
        }
        return SegmentIntersection::None;
    }

    if d2_len < EPS {
        // Second segment is a point, check if it's on first segment
        if is_point_on_segment(p3, p1, p2) {
            return SegmentIntersection::Point(p3);
        }
        return SegmentIntersection::None;
    }

    let cross_d1_d2 = d1.cross(&d2);
    let cross_len_sq = cross_d1_d2.dot(&cross_d1_d2);

    // Check if segments are parallel
    if cross_len_sq < EPS * EPS {
        // Parallel segments - check if collinear
        let cross_r_d1 = r.cross(&d1);
        if cross_r_d1.length() < EPS {
            // Collinear - check for overlap
            return check_collinear_overlap(p1, p2, p3, p4, &d1);
        }
        return SegmentIntersection::Parallel;
    }

    // Non-parallel segments - find closest points
    // Using the formula for closest points on two lines
    // Line 1: p1 + t * d1
    // Line 2: p3 + s * d2

    let cross_r_d2 = r.cross(&d2);
    let cross_r_d1 = r.cross(&d1);

    let t = cross_r_d2.dot(&cross_d1_d2) / cross_len_sq;
    let s = cross_r_d1.dot(&cross_d1_d2) / cross_len_sq;

    // Points on the infinite lines
    let point_on_line1 = Point::new(p1.x + t * d1.dx, p1.y + t * d1.dy, p1.z + t * d1.dz);
    let point_on_line2 = Point::new(p3.x + s * d2.dx, p3.y + s * d2.dy, p3.z + s * d2.dz);

    // Check if the closest points are actually the same (lines intersect)
    let distance = (point_on_line1 - point_on_line2).length();

    if distance < EPS {
        // Lines intersect - check if within segment bounds
        if (-EPS..=1.0 + EPS).contains(&t) && (-EPS..=1.0 + EPS).contains(&s) {
            return SegmentIntersection::Point(point_on_line1);
        }
        // Intersection point is outside segment bounds
        return SegmentIntersection::None;
    }

    // Lines are skew - find closest points on actual segments
    let t_clamped = t.clamp(0.0, 1.0);
    let s_clamped = s.clamp(0.0, 1.0);

    let closest1 = Point::new(
        p1.x + t_clamped * d1.dx,
        p1.y + t_clamped * d1.dy,
        p1.z + t_clamped * d1.dz,
    );
    let closest2 = Point::new(
        p3.x + s_clamped * d2.dx,
        p3.y + s_clamped * d2.dy,
        p3.z + s_clamped * d2.dz,
    );

    SegmentIntersection::Skew(closest1, closest2)
}

/// Helper to check overlap of collinear segments.
fn check_collinear_overlap(
    p1: Point,
    _p2: Point,
    p3: Point,
    p4: Point,
    d1: &Vector,
) -> SegmentIntersection {
    // Project all points onto the line direction
    let d1_len_sq = d1.dot(d1);
    if d1_len_sq < EPS * EPS {
        return SegmentIntersection::None;
    }

    let t1: f64 = 0.0;
    let t2: f64 = 1.0;
    let t3 = (p3 - p1).dot(d1) / d1_len_sq;
    let t4 = (p4 - p1).dot(d1) / d1_len_sq;

    // Sort the t values for segment 2
    let (t3_min, t3_max) = if t3 < t4 { (t3, t4) } else { (t4, t3) };

    // Find overlap
    let overlap_start = t1.max(t3_min);
    let overlap_end = t2.min(t3_max);

    if overlap_start > overlap_end + EPS {
        return SegmentIntersection::Parallel; // Collinear but no overlap
    }

    if (overlap_end - overlap_start).abs() < EPS {
        // Single point overlap
        let pt = Point::new(
            p1.x + overlap_start * d1.dx,
            p1.y + overlap_start * d1.dy,
            p1.z + overlap_start * d1.dz,
        );
        return SegmentIntersection::Point(pt);
    }

    // Segment overlap
    let start_pt = Point::new(
        p1.x + overlap_start * d1.dx,
        p1.y + overlap_start * d1.dy,
        p1.z + overlap_start * d1.dz,
    );
    let end_pt = Point::new(
        p1.x + overlap_end * d1.dx,
        p1.y + overlap_end * d1.dy,
        p1.z + overlap_end * d1.dz,
    );

    SegmentIntersection::Collinear(start_pt, end_pt)
}

/// Checks if a point lies on a line segment.
fn is_point_on_segment(pt: Point, p1: Point, p2: Point) -> bool {
    pt.is_on_segment(p1, p2)
}

/// Checks if a line segment crosses a polygon.
///
/// A segment crosses a polygon if:
/// 1. The segment intersects the polygon's plane
/// 2. The intersection point lies inside the polygon
///
/// Returns `Some(Point)` with the intersection point if crossing occurs,
/// `None` otherwise.
pub fn segment_crosses_polygon(
    seg_start: Point,
    seg_end: Point,
    polygon: &Polygon,
) -> Option<Point> {
    let (a, b, c, d) = polygon.plane_coefficients();
    let vn = &polygon.vn;

    // Calculate signed distances from segment endpoints to plane
    let dist_start = a * seg_start.x + b * seg_start.y + c * seg_start.z + d;
    let dist_end = a * seg_end.x + b * seg_end.y + c * seg_end.z + d;

    // Check if segment is entirely on one side of the plane
    if dist_start > EPS && dist_end > EPS {
        return None; // Both points above plane
    }
    if dist_start < -EPS && dist_end < -EPS {
        return None; // Both points below plane
    }

    // Check if segment lies in the plane
    if dist_start.abs() < EPS && dist_end.abs() < EPS {
        // Segment is in the plane - check if either endpoint is inside
        if polygon.is_point_inside(seg_start, true) {
            return Some(seg_start);
        }
        if polygon.is_point_inside(seg_end, true) {
            return Some(seg_end);
        }
        // Check if segment crosses polygon edges
        // This is a more complex case - for now, return None
        // TODO: Handle coplanar segment-polygon intersection
        return None;
    }

    // Segment crosses the plane - find intersection point
    // Parametric form: P = seg_start + t * (seg_end - seg_start)
    // Plane equation: a*Px + b*Py + c*Pz + d = 0
    let seg_vec = seg_end - seg_start;
    let denom = vn.dx * seg_vec.dx + vn.dy * seg_vec.dy + vn.dz * seg_vec.dz;

    if denom.abs() < EPS {
        // Segment is parallel to plane (shouldn't happen given earlier checks)
        return None;
    }

    let t = -dist_start / denom;

    // Check if intersection is within segment bounds
    if !(-EPS..=1.0 + EPS).contains(&t) {
        return None;
    }

    let intersection = Point::new(
        seg_start.x + t * seg_vec.dx,
        seg_start.y + t * seg_vec.dy,
        seg_start.z + t * seg_vec.dz,
    );

    // Check if intersection point is inside the polygon
    if polygon.is_point_inside(intersection, true) {
        Some(intersection)
    } else {
        None
    }
}

/// Calculates the distance between a point and a line segment.
///
/// Returns the minimum distance from the point to any point on the segment.
pub fn distance_point_to_segment(pt: Point, p1: Point, p2: Point) -> f64 {
    let seg_vec = p2 - p1;
    let pt_vec = pt - p1;

    let seg_len_sq = seg_vec.dot(&seg_vec);

    if seg_len_sq < EPS * EPS {
        // Segment is a point
        return pt_vec.length();
    }

    // Project pt onto the line, clamped to segment
    let t = (pt_vec.dot(&seg_vec) / seg_len_sq).clamp(0.0, 1.0);

    let closest = Point::new(
        p1.x + t * seg_vec.dx,
        p1.y + t * seg_vec.dy,
        p1.z + t * seg_vec.dz,
    );

    (pt - closest).length()
}

/// Calculates the distance between a point and an infinite line.
///
/// The line is defined by two points p1 and p2.
pub fn distance_point_to_line(pt: Point, p1: Point, p2: Point) -> f64 {
    let line_vec = p2 - p1;
    let pt_vec = pt - p1;

    let line_len_sq = line_vec.dot(&line_vec);

    if line_len_sq < EPS * EPS {
        // Line is a point
        return pt_vec.length();
    }

    // Cross product gives area of parallelogram
    let cross = pt_vec.cross(&line_vec);

    // Distance = area / base = |cross| / |line_vec|
    cross.length() / line_len_sq.sqrt()
}

/// Finds the closest point on a line to a given point.
///
/// The line is defined by two points p1 and p2.
pub fn closest_point_on_line(pt: Point, p1: Point, p2: Point) -> Point {
    let line_vec = p2 - p1;
    let pt_vec = pt - p1;

    let line_len_sq = line_vec.dot(&line_vec);

    if line_len_sq < EPS * EPS {
        // Line is a point
        return p1;
    }

    let t = pt_vec.dot(&line_vec) / line_len_sq;

    Point::new(
        p1.x + t * line_vec.dx,
        p1.y + t * line_vec.dy,
        p1.z + t * line_vec.dz,
    )
}

/// Finds the closest point on a segment to a given point.
///
/// The segment is defined by two endpoints p1 and p2.
pub fn closest_point_on_segment(pt: Point, p1: Point, p2: Point) -> Point {
    let seg_vec = p2 - p1;
    let pt_vec = pt - p1;

    let seg_len_sq = seg_vec.dot(&seg_vec);

    if seg_len_sq < EPS * EPS {
        // Segment is a point
        return p1;
    }

    let t = (pt_vec.dot(&seg_vec) / seg_len_sq).clamp(0.0, 1.0);

    Point::new(
        p1.x + t * seg_vec.dx,
        p1.y + t * seg_vec.dy,
        p1.z + t * seg_vec.dz,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::IsClose;

    #[test]
    fn test_parallel_segments() {
        // Parallel segments along X axis
        let p1 = Point::new(0.0, 0.0, 0.0);
        let p2 = Point::new(1.0, 0.0, 0.0);
        let p3 = Point::new(0.0, 1.0, 0.0);
        let p4 = Point::new(1.0, 1.0, 0.0);

        assert!(are_segments_parallel(p1, p2, p3, p4));
    }

    #[test]
    fn test_non_parallel_segments() {
        // Perpendicular segments
        let p1 = Point::new(0.0, 0.0, 0.0);
        let p2 = Point::new(1.0, 0.0, 0.0);
        let p3 = Point::new(0.5, -0.5, 0.0);
        let p4 = Point::new(0.5, 0.5, 0.0);

        assert!(!are_segments_parallel(p1, p2, p3, p4));
    }

    #[test]
    fn test_segment_intersection_crossing() {
        // Two segments that cross at (0.5, 0.5, 0)
        let p1 = Point::new(0.0, 0.0, 0.0);
        let p2 = Point::new(1.0, 1.0, 0.0);
        let p3 = Point::new(0.0, 1.0, 0.0);
        let p4 = Point::new(1.0, 0.0, 0.0);

        let result = segment_intersection(p1, p2, p3, p4);
        match result {
            SegmentIntersection::Point(pt) => {
                assert!(pt.is_close(&Point::new(0.5, 0.5, 0.0)));
            }
            _ => panic!("Expected Point intersection, got {:?}", result),
        }
    }

    #[test]
    fn test_segment_intersection_no_crossing() {
        // Two segments that don't intersect
        let p1 = Point::new(0.0, 0.0, 0.0);
        let p2 = Point::new(1.0, 0.0, 0.0);
        let p3 = Point::new(0.0, 1.0, 0.0);
        let p4 = Point::new(1.0, 1.0, 0.0);

        let result = segment_intersection(p1, p2, p3, p4);
        assert!(matches!(result, SegmentIntersection::Parallel));
    }

    #[test]
    fn test_segment_intersection_collinear_overlap() {
        // Two collinear overlapping segments
        let p1 = Point::new(0.0, 0.0, 0.0);
        let p2 = Point::new(2.0, 0.0, 0.0);
        let p3 = Point::new(1.0, 0.0, 0.0);
        let p4 = Point::new(3.0, 0.0, 0.0);

        let result = segment_intersection(p1, p2, p3, p4);
        match result {
            SegmentIntersection::Collinear(start, end) => {
                assert!(start.is_close(&Point::new(1.0, 0.0, 0.0)));
                assert!(end.is_close(&Point::new(2.0, 0.0, 0.0)));
            }
            _ => panic!("Expected Collinear intersection, got {:?}", result),
        }
    }

    #[test]
    fn test_segment_intersection_collinear_no_overlap() {
        // Two collinear non-overlapping segments
        let p1 = Point::new(0.0, 0.0, 0.0);
        let p2 = Point::new(1.0, 0.0, 0.0);
        let p3 = Point::new(2.0, 0.0, 0.0);
        let p4 = Point::new(3.0, 0.0, 0.0);

        let result = segment_intersection(p1, p2, p3, p4);
        assert!(matches!(result, SegmentIntersection::Parallel));
    }

    #[test]
    fn test_segment_intersection_skew() {
        // Two skew segments in 3D
        let p1 = Point::new(0.0, 0.0, 0.0);
        let p2 = Point::new(1.0, 0.0, 0.0);
        let p3 = Point::new(0.5, 0.5, 0.5);
        let p4 = Point::new(0.5, 0.5, 1.5);

        let result = segment_intersection(p1, p2, p3, p4);
        match result {
            SegmentIntersection::Skew(_, _) => {}
            _ => panic!("Expected Skew intersection, got {:?}", result),
        }
    }

    #[test]
    fn test_segment_intersection_t_shape() {
        // T-shaped intersection: one segment ends at the middle of another
        let p1 = Point::new(0.0, 0.0, 0.0);
        let p2 = Point::new(2.0, 0.0, 0.0);
        let p3 = Point::new(1.0, 0.0, 0.0);
        let p4 = Point::new(1.0, 1.0, 0.0);

        let result = segment_intersection(p1, p2, p3, p4);
        match result {
            SegmentIntersection::Point(pt) => {
                assert!(pt.is_close(&Point::new(1.0, 0.0, 0.0)));
            }
            _ => panic!("Expected Point intersection, got {:?}", result),
        }
    }

    #[test]
    fn test_distance_point_to_segment() {
        let p1 = Point::new(0.0, 0.0, 0.0);
        let p2 = Point::new(2.0, 0.0, 0.0);

        // Point directly above middle of segment
        let pt = Point::new(1.0, 1.0, 0.0);
        let dist = distance_point_to_segment(pt, p1, p2);
        assert!(dist.is_close(1.0));

        // Point beyond segment end
        let pt = Point::new(3.0, 0.0, 0.0);
        let dist = distance_point_to_segment(pt, p1, p2);
        assert!(dist.is_close(1.0));

        // Point at segment start
        let pt = Point::new(0.0, 0.0, 0.0);
        let dist = distance_point_to_segment(pt, p1, p2);
        assert!(dist.is_close(0.0));
    }

    #[test]
    fn test_distance_point_to_line() {
        let p1 = Point::new(0.0, 0.0, 0.0);
        let p2 = Point::new(2.0, 0.0, 0.0);

        // Point above line
        let pt = Point::new(1.0, 1.0, 0.0);
        let dist = distance_point_to_line(pt, p1, p2);
        assert!(dist.is_close(1.0));

        // Point beyond segment but on extended line
        let pt = Point::new(5.0, 0.0, 0.0);
        let dist = distance_point_to_line(pt, p1, p2);
        assert!(dist.is_close(0.0));
    }

    #[test]
    fn test_closest_point_on_line() {
        let p1 = Point::new(0.0, 0.0, 0.0);
        let p2 = Point::new(2.0, 0.0, 0.0);

        // Point above line
        let pt = Point::new(1.0, 1.0, 0.0);
        let closest = closest_point_on_line(pt, p1, p2);
        assert!(closest.is_close(&Point::new(1.0, 0.0, 0.0)));

        // Point beyond segment
        let pt = Point::new(5.0, 1.0, 0.0);
        let closest = closest_point_on_line(pt, p1, p2);
        assert!(closest.is_close(&Point::new(5.0, 0.0, 0.0)));
    }

    #[test]
    fn test_closest_point_on_segment() {
        let p1 = Point::new(0.0, 0.0, 0.0);
        let p2 = Point::new(2.0, 0.0, 0.0);

        // Point above segment
        let pt = Point::new(1.0, 1.0, 0.0);
        let closest = closest_point_on_segment(pt, p1, p2);
        assert!(closest.is_close(&Point::new(1.0, 0.0, 0.0)));

        // Point beyond segment - should clamp to endpoint
        let pt = Point::new(5.0, 1.0, 0.0);
        let closest = closest_point_on_segment(pt, p1, p2);
        assert!(closest.is_close(&Point::new(2.0, 0.0, 0.0)));
    }

    #[test]
    fn test_segment_crosses_polygon() -> anyhow::Result<()> {
        // Create a square polygon in the XY plane at z=0
        let pts = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
        ];
        let polygon = Polygon::new("square", pts, None)?;

        // Segment that crosses through the center
        let seg_start = Point::new(0.5, 0.5, -1.0);
        let seg_end = Point::new(0.5, 0.5, 1.0);
        let result = segment_crosses_polygon(seg_start, seg_end, &polygon);
        assert!(result.is_some());
        let intersection = result.unwrap();
        assert!(intersection.is_close(&Point::new(0.5, 0.5, 0.0)));

        Ok(())
    }

    #[test]
    fn test_segment_misses_polygon() -> anyhow::Result<()> {
        // Create a square polygon in the XY plane at z=0
        let pts = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
        ];
        let polygon = Polygon::new("square", pts, None)?;

        // Segment that misses the polygon (outside bounds)
        let seg_start = Point::new(2.0, 2.0, -1.0);
        let seg_end = Point::new(2.0, 2.0, 1.0);
        let result = segment_crosses_polygon(seg_start, seg_end, &polygon);
        assert!(result.is_none());

        // Segment parallel to polygon plane (above it)
        let seg_start = Point::new(0.0, 0.5, 1.0);
        let seg_end = Point::new(1.0, 0.5, 1.0);
        let result = segment_crosses_polygon(seg_start, seg_end, &polygon);
        assert!(result.is_none());

        Ok(())
    }

    #[test]
    fn test_segment_touches_polygon_edge() -> anyhow::Result<()> {
        // Create a square polygon in the XY plane at z=0
        let pts = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
        ];
        let polygon = Polygon::new("square", pts, None)?;

        // Segment that crosses at edge
        let seg_start = Point::new(0.5, 0.0, -1.0);
        let seg_end = Point::new(0.5, 0.0, 1.0);
        let result = segment_crosses_polygon(seg_start, seg_end, &polygon);
        assert!(result.is_some());

        Ok(())
    }
}
