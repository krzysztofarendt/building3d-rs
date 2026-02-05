//! Polygon slicing operations.
//!
//! This module provides functions for cutting polygons along a line,
//! producing two new polygons.

use crate::geom::segment::{SegmentIntersection, segment_intersection};
use crate::{Point, Polygon};
use anyhow::{Result, anyhow};

/// Result of slicing a polygon.
#[derive(Debug, Clone)]
pub struct SliceResult {
    /// The first polygon (on one side of the slice line)
    pub poly1: Polygon,
    /// The second polygon (on the other side of the slice line)
    pub poly2: Polygon,
}

/// Slices a polygon along a line defined by two points.
///
/// The slice line must intersect the polygon at exactly two points
/// (either on edges or at vertices). The function returns two new
/// polygons created by the slice.
///
/// # Arguments
/// * `polygon` - The polygon to slice
/// * `slice_start` - First point defining the slice line
/// * `slice_end` - Second point defining the slice line
///
/// # Returns
/// A `SliceResult` containing two new polygons, or an error if the
/// slice is invalid (e.g., doesn't intersect polygon at exactly 2 points).
///
/// # Example
/// ```
/// use building3d::{Point, Polygon};
/// use building3d::geom::polygon::slice::slice_polygon;
///
/// let pts = vec![
///     Point::new(0.0, 0.0, 0.0),
///     Point::new(2.0, 0.0, 0.0),
///     Point::new(2.0, 2.0, 0.0),
///     Point::new(0.0, 2.0, 0.0),
/// ];
/// let poly = Polygon::new("square", pts, None).unwrap();
///
/// // Slice horizontally through the middle
/// let result = slice_polygon(&poly, Point::new(-1.0, 1.0, 0.0), Point::new(3.0, 1.0, 0.0)).unwrap();
/// // result.poly1 and result.poly2 are the two halves
/// ```
pub fn slice_polygon(
    polygon: &Polygon,
    slice_start: Point,
    slice_end: Point,
) -> Result<SliceResult> {
    let vertices = polygon.vertices();
    let n = vertices.len();

    if n < 3 {
        return Err(anyhow!("Polygon must have at least 3 vertices"));
    }

    // Find intersection points of the slice line with polygon edges
    let mut intersections: Vec<(usize, Point)> = Vec::new();

    for i in 0..n {
        let p1 = vertices[i];
        let p2 = vertices[(i + 1) % n];

        match segment_intersection(slice_start, slice_end, p1, p2) {
            SegmentIntersection::Point(pt) => {
                // Check if this intersection is actually on the polygon edge
                // (segment_intersection checks if it's on the infinite line)
                if is_point_on_segment(pt, p1, p2) {
                    // Avoid duplicate intersections at vertices
                    let is_duplicate = intersections
                        .iter()
                        .any(|(_, existing)| existing.is_close(&pt));
                    if !is_duplicate {
                        intersections.push((i, pt));
                    }
                }
            }
            SegmentIntersection::Collinear(start, end) => {
                // Edge is collinear with slice line - use endpoints
                let is_dup1 = intersections
                    .iter()
                    .any(|(_, existing)| existing.is_close(&start));
                let is_dup2 = intersections
                    .iter()
                    .any(|(_, existing)| existing.is_close(&end));
                if !is_dup1 {
                    intersections.push((i, start));
                }
                if !is_dup2 {
                    intersections.push((i, end));
                }
            }
            _ => {}
        }
    }

    // We need exactly 2 intersection points to create a valid slice
    if intersections.len() != 2 {
        return Err(anyhow!(
            "Slice line must intersect polygon at exactly 2 points, found {}",
            intersections.len()
        ));
    }

    // Sort intersections by edge index
    intersections.sort_by_key(|(idx, _)| *idx);

    let (idx1, pt1) = intersections[0];
    let (idx2, pt2) = intersections[1];

    // Build the two new polygons
    // Polygon 1: from pt1, along edges to pt2
    // Polygon 2: from pt2, along edges back to pt1

    let mut poly1_pts: Vec<Point> = Vec::new();
    let mut poly2_pts: Vec<Point> = Vec::new();

    // Start polygon 1 with first intersection point
    poly1_pts.push(pt1);

    // Add vertices from idx1+1 to idx2 (inclusive of vertices after idx1, up to idx2)
    let mut i = (idx1 + 1) % n;
    while i != (idx2 + 1) % n {
        // Skip if vertex is same as intersection point
        if !vertices[i].is_close(&pt1) && !vertices[i].is_close(&pt2) {
            poly1_pts.push(vertices[i]);
        }
        i = (i + 1) % n;
    }

    // Add second intersection point if not already at a vertex
    if !poly1_pts.last().is_some_and(|p| p.is_close(&pt2)) {
        poly1_pts.push(pt2);
    }

    // Start polygon 2 with second intersection point
    poly2_pts.push(pt2);

    // Add vertices from idx2+1 back around to idx1
    let mut i = (idx2 + 1) % n;
    while i != (idx1 + 1) % n {
        // Skip if vertex is same as intersection point
        if !vertices[i].is_close(&pt1) && !vertices[i].is_close(&pt2) {
            poly2_pts.push(vertices[i]);
        }
        i = (i + 1) % n;
    }

    // Add first intersection point if not already at a vertex
    if !poly2_pts.last().is_some_and(|p| p.is_close(&pt1)) {
        poly2_pts.push(pt1);
    }

    // Remove duplicate consecutive points
    poly1_pts = remove_consecutive_duplicates(poly1_pts);
    poly2_pts = remove_consecutive_duplicates(poly2_pts);

    // Validate we have enough points for polygons
    if poly1_pts.len() < 3 {
        return Err(anyhow!("Resulting polygon 1 has fewer than 3 vertices"));
    }
    if poly2_pts.len() < 3 {
        return Err(anyhow!("Resulting polygon 2 has fewer than 3 vertices"));
    }

    // Create the new polygons
    let poly1 = Polygon::new(&format!("{}_1", polygon.name), poly1_pts, Some(polygon.vn))?;
    let poly2 = Polygon::new(&format!("{}_2", polygon.name), poly2_pts, Some(polygon.vn))?;

    Ok(SliceResult { poly1, poly2 })
}

/// Checks if a point lies on a line segment (within tolerance).
fn is_point_on_segment(pt: Point, seg_start: Point, seg_end: Point) -> bool {
    let seg_vec = seg_end - seg_start;
    let pt_vec = pt - seg_start;

    let seg_len_sq = seg_vec.dot(&seg_vec);
    if seg_len_sq < 1e-20 {
        // Degenerate segment
        return pt.is_close(&seg_start);
    }

    let t = pt_vec.dot(&seg_vec) / seg_len_sq;

    // Check if t is in [0, 1] with some tolerance
    if !(-1e-10..=1.0 + 1e-10).contains(&t) {
        return false;
    }

    // Check if point is actually on the line
    let projected = seg_start + seg_vec * t.clamp(0.0, 1.0);
    pt.is_close(&projected)
}

/// Removes consecutive duplicate points from a list.
fn remove_consecutive_duplicates(mut pts: Vec<Point>) -> Vec<Point> {
    if pts.len() < 2 {
        return pts;
    }

    let mut result: Vec<Point> = Vec::with_capacity(pts.len());
    result.push(pts[0]);

    for pt in pts.drain(1..) {
        if !result.last().unwrap().is_close(&pt) {
            result.push(pt);
        }
    }

    // Also check if first and last are the same
    if result.len() > 1 && result.first().unwrap().is_close(result.last().unwrap()) {
        result.pop();
    }

    result
}

/// Slices a polygon at specific points that lie on its edges.
///
/// This is a convenience function when you already know the exact
/// intersection points on the polygon boundary.
///
/// # Arguments
/// * `polygon` - The polygon to slice
/// * `pt1` - First slice point (must be on polygon boundary)
/// * `pt2` - Second slice point (must be on polygon boundary)
pub fn slice_polygon_at_points(polygon: &Polygon, pt1: Point, pt2: Point) -> Result<SliceResult> {
    // Use the points directly as the slice line
    slice_polygon(polygon, pt1, pt2)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_square(size: f64) -> Result<Polygon> {
        let pts = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(size, 0.0, 0.0),
            Point::new(size, size, 0.0),
            Point::new(0.0, size, 0.0),
        ];
        Polygon::new("square", pts, None)
    }

    #[test]
    fn test_slice_square_horizontally() -> Result<()> {
        let square = make_square(2.0)?;

        // Slice horizontally through the middle (y = 1)
        let result = slice_polygon(
            &square,
            Point::new(-1.0, 1.0, 0.0),
            Point::new(3.0, 1.0, 0.0),
        )?;

        // Both polygons should be valid
        assert!(result.poly1.vertices().len() >= 3);
        assert!(result.poly2.vertices().len() >= 3);

        // Combined area should equal original
        let original_area = square.area();
        let combined_area = result.poly1.area() + result.poly2.area();
        assert!(
            (original_area - combined_area).abs() < 0.01,
            "Area mismatch: {} vs {}",
            original_area,
            combined_area
        );

        // Each half should be approximately half the original area
        assert!((result.poly1.area() - 2.0).abs() < 0.1);
        assert!((result.poly2.area() - 2.0).abs() < 0.1);

        Ok(())
    }

    #[test]
    fn test_slice_square_vertically() -> Result<()> {
        let square = make_square(2.0)?;

        // Slice vertically through the middle (x = 1)
        let result = slice_polygon(
            &square,
            Point::new(1.0, -1.0, 0.0),
            Point::new(1.0, 3.0, 0.0),
        )?;

        // Both polygons should be valid rectangles
        assert!(result.poly1.vertices().len() >= 3);
        assert!(result.poly2.vertices().len() >= 3);

        // Combined area should equal original
        let original_area = square.area();
        let combined_area = result.poly1.area() + result.poly2.area();
        assert!(
            (original_area - combined_area).abs() < 0.01,
            "Area mismatch: {} vs {}",
            original_area,
            combined_area
        );

        Ok(())
    }

    #[test]
    fn test_slice_square_diagonally() -> Result<()> {
        let square = make_square(2.0)?;

        // Slice diagonally from corner to corner
        let result = slice_polygon(
            &square,
            Point::new(-1.0, -1.0, 0.0),
            Point::new(3.0, 3.0, 0.0),
        )?;

        // Both polygons should be triangles
        assert_eq!(result.poly1.vertices().len(), 3);
        assert_eq!(result.poly2.vertices().len(), 3);

        // Each triangle should be half the original area
        let half_area = square.area() / 2.0;
        assert!((result.poly1.area() - half_area).abs() < 0.1);
        assert!((result.poly2.area() - half_area).abs() < 0.1);

        Ok(())
    }

    #[test]
    fn test_slice_through_vertex() -> Result<()> {
        let square = make_square(2.0)?;

        // Slice from vertex (0,0) to midpoint of opposite edge (2, 1)
        let result = slice_polygon(
            &square,
            Point::new(0.0, 0.0, 0.0),
            Point::new(2.0, 1.0, 0.0),
        )?;

        // Should produce two valid polygons
        assert!(result.poly1.vertices().len() >= 3);
        assert!(result.poly2.vertices().len() >= 3);

        // Combined area should equal original
        let original_area = square.area();
        let combined_area = result.poly1.area() + result.poly2.area();
        assert!(
            (original_area - combined_area).abs() < 0.01,
            "Area mismatch: {} vs {}",
            original_area,
            combined_area
        );

        Ok(())
    }

    #[test]
    fn test_slice_complex_polygon() -> Result<()> {
        // L-shaped polygon
        let pts = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(2.0, 0.0, 0.0),
            Point::new(2.0, 1.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
            Point::new(1.0, 2.0, 0.0),
            Point::new(0.0, 2.0, 0.0),
        ];
        let l_shape = Polygon::new("l_shape", pts, None)?;
        let original_area = l_shape.area();

        // Slice horizontally at y = 0.5
        let result = slice_polygon(
            &l_shape,
            Point::new(-1.0, 0.5, 0.0),
            Point::new(3.0, 0.5, 0.0),
        )?;

        // Both parts should be valid
        assert!(result.poly1.vertices().len() >= 3);
        assert!(result.poly2.vertices().len() >= 3);

        // Combined area should equal original
        let combined_area = result.poly1.area() + result.poly2.area();
        assert!(
            (original_area - combined_area).abs() < 0.01,
            "Area mismatch: {} vs {}",
            original_area,
            combined_area
        );

        Ok(())
    }

    #[test]
    fn test_slice_invalid_no_intersection() {
        let square = make_square(2.0).unwrap();

        // Slice line that doesn't intersect the polygon
        let result = slice_polygon(
            &square,
            Point::new(5.0, 0.0, 0.0),
            Point::new(5.0, 2.0, 0.0),
        );

        assert!(result.is_err());
    }

    #[test]
    fn test_slice_invalid_single_intersection() {
        let square = make_square(2.0).unwrap();

        // Slice line that only touches one vertex
        let result = slice_polygon(
            &square,
            Point::new(0.0, 0.0, 0.0),
            Point::new(-1.0, -1.0, 0.0),
        );

        // Should fail - only 1 intersection point
        assert!(result.is_err());
    }

    #[test]
    fn test_slice_preserves_normal() -> Result<()> {
        let square = make_square(2.0)?;
        let original_normal = square.vn;

        let result = slice_polygon(
            &square,
            Point::new(-1.0, 1.0, 0.0),
            Point::new(3.0, 1.0, 0.0),
        )?;

        // Both new polygons should have the same normal as the original
        assert!(result.poly1.vn.is_close(&original_normal));
        assert!(result.poly2.vn.is_close(&original_normal));

        Ok(())
    }
}
