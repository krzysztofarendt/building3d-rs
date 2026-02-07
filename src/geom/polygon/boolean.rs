//! Boolean operations on polygons.
//!
//! This module provides functions for computing boolean operations
//! (intersection, union, difference) on coplanar polygons.
//!
//! These operations are essential for:
//! - Creating windows/openings (subtraction)
//! - Merging adjacent polygons (union)
//! - Finding overlap areas (intersection)

use crate::geom::projection::PlaneBasis;
use crate::{Point, Polygon, Vector};
use anyhow::{Result, anyhow};

/// Result of a boolean operation that may produce multiple polygons.
#[derive(Debug, Clone)]
pub struct BooleanResult {
    /// Resulting polygons from the operation
    pub polygons: Vec<Polygon>,
}

impl BooleanResult {
    /// Creates a new boolean result.
    pub fn new(polygons: Vec<Polygon>) -> Self {
        Self { polygons }
    }

    /// Returns true if the result is empty (no polygons).
    pub fn is_empty(&self) -> bool {
        self.polygons.is_empty()
    }

    /// Returns the number of resulting polygons.
    pub fn len(&self) -> usize {
        self.polygons.len()
    }

    /// Returns the total area of all resulting polygons.
    pub fn total_area(&self) -> f64 {
        self.polygons.iter().map(|p| p.area()).sum()
    }
}

/// Computes the intersection of two coplanar polygons.
///
/// Returns the polygon representing the area common to both input polygons.
/// Uses the Sutherland-Hodgman algorithm for polygon clipping.
///
/// # Arguments
/// * `poly1` - First polygon
/// * `poly2` - Second polygon (used as the clipping polygon)
///
/// # Returns
/// A `BooleanResult` containing the intersection polygon(s), or an empty
/// result if the polygons don't intersect.
pub fn polygon_intersection(poly1: &Polygon, poly2: &Polygon) -> Result<BooleanResult> {
    // Verify polygons are coplanar
    if !are_polygons_coplanar(poly1, poly2) {
        return Err(anyhow!("Polygons must be coplanar for boolean operations"));
    }

    // Get vertices
    let subject = poly1.vertices();
    let clip_poly = poly2.vertices();

    // Sutherland-Hodgman clipping
    let result_pts = sutherland_hodgman(subject, clip_poly);

    if result_pts.len() < 3 {
        return Ok(BooleanResult::new(vec![]));
    }

    // Create result polygon
    let result_poly = Polygon::new(
        &format!("{}_intersect_{}", poly1.name, poly2.name),
        result_pts,
        Some(poly1.vn),
    )?;

    Ok(BooleanResult::new(vec![result_poly]))
}

/// Computes the difference of two coplanar polygons (poly1 - poly2).
///
/// Returns the part of poly1 that is NOT inside poly2.
/// This is useful for creating holes/openings in polygons.
///
/// Note: This is a simplified implementation that works best when poly2
/// is entirely inside poly1 or they don't overlap at all.
///
/// # Arguments
/// * `poly1` - The polygon to subtract from
/// * `poly2` - The polygon to subtract (the "hole")
///
/// # Returns
/// A `BooleanResult` containing the resulting polygon(s).
pub fn polygon_difference(poly1: &Polygon, poly2: &Polygon) -> Result<BooleanResult> {
    // Verify polygons are coplanar
    if !are_polygons_coplanar(poly1, poly2) {
        return Err(anyhow!("Polygons must be coplanar for boolean operations"));
    }

    // Check if there's any intersection
    let intersection = polygon_intersection(poly1, poly2)?;

    if intersection.is_empty() {
        // No overlap - return original polygon
        return Ok(BooleanResult::new(vec![poly1.clone()]));
    }

    // For a simple implementation, we'll check if poly2 is entirely inside poly1
    // If so, this creates a "hole" which we represent as the boundary vertices
    // connected to the hole vertices

    let poly2_inside = poly2
        .vertices()
        .iter()
        .all(|p| poly1.is_point_inside(*p, true));

    if poly2_inside {
        // Create a polygon with a "bridge" to the hole
        // This is a simplified approach - true polygon-with-holes requires
        // a more complex data structure
        let result = create_polygon_with_hole(poly1, poly2)?;
        return Ok(BooleanResult::new(vec![result]));
    }

    // For partial overlap, return the original minus the intersection
    // This is a simplified approximation
    let remaining_area = poly1.area() - intersection.total_area();

    if remaining_area < 1e-10 {
        // poly1 is entirely contained in poly2
        return Ok(BooleanResult::new(vec![]));
    }

    // For complex cases, return the original with a note that this is approximate
    // A full implementation would require proper polygon clipping
    Ok(BooleanResult::new(vec![poly1.clone()]))
}

/// Sutherland-Hodgman polygon clipping algorithm.
/// Works by projecting to 2D, clipping, then back to 3D.
fn sutherland_hodgman(subject: &[Point], clip_poly: &[Point]) -> Vec<Point> {
    if subject.is_empty() || clip_poly.len() < 3 {
        return vec![];
    }

    let Some(basis) = PlaneBasis::from_polygon(clip_poly) else {
        return vec![];
    };

    // Project to 2D
    let subject_2d: Vec<(f64, f64)> = subject.iter().map(|p| basis.project(*p)).collect();
    let clip_2d: Vec<(f64, f64)> = clip_poly.iter().map(|p| basis.project(*p)).collect();

    // Run 2D Sutherland-Hodgman
    let result_2d = sutherland_hodgman_2d(&subject_2d, &clip_2d);

    if result_2d.len() < 3 {
        return vec![];
    }

    // Project back to 3D
    result_2d
        .iter()
        .map(|&(u, v)| basis.unproject(u, v))
        .collect()
}

/// 2D Sutherland-Hodgman polygon clipping algorithm.
fn sutherland_hodgman_2d(subject: &[(f64, f64)], clip_poly: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let mut output = subject.to_vec();

    for i in 0..clip_poly.len() {
        if output.is_empty() {
            break;
        }

        let edge_start = clip_poly[i];
        let edge_end = clip_poly[(i + 1) % clip_poly.len()];

        let input = output;
        output = Vec::new();

        for j in 0..input.len() {
            let current = input[j];
            let previous = input[(j + input.len() - 1) % input.len()];

            let curr_inside = is_inside_edge_2d(current, edge_start, edge_end);
            let prev_inside = is_inside_edge_2d(previous, edge_start, edge_end);

            if curr_inside {
                if !prev_inside
                    && let Some(intersection) =
                        line_intersection_2d(previous, current, edge_start, edge_end)
                {
                    output.push(intersection);
                }
                output.push(current);
            } else if prev_inside
                && let Some(intersection) =
                    line_intersection_2d(previous, current, edge_start, edge_end)
            {
                output.push(intersection);
            }
        }
    }

    // Remove duplicates
    let mut result = Vec::with_capacity(output.len());
    for pt in output {
        if result.is_empty() || !is_close_2d(result.last().unwrap(), &pt) {
            result.push(pt);
        }
    }
    if result.len() > 1 && is_close_2d(result.first().unwrap(), result.last().unwrap()) {
        result.pop();
    }

    result
}

/// Checks if a 2D point is on the inside of an edge (left side for CCW).
fn is_inside_edge_2d(point: (f64, f64), edge_start: (f64, f64), edge_end: (f64, f64)) -> bool {
    // Cross product z-component: (edge) x (point - edge_start)
    let edge_x = edge_end.0 - edge_start.0;
    let edge_y = edge_end.1 - edge_start.1;
    let to_point_x = point.0 - edge_start.0;
    let to_point_y = point.1 - edge_start.1;

    let cross_z = edge_x * to_point_y - edge_y * to_point_x;

    // For CCW polygons, inside is on the left (cross >= 0)
    cross_z >= -1e-10
}

/// Computes intersection of line (p1,p2) with line (p3,p4).
/// For Sutherland-Hodgman, we need intersection with the infinite edge line.
fn line_intersection_2d(
    p1: (f64, f64),
    p2: (f64, f64), // subject edge (segment)
    p3: (f64, f64),
    p4: (f64, f64), // clip edge (treated as infinite line)
) -> Option<(f64, f64)> {
    let d1x = p2.0 - p1.0;
    let d1y = p2.1 - p1.1;
    let d2x = p4.0 - p3.0;
    let d2y = p4.1 - p3.1;

    let cross = d1x * d2y - d1y * d2x;

    if cross.abs() < 1e-10 {
        return None; // Parallel
    }

    // Compute t for the subject segment (p1->p2)
    let d3x = p3.0 - p1.0;
    let d3y = p3.1 - p1.1;

    let t = (d3x * d2y - d3y * d2x) / cross;

    // For Sutherland-Hodgman, we want intersection with the clip edge line,
    // but t should be in [0, 1] for the subject segment
    if (-1e-10..=1.0 + 1e-10).contains(&t) {
        Some((p1.0 + t * d1x, p1.1 + t * d1y))
    } else {
        None
    }
}

/// Checks if two 2D points are close.
fn is_close_2d(a: &(f64, f64), b: &(f64, f64)) -> bool {
    (a.0 - b.0).abs() < 1e-10 && (a.1 - b.1).abs() < 1e-10
}

/// Checks if two polygons are coplanar.
fn are_polygons_coplanar(poly1: &Polygon, poly2: &Polygon) -> bool {
    // Check if normals are parallel
    let n1 = poly1.vn;
    let n2 = poly2.vn;

    let cross = n1.cross(&n2);
    if cross.length() > 0.01 {
        return false;
    }

    // Check if a point from poly2 lies on poly1's plane
    let (a, b, c, d) = poly1.plane_coefficients();
    let p = poly2.vertices()[0];

    let dist = (a * p.x + b * p.y + c * p.z + d).abs();
    dist < 0.01
}

/// Creates a polygon with a hole by merging the hole into the outer boundary.
fn create_polygon_with_hole(outer: &Polygon, hole: &Polygon) -> Result<Polygon> {
    let outer_verts = outer.vertices().to_vec();
    let mut hole_verts = hole.vertices().to_vec();

    if outer_verts.is_empty() || hole_verts.is_empty() {
        return Err(anyhow!("Cannot create polygon with hole: empty vertices"));
    }

    // `Polygon::with_holes()` requires hole winding opposite to the outer boundary
    // (CW vs outer CCW w.r.t the polygon normal).
    let outer_signed = ring_signed_area(&outer_verts, &outer.vn);
    let hole_signed = ring_signed_area(&hole_verts, &outer.vn);
    if outer_signed.abs() < crate::geom::EPS || hole_signed.abs() < crate::geom::EPS {
        return Err(anyhow!("Cannot create polygon with hole: degenerate ring"));
    }
    if outer_signed.signum() == hole_signed.signum() {
        hole_verts.reverse();
    }

    Polygon::with_holes(
        &format!("{}_with_hole", outer.name),
        outer_verts,
        vec![hole_verts],
        Some(outer.vn),
    )
}

fn ring_signed_area(pts: &[Point], vn: &Vector) -> f64 {
    let n = pts.len();
    if n < 3 {
        return 0.0;
    }
    let mut cross_sum = Vector::new(0.0, 0.0, 0.0);
    for i in 0..n {
        let v1 = Vector::from_a_point(pts[i]);
        let v2 = Vector::from_a_point(pts[(i + 1) % n]);
        cross_sum = cross_sum + v1.cross(&v2);
    }
    0.5 * cross_sum.dot(vn)
}

/// Computes overlap area between two coplanar polygons.
///
/// This is a convenience function that returns just the area
/// of the intersection.
pub fn polygon_overlap_area(poly1: &Polygon, poly2: &Polygon) -> Result<f64> {
    let intersection = polygon_intersection(poly1, poly2)?;
    Ok(intersection.total_area())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_square(size: f64, origin: (f64, f64)) -> Result<Polygon> {
        let pts = vec![
            Point::new(origin.0, origin.1, 0.0),
            Point::new(origin.0 + size, origin.1, 0.0),
            Point::new(origin.0 + size, origin.1 + size, 0.0),
            Point::new(origin.0, origin.1 + size, 0.0),
        ];
        Polygon::new("square", pts, None)
    }

    #[test]
    fn test_polygon_intersection_overlapping() -> Result<()> {
        // Two overlapping squares
        let poly1 = make_square(2.0, (0.0, 0.0))?;
        let poly2 = make_square(2.0, (1.0, 1.0))?;

        let result = polygon_intersection(&poly1, &poly2)?;

        // Should have one intersection polygon
        assert_eq!(result.len(), 1, "Should have 1 intersection polygon");

        // Intersection area should be 1.0 (1x1 square overlap)
        let area = result.total_area();
        assert!(
            (area - 1.0).abs() < 0.1,
            "Intersection area should be ~1.0, got {}",
            area
        );

        Ok(())
    }

    #[test]
    fn test_polygon_intersection_no_overlap() -> Result<()> {
        // Two non-overlapping squares
        let poly1 = make_square(1.0, (0.0, 0.0))?;
        let poly2 = make_square(1.0, (5.0, 5.0))?;

        let result = polygon_intersection(&poly1, &poly2)?;

        // Should have no intersection
        assert!(result.is_empty());

        Ok(())
    }

    #[test]
    fn test_polygon_intersection_contained() -> Result<()> {
        // Small square inside large square
        let large = make_square(4.0, (0.0, 0.0))?;
        let small = make_square(1.0, (1.0, 1.0))?;

        let result = polygon_intersection(&large, &small)?;

        // Intersection should be the small square
        assert_eq!(result.len(), 1);
        let area = result.total_area();
        assert!(
            (area - 1.0).abs() < 0.1,
            "Intersection area should be ~1.0, got {}",
            area
        );

        Ok(())
    }

    #[test]
    fn test_polygon_difference_no_overlap() -> Result<()> {
        let poly1 = make_square(2.0, (0.0, 0.0))?;
        let poly2 = make_square(1.0, (5.0, 5.0))?;

        let result = polygon_difference(&poly1, &poly2)?;

        // Should return original polygon
        assert_eq!(result.len(), 1);
        let area = result.total_area();
        assert!(
            (area - 4.0).abs() < 0.1,
            "Should return original area ~4.0, got {}",
            area
        );

        Ok(())
    }

    #[test]
    fn test_polygon_difference_creates_hole_with_correct_winding() -> Result<()> {
        let poly1 = make_square(4.0, (0.0, 0.0))?;
        let poly2 = make_square(1.0, (1.0, 1.0))?;

        let result = polygon_difference(&poly1, &poly2)?;
        assert_eq!(result.len(), 1);

        let p = &result.polygons[0];
        assert!(p.has_holes());
        assert_eq!(p.holes().len(), 1);
        assert!(
            (p.area() - 15.0).abs() < 0.1,
            "Area should be ~15.0, got {}",
            p.area()
        );

        // Point inside hole should be outside
        assert!(!p.is_point_inside(Point::new(1.5, 1.5, 0.0), false));
        // Point in solid part should be inside
        assert!(p.is_point_inside(Point::new(0.5, 0.5, 0.0), false));

        Ok(())
    }

    #[test]
    fn test_polygon_overlap_area() -> Result<()> {
        let poly1 = make_square(2.0, (0.0, 0.0))?;
        let poly2 = make_square(2.0, (1.0, 1.0))?;

        let area = polygon_overlap_area(&poly1, &poly2)?;

        assert!(
            (area - 1.0).abs() < 0.1,
            "Overlap area should be ~1.0, got {}",
            area
        );

        Ok(())
    }

    #[test]
    fn test_non_coplanar_error() {
        // Create two non-coplanar polygons
        let pts1 = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
        ];
        let pts2 = vec![
            Point::new(0.0, 0.0, 1.0),
            Point::new(1.0, 0.0, 1.0),
            Point::new(1.0, 0.0, 2.0),
            Point::new(0.0, 0.0, 2.0),
        ];

        let poly1 = Polygon::new("p1", pts1, None).unwrap();
        let poly2 = Polygon::new("p2", pts2, None).unwrap();

        let result = polygon_intersection(&poly1, &poly2);
        assert!(result.is_err());
    }

    #[test]
    fn test_boolean_result_methods() {
        let result = BooleanResult::new(vec![]);
        assert!(result.is_empty());
        assert_eq!(result.len(), 0);
        assert_eq!(result.total_area(), 0.0);
    }
}
