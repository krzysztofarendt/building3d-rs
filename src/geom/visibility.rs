//! Visibility analysis.
//!
//! This module provides functions for determining visibility between points
//! in the presence of blocking geometry (polygons).

use crate::geom::ray::Ray;
use crate::{Point, Polygon};

/// Checks if two points are visible to each other.
///
/// Two points are visible if the line segment between them does not
/// intersect any of the blocking polygons.
///
/// # Arguments
/// * `p1` - First point
/// * `p2` - Second point
/// * `blockers` - Polygons that may block visibility
///
/// # Returns
/// `true` if the points can "see" each other, `false` if blocked.
pub fn are_points_visible(p1: Point, p2: Point, blockers: &[&Polygon]) -> bool {
    // Create ray from p1 towards p2
    let ray = match Ray::from_points(p1, p2) {
        Some(r) => r,
        None => return true, // Points are the same, trivially visible
    };

    // Calculate distance from p1 to p2
    let distance = (p2 - p1).length();

    // Check if any blocker intersects the ray between p1 and p2
    for polygon in blockers {
        if let Some((t, _)) = ray.intersect_polygon(polygon) {
            // Check if intersection is between p1 and p2 (not beyond p2)
            // t is the ray parameter, and ray direction is normalized,
            // so t represents distance from p1
            if t < distance - 1e-6 {
                return false; // Blocked
            }
        }
    }

    true
}

/// Computes a visibility matrix between multiple points.
///
/// Returns a 2D vector where `result[i][j]` is `true` if point `i`
/// can see point `j`.
///
/// The matrix is symmetric: `result[i][j] == result[j][i]`.
/// Diagonal entries are always `true` (a point can see itself).
///
/// # Arguments
/// * `points` - List of points to check visibility between
/// * `blockers` - Polygons that may block visibility
///
/// # Returns
/// A 2D vector of booleans representing the visibility matrix.
pub fn visibility_matrix(points: &[Point], blockers: &[&Polygon]) -> Vec<Vec<bool>> {
    let n = points.len();
    let mut matrix = vec![vec![true; n]; n];

    for i in 0..n {
        for j in (i + 1)..n {
            let visible = are_points_visible(points[i], points[j], blockers);
            matrix[i][j] = visible;
            matrix[j][i] = visible;
        }
    }

    matrix
}

/// Finds all points visible from a given point.
///
/// # Arguments
/// * `from` - The point to check visibility from
/// * `targets` - Points to check visibility to
/// * `blockers` - Polygons that may block visibility
///
/// # Returns
/// Indices of visible target points.
pub fn visible_points(from: Point, targets: &[Point], blockers: &[&Polygon]) -> Vec<usize> {
    targets
        .iter()
        .enumerate()
        .filter(|(_, target)| are_points_visible(from, **target, blockers))
        .map(|(idx, _)| idx)
        .collect()
}

/// Counts how many targets are visible from a point.
///
/// # Arguments
/// * `from` - The point to check visibility from
/// * `targets` - Points to check visibility to
/// * `blockers` - Polygons that may block visibility
///
/// # Returns
/// Number of visible targets.
pub fn visibility_count(from: Point, targets: &[Point], blockers: &[&Polygon]) -> usize {
    targets
        .iter()
        .filter(|&&target| are_points_visible(from, target, blockers))
        .count()
}

/// Checks if a point is visible from any of the source points.
///
/// # Arguments
/// * `target` - The point to check visibility to
/// * `sources` - Points to check visibility from
/// * `blockers` - Polygons that may block visibility
///
/// # Returns
/// `true` if visible from at least one source.
pub fn is_visible_from_any(target: Point, sources: &[Point], blockers: &[&Polygon]) -> bool {
    sources
        .iter()
        .any(|&source| are_points_visible(source, target, blockers))
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::Result;

    fn make_blocker_at_z(z: f64) -> Result<Polygon> {
        let pts = vec![
            Point::new(-10.0, -10.0, z),
            Point::new(10.0, -10.0, z),
            Point::new(10.0, 10.0, z),
            Point::new(0.0, 10.0, z),
        ];
        Polygon::new("blocker", pts, None)
    }

    fn make_small_blocker() -> Result<Polygon> {
        // Small square blocker at z=5, centered at (0,0)
        let pts = vec![
            Point::new(-1.0, -1.0, 5.0),
            Point::new(1.0, -1.0, 5.0),
            Point::new(1.0, 1.0, 5.0),
            Point::new(-1.0, 1.0, 5.0),
        ];
        Polygon::new("small_blocker", pts, None)
    }

    #[test]
    fn test_points_visible_no_blockers() {
        let p1 = Point::new(0.0, 0.0, 0.0);
        let p2 = Point::new(10.0, 10.0, 10.0);
        let blockers: Vec<&Polygon> = vec![];

        assert!(are_points_visible(p1, p2, &blockers));
    }

    #[test]
    fn test_points_blocked() -> Result<()> {
        let p1 = Point::new(0.0, 0.0, 0.0);
        let p2 = Point::new(0.0, 0.0, 10.0);
        let blocker = make_blocker_at_z(5.0)?;
        let blockers: Vec<&Polygon> = vec![&blocker];

        assert!(!are_points_visible(p1, p2, &blockers));
        Ok(())
    }

    #[test]
    fn test_points_visible_blocker_not_in_path() -> Result<()> {
        let p1 = Point::new(0.0, 0.0, 0.0);
        let p2 = Point::new(0.0, 0.0, 3.0);
        let blocker = make_blocker_at_z(5.0)?; // Blocker is beyond p2
        let blockers: Vec<&Polygon> = vec![&blocker];

        assert!(are_points_visible(p1, p2, &blockers));
        Ok(())
    }

    #[test]
    fn test_points_visible_around_blocker() -> Result<()> {
        let p1 = Point::new(0.0, 0.0, 0.0);
        let p2 = Point::new(5.0, 5.0, 10.0); // Path goes around the small blocker
        let blocker = make_small_blocker()?;
        let blockers: Vec<&Polygon> = vec![&blocker];

        assert!(are_points_visible(p1, p2, &blockers));
        Ok(())
    }

    #[test]
    fn test_visibility_matrix_no_blockers() {
        let points = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(2.0, 0.0, 0.0),
        ];
        let blockers: Vec<&Polygon> = vec![];

        let matrix = visibility_matrix(&points, &blockers);

        // All points should see all other points
        assert_eq!(matrix.len(), 3);
        for row in &matrix {
            for &visible in row {
                assert!(visible);
            }
        }
    }

    #[test]
    fn test_visibility_matrix_with_blocker() -> Result<()> {
        let points = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(0.0, 0.0, 10.0),
            Point::new(5.0, 5.0, 5.0), // Off to the side
        ];
        let blocker = make_blocker_at_z(5.0)?;
        let blockers: Vec<&Polygon> = vec![&blocker];

        let matrix = visibility_matrix(&points, &blockers);

        // Point 0 and 1 should not see each other (blocker in between)
        assert!(!matrix[0][1]);
        assert!(!matrix[1][0]);

        // Diagonal should be true
        assert!(matrix[0][0]);
        assert!(matrix[1][1]);
        assert!(matrix[2][2]);

        Ok(())
    }

    #[test]
    fn test_visible_points() -> Result<()> {
        let from = Point::new(0.0, 0.0, 0.0);
        let targets = vec![
            Point::new(0.0, 0.0, 10.0), // Blocked
            Point::new(5.0, 5.0, 5.0),   // Visible (around blocker)
            Point::new(0.0, 0.0, 3.0),   // Visible (before blocker)
        ];
        let blocker = make_blocker_at_z(5.0)?;
        let blockers: Vec<&Polygon> = vec![&blocker];

        let visible = visible_points(from, &targets, &blockers);

        // Targets 1 and 2 should be visible
        assert!(!visible.contains(&0)); // Blocked
        assert!(visible.contains(&1));  // Around blocker
        assert!(visible.contains(&2));  // Before blocker

        Ok(())
    }

    #[test]
    fn test_visibility_count() -> Result<()> {
        let from = Point::new(0.0, 0.0, 0.0);
        let targets = vec![
            Point::new(0.0, 0.0, 10.0), // Blocked
            Point::new(5.0, 5.0, 5.0),   // Visible
            Point::new(0.0, 0.0, 3.0),   // Visible
        ];
        let blocker = make_blocker_at_z(5.0)?;
        let blockers: Vec<&Polygon> = vec![&blocker];

        let count = visibility_count(from, &targets, &blockers);
        assert_eq!(count, 2);

        Ok(())
    }

    #[test]
    fn test_is_visible_from_any() -> Result<()> {
        let target = Point::new(0.0, 0.0, 10.0);
        let sources = vec![
            Point::new(0.0, 0.0, 0.0),  // Blocked by blocker
            Point::new(0.0, 0.0, 12.0), // Above target, visible
        ];
        let blocker = make_blocker_at_z(5.0)?;
        let blockers: Vec<&Polygon> = vec![&blocker];

        // Should be visible from source[1] which is above the blocker
        assert!(is_visible_from_any(target, &sources, &blockers));

        Ok(())
    }

    #[test]
    fn test_same_point_visibility() {
        let p = Point::new(1.0, 2.0, 3.0);
        let blockers: Vec<&Polygon> = vec![];

        // A point should be visible to itself
        assert!(are_points_visible(p, p, &blockers));
    }
}
