//! Point-in-solid containment testing using ray casting.
//!
//! This module implements the ray casting algorithm to determine if a point
//! lies inside a closed solid. The algorithm casts a ray from the test point
//! and counts the number of polygon intersections - an odd count means inside.

use crate::Point;
use crate::geom::EPS;
use crate::geom::bboxes::bounding_box;
use crate::geom::segment::segment_crosses_polygon;
use crate::geom::solid::Solid;

/// Checks if a point lies inside a solid using the ray casting algorithm.
///
/// The algorithm casts a ray from the test point in a random direction and
/// counts how many times it crosses the solid's surface. An odd number of
/// crossings means the point is inside.
///
/// # Arguments
/// * `solid` - The solid to test against
/// * `ptest` - The point to test
///
/// # Returns
/// `true` if the point is inside the solid, `false` otherwise.
/// Points exactly on the boundary may return either value due to numerical precision.
pub fn is_point_inside_solid(solid: &Solid, ptest: Point) -> bool {
    // Quick rejection: check bounding box first
    let mesh = solid.copy_mesh();
    let (bbox_min, bbox_max) = bounding_box(&mesh.vertices);

    // Add small tolerance to bbox check
    if ptest.x < bbox_min.x - EPS
        || ptest.x > bbox_max.x + EPS
        || ptest.y < bbox_min.y - EPS
        || ptest.y > bbox_max.y + EPS
        || ptest.z < bbox_min.z - EPS
        || ptest.z > bbox_max.z + EPS
    {
        return false;
    }

    // Cast ray in multiple directions to handle edge cases
    // If results are inconsistent, the point is likely on the boundary
    let directions = [
        (1.0, 0.0, 0.0), // +X
        (0.0, 1.0, 0.0), // +Y
        (0.0, 0.0, 1.0), // +Z
        (1.0, 1.0, 1.0), // Diagonal
    ];

    let mut inside_count = 0;
    let mut outside_count = 0;

    for (dx, dy, dz) in directions {
        let result = cast_ray(solid, ptest, dx, dy, dz, &bbox_min, &bbox_max);
        if result {
            inside_count += 1;
        } else {
            outside_count += 1;
        }
    }

    // Use majority vote
    inside_count > outside_count
}

/// Casts a ray from the test point and counts polygon crossings.
fn cast_ray(
    solid: &Solid,
    ptest: Point,
    dx: f64,
    dy: f64,
    dz: f64,
    bbox_min: &Point,
    bbox_max: &Point,
) -> bool {
    // Normalize direction
    let len = (dx * dx + dy * dy + dz * dz).sqrt();
    let dx = dx / len;
    let dy = dy / len;
    let dz = dz / len;

    // Calculate ray endpoint far enough to exit the solid.
    // Use the full bbox diagonal length so this remains valid even if `ptest`
    // is near one corner and the ray points toward the opposite corner.
    let bbox_diag = ((bbox_max.x - bbox_min.x).powi(2)
        + (bbox_max.y - bbox_min.y).powi(2)
        + (bbox_max.z - bbox_min.z).powi(2))
    .sqrt();
    let ray_length = bbox_diag * 2.0 + 10.0;

    let ray_end = Point::new(
        ptest.x + dx * ray_length,
        ptest.y + dy * ray_length,
        ptest.z + dz * ray_length,
    );

    // Count crossings
    let mut crossing_count = 0;

    for polygon in solid.polygons() {
        if let Some(_intersection) = segment_crosses_polygon(ptest, ray_end, polygon) {
            crossing_count += 1;
        }
    }

    // Odd number of crossings = inside
    crossing_count % 2 == 1
}

/// Checks if a point lies on the boundary (surface) of a solid.
///
/// A point is on the boundary if it lies on any of the solid's polygon surfaces.
pub fn is_point_at_boundary(solid: &Solid, ptest: Point) -> bool {
    for polygon in solid.polygons() {
        if polygon.is_point_inside(ptest, true) {
            return true;
        }
    }
    false
}

/// Checks if a point lies strictly inside a solid (not on boundary).
pub fn is_point_strictly_inside(solid: &Solid, ptest: Point) -> bool {
    if is_point_at_boundary(solid, ptest) {
        return false;
    }
    is_point_inside_solid(solid, ptest)
}

use crate::HasMesh;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::bboxes::are_bboxes_overlapping;

    #[test]
    fn test_point_inside_box() {
        let solid = Solid::from_box(2.0, 2.0, 2.0, Some((0.0, 0.0, 0.0)), "box").unwrap();

        // Center point - should be inside
        let center = Point::new(1.0, 1.0, 1.0);
        assert!(is_point_inside_solid(&solid, center));

        // Off-center but still inside
        let inside = Point::new(0.5, 0.5, 0.5);
        assert!(is_point_inside_solid(&solid, inside));

        let inside = Point::new(1.5, 1.5, 1.5);
        assert!(is_point_inside_solid(&solid, inside));
    }

    #[test]
    fn test_point_outside_box() {
        let solid = Solid::from_box(2.0, 2.0, 2.0, Some((0.0, 0.0, 0.0)), "box").unwrap();

        // Outside on each axis
        let outside = Point::new(-1.0, 1.0, 1.0);
        assert!(!is_point_inside_solid(&solid, outside));

        let outside = Point::new(3.0, 1.0, 1.0);
        assert!(!is_point_inside_solid(&solid, outside));

        let outside = Point::new(1.0, -1.0, 1.0);
        assert!(!is_point_inside_solid(&solid, outside));

        let outside = Point::new(1.0, 1.0, 3.0);
        assert!(!is_point_inside_solid(&solid, outside));

        // Far outside
        let outside = Point::new(10.0, 10.0, 10.0);
        assert!(!is_point_inside_solid(&solid, outside));
    }

    #[test]
    fn test_point_on_face() {
        let solid = Solid::from_box(2.0, 2.0, 2.0, Some((0.0, 0.0, 0.0)), "box").unwrap();

        // Point on face - should be at boundary
        let on_face = Point::new(1.0, 1.0, 0.0); // On floor
        assert!(is_point_at_boundary(&solid, on_face));

        let on_face = Point::new(1.0, 1.0, 2.0); // On ceiling
        assert!(is_point_at_boundary(&solid, on_face));

        let on_face = Point::new(0.0, 1.0, 1.0); // On side wall
        assert!(is_point_at_boundary(&solid, on_face));
    }

    #[test]
    fn test_point_strictly_inside() {
        let solid = Solid::from_box(2.0, 2.0, 2.0, Some((0.0, 0.0, 0.0)), "box").unwrap();

        // Center - strictly inside
        let center = Point::new(1.0, 1.0, 1.0);
        assert!(is_point_strictly_inside(&solid, center));

        // On face - not strictly inside
        let on_face = Point::new(1.0, 1.0, 0.0);
        assert!(!is_point_strictly_inside(&solid, on_face));
    }

    #[test]
    fn test_point_at_corner() {
        let solid = Solid::from_box(2.0, 2.0, 2.0, Some((0.0, 0.0, 0.0)), "box").unwrap();

        // Corner vertex - should be at boundary
        let corner = Point::new(0.0, 0.0, 0.0);
        assert!(is_point_at_boundary(&solid, corner));

        let corner = Point::new(2.0, 2.0, 2.0);
        assert!(is_point_at_boundary(&solid, corner));
    }

    #[test]
    fn test_point_on_edge() {
        let solid = Solid::from_box(2.0, 2.0, 2.0, Some((0.0, 0.0, 0.0)), "box").unwrap();

        // Point on edge - should be at boundary
        let on_edge = Point::new(1.0, 0.0, 0.0); // Bottom front edge
        assert!(is_point_at_boundary(&solid, on_edge));
    }

    #[test]
    fn test_translated_solid() {
        let mut solid = Solid::from_box(1.0, 1.0, 1.0, None, "box").unwrap();
        solid.translate(&crate::Vector::new(5.0, 5.0, 5.0));

        // Original center (0.5, 0.5, 0.5) is now outside
        let old_center = Point::new(0.5, 0.5, 0.5);
        assert!(!is_point_inside_solid(&solid, old_center));

        // New center (5.5, 5.5, 5.5) is inside
        let new_center = Point::new(5.5, 5.5, 5.5);
        assert!(is_point_inside_solid(&solid, new_center));
    }

    #[test]
    fn test_bboxes_overlapping() {
        let min1 = Point::new(0.0, 0.0, 0.0);
        let max1 = Point::new(2.0, 2.0, 2.0);

        // Overlapping box
        let min2 = Point::new(1.0, 1.0, 1.0);
        let max2 = Point::new(3.0, 3.0, 3.0);
        assert!(are_bboxes_overlapping(min1, max1, min2, max2));

        // Touching box
        let min3 = Point::new(2.0, 0.0, 0.0);
        let max3 = Point::new(4.0, 2.0, 2.0);
        assert!(are_bboxes_overlapping(min1, max1, min3, max3));

        // Separated box
        let min4 = Point::new(5.0, 5.0, 5.0);
        let max4 = Point::new(6.0, 6.0, 6.0);
        assert!(!are_bboxes_overlapping(min1, max1, min4, max4));
    }
}
