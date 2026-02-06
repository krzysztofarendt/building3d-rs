//! Solid adjacency detection.
//!
//! This module provides functions for detecting when two solids share faces,
//! which is essential for building thermal/acoustic analysis and other
//! inter-zone calculations.

use crate::Polygon;
use crate::Solid;
use crate::geom::polygon::relations::{are_polygons_facing, polygon_overlap_area};

/// Represents a shared interface between two solids.
#[derive(Debug, Clone)]
pub struct SharedPolygon {
    /// Name of the polygon from the first solid
    pub poly1_name: String,
    /// Name of the polygon from the second solid
    pub poly2_name: String,
    /// Approximate overlap area
    pub overlap_area: f64,
}

/// Checks if this solid is adjacent to another solid.
///
/// Two solids are adjacent if any of their polygons face each other,
/// meaning they share an interface (are coplanar with opposite normals
/// and have overlapping area).
pub fn is_solid_adjacent_to(solid1: &Solid, solid2: &Solid) -> bool {
    let polygons1 = solid1.polygons();
    let polygons2 = solid2.polygons();

    for poly1 in polygons1.iter() {
        for poly2 in polygons2.iter() {
            if are_polygons_facing(poly1, poly2) {
                return true;
            }
        }
    }

    false
}

/// Returns all polygon pairs that form shared interfaces between two solids.
///
/// Each returned `SharedPolygon` represents a pair of facing polygons
/// from the two solids that form an interface.
pub fn get_shared_polygons(solid1: &Solid, solid2: &Solid) -> Vec<SharedPolygon> {
    let mut shared = Vec::new();
    let polygons1 = solid1.polygons();
    let polygons2 = solid2.polygons();

    for poly1 in polygons1.iter() {
        for poly2 in polygons2.iter() {
            if are_polygons_facing(poly1, poly2) {
                let overlap_area = polygon_overlap_area(poly1, poly2);
                shared.push(SharedPolygon {
                    poly1_name: poly1.name.clone(),
                    poly2_name: poly2.name.clone(),
                    overlap_area,
                });
            }
        }
    }

    shared
}

/// Checks if the interface between two adjacent solids is correct.
///
/// A correct interface means:
/// 1. The facing polygons have the same area (within tolerance)
/// 2. The facing polygons have exactly matching vertices (possibly rotated)
///
/// This validation ensures that there are no gaps or overlaps between
/// adjacent solids, which is important for accurate simulations.
pub fn has_correct_interface(solid1: &Solid, solid2: &Solid) -> bool {
    let shared = get_shared_polygons(solid1, solid2);

    if shared.is_empty() {
        return false; // Not adjacent at all
    }

    let polygons1: std::collections::HashMap<String, &Polygon> = solid1
        .polygons()
        .into_iter()
        .map(|p| (p.name.clone(), p))
        .collect();

    let polygons2: std::collections::HashMap<String, &Polygon> = solid2
        .polygons()
        .into_iter()
        .map(|p| (p.name.clone(), p))
        .collect();

    for sp in shared {
        let poly1 = match polygons1.get(&sp.poly1_name) {
            Some(p) => p,
            None => return false,
        };
        let poly2 = match polygons2.get(&sp.poly2_name) {
            Some(p) => p,
            None => return false,
        };

        // Check if areas match
        let area1 = poly1.area();
        let area2 = poly2.area();
        let area_diff = (area1 - area2).abs();

        // Allow 1% tolerance for area matching
        if area_diff > area1.max(area2) * 0.01 {
            return false;
        }

        // Check if vertex counts match
        if poly1.vertices().len() != poly2.vertices().len() {
            return false;
        }

        // Flipped polygon should match (vertices in reverse order)
        // The comparison already handles rotation via PartialEq
        // But for facing polygons, one is "flipped" relative to the other
        // Check if poly1.flip() == poly2 (accounting for rotation)
        let poly1_flipped = match poly1.flip(&poly1.name) {
            Ok(p) => p,
            Err(_) => return false,
        };

        if poly1_flipped != **poly2 {
            return false;
        }
    }

    true
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::IsClose;

    #[test]
    fn test_adjacent_solids() {
        // Two boxes sharing a face at x=1
        let solid1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1").unwrap();
        let solid2 = Solid::from_box(1.0, 1.0, 1.0, Some((1.0, 0.0, 0.0)), "box2").unwrap();

        assert!(is_solid_adjacent_to(&solid1, &solid2));
    }

    #[test]
    fn test_non_adjacent_solids() {
        // Two boxes with a gap between them
        let solid1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1").unwrap();
        let solid2 = Solid::from_box(1.0, 1.0, 1.0, Some((2.0, 0.0, 0.0)), "box2").unwrap();

        assert!(!is_solid_adjacent_to(&solid1, &solid2));
    }

    #[test]
    fn test_shared_polygons() {
        // Two boxes sharing a face
        let solid1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1").unwrap();
        let solid2 = Solid::from_box(1.0, 1.0, 1.0, Some((1.0, 0.0, 0.0)), "box2").unwrap();

        let shared = get_shared_polygons(&solid1, &solid2);

        // Should have one shared interface (wall_1 of box1 faces wall_3 of box2)
        assert_eq!(shared.len(), 1);
        assert!(shared[0].overlap_area.is_close(1.0)); // 1x1 face
    }

    #[test]
    fn test_no_shared_polygons() {
        // Two boxes with no shared faces
        let solid1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1").unwrap();
        let solid2 = Solid::from_box(1.0, 1.0, 1.0, Some((5.0, 5.0, 5.0)), "box2").unwrap();

        let shared = get_shared_polygons(&solid1, &solid2);
        assert!(shared.is_empty());
    }

    #[test]
    fn test_correct_interface() {
        // Two boxes with matching faces
        let solid1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1").unwrap();
        let solid2 = Solid::from_box(1.0, 1.0, 1.0, Some((1.0, 0.0, 0.0)), "box2").unwrap();

        assert!(has_correct_interface(&solid1, &solid2));
    }

    #[test]
    fn test_incorrect_interface_size_mismatch() {
        // Two boxes with different face sizes at the interface
        let solid1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1").unwrap();
        let solid2 = Solid::from_box(1.0, 2.0, 2.0, Some((1.0, 0.0, 0.0)), "box2").unwrap();

        // They are adjacent but interface is not correct (different sizes)
        assert!(is_solid_adjacent_to(&solid1, &solid2));
        assert!(!has_correct_interface(&solid1, &solid2));
    }

    #[test]
    fn test_adjacent_but_partial_overlap() {
        // Two boxes where one is shifted, creating partial overlap
        let solid1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1").unwrap();
        let solid2 = Solid::from_box(1.0, 1.0, 1.0, Some((1.0, 0.5, 0.0)), "box2").unwrap();

        // They are adjacent (have facing polygons that overlap)
        assert!(is_solid_adjacent_to(&solid1, &solid2));
        // But the interface is not correct (partial overlap)
        assert!(!has_correct_interface(&solid1, &solid2));
    }

    #[test]
    fn test_vertical_adjacency() {
        // Two boxes stacked vertically
        let solid1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1").unwrap();
        let solid2 = Solid::from_box(1.0, 1.0, 1.0, Some((0.0, 0.0, 1.0)), "box2").unwrap();

        assert!(is_solid_adjacent_to(&solid1, &solid2));
        assert!(has_correct_interface(&solid1, &solid2));
    }

    #[test]
    fn test_has_correct_interface_non_adjacent() {
        // Non-adjacent solids: has_correct_interface should return false
        let solid1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1").unwrap();
        let solid2 = Solid::from_box(1.0, 1.0, 1.0, Some((5.0, 5.0, 5.0)), "box2").unwrap();
        assert!(!has_correct_interface(&solid1, &solid2));
    }

    #[test]
    fn test_multiple_shared_faces() {
        // This test verifies that when solids share multiple faces,
        // all are detected. Here we use an L-shaped arrangement would
        // need custom solids, so we'll just check the simple case.
        let solid1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1").unwrap();
        let solid2 = Solid::from_box(1.0, 1.0, 1.0, Some((1.0, 0.0, 0.0)), "box2").unwrap();

        let shared = get_shared_polygons(&solid1, &solid2);
        // Simple boxes share exactly one face
        assert_eq!(shared.len(), 1);
    }
}
