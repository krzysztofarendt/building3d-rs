//! Delaunay tetrahedralization via the Bowyer-Watson incremental insertion algorithm.
//!
//! This module provides a proper Delaunay tetrahedralization that works for both
//! convex and concave meshes. For concave meshes, a post-processing culling step
//! removes tetrahedra whose centroids lie outside the original surface.

use std::collections::HashMap;

use crate::Point;
use crate::Vector;
use crate::geom::bboxes::bounding_box;
use crate::geom::point::check::are_points_coplanar;
use crate::geom::ray::moller_trumbore;
use crate::geom::tetrahedron::{circumsphere, tetrahedron_centroid, tetrahedron_volume};
use crate::geom::triangles::TriangleIndex;

/// Internal tetrahedron representation with cached circumsphere data.
struct BwTet {
    v: [usize; 4],
    center: Point,
    radius_sq: f64,
}

/// A face key with sorted vertex indices for hashing.
#[derive(Hash, Eq, PartialEq)]
struct FaceKey([usize; 3]);

impl FaceKey {
    fn new(a: usize, b: usize, c: usize) -> Self {
        let mut arr = [a, b, c];
        arr.sort();
        FaceKey(arr)
    }
}

/// Creates a super-tetrahedron enclosing all given points.
///
/// Returns four new points that form a tetrahedron large enough to contain
/// all input points within its circumsphere.
fn super_tetrahedron(points: &[Point]) -> [Point; 4] {
    let (pmin, pmax) = bounding_box(points);

    let cx = (pmin.x + pmax.x) * 0.5;
    let cy = (pmin.y + pmax.y) * 0.5;
    let cz = (pmin.z + pmax.z) * 0.5;

    let dx = (pmax.x - pmin.x).max(1e-6);
    let dy = (pmax.y - pmin.y).max(1e-6);
    let dz = (pmax.z - pmin.z).max(1e-6);

    // Half-diagonal of the bounding box
    let half_diag = (dx * dx + dy * dy + dz * dz).sqrt() * 0.5;
    // Scale factor to ensure all points are well inside the super-tet
    let scale = 4.0 * half_diag;

    // Regular tetrahedron vertices centered at (cx, cy, cz)
    // Using a regular tetrahedron inscribed in a sphere of radius `scale`
    [
        Point::new(cx + scale * 1.0, cy + scale * 0.0, cz - scale * 0.707),
        Point::new(cx - scale * 1.0, cy + scale * 0.0, cz - scale * 0.707),
        Point::new(cx + scale * 0.0, cy + scale * 1.0, cz + scale * 0.707),
        Point::new(cx + scale * 0.0, cy - scale * 1.0, cz + scale * 0.707),
    ]
}

/// The four faces of a tetrahedron as triples of vertex indices.
fn tet_faces(v: &[usize; 4]) -> [[usize; 3]; 4] {
    [
        [v[0], v[1], v[2]],
        [v[0], v[1], v[3]],
        [v[0], v[2], v[3]],
        [v[1], v[2], v[3]],
    ]
}

/// Bowyer-Watson incremental Delaunay tetrahedralization.
///
/// Returns `None` if fewer than 4 points are provided or all points are coplanar.
/// The returned tetrahedra reference indices into the original `points` slice.
pub(crate) fn bowyer_watson(points: &[Point]) -> Option<Vec<[usize; 4]>> {
    let n = points.len();
    if n < 4 {
        return None;
    }
    if are_points_coplanar(points) {
        return None;
    }

    // Build combined vertex array: original points + super-tet vertices
    let super_pts = super_tetrahedron(points);
    let mut all_points: Vec<Point> = points.to_vec();
    all_points.extend_from_slice(&super_pts);

    // Super-tet vertex indices
    let si = [n, n + 1, n + 2, n + 3];

    // Initialize with super-tetrahedron
    let (center, radius_sq) = circumsphere(
        all_points[si[0]],
        all_points[si[1]],
        all_points[si[2]],
        all_points[si[3]],
    )?;

    let mut tets: Vec<BwTet> = vec![BwTet {
        v: si,
        center,
        radius_sq,
    }];

    // Insert points one at a time
    for i in 0..n {
        let pt = all_points[i];

        // Find bad tets: those whose circumsphere contains the new point
        let mut bad_indices: Vec<usize> = Vec::new();
        for (ti, tet) in tets.iter().enumerate() {
            let dx = tet.center.x - pt.x;
            let dy = tet.center.y - pt.y;
            let dz = tet.center.z - pt.z;
            let dist_sq = dx * dx + dy * dy + dz * dz;
            if dist_sq < tet.radius_sq + 1e-10 {
                bad_indices.push(ti);
            }
        }

        if bad_indices.is_empty() {
            continue;
        }

        // Find cavity boundary: faces shared by exactly one bad tet
        let mut face_count: HashMap<FaceKey, (usize, [usize; 3])> = HashMap::new();
        for &bi in &bad_indices {
            let faces = tet_faces(&tets[bi].v);
            for face in &faces {
                let key = FaceKey::new(face[0], face[1], face[2]);
                face_count
                    .entry(key)
                    .and_modify(|(count, _)| *count += 1)
                    .or_insert((1, *face));
            }
        }

        // Boundary faces are those with count == 1
        let boundary_faces: Vec<[usize; 3]> = face_count
            .into_values()
            .filter(|(count, _)| *count == 1)
            .map(|(_, face)| face)
            .collect();

        // Remove bad tets (in reverse order to preserve indices)
        bad_indices.sort_unstable();
        for &bi in bad_indices.iter().rev() {
            tets.swap_remove(bi);
        }

        // Create new tets from boundary faces + new point
        for face in &boundary_faces {
            let v = [face[0], face[1], face[2], i];
            if let Some((center, radius_sq)) = circumsphere(
                all_points[v[0]],
                all_points[v[1]],
                all_points[v[2]],
                all_points[v[3]],
            ) {
                tets.push(BwTet {
                    v,
                    center,
                    radius_sq,
                });
            }
        }
    }

    // Remove tets referencing super-tet vertices
    tets.retain(|t| t.v.iter().all(|&vi| vi < n));

    // Filter out near-zero volume tets
    let result: Vec<[usize; 4]> = tets
        .into_iter()
        .filter(|t| {
            let vol = tetrahedron_volume(
                all_points[t.v[0]],
                all_points[t.v[1]],
                all_points[t.v[2]],
                all_points[t.v[3]],
            );
            vol > 1e-15
        })
        .map(|t| t.v)
        .collect();

    if result.is_empty() {
        None
    } else {
        Some(result)
    }
}

/// Cull tetrahedra whose centroids lie outside the original surface mesh.
///
/// Uses ray-casting with majority voting to determine inside/outside status.
pub(crate) fn cull_exterior_tets(
    vertices: &[Point],
    tets: Vec<[usize; 4]>,
    surface_faces: &[TriangleIndex],
) -> Vec<[usize; 4]> {
    if surface_faces.is_empty() {
        return tets;
    }

    // Use irrational-ish directions to avoid hitting edges/vertices of axis-aligned geometry
    let ray_dirs = [
        Vector::new(1.0, 0.3141592653, 0.2718281828),
        Vector::new(0.1732050808, 1.0, 0.1414213562),
        Vector::new(0.2236067977, 0.1618033988, 1.0),
        Vector::new(0.7071167811, 0.5773502691, 0.4142135623),
        Vector::new(-1.0, 0.2718281828, -0.3141592653),
    ];

    tets.into_iter()
        .filter(|t| {
            let centroid = tetrahedron_centroid(
                vertices[t[0]],
                vertices[t[1]],
                vertices[t[2]],
                vertices[t[3]],
            );

            let mut inside_votes = 0;
            for dir in &ray_dirs {
                let mut count = 0usize;
                for face in surface_faces {
                    let p0 = vertices[face.0];
                    let p1 = vertices[face.1];
                    let p2 = vertices[face.2];
                    // Use small negative tolerance to avoid double-counting edge hits
                    if moller_trumbore(centroid, *dir, p0, p1, p2, -1e-6).is_some() {
                        count += 1;
                    }
                }
                if count % 2 == 1 {
                    inside_votes += 1;
                }
            }

            // Majority vote: at least 3 of 5 rays say inside
            inside_votes >= 3
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bowyer_watson_single_tet() {
        let points = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
            Point::new(0.0, 0.0, 1.0),
        ];
        let tets = bowyer_watson(&points).unwrap();
        assert_eq!(tets.len(), 1);
    }

    #[test]
    fn test_bowyer_watson_cube() {
        let points = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
            Point::new(0.0, 0.0, 1.0),
            Point::new(1.0, 0.0, 1.0),
            Point::new(1.0, 1.0, 1.0),
            Point::new(0.0, 1.0, 1.0),
        ];
        let tets = bowyer_watson(&points).unwrap();

        // A cube should be decomposed into 5 or 6 tets
        assert!(
            tets.len() >= 5 && tets.len() <= 12,
            "Expected 5-12 tets, got {}",
            tets.len()
        );

        // Total volume should be 1.0
        let total_vol: f64 = tets
            .iter()
            .map(|t| tetrahedron_volume(points[t[0]], points[t[1]], points[t[2]], points[t[3]]))
            .sum();
        assert!(
            (total_vol - 1.0).abs() < 1e-10,
            "Volume should be 1.0, got {}",
            total_vol
        );
    }

    #[test]
    fn test_bowyer_watson_delaunay_property() {
        let points = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
            Point::new(0.0, 0.0, 1.0),
            Point::new(1.0, 1.0, 1.0),
        ];
        let tets = bowyer_watson(&points).unwrap();

        // No point should lie strictly inside any tet's circumsphere
        for t in &tets {
            let p0 = points[t[0]];
            let p1 = points[t[1]];
            let p2 = points[t[2]];
            let p3 = points[t[3]];

            if let Some((center, radius_sq)) = circumsphere(p0, p1, p2, p3) {
                for (pi, pt) in points.iter().enumerate() {
                    if pi == t[0] || pi == t[1] || pi == t[2] || pi == t[3] {
                        continue;
                    }
                    let dx = center.x - pt.x;
                    let dy = center.y - pt.y;
                    let dz = center.z - pt.z;
                    let dist_sq = dx * dx + dy * dy + dz * dz;
                    assert!(
                        dist_sq >= radius_sq - 1e-8,
                        "Point {} is inside circumsphere of tet {:?}",
                        pi,
                        t
                    );
                }
            }
        }
    }

    #[test]
    fn test_bowyer_watson_coplanar_returns_none() {
        let points = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
        ];
        assert!(bowyer_watson(&points).is_none());
    }

    #[test]
    fn test_cull_exterior_convex() {
        // For a convex box, no tets should be culled
        let points = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
            Point::new(0.0, 0.0, 1.0),
            Point::new(1.0, 0.0, 1.0),
            Point::new(1.0, 1.0, 1.0),
            Point::new(0.0, 1.0, 1.0),
        ];

        // Surface triangles of the cube (2 triangles per face, 6 faces)
        let faces = vec![
            // Bottom (z=0)
            TriangleIndex(0, 1, 2),
            TriangleIndex(0, 2, 3),
            // Top (z=1)
            TriangleIndex(4, 6, 5),
            TriangleIndex(4, 7, 6),
            // Front (y=0)
            TriangleIndex(0, 5, 1),
            TriangleIndex(0, 4, 5),
            // Back (y=1)
            TriangleIndex(2, 7, 3),
            TriangleIndex(2, 6, 7),
            // Left (x=0)
            TriangleIndex(0, 7, 4),
            TriangleIndex(0, 3, 7),
            // Right (x=1)
            TriangleIndex(1, 6, 2),
            TriangleIndex(1, 5, 6),
        ];

        let tets = bowyer_watson(&points).unwrap();
        let n_before = tets.len();
        let culled = cull_exterior_tets(&points, tets, &faces);

        assert_eq!(
            culled.len(),
            n_before,
            "No tets should be culled for a convex box"
        );
    }
}
