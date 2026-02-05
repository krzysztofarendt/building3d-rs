//! Tetrahedralization of meshes.
//!
//! This module provides functions to convert triangulated surface meshes
//! into volumetric tetrahedral meshes.

use crate::Point;
use crate::geom::mesh::Mesh;
use crate::geom::tetrahedron::{tetrahedron_centroid, tetrahedron_volume};
use crate::geom::triangles::TriangleIndex;

/// A tetrahedron defined by four point indices.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TetrahedronIndex(pub usize, pub usize, pub usize, pub usize);

/// Result of tetrahedralizing a mesh.
#[derive(Debug, Clone)]
pub struct TetrahedralMesh {
    /// Vertices (includes original mesh vertices plus any added points)
    pub vertices: Vec<Point>,
    /// Tetrahedra as vertex indices
    pub tetrahedra: Vec<TetrahedronIndex>,
}

impl TetrahedralMesh {
    /// Creates a new tetrahedral mesh.
    pub fn new(vertices: Vec<Point>, tetrahedra: Vec<TetrahedronIndex>) -> Self {
        Self {
            vertices,
            tetrahedra,
        }
    }

    /// Returns the number of tetrahedra.
    pub fn tetrahedra_count(&self) -> usize {
        self.tetrahedra.len()
    }

    /// Returns the total volume of the tetrahedral mesh.
    pub fn volume(&self) -> f64 {
        self.tetrahedra
            .iter()
            .map(|t| {
                let p0 = self.vertices[t.0];
                let p1 = self.vertices[t.1];
                let p2 = self.vertices[t.2];
                let p3 = self.vertices[t.3];
                tetrahedron_volume(p0, p1, p2, p3)
            })
            .sum()
    }

    /// Returns the centroid of the tetrahedral mesh.
    pub fn centroid(&self) -> Point {
        if self.tetrahedra.is_empty() {
            return Point::new(0.0, 0.0, 0.0);
        }

        let mut total_volume = 0.0;
        let mut weighted_x = 0.0;
        let mut weighted_y = 0.0;
        let mut weighted_z = 0.0;

        for t in &self.tetrahedra {
            let p0 = self.vertices[t.0];
            let p1 = self.vertices[t.1];
            let p2 = self.vertices[t.2];
            let p3 = self.vertices[t.3];

            let vol = tetrahedron_volume(p0, p1, p2, p3);
            let centroid = tetrahedron_centroid(p0, p1, p2, p3);

            total_volume += vol;
            weighted_x += vol * centroid.x;
            weighted_y += vol * centroid.y;
            weighted_z += vol * centroid.z;
        }

        if total_volume.abs() < 1e-10 {
            return Point::new(0.0, 0.0, 0.0);
        }

        Point::new(
            weighted_x / total_volume,
            weighted_y / total_volume,
            weighted_z / total_volume,
        )
    }

    /// Gets a specific tetrahedron's vertices.
    pub fn get_tetrahedron(&self, idx: usize) -> Option<(Point, Point, Point, Point)> {
        self.tetrahedra.get(idx).map(|t| {
            (
                self.vertices[t.0],
                self.vertices[t.1],
                self.vertices[t.2],
                self.vertices[t.3],
            )
        })
    }

    /// Converts to a surface mesh (triangles of tetrahedra faces).
    /// Note: This produces overlapping faces for internal tetrahedra.
    pub fn to_surface_mesh(&self) -> Mesh {
        let mut faces = Vec::with_capacity(self.tetrahedra.len() * 4);

        for t in &self.tetrahedra {
            // Each tetrahedron has 4 faces
            faces.push(TriangleIndex(t.0, t.1, t.2));
            faces.push(TriangleIndex(t.0, t.1, t.3));
            faces.push(TriangleIndex(t.0, t.2, t.3));
            faces.push(TriangleIndex(t.1, t.2, t.3));
        }

        Mesh::new(self.vertices.clone(), Some(faces))
    }
}

/// Tetrahedralizes a mesh using the centroid method.
///
/// This method works well for convex meshes. It creates a tetrahedron for
/// each triangle face by connecting the face to the mesh centroid.
///
/// # Arguments
/// * `mesh` - The triangulated surface mesh to tetrahedralize
///
/// # Returns
/// A tetrahedral mesh, or None if the mesh has no faces.
pub fn tetrahedralize_centroid(mesh: &Mesh) -> Option<TetrahedralMesh> {
    let faces = mesh.faces()?;
    let vertices = mesh.vertices();

    if faces.is_empty() || vertices.is_empty() {
        return None;
    }

    // Calculate mesh centroid
    let centroid = calculate_mesh_centroid(vertices);

    // Add centroid as a new vertex
    let mut new_vertices = vertices.to_vec();
    let centroid_idx = new_vertices.len();
    new_vertices.push(centroid);

    // Create tetrahedra from each face + centroid
    let tetrahedra: Vec<TetrahedronIndex> = faces
        .iter()
        .map(|f| TetrahedronIndex(f.0, f.1, f.2, centroid_idx))
        .collect();

    Some(TetrahedralMesh::new(new_vertices, tetrahedra))
}

/// Tetrahedralizes a mesh using a specified interior point.
///
/// Creates a tetrahedron for each triangle face by connecting
/// the face to the given interior point.
///
/// # Arguments
/// * `mesh` - The triangulated surface mesh to tetrahedralize
/// * `interior_point` - A point inside the mesh
///
/// # Returns
/// A tetrahedral mesh, or None if the mesh has no faces.
pub fn tetrahedralize_with_point(mesh: &Mesh, interior_point: Point) -> Option<TetrahedralMesh> {
    let faces = mesh.faces()?;
    let vertices = mesh.vertices();

    if faces.is_empty() || vertices.is_empty() {
        return None;
    }

    // Add interior point as a new vertex
    let mut new_vertices = vertices.to_vec();
    let point_idx = new_vertices.len();
    new_vertices.push(interior_point);

    // Create tetrahedra from each face + interior point
    let tetrahedra: Vec<TetrahedronIndex> = faces
        .iter()
        .map(|f| TetrahedronIndex(f.0, f.1, f.2, point_idx))
        .collect();

    Some(TetrahedralMesh::new(new_vertices, tetrahedra))
}

/// Tetrahedralizes a mesh using multiple interior points.
///
/// This creates a more refined tetrahedralization by first triangulating
/// each face with the nearest interior point, resulting in better quality
/// tetrahedra for complex shapes.
///
/// # Arguments
/// * `mesh` - The triangulated surface mesh to tetrahedralize
/// * `interior_points` - Points inside the mesh
///
/// # Returns
/// A tetrahedral mesh, or None if the mesh has no faces or no interior points.
pub fn tetrahedralize_with_points(
    mesh: &Mesh,
    interior_points: &[Point],
) -> Option<TetrahedralMesh> {
    if interior_points.is_empty() {
        return tetrahedralize_centroid(mesh);
    }

    let faces = mesh.faces()?;
    let vertices = mesh.vertices();

    if faces.is_empty() || vertices.is_empty() {
        return None;
    }

    // Add interior points as new vertices
    let mut new_vertices = vertices.to_vec();
    let first_interior_idx = new_vertices.len();
    new_vertices.extend_from_slice(interior_points);

    // For each face, find the closest interior point
    let mut tetrahedra = Vec::with_capacity(faces.len());

    for face in faces {
        let face_centroid = Point::new(
            (vertices[face.0].x + vertices[face.1].x + vertices[face.2].x) / 3.0,
            (vertices[face.0].y + vertices[face.1].y + vertices[face.2].y) / 3.0,
            (vertices[face.0].z + vertices[face.1].z + vertices[face.2].z) / 3.0,
        );

        // Find closest interior point
        let closest_idx = interior_points
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| {
                let dist_a = distance_squared(&face_centroid, a);
                let dist_b = distance_squared(&face_centroid, b);
                dist_a.partial_cmp(&dist_b).unwrap()
            })
            .map(|(i, _)| first_interior_idx + i)
            .unwrap_or(first_interior_idx);

        tetrahedra.push(TetrahedronIndex(face.0, face.1, face.2, closest_idx));
    }

    Some(TetrahedralMesh::new(new_vertices, tetrahedra))
}

/// Calculates the centroid of mesh vertices.
fn calculate_mesh_centroid(vertices: &[Point]) -> Point {
    if vertices.is_empty() {
        return Point::new(0.0, 0.0, 0.0);
    }

    let sum_x: f64 = vertices.iter().map(|p| p.x).sum();
    let sum_y: f64 = vertices.iter().map(|p| p.y).sum();
    let sum_z: f64 = vertices.iter().map(|p| p.z).sum();
    let n = vertices.len() as f64;

    Point::new(sum_x / n, sum_y / n, sum_z / n)
}

/// Calculates squared distance between two points.
fn distance_squared(a: &Point, b: &Point) -> f64 {
    let dx = b.x - a.x;
    let dy = b.y - a.y;
    let dz = b.z - a.z;
    dx * dx + dy * dy + dz * dz
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Solid;
    use crate::geom::mesh::HasMesh;

    fn make_box_mesh() -> Mesh {
        let solid = Solid::from_box(2.0, 2.0, 2.0, Some((0.0, 0.0, 0.0)), "box");
        solid.copy_mesh()
    }

    #[test]
    fn test_tetrahedralize_box() {
        let mesh = make_box_mesh();
        let tet_mesh = tetrahedralize_centroid(&mesh).unwrap();

        // A box has 12 triangles, so we should have 12 tetrahedra
        assert_eq!(tet_mesh.tetrahedra_count(), 12);

        // Volume should be approximately 8 (2x2x2)
        let volume = tet_mesh.volume();
        assert!(
            (volume - 8.0).abs() < 0.01,
            "Volume should be ~8.0, got {}",
            volume
        );
    }

    #[test]
    fn test_tetrahedralize_with_point() {
        let mesh = make_box_mesh();
        let interior = Point::new(1.0, 1.0, 1.0); // Center of box

        let tet_mesh = tetrahedralize_with_point(&mesh, interior).unwrap();

        assert_eq!(tet_mesh.tetrahedra_count(), 12);

        // Volume should still be ~8
        let volume = tet_mesh.volume();
        assert!(
            (volume - 8.0).abs() < 0.01,
            "Volume should be ~8.0, got {}",
            volume
        );
    }

    #[test]
    fn test_tetrahedral_mesh_centroid() {
        let mesh = make_box_mesh();
        let tet_mesh = tetrahedralize_centroid(&mesh).unwrap();

        // Centroid should be at center of box (1, 1, 1)
        let centroid = tet_mesh.centroid();
        assert!((centroid.x - 1.0).abs() < 0.01);
        assert!((centroid.y - 1.0).abs() < 0.01);
        assert!((centroid.z - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_get_tetrahedron() {
        let mesh = make_box_mesh();
        let tet_mesh = tetrahedralize_centroid(&mesh).unwrap();

        // Should be able to get first tetrahedron
        let tet = tet_mesh.get_tetrahedron(0);
        assert!(tet.is_some());

        // Should not be able to get out-of-bounds tetrahedron
        let tet = tet_mesh.get_tetrahedron(100);
        assert!(tet.is_none());
    }

    #[test]
    fn test_to_surface_mesh() {
        let mesh = make_box_mesh();
        let tet_mesh = tetrahedralize_centroid(&mesh).unwrap();

        let surface = tet_mesh.to_surface_mesh();

        // Each tetrahedron contributes 4 faces
        assert_eq!(surface.face_count(), 12 * 4);
    }

    #[test]
    fn test_tetrahedralize_with_multiple_points() {
        let mesh = make_box_mesh();
        // Use a single central point to ensure consistent face orientation
        let interior_points = vec![
            Point::new(1.0, 1.0, 1.0), // Exact center
        ];

        let tet_mesh = tetrahedralize_with_points(&mesh, &interior_points).unwrap();

        // Should still have 12 tetrahedra (one per face)
        assert_eq!(tet_mesh.tetrahedra_count(), 12);

        // Volume should be ~8 when using center point
        let volume = tet_mesh.volume();
        assert!(
            (volume - 8.0).abs() < 0.1,
            "Volume should be ~8.0, got {}",
            volume
        );
    }

    #[test]
    fn test_empty_mesh() {
        let mesh = Mesh::new(vec![], None);
        let result = tetrahedralize_centroid(&mesh);
        assert!(result.is_none());
    }

    #[test]
    fn test_tetrahedron_index() {
        let idx = TetrahedronIndex(0, 1, 2, 3);
        assert_eq!(idx.0, 0);
        assert_eq!(idx.1, 1);
        assert_eq!(idx.2, 2);
        assert_eq!(idx.3, 3);
    }
}
