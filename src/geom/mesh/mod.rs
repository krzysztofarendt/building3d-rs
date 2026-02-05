//! Mesh representation and quality analysis.

pub mod quality;
pub mod tetrahedralize;

use crate::Point;
use crate::TriangleIndex;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// A triangle mesh defined by vertices and optional face indices.
///
/// When `faces` is `None` the mesh represents a point cloud only.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Mesh {
    pub vertices: Vec<Point>,
    pub faces: Option<Vec<TriangleIndex>>,
}

impl Mesh {
    /// Creates a new mesh with the given vertices and optional faces.
    pub fn new(vertices: Vec<Point>, faces: Option<Vec<TriangleIndex>>) -> Self {
        Self { vertices, faces }
    }

    /// Returns a reference to the vertices.
    pub fn vertices(&self) -> &[Point] {
        &self.vertices
    }

    /// Returns a reference to the faces if present.
    pub fn faces(&self) -> Option<&[TriangleIndex]> {
        self.faces.as_deref()
    }

    /// Returns the number of vertices.
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }

    /// Returns the number of faces (triangles).
    pub fn face_count(&self) -> usize {
        self.faces.as_ref().map_or(0, |f| f.len())
    }

    /// Returns a new mesh with duplicate vertices merged.
    ///
    /// Vertices are considered identical when they quantize to the same
    /// `(i64, i64, i64)` key at 1e9 scale (≈ 1 nm precision). Face indices
    /// are remapped accordingly. If no faces are present, returns self
    /// unchanged.
    pub fn deduplicate_vertices(self) -> Self {
        let faces = match self.faces {
            Some(ref f) => f,
            None => return self,
        };

        const SCALE: f64 = 1e9;

        let mut key_map: HashMap<(i64, i64, i64), usize> = HashMap::new();
        let mut new_vertices: Vec<Point> = Vec::new();
        let mut old_to_new: Vec<usize> = Vec::with_capacity(self.vertices.len());

        for p in &self.vertices {
            let key = (
                (p.x * SCALE).round() as i64,
                (p.y * SCALE).round() as i64,
                (p.z * SCALE).round() as i64,
            );
            let new_idx = match key_map.get(&key) {
                Some(&idx) => idx,
                None => {
                    let idx = new_vertices.len();
                    new_vertices.push(*p);
                    key_map.insert(key, idx);
                    idx
                }
            };
            old_to_new.push(new_idx);
        }

        let new_faces: Vec<TriangleIndex> = faces
            .iter()
            .map(|t| TriangleIndex(old_to_new[t.0], old_to_new[t.1], old_to_new[t.2]))
            .collect();

        Self {
            vertices: new_vertices,
            faces: Some(new_faces),
        }
    }
}

/// Trait for types that can produce a triangulated [`Mesh`].
///
/// Implemented by every level of the building hierarchy so that drawing and
/// I/O code can work polymorphically.
pub trait HasMesh {
    /// Returns a deep copy of the mesh for this entity.
    fn copy_mesh(&self) -> Mesh;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{HasMesh, Solid};

    #[test]
    fn test_cube_dedup() {
        // A box solid produces 6 quads → 12 triangles → 24 vertices (4 per face),
        // but only 8 unique corners.
        let solid = Solid::from_box(1.0, 1.0, 1.0, None, "cube").unwrap();
        let mesh = solid.copy_mesh();
        assert_eq!(mesh.vertex_count(), 24);
        assert_eq!(mesh.face_count(), 12);

        let deduped = mesh.deduplicate_vertices();
        assert_eq!(deduped.vertex_count(), 8);
        assert_eq!(deduped.face_count(), 12);
    }

    #[test]
    fn test_no_faces_unchanged() {
        let mesh = Mesh {
            vertices: vec![Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0)],
            faces: None,
        };
        let deduped = mesh.deduplicate_vertices();
        assert_eq!(deduped.vertex_count(), 2);
        assert!(deduped.faces.is_none());
    }

    #[test]
    fn test_dedup_index_validity() {
        let solid = Solid::from_box(2.0, 3.0, 4.0, None, "box").unwrap();
        let deduped = solid.copy_mesh().deduplicate_vertices();
        let vc = deduped.vertex_count();
        for tri in deduped.faces.unwrap() {
            assert!(tri.0 < vc);
            assert!(tri.1 < vc);
            assert!(tri.2 < vc);
        }
    }
}
