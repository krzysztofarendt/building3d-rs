//! Mesh representation and quality analysis.

pub mod quality;
pub mod tetrahedralize;

use crate::Point;
use crate::TriangleIndex;
use serde::{Deserialize, Serialize};

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
}

pub trait HasMesh {
    fn copy_mesh(&self) -> Mesh;
}
