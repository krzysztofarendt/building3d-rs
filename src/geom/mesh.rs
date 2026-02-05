use crate::Point;
use crate::TriangleIndex;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Mesh {
    pub vertices: Vec<Point>,
    pub faces: Option<Vec<TriangleIndex>>,
}

pub trait HasMesh {
    fn copy_mesh(&self) -> Mesh;
}
