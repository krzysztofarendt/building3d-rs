use crate::Point;
use crate::TriangleIndex;

#[derive(Debug, Clone)]
pub struct Mesh {
    pub vertices: Vec<Point>,
    pub faces: Option<Vec<TriangleIndex>>,
}

pub trait HasMesh {
    fn copy_mesh(&self) -> Mesh;
}
