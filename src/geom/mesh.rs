use crate::Point;
use crate::TriangleIndex;

#[derive(Debug)]
pub struct Mesh {
    pub vertices: Vec<Point>,
    pub faces: Option<Vec<TriangleIndex>>,
}

pub trait HasMesh {
    fn get_mesh(&self) -> Mesh;
}
