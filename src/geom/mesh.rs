use crate::Point;
use crate::TriangleIndex;

#[derive(Debug)]
pub struct Mesh {
    pub vertices: Vec<Point>,
    pub triangles: Vec<TriangleIndex>,
}

pub trait GetMesh {
    fn get_mesh(&self) -> Mesh;
}
