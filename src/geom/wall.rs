use crate::geom::polygon::Polygon;
use crate::random_id;
// use anyhow::Result;

#[derive(Debug, Clone)]
pub struct Wall {
    pub name: String,
    pub polygons: Vec<Polygon>,
    pub uid: String,
}

impl Wall {
    pub fn new(name: String, mut polygons: Vec<Polygon>) -> Self {
        let uid = random_id();
        for p in polygons.iter_mut() {
            p.parent = Some(uid.clone());
        }
        Self {
            name,
            polygons,
            uid,
        }
    }

    pub fn polygons(&self) -> &[Polygon] {
        &self.polygons
    }
}
