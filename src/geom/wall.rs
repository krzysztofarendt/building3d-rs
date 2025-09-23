use crate::geom::polygon::Polygon;
use crate::random_id;
use anyhow::{Result, anyhow};
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct Wall {
    pub name: String,
    pub uid: String,
    pub parent: Option<String>,
    polygons: HashMap<String, Polygon>,
}

impl Wall {
    pub fn new(name: String, mut polygons: Vec<Polygon>) -> Self {
        let uid = random_id();
        for p in polygons.iter_mut() {
            p.parent = Some(uid.clone());
        }
        let parent = None;
        let polygons: HashMap<String, Polygon> =
            polygons.into_iter().map(|x| (x.name.clone(), x)).collect();
        Self {
            name,
            uid,
            parent,
            polygons,
        }
    }

    pub fn polygons(&self) -> Vec<&Polygon> {
        self.polygons.values().collect()
    }

    pub fn add_polygon(&mut self, mut polygon: Polygon) -> Result<()> {
        if self.polygons.contains_key(&polygon.name) {
            return Err(anyhow!("Polygon is already present: {}", &polygon.name));
        }
        polygon.parent = Some(self.uid.clone());
        self.polygons.insert(polygon.name.clone(), polygon);

        Ok(())
    }
}
