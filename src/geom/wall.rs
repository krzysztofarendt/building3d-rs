use crate::Point;
use crate::TriangleIndex;
use crate::Vector;
use crate::geom::polygon::Polygon;
use crate::random_id;
use crate::sortbyname::{HasName, SortByName};
use crate::{GetMesh, Mesh};
use anyhow::{Result, anyhow};
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct Wall {
    pub name: String,
    pub uid: String,
    pub parent: Option<String>,
    polygons: HashMap<String, Polygon>,
}

impl HasName for Wall {
    fn name(&self) -> &str {
        &self.name
    }
}

impl GetMesh for Wall {
    fn get_mesh(&self) -> Mesh {
        let polygons = self.polygons();
        let vertices: Vec<Point> = polygons.iter().flat_map(|&p| p.pts.clone()).collect();
        let mut triangles: Vec<TriangleIndex> = Vec::new();
        let mut num_vertices = 0;

        for &poly in polygons.iter() {
            let mut tri: Vec<TriangleIndex> = poly.tri.clone();
            tri = tri
                .into_iter()
                .map(|t| TriangleIndex(t.0 + num_vertices, t.1 + num_vertices, t.2 + num_vertices))
                .collect();
            triangles.extend(tri.into_iter());
            num_vertices += poly.pts.len();
        }

        Mesh {
            vertices,
            triangles,
        }
    }
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
        let mut polygons: Vec<&Polygon> = self.polygons.values().collect();
        polygons.as_mut_slice().sort_by_name();

        polygons
    }

    pub fn rotate(&mut self, angle: f64, rot_vec: &Vector) {
        let mut polygons: Vec<&mut Polygon> = self.polygons.values_mut().collect();
        for poly in polygons.iter_mut() {
            poly.rotate(angle, rot_vec);
        }
    }

    pub fn translate(&mut self, vec: &Vector) {
        let mut polygons: Vec<&mut Polygon> = self.polygons.values_mut().collect();
        for poly in polygons.iter_mut() {
            poly.translate(vec);
        }
    }

    pub fn add_polygon(&mut self, mut polygon: Polygon) -> Result<()> {
        if self.polygons.contains_key(&polygon.name) {
            return Err(anyhow!("Polygon is already present: {}", &polygon.name));
        }
        polygon.parent = Some(self.uid.clone());
        self.polygons.insert(polygon.name.clone(), polygon);

        Ok(())
    }

    pub fn replace_polygon(&mut self, old_name: &str, new_poly: Vec<Polygon>) -> Result<()> {
        let removed: Option<_> = self.polygons.remove(old_name);
        if removed.is_none() {
            return Err(anyhow!(
                "No such polygon ({}) in this wall ({})",
                old_name,
                &self.name
            ));
        }
        for pl in new_poly.into_iter() {
            self.add_polygon(pl)?; // TODO: Add some check, e.g. compare surface area
        }

        Ok(())
    }
}
