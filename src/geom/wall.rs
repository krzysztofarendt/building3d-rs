use crate::Point;
use crate::TriangleIndex;
use crate::UID;
use crate::Vector;
use crate::geom;
use crate::geom::polygon::Polygon;
use crate::{HasMesh, Mesh};
use crate::{HasName, SortByName};
use anyhow::{Result, anyhow};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Wall {
    pub name: String,
    pub uid: UID,
    pub parent: Option<UID>,
    polygons: HashMap<String, Polygon>,
}

impl HasName for Wall {
    fn get_name(&self) -> &str {
        &self.name
    }
}

impl HasMesh for Wall {
    fn copy_mesh(&self) -> Mesh {
        let polygons = self.polygons();
        let mut vertices: Vec<Point> = Vec::new();
        let mut triangles: Vec<TriangleIndex> = Vec::new();
        let mut num_vertices = 0;

        for &poly in polygons.iter() {
            let mesh = poly.copy_mesh();
            vertices.extend(mesh.vertices);
            if let Some(faces) = mesh.faces {
                let tri = faces
                    .into_iter()
                    .map(|t| {
                        TriangleIndex(t.0 + num_vertices, t.1 + num_vertices, t.2 + num_vertices)
                    })
                    .collect::<Vec<_>>();
                triangles.extend(tri);
            }
            num_vertices += poly.mesh_ref().vertices.len();
        }

        Mesh {
            vertices,
            faces: Some(triangles),
        }
    }
}

impl Wall {
    pub fn new(name: &str, mut polygons: Vec<Polygon>) -> Self {
        let name = geom::validate_name(name).expect("Invalid wall name");
        let uid = UID::new();
        for p in polygons.iter_mut() {
            p.parent = Some(uid.clone());
        }
        let parent = None;
        let mut map: HashMap<String, Polygon> = HashMap::new();
        for poly in polygons {
            if map.contains_key(&poly.name) {
                panic!("Polygon is already present in Wall::new(): {}", &poly.name);
            }
            map.insert(poly.name.clone(), poly);
        }
        Self {
            name: name.to_string(),
            uid,
            parent,
            polygons: map,
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
        polygon.name = geom::validate_name(&polygon.name)?;
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

    pub(crate) fn repair_parents(&mut self) {
        for poly in self.polygons.values_mut() {
            poly.parent = Some(self.uid.clone());
        }
    }
}
