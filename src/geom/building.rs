use crate::random_id;
use crate::Point;
use crate::Polygon;
use crate::Solid;
use crate::TriangleIndex;
use crate::Vector;
use crate::Wall;
use crate::{HasMesh, Mesh};
use crate::{HasName, SortByName};
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct Building {
    pub name: String,
    pub uid: String,
    pub parent: Option<String>,
    solids: HashMap<String, Solid>,
}

impl HasName for Building {
    fn get_name(&self) -> &str {
        &self.name
    }
}

impl HasMesh for Building {
    fn copy_mesh(&self) -> Mesh {
        let polygons = self.polygons();
        let mut vertices: Vec<Point> = Vec::new();
        let mut triangles: Vec<TriangleIndex> = Vec::new();
        let mut num_vertices = 0;

        for &poly in polygons.iter() {
            let mesh = poly.copy_mesh();
            vertices.extend(mesh.vertices);
            let mut tri: Vec<TriangleIndex> = mesh.faces.unwrap();
            tri = tri
                .into_iter()
                .map(|t| TriangleIndex(t.0 + num_vertices, t.1 + num_vertices, t.2 + num_vertices))
                .collect();
            triangles.extend(tri.into_iter());
            num_vertices += poly.mesh_ref().vertices.len();
        }

        Mesh {
            vertices,
            faces: Some(triangles),
        }
    }
}

impl Building {
    pub fn new(name: &str, mut solids: Vec<Solid>) -> Self {
        let uid = random_id();
        for s in solids.iter_mut() {
            s.parent = Some(uid.clone());
        }
        let parent = None;
        let solids: HashMap<String, Solid> =
            solids.into_iter().map(|x| (x.name.clone(), x)).collect();

        Self {
            name: name.to_string(),
            uid,
            parent,
            solids,
        }
    }

    /// Return solids sorted by name
    pub fn solids(&self) -> Vec<&Solid> {
        let mut solids: Vec<&Solid> = self.solids.values().collect();
        solids.as_mut_slice().sort_by_name();

        solids
    }

    /// Return walls sorted by name
    pub fn walls(&self) -> Vec<&Wall> {
        let solids = self.solids();
        let walls: Vec<&Wall> = solids.iter().flat_map(|s| s.walls()).collect();
        // NOTE: Walls are already sorted

        walls
    }

    pub fn polygons(&self) -> Vec<&Polygon> {
        let walls = self.walls();
        let polygons: Vec<&Polygon> = walls.iter().flat_map(|w| w.polygons()).collect();
        // NOTE: Polygons are already sorted

        polygons
    }

    pub fn rotate(&mut self, angle: f64, rot_vec: &Vector) {
        let mut solids: Vec<&mut Solid> = self.solids.values_mut().collect();
        for sld in solids.iter_mut() {
            sld.rotate(angle, rot_vec);
        }
    }

    pub fn translate(&mut self, vec: &Vector) {
        let mut solids: Vec<&mut Solid> = self.solids.values_mut().collect();
        for sld in solids.iter_mut() {
            sld.translate(vec);
        }
    }

    pub fn volume(&self) -> f64 {
        let mut v = 0.;
        for sld in self.solids() {
            v += sld.volume();
        }
        v
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_volume() {
        let s0 = Solid::from_box(1., 1., 1., None, None);
        let s1 = Solid::from_box(1., 2., 3., None, None);
        let bdg = Building::new("building", vec![s0, s1]);
        let expected_vol = 1. * 2. * 3. + 1.;
        assert!((bdg.volume() - expected_vol).abs() < 1e-4);
    }
}
