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
    pub fn new(name: &str, mut polygons: Vec<Polygon>) -> Result<Self> {
        let name = geom::validate_name(name)?;
        let uid = UID::new();
        for p in polygons.iter_mut() {
            p.parent = Some(uid.clone());
        }
        let parent = None;
        let mut map: HashMap<String, Polygon> = HashMap::new();
        for poly in polygons {
            if map.contains_key(&poly.name) {
                return Err(anyhow!(
                    "Polygon is already present in Wall::new(): {}",
                    &poly.name
                ));
            }
            map.insert(poly.name.clone(), poly);
        }
        Ok(Self {
            name: name.to_string(),
            uid,
            parent,
            polygons: map,
        })
    }

    pub fn polygons(&self) -> Vec<&Polygon> {
        let mut polygons: Vec<&Polygon> = self.polygons.values().collect();
        polygons.as_mut_slice().sort_by_name();

        polygons
    }

    /// Gets a polygon by name via direct HashMap lookup (O(1)).
    pub fn get_polygon(&self, name: &str) -> Option<&Polygon> {
        self.polygons.get(name)
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::IsClose;

    fn make_triangle(name: &str, z: f64) -> Result<Polygon> {
        Polygon::new(
            name,
            vec![
                Point::new(0., 0., z),
                Point::new(1., 0., z),
                Point::new(0.5, 1., z),
            ],
            None,
        )
    }

    #[test]
    fn test_add_polygon_happy_path() -> Result<()> {
        let poly1 = make_triangle("p1", 0.)?;
        let mut wall = Wall::new("w", vec![poly1])?;
        assert_eq!(wall.polygons().len(), 1);

        let poly2 = make_triangle("p2", 1.)?;
        wall.add_polygon(poly2)?;
        assert_eq!(wall.polygons().len(), 2);
        Ok(())
    }

    #[test]
    fn test_add_polygon_duplicate_error() -> Result<()> {
        let poly = make_triangle("p1", 0.)?;
        let mut wall = Wall::new("w", vec![poly])?;
        let dup = make_triangle("p1", 1.)?;
        assert!(wall.add_polygon(dup).is_err());
        Ok(())
    }

    #[test]
    fn test_replace_polygon_happy_path() -> Result<()> {
        let poly = make_triangle("p1", 0.)?;
        let mut wall = Wall::new("w", vec![poly])?;

        let replacement = make_triangle("p1_new", 2.)?;
        wall.replace_polygon("p1", vec![replacement])?;
        assert_eq!(wall.polygons().len(), 1);
        assert!(wall.get_polygon("p1_new").is_some());
        assert!(wall.get_polygon("p1").is_none());
        Ok(())
    }

    #[test]
    fn test_replace_polygon_missing_error() -> Result<()> {
        let poly = make_triangle("p1", 0.)?;
        let mut wall = Wall::new("w", vec![poly])?;
        let replacement = make_triangle("p2", 1.)?;
        assert!(wall.replace_polygon("nonexistent", vec![replacement]).is_err());
        Ok(())
    }

    #[test]
    fn test_get_polygon() -> Result<()> {
        let poly = make_triangle("p1", 0.)?;
        let wall = Wall::new("w", vec![poly])?;
        assert!(wall.get_polygon("p1").is_some());
        assert!(wall.get_polygon("missing").is_none());
        Ok(())
    }

    #[test]
    fn test_new_duplicate_polygon_error() -> Result<()> {
        let p1 = make_triangle("dup", 0.)?;
        let p2 = make_triangle("dup", 1.)?;
        let result = Wall::new("w", vec![p1, p2]);
        assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn test_copy_mesh() -> Result<()> {
        let poly = make_triangle("p1", 0.)?;
        let wall = Wall::new("w", vec![poly])?;
        let mesh = wall.copy_mesh();
        assert_eq!(mesh.vertices.len(), 3);
        assert!(mesh.faces.is_some());
        assert_eq!(mesh.faces.as_ref().unwrap().len(), 1);
        Ok(())
    }

    #[test]
    fn test_rotate_and_translate() -> Result<()> {
        let poly = Polygon::new(
            "p1",
            vec![
                Point::new(0., 0., 0.),
                Point::new(1., 0., 0.),
                Point::new(1., 1., 0.),
                Point::new(0., 1., 0.),
            ],
            None,
        )?;
        let mut wall = Wall::new("w", vec![poly])?;

        // Translate
        wall.translate(&Vector::new(10., 0., 0.));
        let p = &wall.polygons()[0];
        assert!(p.vertices()[0].x.is_close(10.));

        // Rotate 90 deg around z axis
        wall.rotate(std::f64::consts::PI / 2., &Vector::new(0., 0., 1.));
        // After rotation, points should still form a valid polygon
        assert_eq!(wall.polygons()[0].vertices().len(), 4);
        Ok(())
    }
}
