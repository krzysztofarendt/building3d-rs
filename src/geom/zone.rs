//! Zone container for grouping solids.
//!
//! A Zone is an intermediate level in the building hierarchy:
//! Building → Zone → Solid → Wall → Polygon

use crate::Point;
use crate::Polygon;
use crate::Solid;
use crate::TriangleIndex;
use crate::UID;
use crate::Vector;
use crate::Wall;
use crate::geom::bboxes::bounding_box;
use crate::{HasMesh, Mesh};
use crate::{HasName, SortByName};
use anyhow::{Result, anyhow};
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct Zone {
    pub name: String,
    pub uid: UID,
    pub parent: Option<UID>,
    solids: HashMap<String, Solid>,
}

impl HasName for Zone {
    fn get_name(&self) -> &str {
        &self.name
    }
}

impl HasMesh for Zone {
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

impl Zone {
    /// Creates a new zone with the given name and solids.
    pub fn new(name: &str, mut solids: Vec<Solid>) -> Self {
        let uid = UID::new();
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

    /// Returns solids sorted by name.
    pub fn solids(&self) -> Vec<&Solid> {
        let mut solids: Vec<&Solid> = self.solids.values().collect();
        solids.as_mut_slice().sort_by_name();
        solids
    }

    /// Returns walls from all solids, sorted by name.
    pub fn walls(&self) -> Vec<&Wall> {
        let solids = self.solids();
        let walls: Vec<&Wall> = solids.iter().flat_map(|s| s.walls()).collect();
        walls
    }

    /// Returns polygons from all solids, sorted by name.
    pub fn polygons(&self) -> Vec<&Polygon> {
        let walls = self.walls();
        let polygons: Vec<&Polygon> = walls.iter().flat_map(|w| w.polygons()).collect();
        polygons
    }

    /// Adds a solid to the zone.
    pub fn add_solid(&mut self, mut solid: Solid) -> Result<()> {
        if self.solids.contains_key(&solid.name) {
            return Err(anyhow!("Solid is already present: {}", &solid.name));
        }
        solid.parent = Some(self.uid.clone());
        self.solids.insert(solid.name.clone(), solid);
        Ok(())
    }

    /// Rotates the zone around a vector by an angle (in place).
    pub fn rotate(&mut self, angle: f64, rot_vec: &Vector) {
        for solid in self.solids.values_mut() {
            solid.rotate(angle, rot_vec);
        }
    }

    /// Translates the zone by a vector (in place).
    pub fn translate(&mut self, vec: &Vector) {
        for solid in self.solids.values_mut() {
            solid.translate(vec);
        }
    }

    /// Calculates the total volume of all solids in the zone.
    pub fn volume(&self) -> f64 {
        self.solids().iter().map(|s| s.volume()).sum()
    }

    /// Returns the bounding box of the zone as (min_point, max_point).
    pub fn bbox(&self) -> (Point, Point) {
        let mesh = self.copy_mesh();
        bounding_box(&mesh.vertices)
    }

    /// Checks if a point lies inside any solid in the zone.
    pub fn is_point_inside(&self, ptest: Point) -> bool {
        for solid in self.solids.values() {
            if solid.is_point_inside(ptest) {
                return true;
            }
        }
        false
    }

    /// Checks if a point lies on the boundary of any solid in the zone.
    pub fn is_point_at_boundary(&self, ptest: Point) -> bool {
        for solid in self.solids.values() {
            if solid.is_point_at_boundary(ptest) {
                return true;
            }
        }
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::IsClose;

    #[test]
    fn test_zone_creation() {
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1");
        let s2 = Solid::from_box(2.0, 2.0, 2.0, Some((2.0, 0.0, 0.0)), "box2");
        let zone = Zone::new("zone1", vec![s1, s2]);

        assert_eq!(zone.name, "zone1");
        assert_eq!(zone.solids().len(), 2);
    }

    #[test]
    fn test_zone_volume() {
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1");
        let s2 = Solid::from_box(1.0, 2.0, 3.0, None, "box2");
        let zone = Zone::new("zone1", vec![s1, s2]);

        let expected = 1.0 + 6.0; // 1x1x1 + 1x2x3
        let actual = zone.volume();
        assert!((actual - expected).abs() < 1e-4);
    }

    #[test]
    fn test_cube_volume() {
        // Test that a 2x2x2 cube works
        let s = Solid::from_box(2.0, 2.0, 2.0, None, "cube");
        let expected = 8.0;
        eprintln!("Cube volume: {}", s.volume());
        assert!((s.volume() - expected).abs() < 1e-4);
    }

    #[test]
    fn test_zone_add_solid() -> Result<()> {
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1");
        let mut zone = Zone::new("zone1", vec![s1]);

        assert_eq!(zone.solids().len(), 1);

        let s2 = Solid::from_box(2.0, 2.0, 2.0, None, "box2");
        zone.add_solid(s2)?;

        assert_eq!(zone.solids().len(), 2);
        Ok(())
    }

    #[test]
    fn test_zone_add_duplicate_solid() {
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1");
        let mut zone = Zone::new("zone1", vec![s1]);

        let s2 = Solid::from_box(2.0, 2.0, 2.0, None, "box1"); // Same name
        let result = zone.add_solid(s2);

        assert!(result.is_err());
    }

    #[test]
    fn test_zone_walls_and_polygons() {
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1");
        let zone = Zone::new("zone1", vec![s1]);

        // A box has 6 walls
        assert_eq!(zone.walls().len(), 6);
        // A box has 6 polygons (one per wall)
        assert_eq!(zone.polygons().len(), 6);
    }

    #[test]
    fn test_zone_translate() {
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1");
        let mut zone = Zone::new("zone1", vec![s1]);

        let (min_before, _) = zone.bbox();
        assert!(min_before.is_close(&Point::new(0.0, 0.0, 0.0)));

        zone.translate(&Vector::new(5.0, 5.0, 5.0));

        let (min_after, _) = zone.bbox();
        assert!(min_after.is_close(&Point::new(5.0, 5.0, 5.0)));
    }

    #[test]
    fn test_zone_point_inside() {
        let s1 = Solid::from_box(2.0, 2.0, 2.0, None, "box1");
        let zone = Zone::new("zone1", vec![s1]);

        // Inside the box
        assert!(zone.is_point_inside(Point::new(1.0, 1.0, 1.0)));

        // Outside the box
        assert!(!zone.is_point_inside(Point::new(5.0, 5.0, 5.0)));
    }

    #[test]
    fn test_zone_bbox() {
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1");
        let s2 = Solid::from_box(1.0, 1.0, 1.0, Some((2.0, 2.0, 2.0)), "box2");
        let zone = Zone::new("zone1", vec![s1, s2]);

        let (min, max) = zone.bbox();
        assert!(min.is_close(&Point::new(0.0, 0.0, 0.0)));
        assert!(max.is_close(&Point::new(3.0, 3.0, 3.0)));
    }
}
