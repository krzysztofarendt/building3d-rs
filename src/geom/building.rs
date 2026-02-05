//! Building container - the top level of the hierarchy.
//!
//! Hierarchy: Building → Zone → Solid → Wall → Polygon

use crate::Point;
use crate::Polygon;
use crate::Solid;
use crate::TriangleIndex;
use crate::UID;
use crate::Vector;
use crate::Wall;
use crate::geom::bboxes::bounding_box;
use crate::geom::zone::Zone;
use crate::{HasMesh, Mesh};
use crate::{HasName, SortByName};
use anyhow::{Result, anyhow};
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct Building {
    pub name: String,
    pub uid: UID,
    pub parent: Option<UID>,
    zones: HashMap<String, Zone>,
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
    /// Creates a new building with the given name and zones.
    pub fn new(name: &str, mut zones: Vec<Zone>) -> Self {
        let uid = UID::new();
        for z in zones.iter_mut() {
            z.parent = Some(uid.clone());
        }
        let parent = None;
        let zones: HashMap<String, Zone> = zones.into_iter().map(|x| (x.name.clone(), x)).collect();

        Self {
            name: name.to_string(),
            uid,
            parent,
            zones,
        }
    }

    /// Creates a new building from solids directly (convenience method).
    ///
    /// Each solid is wrapped in its own zone with the same name.
    pub fn from_solids(name: &str, solids: Vec<Solid>) -> Self {
        let zones: Vec<Zone> = solids
            .into_iter()
            .map(|s| {
                let zone_name = s.name.clone();
                Zone::new(&zone_name, vec![s])
            })
            .collect();
        Self::new(name, zones)
    }

    /// Returns zones sorted by name.
    pub fn zones(&self) -> Vec<&Zone> {
        let mut zones: Vec<&Zone> = self.zones.values().collect();
        zones.as_mut_slice().sort_by_name();
        zones
    }

    /// Returns solids from all zones, sorted by name.
    pub fn solids(&self) -> Vec<&Solid> {
        let zones = self.zones();
        let solids: Vec<&Solid> = zones.iter().flat_map(|z| z.solids()).collect();
        solids
    }

    /// Returns walls from all zones, sorted by name.
    pub fn walls(&self) -> Vec<&Wall> {
        let solids = self.solids();
        let walls: Vec<&Wall> = solids.iter().flat_map(|s| s.walls()).collect();
        walls
    }

    /// Returns polygons from all zones, sorted by name.
    pub fn polygons(&self) -> Vec<&Polygon> {
        let walls = self.walls();
        let polygons: Vec<&Polygon> = walls.iter().flat_map(|w| w.polygons()).collect();
        polygons
    }

    /// Adds a zone to the building.
    pub fn add_zone(&mut self, mut zone: Zone) -> Result<()> {
        if self.zones.contains_key(&zone.name) {
            return Err(anyhow!("Zone is already present: {}", &zone.name));
        }
        zone.parent = Some(self.uid.clone());
        self.zones.insert(zone.name.clone(), zone);
        Ok(())
    }

    /// Rotates the building around a vector by an angle (in place).
    pub fn rotate(&mut self, angle: f64, rot_vec: &Vector) {
        for zone in self.zones.values_mut() {
            zone.rotate(angle, rot_vec);
        }
    }

    /// Translates the building by a vector (in place).
    pub fn translate(&mut self, vec: &Vector) {
        for zone in self.zones.values_mut() {
            zone.translate(vec);
        }
    }

    /// Calculates the total volume of all solids in the building.
    pub fn volume(&self) -> f64 {
        self.zones().iter().map(|z| z.volume()).sum()
    }

    /// Returns the bounding box of the building as (min_point, max_point).
    pub fn bbox(&self) -> (Point, Point) {
        let mesh = self.copy_mesh();
        bounding_box(&mesh.vertices)
    }

    /// Checks if a point lies inside any solid in the building.
    pub fn is_point_inside(&self, ptest: Point) -> bool {
        for zone in self.zones.values() {
            if zone.is_point_inside(ptest) {
                return true;
            }
        }
        false
    }

    /// Checks if a point lies on the boundary of any solid in the building.
    pub fn is_point_at_boundary(&self, ptest: Point) -> bool {
        for zone in self.zones.values() {
            if zone.is_point_at_boundary(ptest) {
                return true;
            }
        }
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_building_from_zones() {
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1");
        let z1 = Zone::new("zone1", vec![s1]);

        let s2 = Solid::from_box(2.0, 2.0, 2.0, None, "box2");
        let z2 = Zone::new("zone2", vec![s2]);

        let bdg = Building::new("building", vec![z1, z2]);

        assert_eq!(bdg.zones().len(), 2);
        assert_eq!(bdg.solids().len(), 2);
    }

    #[test]
    fn test_building_from_solids() {
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1");
        let s2 = Solid::from_box(2.0, 2.0, 2.0, None, "box2");

        let bdg = Building::from_solids("building", vec![s1, s2]);

        assert_eq!(bdg.zones().len(), 2);
        assert_eq!(bdg.solids().len(), 2);
    }

    #[test]
    fn test_volume() {
        let s0 = Solid::from_box(1., 1., 1., None, "box_0");
        let s1 = Solid::from_box(1., 2., 3., None, "box_1");
        let bdg = Building::from_solids("building", vec![s0, s1]);
        let expected_vol = 1. * 2. * 3. + 1.;
        assert!((bdg.volume() - expected_vol).abs() < 1e-4);
    }

    #[test]
    fn test_building_add_zone() -> Result<()> {
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1");
        let z1 = Zone::new("zone1", vec![s1]);
        let mut bdg = Building::new("building", vec![z1]);

        assert_eq!(bdg.zones().len(), 1);

        let s2 = Solid::from_box(2.0, 2.0, 2.0, None, "box2");
        let z2 = Zone::new("zone2", vec![s2]);
        bdg.add_zone(z2)?;

        assert_eq!(bdg.zones().len(), 2);
        Ok(())
    }

    #[test]
    fn test_building_translate() {
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1");
        let z1 = Zone::new("zone1", vec![s1]);
        let mut bdg = Building::new("building", vec![z1]);

        let (min_before, _) = bdg.bbox();
        assert!(min_before.is_close(&Point::new(0.0, 0.0, 0.0)));

        bdg.translate(&Vector::new(10.0, 10.0, 10.0));

        let (min_after, _) = bdg.bbox();
        assert!(min_after.is_close(&Point::new(10.0, 10.0, 10.0)));
    }

    #[test]
    fn test_building_point_inside() {
        let s1 = Solid::from_box(2.0, 2.0, 2.0, None, "box1");
        let z1 = Zone::new("zone1", vec![s1]);
        let bdg = Building::new("building", vec![z1]);

        // Inside
        assert!(bdg.is_point_inside(Point::new(1.0, 1.0, 1.0)));

        // Outside
        assert!(!bdg.is_point_inside(Point::new(5.0, 5.0, 5.0)));
    }
}
