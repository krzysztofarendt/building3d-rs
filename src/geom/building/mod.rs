//! Building container - the top level of the hierarchy.
//!
//! Hierarchy: Building → Zone → Solid → Wall → Polygon

pub mod graph;

use crate::Point;
use crate::Polygon;
use crate::Solid;
use crate::TriangleIndex;
use crate::UID;
use crate::Vector;
use crate::Wall;
use crate::geom;
use crate::geom::bboxes::bounding_box;
use crate::geom::zone::Zone;
use crate::{HasMesh, Mesh};
use crate::{HasName, SortByName};
use anyhow::{Result, anyhow};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, Serialize, Deserialize)]
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

impl Building {
    /// Creates a new building with the given name and zones.
    pub fn new(name: &str, mut zones: Vec<Zone>) -> Result<Self> {
        let name = geom::validate_name(name)?;
        let uid = UID::new();
        for z in zones.iter_mut() {
            z.parent = Some(uid.clone());
        }
        let parent = None;
        let mut map: HashMap<String, Zone> = HashMap::new();
        for zone in zones {
            if map.contains_key(&zone.name) {
                return Err(anyhow!(
                    "Zone is already present in Building::new(): {}",
                    &zone.name
                ));
            }
            map.insert(zone.name.clone(), zone);
        }

        Ok(Self {
            name: name.to_string(),
            uid,
            parent,
            zones: map,
        })
    }

    /// Creates a new building from solids directly (convenience method).
    ///
    /// Each solid is wrapped in its own zone with the same name.
    pub fn from_solids(name: &str, solids: Vec<Solid>) -> Result<Self> {
        let mut zones: Vec<Zone> = Vec::new();
        for s in solids {
            let zone_name = s.name.clone();
            zones.push(Zone::new(&zone_name, vec![s])?);
        }
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
        zone.name = geom::validate_name(&zone.name)?;
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

    /// Builds an adjacency graph for the building.
    ///
    /// Returns a HashMap mapping entity paths to lists of adjacent entity paths.
    /// Use GraphParams to control the level (polygon/wall/solid/zone) and
    /// relationship types (facing/touching).
    pub fn get_graph(&self, params: graph::GraphParams) -> HashMap<String, Vec<String>> {
        graph::get_graph(self, params)
    }

    /// Returns all edges in the building graph as a flat list.
    pub fn get_graph_edges(&self, params: graph::GraphParams) -> Vec<graph::GraphEdge> {
        graph::get_graph_edges(self, params)
    }

    /// Finds all interfaces between adjacent solids.
    ///
    /// Returns detailed information about each interface, including whether
    /// the interface is correct (faces match exactly).
    pub fn stitch_solids(&self) -> Vec<graph::StitchInfo> {
        graph::stitch_solids(self)
    }

    /// Gets a zone by path (zone_name).
    pub fn get_zone(&self, path: &str) -> Option<&Zone> {
        self.zones.get(path)
    }

    /// Gets a solid by path (zone_name/solid_name).
    pub fn get_solid(&self, path: &str) -> Option<&Solid> {
        let parts: Vec<&str> = path.split('/').collect();
        if parts.len() != 2 {
            return None;
        }
        self.zones.get(parts[0])?.get_solid(parts[1])
    }

    /// Gets a wall by path (zone_name/solid_name/wall_name).
    pub fn get_wall(&self, path: &str) -> Option<&Wall> {
        let parts: Vec<&str> = path.split('/').collect();
        if parts.len() != 3 {
            return None;
        }
        self.zones
            .get(parts[0])?
            .get_solid(parts[1])?
            .get_wall(parts[2])
    }

    /// Gets a polygon by path (zone_name/solid_name/wall_name/polygon_name).
    pub fn get_polygon(&self, path: &str) -> Option<&Polygon> {
        let parts: Vec<&str> = path.split('/').collect();
        if parts.len() != 4 {
            return None;
        }
        self.zones
            .get(parts[0])?
            .get_solid(parts[1])?
            .get_wall(parts[2])?
            .get_polygon(parts[3])
    }

    /// Validates the structural integrity of the building.
    ///
    /// Checks for:
    /// - Duplicate UIDs across all entities
    /// - Empty containers (zones with no solids, solids with no walls, walls with no polygons)
    /// - Mesh face indices that are out of bounds
    pub fn validate(&self) -> Result<()> {
        let mut uids: HashSet<&str> = HashSet::new();

        // Check building UID
        if !uids.insert(self.uid.as_str()) {
            return Err(anyhow!("Duplicate UID: {}", self.uid.as_str()));
        }

        for zone in self.zones.values() {
            if !uids.insert(zone.uid.as_str()) {
                return Err(anyhow!("Duplicate UID: {}", zone.uid.as_str()));
            }
            if zone.solids().is_empty() {
                return Err(anyhow!("Zone '{}' has no solids", zone.name));
            }

            for solid in zone.solids() {
                if !uids.insert(solid.uid.as_str()) {
                    return Err(anyhow!("Duplicate UID: {}", solid.uid.as_str()));
                }
                if solid.walls().is_empty() {
                    return Err(anyhow!("Solid '{}' has no walls", solid.name));
                }

                for wall in solid.walls() {
                    if !uids.insert(wall.uid.as_str()) {
                        return Err(anyhow!("Duplicate UID: {}", wall.uid.as_str()));
                    }
                    if wall.polygons().is_empty() {
                        return Err(anyhow!("Wall '{}' has no polygons", wall.name));
                    }

                    for poly in wall.polygons() {
                        if !uids.insert(poly.uid.as_str()) {
                            return Err(anyhow!("Duplicate UID: {}", poly.uid.as_str()));
                        }
                        let mesh = poly.mesh_ref();
                        let vc = mesh.vertex_count();
                        if let Some(faces) = &mesh.faces {
                            for tri in faces {
                                if tri.0 >= vc || tri.1 >= vc || tri.2 >= vc {
                                    return Err(anyhow!(
                                        "Face index out of bounds in polygon '{}': ({}, {}, {}) >= {}",
                                        poly.name,
                                        tri.0,
                                        tri.1,
                                        tri.2,
                                        vc
                                    ));
                                }
                            }
                        }
                    }
                }
            }
        }

        Ok(())
    }

    pub fn repair_parents(&mut self) {
        for zone in self.zones.values_mut() {
            zone.parent = Some(self.uid.clone());
            zone.repair_parents();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_building_from_zones() -> Result<()> {
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1")?;
        let z1 = Zone::new("zone1", vec![s1])?;

        let s2 = Solid::from_box(2.0, 2.0, 2.0, None, "box2")?;
        let z2 = Zone::new("zone2", vec![s2])?;

        let bdg = Building::new("building", vec![z1, z2])?;

        assert_eq!(bdg.zones().len(), 2);
        assert_eq!(bdg.solids().len(), 2);
        Ok(())
    }

    #[test]
    fn test_building_from_solids() -> Result<()> {
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1")?;
        let s2 = Solid::from_box(2.0, 2.0, 2.0, None, "box2")?;

        let bdg = Building::from_solids("building", vec![s1, s2])?;

        assert_eq!(bdg.zones().len(), 2);
        assert_eq!(bdg.solids().len(), 2);
        Ok(())
    }

    #[test]
    fn test_volume() -> Result<()> {
        let s0 = Solid::from_box(1., 1., 1., None, "box_0")?;
        let s1 = Solid::from_box(1., 2., 3., None, "box_1")?;
        let bdg = Building::from_solids("building", vec![s0, s1])?;
        let expected_vol = 1. * 2. * 3. + 1.;
        assert!((bdg.volume() - expected_vol).abs() < 1e-4);
        Ok(())
    }

    #[test]
    fn test_building_add_zone() -> Result<()> {
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1")?;
        let z1 = Zone::new("zone1", vec![s1])?;
        let mut bdg = Building::new("building", vec![z1])?;

        assert_eq!(bdg.zones().len(), 1);

        let s2 = Solid::from_box(2.0, 2.0, 2.0, None, "box2")?;
        let z2 = Zone::new("zone2", vec![s2])?;
        bdg.add_zone(z2)?;

        assert_eq!(bdg.zones().len(), 2);
        Ok(())
    }

    #[test]
    fn test_building_translate() -> Result<()> {
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1")?;
        let z1 = Zone::new("zone1", vec![s1])?;
        let mut bdg = Building::new("building", vec![z1])?;

        let (min_before, _) = bdg.bbox();
        assert!(min_before.is_close(&Point::new(0.0, 0.0, 0.0)));

        bdg.translate(&Vector::new(10.0, 10.0, 10.0));

        let (min_after, _) = bdg.bbox();
        assert!(min_after.is_close(&Point::new(10.0, 10.0, 10.0)));
        Ok(())
    }

    #[test]
    fn test_building_point_inside() -> Result<()> {
        let s1 = Solid::from_box(2.0, 2.0, 2.0, None, "box1")?;
        let z1 = Zone::new("zone1", vec![s1])?;
        let bdg = Building::new("building", vec![z1])?;

        // Inside
        assert!(bdg.is_point_inside(Point::new(1.0, 1.0, 1.0)));

        // Outside
        assert!(!bdg.is_point_inside(Point::new(5.0, 5.0, 5.0)));
        Ok(())
    }

    #[test]
    fn test_path_access_zone() -> Result<()> {
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1")?;
        let z1 = Zone::new("zone1", vec![s1])?;
        let bdg = Building::new("building", vec![z1])?;

        // Valid path
        let zone = bdg.get_zone("zone1");
        assert!(zone.is_some());
        assert_eq!(zone.unwrap().name, "zone1");

        // Invalid path
        assert!(bdg.get_zone("nonexistent").is_none());
        Ok(())
    }

    #[test]
    fn test_path_access_solid() -> Result<()> {
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1")?;
        let z1 = Zone::new("zone1", vec![s1])?;
        let bdg = Building::new("building", vec![z1])?;

        // Valid path
        let solid = bdg.get_solid("zone1/box1");
        assert!(solid.is_some());
        assert_eq!(solid.unwrap().name, "box1");

        // Invalid paths
        assert!(bdg.get_solid("zone1").is_none()); // Wrong format
        assert!(bdg.get_solid("zone1/nonexistent").is_none());
        assert!(bdg.get_solid("nonexistent/box1").is_none());
        Ok(())
    }

    #[test]
    fn test_path_access_wall() -> Result<()> {
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1")?;
        let z1 = Zone::new("zone1", vec![s1])?;
        let bdg = Building::new("building", vec![z1])?;

        // Valid path - box has wall_0, wall_1, wall_2, wall_3, floor, ceiling
        let wall = bdg.get_wall("zone1/box1/wall_0");
        assert!(wall.is_some());
        assert_eq!(wall.unwrap().name, "wall_0");

        // Invalid paths
        assert!(bdg.get_wall("zone1/box1").is_none()); // Wrong format
        assert!(bdg.get_wall("zone1/box1/nonexistent").is_none());
        Ok(())
    }

    #[test]
    fn test_validate_ok() -> Result<()> {
        let s = Solid::from_box(1.0, 1.0, 1.0, None, "box1")?;
        let bdg = Building::from_solids("building", vec![s])?;
        bdg.validate()?;
        Ok(())
    }

    #[test]
    fn test_path_access_polygon() -> Result<()> {
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1")?;
        let z1 = Zone::new("zone1", vec![s1])?;
        let bdg = Building::new("building", vec![z1])?;

        // Valid path - floor wall has floor polygon
        let poly = bdg.get_polygon("zone1/box1/floor/floor");
        assert!(poly.is_some());
        assert_eq!(poly.unwrap().name, "floor");

        // Invalid paths
        assert!(bdg.get_polygon("zone1/box1/floor").is_none()); // Wrong format
        assert!(bdg.get_polygon("zone1/box1/floor/nonexistent").is_none());
        Ok(())
    }

    #[test]
    fn test_building_duplicate_zone_error() -> Result<()> {
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1")?;
        let z1 = Zone::new("zone1", vec![s1])?;
        let s2 = Solid::from_box(1.0, 1.0, 1.0, None, "box2")?;
        let z2 = Zone::new("zone1", vec![s2])?; // Same name
        let result = Building::new("building", vec![z1, z2]);
        assert!(result.is_err());
        Ok(())
    }

    #[test]
    fn test_building_add_duplicate_zone_error() -> Result<()> {
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1")?;
        let z1 = Zone::new("zone1", vec![s1])?;
        let mut bdg = Building::new("building", vec![z1])?;

        let s2 = Solid::from_box(1.0, 1.0, 1.0, None, "box2")?;
        let z2 = Zone::new("zone1", vec![s2])?; // Same name
        assert!(bdg.add_zone(z2).is_err());
        Ok(())
    }

    #[test]
    fn test_building_rotate() -> Result<()> {
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1")?;
        let z1 = Zone::new("zone1", vec![s1])?;
        let mut bdg = Building::new("building", vec![z1])?;

        bdg.rotate(std::f64::consts::PI / 2., &Vector::new(0., 0., 1.));
        assert_eq!(bdg.zones().len(), 1);
        assert_eq!(bdg.solids().len(), 1);
        Ok(())
    }

    #[test]
    fn test_building_is_point_at_boundary() -> Result<()> {
        let s1 = Solid::from_box(2.0, 2.0, 2.0, None, "box1")?;
        let z1 = Zone::new("zone1", vec![s1])?;
        let bdg = Building::new("building", vec![z1])?;

        // On boundary
        assert!(bdg.is_point_at_boundary(Point::new(0.0, 1.0, 1.0)));
        // Outside
        assert!(!bdg.is_point_at_boundary(Point::new(5.0, 5.0, 5.0)));
        Ok(())
    }

    #[test]
    fn test_building_get_graph() -> Result<()> {
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1")?;
        let s2 = Solid::from_box(1.0, 1.0, 1.0, Some((1.0, 0.0, 0.0)), "box2")?;
        let bdg = Building::from_solids("building", vec![s1, s2])?;

        let params = graph::GraphParams::default();
        let graph = bdg.get_graph(params);
        assert_eq!(graph.len(), 2);
        Ok(())
    }

    #[test]
    fn test_building_get_graph_edges() -> Result<()> {
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1")?;
        let s2 = Solid::from_box(1.0, 1.0, 1.0, Some((1.0, 0.0, 0.0)), "box2")?;
        let bdg = Building::from_solids("building", vec![s1, s2])?;

        let params = graph::GraphParams::default();
        let edges = bdg.get_graph_edges(params);
        assert_eq!(edges.len(), 1);
        Ok(())
    }

    #[test]
    fn test_building_stitch_solids() -> Result<()> {
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1")?;
        let s2 = Solid::from_box(1.0, 1.0, 1.0, Some((1.0, 0.0, 0.0)), "box2")?;
        let bdg = Building::from_solids("building", vec![s1, s2])?;

        let stitches = bdg.stitch_solids();
        assert_eq!(stitches.len(), 1);
        Ok(())
    }

    #[test]
    fn test_building_walls_and_polygons() -> Result<()> {
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1")?;
        let bdg = Building::from_solids("building", vec![s1])?;

        assert_eq!(bdg.walls().len(), 6);
        assert_eq!(bdg.polygons().len(), 6);
        Ok(())
    }

    #[test]
    fn test_building_has_name() -> Result<()> {
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1")?;
        let bdg = Building::from_solids("mybuilding", vec![s1])?;
        assert_eq!(bdg.get_name(), "mybuilding");
        Ok(())
    }

    #[test]
    fn test_building_repair_parents() -> Result<()> {
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1")?;
        let mut bdg = Building::from_solids("building", vec![s1])?;
        bdg.repair_parents();
        // Verify structure is still valid
        bdg.validate()?;
        Ok(())
    }
}
