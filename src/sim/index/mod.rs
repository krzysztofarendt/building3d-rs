//! Building geometry indexing for simulations.
//!
//! This module provides a stable, domain-agnostic index of polygon surfaces.
//! Simulations attach their own metadata overlays keyed by [`UID`], rather than
//! storing domain-specific fields on geometry types.

use std::collections::HashMap;

use crate::{Building, UID};

/// Reference information for a polygon surface in the building hierarchy.
#[derive(Debug, Clone)]
pub struct SurfaceRef {
    /// Full path: `zone/solid/wall/polygon`.
    pub path: String,
    /// Polygon UID.
    pub polygon_uid: UID,
    /// Owning zone UID.
    pub zone_uid: UID,
    /// Zone name (path component).
    pub zone_name: String,
}

/// Domain-agnostic lookup table for polygon surfaces.
#[derive(Debug, Clone)]
pub struct SurfaceIndex {
    pub surfaces: Vec<SurfaceRef>,
    path_to_polygon_uid: HashMap<String, UID>,
    polygon_uid_to_zone_uid: HashMap<UID, UID>,
    polygon_uid_to_path: HashMap<UID, String>,
}

impl SurfaceIndex {
    pub fn new(building: &Building) -> Self {
        let mut surfaces = Vec::new();
        let mut path_to_polygon_uid = HashMap::new();
        let mut polygon_uid_to_zone_uid = HashMap::new();
        let mut polygon_uid_to_path = HashMap::new();

        for zone in building.zones() {
            for solid in zone.solids() {
                for wall in solid.walls() {
                    for polygon in wall.polygons() {
                        let path = format!(
                            "{}/{}/{}/{}",
                            zone.name, solid.name, wall.name, polygon.name
                        );
                        let polygon_uid = polygon.uid.clone();
                        let zone_uid = zone.uid.clone();

                        path_to_polygon_uid.insert(path.clone(), polygon_uid.clone());
                        polygon_uid_to_zone_uid.insert(polygon_uid.clone(), zone_uid.clone());
                        polygon_uid_to_path.insert(polygon_uid.clone(), path.clone());

                        surfaces.push(SurfaceRef {
                            path,
                            polygon_uid,
                            zone_uid,
                            zone_name: zone.name.clone(),
                        });
                    }
                }
            }
        }

        Self {
            surfaces,
            path_to_polygon_uid,
            polygon_uid_to_zone_uid,
            polygon_uid_to_path,
        }
    }

    pub fn polygon_uid_by_path(&self, path: &str) -> Option<&UID> {
        self.path_to_polygon_uid.get(path)
    }

    pub fn zone_uid_by_polygon_uid(&self, polygon_uid: &UID) -> Option<&UID> {
        self.polygon_uid_to_zone_uid.get(polygon_uid)
    }

    pub fn path_by_polygon_uid(&self, polygon_uid: &UID) -> Option<&str> {
        self.polygon_uid_to_path
            .get(polygon_uid)
            .map(|s| s.as_str())
    }
}
