//! Thermal boundary classification for energy simulations.
//!
//! This module imposes thermal semantics ("exterior", "interface") onto polygon
//! surfaces as an overlay keyed by [`UID`]. No domain-specific fields are added
//! to geometry types.

use std::collections::HashMap;

use crate::geom::building::graph::{GraphLevel, GraphParams, get_graph_edges};
use crate::sim::index::SurfaceIndex;
use crate::{Building, UID};

/// Information about a polygon interface to other polygon(s).
#[derive(Debug, Clone, Default)]
pub struct ThermalInterface {
    /// True if the polygon faces another polygon in the same zone (internal partition
    /// between solids in the same thermal zone).
    pub same_zone: bool,
    /// Zone UIDs of neighboring zones that this polygon interfaces with.
    pub other_zones: Vec<UID>,
}

/// Thermal boundary overlay for all polygon surfaces in a building.
#[derive(Debug, Clone, Default)]
pub struct ThermalBoundaries {
    /// Interface information for non-exterior polygons, keyed by polygon UID.
    ///
    /// Any polygon not present in this map is treated as `Exterior`.
    pub interfaces: HashMap<UID, ThermalInterface>,
}

impl ThermalBoundaries {
    /// Classifies polygon surfaces using the polygon-facing graph:
    /// - Polygons that face another polygon are interfaces (same-zone or inter-zone).
    /// - Polygons with no facing neighbor are treated as exterior surfaces.
    pub fn classify(building: &Building, index: &SurfaceIndex) -> Self {
        let params = GraphParams {
            level: GraphLevel::Polygon,
            facing: true,
            touching: false,
        };
        let edges = get_graph_edges(building, params);

        let mut interfaces: HashMap<UID, ThermalInterface> = HashMap::new();

        for edge in edges {
            if edge.relationship != "facing" {
                continue;
            }

            let Some(uid1) = index.polygon_uid_by_path(&edge.path1).cloned() else {
                continue;
            };
            let Some(uid2) = index.polygon_uid_by_path(&edge.path2).cloned() else {
                continue;
            };

            let zone1 = index.zone_uid_by_polygon_uid(&uid1).cloned();
            let zone2 = index.zone_uid_by_polygon_uid(&uid2).cloned();

            match (zone1, zone2) {
                (Some(z1), Some(z2)) if z1 == z2 => {
                    interfaces.entry(uid1).or_default().same_zone = true;
                    interfaces.entry(uid2).or_default().same_zone = true;
                }
                (Some(z1), Some(z2)) => {
                    interfaces
                        .entry(uid1)
                        .or_default()
                        .other_zones
                        .push(z2.clone());
                    interfaces
                        .entry(uid2)
                        .or_default()
                        .other_zones
                        .push(z1.clone());
                }
                _ => {
                    // If we cannot resolve zones for some reason, still mark as non-exterior.
                    interfaces.entry(uid1).or_default();
                    interfaces.entry(uid2).or_default();
                }
            }
        }

        Self { interfaces }
    }

    /// Returns true if the polygon is considered exterior for thermal transmission.
    pub fn is_exterior(&self, polygon_uid: &UID) -> bool {
        !self.interfaces.contains_key(polygon_uid)
    }

    pub fn interface(&self, polygon_uid: &UID) -> Option<&ThermalInterface> {
        self.interfaces.get(polygon_uid)
    }
}
