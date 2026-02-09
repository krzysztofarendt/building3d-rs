//! Cross-domain surface semantics overlay.
//!
//! This module classifies polygon surfaces using the polygon-facing graph and
//! provides a shared, UID-keyed "surface kind" that multiple physics domains
//! can use consistently (acoustics, lighting/shortwave, thermal).
//!
//! The key invariant is that **same-zone interfaces** (two facing polygons
//! belonging to the same zone) represent *internal* faces created by splitting
//! a zone volume into multiple solids. Simulations should treat these as
//! transparent/ignored by default so that results do not change when a zone is
//! decomposed into multiple solids.

use std::collections::HashMap;

use crate::geom::building::graph::{GraphLevel, GraphParams, get_graph_edges};
use crate::sim::index::SurfaceIndex;
use crate::{Building, UID};

/// Information about a polygon interface to other polygon(s).
#[derive(Debug, Clone, Default)]
pub struct SurfaceInterface {
    /// True if the polygon faces another polygon in the same zone.
    pub same_zone: bool,
    /// Zone UIDs of neighboring zones that this polygon interfaces with.
    pub other_zones: Vec<UID>,
}

/// Classification of a polygon surface.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SurfaceKind {
    /// No facing polygon neighbor exists (exterior envelope surface).
    Exterior,
    /// Faces another polygon in the same zone.
    SameZoneInterface,
    /// Faces another polygon in a different zone.
    InterZoneInterface,
    /// Faces another polygon, but zones could not be resolved.
    UnknownInterface,
}

/// Surface semantics overlay for all polygon surfaces in a building.
#[derive(Debug, Clone, Default)]
pub struct SurfaceSemantics {
    /// Interface information for non-exterior polygons, keyed by polygon UID.
    ///
    /// Any polygon not present in this map is treated as `Exterior`.
    pub interfaces: HashMap<UID, SurfaceInterface>,
    /// Cached list of facing polygon UID pairs, as reported by the polygon-facing graph.
    ///
    /// Each pair appears at most once.
    pub facing_pairs: Vec<(UID, UID)>,
}

impl SurfaceSemantics {
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

        let mut interfaces: HashMap<UID, SurfaceInterface> = HashMap::new();
        let mut facing_pairs = Vec::new();

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

            facing_pairs.push((uid1.clone(), uid2.clone()));

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

        Self {
            interfaces,
            facing_pairs,
        }
    }

    /// Returns true if the polygon is considered exterior.
    pub fn is_exterior(&self, polygon_uid: &UID) -> bool {
        !self.interfaces.contains_key(polygon_uid)
    }

    /// Returns true if the polygon is a same-zone interface.
    pub fn is_same_zone_interface(&self, polygon_uid: &UID) -> bool {
        self.interfaces
            .get(polygon_uid)
            .map(|i| i.same_zone)
            .unwrap_or(false)
    }

    pub fn interface(&self, polygon_uid: &UID) -> Option<&SurfaceInterface> {
        self.interfaces.get(polygon_uid)
    }

    pub fn kind(&self, polygon_uid: &UID) -> SurfaceKind {
        let Some(i) = self.interfaces.get(polygon_uid) else {
            return SurfaceKind::Exterior;
        };
        if i.same_zone {
            return SurfaceKind::SameZoneInterface;
        }
        if !i.other_zones.is_empty() {
            return SurfaceKind::InterZoneInterface;
        }
        SurfaceKind::UnknownInterface
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Solid, Zone};

    #[test]
    fn test_same_zone_interfaces_classified() {
        // Two adjacent boxes in the same zone — their shared face should be SameZoneInterface.
        let s0 = Solid::from_box(1.0, 1.0, 1.0, None, "s0").unwrap();
        let s1 = Solid::from_box(1.0, 1.0, 1.0, Some((1.0, 0.0, 0.0)), "s1").unwrap();
        let zone = Zone::new("z", vec![s0, s1]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let index = SurfaceIndex::new(&building);
        let sem = SurfaceSemantics::classify(&building, &index);

        let same_zone_count = index
            .surfaces
            .iter()
            .filter(|s| sem.is_same_zone_interface(&s.polygon_uid))
            .count();
        assert!(same_zone_count > 0);

        // Adjacent solids share exactly one interface: two facing polygons.
        assert_eq!(same_zone_count, 2);
    }

    #[test]
    fn test_different_zone_interfaces_not_same_zone() {
        let s0 = Solid::from_box(1.0, 1.0, 1.0, None, "s0").unwrap();
        let s1 = Solid::from_box(1.0, 1.0, 1.0, Some((1.0, 0.0, 0.0)), "s1").unwrap();
        let z0 = Zone::new("z0", vec![s0]).unwrap();
        let z1 = Zone::new("z1", vec![s1]).unwrap();
        let building = Building::new("b", vec![z0, z1]).unwrap();

        let index = SurfaceIndex::new(&building);
        let sem = SurfaceSemantics::classify(&building, &index);

        let same_zone_count = index
            .surfaces
            .iter()
            .filter(|s| sem.is_same_zone_interface(&s.polygon_uid))
            .count();
        assert_eq!(same_zone_count, 0);
    }
}
