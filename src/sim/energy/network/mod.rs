//! Thermal network representation (foundation for pluggable heat-transfer models).
//!
//! The goal is to represent heat transfer as a graph of nodes (temperatures)
//! connected by components (conductances/capacitances/sources). This provides a
//! stable integration surface for personalized simulations: wall models can be
//! swapped (steady-state U, 2R1C, FD, 3D FEM, ...) without changing the rest of
//! the simulation pipeline.
//!
//! Important: thermal semantics like "exterior" or "inter-zone" are imposed via
//! overlays (e.g. [`ThermalBoundaries`]) keyed by [`UID`], not stored on geometry.

mod multizone;
mod multizone_envelope;
pub(crate) mod solve;

use std::collections::{HashMap, HashSet};

use crate::geom::polygon::relations::polygon_overlap_area;
use crate::sim::energy::boundary::ThermalBoundaries;
use crate::sim::energy::config::ThermalConfig;
use crate::sim::index::SurfaceIndex;
use crate::{Building, UID};

pub use multizone::{MultiZoneAirModel, MultiZoneStepResult};
pub use multizone_envelope::MultiZoneEnvelopeRcModel;

/// Inter-zone conductance between two thermal zones (W/K).
#[derive(Debug, Clone)]
pub struct InterZoneConductance {
    pub zone_a: UID,
    pub zone_b: UID,
    pub conductance_w_per_k: f64,
}

/// Minimal thermal network summary for zone-air-only models.
///
/// This intentionally only covers "conductance to outside" and "conductance
/// between zones". Future steps can add node capacities, sources, and additional
/// node types (surface temperatures, radiant node, etc.).
#[derive(Debug, Clone)]
pub struct ThermalNetwork {
    exterior_conductance_w_per_k: HashMap<UID, f64>,
    interzone_conductance_w_per_k: HashMap<(UID, UID), f64>,
}

impl ThermalNetwork {
    pub fn build(
        building: &Building,
        config: &ThermalConfig,
        index: &SurfaceIndex,
        boundaries: &ThermalBoundaries,
    ) -> Self {
        Self::build_internal(building, config, index, boundaries, None)
    }

    pub fn build_with_ignored_exterior_polygons(
        building: &Building,
        config: &ThermalConfig,
        index: &SurfaceIndex,
        boundaries: &ThermalBoundaries,
        ignored_exterior_polygon_uids: &HashSet<UID>,
    ) -> Self {
        Self::build_internal(
            building,
            config,
            index,
            boundaries,
            Some(ignored_exterior_polygon_uids),
        )
    }

    fn build_internal(
        building: &Building,
        config: &ThermalConfig,
        index: &SurfaceIndex,
        boundaries: &ThermalBoundaries,
        ignored_exterior_polygon_uids: Option<&HashSet<UID>>,
    ) -> Self {
        let mut exterior_conductance_w_per_k: HashMap<UID, f64> = HashMap::new();

        // Exterior conductance: sum(U*A) over exterior surfaces.
        for surface in &index.surfaces {
            if !boundaries.is_exterior(&surface.polygon_uid) {
                continue;
            }
            if let Some(ignored) = ignored_exterior_polygon_uids
                && ignored.contains(&surface.polygon_uid)
            {
                continue;
            }
            let u = config.resolve_u_value_for_surface(&surface.polygon_uid, &surface.path);
            *exterior_conductance_w_per_k
                .entry(surface.zone_uid.clone())
                .or_insert(0.0) += u * surface.area_m2;
        }

        // Inter-zone conductance from facing polygon pairs.
        let mut interzone_conductance_w_per_k: HashMap<(UID, UID), f64> = HashMap::new();
        for (poly1_uid, poly2_uid) in &boundaries.facing_pairs {
            let (Some(zone1_uid), Some(zone2_uid)) = (
                index.zone_uid_by_polygon_uid(poly1_uid).cloned(),
                index.zone_uid_by_polygon_uid(poly2_uid).cloned(),
            ) else {
                continue;
            };

            if zone1_uid == zone2_uid {
                continue;
            }

            let (Some(path1), Some(path2)) = (
                index.path_by_polygon_uid(poly1_uid),
                index.path_by_polygon_uid(poly2_uid),
            ) else {
                continue;
            };

            let (Some(poly1), Some(poly2)) =
                (building.get_polygon(path1), building.get_polygon(path2))
            else {
                continue;
            };

            let overlap_area = polygon_overlap_area(poly1, poly2);
            if overlap_area <= 0.0 {
                continue;
            }

            let u1 = config.resolve_u_value_for_surface(poly1_uid, path1);
            let u2 = config.resolve_u_value_for_surface(poly2_uid, path2);

            let k = config.interzone_conductance_w_per_k(u1, u2, overlap_area);
            if k <= 0.0 {
                continue;
            }

            let key = canonical_zone_pair(&zone1_uid, &zone2_uid);
            *interzone_conductance_w_per_k.entry(key).or_insert(0.0) += k;
        }

        Self {
            exterior_conductance_w_per_k,
            interzone_conductance_w_per_k,
        }
    }

    /// Total exterior conductance UA for the building (W/K).
    pub fn exterior_conductance_total_w_per_k(&self) -> f64 {
        self.exterior_conductance_w_per_k.values().sum()
    }

    pub fn exterior_conductance_by_zone_w_per_k(&self, zone_uid: &UID) -> f64 {
        self.exterior_conductance_w_per_k
            .get(zone_uid)
            .cloned()
            .unwrap_or(0.0)
    }

    /// Returns all inter-zone conductances as a flat list.
    pub fn interzone_conductances(&self) -> Vec<InterZoneConductance> {
        self.interzone_conductance_w_per_k
            .iter()
            .map(|((a, b), k)| InterZoneConductance {
                zone_a: a.clone(),
                zone_b: b.clone(),
                conductance_w_per_k: *k,
            })
            .collect()
    }
}

fn canonical_zone_pair(a: &UID, b: &UID) -> (UID, UID) {
    if a.as_str() <= b.as_str() {
        (a.clone(), b.clone())
    } else {
        (b.clone(), a.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::energy::zone::calculate_heat_balance_with_boundaries;
    use crate::{Solid, Zone};

    #[test]
    fn test_network_matches_exterior_ua_simple_box() {
        let s = Solid::from_box(5.0, 5.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let index = SurfaceIndex::new(&building);
        let boundaries = ThermalBoundaries::classify(&building, &index);

        let config = ThermalConfig::new();
        let network = ThermalNetwork::build(&building, &config, &index, &boundaries);

        let dt = config.indoor_temperature - config.outdoor_temperature;
        let steady = calculate_heat_balance_with_boundaries(&building, &config, &boundaries);
        let ua_from_steady = steady.transmission_loss / dt;

        assert!(
            (network.exterior_conductance_total_w_per_k() - ua_from_steady).abs() < 1e-10,
            "Network UA should match steady-state UA"
        );
    }

    #[test]
    fn test_network_excludes_internal_interfaces_same_zone() {
        let s0 = Solid::from_box(1.0, 1.0, 1.0, None, "s0").unwrap();
        let s1 = Solid::from_box(1.0, 1.0, 1.0, Some((1.0, 0.0, 0.0)), "s1").unwrap();
        let zone = Zone::new("z", vec![s0, s1]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let index = SurfaceIndex::new(&building);
        let boundaries = ThermalBoundaries::classify(&building, &index);

        let config = ThermalConfig::new();
        let network = ThermalNetwork::build(&building, &config, &index, &boundaries);

        // Exterior area of 2x1x1 prism is 10 m².
        let expected_ua = config.default_u_value * 10.0;
        assert!(
            (network.exterior_conductance_total_w_per_k() - expected_ua).abs() < 1e-10,
            "Expected UA={expected_ua}, got {}",
            network.exterior_conductance_total_w_per_k()
        );
    }

    #[test]
    fn test_network_interzone_conductance_detected() {
        let s0 = Solid::from_box(1.0, 1.0, 1.0, None, "s0").unwrap();
        let s1 = Solid::from_box(1.0, 1.0, 1.0, Some((1.0, 0.0, 0.0)), "s1").unwrap();
        let z0 = Zone::new("z0", vec![s0]).unwrap();
        let z1 = Zone::new("z1", vec![s1]).unwrap();
        let building = Building::new("b", vec![z0, z1]).unwrap();

        let index = SurfaceIndex::new(&building);
        let boundaries = ThermalBoundaries::classify(&building, &index);

        let config = ThermalConfig::new();
        let network = ThermalNetwork::build(&building, &config, &index, &boundaries);

        // Each zone has 5 m² exterior (one face becomes inter-zone interface).
        let expected_exterior_ua_total = config.default_u_value * (5.0 + 5.0);
        assert!(
            (network.exterior_conductance_total_w_per_k() - expected_exterior_ua_total).abs()
                < 1e-10
        );

        let interzone = network.interzone_conductances();
        assert_eq!(interzone.len(), 1, "Should detect one inter-zone interface");

        // Default U=2.0, overlap area=1.0 => conductance=2.0 W/K.
        assert!(
            (interzone[0].conductance_w_per_k - 2.0).abs() < 1e-10,
            "Expected inter-zone conductance ~2.0, got {}",
            interzone[0].conductance_w_per_k
        );
    }
}
