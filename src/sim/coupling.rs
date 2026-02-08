//! Cross-domain coupling payloads for composed simulations.
//!
//! These types define small, stable "contracts" that independent simulation
//! modules can exchange via the [`crate::sim::framework::Bus`].
//!
//! The core guideline is:
//! - key cross-module outputs by stable `UID` (not by path strings),
//! - keep core values in radiometric SI units (W, W/m²),
//! - convert to photometric units only at output/reporting boundaries.

use std::collections::HashMap;

use crate::UID;

/// Shortwave radiant power absorbed at polygon surfaces, keyed by polygon `UID` [W].
///
/// Note: "shortwave" usually means solar radiation. If a lighting simulation includes
/// artificial light sources, callers should avoid treating this as solar gain unless
/// sources are restricted appropriately.
#[derive(Debug, Clone, Default)]
pub struct ShortwaveAbsorbedWPerPolygon {
    pub watts_by_polygon_uid: HashMap<UID, f64>,
}

/// Shortwave radiant power transmitted into each zone, keyed by zone `UID` [W].
#[derive(Debug, Clone, Default)]
pub struct ShortwaveTransmittedWPerZone {
    pub watts_by_zone_uid: HashMap<UID, f64>,
}
