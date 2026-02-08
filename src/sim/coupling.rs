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

/// Outdoor air (dry-bulb) temperature [°C].
///
/// This is a minimal weather input for step-based thermal simulations.
#[derive(Debug, Clone, Copy)]
pub struct OutdoorAirTemperatureC(pub f64);

/// Internal heat gains per zone (people + equipment + lighting), keyed by zone `UID` [W].
#[derive(Debug, Clone, Default)]
pub struct InternalGainsWPerZone {
    pub watts_by_zone_uid: HashMap<UID, f64>,
}

/// Total internal heat gains for the whole building [W].
///
/// When only a building-wide number is available, thermal simulations may distribute
/// this across zones by volume as a fallback.
#[derive(Debug, Clone, Copy)]
pub struct InternalGainsWTotal(pub f64);

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
