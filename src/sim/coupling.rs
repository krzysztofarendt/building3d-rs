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

/// Index of the current weather timestep (hour-of-year) in the active weather dataset.
///
/// This is typically an EPW record index (0..8759), but may differ if the dataset has
/// leap-year hours (8784) or if a shorter custom series is used.
#[derive(Debug, Clone, Copy)]
pub struct WeatherHourIndex(pub usize);

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

/// Outdoor wind speed [m/s].
///
/// Used by dynamic exterior convection models (DOE-2) to compute wind-driven
/// forced convection coefficients.
#[derive(Debug, Clone, Copy)]
pub struct OutdoorWindSpeedMPerS(pub f64);

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

/// Shortwave radiant power transmitted and distributed to interior polygon surfaces [W].
///
/// This is a refinement of [`ShortwaveTransmittedWPerZone`]: instead of a bulk
/// per-zone number, the power is allocated to specific interior surfaces (e.g.
/// beam solar concentrated on floors, diffuse spread over all interiors).
#[derive(Debug, Clone, Default)]
pub struct ShortwaveTransmittedWPerPolygon {
    pub watts_by_polygon_uid: HashMap<UID, f64>,
}
