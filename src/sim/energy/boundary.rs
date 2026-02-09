//! Thermal boundary classification for energy simulations.
//!
//! This module historically exposed `ThermalBoundaries` as a thermal-specific overlay.
//! As part of the cross-domain unification plan (see `PHYSICS.md`), the underlying
//! implementation is now shared across domains in `sim::surfaces`.

pub use crate::sim::surfaces::{
    SurfaceInterface as ThermalInterface, SurfaceKind as ThermalSurfaceKind,
    SurfaceSemantics as ThermalBoundaries,
};
