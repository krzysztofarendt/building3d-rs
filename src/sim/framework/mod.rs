//! Generic simulation framework.
//!
//! This module is intentionally domain-agnostic: it provides a small runtime
//! for composing simulation modules (energy, acoustics, lighting, etc.) without
//! hardcoding domain-specific metadata onto geometry types.

pub mod bus;
pub mod context;
pub mod module;
pub mod pipeline;

pub use bus::Bus;
pub use context::SimContext;
pub use module::SimModule;
pub use pipeline::Pipeline;
