use anyhow::Result;

use super::{Bus, SimContext};

/// A composable simulation module.
///
/// Modules can run standalone (acoustics-only, lighting-only) or as part of a
/// personalized multi-physics pipeline. Communication is done via the [`Bus`].
pub trait SimModule {
    /// Human-readable identifier for debugging / telemetry.
    fn name(&self) -> &'static str;

    /// Optional one-time initialization hook.
    fn init(&mut self, _ctx: &SimContext, _bus: &mut Bus) -> Result<()> {
        Ok(())
    }

    /// Advances the module by one step.
    fn step(&mut self, ctx: &SimContext, bus: &mut Bus) -> Result<()>;
}
