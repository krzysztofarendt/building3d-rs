use anyhow::Result;

use crate::sim::framework::{Bus, SimContext, SimModule};

use super::config::LightingConfig;
use super::simulation::LightingSimulation;

/// Pipeline wrapper for running lighting simulations in a composed workflow.
///
/// On the first `step()` call, this module:
/// - runs the lighting simulation (forward tracer),
/// - publishes `LightingResult`.
///
/// Cross-domain coupling payloads (e.g. shortwave solar gains) should be produced by a
/// dedicated module to keep “single source of truth” behavior explicit in composed
/// pipelines.
pub struct LightingModule {
    config: LightingConfig,
    has_run: bool,
    simulation: Option<LightingSimulation>,
}

impl LightingModule {
    pub fn new(config: LightingConfig) -> Self {
        Self {
            config,
            has_run: false,
            simulation: None,
        }
    }
}

impl SimModule for LightingModule {
    fn name(&self) -> &'static str {
        "lighting"
    }

    fn init(&mut self, ctx: &SimContext, _bus: &mut Bus) -> Result<()> {
        self.simulation = Some(LightingSimulation::new(ctx.building, self.config.clone())?);
        Ok(())
    }

    fn step(&mut self, ctx: &SimContext, bus: &mut Bus) -> Result<()> {
        if self.has_run {
            return Ok(());
        }

        let Some(sim) = self.simulation.as_ref() else {
            self.simulation = Some(LightingSimulation::new(ctx.building, self.config.clone())?);
            return self.step(ctx, bus);
        };

        let result = sim.run();
        bus.put(result);

        self.has_run = true;
        Ok(())
    }
}
