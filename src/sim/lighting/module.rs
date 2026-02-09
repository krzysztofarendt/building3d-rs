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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::framework::Pipeline;
    use crate::sim::index::SurfaceIndex;
    use crate::sim::lighting::result::LightingResult;
    use crate::sim::lighting::sources::PointLight;
    use crate::{Building, Point, Solid, Zone};

    #[test]
    fn test_lighting_module_runs_once_and_publishes_result() {
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();
        let index = SurfaceIndex::new(&building);
        let ctx = SimContext::new(&building, &index);

        let mut config = LightingConfig::new();
        config.num_rays = 200;
        config.max_bounces = 1;
        config
            .point_lights
            .push(PointLight::white(Point::new(1.0, 1.0, 1.5), 1000.0));

        let mut bus = Bus::new();
        let mut module = LightingModule::new(config);

        // Exercise lazy-init branch by calling step() without init().
        module.step(&ctx, &mut bus).unwrap();
        assert!(bus.get::<LightingResult>().is_some());

        // Second step should be a no-op.
        let hit_count_before = bus.get::<LightingResult>().unwrap().hit_count.len();
        module.step(&ctx, &mut bus).unwrap();
        let hit_count_after = bus.get::<LightingResult>().unwrap().hit_count.len();
        assert_eq!(hit_count_before, hit_count_after);
    }

    #[test]
    fn test_module_in_pipeline() {
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();
        let index = SurfaceIndex::new(&building);
        let ctx = SimContext::new(&building, &index);

        let mut config = LightingConfig::new();
        config.num_rays = 50;
        config.max_bounces = 1;

        let mut pipeline = Pipeline::new().with_module(LightingModule::new(config));
        let mut bus = Bus::new();
        pipeline.init(&ctx, &mut bus).unwrap();
        pipeline.step(&ctx, &mut bus).unwrap();

        assert!(bus.get::<LightingResult>().is_some());
    }
}
