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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::index::SurfaceIndex;
    use crate::{Building, Solid, Zone};

    struct Dummy;
    impl SimModule for Dummy {
        fn name(&self) -> &'static str {
            "dummy"
        }

        fn step(&mut self, _ctx: &SimContext, bus: &mut Bus) -> Result<()> {
            bus.put(42_u32);
            Ok(())
        }
    }

    #[test]
    fn test_default_init_is_ok() {
        let s0 = Solid::from_box(1.0, 1.0, 1.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();
        let index = SurfaceIndex::new(&building);
        let ctx = SimContext::new(&building, &index);

        let mut bus = Bus::new();
        let mut m = Dummy;

        m.init(&ctx, &mut bus).unwrap();
        m.step(&ctx, &mut bus).unwrap();
        assert_eq!(bus.get::<u32>(), Some(&42));
    }
}
