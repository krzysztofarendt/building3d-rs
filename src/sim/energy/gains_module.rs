use std::sync::Arc;

use anyhow::Result;

use crate::sim::coupling::{InternalGainsWTotal, WeatherHourIndex};
use crate::sim::framework::{Bus, SimContext, SimModule};

use super::schedule::InternalGainsProfile;

/// Publishes time-varying internal gains based on an [`InternalGainsProfile`].
///
/// This module consumes [`WeatherHourIndex`] as the shared step counter.
#[derive(Clone)]
pub struct InternalGainsModuleConfig {
    pub profile: Arc<InternalGainsProfile>,
    /// Optional multiplier for scenario scaling (e.g. partial occupancy).
    pub multiplier: f64,
}

impl InternalGainsModuleConfig {
    pub fn new(profile: Arc<InternalGainsProfile>) -> Self {
        Self {
            profile,
            multiplier: 1.0,
        }
    }
}

pub struct InternalGainsModule {
    config: InternalGainsModuleConfig,
}

impl InternalGainsModule {
    pub fn new(config: InternalGainsModuleConfig) -> Self {
        Self { config }
    }
}

impl SimModule for InternalGainsModule {
    fn name(&self) -> &'static str {
        "internal_gains"
    }

    fn step(&mut self, _ctx: &SimContext, bus: &mut Bus) -> Result<()> {
        let Some(hour_idx) = bus.get::<WeatherHourIndex>().map(|i| i.0) else {
            anyhow::bail!("InternalGainsModule requires WeatherHourIndex on the Bus");
        };
        let gains = self.config.profile.gains_at(hour_idx) * self.config.multiplier;
        bus.put(InternalGainsWTotal(gains));
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::index::SurfaceIndex;
    use crate::{Building, Solid, Zone};

    #[test]
    fn test_internal_gains_module_publishes_total() -> Result<()> {
        let s = Solid::from_box(1.0, 1.0, 1.0, None, "s").unwrap();
        let z = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![z]).unwrap();
        let index = SurfaceIndex::new(&building);
        let ctx = SimContext::new(&building, &index);

        let profile = Arc::new(InternalGainsProfile::office(100.0));
        let mut module = InternalGainsModule::new(InternalGainsModuleConfig::new(profile));

        let mut bus = Bus::new();
        bus.put(WeatherHourIndex(10)); // Mon 10am: occupied

        module.step(&ctx, &mut bus)?;
        let g = bus.get::<InternalGainsWTotal>().unwrap().0;
        assert!(g > 1000.0);
        Ok(())
    }
}
