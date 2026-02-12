use std::sync::Arc;

use anyhow::Result;

use crate::sim::coupling::{OutdoorAirTemperatureC, OutdoorWindSpeedMPerS, WeatherHourIndex};
use crate::sim::framework::{Bus, SimContext, SimModule};

use super::weather::WeatherData;

/// Step-based weather publisher backed by EPW-derived [`WeatherData`].
///
/// This module intentionally publishes *minimal* coupling payloads (e.g. outdoor temperature
/// and the hour index). Other modules can share the same `Arc<WeatherData>` if they need
/// additional fields (DNI/DHI, timestamps, etc.).
#[derive(Clone)]
pub struct WeatherModuleConfig {
    pub weather: Arc<WeatherData>,
}

impl WeatherModuleConfig {
    pub fn new(weather: Arc<WeatherData>) -> Self {
        Self { weather }
    }
}

pub struct WeatherModule {
    config: WeatherModuleConfig,
    hour_index: usize,
}

impl WeatherModule {
    pub fn new(config: WeatherModuleConfig) -> Self {
        Self {
            config,
            hour_index: 0,
        }
    }

    fn publish(&self, bus: &mut Bus, hour_index: usize) -> Result<()> {
        anyhow::ensure!(
            hour_index < self.config.weather.records.len(),
            "WeatherModule: hour index {hour_index} out of range (len={})",
            self.config.weather.records.len()
        );
        let record = &self.config.weather.records[hour_index];
        bus.put(WeatherHourIndex(hour_index));
        bus.put(OutdoorAirTemperatureC(record.dry_bulb_temperature));
        bus.put(OutdoorWindSpeedMPerS(record.wind_speed));
        Ok(())
    }
}

impl SimModule for WeatherModule {
    fn name(&self) -> &'static str {
        "weather"
    }

    fn init(&mut self, _ctx: &SimContext, bus: &mut Bus) -> Result<()> {
        // Make the dataset discoverable to downstream modules.
        bus.put(self.config.weather.clone());
        Ok(())
    }

    fn step(&mut self, _ctx: &SimContext, bus: &mut Bus) -> Result<()> {
        self.publish(bus, self.hour_index)?;
        self.hour_index += 1;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::index::SurfaceIndex;
    use crate::{Building, Solid, Zone};

    #[test]
    fn test_weather_module_publishes_temperature_and_index() -> Result<()> {
        let s = Solid::from_box(1.0, 1.0, 1.0, None, "s").unwrap();
        let z = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![z]).unwrap();
        let index = SurfaceIndex::new(&building);
        let ctx = SimContext::new(&building, &index);

        let weather = Arc::new(WeatherData::synthetic("X", 0.0, 0.0, 10.0, 0.0));
        let mut module = WeatherModule::new(WeatherModuleConfig::new(weather));
        let mut bus = Bus::new();

        module.init(&ctx, &mut bus)?;
        module.step(&ctx, &mut bus)?;
        let idx0 = bus.get::<WeatherHourIndex>().unwrap().0;
        assert_eq!(idx0, 0);

        module.step(&ctx, &mut bus)?;
        let idx1 = bus.get::<WeatherHourIndex>().unwrap().0;
        assert_eq!(idx1, 1);

        let t = bus.get::<OutdoorAirTemperatureC>().unwrap().0;
        assert!(t.is_finite());
        Ok(())
    }
}
