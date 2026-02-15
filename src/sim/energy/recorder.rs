use anyhow::Result;

use crate::sim::coupling::WeatherHourIndex;
use crate::sim::framework::{Bus, SimContext, SimModule};

use super::network::MultiZoneStepResult;
use super::simulation::{AnnualResult, MultiZoneAnnualResult};
use super::weather::WeatherData;

/// Mutable, Bus-stored recording buffer for multi-zone thermal simulations.
///
/// The intended workflow is:
/// 1) `MultiZoneRecorderModule` initializes this on the Bus,
/// 2) the recorder appends one sample per step,
/// 3) the caller takes the data from the Bus and finalizes it into a
///    [`MultiZoneAnnualResult`].
#[derive(Debug, Default)]
pub struct MultiZoneRecorderData {
    zone_uids: Vec<crate::UID>,
    zone_names: Vec<String>,
    hourly_zone_temperatures_c: Vec<Vec<f64>>,
    hourly_zone_heating_w: Vec<Vec<f64>>,
    hourly_zone_cooling_w: Vec<Vec<f64>>,
    hourly_heating_w: Vec<f64>,
    hourly_cooling_w: Vec<f64>,
    // Month-of-year per step (1-12) when available.
    month_by_step: Vec<u8>,
}

impl MultiZoneRecorderData {
    fn ensure_zones_initialized(&mut self, step: &MultiZoneStepResult) {
        if !self.zone_uids.is_empty() {
            return;
        }

        self.zone_uids = step.zone_uids.clone();
        self.zone_names = step.zone_names.clone();

        let n = self.zone_uids.len();
        self.hourly_zone_temperatures_c = vec![Vec::new(); n];
        self.hourly_zone_heating_w = vec![Vec::new(); n];
        self.hourly_zone_cooling_w = vec![Vec::new(); n];
    }

    fn push_step(&mut self, step: &MultiZoneStepResult, month: Option<u8>) -> anyhow::Result<()> {
        self.ensure_zones_initialized(step);

        let n = self.zone_uids.len();
        anyhow::ensure!(
            step.zone_temperatures_c.len() == n
                && step.zone_heating_w.len() == n
                && step.zone_cooling_w.len() == n,
            "MultiZoneRecorderData::push_step: inconsistent zone count \
             (expected {n}, got temps={}, heating={}, cooling={})",
            step.zone_temperatures_c.len(),
            step.zone_heating_w.len(),
            step.zone_cooling_w.len(),
        );

        let hour_heating: f64 = step.zone_heating_w.iter().sum();
        let hour_cooling: f64 = step.zone_cooling_w.iter().sum();
        self.hourly_heating_w.push(hour_heating);
        self.hourly_cooling_w.push(hour_cooling);

        for z in 0..n {
            self.hourly_zone_temperatures_c[z].push(step.zone_temperatures_c[z]);
            self.hourly_zone_heating_w[z].push(step.zone_heating_w[z]);
            self.hourly_zone_cooling_w[z].push(step.zone_cooling_w[z]);
        }

        let m = month.unwrap_or(1).clamp(1, 12);
        self.month_by_step.push(m);

        Ok(())
    }

    pub fn finalize(self, dt_s: f64) -> MultiZoneAnnualResult {
        let dt_h = (dt_s / 3600.0).max(0.0);
        let to_kwh = dt_h / 1000.0;

        let num_steps = self.hourly_heating_w.len();
        let mut annual_heating_wh_equiv = 0.0;
        let mut annual_cooling_wh_equiv = 0.0;
        let mut peak_heating = 0.0_f64;
        let mut peak_cooling = 0.0_f64;
        let mut monthly_heating = [0.0; 12];
        let mut monthly_cooling = [0.0; 12];

        for i in 0..num_steps {
            let h = self.hourly_heating_w[i];
            let c = self.hourly_cooling_w[i];
            annual_heating_wh_equiv += h;
            annual_cooling_wh_equiv += c;
            peak_heating = peak_heating.max(h);
            peak_cooling = peak_cooling.max(c);

            let month = self.month_by_step.get(i).cloned().unwrap_or(1);
            let month_idx = (month as usize) - 1;
            monthly_heating[month_idx] += h;
            monthly_cooling[month_idx] += c;
        }

        let annual = AnnualResult {
            hourly_heating: self.hourly_heating_w,
            hourly_cooling: self.hourly_cooling_w,
            annual_heating_kwh: annual_heating_wh_equiv * to_kwh,
            annual_cooling_kwh: annual_cooling_wh_equiv * to_kwh,
            peak_heating,
            peak_cooling,
            monthly_heating_kwh: monthly_heating.map(|v| v * to_kwh),
            monthly_cooling_kwh: monthly_cooling.map(|v| v * to_kwh),
            min_zone_temp_c: 0.0,
            max_zone_temp_c: 0.0,
            hourly_zone_temp_c: Vec::new(),
            num_zones: 1,
            per_zone_hourly_temp_c: Vec::new(),
            per_zone_min_temp_c: Vec::new(),
            per_zone_max_temp_c: Vec::new(),
        };

        MultiZoneAnnualResult {
            zone_uids: self.zone_uids,
            zone_names: self.zone_names,
            hourly_zone_temperatures_c: self.hourly_zone_temperatures_c,
            hourly_zone_heating_w: self.hourly_zone_heating_w,
            hourly_zone_cooling_w: self.hourly_zone_cooling_w,
            annual,
        }
    }
}

/// Records `MultiZoneStepResult` samples into a [`MultiZoneRecorderData`] stored on the Bus.
pub struct MultiZoneRecorderModule {
    dt_s: f64,
}

impl MultiZoneRecorderModule {
    pub fn new(dt_s: f64) -> Self {
        Self { dt_s }
    }

    /// Removes the recorder data from the Bus and converts it into a finalized result.
    pub fn take_result(bus: &mut Bus, dt_s: f64) -> Result<MultiZoneAnnualResult> {
        let Some(data) = bus.take::<MultiZoneRecorderData>() else {
            anyhow::bail!("MultiZoneRecorderData not found on Bus");
        };
        Ok(data.finalize(dt_s))
    }
}

impl SimModule for MultiZoneRecorderModule {
    fn name(&self) -> &'static str {
        "energy_recorder"
    }

    fn init(&mut self, _ctx: &SimContext, bus: &mut Bus) -> Result<()> {
        if bus.get::<MultiZoneRecorderData>().is_none() {
            bus.put(MultiZoneRecorderData::default());
        }
        Ok(())
    }

    fn step(&mut self, _ctx: &SimContext, bus: &mut Bus) -> Result<()> {
        let Some(step) = bus.get::<MultiZoneStepResult>().cloned() else {
            anyhow::bail!("MultiZoneRecorderModule requires MultiZoneStepResult on the Bus");
        };

        let month = month_from_weather(bus);

        let Some(data) = bus.get_mut::<MultiZoneRecorderData>() else {
            anyhow::bail!("MultiZoneRecorderData not initialized on Bus");
        };
        data.push_step(&step, month)?;

        // Keep dt_s here for future: might store it in RecorderData if needed.
        let _ = self.dt_s;
        Ok(())
    }
}

fn month_from_weather(bus: &Bus) -> Option<u8> {
    let hour_index = bus.get::<WeatherHourIndex>()?.0;
    let weather = bus.get::<std::sync::Arc<WeatherData>>()?;
    weather.records.get(hour_index).map(|r| r.month)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::coupling::OutdoorAirTemperatureC;
    use crate::sim::energy::config::ThermalConfig;
    use crate::sim::energy::hvac::HvacIdealLoads;
    use crate::sim::energy::module::{EnergyModule, EnergyModuleConfig};
    use crate::sim::energy::weather_module::{WeatherModule, WeatherModuleConfig};
    use crate::sim::framework::Pipeline;
    use crate::sim::index::SurfaceIndex;
    use crate::{Building, Solid, Zone};
    use std::sync::Arc;

    #[test]
    fn test_multizone_recorder_pipeline_smoke() -> Result<()> {
        let s = Solid::from_box(3.0, 3.0, 3.0, None, "s").unwrap();
        let z = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![z]).unwrap();
        let index = SurfaceIndex::new(&building);
        let ctx = SimContext::new(&building, &index);

        let mut thermal = ThermalConfig::new();
        thermal.default_u_value = 0.0;
        // Ensure the thermal system is well-posed even with `default_u_value = 0.0`.
        thermal.infiltration_ach = 0.1;
        thermal.indoor_temperature = 20.0;
        thermal.outdoor_temperature = 0.0;
        thermal.thermal_capacity_j_per_m3_k = 0.0; // steady

        let hvac = HvacIdealLoads::with_setpoints(20.0, 20.0);

        let weather = Arc::new(WeatherData::synthetic("X", 52.0, 13.0, 0.0, 0.0));

        let mut bus = Bus::new();
        // Provide a fallback temperature too (should be overwritten by WeatherModule each step).
        bus.put(OutdoorAirTemperatureC(0.0));

        let mut pipeline = Pipeline::new()
            .with_module(WeatherModule::new(WeatherModuleConfig::new(weather)))
            .with_module(EnergyModule::new(EnergyModuleConfig {
                thermal,
                hvac,
                dt_s: 3600.0,
                steady_state: true,
                model_kind: crate::sim::energy::module::EnergyModelKind::AirOnly,
                material_library: None,
                default_envelope_capacity_j_per_m2_k: 0.0,
            }))
            .with_module(MultiZoneRecorderModule::new(3600.0));

        pipeline.init(&ctx, &mut bus)?;
        for _ in 0..3 {
            pipeline.step(&ctx, &mut bus)?;
        }

        let result = MultiZoneRecorderModule::take_result(&mut bus, 3600.0)?;
        assert_eq!(result.annual.hourly_heating.len(), 3);
        Ok(())
    }
}
