use anyhow::Result;

use crate::sim::coupling::{
    InternalGainsWPerZone, InternalGainsWTotal, OutdoorAirTemperatureC,
    ShortwaveTransmittedWPerZone,
};
use crate::sim::framework::{Bus, SimContext, SimModule};

use super::boundary::ThermalBoundaries;
use super::config::ThermalConfig;
use super::hvac::HvacIdealLoads;
use super::network::{MultiZoneAirModel, ThermalNetwork};

/// Configuration for step-based multi-zone thermal simulation.
#[derive(Debug, Clone)]
pub struct EnergyModuleConfig {
    /// Thermal envelope and model parameters.
    pub thermal: ThermalConfig,
    /// Ideal HVAC controller.
    pub hvac: HvacIdealLoads,
    /// Timestep in seconds.
    pub dt_s: f64,
    /// If true, run in steady-state mode (equivalent to zero thermal capacity).
    pub steady_state: bool,
}

impl Default for EnergyModuleConfig {
    fn default() -> Self {
        Self {
            thermal: ThermalConfig::new(),
            hvac: HvacIdealLoads::new(),
            dt_s: 3600.0,
            steady_state: false,
        }
    }
}

/// Step-based multi-zone thermal simulation module.
///
/// Inputs (via [`Bus`]):
/// - [`OutdoorAirTemperatureC`] (optional; falls back to `thermal.outdoor_temperature`)
/// - [`InternalGainsWPerZone`] or [`InternalGainsWTotal`] (optional; default 0)
/// - [`ShortwaveTransmittedWPerZone`] (optional; treated as solar shortwave gains)
///
/// Outputs (via [`Bus`]):
/// - [`MultiZoneStepResult`] for the latest step
pub struct EnergyModule {
    config: EnergyModuleConfig,
    model: Option<MultiZoneAirModel>,
    step_index: usize,
    zone_volumes_m3: Vec<f64>,
    total_volume_m3: f64,
    zone_uids: Vec<crate::UID>,
}

impl EnergyModule {
    pub fn new(config: EnergyModuleConfig) -> Self {
        Self {
            config,
            model: None,
            step_index: 0,
            zone_volumes_m3: vec![],
            total_volume_m3: 0.0,
            zone_uids: vec![],
        }
    }

    fn gains_by_zone(&self, bus: &Bus) -> Vec<f64> {
        let n = self.zone_uids.len();
        let mut gains = vec![0.0; n];

        let internal_by_zone = bus
            .get::<InternalGainsWPerZone>()
            .map(|g| &g.watts_by_zone_uid);

        let internal_total = bus.get::<InternalGainsWTotal>().map(|g| g.0).unwrap_or(0.0);

        let solar_by_zone = bus
            .get::<ShortwaveTransmittedWPerZone>()
            .map(|g| &g.watts_by_zone_uid);

        for i in 0..n {
            let uid = &self.zone_uids[i];

            let internal_i = if let Some(map) = internal_by_zone {
                map.get(uid).cloned().unwrap_or(0.0)
            } else if self.total_volume_m3 > 1e-14 {
                internal_total * (self.zone_volumes_m3[i] / self.total_volume_m3)
            } else {
                0.0
            };

            let solar_i = solar_by_zone
                .and_then(|map| map.get(uid).cloned())
                .unwrap_or(0.0);

            gains[i] = internal_i + solar_i;
        }

        gains
    }
}

impl SimModule for EnergyModule {
    fn name(&self) -> &'static str {
        "energy"
    }

    fn init(&mut self, ctx: &SimContext, _bus: &mut Bus) -> Result<()> {
        let boundaries = ThermalBoundaries::classify(ctx.building, ctx.surface_index);
        let network = ThermalNetwork::build(
            ctx.building,
            &self.config.thermal,
            ctx.surface_index,
            &boundaries,
        );

        let zones = ctx.building.zones();
        self.zone_volumes_m3 = zones.iter().map(|z| z.volume()).collect();
        self.total_volume_m3 = self.zone_volumes_m3.iter().sum();
        self.zone_uids = zones.iter().map(|z| z.uid.clone()).collect();

        let cap = if self.config.steady_state {
            0.0
        } else {
            self.config.thermal.thermal_capacity_j_per_m3_k
        };

        self.model = Some(MultiZoneAirModel::new(
            ctx.building,
            &network,
            self.config.thermal.infiltration_ach,
            cap,
            self.config.thermal.indoor_temperature,
        ));
        Ok(())
    }

    fn step(&mut self, _ctx: &SimContext, bus: &mut Bus) -> Result<()> {
        let outdoor_temp_c = bus
            .get::<OutdoorAirTemperatureC>()
            .map(|t| t.0)
            .unwrap_or(self.config.thermal.outdoor_temperature);

        let gains = self.gains_by_zone(bus);

        let Some(model) = self.model.as_mut() else {
            anyhow::bail!("EnergyModule not initialized");
        };

        let result = model.step(outdoor_temp_c, &gains, &self.config.hvac, self.config.dt_s)?;
        bus.put(result.clone());

        self.step_index += 1;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::energy::network::MultiZoneStepResult;
    use crate::sim::framework::SimContext;
    use crate::sim::index::SurfaceIndex;
    use crate::{Building, Solid, Zone};

    #[test]
    fn test_energy_module_consumes_shortwave_transmitted_per_zone() {
        // Two zones with no envelope loss and no coupling; shortwave in z0 should only heat z0.
        let s0 = Solid::from_box(1.0, 1.0, 1.0, None, "s0").unwrap();
        let s1 = Solid::from_box(1.0, 1.0, 1.0, Some((3.0, 0.0, 0.0)), "s1").unwrap();
        let z0 = Zone::new("z0", vec![s0]).unwrap();
        let z1 = Zone::new("z1", vec![s1]).unwrap();
        let building = Building::new("b", vec![z0, z1]).unwrap();

        let surface_index = SurfaceIndex::new(&building);
        let ctx = SimContext::new(&building, &surface_index);
        let mut bus = Bus::new();

        let mut thermal = ThermalConfig::new();
        thermal.default_u_value = 0.0;
        thermal.infiltration_ach = 0.0;
        thermal.indoor_temperature = 20.0;
        thermal.outdoor_temperature = 20.0;
        thermal.thermal_capacity_j_per_m3_k = 50_000.0;

        let hvac = HvacIdealLoads::with_setpoints(-1e9, 1e9);

        let mut module = EnergyModule::new(EnergyModuleConfig {
            thermal,
            hvac,
            dt_s: 3600.0,
            steady_state: false,
        });
        module.init(&ctx, &mut bus).unwrap();

        bus.put(OutdoorAirTemperatureC(20.0));

        let mut sw = ShortwaveTransmittedWPerZone::default();
        sw.watts_by_zone_uid
            .insert(module.zone_uids[0].clone(), 1000.0);
        bus.put(sw);

        module.step(&ctx, &mut bus).unwrap();

        let out = bus.get::<MultiZoneStepResult>().unwrap();
        assert!(out.zone_temperatures_c[0] > 20.0);
        assert!((out.zone_temperatures_c[1] - 20.0).abs() < 1e-6);
    }
}
