use anyhow::Result;

use crate::sim::coupling::{
    InternalGainsWPerZone, InternalGainsWTotal, OutdoorAirTemperatureC,
    ShortwaveAbsorbedWPerPolygon, ShortwaveTransmittedWPerZone,
};
use crate::sim::framework::{Bus, SimContext, SimModule};

use super::boundary::ThermalBoundaries;
use super::config::ThermalConfig;
use super::hvac::HvacIdealLoads;
use super::network::{MultiZoneAirModel, MultiZoneEnvelopeRcModel, ThermalNetwork};
use crate::sim::materials::MaterialLibrary;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnergyModelKind {
    /// Zone air nodes only (current baseline).
    AirOnly,
    /// Zone air nodes + one aggregated 2R1C envelope node per zone.
    EnvelopeRc2R1C,
}

impl Default for EnergyModelKind {
    fn default() -> Self {
        Self::AirOnly
    }
}

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
    /// Selects the internal thermal network model.
    pub model_kind: EnergyModelKind,
    /// Optional material library (legacy convenience wiring).
    ///
    /// If `thermal.material_library` is not set, this is copied into it during `init()`.
    /// Prefer setting `thermal.material_library` directly so U-values and capacities are
    /// resolved consistently from one source of truth.
    pub material_library: Option<MaterialLibrary>,
    /// Fallback envelope capacity per exterior area [J/(m²·K)] for surfaces that
    /// do not have `ThermalMaterial.thermal_capacity` in the material library.
    pub default_envelope_capacity_j_per_m2_k: f64,
}

impl Default for EnergyModuleConfig {
    fn default() -> Self {
        Self {
            thermal: ThermalConfig::new(),
            hvac: HvacIdealLoads::new(),
            dt_s: 3600.0,
            steady_state: false,
            model_kind: EnergyModelKind::AirOnly,
            material_library: None,
            default_envelope_capacity_j_per_m2_k: 0.0,
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
    model: Option<EnergyModel>,
    step_index: usize,
    zone_volumes_m3: Vec<f64>,
    total_volume_m3: f64,
    zone_uids: Vec<crate::UID>,
    polygon_uid_to_zone_idx: std::collections::HashMap<crate::UID, usize>,
    boundaries: Option<ThermalBoundaries>,
}

enum EnergyModel {
    Air(MultiZoneAirModel),
    EnvelopeRc(MultiZoneEnvelopeRcModel),
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
            polygon_uid_to_zone_idx: std::collections::HashMap::new(),
            boundaries: None,
        }
    }

    fn gains_by_zone_split(&self, bus: &Bus) -> (Vec<f64>, Vec<f64>) {
        let n = self.zone_uids.len();
        let mut air_gains = vec![0.0; n];
        let mut env_gains = vec![0.0; n];

        let internal_by_zone = bus
            .get::<InternalGainsWPerZone>()
            .map(|g| &g.watts_by_zone_uid);

        let internal_total = bus.get::<InternalGainsWTotal>().map(|g| g.0).unwrap_or(0.0);

        let solar_by_zone = bus
            .get::<ShortwaveTransmittedWPerZone>()
            .map(|g| &g.watts_by_zone_uid);

        // Apply internal + transmitted shortwave to the air node.
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

            air_gains[i] += internal_i + solar_i;
        }

        // Per-polygon absorbed shortwave: split between envelope (exterior) and air (non-exterior).
        let is_exterior = |uid: &crate::UID| {
            self.boundaries
                .as_ref()
                .map(|b| b.is_exterior(uid))
                .unwrap_or(false)
        };

        if let Some(abs) = bus.get::<ShortwaveAbsorbedWPerPolygon>() {
            for (polygon_uid, w) in &abs.watts_by_polygon_uid {
                if *w == 0.0 {
                    continue;
                }
                if let Some(&zone_idx) = self.polygon_uid_to_zone_idx.get(polygon_uid) {
                    if is_exterior(polygon_uid) {
                        env_gains[zone_idx] += *w;
                    } else {
                        air_gains[zone_idx] += *w;
                    }
                }
            }
        }

        (air_gains, env_gains)
    }
}

impl SimModule for EnergyModule {
    fn name(&self) -> &'static str {
        "energy"
    }

    fn init(&mut self, ctx: &SimContext, _bus: &mut Bus) -> Result<()> {
        // Backward compatible wiring: allow the module-level material library field to
        // populate `ThermalConfig.material_library` when the latter is unset.
        if self.config.thermal.material_library.is_none() {
            self.config.thermal.material_library = self.config.material_library.clone();
        }

        let boundaries = ThermalBoundaries::classify(ctx.building, ctx.surface_index);
        self.boundaries = Some(boundaries.clone());
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
        self.polygon_uid_to_zone_idx.clear();
        let uid_to_idx: std::collections::HashMap<&str, usize> = self
            .zone_uids
            .iter()
            .enumerate()
            .map(|(i, uid)| (uid.as_str(), i))
            .collect();
        for s in &ctx.surface_index.surfaces {
            if let Some(&idx) = uid_to_idx.get(s.zone_uid.as_str()) {
                self.polygon_uid_to_zone_idx
                    .insert(s.polygon_uid.clone(), idx);
            }
        }

        let cap = if self.config.steady_state {
            0.0
        } else {
            self.config.thermal.thermal_capacity_j_per_m3_k
        };

        self.model = Some(match self.config.model_kind {
            EnergyModelKind::AirOnly => EnergyModel::Air(MultiZoneAirModel::new(
                ctx.building,
                &network,
                self.config.thermal.infiltration_ach,
                cap,
                self.config.thermal.indoor_temperature,
            )),
            EnergyModelKind::EnvelopeRc2R1C => {
                let default_env_cap = if self.config.steady_state {
                    0.0
                } else {
                    self.config.default_envelope_capacity_j_per_m2_k
                };

                EnergyModel::EnvelopeRc(MultiZoneEnvelopeRcModel::new_with_thermal_config(
                    ctx.building,
                    &network,
                    ctx.surface_index,
                    &boundaries,
                    &self.config.thermal,
                    self.config.thermal.infiltration_ach,
                    cap,
                    default_env_cap,
                    self.config.thermal.indoor_temperature,
                ))
            }
        });
        Ok(())
    }

    fn step(&mut self, _ctx: &SimContext, bus: &mut Bus) -> Result<()> {
        let outdoor_temp_c = bus
            .get::<OutdoorAirTemperatureC>()
            .map(|t| t.0)
            .unwrap_or(self.config.thermal.outdoor_temperature);

        let (air_gains, env_gains) = self.gains_by_zone_split(bus);

        let Some(model) = self.model.as_mut() else {
            anyhow::bail!("EnergyModule not initialized");
        };

        let result = match model {
            EnergyModel::Air(m) => {
                let mut gains = air_gains;
                for (i, v) in env_gains.into_iter().enumerate() {
                    gains[i] += v;
                }
                m.step(outdoor_temp_c, &gains, &self.config.hvac, self.config.dt_s)?
            }
            EnergyModel::EnvelopeRc(m) => m.step(
                outdoor_temp_c,
                &air_gains,
                &env_gains,
                &self.config.hvac,
                self.config.dt_s,
            )?,
        };
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
    fn test_energy_module_name() {
        let module = EnergyModule::new(EnergyModuleConfig::default());
        assert_eq!(module.name(), "energy");
    }

    #[test]
    fn test_energy_module_config_default_trait() {
        let cfg: EnergyModuleConfig = Default::default();
        assert_eq!(cfg.model_kind, EnergyModelKind::AirOnly);
        assert!((cfg.dt_s - 3600.0).abs() < 1e-12);
    }

    #[test]
    fn test_energy_module_step_without_init_errors() {
        let s = Solid::from_box(1.0, 1.0, 1.0, None, "s").unwrap();
        let z = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![z]).unwrap();
        let surface_index = SurfaceIndex::new(&building);
        let ctx = SimContext::new(&building, &surface_index);
        let mut bus = Bus::new();

        let module_cfg: EnergyModuleConfig = Default::default();
        let mut module = EnergyModule::new(module_cfg);
        let err = module.step(&ctx, &mut bus).unwrap_err();
        assert!(err.to_string().contains("not initialized"));
    }

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
            model_kind: EnergyModelKind::AirOnly,
            material_library: None,
            default_envelope_capacity_j_per_m2_k: 0.0,
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

    #[test]
    fn test_gains_by_zone_split_routes_absorbed_shortwave_by_exteriority() {
        // Two adjacent solids in the same zone => shared face is not exterior.
        let s0 = Solid::from_box(1.0, 1.0, 1.0, None, "s0").unwrap();
        let s1 = Solid::from_box(1.0, 1.0, 1.0, Some((1.0, 0.0, 0.0)), "s1").unwrap();
        let z = Zone::new("z", vec![s0, s1]).unwrap();
        let building = Building::new("b", vec![z]).unwrap();

        let surface_index = SurfaceIndex::new(&building);
        let ctx = SimContext::new(&building, &surface_index);
        let mut bus = Bus::new();

        let mut thermal = ThermalConfig::new();
        thermal.default_u_value = 0.0;
        thermal.infiltration_ach = 0.0;
        thermal.indoor_temperature = 20.0;
        thermal.outdoor_temperature = 20.0;

        let hvac = HvacIdealLoads::with_setpoints(-1e9, 1e9);
        let mut module = EnergyModule::new(EnergyModuleConfig {
            thermal,
            hvac,
            dt_s: 3600.0,
            steady_state: false,
            model_kind: EnergyModelKind::AirOnly,
            material_library: None,
            default_envelope_capacity_j_per_m2_k: 0.0,
        });
        module.init(&ctx, &mut bus).unwrap();

        let interior = surface_index
            .polygon_uid_by_path("z/s0/wall_1/poly_1")
            .unwrap()
            .clone();
        let exterior = surface_index
            .polygon_uid_by_path("z/s0/floor/floor")
            .unwrap()
            .clone();
        let boundaries = module.boundaries.as_ref().unwrap();
        assert!(!boundaries.is_exterior(&interior));
        assert!(boundaries.is_exterior(&exterior));

        let mut internal = InternalGainsWPerZone::default();
        internal
            .watts_by_zone_uid
            .insert(module.zone_uids[0].clone(), 50.0);
        bus.put(internal);

        let mut transmitted = ShortwaveTransmittedWPerZone::default();
        transmitted
            .watts_by_zone_uid
            .insert(module.zone_uids[0].clone(), 25.0);
        bus.put(transmitted);

        let mut absorbed = ShortwaveAbsorbedWPerPolygon::default();
        absorbed.watts_by_polygon_uid.insert(interior, 100.0);
        absorbed.watts_by_polygon_uid.insert(exterior, 200.0);
        bus.put(absorbed);

        let (air_gains, env_gains) = module.gains_by_zone_split(&bus);
        assert_eq!(air_gains.len(), 1);
        assert_eq!(env_gains.len(), 1);
        assert!((air_gains[0] - (50.0 + 25.0 + 100.0)).abs() < 1e-12);
        assert!((env_gains[0] - 200.0).abs() < 1e-12);
    }

    #[test]
    fn test_gains_by_zone_split_internal_total_with_zero_total_volume_returns_zero() {
        // Cover the defensive branch when `total_volume_m3` is effectively zero.
        let s0 = Solid::from_box(1.0, 1.0, 1.0, None, "s0").unwrap();
        let s1 = Solid::from_box(1.0, 1.0, 1.0, Some((3.0, 0.0, 0.0)), "s1").unwrap();
        let z0 = Zone::new("z0", vec![s0]).unwrap();
        let z1 = Zone::new("z1", vec![s1]).unwrap();
        let building = Building::new("b", vec![z0, z1]).unwrap();

        let surface_index = SurfaceIndex::new(&building);
        let ctx = SimContext::new(&building, &surface_index);
        let mut bus = Bus::new();

        let mut module = EnergyModule::new(EnergyModuleConfig::default());
        module.init(&ctx, &mut bus).unwrap();
        assert!(module.total_volume_m3 > 1e-6);

        // Force the "else { 0.0 }" internal gains branch.
        module.total_volume_m3 = 0.0;
        bus.put(InternalGainsWTotal(123.0));

        let (air_gains, env_gains) = module.gains_by_zone_split(&bus);
        assert_eq!(air_gains, vec![0.0, 0.0]);
        assert_eq!(env_gains, vec![0.0, 0.0]);
    }

    #[test]
    fn test_energy_module_envelope_rc_steady_state_init_path() {
        let s = Solid::from_box(2.0, 2.0, 2.0, None, "s").unwrap();
        let z = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![z]).unwrap();

        let surface_index = SurfaceIndex::new(&building);
        let ctx = SimContext::new(&building, &surface_index);
        let mut bus = Bus::new();

        let cfg = EnergyModuleConfig {
            model_kind: EnergyModelKind::EnvelopeRc2R1C,
            steady_state: true,
            ..EnergyModuleConfig::default()
        };
        let mut module = EnergyModule::new(cfg);
        module.init(&ctx, &mut bus).unwrap();

        match module.model.as_ref().unwrap() {
            EnergyModel::EnvelopeRc(_) => {}
            _ => panic!("expected EnvelopeRc model"),
        }
    }

    #[test]
    fn test_energy_module_envelope_rc_has_thermal_lag() {
        // No gains, HVAC off, outdoor drops: air-only steady-state should jump to outdoor,
        // while RC envelope retains heat for at least one step.
        let s = Solid::from_box(3.0, 3.0, 3.0, None, "s").unwrap();
        let z = Zone::new("z", vec![s]).unwrap();
        let building = Building::new("b", vec![z]).unwrap();

        let surface_index = SurfaceIndex::new(&building);
        let ctx = SimContext::new(&building, &surface_index);

        let mut thermal = ThermalConfig::new();
        thermal.default_u_value = 1.0;
        thermal.infiltration_ach = 0.0;
        thermal.thermal_capacity_j_per_m3_k = 0.0; // isolate envelope capacity effect
        thermal.indoor_temperature = 20.0;
        thermal.outdoor_temperature = 20.0;

        let hvac = HvacIdealLoads::with_setpoints(-1e9, 1e9);

        let mut bus_air = Bus::new();
        let mut air_only = EnergyModule::new(EnergyModuleConfig {
            thermal: thermal.clone(),
            hvac: hvac.clone(),
            dt_s: 3600.0,
            steady_state: true,
            model_kind: EnergyModelKind::AirOnly,
            material_library: None,
            default_envelope_capacity_j_per_m2_k: 0.0,
        });
        air_only.init(&ctx, &mut bus_air).unwrap();
        bus_air.put(OutdoorAirTemperatureC(0.0));
        air_only.step(&ctx, &mut bus_air).unwrap();
        let out_air = bus_air.get::<MultiZoneStepResult>().unwrap();
        let t_air = out_air.zone_temperatures_c[0];

        let mut bus_rc = Bus::new();
        let mut rc = EnergyModule::new(EnergyModuleConfig {
            thermal,
            hvac,
            dt_s: 3600.0,
            steady_state: false,
            model_kind: EnergyModelKind::EnvelopeRc2R1C,
            material_library: None,
            default_envelope_capacity_j_per_m2_k: 200_000.0,
        });
        rc.init(&ctx, &mut bus_rc).unwrap();
        bus_rc.put(OutdoorAirTemperatureC(0.0));
        rc.step(&ctx, &mut bus_rc).unwrap();
        let out_rc = bus_rc.get::<MultiZoneStepResult>().unwrap();
        let t_rc = out_rc.zone_temperatures_c[0];

        assert!(t_air < 1e-6, "air-only steady-state should go to outdoor");
        assert!(t_rc > 0.1, "RC envelope should retain heat at first step");
    }

    #[test]
    fn test_energy_module_consumes_shortwave_absorbed_per_polygon() {
        // Two zones with no envelope loss and no coupling; absorbed shortwave on a polygon
        // in z0 should only heat z0 (currently treated as a zone gain).
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
            model_kind: EnergyModelKind::AirOnly,
            material_library: None,
            default_envelope_capacity_j_per_m2_k: 0.0,
        });
        module.init(&ctx, &mut bus).unwrap();

        bus.put(OutdoorAirTemperatureC(20.0));

        let z0_uid = ctx.building.zones().first().unwrap().uid.clone();
        let polygon_uid_z0 = surface_index
            .surfaces
            .iter()
            .find(|s| s.zone_uid == z0_uid)
            .unwrap()
            .polygon_uid
            .clone();

        let mut abs = ShortwaveAbsorbedWPerPolygon::default();
        abs.watts_by_polygon_uid.insert(polygon_uid_z0, 1000.0);
        bus.put(abs);

        module.step(&ctx, &mut bus).unwrap();
        let out = bus.get::<MultiZoneStepResult>().unwrap();
        assert!(out.zone_temperatures_c[0] > 20.0);
        assert!((out.zone_temperatures_c[1] - 20.0).abs() < 1e-6);
    }
}
