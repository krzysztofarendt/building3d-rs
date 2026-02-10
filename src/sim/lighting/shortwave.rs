use anyhow::Result;

use std::collections::HashMap;
use std::sync::Arc;

use crate::sim::coupling::OutdoorAirTemperatureC;
use crate::sim::coupling::WeatherHourIndex;
use crate::sim::coupling::{ShortwaveAbsorbedWPerPolygon, ShortwaveTransmittedWPerZone};
use crate::sim::energy::boundary::ThermalBoundaries;
use crate::sim::energy::solar_bridge::{
    SolarGainConfig, SolarHourParams, compute_solar_gains_per_zone_with_materials,
};
use crate::sim::energy::weather::WeatherData;
use crate::sim::engine::FlatScene;
use crate::sim::framework::{Bus, SimContext, SimModule};
use crate::sim::lighting::result::LightingResult;
use crate::sim::lighting::solar::SolarPosition;
use crate::sim::materials::{MaterialLibrary, OpticalMaterial};
use crate::{Point, Vector};

use super::sources::Rgb;

const DEFAULT_OPAQUE_ABSORPTANCE: f64 = 0.7;

fn opaque_absorptance_for_path(path: &str, material_library: Option<&MaterialLibrary>) -> f64 {
    let Some(lib) = material_library else {
        return DEFAULT_OPAQUE_ABSORPTANCE;
    };
    let Some(mat) = lib.lookup(path) else {
        return DEFAULT_OPAQUE_ABSORPTANCE;
    };
    let Some(opt) = mat.optical.as_ref() else {
        return DEFAULT_OPAQUE_ABSORPTANCE;
    };

    // Convert RGB optical props into a scalar absorptance estimate.
    let mut a = 0.0;
    for c in 0..3 {
        let absorb =
            (1.0 - opt.diffuse_reflectance[c] - opt.specular_reflectance[c] - opt.transmittance[c])
                .clamp(0.0, 1.0);
        a += absorb;
    }
    (a / 3.0).clamp(0.0, 1.0)
}

/// Produces UID-keyed shortwave coupling payloads from an existing `LightingResult`.
///
/// This keeps "single source of truth" explicit: a composed pipeline can run
/// (a) a ray-based lighting simulation that writes `LightingResult`, then
/// (b) this module to convert it into shortwave payloads for thermal.
#[derive(Clone)]
pub struct LightingToShortwaveConfig {
    pub material_library: Option<MaterialLibrary>,
    pub default_reflectance: Rgb,
}

impl LightingToShortwaveConfig {
    pub fn new() -> Self {
        Self {
            material_library: None,
            default_reflectance: [0.5, 0.5, 0.5],
        }
    }
}

impl Default for LightingToShortwaveConfig {
    fn default() -> Self {
        Self::new()
    }
}

pub struct LightingToShortwaveModule {
    config: LightingToShortwaveConfig,
    has_run: bool,
}

impl LightingToShortwaveModule {
    pub fn new(config: LightingToShortwaveConfig) -> Self {
        Self {
            config,
            has_run: false,
        }
    }
}

impl SimModule for LightingToShortwaveModule {
    fn name(&self) -> &'static str {
        "lighting_to_shortwave"
    }

    fn step(&mut self, ctx: &SimContext, bus: &mut Bus) -> Result<()> {
        if self.has_run {
            return Ok(());
        }

        let Some(result) = bus.get::<LightingResult>() else {
            anyhow::bail!("LightingToShortwaveModule requires LightingResult on the Bus");
        };

        let (absorbed, transmitted) = shortwave_payloads_from_lighting(result, ctx, &self.config);
        bus.put(absorbed);
        bus.put(transmitted);

        self.has_run = true;
        Ok(())
    }
}

fn shortwave_payloads_from_lighting(
    result: &LightingResult,
    ctx: &SimContext,
    config: &LightingToShortwaveConfig,
) -> (ShortwaveAbsorbedWPerPolygon, ShortwaveTransmittedWPerZone) {
    let mut absorbed = ShortwaveAbsorbedWPerPolygon::default();
    let mut transmitted = ShortwaveTransmittedWPerZone::default();

    for (polygon_uid, flux_rgb) in &result.incident_flux {
        let Some(path) = ctx.surface_index.path_by_polygon_uid(polygon_uid) else {
            continue;
        };

        let optical = config
            .material_library
            .as_ref()
            .and_then(|lib| lib.lookup(path))
            .and_then(|mat| mat.optical.as_ref());

        let (diff, spec, trans) = optical_props_or_default(optical, config.default_reflectance);

        let absorbed_w = absorbed_w_from_flux(flux_rgb, diff, spec, trans);
        if absorbed_w > 0.0 {
            *absorbed
                .watts_by_polygon_uid
                .entry(polygon_uid.clone())
                .or_insert(0.0) += absorbed_w;
        }

        let transmitted_w = transmitted_w_from_flux(flux_rgb, trans);
        if transmitted_w > 0.0 {
            let Some(zone_uid) = ctx.surface_index.zone_uid_by_polygon_uid(polygon_uid) else {
                continue;
            };
            *transmitted
                .watts_by_zone_uid
                .entry(zone_uid.clone())
                .or_insert(0.0) += transmitted_w;
        }
    }

    (absorbed, transmitted)
}

fn optical_props_or_default(
    optical: Option<&OpticalMaterial>,
    default_diffuse: Rgb,
) -> (Rgb, Rgb, Rgb) {
    match optical {
        Some(opt) => (
            opt.diffuse_reflectance,
            opt.specular_reflectance,
            opt.transmittance,
        ),
        None => (default_diffuse, [0.0; 3], [0.0; 3]),
    }
}

fn absorbed_w_from_flux(flux: &Rgb, diff: Rgb, spec: Rgb, trans: Rgb) -> f64 {
    let mut absorbed = 0.0;
    for c in 0..3 {
        let absorbance = (1.0 - diff[c] - spec[c] - trans[c]).clamp(0.0, 1.0);
        absorbed += flux[c] * absorbance;
    }
    absorbed
}

fn transmitted_w_from_flux(flux: &Rgb, trans: Rgb) -> f64 {
    flux[0] * trans[0] + flux[1] * trans[1] + flux[2] * trans[2]
}

/// Deterministic (non-ray) solar shortwave producer based on EPW radiation inputs.
///
/// This module uses the same SHGC-based approximation as the energy module's
/// `compute_solar_gains_per_zone_with_materials()` and publishes the resulting gains as
/// `ShortwaveTransmittedWPerZone`.
///
/// It is intended as a "first producer" for composed pipelines before adding
/// ray-traced shading / complex glazing.
#[derive(Clone)]
pub struct SolarShortwaveConfig {
    pub solar_params: SolarHourParams,
    pub gain_config: SolarGainConfig,
    pub material_library: Option<MaterialLibrary>,
}

pub struct SolarShortwaveModule {
    config: SolarShortwaveConfig,
    has_run: bool,
}

impl SolarShortwaveModule {
    pub fn new(config: SolarShortwaveConfig) -> Self {
        Self {
            config,
            has_run: false,
        }
    }
}

impl SimModule for SolarShortwaveModule {
    fn name(&self) -> &'static str {
        "solar_shortwave"
    }

    fn step(&mut self, ctx: &SimContext, bus: &mut Bus) -> Result<()> {
        if self.has_run {
            return Ok(());
        }

        let gains_by_zone = compute_solar_gains_per_zone_with_materials(
            ctx.building,
            &self.config.solar_params,
            &self.config.gain_config,
            self.config.material_library.as_ref(),
        );

        let transmitted = ShortwaveTransmittedWPerZone {
            watts_by_zone_uid: gains_by_zone,
        };

        let boundaries = ThermalBoundaries::classify(ctx.building, ctx.surface_index);
        let absorbed = compute_unshaded_opaque_absorbed(
            ctx,
            &boundaries,
            &self.config.solar_params,
            &self.config.gain_config,
            self.config.material_library.as_ref(),
        );

        bus.put(transmitted);
        bus.put(absorbed);

        self.has_run = true;
        Ok(())
    }
}

/// EPW-driven (time-series) solar shortwave producer for composed pipelines.
///
/// This module consumes [`WeatherHourIndex`] (typically published by `energy::WeatherModule`)
/// and uses the referenced EPW record (DNI/DHI + timestamp) to compute per-zone solar gains.
#[derive(Clone)]
pub struct SolarEpwConfig {
    pub weather: Arc<WeatherData>,
    pub gain_config: SolarGainConfig,
    pub material_library: Option<MaterialLibrary>,
}

pub struct SolarEpwModule {
    config: SolarEpwConfig,
    boundaries: Option<ThermalBoundaries>,
}

impl SolarEpwModule {
    pub fn new(config: SolarEpwConfig) -> Self {
        Self {
            config,
            boundaries: None,
        }
    }
}

impl SimModule for SolarEpwModule {
    fn name(&self) -> &'static str {
        "solar_epw_shortwave"
    }

    fn init(&mut self, ctx: &SimContext, _bus: &mut Bus) -> Result<()> {
        if self.boundaries.is_none() {
            self.boundaries = Some(ThermalBoundaries::classify(ctx.building, ctx.surface_index));
        }
        Ok(())
    }

    fn step(&mut self, ctx: &SimContext, bus: &mut Bus) -> Result<()> {
        if self.boundaries.is_none() {
            self.boundaries = Some(ThermalBoundaries::classify(ctx.building, ctx.surface_index));
        }

        let Some(hour_index) = bus.get::<WeatherHourIndex>().map(|i| i.0) else {
            anyhow::bail!("SolarEpwModule requires WeatherHourIndex on the Bus");
        };
        anyhow::ensure!(
            hour_index < self.config.weather.records.len(),
            "SolarEpwModule: hour index {hour_index} out of range (len={})",
            self.config.weather.records.len()
        );
        let record = &self.config.weather.records[hour_index];

        let params = SolarHourParams {
            global_horizontal_irradiance: record.global_horizontal_radiation,
            direct_normal_irradiance: record.direct_normal_radiation,
            diffuse_horizontal_irradiance: record.diffuse_horizontal_radiation,
            day_of_year: day_of_year(record.month, record.day),
            local_time_hours: record.hour as f64 - 0.5,
            latitude: self.config.weather.latitude,
            longitude: self.config.weather.longitude,
            timezone: self.config.weather.timezone,
        };

        let gains_by_zone = compute_solar_gains_per_zone_with_materials(
            ctx.building,
            &params,
            &self.config.gain_config,
            self.config.material_library.as_ref(),
        );

        let transmitted = ShortwaveTransmittedWPerZone {
            watts_by_zone_uid: gains_by_zone,
        };

        let boundaries = self.boundaries.as_ref().expect("set above");
        let absorbed = compute_unshaded_opaque_absorbed(
            ctx,
            boundaries,
            &params,
            &self.config.gain_config,
            self.config.material_library.as_ref(),
        );

        bus.put(transmitted);
        bus.put(absorbed);
        Ok(())
    }
}

/// EPW-driven (time-series) solar shortwave producer that reads `Arc<WeatherData>` from the `Bus`.
///
/// This is a convenience wrapper for composed pipelines that already include `energy::WeatherModule`.
/// It avoids duplicating EPW dataset plumbing across module configs.
#[derive(Clone)]
pub struct SolarEpwBusConfig {
    pub gain_config: SolarGainConfig,
    pub material_library: Option<MaterialLibrary>,
}

impl SolarEpwBusConfig {
    pub fn new(gain_config: SolarGainConfig) -> Self {
        Self {
            gain_config,
            material_library: None,
        }
    }
}

pub struct SolarEpwBusModule {
    config: SolarEpwBusConfig,
    inner: Option<SolarEpwModule>,
}

impl SolarEpwBusModule {
    pub fn new(config: SolarEpwBusConfig) -> Self {
        Self {
            config,
            inner: None,
        }
    }

    fn ensure_initialized(&mut self, ctx: &SimContext, bus: &mut Bus) -> Result<()> {
        if self.inner.is_some() {
            return Ok(());
        }

        let Some(weather) = bus.get::<Arc<WeatherData>>().cloned() else {
            anyhow::bail!(
                "SolarEpwBusModule requires Arc<WeatherData> on the Bus (published by WeatherModule::init)"
            );
        };

        let mut inner = SolarEpwModule::new(SolarEpwConfig {
            weather,
            gain_config: self.config.gain_config.clone(),
            material_library: self.config.material_library.clone(),
        });
        inner.init(ctx, bus)?;
        self.inner = Some(inner);
        Ok(())
    }
}

impl SimModule for SolarEpwBusModule {
    fn name(&self) -> &'static str {
        "solar_epw_shortwave_bus"
    }

    fn init(&mut self, ctx: &SimContext, bus: &mut Bus) -> Result<()> {
        if bus.get::<Arc<WeatherData>>().is_some() {
            self.ensure_initialized(ctx, bus)?;
        }
        Ok(())
    }

    fn step(&mut self, ctx: &SimContext, bus: &mut Bus) -> Result<()> {
        self.ensure_initialized(ctx, bus)?;
        self.inner
            .as_mut()
            .expect("initialized above")
            .step(ctx, bus)
    }
}

/// EPW-driven (time-series) solar shortwave producer with direct-sun occlusion.
///
/// This is the `WeatherHourIndex`-based counterpart of [`SolarShortwaveShadedStepModule`].
/// It is intended for composed pipelines where a separate weather module publishes the
/// current hour index and outdoor temperature.
///
/// Each `step()` publishes:
/// - [`ShortwaveTransmittedWPerZone`] (per zone, shaded direct + unshaded diffuse)
/// - [`ShortwaveAbsorbedWPerPolygon`] (currently empty placeholder)
#[derive(Clone)]
pub struct SolarEpwShadedConfig {
    pub weather: Arc<WeatherData>,
    pub gain_config: SolarGainConfig,
    pub material_library: Option<MaterialLibrary>,
    /// Voxel size for `FlatScene` acceleration.
    pub voxel_size: f64,
}

impl SolarEpwShadedConfig {
    pub fn new(weather: Arc<WeatherData>, gain_config: SolarGainConfig) -> Self {
        Self {
            weather,
            gain_config,
            material_library: None,
            voxel_size: 0.5,
        }
    }
}

pub struct SolarEpwShadedModule {
    config: SolarEpwShadedConfig,
    scene: Option<FlatScene>,
    polygon_idx_by_uid: HashMap<crate::UID, usize>,
    boundaries: Option<ThermalBoundaries>,
}

impl SolarEpwShadedModule {
    pub fn new(config: SolarEpwShadedConfig) -> Self {
        Self {
            config,
            scene: None,
            polygon_idx_by_uid: HashMap::new(),
            boundaries: None,
        }
    }

    fn ensure_initialized(&mut self, ctx: &SimContext) {
        if self.scene.is_some() {
            return;
        }

        let scene = FlatScene::new(ctx.building, self.config.voxel_size, true);
        self.polygon_idx_by_uid = scene
            .polygons
            .iter()
            .enumerate()
            .map(|(idx, p)| (p.uid.clone(), idx))
            .collect();
        self.scene = Some(scene);
        self.boundaries = Some(ThermalBoundaries::classify(ctx.building, ctx.surface_index));
    }
}

impl SimModule for SolarEpwShadedModule {
    fn name(&self) -> &'static str {
        "solar_epw_shaded_shortwave"
    }

    fn init(&mut self, ctx: &SimContext, _bus: &mut Bus) -> Result<()> {
        self.ensure_initialized(ctx);
        Ok(())
    }

    fn step(&mut self, ctx: &SimContext, bus: &mut Bus) -> Result<()> {
        self.ensure_initialized(ctx);

        let Some(hour_idx) = bus.get::<WeatherHourIndex>().map(|i| i.0) else {
            anyhow::bail!("SolarEpwShadedModule requires WeatherHourIndex on the Bus");
        };
        anyhow::ensure!(
            hour_idx < self.config.weather.records.len(),
            "SolarEpwShadedModule: hour_idx {hour_idx} out of range (len={})",
            self.config.weather.records.len()
        );

        let scene = self.scene.as_ref().expect("initialized");
        let boundaries = self.boundaries.as_ref().expect("initialized");

        let record = &self.config.weather.records[hour_idx];
        let solar_pos = SolarPosition::calculate_from_local_time(
            self.config.weather.latitude,
            self.config.weather.longitude,
            self.config.weather.timezone,
            day_of_year(record.month, record.day),
            record.hour as f64 - 0.5,
        );

        let (transmitted, absorbed) = compute_shaded_shortwave(
            ctx,
            scene,
            &self.polygon_idx_by_uid,
            boundaries,
            solar_pos.to_direction(),
            solar_pos.is_above_horizon(),
            record.direct_normal_radiation.max(0.0),
            record.diffuse_horizontal_radiation.max(0.0),
            &self.config.gain_config,
            self.config.material_library.as_ref(),
        );

        bus.put(transmitted);
        bus.put(absorbed);
        Ok(())
    }
}

/// Shaded EPW-driven (time-series) shortwave producer that reads `Arc<WeatherData>` from the `Bus`.
///
/// Convenience wrapper for composed pipelines that already include `energy::WeatherModule`.
#[derive(Clone)]
pub struct SolarEpwShadedBusConfig {
    pub gain_config: SolarGainConfig,
    pub material_library: Option<MaterialLibrary>,
    pub voxel_size: f64,
}

impl SolarEpwShadedBusConfig {
    pub fn new(gain_config: SolarGainConfig) -> Self {
        Self {
            gain_config,
            material_library: None,
            voxel_size: 0.5,
        }
    }
}

pub struct SolarEpwShadedBusModule {
    config: SolarEpwShadedBusConfig,
    inner: Option<SolarEpwShadedModule>,
}

impl SolarEpwShadedBusModule {
    pub fn new(config: SolarEpwShadedBusConfig) -> Self {
        Self {
            config,
            inner: None,
        }
    }

    fn ensure_initialized(&mut self, ctx: &SimContext, bus: &mut Bus) -> Result<()> {
        if self.inner.is_some() {
            return Ok(());
        }

        let Some(weather) = bus.get::<Arc<WeatherData>>().cloned() else {
            anyhow::bail!(
                "SolarEpwShadedBusModule requires Arc<WeatherData> on the Bus (published by WeatherModule::init)"
            );
        };

        let mut inner = SolarEpwShadedModule::new(SolarEpwShadedConfig {
            weather,
            gain_config: self.config.gain_config.clone(),
            material_library: self.config.material_library.clone(),
            voxel_size: self.config.voxel_size,
        });
        inner.init(ctx, bus)?;
        self.inner = Some(inner);
        Ok(())
    }
}

impl SimModule for SolarEpwShadedBusModule {
    fn name(&self) -> &'static str {
        "solar_epw_shaded_shortwave_bus"
    }

    fn init(&mut self, ctx: &SimContext, bus: &mut Bus) -> Result<()> {
        if bus.get::<Arc<WeatherData>>().is_some() {
            self.ensure_initialized(ctx, bus)?;
        }
        Ok(())
    }

    fn step(&mut self, ctx: &SimContext, bus: &mut Bus) -> Result<()> {
        self.ensure_initialized(ctx, bus)?;
        self.inner
            .as_mut()
            .expect("initialized above")
            .step(ctx, bus)
    }
}

/// Step-based EPW-driven solar shortwave producer for composed pipelines.
///
/// Each `step()` publishes:
/// - [`OutdoorAirTemperatureC`] (from the weather record)
/// - [`ShortwaveTransmittedWPerZone`] (SHGC-based approximation, per zone)
/// - [`ShortwaveAbsorbedWPerPolygon`] (currently empty placeholder)
///
/// This is intended to pair with `EnergyModule` (`sim::energy::module`) in a pipeline where the
/// caller advances time by repeatedly calling `Pipeline::step()`.
#[derive(Clone)]
pub struct SolarShortwaveStepConfig {
    pub weather: WeatherData,
    pub gain_config: SolarGainConfig,
    pub material_library: Option<MaterialLibrary>,
    /// Starting weather record index (0-based).
    pub start_hour_idx: usize,
}

impl SolarShortwaveStepConfig {
    pub fn new(weather: WeatherData, gain_config: SolarGainConfig) -> Self {
        Self {
            weather,
            gain_config,
            material_library: None,
            start_hour_idx: 0,
        }
    }
}

pub struct SolarShortwaveStepModule {
    config: SolarShortwaveStepConfig,
    hour_idx: usize,
    boundaries: Option<ThermalBoundaries>,
}

impl SolarShortwaveStepModule {
    pub fn new(config: SolarShortwaveStepConfig) -> Self {
        let hour_idx = config.start_hour_idx;
        Self {
            config,
            hour_idx,
            boundaries: None,
        }
    }

    pub fn hour_idx(&self) -> usize {
        self.hour_idx
    }
}

impl SimModule for SolarShortwaveStepModule {
    fn name(&self) -> &'static str {
        "solar_shortwave_step"
    }

    fn init(&mut self, ctx: &SimContext, _bus: &mut Bus) -> Result<()> {
        if self.boundaries.is_none() {
            self.boundaries = Some(ThermalBoundaries::classify(ctx.building, ctx.surface_index));
        }
        Ok(())
    }

    fn step(&mut self, ctx: &SimContext, bus: &mut Bus) -> Result<()> {
        if self.boundaries.is_none() {
            self.boundaries = Some(ThermalBoundaries::classify(ctx.building, ctx.surface_index));
        }

        if self.hour_idx >= self.config.weather.records.len() {
            anyhow::bail!(
                "SolarShortwaveStepModule: hour_idx {} out of range (len={})",
                self.hour_idx,
                self.config.weather.records.len()
            );
        }

        let record = &self.config.weather.records[self.hour_idx];
        bus.put(OutdoorAirTemperatureC(record.dry_bulb_temperature));

        let params = SolarHourParams {
            global_horizontal_irradiance: record.global_horizontal_radiation,
            direct_normal_irradiance: record.direct_normal_radiation,
            diffuse_horizontal_irradiance: record.diffuse_horizontal_radiation,
            day_of_year: day_of_year(record.month, record.day),
            local_time_hours: record.hour as f64 - 0.5,
            latitude: self.config.weather.latitude,
            longitude: self.config.weather.longitude,
            timezone: self.config.weather.timezone,
        };

        let gains_by_zone = compute_solar_gains_per_zone_with_materials(
            ctx.building,
            &params,
            &self.config.gain_config,
            self.config.material_library.as_ref(),
        );

        let transmitted = ShortwaveTransmittedWPerZone {
            watts_by_zone_uid: gains_by_zone,
        };

        let boundaries = self.boundaries.as_ref().expect("set above");
        let absorbed = compute_unshaded_opaque_absorbed(
            ctx,
            boundaries,
            &params,
            &self.config.gain_config,
            self.config.material_library.as_ref(),
        );

        bus.put(transmitted);
        bus.put(absorbed);

        self.hour_idx += 1;
        Ok(())
    }
}

/// Step-based EPW-driven solar shortwave producer with direct-sun occlusion.
///
/// This is an incremental bridge between the deterministic EPW+SHGC model and a
/// fully ray-traced lighting-driven producer:
/// - Direct component (DNI) is reduced by a hard-shadow visibility fraction computed
///   via `FlatScene` ray casting.
/// - Diffuse component (DHI) remains the isotropic-sky approximation (no occlusion yet).
///
/// Each `step()` publishes:
/// - [`OutdoorAirTemperatureC`] (from the weather record)
/// - [`ShortwaveTransmittedWPerZone`] (per zone, shaded direct + unshaded diffuse)
/// - [`ShortwaveAbsorbedWPerPolygon`] (currently empty placeholder)
#[derive(Clone)]
pub struct SolarShortwaveShadedStepConfig {
    pub weather: WeatherData,
    pub gain_config: SolarGainConfig,
    pub material_library: Option<MaterialLibrary>,
    /// Starting weather record index (0-based).
    pub start_hour_idx: usize,
    /// Voxel size for `FlatScene` acceleration.
    pub voxel_size: f64,
}

impl SolarShortwaveShadedStepConfig {
    pub fn new(weather: WeatherData, gain_config: SolarGainConfig) -> Self {
        Self {
            weather,
            gain_config,
            material_library: None,
            start_hour_idx: 0,
            voxel_size: 0.5,
        }
    }
}

pub struct SolarShortwaveShadedStepModule {
    config: SolarShortwaveShadedStepConfig,
    hour_idx: usize,
    scene: Option<FlatScene>,
    polygon_idx_by_uid: HashMap<crate::UID, usize>,
    boundaries: Option<ThermalBoundaries>,
}

impl SolarShortwaveShadedStepModule {
    pub fn new(config: SolarShortwaveShadedStepConfig) -> Self {
        let hour_idx = config.start_hour_idx;
        Self {
            config,
            hour_idx,
            scene: None,
            polygon_idx_by_uid: HashMap::new(),
            boundaries: None,
        }
    }

    pub fn hour_idx(&self) -> usize {
        self.hour_idx
    }

    fn ensure_initialized(&mut self, ctx: &SimContext) {
        if self.scene.is_some() {
            return;
        }

        let scene = FlatScene::new(ctx.building, self.config.voxel_size, true);
        self.polygon_idx_by_uid = scene
            .polygons
            .iter()
            .enumerate()
            .map(|(idx, p)| (p.uid.clone(), idx))
            .collect();
        self.scene = Some(scene);
        self.boundaries = Some(ThermalBoundaries::classify(ctx.building, ctx.surface_index));
    }
}

impl SimModule for SolarShortwaveShadedStepModule {
    fn name(&self) -> &'static str {
        "solar_shortwave_shaded_step"
    }

    fn init(&mut self, ctx: &SimContext, _bus: &mut Bus) -> Result<()> {
        self.ensure_initialized(ctx);
        Ok(())
    }

    fn step(&mut self, ctx: &SimContext, bus: &mut Bus) -> Result<()> {
        self.ensure_initialized(ctx);

        if self.hour_idx >= self.config.weather.records.len() {
            anyhow::bail!(
                "SolarShortwaveShadedStepModule: hour_idx {} out of range (len={})",
                self.hour_idx,
                self.config.weather.records.len()
            );
        }

        let scene = self.scene.as_ref().expect("initialized");
        let boundaries = self.boundaries.as_ref().expect("initialized");

        let record = &self.config.weather.records[self.hour_idx];
        bus.put(OutdoorAirTemperatureC(record.dry_bulb_temperature));

        let solar_pos = SolarPosition::calculate_from_local_time(
            self.config.weather.latitude,
            self.config.weather.longitude,
            self.config.weather.timezone,
            day_of_year(record.month, record.day),
            record.hour as f64 - 0.5,
        );

        let (transmitted, absorbed) = compute_shaded_shortwave(
            ctx,
            scene,
            &self.polygon_idx_by_uid,
            boundaries,
            solar_pos.to_direction(),
            solar_pos.is_above_horizon(),
            record.direct_normal_radiation.max(0.0),
            record.diffuse_horizontal_radiation.max(0.0),
            &self.config.gain_config,
            self.config.material_library.as_ref(),
        );

        bus.put(transmitted);
        bus.put(absorbed);

        self.hour_idx += 1;
        Ok(())
    }
}

/// Computes absorbed shortwave on exterior opaque surfaces (unshaded).
///
/// Shared helper for `SolarShortwaveModule`, `SolarEpwModule`, and
/// `SolarShortwaveStepModule`.  Glazing surfaces (identified via SHGC) are
/// skipped — their energy is accounted for via the transmitted-to-zone path.
fn compute_unshaded_opaque_absorbed(
    ctx: &SimContext,
    boundaries: &ThermalBoundaries,
    params: &SolarHourParams,
    gain_config: &SolarGainConfig,
    material_library: Option<&MaterialLibrary>,
) -> ShortwaveAbsorbedWPerPolygon {
    let solar_pos = SolarPosition::calculate_from_local_time(
        params.latitude,
        params.longitude,
        params.timezone,
        params.day_of_year,
        params.local_time_hours,
    );
    let sun_dir = solar_pos.to_direction();
    let sun_above = solar_pos.is_above_horizon();

    let mut absorbed = ShortwaveAbsorbedWPerPolygon::default();

    for surface in &ctx.surface_index.surfaces {
        if !boundaries.is_exterior(&surface.polygon_uid) {
            continue;
        }
        if gain_config
            .resolve_shgc(&surface.path, material_library)
            .is_some()
        {
            continue; // glazing handled via transmitted-to-zone
        }

        let Some(poly) = ctx.building.get_polygon(&surface.path) else {
            continue;
        };
        let area = surface.area_m2;
        if area <= 0.0 {
            continue;
        }

        let normal = poly.vn;
        let mut incident = 0.0;
        let sky_view = 0.5 * (1.0 + normal.dz.max(0.0));
        incident += params.diffuse_horizontal_irradiance.max(0.0) * sky_view;

        if sun_above && params.direct_normal_irradiance > 0.0 {
            let cos_incidence = sun_dir.dot(&normal).max(0.0);
            incident += params.direct_normal_irradiance.max(0.0) * cos_incidence;
        }

        if incident <= 0.0 {
            continue;
        }

        let absorptance = opaque_absorptance_for_path(&surface.path, material_library);
        let q_abs = incident * area * absorptance;
        if q_abs > 0.0 {
            absorbed
                .watts_by_polygon_uid
                .insert(surface.polygon_uid.clone(), q_abs);
        }
    }

    absorbed
}

/// Computes both transmitted and absorbed shortwave with direct-sun occlusion.
///
/// Shared helper for `SolarEpwShadedModule` and `SolarShortwaveShadedStepModule`.
/// Direct component is reduced by a visibility fraction computed via `FlatScene`
/// ray casting; diffuse component uses the isotropic sky approximation.
#[allow(clippy::too_many_arguments)]
fn compute_shaded_shortwave(
    ctx: &SimContext,
    scene: &FlatScene,
    polygon_idx_by_uid: &HashMap<crate::UID, usize>,
    boundaries: &ThermalBoundaries,
    sun_dir: Vector,
    sun_above: bool,
    dni: f64,
    dhi: f64,
    gain_config: &SolarGainConfig,
    material_library: Option<&MaterialLibrary>,
) -> (ShortwaveTransmittedWPerZone, ShortwaveAbsorbedWPerPolygon) {
    let mut gains_by_zone: HashMap<crate::UID, f64> = HashMap::new();
    let mut absorbed = ShortwaveAbsorbedWPerPolygon::default();

    for surface in &ctx.surface_index.surfaces {
        if !boundaries.is_exterior(&surface.polygon_uid) {
            continue;
        }

        let shgc = gain_config
            .resolve_shgc(&surface.path, material_library)
            .map(|v| v.clamp(0.0, 1.0));

        let Some(&poly_idx) = polygon_idx_by_uid.get(&surface.polygon_uid) else {
            continue;
        };
        let poly = &scene.polygons[poly_idx];
        let area = poly.area();
        if area <= 0.0 {
            continue;
        }

        let normal = poly.vn;
        let mut q = 0.0;

        // Diffuse component: isotropic sky view factor.
        let sky_view = 0.5 * (1.0 + normal.dz.max(0.0));
        q += dhi * sky_view;

        // Direct component: DNI * cos(incidence) * visibility.
        if sun_above && dni > 0.0 {
            let cos_incidence = sun_dir.dot(&normal).max(0.0);
            if cos_incidence > 0.0 {
                let vis = visibility_fraction(scene, poly_idx, sun_dir);
                q += dni * cos_incidence * vis;
            }
        }

        if q <= 0.0 {
            continue;
        }

        if let Some(shgc) = shgc {
            let q_trans = q * area * shgc;
            if q_trans != 0.0 {
                *gains_by_zone.entry(surface.zone_uid.clone()).or_insert(0.0) += q_trans;
            }
        } else {
            let absorptance = opaque_absorptance_for_path(&surface.path, material_library);
            let q_abs = q * area * absorptance;
            if q_abs > 0.0 {
                absorbed
                    .watts_by_polygon_uid
                    .insert(surface.polygon_uid.clone(), q_abs);
            }
        }
    }

    let transmitted = ShortwaveTransmittedWPerZone {
        watts_by_zone_uid: gains_by_zone,
    };
    (transmitted, absorbed)
}

fn day_of_year(month: u8, day: u8) -> u16 {
    const DAYS_BEFORE_MONTH: [u16; 12] = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];
    let m = (month as usize).saturating_sub(1).min(11);
    DAYS_BEFORE_MONTH[m] + day as u16
}

fn visibility_fraction(scene: &FlatScene, polygon_idx: usize, sun_dir: Vector) -> f64 {
    let Ok(dir) = sun_dir.normalize() else {
        return 0.0;
    };

    let poly = &scene.polygons[polygon_idx];
    let sample_points = triangle_centroids(poly);
    if sample_points.is_empty() {
        return 0.0;
    }

    // Offset along the ray direction to avoid immediate self-intersection at t≈0.
    const ORIGIN_EPS: f64 = 1e-4;

    let mut visible = 0usize;
    for p in &sample_points {
        let origin = *p + dir * ORIGIN_EPS;
        match scene.find_target_surface(origin, dir) {
            None => visible += 1,
            Some((hit_idx, _t)) if hit_idx == polygon_idx => visible += 1,
            _ => {}
        }
    }

    visible as f64 / sample_points.len() as f64
}

fn triangle_centroids(poly: &crate::Polygon) -> Vec<Point> {
    let verts = poly.vertices();
    let Some(tris) = poly.triangles() else {
        return vec![];
    };

    let mut pts = Vec::with_capacity(tris.len());
    for tri in tris {
        let i0 = tri.0;
        let i1 = tri.1;
        let i2 = tri.2;
        if i0 >= verts.len() || i1 >= verts.len() || i2 >= verts.len() {
            continue;
        }
        let p0 = verts[i0];
        let p1 = verts[i1];
        let p2 = verts[i2];
        pts.push(Point::new(
            (p0.x + p1.x + p2.x) / 3.0,
            (p0.y + p1.y + p2.y) / 3.0,
            (p0.z + p1.z + p2.z) / 3.0,
        ));
    }
    pts
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::energy::weather::HourlyRecord;
    use crate::sim::index::SurfaceIndex;
    use crate::{Building, Polygon, Solid, Vector, Wall, Zone};

    #[test]
    fn test_shortwave_payloads_from_lighting_uid_keyed() -> Result<()> {
        let poly = Polygon::new(
            "glass",
            vec![
                Point::new(0.0, 0.0, 0.0),
                Point::new(1.0, 0.0, 0.0),
                Point::new(1.0, 1.0, 0.0),
                Point::new(0.0, 1.0, 0.0),
            ],
            None,
        )?;
        let wall = Wall::new("window", vec![poly])?;
        let solid = Solid::new("room", vec![wall])?;
        let zone = Zone::new("z", vec![solid])?;
        let building = Building::new("b", vec![zone])?;

        let index = SurfaceIndex::new(&building);
        let ctx = SimContext::new(&building, &index);

        let polygon_uid = index.surfaces.first().unwrap().polygon_uid.clone();

        let mut lighting = LightingResult::new();
        lighting
            .incident_flux
            .insert(polygon_uid.clone(), [10.0, 10.0, 10.0]);

        let mut bus = Bus::new();
        bus.put(lighting);

        let mut config = LightingToShortwaveConfig::new();
        config.default_reflectance = [0.5, 0.5, 0.5]; // absorbance 0.5

        let mut module = LightingToShortwaveModule::new(config);
        module.step(&ctx, &mut bus)?;

        let absorbed = bus.get::<ShortwaveAbsorbedWPerPolygon>().unwrap();
        let a = absorbed.watts_by_polygon_uid.get(&polygon_uid).unwrap();
        assert!((a - 15.0).abs() < 1e-10); // sum(10 * 0.5)

        Ok(())
    }

    #[test]
    fn test_solar_shortwave_module_produces_zone_gains() -> Result<()> {
        // Single zone with a 1m² "window" polygon on a horizontal plane.
        let poly = Polygon::new(
            "glass",
            vec![
                Point::new(0.0, 0.0, 0.0),
                Point::new(1.0, 0.0, 0.0),
                Point::new(1.0, 1.0, 0.0),
                Point::new(0.0, 1.0, 0.0),
            ],
            None,
        )?;
        let wall = Wall::new("window", vec![poly])?;
        let solid = Solid::new("room", vec![wall])?;
        let zone = Zone::new("zone", vec![solid])?;
        let building = Building::new("b", vec![zone])?;

        let index = SurfaceIndex::new(&building);
        let ctx = SimContext::new(&building, &index);

        let params = SolarHourParams {
            global_horizontal_irradiance: 900.0,
            direct_normal_irradiance: 500.0,
            diffuse_horizontal_irradiance: 200.0,
            day_of_year: 80, // ~equinox
            local_time_hours: 12.0,
            latitude: 0.0,
            longitude: 0.0,
            timezone: 0.0,
        };

        let config = SolarShortwaveConfig {
            solar_params: params,
            gain_config: SolarGainConfig::new(), // default patterns include "window"/"glass"
            material_library: None,
        };

        let mut bus = Bus::new();
        let mut module = SolarShortwaveModule::new(config);
        module.step(&ctx, &mut bus)?;

        let transmitted = bus.get::<ShortwaveTransmittedWPerZone>().unwrap();
        let zone_uid = ctx.building.zones().first().unwrap().uid.clone();
        let q = transmitted
            .watts_by_zone_uid
            .get(&zone_uid)
            .cloned()
            .unwrap_or(0.0);
        // For equinox at equator at noon, cos_incidence should be near 1 for a +Z surface.
        // Expected approx: (500 + 200) * area(1) * SHGC(0.6) = 420 W.
        assert!((q - 420.0).abs() < 60.0, "q={q}");
        Ok(())
    }

    #[test]
    fn test_same_zone_split_invariant_unshaded_opaque_absorbed_total() {
        let one = {
            let s = Solid::from_box(2.0, 1.0, 1.0, None, "s").unwrap();
            let z = Zone::new("z", vec![s]).unwrap();
            Building::new("b", vec![z]).unwrap()
        };
        let split = {
            let s0 = Solid::from_box(1.0, 1.0, 1.0, None, "s0").unwrap();
            let s1 = Solid::from_box(1.0, 1.0, 1.0, Some((1.0, 0.0, 0.0)), "s1").unwrap();
            let z = Zone::new("z", vec![s0, s1]).unwrap();
            Building::new("b", vec![z]).unwrap()
        };

        let params = SolarHourParams {
            global_horizontal_irradiance: 1000.0,
            direct_normal_irradiance: 800.0,
            diffuse_horizontal_irradiance: 200.0,
            day_of_year: 80,
            local_time_hours: 12.0,
            latitude: 0.0,
            longitude: 0.0,
            timezone: 0.0,
        };
        let gain_cfg = SolarGainConfig::new();

        let index_one = SurfaceIndex::new(&one);
        let ctx_one = SimContext::new(&one, &index_one);
        let b_one = ThermalBoundaries::classify(&one, &index_one);
        let absorbed_one =
            compute_unshaded_opaque_absorbed(&ctx_one, &b_one, &params, &gain_cfg, None);
        let total_one: f64 = absorbed_one.watts_by_polygon_uid.values().sum();

        let index_split = SurfaceIndex::new(&split);
        let ctx_split = SimContext::new(&split, &index_split);
        let b_split = ThermalBoundaries::classify(&split, &index_split);
        let absorbed_split =
            compute_unshaded_opaque_absorbed(&ctx_split, &b_split, &params, &gain_cfg, None);
        let total_split: f64 = absorbed_split.watts_by_polygon_uid.values().sum();

        assert!(
            (total_one - total_split).abs() < 1e-8,
            "Expected same total absorbed shortwave. one={total_one:.6}, split={total_split:.6}"
        );
    }

    #[test]
    fn test_solar_epw_module_consumes_weather_hour_index() -> Result<()> {
        // Single zone with a 1m² "glass" polygon facing up: diffuse-only should add gains.
        let poly = Polygon::new(
            "glass",
            vec![
                Point::new(0.0, 0.0, 0.0),
                Point::new(1.0, 0.0, 0.0),
                Point::new(1.0, 1.0, 0.0),
                Point::new(0.0, 1.0, 0.0),
            ],
            None,
        )?;
        let wall = Wall::new("window", vec![poly])?;
        let solid = Solid::new("room", vec![wall])?;
        let zone = Zone::new("z", vec![solid])?;
        let building = Building::new("b", vec![zone])?;

        let index = SurfaceIndex::new(&building);
        let ctx = SimContext::new(&building, &index);

        let mut weather = WeatherData::synthetic("X", 52.0, 13.0, 10.0, 0.0);
        weather.records[0].diffuse_horizontal_radiation = 100.0;
        weather.records[0].direct_normal_radiation = 0.0;
        let weather = Arc::new(weather);

        let mut bus = Bus::new();
        bus.put(WeatherHourIndex(0));

        let mut module = SolarEpwModule::new(SolarEpwConfig {
            weather,
            gain_config: SolarGainConfig::new(),
            material_library: None,
        });
        module.step(&ctx, &mut bus)?;

        let transmitted = bus.get::<ShortwaveTransmittedWPerZone>().unwrap();
        let zone_uid = ctx.building.zones().first().unwrap().uid.clone();
        let g = transmitted
            .watts_by_zone_uid
            .get(&zone_uid)
            .cloned()
            .unwrap_or(0.0);
        assert!(g > 0.0);
        Ok(())
    }

    #[test]
    fn test_solar_epw_bus_module_consumes_weather_from_bus() -> Result<()> {
        let poly = Polygon::new(
            "glass",
            vec![
                Point::new(0.0, 0.0, 0.0),
                Point::new(1.0, 0.0, 0.0),
                Point::new(1.0, 1.0, 0.0),
                Point::new(0.0, 1.0, 0.0),
            ],
            None,
        )?;
        let wall = Wall::new("window", vec![poly])?;
        let solid = Solid::new("room", vec![wall])?;
        let zone = Zone::new("z", vec![solid])?;
        let building = Building::new("b", vec![zone])?;

        let index = SurfaceIndex::new(&building);
        let ctx = SimContext::new(&building, &index);

        let mut weather = WeatherData::synthetic("X", 52.0, 13.0, 10.0, 0.0);
        weather.records[0].diffuse_horizontal_radiation = 100.0;
        weather.records[0].direct_normal_radiation = 0.0;
        let weather = Arc::new(weather);

        let mut bus = Bus::new();
        bus.put(weather);
        bus.put(WeatherHourIndex(0));

        let mut module = SolarEpwBusModule::new(SolarEpwBusConfig::new(SolarGainConfig::new()));
        module.init(&ctx, &mut bus)?;
        module.step(&ctx, &mut bus)?;

        let transmitted = bus.get::<ShortwaveTransmittedWPerZone>().unwrap();
        let zone_uid = ctx.building.zones().first().unwrap().uid.clone();
        let g = transmitted
            .watts_by_zone_uid
            .get(&zone_uid)
            .cloned()
            .unwrap_or(0.0);
        assert!(g > 0.0);
        Ok(())
    }

    #[test]
    fn test_solar_epw_shaded_bus_module_consumes_weather_from_bus() -> Result<()> {
        let win = Polygon::new(
            "glass",
            vec![
                Point::new(0.0, 0.0, 0.0),
                Point::new(0.0, 1.0, 0.0),
                Point::new(0.0, 1.0, 1.0),
                Point::new(0.0, 0.0, 1.0),
            ],
            Some(Vector::new(1.0, 0.0, 0.0)),
        )?;
        let wall = Wall::new("window", vec![win])?;
        let solid = Solid::new("room", vec![wall])?;
        let zone = Zone::new("z", vec![solid])?;
        let building = Building::new("b", vec![zone])?;

        let index = SurfaceIndex::new(&building);
        let ctx = SimContext::new(&building, &index);

        let mut weather = WeatherData::synthetic("X", 0.0, 0.0, 10.0, 0.0);
        weather.records[0].month = 3;
        weather.records[0].day = 21;
        weather.records[0].hour = 7;
        weather.records[0].direct_normal_radiation = 1000.0;
        weather.records[0].diffuse_horizontal_radiation = 0.0;
        let weather = Arc::new(weather);

        let mut bus = Bus::new();
        bus.put(weather);
        bus.put(WeatherHourIndex(0));

        let mut module = SolarEpwShadedBusModule::new(SolarEpwShadedBusConfig {
            gain_config: SolarGainConfig::new(),
            material_library: None,
            voxel_size: 0.25,
        });
        module.init(&ctx, &mut bus)?;
        module.step(&ctx, &mut bus)?;

        let transmitted = bus.get::<ShortwaveTransmittedWPerZone>().unwrap();
        let zone_uid = ctx.building.zones().first().unwrap().uid.clone();
        let g = transmitted
            .watts_by_zone_uid
            .get(&zone_uid)
            .cloned()
            .unwrap_or(0.0);
        assert!(g > 0.0);
        Ok(())
    }

    #[test]
    fn test_visibility_fraction_unoccluded() -> Result<()> {
        let win = Polygon::new(
            "win",
            vec![
                Point::new(0.0, 0.0, 0.0),
                Point::new(0.0, 1.0, 0.0),
                Point::new(0.0, 1.0, 1.0),
                Point::new(0.0, 0.0, 1.0),
            ],
            Some(Vector::new(1.0, 0.0, 0.0)),
        )?;
        let wall = Wall::new("window", vec![win])?;
        let solid = Solid::new("room", vec![wall])?;
        let zone = Zone::new("z", vec![solid])?;
        let building = Building::new("b", vec![zone])?;

        let scene = FlatScene::new(&building, 0.25, true);
        let idx = scene
            .paths
            .iter()
            .position(|p| p.ends_with("/win"))
            .unwrap();
        let vis = visibility_fraction(&scene, idx, Vector::new(1.0, 0.0, 0.0));
        assert!((vis - 1.0).abs() < 1e-12, "vis={vis}");
        Ok(())
    }

    #[test]
    fn test_visibility_fraction_occluded() -> Result<()> {
        let win = Polygon::new(
            "win",
            vec![
                Point::new(0.0, 0.0, 0.0),
                Point::new(0.0, 1.0, 0.0),
                Point::new(0.0, 1.0, 1.0),
                Point::new(0.0, 0.0, 1.0),
            ],
            Some(Vector::new(1.0, 0.0, 0.0)),
        )?;
        let shade = Polygon::new(
            "shade",
            vec![
                Point::new(0.25, -2.0, -2.0),
                Point::new(0.25, 3.0, -2.0),
                Point::new(0.25, 3.0, 3.0),
                Point::new(0.25, -2.0, 3.0),
            ],
            Some(Vector::new(-1.0, 0.0, 0.0)),
        )?;

        let wall = Wall::new("window", vec![win])?;
        let wall2 = Wall::new("shade", vec![shade])?;
        let solid = Solid::new("room", vec![wall, wall2])?;
        let zone = Zone::new("z", vec![solid])?;
        let building = Building::new("b", vec![zone])?;

        let scene = FlatScene::new(&building, 0.25, true);
        let idx = scene
            .paths
            .iter()
            .position(|p| p.ends_with("/win"))
            .unwrap();
        let vis = visibility_fraction(&scene, idx, Vector::new(1.0, 0.0, 0.0));
        assert!(vis < 0.01, "vis={vis}");
        Ok(())
    }

    #[test]
    fn test_solar_shortwave_shaded_step_module_blocks_direct_component() -> Result<()> {
        fn weather_one_hour(dni: f64, dhi: f64) -> WeatherData {
            WeatherData {
                location: "x".to_string(),
                latitude: 0.0,
                longitude: 0.0,
                timezone: 0.0,
                elevation: 0.0,
                records: vec![HourlyRecord {
                    month: 3,
                    day: 21,
                    hour: 7, // near sunrise at equator: sun dir ~ +X
                    dry_bulb_temperature: 20.0,
                    relative_humidity: 50.0,
                    global_horizontal_radiation: dni + dhi,
                    direct_normal_radiation: dni,
                    diffuse_horizontal_radiation: dhi,
                    wind_speed: 0.0,
                    wind_direction: 0.0,
                }],
            }
        }

        fn make_building(with_shade: bool) -> Result<Building> {
            let win = Polygon::new(
                "glass",
                vec![
                    Point::new(0.0, 0.0, 0.0),
                    Point::new(0.0, 1.0, 0.0),
                    Point::new(0.0, 1.0, 1.0),
                    Point::new(0.0, 0.0, 1.0),
                ],
                Some(Vector::new(1.0, 0.0, 0.0)),
            )?;

            let mut walls = vec![Wall::new("window", vec![win])?];
            if with_shade {
                let shade = Polygon::new(
                    "panel",
                    vec![
                        Point::new(0.25, -2.0, -2.0),
                        Point::new(0.25, 3.0, -2.0),
                        Point::new(0.25, 3.0, 3.0),
                        Point::new(0.25, -2.0, 3.0),
                    ],
                    Some(Vector::new(-1.0, 0.0, 0.0)),
                )?;
                walls.push(Wall::new("shade", vec![shade])?);
            }

            let solid = Solid::new("room", walls)?;
            let zone = Zone::new("zone", vec![solid])?;
            Building::new("b", vec![zone])
        }

        let weather = weather_one_hour(1000.0, 0.0);
        let gain_config = SolarGainConfig::new();

        let b_clear = make_building(false)?;
        let idx_clear = SurfaceIndex::new(&b_clear);
        let ctx_clear = SimContext::new(&b_clear, &idx_clear);
        let mut bus_clear = Bus::new();
        let mut m_clear = SolarShortwaveShadedStepModule::new(SolarShortwaveShadedStepConfig {
            weather: weather.clone(),
            gain_config: gain_config.clone(),
            material_library: None,
            start_hour_idx: 0,
            voxel_size: 0.25,
        });
        m_clear.init(&ctx_clear, &mut bus_clear)?;
        m_clear.step(&ctx_clear, &mut bus_clear)?;
        let q_clear = bus_clear
            .get::<ShortwaveTransmittedWPerZone>()
            .unwrap()
            .watts_by_zone_uid
            .values()
            .next()
            .cloned()
            .unwrap_or(0.0);

        let b_shaded = make_building(true)?;
        let idx_shaded = SurfaceIndex::new(&b_shaded);
        let ctx_shaded = SimContext::new(&b_shaded, &idx_shaded);
        let mut bus_shaded = Bus::new();
        let mut m_shaded = SolarShortwaveShadedStepModule::new(SolarShortwaveShadedStepConfig {
            weather,
            gain_config,
            material_library: None,
            start_hour_idx: 0,
            voxel_size: 0.25,
        });
        m_shaded.init(&ctx_shaded, &mut bus_shaded)?;
        m_shaded.step(&ctx_shaded, &mut bus_shaded)?;
        let q_shaded = bus_shaded
            .get::<ShortwaveTransmittedWPerZone>()
            .unwrap()
            .watts_by_zone_uid
            .values()
            .next()
            .cloned()
            .unwrap_or(0.0);

        assert!(q_clear > 100.0, "q_clear={q_clear}");
        assert!(
            q_shaded < q_clear * 0.1,
            "q_shaded={q_shaded} q_clear={q_clear}"
        );
        Ok(())
    }

    #[test]
    fn test_solar_epw_shaded_module_blocks_direct_component() -> Result<()> {
        fn weather_one_hour(dni: f64, dhi: f64) -> WeatherData {
            WeatherData {
                location: "x".to_string(),
                latitude: 0.0,
                longitude: 0.0,
                timezone: 0.0,
                elevation: 0.0,
                records: vec![HourlyRecord {
                    month: 3,
                    day: 21,
                    hour: 7, // near sunrise at equator: sun dir ~ +X
                    dry_bulb_temperature: 20.0,
                    relative_humidity: 50.0,
                    global_horizontal_radiation: dni + dhi,
                    direct_normal_radiation: dni,
                    diffuse_horizontal_radiation: dhi,
                    wind_speed: 0.0,
                    wind_direction: 0.0,
                }],
            }
        }

        fn make_building(with_shade: bool) -> Result<Building> {
            let win = Polygon::new(
                "glass",
                vec![
                    Point::new(0.0, 0.0, 0.0),
                    Point::new(0.0, 1.0, 0.0),
                    Point::new(0.0, 1.0, 1.0),
                    Point::new(0.0, 0.0, 1.0),
                ],
                Some(Vector::new(1.0, 0.0, 0.0)),
            )?;

            let mut walls = vec![Wall::new("window", vec![win])?];
            if with_shade {
                let shade = Polygon::new(
                    "panel",
                    vec![
                        Point::new(0.25, -2.0, -2.0),
                        Point::new(0.25, 3.0, -2.0),
                        Point::new(0.25, 3.0, 3.0),
                        Point::new(0.25, -2.0, 3.0),
                    ],
                    Some(Vector::new(-1.0, 0.0, 0.0)),
                )?;
                walls.push(Wall::new("shade", vec![shade])?);
            }

            let solid = Solid::new("room", walls)?;
            let zone = Zone::new("zone", vec![solid])?;
            Building::new("b", vec![zone])
        }

        let weather = Arc::new(weather_one_hour(1000.0, 0.0));
        let gain_config = SolarGainConfig::new();

        let b_clear = make_building(false)?;
        let idx_clear = SurfaceIndex::new(&b_clear);
        let ctx_clear = SimContext::new(&b_clear, &idx_clear);
        let mut bus_clear = Bus::new();
        bus_clear.put(WeatherHourIndex(0));
        let mut m_clear = SolarEpwShadedModule::new(SolarEpwShadedConfig {
            weather: weather.clone(),
            gain_config: gain_config.clone(),
            material_library: None,
            voxel_size: 0.25,
        });
        m_clear.init(&ctx_clear, &mut bus_clear)?;
        m_clear.step(&ctx_clear, &mut bus_clear)?;
        let q_clear = bus_clear
            .get::<ShortwaveTransmittedWPerZone>()
            .unwrap()
            .watts_by_zone_uid
            .values()
            .next()
            .cloned()
            .unwrap_or(0.0);

        let b_shaded = make_building(true)?;
        let idx_shaded = SurfaceIndex::new(&b_shaded);
        let ctx_shaded = SimContext::new(&b_shaded, &idx_shaded);
        let mut bus_shaded = Bus::new();
        bus_shaded.put(WeatherHourIndex(0));
        let mut m_shaded = SolarEpwShadedModule::new(SolarEpwShadedConfig {
            weather,
            gain_config,
            material_library: None,
            voxel_size: 0.25,
        });
        m_shaded.init(&ctx_shaded, &mut bus_shaded)?;
        m_shaded.step(&ctx_shaded, &mut bus_shaded)?;
        let q_shaded = bus_shaded
            .get::<ShortwaveTransmittedWPerZone>()
            .unwrap()
            .watts_by_zone_uid
            .values()
            .next()
            .cloned()
            .unwrap_or(0.0);

        assert!(q_clear > 100.0, "q_clear={q_clear}");
        assert!(
            q_shaded < q_clear * 0.1,
            "q_shaded={q_shaded} q_clear={q_clear}"
        );
        Ok(())
    }

    #[test]
    fn test_solar_epw_shaded_module_produces_opaque_absorbed_per_polygon() -> Result<()> {
        // A single exterior opaque surface should receive absorbed shortwave.
        fn weather_one_hour(dni: f64, dhi: f64) -> WeatherData {
            WeatherData {
                location: "x".to_string(),
                latitude: 0.0,
                longitude: 0.0,
                timezone: 0.0,
                elevation: 0.0,
                records: vec![HourlyRecord {
                    month: 3,
                    day: 21,
                    hour: 12,
                    dry_bulb_temperature: 20.0,
                    relative_humidity: 50.0,
                    global_horizontal_radiation: dni + dhi,
                    direct_normal_radiation: dni,
                    diffuse_horizontal_radiation: dhi,
                    wind_speed: 0.0,
                    wind_direction: 0.0,
                }],
            }
        }

        // Opaque roof polygon (not glazing by name).
        let roof = Polygon::new(
            "roof",
            vec![
                Point::new(0.0, 0.0, 1.0),
                Point::new(1.0, 0.0, 1.0),
                Point::new(1.0, 1.0, 1.0),
                Point::new(0.0, 1.0, 1.0),
            ],
            Some(Vector::new(0.0, 0.0, 1.0)),
        )?;
        let wall = Wall::new("opaque", vec![roof])?;
        let solid = Solid::new("room", vec![wall])?;
        let zone = Zone::new("zone", vec![solid])?;
        let building = Building::new("b", vec![zone])?;

        let index = SurfaceIndex::new(&building);
        let ctx = SimContext::new(&building, &index);
        let weather = Arc::new(weather_one_hour(0.0, 200.0));

        let mut bus = Bus::new();
        bus.put(WeatherHourIndex(0));

        let mut module = SolarEpwShadedModule::new(SolarEpwShadedConfig {
            weather,
            gain_config: SolarGainConfig::new(),
            material_library: None,
            voxel_size: 0.25,
        });
        module.init(&ctx, &mut bus)?;
        module.step(&ctx, &mut bus)?;

        let absorbed = bus.get::<ShortwaveAbsorbedWPerPolygon>().unwrap();
        assert!(
            !absorbed.watts_by_polygon_uid.is_empty(),
            "expected absorbed shortwave on opaque exterior surface"
        );
        Ok(())
    }
}
