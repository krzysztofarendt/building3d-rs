use anyhow::Result;
use std::sync::Arc;

use crate::sim::coupling::WeatherHourIndex;
use crate::sim::coupling::{ShortwaveAbsorbedWPerPolygon, ShortwaveTransmittedWPerZone};
use crate::sim::energy::solar_bridge::{
    SolarGainConfig, SolarHourParams, compute_solar_gains_per_zone_with_materials,
};
use crate::sim::energy::weather::WeatherData;
use crate::sim::framework::{Bus, SimContext, SimModule};
use crate::sim::lighting::result::LightingResult;
use crate::sim::materials::{MaterialLibrary, OpticalMaterial};

use super::sources::Rgb;

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

        let mut transmitted = ShortwaveTransmittedWPerZone::default();
        transmitted.watts_by_zone_uid = gains_by_zone;

        // For now, this module does not attempt to apportion absorbed shortwave to surfaces.
        let absorbed = ShortwaveAbsorbedWPerPolygon::default();

        bus.put(transmitted);
        bus.put(absorbed);

        self.has_run = true;
        Ok(())
    }
}

/// EPW-driven (time-series) solar shortwave producer.
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
}

impl SolarEpwModule {
    pub fn new(config: SolarEpwConfig) -> Self {
        Self { config }
    }
}

impl SimModule for SolarEpwModule {
    fn name(&self) -> &'static str {
        "solar_epw_shortwave"
    }

    fn step(&mut self, ctx: &SimContext, bus: &mut Bus) -> Result<()> {
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
            direct_normal_irradiance: record.direct_normal_radiation,
            diffuse_horizontal_irradiance: record.diffuse_horizontal_radiation,
            day_of_year: day_of_year(record.month, record.day),
            hour: record.hour as f64,
            latitude: self.config.weather.latitude,
            longitude: self.config.weather.longitude,
        };

        let gains_by_zone = compute_solar_gains_per_zone_with_materials(
            ctx.building,
            &params,
            &self.config.gain_config,
            self.config.material_library.as_ref(),
        );

        let mut transmitted = ShortwaveTransmittedWPerZone::default();
        transmitted.watts_by_zone_uid = gains_by_zone;

        // For now, this module does not attempt to apportion absorbed shortwave to surfaces.
        let absorbed = ShortwaveAbsorbedWPerPolygon::default();

        bus.put(transmitted);
        bus.put(absorbed);
        Ok(())
    }
}

fn day_of_year(month: u8, day: u8) -> u16 {
    const DAYS_BEFORE_MONTH: [u16; 12] = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];
    let m = (month as usize).saturating_sub(1).min(11);
    DAYS_BEFORE_MONTH[m] + day as u16
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::index::SurfaceIndex;
    use crate::{Building, Point, Polygon, Solid, Wall, Zone};

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
            direct_normal_irradiance: 500.0,
            diffuse_horizontal_irradiance: 200.0,
            day_of_year: 80, // ~equinox
            hour: 12.0,
            latitude: 0.0,
            longitude: 0.0,
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
}
