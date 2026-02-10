use std::sync::Arc;

use anyhow::Result;

use building3d::sim::coupling::ShortwaveTransmittedWPerZone;
use building3d::sim::energy::gains_module::{InternalGainsModule, InternalGainsModuleConfig};
use building3d::sim::energy::module::{EnergyModule, EnergyModuleConfig};
use building3d::sim::energy::schedule::InternalGainsProfile;
use building3d::sim::energy::weather::HourlyRecord;
use building3d::sim::energy::weather_module::{WeatherModule, WeatherModuleConfig};
use building3d::sim::framework::{Bus, Pipeline, SimContext};
use building3d::sim::index::SurfaceIndex;
use building3d::sim::lighting::shortwave::{
    SolarEpwBusConfig, SolarEpwBusModule, SolarEpwShadedBusConfig, SolarEpwShadedBusModule,
};
use building3d::sim::materials::MaterialLibrary;
use building3d::{Building, Solid, Zone};

fn main() -> Result<()> {
    // A small "room" plus a thin exterior shading slab in front of the east wall.
    let room = Solid::from_box(2.0, 2.0, 2.0, None, "room")?;
    let shade = Solid::from_box(0.02, 4.0, 4.0, Some((2.20, -1.0, -1.0)), "shade")?;
    let zone = Zone::new("z0", vec![room, shade])?;
    let building = Building::new("b", vec![zone])?;

    let surface_index = SurfaceIndex::new(&building);
    let ctx = SimContext::new(&building, &surface_index);

    // Mark the east wall polygon of the room (`room/wall_1/poly_1`) as glazing.
    let mut materials = MaterialLibrary::with_presets();
    materials.assign("room/wall_1/poly_1", "glass");
    let materials = Some(materials);

    // One EPW-like record: equator, equinox-ish, early morning => sun direction ~ +X.
    let weather = Arc::new(building3d::sim::energy::weather::WeatherData {
        location: "synthetic".to_string(),
        latitude: 0.0,
        longitude: 0.0,
        timezone: 0.0,
        elevation: 0.0,
        records: vec![HourlyRecord {
            month: 3,
            day: 21,
            hour: 7,
            dry_bulb_temperature: 20.0,
            relative_humidity: 50.0,
            global_horizontal_radiation: 1000.0,
            direct_normal_radiation: 1000.0,
            diffuse_horizontal_radiation: 0.0,
            horizontal_infrared_radiation: 300.0,
            wind_speed: 0.0,
            wind_direction: 0.0,
        }],
    });

    let gains_profile = Arc::new(InternalGainsProfile::office(0.0));

    let mut energy_cfg = EnergyModuleConfig::default();
    // For this demo, isolate the solar coupling effect.
    energy_cfg.thermal.default_u_value = 0.0;
    energy_cfg.thermal.infiltration_ach = 0.0;
    energy_cfg.hvac = building3d::sim::energy::hvac::HvacIdealLoads::with_setpoints(-1e9, 1e9);

    let mut bus_unshaded = Bus::new();
    let mut pipeline_unshaded = Pipeline::new()
        .with_module(WeatherModule::new(WeatherModuleConfig::new(
            weather.clone(),
        )))
        .with_module(InternalGainsModule::new(InternalGainsModuleConfig::new(
            gains_profile.clone(),
        )))
        .with_module(SolarEpwBusModule::new(SolarEpwBusConfig {
            gain_config: building3d::sim::energy::solar_bridge::SolarGainConfig::new(),
            material_library: materials.clone(),
        }))
        .with_module(EnergyModule::new(energy_cfg.clone()));
    pipeline_unshaded.init(&ctx, &mut bus_unshaded)?;
    pipeline_unshaded.step(&ctx, &mut bus_unshaded)?;

    let q_unshaded = bus_unshaded
        .get::<ShortwaveTransmittedWPerZone>()
        .and_then(|q| q.watts_by_zone_uid.values().next().cloned())
        .unwrap_or(0.0);

    let mut bus_shaded = Bus::new();
    let mut pipeline_shaded = Pipeline::new()
        .with_module(WeatherModule::new(WeatherModuleConfig::new(
            weather.clone(),
        )))
        .with_module(InternalGainsModule::new(InternalGainsModuleConfig::new(
            gains_profile,
        )))
        .with_module(SolarEpwShadedBusModule::new(SolarEpwShadedBusConfig {
            gain_config: building3d::sim::energy::solar_bridge::SolarGainConfig::new(),
            material_library: materials,
            voxel_size: 0.25,
        }))
        .with_module(EnergyModule::new(energy_cfg));
    pipeline_shaded.init(&ctx, &mut bus_shaded)?;
    pipeline_shaded.step(&ctx, &mut bus_shaded)?;

    let q_shaded = bus_shaded
        .get::<ShortwaveTransmittedWPerZone>()
        .and_then(|q| q.watts_by_zone_uid.values().next().cloned())
        .unwrap_or(0.0);

    println!("unshaded shortwave transmitted: {:.1} W", q_unshaded);
    println!("shaded shortwave transmitted:   {:.1} W", q_shaded);
    Ok(())
}
