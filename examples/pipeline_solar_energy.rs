use anyhow::Result;

use building3d::sim::energy::module::{EnergyModule, EnergyModuleConfig};
use building3d::sim::energy::weather::WeatherData;
use building3d::sim::framework::{Bus, Pipeline, SimContext};
use building3d::sim::index::SurfaceIndex;
use building3d::sim::lighting::shortwave::{SolarShortwaveStepConfig, SolarShortwaveStepModule};
use building3d::sim::materials::MaterialLibrary;
use building3d::{Building, Solid, Zone};

fn main() -> Result<()> {
    // Two adjacent zones to demonstrate inter-zone coupling.
    let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "s0")?;
    let s1 = Solid::from_box(2.0, 2.0, 2.0, Some((2.0, 0.0, 0.0)), "s1")?;
    let z0 = Zone::new("z0", vec![s0])?;
    let z1 = Zone::new("z1", vec![s1])?;
    let building = Building::new("b", vec![z0, z1])?;

    let surface_index = SurfaceIndex::new(&building);
    let ctx = SimContext::new(&building, &surface_index);
    let mut bus = Bus::new();

    // Synthetic weather with a simple solar profile.
    let weather = WeatherData::synthetic("synthetic", 0.0, 0.0, 15.0, 0.0);

    // Mark all surfaces in zone z0 as "glass" so the SHGC approximation produces non-zero gains.
    let mut materials = MaterialLibrary::with_presets();
    materials.assign("z0", "glass");

    let solar_cfg = building3d::sim::energy::solar_bridge::SolarGainConfig::new();
    let mut solar_step_cfg = SolarShortwaveStepConfig::new(weather, solar_cfg);
    solar_step_cfg.material_library = Some(materials);
    solar_step_cfg.start_hour_idx = 0;

    let solar = SolarShortwaveStepModule::new(solar_step_cfg);

    let mut energy_cfg = EnergyModuleConfig::default();
    energy_cfg.thermal.outdoor_temperature = 15.0;

    let energy = EnergyModule::new(energy_cfg);

    let mut pipeline = Pipeline::new().with_module(solar).with_module(energy);
    pipeline.init(&ctx, &mut bus)?;

    // Advance 24 hours.
    for step in 0..24 {
        pipeline.step(&ctx, &mut bus)?;
        let out = bus
            .get::<building3d::sim::energy::network::MultiZoneStepResult>()
            .unwrap();
        println!(
            "step {:02}: T = [{:.2}, {:.2}] °C, HVAC heat=[{:.0},{:.0}] W, cool=[{:.0},{:.0}] W",
            step,
            out.zone_temperatures_c[0],
            out.zone_temperatures_c[1],
            out.zone_heating_w[0],
            out.zone_heating_w[1],
            out.zone_cooling_w[0],
            out.zone_cooling_w[1],
        );
    }

    Ok(())
}
