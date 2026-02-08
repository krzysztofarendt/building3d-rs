use anyhow::Result;

use crate::sim::coupling::{ShortwaveAbsorbedWPerPolygon, ShortwaveTransmittedWPerZone};
use crate::sim::framework::{Bus, SimContext, SimModule};
use crate::sim::materials::OpticalMaterial;

use super::config::LightingConfig;
use super::result::LightingResult;
use super::simulation::LightingSimulation;
use super::sources::Rgb;

/// Pipeline wrapper for running lighting simulations in a composed workflow.
///
/// On the first `step()` call, this module:
/// - runs the lighting simulation (forward tracer),
/// - publishes `LightingResult`,
/// - publishes shortwave coupling payloads keyed by `UID`.
///
/// Notes:
/// - The coupling payloads are computed from the simulation’s total incident flux and
///   the assigned `OpticalMaterial` reflectance/transmittance fractions.
/// - Callers should treat the produced payloads as "solar shortwave" only when the
///   configured lighting sources represent solar radiation (e.g., directional sun/sky),
///   not arbitrary artificial lights.
pub struct LightingModule {
    config: LightingConfig,
    has_run: bool,
    simulation: Option<LightingSimulation>,
}

impl LightingModule {
    pub fn new(config: LightingConfig) -> Self {
        Self {
            config,
            has_run: false,
            simulation: None,
        }
    }
}

impl SimModule for LightingModule {
    fn name(&self) -> &'static str {
        "lighting"
    }

    fn init(&mut self, ctx: &SimContext, _bus: &mut Bus) -> Result<()> {
        self.simulation = Some(LightingSimulation::new(ctx.building, self.config.clone())?);
        Ok(())
    }

    fn step(&mut self, ctx: &SimContext, bus: &mut Bus) -> Result<()> {
        if self.has_run {
            return Ok(());
        }

        let Some(sim) = self.simulation.as_ref() else {
            self.simulation = Some(LightingSimulation::new(ctx.building, self.config.clone())?);
            return self.step(ctx, bus);
        };

        let result = sim.run();
        let (absorbed, transmitted) = shortwave_payloads_from_result(&result, ctx, &self.config);

        bus.put(result);
        bus.put(absorbed);
        bus.put(transmitted);

        self.has_run = true;
        Ok(())
    }
}

fn shortwave_payloads_from_result(
    result: &LightingResult,
    ctx: &SimContext,
    config: &LightingConfig,
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::index::SurfaceIndex;
    use crate::{Building, Point, Polygon, Solid, Wall, Zone};

    #[test]
    fn test_shortwave_payloads_keyed_by_uid() -> Result<()> {
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

        // Get polygon UID for the single surface.
        let surface = index.surfaces.first().unwrap();
        let polygon_uid = surface.polygon_uid.clone();

        let mut result = LightingResult::new();
        result
            .incident_flux
            .insert(polygon_uid.clone(), [10.0, 10.0, 10.0]);

        let mut config = LightingConfig::new();
        // Diffuse reflectance 0.5 => absorbance 0.5 (no spec/trans)
        config.default_reflectance = [0.5, 0.5, 0.5];

        let (absorbed, transmitted) = shortwave_payloads_from_result(&result, &ctx, &config);

        assert_eq!(absorbed.watts_by_polygon_uid.len(), 1);
        let a = absorbed.watts_by_polygon_uid.get(&polygon_uid).unwrap();
        // absorbed = sum(10 * 0.5) = 15 W across RGB
        assert!((a - 15.0).abs() < 1e-10);

        assert!(transmitted.watts_by_zone_uid.is_empty());
        Ok(())
    }
}
