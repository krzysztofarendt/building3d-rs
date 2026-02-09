use std::collections::HashMap;

use anyhow::Result;
use rand::Rng;

use crate::sim::engine::FlatScene;
use crate::sim::engine::reflection::{Diffuse, ReflectionModel, Specular};
use crate::{Building, Point, Vector};

use super::config::LightingConfig;
use super::result::LightingResult;
use super::sensor::SensorGrid;
use super::sources::{LightSource, Rgb};

/// Forward ray tracing lighting simulation.
pub struct LightingSimulation {
    config: LightingConfig,
    scene: FlatScene,
    diffuse_reflectance: Vec<Rgb>,
    specular_reflectance: Vec<Rgb>,
    transmittance: Vec<Rgb>,
    polygon_uids: Vec<crate::UID>,
    areas_m2: Vec<f64>,
    /// Maps polygon index → sensor grid index (only for sensor-equipped polygons).
    sensor_map: HashMap<usize, usize>,
    /// Sensor grids generated for matching polygons.
    sensor_grids: Vec<SensorGrid>,
    /// Area per sensor cell (spacing^2).
    sensor_area: f64,
}

impl LightingSimulation {
    pub fn new(building: &Building, config: LightingConfig) -> Result<Self> {
        let scene = FlatScene::new(building, config.voxel_size, false);

        // Resolve optical properties per polygon
        let mut diffuse_reflectance: Vec<Rgb> = Vec::with_capacity(scene.paths.len());
        let mut specular_reflectance: Vec<Rgb> = Vec::with_capacity(scene.paths.len());
        let mut transmittance: Vec<Rgb> = Vec::with_capacity(scene.paths.len());

        for path in &scene.paths {
            let optical = config
                .material_library
                .as_ref()
                .and_then(|lib| lib.lookup(path))
                .and_then(|mat| mat.optical.as_ref());

            if let Some(opt) = optical {
                diffuse_reflectance.push(opt.diffuse_reflectance);
                specular_reflectance.push(opt.specular_reflectance);
                transmittance.push(opt.transmittance);
            } else {
                diffuse_reflectance.push(config.default_reflectance);
                specular_reflectance.push([0.0; 3]);
                transmittance.push([0.0; 3]);
            }
        }

        let polygon_uids: Vec<crate::UID> = scene.polygons.iter().map(|p| p.uid.clone()).collect();
        let areas_m2: Vec<f64> = scene.polygons.iter().map(|p| p.area()).collect();

        // Generate sensor grids if configured
        let mut sensor_map = HashMap::new();
        let mut sensor_grids = Vec::new();
        let sensor_area = config.sensor_spacing.map(|s| s * s).unwrap_or(1.0);

        if let Some(spacing) = config.sensor_spacing {
            for (i, poly) in scene.polygons.iter().enumerate() {
                let path = &scene.paths[i];
                let matches = config.sensor_patterns.is_empty()
                    || config
                        .sensor_patterns
                        .iter()
                        .any(|pat| path.contains(pat.as_str()));
                if matches {
                    let grid = SensorGrid::generate(poly, spacing, path);
                    if !grid.sensors.is_empty() {
                        sensor_map.insert(i, sensor_grids.len());
                        sensor_grids.push(grid);
                    }
                }
            }
        }

        Ok(Self {
            config,
            scene,
            diffuse_reflectance,
            specular_reflectance,
            transmittance,
            polygon_uids,
            areas_m2,
            sensor_map,
            sensor_grids,
            sensor_area,
        })
    }

    /// Runs the forward ray tracing simulation.
    pub fn run(&self) -> LightingResult {
        let mut result = LightingResult::new();
        result.sensor_grids = self.sensor_grids.clone();
        let mut rng = rand::thread_rng();

        // Trace rays from each point light
        for light in &self.config.point_lights {
            let rays_per_light = self.config.num_rays;
            // Each ray carries an equal share of the total flux
            let scale = 1.0 / rays_per_light as f64;

            for _ in 0..rays_per_light {
                let dir = random_unit_vector(&mut rng);
                let intensity = light.intensity(dir);
                let energy = [
                    intensity[0] * scale,
                    intensity[1] * scale,
                    intensity[2] * scale,
                ];

                self.trace_ray(light.position(), dir, energy, &mut result, &mut rng);
            }
        }

        // Trace rays from directional lights
        for light in &self.config.directional_lights {
            let rays_per_light = self.config.num_rays;
            let irradiance = light.intensity(light.direction);

            // Spawn rays from the top of the bounding box in the light direction
            let dir_norm = match light.direction.normalize() {
                Ok(v) => v,
                Err(_) => continue,
            };

            for _ in 0..rays_per_light {
                // Random point on the scene bounding box face
                let origin = random_point_on_bbox_face(&self.scene, dir_norm, &mut rng);
                let energy_per_ray = [
                    irradiance[0] / rays_per_light as f64,
                    irradiance[1] / rays_per_light as f64,
                    irradiance[2] / rays_per_light as f64,
                ];

                self.trace_ray(origin, dir_norm, energy_per_ray, &mut result, &mut rng);
            }
        }

        result.finalize(&self.polygon_uids, &self.areas_m2);
        result
    }

    fn trace_ray(
        &self,
        origin: Point,
        direction: Vector,
        energy: Rgb,
        result: &mut LightingResult,
        rng: &mut impl Rng,
    ) {
        let mut pos = origin;
        let mut dir = direction;
        let mut e = energy;

        for _bounce in 0..self.config.max_bounces {
            let total_e = e[0] + e[1] + e[2];
            if total_e < self.config.min_energy {
                break;
            }

            if let Some((idx, dist)) = self.scene.find_target_surface_global(pos, dir) {
                // Record hit
                result.record_hit_index(idx, e);

                // Record sensor hit if this polygon has a sensor grid
                if let Some(&grid_idx) = self.sensor_map.get(&idx) {
                    let hit_pos = pos + dir * dist;
                    result.record_sensor_hit(grid_idx, hit_pos, e, self.sensor_area);
                }

                // Look up per-polygon optical properties
                let diff = self.diffuse_reflectance[idx];
                let spec = self.specular_reflectance[idx];
                let trans = self.transmittance[idx];

                // Compute selection probabilities from max channel values
                let p_diff = diff[0].max(diff[1]).max(diff[2]);
                let p_spec = spec[0].max(spec[1]).max(spec[2]);
                let p_trans = trans[0].max(trans[1]).max(trans[2]);
                let total = p_diff + p_spec + p_trans;

                if total < 1e-10 {
                    break; // Fully absorbed
                }

                // Random roll to pick interaction type
                let roll: f64 = rng.gen_range(0.0..total);
                let normal = self.scene.polygons[idx].vn;
                pos = pos + dir * dist;

                if roll < p_diff {
                    // Diffuse reflection
                    e = [
                        e[0] * diff[0] / p_diff,
                        e[1] * diff[1] / p_diff,
                        e[2] * diff[2] / p_diff,
                    ];
                    dir = Diffuse.reflect(dir, normal);
                } else if roll < p_diff + p_spec {
                    // Specular reflection
                    e = [
                        e[0] * spec[0] / p_spec,
                        e[1] * spec[1] / p_spec,
                        e[2] * spec[2] / p_spec,
                    ];
                    dir = Specular.reflect(dir, normal);
                } else {
                    // Transmission: continue in same direction, offset past surface
                    e = [
                        e[0] * trans[0] / p_trans,
                        e[1] * trans[1] / p_trans,
                        e[2] * trans[2] / p_trans,
                    ];
                    pos = pos + dir * 1e-6;
                    // dir unchanged
                }

                // Russian roulette for path termination
                let max_channel = e[0].max(e[1]).max(e[2]);
                if max_channel < 1e-10 {
                    break;
                }
                let survival = max_channel.min(1.0);
                let rr: f64 = rng.gen_range(0.0..1.0);
                if rr > survival {
                    break;
                }
                e = [e[0] / survival, e[1] / survival, e[2] / survival];
            } else {
                break; // Ray escapes scene
            }
        }
    }
}

/// Generate a random unit vector.
fn random_unit_vector(rng: &mut impl Rng) -> Vector {
    loop {
        let x: f64 = rng.gen_range(-1.0..1.0);
        let y: f64 = rng.gen_range(-1.0..1.0);
        let z: f64 = rng.gen_range(-1.0..1.0);
        let len2 = x * x + y * y + z * z;
        if len2 > 1e-6 && len2 <= 1.0 {
            let len = len2.sqrt();
            return Vector::new(x / len, y / len, z / len);
        }
    }
}

/// Generate a random point on the bounding box face facing the given direction.
fn random_point_on_bbox_face(scene: &FlatScene, direction: Vector, rng: &mut impl Rng) -> Point {
    let min = scene.bbox_min;
    let max = scene.bbox_max;
    let margin = 0.01;

    // Choose the face opposite to the light direction
    let abs_x = direction.dx.abs();
    let abs_y = direction.dy.abs();
    let abs_z = direction.dz.abs();

    if abs_z >= abs_x && abs_z >= abs_y {
        // Z-dominant: spawn on top or bottom
        let z = if direction.dz < 0.0 {
            max.z + margin
        } else {
            min.z - margin
        };
        let x = rng.gen_range(min.x..max.x);
        let y = rng.gen_range(min.y..max.y);
        Point::new(x, y, z)
    } else if abs_y >= abs_x {
        let y = if direction.dy < 0.0 {
            max.y + margin
        } else {
            min.y - margin
        };
        let x = rng.gen_range(min.x..max.x);
        let z = rng.gen_range(min.z..max.z);
        Point::new(x, y, z)
    } else {
        let x = if direction.dx < 0.0 {
            max.x + margin
        } else {
            min.x - margin
        };
        let y = rng.gen_range(min.y..max.y);
        let z = rng.gen_range(min.z..max.z);
        Point::new(x, y, z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::index::SurfaceIndex;
    use crate::sim::lighting::sources::DirectionalLight;
    use crate::sim::lighting::sources::PointLight;
    use crate::{Solid, Zone};
    use rand::SeedableRng;

    #[test]
    fn test_lighting_simulation_basic() {
        let s0 = Solid::from_box(4.0, 4.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config = LightingConfig::new();
        config.num_rays = 1000;
        config.max_bounces = 3;
        config
            .point_lights
            .push(PointLight::white(Point::new(2.0, 2.0, 2.5), 1000.0));

        let sim = LightingSimulation::new(&building, config).unwrap();
        let result = sim.run();

        // Some surfaces should have been hit
        assert!(!result.hit_count.is_empty(), "Some surfaces should be hit");
        assert!(
            !result.illuminance.is_empty(),
            "Illuminance should be computed"
        );
    }

    #[test]
    fn test_inverse_square_falloff() {
        // Verify 1/r^2 by comparing total power hitting opposing walls at
        // different distances from a point light.
        //
        // In a 4x4x10 box with light at (2, 2, 2), the floor (z=0) is at
        // distance 2 from the light, and the ceiling (z=10) is at distance 8.
        // Both are 4x4 = 16 m^2.
        //
        // For a centered point light hitting identical-size parallel surfaces:
        //   P_floor / P_ceiling ≈ (r_ceiling / r_floor)^2 = (8/2)^2 = 16
        //
        // This is approximate because the surfaces are finite, not infinitesimal.
        let s0 = Solid::from_box(4.0, 4.0, 10.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config = LightingConfig::new();
        config.num_rays = 500_000;
        config.max_bounces = 1; // first bounce only
        config.default_reflectance = [0.0, 0.0, 0.0]; // absorb all light
        config
            .point_lights
            .push(PointLight::white(Point::new(2.0, 2.0, 2.0), 3000.0));

        let sim = LightingSimulation::new(&building, config).unwrap();
        let result = sim.run();

        // Compare total flux (not irradiance) on floor vs ceiling
        let mut floor_flux = 0.0;
        let mut ceiling_flux = 0.0;

        let index = SurfaceIndex::new(&building);
        for (polygon_uid, flux) in &result.incident_flux {
            let path = index.path_by_polygon_uid(polygon_uid).unwrap_or("");
            let total = flux[0] + flux[1] + flux[2];
            if path.contains("floor") {
                floor_flux += total;
            } else if path.contains("ceiling") {
                ceiling_flux += total;
            }
        }

        // Both surfaces receive some flux
        assert!(floor_flux > 0.0, "Floor should receive flux");
        assert!(ceiling_flux > 0.0, "Ceiling should receive flux");

        // The floor (r=2) should receive more flux than ceiling (r=8)
        let ratio = floor_flux / ceiling_flux;

        // Expected ratio ≈ 16 for point source. Due to finite surface effects,
        // the actual ratio will be lower. Accept range 5-50.
        assert!(
            ratio > 5.0,
            "Floor should receive much more flux than ceiling. \
             ratio={ratio:.2} (floor={floor_flux:.2}, ceiling={ceiling_flux:.2})"
        );
    }

    #[test]
    fn test_lighting_no_lights() {
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let config = LightingConfig::new();
        let sim = LightingSimulation::new(&building, config).unwrap();
        let result = sim.run();

        assert!(result.hit_count.is_empty(), "No lights means no hits");
    }

    #[test]
    fn test_random_unit_vector_is_unit_length() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(123);
        let v = random_unit_vector(&mut rng);
        let len = v.length();
        assert!((len - 1.0).abs() < 1e-10, "len={len}");
    }

    #[test]
    fn test_random_point_on_bbox_face_spawns_on_expected_face() {
        let s0 = Solid::from_box(2.0, 3.0, 4.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();
        let scene = FlatScene::new(&building, 0.5, false);

        let mut rng = rand::rngs::StdRng::seed_from_u64(999);

        // Z-dominant: negative dz -> spawn on top (max.z + margin).
        let p = random_point_on_bbox_face(&scene, Vector::new(0.0, 0.0, -1.0), &mut rng);
        assert!(p.z > scene.bbox_max.z);

        // Y-dominant: positive dy -> spawn on min.y - margin.
        let p = random_point_on_bbox_face(&scene, Vector::new(0.0, 2.0, 0.0), &mut rng);
        assert!(p.y < scene.bbox_min.y);

        // X-dominant: negative dx -> spawn on max.x + margin.
        let p = random_point_on_bbox_face(&scene, Vector::new(-3.0, 0.0, 0.0), &mut rng);
        assert!(p.x > scene.bbox_max.x);
    }

    #[test]
    fn test_directional_light_zero_direction_is_skipped() {
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config = LightingConfig::new();
        config.num_rays = 100;
        config.max_bounces = 1;
        config.directional_lights.push(DirectionalLight::new(
            Vector::new(0.0, 0.0, 0.0),
            [100.0; 3],
        ));

        let sim = LightingSimulation::new(&building, config).unwrap();
        let result = sim.run();
        assert!(result.hit_count.is_empty());
    }

    #[test]
    fn test_sensor_grids_receive_hits_from_directional_light() {
        let s0 = Solid::from_box(4.0, 4.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config = LightingConfig::new();
        config.num_rays = 500;
        config.max_bounces = 1;
        config.default_reflectance = [0.0; 3];
        config.sensor_spacing = Some(1.0);
        config.sensor_patterns = vec!["floor".to_string()];
        config.directional_lights.push(DirectionalLight::new(
            Vector::new(0.0, 0.0, 1.0),
            [300.0; 3],
        ));

        let sim = LightingSimulation::new(&building, config).unwrap();
        let result = sim.run();

        assert!(
            !result.sensor_grids.is_empty(),
            "Expected at least one generated sensor grid"
        );
        let total_sensor_hits: usize = result
            .sensor_grids
            .iter()
            .flat_map(|g| g.sensors.iter())
            .map(|s| s.hit_count)
            .sum();
        assert!(
            total_sensor_hits > 0,
            "Expected at least one sensor hit from directional rays"
        );
    }

    // ── Physics verification tests ──────────────────────────────────────

    #[test]
    fn test_energy_conservation_first_bounce() {
        // With fully absorbing walls (reflectance=0) and one bounce,
        // all rays hit exactly one surface. In a closed box, total collected
        // flux should be consistent across two runs with the same setup.
        //
        // The PointLight emits intensity in W/sr; ray energy = intensity/N_rays.
        // Total collected = source_flux / (4*pi) for an omnidirectional source.
        let s0 = Solid::from_box(4.0, 4.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let source_flux = 1000.0;
        let expected = source_flux / (4.0 * std::f64::consts::PI);
        let mut config = LightingConfig::new();
        config.num_rays = 200_000;
        config.max_bounces = 1;
        config.default_reflectance = [0.0, 0.0, 0.0]; // fully absorbing
        config
            .point_lights
            .push(PointLight::white(Point::new(2.0, 2.0, 1.5), source_flux));

        let sim = LightingSimulation::new(&building, config).unwrap();
        let result = sim.run();

        let total_collected: f64 = result
            .incident_flux
            .values()
            .map(|f| f[0] + f[1] + f[2])
            .sum();

        // Total should be close to source_flux/(4*pi) ≈ 79.6 W.
        // Allow 5% tolerance for ray sampling variance.
        let ratio = total_collected / expected;
        assert!(
            ratio > 0.95 && ratio < 1.05,
            "Total collected flux should approximate source_flux/(4*pi). \
             ratio={ratio:.3} (collected={total_collected:.1}, expected={expected:.1})"
        );
    }

    #[test]
    fn test_more_bounces_more_total_flux() {
        // With reflective walls, more bounces should lead to more total
        // flux being deposited across surfaces (energy gets redistributed).
        let s0 = Solid::from_box(3.0, 3.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let source_flux = 1000.0;

        let run_with_bounces = |max_bounces: usize| -> f64 {
            let mut config = LightingConfig::new();
            config.num_rays = 100_000;
            config.max_bounces = max_bounces;
            config.default_reflectance = [0.5, 0.5, 0.5];
            config
                .point_lights
                .push(PointLight::white(Point::new(1.5, 1.5, 1.5), source_flux));

            let sim = LightingSimulation::new(&building, config).unwrap();
            let result = sim.run();
            result
                .incident_flux
                .values()
                .map(|f| f[0] + f[1] + f[2])
                .sum()
        };

        let flux_1 = run_with_bounces(1);
        let flux_5 = run_with_bounces(5);

        assert!(
            flux_5 > flux_1,
            "More bounces should deposit more total flux: 5-bounce={flux_5:.1}, 1-bounce={flux_1:.1}"
        );
    }

    #[test]
    fn test_double_source_power_doubles_illuminance() {
        // Illuminance should scale linearly with source power.
        let s0 = Solid::from_box(3.0, 3.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let run_with_flux = |flux: f64| -> f64 {
            let mut config = LightingConfig::new();
            config.num_rays = 200_000;
            config.max_bounces = 1;
            config.default_reflectance = [0.0, 0.0, 0.0];
            config
                .point_lights
                .push(PointLight::white(Point::new(1.5, 1.5, 1.5), flux));

            let sim = LightingSimulation::new(&building, config).unwrap();
            let result = sim.run();
            result
                .incident_flux
                .values()
                .map(|f| f[0] + f[1] + f[2])
                .sum()
        };

        let flux_500 = run_with_flux(500.0);
        let flux_1000 = run_with_flux(1000.0);

        let ratio = flux_1000 / flux_500;
        assert!(
            (ratio - 2.0).abs() < 0.2,
            "Doubling source power should double collected flux. \
             ratio={ratio:.3} (1000W={flux_1000:.1}, 500W={flux_500:.1})"
        );
    }

    #[test]
    fn test_reflective_walls_increase_total_flux() {
        // Reflective walls should increase total deposited flux compared
        // to fully absorbing walls (energy bounces and gets recorded multiple times).
        let s0 = Solid::from_box(3.0, 3.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let source_flux = 1000.0;

        let run_with_reflectance = |refl: f64| -> f64 {
            let mut config = LightingConfig::new();
            config.num_rays = 100_000;
            config.max_bounces = 10;
            config.default_reflectance = [refl, refl, refl];
            config
                .point_lights
                .push(PointLight::white(Point::new(1.5, 1.5, 1.5), source_flux));

            let sim = LightingSimulation::new(&building, config).unwrap();
            let result = sim.run();
            result
                .incident_flux
                .values()
                .map(|f| f[0] + f[1] + f[2])
                .sum()
        };

        let flux_absorbing = run_with_reflectance(0.0);
        let flux_reflective = run_with_reflectance(0.8);

        assert!(
            flux_reflective > flux_absorbing * 1.5,
            "Reflective walls should deposit significantly more total flux. \
             reflective={flux_reflective:.1}, absorbing={flux_absorbing:.1}"
        );
    }
}
