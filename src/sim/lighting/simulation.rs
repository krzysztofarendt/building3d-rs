use std::collections::HashMap;

use anyhow::Result;
use rand::Rng;

use crate::sim::engine::FlatScene;
use crate::sim::engine::reflection::Diffuse;
use crate::sim::engine::reflection::ReflectionModel;
use crate::{Building, Point, Vector};

use super::config::LightingConfig;
use super::result::LightingResult;
use super::sources::{LightSource, Rgb};

/// Forward ray tracing lighting simulation.
pub struct LightingSimulation {
    config: LightingConfig,
    scene: FlatScene,
    reflectance: Vec<Rgb>,
    areas: HashMap<String, f64>,
}

impl LightingSimulation {
    pub fn new(building: &Building, config: LightingConfig) -> Result<Self> {
        let scene = FlatScene::new(building, config.voxel_size, false);

        // Resolve reflectance per polygon
        let reflectance: Vec<Rgb> = scene
            .paths
            .iter()
            .map(|path| {
                if let Some(optical) = config
                    .material_library
                    .as_ref()
                    .and_then(|lib| lib.lookup(path))
                    .and_then(|mat| mat.optical.as_ref())
                {
                    return optical.diffuse_reflectance;
                }
                config.default_reflectance
            })
            .collect();

        // Compute polygon areas
        let mut areas = HashMap::new();
        for (i, poly) in scene.polygons.iter().enumerate() {
            areas.insert(scene.paths[i].clone(), poly.area());
        }

        Ok(Self {
            config,
            scene,
            reflectance,
            areas,
        })
    }

    /// Runs the forward ray tracing simulation.
    pub fn run(&self) -> LightingResult {
        let mut result = LightingResult::new();
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

        result.compute_illuminance(&self.areas);
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
        let reflection_model = Diffuse;

        for _bounce in 0..self.config.max_bounces {
            let total_e = e[0] + e[1] + e[2];
            if total_e < self.config.min_energy {
                break;
            }

            if let Some((idx, dist)) = self.scene.find_target_surface_global(pos, dir) {
                // Record hit
                result.record_hit(&self.scene.paths[idx], e);

                // Reflect with diffuse reflectance
                let refl = self.reflectance[idx];
                e = [e[0] * refl[0], e[1] * refl[1], e[2] * refl[2]];

                // Russian roulette for path termination
                let max_refl = refl[0].max(refl[1]).max(refl[2]);
                if max_refl < 1e-10 {
                    break;
                }
                let rr: f64 = rng.gen_range(0.0..1.0);
                if rr > max_refl {
                    break;
                }
                // Compensate for Russian roulette
                e = [e[0] / max_refl, e[1] / max_refl, e[2] / max_refl];

                // Move to hit point and reflect
                let normal = self.scene.polygons[idx].vn;
                pos = pos + dir * dist;
                dir = reflection_model.reflect(dir, normal);
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
    use crate::sim::lighting::sources::PointLight;
    use crate::{Solid, Zone};

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
    fn test_lighting_no_lights() {
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let config = LightingConfig::new();
        let sim = LightingSimulation::new(&building, config).unwrap();
        let result = sim.run();

        assert!(result.hit_count.is_empty(), "No lights means no hits");
    }
}
