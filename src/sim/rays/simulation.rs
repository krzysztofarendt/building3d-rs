use anyhow::Result;
use rayon::prelude::*;

use crate::sim::engine::absorption::{
    AbsorptionModel, AirAbsorption, FrequencyDependentAbsorption, ScalarAbsorption,
};
use crate::sim::engine::propagation::{FixedTimeStep, PropagationModel};
use crate::sim::engine::reflection::{Hybrid, ReflectionModel, Specular};
use crate::sim::engine::{FlatScene, RayBatch};
use crate::sim::materials::NUM_OCTAVE_BANDS;
use crate::{Building, Point, Vector};

use super::config::{AcousticMode, SimulationConfig};

/// Result of a ray tracing simulation.
pub struct SimulationResult {
    /// Ray positions per time step: positions[step][ray]
    pub positions: Vec<Vec<Point>>,
    /// Ray energies per time step: energies[step][ray] (scalar: sum across bands)
    pub energies: Vec<Vec<f64>>,
    /// Frequency-band energies per time step (if frequency-dependent mode).
    /// band_energies[step][ray] is [f64; NUM_OCTAVE_BANDS].
    pub band_energies: Option<Vec<Vec<[f64; NUM_OCTAVE_BANDS]>>>,
    /// Absorber hits per time step: hits[step][absorber]
    pub hits: Vec<Vec<f64>>,
    /// Frequency-band absorber hits (if frequency-dependent mode).
    /// band_hits[step][absorber] is [f64; NUM_OCTAVE_BANDS].
    pub band_hits: Option<Vec<Vec<[f64; NUM_OCTAVE_BANDS]>>>,
    /// Configuration used for this simulation
    pub config: SimulationConfig,
}

pub struct Simulation {
    config: SimulationConfig,
    scene: FlatScene,
}

impl Simulation {
    pub fn new(building: &Building, config: SimulationConfig) -> Result<Self> {
        let scene = FlatScene::new(building, config.voxel_size, config.search_transparent);
        Ok(Self { config, scene })
    }

    pub fn run(self) -> SimulationResult {
        match self.config.acoustic_mode {
            AcousticMode::Scalar => self.run_scalar(),
            AcousticMode::FrequencyDependent => self.run_frequency_dependent(),
        }
    }

    fn run_scalar(self) -> SimulationResult {
        let num_rays = self.config.num_rays;
        let num_steps = self.config.num_steps;
        let num_absorbers = self.config.absorbers.len();
        let dt = self.config.time_step;
        let speed = self.config.ray_speed;
        let absorber_r2 = self.config.absorber_radius * self.config.absorber_radius;

        // Build scalar absorption
        let coefficients: Vec<f64> = self
            .scene
            .paths
            .iter()
            .map(|path| {
                self.config
                    .absorption
                    .get(path)
                    .copied()
                    .unwrap_or(self.config.default_absorption)
            })
            .collect();
        let absorption = ScalarAbsorption::new(coefficients);

        let reflection_model = Specular;
        let propagation_model = FixedTimeStep;
        let reflection_dist = propagation_model.reflection_distance(speed, dt);

        let batch = RayBatch::new(self.config.source, speed, num_rays);
        let mut positions: Vec<Point> = batch.rays.iter().map(|r| r.position).collect();
        let mut velocities: Vec<Vector> = batch.rays.iter().map(|r| r.velocity).collect();
        let mut energies: Vec<f64> = batch.rays.iter().map(|r| r.energy).collect();

        let store_history = self.config.store_ray_history;
        let mut all_positions: Vec<Vec<Point>> = Vec::with_capacity(if store_history {
            num_steps
        } else {
            0
        });
        let mut all_energies: Vec<Vec<f64>> = Vec::with_capacity(if store_history {
            num_steps
        } else {
            0
        });
        let mut all_hits: Vec<Vec<f64>> = Vec::with_capacity(num_steps);

        let eps = 1e-10;
        let bbox_margin = 1e-6;
        let min_alive_fraction = self.config.min_alive_fraction;

        // Track whether each ray was inside each absorber on the previous step
        let mut was_inside: Vec<Vec<bool>> = vec![vec![false; num_absorbers]; num_rays];

        for _step in 0..num_steps {
            let mut step_hits = vec![0.0; num_absorbers];

            for ray in 0..num_rays {
                if energies[ray] <= eps {
                    continue;
                }
                for (ai, absorber) in self.config.absorbers.iter().enumerate() {
                    let dx = positions[ray].x - absorber.x;
                    let dy = positions[ray].y - absorber.y;
                    let dz = positions[ray].z - absorber.z;
                    let dist2 = dx * dx + dy * dy + dz * dz;
                    if dist2 <= absorber_r2 {
                        if !was_inside[ray][ai] {
                            // First entry: record the hit
                            step_hits[ai] += energies[ray];
                            was_inside[ray][ai] = true;
                        }
                        // Ray continues (NOT killed)
                    } else {
                        was_inside[ray][ai] = false;
                    }
                }
            }

            positions
                .par_iter_mut()
                .zip(velocities.par_iter_mut())
                .zip(energies.par_iter_mut())
                .for_each(|((pos, vel), energy)| {
                    if *energy <= eps {
                        return;
                    }

                    let orig_pos = *pos;
                    if !self.scene.is_in_bounds(orig_pos, bbox_margin) {
                        *energy = 0.0;
                        return;
                    }

                    let vel_norm = match vel.normalize() {
                        Ok(v) => v,
                        Err(_) => {
                            *energy = 0.0;
                            return;
                        }
                    };

                    if let Some((target_idx, target_dist)) =
                        self.scene.find_target_surface(orig_pos, vel_norm)
                    {
                        if target_dist <= reflection_dist {
                            let mut current_target = Some((target_idx, target_dist));
                            let mut search_pos = orig_pos;

                            while let Some((tidx, _tdist)) = current_target {
                                if *energy <= eps {
                                    break;
                                }

                                *energy = absorption.apply(*energy, tidx);
                                if *energy <= eps {
                                    *energy = 0.0;
                                    break;
                                }

                                let normal = self.scene.polygons[tidx].vn;
                                *vel = reflection_model.reflect(*vel, normal);

                                let new_vel_norm = match vel.normalize() {
                                    Ok(v) => v,
                                    Err(_) => break,
                                };

                                current_target =
                                    self.scene.find_target_surface(search_pos, new_vel_norm);

                                if let Some((_, next_dist)) = current_target {
                                    if next_dist > reflection_dist {
                                        current_target = None;
                                    } else {
                                        search_pos =
                                            search_pos + new_vel_norm * (next_dist * 0.5);
                                    }
                                }
                            }

                            if *energy > eps {
                                *pos = propagation_model.advance(orig_pos, *vel, dt);
                            }
                        } else {
                            *pos = propagation_model.advance(orig_pos, *vel, dt);
                        }
                    } else {
                        *pos = propagation_model.advance(orig_pos, *vel, dt);
                    }
                });

            if store_history {
                all_positions.push(positions.clone());
                all_energies.push(energies.clone());
            }
            all_hits.push(step_hits);

            // Early termination if too few rays are alive
            if min_alive_fraction > 0.0 {
                let alive = energies.iter().filter(|&&e| e > eps).count();
                if (alive as f64) / (num_rays as f64) < min_alive_fraction {
                    break;
                }
            }
        }

        SimulationResult {
            positions: all_positions,
            energies: all_energies,
            band_energies: None,
            hits: all_hits,
            band_hits: None,
            config: self.config,
        }
    }

    fn run_frequency_dependent(self) -> SimulationResult {
        let num_rays = self.config.num_rays;
        let num_steps = self.config.num_steps;
        let num_absorbers = self.config.absorbers.len();
        let dt = self.config.time_step;
        let speed = self.config.ray_speed;
        let absorber_r2 = self.config.absorber_radius * self.config.absorber_radius;

        // Build frequency-dependent absorption
        let band_coeffs = self.config.resolve_band_absorption(&self.scene.paths);
        let absorption = FrequencyDependentAbsorption::new(band_coeffs);

        // Build scattering coefficients for hybrid reflection
        let scattering_coeffs = self.config.resolve_scattering(&self.scene.paths);

        // Air absorption (optional)
        let air_absorption = if self.config.enable_air_absorption {
            Some(AirAbsorption::standard())
        } else {
            None
        };

        let propagation_model = FixedTimeStep;
        let reflection_dist = propagation_model.reflection_distance(speed, dt);

        let batch = RayBatch::new(self.config.source, speed, num_rays);
        let mut positions: Vec<Point> = batch.rays.iter().map(|r| r.position).collect();
        let mut velocities: Vec<Vector> = batch.rays.iter().map(|r| r.velocity).collect();
        let mut band_energies: Vec<[f64; NUM_OCTAVE_BANDS]> =
            vec![[1.0; NUM_OCTAVE_BANDS]; num_rays];

        let store_history = self.config.store_ray_history;
        let hist_cap = if store_history { num_steps } else { 0 };
        let mut all_positions: Vec<Vec<Point>> = Vec::with_capacity(hist_cap);
        let mut all_energies: Vec<Vec<f64>> = Vec::with_capacity(hist_cap);
        let mut all_band_energies: Vec<Vec<[f64; NUM_OCTAVE_BANDS]>> =
            Vec::with_capacity(hist_cap);
        let mut all_hits: Vec<Vec<f64>> = Vec::with_capacity(num_steps);
        let mut all_band_hits: Vec<Vec<[f64; NUM_OCTAVE_BANDS]>> = Vec::with_capacity(num_steps);

        let eps = 1e-10;
        let bbox_margin = 1e-6;
        let min_alive_fraction = self.config.min_alive_fraction;

        // Track whether each ray was inside each absorber on the previous step
        let mut was_inside: Vec<Vec<bool>> = vec![vec![false; num_absorbers]; num_rays];

        for _step in 0..num_steps {
            let mut step_hits = vec![0.0; num_absorbers];
            let mut step_band_hits = vec![[0.0; NUM_OCTAVE_BANDS]; num_absorbers];

            // Check absorbers
            for ray in 0..num_rays {
                let total_energy: f64 = band_energies[ray].iter().sum();
                if total_energy <= eps {
                    continue;
                }
                for (ai, absorber) in self.config.absorbers.iter().enumerate() {
                    let dx = positions[ray].x - absorber.x;
                    let dy = positions[ray].y - absorber.y;
                    let dz = positions[ray].z - absorber.z;
                    let dist2 = dx * dx + dy * dy + dz * dz;
                    if dist2 <= absorber_r2 {
                        if !was_inside[ray][ai] {
                            // First entry: record the hit
                            step_hits[ai] += total_energy;
                            for b in 0..NUM_OCTAVE_BANDS {
                                step_band_hits[ai][b] += band_energies[ray][b];
                            }
                            was_inside[ray][ai] = true;
                        }
                        // Ray continues (NOT killed)
                    } else {
                        was_inside[ray][ai] = false;
                    }
                }
            }

            // Propagate rays (parallel)
            positions
                .par_iter_mut()
                .zip(velocities.par_iter_mut())
                .zip(band_energies.par_iter_mut())
                .for_each(|((pos, vel), be)| {
                    let total_energy: f64 = be.iter().sum();
                    if total_energy <= eps {
                        return;
                    }

                    let orig_pos = *pos;
                    if !self.scene.is_in_bounds(orig_pos, bbox_margin) {
                        *be = [0.0; NUM_OCTAVE_BANDS];
                        return;
                    }

                    let vel_norm = match vel.normalize() {
                        Ok(v) => v,
                        Err(_) => {
                            *be = [0.0; NUM_OCTAVE_BANDS];
                            return;
                        }
                    };

                    if let Some((target_idx, target_dist)) =
                        self.scene.find_target_surface(orig_pos, vel_norm)
                    {
                        if target_dist <= reflection_dist {
                            let mut current_target = Some((target_idx, target_dist));
                            let mut search_pos = orig_pos;

                            while let Some((tidx, _tdist)) = current_target {
                                let total_e: f64 = be.iter().sum();
                                if total_e <= eps {
                                    break;
                                }

                                // Apply frequency-dependent absorption
                                *be = absorption.apply_bands(be, tidx);

                                let total_e: f64 = be.iter().sum();
                                if total_e <= eps {
                                    *be = [0.0; NUM_OCTAVE_BANDS];
                                    break;
                                }

                                // Use hybrid reflection with mean scattering
                                let mean_scattering: f64 =
                                    scattering_coeffs[tidx].iter().sum::<f64>()
                                        / NUM_OCTAVE_BANDS as f64;
                                let reflection_model = Hybrid::new(mean_scattering);

                                let normal = self.scene.polygons[tidx].vn;
                                *vel = reflection_model.reflect(*vel, normal);

                                let new_vel_norm = match vel.normalize() {
                                    Ok(v) => v,
                                    Err(_) => break,
                                };

                                current_target =
                                    self.scene.find_target_surface(search_pos, new_vel_norm);

                                if let Some((_, next_dist)) = current_target {
                                    if next_dist > reflection_dist {
                                        current_target = None;
                                    } else {
                                        search_pos =
                                            search_pos + new_vel_norm * (next_dist * 0.5);
                                    }
                                }
                            }

                            // Apply air absorption for this step
                            if let Some(ref air) = air_absorption {
                                let step_distance = speed * dt;
                                let factors = air.apply_distance(step_distance);
                                for (b, f) in be.iter_mut().zip(factors.iter()) {
                                    *b *= f;
                                }
                            }

                            let total_e: f64 = be.iter().sum();
                            if total_e > eps {
                                *pos = propagation_model.advance(orig_pos, *vel, dt);
                            }
                        } else {
                            // Apply air absorption for free-flight step
                            if let Some(ref air) = air_absorption {
                                let step_distance = speed * dt;
                                let factors = air.apply_distance(step_distance);
                                for (b, f) in be.iter_mut().zip(factors.iter()) {
                                    *b *= f;
                                }
                            }
                            *pos = propagation_model.advance(orig_pos, *vel, dt);
                        }
                    } else {
                        if let Some(ref air) = air_absorption {
                            let step_distance = speed * dt;
                            let factors = air.apply_distance(step_distance);
                            for (b, f) in be.iter_mut().zip(factors.iter()) {
                                *b *= f;
                            }
                        }
                        *pos = propagation_model.advance(orig_pos, *vel, dt);
                    }
                });

            // Record state
            if store_history {
                let scalar_energies: Vec<f64> =
                    band_energies.iter().map(|be| be.iter().sum()).collect();
                all_positions.push(positions.clone());
                all_energies.push(scalar_energies);
                all_band_energies.push(band_energies.clone());
            }
            all_hits.push(step_hits);
            all_band_hits.push(step_band_hits);

            // Early termination if too few rays are alive
            if min_alive_fraction > 0.0 {
                let alive = band_energies
                    .iter()
                    .filter(|be| be.iter().sum::<f64>() > eps)
                    .count();
                if (alive as f64) / (num_rays as f64) < min_alive_fraction {
                    break;
                }
            }
        }

        SimulationResult {
            positions: all_positions,
            energies: all_energies,
            band_energies: Some(all_band_energies),
            hits: all_hits,
            band_hits: Some(all_band_hits),
            config: self.config,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::materials::{AcousticMaterial, Material, MaterialLibrary};
    use crate::{Solid, Zone};

    #[test]
    fn test_simulation_basic() {
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "s0").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config = SimulationConfig::new();
        config.num_steps = 10;
        config.num_rays = 5;
        config.source = Point::new(1.0, 1.0, 1.0);
        config.store_ray_history = true;

        let sim = Simulation::new(&building, config).unwrap();
        let result = sim.run();

        assert_eq!(result.positions.len(), 10);
        assert_eq!(result.energies.len(), 10);
        assert_eq!(result.positions[0].len(), 5);
        assert!(result.band_energies.is_none());
    }

    #[test]
    fn test_simulation_with_absorber() {
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "s0").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config = SimulationConfig::new();
        config.num_steps = 200;
        config.num_rays = 200;
        config.source = Point::new(1.0, 1.0, 1.0);
        config.absorbers = vec![Point::new(0.5, 0.5, 0.5)];
        config.absorber_radius = 0.3;

        let sim = Simulation::new(&building, config).unwrap();
        let result = sim.run();

        let total_hit: f64 = result.hits.iter().map(|h| h[0]).sum();
        assert!(total_hit > 0.0, "Some rays should hit the absorber");
    }

    #[test]
    fn test_simulation_transparent() {
        // Two adjacent boxes in the same zone — rays should pass through the shared face
        let s0 = Solid::from_box(1.0, 1.0, 1.0, None, "s0").unwrap();
        let s1 = Solid::from_box(1.0, 1.0, 1.0, Some((1.0, 0.0, 0.0)), "s1").unwrap();
        let zone = Zone::new("z", vec![s0, s1]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config = SimulationConfig::new();
        config.num_steps = 200;
        config.num_rays = 200;
        config.source = Point::new(0.5, 0.5, 0.5);
        config.absorbers = vec![Point::new(1.5, 0.5, 0.5)];
        config.absorber_radius = 0.4;

        let sim = Simulation::new(&building, config).unwrap();
        let result = sim.run();

        let total_hit: f64 = result.hits.iter().map(|h| h[0]).sum();
        assert!(
            total_hit > 0.0,
            "Rays should pass through transparent face and hit absorber"
        );
    }

    #[test]
    fn test_frequency_dependent_simulation() {
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "s0").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut lib = MaterialLibrary::new();
        lib.add(
            Material::new("absorber").with_acoustic(AcousticMaterial::new(
                "absorber",
                [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
                [0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
            )),
        );
        // Assign to all surfaces
        lib.assign("/", "absorber");

        let mut config = SimulationConfig::new();
        config.num_steps = 50;
        config.num_rays = 20;
        config.source = Point::new(1.0, 1.0, 1.0);
        config.acoustic_mode = AcousticMode::FrequencyDependent;
        config.material_library = Some(lib);
        config.store_ray_history = true;

        let sim = Simulation::new(&building, config).unwrap();
        let result = sim.run();

        assert!(result.band_energies.is_some());
        let band_e = result.band_energies.as_ref().unwrap();
        assert_eq!(band_e.len(), 50);
        assert_eq!(band_e[0].len(), 20);

        // High-frequency bands should lose more energy than low-frequency
        let last_step = &band_e[band_e.len() - 1];
        let avg_low: f64 = last_step.iter().map(|be| be[0]).sum::<f64>() / 20.0;
        let avg_high: f64 = last_step.iter().map(|be| be[5]).sum::<f64>() / 20.0;
        assert!(
            avg_low >= avg_high,
            "Low frequency should retain more energy than high frequency"
        );
    }

    #[test]
    fn test_frequency_dependent_with_air_absorption() {
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "s0").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config = SimulationConfig::new();
        config.num_steps = 50;
        config.num_rays = 10;
        config.source = Point::new(1.0, 1.0, 1.0);
        config.acoustic_mode = AcousticMode::FrequencyDependent;
        config.enable_air_absorption = true;

        let sim = Simulation::new(&building, config).unwrap();
        let result = sim.run();

        assert!(result.band_energies.is_some());
        assert!(result.band_hits.is_some());
    }
}
