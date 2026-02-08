use anyhow::Result;

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
        let mut all_positions: Vec<Vec<Point>> =
            Vec::with_capacity(if store_history { num_steps } else { 0 });
        let mut all_energies: Vec<Vec<f64>> =
            Vec::with_capacity(if store_history { num_steps } else { 0 });
        let mut all_hits: Vec<Vec<f64>> = Vec::with_capacity(num_steps);

        let eps = 1e-10;
        let bbox_margin = 1e-6;

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
                        step_hits[ai] += energies[ray];
                        energies[ray] = 0.0;
                        break;
                    }
                }
            }

            for ray in 0..num_rays {
                if energies[ray] <= eps {
                    continue;
                }

                let pos = positions[ray];
                if !self.scene.is_in_bounds(pos, bbox_margin) {
                    energies[ray] = 0.0;
                    continue;
                }

                let vel_norm = match velocities[ray].normalize() {
                    Ok(v) => v,
                    Err(_) => {
                        energies[ray] = 0.0;
                        continue;
                    }
                };

                if let Some((target_idx, target_dist)) =
                    self.scene.find_target_surface(pos, vel_norm)
                {
                    if target_dist <= reflection_dist {
                        let mut current_target = Some((target_idx, target_dist));
                        let mut search_pos = pos;

                        while let Some((tidx, _tdist)) = current_target {
                            if energies[ray] <= eps {
                                break;
                            }

                            energies[ray] = absorption.apply(energies[ray], tidx);
                            if energies[ray] <= eps {
                                energies[ray] = 0.0;
                                break;
                            }

                            let normal = self.scene.polygons[tidx].vn;
                            velocities[ray] = reflection_model.reflect(velocities[ray], normal);

                            let new_vel_norm = match velocities[ray].normalize() {
                                Ok(v) => v,
                                Err(_) => break,
                            };

                            current_target =
                                self.scene.find_target_surface(search_pos, new_vel_norm);

                            if let Some((_, next_dist)) = current_target {
                                if next_dist > reflection_dist {
                                    current_target = None;
                                } else {
                                    search_pos = search_pos + new_vel_norm * (next_dist * 0.5);
                                }
                            }
                        }

                        if energies[ray] > eps {
                            positions[ray] = propagation_model.advance(pos, velocities[ray], dt);
                        }
                    } else {
                        positions[ray] = propagation_model.advance(pos, velocities[ray], dt);
                    }
                } else {
                    positions[ray] = propagation_model.advance(pos, velocities[ray], dt);
                }
            }

            if store_history {
                all_positions.push(positions.clone());
                all_energies.push(energies.clone());
            }
            all_hits.push(step_hits);
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
        let mut all_band_energies: Vec<Vec<[f64; NUM_OCTAVE_BANDS]>> = Vec::with_capacity(hist_cap);
        let mut all_hits: Vec<Vec<f64>> = Vec::with_capacity(num_steps);
        let mut all_band_hits: Vec<Vec<[f64; NUM_OCTAVE_BANDS]>> = Vec::with_capacity(num_steps);

        let eps = 1e-10;
        let bbox_margin = 1e-6;

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
                        step_hits[ai] += total_energy;
                        for b in 0..NUM_OCTAVE_BANDS {
                            step_band_hits[ai][b] += band_energies[ray][b];
                        }
                        band_energies[ray] = [0.0; NUM_OCTAVE_BANDS];
                        break;
                    }
                }
            }

            // Propagate rays
            for ray in 0..num_rays {
                let total_energy: f64 = band_energies[ray].iter().sum();
                if total_energy <= eps {
                    continue;
                }

                let pos = positions[ray];
                if !self.scene.is_in_bounds(pos, bbox_margin) {
                    band_energies[ray] = [0.0; NUM_OCTAVE_BANDS];
                    continue;
                }

                let vel_norm = match velocities[ray].normalize() {
                    Ok(v) => v,
                    Err(_) => {
                        band_energies[ray] = [0.0; NUM_OCTAVE_BANDS];
                        continue;
                    }
                };

                if let Some((target_idx, target_dist)) =
                    self.scene.find_target_surface(pos, vel_norm)
                {
                    if target_dist <= reflection_dist {
                        let mut current_target = Some((target_idx, target_dist));
                        let mut search_pos = pos;

                        while let Some((tidx, _tdist)) = current_target {
                            let total_e: f64 = band_energies[ray].iter().sum();
                            if total_e <= eps {
                                break;
                            }

                            // Apply frequency-dependent absorption
                            band_energies[ray] = absorption.apply_bands(&band_energies[ray], tidx);

                            let total_e: f64 = band_energies[ray].iter().sum();
                            if total_e <= eps {
                                band_energies[ray] = [0.0; NUM_OCTAVE_BANDS];
                                break;
                            }

                            // Use hybrid reflection with mean scattering
                            let mean_scattering: f64 = scattering_coeffs[tidx].iter().sum::<f64>()
                                / NUM_OCTAVE_BANDS as f64;
                            let reflection_model = Hybrid::new(mean_scattering);

                            let normal = self.scene.polygons[tidx].vn;
                            velocities[ray] = reflection_model.reflect(velocities[ray], normal);

                            let new_vel_norm = match velocities[ray].normalize() {
                                Ok(v) => v,
                                Err(_) => break,
                            };

                            current_target =
                                self.scene.find_target_surface(search_pos, new_vel_norm);

                            if let Some((_, next_dist)) = current_target {
                                if next_dist > reflection_dist {
                                    current_target = None;
                                } else {
                                    search_pos = search_pos + new_vel_norm * (next_dist * 0.5);
                                }
                            }
                        }

                        // Apply air absorption for this step
                        if let Some(ref air) = air_absorption {
                            let step_distance = speed * dt;
                            let factors = air.apply_distance(step_distance);
                            for (be, f) in band_energies[ray].iter_mut().zip(factors.iter()) {
                                *be *= f;
                            }
                        }

                        let total_e: f64 = band_energies[ray].iter().sum();
                        if total_e > eps {
                            positions[ray] = propagation_model.advance(pos, velocities[ray], dt);
                        }
                    } else {
                        // Apply air absorption for free-flight step
                        if let Some(ref air) = air_absorption {
                            let step_distance = speed * dt;
                            let factors = air.apply_distance(step_distance);
                            for (be, f) in band_energies[ray].iter_mut().zip(factors.iter()) {
                                *be *= f;
                            }
                        }
                        positions[ray] = propagation_model.advance(pos, velocities[ray], dt);
                    }
                } else {
                    if let Some(ref air) = air_absorption {
                        let step_distance = speed * dt;
                        let factors = air.apply_distance(step_distance);
                        for (be, f) in band_energies[ray].iter_mut().zip(factors.iter()) {
                            *be *= f;
                        }
                    }
                    positions[ray] = propagation_model.advance(pos, velocities[ray], dt);
                }
            }

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

    // ── Physics verification tests ──────────────────────────────────────

    #[test]
    fn test_energy_decay_monotonic() {
        // Total ray energy should decrease (or stay equal) over time
        // as energy is absorbed at each reflection.
        //
        // Use a small room (2m) with high absorption (0.5) and many time
        // steps (2000) so rays undergo many reflections with significant loss.
        // Mean free path = 4*V/S = 4*8/24 = 1.33m, distance per step = 343*2.5e-5 = 0.00858m.
        // Over 2000 steps distance = 17.2m → ~13 reflections → (0.5)^13 ≈ 0.0001.
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config = SimulationConfig::new();
        config.num_steps = 2000;
        config.num_rays = 100;
        config.source = Point::new(1.0, 1.0, 1.0);
        config.default_absorption = 0.5;
        config.store_ray_history = true;

        let sim = Simulation::new(&building, config).unwrap();
        let result = sim.run();

        // Compute total energy at each step
        let total_per_step: Vec<f64> = result
            .energies
            .iter()
            .map(|step| step.iter().sum::<f64>())
            .collect();

        // Energy should be non-increasing (allowing for small numerical noise)
        for i in 1..total_per_step.len() {
            assert!(
                total_per_step[i] <= total_per_step[i - 1] + 1e-10,
                "Energy should not increase: step {i} has {:.6} > step {} had {:.6}",
                total_per_step[i],
                i - 1,
                total_per_step[i - 1]
            );
        }

        // Energy should actually decrease significantly over many reflections
        let first = total_per_step[0];
        let last = *total_per_step.last().unwrap();
        assert!(
            last < first * 0.1,
            "Energy should decay to <10% of initial: first={first:.3}, last={last:.3}"
        );
    }

    #[test]
    fn test_higher_absorption_faster_decay() {
        // Higher absorption coefficients should cause energy to decay faster.
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let run_with_absorption = |alpha: f64| -> f64 {
            let mut config = SimulationConfig::new();
            config.num_steps = 2000;
            config.num_rays = 200;
            config.source = Point::new(1.0, 1.0, 1.0);
            config.default_absorption = alpha;
            config.store_ray_history = true;

            let sim = Simulation::new(&building, config).unwrap();
            let result = sim.run();

            // Return total energy at last step
            result.energies.last().unwrap().iter().sum()
        };

        let energy_low_abs = run_with_absorption(0.1);
        let energy_high_abs = run_with_absorption(0.5);

        assert!(
            energy_high_abs < energy_low_abs,
            "Higher absorption should result in less remaining energy. \
             α=0.1 → {energy_low_abs:.4}, α=0.5 → {energy_high_abs:.4}"
        );
    }

    #[test]
    fn test_larger_room_slower_decay() {
        // In a larger room, rays travel longer between reflections so energy
        // decays more slowly per time step (Sabine: RT ∝ V/A).
        let run_for_room_size = |size: f64| -> f64 {
            let s0 = Solid::from_box(size, size, size, None, "room").unwrap();
            let zone = Zone::new("z", vec![s0]).unwrap();
            let building = Building::new("b", vec![zone]).unwrap();

            let mut config = SimulationConfig::new();
            config.num_steps = 2000;
            config.num_rays = 200;
            config.source = Point::new(size / 2.0, size / 2.0, size / 2.0);
            config.default_absorption = 0.3;
            config.store_ray_history = true;

            let sim = Simulation::new(&building, config).unwrap();
            let result = sim.run();

            result.energies.last().unwrap().iter().sum()
        };

        let energy_small = run_for_room_size(2.0); // 2×2×2 m
        let energy_large = run_for_room_size(6.0); // 6×6×6 m

        assert!(
            energy_large > energy_small,
            "Larger room should have more remaining energy (slower decay). \
             small={energy_small:.4}, large={energy_large:.4}"
        );
    }

    #[test]
    fn test_absorber_collects_energy_proportional_to_solid_angle() {
        // A larger absorber should collect more energy.
        let s0 = Solid::from_box(4.0, 4.0, 4.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let run_with_radius = |radius: f64| -> f64 {
            let mut config = SimulationConfig::new();
            config.num_steps = 500;
            config.num_rays = 1000;
            config.source = Point::new(2.0, 2.0, 2.0);
            config.absorbers = vec![Point::new(1.0, 1.0, 1.0)];
            config.absorber_radius = radius;

            let sim = Simulation::new(&building, config).unwrap();
            let result = sim.run();
            result.hits.iter().map(|h| h[0]).sum()
        };

        let energy_small = run_with_radius(0.2);
        let energy_large = run_with_radius(0.5);

        assert!(
            energy_large > energy_small,
            "Larger absorber should collect more energy. \
             r=0.2 → {energy_small:.4}, r=0.5 → {energy_large:.4}"
        );
    }

    #[test]
    fn test_frequency_dependent_decay_ordering() {
        // With frequency-dependent absorption [0.05, 0.10, 0.20, 0.35, 0.50, 0.70],
        // higher bands should have less remaining energy after many steps.
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut lib = MaterialLibrary::new();
        lib.add(
            Material::new("walls").with_acoustic(AcousticMaterial::new(
                "walls",
                [0.05, 0.10, 0.20, 0.35, 0.50, 0.70],
                [0.1; 6],
            )),
        );
        lib.assign("/", "walls");

        let mut config = SimulationConfig::new();
        config.num_steps = 2000;
        config.num_rays = 200;
        config.source = Point::new(1.0, 1.0, 1.0);
        config.acoustic_mode = AcousticMode::FrequencyDependent;
        config.material_library = Some(lib);
        config.store_ray_history = true;

        let sim = Simulation::new(&building, config).unwrap();
        let result = sim.run();

        let band_e = result.band_energies.as_ref().unwrap();
        let last_step = band_e.last().unwrap();

        // Sum remaining energy per band across all rays
        let mut band_totals = [0.0_f64; 6];
        for ray_bands in last_step {
            for (b, &e) in ray_bands.iter().enumerate() {
                band_totals[b] += e;
            }
        }

        // Each band should have less energy than the previous (lower absorption)
        for b in 1..6 {
            assert!(
                band_totals[b] <= band_totals[b - 1] + 1e-6,
                "Band {b} ({:.4}) should have ≤ energy than band {} ({:.4})",
                band_totals[b],
                b - 1,
                band_totals[b - 1]
            );
        }
    }
}
