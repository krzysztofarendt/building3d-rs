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

/// Per-ray energy below this value is treated as "dead".
pub const ENERGY_EPS: f64 = 1e-10;

#[derive(Debug, Clone, Copy)]
pub struct SimulationProgress {
    /// Number of completed steps (0..=num_steps).
    pub steps_done: usize,
    /// Target number of steps from the configuration (may stop early).
    pub num_steps: usize,
    /// Simulation time step (seconds).
    pub dt_s: f64,
    /// Simulated time elapsed (seconds).
    pub sim_time_s: f64,
    /// Total number of rays.
    pub num_rays: usize,
    /// Rays with energy > `ENERGY_EPS`.
    pub alive_rays: usize,
    /// Total scalar energy remaining (frequency-dependent mode sums bands).
    pub total_energy: f64,
    /// Total scalar energy at start (for fractions).
    pub initial_total_energy: f64,
}

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

trait ProgressReporter {
    fn every_steps(&self) -> usize;
    fn report(&mut self, progress: &SimulationProgress);
}

struct NoProgress;
impl ProgressReporter for NoProgress {
    fn every_steps(&self) -> usize {
        0
    }
    fn report(&mut self, _progress: &SimulationProgress) {}
}

struct FnProgress<F> {
    every_steps: usize,
    f: F,
}
impl<F> ProgressReporter for FnProgress<F>
where
    F: FnMut(&SimulationProgress),
{
    fn every_steps(&self) -> usize {
        self.every_steps
    }
    fn report(&mut self, progress: &SimulationProgress) {
        (self.f)(progress);
    }
}

impl Simulation {
    pub fn new(building: &Building, config: SimulationConfig) -> Result<Self> {
        let scene = FlatScene::new(building, config.voxel_size, config.search_transparent);
        Ok(Self { config, scene })
    }

    pub fn run(self) -> SimulationResult {
        match self.config.acoustic_mode {
            AcousticMode::Scalar => self.run_scalar_with_progress(NoProgress),
            AcousticMode::FrequencyDependent => {
                self.run_frequency_dependent_with_progress(NoProgress)
            }
        }
    }

    /// Runs the simulation while periodically reporting progress.
    ///
    /// - `every_steps=0` disables progress reporting.
    /// - The reporter is called once at start (`steps_done=0`) and then every `every_steps`,
    ///   plus once at the end (or early-termination step).
    pub fn run_with_progress<F>(self, every_steps: usize, report: F) -> SimulationResult
    where
        F: FnMut(&SimulationProgress),
    {
        let reporter = FnProgress {
            every_steps,
            f: report,
        };
        match self.config.acoustic_mode {
            AcousticMode::Scalar => self.run_scalar_with_progress(reporter),
            AcousticMode::FrequencyDependent => {
                self.run_frequency_dependent_with_progress(reporter)
            }
        }
    }

    fn run_scalar_with_progress<R: ProgressReporter>(self, mut reporter: R) -> SimulationResult {
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

        let eps = ENERGY_EPS;
        let bbox_margin = 1e-6;
        let min_alive_fraction = self.config.min_alive_fraction;

        // Track whether each ray was inside each absorber on the previous step
        let mut was_inside: Vec<Vec<bool>> = vec![vec![false; num_absorbers]; num_rays];

        let initial_total_energy: f64 = num_rays as f64; // energies start at 1.0
        let report_every = reporter.every_steps();
        if report_every > 0 {
            reporter.report(&SimulationProgress {
                steps_done: 0,
                num_steps,
                dt_s: dt,
                sim_time_s: 0.0,
                num_rays,
                alive_rays: num_rays,
                total_energy: initial_total_energy,
                initial_total_energy,
            });
        }

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
                                        search_pos = search_pos + new_vel_norm * (next_dist * 0.5);
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

            let steps_done = all_hits.len();
            let should_report = report_every > 0
                && (steps_done.is_multiple_of(report_every) || steps_done == num_steps);
            let needs_stats = should_report || min_alive_fraction > 0.0;

            if needs_stats {
                let mut alive: usize = 0;
                let mut total_energy: f64 = 0.0;
                for &e in energies.iter() {
                    if e > eps {
                        alive += 1;
                    }
                    total_energy += e;
                }

                let alive_fraction = (alive as f64) / (num_rays as f64);
                let would_terminate =
                    min_alive_fraction > 0.0 && alive_fraction < min_alive_fraction;

                if (should_report || would_terminate) && report_every > 0 {
                    reporter.report(&SimulationProgress {
                        steps_done,
                        num_steps,
                        dt_s: dt,
                        sim_time_s: steps_done as f64 * dt,
                        num_rays,
                        alive_rays: alive,
                        total_energy,
                        initial_total_energy,
                    });
                }

                if would_terminate {
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

    fn run_frequency_dependent_with_progress<R: ProgressReporter>(
        self,
        mut reporter: R,
    ) -> SimulationResult {
        if let Some(band) = self.config.single_band_index {
            return self.run_frequency_single_band_with_progress(reporter, band);
        }

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
        let store_band_history = store_history && self.config.store_ray_band_history;
        let hist_cap = if store_history { num_steps } else { 0 };
        let mut all_positions: Vec<Vec<Point>> = Vec::with_capacity(hist_cap);
        let mut all_energies: Vec<Vec<f64>> = Vec::with_capacity(hist_cap);
        let mut all_band_energies: Option<Vec<Vec<[f64; NUM_OCTAVE_BANDS]>>> =
            store_band_history.then(|| Vec::with_capacity(hist_cap));
        let mut all_hits: Vec<Vec<f64>> = Vec::with_capacity(num_steps);
        let mut all_band_hits: Vec<Vec<[f64; NUM_OCTAVE_BANDS]>> = Vec::with_capacity(num_steps);

        let eps = ENERGY_EPS;
        let bbox_margin = 1e-6;
        let min_alive_fraction = self.config.min_alive_fraction;

        // Track whether each ray was inside each absorber on the previous step
        let mut was_inside: Vec<Vec<bool>> = vec![vec![false; num_absorbers]; num_rays];

        let initial_total_energy: f64 = (num_rays as f64) * (NUM_OCTAVE_BANDS as f64);
        let report_every = reporter.every_steps();
        if report_every > 0 {
            reporter.report(&SimulationProgress {
                steps_done: 0,
                num_steps,
                dt_s: dt,
                sim_time_s: 0.0,
                num_rays,
                alive_rays: num_rays,
                total_energy: initial_total_energy,
                initial_total_energy,
            });
        }

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
                                        search_pos = search_pos + new_vel_norm * (next_dist * 0.5);
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
                if let Some(ref mut v) = all_band_energies {
                    v.push(band_energies.clone());
                }
            }
            all_hits.push(step_hits);
            all_band_hits.push(step_band_hits);

            let steps_done = all_hits.len();
            let should_report = report_every > 0
                && (steps_done.is_multiple_of(report_every) || steps_done == num_steps);
            let needs_stats = should_report || min_alive_fraction > 0.0;

            if needs_stats {
                let mut alive: usize = 0;
                let mut total_energy: f64 = 0.0;
                for be in band_energies.iter() {
                    let e = be.iter().sum::<f64>();
                    if e > eps {
                        alive += 1;
                    }
                    total_energy += e;
                }

                let alive_fraction = (alive as f64) / (num_rays as f64);
                let would_terminate =
                    min_alive_fraction > 0.0 && alive_fraction < min_alive_fraction;

                if (should_report || would_terminate) && report_every > 0 {
                    reporter.report(&SimulationProgress {
                        steps_done,
                        num_steps,
                        dt_s: dt,
                        sim_time_s: steps_done as f64 * dt,
                        num_rays,
                        alive_rays: alive,
                        total_energy,
                        initial_total_energy,
                    });
                }

                if would_terminate {
                    break;
                }
            }
        }

        SimulationResult {
            positions: all_positions,
            energies: all_energies,
            band_energies: all_band_energies,
            hits: all_hits,
            band_hits: Some(all_band_hits),
            config: self.config,
        }
    }

    fn run_frequency_single_band_with_progress<R: ProgressReporter>(
        self,
        mut reporter: R,
        band: usize,
    ) -> SimulationResult {
        let num_rays = self.config.num_rays;
        let num_steps = self.config.num_steps;
        let num_absorbers = self.config.absorbers.len();
        let dt = self.config.time_step;
        let speed = self.config.ray_speed;
        let absorber_r2 = self.config.absorber_radius * self.config.absorber_radius;

        let band = band.min(NUM_OCTAVE_BANDS.saturating_sub(1));

        let band_absorption = self
            .config
            .resolve_band_absorption(&self.scene.paths)
            .into_iter()
            .map(|a| a[band])
            .collect::<Vec<f64>>();
        let band_scattering = self
            .config
            .resolve_scattering(&self.scene.paths)
            .into_iter()
            .map(|s| s[band])
            .collect::<Vec<f64>>();

        // Air absorption (optional)
        let air_factor = if self.config.enable_air_absorption {
            let air = AirAbsorption::standard();
            let step_distance = speed * dt;
            air.apply_distance(step_distance)[band]
        } else {
            1.0
        };

        let propagation_model = FixedTimeStep;
        let reflection_dist = propagation_model.reflection_distance(speed, dt);

        let batch = RayBatch::new(self.config.source, speed, num_rays);
        let mut positions: Vec<Point> = batch.rays.iter().map(|r| r.position).collect();
        let mut velocities: Vec<Vector> = batch.rays.iter().map(|r| r.velocity).collect();
        let mut energies: Vec<f64> = vec![1.0; num_rays];

        let store_history = self.config.store_ray_history;
        let store_band_history = store_history && self.config.store_ray_band_history;
        let hist_cap = if store_history { num_steps } else { 0 };
        let mut all_positions: Vec<Vec<Point>> = Vec::with_capacity(hist_cap);
        let mut all_energies: Vec<Vec<f64>> = Vec::with_capacity(hist_cap);
        let mut all_band_energies: Option<Vec<Vec<[f64; NUM_OCTAVE_BANDS]>>> =
            store_band_history.then(|| Vec::with_capacity(hist_cap));
        let mut all_hits: Vec<Vec<f64>> = Vec::with_capacity(num_steps);
        let mut all_band_hits: Vec<Vec<[f64; NUM_OCTAVE_BANDS]>> = Vec::with_capacity(num_steps);

        let eps = ENERGY_EPS;
        let bbox_margin = 1e-6;
        let min_alive_fraction = self.config.min_alive_fraction;

        // Track whether each ray was inside each absorber on the previous step
        let mut was_inside: Vec<Vec<bool>> = vec![vec![false; num_absorbers]; num_rays];

        let initial_total_energy: f64 = num_rays as f64;
        let report_every = reporter.every_steps();
        if report_every > 0 {
            reporter.report(&SimulationProgress {
                steps_done: 0,
                num_steps,
                dt_s: dt,
                sim_time_s: 0.0,
                num_rays,
                alive_rays: num_rays,
                total_energy: initial_total_energy,
                initial_total_energy,
            });
        }

        for _step in 0..num_steps {
            let mut step_hits = vec![0.0; num_absorbers];
            let mut step_band_hits = vec![[0.0; NUM_OCTAVE_BANDS]; num_absorbers];

            // Check absorbers
            for ray in 0..num_rays {
                let e = energies[ray];
                if e <= eps {
                    continue;
                }
                for (ai, absorber) in self.config.absorbers.iter().enumerate() {
                    let dx = positions[ray].x - absorber.x;
                    let dy = positions[ray].y - absorber.y;
                    let dz = positions[ray].z - absorber.z;
                    let dist2 = dx * dx + dy * dy + dz * dz;
                    if dist2 <= absorber_r2 {
                        if !was_inside[ray][ai] {
                            step_hits[ai] += e;
                            step_band_hits[ai][band] += e;
                            was_inside[ray][ai] = true;
                        }
                    } else {
                        was_inside[ray][ai] = false;
                    }
                }
            }

            // Propagate rays (parallel)
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

                                let alpha = band_absorption.get(tidx).copied().unwrap_or(0.0);
                                *energy *= 1.0 - alpha;
                                if *energy <= eps {
                                    *energy = 0.0;
                                    break;
                                }

                                let scattering = band_scattering.get(tidx).copied().unwrap_or(0.0);
                                let reflection_model = Hybrid::new(scattering);

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
                                        search_pos = search_pos + new_vel_norm * (next_dist * 0.5);
                                    }
                                }
                            }

                            if *energy > eps {
                                *energy *= air_factor;
                                *pos = propagation_model.advance(orig_pos, *vel, dt);
                            }
                        } else {
                            *energy *= air_factor;
                            *pos = propagation_model.advance(orig_pos, *vel, dt);
                        }
                    } else {
                        *energy *= air_factor;
                        *pos = propagation_model.advance(orig_pos, *vel, dt);
                    }
                });

            if store_history {
                all_positions.push(positions.clone());
                all_energies.push(energies.clone());
                if let Some(ref mut v) = all_band_energies {
                    let bands: Vec<[f64; NUM_OCTAVE_BANDS]> = energies
                        .iter()
                        .map(|&e| {
                            let mut arr = [0.0; NUM_OCTAVE_BANDS];
                            arr[band] = e;
                            arr
                        })
                        .collect();
                    v.push(bands);
                }
            }

            all_hits.push(step_hits);
            all_band_hits.push(step_band_hits);

            let steps_done = all_hits.len();
            let should_report = report_every > 0
                && (steps_done.is_multiple_of(report_every) || steps_done == num_steps);
            let needs_stats = should_report || min_alive_fraction > 0.0;

            if needs_stats {
                let mut alive: usize = 0;
                let mut total_energy: f64 = 0.0;
                for &e in energies.iter() {
                    if e > eps {
                        alive += 1;
                    }
                    total_energy += e;
                }

                let alive_fraction = (alive as f64) / (num_rays as f64);
                let would_terminate =
                    min_alive_fraction > 0.0 && alive_fraction < min_alive_fraction;

                if (should_report || would_terminate) && report_every > 0 {
                    reporter.report(&SimulationProgress {
                        steps_done,
                        num_steps,
                        dt_s: dt,
                        sim_time_s: steps_done as f64 * dt,
                        num_rays,
                        alive_rays: alive,
                        total_energy,
                        initial_total_energy,
                    });
                }

                if would_terminate {
                    break;
                }
            }
        }

        SimulationResult {
            positions: all_positions,
            energies: all_energies,
            band_energies: all_band_energies,
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
    fn test_progress_reporter_is_called() {
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "s0").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config = SimulationConfig::new();
        config.num_steps = 10;
        config.num_rays = 5;
        config.source = Point::new(1.0, 1.0, 1.0);
        config.default_absorption = 0.0;
        config.store_ray_history = false;

        let mut calls: usize = 0;
        let mut last_steps_done: usize = 999;
        let sim = Simulation::new(&building, config).unwrap();
        let _result = sim.run_with_progress(3, |p| {
            calls += 1;
            last_steps_done = p.steps_done;
        });

        // Called at start (0), then at steps 3/6/9, and at end (10).
        assert_eq!(calls, 5);
        assert_eq!(last_steps_done, 10);
    }

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
        config.store_ray_band_history = true;

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
    fn test_frequency_dependent_band_history_is_opt_in() {
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "s0").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config = SimulationConfig::new();
        config.num_steps = 5;
        config.num_rays = 3;
        config.source = Point::new(1.0, 1.0, 1.0);
        config.acoustic_mode = AcousticMode::FrequencyDependent;
        config.store_ray_history = true;
        // store_ray_band_history intentionally left as default (false).

        let sim = Simulation::new(&building, config).unwrap();
        let result = sim.run();

        assert!(result.band_hits.is_some());
        assert!(result.band_energies.is_none());
        assert_eq!(result.positions.len(), 5);
        assert_eq!(result.energies.len(), 5);
    }

    #[test]
    fn test_frequency_dependent_single_band_runs() {
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "s0").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config = SimulationConfig::new();
        config.num_steps = 10;
        config.num_rays = 20;
        config.source = Point::new(1.0, 1.0, 1.0);
        config.acoustic_mode = AcousticMode::FrequencyDependent;
        config.single_band_index = Some(3);
        config.store_ray_history = true;

        let sim = Simulation::new(&building, config).unwrap();
        let result = sim.run();

        assert_eq!(result.positions.len(), 10);
        assert_eq!(result.energies.len(), 10);
        assert!(result.band_hits.is_some());
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
        config.store_ray_history = true;
        config.store_ray_band_history = true;

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
    fn test_scalar_min_alive_fraction_early_termination() {
        // With very high absorption (0.99) and min_alive_fraction = 0.99,
        // simulation should terminate early because rays die quickly.
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config = SimulationConfig::new();
        config.num_steps = 10_000;
        config.num_rays = 50;
        config.source = Point::new(1.0, 1.0, 1.0);
        config.default_absorption = 0.99;
        config.min_alive_fraction = 0.99;
        config.store_ray_history = false;

        let sim = Simulation::new(&building, config).unwrap();
        let result = sim.run();

        // Should have terminated early (well before 10_000 steps)
        assert!(
            result.hits.len() < 10_000,
            "Should terminate early: ran {} steps",
            result.hits.len()
        );
    }

    #[test]
    fn test_frequency_dependent_with_progress_callback() {
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "s0").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut lib = MaterialLibrary::new();
        lib.add(Material::new("mat").with_acoustic(AcousticMaterial::new(
            "mat",
            [0.1, 0.2, 0.3, 0.4, 0.5, 0.6],
            [0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
        )));
        lib.assign("/", "mat");

        let mut config = SimulationConfig::new();
        config.num_steps = 10;
        config.num_rays = 5;
        config.source = Point::new(1.0, 1.0, 1.0);
        config.acoustic_mode = AcousticMode::FrequencyDependent;
        config.material_library = Some(lib);
        config.store_ray_history = false;

        let mut calls: usize = 0;
        let sim = Simulation::new(&building, config).unwrap();
        let _result = sim.run_with_progress(5, |_p| {
            calls += 1;
        });

        // Called at start (0), then at steps 5 and 10.
        assert_eq!(calls, 3);
    }

    #[test]
    fn test_frequency_dependent_min_alive_fraction_early_termination() {
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut lib = MaterialLibrary::new();
        lib.add(Material::new("hi_abs").with_acoustic(AcousticMaterial::new(
            "hi_abs",
            [0.99, 0.99, 0.99, 0.99, 0.99, 0.99],
            [0.0; 6],
        )));
        lib.assign("/", "hi_abs");

        let mut config = SimulationConfig::new();
        config.num_steps = 10_000;
        config.num_rays = 50;
        config.source = Point::new(1.0, 1.0, 1.0);
        config.acoustic_mode = AcousticMode::FrequencyDependent;
        config.material_library = Some(lib);
        config.min_alive_fraction = 0.99;
        config.store_ray_history = false;

        let sim = Simulation::new(&building, config).unwrap();
        let result = sim.run();

        assert!(
            result.hits.len() < 10_000,
            "Should terminate early: ran {} steps",
            result.hits.len()
        );
        assert!(result.band_hits.is_some());
    }

    #[test]
    fn test_single_band_with_air_absorption() {
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config_no_air = SimulationConfig::new();
        config_no_air.num_steps = 100;
        config_no_air.num_rays = 20;
        config_no_air.source = Point::new(1.0, 1.0, 1.0);
        config_no_air.acoustic_mode = AcousticMode::FrequencyDependent;
        config_no_air.single_band_index = Some(5);
        config_no_air.enable_air_absorption = false;
        config_no_air.store_ray_history = true;

        let mut config_air = SimulationConfig::new();
        config_air.num_steps = 100;
        config_air.num_rays = 20;
        config_air.source = Point::new(1.0, 1.0, 1.0);
        config_air.acoustic_mode = AcousticMode::FrequencyDependent;
        config_air.single_band_index = Some(5);
        config_air.enable_air_absorption = true;
        config_air.store_ray_history = true;

        let sim_no_air = Simulation::new(&building, config_no_air).unwrap();
        let result_no_air = sim_no_air.run();

        let sim_air = Simulation::new(&building, config_air).unwrap();
        let result_air = sim_air.run();

        let energy_no_air: f64 = result_no_air.energies.last().unwrap().iter().sum();
        let energy_air: f64 = result_air.energies.last().unwrap().iter().sum();

        assert!(
            energy_air <= energy_no_air + 1e-10,
            "Air absorption should reduce or maintain energy: \
             no_air={energy_no_air:.6}, air={energy_air:.6}"
        );
    }

    #[test]
    fn test_single_band_with_progress_and_early_termination() {
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut lib = MaterialLibrary::new();
        lib.add(Material::new("hi_abs").with_acoustic(AcousticMaterial::new(
            "hi_abs",
            [0.99, 0.99, 0.99, 0.99, 0.99, 0.99],
            [0.0; 6],
        )));
        lib.assign("/", "hi_abs");

        let mut config = SimulationConfig::new();
        config.num_steps = 10_000;
        config.num_rays = 50;
        config.source = Point::new(1.0, 1.0, 1.0);
        config.acoustic_mode = AcousticMode::FrequencyDependent;
        config.single_band_index = Some(3);
        config.material_library = Some(lib);
        config.min_alive_fraction = 0.99;
        config.store_ray_history = true;
        config.store_ray_band_history = true;

        let mut report_count: usize = 0;
        let sim = Simulation::new(&building, config).unwrap();
        let result = sim.run_with_progress(5, |_p| {
            report_count += 1;
        });

        // Should terminate early
        assert!(
            result.hits.len() < 10_000,
            "Should terminate early: ran {} steps",
            result.hits.len()
        );
        assert!(report_count > 0, "Should have reported at least once");
        assert!(result.band_hits.is_some());
        assert!(result.band_energies.is_some());
    }

    #[test]
    fn test_single_band_absorber_records_only_target_band() {
        let s0 = Solid::from_box(3.0, 3.0, 3.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config = SimulationConfig::new();
        config.num_steps = 500;
        config.num_rays = 500;
        config.source = Point::new(1.5, 1.5, 1.5);
        config.absorbers = vec![Point::new(0.5, 0.5, 0.5)];
        config.absorber_radius = 0.8;
        config.acoustic_mode = AcousticMode::FrequencyDependent;
        config.single_band_index = Some(2);
        config.store_ray_history = false;

        let sim = Simulation::new(&building, config).unwrap();
        let result = sim.run();

        let band_hits = result.band_hits.as_ref().unwrap();
        let total_hit: f64 = result.hits.iter().map(|h| h[0]).sum();
        assert!(total_hit > 0.0, "Should have absorber hits");

        // Only band 2 should have non-zero hits
        for step_hits in band_hits {
            for absorber_bh in step_hits {
                for (b, &h) in absorber_bh.iter().enumerate() {
                    if b != 2 {
                        assert!(h.abs() < 1e-15, "Band {b} should have zero hits, got {h}");
                    }
                }
            }
        }
    }

    #[test]
    fn test_frequency_dependent_decay_ordering() {
        // With frequency-dependent absorption [0.05, 0.10, 0.20, 0.35, 0.50, 0.70],
        // higher bands should have less remaining energy after many steps.
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut lib = MaterialLibrary::new();
        lib.add(Material::new("walls").with_acoustic(AcousticMaterial::new(
            "walls",
            [0.05, 0.10, 0.20, 0.35, 0.50, 0.70],
            [0.1; 6],
        )));
        lib.assign("/", "walls");

        let mut config = SimulationConfig::new();
        config.num_steps = 2000;
        config.num_rays = 200;
        config.source = Point::new(1.0, 1.0, 1.0);
        config.acoustic_mode = AcousticMode::FrequencyDependent;
        config.material_library = Some(lib);
        config.store_ray_history = true;
        config.store_ray_band_history = true;

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
