use std::collections::HashMap;

use crate::sim::materials::{AcousticMaterial, MaterialLibrary, NUM_OCTAVE_BANDS};
use crate::{Building, Point};

/// Acoustic simulation mode.
#[derive(Debug, Clone)]
pub enum AcousticMode {
    /// Scalar absorption: single coefficient per surface.
    Scalar,
    /// Frequency-dependent absorption: per-octave-band coefficients.
    FrequencyDependent,
}

pub struct SimulationConfig {
    // Engine
    pub time_step: f64,
    pub num_steps: usize,
    pub batch_size: usize,
    pub voxel_size: f64,
    pub search_transparent: bool,

    // Rays
    pub num_rays: usize,
    pub ray_speed: f64,
    pub source: Point,
    pub absorbers: Vec<Point>,
    pub absorber_radius: f64,

    // Surfaces (scalar mode)
    pub default_absorption: f64,
    pub absorption: HashMap<String, f64>,

    // Surfaces (frequency-dependent mode)
    pub acoustic_mode: AcousticMode,
    pub material_library: Option<MaterialLibrary>,
    pub default_acoustic_material: Option<AcousticMaterial>,
    /// If set (0..NUM_OCTAVE_BANDS), simulate a single octave band only.
    ///
    /// Intended mainly for visualization/debugging to reduce compute.
    pub single_band_index: Option<usize>,

    // Air absorption
    pub enable_air_absorption: bool,

    // Early termination
    /// Minimum fraction of rays that must still be alive (energy > eps) to continue.
    /// When the alive fraction drops below this threshold, the simulation ends early.
    /// Set to 0.0 to disable early termination. Default: 0.0 (disabled).
    pub min_alive_fraction: f64,

    // Memory
    /// If `true`, store per-step ray positions and energies in the result.
    /// Required for `draw_simulation()` but uses ~800 MB for 5000 rays × 2000 steps.
    /// Set to `false` when you only need absorber hits (e.g. for acoustic metrics).
    pub store_ray_history: bool,
    /// If `true` (and `store_ray_history=true`), store per-ray per-band energies
    /// for frequency-dependent simulations.
    ///
    /// This is expensive and usually unnecessary for visualization (which only
    /// needs scalar energies). Default: `false`.
    pub store_ray_band_history: bool,
}

impl SimulationConfig {
    pub fn new() -> Self {
        Self {
            time_step: 2.5e-5,
            num_steps: 1000,
            batch_size: 100,
            voxel_size: 0.1,
            search_transparent: true,
            num_rays: 1000,
            ray_speed: 343.0,
            source: Point::new(0.0, 0.0, 0.0),
            absorbers: Vec::new(),
            absorber_radius: 0.1,
            default_absorption: 0.2,
            absorption: HashMap::new(),
            acoustic_mode: AcousticMode::Scalar,
            material_library: None,
            default_acoustic_material: None,
            single_band_index: None,
            enable_air_absorption: false,
            min_alive_fraction: 0.0,
            store_ray_history: false,
            store_ray_band_history: false,
        }
    }

    /// Sets absorption for surfaces matching the given path pattern.
    ///
    /// Supports partial paths: if `path` is a substring of a polygon's full path
    /// (zone/solid/wall/polygon), it will match. For example, "wall_0" matches
    /// all polygons in wall_0 across all zones and solids.
    pub fn set_absorption(&mut self, path: &str, value: f64, building: &Building) {
        for zone in building.zones() {
            for solid in zone.solids() {
                for wall in solid.walls() {
                    for poly in wall.polygons() {
                        let full_path =
                            format!("{}/{}/{}/{}", zone.name, solid.name, wall.name, poly.name);
                        if full_path.contains(path) {
                            self.absorption.insert(full_path, value);
                        }
                    }
                }
            }
        }
    }

    /// Resolves frequency-dependent absorption coefficients for each polygon path.
    ///
    /// Returns a vector of `[f64; NUM_OCTAVE_BANDS]` indexed by polygon order.
    pub fn resolve_band_absorption(&self, paths: &[String]) -> Vec<[f64; NUM_OCTAVE_BANDS]> {
        let default = self
            .default_acoustic_material
            .as_ref()
            .map(|m| m.absorption)
            .unwrap_or([self.default_absorption; NUM_OCTAVE_BANDS]);

        paths
            .iter()
            .map(|path| {
                if let Some(ref lib) = self.material_library
                    && let Some(mat) = lib.lookup(path)
                    && let Some(ref acoustic) = mat.acoustic
                {
                    return acoustic.absorption;
                }
                default
            })
            .collect()
    }

    /// Resolves scattering coefficients for each polygon path.
    pub fn resolve_scattering(&self, paths: &[String]) -> Vec<[f64; NUM_OCTAVE_BANDS]> {
        let default = self
            .default_acoustic_material
            .as_ref()
            .map(|m| m.scattering)
            .unwrap_or([0.0; NUM_OCTAVE_BANDS]);

        paths
            .iter()
            .map(|path| {
                if let Some(ref lib) = self.material_library
                    && let Some(mat) = lib.lookup(path)
                    && let Some(ref acoustic) = mat.acoustic
                {
                    return acoustic.scattering;
                }
                default
            })
            .collect()
    }
}

impl Default for SimulationConfig {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::materials::{AcousticMaterial, Material, MaterialLibrary};
    use crate::{Solid, Zone};

    #[test]
    fn test_config_defaults() {
        let config = SimulationConfig::new();
        assert_eq!(config.num_steps, 1000);
        assert_eq!(config.num_rays, 1000);
        assert!((config.ray_speed - 343.0).abs() < 1e-10);
        assert!((config.default_absorption - 0.2).abs() < 1e-10);
    }

    #[test]
    fn test_config_default_trait() {
        let config: SimulationConfig = Default::default();
        assert_eq!(config.num_steps, 1000);
    }

    #[test]
    fn test_set_absorption() {
        let s = Solid::from_box(1.0, 1.0, 1.0, None, "box1").unwrap();
        let z = Zone::new("zone1", vec![s]).unwrap();
        let b = Building::new("building", vec![z]).unwrap();

        let mut config = SimulationConfig::new();
        config.set_absorption("floor", 0.9, &b);

        // Should have matched floor/floor polygon
        assert!(!config.absorption.is_empty());
        for (path, val) in &config.absorption {
            assert!(path.contains("floor"));
            assert!((*val - 0.9).abs() < 1e-10);
        }
    }

    #[test]
    fn test_resolve_band_absorption_with_library() {
        let mut lib = MaterialLibrary::new();
        lib.add(Material::new("carpet").with_acoustic(AcousticMaterial::new(
            "carpet",
            [0.02, 0.06, 0.14, 0.37, 0.60, 0.65],
            [0.40, 0.40, 0.40, 0.50, 0.50, 0.50],
        )));
        lib.assign("floor", "carpet");

        let mut config = SimulationConfig::new();
        config.acoustic_mode = AcousticMode::FrequencyDependent;
        config.material_library = Some(lib);

        let paths = vec![
            "zone/room/floor/floor_0".to_string(),
            "zone/room/wall/north".to_string(),
        ];
        let absorption = config.resolve_band_absorption(&paths);

        // Floor should get carpet absorption
        assert!((absorption[0][0] - 0.02).abs() < 1e-10);
        assert!((absorption[0][5] - 0.65).abs() < 1e-10);

        // Wall should get default
        assert!((absorption[1][0] - 0.2).abs() < 1e-10);
    }
}
