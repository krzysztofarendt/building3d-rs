use crate::sim::materials::MaterialLibrary;

use super::sources::{DirectionalLight, PointLight, Rgb};

/// Configuration for a lighting simulation.
pub struct LightingConfig {
    /// Point light sources.
    pub point_lights: Vec<PointLight>,
    /// Directional light sources.
    pub directional_lights: Vec<DirectionalLight>,
    /// Number of rays per light source.
    pub num_rays: usize,
    /// Maximum number of bounces.
    pub max_bounces: usize,
    /// Minimum ray energy to continue tracing (sum of RGB).
    pub min_energy: f64,
    /// Voxel grid size for spatial acceleration.
    pub voxel_size: f64,
    /// Material library for optical properties.
    pub material_library: Option<MaterialLibrary>,
    /// Default diffuse reflectance if no material is assigned.
    pub default_reflectance: Rgb,
}

impl LightingConfig {
    pub fn new() -> Self {
        Self {
            point_lights: Vec::new(),
            directional_lights: Vec::new(),
            num_rays: 10000,
            max_bounces: 5,
            min_energy: 1e-6,
            voxel_size: 0.1,
            material_library: None,
            default_reflectance: [0.5, 0.5, 0.5],
        }
    }
}

impl Default for LightingConfig {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_config_defaults() {
        let config = LightingConfig::new();
        assert_eq!(config.num_rays, 10000);
        assert_eq!(config.max_bounces, 5);
        assert!(config.point_lights.is_empty());
    }
}
