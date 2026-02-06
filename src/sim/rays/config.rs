use std::collections::HashMap;

use crate::{Building, Point};

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

    // Surfaces
    pub default_absorption: f64,
    pub absorption: HashMap<String, f64>,
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
}

impl Default for SimulationConfig {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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
}
