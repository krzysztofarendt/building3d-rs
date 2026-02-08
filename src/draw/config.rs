/// RGBA color tuple (red, green, blue, alpha) with values in 0.0..=1.0.
pub type Rgba = (f32, f32, f32, f32);

/// Configuration for Rerun visualization sessions.
///
/// Controls session naming, entity prefixes, and default colors/sizes
/// for drawing functions and simulation visualization.
pub struct RerunConfig {
    // Labels
    pub session_name: String,
    pub entity_prefix: String,

    // General drawing defaults
    pub face_color: Rgba,
    pub edge_color: Rgba,
    pub edge_radius: f32,
    pub point_color: Rgba,
    pub point_radius: f32,

    // Simulation drawing
    pub sim_building_color: Rgba,
    pub sim_absorber_color: Rgba,
    pub sim_source_color: Rgba,
    pub sim_source_radius: f32,
    pub sim_ray_color_high: Rgba,
    pub sim_ray_color_low: Rgba,
    pub sim_ray_radius: f32,
    pub sim_ray_energy_threshold: f64,
}

impl RerunConfig {
    pub fn new() -> Self {
        Self {
            session_name: "building3d".to_string(),
            entity_prefix: "Building3d".to_string(),

            face_color: (1.0, 1.0, 1.0, 0.2),
            edge_color: (0.0, 0.0, 1.0, 0.5),
            edge_radius: 0.01,
            point_color: (0.0, 1.0, 0.0, 1.0),
            point_radius: 0.05,

            sim_building_color: (0.8, 0.8, 0.8, 0.15),
            sim_absorber_color: (1.0, 0.0, 0.0, 0.5),
            sim_source_color: (0.0, 1.0, 0.0, 1.0),
            sim_source_radius: 0.02,
            sim_ray_color_high: (1.0, 0.0, 0.0, 0.8),
            sim_ray_color_low: (1.0, 1.0, 1.0, 0.1),
            sim_ray_radius: 0.04,
            sim_ray_energy_threshold: 1e-10,
        }
    }
}

impl Default for RerunConfig {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_defaults() {
        let config = RerunConfig::new();
        assert_eq!(config.session_name, "building3d");
        assert_eq!(config.entity_prefix, "Building3d");
        assert_eq!(config.face_color, (1.0, 1.0, 1.0, 0.2));
        assert_eq!(config.sim_ray_color_high, (1.0, 0.0, 0.0, 0.8));
        assert_eq!(config.sim_ray_color_low, (1.0, 1.0, 1.0, 0.1));
        assert_eq!(config.sim_ray_energy_threshold, 1e-10);
    }

    #[test]
    fn test_default_trait() {
        let config = RerunConfig::default();
        assert_eq!(config.session_name, "building3d");
    }

    #[test]
    fn test_custom_values() {
        let mut config = RerunConfig::new();
        config.session_name = "my_session".to_string();
        config.entity_prefix = "MyPrefix".to_string();
        config.sim_ray_radius = 0.1;
        assert_eq!(config.session_name, "my_session");
        assert_eq!(config.entity_prefix, "MyPrefix");
        assert_eq!(config.sim_ray_radius, 0.1);
    }
}
