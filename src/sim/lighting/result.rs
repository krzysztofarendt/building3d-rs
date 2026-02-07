use std::collections::HashMap;

use super::sensor::SensorGrid;
use super::sources::Rgb;

/// Result of a lighting simulation.
pub struct LightingResult {
    /// Irradiance per polygon path [W/m^2] (average per channel on the surface).
    pub illuminance: HashMap<String, Rgb>,
    /// Total incident radiant power [W] hitting each polygon (per channel).
    pub incident_flux: HashMap<String, Rgb>,
    /// Number of ray hits per polygon.
    pub hit_count: HashMap<String, usize>,
    /// Per-point sensor grids (populated when sensor_spacing is set).
    pub sensor_grids: Vec<SensorGrid>,
}

impl LightingResult {
    pub fn new() -> Self {
        Self {
            illuminance: HashMap::new(),
            incident_flux: HashMap::new(),
            hit_count: HashMap::new(),
            sensor_grids: Vec::new(),
        }
    }

    /// Records a ray hit on a polygon.
    pub fn record_hit(&mut self, path: &str, energy: Rgb) {
        let flux = self
            .incident_flux
            .entry(path.to_string())
            .or_insert([0.0; 3]);
        flux[0] += energy[0];
        flux[1] += energy[1];
        flux[2] += energy[2];

        *self.hit_count.entry(path.to_string()).or_insert(0) += 1;
    }

    /// Records a ray hit on a sensor grid, accumulating energy on the nearest sensor.
    pub fn record_sensor_hit(
        &mut self,
        grid_idx: usize,
        hit_pos: crate::Point,
        energy: Rgb,
        sensor_area: f64,
    ) {
        if grid_idx >= self.sensor_grids.len() {
            return;
        }
        let grid = &mut self.sensor_grids[grid_idx];
        if grid.sensors.is_empty() || sensor_area <= 0.0 {
            return;
        }
        // Find nearest sensor
        let mut best_i = 0;
        let mut best_d2 = f64::INFINITY;
        for (i, s) in grid.sensors.iter().enumerate() {
            let dx = s.position.x - hit_pos.x;
            let dy = s.position.y - hit_pos.y;
            let dz = s.position.z - hit_pos.z;
            let d2 = dx * dx + dy * dy + dz * dz;
            if d2 < best_d2 {
                best_d2 = d2;
                best_i = i;
            }
        }
        let s = &mut grid.sensors[best_i];
        s.illuminance[0] += energy[0] / sensor_area;
        s.illuminance[1] += energy[1] / sensor_area;
        s.illuminance[2] += energy[2] / sensor_area;
        s.hit_count += 1;
    }

    /// Computes irradiance from accumulated flux and polygon areas.
    pub fn compute_illuminance(&mut self, areas: &HashMap<String, f64>) {
        for (path, flux) in &self.incident_flux {
            if let Some(&area) = areas.get(path).filter(|&&a| a > 0.0) {
                self.illuminance.insert(
                    path.clone(),
                    [flux[0] / area, flux[1] / area, flux[2] / area],
                );
            }
        }
    }

    /// Returns the average irradiance [W/m^2] across all surfaces.
    pub fn average_illuminance(&self) -> Rgb {
        if self.illuminance.is_empty() {
            return [0.0; 3];
        }
        let mut sum = [0.0; 3];
        for ill in self.illuminance.values() {
            sum[0] += ill[0];
            sum[1] += ill[1];
            sum[2] += ill[2];
        }
        let n = self.illuminance.len() as f64;
        [sum[0] / n, sum[1] / n, sum[2] / n]
    }
}

impl Default for LightingResult {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_record_and_compute() {
        let mut result = LightingResult::new();
        result.record_hit("zone/room/floor/f0", [100.0, 100.0, 100.0]);
        result.record_hit("zone/room/floor/f0", [50.0, 50.0, 50.0]);

        let mut areas = HashMap::new();
        areas.insert("zone/room/floor/f0".to_string(), 10.0);
        result.compute_illuminance(&areas);

        let ill = result.illuminance.get("zone/room/floor/f0").unwrap();
        assert!((ill[0] - 15.0).abs() < 1e-10); // 150 / 10
    }
}
