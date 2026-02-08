use std::collections::HashMap;

use super::sensor::SensorGrid;
use super::sources::Rgb;
use crate::Point;

/// Result of a lighting simulation.
pub struct LightingResult {
    /// Irradiance per polygon path [W/m^2] (average per channel on the surface).
    pub illuminance: HashMap<String, Rgb>,
    /// Total incident radiant power [W] hitting each polygon (per channel).
    pub incident_flux: HashMap<String, Rgb>,
    /// Number of ray hits per polygon.
    pub hit_count: HashMap<String, usize>,
    /// Sensor grids with recorded illuminance data.
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

    /// Records a ray hit on a sensor grid, distributing energy to the nearest sensor.
    pub fn record_sensor_hit(
        &mut self,
        grid_idx: usize,
        hit_pos: Point,
        energy: Rgb,
        sensor_area: f64,
    ) {
        if let Some(grid) = self.sensor_grids.get_mut(grid_idx) {
            // Find nearest sensor
            let mut best_idx = None;
            let mut best_dist2 = f64::INFINITY;
            for (i, sensor) in grid.sensors.iter().enumerate() {
                let d = hit_pos - sensor.position;
                let dist2 = d.dx * d.dx + d.dy * d.dy + d.dz * d.dz;
                if dist2 < best_dist2 {
                    best_dist2 = dist2;
                    best_idx = Some(i);
                }
            }
            if let Some(idx) = best_idx {
                let sensor = &mut grid.sensors[idx];
                // Convert flux to irradiance: E = flux / area
                sensor.illuminance[0] += energy[0] / sensor_area;
                sensor.illuminance[1] += energy[1] / sensor_area;
                sensor.illuminance[2] += energy[2] / sensor_area;
                sensor.hit_count += 1;
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
