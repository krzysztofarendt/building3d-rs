use std::collections::HashMap;

use crate::UID;

use super::sensor::SensorGrid;
use super::sources::Rgb;
use crate::Point;

/// Result of a lighting simulation.
///
/// Cross-module-friendly outputs are keyed by polygon `UID`, with optional path
/// strings available via [`crate::sim::index::SurfaceIndex`] for reporting.
pub struct LightingResult {
    /// Irradiance per polygon [W/m^2] (average per channel on the surface).
    pub illuminance: HashMap<UID, Rgb>,
    /// Total incident radiant power [W] hitting each polygon (per channel).
    pub incident_flux: HashMap<UID, Rgb>,
    /// Number of ray hits per polygon.
    pub hit_count: HashMap<UID, usize>,
    /// Sensor grids with recorded illuminance data.
    pub sensor_grids: Vec<SensorGrid>,

    // Accumulators keyed by polygon index (to avoid cloning UID per ray hit).
    incident_flux_by_index: HashMap<usize, Rgb>,
    hit_count_by_index: HashMap<usize, usize>,
}

impl LightingResult {
    pub fn new() -> Self {
        Self {
            illuminance: HashMap::new(),
            incident_flux: HashMap::new(),
            hit_count: HashMap::new(),
            sensor_grids: Vec::new(),
            incident_flux_by_index: HashMap::new(),
            hit_count_by_index: HashMap::new(),
        }
    }

    /// Records a ray hit on the polygon at `polygon_index`.
    pub fn record_hit_index(&mut self, polygon_index: usize, energy: Rgb) {
        let flux = self
            .incident_flux_by_index
            .entry(polygon_index)
            .or_insert([0.0; 3]);
        flux[0] += energy[0];
        flux[1] += energy[1];
        flux[2] += energy[2];

        *self.hit_count_by_index.entry(polygon_index).or_insert(0) += 1;
    }

    /// Finalizes UID-keyed maps from the internal accumulators.
    ///
    /// - `polygon_uids_by_index` must be aligned with the scene polygon list
    /// - `areas_m2_by_index` is used to compute irradiance
    pub fn finalize(&mut self, polygon_uids_by_index: &[UID], areas_m2_by_index: &[f64]) {
        self.incident_flux.clear();
        self.hit_count.clear();
        self.illuminance.clear();

        for (idx, flux) in &self.incident_flux_by_index {
            let Some(uid) = polygon_uids_by_index.get(*idx) else {
                continue;
            };
            self.incident_flux.insert(uid.clone(), *flux);

            if let Some(hits) = self.hit_count_by_index.get(idx) {
                self.hit_count.insert(uid.clone(), *hits);
            }

            if let Some(&area) = areas_m2_by_index
                .get(*idx)
                .filter(|a| a.is_finite() && **a > 0.0)
            {
                self.illuminance.insert(
                    uid.clone(),
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
                sensor.illuminance[0] += energy[0] / sensor_area;
                sensor.illuminance[1] += energy[1] / sensor_area;
                sensor.illuminance[2] += energy[2] / sensor_area;
                sensor.hit_count += 1;
            }
        }
    }

    /// Returns the average irradiance [W/m^2] across all recorded polygons.
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
    use crate::Vector;
    use crate::sim::lighting::sensor::{Sensor, SensorGrid};

    #[test]
    fn test_default_trait() {
        let result: LightingResult = Default::default();
        assert!(result.illuminance.is_empty());
        assert!(result.incident_flux.is_empty());
        assert!(result.hit_count.is_empty());
    }

    #[test]
    fn test_record_and_finalize() {
        let u0 = UID::from("poly-0");
        let mut result = LightingResult::new();

        result.record_hit_index(0, [100.0, 100.0, 100.0]);
        result.record_hit_index(0, [50.0, 50.0, 50.0]);
        result.finalize(&[u0.clone()], &[10.0]);

        let flux = result.incident_flux.get(&u0).unwrap();
        assert!((flux[0] - 150.0).abs() < 1e-10);

        let ill = result.illuminance.get(&u0).unwrap();
        assert!((ill[0] - 15.0).abs() < 1e-10); // 150 / 10

        let hits = result.hit_count.get(&u0).unwrap();
        assert_eq!(*hits, 2);
    }

    #[test]
    fn test_record_sensor_hit_updates_nearest_sensor() {
        let mut result = LightingResult::new();
        result.sensor_grids.push(SensorGrid {
            polygon_path: "zone/room/floor/floor".to_string(),
            sensors: vec![
                Sensor {
                    position: Point::new(0.0, 0.0, 0.0),
                    normal: Vector::new(0.0, 0.0, 1.0),
                    illuminance: [0.0; 3],
                    hit_count: 0,
                },
                Sensor {
                    position: Point::new(10.0, 0.0, 0.0),
                    normal: Vector::new(0.0, 0.0, 1.0),
                    illuminance: [0.0; 3],
                    hit_count: 0,
                },
            ],
        });

        // Hit near the first sensor.
        result.record_sensor_hit(0, Point::new(0.1, 0.0, 0.0), [3.0, 6.0, 9.0], 3.0);
        let g = &result.sensor_grids[0];
        assert_eq!(g.sensors[0].hit_count, 1);
        assert!(g.sensors[1].hit_count == 0);
        assert!((g.sensors[0].illuminance[0] - 1.0).abs() < 1e-12);
        assert!((g.sensors[0].illuminance[2] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn test_average_illuminance_empty_and_nonempty() {
        let result = LightingResult::new();
        assert_eq!(result.average_illuminance(), [0.0; 3]);

        let u0 = UID::from("poly-0");
        let u1 = UID::from("poly-1");
        let mut result = LightingResult::new();
        result.illuminance.insert(u0, [1.0, 2.0, 3.0]);
        result.illuminance.insert(u1, [3.0, 4.0, 5.0]);
        let avg = result.average_illuminance();
        assert!((avg[0] - 2.0).abs() < 1e-12);
        assert!((avg[2] - 4.0).abs() < 1e-12);
    }
}
