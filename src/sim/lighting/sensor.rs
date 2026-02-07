use crate::{Point, Polygon, Vector};

use super::sources::Rgb;

/// A single sensor point on a surface.
#[derive(Debug, Clone)]
pub struct Sensor {
    pub position: Point,
    pub normal: Vector,
    /// Accumulated irradiance [R, G, B] in W/m^2 per channel.
    pub illuminance: Rgb,
    /// Number of ray hits.
    pub hit_count: usize,
}

/// A grid of sensors on a polygon surface.
pub struct SensorGrid {
    pub sensors: Vec<Sensor>,
    /// Polygon path this grid belongs to.
    pub polygon_path: String,
}

impl SensorGrid {
    /// Generates a sensor grid on a polygon with the given spacing.
    ///
    /// Places sensors in a rectangular grid aligned with the polygon's local
    /// coordinate system. Only sensors that fall inside the polygon are kept.
    pub fn generate(polygon: &Polygon, spacing: f64, path: &str) -> Self {
        let vertices = polygon.vertices();
        if vertices.len() < 3 {
            return Self {
                sensors: Vec::new(),
                polygon_path: path.to_string(),
            };
        }

        let normal = polygon.vn;

        // Build local coordinate system on the polygon plane
        let edge = vertices[1] - vertices[0];
        let u_axis = match edge.normalize() {
            Ok(v) => v,
            Err(_) => {
                return Self {
                    sensors: Vec::new(),
                    polygon_path: path.to_string(),
                };
            }
        };
        let v_axis = normal.cross(&u_axis);

        // Project vertices to local 2D coordinates
        let origin = vertices[0];
        let local_pts: Vec<(f64, f64)> = vertices
            .iter()
            .map(|p| {
                let d = *p - origin;
                (d.dot(&u_axis), d.dot(&v_axis))
            })
            .collect();

        // Find bounding box in local coordinates
        let u_min = local_pts.iter().map(|p| p.0).fold(f64::INFINITY, f64::min);
        let u_max = local_pts
            .iter()
            .map(|p| p.0)
            .fold(f64::NEG_INFINITY, f64::max);
        let v_min = local_pts.iter().map(|p| p.1).fold(f64::INFINITY, f64::min);
        let v_max = local_pts
            .iter()
            .map(|p| p.1)
            .fold(f64::NEG_INFINITY, f64::max);

        let mut sensors = Vec::new();

        let mut u = u_min + spacing * 0.5;
        while u <= u_max {
            let mut v = v_min + spacing * 0.5;
            while v <= v_max {
                // Convert back to 3D
                let point = origin + u_axis * u + v_axis * v;

                // Check if the point is inside the polygon
                if polygon.is_point_inside(point, true) {
                    sensors.push(Sensor {
                        position: point,
                        normal,
                        illuminance: [0.0; 3],
                        hit_count: 0,
                    });
                }
                v += spacing;
            }
            u += spacing;
        }

        Self {
            sensors,
            polygon_path: path.to_string(),
        }
    }

    /// Returns the average irradiance [W/m^2] across all sensors.
    pub fn average_illuminance(&self) -> Rgb {
        if self.sensors.is_empty() {
            return [0.0; 3];
        }
        let mut sum = [0.0; 3];
        for s in &self.sensors {
            sum[0] += s.illuminance[0];
            sum[1] += s.illuminance[1];
            sum[2] += s.illuminance[2];
        }
        let n = self.sensors.len() as f64;
        [sum[0] / n, sum[1] / n, sum[2] / n]
    }

    /// Computes the daylight factor for each sensor.
    ///
    /// Daylight factor = interior irradiance / unobstructed exterior irradiance.
    pub fn daylight_factors(&self, exterior_illuminance: f64) -> Vec<f64> {
        if exterior_illuminance <= 0.0 {
            return vec![0.0; self.sensors.len()];
        }
        self.sensors
            .iter()
            .map(|s| {
                let total = s.illuminance[0] + s.illuminance[1] + s.illuminance[2];
                total / (3.0 * exterior_illuminance)
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::Result;

    #[test]
    fn test_sensor_grid_generation() -> Result<()> {
        // Create a simple floor polygon (2x2 at z=0)
        let pts = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(2.0, 0.0, 0.0),
            Point::new(2.0, 2.0, 0.0),
            Point::new(0.0, 2.0, 0.0),
        ];
        let poly = Polygon::new("floor", pts, None)?;

        let grid = SensorGrid::generate(&poly, 0.5, "zone/room/floor/floor");

        // With 0.5 spacing on 2x2 surface, expect roughly 4x4 = 16 sensors
        // (some at edges might not be inside)
        assert!(
            !grid.sensors.is_empty(),
            "Should generate some sensor points"
        );
        assert!(
            grid.sensors.len() >= 4,
            "Should have at least a few sensors"
        );

        // All sensors should have z ≈ 0
        for s in &grid.sensors {
            assert!((s.position.z - 0.0).abs() < 1e-6);
        }

        Ok(())
    }

    #[test]
    fn test_daylight_factor() -> Result<()> {
        let pts = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
        ];
        let poly = Polygon::new("floor", pts, None)?;

        let mut grid = SensorGrid::generate(&poly, 0.5, "zone/room/floor/floor");

        // Set some illuminance
        for s in &mut grid.sensors {
            s.illuminance = [100.0, 100.0, 100.0];
        }

        let df = grid.daylight_factors(10000.0);
        for &d in &df {
            assert!((d - 0.01).abs() < 1e-6, "DF should be 100/10000 = 0.01");
        }

        Ok(())
    }

    #[test]
    fn test_average_illuminance() -> Result<()> {
        let pts = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
        ];
        let poly = Polygon::new("floor", pts, None)?;

        let mut grid = SensorGrid::generate(&poly, 0.5, "test");
        for s in &mut grid.sensors {
            s.illuminance = [200.0, 300.0, 400.0];
        }
        let avg = grid.average_illuminance();
        assert!((avg[0] - 200.0).abs() < 1e-10);
        assert!((avg[1] - 300.0).abs() < 1e-10);

        Ok(())
    }
}
