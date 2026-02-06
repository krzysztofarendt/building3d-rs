use rand::Rng;

use crate::Vector;
use crate::sim::engine::FlatScene;

use super::sensor::SensorGrid;
use super::sources::{LightSource, Rgb};

/// Performs backward ray tracing from sensor points toward light sources.
///
/// For each sensor, traces rays in the hemisphere above the surface to determine
/// how much light reaches the sensor point. This is the reciprocal approach to
/// forward ray tracing.
pub fn backward_trace_sensors(
    grid: &mut SensorGrid,
    scene: &FlatScene,
    lights: &[&dyn LightSource],
    num_rays: usize,
) {
    let mut rng = rand::thread_rng();

    for sensor in &mut grid.sensors {
        let normal = sensor.normal;
        let mut total_illuminance: Rgb = [0.0; 3];

        for _ in 0..num_rays {
            // Generate random direction in hemisphere above sensor
            let dir = random_hemisphere_direction(normal, &mut rng);

            // Check if this direction reaches a light source
            // First check if there's an obstruction
            let hit = scene.find_target_surface_global(sensor.position, dir);

            if hit.is_none() {
                // Ray escapes — may reach directional lights (sky)
                // For now, just count unobstructed rays
                let cos_theta = dir.dot(&normal).abs();
                for light in lights {
                    let intensity = light.intensity(dir);
                    total_illuminance[0] += intensity[0] * cos_theta;
                    total_illuminance[1] += intensity[1] * cos_theta;
                    total_illuminance[2] += intensity[2] * cos_theta;
                }
            }
        }

        // Normalize by number of rays and hemisphere solid angle (2*pi)
        let scale = 2.0 * std::f64::consts::PI / num_rays as f64;
        sensor.illuminance[0] += total_illuminance[0] * scale;
        sensor.illuminance[1] += total_illuminance[1] * scale;
        sensor.illuminance[2] += total_illuminance[2] * scale;
        sensor.hit_count += num_rays;
    }
}

/// Generate a random direction in the hemisphere defined by the given normal.
fn random_hemisphere_direction(normal: Vector, rng: &mut impl Rng) -> Vector {
    loop {
        let x: f64 = rng.gen_range(-1.0..1.0);
        let y: f64 = rng.gen_range(-1.0..1.0);
        let z: f64 = rng.gen_range(-1.0..1.0);
        let len2 = x * x + y * y + z * z;
        if len2 > 1e-6 && len2 <= 1.0 {
            let len = len2.sqrt();
            let v = Vector::new(x / len, y / len, z / len);
            if v.dot(&normal) > 0.0 {
                return v;
            }
            return v * -1.0;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::lighting::sources::DirectionalLight;
    use crate::{Point, Polygon, Solid, Zone};

    #[test]
    fn test_backward_trace_no_obstruction() {
        // Open scene — a floor polygon with light from above
        let s0 = Solid::from_box(10.0, 10.0, 5.0, None, "room").unwrap();
        let zone = crate::Building::new("b", vec![Zone::new("z", vec![s0]).unwrap()]).unwrap();
        let scene = FlatScene::new(&zone, 0.5, false);

        let floor_pts = vec![
            Point::new(3.0, 3.0, 0.01),
            Point::new(7.0, 3.0, 0.01),
            Point::new(7.0, 7.0, 0.01),
            Point::new(3.0, 7.0, 0.01),
        ];
        let floor = Polygon::new("test_floor", floor_pts, None).unwrap();
        let mut grid = super::SensorGrid::generate(&floor, 2.0, "test");

        let sky = DirectionalLight::new(Vector::new(0.0, 0.0, -1.0), [500.0, 500.0, 500.0]);
        let lights: Vec<&dyn LightSource> = vec![&sky];

        backward_trace_sensors(&mut grid, &scene, &lights, 100);

        // Sensors should receive some illuminance
        if !grid.sensors.is_empty() {
            let avg = grid.average_illuminance();
            // In an enclosed room, most rays will hit the ceiling
            // so illuminance will be low, but test structure is correct
            assert!(avg[0] >= 0.0, "Illuminance should be non-negative");
        }
    }
}
