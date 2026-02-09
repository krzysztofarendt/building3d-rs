use anyhow::Result;
use rerun as rr;

use crate::sim::lighting::result::LightingResult;
use crate::sim::lighting::sensor::SensorGrid;
use crate::{Building, HasMesh};

use super::config::RerunConfig;
use super::rerun::{draw_faces, start_session};

const SESSION_NAME: &str = "Building3d";

/// Draws a lighting simulation result as a heatmap on the building surfaces.
///
/// Each polygon is colored based on its illuminance level:
/// blue (low) → green → yellow (high).
pub fn draw_illuminance_heatmap(
    session: &rr::RecordingStream,
    result: &LightingResult,
    building: &Building,
) -> Result<()> {
    let config = RerunConfig::default();

    // Draw building mesh as wireframe background
    draw_faces(session, building, (0.9, 0.9, 0.9, 0.1), &config)?;

    // Find illuminance range for normalization
    let mut max_lux = 0.0_f64;
    for ill in result.illuminance.values() {
        let total = ill[0] + ill[1] + ill[2];
        max_lux = max_lux.max(total / 3.0);
    }
    if max_lux < 1e-10 {
        max_lux = 1.0;
    }

    // Draw each polygon with color based on illuminance
    for zone in building.zones() {
        for solid in zone.solids() {
            for wall in solid.walls() {
                for polygon in wall.polygons() {
                    let path = format!(
                        "{}/{}/{}/{}",
                        zone.name, solid.name, wall.name, polygon.name
                    );

                    let lux = result
                        .illuminance
                        .get(&polygon.uid)
                        .map(|ill| (ill[0] + ill[1] + ill[2]) / 3.0)
                        .unwrap_or(0.0);

                    let t = (lux / max_lux).clamp(0.0, 1.0) as f32;
                    let (r, g, b) = lux_to_color(t);
                    let a = if polygon.name.contains("glass") {
                        0.15
                    } else {
                        0.8
                    };

                    let mesh = polygon.copy_mesh();
                    let vertices = mesh.vertices;
                    let triangles = mesh.faces.unwrap_or_default();

                    let name = format!("{}/illuminance/{}", SESSION_NAME, path);
                    session.log_static(
                        name,
                        &rr::Mesh3D::new(vertices)
                            .with_triangle_indices(triangles)
                            .with_albedo_factor(rr::Rgba32::from_linear_unmultiplied_rgba_f32(
                                r, g, b, a,
                            )),
                    )?;
                }
            }
        }
    }

    Ok(())
}

/// Draws a sensor grid as colored points.
pub fn draw_sensor_grid(
    session: &rr::RecordingStream,
    grid: &SensorGrid,
    max_lux: f64,
) -> Result<()> {
    if grid.sensors.is_empty() {
        return Ok(());
    }

    let max_lux = if max_lux < 1e-10 { 1.0 } else { max_lux };

    let mut pts = Vec::new();
    let mut colors = Vec::new();
    let mut radii = Vec::new();

    // Use log scale: map [0, max_lux] → [0, 1] via ln(1+x)/ln(1+max)
    let log_max = (1.0 + max_lux).ln();

    for sensor in &grid.sensors {
        pts.push(sensor.position);
        let lux = (sensor.illuminance[0] + sensor.illuminance[1] + sensor.illuminance[2]) / 3.0;
        let t = ((1.0 + lux).ln() / log_max).clamp(0.0, 1.0) as f32;
        let (r, g, b) = lux_to_color(t);
        colors.push(rr::Color(rr::Rgba32::from_linear_unmultiplied_rgba_f32(
            r, g, b, 1.0,
        )));
        radii.push(0.05_f32);
    }

    let name = format!("{}/sensors/{}", SESSION_NAME, grid.polygon_path);
    session.log_static(
        name,
        &rr::Points3D::new(pts).with_radii(radii).with_colors(colors),
    )?;

    Ok(())
}

/// Maps a normalized value (0-1) to a blue-green-yellow color ramp.
fn lux_to_color(t: f32) -> (f32, f32, f32) {
    if t < 0.5 {
        // Blue to green
        let s = t * 2.0;
        (0.0, s, 1.0 - s)
    } else {
        // Green to yellow
        let s = (t - 0.5) * 2.0;
        (s, 1.0, 0.0)
    }
}

/// Convenience: creates a session and draws the illuminance heatmap.
pub fn show_illuminance(result: &LightingResult, building: &Building) -> Result<()> {
    let config = RerunConfig::default();
    let session = start_session(&config)?;
    draw_illuminance_heatmap(&session, result, building)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::lighting::result::LightingResult;
    use crate::{Polygon, Solid, Wall, Zone};

    fn buffered_session() -> rr::RecordingStream {
        rr::RecordingStreamBuilder::new("test").buffered().unwrap()
    }

    #[test]
    fn test_lux_to_color_endpoints() {
        assert_eq!(lux_to_color(0.0), (0.0, 0.0, 1.0));
        assert_eq!(lux_to_color(0.5), (0.0, 1.0, 0.0));
        assert_eq!(lux_to_color(1.0), (1.0, 1.0, 0.0));
    }

    #[test]
    fn test_draw_sensor_grid_empty_is_ok() {
        let session = buffered_session();
        let grid = SensorGrid {
            sensors: Vec::new(),
            polygon_path: "zone/room/floor/floor".to_string(),
        };
        draw_sensor_grid(&session, &grid, 100.0).unwrap();
    }

    #[test]
    fn test_draw_sensor_grid_nonempty_is_ok() {
        let session = buffered_session();

        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        // Generate from an actual polygon so normals/positions are sensible.
        let poly = building.zones()[0].solids()[0].walls()[0].polygons()[0];
        let grid = SensorGrid::generate(poly, 0.5, "zone/room/floor/floor");

        draw_sensor_grid(&session, &grid, 0.0).unwrap(); // exercises max_lux floor
    }

    #[test]
    fn test_draw_illuminance_heatmap_is_ok() {
        let session = buffered_session();

        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let poly = building.zones()[0].solids()[0].walls()[0].polygons()[0];

        let mut result = LightingResult::new();
        result
            .illuminance
            .insert(poly.uid.clone(), [10.0, 20.0, 30.0]);
        draw_illuminance_heatmap(&session, &result, &building).unwrap();

        // Also exercise the max_lux floor branch (all zeros / empty map).
        let result = LightingResult::new();
        draw_illuminance_heatmap(&session, &result, &building).unwrap();
    }

    #[test]
    fn test_draw_illuminance_heatmap_glass_alpha_branch() {
        let session = buffered_session();

        let glass = Polygon::new(
            "glass_panel",
            vec![
                crate::Point::new(0.0, 0.0, 0.0),
                crate::Point::new(1.0, 0.0, 0.0),
                crate::Point::new(1.0, 1.0, 0.0),
                crate::Point::new(0.0, 1.0, 0.0),
            ],
            None,
        )
        .unwrap();
        let opaque = Polygon::new(
            "opaque_panel",
            vec![
                crate::Point::new(0.0, 0.0, 1.0),
                crate::Point::new(1.0, 0.0, 1.0),
                crate::Point::new(1.0, 1.0, 1.0),
                crate::Point::new(0.0, 1.0, 1.0),
            ],
            None,
        )
        .unwrap();

        let wall = Wall::new("w", vec![glass.clone(), opaque.clone()]).unwrap();
        let solid = Solid::new("s", vec![wall]).unwrap();
        let zone = Zone::new("z", vec![solid]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut result = LightingResult::new();
        result
            .illuminance
            .insert(glass.uid.clone(), [10.0, 10.0, 10.0]);
        result
            .illuminance
            .insert(opaque.uid.clone(), [10.0, 10.0, 10.0]);

        draw_illuminance_heatmap(&session, &result, &building).unwrap();
    }

    #[test]
    fn test_show_illuminance_is_ok_when_rerun_disabled() {
        let prev = std::env::var("RERUN").ok();
        unsafe { std::env::set_var("RERUN", "0") };

        let s0 = Solid::from_box(1.0, 1.0, 1.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();
        let result = LightingResult::new();

        show_illuminance(&result, &building).unwrap();

        match prev {
            Some(v) => unsafe { std::env::set_var("RERUN", v) },
            None => unsafe { std::env::remove_var("RERUN") },
        }
    }
}
