use anyhow::Result;
use rerun as rr;

use crate::sim::lighting::result::LightingResult;
use crate::sim::lighting::sensor::SensorGrid;
use crate::{Building, HasMesh};

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
    // Draw building mesh as wireframe background
    draw_faces(session, building, (0.9, 0.9, 0.9, 0.1))?;

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
                        .get(&path)
                        .map(|ill| (ill[0] + ill[1] + ill[2]) / 3.0)
                        .unwrap_or(0.0);

                    let t = (lux / max_lux).clamp(0.0, 1.0) as f32;
                    let (r, g, b) = lux_to_color(t);

                    let mesh = polygon.copy_mesh();
                    let vertices = mesh.vertices;
                    let triangles = mesh.faces.unwrap_or_default();

                    let name = format!("{}/illuminance/{}", SESSION_NAME, path);
                    session.log_static(
                        name,
                        &rr::Mesh3D::new(vertices)
                            .with_triangle_indices(triangles)
                            .with_albedo_factor(rr::Rgba32::from_linear_unmultiplied_rgba_f32(
                                r, g, b, 0.8,
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

    for sensor in &grid.sensors {
        pts.push(sensor.position);
        let lux = (sensor.illuminance[0] + sensor.illuminance[1] + sensor.illuminance[2]) / 3.0;
        let t = (lux / max_lux).clamp(0.0, 1.0) as f32;
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
    let session = start_session()?;
    draw_illuminance_heatmap(&session, result, building)?;
    Ok(())
}
