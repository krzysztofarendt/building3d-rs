use anyhow::Result;
use rerun as rr;

use crate::sim::energy::result::ThermalResult;
use crate::sim::energy::simulation::AnnualResult;
use crate::{Building, HasMesh};

use super::rerun::draw_faces;

const SESSION_NAME: &str = "Building3d";

/// Draws a thermal simulation result as a heatmap on the building surfaces.
///
/// Each polygon is colored based on its heat loss:
/// blue (low loss) → red (high loss).
pub fn draw_heat_loss_heatmap(
    session: &rr::RecordingStream,
    result: &ThermalResult,
    building: &Building,
) -> Result<()> {
    // Draw building mesh as transparent background
    draw_faces(session, building, (0.9, 0.9, 0.9, 0.1))?;

    // Find max heat loss for normalization
    let mut max_loss = 0.0_f64;
    for loss in result.surface_heat_loss.values() {
        max_loss = max_loss.max(loss.abs());
    }
    if max_loss < 1e-10 {
        max_loss = 1.0;
    }

    // Draw each polygon colored by heat loss
    for zone in building.zones() {
        for solid in zone.solids() {
            for wall in solid.walls() {
                for polygon in wall.polygons() {
                    let path = format!(
                        "{}/{}/{}/{}",
                        zone.name, solid.name, wall.name, polygon.name
                    );

                    let loss = result.surface_heat_loss.get(&path).copied().unwrap_or(0.0);
                    let t = (loss.abs() / max_loss).clamp(0.0, 1.0) as f32;
                    let (r, g, b) = heat_loss_color(t);

                    let mesh = polygon.copy_mesh();
                    let vertices = mesh.vertices;
                    let triangles = mesh.faces.unwrap_or_default();

                    let name = format!("{}/heat_loss/{}", SESSION_NAME, path);
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

/// Draws annual simulation results as an animated temperature timeline.
///
/// Logs hourly heating/cooling demands as scalar values to Rerun's timeline.
pub fn draw_annual_timeline(session: &rr::RecordingStream, result: &AnnualResult) -> Result<()> {
    for (hour, (heating, cooling)) in result
        .hourly_heating
        .iter()
        .zip(result.hourly_cooling.iter())
        .enumerate()
    {
        session.set_time_sequence("hour", hour as i64);

        session.log(
            format!("{}/energy/heating_demand", SESSION_NAME),
            &rr::Scalars::new([*heating]),
        )?;
        session.log(
            format!("{}/energy/cooling_demand", SESSION_NAME),
            &rr::Scalars::new([*cooling]),
        )?;
    }

    Ok(())
}

/// Maps a normalized heat loss value (0-1) to a blue-red color ramp.
fn heat_loss_color(t: f32) -> (f32, f32, f32) {
    // Blue (low loss) → White (medium) → Red (high loss)
    if t < 0.5 {
        let s = t * 2.0;
        (s, s, 1.0)
    } else {
        let s = (t - 0.5) * 2.0;
        (1.0, 1.0 - s, 1.0 - s)
    }
}
