use anyhow::Result;
use rerun as rr;

use crate::sim::energy::result::ThermalResult;
use crate::sim::energy::simulation::AnnualResult;
use crate::{Building, HasMesh};

use super::config::RerunConfig;
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
    let config = RerunConfig::default();

    // Draw building mesh as transparent background
    draw_faces(session, building, (0.9, 0.9, 0.9, 0.1), &config)?;

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

/// Draws annual simulation results as a time-series chart.
///
/// Logs hourly heating/cooling demands as scalar values on the "hour" timeline.
/// Uses SeriesLines to configure line colors and labels.
pub fn draw_annual_timeline(session: &rr::RecordingStream, result: &AnnualResult) -> Result<()> {
    // Use top-level paths (not under Building3d/) so Rerun auto-generates
    // a dedicated time series view separate from the 3D spatial view.
    let heating_path = "energy/heating_demand";
    let cooling_path = "energy/cooling_demand";

    // Configure series appearance (static, logged once)
    session.log_static(
        heating_path,
        &rr::SeriesLines::new()
            .with_colors([[255, 80, 80]])
            .with_names(["Heating (W)"])
            .with_widths([1.5]),
    )?;
    session.log_static(
        cooling_path,
        &rr::SeriesLines::new()
            .with_colors([[80, 80, 255]])
            .with_names(["Cooling (W)"])
            .with_widths([1.5]),
    )?;

    // Log hourly values
    for (hour, (heating, cooling)) in result
        .hourly_heating
        .iter()
        .zip(result.hourly_cooling.iter())
        .enumerate()
    {
        session.set_time_sequence("hour", hour as i64);
        session.log(heating_path, &rr::Scalars::single(*heating))?;
        session.log(cooling_path, &rr::Scalars::single(*cooling))?;
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Solid, Zone};

    fn buffered_session() -> rr::RecordingStream {
        rr::RecordingStreamBuilder::new("test").buffered().unwrap()
    }

    #[test]
    fn test_heat_loss_color_endpoints() {
        assert_eq!(heat_loss_color(0.0), (0.0, 0.0, 1.0));
        assert_eq!(heat_loss_color(0.5), (1.0, 1.0, 1.0));
        assert_eq!(heat_loss_color(1.0), (1.0, 0.0, 0.0));
    }

    #[test]
    fn test_draw_heat_loss_heatmap_is_ok() {
        let session = buffered_session();

        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "room").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        // Build a result with per-surface path keys.
        let mut result = ThermalResult::new();
        for zone in building.zones() {
            for solid in zone.solids() {
                for wall in solid.walls() {
                    for polygon in wall.polygons() {
                        let path = format!(
                            "{}/{}/{}/{}",
                            zone.name, solid.name, wall.name, polygon.name
                        );
                        result.surface_heat_loss.insert(path, -10.0);
                    }
                }
            }
        }

        draw_heat_loss_heatmap(&session, &result, &building).unwrap();

        // Also exercise max_loss floor branch.
        let result = ThermalResult::new();
        draw_heat_loss_heatmap(&session, &result, &building).unwrap();
    }

    #[test]
    fn test_draw_annual_timeline_is_ok() {
        let session = buffered_session();
        let result = AnnualResult {
            hourly_heating: vec![0.0, 100.0, 200.0],
            hourly_cooling: vec![50.0, 0.0, 25.0],
            annual_heating_kwh: 0.0,
            annual_cooling_kwh: 0.0,
            peak_heating: 200.0,
            peak_cooling: 50.0,
            monthly_heating_kwh: [0.0; 12],
            monthly_cooling_kwh: [0.0; 12],
            min_zone_temp_c: 0.0,
            max_zone_temp_c: 0.0,
            hourly_zone_temp_c: Vec::new(),
            num_zones: 1,
            per_zone_hourly_temp_c: Vec::new(),
            per_zone_min_temp_c: Vec::new(),
            per_zone_max_temp_c: Vec::new(),
        };
        draw_annual_timeline(&session, &result).unwrap();
    }
}
