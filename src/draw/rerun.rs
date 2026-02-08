use crate::HasName;
use crate::Point;
use crate::geom::triangles::TriangleIndex;
use crate::sim::materials::NUM_OCTAVE_BANDS;
use crate::sim::rays::AcousticMode;
use crate::sim::rays::SimulationResult;
use crate::{Building, HasMesh, Mesh};
use anyhow::Result;
use rerun as rr;

use super::config::{RerunConfig, Rgba, SimRayColormap, SimRayEnergyScale};

/// Converts Point to native format of Rerun
impl From<Point> for rr::Vec3D {
    fn from(val: Point) -> Self {
        rr::Vec3D([val.x as f32, val.y as f32, val.z as f32])
    }
}

/// Converts TriangleIndex to native format of Rerun
impl From<TriangleIndex> for rr::TriangleIndices {
    fn from(val: TriangleIndex) -> Self {
        rr::TriangleIndices(rr::datatypes::UVec3D([
            val.0 as u32,
            val.1 as u32,
            val.2 as u32,
        ]))
    }
}

fn color(rgba: Rgba) -> rr::Color {
    let (r, g, b, a) = rgba;
    rr::Color(rr::Rgba32::from_linear_unmultiplied_rgba_f32(r, g, b, a))
}

/// Linearly interpolate between two RGBA colors.
fn lerp_color(low: Rgba, high: Rgba, t: f32) -> Rgba {
    let t = t.clamp(0.0, 1.0);
    (
        low.0 + (high.0 - low.0) * t,
        low.1 + (high.1 - low.1) * t,
        low.2 + (high.2 - low.2) * t,
        low.3 + (high.3 - low.3) * t,
    )
}

fn hsv_to_rgba(h_deg: f32, s: f32, v: f32, a: f32) -> Rgba {
    let h = (h_deg.rem_euclid(360.0) / 60.0).clamp(0.0, 6.0);
    let c = v * s;
    let x = c * (1.0 - ((h % 2.0) - 1.0).abs());
    let (r1, g1, b1) = match h as u32 {
        0 => (c, x, 0.0),
        1 => (x, c, 0.0),
        2 => (0.0, c, x),
        3 => (0.0, x, c),
        4 => (x, 0.0, c),
        _ => (c, 0.0, x),
    };
    let m = v - c;
    (r1 + m, g1 + m, b1 + m, a)
}

fn energy_to_color_norm(energy: f64, config: &RerunConfig) -> f64 {
    let e = energy.clamp(0.0, 1.0);
    match config.sim_ray_energy_scale {
        SimRayEnergyScale::Linear => e,
        SimRayEnergyScale::Log => {
            let min = config
                .sim_ray_color_energy_min
                .max(config.sim_ray_energy_threshold)
                .clamp(1e-300, 1.0);
            let e = e.clamp(min, 1.0);
            let denom = 0.0 - min.ln();
            if denom <= 0.0 {
                return e;
            }

            // u: 1.0 at full energy, 0.0 at min energy.
            let u = ((e.ln() - min.ln()) / denom).clamp(0.0, 1.0);
            // Curve in "loss" space to make early changes more visible.
            let gamma = config.sim_ray_color_gamma.clamp(1e-6, 100.0);
            let loss = (1.0 - u).powf(gamma);
            (1.0 - loss).clamp(0.0, 1.0)
        }
    }
}

fn ray_color(energy: f64, config: &RerunConfig) -> rr::Color {
    let t = energy_to_color_norm(energy, config) as f32; // 1.0 = high energy
    match config.sim_ray_colormap {
        SimRayColormap::Lerp => {
            let c = lerp_color(config.sim_ray_color_low, config.sim_ray_color_high, t);
            color(c)
        }
        SimRayColormap::Rainbow => {
            // Red (high energy) -> ... -> Blue (low energy)
            let hue = (1.0 - t) * 240.0;
            // Use the "high" alpha as the common alpha.
            let a = config.sim_ray_color_high.3;
            color(hsv_to_rgba(hue, 1.0, 1.0, a))
        }
    }
}

fn energy_full_scale(result: &SimulationResult) -> f64 {
    match result.config.acoustic_mode {
        AcousticMode::Scalar => 1.0,
        AcousticMode::FrequencyDependent => {
            if result.config.single_band_index.is_some() {
                1.0
            } else {
                NUM_OCTAVE_BANDS as f64
            }
        }
    }
}

pub fn start_session(config: &RerunConfig) -> Result<rr::RecordingStream> {
    // Connect to the Rerun gRPC server using the default address and port: localhost:9876
    let session = rr::RecordingStreamBuilder::new(config.session_name.as_str()).spawn()?;

    Ok(session)
}

pub fn draw_faces<T: HasMesh + HasName>(
    session: &rr::RecordingStream,
    model: &T,
    rgba: Rgba,
    config: &RerunConfig,
) -> Result<()> {
    let mesh: Mesh = model.copy_mesh();
    let vertices: Vec<Point> = mesh.vertices;
    let triangles: Vec<TriangleIndex> = mesh.faces.unwrap_or_default();

    let (r, g, b, a) = rgba;

    let name = format!("{}/{}", config.entity_prefix, model.get_name());
    session.log_static(
        name,
        &rr::Mesh3D::new(vertices)
            .with_triangle_indices(triangles)
            .with_albedo_factor(rr::Rgba32::from_linear_unmultiplied_rgba_f32(r, g, b, a)),
    )?;

    Ok(())
}

pub fn draw_edges<T: HasMesh + HasName>(
    session: &rr::RecordingStream,
    model: &T,
    radius: f32,
    rgba: Rgba,
    config: &RerunConfig,
) -> Result<()> {
    let mesh: Mesh = model.copy_mesh();
    let vertices: Vec<Point> = mesh.vertices;
    let triangles: Vec<TriangleIndex> = mesh.faces.unwrap_or_default();

    let mut lines: Vec<Vec<rr::Vec3D>> = Vec::new();
    let mut radii: Vec<f32> = Vec::new();
    let mut colors: Vec<rr::Color> = Vec::new();

    for t in triangles.iter() {
        lines.push(Vec::new());
        let index = lines.len() - 1;
        lines[index].push(rr::Vec3D::from(vertices[t.0]));
        lines[index].push(rr::Vec3D::from(vertices[t.1]));
        lines[index].push(rr::Vec3D::from(vertices[t.2]));
        lines[index].push(rr::Vec3D::from(vertices[t.0]));
        radii.push(radius);
        colors.push(color(rgba));
    }

    let name = format!("{}/{}", config.entity_prefix, model.get_name());
    session.log_static(
        name,
        &rr::LineStrips3D::new(lines)
            .with_radii(radii)
            .with_colors(colors),
    )?;

    Ok(())
}

pub fn draw_points<T: HasMesh + HasName>(
    session: &rr::RecordingStream,
    model: &T,
    radius: f32,
    rgba: Rgba,
    config: &RerunConfig,
) -> Result<()> {
    let mesh: Mesh = model.copy_mesh();
    let vertices: Vec<Point> = mesh.vertices;

    let mut radii: Vec<f32> = Vec::new();
    let mut colors: Vec<rr::Color> = Vec::new();

    for _ in 0..vertices.len() {
        radii.push(radius);
        colors.push(color(rgba));
    }

    let name = format!("{}/{}", config.entity_prefix, model.get_name());
    session.log_static(
        name,
        &rr::Points3D::new(vertices)
            .with_radii(radii)
            .with_colors(colors),
    )?;

    Ok(())
}

/// Draws a ray tracing simulation result with animated ray positions.
///
/// The building mesh is drawn as semi-transparent static geometry.
/// Ray positions are logged per time step using Rerun's timeline for animation.
/// Ray color intensity reflects remaining energy.
pub fn draw_simulation(
    session: &rr::RecordingStream,
    result: &SimulationResult,
    building: &Building,
    config: &RerunConfig,
) -> Result<()> {
    // Draw building mesh (semi-transparent)
    draw_faces(session, building, config.sim_building_color, config)?;

    // Draw absorbers as static red spheres
    for (i, absorber) in result.config.absorbers.iter().enumerate() {
        let name = format!("{}/absorbers/{}", config.entity_prefix, i);
        session.log_static(
            name,
            &rr::Points3D::new([*absorber])
                .with_radii([result.config.absorber_radius as f32])
                .with_colors([color(config.sim_absorber_color)]),
        )?;
    }

    // Draw source as a static green sphere
    session.log_static(
        format!("{}/source", config.entity_prefix),
        &rr::Points3D::new([result.config.source])
            .with_radii([config.sim_source_radius])
            .with_colors([color(config.sim_source_color)]),
    )?;

    // Animate ray positions over time steps
    let energy_scale = energy_full_scale(result).max(1e-300);
    for (step, (positions, energies)) in result
        .positions
        .iter()
        .zip(result.energies.iter())
        .enumerate()
    {
        if config.sim_playback_dt_s > 0.0 {
            session.set_duration_secs("step", (step as f64) * config.sim_playback_dt_s.max(0.0));
        } else {
            session.set_time_sequence("step", step as i64);
        }

        // Collect alive rays (energy > threshold)
        let mut pts: Vec<Point> = Vec::new();
        let mut colors_vec: Vec<rr::Color> = Vec::new();
        let mut radii_vec: Vec<f32> = Vec::new();

        for (pos, &energy) in positions.iter().zip(energies.iter()) {
            if energy > config.sim_ray_energy_threshold {
                pts.push(*pos);
                colors_vec.push(ray_color(energy / energy_scale, config));
                radii_vec.push(config.sim_ray_radius);
            }
        }

        if !pts.is_empty() {
            session.log(
                format!("{}/rays", config.entity_prefix),
                &rr::Points3D::new(pts)
                    .with_radii(radii_vec)
                    .with_colors(colors_vec),
            )?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lerp_color_zero() {
        let low = (0.0, 0.0, 0.0, 0.0);
        let high = (1.0, 1.0, 1.0, 1.0);
        let result = lerp_color(low, high, 0.0);
        assert_eq!(result, (0.0, 0.0, 0.0, 0.0));
    }

    #[test]
    fn test_lerp_color_one() {
        let low = (0.0, 0.0, 0.0, 0.0);
        let high = (1.0, 1.0, 1.0, 1.0);
        let result = lerp_color(low, high, 1.0);
        assert_eq!(result, (1.0, 1.0, 1.0, 1.0));
    }

    #[test]
    fn test_lerp_color_mid() {
        let low = (0.0, 0.0, 0.0, 0.0);
        let high = (1.0, 1.0, 1.0, 1.0);
        let result = lerp_color(low, high, 0.5);
        assert_eq!(result, (0.5, 0.5, 0.5, 0.5));
    }

    #[test]
    fn test_lerp_color_clamps() {
        let low = (0.0, 0.0, 0.0, 0.0);
        let high = (1.0, 1.0, 1.0, 1.0);
        assert_eq!(lerp_color(low, high, -1.0), (0.0, 0.0, 0.0, 0.0));
        assert_eq!(lerp_color(low, high, 2.0), (1.0, 1.0, 1.0, 1.0));
    }
}
