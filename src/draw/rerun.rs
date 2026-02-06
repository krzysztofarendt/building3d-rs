use crate::HasName;
use crate::Point;
use crate::geom::triangles::TriangleIndex;
use crate::sim::rays::SimulationResult;
use crate::{Building, HasMesh, Mesh};
use anyhow::Result;
use rerun as rr;

const SESSION_NAME: &str = "Building3d";

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

fn color(rgba: (f32, f32, f32, f32)) -> rr::Color {
    let (r, g, b, a) = rgba;
    rr::Color(rr::Rgba32::from_linear_unmultiplied_rgba_f32(r, g, b, a))
}

pub fn start_session() -> Result<rr::RecordingStream> {
    // Connect to the Rerun gRPC server using the default address and port: localhost:9876
    let session = rr::RecordingStreamBuilder::new("building3d").spawn()?;

    Ok(session)
}

pub fn draw_faces<T: HasMesh + HasName>(
    session: &rr::RecordingStream,
    model: &T,
    rgba: (f32, f32, f32, f32),
) -> Result<()> {
    let mesh: Mesh = model.copy_mesh();
    let vertices: Vec<Point> = mesh.vertices;
    let triangles: Vec<TriangleIndex> = mesh.faces.unwrap_or_default();

    let (r, g, b, a) = rgba;

    let name = format!("{}/{}", SESSION_NAME, model.get_name());
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
    rgba: (f32, f32, f32, f32),
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

    let name = format!("{}/{}", SESSION_NAME, model.get_name());
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
    rgba: (f32, f32, f32, f32),
) -> Result<()> {
    let mesh: Mesh = model.copy_mesh();
    let vertices: Vec<Point> = mesh.vertices;

    let mut radii: Vec<f32> = Vec::new();
    let mut colors: Vec<rr::Color> = Vec::new();

    for _ in 0..vertices.len() {
        radii.push(radius);
        colors.push(color(rgba));
    }

    let name = format!("{}/{}", SESSION_NAME, model.get_name());
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
) -> Result<()> {
    // Draw building mesh (semi-transparent)
    draw_faces(session, building, (0.8, 0.8, 0.8, 0.15))?;

    // Draw absorbers as static red spheres
    for (i, absorber) in result.config.absorbers.iter().enumerate() {
        let name = format!("{}/absorbers/{}", SESSION_NAME, i);
        session.log_static(
            name,
            &rr::Points3D::new([*absorber])
                .with_radii([result.config.absorber_radius as f32])
                .with_colors([color((1.0, 0.0, 0.0, 0.5))]),
        )?;
    }

    // Draw source as a static green sphere
    session.log_static(
        format!("{}/source", SESSION_NAME),
        &rr::Points3D::new([result.config.source])
            .with_radii([0.02f32])
            .with_colors([color((0.0, 1.0, 0.0, 1.0))]),
    )?;

    // Animate ray positions over time steps
    for (step, (positions, energies)) in result
        .positions
        .iter()
        .zip(result.energies.iter())
        .enumerate()
    {
        session.set_time_sequence("step", step as i64);

        // Collect alive rays (energy > 0)
        let mut pts: Vec<Point> = Vec::new();
        let mut colors_vec: Vec<rr::Color> = Vec::new();
        let mut radii_vec: Vec<f32> = Vec::new();

        for (pos, &energy) in positions.iter().zip(energies.iter()) {
            if energy > 1e-10 {
                pts.push(*pos);
                // Color from bright red (high energy) to dark red (low energy)
                let e = energy.clamp(0.0, 1.0) as f32;
                colors_vec.push(color((e, 0.0, 0.0, 0.8)));
                radii_vec.push(0.04);
            }
        }

        if !pts.is_empty() {
            session.log(
                format!("{}/rays", SESSION_NAME),
                &rr::Points3D::new(pts)
                    .with_radii(radii_vec)
                    .with_colors(colors_vec),
            )?;
        }
    }

    Ok(())
}
