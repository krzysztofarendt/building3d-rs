use crate::HasName;
use crate::Point;
use crate::geom::triangles::TriangleIndex;
use crate::{HasMesh, Mesh};
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
    let mesh: Mesh = model.get_mesh();
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
    let mesh: Mesh = model.get_mesh();
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
    let mesh: Mesh = model.get_mesh();
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
        &rr::Points3D::new(vertices).with_radii(radii).with_colors(colors),
    )?;

    Ok(())
}

