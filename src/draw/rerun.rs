use crate::Building;
use crate::Point;
use crate::geom::triangles::TriangleIndex;
use anyhow::Result;
use rerun as rr;

impl From<Point> for rr::Position3D {
    fn from(val: Point) -> Self {
        rr::Position3D(rr::Vec3D([val.x as f32, val.y as f32, val.z as f32]))
    }
}

impl From<TriangleIndex> for rr::TriangleIndices {
    fn from(val: TriangleIndex) -> Self {
        rr::TriangleIndices(rr::datatypes::UVec3D([
            val.0 as u32,
            val.1 as u32,
            val.2 as u32,
        ]))
    }
}

pub fn show(building: &Building) -> Result<()> {
    let mesh = building.mesh();
    let vertices: Vec<Point> = mesh.vertices;
    let triangles: Vec<TriangleIndex> = mesh.triangles;

    let alpha: f32 = 0.2;
    let red: f32 = 1.;
    let green: f32 = 1.;
    let blue: f32 = 1.;

    // Connect to the Rerun gRPC server using the default address and port: localhost:9876
    let rec = rr::RecordingStreamBuilder::new("building3d").spawn()?;
    rec.log_static(
        "model",
        &rr::Mesh3D::new(vertices)
            .with_triangle_indices(triangles)
            .with_albedo_factor(rr::Rgba32::from_linear_unmultiplied_rgba_f32(
                red, green, blue, alpha,
            )),
    )?;

    Ok(())
}
