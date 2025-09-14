// TODO: Triangulation first

// use three_d::{Window, WindowSettings};
// use three_d::Camera;
// use three_d::vec3;
// use three_d::degrees;
// use three_d::control::OrbitControl;
// use three_d::CpuMesh;
// use three_d::Positions;
// use three_d::Indices;
// use three_d::Mesh;

// use crate::geom::polygon::Polygon;
// use anyhow::Result;


// pub fn draw_polygon(poly: Polygon) -> Result<()> {

//     let vertices = poly.pts;

//     // Window & GL
//     let window = Window::new(WindowSettings {
//             title: "Polygon".into(),
//             ..Default::default()
//     })?;
//     let context = window.gl();

//     // Camera & control
//     // TODO: Camera position should automatically adapt to the size of mesh
//     let mut camera = Camera::new_perspective(
//         window.viewport(),
//         vec3(0.0, 1.0, 8.0),
//         vec3(0.0, 0.0, 0.0),
//         vec3(0.0, 1.0, 0.0),
//         degrees(45.0),
//         0.1,
//         1000.0,
//     );
//     let mut control = OrbitControl::new(camera.target(), 0.5, 50.0);

//     // Build mesh
//     let mut cpu = CpuMesh {
//         positions: Positions::F32(vertices),
//         indices: Indices::U32(triangles),
//         ..Default::default()
//     };
//     cpu.compute_normals();
//     let mesh = Mesh::new(&context, &cpu);

//     Ok(())
// }


