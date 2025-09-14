use anyhow::Result;
use three_d::Blend;
use three_d::Camera;
use three_d::ColorMaterial;
use three_d::CpuMaterial;
use three_d::CpuMesh;
use three_d::Cull;
use three_d::Gm;
use three_d::Indices;
use three_d::Mesh;
use three_d::PhysicalMaterial;
use three_d::Positions;
use three_d::RenderStates;
use three_d::Srgba;
use three_d::WriteMask;
use three_d::control::OrbitControl;
use three_d::degrees;
use three_d::vec3;
use three_d::{
    ClearState, DepthTest, FrameOutput, InnerSpace, InstancedMesh, Instances, Mat4, Quat, Vec3,
};
use three_d::{Window, WindowSettings};

use crate::geom::point::Point;
use crate::geom::polygon::Polygon;
use crate::geom::triangles::TriangleIndex;

const MAX_DISTANCE: f32 = 1000.0;

fn points_to_positions(pts: &[Point]) -> Positions {
    Positions::F64(pts.iter().map(|p| vec3(p.x, p.y, p.z)).collect())
}

fn triangles_to_indices(tri: &[TriangleIndex]) -> Indices {
    Indices::U32(
        tri.iter()
            .flat_map(|x| [x.0 as u32, x.1 as u32, x.2 as u32])
            .collect(),
    )
}

pub fn draw_polygon(poly: &Polygon) -> Result<()> {
    // Window & GL
    let window = Window::new(WindowSettings {
        title: "Polygon".into(),
        ..Default::default()
    })?;
    let context = window.gl();

    // Build mesh
    let mut cpu = CpuMesh {
        positions: points_to_positions(&poly.pts),
        indices: triangles_to_indices(&poly.tri),
        ..Default::default()
    };
    cpu.compute_normals();
    let mesh = Mesh::new(&context, &cpu);

    // Translucent fill
    let fill = Gm::new(
        mesh,
        ColorMaterial {
            color: Srgba::new(0, 90, 255, 160),
            render_states: RenderStates {
                depth_test: DepthTest::Always,
                write_mask: WriteMask::COLOR,
                blend: Blend::TRANSPARENCY,
                ..Default::default()
            },
            is_transparent: true,
            ..Default::default()
        },
    );

    // Wireframe (cylinders)
    let mut wire_mat = PhysicalMaterial::new_opaque(
        &context,
        &CpuMaterial {
            albedo: Srgba::new_opaque(0, 0, 0),
            roughness: 0.7,
            metallic: 0.0,
            ..Default::default()
        },
    );
    wire_mat.render_states.cull = Cull::Back;

    // Compute positions as Vec3<f32> for transformations
    let positions = poly
        .pts
        .iter()
        .map(|p| vec3(p.x as f32, p.y as f32, p.z as f32))
        .collect::<Vec<Vec3>>();

    // Compute center and radius for camera framing and scaling
    let center = positions.iter().copied().reduce(|a, b| a + b).unwrap() / positions.len() as f32;
    let radius = positions
        .iter()
        .map(|p| (p - center).magnitude())
        .fold(0.0, f32::max);

    // Camera & control adapt to mesh size
    let mut camera = Camera::new_perspective(
        window.viewport(),
        center + vec3(1.0, 1.0, 1.0).normalize() * (radius * 2.0),
        center,
        vec3(0.0, 1.0, 0.0),
        degrees(45.0),
        0.1,
        radius * MAX_DISTANCE,
    );
    let mut control = OrbitControl::new(center, radius * 0.5, radius * MAX_DISTANCE);

    // Prepare instanced meshes for vertices and edges
    let mut sphere_cpu = CpuMesh::sphere(16);
    sphere_cpu.transform(Mat4::from_scale(radius * 0.03))?;
    let vertex_instances = Instances {
        transformations: positions
            .iter()
            .copied()
            .map(Mat4::from_translation)
            .collect(),
        ..Default::default()
    };
    let vertex_gm = Gm::new(
        InstancedMesh::new(&context, &vertex_instances, &sphere_cpu),
        ColorMaterial {
            color: Srgba::new_opaque(255, 0, 0),
            render_states: RenderStates {
                depth_test: DepthTest::Always,
                write_mask: WriteMask::COLOR,
                ..Default::default()
            },
            ..Default::default()
        },
    );

    let base_cylinder = CpuMesh::cylinder(12);

    // Outer polygon boundary edges (gray, slightly thicker)
    let mut cyl_bound = base_cylinder.clone();
    cyl_bound.transform(Mat4::from_nonuniform_scale(
        1.0,
        radius * 0.005,
        radius * 0.005,
    ))?;

    // Triangle edges (light gray, thinner)
    let mut cyl_tri = base_cylinder;
    cyl_tri.transform(Mat4::from_nonuniform_scale(
        1.0,
        radius * 0.002,
        radius * 0.002,
    ))?;

    // Function to compute a cylinder transform between two points
    fn edge_transform(p1: Vec3, p2: Vec3) -> Mat4 {
        Mat4::from_translation(p1)
            * Into::<Mat4>::into(Quat::from_arc(
                vec3(1.0, 0.0, 0.0),
                (p2 - p1).normalize(),
                None,
            ))
            * Mat4::from_nonuniform_scale((p2 - p1).magnitude(), 1.0, 1.0)
    }

    // Triangle edge instances
    let tri_transforms = poly
        .tri
        .iter()
        .flat_map(|t| {
            let [i1, i2, i3] = [t.0, t.1, t.2];
            [
                edge_transform(positions[i1], positions[i2]),
                edge_transform(positions[i2], positions[i3]),
                edge_transform(positions[i3], positions[i1]),
            ]
        })
        .collect::<Vec<Mat4>>();
    let tri_instances = Instances {
        transformations: tri_transforms,
        ..Default::default()
    };
    let tri_gm = Gm::new(
        InstancedMesh::new(&context, &tri_instances, &cyl_tri),
        ColorMaterial {
            color: Srgba::new_opaque(200, 200, 200),
            render_states: RenderStates {
                depth_test: DepthTest::Always,
                write_mask: WriteMask::COLOR,
                ..Default::default()
            },
            ..Default::default()
        },
    );

    // Boundary edge instances
    let mut edge_counts = std::collections::HashMap::new();
    for t in &poly.tri {
        let [i1, i2, i3] = [t.0, t.1, t.2];
        for &(a, b) in &[(i1, i2), (i2, i3), (i3, i1)] {
            let key = if a < b { (a, b) } else { (b, a) };
            *edge_counts.entry(key).or_insert(0) += 1;
        }
    }
    let bound_transforms = edge_counts
        .into_iter()
        .filter_map(|((a, b), count)| {
            if count == 1 {
                Some(edge_transform(positions[a], positions[b]))
            } else {
                None
            }
        })
        .collect::<Vec<Mat4>>();
    let bound_instances = Instances {
        transformations: bound_transforms,
        ..Default::default()
    };
    let bound_gm = Gm::new(
        InstancedMesh::new(&context, &bound_instances, &cyl_bound),
        ColorMaterial {
            color: Srgba::new_opaque(100, 100, 100),
            render_states: RenderStates {
                depth_test: DepthTest::Always,
                write_mask: WriteMask::COLOR,
                ..Default::default()
            },
            ..Default::default()
        },
    );

    // Render loop
    window.render_loop(move |mut frame_input| {
        camera.set_viewport(frame_input.viewport);
        control.handle_events(&mut camera, &mut frame_input.events);

        frame_input
            .screen()
            .clear(ClearState::color_and_depth(0.9, 0.9, 0.9, 1.0, 1.0))
            .render(
                &camera,
                fill.into_iter()
                    .chain(&tri_gm)
                    .chain(&bound_gm)
                    .chain(&vertex_gm),
                &[],
            );

        FrameOutput::default()
    });
    Ok(())
}
