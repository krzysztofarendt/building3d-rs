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
const VERTEX_SPHERE_SIZE: f32 = 1.0; // percent of bounding radius for vertex sphere scale (and mesh resolution)

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

pub fn draw_polygons(polygons: &[Polygon]) -> Result<()> {
    // Window & GL
    let window = Window::new(WindowSettings {
        title: "Polygon".into(),
        ..Default::default()
    })?;
    let context = window.gl();

    // Translucent fill for each polygon
    let mut fill_gms = Vec::new();
    for poly in polygons {
        let mut cpu = CpuMesh {
            positions: points_to_positions(&poly.pts),
            indices: triangles_to_indices(&poly.tri),
            ..Default::default()
        };
        cpu.compute_normals();
        let mesh = Mesh::new(&context, &cpu);
        let mat = ColorMaterial {
            color: Srgba::new(0, 90, 255, 160),
            render_states: RenderStates {
                depth_test: DepthTest::Always,
                write_mask: WriteMask::COLOR,
                blend: Blend::TRANSPARENCY,
                ..Default::default()
            },
            is_transparent: true,
            ..Default::default()
        };
        fill_gms.push(Gm::new(mesh, mat));
    }

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

    // Gather all vertices, triangle-edges, and boundary-edges across polygons
    let mut flat_positions = Vec::new();
    let mut vertex_transforms = Vec::new();
    let mut tri_transforms = Vec::new();
    let mut bound_transforms = Vec::new();
    for poly in polygons {
        let positions = poly
            .pts
            .iter()
            .map(|p| vec3(p.x as f32, p.y as f32, p.z as f32))
            .collect::<Vec<Vec3>>();
        for &p in &positions {
            flat_positions.push(p);
            vertex_transforms.push(Mat4::from_translation(p));
        }
        let edge_transform = |p1: Vec3, p2: Vec3| {
            Mat4::from_translation(p1)
                * Into::<Mat4>::into(Quat::from_arc(
                    vec3(1.0, 0.0, 0.0),
                    (p2 - p1).normalize(),
                    None,
                ))
                * Mat4::from_nonuniform_scale((p2 - p1).magnitude(), 1.0, 1.0)
        };
        for t in &poly.tri {
            let [i1, i2, i3] = [t.0, t.1, t.2];
            tri_transforms.push(edge_transform(positions[i1], positions[i2]));
            tri_transforms.push(edge_transform(positions[i2], positions[i3]));
            tri_transforms.push(edge_transform(positions[i3], positions[i1]));
        }
        let mut counts = std::collections::HashMap::new();
        for t in &poly.tri {
            let [i1, i2, i3] = [t.0, t.1, t.2];
            for &(a, b) in &[(i1, i2), (i2, i3), (i3, i1)] {
                let key = if a < b { (a, b) } else { (b, a) };
                *counts.entry(key).or_insert(0) += 1;
            }
        }
        for ((a, b), count) in counts {
            if count == 1 {
                let p1 = positions[a];
                let p2 = positions[b];
                bound_transforms.push(edge_transform(p1, p2));
            }
        }
    }

    // Compute camera frame from all vertices
    let center =
        flat_positions.iter().copied().reduce(|a, b| a + b).unwrap() / flat_positions.len() as f32;
    let radius = flat_positions
        .iter()
        .map(|p| (p - center).magnitude())
        .fold(0.0, f32::max);
    // Camera & control adapt to mesh size, but start no further than MAX_DISTANCE
    let start_dist = (radius * 2.0).min(MAX_DISTANCE);
    let far_dist = radius * MAX_DISTANCE;
    let mut camera = Camera::new_perspective(
        window.viewport(),
        center + vec3(1.0, 1.0, 1.0).normalize() * start_dist,
        center,
        vec3(0.0, 1.0, 0.0),
        degrees(45.0),
        0.1,
        far_dist,
    );
    let mut control = OrbitControl::new(center, radius * 0.5, far_dist);

    // Build instanced meshes for vertices, triangle-edges, and boundary-edges
    let mut sphere_cpu = CpuMesh::sphere(16);
    // scale sphere radius as percentage of bounding radius
    let sphere_scale = radius * (VERTEX_SPHERE_SIZE / 100.0);
    sphere_cpu.transform(Mat4::from_scale(sphere_scale))?;
    let vertex_gm = Gm::new(
        InstancedMesh::new(
            &context,
            &Instances {
                transformations: vertex_transforms,
                ..Default::default()
            },
            &sphere_cpu,
        ),
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
    let mut cyl_bound = CpuMesh::cylinder(12);
    cyl_bound.transform(Mat4::from_nonuniform_scale(
        1.0,
        radius * 0.005,
        radius * 0.005,
    ))?;
    let bound_gm = Gm::new(
        InstancedMesh::new(
            &context,
            &Instances {
                transformations: bound_transforms,
                ..Default::default()
            },
            &cyl_bound,
        ),
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
    let mut cyl_tri = CpuMesh::cylinder(12);
    cyl_tri.transform(Mat4::from_nonuniform_scale(
        1.0,
        radius * 0.002,
        radius * 0.002,
    ))?;
    let tri_gm = Gm::new(
        InstancedMesh::new(
            &context,
            &Instances {
                transformations: tri_transforms,
                ..Default::default()
            },
            &cyl_tri,
        ),
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

    // Render loop
    window.render_loop(move |mut frame_input| {
        camera.set_viewport(frame_input.viewport);
        control.handle_events(&mut camera, &mut frame_input.events);

        let objects = fill_gms
            .iter()
            .flat_map(|g| g.into_iter())
            .chain(&tri_gm)
            .chain(&bound_gm)
            .chain(&vertex_gm);
        frame_input
            .screen()
            .clear(ClearState::color_and_depth(0.9, 0.9, 0.9, 1.0, 1.0))
            .render(&camera, objects, &[]);

        FrameOutput::default()
    });
    Ok(())
}
