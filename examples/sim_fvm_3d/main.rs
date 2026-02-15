//! 3D FVM heat transfer example: aluminum block with volumetric heat source.
//!
//! Demonstrates the full 3D pipeline:
//!   Solid -> TetrahedralMesh -> FvmMesh -> FvmSparseSolver
//!
//! A 1×1×1 m aluminum block is uniformly cooled on all surfaces (h=10, T_fluid=20°C)
//! while cells near the center receive a volumetric heat source.  The simulation
//! watches the hot spot develop and equilibrate over 10 minutes.
//!
//! Run with a Rerun viewer open: `cargo run --example sim_fvm_3d`

use anyhow::Result;
use building3d::draw::rerun::{draw_edges, start_session};
use building3d::geom::bboxes::bounding_box;
use building3d::geom::mesh::tetrahedralize::tetrahedralize_delaunay_refined;
use building3d::geom::solid::containment::is_point_inside_solid;
use building3d::geom::tetrahedron::tetrahedron_centroid;
use building3d::sim::heat_transfer::boundary::BoundaryCondition;
use building3d::sim::heat_transfer::mesh::BOUNDARY;
use building3d::sim::heat_transfer::mesh_3d::build_3d_mesh_uniform;
use building3d::sim::heat_transfer::solver_sparse::FvmSparseSolver;
use building3d::{HasMesh, HasName, Mesh, Point, RerunConfig, Solid};
use rerun as rr;

// ---------------------------------------------------------------------------
// Helper: NamedMesh wrapper so raw Mesh can be passed to draw_edges
// ---------------------------------------------------------------------------
struct NamedMesh {
    mesh: Mesh,
    name: String,
}
impl HasMesh for NamedMesh {
    fn copy_mesh(&self) -> Mesh {
        self.mesh.clone()
    }
}
impl HasName for NamedMesh {
    fn get_name(&self) -> &str {
        &self.name
    }
}

// ---------------------------------------------------------------------------
// Helpers copied from draw_delaunay_refined.rs
// ---------------------------------------------------------------------------

fn dist(a: Point, b: Point) -> f64 {
    let dx = b.x - a.x;
    let dy = b.y - a.y;
    let dz = b.z - a.z;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

fn generate_surface_points(mesh: &Mesh, spacing: f64) -> Vec<Point> {
    let verts = mesh.vertices();
    let faces = match mesh.faces() {
        Some(f) => f,
        None => return Vec::new(),
    };
    let mut points = Vec::new();
    for face in faces {
        let a = verts[face.0];
        let b = verts[face.1];
        let c = verts[face.2];
        let max_edge = dist(a, b).max(dist(b, c)).max(dist(c, a));
        let n = (max_edge / spacing).ceil() as usize;
        if n < 2 {
            continue;
        }
        for i in 0..=n {
            for j in 0..=(n - i) {
                let k = n - i - j;
                if (i == n && j == 0 && k == 0)
                    || (i == 0 && j == n && k == 0)
                    || (i == 0 && j == 0 && k == n)
                {
                    continue;
                }
                let u = i as f64 / n as f64;
                let v = j as f64 / n as f64;
                let w = k as f64 / n as f64;
                points.push(Point::new(
                    u * a.x + v * b.x + w * c.x,
                    u * a.y + v * b.y + w * c.y,
                    u * a.z + v * b.z + w * c.z,
                ));
            }
        }
    }
    points
}

fn deduplicate_points(points: Vec<Point>, eps: f64) -> Vec<Point> {
    let eps_sq = eps * eps;
    let mut unique = Vec::with_capacity(points.len());
    for p in points {
        let is_dup = unique.iter().any(|q: &Point| {
            let dx = p.x - q.x;
            let dy = p.y - q.y;
            let dz = p.z - q.z;
            dx * dx + dy * dy + dz * dz < eps_sq
        });
        if !is_dup {
            unique.push(p);
        }
    }
    unique
}

fn generate_interior_points(solid: &Solid, spacing: f64) -> Vec<Point> {
    let mesh = solid.copy_mesh();
    let (bmin, bmax) = bounding_box(&mesh.vertices);
    let mut points = Vec::new();
    let mut x = bmin.x + spacing;
    while x < bmax.x {
        let mut y = bmin.y + spacing;
        while y < bmax.y {
            let mut z = bmin.z + spacing;
            while z < bmax.z {
                let p = Point::new(x, y, z);
                if is_point_inside_solid(solid, p) {
                    points.push(p);
                }
                z += spacing;
            }
            y += spacing;
        }
        x += spacing;
    }
    points
}

// ---------------------------------------------------------------------------
// Color map: blue (cold) -> white -> red (hot)
// ---------------------------------------------------------------------------
fn heat_color(t: f32) -> [u8; 4] {
    let (r, g, b) = if t < 0.5 {
        let s = t * 2.0;
        (s, s, 1.0)
    } else {
        let s = (t - 0.5) * 2.0;
        (1.0, 1.0 - s, 1.0 - s)
    };
    [(r * 255.0) as u8, (g * 255.0) as u8, (b * 255.0) as u8, 255]
}

fn main() -> Result<()> {
    // --- 1. Geometry -> Tet Mesh ---
    let spacing = 0.2;
    let solid = Solid::from_box(1.0, 1.0, 1.0, Some((0.0, 0.0, 0.0)), "block")?;
    let mesh = solid.copy_mesh();

    let surface_pts = deduplicate_points(generate_surface_points(&mesh, spacing), spacing * 0.1);
    let interior_pts = generate_interior_points(&solid, spacing);
    let mut extra = surface_pts;
    extra.extend(interior_pts.iter());

    let tet_mesh =
        tetrahedralize_delaunay_refined(&mesh, &extra).expect("tetrahedralization failed");
    println!(
        "Tet mesh: {} vertices, {} tetrahedra, volume = {:.4}",
        tet_mesh.vertices.len(),
        tet_mesh.tetrahedra_count(),
        tet_mesh.volume()
    );

    // --- 2. Tet Mesh -> FVM Mesh (aluminum: k=200, rho=2700, cp=900) ---
    let fvm_mesh = build_3d_mesh_uniform(&tet_mesh, 200.0, 2700.0, 900.0);
    let num_cells = fvm_mesh.cells.len();
    let num_faces = fvm_mesh.faces.len();
    println!("FVM mesh: {} cells, {} faces", num_cells, num_faces);

    // --- 3. Compute cell centroids and identify source cells ---
    let box_center = Point::new(0.5, 0.5, 0.5);
    let source_radius = 0.15;
    let total_source_w = 5000.0;

    let centroids: Vec<Point> = tet_mesh
        .tetrahedra
        .iter()
        .map(|t| {
            tetrahedron_centroid(
                tet_mesh.vertices[t.0],
                tet_mesh.vertices[t.1],
                tet_mesh.vertices[t.2],
                tet_mesh.vertices[t.3],
            )
        })
        .collect();

    let source_cells: Vec<usize> = centroids
        .iter()
        .enumerate()
        .filter(|(_, c)| dist(**c, box_center) < source_radius)
        .map(|(i, _)| i)
        .collect();

    let per_cell_source = if source_cells.is_empty() {
        0.0
    } else {
        total_source_w / source_cells.len() as f64
    };
    println!(
        "Source: {} cells within {:.2}m of center, {:.1} W each ({:.0} W total)",
        source_cells.len(),
        source_radius,
        per_cell_source,
        per_cell_source * source_cells.len() as f64
    );

    let mut sources = vec![0.0; num_cells];
    for &ci in &source_cells {
        sources[ci] = per_cell_source;
    }

    // --- 4. Boundary conditions: convective on all boundary faces ---
    let bc = BoundaryCondition::Convective {
        h: 10.0,
        t_fluid: 20.0,
    };
    let boundary_conditions: Vec<(usize, BoundaryCondition)> = fvm_mesh
        .faces
        .iter()
        .enumerate()
        .filter(|(_, f)| f.cell_left == BOUNDARY || f.cell_right == BOUNDARY)
        .map(|(idx, _)| (idx, bc))
        .collect();
    println!("Boundary faces: {}", boundary_conditions.len());

    // --- 5. Solver ---
    let mut solver = FvmSparseSolver::new(fvm_mesh, 20.0);

    // --- 6. Rerun visualization setup ---
    let config = RerunConfig::new();
    let session = start_session(&config)?;

    // Draw box wireframe (static)
    let box_named = NamedMesh {
        mesh: solid.copy_mesh(),
        name: "block_wireframe".to_string(),
    };
    draw_edges(&session, &box_named, 0.005, (0.5, 0.5, 0.5, 1.0), &config)?;

    // Configure time-series appearance
    session.log_static(
        "temperature/max",
        &rr::SeriesLines::new()
            .with_colors([[255, 50, 50]])
            .with_names(["T_max (°C)"])
            .with_widths([2.0]),
    )?;
    session.log_static(
        "temperature/mean",
        &rr::SeriesLines::new()
            .with_colors([[50, 200, 50]])
            .with_names(["T_mean (°C)"])
            .with_widths([2.0]),
    )?;
    session.log_static(
        "temperature/min",
        &rr::SeriesLines::new()
            .with_colors([[50, 50, 255]])
            .with_names(["T_min (°C)"])
            .with_widths([2.0]),
    )?;

    // Precompute cross-section filter: tets whose centroid.y is in [0.4, 0.6]
    let slice_cells: Vec<usize> = centroids
        .iter()
        .enumerate()
        .filter(|(_, c)| c.y >= 0.4 && c.y <= 0.6)
        .map(|(i, _)| i)
        .collect();
    println!("Cross-section cells: {}", slice_cells.len());

    // --- 7. Time loop ---
    let dt = 1.0;
    let num_steps = 600;
    let log_every = 10;

    println!(
        "Running {} steps (dt={} s, total={} s)...",
        num_steps,
        dt,
        num_steps as f64 * dt
    );

    for step in 0..=num_steps {
        if step > 0 {
            solver.step(dt, &boundary_conditions, &sources);
        }

        let temps = solver.temperatures();

        if step % log_every == 0 {
            // Temperature statistics
            let t_min = temps.iter().cloned().fold(f64::INFINITY, f64::min);
            let t_max = temps.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            let t_mean = temps.iter().sum::<f64>() / temps.len() as f64;

            session.set_time_sequence("step", step as i64);

            // Log scalars
            session.log("temperature/max", &rr::Scalars::single(t_max))?;
            session.log("temperature/mean", &rr::Scalars::single(t_mean))?;
            session.log("temperature/min", &rr::Scalars::single(t_min))?;

            // Build cross-section Mesh3D with per-vertex colors
            let mut vertices: Vec<rr::Vec3D> = Vec::new();
            let mut triangles: Vec<rr::TriangleIndices> = Vec::new();
            let mut colors: Vec<rr::Color> = Vec::new();

            for &ci in &slice_cells {
                let tet = &tet_mesh.tetrahedra[ci];
                let t_norm = if t_max > t_min {
                    ((temps[ci] - t_min) / (t_max - t_min)).clamp(0.0, 1.0) as f32
                } else {
                    0.5
                };
                let col = heat_color(t_norm);
                let rgba = rr::Color::from_unmultiplied_rgba(col[0], col[1], col[2], col[3]);

                // 4 vertices per tet (duplicated so each tet gets uniform color)
                let base = vertices.len() as u32;
                vertices.push(tet_mesh.vertices[tet.0].into());
                vertices.push(tet_mesh.vertices[tet.1].into());
                vertices.push(tet_mesh.vertices[tet.2].into());
                vertices.push(tet_mesh.vertices[tet.3].into());

                for _ in 0..4 {
                    colors.push(rgba);
                }

                // 4 triangular faces of the tetrahedron
                let tri =
                    |a: u32, b: u32, c: u32| rr::TriangleIndices(rr::datatypes::UVec3D([a, b, c]));
                triangles.push(tri(base, base + 1, base + 2));
                triangles.push(tri(base, base + 1, base + 3));
                triangles.push(tri(base, base + 2, base + 3));
                triangles.push(tri(base + 1, base + 2, base + 3));
            }

            session.log(
                "Building3d/cross_section",
                &rr::Mesh3D::new(vertices)
                    .with_triangle_indices(triangles)
                    .with_vertex_colors(colors),
            )?;

            if step % 100 == 0 {
                println!(
                    "  step {:>4}: T_min={:.2}, T_mean={:.2}, T_max={:.2}",
                    step, t_min, t_mean, t_max
                );
            }
        }
    }

    println!("Done.");
    Ok(())
}
