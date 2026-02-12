use anyhow::Result;
use building3d::draw::rerun::{draw_edges, draw_faces, start_session};
use building3d::geom::bboxes::bounding_box;
use building3d::geom::mesh::tetrahedralize::tetrahedralize_delaunay_refined;
use building3d::geom::solid::containment::is_point_inside_solid;
use building3d::{FloorPlan, HasMesh, HasName, Mesh, Point, RerunConfig, Solid, Vector};

/// Lightweight wrapper so raw `Mesh` values can be passed to `draw_*`.
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

/// Distance between two points.
fn dist(a: Point, b: Point) -> f64 {
    let dx = b.x - a.x;
    let dy = b.y - a.y;
    let dz = b.z - a.z;
    (dx * dx + dy * dy + dz * dz).sqrt()
}

/// Generate points on the surface of a mesh by subdividing each triangle face.
///
/// The number of subdivisions per triangle is derived from its longest edge
/// divided by `spacing`, so that surface element sizes match the interior grid.
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

        // Compute subdivisions from longest edge
        let max_edge = dist(a, b).max(dist(b, c)).max(dist(c, a));
        let n = (max_edge / spacing).ceil() as usize;
        if n < 2 {
            continue; // triangle is already smaller than spacing
        }

        // Sample points using barycentric coords, skipping the 3 corner vertices.
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

/// Remove near-duplicate points (within epsilon tolerance).
fn deduplicate_points(points: Vec<Point>, eps: f64) -> Vec<Point> {
    let eps_sq = eps * eps;
    let mut unique = Vec::with_capacity(points.len());
    for p in points {
        let is_dup = unique
            .iter()
            .any(|q: &Point| {
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

/// Generate interior points inside a solid on a regular 3D grid.
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

fn main() -> Result<()> {
    let spacing = 0.25;

    // --- Convex box ---
    let box_solid = Solid::from_box(2.0, 2.0, 2.0, Some((0.0, 0.0, 0.0)), "box")?;
    let box_mesh = box_solid.copy_mesh();
    let box_surface = deduplicate_points(generate_surface_points(&box_mesh, spacing), spacing * 0.1);
    let box_interior = generate_interior_points(&box_solid, spacing);
    let mut box_extra = box_surface;
    box_extra.extend(box_interior.iter());
    let box_tets = tetrahedralize_delaunay_refined(&box_mesh, &box_extra)
        .expect("box tetrahedralization failed");
    println!(
        "Box: {} surface + {} interior = {} extra pts, {} tetrahedra, volume = {:.3}",
        box_extra.len() - box_interior.len(),
        box_interior.len(),
        box_extra.len(),
        box_tets.tetrahedra_count(),
        box_tets.volume()
    );

    // --- Concave L-shape ---
    let fp = FloorPlan {
        plan: vec![
            (0.0, 0.0),
            (3.0, 0.0),
            (3.0, 1.0),
            (1.0, 1.0),
            (1.0, 2.0),
            (0.0, 2.0),
        ],
        height: 1.0,
        name: "l_shape".to_string(),
        ..Default::default()
    };
    let mut l_solid = Solid::from_floor_plan(fp)?;
    l_solid.translate(&Vector::new(4.0, 0.0, 0.0));
    let l_mesh = l_solid.copy_mesh();
    let l_surface = deduplicate_points(generate_surface_points(&l_mesh, spacing), spacing * 0.1);
    let l_interior = generate_interior_points(&l_solid, spacing);
    let mut l_extra = l_surface;
    l_extra.extend(l_interior.iter());
    let l_tets = tetrahedralize_delaunay_refined(&l_mesh, &l_extra)
        .expect("L-shape tetrahedralization failed");
    println!(
        "L-shape: {} surface + {} interior = {} extra pts, {} tetrahedra, volume = {:.3}",
        l_extra.len() - l_interior.len(),
        l_interior.len(),
        l_extra.len(),
        l_tets.tetrahedra_count(),
        l_tets.volume()
    );

    // Start Rerun session
    let config = RerunConfig::new();
    let session = start_session(&config)?;

    // Draw original surfaces (semi-transparent)
    let surface_color = (1.0, 1.0, 1.0, 0.15);
    draw_faces(&session, &box_solid, surface_color, &config)?;
    draw_faces(&session, &l_solid, surface_color, &config)?;

    // Draw tet edges
    let box_tet_surface = NamedMesh {
        mesh: box_tets.to_surface_mesh(),
        name: "box_delaunay_refined".to_string(),
    };
    let l_tet_surface = NamedMesh {
        mesh: l_tets.to_surface_mesh(),
        name: "l_shape_delaunay_refined".to_string(),
    };

    let edge_radius = 0.005;
    let edge_color = (0.2, 0.8, 1.0, 0.6);
    draw_edges(&session, &box_tet_surface, edge_radius, edge_color, &config)?;
    draw_edges(&session, &l_tet_surface, edge_radius, edge_color, &config)?;

    Ok(())
}
