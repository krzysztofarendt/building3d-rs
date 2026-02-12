use anyhow::Result;
use building3d::draw::rerun::{draw_edges, draw_faces, start_session};
use building3d::geom::mesh::tetrahedralize::tetrahedralize_delaunay;
use building3d::{FloorPlan, HasMesh, HasName, Mesh, RerunConfig, Solid};

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

fn main() -> Result<()> {
    // --- Convex box ---
    let box_solid = Solid::from_box(2.0, 2.0, 2.0, Some((0.0, 0.0, 0.0)), "box")?;

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
    let l_solid = Solid::from_floor_plan(fp)?;

    // Tetrahedralize both
    let box_mesh = box_solid.copy_mesh();
    let box_tets = tetrahedralize_delaunay(&box_mesh).expect("box tetrahedralization failed");
    println!(
        "Box: {} tetrahedra, volume = {:.3}",
        box_tets.tetrahedra_count(),
        box_tets.volume()
    );

    let l_mesh = l_solid.copy_mesh();
    let l_tets = tetrahedralize_delaunay(&l_mesh).expect("L-shape tetrahedralization failed");
    println!(
        "L-shape: {} tetrahedra, volume = {:.3}",
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

    // Draw tet edges by converting TetrahedralMesh → surface Mesh (all tet faces)
    let box_tet_surface = NamedMesh {
        mesh: box_tets.to_surface_mesh(),
        name: "box_delaunay".to_string(),
    };
    let l_tet_surface = NamedMesh {
        mesh: l_tets.to_surface_mesh(),
        name: "l_shape_delaunay".to_string(),
    };

    let edge_radius = 0.005;
    let edge_color = (0.2, 0.8, 1.0, 0.6);
    draw_edges(&session, &box_tet_surface, edge_radius, edge_color, &config)?;
    draw_edges(&session, &l_tet_surface, edge_radius, edge_color, &config)?;

    Ok(())
}
