use crate::Point;
use crate::geom::mesh::tetrahedralize::{TetrahedralMesh, TetrahedronIndex};
use crate::geom::tetrahedron::{tetrahedron_centroid, tetrahedron_volume};
use crate::sim::heat_transfer::mesh::{BOUNDARY, FvmCell, FvmFace, FvmMesh};
use std::collections::{HashMap, HashSet};

/// Thermal properties assigned to one tetrahedral FVM cell.
#[derive(Debug, Clone, Copy)]
pub struct CellThermalProperties {
    /// Thermal conductivity in W/(m*K).
    pub conductivity: f64,
    /// Density in kg/m^3.
    pub density: f64,
    /// Specific heat in J/(kg*K).
    pub specific_heat: f64,
}

impl CellThermalProperties {
    pub fn new(conductivity: f64, density: f64, specific_heat: f64) -> Self {
        Self {
            conductivity,
            density,
            specific_heat,
        }
    }
}

#[derive(Debug, Clone)]
struct PendingFace {
    cell_idx: usize,
    area: f64,
    distance_to_face: f64,
}

/// Build a 3D finite-volume mesh from a tetrahedral mesh using one uniform
/// material for all tetrahedra.
pub fn build_3d_mesh_uniform(
    tet_mesh: &TetrahedralMesh,
    conductivity: f64,
    density: f64,
    specific_heat: f64,
) -> FvmMesh {
    let props = vec![
        CellThermalProperties {
            conductivity,
            density,
            specific_heat,
        };
        tet_mesh.tetrahedra.len()
    ];
    build_3d_mesh(tet_mesh, &props)
}

/// Build a 3D finite-volume mesh from a tetrahedral mesh.
///
/// - One [`FvmCell`] is created per tetrahedron.
/// - Shared triangular faces become interior [`FvmFace`] entries.
/// - Unmatched triangular faces become boundary [`FvmFace`] entries.
///
/// Interior-face conductance uses a series-resistance form:
/// `K = A / (d_l / k_l + d_r / k_r)`,
/// where `d_l` and `d_r` are centroid-to-face distances for the left/right
/// tetrahedra, and `k_l`, `k_r` are their conductivities.
///
/// Boundary-face conductance is `K = k * A / d`, where `d` is the
/// centroid-to-face distance of the adjacent tetrahedron.
pub fn build_3d_mesh(
    tet_mesh: &TetrahedralMesh,
    cell_properties: &[CellThermalProperties],
) -> FvmMesh {
    assert_eq!(
        tet_mesh.tetrahedra.len(),
        cell_properties.len(),
        "cell properties length ({}) must match tetrahedra count ({})",
        cell_properties.len(),
        tet_mesh.tetrahedra.len()
    );

    let mut cells = Vec::with_capacity(tet_mesh.tetrahedra.len());
    let mut centroids = Vec::with_capacity(tet_mesh.tetrahedra.len());

    for (i, tet) in tet_mesh.tetrahedra.iter().enumerate() {
        let (p0, p1, p2, p3) = tetra_points(&tet_mesh.vertices, tet, i);
        let volume = tetrahedron_volume(p0, p1, p2, p3);
        assert!(volume > 0.0, "tetrahedron {i} has zero/negative volume");

        let props = cell_properties[i];
        cells.push(FvmCell {
            volume,
            conductivity: props.conductivity,
            density: props.density,
            specific_heat: props.specific_heat,
        });
        centroids.push(tetrahedron_centroid(p0, p1, p2, p3));
    }

    let mut faces: Vec<FvmFace> = Vec::new();
    let mut pending: HashMap<[usize; 3], PendingFace> = HashMap::new();
    let mut paired: HashSet<[usize; 3]> = HashSet::new();

    for (cell_idx, tet) in tet_mesh.tetrahedra.iter().enumerate() {
        let (p0, p1, p2, p3) = tetra_points(&tet_mesh.vertices, tet, cell_idx);
        let tri_faces = [
            ([tet.0, tet.1, tet.2], (p0, p1, p2)),
            ([tet.0, tet.1, tet.3], (p0, p1, p3)),
            ([tet.0, tet.2, tet.3], (p0, p2, p3)),
            ([tet.1, tet.2, tet.3], (p1, p2, p3)),
        ];

        for (idxs, (a, b, c)) in tri_faces {
            let key = sorted_face_key(idxs);
            assert!(
                !paired.contains(&key),
                "Non-manifold tetrahedral mesh: face {:?} is shared by more than 2 cells",
                key
            );

            let area = triangle_area(a, b, c);
            assert!(
                area > 0.0,
                "Degenerate triangle face in tetrahedron {cell_idx}"
            );

            let d_this = point_to_plane_distance(centroids[cell_idx], a, b, c);
            assert!(
                d_this > 0.0,
                "Invalid centroid-to-face distance in tetrahedron {cell_idx}"
            );

            if let Some(other) = pending.remove(&key) {
                let k_l = cells[other.cell_idx].conductivity;
                let k_r = cells[cell_idx].conductivity;
                let resistance = other.distance_to_face / k_l + d_this / k_r;
                let conductance = other.area / resistance;

                faces.push(FvmFace {
                    cell_left: other.cell_idx,
                    cell_right: cell_idx,
                    area: other.area,
                    distance: other.distance_to_face + d_this,
                    conductance,
                });
                paired.insert(key);
            } else {
                pending.insert(
                    key,
                    PendingFace {
                        cell_idx,
                        area,
                        distance_to_face: d_this,
                    },
                );
            }
        }
    }

    let mut boundary_keys: Vec<[usize; 3]> = pending.keys().copied().collect();
    boundary_keys.sort_unstable();
    for key in boundary_keys {
        let entry = pending
            .get(&key)
            .expect("boundary face key disappeared unexpectedly");
        let k = cells[entry.cell_idx].conductivity;
        faces.push(FvmFace {
            cell_left: entry.cell_idx,
            cell_right: BOUNDARY,
            area: entry.area,
            distance: entry.distance_to_face,
            conductance: k * entry.area / entry.distance_to_face,
        });
    }

    FvmMesh { cells, faces }
}

fn tetra_points(
    vertices: &[Point],
    tet: &TetrahedronIndex,
    tet_idx: usize,
) -> (Point, Point, Point, Point) {
    assert!(
        tet.0 < vertices.len()
            && tet.1 < vertices.len()
            && tet.2 < vertices.len()
            && tet.3 < vertices.len(),
        "tetrahedron {tet_idx} has an out-of-bounds vertex index"
    );
    (
        vertices[tet.0],
        vertices[tet.1],
        vertices[tet.2],
        vertices[tet.3],
    )
}

fn sorted_face_key(mut face: [usize; 3]) -> [usize; 3] {
    face.sort_unstable();
    face
}

fn triangle_area(a: Point, b: Point, c: Point) -> f64 {
    let ab = b - a;
    let ac = c - a;
    0.5 * ab.cross(&ac).length()
}

fn point_to_plane_distance(p: Point, a: Point, b: Point, c: Point) -> f64 {
    let ab = b - a;
    let ac = c - a;
    let n = ab.cross(&ac);
    let n_len = n.length();
    if n_len <= 1e-15 {
        return 0.0;
    }
    ((p - a).dot(&n)).abs() / n_len
}

#[cfg(test)]
mod tests {
    use super::*;

    fn base_two_tetra_mesh() -> TetrahedralMesh {
        // Two tetrahedra sharing face (0,1,2) in z=0 plane.
        let vertices = vec![
            Point::new(0.0, 0.0, 0.0),  // 0
            Point::new(1.0, 0.0, 0.0),  // 1
            Point::new(0.0, 1.0, 0.0),  // 2
            Point::new(0.0, 0.0, 1.0),  // 3
            Point::new(0.0, 0.0, -1.0), // 4
        ];
        let tetrahedra = vec![TetrahedronIndex(0, 1, 2, 3), TetrahedronIndex(0, 1, 2, 4)];
        TetrahedralMesh::new(vertices, tetrahedra)
    }

    #[test]
    fn test_single_tetra_builds_boundary_faces() {
        let tet_mesh = TetrahedralMesh::new(
            vec![
                Point::new(0.0, 0.0, 0.0),
                Point::new(1.0, 0.0, 0.0),
                Point::new(0.0, 1.0, 0.0),
                Point::new(0.0, 0.0, 1.0),
            ],
            vec![TetrahedronIndex(0, 1, 2, 3)],
        );

        let mesh = build_3d_mesh_uniform(&tet_mesh, 2.0, 1000.0, 1000.0);
        assert_eq!(mesh.cells.len(), 1);
        assert_eq!(mesh.faces.len(), 4);

        // Volume of this tetrahedron is 1/6.
        assert!((mesh.cells[0].volume - (1.0 / 6.0)).abs() < 1e-12);

        let boundary_faces = mesh
            .faces
            .iter()
            .filter(|f| f.cell_right == BOUNDARY || f.cell_left == BOUNDARY)
            .count();
        assert_eq!(boundary_faces, 4);

        for face in &mesh.faces {
            assert!(face.distance > 0.0);
            assert!(face.area > 0.0);
            assert!(face.conductance > 0.0);
        }
    }

    #[test]
    fn test_two_tetra_one_interior_face() {
        let tet_mesh = base_two_tetra_mesh();
        let mesh = build_3d_mesh_uniform(&tet_mesh, 2.0, 1000.0, 1000.0);

        assert_eq!(mesh.cells.len(), 2);
        // 2 tetra => 8 local faces, one shared pair -> 7 unique faces.
        assert_eq!(mesh.faces.len(), 7);

        let interior_faces: Vec<&FvmFace> = mesh
            .faces
            .iter()
            .filter(|f| f.cell_left != BOUNDARY && f.cell_right != BOUNDARY)
            .collect();
        assert_eq!(interior_faces.len(), 1);
        let f = interior_faces[0];

        // Shared face is triangle (0,1,2): area = 0.5.
        assert!((f.area - 0.5).abs() < 1e-12, "area={}", f.area);
        // Centroid distances are 0.25 on each side, so total is 0.5.
        assert!((f.distance - 0.5).abs() < 1e-12, "distance={}", f.distance);
        // Uniform k=2: K = A / (0.25/2 + 0.25/2) = 2.0 W/K.
        assert!((f.conductance - 2.0).abs() < 1e-12, "K={}", f.conductance);
    }

    #[test]
    fn test_two_tetra_heterogeneous_conductance() {
        let tet_mesh = base_two_tetra_mesh();
        let props = vec![
            CellThermalProperties::new(2.0, 1000.0, 1000.0),
            CellThermalProperties::new(4.0, 1000.0, 1000.0),
        ];
        let mesh = build_3d_mesh(&tet_mesh, &props);

        let interior = mesh
            .faces
            .iter()
            .find(|f| f.cell_left != BOUNDARY && f.cell_right != BOUNDARY)
            .expect("expected one interior face");

        // K = A / (d1/k1 + d2/k2) = 0.5 / (0.25/2 + 0.25/4) = 2.666...
        let expected = 0.5 / (0.25 / 2.0 + 0.25 / 4.0);
        assert!(
            (interior.conductance - expected).abs() < 1e-12,
            "K={}, expected={}",
            interior.conductance,
            expected
        );
    }

    #[test]
    #[should_panic(expected = "cell properties length")]
    fn test_property_count_must_match_tet_count() {
        let tet_mesh = base_two_tetra_mesh();
        let _ = build_3d_mesh(
            &tet_mesh,
            &[CellThermalProperties::new(1.0, 1000.0, 1000.0)],
        );
    }

    #[test]
    #[should_panic(expected = "Non-manifold tetrahedral mesh")]
    fn test_non_manifold_face_panics() {
        // Three tetrahedra sharing face (0,1,2) => non-manifold.
        let vertices = vec![
            Point::new(0.0, 0.0, 0.0),  // 0
            Point::new(1.0, 0.0, 0.0),  // 1
            Point::new(0.0, 1.0, 0.0),  // 2
            Point::new(0.0, 0.0, 1.0),  // 3
            Point::new(0.0, 0.0, -1.0), // 4
            Point::new(0.0, 0.0, 2.0),  // 5
        ];
        let tet_mesh = TetrahedralMesh::new(
            vertices,
            vec![
                TetrahedronIndex(0, 1, 2, 3),
                TetrahedronIndex(0, 1, 2, 4),
                TetrahedronIndex(0, 1, 2, 5),
            ],
        );
        let _ = build_3d_mesh_uniform(&tet_mesh, 2.0, 1000.0, 1000.0);
    }
}
