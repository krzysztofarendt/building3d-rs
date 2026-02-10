use crate::sim::energy::construction::WallConstruction;
use crate::sim::heat_transfer::mesh::{FvmCell, FvmFace, FvmMesh, BOUNDARY};

/// Maximum thickness (m) of a single cell.  Layers thicker than this are
/// subdivided for accuracy.
const MAX_CELL_THICKNESS: f64 = 0.05;

/// Build a 1D finite-volume mesh from a [`WallConstruction`].
///
/// Layers are ordered **exterior → interior** (matching `WallConstruction`).
/// Each layer produces at least one cell; thick layers are subdivided so that
/// no cell exceeds [`MAX_CELL_THICKNESS`].
///
/// `wall_area` is the polygon area in m^2 — it sets the face area for all
/// faces in the 1D mesh and scales cell volumes.
pub fn build_1d_mesh(construction: &WallConstruction, wall_area: f64) -> FvmMesh {
    let mut cells = Vec::new();

    // Subdivide layers into cells
    for layer in &construction.layers {
        let n = (layer.thickness / MAX_CELL_THICKNESS).ceil().max(1.0) as usize;
        let dx = layer.thickness / n as f64;
        for _ in 0..n {
            cells.push(FvmCell {
                volume: wall_area * dx,
                conductivity: layer.conductivity,
                density: layer.density,
                specific_heat: layer.specific_heat,
            });
        }
    }

    let n = cells.len();
    let mut faces = Vec::with_capacity(n + 1);

    // Exterior boundary face (BOUNDARY | cell_0)
    if n > 0 {
        let half_dx = cell_thickness(&cells[0], wall_area) / 2.0;
        faces.push(FvmFace {
            cell_left: BOUNDARY,
            cell_right: 0,
            area: wall_area,
            distance: half_dx,
            conductance: cells[0].conductivity * wall_area / half_dx,
        });
    }

    // Interior faces between adjacent cells.
    // The conductance uses the series-resistance formula:
    //   K = A / (half_dx_L / k_L + half_dx_R / k_R)
    // This reduces to the harmonic-mean formula when half-distances are equal,
    // and correctly handles material interfaces with different cell sizes.
    for i in 0..n.saturating_sub(1) {
        let dx_l = cell_thickness(&cells[i], wall_area) / 2.0;
        let dx_r = cell_thickness(&cells[i + 1], wall_area) / 2.0;
        let distance = dx_l + dx_r;
        let r_face = dx_l / cells[i].conductivity + dx_r / cells[i + 1].conductivity;
        faces.push(FvmFace {
            cell_left: i,
            cell_right: i + 1,
            area: wall_area,
            distance,
            conductance: wall_area / r_face,
        });
    }

    // Interior boundary face (cell_N | BOUNDARY)
    if n > 0 {
        let last = n - 1;
        let half_dx = cell_thickness(&cells[last], wall_area) / 2.0;
        faces.push(FvmFace {
            cell_left: last,
            cell_right: BOUNDARY,
            area: wall_area,
            distance: half_dx,
            conductance: cells[last].conductivity * wall_area / half_dx,
        });
    }

    FvmMesh { cells, faces }
}

/// Recover cell thickness from volume and wall area.
fn cell_thickness(cell: &FvmCell, wall_area: f64) -> f64 {
    cell.volume / wall_area
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::materials::Layer;

    #[test]
    fn test_single_layer_mesh() {
        let construction = WallConstruction::new(
            "test",
            vec![Layer {
                name: "concrete".into(),
                thickness: 0.20,
                conductivity: 1.4,
                density: 2300.0,
                specific_heat: 880.0,
            }],
        );
        let mesh = build_1d_mesh(&construction, 10.0);

        // 0.20 m / 0.05 max = 4 cells
        assert_eq!(mesh.cells.len(), 4);
        // 4 cells => 3 interior faces + 2 boundary faces = 5
        assert_eq!(mesh.faces.len(), 5);

        // Total thickness
        let total: f64 = mesh
            .cells
            .iter()
            .map(|c| cell_thickness(c, 10.0))
            .sum();
        assert!((total - 0.20).abs() < 1e-12, "total thickness = {total}");
    }

    #[test]
    fn test_three_layer_mesh() {
        let construction = WallConstruction::new(
            "insulated",
            vec![
                Layer {
                    name: "plaster".into(),
                    thickness: 0.02,
                    conductivity: 0.87,
                    density: 1800.0,
                    specific_heat: 840.0,
                },
                Layer {
                    name: "insulation".into(),
                    thickness: 0.10,
                    conductivity: 0.04,
                    density: 30.0,
                    specific_heat: 1030.0,
                },
                Layer {
                    name: "concrete".into(),
                    thickness: 0.15,
                    conductivity: 1.4,
                    density: 2300.0,
                    specific_heat: 880.0,
                },
            ],
        );
        let area = 12.0;
        let mesh = build_1d_mesh(&construction, area);

        // plaster: 0.02/0.05 => 1 cell
        // insulation: 0.10/0.05 => 2 cells
        // concrete: 0.15/0.05 => 3 cells
        // total: 6 cells
        assert_eq!(mesh.cells.len(), 6);
        // 6 cells => 5 interior + 2 boundary = 7 faces
        assert_eq!(mesh.faces.len(), 7);

        // Total thickness
        let total: f64 = mesh
            .cells
            .iter()
            .map(|c| cell_thickness(c, area))
            .sum();
        assert!(
            (total - 0.27).abs() < 1e-12,
            "total thickness = {total}, expected 0.27"
        );

        // Boundary faces
        assert_eq!(
            mesh.exterior_boundary_face(),
            Some(0),
            "exterior boundary should be face 0"
        );
        assert_eq!(
            mesh.interior_boundary_face(),
            Some(6),
            "interior boundary should be last face"
        );
    }

    #[test]
    fn test_thin_layer_gets_one_cell() {
        let construction = WallConstruction::new(
            "thin",
            vec![Layer {
                name: "membrane".into(),
                thickness: 0.001,
                conductivity: 0.2,
                density: 1000.0,
                specific_heat: 1000.0,
            }],
        );
        let mesh = build_1d_mesh(&construction, 1.0);
        assert_eq!(mesh.cells.len(), 1);
        assert_eq!(mesh.faces.len(), 2); // 2 boundary faces
    }

    #[test]
    fn test_face_conductance_same_material() {
        // Two cells of same material: series resistance = harmonic mean
        let construction = WallConstruction::new(
            "test",
            vec![Layer {
                name: "concrete".into(),
                thickness: 0.10,
                conductivity: 1.4,
                density: 2300.0,
                specific_heat: 880.0,
            }],
        );
        let mesh = build_1d_mesh(&construction, 1.0);
        // 0.10m / 0.05 = 2 cells, 1 interior face
        let interior_face = &mesh.faces[1]; // face between cell 0 and cell 1
        // K = A / (0.025/1.4 + 0.025/1.4) = 1.0 / (0.03571) = 28.0
        assert!(
            (interior_face.conductance - 28.0).abs() < 0.01,
            "got {}",
            interior_face.conductance
        );
    }
}
