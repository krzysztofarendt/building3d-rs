use crate::sim::heat_transfer::boundary::BoundaryCondition;
use crate::sim::heat_transfer::mesh::{BOUNDARY, FvmFace, FvmMesh};
use std::collections::HashSet;

/// Configuration for the sparse FVM solver.
#[derive(Debug, Clone, Copy)]
pub struct SparseSolverConfig {
    /// Maximum number of PCG iterations per time step.
    pub max_iterations: usize,
    /// Relative residual tolerance.
    pub rel_tolerance: f64,
    /// Absolute residual tolerance.
    pub abs_tolerance: f64,
}

impl Default for SparseSolverConfig {
    fn default() -> Self {
        Self {
            max_iterations: 2000,
            rel_tolerance: 1e-9,
            abs_tolerance: 1e-12,
        }
    }
}

/// Generic sparse FVM solver for dimension-agnostic meshes.
///
/// Unlike [`crate::sim::heat_transfer::solver::FvmWallSolver`], this solver
/// does not assume a tridiagonal topology and uses preconditioned conjugate
/// gradients (Jacobi preconditioner) on the assembled sparse system.
#[derive(Debug, Clone)]
pub struct FvmSparseSolver {
    mesh: FvmMesh,
    temperatures: Vec<f64>,
    config: SparseSolverConfig,
}

impl FvmSparseSolver {
    /// Create a solver with default PCG settings.
    pub fn new(mesh: FvmMesh, initial_temperature: f64) -> Self {
        Self::with_config(mesh, initial_temperature, SparseSolverConfig::default())
    }

    /// Create a solver with custom PCG settings.
    pub fn with_config(
        mesh: FvmMesh,
        initial_temperature: f64,
        config: SparseSolverConfig,
    ) -> Self {
        let n = mesh.cells.len();
        Self {
            mesh,
            temperatures: vec![initial_temperature; n],
            config,
        }
    }

    /// Advance one implicit FVM step on a generic sparse mesh.
    ///
    /// `boundary_conditions` maps boundary face indices to boundary conditions.
    /// Boundary faces omitted from this list are treated as adiabatic.
    ///
    /// `sources` are per-cell source terms in watts.
    pub fn step(
        &mut self,
        dt: f64,
        boundary_conditions: &[(usize, BoundaryCondition)],
        sources: &[f64],
    ) -> &[f64] {
        assert!(dt > 0.0, "time step must be > 0");
        let n = self.mesh.cells.len();
        if n == 0 {
            return &self.temperatures;
        }

        // Sparse matrix A = C/dt + K and RHS b.
        let mut diag = vec![0.0; n];
        let mut off = vec![Vec::<(usize, f64)>::new(); n];
        let mut rhs = vec![0.0; n];

        for i in 0..n {
            let cap = self.mesh.cells[i].capacity();
            diag[i] += cap / dt;
            rhs[i] += (cap / dt) * self.temperatures[i];
            if i < sources.len() {
                rhs[i] += sources[i];
            }
        }

        for face in &self.mesh.faces {
            match boundary_adjacent_cell(face) {
                Some(_) => {
                    // BC contributions are applied explicitly from boundary_conditions.
                }
                None => {
                    // Interior coupling.
                    let i = face.cell_left;
                    let j = face.cell_right;
                    let k = face.conductance;
                    diag[i] += k;
                    diag[j] += k;
                    off[i].push((j, -k));
                    off[j].push((i, -k));
                }
            }
        }

        let mut seen = HashSet::new();
        for (face_idx, bc) in boundary_conditions {
            assert!(
                *face_idx < self.mesh.faces.len(),
                "boundary face index {} out of range ({} faces)",
                face_idx,
                self.mesh.faces.len()
            );
            assert!(
                seen.insert(*face_idx),
                "duplicate boundary condition for face index {}",
                face_idx
            );

            let face = &self.mesh.faces[*face_idx];
            let cell_idx = boundary_adjacent_cell(face)
                .expect("boundary condition can only be applied to boundary faces");
            apply_boundary_contribution(*bc, face, cell_idx, &mut diag, &mut rhs);
        }

        let x0 = self.temperatures.clone();
        self.temperatures = pcg_solve(&diag, &off, &rhs, &x0, self.config);
        &self.temperatures
    }

    /// Access current cell temperatures.
    pub fn temperatures(&self) -> &[f64] {
        &self.temperatures
    }
}

fn boundary_adjacent_cell(face: &FvmFace) -> Option<usize> {
    let left_boundary = face.cell_left == BOUNDARY;
    let right_boundary = face.cell_right == BOUNDARY;
    match (left_boundary, right_boundary) {
        (true, true) => panic!("invalid face: both sides are boundary"),
        (true, false) => Some(face.cell_right),
        (false, true) => Some(face.cell_left),
        (false, false) => None,
    }
}

fn apply_boundary_contribution(
    bc: BoundaryCondition,
    face: &FvmFace,
    cell_idx: usize,
    diag: &mut [f64],
    rhs: &mut [f64],
) {
    match bc {
        BoundaryCondition::Dirichlet { temperature } => {
            diag[cell_idx] += face.conductance;
            rhs[cell_idx] += face.conductance * temperature;
        }
        BoundaryCondition::Neumann { heat_flux } => {
            rhs[cell_idx] += heat_flux * face.area;
        }
        BoundaryCondition::Convective { h, t_fluid } => {
            let h_a = h * face.area;
            if h_a <= 0.0 {
                return;
            }
            let k_eff = if face.conductance > 0.0 {
                1.0 / (1.0 / h_a + 1.0 / face.conductance)
            } else {
                h_a
            };
            diag[cell_idx] += k_eff;
            rhs[cell_idx] += k_eff * t_fluid;
        }
        BoundaryCondition::ConvectiveWithFlux {
            h,
            t_fluid,
            heat_flux,
        } => {
            let h_a = h * face.area;
            if h_a <= 0.0 {
                rhs[cell_idx] += heat_flux * face.area;
                return;
            }
            let k_eff = if face.conductance > 0.0 {
                1.0 / (1.0 / h_a + 1.0 / face.conductance)
            } else {
                h_a
            };
            diag[cell_idx] += k_eff;
            rhs[cell_idx] += k_eff * t_fluid;
            let alpha = if face.conductance > 0.0 {
                face.conductance / (face.conductance + h_a)
            } else {
                0.0
            };
            rhs[cell_idx] += alpha * heat_flux * face.area;
        }
    }
}

fn pcg_solve(
    diag: &[f64],
    off: &[Vec<(usize, f64)>],
    b: &[f64],
    x0: &[f64],
    config: SparseSolverConfig,
) -> Vec<f64> {
    let n = b.len();
    if n == 0 {
        return Vec::new();
    }

    let mut x = x0.to_vec();
    let mut ax = vec![0.0; n];
    apply_matrix(diag, off, &x, &mut ax);

    let mut r = vec![0.0; n];
    for i in 0..n {
        r[i] = b[i] - ax[i];
    }

    let b_norm = l2_norm(b).max(1.0);
    let tol = config.abs_tolerance.max(config.rel_tolerance * b_norm);
    if l2_norm(&r) <= tol {
        return x;
    }

    let mut z = vec![0.0; n];
    for i in 0..n {
        z[i] = if diag[i].abs() > 1e-30 {
            r[i] / diag[i]
        } else {
            r[i]
        };
    }
    let mut p = z.clone();
    let mut rz_old = dot(&r, &z);
    if rz_old.abs() < 1e-30 {
        return x;
    }

    let mut ap = vec![0.0; n];
    for _ in 0..config.max_iterations {
        apply_matrix(diag, off, &p, &mut ap);
        let denom = dot(&p, &ap);
        if denom.abs() < 1e-30 {
            break;
        }

        let alpha = rz_old / denom;
        for i in 0..n {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }

        if l2_norm(&r) <= tol {
            break;
        }

        for i in 0..n {
            z[i] = if diag[i].abs() > 1e-30 {
                r[i] / diag[i]
            } else {
                r[i]
            };
        }
        let rz_new = dot(&r, &z);
        if rz_old.abs() < 1e-30 {
            break;
        }
        let beta = rz_new / rz_old;
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
        rz_old = rz_new;
    }

    x
}

fn apply_matrix(diag: &[f64], off: &[Vec<(usize, f64)>], x: &[f64], y: &mut [f64]) {
    for i in 0..diag.len() {
        let mut sum = diag[i] * x[i];
        for (j, a_ij) in &off[i] {
            sum += a_ij * x[*j];
        }
        y[i] = sum;
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn l2_norm(a: &[f64]) -> f64 {
    dot(a, a).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Point;
    use crate::geom::mesh::tetrahedralize::{TetrahedralMesh, TetrahedronIndex};
    use crate::sim::heat_transfer::mesh::{FvmCell, FvmFace};
    use crate::sim::heat_transfer::mesh_3d::build_3d_mesh_uniform;

    #[test]
    fn test_sparse_solver_two_cell_steady_state_dirichlet() {
        // Two-cell chain with one interior and two boundary faces.
        let mesh = FvmMesh {
            cells: vec![
                FvmCell {
                    volume: 1.0,
                    conductivity: 1.0,
                    density: 1.0,
                    specific_heat: 1.0,
                },
                FvmCell {
                    volume: 1.0,
                    conductivity: 1.0,
                    density: 1.0,
                    specific_heat: 1.0,
                },
            ],
            faces: vec![
                FvmFace {
                    cell_left: BOUNDARY,
                    cell_right: 0,
                    area: 1.0,
                    distance: 0.5,
                    conductance: 2.0,
                },
                FvmFace {
                    cell_left: 0,
                    cell_right: 1,
                    area: 1.0,
                    distance: 1.0,
                    conductance: 1.0,
                },
                FvmFace {
                    cell_left: 1,
                    cell_right: BOUNDARY,
                    area: 1.0,
                    distance: 0.5,
                    conductance: 2.0,
                },
            ],
        };

        let mut solver = FvmSparseSolver::new(mesh, 10.0);
        let bcs = vec![
            (0, BoundaryCondition::Dirichlet { temperature: 0.0 }),
            (2, BoundaryCondition::Dirichlet { temperature: 30.0 }),
        ];
        let sources = vec![0.0, 0.0];

        for _ in 0..4 {
            solver.step(1e8, &bcs, &sources);
        }

        let t = solver.temperatures();
        assert!((t[0] - 7.5).abs() < 1e-6, "t0={}", t[0]);
        assert!((t[1] - 22.5).abs() < 1e-6, "t1={}", t[1]);
    }

    #[test]
    fn test_sparse_solver_works_with_3d_mesh_builder_output() {
        let tet_mesh = TetrahedralMesh::new(
            vec![
                Point::new(0.0, 0.0, 0.0),
                Point::new(1.0, 0.0, 0.0),
                Point::new(0.0, 1.0, 0.0),
                Point::new(0.0, 0.0, 1.0),
                Point::new(0.0, 0.0, -1.0),
            ],
            vec![TetrahedronIndex(0, 1, 2, 3), TetrahedronIndex(0, 1, 2, 4)],
        );
        let mesh = build_3d_mesh_uniform(&tet_mesh, 2.0, 1.0, 1.0);

        let mut solver = FvmSparseSolver::new(mesh.clone(), 0.0);
        let bcs: Vec<(usize, BoundaryCondition)> = mesh
            .faces
            .iter()
            .enumerate()
            .filter(|(_, f)| f.cell_left == BOUNDARY || f.cell_right == BOUNDARY)
            .map(|(idx, _)| (idx, BoundaryCondition::Dirichlet { temperature: 20.0 }))
            .collect();
        let sources = vec![0.0; mesh.cells.len()];

        for _ in 0..3 {
            solver.step(1e8, &bcs, &sources);
        }
        for &ti in solver.temperatures() {
            assert!((ti - 20.0).abs() < 1e-6, "T={ti}");
        }
    }

    #[test]
    #[should_panic(expected = "boundary condition can only be applied to boundary faces")]
    fn test_sparse_solver_rejects_bc_on_interior_face() {
        let mesh = FvmMesh {
            cells: vec![
                FvmCell {
                    volume: 1.0,
                    conductivity: 1.0,
                    density: 1.0,
                    specific_heat: 1.0,
                },
                FvmCell {
                    volume: 1.0,
                    conductivity: 1.0,
                    density: 1.0,
                    specific_heat: 1.0,
                },
            ],
            faces: vec![FvmFace {
                cell_left: 0,
                cell_right: 1,
                area: 1.0,
                distance: 1.0,
                conductance: 1.0,
            }],
        };
        let mut solver = FvmSparseSolver::new(mesh, 10.0);
        let _ = solver.step(
            60.0,
            &[(0, BoundaryCondition::Dirichlet { temperature: 0.0 })],
            &[0.0, 0.0],
        );
    }
}
