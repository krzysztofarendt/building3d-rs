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

    /// Heat flux at a boundary face [W/m²].
    /// Positive = heat flows from domain outward through this face.
    pub fn boundary_heat_flux(&self, face_idx: usize, bc: &BoundaryCondition) -> f64 {
        let (cell_idx, face) = get_boundary_face_info(&self.mesh, face_idx);
        let t_centroid = self.temperatures[cell_idx];
        match *bc {
            BoundaryCondition::Dirichlet { temperature } => {
                let k_per_area = face.conductance / face.area;
                k_per_area * (t_centroid - temperature)
            }
            BoundaryCondition::Neumann { heat_flux } => -heat_flux,
            BoundaryCondition::Convective { h, t_fluid }
            | BoundaryCondition::ConvectiveWithFlux {
                h,
                t_fluid,
                heat_flux: _,
            }
            | BoundaryCondition::ConvectiveWithFluxToDomain {
                h,
                t_fluid,
                heat_flux: _,
            } => {
                if h <= 0.0 {
                    return 0.0;
                }
                let k_cond = face.conductance / face.area;
                let k_combined = 1.0 / (1.0 / k_cond + 1.0 / h);
                k_combined * (t_centroid - t_fluid)
            }
        }
    }

    /// Reconstructed surface temperature at a boundary face [°C].
    /// Accounts for half-cell conduction resistance between centroid and surface.
    pub fn boundary_surface_temperature(&self, face_idx: usize, bc: &BoundaryCondition) -> f64 {
        let (cell_idx, face) = get_boundary_face_info(&self.mesh, face_idx);
        let t_centroid = self.temperatures[cell_idx];
        let k_face_per_area = (face.conductance / face.area).max(0.0);
        match *bc {
            BoundaryCondition::Dirichlet { temperature } => temperature,
            BoundaryCondition::Neumann { .. } => t_centroid,
            BoundaryCondition::Convective { h, t_fluid } => {
                if h <= 0.0 {
                    return t_centroid;
                }
                (k_face_per_area * t_centroid + h * t_fluid) / (k_face_per_area + h)
            }
            BoundaryCondition::ConvectiveWithFlux {
                h,
                t_fluid,
                heat_flux,
            } => {
                if h <= 0.0 {
                    return t_centroid;
                }
                let alpha = self.boundary_conduction_fraction(face_idx, h);
                let q_into_wall = alpha * heat_flux;
                (k_face_per_area * t_centroid + h * t_fluid + q_into_wall)
                    / (k_face_per_area + h)
            }
            BoundaryCondition::ConvectiveWithFluxToDomain {
                h,
                t_fluid,
                heat_flux: _,
            } => {
                if h <= 0.0 {
                    return t_centroid;
                }
                (k_face_per_area * t_centroid + h * t_fluid) / (k_face_per_area + h)
            }
        }
    }

    /// Total thermal capacity of all cells [J/K].
    pub fn total_capacity_j_per_k(&self) -> f64 {
        self.mesh.cells.iter().map(|c| c.capacity()).sum()
    }

    /// Fraction (0..1) of a surface source that enters the domain
    /// when the boundary has convective film coefficient h [W/(m²·K)].
    pub fn boundary_conduction_fraction(&self, face_idx: usize, h: f64) -> f64 {
        let (_cell_idx, face) = get_boundary_face_info(&self.mesh, face_idx);
        let h_a = h.max(0.0) * face.area;
        if h_a <= 0.0 {
            return 1.0;
        }
        (face.conductance / (face.conductance + h_a)).clamp(0.0, 1.0)
    }
}

/// Returns `(cell_idx, &FvmFace)` for a boundary face, or panics with a clear message.
fn get_boundary_face_info(mesh: &FvmMesh, face_idx: usize) -> (usize, &FvmFace) {
    assert!(
        face_idx < mesh.faces.len(),
        "face index {} out of range ({} faces)",
        face_idx,
        mesh.faces.len()
    );
    let face = &mesh.faces[face_idx];
    let cell_idx = boundary_adjacent_cell(face).unwrap_or_else(|| {
        panic!(
            "face {} is not a boundary face (cell_left={}, cell_right={})",
            face_idx, face.cell_left, face.cell_right
        )
    });
    (cell_idx, face)
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
        BoundaryCondition::ConvectiveWithFluxToDomain {
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
            // Full flux enters the domain (no convective split)
            rhs[cell_idx] += heat_flux * face.area;
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

    /// Helper: build a two-cell chain mesh for boundary method tests.
    /// Layout: [boundary|--cell0--|--cell1--|boundary]
    /// Each cell has volume = area * dx, conductance at boundary = k*A/half_dx.
    fn two_cell_chain(k: f64, rho: f64, cp: f64, dx: f64, area: f64) -> FvmMesh {
        let half_dx = dx / 2.0;
        FvmMesh {
            cells: vec![
                FvmCell {
                    volume: area * dx,
                    conductivity: k,
                    density: rho,
                    specific_heat: cp,
                },
                FvmCell {
                    volume: area * dx,
                    conductivity: k,
                    density: rho,
                    specific_heat: cp,
                },
            ],
            faces: vec![
                // face 0: left boundary -> cell 0
                FvmFace {
                    cell_left: BOUNDARY,
                    cell_right: 0,
                    area,
                    distance: half_dx,
                    conductance: k * area / half_dx,
                },
                // face 1: cell 0 -> cell 1 (interior)
                FvmFace {
                    cell_left: 0,
                    cell_right: 1,
                    area,
                    distance: dx,
                    conductance: k * area / dx,
                },
                // face 2: cell 1 -> right boundary
                FvmFace {
                    cell_left: 1,
                    cell_right: BOUNDARY,
                    area,
                    distance: half_dx,
                    conductance: k * area / half_dx,
                },
            ],
        }
    }

    #[test]
    fn test_sparse_boundary_heat_flux_dirichlet() {
        let k = 1.4;
        let dx = 0.1;
        let area = 1.0;
        let mesh = two_cell_chain(k, 2300.0, 880.0, dx, area);
        let mut solver = FvmSparseSolver::new(mesh, 10.0);

        let t_left = 0.0;
        let t_right = 30.0;
        let bc_left = BoundaryCondition::Dirichlet { temperature: t_left };
        let bc_right = BoundaryCondition::Dirichlet { temperature: t_right };
        let bcs = vec![(0, bc_left), (2, bc_right)];
        let sources = vec![0.0, 0.0];

        // Run to steady state
        for _ in 0..10 {
            solver.step(1e8, &bcs, &sources);
        }

        // Steady-state flux: q = k * dT / L = 1.4 * 30 / 0.2 = 210 W/m²
        let expected_flux = k * (t_right - t_left) / (2.0 * dx);

        // Left boundary: heat flows from domain outward (domain is warmer than 0°C)
        // so flux should be positive-ish at left... actually let's think:
        // cell0 centroid T ~ 7.5, bc T = 0 => flux = k_per_area * (7.5 - 0) > 0 => out
        // But the expected uniform flux is q = k * 30 / 0.2 = 210 W/m²
        let flux_left = solver.boundary_heat_flux(0, &bc_left);
        assert!(
            (flux_left - expected_flux).abs() < 1.0,
            "left flux = {flux_left}, expected {expected_flux}"
        );

        // Right boundary: heat flows into domain from the hot side
        // cell1 centroid T ~ 22.5, bc T = 30 => flux = k_per_area * (22.5 - 30) < 0 => into domain
        let flux_right = solver.boundary_heat_flux(2, &bc_right);
        assert!(
            (flux_right + expected_flux).abs() < 1.0,
            "right flux = {flux_right}, expected {}",
            -expected_flux
        );
    }

    #[test]
    fn test_sparse_boundary_heat_flux_convective() {
        let k = 1.4;
        let dx = 0.1;
        let area = 1.0;
        let h_left = 25.0;
        let h_right = 7.7;
        let t_left = -10.0;
        let t_right = 20.0;

        let mesh = two_cell_chain(k, 2300.0, 880.0, dx, area);
        let mut solver = FvmSparseSolver::new(mesh, 10.0);

        let bc_left = BoundaryCondition::Convective {
            h: h_left,
            t_fluid: t_left,
        };
        let bc_right = BoundaryCondition::Convective {
            h: h_right,
            t_fluid: t_right,
        };
        let bcs = vec![(0, bc_left), (2, bc_right)];
        let sources = vec![0.0, 0.0];

        for _ in 0..20 {
            solver.step(1e8, &bcs, &sources);
        }

        // R_total = 1/h_left + L/k + 1/h_right
        let l_total = 2.0 * dx;
        let r_total = 1.0 / h_left + l_total / k + 1.0 / h_right;
        let expected_flux = (t_left - t_right) / r_total; // negative (heat flows right to left)

        let flux_right = solver.boundary_heat_flux(2, &bc_right);
        // flux_right should be negative (heat enters domain from right)
        assert!(
            (flux_right - expected_flux).abs() < expected_flux.abs() * 0.01,
            "right flux = {flux_right}, expected {expected_flux}"
        );
    }

    #[test]
    fn test_sparse_boundary_surface_temp_dirichlet() {
        let mesh = two_cell_chain(1.4, 2300.0, 880.0, 0.1, 1.0);
        let mut solver = FvmSparseSolver::new(mesh, 10.0);

        let bc = BoundaryCondition::Dirichlet { temperature: 42.0 };
        let bcs = vec![
            (0, bc),
            (2, BoundaryCondition::Dirichlet { temperature: 10.0 }),
        ];
        solver.step(1e8, &bcs, &[0.0, 0.0]);

        let t_surf = solver.boundary_surface_temperature(0, &bc);
        assert!(
            (t_surf - 42.0).abs() < 1e-12,
            "Dirichlet surface temp should be prescribed: got {t_surf}"
        );
    }

    #[test]
    fn test_sparse_boundary_surface_temp_convective() {
        let k = 1.4;
        let dx = 0.1;
        let area = 1.0;
        let h = 7.7;
        let t_fluid = 20.0;

        let mesh = two_cell_chain(k, 2300.0, 880.0, dx, area);
        let mut solver = FvmSparseSolver::new(mesh, 10.0);

        let bc_left = BoundaryCondition::Dirichlet { temperature: -10.0 };
        let bc_right = BoundaryCondition::Convective { h, t_fluid };
        let bcs = vec![(0, bc_left), (2, bc_right)];

        for _ in 0..20 {
            solver.step(1e8, &bcs, &[0.0, 0.0]);
        }

        let t_centroid = solver.temperatures()[1];
        let t_surf = solver.boundary_surface_temperature(2, &bc_right);

        // Surface temp should be between centroid and fluid temp
        assert!(
            t_surf > t_centroid.min(t_fluid) - 1e-6
                && t_surf < t_centroid.max(t_fluid) + 1e-6,
            "surface temp {t_surf} should be between centroid {t_centroid} and fluid {t_fluid}"
        );

        // Verify formula: T_s = (k_face * T_centroid + h * T_fluid) / (k_face + h)
        let k_face_per_area = k / (dx / 2.0); // conductance / area
        let expected = (k_face_per_area * t_centroid + h * t_fluid) / (k_face_per_area + h);
        assert!(
            (t_surf - expected).abs() < 1e-6,
            "surface temp {t_surf} != expected {expected}"
        );
    }

    #[test]
    fn test_sparse_total_capacity() {
        let k = 1.4;
        let rho = 2300.0;
        let cp = 880.0;
        let dx = 0.1;
        let area = 1.0;
        let mesh = two_cell_chain(k, rho, cp, dx, area);
        let solver = FvmSparseSolver::new(mesh, 10.0);

        let expected = 2.0 * rho * cp * (area * dx);
        let actual = solver.total_capacity_j_per_k();
        assert!(
            (actual - expected).abs() < 1e-6,
            "capacity = {actual}, expected {expected}"
        );
    }

    #[test]
    fn test_sparse_boundary_conduction_fraction() {
        let k = 1.4;
        let dx = 0.1;
        let area = 1.0;
        let h = 25.0;

        let mesh = two_cell_chain(k, 2300.0, 880.0, dx, area);
        let solver = FvmSparseSolver::new(mesh, 10.0);

        // face 0: boundary face, conductance = k * area / half_dx = 1.4 * 1.0 / 0.05 = 28.0
        // h_a = 25.0 * 1.0 = 25.0
        // fraction = 28.0 / (28.0 + 25.0) = 28/53
        let expected = 28.0 / (28.0 + 25.0);
        let actual = solver.boundary_conduction_fraction(0, h);
        assert!(
            (actual - expected).abs() < 1e-10,
            "fraction = {actual}, expected {expected}"
        );

        // With h=0, all flux enters the domain
        assert!(
            (solver.boundary_conduction_fraction(0, 0.0) - 1.0).abs() < 1e-10,
            "fraction with h=0 should be 1.0"
        );
    }

    #[test]
    fn test_sparse_neumann_flux_steady_state() {
        // Two-cell chain: Neumann on left (flux into domain), Dirichlet on right.
        // At steady state the imposed flux flows through both cells to the fixed boundary.
        let k = 1.4;
        let dx = 0.1;
        let area = 1.0;
        let q_in = 100.0; // W/m² into domain at left boundary
        let t_right = 20.0;

        let mesh = two_cell_chain(k, 2300.0, 880.0, dx, area);
        let mut solver = FvmSparseSolver::new(mesh, 10.0);

        let bcs = vec![
            (0, BoundaryCondition::Neumann { heat_flux: q_in }),
            (
                2,
                BoundaryCondition::Dirichlet {
                    temperature: t_right,
                },
            ),
        ];
        let sources = vec![0.0, 0.0];

        for _ in 0..20 {
            solver.step(1e8, &bcs, &sources);
        }

        // At steady state, the flux q passes through both cells.
        // T(x) = T_right + q * (L - x) / k, where L = 2*dx, x measured from left.
        // Cell 0 centroid at x = dx/2:  T0 = 20 + 100 * (0.2 - 0.05) / 1.4
        // Cell 1 centroid at x = 3*dx/2: T1 = 20 + 100 * (0.2 - 0.15) / 1.4
        let l_total = 2.0 * dx;
        let t0_expected = t_right + q_in * (l_total - dx / 2.0) / k;
        let t1_expected = t_right + q_in * (l_total - 3.0 * dx / 2.0) / k;

        let t = solver.temperatures();
        assert!(
            (t[0] - t0_expected).abs() < 0.5,
            "t0={}, expected {}",
            t[0],
            t0_expected
        );
        assert!(
            (t[1] - t1_expected).abs() < 0.5,
            "t1={}, expected {}",
            t[1],
            t1_expected
        );

        // Verify reported heat flux at left boundary (positive = outward = -q_in)
        let flux_left = solver.boundary_heat_flux(
            0,
            &BoundaryCondition::Neumann { heat_flux: q_in },
        );
        assert!(
            (flux_left + q_in).abs() < 1e-10,
            "left flux = {flux_left}, expected {}",
            -q_in
        );
    }

    /// ConvectiveWithFlux: surface source splits between convection and conduction.
    /// At steady state, ConvectiveWithFlux with identical convective parameters should
    /// produce a warmer wall than plain Convective (no flux), colder than
    /// ConvectiveWithFluxToDomain (full flux into domain).
    #[test]
    fn test_sparse_convective_with_flux_ordering() {
        let k = 1.4;
        let dx = 0.1;
        let area = 1.0;
        let h = 25.0;
        let t_fluid = 0.0;
        let t_right = 20.0;
        let q_solar = 200.0; // W/m² surface source

        let mesh = two_cell_chain(k, 2300.0, 880.0, dx, area);

        // 1. Plain Convective (no solar)
        let mut solver_plain = FvmSparseSolver::new(mesh.clone(), 10.0);
        let bcs_plain = vec![
            (
                0,
                BoundaryCondition::Convective {
                    h,
                    t_fluid,
                },
            ),
            (
                2,
                BoundaryCondition::Dirichlet {
                    temperature: t_right,
                },
            ),
        ];
        let sources = vec![0.0, 0.0];
        for _ in 0..30 {
            solver_plain.step(1e8, &bcs_plain, &sources);
        }

        // 2. ConvectiveWithFlux (alpha fraction enters wall)
        let mut solver_split = FvmSparseSolver::new(mesh.clone(), 10.0);
        let bcs_split = vec![
            (
                0,
                BoundaryCondition::ConvectiveWithFlux {
                    h,
                    t_fluid,
                    heat_flux: q_solar,
                },
            ),
            (
                2,
                BoundaryCondition::Dirichlet {
                    temperature: t_right,
                },
            ),
        ];
        for _ in 0..30 {
            solver_split.step(1e8, &bcs_split, &sources);
        }

        // 3. ConvectiveWithFluxToDomain (100% enters wall)
        let mut solver_full = FvmSparseSolver::new(mesh, 10.0);
        let bcs_full = vec![
            (
                0,
                BoundaryCondition::ConvectiveWithFluxToDomain {
                    h,
                    t_fluid,
                    heat_flux: q_solar,
                },
            ),
            (
                2,
                BoundaryCondition::Dirichlet {
                    temperature: t_right,
                },
            ),
        ];
        for _ in 0..30 {
            solver_full.step(1e8, &bcs_full, &sources);
        }

        // Ordering: plain < split < full for all cell temperatures
        let t_plain = solver_plain.temperatures();
        let t_split = solver_split.temperatures();
        let t_full = solver_full.temperatures();

        for i in 0..2 {
            assert!(
                t_split[i] > t_plain[i],
                "cell {i}: split={} should be > plain={}",
                t_split[i],
                t_plain[i]
            );
            assert!(
                t_full[i] > t_split[i],
                "cell {i}: full={} should be > split={}",
                t_full[i],
                t_split[i]
            );
        }
    }

    /// ConvectiveWithFlux with h=0 should behave identically to pure Neumann.
    #[test]
    fn test_sparse_convective_with_flux_zero_h_is_neumann() {
        let k = 1.4;
        let dx = 0.1;
        let area = 1.0;
        let q_flux = 50.0;
        let t_right = 20.0;

        let mesh = two_cell_chain(k, 2300.0, 880.0, dx, area);

        // ConvectiveWithFlux, h=0
        let mut solver_cwf = FvmSparseSolver::new(mesh.clone(), 10.0);
        let bcs_cwf = vec![
            (
                0,
                BoundaryCondition::ConvectiveWithFlux {
                    h: 0.0,
                    t_fluid: 0.0,
                    heat_flux: q_flux,
                },
            ),
            (
                2,
                BoundaryCondition::Dirichlet {
                    temperature: t_right,
                },
            ),
        ];
        let sources = vec![0.0, 0.0];
        for _ in 0..20 {
            solver_cwf.step(1e8, &bcs_cwf, &sources);
        }

        // Pure Neumann
        let mut solver_nm = FvmSparseSolver::new(mesh, 10.0);
        let bcs_nm = vec![
            (0, BoundaryCondition::Neumann { heat_flux: q_flux }),
            (
                2,
                BoundaryCondition::Dirichlet {
                    temperature: t_right,
                },
            ),
        ];
        for _ in 0..20 {
            solver_nm.step(1e8, &bcs_nm, &sources);
        }

        let t_cwf = solver_cwf.temperatures();
        let t_nm = solver_nm.temperatures();
        for i in 0..2 {
            assert!(
                (t_cwf[i] - t_nm[i]).abs() < 1e-10,
                "cell {i}: cwf={} vs neumann={}",
                t_cwf[i],
                t_nm[i]
            );
        }
    }

    /// ConvectiveWithFluxToDomain with h=0 should also degenerate to pure Neumann.
    #[test]
    fn test_sparse_convective_with_flux_to_domain_zero_h_is_neumann() {
        let k = 1.4;
        let dx = 0.1;
        let area = 1.0;
        let q_flux = 50.0;
        let t_right = 20.0;

        let mesh = two_cell_chain(k, 2300.0, 880.0, dx, area);

        // ConvectiveWithFluxToDomain, h=0
        let mut solver_td = FvmSparseSolver::new(mesh.clone(), 10.0);
        let bcs_td = vec![
            (
                0,
                BoundaryCondition::ConvectiveWithFluxToDomain {
                    h: 0.0,
                    t_fluid: 0.0,
                    heat_flux: q_flux,
                },
            ),
            (
                2,
                BoundaryCondition::Dirichlet {
                    temperature: t_right,
                },
            ),
        ];
        let sources = vec![0.0, 0.0];
        for _ in 0..20 {
            solver_td.step(1e8, &bcs_td, &sources);
        }

        // Pure Neumann
        let mut solver_nm = FvmSparseSolver::new(mesh, 10.0);
        let bcs_nm = vec![
            (0, BoundaryCondition::Neumann { heat_flux: q_flux }),
            (
                2,
                BoundaryCondition::Dirichlet {
                    temperature: t_right,
                },
            ),
        ];
        for _ in 0..20 {
            solver_nm.step(1e8, &bcs_nm, &sources);
        }

        let t_td = solver_td.temperatures();
        let t_nm = solver_nm.temperatures();
        for i in 0..2 {
            assert!(
                (t_td[i] - t_nm[i]).abs() < 1e-10,
                "cell {i}: to_domain={} vs neumann={}",
                t_td[i],
                t_nm[i]
            );
        }
    }

    /// Verify surface temperature reconstruction for ConvectiveWithFlux.
    #[test]
    fn test_sparse_surface_temp_convective_with_flux() {
        let k = 1.4;
        let dx = 0.1;
        let area = 1.0;
        let h = 25.0;
        let t_fluid = 0.0;
        let q_solar = 200.0;

        let mesh = two_cell_chain(k, 2300.0, 880.0, dx, area);
        let mut solver = FvmSparseSolver::new(mesh, 10.0);

        let bc_left = BoundaryCondition::ConvectiveWithFlux {
            h,
            t_fluid,
            heat_flux: q_solar,
        };
        let bcs = vec![
            (0, bc_left),
            (
                2,
                BoundaryCondition::Dirichlet { temperature: 20.0 },
            ),
        ];
        let sources = vec![0.0, 0.0];

        for _ in 0..30 {
            solver.step(1e8, &bcs, &sources);
        }

        let t_centroid = solver.temperatures()[0];
        let t_surf = solver.boundary_surface_temperature(0, &bc_left);

        // Surface temp should be between centroid and fluid, but boosted by the
        // solar source fraction
        assert!(
            t_surf > t_fluid,
            "surface temp {t_surf} should exceed fluid temp {t_fluid} due to solar"
        );
        // Surface temp should be between centroid and fluid (in the general sense
        // that it lies on the conduction-convection path)
        assert!(
            t_surf != t_centroid,
            "surface temp should differ from centroid due to half-cell resistance"
        );
    }
}
