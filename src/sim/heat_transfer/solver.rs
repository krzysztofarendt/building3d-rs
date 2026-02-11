use crate::sim::heat_transfer::boundary::BoundaryCondition;
use crate::sim::heat_transfer::mesh::{BOUNDARY, FvmMesh};

/// 1D FVM wall solver implementing implicit (backward Euler) time stepping.
///
/// Each exterior polygon with a `WallConstruction` gets its own solver instance.
/// The solver assembles a tridiagonal system and solves it with the Thomas
/// algorithm every time step.
pub struct FvmWallSolver {
    mesh: FvmMesh,
    /// Current cell temperatures [C].
    temperatures: Vec<f64>,
    /// Wall polygon area [m^2], needed for BC flux scaling.
    wall_area: f64,
}

impl Clone for FvmWallSolver {
    fn clone(&self) -> Self {
        Self {
            mesh: self.mesh.clone(),
            temperatures: self.temperatures.clone(),
            wall_area: self.wall_area,
        }
    }
}

impl FvmWallSolver {
    /// Create a new solver from a mesh and an initial uniform temperature.
    pub fn new(mesh: FvmMesh, initial_temperature: f64) -> Self {
        let n = mesh.cells.len();
        Self {
            wall_area: if n > 0 {
                // Recover area from any boundary face
                mesh.faces[0].area
            } else {
                1.0
            },
            temperatures: vec![initial_temperature; n],
            mesh,
        }
    }

    /// Advance one time step.
    ///
    /// `dt` is the time step in seconds.
    /// `bc_exterior` applies to the left (exterior) boundary.
    /// `bc_interior` applies to the right (interior) boundary.
    /// `sources` is a per-cell volumetric heat source in watts.
    ///
    /// Returns the updated cell temperatures.
    pub fn step(
        &mut self,
        dt: f64,
        bc_exterior: &BoundaryCondition,
        bc_interior: &BoundaryCondition,
        sources: &[f64],
    ) -> &[f64] {
        let n = self.mesh.cells.len();
        if n == 0 {
            return &self.temperatures;
        }

        // Tridiagonal coefficients: a[i]*T[i-1] + b[i]*T[i] + c[i]*T[i+1] = d[i]
        let mut a = vec![0.0; n]; // sub-diagonal
        let mut b = vec![0.0; n]; // diagonal
        let mut c = vec![0.0; n]; // super-diagonal
        let mut d = vec![0.0; n]; // RHS

        // Capacity contribution: C_i / dt
        for i in 0..n {
            let cap = self.mesh.cells[i].capacity();
            b[i] += cap / dt;
            d[i] += (cap / dt) * self.temperatures[i];
        }

        // Source terms
        for i in 0..n {
            if i < sources.len() {
                d[i] += sources[i];
            }
        }

        // Interior face conductances (between adjacent cells)
        for face in &self.mesh.faces {
            if face.cell_left == BOUNDARY || face.cell_right == BOUNDARY {
                continue;
            }
            let il = face.cell_left;
            let ir = face.cell_right;
            let k = face.conductance;

            // Adds -K to off-diagonal, +K to diagonal (implicit coupling)
            b[il] += k;
            b[ir] += k;
            // il equation: ... - K * T[ir]
            // ir equation: ... - K * T[il]
            if ir == il + 1 {
                c[il] -= k;
                a[ir] -= k;
            }
        }

        // Exterior boundary condition (face with cell_left == BOUNDARY)
        if let Some(fi) = self.mesh.exterior_boundary_face() {
            let face = &self.mesh.faces[fi];
            let cell_idx = face.cell_right;
            apply_bc(
                bc_exterior,
                cell_idx,
                face.conductance,
                self.wall_area,
                &mut b,
                &mut d,
            );
        }

        // Interior boundary condition (face with cell_right == BOUNDARY)
        if let Some(fi) = self.mesh.interior_boundary_face() {
            let face = &self.mesh.faces[fi];
            let cell_idx = face.cell_left;
            apply_bc(
                bc_interior,
                cell_idx,
                face.conductance,
                self.wall_area,
                &mut b,
                &mut d,
            );
        }

        // Solve tridiagonal system with Thomas algorithm
        thomas_solve(&a, &mut b, &c, &mut d);
        self.temperatures.copy_from_slice(&d);

        &self.temperatures
    }

    /// Temperature of the exterior-most cell (first cell).
    pub fn exterior_surface_temp(&self) -> f64 {
        self.temperatures.first().copied().unwrap_or(0.0)
    }

    /// Temperature of the interior-most cell (last cell).
    pub fn interior_surface_temp(&self) -> f64 {
        self.temperatures.last().copied().unwrap_or(0.0)
    }

    /// Heat flux through the interior surface into the zone [W/m^2].
    ///
    /// Positive means heat flows from wall into the zone.
    pub fn interior_heat_flux(&self, bc_interior: &BoundaryCondition) -> f64 {
        let n = self.mesh.cells.len();
        if n == 0 {
            return 0.0;
        }
        let t_centroid = self.temperatures[n - 1];
        self.boundary_flux(t_centroid, bc_interior, self.mesh.interior_boundary_face())
    }

    /// Heat flux through the exterior surface [W/m^2].
    ///
    /// Positive means heat flows from the wall to the outside.
    pub fn exterior_heat_flux(&self, bc_exterior: &BoundaryCondition) -> f64 {
        if self.mesh.cells.is_empty() {
            return 0.0;
        }
        let t_centroid = self.temperatures[0];
        self.boundary_flux(t_centroid, bc_exterior, self.mesh.exterior_boundary_face())
    }

    /// Conductance (W/K) of the exterior boundary face (includes wall area).
    pub fn exterior_boundary_face_conductance_w_per_k(&self) -> Option<f64> {
        self.mesh
            .exterior_boundary_face()
            .map(|fi| self.mesh.faces[fi].conductance)
    }

    /// Conductance (W/K) of the interior boundary face (includes wall area).
    pub fn interior_boundary_face_conductance_w_per_k(&self) -> Option<f64> {
        self.mesh
            .interior_boundary_face()
            .map(|fi| self.mesh.faces[fi].conductance)
    }

    /// Fraction (0..1) of a surface source that enters the wall domain when the boundary is
    /// modeled as convection with film coefficient `h` (W/(m²·K)).
    ///
    /// This mirrors the split used by [`BoundaryCondition::ConvectiveWithFlux`] in `apply_bc()`.
    pub fn boundary_conduction_fraction(&self, h: f64, is_interior: bool) -> f64 {
        let face_conductance = if is_interior {
            self.interior_boundary_face_conductance_w_per_k()
        } else {
            self.exterior_boundary_face_conductance_w_per_k()
        };
        let Some(k_face) = face_conductance else {
            return 0.0;
        };
        let h_a = h.max(0.0) * self.wall_area;
        if h_a <= 0.0 {
            return 1.0;
        }
        (k_face / (k_face + h_a)).clamp(0.0, 1.0)
    }

    /// Compute heat flux at a boundary surface.
    ///
    /// For convective BCs, the actual surface temperature is computed from the
    /// cell centroid temperature by accounting for the half-cell conduction
    /// resistance between centroid and wall surface.
    fn boundary_flux(
        &self,
        t_centroid: f64,
        bc: &BoundaryCondition,
        face_idx: Option<usize>,
    ) -> f64 {
        match bc {
            BoundaryCondition::Convective { h, t_fluid } => {
                if *h <= 0.0 {
                    return 0.0;
                }
                // Series resistance: half-cell conduction + convective film.
                // q = (T_centroid - T_fluid) / (half_dx/k + 1/h)
                if let Some(fi) = face_idx {
                    let k_cond = self.mesh.faces[fi].conductance / self.wall_area;
                    let k_combined = 1.0 / (1.0 / k_cond + 1.0 / h);
                    k_combined * (t_centroid - t_fluid)
                } else {
                    h * (t_centroid - t_fluid)
                }
            }
            BoundaryCondition::ConvectiveWithFlux {
                h,
                t_fluid,
                heat_flux,
            }
            | BoundaryCondition::ConvectiveWithFluxToDomain {
                h,
                t_fluid,
                heat_flux,
            } => {
                // Report only the convective component (wall ↔ fluid heat transfer).
                //
                // The imposed surface source affects the wall temperatures and thus the
                // convective flux indirectly, but it is not itself a wall→fluid heat flux term.
                let _ = heat_flux;
                if *h <= 0.0 {
                    return 0.0;
                }
                if let Some(fi) = face_idx {
                    let k_cond = self.mesh.faces[fi].conductance / self.wall_area;
                    let k_combined = 1.0 / (1.0 / k_cond + 1.0 / h);
                    k_combined * (t_centroid - t_fluid)
                } else {
                    h * (t_centroid - t_fluid)
                }
            }
            BoundaryCondition::Dirichlet { temperature } => {
                if let Some(fi) = face_idx {
                    let k_per_area = self.mesh.faces[fi].conductance / self.wall_area;
                    k_per_area * (t_centroid - temperature)
                } else {
                    0.0
                }
            }
            BoundaryCondition::Neumann { heat_flux } => {
                // Neumann BC is defined as positive into the wall domain.
                // Reported heat fluxes are positive out of the wall.
                -*heat_flux
            }
        }
    }

    /// Access current cell temperatures.
    pub fn temperatures(&self) -> &[f64] {
        &self.temperatures
    }

    /// Total thermal capacity of the wall mesh (sum of cell capacities) [J/K].
    pub fn total_capacity_j_per_k(&self) -> f64 {
        self.mesh.cells.iter().map(|c| c.capacity()).sum()
    }

    /// Estimate the actual boundary **surface** temperature (not the cell centroid).
    ///
    /// For convective boundary conditions, the solver eliminates the surface node and
    /// applies an effective conductance between the fluid and the boundary-adjacent cell.
    /// This helper reconstructs the corresponding surface temperature from the current
    /// boundary-adjacent cell temperature, the half-cell conductance, and any imposed
    /// surface source term.
    fn boundary_surface_temperature(
        &self,
        t_centroid: f64,
        bc: &BoundaryCondition,
        face_idx: Option<usize>,
        is_interior: bool,
    ) -> f64 {
        let Some(fi) = face_idx else {
            return t_centroid;
        };
        if self.wall_area <= 0.0 {
            return t_centroid;
        }
        let k_face_per_area = (self.mesh.faces[fi].conductance / self.wall_area).max(0.0);
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
                let alpha = self.boundary_conduction_fraction(h, is_interior);
                let q_into_wall = alpha * heat_flux;
                (k_face_per_area * t_centroid + h * t_fluid + q_into_wall) / (k_face_per_area + h)
            }
            BoundaryCondition::ConvectiveWithFluxToDomain {
                h,
                t_fluid,
                heat_flux,
            } => {
                if h <= 0.0 {
                    return t_centroid;
                }
                // `ConvectiveWithFluxToDomain` represents a source injected into the wall domain
                // (added to the boundary-adjacent control volume), not a surface energy-balance
                // term. The effect of `heat_flux` is already reflected in `t_centroid`, so we do
                // not add it again when reconstructing the surface temperature.
                let _ = heat_flux;
                (k_face_per_area * t_centroid + h * t_fluid) / (k_face_per_area + h)
            }
        }
    }

    /// Estimate the interior boundary **surface** temperature (°C).
    pub fn interior_surface_temperature(&self, bc_interior: &BoundaryCondition) -> f64 {
        let n = self.mesh.cells.len();
        if n == 0 {
            return 0.0;
        }
        let t_centroid = self.temperatures[n - 1];
        self.boundary_surface_temperature(
            t_centroid,
            bc_interior,
            self.mesh.interior_boundary_face(),
            true,
        )
    }

    /// Estimate the exterior boundary **surface** temperature (°C).
    pub fn exterior_surface_temperature(&self, bc_exterior: &BoundaryCondition) -> f64 {
        if self.mesh.cells.is_empty() {
            return 0.0;
        }
        let t_centroid = self.temperatures[0];
        self.boundary_surface_temperature(
            t_centroid,
            bc_exterior,
            self.mesh.exterior_boundary_face(),
            false,
        )
    }
}

/// Apply a boundary condition to the tridiagonal system.
///
/// `cell_idx` is the cell adjacent to the boundary.
/// `face_conductance` is the precomputed k*A/d for the boundary face.
fn apply_bc(
    bc: &BoundaryCondition,
    cell_idx: usize,
    face_conductance: f64,
    wall_area: f64,
    diag: &mut [f64],
    rhs: &mut [f64],
) {
    match bc {
        BoundaryCondition::Dirichlet { temperature } => {
            // The boundary face has conductance K = k*A/half_dx.
            // Contribution: K * (T_boundary - T_cell) added implicitly.
            diag[cell_idx] += face_conductance;
            rhs[cell_idx] += face_conductance * temperature;
        }
        BoundaryCondition::Neumann { heat_flux } => {
            // Prescribed flux: Q = heat_flux * A  (positive into domain)
            rhs[cell_idx] += heat_flux * wall_area;
        }
        BoundaryCondition::Convective { h, t_fluid } => {
            // The path from fluid to cell centroid has two resistances in series:
            //   1/(h*A)        — convective film
            //   1/K_face       — half-cell conduction (K_face = k*A/half_dx)
            // Effective conductance: K_eff = 1 / (1/(h*A) + 1/K_face)
            let h_a = h * wall_area;
            if h_a <= 0.0 {
                return;
            }
            let k_eff = if face_conductance > 0.0 {
                1.0 / (1.0 / h_a + 1.0 / face_conductance)
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
            let h_a = h * wall_area;
            if h_a <= 0.0 {
                // Pure imposed flux into the domain.
                rhs[cell_idx] += heat_flux * wall_area;
                return;
            }
            let k_eff = if face_conductance > 0.0 {
                1.0 / (1.0 / h_a + 1.0 / face_conductance)
            } else {
                h_a
            };
            diag[cell_idx] += k_eff;
            rhs[cell_idx] += k_eff * t_fluid;
            // Split the imposed surface flux between convection to the fluid and conduction into
            // the wall. For a boundary face with half-cell conductance `K_face`, the fraction that
            // enters the wall control volume is:
            //   alpha = K_face / (K_face + h*A)
            // so that the added RHS term is `alpha * q_solar * A`.
            //
            // This prevents incorrectly injecting 100% of absorbed shortwave into the wall when
            // there is a strong convective loss path to outdoors.
            let alpha = if face_conductance > 0.0 {
                face_conductance / (face_conductance + h_a)
            } else {
                0.0
            };
            rhs[cell_idx] += alpha * heat_flux * wall_area;
        }
        BoundaryCondition::ConvectiveWithFluxToDomain {
            h,
            t_fluid,
            heat_flux,
        } => {
            let h_a = h * wall_area;
            if h_a <= 0.0 {
                rhs[cell_idx] += heat_flux * wall_area;
                return;
            }
            let k_eff = if face_conductance > 0.0 {
                1.0 / (1.0 / h_a + 1.0 / face_conductance)
            } else {
                h_a
            };
            diag[cell_idx] += k_eff;
            rhs[cell_idx] += k_eff * t_fluid;
            // Apply 100% of the imposed surface source into the domain.
            rhs[cell_idx] += heat_flux * wall_area;
        }
    }
}

/// Thomas algorithm for solving a tridiagonal system.
///
/// On entry:
/// - `a[i]` is the sub-diagonal (a[0] unused)
/// - `b[i]` is the diagonal
/// - `c[i]` is the super-diagonal (c[n-1] unused)
/// - `d[i]` is the RHS
///
/// On exit `d` contains the solution.
fn thomas_solve(a: &[f64], b: &mut [f64], c: &[f64], d: &mut [f64]) {
    let n = b.len();
    if n == 0 {
        return;
    }

    // Forward sweep
    for i in 1..n {
        if b[i - 1].abs() < 1e-30 {
            continue;
        }
        let w = a[i] / b[i - 1];
        b[i] -= w * c[i - 1];
        d[i] -= w * d[i - 1];
    }

    // Back substitution
    d[n - 1] /= b[n - 1];
    for i in (0..n - 1).rev() {
        d[i] = (d[i] - c[i] * d[i + 1]) / b[i];
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::energy::construction::WallConstruction;
    use crate::sim::heat_transfer::mesh_1d::build_1d_mesh;
    use crate::sim::materials::Layer;

    /// Steady-state single-layer wall with Dirichlet BCs.
    /// T_ext = 0, T_int = 20 => linear profile, q = k * dT / L.
    #[test]
    fn test_steady_state_dirichlet() {
        let k = 1.4;
        let thickness = 0.20;
        let area = 1.0;
        let construction = WallConstruction::new(
            "test",
            vec![Layer {
                name: "concrete".into(),
                thickness,
                conductivity: k,
                density: 2300.0,
                specific_heat: 880.0,
            }],
        );
        let mesh = build_1d_mesh(&construction, area);
        let mut solver = FvmWallSolver::new(mesh, 10.0);

        let bc_ext = BoundaryCondition::Dirichlet { temperature: 0.0 };
        let bc_int = BoundaryCondition::Dirichlet { temperature: 20.0 };
        let sources = vec![0.0; solver.temperatures().len()];

        // Run to steady state (large dt effectively makes capacity negligible)
        let dt = 1e8;
        for _ in 0..10 {
            solver.step(dt, &bc_ext, &bc_int, &sources);
        }

        let temps = solver.temperatures();
        let n = temps.len();

        // Check linear profile
        for i in 0..n {
            let x_frac = (i as f64 + 0.5) / n as f64;
            let expected = 0.0 + 20.0 * x_frac;
            assert!(
                (temps[i] - expected).abs() < 0.01,
                "cell {i}: got {}, expected {}",
                temps[i],
                expected
            );
        }

        // Heat flux into zone: negative because heat flows from zone (warm) to
        // exterior (cold).  Magnitude = k * dT / L = 1.4 * 20 / 0.20 = 140 W/m^2.
        let expected_flux = -(k * 20.0 / thickness);
        let flux = solver.interior_heat_flux(&bc_int);
        assert!(
            (flux - expected_flux).abs() < 1.0,
            "flux = {flux}, expected {expected_flux}"
        );
    }

    /// Verify that the Thomas algorithm solves a simple 3x3 system correctly.
    #[test]
    fn test_thomas_algorithm() {
        // System: [2 -1 0; -1 2 -1; 0 -1 2] * x = [1; 0; 1]
        // Solution: x = [1; 1; 1]
        let a = vec![0.0, -1.0, -1.0];
        let mut b = vec![2.0, 2.0, 2.0];
        let c = vec![-1.0, -1.0, 0.0];
        let mut d = vec![1.0, 0.0, 1.0];

        thomas_solve(&a, &mut b, &c, &mut d);

        assert!((d[0] - 1.0).abs() < 1e-10);
        assert!((d[1] - 1.0).abs() < 1e-10);
        assert!((d[2] - 1.0).abs() < 1e-10);
    }

    /// Multi-layer wall with Dirichlet BCs at steady state.
    /// Verify total heat flux matches U*dT.
    #[test]
    fn test_steady_state_multilayer() {
        let construction = WallConstruction::new(
            "insulated",
            vec![
                Layer {
                    name: "plaster_ext".into(),
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
        let area = 1.0;
        let mesh = build_1d_mesh(&construction, area);
        let mut solver = FvmWallSolver::new(mesh, 10.0);

        let t_ext = -5.0;
        let t_int = 20.0;
        let bc_ext = BoundaryCondition::Dirichlet { temperature: t_ext };
        let bc_int = BoundaryCondition::Dirichlet { temperature: t_int };
        let sources = vec![0.0; solver.temperatures().len()];

        // Run to steady state
        let dt = 1e8;
        for _ in 0..20 {
            solver.step(dt, &bc_ext, &bc_int, &sources);
        }

        // Total R = 0.02/0.87 + 0.10/0.04 + 0.15/1.4
        // Flux is negative (heat flows from zone to exterior).
        let r_total = 0.02 / 0.87 + 0.10 / 0.04 + 0.15 / 1.4;
        let expected_flux = -(t_int - t_ext) / r_total;
        let flux = solver.interior_heat_flux(&bc_int);

        assert!(
            (flux - expected_flux).abs() < 0.5,
            "flux = {flux}, expected {expected_flux} (R = {r_total})"
        );
    }

    /// Convective BCs on both sides, steady-state.
    /// q = (T_out - T_in) / (1/h_out + R_wall + 1/h_in)
    #[test]
    fn test_steady_state_convective() {
        let construction = WallConstruction::new(
            "concrete",
            vec![Layer {
                name: "concrete".into(),
                thickness: 0.20,
                conductivity: 1.4,
                density: 2300.0,
                specific_heat: 880.0,
            }],
        );
        let area = 1.0;
        let mesh = build_1d_mesh(&construction, area);
        let mut solver = FvmWallSolver::new(mesh, 10.0);

        let h_out = 25.0;
        let h_in = 7.7;
        let t_out = -10.0;
        let t_in = 20.0;

        let bc_ext = BoundaryCondition::Convective {
            h: h_out,
            t_fluid: t_out,
        };
        let bc_int = BoundaryCondition::Convective {
            h: h_in,
            t_fluid: t_in,
        };
        let sources = vec![0.0; solver.temperatures().len()];

        // Run to steady state
        let dt = 1e8;
        for _ in 0..20 {
            solver.step(dt, &bc_ext, &bc_int, &sources);
        }

        // R_total = 1/h_out + L/k + 1/h_in
        let r_total = 1.0 / h_out + 0.20 / 1.4 + 1.0 / h_in;
        let expected_flux = (t_out - t_in) / r_total;
        let flux = solver.interior_heat_flux(&bc_int);

        // Allow ~1% error from FVM spatial discretization
        assert!(
            (flux - expected_flux).abs() < expected_flux.abs() * 0.01,
            "flux = {flux}, expected {expected_flux} (R = {r_total})"
        );
    }

    /// Multi-layer wall with convective BCs: verify interface temperatures.
    #[test]
    fn test_convective_multilayer_interface_temps() {
        let construction = WallConstruction::new(
            "insulated",
            vec![
                Layer {
                    name: "plaster_ext".into(),
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
        let area = 1.0;
        let h_out = 25.0;
        let h_in = 7.7;
        let t_out = -10.0;
        let t_in = 20.0;

        let mesh = build_1d_mesh(&construction, area);
        let mut solver = FvmWallSolver::new(mesh, 10.0);

        let bc_ext = BoundaryCondition::Convective {
            h: h_out,
            t_fluid: t_out,
        };
        let bc_int = BoundaryCondition::Convective {
            h: h_in,
            t_fluid: t_in,
        };
        let sources = vec![0.0; solver.temperatures().len()];

        let dt = 1e8;
        for _ in 0..20 {
            solver.step(dt, &bc_ext, &bc_int, &sources);
        }

        // R_total = 1/h_out + 0.02/0.87 + 0.10/0.04 + 0.15/1.4 + 1/h_in
        let r_total = 1.0 / h_out + 0.02 / 0.87 + 0.10 / 0.04 + 0.15 / 1.4 + 1.0 / h_in;
        let q = (t_out - t_in) / r_total; // W/m^2 (negative)

        // Interior surface temperature: T_in - q / h_in
        // (q is negative, so T_surface < T_in)
        let t_int_surface = t_in + q / h_in;

        // Exterior surface temperature: T_out - q / h_out
        let t_ext_surface = t_out - q / h_out;

        let actual_t_int = solver.interior_surface_temp();
        let actual_t_ext = solver.exterior_surface_temp();

        // Cell temps are at centroids, not surfaces, so allow ~1 degree tolerance
        // for the half-cell offset.
        assert!(
            (actual_t_int - t_int_surface).abs() < 1.5,
            "interior: got {actual_t_int}, expected ~{t_int_surface}"
        );
        assert!(
            (actual_t_ext - t_ext_surface).abs() < 1.5,
            "exterior: got {actual_t_ext}, expected ~{t_ext_surface}"
        );
    }

    /// Transient step response: semi-infinite solid with sudden surface
    /// temperature change.  Compare to erfc analytical solution at selected
    /// points.
    #[test]
    fn test_transient_step_response() {
        // Thick concrete slab, only test interior cells (away from far boundary)
        let k = 1.4;
        let rho = 2300.0;
        let cp = 880.0;
        let alpha = k / (rho * cp); // thermal diffusivity m^2/s
        let thickness = 1.0; // 1m thick slab
        let area = 1.0;

        let construction = WallConstruction::new(
            "thick_concrete",
            vec![Layer {
                name: "concrete".into(),
                thickness,
                conductivity: k,
                density: rho,
                specific_heat: cp,
            }],
        );
        let mesh = build_1d_mesh(&construction, area);
        // 1.0m / 0.05 = 20 cells
        let n = mesh.cells.len();
        assert_eq!(n, 20);

        let t_initial = 20.0;
        let t_surface = 0.0; // sudden step at exterior surface
        let mut solver = FvmWallSolver::new(mesh, t_initial);

        let bc_ext = BoundaryCondition::Dirichlet {
            temperature: t_surface,
        };
        // Far boundary: adiabatic (Neumann zero flux)
        let bc_int = BoundaryCondition::Neumann { heat_flux: 0.0 };
        let sources = vec![0.0; n];

        // Simulate 1 hour with 60s time steps
        let dt = 60.0;
        let total_time = 3600.0;
        let steps = (total_time / dt) as usize;
        for _ in 0..steps {
            solver.step(dt, &bc_ext, &bc_int, &sources);
        }

        // Analytical solution: T(x,t) = T_initial + (T_surface - T_initial) * erfc(x / (2*sqrt(alpha*t)))
        // At t = 3600s, sqrt(alpha*t) = sqrt(6.917e-7 * 3600) = sqrt(0.002490) = 0.04990
        let sqrt_alpha_t = (alpha * total_time).sqrt();

        // Check a few cells near the exterior surface
        let temps = solver.temperatures();
        for i in 0..5 {
            let dx = thickness / n as f64;
            let x = (i as f64 + 0.5) * dx; // cell centroid
            let xi = x / (2.0 * sqrt_alpha_t);
            let expected = t_initial + (t_surface - t_initial) * erfc(xi);

            // Allow 15% relative error or 0.5 degree absolute (FVM discretization error)
            let abs_err = (temps[i] - expected).abs();
            let rel_err = abs_err / (t_initial - t_surface).abs();
            assert!(
                abs_err < 0.5 || rel_err < 0.15,
                "cell {i} (x={x:.3}m): got {:.2}, expected {:.2} (err={abs_err:.2}, rel={rel_err:.3})",
                temps[i],
                expected,
            );
        }
    }

    #[test]
    fn test_neumann_flux_reporting_sign() {
        let construction = WallConstruction::new(
            "test",
            vec![Layer {
                name: "concrete".into(),
                thickness: 0.05,
                conductivity: 1.4,
                density: 2300.0,
                specific_heat: 880.0,
            }],
        );
        let mesh = build_1d_mesh(&construction, 1.0);
        let solver = FvmWallSolver::new(mesh, 20.0);

        // Neumann is defined as positive into the wall (domain).
        // Reported fluxes are positive out of the wall.
        let q_in = 12.3;
        let bc = BoundaryCondition::Neumann { heat_flux: q_in };

        assert!((solver.interior_heat_flux(&bc) + q_in).abs() < 1e-12);
        assert!((solver.exterior_heat_flux(&bc) + q_in).abs() < 1e-12);
    }

    /// Complementary error function approximation (Abramowitz & Stegun 7.1.26).
    fn erfc(x: f64) -> f64 {
        if x < 0.0 {
            return 2.0 - erfc(-x);
        }
        let t = 1.0 / (1.0 + 0.3275911 * x);
        let poly = t
            * (0.254829592
                + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
        poly * (-x * x).exp()
    }
}
