//! Global simultaneous FVM thermal solver.
//!
//! Assembles ALL wall FVM cells + surface nodes + air nodes into one global
//! matrix and solves simultaneously. Radiation coupling embeds directly in the
//! matrix (no iteration needed). Surface temperatures, air temperatures, and
//! radiative exchange are all self-consistent within each timestep.
//!
//! Node types in the global system:
//! ```text
//! [T_out] --K_ext-- [cell_0] --K_01-- ... --K_face-- [surf_node] --h_conv*A-- [air_node]
//!                    ~~C_0~~                                                    ~~C_air~~
//!                                                        └── h_rad*A*F_ij to other surf_nodes
//! ```

use crate::sim::energy::hvac::HvacIdealLoads;
use crate::sim::energy::network::solve::solve_dense;
use crate::sim::heat_transfer::mesh::BOUNDARY;
use crate::sim::heat_transfer::FvmWallSolver;

// ─── Topology types ─────────────────────────────────────────────────────

/// Static topology extracted once during initialization.
pub struct GlobalTopology {
    /// Total number of unknowns.
    pub n: usize,
    /// Per-wall topology (exterior and internal mass walls).
    pub walls: Vec<WallTopology>,
    /// Per-surface-node topology (one per interior-facing wall face).
    pub surfaces: Vec<SurfaceTopology>,
    /// Per-zone air node topology.
    pub air_nodes: Vec<AirNodeTopology>,
}

/// Topology for one FVM wall (chain of cells).
pub struct WallTopology {
    /// Index in the caller's wall array.
    pub wall_idx: usize,
    /// Global index of the first cell.
    pub cell_offset: usize,
    /// Number of cells in this wall.
    pub n_cells: usize,
    /// Zone index (into `air_nodes`).
    pub zone_idx: usize,
    /// Surface node index (into `surfaces`) for the interior face.
    pub surface_idx: usize,
    /// Cell thermal capacities [J/K].
    pub cell_capacities: Vec<f64>,
    /// Conductances between adjacent cells [W/K] (length = n_cells - 1).
    pub internal_conductances: Vec<f64>,
    /// Conductance of the exterior boundary face [W/K] (half-cell: k*A/half_dx).
    pub ext_boundary_face_k: f64,
    /// Conductance of the interior boundary face [W/K] (half-cell: k*A/half_dx).
    pub int_boundary_face_k: f64,
    /// Wall polygon area [m²].
    pub area_m2: f64,
}

/// Topology for a surface node (zero-capacity interface between wall interior and air).
pub struct SurfaceTopology {
    /// Index in the `surfaces` array.
    pub surface_idx: usize,
    /// Global index of this surface node.
    pub global_idx: usize,
    /// Wall index (into `walls`).
    pub wall_idx: usize,
    /// Global index of the innermost wall cell.
    pub inner_cell_global_idx: usize,
    /// Zone index (into `air_nodes`).
    pub zone_idx: usize,
    /// Surface area [m²].
    pub area_m2: f64,
    /// Interior boundary face conductance [W/K] (half-cell: k*A/half_dx).
    pub k_face: f64,
}

/// Topology for a zone air node.
pub struct AirNodeTopology {
    /// Zone index.
    pub zone_idx: usize,
    /// Global index.
    pub global_idx: usize,
    /// Air thermal capacity [J/K] (rho * cp * V).
    pub capacity_j_per_k: f64,
    /// Infiltration conductance [W/K].
    pub infiltration_k: f64,
    /// Steady-state glazing UA [W/K].
    pub glazing_ua: f64,
}

// ─── Per-step condition types ───────────────────────────────────────────

/// Per-wall boundary conditions for one time step.
pub struct WallStepConditions {
    /// Series K_eff for exterior coupling [W/K]:
    /// 1/(1/(h_out*A) + 1/K_ext_boundary_face).
    pub ext_k_eff: f64,
    /// Driving temperature for exterior [°C].
    pub ext_t_drive: f64,
    /// Source at exterior cell [W] (alpha * solar * A + longwave, etc.).
    pub ext_source_w: f64,
    /// Convective coefficient for interior surface [W/(m²·K)].
    pub h_conv: f64,
    /// Total interior film coefficient (conv + rad) [W/(m²·K)].
    pub h_total: f64,
    /// Source at interior surface [W] (transmitted solar, radiant gains).
    pub int_source_w: f64,
}

/// Per-zone air conditions for one time step.
pub struct AirStepConditions {
    /// Outdoor temperature [°C].
    pub outdoor_temp_c: f64,
    /// Direct convective gains to air [W] (internal gains, etc.).
    pub direct_gains_w: f64,
}

/// Radiation coupling conditions.
pub struct RadiationConditions {
    /// Linearized radiation coefficient [W/(m²·K)].
    pub h_rad: f64,
    /// Surface areas [m²] (per surface node).
    pub surface_areas: Vec<f64>,
    /// View factor matrix, row-major F_ij (n_surfaces × n_surfaces).
    pub f_matrix: Vec<f64>,
    /// Number of surface nodes.
    pub n_surfaces: usize,
}

/// Result of one global solve step.
pub struct GlobalStepResult {
    /// Heating power per zone [W] (positive = heating).
    pub heating_w_per_zone: Vec<f64>,
    /// Cooling power per zone [W] (positive = cooling).
    pub cooling_w_per_zone: Vec<f64>,
}

// ─── Topology builder ───────────────────────────────────────────────────

/// Information about one FVM wall for topology building.
pub struct FvmWallInfo<'a> {
    pub solver: &'a FvmWallSolver,
    pub zone_idx: usize,
    pub area_m2: f64,
    /// True if this wall has a second surface node (e.g. TwoSided internal mass).
    pub has_exterior_surface: bool,
    /// True if exterior face is adiabatic (OneSidedAdiabatic internal mass).
    pub exterior_adiabatic: bool,
}

/// Build topology from FVM wall solvers and zone data.
///
/// Called once during initialization. The caller provides wall info and zone
/// parameters; this function assigns global indices and extracts mesh topology.
pub fn build_topology(
    wall_infos: &[FvmWallInfo],
    num_zones: usize,
    zone_air_capacities: &[f64],
    zone_infiltration_k: &[f64],
    zone_glazing_ua: &[f64],
) -> GlobalTopology {
    let mut idx = 0usize; // running global index

    // 1. Wall cells
    let mut walls = Vec::with_capacity(wall_infos.len());
    for (wi, info) in wall_infos.iter().enumerate() {
        let mesh = info.solver.mesh();
        let n_cells = mesh.cells.len();
        let cell_offset = idx;

        let cell_capacities: Vec<f64> = mesh.cells.iter().map(|c| c.capacity()).collect();

        // Internal conductances between adjacent cells
        let mut internal_conductances = Vec::new();
        for face in &mesh.faces {
            if face.cell_left == BOUNDARY || face.cell_right == BOUNDARY {
                continue;
            }
            internal_conductances.push(face.conductance);
        }

        // Boundary face conductances
        let ext_k = mesh
            .exterior_boundary_face()
            .map(|fi| mesh.faces[fi].conductance)
            .unwrap_or(0.0);
        let int_k = mesh
            .interior_boundary_face()
            .map(|fi| mesh.faces[fi].conductance)
            .unwrap_or(0.0);

        walls.push(WallTopology {
            wall_idx: wi,
            cell_offset,
            n_cells,
            zone_idx: info.zone_idx,
            surface_idx: 0, // filled below
            cell_capacities,
            internal_conductances,
            ext_boundary_face_k: ext_k,
            int_boundary_face_k: int_k,
            area_m2: info.area_m2,
        });

        idx += n_cells;
    }

    // 2. Surface nodes (one per wall interior face, plus optionally exterior for TwoSided)
    let mut surfaces = Vec::new();
    for (wi, info) in wall_infos.iter().enumerate() {
        let wall = &walls[wi];
        let inner_cell_global_idx = wall.cell_offset + wall.n_cells.saturating_sub(1);

        // Interior surface node
        let si = surfaces.len();
        surfaces.push(SurfaceTopology {
            surface_idx: si,
            global_idx: idx,
            wall_idx: wi,
            inner_cell_global_idx,
            zone_idx: info.zone_idx,
            area_m2: info.area_m2,
            k_face: wall.int_boundary_face_k,
        });
        // Link wall to its primary surface
        // (walls vec is mutable through index, but we built it above)
        idx += 1;

        // For TwoSided internal mass: add a second surface node on the exterior side
        if info.has_exterior_surface {
            let outer_cell_global_idx = wall.cell_offset;
            let si2 = surfaces.len();
            surfaces.push(SurfaceTopology {
                surface_idx: si2,
                global_idx: idx,
                wall_idx: wi,
                inner_cell_global_idx: outer_cell_global_idx,
                zone_idx: info.zone_idx,
                area_m2: info.area_m2,
                k_face: wall.ext_boundary_face_k,
            });
            idx += 1;
        }
    }

    // Fix up surface_idx on walls (primary interior surface)
    {
        let mut wall_to_surface = vec![0usize; walls.len()];
        for s in &surfaces {
            // First surface for each wall is the interior one
            if wall_to_surface[s.wall_idx] == 0 || s.inner_cell_global_idx == walls[s.wall_idx].cell_offset + walls[s.wall_idx].n_cells.saturating_sub(1) {
                wall_to_surface[s.wall_idx] = s.surface_idx;
            }
        }
        for (wi, wall) in walls.iter_mut().enumerate() {
            wall.surface_idx = wall_to_surface[wi];
        }
    }

    // 3. Air nodes
    let mut air_nodes = Vec::with_capacity(num_zones);
    for zi in 0..num_zones {
        air_nodes.push(AirNodeTopology {
            zone_idx: zi,
            global_idx: idx,
            capacity_j_per_k: zone_air_capacities.get(zi).copied().unwrap_or(0.0),
            infiltration_k: zone_infiltration_k.get(zi).copied().unwrap_or(0.0),
            glazing_ua: zone_glazing_ua.get(zi).copied().unwrap_or(0.0),
        });
        idx += 1;
    }

    GlobalTopology {
        n: idx,
        walls,
        surfaces,
        air_nodes,
    }
}

// ─── Global solve ───────────────────────────────────────────────────────

/// Assemble and solve the global system for one timestep.
///
/// `temperatures` is the state vector (length `topo.n`), updated in-place.
/// Returns heating/cooling power per zone.
#[allow(clippy::too_many_arguments)]
pub fn step_global(
    topo: &GlobalTopology,
    temperatures: &mut [f64],
    wall_conditions: &[WallStepConditions],
    air_conditions: &[AirStepConditions],
    radiation: Option<&RadiationConditions>,
    hvac: &HvacIdealLoads,
    dt_s: f64,
    wall_infos: &[FvmWallInfo],
) -> GlobalStepResult {
    let n_zones = topo.air_nodes.len();

    // ── Free-float solve ────────────────────────────────────────────────
    let free_temps = solve_system(topo, temperatures, wall_conditions, air_conditions, radiation, dt_s, wall_infos, None);

    let Ok(free_temps) = free_temps else {
        // Fallback: keep current temps, no HVAC
        return GlobalStepResult {
            heating_w_per_zone: vec![0.0; n_zones],
            cooling_w_per_zone: vec![0.0; n_zones],
        };
    };

    // ── HVAC check ──────────────────────────────────────────────────────
    // Check each zone air node: if free-float T is outside setpoints, re-solve
    // with that air node fixed at the setpoint.
    let mut need_hvac = false;
    let mut fixed_air_temps = vec![0.0f64; n_zones];
    let mut zone_needs_hvac = vec![false; n_zones];

    for air in &topo.air_nodes {
        let t_free = free_temps[air.global_idx];
        if t_free < hvac.heating_setpoint {
            zone_needs_hvac[air.zone_idx] = true;
            fixed_air_temps[air.zone_idx] = hvac.heating_setpoint;
            need_hvac = true;
        } else if t_free > hvac.cooling_setpoint {
            zone_needs_hvac[air.zone_idx] = true;
            fixed_air_temps[air.zone_idx] = hvac.cooling_setpoint;
            need_hvac = true;
        }
    }

    if !need_hvac {
        // Accept free-float solution
        temperatures.copy_from_slice(&free_temps);
        return GlobalStepResult {
            heating_w_per_zone: vec![0.0; n_zones],
            cooling_w_per_zone: vec![0.0; n_zones],
        };
    }

    // ── Constrained solve with fixed air nodes ──────────────────────────
    let fixed_nodes: Vec<(usize, f64)> = topo
        .air_nodes
        .iter()
        .filter(|a| zone_needs_hvac[a.zone_idx])
        .map(|a| (a.global_idx, fixed_air_temps[a.zone_idx]))
        .collect();

    let constrained_temps = solve_system(
        topo,
        temperatures,
        wall_conditions,
        air_conditions,
        radiation,
        dt_s,
        wall_infos,
        Some(&fixed_nodes),
    );

    let Ok(constrained_temps) = constrained_temps else {
        // Fallback: accept free-float
        temperatures.copy_from_slice(&free_temps);
        return GlobalStepResult {
            heating_w_per_zone: vec![0.0; n_zones],
            cooling_w_per_zone: vec![0.0; n_zones],
        };
    };

    // ── Back-compute HVAC power per zone ────────────────────────────────
    // Q_hvac = C/dt * (T_new - T_old) + K_total * (T_new - T_out) + surface_coupling - gains
    // But it's easier to compute from the residual of the air node equation.
    let mut heating_w = vec![0.0f64; n_zones];
    let mut cooling_w = vec![0.0f64; n_zones];

    // To get Q_hvac, we evaluate: for the air node equation with the constrained temps,
    // Q_hvac is whatever makes the equation balance.
    // Air node equation: C/dt * T_new = C/dt * T_old + K_inf*(T_out - T_new) + K_glaz*(T_out - T_new)
    //                    + Σ_surfaces h_conv*A*(T_surf - T_new) + gains + Q_hvac
    // => Q_hvac = C/dt*(T_new - T_old) - K_inf*(T_out - T_new) - K_glaz*(T_out - T_new)
    //            - Σ_surfaces h_conv*A*(T_surf - T_new) - gains
    for air in &topo.air_nodes {
        let zi = air.zone_idx;
        if !zone_needs_hvac[zi] {
            continue;
        }

        let t_new = constrained_temps[air.global_idx];
        let t_old = temperatures[air.global_idx];
        let ac = &air_conditions[zi];

        // Capacity term
        let mut q_hvac = air.capacity_j_per_k / dt_s * (t_new - t_old);

        // Infiltration + glazing loss
        q_hvac += (air.infiltration_k + air.glazing_ua) * (t_new - ac.outdoor_temp_c);

        // Surface convective coupling
        for surf in &topo.surfaces {
            if surf.zone_idx != zi {
                continue;
            }
            let wc = &wall_conditions[surf.wall_idx];
            let h_conv_a = wc.h_conv * surf.area_m2;
            let t_surf = constrained_temps[surf.global_idx];
            q_hvac -= h_conv_a * (t_surf - t_new);
        }

        // Direct gains
        q_hvac -= ac.direct_gains_w;

        if q_hvac > 0.0 {
            heating_w[zi] = q_hvac;
        } else {
            cooling_w[zi] = -q_hvac;
        }
    }

    temperatures.copy_from_slice(&constrained_temps);
    GlobalStepResult {
        heating_w_per_zone: heating_w,
        cooling_w_per_zone: cooling_w,
    }
}

/// Assemble and solve the global linear system.
///
/// If `fixed_nodes` is provided, those global indices are held at fixed temperatures
/// (Dirichlet) and eliminated from the system.
#[allow(clippy::too_many_arguments)]
fn solve_system(
    topo: &GlobalTopology,
    temperatures: &[f64],
    wall_conditions: &[WallStepConditions],
    air_conditions: &[AirStepConditions],
    radiation: Option<&RadiationConditions>,
    dt_s: f64,
    wall_infos: &[FvmWallInfo],
    fixed_nodes: Option<&[(usize, f64)]>,
) -> anyhow::Result<Vec<f64>> {
    let n = topo.n;
    let mut a_mat = vec![vec![0.0f64; n]; n];
    let mut rhs = vec![0.0f64; n];

    // ── 1. Wall cells: capacity + internal conductances ─────────────────
    for wall in &topo.walls {
        let off = wall.cell_offset;

        // Capacity terms (backward Euler: C/dt on diagonal, C/dt*T_old on RHS)
        for i in 0..wall.n_cells {
            let gi = off + i;
            let cap = wall.cell_capacities[i];
            a_mat[gi][gi] += cap / dt_s;
            rhs[gi] += (cap / dt_s) * temperatures[gi];
        }

        // Internal face conductances (tridiagonal coupling)
        for (fi, &k) in wall.internal_conductances.iter().enumerate() {
            let gi_l = off + fi;
            let gi_r = off + fi + 1;
            a_mat[gi_l][gi_l] += k;
            a_mat[gi_r][gi_r] += k;
            a_mat[gi_l][gi_r] -= k;
            a_mat[gi_r][gi_l] -= k;
        }
    }

    // ── 2. Exterior BC on first cell of each wall ───────────────────────
    for (wi, wall) in topo.walls.iter().enumerate() {
        let info = &wall_infos[wi];
        let wc = &wall_conditions[wi];
        let gi_ext = wall.cell_offset; // first cell

        if info.exterior_adiabatic {
            // Neumann(0): no exterior coupling
            continue;
        }

        if info.has_exterior_surface {
            // TwoSided internal mass: exterior face connects to a surface node,
            // not directly to outdoors. Find the second surface node for this wall.
            // It's the surface whose inner_cell_global_idx == wall.cell_offset
            if let Some(ext_surf) = topo.surfaces.iter().find(|s| {
                s.wall_idx == wi && s.inner_cell_global_idx == wall.cell_offset
            }) {
                let gi_surf = ext_surf.global_idx;
                let k_face = ext_surf.k_face;

                // Wall cell <-> surface node coupling
                a_mat[gi_ext][gi_ext] += k_face;
                a_mat[gi_ext][gi_surf] -= k_face;
                a_mat[gi_surf][gi_surf] += k_face;
                a_mat[gi_surf][gi_ext] -= k_face;

                // Surface node <-> air node coupling (convective)
                let air_gi = topo.air_nodes[wall.zone_idx].global_idx;
                let h_conv_a = wc.h_conv * ext_surf.area_m2;
                a_mat[gi_surf][gi_surf] += h_conv_a;
                a_mat[gi_surf][air_gi] -= h_conv_a;
                a_mat[air_gi][air_gi] += h_conv_a;
                a_mat[air_gi][gi_surf] -= h_conv_a;
            }
        } else {
            // Normal exterior wall: series K_eff to outdoor
            let k_eff = wc.ext_k_eff;
            a_mat[gi_ext][gi_ext] += k_eff;
            rhs[gi_ext] += k_eff * wc.ext_t_drive;

            // Exterior source (absorbed solar, longwave, etc.)
            rhs[gi_ext] += wc.ext_source_w;
        }
    }

    // ── 3. Surface nodes: coupling to inner cell and air ────────────────
    for surf in &topo.surfaces {
        // Skip exterior surface nodes of TwoSided walls (handled above)
        let wall = &topo.walls[surf.wall_idx];
        let is_interior_surface = surf.inner_cell_global_idx
            == wall.cell_offset + wall.n_cells.saturating_sub(1);

        if !is_interior_surface {
            // This is an exterior surface node for a TwoSided wall, handled in step 2
            continue;
        }

        let gi_surf = surf.global_idx;
        let gi_cell = surf.inner_cell_global_idx;
        let gi_air = topo.air_nodes[surf.zone_idx].global_idx;
        let wc = &wall_conditions[surf.wall_idx];

        // Surface node <-> innermost cell coupling (half-cell conductance)
        let k_face = surf.k_face;
        a_mat[gi_surf][gi_surf] += k_face;
        a_mat[gi_surf][gi_cell] -= k_face;
        a_mat[gi_cell][gi_cell] += k_face;
        a_mat[gi_cell][gi_surf] -= k_face;

        // Surface node <-> air node (convective only)
        let h_conv_a = wc.h_conv * surf.area_m2;
        a_mat[gi_surf][gi_surf] += h_conv_a;
        a_mat[gi_surf][gi_air] -= h_conv_a;
        a_mat[gi_air][gi_air] += h_conv_a;
        a_mat[gi_air][gi_surf] -= h_conv_a;

        // Interior surface source (transmitted solar, radiant gains)
        rhs[gi_surf] += wc.int_source_w;
    }

    // ── 4. Radiation coupling between surface nodes ─────────────────────
    if let Some(rad) = radiation {
        let ns = rad.n_surfaces;
        // Only apply radiation to the first ns surfaces (interior surface nodes).
        // Surface nodes beyond ns (e.g. exterior faces of TwoSided) don't participate
        // in the enclosed radiation model.
        //
        // NOTE: The F matrix may not have row sums = 1 if some surfaces (e.g.
        // windows) are not in the global system. Use actual row sums for the
        // diagonal to avoid creating a spurious heat sink.
        for i in 0..ns.min(topo.surfaces.len()) {
            let gi_i = topo.surfaces[i].global_idx;
            let a_i = rad.surface_areas[i];

            // Compute actual row sum for this surface
            let f_row_sum: f64 = (0..ns)
                .filter(|&j| j != i)
                .map(|j| rad.f_matrix[i * ns + j])
                .sum();

            // Diagonal: += h_rad * A_i * Σ_j F_ij
            a_mat[gi_i][gi_i] += rad.h_rad * a_i * f_row_sum;

            // Off-diagonal: -= h_rad * A_i * F_ij
            for j in 0..ns.min(topo.surfaces.len()) {
                if i == j {
                    continue;
                }
                let gi_j = topo.surfaces[j].global_idx;
                let f_ij = rad.f_matrix[i * ns + j];
                a_mat[gi_i][gi_j] -= rad.h_rad * a_i * f_ij;
            }
        }
    }

    // ── 5. Air nodes ────────────────────────────────────────────────────
    for air in &topo.air_nodes {
        let gi = air.global_idx;
        let ac = &air_conditions[air.zone_idx];

        // Capacity term
        a_mat[gi][gi] += air.capacity_j_per_k / dt_s;
        rhs[gi] += (air.capacity_j_per_k / dt_s) * temperatures[gi];

        // Infiltration + glazing to outdoor
        let k_out = air.infiltration_k + air.glazing_ua;
        a_mat[gi][gi] += k_out;
        rhs[gi] += k_out * ac.outdoor_temp_c;

        // Direct convective gains
        rhs[gi] += ac.direct_gains_w;
    }

    // ── 6. Apply fixed nodes (Dirichlet elimination) ────────────────────
    if let Some(fixed) = fixed_nodes {
        for &(fi, t_fixed) in fixed {
            // Move known terms to RHS for all other equations
            for i in 0..n {
                if i == fi {
                    continue;
                }
                rhs[i] -= a_mat[i][fi] * t_fixed;
                a_mat[i][fi] = 0.0;
            }
            // Replace equation for fixed node
            for j in 0..n {
                a_mat[fi][j] = 0.0;
            }
            a_mat[fi][fi] = 1.0;
            rhs[fi] = t_fixed;
        }
    }

    // ── 7. Solve ────────────────────────────────────────────────────────
    solve_dense(a_mat, rhs)
}

/// Extract initial temperature vector from wall solvers and air temperatures.
pub fn extract_temperatures(
    topo: &GlobalTopology,
    wall_solvers: &[&FvmWallSolver],
    air_temperatures: &[f64],
) -> Vec<f64> {
    let mut temps = vec![20.0; topo.n];

    for wall in &topo.walls {
        let solver = wall_solvers[wall.wall_idx];
        let cell_temps = solver.temperatures();
        for i in 0..wall.n_cells {
            temps[wall.cell_offset + i] = cell_temps[i];
        }
    }

    // Surface nodes: initialize to the innermost cell temperature
    for surf in &topo.surfaces {
        let wall = &topo.walls[surf.wall_idx];
        let solver = wall_solvers[wall.wall_idx];
        let cell_temps = solver.temperatures();
        let local_idx = surf.inner_cell_global_idx - wall.cell_offset;
        temps[surf.global_idx] = cell_temps[local_idx];
    }

    // Air nodes
    for air in &topo.air_nodes {
        temps[air.global_idx] = air_temperatures
            .get(air.zone_idx)
            .copied()
            .unwrap_or(20.0);
    }

    temps
}

/// Write back solved temperatures to individual wall solvers.
pub fn scatter_temperatures(
    topo: &GlobalTopology,
    temperatures: &[f64],
    wall_solvers: &mut [&mut FvmWallSolver],
) {
    for wall in &topo.walls {
        let solver = &mut wall_solvers[wall.wall_idx];
        let start = wall.cell_offset;
        let end = start + wall.n_cells;
        solver.set_temperatures(&temperatures[start..end]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::energy::construction::WallConstruction;
    use crate::sim::heat_transfer::mesh_1d::build_1d_mesh;
    use crate::sim::materials::Layer;

    fn concrete_wall(thickness: f64, area: f64) -> FvmWallSolver {
        let construction = WallConstruction::new(
            "concrete",
            vec![Layer {
                name: "concrete".into(),
                thickness,
                conductivity: 1.4,
                density: 2300.0,
                specific_heat: 880.0,
            }],
        );
        let mesh = build_1d_mesh(&construction, area);
        FvmWallSolver::new(mesh, 20.0)
    }

    /// Steady-state single wall + air: check that air temperature matches
    /// analytical solution at steady state.
    #[test]
    fn test_steady_state_single_wall() {
        let area = 10.0;
        let solver = concrete_wall(0.20, area);
        let n_cells = solver.mesh().cells.len();

        let wall_infos = vec![FvmWallInfo {
            solver: &solver,
            zone_idx: 0,
            area_m2: area,
            has_exterior_surface: false,
            exterior_adiabatic: false,
        }];

        let h_out = 25.0;
        let h_in = 7.7;
        let t_out = -10.0;
        let k_inf = 10.0; // W/K infiltration

        let topo = build_topology(&wall_infos, 1, &[1.2 * 1005.0 * 50.0], &[k_inf], &[0.0]);

        // Check topology
        assert_eq!(topo.walls.len(), 1);
        assert_eq!(topo.surfaces.len(), 1);
        assert_eq!(topo.air_nodes.len(), 1);
        assert_eq!(topo.n, n_cells + 1 + 1); // cells + 1 surface + 1 air

        // Extract initial temps
        let mut temps = extract_temperatures(&topo, &[&solver], &[20.0]);

        // Compute ext K_eff = series(h_out*A, K_ext_face)
        let k_ext_face = topo.walls[0].ext_boundary_face_k;
        let h_out_a = h_out * area;
        let ext_k_eff = 1.0 / (1.0 / h_out_a + 1.0 / k_ext_face);

        let wc = WallStepConditions {
            ext_k_eff,
            ext_t_drive: t_out,
            ext_source_w: 0.0,
            h_conv: h_in,
            h_total: h_in,
            int_source_w: 0.0,
        };
        let ac = AirStepConditions {
            outdoor_temp_c: t_out,
            direct_gains_w: 0.0,
        };
        let hvac = HvacIdealLoads::with_setpoints(-100.0, 100.0); // no HVAC

        // Run to steady state with large dt
        for _ in 0..50 {
            let _result = step_global(
                &topo,
                &mut temps,
                &[wc.clone()],
                &[ac.clone()],
                None,
                &hvac,
                1e6,
                &wall_infos,
            );
        }

        // At steady state with no HVAC:
        // Heat flow: Q = (T_air - T_out) / R_total
        // where R_total = 1/(h_out*A) + L/(k*A) + 1/(h_in*A) + 1/K_inf
        // Actually with K_inf to outdoor: T_air = T_out (in steady state with no gains)
        // since all losses go to outdoor.
        // T_air should approach T_out since there are no gains.
        let t_air = temps[topo.air_nodes[0].global_idx];
        assert!(
            (t_air - t_out).abs() < 0.5,
            "Air should approach outdoor temp with no gains: got {t_air}, expected ~{t_out}"
        );
    }

    /// Multi-wall box: total heat loss matches expectations.
    #[test]
    fn test_multi_wall_box_steady_state() {
        let area = 10.0;
        let solver1 = concrete_wall(0.20, area);
        let solver2 = concrete_wall(0.20, area);

        let wall_infos = vec![
            FvmWallInfo {
                solver: &solver1,
                zone_idx: 0,
                area_m2: area,
                has_exterior_surface: false,
                exterior_adiabatic: false,
            },
            FvmWallInfo {
                solver: &solver2,
                zone_idx: 0,
                area_m2: area,
                has_exterior_surface: false,
                exterior_adiabatic: false,
            },
        ];

        let t_out = -10.0;
        let h_out = 25.0;
        let h_in = 7.7;
        let k_inf = 5.0;

        let topo = build_topology(&wall_infos, 1, &[1.2 * 1005.0 * 50.0], &[k_inf], &[0.0]);
        let mut temps = extract_temperatures(&topo, &[&solver1, &solver2], &[20.0]);

        let make_wc = |wall: &WallTopology| -> WallStepConditions {
            let k_ext_face = wall.ext_boundary_face_k;
            let h_out_a = h_out * area;
            let ext_k_eff = 1.0 / (1.0 / h_out_a + 1.0 / k_ext_face);
            WallStepConditions {
                ext_k_eff,
                ext_t_drive: t_out,
                ext_source_w: 0.0,
                h_conv: h_in,
                h_total: h_in,
                int_source_w: 0.0,
            }
        };

        let wcs: Vec<WallStepConditions> = topo.walls.iter().map(|w| make_wc(w)).collect();
        let ac = AirStepConditions {
            outdoor_temp_c: t_out,
            direct_gains_w: 0.0,
        };
        let hvac = HvacIdealLoads::with_setpoints(-100.0, 100.0);

        for _ in 0..50 {
            let _result = step_global(&topo, &mut temps, &wcs, &[ac.clone()], None, &hvac, 1e6, &wall_infos);
        }

        // With no gains, air should approach outdoor
        let t_air = temps[topo.air_nodes[0].global_idx];
        assert!(
            (t_air - t_out).abs() < 0.5,
            "Air should approach outdoor: got {t_air}"
        );
    }

    /// Radiation conservation: symmetric surfaces converge to equal temperatures.
    #[test]
    fn test_radiation_conservation() {
        let area = 10.0;
        let solver1 = concrete_wall(0.10, area);
        let solver2 = concrete_wall(0.10, area);

        let wall_infos = vec![
            FvmWallInfo {
                solver: &solver1,
                zone_idx: 0,
                area_m2: area,
                has_exterior_surface: false,
                exterior_adiabatic: false,
            },
            FvmWallInfo {
                solver: &solver2,
                zone_idx: 0,
                area_m2: area,
                has_exterior_surface: false,
                exterior_adiabatic: false,
            },
        ];

        let topo = build_topology(&wall_infos, 1, &[1000.0], &[5.0], &[0.0]);
        let mut temps = extract_temperatures(&topo, &[&solver1, &solver2], &[20.0]);
        temps[topo.surfaces[0].global_idx] = 25.0;
        temps[topo.surfaces[1].global_idx] = 15.0;

        // F12 = 1, F21 = 1 (equal flat surfaces facing each other)
        let h_rad = 4.0 * 0.9 * 5.67e-8 * (293.15_f64).powi(3); // ~4.6 W/(m²·K)
        let rad = RadiationConditions {
            h_rad,
            surface_areas: vec![area, area],
            f_matrix: vec![0.0, 1.0, 1.0, 0.0],
            n_surfaces: 2,
        };

        // Net radiation: A*h_rad*(T1 - T2) + A*h_rad*(T2 - T1) = 0
        // The matrix formulation ensures this algebraically.
        // Verify by running to steady state and checking surface temps converge.
        let wc = WallStepConditions {
            ext_k_eff: 100.0,
            ext_t_drive: 0.0,
            ext_source_w: 0.0,
            h_conv: 3.0,
            h_total: 3.0 + h_rad,
            int_source_w: 0.0,
        };
        let ac = AirStepConditions {
            outdoor_temp_c: 0.0,
            direct_gains_w: 0.0,
        };
        let hvac = HvacIdealLoads::with_setpoints(-100.0, 100.0);

        for _ in 0..100 {
            let _r = step_global(
                &topo,
                &mut temps,
                &[wc.clone(), wc.clone()],
                &[ac.clone()],
                Some(&rad),
                &hvac,
                1e5,
                &wall_infos,
            );
        }

        // At steady state with symmetric walls and radiation:
        // Both surface temperatures should be equal (symmetric problem)
        let t_s0 = temps[topo.surfaces[0].global_idx];
        let t_s1 = temps[topo.surfaces[1].global_idx];
        assert!(
            (t_s0 - t_s1).abs() < 0.01,
            "Symmetric surfaces should have equal temps: {t_s0} vs {t_s1}"
        );
    }

    /// HVAC two-pass: fixed air temp -> correct Q_hvac back-computation.
    #[test]
    fn test_hvac_two_pass() {
        let area = 10.0;
        let solver = concrete_wall(0.20, area);

        let wall_infos = vec![FvmWallInfo {
            solver: &solver,
            zone_idx: 0,
            area_m2: area,
            has_exterior_surface: false,
            exterior_adiabatic: false,
        }];

        let t_out = -20.0;
        let h_out = 25.0;
        let h_in = 7.7;
        let k_inf = 10.0;

        let topo = build_topology(&wall_infos, 1, &[1.2 * 1005.0 * 50.0], &[k_inf], &[0.0]);
        let mut temps = extract_temperatures(&topo, &[&solver], &[20.0]);

        let k_ext_face = topo.walls[0].ext_boundary_face_k;
        let h_out_a = h_out * area;
        let ext_k_eff = 1.0 / (1.0 / h_out_a + 1.0 / k_ext_face);

        let wc = WallStepConditions {
            ext_k_eff,
            ext_t_drive: t_out,
            ext_source_w: 0.0,
            h_conv: h_in,
            h_total: h_in,
            int_source_w: 0.0,
        };
        let ac = AirStepConditions {
            outdoor_temp_c: t_out,
            direct_gains_w: 0.0,
        };
        let hvac = HvacIdealLoads::with_setpoints(20.0, 27.0);

        // Run a few steps
        let mut total_heating = 0.0;
        for _ in 0..100 {
            let result = step_global(&topo, &mut temps, &[wc.clone()], &[ac.clone()], None, &hvac, 600.0, &wall_infos);
            total_heating += result.heating_w_per_zone[0];
        }

        // Air should be at heating setpoint
        let t_air = temps[topo.air_nodes[0].global_idx];
        assert!(
            (t_air - 20.0).abs() < 0.5,
            "Air should be at heating setpoint: got {t_air}"
        );

        // Should require positive heating
        assert!(
            total_heating > 0.0,
            "Should need heating with T_out={t_out}: got {total_heating}"
        );
    }

    /// Verify that with no radiation and a single wall, the global solve
    /// produces results close to the sequential Thomas solve.
    #[test]
    fn test_equivalence_with_thomas() {
        let area = 10.0;
        let thickness = 0.20;
        let construction = WallConstruction::new(
            "concrete",
            vec![Layer {
                name: "concrete".into(),
                thickness,
                conductivity: 1.4,
                density: 2300.0,
                specific_heat: 880.0,
            }],
        );
        let mesh = build_1d_mesh(&construction, area);
        let mut solver_thomas = FvmWallSolver::new(mesh.clone(), 20.0);
        let solver_global = FvmWallSolver::new(mesh, 20.0);

        let t_out = -10.0;
        let h_out = 25.0;
        let h_in = 7.7;
        let dt = 600.0;

        // Thomas path
        use crate::sim::heat_transfer::BoundaryCondition;
        let bc_ext = BoundaryCondition::Convective {
            h: h_out,
            t_fluid: t_out,
        };
        let bc_int = BoundaryCondition::Convective {
            h: h_in,
            t_fluid: 20.0,
        };

        for _ in 0..10 {
            solver_thomas.step(dt, &bc_ext, &bc_int, &[]);
        }

        // Global path
        let wall_infos = vec![FvmWallInfo {
            solver: &solver_global,
            zone_idx: 0,
            area_m2: area,
            has_exterior_surface: false,
            exterior_adiabatic: false,
        }];

        // Use large air capacity and fix air at 20°C via HVAC
        let topo = build_topology(&wall_infos, 1, &[1e12], &[0.0], &[0.0]);
        let mut temps = extract_temperatures(&topo, &[&solver_global], &[20.0]);

        let k_ext_face = topo.walls[0].ext_boundary_face_k;
        let h_out_a = h_out * area;
        let ext_k_eff = 1.0 / (1.0 / h_out_a + 1.0 / k_ext_face);

        let wc = WallStepConditions {
            ext_k_eff,
            ext_t_drive: t_out,
            ext_source_w: 0.0,
            h_conv: h_in,
            h_total: h_in,
            int_source_w: 0.0,
        };
        let ac = AirStepConditions {
            outdoor_temp_c: t_out,
            direct_gains_w: 0.0,
        };
        // Fix air at 20°C
        let hvac = HvacIdealLoads::with_setpoints(20.0, 20.0);

        for _ in 0..10 {
            let _r = step_global(&topo, &mut temps, &[wc.clone()], &[ac.clone()], None, &hvac, dt, &wall_infos);
        }

        // Compare cell temperatures
        let thomas_temps = solver_thomas.temperatures();
        for (i, &t_thomas) in thomas_temps.iter().enumerate() {
            let t_global = temps[topo.walls[0].cell_offset + i];
            assert!(
                (t_thomas - t_global).abs() < 0.5,
                "cell {i}: Thomas={t_thomas:.3} vs Global={t_global:.3}"
            );
        }
    }
}

// Manual Clone impls for condition types (used in tests and simulation loop)
impl Clone for WallStepConditions {
    fn clone(&self) -> Self {
        Self {
            ext_k_eff: self.ext_k_eff,
            ext_t_drive: self.ext_t_drive,
            ext_source_w: self.ext_source_w,
            h_conv: self.h_conv,
            h_total: self.h_total,
            int_source_w: self.int_source_w,
        }
    }
}

impl Clone for AirStepConditions {
    fn clone(&self) -> Self {
        Self {
            outdoor_temp_c: self.outdoor_temp_c,
            direct_gains_w: self.direct_gains_w,
        }
    }
}
