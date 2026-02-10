/// Boundary condition applied at a mesh boundary face.
#[derive(Debug, Clone, Copy)]
pub enum BoundaryCondition {
    /// Fixed temperature: T_surface = temperature.
    Dirichlet { temperature: f64 },
    /// Fixed heat flux: q = heat_flux  [W/m^2], positive into the domain.
    Neumann { heat_flux: f64 },
    /// Convective: q = h * (T_fluid - T_surface).
    ///
    /// This is the primary mode in building simulation.
    /// Exterior: h ~ 10-25 W/(m^2*K), interior: h ~ 3-8 W/(m^2*K).
    Convective { h: f64, t_fluid: f64 },
}
