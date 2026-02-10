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
    /// Convective with an additional imposed heat flux into the domain.
    ///
    /// This represents an exterior surface balance term like absorbed shortwave:
    /// `q_total = h*(T_fluid - T_surface) + heat_flux`.
    ///
    /// Important: `heat_flux` is a *surface source* (e.g. absorbed shortwave) in W/m²,
    /// positive into the wall domain. When `h > 0`, only a fraction of this source
    /// enters the wall; the remainder is lost immediately via convection to `t_fluid`
    /// (the split depends on the half-cell conductance at the boundary).
    ConvectiveWithFlux {
        h: f64,
        t_fluid: f64,
        heat_flux: f64,
    },
}
