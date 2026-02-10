/// Sentinel index indicating a boundary (no neighbor cell).
pub const BOUNDARY: usize = usize::MAX;

/// A single finite-volume cell.
///
/// In 1D this represents a slab of wall material; in 3D it would be a
/// tetrahedron.  The solver sees only volume and material properties.
#[derive(Debug, Clone)]
pub struct FvmCell {
    /// Cell volume in m^3 (1D: wall_area * dx).
    pub volume: f64,
    /// Thermal conductivity in W/(m*K).
    pub conductivity: f64,
    /// Density in kg/m^3.
    pub density: f64,
    /// Specific heat capacity in J/(kg*K).
    pub specific_heat: f64,
}

impl FvmCell {
    /// Thermal capacity of this cell: rho * c_p * V  [J/K].
    pub fn capacity(&self) -> f64 {
        self.density * self.specific_heat * self.volume
    }
}

/// A face shared between two cells (or between a cell and a boundary).
///
/// `cell_left` and `cell_right` are indices into `FvmMesh::cells`.
/// If one side is a boundary, that index is [`BOUNDARY`].
#[derive(Debug, Clone)]
pub struct FvmFace {
    /// Left cell index (exterior side in 1D), or [`BOUNDARY`].
    pub cell_left: usize,
    /// Right cell index (interior side in 1D), or [`BOUNDARY`].
    pub cell_right: usize,
    /// Face area in m^2 (1D: wall polygon area).
    pub area: f64,
    /// Centroid-to-centroid distance across this face in m.
    pub distance: f64,
    /// Precomputed conductance: k_face * area / distance  [W/K].
    pub conductance: f64,
}

/// Dimension-agnostic finite-volume mesh.
///
/// The solver operates on this mesh without knowing whether it came from a 1D
/// wall layering or a 3D tetrahedral decomposition.
#[derive(Debug, Clone)]
pub struct FvmMesh {
    pub cells: Vec<FvmCell>,
    pub faces: Vec<FvmFace>,
}

impl FvmMesh {
    /// Returns the index of the exterior boundary face (first face whose
    /// `cell_left` is [`BOUNDARY`]).
    pub fn exterior_boundary_face(&self) -> Option<usize> {
        self.faces.iter().position(|f| f.cell_left == BOUNDARY)
    }

    /// Returns the index of the interior boundary face (last face whose
    /// `cell_right` is [`BOUNDARY`]).
    pub fn interior_boundary_face(&self) -> Option<usize> {
        self.faces.iter().rposition(|f| f.cell_right == BOUNDARY)
    }
}
