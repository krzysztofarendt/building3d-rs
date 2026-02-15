//! Finite Volume Method (FVM) heat transfer solver.
//!
//! Provides 1D conduction through wall layers using a dimension-agnostic FVM
//! mesh.  The same mesh/solver design extends to 3D volumetric conduction on
//! tetrahedral meshes.
//!
//! # Architecture
//!
//! ```text
//! WallConstruction ──► build_1d_mesh() ──► FvmMesh ──► FvmWallSolver
//!                                                        │
//!                                               step() / interior_heat_flux()
//! ```
//!
//! The solver sees only cells and faces (never coordinates), so the same
//! assembly and solve code works for any dimensionality.

pub mod boundary;
pub mod mesh;
pub mod mesh_1d;
pub mod mesh_3d;
pub mod solver;
pub mod solver_sparse;

pub use boundary::BoundaryCondition;
pub use mesh::{FvmCell, FvmFace, FvmMesh};
pub use mesh_1d::build_1d_mesh;
pub use mesh_3d::{CellThermalProperties, build_3d_mesh, build_3d_mesh_uniform};
pub use solver::{FvmSolverSnapshot, FvmWallSolver};
pub use solver_sparse::{FvmSparseSolver, SparseSolverConfig};
