mod config;
mod simulation;

pub use config::{AcousticMode, SimulationConfig};
pub use simulation::{ENERGY_EPS, Simulation, SimulationProgress, SimulationResult};

// Re-export VoxelGrid from engine for backward compatibility
pub use crate::sim::engine::voxel_grid::VoxelGrid;
