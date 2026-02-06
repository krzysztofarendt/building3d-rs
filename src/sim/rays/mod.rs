mod config;
mod simulation;

pub use config::SimulationConfig;
pub use simulation::{Simulation, SimulationResult};

// Re-export VoxelGrid from engine for backward compatibility
pub use crate::sim::engine::voxel_grid::VoxelGrid;
