//! File I/O for building models.
//!
//! This module provides functions for reading and writing building models
//! in various formats.

pub mod b3d;
pub mod stl;

pub use b3d::{read_b3d, write_b3d};
pub use stl::{read_stl, write_stl, StlFormat};
