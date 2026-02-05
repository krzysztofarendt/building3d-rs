//! File I/O for building models.
//!
//! This module provides functions for reading and writing building models
//! in various formats:
//! - B3D: Native JSON format preserving full hierarchy
//! - STL: Triangulated mesh format (ASCII/Binary)
//! - BIM: dotbim format for BIM interoperability

pub mod b3d;
pub mod bim;
pub mod stl;

pub use b3d::{read_b3d, write_b3d};
pub use bim::{read_bim, write_bim};
pub use stl::{StlFormat, read_stl, write_stl};
