use anyhow::{Result, anyhow};

pub mod point;
pub mod polygon;
pub mod vector;
pub mod bboxes;
pub mod triangles;
pub mod wall;

/// Path separator (building/zone/solid/wall/polygon)
const SEP: char = '/';

/// Validates object name w.r.t. to the presence of the path separator.
fn validate_name(name: &str) -> Result<String> {
    if name.contains(SEP) {
        Err(anyhow!("Name {} contains forbidden character '{}'", name, SEP))
    } else {
        Ok(name.to_string())
    }
}

/// Geometric precision
const EPS: f64 = 1e-13;

/// Trait enabling to check if two f64 floats are almost equal.
trait IsClose {
    /// Checks if this float is almost equal to the other.
    fn is_close(&self, other: Self) -> bool;
}

impl IsClose for f64 {
    fn is_close(&self, other: Self) -> bool {
        let remainder = (self - other).abs();
        remainder <= EPS
    }
}

