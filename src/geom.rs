use anyhow::{Result, anyhow};

pub mod bboxes;
pub mod building;
pub mod mesh;
pub mod point;
pub mod polygon;
pub mod ray;
pub mod rotation;
pub mod segment;
pub mod solid;
pub mod tetrahedron;
pub mod triangles;
pub mod vector;
pub mod visibility;
pub mod wall;
pub mod zone;

/// Path separator (building/zone/solid/wall/polygon)
const SEP: char = '/';

/// Validates object name w.r.t. to the presence of the path separator.
fn validate_name(name: &str) -> Result<String> {
    if name.contains(SEP) {
        Err(anyhow!(
            "Name {} contains forbidden character '{}'",
            name,
            SEP
        ))
    } else {
        Ok(name.to_string())
    }
}

/// Geometric precision
const EPS: f64 = 1e-13;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_validate_name_with_separator() {
        let result = validate_name("zone/solid");
        assert!(result.is_err());
    }

    #[test]
    fn test_validate_name_ok() {
        let result = validate_name("valid_name");
        assert!(result.is_ok());
        assert_eq!(result.unwrap(), "valid_name");
    }
}

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
