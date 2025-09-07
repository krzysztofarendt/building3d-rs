pub mod point;
pub mod polygon;
pub mod vector;

/// Geometric precision
const EPS: f64 = 1e-13;

/// Trait enabling to check if two f64 floats are almost equal.
trait IsClose {
    /// Checks if this float is almost equal to the other.
    fn is_close(&self, other: &Self) -> bool;
}

impl IsClose for f64 {
    fn is_close(&self, other: &Self) -> bool {
        let remainder = (*self - *other).abs();
        remainder <= EPS
    }
}
