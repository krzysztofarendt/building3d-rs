//! Vector utility functions like min(), max(), roll()
pub fn max(vec: &[f64]) -> f64 {
    vec.iter().cloned().max_by(f64::total_cmp).unwrap()
}

pub fn min(vec: &[f64]) -> f64 {
    vec.iter().cloned().min_by(f64::total_cmp).unwrap()
}

/// Roll vector or array forward by k. It modifies the collection in-place.
pub fn roll<T>(v: &mut [T], k: usize) {
    v.rotate_right(k);
}

/// Reverses `v` and returns it as a new vector.
pub fn flip<T: Clone>(v: &[T]) -> Vec<T> {
    let mut reversed: Vec<T> = Vec::new();
    for i in (0..v.len()).rev() {
        reversed.push(v[i].clone());
    }
    reversed
}

/// Checks if two arrays or vectors are almost equal.
///
/// Elements in both containers must be in the same order.
pub fn almost_equal(a: &[f64], b: &[f64], eps: f64) -> bool {
    if a.len() != b.len() {
        return false;
    }
    a.iter().zip(b.iter()).all(|(&x, &y)| (x - y).abs() <= eps)
}
