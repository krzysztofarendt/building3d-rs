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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_max() {
        assert_eq!(max(&[1.0, 3.0, 2.0]), 3.0);
        assert_eq!(max(&[-5.0, -1.0, -3.0]), -1.0);
        assert_eq!(max(&[42.0]), 42.0);
    }

    #[test]
    fn test_min() {
        assert_eq!(min(&[1.0, 3.0, 2.0]), 1.0);
        assert_eq!(min(&[-5.0, -1.0, -3.0]), -5.0);
        assert_eq!(min(&[42.0]), 42.0);
    }

    #[test]
    fn test_flip() {
        let v = vec![1, 2, 3, 4];
        assert_eq!(flip(&v), vec![4, 3, 2, 1]);
        let empty: Vec<i32> = vec![];
        assert_eq!(flip(&empty), Vec::<i32>::new());
    }

    #[test]
    fn test_roll() {
        let mut v = vec![1, 2, 3, 4];
        roll(&mut v, 1);
        assert_eq!(v, vec![4, 1, 2, 3]);
    }

    #[test]
    fn test_almost_equal_true() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![1.0, 2.0, 3.0];
        assert!(almost_equal(&a, &b, 1e-10));
    }

    #[test]
    fn test_almost_equal_false() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![1.0, 2.0, 4.0];
        assert!(!almost_equal(&a, &b, 1e-10));
    }

    #[test]
    fn test_almost_equal_different_lengths() {
        let a = vec![1.0, 2.0];
        let b = vec![1.0, 2.0, 3.0];
        assert!(!almost_equal(&a, &b, 1e-10));
    }
}
