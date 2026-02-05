use super::*;

/// Checks if all vectors are (almost) equal
pub fn are_vectors_close(vectors: &[Vector]) -> bool {
    if vectors.is_empty() {
        return true;
    }
    let mut all_close = true;
    let v0 = vectors[0];
    for v in vectors.iter().skip(1) {
        if !v.is_close(&v0) {
            all_close = false;
            break;
        }
    }
    all_close
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_are_vectors_close() {
        let v0 = Vector::new(1.1, 2.2, 3.3);
        let v1 = Vector::new(1.1, 2.2, 3.3);
        let v2 = Vector::new(1.1, 2.2, 3.3);
        let v3 = Vector::new(1.1, 2.2, 3.3);
        assert!(are_vectors_close(&[v0, v1, v2, v3]));
    }
}
