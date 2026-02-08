use anyhow::{Context, Result};

/// Solves a dense linear system `A * x = b` using Gaussian elimination with
/// partial pivoting.
///
/// Intended for small systems (e.g. multi-zone air node solves). This keeps the
/// dependency surface small and deterministic.
pub fn solve_dense(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Result<Vec<f64>> {
    let n = a.len();
    if n == 0 {
        return Ok(vec![]);
    }
    anyhow::ensure!(b.len() == n, "b length mismatch");
    for (i, row) in a.iter().enumerate() {
        anyhow::ensure!(row.len() == n, "A row {i} length mismatch");
    }

    // Forward elimination
    for col in 0..n {
        // Pivot selection
        let mut pivot_row = col;
        let mut pivot_val = a[col][col].abs();
        for r in (col + 1)..n {
            let v = a[r][col].abs();
            if v > pivot_val {
                pivot_val = v;
                pivot_row = r;
            }
        }

        anyhow::ensure!(
            pivot_val > 1e-14,
            "Singular matrix (pivot too small) at column {col}"
        );

        if pivot_row != col {
            a.swap(pivot_row, col);
            b.swap(pivot_row, col);
        }

        let pivot = a[col][col];
        for r in (col + 1)..n {
            let factor = a[r][col] / pivot;
            if factor == 0.0 {
                continue;
            }
            a[r][col] = 0.0;
            for c in (col + 1)..n {
                a[r][c] -= factor * a[col][c];
            }
            b[r] -= factor * b[col];
        }
    }

    // Back substitution
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut rhs = b[i];
        for j in (i + 1)..n {
            rhs -= a[i][j] * x[j];
        }
        x[i] = rhs / a[i][i];
    }

    // Basic sanity: reject NaNs/Infs early.
    for (i, xi) in x.iter().enumerate() {
        xi.is_finite()
            .then_some(())
            .with_context(|| format!("Non-finite solution at index {i}: {xi}"))?;
    }

    Ok(x)
}
