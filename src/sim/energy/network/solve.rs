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

/// Solves an implicit-conduction-style linear system with optional fixed nodes.
///
/// This helper is intended for small thermal networks where the matrix is naturally expressed as:
/// - `diag[i]`: diagonal term for node i (e.g., `C/dt + ΣK`)
/// - `neighbors[i]`: list of `(j, k_ij)` conductances to other unknown nodes
/// - `rhs_base[i]`: RHS term excluding neighbor unknowns (e.g., `C/dt*T_prev + K_boundary*T_boundary + Q_sources`)
///
/// Fixed nodes (Dirichlet temperatures) are eliminated into the RHS.
pub fn solve_with_fixed_nodes(
    diag: &[f64],
    rhs_base: &[f64],
    neighbors: &[Vec<(usize, f64)>],
    is_fixed: &[bool],
    fixed_value: &[f64],
) -> Result<Vec<f64>> {
    let n = diag.len();
    anyhow::ensure!(rhs_base.len() == n, "rhs length mismatch");
    anyhow::ensure!(neighbors.len() == n, "neighbors length mismatch");
    anyhow::ensure!(is_fixed.len() == n, "is_fixed length mismatch");
    anyhow::ensure!(fixed_value.len() == n, "fixed_value length mismatch");

    let mut unknown = Vec::new();
    let mut pos = vec![usize::MAX; n];
    for i in 0..n {
        if !is_fixed[i] {
            pos[i] = unknown.len();
            unknown.push(i);
        }
    }

    let mut temps = vec![0.0; n];
    for i in 0..n {
        if is_fixed[i] {
            temps[i] = fixed_value[i];
        }
    }

    if unknown.is_empty() {
        return Ok(temps);
    }

    let m = unknown.len();
    let mut a = vec![vec![0.0; m]; m];
    let mut b = vec![0.0; m];

    for (row_idx, &i) in unknown.iter().enumerate() {
        a[row_idx][row_idx] = diag[i];
        let mut rhs = rhs_base[i];
        for &(j, k) in &neighbors[i] {
            anyhow::ensure!(j < n, "neighbor index out of bounds");
            if is_fixed[j] {
                rhs += k * fixed_value[j];
            } else {
                let col_idx = pos[j];
                a[row_idx][col_idx] -= k;
            }
        }
        b[row_idx] = rhs;
    }

    let x = solve_dense(a, b)?;
    for (row_idx, &i) in unknown.iter().enumerate() {
        temps[i] = x[row_idx];
    }

    Ok(temps)
}
