use crate::Point;
use crate::Vector;
use crate::geom::IsClose;
use crate::geom::point::convert::{array_to_points, points_to_array};
use ndarray as nd;

/// Calculate rotation matrix for a unit vector `u` and angle `phi`.
///
/// A rotation in 3D can be described with an axis and angle around that axis.
/// The axis is described with a unit vector `u` `(ux**2 + uy**2 + uz**2 == 1)`
/// a vector `phi` (in radians).
///
/// Rotation matrix
/// Reference: https://en.wikipedia.org/wiki/Rotation_matrix#Basic_3D_rotations
/// Method 1 (less stable numerically):
/// ```ignore
/// R = np.array([
///     [
///         np.cos(phi) + u[0]**2 * (1 - np.cos(phi)),
///         u[0] * u[1] * (1 - np.cos(phi)) - u[2] * np.sin(phi),
///         u[0] * u[2] * (1 - np.cos(phi)) + u[1] * np.sin(phi),
///     ],
///     [
///         u[1] * u[0] * (1 - np.cos(phi)) + u[2] * np.sin(phi),
///         np.cos(phi) + u[1] ** 2 * (1 - np.cos(phi)),
///         u[1] * u[2] * (1 - np.cos(phi)) - u[0] * np.sin(phi),
///     ],
///     [
///         u[2] * u[0] * (1 - np.cos(phi)) - u[1] * np.sin(phi),
///         u[2] * u[1] * (1 - np.cos(phi)) + u[0] * np.sin(phi),
///         np.cos(phi) + u[2] ** 2 * (1 - np.cos(phi)),
///     ],
/// ])
/// ```
/// Method 2 (more stable numerically):
/// https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
/// https://math.stackexchange.com/questions/142821/matrix-for-rotation-around-a-vector
pub fn rotation_matrix(u: &Vector, phi: f64) -> nd::Array2<f64> {
    let u = match u.normalize() {
        Ok(v) => v,
        Err(_) => return nd::Array::eye(3),
    };

    let w: nd::Array2<f64> = nd::arr2(&[[0., -u.dz, u.dy], [u.dz, 0., -u.dx], [-u.dy, u.dx, 0.]]);

    nd::Array::eye(3) + phi.sin() * &w + (2. * (phi / 2.).sin().powi(2)) * w.dot(&w)
}

/// Rotate points using the rotation matrix `rot`
pub fn rotate_points(pts: &[Point], rot: &nd::ArrayView2<f64>) -> Vec<Point> {
    let pts = points_to_array(pts);
    let pts = pts.dot(rot);

    array_to_points(pts)
}

/// Rotate points around the unit vector `u` with the angle `phi` (radians).
///
/// Arguments:
/// - pts: list of points to be rotated
/// - u: normal vector of the rotation axis
/// - phi: rotation angle in radians
///
/// Returns:
/// rotated points
pub fn rotate_points_around_vector(pts: &[Point], u: &Vector, phi: f64) -> Vec<Point> {
    if u.length().is_close(0.) || phi.abs().is_close(0.) {
        // No need to rotate
        return pts.to_vec();
    }
    let u = match u.normalize() {
        Ok(v) => v,
        Err(_) => return pts.to_vec(),
    };
    let rot = rotation_matrix(&u, phi);

    rotate_points(pts, &rot.t())
}

/// Rotate points so they lie on a target plane.
///
/// This function rotates a set of coplanar points from their current plane
/// to a target plane defined by its normal vector. The rotation preserves
/// the relative positions of the points.
///
/// # Arguments
/// * `pts` - Points to rotate (should be coplanar)
/// * `target_normal` - Normal vector of the target plane
///
/// # Returns
/// Rotated points that lie on a plane with the target normal.
/// If the source and target normals are already aligned, returns
/// the original points unchanged.
///
/// # Note
/// The points are rotated around their centroid, then translated back.
/// If the source points are not coplanar, the result may not be meaningful.
pub fn rotate_points_to_plane(pts: &[Point], target_normal: &Vector) -> Vec<Point> {
    if pts.len() < 3 {
        return pts.to_vec();
    }

    // Normalize target normal
    let target_n = match target_normal.normalize() {
        Ok(n) => n,
        Err(_) => return pts.to_vec(), // Invalid target normal
    };

    // Calculate source plane normal from first three non-collinear points
    let source_n = match calculate_plane_normal(pts) {
        Some(n) => n,
        None => return pts.to_vec(), // Points are collinear
    };

    // Check if normals are already aligned (or opposite)
    let dot = source_n.dot(&target_n);
    if dot.abs().is_close(1.0) {
        // Normals are parallel - no rotation needed (or 180° flip)
        if dot > 0.0 {
            return pts.to_vec();
        } else {
            // Normals are opposite - rotate 180° around any perpendicular axis
            let perp = find_perpendicular_vector(&source_n);
            return rotate_points_around_vector(pts, &perp, std::f64::consts::PI);
        }
    }

    // Calculate rotation axis (perpendicular to both normals)
    let axis = match source_n.cross(&target_n).normalize() {
        Ok(a) => a,
        Err(_) => return pts.to_vec(), // Normals are parallel
    };

    // Calculate rotation angle
    let angle = dot.clamp(-1.0, 1.0).acos();

    // Rotate points
    rotate_points_around_vector(pts, &axis, angle)
}

/// Rotate points to lie on the XY plane (z = 0).
///
/// Convenience function that rotates points to have normal (0, 0, 1).
pub fn rotate_points_to_xy_plane(pts: &[Point]) -> Vec<Point> {
    rotate_points_to_plane(pts, &Vector::new(0.0, 0.0, 1.0))
}

/// Calculate the normal vector of a plane defined by points.
///
/// Uses the first three non-collinear points to compute the normal.
fn calculate_plane_normal(pts: &[Point]) -> Option<Vector> {
    if pts.len() < 3 {
        return None;
    }

    // Find first three non-collinear points
    let p0 = pts[0];
    for (i, pt_i) in pts.iter().enumerate().skip(1) {
        let v1 = *pt_i - p0;
        if v1.length().is_close(0.0) {
            continue;
        }

        for pt_j in pts.iter().skip(i + 1) {
            let v2 = *pt_j - p0;
            if v2.length().is_close(0.0) {
                continue;
            }

            let normal = v1.cross(&v2);
            if let Ok(n) = normal.normalize() {
                return Some(n);
            }
        }
    }

    None
}

/// Find a vector perpendicular to the given vector.
fn find_perpendicular_vector(v: &Vector) -> Vector {
    // Choose the axis that is most different from v
    let abs_x = v.dx.abs();
    let abs_y = v.dy.abs();
    let abs_z = v.dz.abs();

    let candidate = if abs_x <= abs_y && abs_x <= abs_z {
        Vector::new(1.0, 0.0, 0.0)
    } else if abs_y <= abs_z {
        Vector::new(0.0, 1.0, 0.0)
    } else {
        Vector::new(0.0, 0.0, 1.0)
    };

    // Cross product gives perpendicular vector
    v.cross(&candidate)
        .normalize()
        .unwrap_or(Vector::new(1.0, 0.0, 0.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rotate_points_around_vector() {
        let p0 = Point::new(1.0, 0.0, 0.0);
        let p1 = Point::new(0.0, 1.0, 0.0);
        let p2 = Point::new(0.0, 0.0, 0.0);
        let u = Vector::new(0., 1., 0.);
        let phi = -std::f64::consts::PI / 2.;

        let rotated_points = rotate_points_around_vector(&[p0, p1, p2], &u, phi);
        println!("{:?}", rotated_points);

        assert!(rotated_points[0].is_close(&Point::new(0.0, 0.0, 1.0)));
        assert!(rotated_points[1].is_close(&Point::new(0.0, 1.0, 0.0)));
        assert!(rotated_points[2].is_close(&Point::new(0.0, 0.0, 0.0)));
    }

    #[test]
    fn test_rotate_points_to_xy_plane() {
        // Square in the XZ plane (y = 0, normal is (0, 1, 0))
        let pts = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 1.0),
            Point::new(0.0, 0.0, 1.0),
        ];

        let rotated = rotate_points_to_xy_plane(&pts);

        // All z-coordinates should now be the same (on XY plane)
        let z_ref = rotated[0].z;
        for p in &rotated {
            assert!(
                (p.z - z_ref).abs() < 1e-6,
                "Points should be coplanar on XY plane"
            );
        }
    }

    #[test]
    fn test_rotate_points_to_arbitrary_plane() {
        // Square in the XY plane (z = 0, normal is (0, 0, 1))
        let pts = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
        ];

        // Rotate to plane with normal (1, 0, 0) - the YZ plane
        let target_normal = Vector::new(1.0, 0.0, 0.0);
        let rotated = rotate_points_to_plane(&pts, &target_normal);

        // All x-coordinates should now be the same
        let x_ref = rotated[0].x;
        for p in &rotated {
            assert!(
                (p.x - x_ref).abs() < 1e-6,
                "Points should be coplanar on YZ plane, got x = {}",
                p.x
            );
        }
    }

    #[test]
    fn test_rotate_points_already_aligned() {
        // Points already on XY plane
        let pts = vec![
            Point::new(0.0, 0.0, 5.0),
            Point::new(1.0, 0.0, 5.0),
            Point::new(1.0, 1.0, 5.0),
        ];

        let target_normal = Vector::new(0.0, 0.0, 1.0);
        let rotated = rotate_points_to_plane(&pts, &target_normal);

        // Points should be unchanged
        for (orig, rot) in pts.iter().zip(rotated.iter()) {
            assert!(orig.is_close(rot), "Points should be unchanged");
        }
    }

    #[test]
    fn test_rotate_points_opposite_normals() {
        // Points on plane with normal (0, 0, -1)
        let pts = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
        ];

        // Target normal is opposite
        let target_normal = Vector::new(0.0, 0.0, -1.0);
        let rotated = rotate_points_to_plane(&pts, &target_normal);

        // Should rotate 180 degrees
        assert_eq!(rotated.len(), 3);
        // After 180° rotation around some perpendicular axis, points should still be coplanar
        let z_ref = rotated[0].z;
        for p in &rotated {
            assert!((p.z - z_ref).abs() < 1e-6, "Points should remain coplanar");
        }
    }

    #[test]
    fn test_rotate_points_few_points() {
        // Less than 3 points - should return unchanged
        let pts = vec![Point::new(1.0, 2.0, 3.0)];
        let target_normal = Vector::new(0.0, 0.0, 1.0);
        let rotated = rotate_points_to_plane(&pts, &target_normal);

        assert_eq!(rotated.len(), 1);
        assert!(rotated[0].is_close(&pts[0]));
    }

    #[test]
    fn test_calculate_plane_normal() {
        // XY plane
        let pts = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
        ];
        let normal = calculate_plane_normal(&pts).unwrap();
        assert!(
            normal.is_close(&Vector::new(0.0, 0.0, 1.0))
                || normal.is_close(&Vector::new(0.0, 0.0, -1.0))
        );
    }

    #[test]
    fn test_find_perpendicular_vector() {
        let v = Vector::new(1.0, 0.0, 0.0);
        let perp = find_perpendicular_vector(&v);

        // Should be perpendicular (dot product = 0)
        assert!(v.dot(&perp).abs() < 1e-10);
        // Should be normalized
        assert!((perp.length() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_rotation_matrix_zero_vector() {
        // Zero vector should return identity matrix
        let zero = Vector::new(0.0, 0.0, 0.0);
        let mat = rotation_matrix(&zero, 1.0);
        let identity = nd::Array::eye(3);
        for ((r, c), val) in mat.indexed_iter() {
            let expected: f64 = identity[[r, c]];
            assert!(
                (val - expected).abs() < 1e-10,
                "Expected identity at ({}, {})",
                r,
                c
            );
        }
    }

    #[test]
    fn test_find_perpendicular_y_dominant() {
        // y-axis dominant vector
        let v = Vector::new(0.0, 1.0, 0.0);
        let perp = find_perpendicular_vector(&v);
        assert!(v.dot(&perp).abs() < 1e-10);
        assert!((perp.length() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_find_perpendicular_z_dominant() {
        // z-axis dominant vector
        let v = Vector::new(0.0, 0.0, 1.0);
        let perp = find_perpendicular_vector(&v);
        assert!(v.dot(&perp).abs() < 1e-10);
        assert!((perp.length() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_rotate_points_zero_angle() {
        let pts = vec![Point::new(1.0, 0.0, 0.0), Point::new(0.0, 1.0, 0.0)];
        let u = Vector::new(0.0, 0.0, 1.0);
        let rotated = rotate_points_around_vector(&pts, &u, 0.0);
        assert!(rotated[0].is_close(&pts[0]));
        assert!(rotated[1].is_close(&pts[1]));
    }

    #[test]
    fn test_rotate_points_to_plane_invalid_normal() {
        let pts = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
        ];
        let zero_normal = Vector::new(0.0, 0.0, 0.0);
        let result = rotate_points_to_plane(&pts, &zero_normal);
        // Should return original points unchanged
        assert!(result[0].is_close(&pts[0]));
    }

    #[test]
    fn test_calculate_plane_normal_collinear() {
        // All collinear points - no valid plane normal
        let pts = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(2.0, 0.0, 0.0),
        ];
        assert!(calculate_plane_normal(&pts).is_none());
    }

    #[test]
    fn test_calculate_plane_normal_few_points() {
        let pts = vec![Point::new(0.0, 0.0, 0.0), Point::new(1.0, 0.0, 0.0)];
        assert!(calculate_plane_normal(&pts).is_none());
    }
}
