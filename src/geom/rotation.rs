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
    if !u.length().is_close(1.) {
        panic!("rotation_matrix() requires u to be a unit vector");
    }

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
    let rot = rotation_matrix(u, phi);

    rotate_points(pts, &rot.t())
}

pub fn rotate_points_to_plane() {
    unimplemented!(); // NOTE: This function is not used anywhere in Building3D
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
}
