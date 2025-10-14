use crate::HasName;
use crate::Point;
use crate::UID;
use crate::Vector;
use crate::geom;
use crate::geom::point::check::are_point_sequences_close_rot;
use crate::geom::point::check::are_points_coplanar;
use crate::geom::rotation::rotate_points_around_vector;
use crate::geom::triangles::triangulate;
use crate::{HasMesh, Mesh};
use anyhow::{Result, anyhow};
use std::fmt;

#[derive(Debug, Clone)]
pub struct Polygon {
    /// Polygon name
    pub name: String,
    /// Polygon mesh
    mesh: Mesh,
    /// Normal vector
    pub vn: Vector,
    /// Unique identifier of this polygon
    pub uid: UID,
    /// Unique identifier of the parent wall
    pub parent: Option<UID>,
}

impl HasName for Polygon {
    fn get_name(&self) -> &str {
        &self.name
    }
}

impl HasMesh for Polygon {
    fn copy_mesh(&self) -> Mesh {
        self.mesh.clone()
    }
}

impl Polygon {
    /// Returns a new polygon.
    ///
    /// The normal vector is optional. If it is provided, its validity isn't checked.
    /// If it isn't provided, the normal will be calculated based on the first corner
    /// defined by points: last (-1), first (0), second (1).
    pub fn new(name: &str, pts: Vec<Point>, normal: Option<Vector>) -> Result<Self> {
        if !are_points_coplanar(&pts) || pts.len() < 3 {
            return Err(anyhow!("Polygon points are invalid."));
        }

        let name = geom::validate_name(name)?;

        // Assign normal vector. If it is provided, take it. If None is passed
        // then calculate it from the points of the first corner: last, 0, 1.
        // If the calculated normal is None, the corner points are collinear so go panic.
        // The first corner must be convex.
        // TODO: Validate the normal vector if it is provided (check if it's orthogonal).
        let vn = match normal {
            Some(v) => v,
            None => {
                let last = pts.len() - 1;
                match Vector::normal(pts[last], pts[0], pts[1]) {
                    Ok(v) => v,
                    Err(_) => return Err(anyhow!("Normal vector invalid.")),
                }
            }
        };

        let (pts, tri) = triangulate(pts, vn, 0)?;

        let mesh = Mesh {
            vertices: pts,
            faces: Some(tri),
        };

        Ok(Self {
            name: name.to_string(),
            mesh,
            vn,
            uid: UID::new(),
            parent: None,
        })
    }

    pub fn uid(&self) -> &str {
        self.uid.as_str()
    }

    pub fn mesh_ref(&self) -> &Mesh {
        &self.mesh
    }

    // Copies, flips points, renames, resets parent
    pub fn flip(&self, new_name: &str) -> Result<Self> {
        let mut vertices = self.mesh.vertices.clone();
        vertices.reverse();

        Self::new(new_name, vertices, None)
    }

    pub fn rotate(&mut self, angle: f64, rot_vec: &Vector) {
        self.mesh.vertices = rotate_points_around_vector(&self.mesh.vertices, rot_vec, angle);
    }

    pub fn translate(&mut self, vec: &Vector) {
        for pt in self.mesh.vertices.iter_mut() {
            *pt += vec;
        }
    }
}

impl PartialEq for Polygon {
    fn eq(&self, other: &Self) -> bool {
        are_point_sequences_close_rot(&self.mesh.vertices, &other.mesh.vertices)
    }
}

impl Eq for Polygon {}

impl fmt::Display for Polygon {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let prec = f.precision().unwrap_or(2); // Default 2 decimals
        write!(f, "Polygon(\"{}\", ", self.name)?;
        for (i, p) in self.mesh.vertices.iter().enumerate() {
            if i > 0 {
                write!(f, ", ")?;
            }
            write!(f, "{:.prec$}", p, prec = prec)?;
        }
        write!(f, ")")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_eq() -> Result<()> {
        let pts_a = vec![
            Point::new(0., 0., 0.),
            Point::new(0.5, 0.5, 0.5),
            Point::new(1.0, 2.0, 3.0),
        ];
        let pts_b = vec![
            Point::new(0., 0., 0.),
            Point::new(0.5, 0.5, 0.5),
            Point::new(1.0, 2.0, 3.0),
        ];
        let pts_c = vec![
            Point::new(0., 0., 1.),
            Point::new(0.5, 0.5, 0.5),
            Point::new(1.0, 2.0, 3.0),
        ];
        let poly_a = Polygon::new("a", pts_a, None)?;
        let poly_b = Polygon::new("b", pts_b, None)?;
        let poly_c = Polygon::new("c", pts_c, None)?;
        assert!(poly_a == poly_b);
        assert!(poly_a != poly_c);

        Ok(())
    }
}
