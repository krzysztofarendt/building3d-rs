use crate::Point;
use crate::Vector;
use crate::geom;
use crate::geom::point::check::are_point_sequences_close_rot;
use crate::geom::point::check::are_points_coplanar;
use crate::geom::triangles::{TriangleIndex, triangulate};
use crate::random_id;
use anyhow::{Result, anyhow};
use std::fmt;

#[derive(Debug, Clone)]
pub struct Polygon {
    /// Polygon name
    pub name: String,
    /// Polygon points
    pub pts: Vec<Point>,
    /// Normal vector
    pub vn: Vector,
    /// Triangles (flat list with triangle indices)
    tri: Vec<TriangleIndex>,
    /// Unique identifier
    uid: String,
    // TODO: Add parent (wall)
}

impl Polygon {
    /// Returns a new polygon.
    ///
    /// The normal vector is optional. If it is provided, its validity isn't checked.
    /// If it isn't provided, the normal will be calculate based on the first corner
    /// defined by points: last (-1), first (0), second (1).
    pub fn new(name: String, pts: Vec<Point>, normal: Option<Vector>) -> Result<Self> {
        if !are_points_coplanar(&pts) || pts.len() < 3 {
            return Err(anyhow!("Polygon points are invalid."));
        }

        let name = geom::validate_name(&name)?;

        // Assign normal vector. If it is provided, take it. If None is passed
        // then calculate it from the points of the first corner: last, 0, 1.
        // If the calculated normal is None, the corner points are collinear so go panic.
        // The first corner must be convex.
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

        Ok(Self {
            name,
            pts,
            vn,
            tri,
            uid: random_id(),
        })
    }
}

impl PartialEq for Polygon {
    fn eq(&self, other: &Self) -> bool {
        are_point_sequences_close_rot(&self.pts, &other.pts)
    }
}

impl Eq for Polygon {}

impl fmt::Display for Polygon {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let prec = f.precision().unwrap_or(2); // Default 2 decimals
        write!(f, "Polygon(\"{}\", ", self.name)?;
        for (i, p) in self.pts.iter().enumerate() {
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
        let poly_a = Polygon::new("a".to_string(), pts_a, None)?;
        let poly_b = Polygon::new("b".to_string(), pts_b, None)?;
        let poly_c = Polygon::new("c".to_string(), pts_c, None)?;
        assert!(poly_a == poly_b);
        assert!(poly_a != poly_c);

        Ok(())
    }
}
