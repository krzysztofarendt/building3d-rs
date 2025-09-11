use crate::Point;
use crate::Vector;
use crate::geom;
use crate::geom::point::check::are_points_coplanar;
use crate::random_id;
use anyhow::{Result, anyhow};
use std::fmt;

#[derive(Debug, Clone, PartialEq)]
pub struct Polygon {
    /// Polygon name
    pub name: String,
    /// Polygon points
    pub pts: Vec<Point>,
    /// Normal vector
    pub vn: Vector,
    /// Triangles (flat list with triangle indices)
    tri: Vec<usize>,
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
                    Some(v) => v,
                    None => return Err(anyhow!("Normal vector invalid.")),
                }
            }
        };
        Ok(Self {
            name,
            pts,
            vn,
            tri: Vec::new(), // TODO: Add triangulation
            uid: random_id(),
        })
    }
}

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
