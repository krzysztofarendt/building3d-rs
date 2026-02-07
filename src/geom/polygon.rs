use crate::HasName;
use crate::Point;
use crate::UID;
use crate::Vector;
use crate::geom;
use crate::geom::point::check::are_point_sequences_close_rot;
use crate::geom::point::check::are_points_collinear;
use crate::geom::point::check::are_points_coplanar;
use crate::geom::rotation::rotate_points_around_vector;
use crate::geom::triangles::{triangulate, triangulate_with_holes};
use crate::{HasMesh, Mesh};
use anyhow::{Result, anyhow};
use serde::{Deserialize, Serialize};
use std::fmt;

pub mod boolean;
pub mod containment;
pub mod relations;
pub mod slice;

#[derive(Debug, Clone, Serialize, Deserialize)]
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
    /// Original outer boundary vertices (before hole merging)
    #[serde(default)]
    outer_boundary: Vec<Point>,
    /// Hole boundaries (opposite winding to outer)
    #[serde(default)]
    holes: Vec<Vec<Point>>,
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
        let vn = match normal {
            Some(v) => v
                .normalize()
                .map_err(|_| anyhow!("Normal vector invalid."))?,
            None => {
                let last = pts.len() - 1;
                match Vector::normal(pts[last], pts[0], pts[1]) {
                    Ok(v) => v,
                    Err(_) => return Err(anyhow!("Normal vector invalid.")),
                }
            }
        };

        let outer_boundary = pts.clone();
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
            outer_boundary,
            holes: vec![],
        })
    }

    /// Creates a polygon with holes.
    ///
    /// `outer` is the outer boundary (counter-clockwise winding w.r.t. normal).
    /// `holes` is a list of hole boundaries (clockwise winding, i.e. opposite to outer).
    pub fn with_holes(
        name: &str,
        mut outer: Vec<Point>,
        holes: Vec<Vec<Point>>,
        normal: Option<Vector>,
    ) -> Result<Self> {
        if !are_points_coplanar(&outer) || outer.len() < 3 {
            return Err(anyhow!("Polygon outer boundary points are invalid."));
        }

        let name = geom::validate_name(name)?;

        let vn = match normal {
            Some(v) => v
                .normalize()
                .map_err(|_| anyhow!("Normal vector invalid."))?,
            None => {
                let last = outer.len() - 1;
                match Vector::normal(outer[last], outer[0], outer[1]) {
                    Ok(v) => v,
                    Err(_) => return Err(anyhow!("Normal vector invalid.")),
                }
            }
        };

        if Self::ring_signed_area(&outer, &vn).abs() < crate::geom::EPS {
            return Err(anyhow!("Polygon outer boundary is degenerate."));
        }

        // Normalize outer winding to CCW w.r.t. normal
        if Self::ring_signed_area(&outer, &vn) < 0.0 {
            outer.reverse();
        }

        // Validate & normalize holes
        let plane_d = -vn.dot(&Vector::from_a_point(outer[0]));
        let mut hole_boundaries: Vec<Vec<Point>> = Vec::with_capacity(holes.len());
        for (i, mut hole) in holes.into_iter().enumerate() {
            if hole.len() < 3 {
                return Err(anyhow!("Polygon hole boundary {} is invalid.", i));
            }
            if !are_points_coplanar(&hole) {
                return Err(anyhow!("Polygon hole boundary {} points are invalid.", i));
            }
            if are_points_collinear(&hole) {
                return Err(anyhow!("Polygon hole boundary {} is degenerate.", i));
            }
            for pt in &hole {
                let dist = (plane_d + vn.dot(&Vector::from_a_point(*pt))).abs();
                if dist > crate::geom::EPS {
                    return Err(anyhow!("Polygon hole boundary {} is not coplanar.", i));
                }
            }

            // Normalize hole winding to CW (opposite of outer) w.r.t. normal
            let signed = Self::ring_signed_area(&hole, &vn);
            if signed.abs() < crate::geom::EPS {
                return Err(anyhow!("Polygon hole boundary {} is degenerate.", i));
            }
            if signed > 0.0 {
                hole.reverse();
            }
            hole_boundaries.push(hole);
        }

        let outer_boundary = outer.clone();

        let (pts, tri) = triangulate_with_holes(outer, hole_boundaries.clone(), vn, 0)?;

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
            outer_boundary,
            holes: hole_boundaries,
        })
    }

    fn ring_signed_area(pts: &[Point], vn: &Vector) -> f64 {
        let n = pts.len();
        if n < 3 {
            return 0.0;
        }
        let mut cross_sum = Vector::new(0.0, 0.0, 0.0);
        for i in 0..n {
            let v1 = Vector::from_a_point(pts[i]);
            let v2 = Vector::from_a_point(pts[(i + 1) % n]);
            cross_sum = cross_sum + v1.cross(&v2);
        }
        0.5 * cross_sum.dot(vn)
    }

    pub fn uid(&self) -> &str {
        self.uid.as_str()
    }

    pub fn mesh_ref(&self) -> &Mesh {
        &self.mesh
    }

    // Copies, flips points, renames, resets parent
    pub fn flip(&self, new_name: &str) -> Result<Self> {
        if self.holes.is_empty() {
            let mut vertices = self.effective_outer_boundary().to_vec();
            vertices.reverse();
            Self::new(new_name, vertices, None)
        } else {
            let mut outer = self.outer_boundary.clone();
            outer.reverse();
            let holes: Vec<Vec<Point>> = self
                .holes
                .iter()
                .map(|h| h.iter().rev().copied().collect())
                .collect();
            Self::with_holes(new_name, outer, holes, None)
        }
    }

    pub fn rotate(&mut self, angle: f64, rot_vec: &Vector) {
        self.mesh.vertices = rotate_points_around_vector(&self.mesh.vertices, rot_vec, angle);
        if !self.outer_boundary.is_empty() {
            self.outer_boundary = rotate_points_around_vector(&self.outer_boundary, rot_vec, angle);
        }
        for hole in &mut self.holes {
            *hole = rotate_points_around_vector(hole, rot_vec, angle);
        }
    }

    pub fn translate(&mut self, vec: &Vector) {
        for pt in self.mesh.vertices.iter_mut() {
            *pt += vec;
        }
        for pt in self.outer_boundary.iter_mut() {
            *pt += vec;
        }
        for hole in &mut self.holes {
            for pt in hole.iter_mut() {
                *pt += vec;
            }
        }
    }

    /// Returns the vertices of the polygon.
    pub fn vertices(&self) -> &[Point] {
        &self.mesh.vertices
    }

    /// Returns the triangulation of the polygon.
    pub fn triangles(&self) -> Option<&Vec<crate::TriangleIndex>> {
        self.mesh.faces.as_ref()
    }

    /// Returns the edges of the outer boundary as pairs of consecutive vertices.
    ///
    /// The edges form a closed loop, so the last edge connects the last vertex
    /// back to the first vertex.
    pub fn edges(&self) -> Vec<(Point, Point)> {
        let pts = self.effective_outer_boundary();
        let n = pts.len();
        if n < 2 {
            return vec![];
        }

        let mut edges = Vec::with_capacity(n);
        for i in 0..n {
            edges.push((pts[i], pts[(i + 1) % n]));
        }
        edges
    }

    /// Calculates the area of the polygon.
    ///
    /// For polygons with holes, this is the outer area minus the sum of hole areas.
    /// Uses the cross product method: sum of cross products of consecutive edges,
    /// dotted with the normal vector, divided by 2.
    pub fn area(&self) -> f64 {
        let outer_area = Self::shoelace_area(self.effective_outer_boundary(), &self.vn);
        let hole_area: f64 = self
            .holes
            .iter()
            .map(|h| Self::shoelace_area(h, &self.vn))
            .sum();
        outer_area - hole_area
    }

    /// Computes the area of a closed polygon ring using the 3D shoelace formula.
    fn shoelace_area(pts: &[Point], vn: &Vector) -> f64 {
        let n = pts.len();
        if n < 3 {
            return 0.0;
        }

        let mut cross_sum = Vector::new(0.0, 0.0, 0.0);

        for i in 0..n {
            let v1 = Vector::from_a_point(pts[i]);
            let v2 = Vector::from_a_point(pts[(i + 1) % n]);
            cross_sum = cross_sum + v1.cross(&v2);
        }

        0.5 * cross_sum.dot(vn).abs()
    }

    /// Calculates the centroid of the polygon.
    ///
    /// For triangulated polygons, this is the weighted average of triangle centroids,
    /// where each triangle's weight is its area.
    pub fn centroid(&self) -> Point {
        let pts = &self.mesh.vertices;
        let tri = match &self.mesh.faces {
            Some(faces) => faces,
            None => {
                // Fallback: average of vertices if no triangulation
                if pts.is_empty() {
                    return Point::new(0.0, 0.0, 0.0);
                }
                let sum_x: f64 = pts.iter().map(|p| p.x).sum();
                let sum_y: f64 = pts.iter().map(|p| p.y).sum();
                let sum_z: f64 = pts.iter().map(|p| p.z).sum();
                let n = pts.len() as f64;
                return Point::new(sum_x / n, sum_y / n, sum_z / n);
            }
        };

        if tri.is_empty() {
            // No triangles - return average of vertices
            if pts.is_empty() {
                return Point::new(0.0, 0.0, 0.0);
            }
            let sum_x: f64 = pts.iter().map(|p| p.x).sum();
            let sum_y: f64 = pts.iter().map(|p| p.y).sum();
            let sum_z: f64 = pts.iter().map(|p| p.z).sum();
            let n = pts.len() as f64;
            return Point::new(sum_x / n, sum_y / n, sum_z / n);
        }

        // Weighted average of triangle centroids
        let mut total_area = 0.0;
        let mut weighted_x = 0.0;
        let mut weighted_y = 0.0;
        let mut weighted_z = 0.0;

        for t in tri {
            let p0 = pts[t.0];
            let p1 = pts[t.1];
            let p2 = pts[t.2];

            // Triangle centroid
            let cx = (p0.x + p1.x + p2.x) / 3.0;
            let cy = (p0.y + p1.y + p2.y) / 3.0;
            let cz = (p0.z + p1.z + p2.z) / 3.0;

            // Triangle area (half the magnitude of cross product)
            let v1 = p1 - p0;
            let v2 = p2 - p0;
            let area = 0.5 * v1.cross(&v2).length();

            total_area += area;
            weighted_x += cx * area;
            weighted_y += cy * area;
            weighted_z += cz * area;
        }

        if total_area < 1e-15 {
            // Degenerate polygon - return average of vertices
            let sum_x: f64 = pts.iter().map(|p| p.x).sum();
            let sum_y: f64 = pts.iter().map(|p| p.y).sum();
            let sum_z: f64 = pts.iter().map(|p| p.z).sum();
            let n = pts.len() as f64;
            return Point::new(sum_x / n, sum_y / n, sum_z / n);
        }

        Point::new(
            weighted_x / total_area,
            weighted_y / total_area,
            weighted_z / total_area,
        )
    }

    /// Returns the plane coefficients (a, b, c, d) for the plane equation ax + by + cz + d = 0.
    ///
    /// The coefficients (a, b, c) are the components of the normal vector.
    /// The coefficient d is computed from any point on the plane.
    pub fn plane_coefficients(&self) -> (f64, f64, f64, f64) {
        let a = self.vn.dx;
        let b = self.vn.dy;
        let c = self.vn.dz;

        // d = -(ax + by + cz) for any point on the plane
        let p0 = if !self.mesh.vertices.is_empty() {
            self.mesh.vertices[0]
        } else {
            Point::new(0.0, 0.0, 0.0)
        };

        let d = -(a * p0.x + b * p0.y + c * p0.z);

        (a, b, c, d)
    }

    /// Returns true if this polygon has holes.
    pub fn has_holes(&self) -> bool {
        !self.holes.is_empty()
    }

    /// Returns a reference to the hole boundaries.
    pub fn holes(&self) -> &[Vec<Point>] {
        &self.holes
    }

    /// Returns edges for each hole as pairs of consecutive vertices.
    pub fn hole_edges(&self) -> Vec<Vec<(Point, Point)>> {
        self.holes
            .iter()
            .map(|h| {
                let n = h.len();
                if n < 2 {
                    return vec![];
                }
                let mut edges = Vec::with_capacity(n);
                for i in 0..n {
                    edges.push((h[i], h[(i + 1) % n]));
                }
                edges
            })
            .collect()
    }

    /// Returns the effective outer boundary vertices.
    ///
    /// If `outer_boundary` is populated (new-style), returns it.
    /// Otherwise falls back to `mesh.vertices` (legacy deserialized data without holes).
    fn effective_outer_boundary(&self) -> &[Point] {
        if !self.outer_boundary.is_empty() {
            &self.outer_boundary
        } else {
            &self.mesh.vertices
        }
    }

    /// Checks if a point lies inside the polygon.
    ///
    /// If `boundary_in` is true, points on the boundary (edges or vertices) are considered inside.
    /// The point must lie on the polygon's plane to be considered inside.
    /// For polygons with holes, a point inside a hole is considered outside.
    pub fn is_point_inside(&self, ptest: Point, boundary_in: bool) -> bool {
        let tri = match &self.mesh.faces {
            Some(faces) => faces,
            None => return false,
        };

        containment::is_point_inside_polygon(
            ptest,
            &self.mesh.vertices,
            tri,
            &self.vn,
            boundary_in,
            self.effective_outer_boundary(),
            &self.holes,
        )
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
        let pts = self.effective_outer_boundary();
        for (i, p) in pts.iter().enumerate() {
            if i > 0 {
                write!(f, ", ")?;
            }
            write!(f, "{:.prec$}", p, prec = prec)?;
        }
        if !self.holes.is_empty() {
            write!(f, ", holes={}", self.holes.len())?;
        }
        write!(f, ")")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::IsClose;

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

    #[test]
    fn test_area_unit_square() -> Result<()> {
        let pts = vec![
            Point::new(0., 0., 0.),
            Point::new(1., 0., 0.),
            Point::new(1., 1., 0.),
            Point::new(0., 1., 0.),
        ];
        let poly = Polygon::new("square", pts, None)?;
        let area = poly.area();
        assert!(area.is_close(1.0));
        Ok(())
    }

    #[test]
    fn test_area_rectangle() -> Result<()> {
        let pts = vec![
            Point::new(0., 0., 0.),
            Point::new(2., 0., 0.),
            Point::new(2., 3., 0.),
            Point::new(0., 3., 0.),
        ];
        let poly = Polygon::new("rect", pts, None)?;
        let area = poly.area();
        assert!(area.is_close(6.0));
        Ok(())
    }

    #[test]
    fn test_area_triangle() -> Result<()> {
        let pts = vec![
            Point::new(0., 0., 0.),
            Point::new(2., 0., 0.),
            Point::new(1., 2., 0.),
        ];
        let poly = Polygon::new("tri", pts, None)?;
        let area = poly.area();
        // Area = 0.5 * base * height = 0.5 * 2 * 2 = 2.0
        assert!(area.is_close(2.0));
        Ok(())
    }

    #[test]
    fn test_centroid_square() -> Result<()> {
        let pts = vec![
            Point::new(0., 0., 0.),
            Point::new(2., 0., 0.),
            Point::new(2., 2., 0.),
            Point::new(0., 2., 0.),
        ];
        let poly = Polygon::new("square", pts, None)?;
        let centroid = poly.centroid();
        assert!(centroid.is_close(&Point::new(1.0, 1.0, 0.0)));
        Ok(())
    }

    #[test]
    fn test_centroid_triangle() -> Result<()> {
        let pts = vec![
            Point::new(0., 0., 0.),
            Point::new(3., 0., 0.),
            Point::new(0., 3., 0.),
        ];
        let poly = Polygon::new("tri", pts, None)?;
        let centroid = poly.centroid();
        // Centroid of triangle is at (sum_x/3, sum_y/3, sum_z/3)
        assert!(centroid.is_close(&Point::new(1.0, 1.0, 0.0)));
        Ok(())
    }

    #[test]
    fn test_edges() -> Result<()> {
        let pts = vec![
            Point::new(0., 0., 0.),
            Point::new(1., 0., 0.),
            Point::new(1., 1., 0.),
        ];
        let poly = Polygon::new("tri", pts, None)?;
        let edges = poly.edges();

        assert_eq!(edges.len(), 3);

        // First edge: (0,0,0) -> (1,0,0)
        assert!(edges[0].0.is_close(&Point::new(0., 0., 0.)));
        assert!(edges[0].1.is_close(&Point::new(1., 0., 0.)));

        // Second edge: (1,0,0) -> (1,1,0)
        assert!(edges[1].0.is_close(&Point::new(1., 0., 0.)));
        assert!(edges[1].1.is_close(&Point::new(1., 1., 0.)));

        // Third edge (closing): (1,1,0) -> (0,0,0)
        assert!(edges[2].0.is_close(&Point::new(1., 1., 0.)));
        assert!(edges[2].1.is_close(&Point::new(0., 0., 0.)));

        Ok(())
    }

    #[test]
    fn test_plane_coefficients() -> Result<()> {
        // Polygon in XY plane (z = 0)
        let pts = vec![
            Point::new(0., 0., 0.),
            Point::new(1., 0., 0.),
            Point::new(1., 1., 0.),
            Point::new(0., 1., 0.),
        ];
        let poly = Polygon::new("square", pts, None)?;
        let (a, b, c, d) = poly.plane_coefficients();

        // Normal should be (0, 0, 1), so a=0, b=0, c=1
        // d should be 0 since plane passes through origin
        assert!(a.is_close(0.0));
        assert!(b.is_close(0.0));
        assert!(c.is_close(1.0));
        assert!(d.is_close(0.0));

        Ok(())
    }

    #[test]
    fn test_plane_coefficients_offset() -> Result<()> {
        // Polygon in XY plane at z = 5
        let pts = vec![
            Point::new(0., 0., 5.),
            Point::new(1., 0., 5.),
            Point::new(1., 1., 5.),
            Point::new(0., 1., 5.),
        ];
        let poly = Polygon::new("square", pts, None)?;
        let (a, b, c, d) = poly.plane_coefficients();

        // Normal should be (0, 0, 1)
        // Plane equation: 0*x + 0*y + 1*z + d = 0, so z = -d
        // Since z = 5, d = -5
        assert!(a.is_close(0.0));
        assert!(b.is_close(0.0));
        assert!(c.is_close(1.0));
        assert!(d.is_close(-5.0));

        Ok(())
    }

    #[test]
    fn test_is_point_inside_square() -> Result<()> {
        let pts = vec![
            Point::new(0., 0., 0.),
            Point::new(1., 0., 0.),
            Point::new(1., 1., 0.),
            Point::new(0., 1., 0.),
        ];
        let poly = Polygon::new("square", pts, None)?;

        // Inside
        assert!(poly.is_point_inside(Point::new(0.5, 0.5, 0.0), true));
        assert!(poly.is_point_inside(Point::new(0.5, 0.5, 0.0), false));

        // Outside
        assert!(!poly.is_point_inside(Point::new(1.5, 0.5, 0.0), true));
        assert!(!poly.is_point_inside(Point::new(1.5, 0.5, 0.0), false));

        // On boundary (vertex)
        assert!(poly.is_point_inside(Point::new(0., 0., 0.0), true));
        assert!(!poly.is_point_inside(Point::new(0., 0., 0.0), false));

        // On boundary (edge)
        assert!(poly.is_point_inside(Point::new(0.5, 0., 0.0), true));
        assert!(!poly.is_point_inside(Point::new(0.5, 0., 0.0), false));

        // Above plane
        assert!(!poly.is_point_inside(Point::new(0.5, 0.5, 1.0), true));

        Ok(())
    }

    #[test]
    fn test_vertices_and_triangles() -> Result<()> {
        let pts = vec![
            Point::new(0., 0., 0.),
            Point::new(1., 0., 0.),
            Point::new(1., 1., 0.),
            Point::new(0., 1., 0.),
        ];
        let poly = Polygon::new("square", pts, None)?;

        assert_eq!(poly.vertices().len(), 4);
        assert!(poly.triangles().is_some());
        assert_eq!(poly.triangles().unwrap().len(), 2); // Square = 2 triangles

        Ok(())
    }

    #[test]
    fn test_polygon_with_holes_area() -> Result<()> {
        // 4x4 outer square (CCW), area = 16
        let outer = vec![
            Point::new(0., 0., 0.),
            Point::new(4., 0., 0.),
            Point::new(4., 4., 0.),
            Point::new(0., 4., 0.),
        ];
        // 1x1 hole (CW), area = 1
        let hole = vec![
            Point::new(1.5, 1.5, 0.),
            Point::new(1.5, 2.5, 0.),
            Point::new(2.5, 2.5, 0.),
            Point::new(2.5, 1.5, 0.),
        ];
        let poly = Polygon::with_holes("holed", outer, vec![hole], None)?;
        let area = poly.area();
        assert!(
            (area - 15.0).abs() < 0.01,
            "Area should be 16 - 1 = 15, got {}",
            area
        );
        Ok(())
    }

    #[test]
    fn test_polygon_with_holes_is_point_inside() -> Result<()> {
        let outer = vec![
            Point::new(0., 0., 0.),
            Point::new(4., 0., 0.),
            Point::new(4., 4., 0.),
            Point::new(0., 4., 0.),
        ];
        let hole = vec![
            Point::new(1.5, 1.5, 0.),
            Point::new(1.5, 2.5, 0.),
            Point::new(2.5, 2.5, 0.),
            Point::new(2.5, 1.5, 0.),
        ];
        let poly = Polygon::with_holes("holed", outer, vec![hole], None)?;

        // Point in solid part -> inside
        assert!(poly.is_point_inside(Point::new(0.5, 0.5, 0.0), false));

        // Point in hole -> outside
        assert!(!poly.is_point_inside(Point::new(2.0, 2.0, 0.0), false));

        // Point outside everything -> outside
        assert!(!poly.is_point_inside(Point::new(5.0, 5.0, 0.0), false));

        Ok(())
    }

    #[test]
    fn test_polygon_with_holes_edges() -> Result<()> {
        let outer = vec![
            Point::new(0., 0., 0.),
            Point::new(4., 0., 0.),
            Point::new(4., 4., 0.),
            Point::new(0., 4., 0.),
        ];
        let hole = vec![
            Point::new(1.5, 1.5, 0.),
            Point::new(1.5, 2.5, 0.),
            Point::new(2.5, 2.5, 0.),
            Point::new(2.5, 1.5, 0.),
        ];
        let poly = Polygon::with_holes("holed", outer, vec![hole], None)?;

        // edges() returns only outer boundary edges (4 for a square)
        assert_eq!(poly.edges().len(), 4);

        // hole_edges() returns hole edges
        let he = poly.hole_edges();
        assert_eq!(he.len(), 1);
        assert_eq!(he[0].len(), 4);

        assert!(poly.has_holes());
        assert_eq!(poly.holes().len(), 1);

        Ok(())
    }

    #[test]
    fn test_polygon_with_holes_serialization() -> Result<()> {
        let outer = vec![
            Point::new(0., 0., 0.),
            Point::new(4., 0., 0.),
            Point::new(4., 4., 0.),
            Point::new(0., 4., 0.),
        ];
        let hole = vec![
            Point::new(1.5, 1.5, 0.),
            Point::new(1.5, 2.5, 0.),
            Point::new(2.5, 2.5, 0.),
            Point::new(2.5, 1.5, 0.),
        ];
        let poly = Polygon::with_holes("holed", outer, vec![hole], None)?;

        let json = serde_json::to_string(&poly).unwrap();
        let restored: Polygon = serde_json::from_str(&json).unwrap();

        assert_eq!(restored.holes().len(), 1);
        assert!(restored.has_holes());
        assert!((restored.area() - 15.0).abs() < 0.01);

        Ok(())
    }

    #[test]
    fn test_polygon_with_holes_rejects_empty_hole() {
        let outer = vec![
            Point::new(0., 0., 0.),
            Point::new(4., 0., 0.),
            Point::new(4., 4., 0.),
            Point::new(0., 4., 0.),
        ];

        let err = Polygon::with_holes("holed", outer, vec![vec![]], None).unwrap_err();
        assert!(err.to_string().to_lowercase().contains("hole"));
    }

    #[test]
    fn test_polygon_with_holes_rejects_degenerate_hole() {
        let outer = vec![
            Point::new(0., 0., 0.),
            Point::new(4., 0., 0.),
            Point::new(4., 4., 0.),
            Point::new(0., 4., 0.),
        ];
        let hole = vec![Point::new(1.0, 1.0, 0.0), Point::new(2.0, 1.0, 0.0)];

        let err = Polygon::with_holes("holed", outer, vec![hole], None).unwrap_err();
        assert!(err.to_string().to_lowercase().contains("hole"));
    }

    #[test]
    fn test_polygon_backward_compat() -> Result<()> {
        // Simulate legacy JSON without outer_boundary and holes fields
        let poly = Polygon::new(
            "square",
            vec![
                Point::new(0., 0., 0.),
                Point::new(1., 0., 0.),
                Point::new(1., 1., 0.),
                Point::new(0., 1., 0.),
            ],
            None,
        )?;

        let json = serde_json::to_string(&poly).unwrap();

        // Remove outer_boundary and holes from JSON to simulate old format
        let v: serde_json::Value = serde_json::from_str(&json).unwrap();
        let mut map = v.as_object().unwrap().clone();
        map.remove("outer_boundary");
        map.remove("holes");
        let old_json = serde_json::to_string(&map).unwrap();

        let restored: Polygon = serde_json::from_str(&old_json).unwrap();

        // Should still work - effective_outer_boundary falls back to mesh.vertices
        assert!(!restored.has_holes());
        assert!(restored.area().is_close(1.0));
        assert_eq!(restored.edges().len(), 4);

        Ok(())
    }

    #[test]
    fn test_display_default_precision() -> Result<()> {
        let pts = vec![
            Point::new(0., 0., 0.),
            Point::new(1., 0., 0.),
            Point::new(0.5, 1., 0.),
        ];
        let poly = Polygon::new("tri", pts, None)?;
        let s = format!("{}", poly);
        assert!(s.starts_with("Polygon(\"tri\", "));
        assert!(s.ends_with(")"));
        Ok(())
    }

    #[test]
    fn test_display_custom_precision() -> Result<()> {
        let pts = vec![
            Point::new(0., 0., 0.),
            Point::new(1., 0., 0.),
            Point::new(0.5, 1., 0.),
        ];
        let poly = Polygon::new("tri", pts, None)?;
        let s = format!("{:.4}", poly);
        assert!(s.starts_with("Polygon(\"tri\", "));
        assert!(s.ends_with(")"));
        Ok(())
    }

    #[test]
    fn test_uid() -> Result<()> {
        let pts = vec![
            Point::new(0., 0., 0.),
            Point::new(1., 0., 0.),
            Point::new(0.5, 1., 0.),
        ];
        let poly = Polygon::new("tri", pts, None)?;
        let uid = poly.uid();
        assert!(!uid.is_empty());
        Ok(())
    }

    #[test]
    fn test_flip() -> Result<()> {
        let pts = vec![
            Point::new(0., 0., 0.),
            Point::new(1., 0., 0.),
            Point::new(1., 1., 0.),
            Point::new(0., 1., 0.),
        ];
        let poly = Polygon::new("square", pts, None)?;
        let flipped = poly.flip("flipped")?;

        assert_eq!(flipped.name, "flipped");
        // Normal should be reversed
        let dot = poly.vn.dot(&flipped.vn);
        assert!(dot.is_close(-1.0));
        // Area should be preserved
        assert!(poly.area().is_close(flipped.area()));

        Ok(())
    }
}
