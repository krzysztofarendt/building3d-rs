//! Ray casting infrastructure.
//!
//! This module provides a Ray struct and ray-geometry intersection tests
//! for visibility analysis and other ray-based operations.

use crate::{Point, Polygon, Vector};

/// A ray defined by an origin point and a direction vector.
#[derive(Debug, Clone, Copy)]
pub struct Ray {
    /// Origin point of the ray
    pub origin: Point,
    /// Direction vector (should be normalized for distance calculations)
    pub direction: Vector,
}

impl Ray {
    /// Creates a new ray from origin point and direction vector.
    ///
    /// The direction vector is automatically normalized.
    pub fn new(origin: Point, direction: Vector) -> Option<Self> {
        let normalized = direction.normalize().ok()?;
        Some(Self {
            origin,
            direction: normalized,
        })
    }

    /// Creates a ray from two points (origin to target).
    pub fn from_points(origin: Point, target: Point) -> Option<Self> {
        let direction = target - origin;
        Self::new(origin, direction)
    }

    /// Returns the point along the ray at parameter t.
    ///
    /// point = origin + t * direction
    pub fn point_at(&self, t: f64) -> Point {
        self.origin + self.direction * t
    }

    /// Calculates the intersection of this ray with a polygon.
    ///
    /// Returns `Some((t, point))` if the ray intersects the polygon,
    /// where `t` is the ray parameter and `point` is the intersection point.
    /// Returns `None` if no intersection exists.
    ///
    /// Only returns intersections where t > 0 (in front of ray origin).
    pub fn intersect_polygon(&self, polygon: &Polygon) -> Option<(f64, Point)> {
        // Get plane coefficients
        let (a, b, c, d) = polygon.plane_coefficients();
        let plane_normal = Vector::new(a, b, c);

        // Check if ray is parallel to plane
        let denom = plane_normal.dot(&self.direction);
        if denom.abs() < 1e-10 {
            return None; // Ray parallel to plane
        }

        // Calculate intersection parameter t
        // Plane equation: a*x + b*y + c*z + d = 0
        // Ray: P = origin + t * direction
        // Substitute: a*(ox + t*dx) + b*(oy + t*dy) + c*(oz + t*dz) + d = 0
        // Solve for t: t = -(a*ox + b*oy + c*oz + d) / (a*dx + b*dy + c*dz)
        let origin_dot = a * self.origin.x + b * self.origin.y + c * self.origin.z + d;
        let t = -origin_dot / denom;

        // Only consider intersections in front of ray (t > 0)
        // Use small epsilon to avoid self-intersection at origin
        if t < 1e-10 {
            return None;
        }

        // Calculate intersection point
        let intersection_point = self.point_at(t);

        // Check if intersection point is inside the polygon
        if polygon.is_point_inside(intersection_point, true) {
            Some((t, intersection_point))
        } else {
            None
        }
    }

    /// Calculates the intersection of this ray with multiple polygons.
    ///
    /// Returns the closest intersection (smallest positive t) if any exists.
    pub fn intersect_polygons(&self, polygons: &[&Polygon]) -> Option<(f64, Point, usize)> {
        let mut closest: Option<(f64, Point, usize)> = None;

        for (idx, polygon) in polygons.iter().enumerate() {
            if let Some((t, point)) = self.intersect_polygon(polygon) {
                match &closest {
                    None => closest = Some((t, point, idx)),
                    Some((closest_t, _, _)) if t < *closest_t => {
                        closest = Some((t, point, idx));
                    }
                    _ => {}
                }
            }
        }

        closest
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use anyhow::Result;

    fn make_xy_square() -> Result<Polygon> {
        let pts = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(2.0, 0.0, 0.0),
            Point::new(2.0, 2.0, 0.0),
            Point::new(0.0, 2.0, 0.0),
        ];
        Polygon::new("square", pts, None)
    }

    #[test]
    fn test_ray_creation() {
        let ray = Ray::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0));
        assert!(ray.is_some());

        // Zero direction should fail
        let ray = Ray::new(Point::new(0.0, 0.0, 0.0), Vector::new(0.0, 0.0, 0.0));
        assert!(ray.is_none());
    }

    #[test]
    fn test_ray_from_points() {
        let ray = Ray::from_points(Point::new(0.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0));
        assert!(ray.is_some());
    }

    #[test]
    fn test_ray_point_at() {
        let ray = Ray::new(Point::new(0.0, 0.0, 0.0), Vector::new(1.0, 0.0, 0.0)).unwrap();

        let p = ray.point_at(5.0);
        assert!(p.is_close(&Point::new(5.0, 0.0, 0.0)));
    }

    #[test]
    fn test_ray_polygon_intersection() -> Result<()> {
        let polygon = make_xy_square()?;

        // Ray pointing at polygon from below
        let ray = Ray::new(Point::new(1.0, 1.0, -5.0), Vector::new(0.0, 0.0, 1.0)).unwrap();

        let result = ray.intersect_polygon(&polygon);
        assert!(result.is_some());

        let (t, point) = result.unwrap();
        assert!((t - 5.0).abs() < 1e-6);
        assert!(point.is_close(&Point::new(1.0, 1.0, 0.0)));

        Ok(())
    }

    #[test]
    fn test_ray_misses_polygon() -> Result<()> {
        let polygon = make_xy_square()?;

        // Ray pointing away from polygon
        let ray = Ray::new(Point::new(1.0, 1.0, -5.0), Vector::new(0.0, 0.0, -1.0)).unwrap();

        let result = ray.intersect_polygon(&polygon);
        assert!(result.is_none());

        Ok(())
    }

    #[test]
    fn test_ray_parallel_to_polygon() -> Result<()> {
        let polygon = make_xy_square()?;

        // Ray parallel to polygon plane
        let ray = Ray::new(Point::new(1.0, 1.0, 1.0), Vector::new(1.0, 0.0, 0.0)).unwrap();

        let result = ray.intersect_polygon(&polygon);
        assert!(result.is_none());

        Ok(())
    }

    #[test]
    fn test_ray_outside_polygon_bounds() -> Result<()> {
        let polygon = make_xy_square()?;

        // Ray hits plane but outside polygon
        let ray = Ray::new(Point::new(10.0, 10.0, -5.0), Vector::new(0.0, 0.0, 1.0)).unwrap();

        let result = ray.intersect_polygon(&polygon);
        assert!(result.is_none());

        Ok(())
    }

    #[test]
    fn test_ray_intersect_multiple_polygons() -> Result<()> {
        // Create two parallel squares at different z levels
        let pts1 = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(2.0, 0.0, 0.0),
            Point::new(2.0, 2.0, 0.0),
            Point::new(0.0, 2.0, 0.0),
        ];
        let poly1 = Polygon::new("z0", pts1, None)?;

        let pts2 = vec![
            Point::new(0.0, 0.0, 5.0),
            Point::new(2.0, 0.0, 5.0),
            Point::new(2.0, 2.0, 5.0),
            Point::new(0.0, 2.0, 5.0),
        ];
        let poly2 = Polygon::new("z5", pts2, None)?;

        let polygons: Vec<&Polygon> = vec![&poly1, &poly2];

        // Ray from below, should hit poly1 first
        let ray = Ray::new(Point::new(1.0, 1.0, -2.0), Vector::new(0.0, 0.0, 1.0)).unwrap();

        let result = ray.intersect_polygons(&polygons);
        assert!(result.is_some());

        let (t, point, idx) = result.unwrap();
        assert_eq!(idx, 0); // Should hit poly1 (z=0) first
        assert!((t - 2.0).abs() < 1e-6);
        assert!(point.is_close(&Point::new(1.0, 1.0, 0.0)));

        Ok(())
    }
}
