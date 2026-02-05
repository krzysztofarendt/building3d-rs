//! Mesh quality analysis.
//!
//! This module provides functions to analyze the quality of triangulated meshes,
//! including metrics like aspect ratio, minimum/maximum angles, and edge lengths.

use crate::{Mesh, Point, Vector};

/// Quality metrics for a single triangle.
#[derive(Debug, Clone, Copy)]
pub struct TriangleQuality {
    /// Aspect ratio (longest edge / shortest edge). Ideal is 1.0.
    pub aspect_ratio: f64,
    /// Minimum interior angle in degrees. Ideal is 60° for equilateral.
    pub min_angle: f64,
    /// Maximum interior angle in degrees. Ideal is 60° for equilateral.
    pub max_angle: f64,
    /// Area of the triangle.
    pub area: f64,
    /// Shortest edge length.
    pub min_edge: f64,
    /// Longest edge length.
    pub max_edge: f64,
}

/// Overall quality statistics for a mesh.
#[derive(Debug, Clone)]
pub struct MeshQuality {
    /// Number of triangles analyzed.
    pub triangle_count: usize,
    /// Minimum aspect ratio across all triangles.
    pub min_aspect_ratio: f64,
    /// Maximum aspect ratio across all triangles.
    pub max_aspect_ratio: f64,
    /// Average aspect ratio.
    pub avg_aspect_ratio: f64,
    /// Minimum angle across all triangles (degrees).
    pub min_angle: f64,
    /// Maximum angle across all triangles (degrees).
    pub max_angle: f64,
    /// Total mesh area.
    pub total_area: f64,
    /// Number of degenerate triangles (near-zero area).
    pub degenerate_count: usize,
    /// Number of triangles with aspect ratio > threshold (default 10).
    pub poor_quality_count: usize,
}

impl Default for MeshQuality {
    fn default() -> Self {
        Self {
            triangle_count: 0,
            min_aspect_ratio: f64::INFINITY,
            max_aspect_ratio: 0.0,
            avg_aspect_ratio: 0.0,
            min_angle: f64::INFINITY,
            max_angle: 0.0,
            total_area: 0.0,
            degenerate_count: 0,
            poor_quality_count: 0,
        }
    }
}

/// Analyzes the quality of a single triangle.
///
/// # Arguments
/// * `p0`, `p1`, `p2` - The three vertices of the triangle
///
/// # Returns
/// Quality metrics for the triangle.
pub fn analyze_triangle(p0: Point, p1: Point, p2: Point) -> TriangleQuality {
    // Calculate edge vectors and lengths
    let e0 = p1 - p0;
    let e1 = p2 - p1;
    let e2 = p0 - p2;

    let len0 = e0.length();
    let len1 = e1.length();
    let len2 = e2.length();

    // Edge length statistics
    let min_edge = len0.min(len1).min(len2);
    let max_edge = len0.max(len1).max(len2);

    // Aspect ratio
    let aspect_ratio = if min_edge > 1e-10 {
        max_edge / min_edge
    } else {
        f64::INFINITY
    };

    // Calculate angles using dot products
    // Angle at p0: between edges e0 and -e2
    // Angle at p1: between edges e1 and -e0
    // Angle at p2: between edges e2 and -e1
    let angle0 = angle_between_vectors(&e0, &(e2 * -1.0));
    let angle1 = angle_between_vectors(&e1, &(e0 * -1.0));
    let angle2 = angle_between_vectors(&e2, &(e1 * -1.0));

    let min_angle = angle0.min(angle1).min(angle2);
    let max_angle = angle0.max(angle1).max(angle2);

    // Calculate area using cross product
    let cross = e0.cross(&(e2 * -1.0));
    let area = cross.length() / 2.0;

    TriangleQuality {
        aspect_ratio,
        min_angle,
        max_angle,
        area,
        min_edge,
        max_edge,
    }
}

/// Calculates the angle between two vectors in degrees.
fn angle_between_vectors(v1: &Vector, v2: &Vector) -> f64 {
    let len1 = v1.length();
    let len2 = v2.length();

    if len1 < 1e-10 || len2 < 1e-10 {
        return 0.0;
    }

    let cos_angle = (v1.dot(v2) / (len1 * len2)).clamp(-1.0, 1.0);
    cos_angle.acos().to_degrees()
}

/// Analyzes the quality of a mesh.
///
/// # Arguments
/// * `mesh` - The mesh to analyze
/// * `poor_quality_threshold` - Aspect ratio threshold for "poor quality" triangles (default: 10.0)
///
/// # Returns
/// Overall quality statistics for the mesh.
pub fn analyze_mesh(mesh: &Mesh, poor_quality_threshold: f64) -> MeshQuality {
    let faces = match mesh.faces() {
        Some(f) => f,
        None => return MeshQuality::default(),
    };

    let vertices = mesh.vertices();
    let mut quality = MeshQuality::default();
    let mut total_aspect_ratio = 0.0;

    for face in faces {
        let p0 = vertices[face.0];
        let p1 = vertices[face.1];
        let p2 = vertices[face.2];

        let tri_quality = analyze_triangle(p0, p1, p2);

        quality.triangle_count += 1;
        quality.total_area += tri_quality.area;

        // Track aspect ratio
        quality.min_aspect_ratio = quality.min_aspect_ratio.min(tri_quality.aspect_ratio);
        quality.max_aspect_ratio = quality.max_aspect_ratio.max(tri_quality.aspect_ratio);
        total_aspect_ratio += tri_quality.aspect_ratio;

        // Track angles
        quality.min_angle = quality.min_angle.min(tri_quality.min_angle);
        quality.max_angle = quality.max_angle.max(tri_quality.max_angle);

        // Count degenerate triangles
        if tri_quality.area < 1e-10 {
            quality.degenerate_count += 1;
        }

        // Count poor quality triangles
        if tri_quality.aspect_ratio > poor_quality_threshold {
            quality.poor_quality_count += 1;
        }
    }

    // Calculate average
    if quality.triangle_count > 0 {
        quality.avg_aspect_ratio = total_aspect_ratio / quality.triangle_count as f64;
    }

    quality
}

/// Analyzes mesh quality with default poor quality threshold of 10.0.
pub fn analyze_mesh_default(mesh: &Mesh) -> MeshQuality {
    analyze_mesh(mesh, 10.0)
}

/// Checks if a mesh has acceptable quality.
///
/// # Arguments
/// * `mesh` - The mesh to check
/// * `max_aspect_ratio` - Maximum acceptable aspect ratio
/// * `min_angle` - Minimum acceptable angle in degrees
///
/// # Returns
/// `true` if all triangles meet the quality criteria.
pub fn is_mesh_quality_acceptable(mesh: &Mesh, max_aspect_ratio: f64, min_angle: f64) -> bool {
    let faces = match mesh.faces() {
        Some(f) => f,
        None => return true, // No faces = trivially acceptable
    };

    let vertices = mesh.vertices();

    for face in faces {
        let p0 = vertices[face.0];
        let p1 = vertices[face.1];
        let p2 = vertices[face.2];

        let quality = analyze_triangle(p0, p1, p2);

        if quality.aspect_ratio > max_aspect_ratio || quality.min_angle < min_angle {
            return false;
        }
    }

    true
}

/// Finds triangles with poor quality in a mesh.
///
/// # Arguments
/// * `mesh` - The mesh to analyze
/// * `max_aspect_ratio` - Maximum acceptable aspect ratio
/// * `min_angle` - Minimum acceptable angle in degrees
///
/// # Returns
/// Indices of triangles that fail the quality criteria.
pub fn find_poor_quality_triangles(
    mesh: &Mesh,
    max_aspect_ratio: f64,
    min_angle: f64,
) -> Vec<usize> {
    let faces = match mesh.faces() {
        Some(f) => f,
        None => return vec![],
    };

    let vertices = mesh.vertices();
    let mut poor_triangles = Vec::new();

    for (idx, face) in faces.iter().enumerate() {
        let p0 = vertices[face.0];
        let p1 = vertices[face.1];
        let p2 = vertices[face.2];

        let quality = analyze_triangle(p0, p1, p2);

        if quality.aspect_ratio > max_aspect_ratio || quality.min_angle < min_angle {
            poor_triangles.push(idx);
        }
    }

    poor_triangles
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geom::triangles::TriangleIndex;

    #[test]
    fn test_analyze_equilateral_triangle() {
        // Equilateral triangle with side length 2
        let p0 = Point::new(0.0, 0.0, 0.0);
        let p1 = Point::new(2.0, 0.0, 0.0);
        let p2 = Point::new(1.0, 3.0_f64.sqrt(), 0.0);

        let quality = analyze_triangle(p0, p1, p2);

        // Aspect ratio should be 1.0 for equilateral
        assert!(
            (quality.aspect_ratio - 1.0).abs() < 0.01,
            "Equilateral triangle should have aspect ratio ~1.0, got {}",
            quality.aspect_ratio
        );

        // All angles should be 60 degrees
        assert!(
            (quality.min_angle - 60.0).abs() < 0.1,
            "Min angle should be ~60°, got {}",
            quality.min_angle
        );
        assert!(
            (quality.max_angle - 60.0).abs() < 0.1,
            "Max angle should be ~60°, got {}",
            quality.max_angle
        );

        // Area of equilateral triangle with side 2: sqrt(3)
        let expected_area = 3.0_f64.sqrt();
        assert!(
            (quality.area - expected_area).abs() < 0.01,
            "Area should be ~{}, got {}",
            expected_area,
            quality.area
        );
    }

    #[test]
    fn test_analyze_right_triangle() {
        // Right triangle: 3-4-5
        let p0 = Point::new(0.0, 0.0, 0.0);
        let p1 = Point::new(3.0, 0.0, 0.0);
        let p2 = Point::new(0.0, 4.0, 0.0);

        let quality = analyze_triangle(p0, p1, p2);

        // Aspect ratio: 5/3 ≈ 1.67
        assert!(
            (quality.aspect_ratio - 5.0 / 3.0).abs() < 0.01,
            "3-4-5 triangle aspect ratio should be ~1.67, got {}",
            quality.aspect_ratio
        );

        // Max angle should be 90 degrees
        assert!(
            (quality.max_angle - 90.0).abs() < 0.1,
            "Max angle should be ~90°, got {}",
            quality.max_angle
        );

        // Area = 0.5 * 3 * 4 = 6
        assert!(
            (quality.area - 6.0).abs() < 0.01,
            "Area should be 6.0, got {}",
            quality.area
        );
    }

    #[test]
    fn test_analyze_thin_triangle() {
        // Create a truly thin triangle with one very short edge
        let p0 = Point::new(0.0, 0.0, 0.0);
        let p1 = Point::new(10.0, 0.0, 0.0);
        let p2 = Point::new(0.0, 0.01, 0.0);

        let quality = analyze_triangle(p0, p1, p2);

        // Edge lengths: p0-p1=10, p1-p2≈10.000005, p2-p0=0.01
        // Aspect ratio: ~10/0.01 = 1000
        assert!(
            quality.aspect_ratio > 100.0,
            "Thin triangle should have high aspect ratio, got {}",
            quality.aspect_ratio
        );

        // Should have small minimum angle (at p0)
        assert!(
            quality.min_angle < 1.0,
            "Thin triangle should have small min angle, got {}",
            quality.min_angle
        );
    }

    #[test]
    fn test_analyze_mesh_quality() {
        // Create a simple mesh with two triangles
        let vertices = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(0.5, 0.866, 0.0), // Approximately equilateral
            Point::new(1.5, 0.866, 0.0),
        ];
        let faces = vec![TriangleIndex(0, 1, 2), TriangleIndex(1, 3, 2)];

        let mesh = Mesh::new(vertices, Some(faces));
        let quality = analyze_mesh_default(&mesh);

        assert_eq!(quality.triangle_count, 2);
        assert!(quality.total_area > 0.0);
        assert!(quality.degenerate_count == 0);
    }

    #[test]
    fn test_is_mesh_quality_acceptable() {
        // Good quality mesh
        let vertices = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(2.0, 0.0, 0.0),
            Point::new(1.0, 3.0_f64.sqrt(), 0.0),
        ];
        let faces = vec![TriangleIndex(0, 1, 2)];

        let mesh = Mesh::new(vertices, Some(faces));

        assert!(is_mesh_quality_acceptable(&mesh, 2.0, 30.0));
        assert!(is_mesh_quality_acceptable(&mesh, 1.5, 50.0));
    }

    #[test]
    fn test_find_poor_quality_triangles() {
        // Mix of good and poor triangles
        let vertices = vec![
            // Equilateral triangle (good)
            Point::new(0.0, 0.0, 0.0),
            Point::new(2.0, 0.0, 0.0),
            Point::new(1.0, 3.0_f64.sqrt(), 0.0),
            // Thin triangle (poor)
            Point::new(10.0, 0.0, 0.0),
            Point::new(20.0, 0.0, 0.0),
            Point::new(15.0, 0.1, 0.0),
        ];
        let faces = vec![
            TriangleIndex(0, 1, 2), // Good
            TriangleIndex(3, 4, 5), // Poor
        ];

        let mesh = Mesh::new(vertices, Some(faces));
        let poor = find_poor_quality_triangles(&mesh, 5.0, 10.0);

        assert_eq!(poor.len(), 1);
        assert_eq!(poor[0], 1);
    }

    #[test]
    fn test_degenerate_triangle() {
        // Collinear points (degenerate triangle)
        let p0 = Point::new(0.0, 0.0, 0.0);
        let p1 = Point::new(1.0, 0.0, 0.0);
        let p2 = Point::new(2.0, 0.0, 0.0);

        let quality = analyze_triangle(p0, p1, p2);

        // Degenerate triangle should have zero area
        assert!(
            quality.area < 1e-10,
            "Degenerate triangle should have zero area, got {}",
            quality.area
        );

        // For collinear points p0-p1-p2:
        // edge p0-p1 = 1, edge p1-p2 = 1, edge p2-p0 = 2
        // Aspect ratio = 2/1 = 2 (not infinite)
        // The "degenerate" quality comes from the zero area, not aspect ratio
        // For truly degenerate (coincident points), min_edge would be 0

        // Test with coincident points
        let p0 = Point::new(0.0, 0.0, 0.0);
        let p1 = Point::new(0.0, 0.0, 0.0);
        let p2 = Point::new(1.0, 0.0, 0.0);

        let quality = analyze_triangle(p0, p1, p2);

        assert!(
            quality.area < 1e-10,
            "Degenerate triangle should have zero area"
        );
        assert!(
            quality.aspect_ratio.is_infinite(),
            "Triangle with coincident points should have infinite aspect ratio, got {}",
            quality.aspect_ratio
        );
    }
}
