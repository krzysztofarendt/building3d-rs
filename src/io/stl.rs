//! STL file format I/O.
//!
//! STL (STereoLithography) is a common 3D mesh format that stores triangulated
//! surfaces. Note that STL files lose the building hierarchy information -
//! they only contain raw triangles with normals.

use crate::{HasMesh, Mesh, Point, Solid, Vector};
use anyhow::{Context, Result, anyhow};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// STL file format variant.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StlFormat {
    /// ASCII text format (human-readable, larger file size)
    Ascii,
    /// Binary format (compact, faster to read/write)
    Binary,
}

/// Writes a mesh to an STL file.
///
/// # Arguments
/// * `path` - Path to the output file
/// * `mesh` - The mesh to export (must have triangulated faces)
/// * `name` - Name to use in the STL header
/// * `format` - ASCII or Binary format
///
/// # Example
/// ```no_run
/// use building3d::{Solid, HasMesh};
/// use building3d::io::stl::{write_stl, StlFormat};
/// use std::path::Path;
///
/// let solid = Solid::from_box(1.0, 1.0, 1.0, None, "box").unwrap();
/// let mesh = solid.copy_mesh();
/// write_stl(Path::new("model.stl"), &mesh, "box", StlFormat::Binary).unwrap();
/// ```
pub fn write_stl(path: &Path, mesh: &Mesh, name: &str, format: StlFormat) -> Result<()> {
    match format {
        StlFormat::Ascii => write_stl_ascii(path, mesh, name),
        StlFormat::Binary => write_stl_binary(path, mesh, name),
    }
}

/// Writes a mesh to an ASCII STL file.
fn write_stl_ascii(path: &Path, mesh: &Mesh, name: &str) -> Result<()> {
    let file =
        File::create(path).with_context(|| format!("Failed to create file: {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    let faces = mesh
        .faces
        .as_ref()
        .ok_or_else(|| anyhow!("Mesh has no triangulated faces"))?;

    writeln!(writer, "solid {}", name)?;

    for tri in faces {
        let p0 = mesh.vertices[tri.0];
        let p1 = mesh.vertices[tri.1];
        let p2 = mesh.vertices[tri.2];

        // Calculate face normal
        let v1 = p1 - p0;
        let v2 = p2 - p0;
        let cross = v1.cross(&v2);
        let normal = cross.normalize().unwrap_or(Vector::new(0.0, 0.0, 1.0));

        writeln!(
            writer,
            "  facet normal {} {} {}",
            normal.dx, normal.dy, normal.dz
        )?;
        writeln!(writer, "    outer loop")?;
        writeln!(writer, "      vertex {} {} {}", p0.x, p0.y, p0.z)?;
        writeln!(writer, "      vertex {} {} {}", p1.x, p1.y, p1.z)?;
        writeln!(writer, "      vertex {} {} {}", p2.x, p2.y, p2.z)?;
        writeln!(writer, "    endloop")?;
        writeln!(writer, "  endfacet")?;
    }

    writeln!(writer, "endsolid {}", name)?;

    Ok(())
}

/// Writes a mesh to a binary STL file.
fn write_stl_binary(path: &Path, mesh: &Mesh, name: &str) -> Result<()> {
    let file =
        File::create(path).with_context(|| format!("Failed to create file: {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    let faces = mesh
        .faces
        .as_ref()
        .ok_or_else(|| anyhow!("Mesh has no triangulated faces"))?;

    // 80-byte header
    let mut header = [0u8; 80];
    let header_str = format!("binary STL - {}", name);
    let bytes = header_str.as_bytes();
    let len = bytes.len().min(80);
    header[..len].copy_from_slice(&bytes[..len]);
    writer.write_all(&header)?;

    // Number of triangles (u32 little-endian)
    let num_triangles = faces.len() as u32;
    writer.write_all(&num_triangles.to_le_bytes())?;

    // Write each triangle
    for tri in faces {
        let p0 = mesh.vertices[tri.0];
        let p1 = mesh.vertices[tri.1];
        let p2 = mesh.vertices[tri.2];

        // Calculate face normal
        let v1 = p1 - p0;
        let v2 = p2 - p0;
        let cross = v1.cross(&v2);
        let normal = cross.normalize().unwrap_or(Vector::new(0.0, 0.0, 1.0));

        // Normal (3 x f32)
        writer.write_all(&(normal.dx as f32).to_le_bytes())?;
        writer.write_all(&(normal.dy as f32).to_le_bytes())?;
        writer.write_all(&(normal.dz as f32).to_le_bytes())?;

        // Vertex 1 (3 x f32)
        writer.write_all(&(p0.x as f32).to_le_bytes())?;
        writer.write_all(&(p0.y as f32).to_le_bytes())?;
        writer.write_all(&(p0.z as f32).to_le_bytes())?;

        // Vertex 2 (3 x f32)
        writer.write_all(&(p1.x as f32).to_le_bytes())?;
        writer.write_all(&(p1.y as f32).to_le_bytes())?;
        writer.write_all(&(p1.z as f32).to_le_bytes())?;

        // Vertex 3 (3 x f32)
        writer.write_all(&(p2.x as f32).to_le_bytes())?;
        writer.write_all(&(p2.y as f32).to_le_bytes())?;
        writer.write_all(&(p2.z as f32).to_le_bytes())?;

        // Attribute byte count (unused, set to 0)
        writer.write_all(&0u16.to_le_bytes())?;
    }

    Ok(())
}

/// Reads triangles from an STL file into a Mesh.
///
/// Note: STL files do not contain hierarchy information, so the result
/// is a flat mesh of triangles.
///
/// # Arguments
/// * `path` - Path to the input STL file
///
/// # Returns
/// A Mesh containing the vertices and triangulated faces.
pub fn read_stl(path: &Path) -> Result<Mesh> {
    let file =
        File::open(path).with_context(|| format!("Failed to open file: {}", path.display()))?;
    let mut reader = BufReader::new(file);

    // Check if binary or ASCII by looking at the first few bytes
    let mut header = [0u8; 80];
    std::io::Read::read_exact(&mut reader, &mut header)?;

    // Reset to beginning
    drop(reader);
    let file =
        File::open(path).with_context(|| format!("Failed to open file: {}", path.display()))?;
    let _reader = BufReader::new(file);

    // Check if it starts with "solid" (ASCII) or not (binary)
    // Note: Some binary files can also start with "solid" in header,
    // so we also check file size
    let header_str = String::from_utf8_lossy(&header);
    let metadata = std::fs::metadata(path)?;
    let file_size = metadata.len();

    // Binary STL: 80 header + 4 bytes count + 50 bytes per triangle
    // If "solid" appears and file size doesn't match binary formula, it's ASCII
    if header_str.trim_start().starts_with("solid") {
        // Try to determine if it's really ASCII by checking for "facet"
        let file =
            File::open(path).with_context(|| format!("Failed to open file: {}", path.display()))?;
        let mut test_reader = BufReader::new(file);
        let mut first_lines = String::new();
        for _ in 0..5 {
            let mut line = String::new();
            if test_reader.read_line(&mut line).is_ok() {
                first_lines.push_str(&line);
            }
        }

        if first_lines.contains("facet") || first_lines.contains("vertex") {
            return read_stl_ascii(path);
        }
    }

    // Default to binary
    read_stl_binary(path, file_size)
}

/// Reads an ASCII STL file.
fn read_stl_ascii(path: &Path) -> Result<Mesh> {
    let file =
        File::open(path).with_context(|| format!("Failed to open file: {}", path.display()))?;
    let reader = BufReader::new(file);

    let mut vertices: Vec<Point> = Vec::new();
    let mut faces: Vec<crate::TriangleIndex> = Vec::new();
    let mut vertex_map: HashMap<(i64, i64, i64), usize> = HashMap::new();

    let mut current_vertices: Vec<Point> = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let trimmed = line.trim();

        if trimmed.starts_with("vertex") {
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.len() >= 4 {
                let x: f64 = parts[1].parse().context("Invalid vertex x")?;
                let y: f64 = parts[2].parse().context("Invalid vertex y")?;
                let z: f64 = parts[3].parse().context("Invalid vertex z")?;
                current_vertices.push(Point::new(x, y, z));
            }
        } else if trimmed.starts_with("endloop") {
            if current_vertices.len() == 3 {
                let i0 = add_dedup_vertex(&mut vertex_map, &mut vertices, current_vertices[0]);
                let i1 = add_dedup_vertex(&mut vertex_map, &mut vertices, current_vertices[1]);
                let i2 = add_dedup_vertex(&mut vertex_map, &mut vertices, current_vertices[2]);
                faces.push(crate::TriangleIndex(i0, i1, i2));
            }
            current_vertices.clear();
        }
    }

    Ok(Mesh {
        vertices,
        faces: Some(faces),
    })
}

/// Reads a binary STL file.
fn read_stl_binary(path: &Path, _file_size: u64) -> Result<Mesh> {
    let file =
        File::open(path).with_context(|| format!("Failed to open file: {}", path.display()))?;
    let mut reader = BufReader::new(file);

    // Skip 80-byte header
    let mut header = [0u8; 80];
    std::io::Read::read_exact(&mut reader, &mut header)?;

    // Read triangle count
    let mut count_bytes = [0u8; 4];
    std::io::Read::read_exact(&mut reader, &mut count_bytes)?;
    let num_triangles = u32::from_le_bytes(count_bytes) as usize;

    let mut vertices: Vec<Point> = Vec::with_capacity(num_triangles * 3);
    let mut faces: Vec<crate::TriangleIndex> = Vec::with_capacity(num_triangles);
    let mut vertex_map: HashMap<(i64, i64, i64), usize> = HashMap::new();

    for _ in 0..num_triangles {
        // Skip normal (3 x f32 = 12 bytes)
        let mut normal_bytes = [0u8; 12];
        std::io::Read::read_exact(&mut reader, &mut normal_bytes)?;

        // Read 3 vertices
        let mut tri_idx: [usize; 3] = [0; 3];
        for slot in tri_idx.iter_mut() {
            let mut v_bytes = [0u8; 12];
            std::io::Read::read_exact(&mut reader, &mut v_bytes)?;

            let x = f32::from_le_bytes([v_bytes[0], v_bytes[1], v_bytes[2], v_bytes[3]]) as f64;
            let y = f32::from_le_bytes([v_bytes[4], v_bytes[5], v_bytes[6], v_bytes[7]]) as f64;
            let z = f32::from_le_bytes([v_bytes[8], v_bytes[9], v_bytes[10], v_bytes[11]]) as f64;

            *slot = add_dedup_vertex(&mut vertex_map, &mut vertices, Point::new(x, y, z));
        }

        // Skip attribute byte count (2 bytes)
        let mut attr_bytes = [0u8; 2];
        std::io::Read::read_exact(&mut reader, &mut attr_bytes)?;

        faces.push(crate::TriangleIndex(tri_idx[0], tri_idx[1], tri_idx[2]));
    }

    Ok(Mesh {
        vertices,
        faces: Some(faces),
    })
}

const STL_DEDUP_SCALE: f64 = 1e9;

fn stl_vertex_key(p: Point) -> (i64, i64, i64) {
    (
        (p.x * STL_DEDUP_SCALE).round() as i64,
        (p.y * STL_DEDUP_SCALE).round() as i64,
        (p.z * STL_DEDUP_SCALE).round() as i64,
    )
}

fn add_dedup_vertex(
    map: &mut HashMap<(i64, i64, i64), usize>,
    vertices: &mut Vec<Point>,
    p: Point,
) -> usize {
    let key = stl_vertex_key(p);
    if let Some(&idx) = map.get(&key) {
        return idx;
    }
    let idx = vertices.len();
    vertices.push(p);
    map.insert(key, idx);
    idx
}

/// Convenience function to write any HasMesh object to STL.
pub fn write_stl_from_mesh<T: HasMesh>(
    path: &Path,
    obj: &T,
    name: &str,
    format: StlFormat,
) -> Result<()> {
    let mesh = obj.copy_mesh();
    write_stl(path, &mesh, name, format)
}

/// Creates a Solid from an STL mesh.
///
/// Note: This creates a simple solid with a single wall containing all triangles
/// as individual polygons. The original structure is lost since STL doesn't
/// preserve hierarchy.
pub fn stl_to_solid(mesh: &Mesh, name: &str) -> Result<Solid> {
    // For now, just create a basic solid from the mesh vertices
    // This is a simplified implementation - a full implementation would
    // need to reconstruct walls from connected triangles
    use crate::{Polygon, Wall};

    let faces = mesh
        .faces
        .as_ref()
        .ok_or_else(|| anyhow!("Mesh has no faces"))?;

    let mut polygons: Vec<Polygon> = Vec::new();

    for (i, tri) in faces.iter().enumerate() {
        let p0 = mesh.vertices[tri.0];
        let p1 = mesh.vertices[tri.1];
        let p2 = mesh.vertices[tri.2];

        let poly_name = format!("tri_{}", i);
        let poly = Polygon::new(&poly_name, vec![p0, p1, p2], None)?;
        polygons.push(poly);
    }

    // Group all polygons into a single wall
    let wall = Wall::new("mesh", polygons)?;
    let solid = Solid::new(name, vec![wall])?;

    Ok(solid)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    fn create_test_mesh() -> Mesh {
        // Simple tetrahedron
        let vertices = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(0.5, 1.0, 0.0),
            Point::new(0.5, 0.5, 1.0),
        ];
        let faces = vec![
            crate::TriangleIndex(0, 1, 2), // Base
            crate::TriangleIndex(0, 1, 3), // Side 1
            crate::TriangleIndex(1, 2, 3), // Side 2
            crate::TriangleIndex(2, 0, 3), // Side 3
        ];
        Mesh {
            vertices,
            faces: Some(faces),
        }
    }

    #[test]
    fn test_write_stl_ascii() -> Result<()> {
        let dir = tempdir()?;
        let path = dir.path().join("test.stl");

        let mesh = create_test_mesh();
        write_stl(&path, &mesh, "test", StlFormat::Ascii)?;

        // Verify file exists and has content
        let content = std::fs::read_to_string(&path)?;
        assert!(content.contains("solid test"));
        assert!(content.contains("facet normal"));
        assert!(content.contains("vertex"));
        assert!(content.contains("endsolid test"));

        Ok(())
    }

    #[test]
    fn test_write_stl_binary() -> Result<()> {
        let dir = tempdir()?;
        let path = dir.path().join("test.stl");

        let mesh = create_test_mesh();
        write_stl(&path, &mesh, "test", StlFormat::Binary)?;

        // Verify file exists and has correct size
        // 80 header + 4 count + 4 triangles * 50 bytes = 284 bytes
        let metadata = std::fs::metadata(&path)?;
        assert_eq!(metadata.len(), 80 + 4 + 4 * 50);

        Ok(())
    }

    #[test]
    fn test_read_stl_ascii() -> Result<()> {
        let dir = tempdir()?;
        let path = dir.path().join("test.stl");

        // Write ASCII STL
        let original = create_test_mesh();
        write_stl(&path, &original, "test", StlFormat::Ascii)?;

        // Read it back
        let loaded = read_stl(&path)?;

        // Should have same number of triangles
        assert_eq!(
            loaded.faces.as_ref().unwrap().len(),
            original.faces.as_ref().unwrap().len()
        );

        // Vertices are deduplicated
        assert_eq!(loaded.vertices.len(), original.vertices.len());

        Ok(())
    }

    #[test]
    fn test_read_stl_binary() -> Result<()> {
        let dir = tempdir()?;
        let path = dir.path().join("test.stl");

        // Write binary STL
        let original = create_test_mesh();
        write_stl(&path, &original, "test", StlFormat::Binary)?;

        // Read it back
        let loaded = read_stl(&path)?;

        // Should have same number of triangles
        assert_eq!(
            loaded.faces.as_ref().unwrap().len(),
            original.faces.as_ref().unwrap().len()
        );

        // Vertices are deduplicated
        assert_eq!(loaded.vertices.len(), original.vertices.len());

        Ok(())
    }

    #[test]
    fn test_stl_roundtrip_box() -> Result<()> {
        let dir = tempdir()?;
        let path = dir.path().join("box.stl");

        // Create a box solid
        let solid = Solid::from_box(2.0, 3.0, 4.0, None, "box")?;
        let original_mesh = solid.copy_mesh();

        // Write to STL
        write_stl(&path, &original_mesh, "box", StlFormat::Binary)?;

        // Read back
        let loaded_mesh = read_stl(&path)?;

        // Should have same number of triangles
        // A box has 6 faces, each face is 2 triangles = 12 triangles
        assert_eq!(loaded_mesh.faces.as_ref().unwrap().len(), 12);

        Ok(())
    }

    #[test]
    fn test_write_stl_from_solid() -> Result<()> {
        let dir = tempdir()?;
        let path = dir.path().join("solid.stl");

        let solid = Solid::from_box(1.0, 1.0, 1.0, None, "cube")?;
        write_stl_from_mesh(&path, &solid, "cube", StlFormat::Ascii)?;

        let content = std::fs::read_to_string(&path)?;
        assert!(content.contains("solid cube"));

        Ok(())
    }

    #[test]
    fn test_stl_to_solid() -> Result<()> {
        let mesh = create_test_mesh();
        let solid = stl_to_solid(&mesh, "tetrahedron")?;

        assert_eq!(solid.name, "tetrahedron");
        // Should have 4 polygons (one per triangle)
        assert_eq!(solid.polygons().len(), 4);

        Ok(())
    }
}
