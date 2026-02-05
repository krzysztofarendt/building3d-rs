//! dotbim (.bim) file format I/O.
//!
//! Implements reading and writing of the dotbim format, a minimalist
//! BIM file format using JSON with triangulated meshes.
//!
//! See: https://dotbim.net/

use anyhow::{Context, Result, anyhow};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;

use crate::geom::triangles::TriangleIndex;
use crate::{Building, HasMesh, Mesh, Point, Solid, Zone};

/// Schema version for dotbim format.
const SCHEMA_VERSION: &str = "1.1.0";

/// Root structure of a dotbim file.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BimFile {
    /// Schema version (should be "1.1.0")
    pub schema_version: String,
    /// Meshes containing geometry
    pub meshes: Vec<BimMesh>,
    /// Elements referencing meshes with metadata
    pub elements: Vec<BimElement>,
    /// File-level metadata
    #[serde(default)]
    pub info: HashMap<String, String>,
}

/// A mesh in dotbim format.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BimMesh {
    /// Unique mesh identifier
    pub mesh_id: usize,
    /// Flat array of coordinates: [x0, y0, z0, x1, y1, z1, ...]
    pub coordinates: Vec<f64>,
    /// Flat array of triangle indices: [f0_v0, f0_v1, f0_v2, f1_v0, ...]
    pub indices: Vec<usize>,
}

/// An element in dotbim format.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BimElement {
    /// Reference to mesh_id
    pub mesh_id: usize,
    /// Element type (e.g., "Wall", "Floor", "Solid")
    #[serde(rename = "type")]
    pub element_type: String,
    /// RGBA color [r, g, b, a] where values are 0-255
    pub color: BimColor,
    /// 4x4 transformation matrix (column-major, 16 values)
    /// Identity matrix if no transformation
    #[serde(default = "default_matrix")]
    pub vector: BimVector,
    /// 3D rotation quaternion [x, y, z, w]
    #[serde(default = "default_rotation")]
    pub rotation: BimRotation,
    /// Unique identifier for the element
    pub guid: String,
    /// Element metadata
    #[serde(default)]
    pub info: HashMap<String, String>,
    /// Optional face colors (per-triangle RGBA, flattened)
    #[serde(skip_serializing_if = "Option::is_none")]
    pub face_colors: Option<Vec<u8>>,
}

/// RGBA color as [r, g, b, a] where values are 0-255.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BimColor {
    pub r: u8,
    pub g: u8,
    pub b: u8,
    pub a: u8,
}

/// Translation vector [x, y, z].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BimVector {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

/// Rotation quaternion [x, y, z, w].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BimRotation {
    pub qx: f64,
    pub qy: f64,
    pub qz: f64,
    pub qw: f64,
}

fn default_matrix() -> BimVector {
    BimVector {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    }
}

fn default_rotation() -> BimRotation {
    BimRotation {
        qx: 0.0,
        qy: 0.0,
        qz: 0.0,
        qw: 1.0,
    }
}

impl Default for BimColor {
    fn default() -> Self {
        Self {
            r: 128,
            g: 128,
            b: 128,
            a: 255,
        }
    }
}

/// Writes a Building to a dotbim (.bim) file.
///
/// # Arguments
/// * `path` - Output file path
/// * `building` - Building to export
///
/// # Example
/// ```ignore
/// use building3d::io::bim::write_bim;
/// use std::path::Path;
///
/// write_bim(Path::new("model.bim"), &building).unwrap();
/// ```
pub fn write_bim(path: &Path, building: &Building) -> Result<()> {
    let bim_file = building_to_bim(building)?;

    let file =
        File::create(path).with_context(|| format!("Failed to create file: {}", path.display()))?;
    let writer = BufWriter::new(file);

    serde_json::to_writer_pretty(writer, &bim_file)
        .with_context(|| format!("Failed to write BIM file: {}", path.display()))?;

    Ok(())
}

/// Reads a dotbim (.bim) file.
///
/// Note: Returns the raw BimFile structure. Converting to Building
/// requires additional processing since dotbim doesn't preserve
/// the full Building hierarchy.
///
/// # Arguments
/// * `path` - Input file path
pub fn read_bim(path: &Path) -> Result<BimFile> {
    let file =
        File::open(path).with_context(|| format!("Failed to open file: {}", path.display()))?;
    let reader = BufReader::new(file);

    let bim_file: BimFile = serde_json::from_reader(reader)
        .with_context(|| format!("Failed to parse BIM file: {}", path.display()))?;

    Ok(bim_file)
}

/// Converts a Building to dotbim format.
pub fn building_to_bim(building: &Building) -> Result<BimFile> {
    let mut meshes = Vec::new();
    let mut elements = Vec::new();
    let mut mesh_id = 0;

    // Export each zone
    for zone in building.zones() {
        for solid in zone.solids() {
            let (bim_mesh, bim_element) = solid_to_bim_element(solid, mesh_id)?;
            meshes.push(bim_mesh);
            elements.push(bim_element);
            mesh_id += 1;
        }
    }

    let mut info = HashMap::new();
    info.insert("name".to_string(), building.name.clone());
    info.insert("exported_by".to_string(), "building3d-rs".to_string());

    Ok(BimFile {
        schema_version: SCHEMA_VERSION.to_string(),
        meshes,
        elements,
        info,
    })
}

/// Converts a Solid to a BimMesh and BimElement.
fn solid_to_bim_element(solid: &Solid, mesh_id: usize) -> Result<(BimMesh, BimElement)> {
    let mesh = solid.copy_mesh();

    let bim_mesh = mesh_to_bim_mesh(&mesh, mesh_id)?;

    let mut info = HashMap::new();
    info.insert("name".to_string(), solid.name.clone());

    let bim_element = BimElement {
        mesh_id,
        element_type: "Solid".to_string(),
        color: BimColor::default(),
        vector: default_matrix(),
        rotation: default_rotation(),
        guid: solid.uid.as_str().to_string(),
        info,
        face_colors: None,
    };

    Ok((bim_mesh, bim_element))
}

/// Converts a Mesh to BimMesh format.
fn mesh_to_bim_mesh(mesh: &Mesh, mesh_id: usize) -> Result<BimMesh> {
    let vertices = mesh.vertices();

    // Flatten coordinates
    let mut coordinates = Vec::with_capacity(vertices.len() * 3);
    for v in vertices {
        coordinates.push(v.x);
        coordinates.push(v.y);
        coordinates.push(v.z);
    }

    // Flatten indices
    let faces = mesh.faces().ok_or_else(|| anyhow!("Mesh has no faces"))?;
    let mut indices = Vec::with_capacity(faces.len() * 3);
    for face in faces {
        indices.push(face.0);
        indices.push(face.1);
        indices.push(face.2);
    }

    Ok(BimMesh {
        mesh_id,
        coordinates,
        indices,
    })
}

/// Converts a BimFile to a collection of Solids.
///
/// Note: dotbim doesn't preserve the full Building hierarchy,
/// so we create a flat structure with one Zone containing all Solids.
pub fn bim_to_solids(bim_file: &BimFile) -> Result<Vec<Solid>> {
    let mut solids = Vec::new();

    for element in &bim_file.elements {
        // Find the corresponding mesh
        let mesh = bim_file
            .meshes
            .iter()
            .find(|m| m.mesh_id == element.mesh_id)
            .ok_or_else(|| anyhow!("Mesh {} not found for element", element.mesh_id))?;

        // Convert to Solid
        let solid = bim_element_to_solid(mesh, element)?;
        solids.push(solid);
    }

    Ok(solids)
}

/// Converts a BimMesh and BimElement to a Solid.
fn bim_element_to_solid(mesh: &BimMesh, element: &BimElement) -> Result<Solid> {
    // Parse coordinates into Points
    if !mesh.coordinates.len().is_multiple_of(3) {
        return Err(anyhow!(
            "Invalid coordinate count: {}",
            mesh.coordinates.len()
        ));
    }

    let mut vertices = Vec::with_capacity(mesh.coordinates.len() / 3);
    for chunk in mesh.coordinates.chunks(3) {
        vertices.push(Point::new(chunk[0], chunk[1], chunk[2]));
    }

    // Parse indices into TriangleIndex
    if !mesh.indices.len().is_multiple_of(3) {
        return Err(anyhow!("Invalid index count: {}", mesh.indices.len()));
    }

    let mut faces = Vec::with_capacity(mesh.indices.len() / 3);
    for chunk in mesh.indices.chunks(3) {
        faces.push(TriangleIndex(chunk[0], chunk[1], chunk[2]));
    }

    // Create mesh
    let bim_mesh = Mesh::new(vertices, Some(faces));

    // Get name from element info or use GUID
    let name = element
        .info
        .get("name")
        .cloned()
        .unwrap_or_else(|| format!("element_{}", element.mesh_id));

    // Create a simple box solid from the mesh
    // Note: This is a simplified conversion - the original wall structure is lost
    crate::io::stl::stl_to_solid(&bim_mesh, &name)
}

/// Converts a BimFile to a Building.
///
/// Creates a Building with a single Zone containing all elements as Solids.
pub fn bim_to_building(bim_file: &BimFile, name: &str) -> Result<Building> {
    let solids = bim_to_solids(bim_file)?;

    let zone = Zone::new("imported", solids)?;
    let building = Building::new(name, vec![zone])?;

    Ok(building)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Solid;
    use std::fs;
    use tempfile::tempdir;

    fn make_test_building() -> Result<Building> {
        let solid = Solid::from_box(2.0, 3.0, 4.0, Some((0.0, 0.0, 0.0)), "box")?;
        let zone = Zone::new("zone1", vec![solid])?;
        let building = Building::new("test_building", vec![zone])?;
        Ok(building)
    }

    #[test]
    fn test_write_bim() -> Result<()> {
        let dir = tempdir()?;
        let path = dir.path().join("test.bim");

        let building = make_test_building()?;
        write_bim(&path, &building)?;

        assert!(path.exists());

        // Read back and verify structure
        let content = fs::read_to_string(&path)?;
        assert!(content.contains("schema_version"));
        assert!(content.contains("1.1.0"));
        assert!(content.contains("meshes"));
        assert!(content.contains("elements"));

        Ok(())
    }

    #[test]
    fn test_read_bim() -> Result<()> {
        let dir = tempdir()?;
        let path = dir.path().join("test.bim");

        let building = make_test_building()?;
        write_bim(&path, &building)?;

        let bim_file = read_bim(&path)?;

        assert_eq!(bim_file.schema_version, "1.1.0");
        assert_eq!(bim_file.meshes.len(), 1);
        assert_eq!(bim_file.elements.len(), 1);

        Ok(())
    }

    #[test]
    fn test_bim_roundtrip() -> Result<()> {
        let dir = tempdir()?;
        let path = dir.path().join("test.bim");

        let building = make_test_building()?;
        let original_solid_count = building.zones().iter().flat_map(|z| z.solids()).count();

        write_bim(&path, &building)?;
        let bim_file = read_bim(&path)?;

        // Convert back to Building
        let imported = bim_to_building(&bim_file, "imported")?;
        let imported_solid_count = imported.zones().iter().flat_map(|z| z.solids()).count();

        assert_eq!(original_solid_count, imported_solid_count);

        Ok(())
    }

    #[test]
    fn test_mesh_to_bim_mesh() -> Result<()> {
        let vertices = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(0.5, 1.0, 0.0),
        ];
        let faces = vec![TriangleIndex(0, 1, 2)];
        let mesh = Mesh::new(vertices, Some(faces));

        let bim_mesh = mesh_to_bim_mesh(&mesh, 0)?;

        assert_eq!(bim_mesh.mesh_id, 0);
        assert_eq!(bim_mesh.coordinates.len(), 9); // 3 vertices * 3 coords
        assert_eq!(bim_mesh.indices.len(), 3); // 1 triangle * 3 indices

        // Verify coordinates
        assert_eq!(bim_mesh.coordinates[0], 0.0);
        assert_eq!(bim_mesh.coordinates[3], 1.0);
        assert_eq!(bim_mesh.coordinates[7], 1.0);

        // Verify indices
        assert_eq!(bim_mesh.indices, vec![0, 1, 2]);

        Ok(())
    }

    #[test]
    fn test_bim_element_types() -> Result<()> {
        let dir = tempdir()?;
        let path = dir.path().join("test.bim");

        let building = make_test_building()?;
        write_bim(&path, &building)?;

        let bim_file = read_bim(&path)?;

        // Check element has correct type
        assert_eq!(bim_file.elements[0].element_type, "Solid");

        // Check element has GUID
        assert!(!bim_file.elements[0].guid.is_empty());

        Ok(())
    }

    #[test]
    fn test_bim_info_metadata() -> Result<()> {
        let dir = tempdir()?;
        let path = dir.path().join("test.bim");

        let building = make_test_building()?;
        write_bim(&path, &building)?;

        let bim_file = read_bim(&path)?;

        // Check file-level metadata
        assert!(bim_file.info.contains_key("name"));
        assert_eq!(
            bim_file.info.get("name"),
            Some(&"test_building".to_string())
        );
        assert!(bim_file.info.contains_key("exported_by"));

        Ok(())
    }
}
