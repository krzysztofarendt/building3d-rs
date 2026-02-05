//! B3D file format I/O.
//!
//! B3D is the native JSON format for building3d models. It preserves the full
//! hierarchy structure including UIDs for reference integrity.

use crate::Building;
use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufReader, BufWriter};
use std::path::Path;

/// Writes a building to a B3D (JSON) file.
///
/// # Arguments
/// * `path` - Path to the output file
/// * `building` - The building to serialize
///
/// # Example
/// ```no_run
/// use building3d::{Building, Solid};
/// use building3d::io::write_b3d;
/// use std::path::Path;
///
/// let solid = Solid::from_box(1.0, 1.0, 1.0, None, "box");
/// let building = Building::from_solids("my_building", vec![solid]);
/// write_b3d(Path::new("model.b3d"), &building).unwrap();
/// ```
pub fn write_b3d(path: &Path, building: &Building) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Failed to create file: {}", path.display()))?;
    let writer = BufWriter::new(file);

    serde_json::to_writer_pretty(writer, building)
        .with_context(|| format!("Failed to serialize building to: {}", path.display()))?;

    Ok(())
}

/// Reads a building from a B3D (JSON) file.
///
/// # Arguments
/// * `path` - Path to the input file
///
/// # Returns
/// The deserialized Building
///
/// # Example
/// ```no_run
/// use building3d::io::read_b3d;
/// use std::path::Path;
///
/// let building = read_b3d(Path::new("model.b3d")).unwrap();
/// println!("Loaded building: {}", building.name);
/// ```
pub fn read_b3d(path: &Path) -> Result<Building> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open file: {}", path.display()))?;
    let reader = BufReader::new(file);

    let building: Building = serde_json::from_reader(reader)
        .with_context(|| format!("Failed to deserialize building from: {}", path.display()))?;

    Ok(building)
}

/// Serializes a building to a B3D JSON string.
///
/// Useful for in-memory operations or network transfer.
pub fn to_b3d_string(building: &Building) -> Result<String> {
    serde_json::to_string_pretty(building).context("Failed to serialize building to string")
}

/// Deserializes a building from a B3D JSON string.
///
/// Useful for in-memory operations or network transfer.
pub fn from_b3d_string(json: &str) -> Result<Building> {
    serde_json::from_str(json).context("Failed to deserialize building from string")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Solid, Zone};
    use tempfile::tempdir;

    #[test]
    fn test_write_and_read_b3d() -> Result<()> {
        let dir = tempdir()?;
        let path = dir.path().join("test.b3d");

        // Create a building
        let s1 = Solid::from_box(1.0, 2.0, 3.0, None, "box1");
        let s2 = Solid::from_box(2.0, 2.0, 2.0, Some((2.0, 0.0, 0.0)), "box2");
        let original = Building::from_solids("test_building", vec![s1, s2]);

        // Write to file
        write_b3d(&path, &original)?;

        // Read back
        let loaded = read_b3d(&path)?;

        // Verify
        assert_eq!(loaded.name, original.name);
        assert_eq!(loaded.zones().len(), original.zones().len());
        assert_eq!(loaded.solids().len(), original.solids().len());

        // Check volume is preserved
        assert!((loaded.volume() - original.volume()).abs() < 1e-10);

        Ok(())
    }

    #[test]
    fn test_b3d_roundtrip_with_zones() -> Result<()> {
        let dir = tempdir()?;
        let path = dir.path().join("zones.b3d");

        // Create building with explicit zones
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1");
        let z1 = Zone::new("zone1", vec![s1]);

        let s2 = Solid::from_box(2.0, 2.0, 2.0, None, "box2");
        let z2 = Zone::new("zone2", vec![s2]);

        let original = Building::new("building_with_zones", vec![z1, z2]);

        // Roundtrip
        write_b3d(&path, &original)?;
        let loaded = read_b3d(&path)?;

        // Verify zones
        assert_eq!(loaded.zones().len(), 2);
        assert!(loaded.get_zone("zone1").is_some());
        assert!(loaded.get_zone("zone2").is_some());

        Ok(())
    }

    #[test]
    fn test_b3d_string_roundtrip() -> Result<()> {
        let solid = Solid::from_box(1.0, 1.0, 1.0, None, "cube");
        let original = Building::from_solids("test", vec![solid]);

        // Serialize to string
        let json = to_b3d_string(&original)?;

        // Verify it's valid JSON
        assert!(json.contains("\"name\":"));
        assert!(json.contains("\"test\""));

        // Deserialize back
        let loaded = from_b3d_string(&json)?;

        assert_eq!(loaded.name, original.name);
        assert!((loaded.volume() - original.volume()).abs() < 1e-10);

        Ok(())
    }

    #[test]
    fn test_b3d_preserves_uids() -> Result<()> {
        let dir = tempdir()?;
        let path = dir.path().join("uids.b3d");

        let solid = Solid::from_box(1.0, 1.0, 1.0, None, "box");
        let original = Building::from_solids("test", vec![solid]);
        let original_uid = original.uid.as_str().to_string();

        write_b3d(&path, &original)?;
        let loaded = read_b3d(&path)?;

        // UIDs should be preserved
        assert_eq!(loaded.uid.as_str(), original_uid);

        Ok(())
    }

    #[test]
    fn test_b3d_preserves_geometry() -> Result<()> {
        let dir = tempdir()?;
        let path = dir.path().join("geom.b3d");

        let solid = Solid::from_box(3.0, 4.0, 5.0, Some((1.0, 2.0, 3.0)), "box");
        let original = Building::from_solids("test", vec![solid]);

        write_b3d(&path, &original)?;
        let loaded = read_b3d(&path)?;

        // Volume should be exactly preserved
        let orig_vol = original.volume();
        let load_vol = loaded.volume();
        assert!(
            (orig_vol - load_vol).abs() < 1e-10,
            "Volume mismatch: {} vs {}",
            orig_vol,
            load_vol
        );

        // Bounding box should be preserved
        let (orig_min, orig_max) = original.bbox();
        let (load_min, load_max) = loaded.bbox();
        assert!(orig_min.is_close(&load_min));
        assert!(orig_max.is_close(&load_max));

        Ok(())
    }

    #[test]
    fn test_read_nonexistent_file() {
        let result = read_b3d(Path::new("/nonexistent/path/file.b3d"));
        assert!(result.is_err());
    }
}
