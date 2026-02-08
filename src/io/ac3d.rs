//! AC3D file format reader.
//!
//! Reads `.ac` files (AC3D format) and converts them into a building3d [`Building`].
//! Each unique material becomes a [`Wall`], each surface becomes a [`Polygon`],
//! and all walls go into a single [`Solid`].

use crate::{Building, Point, Polygon, Solid, Vector, Wall};
use anyhow::{Context, Result, anyhow};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Coordinate system convention for the AC3D file.
pub enum Ac3dCoordSystem {
    /// Y-up (AC3D/RAVEN native) — no transform.
    YUp,
    /// Z-up (building3d convention) — applies (x, y, z) → (x, −z, y).
    ZUp,
}

struct Ac3dSurface {
    material_index: usize,
    vertex_indices: Vec<usize>,
}

struct Ac3dObject {
    vertices: Vec<Point>,
    surfaces: Vec<Ac3dSurface>,
}

/// Reads an AC3D file and returns (material_names, Building).
///
/// Each material becomes a [`Wall`]; each surface becomes a [`Polygon`].
/// All walls go into a single [`Solid`] named `solid_name`.
///
/// # Arguments
/// * `path` - Path to the `.ac` file
/// * `solid_name` - Name for the resulting solid
/// * `coord_system` - Coordinate system convention
pub fn read_ac3d(
    path: &Path,
    solid_name: &str,
    coord_system: Ac3dCoordSystem,
) -> Result<(Vec<String>, Building)> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open AC3D file: {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut lines = reader.lines();

    // Verify header
    let header = lines.next().ok_or_else(|| anyhow!("Empty AC3D file"))??;
    if !header.starts_with("AC3D") {
        return Err(anyhow!("Not an AC3D file: missing AC3D header"));
    }

    let mut material_names: Vec<String> = Vec::new();
    let mut objects: Vec<Ac3dObject> = Vec::new();
    let mut current_object: Option<Ac3dObject> = None;

    while let Some(line_result) = lines.next() {
        let line = line_result?;
        let line = line.trim();

        if line.starts_with("MATERIAL ") {
            if let Some(name) = extract_material_name(line) {
                material_names.push(name);
            }
        } else if line == "OBJECT poly" {
            current_object = Some(Ac3dObject {
                vertices: Vec::new(),
                surfaces: Vec::new(),
            });
        } else if let Some(rest) = line.strip_prefix("numvert ") {
            let n: usize = rest
                .trim()
                .parse()
                .with_context(|| format!("Invalid numvert: {line}"))?;
            if let Some(ref mut obj) = current_object {
                for _ in 0..n {
                    let vline = lines
                        .next()
                        .ok_or_else(|| anyhow!("Unexpected EOF reading vertices"))??;
                    let coords: Vec<f64> = vline
                        .split_whitespace()
                        .map(|s| s.parse::<f64>())
                        .collect::<std::result::Result<Vec<_>, _>>()
                        .with_context(|| format!("Invalid vertex line: {vline}"))?;
                    if coords.len() != 3 {
                        return Err(anyhow!("Expected 3 coordinates, got {}", coords.len()));
                    }
                    let pt = match coord_system {
                        Ac3dCoordSystem::YUp => Point::new(coords[0], coords[1], coords[2]),
                        Ac3dCoordSystem::ZUp => Point::new(coords[0], -coords[2], coords[1]),
                    };
                    obj.vertices.push(pt);
                }
            }
        } else if let Some(rest) = line.strip_prefix("numsurf ") {
            let n: usize = rest
                .trim()
                .parse()
                .with_context(|| format!("Invalid numsurf: {line}"))?;
            if let Some(ref mut obj) = current_object {
                for _ in 0..n {
                    let surf = parse_surface(&mut lines)?;
                    obj.surfaces.push(surf);
                }
            }
        } else if line == "kids 0"
            && let Some(obj) = current_object.take()
        {
            objects.push(obj);
        }
    }

    // Group all surfaces by material index across all objects
    let mut material_surfaces: HashMap<usize, Vec<(Vec<Point>, usize)>> = HashMap::new();
    for obj in &objects {
        for surf in &obj.surfaces {
            if surf.vertex_indices.len() < 3 {
                eprintln!(
                    "Warning: skipping surface with {} vertices (need >= 3)",
                    surf.vertex_indices.len()
                );
                continue;
            }
            let pts: Vec<Point> = surf
                .vertex_indices
                .iter()
                .map(|&i| obj.vertices[i])
                .collect();
            material_surfaces
                .entry(surf.material_index)
                .or_default()
                .push((pts, surf.material_index));
        }
    }

    // Build walls from material groups
    let mut walls: Vec<Wall> = Vec::new();
    let mut mat_indices: Vec<usize> = material_surfaces.keys().copied().collect();
    mat_indices.sort();

    for mat_idx in mat_indices {
        let surfaces = material_surfaces.remove(&mat_idx).unwrap();
        let mat_name = if mat_idx < material_names.len() {
            strip_material_prefix(&material_names[mat_idx])
        } else {
            format!("material_{mat_idx}")
        };

        let mut polygons: Vec<Polygon> = Vec::new();
        for (i, (pts, _)) in surfaces.iter().enumerate() {
            let poly_name = format!("{mat_name}_{i:03}");
            let normal = compute_polygon_normal(pts);
            match Polygon::new(&poly_name, pts.clone(), normal) {
                Ok(poly) => polygons.push(poly),
                Err(e) => {
                    eprintln!("Warning: skipping polygon {poly_name}: {e}");
                }
            }
        }

        if !polygons.is_empty() {
            walls.push(Wall::new(&mat_name, polygons)?);
        }
    }

    let solid = Solid::new(solid_name, walls)?;
    let building = Building::from_solids(solid_name, vec![solid])?;

    Ok((material_names, building))
}

/// Parse a single surface block (SURF, mat, refs lines).
fn parse_surface(lines: &mut std::io::Lines<BufReader<File>>) -> Result<Ac3dSurface> {
    let mut material_index: usize = 0;
    let mut vertex_indices: Vec<usize> = Vec::new();

    // Read lines until we get the refs block
    loop {
        let line = lines
            .next()
            .ok_or_else(|| anyhow!("Unexpected EOF in surface block"))??;
        let line = line.trim();

        if line.starts_with("SURF ") {
            // Surface type flag — ignored
        } else if let Some(rest) = line.strip_prefix("mat ") {
            material_index = rest
                .trim()
                .parse()
                .with_context(|| format!("Invalid mat index: {line}"))?;
        } else if let Some(rest) = line.strip_prefix("refs ") {
            let n: usize = rest
                .trim()
                .parse()
                .with_context(|| format!("Invalid refs count: {line}"))?;
            for _ in 0..n {
                let rline = lines
                    .next()
                    .ok_or_else(|| anyhow!("Unexpected EOF reading refs"))??;
                let idx: usize = rline
                    .split_whitespace()
                    .next()
                    .ok_or_else(|| anyhow!("Empty ref line"))?
                    .parse()
                    .with_context(|| format!("Invalid ref index: {rline}"))?;
                vertex_indices.push(idx);
            }
            break;
        }
    }

    Ok(Ac3dSurface {
        material_index,
        vertex_indices,
    })
}

/// Extract the quoted material name from a MATERIAL line.
fn extract_material_name(line: &str) -> Option<String> {
    let start = line.find('"')? + 1;
    let end = line[start..].find('"')? + start;
    Some(line[start..end].to_string())
}

/// Compute a normal vector for a polygon by trying successive vertex triples.
/// Returns `None` if no valid normal can be found (degenerate polygon).
fn compute_polygon_normal(pts: &[Point]) -> Option<Vector> {
    let n = pts.len();
    if n < 3 {
        return None;
    }
    for i in 0..n {
        let a = pts[i];
        let b = pts[(i + 1) % n];
        let c = pts[(i + 2) % n];
        let v1 = Vector::from_points(b, a);
        let v2 = Vector::from_points(b, c);
        let cross = v1.cross(&v2);
        if cross.length() > 1e-12 {
            return cross.normalize().ok();
        }
    }
    None
}

/// Strip common prefixes like "mat_scene09_" from material names.
fn strip_material_prefix(name: &str) -> String {
    // Pattern: "mat_sceneNN_actual_name"
    if let Some(rest) = name.strip_prefix("mat_")
        && let Some(pos) = rest.find('_')
    {
        return rest[pos + 1..].to_string();
    }
    name.to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::HasName;

    #[test]
    fn test_extract_material_name() {
        assert_eq!(
            extract_material_name(r#"MATERIAL "mat_scene09_concrete" rgb 0.1 0.1 0.1"#),
            Some("mat_scene09_concrete".to_string())
        );
    }

    #[test]
    fn test_strip_material_prefix() {
        assert_eq!(strip_material_prefix("mat_scene09_concrete"), "concrete");
        assert_eq!(strip_material_prefix("mat_scene09_windows"), "windows");
        assert_eq!(strip_material_prefix("plain_name"), "plain_name");
    }

    #[test]
    fn test_read_ac3d_scene9() {
        let path = Path::new("validation/bras/informed_sim/RavenModels/scene9.ac");
        if !path.exists() {
            eprintln!("Skipping test: scene9.ac not found");
            return;
        }

        let (materials, building) = read_ac3d(path, "seminar_room", Ac3dCoordSystem::ZUp).unwrap();

        // 5 materials
        assert_eq!(materials.len(), 5);
        assert_eq!(materials[0], "mat_scene09_concrete");
        assert_eq!(materials[1], "mat_scene09_windows");
        assert_eq!(materials[2], "mat_scene09_ceiling");
        assert_eq!(materials[3], "mat_scene09_plaster");
        assert_eq!(materials[4], "mat_scene09_floor");

        // 1 zone, 1 solid, 5 walls
        let zones = building.zones();
        assert_eq!(zones.len(), 1);
        let solids = building.solids();
        assert_eq!(solids.len(), 1);
        let walls = solids[0].walls();
        assert_eq!(walls.len(), 5);

        // Check polygon counts per wall (walls are sorted by name)
        let wall_map: HashMap<&str, usize> = walls
            .iter()
            .map(|w| (w.get_name(), w.polygons().len()))
            .collect();
        assert_eq!(wall_map["concrete"], 147);
        assert_eq!(wall_map["windows"], 3);
        // 4 ceiling surfaces are degenerate (collinear vertices), so 130 instead of 134
        assert_eq!(wall_map["ceiling"], 130);
        assert_eq!(wall_map["plaster"], 17);
        assert_eq!(wall_map["floor"], 29);

        // Total polygons (330 surfaces in file, minus 4 degenerate)
        let total: usize = wall_map.values().sum();
        assert_eq!(total, 326);
    }
}
