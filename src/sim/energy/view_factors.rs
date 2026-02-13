//! View-factor based interior longwave radiation exchange.
//!
//! Provides per-surface mean radiant temperature (MRT) using Monte Carlo
//! view factors instead of area-weighted MRT. This approach is energy-conserving
//! by construction via reciprocity (`A_i * F_ij = A_j * F_ji`), which fixes the
//! energy loss observed with non-uniform `h_rad` in the star-network approach.
//!
//! Physics:
//! - `h_conv = TARP(dT, tilt)` — convection-only
//! - `h_rad = eps * 4 * sigma * T_mean^3` — **uniform** linearized radiation
//! - `T_mrt_i = sum_j(F_ij * T_j)` — per-surface from view factors
//! - `h_total = h_conv + h_rad`
//! - `t_eff = (h_conv * T_air + h_rad * T_mrt_i) / h_total`

use std::collections::HashMap;

use crate::geom::bboxes::bounding_box;
use crate::sim::engine::voxel_grid::VoxelGrid;
use crate::sim::engine::FlatScene;
use crate::sim::index::SurfaceIndex;
use crate::sim::surfaces::SurfaceSemantics;
use crate::{Building, Point, Polygon, UID, Vector};

/// Stefan–Boltzmann constant [W/(m² K⁴)].
const SIGMA: f64 = 5.670_374_419e-8;

// ─── Surface handle ─────────────────────────────────────────────────────

/// Identifies a surface participating in the view-factor enclosure.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum SurfaceHandle {
    /// A geometric polygon in the building (FVM wall, window, etc.).
    Polygon(UID),
    /// A non-geometric internal mass slab (config-defined).
    InternalMass { index: usize },
}

// ─── Per-zone data ──────────────────────────────────────────────────────

struct ViewFactorSurface {
    handle: SurfaceHandle,
    area_m2: f64,
    #[allow(dead_code)]
    cos_tilt: f64,
}

/// View factor matrix and surface metadata for one zone.
pub struct ZoneViewFactors {
    surfaces: Vec<ViewFactorSurface>,
    /// Row-major F_ij matrix [i*n + j].
    f_matrix: Vec<f64>,
    n: usize,
}

// ─── Building-level data ────────────────────────────────────────────────

/// Precomputed view factor data for all zones.
pub struct ViewFactorData {
    zones: Vec<ZoneViewFactors>,
    /// Maps each surface handle to (zone_index, surface_index_within_zone).
    /// Used during construction for internal mass matching; kept for future queries.
    #[allow(dead_code)]
    handle_to_index: HashMap<SurfaceHandle, (usize, usize)>,
}

// ─── Sampling utilities ─────────────────────────────────────────────────

/// Sample a random point uniformly on a polygon's surface.
///
/// Decomposes the polygon into triangles (reusing existing triangulation),
/// selects a triangle proportional to area, then samples via barycentric
/// coordinates.
fn sample_point_on_polygon(polygon: &Polygon, rng: &mut impl rand::Rng) -> Option<Point> {
    let verts = polygon.vertices();
    let tris = polygon.triangles()?;
    if tris.is_empty() {
        return None;
    }

    // Compute triangle areas for weighted selection.
    let areas: Vec<f64> = tris
        .iter()
        .map(|tri| {
            let a = verts[tri.0];
            let b = verts[tri.1];
            let c = verts[tri.2];
            let ab = Vector::new(b.x - a.x, b.y - a.y, b.z - a.z);
            let ac = Vector::new(c.x - a.x, c.y - a.y, c.z - a.z);
            0.5 * ab.cross(&ac).length()
        })
        .collect();

    let total_area: f64 = areas.iter().sum();
    if total_area <= 0.0 {
        return None;
    }

    // Select triangle proportional to area.
    let r: f64 = rng.r#gen::<f64>() * total_area;
    let mut cumulative = 0.0;
    let mut selected = tris.len() - 1;
    for (i, &a) in areas.iter().enumerate() {
        cumulative += a;
        if r <= cumulative {
            selected = i;
            break;
        }
    }

    let tri = &tris[selected];
    let a = verts[tri.0];
    let b = verts[tri.1];
    let c = verts[tri.2];

    // Uniform barycentric sampling.
    let u: f64 = rng.r#gen();
    let v: f64 = rng.r#gen();
    let (su, sv) = if u + v > 1.0 {
        (1.0 - u, 1.0 - v)
    } else {
        (u, v)
    };
    let w = 1.0 - su - sv;

    Some(Point::new(
        a.x * w + b.x * su + c.x * sv,
        a.y * w + b.y * su + c.y * sv,
        a.z * w + b.z * su + c.z * sv,
    ))
}

/// Generate a cosine-weighted hemisphere direction (Malley's method).
///
/// Pattern from `src/sim/engine/reflection.rs` Diffuse::reflect.
fn cosine_weighted_hemisphere_dir(normal: Vector, rng: &mut impl rand::Rng) -> Option<Vector> {
    let n = normal.normalize().ok()?;

    // Build orthonormal basis around the normal.
    let arbitrary = if n.dx.abs() < 0.9 {
        Vector::new(1.0, 0.0, 0.0)
    } else {
        Vector::new(0.0, 1.0, 0.0)
    };
    let tangent = n
        .cross(&arbitrary)
        .normalize()
        .unwrap_or(Vector::new(1.0, 0.0, 0.0));
    let bitangent = n.cross(&tangent);

    // Malley's method: sample uniformly on disk, project onto hemisphere.
    let u1: f64 = rng.r#gen();
    let u2: f64 = rng.r#gen();
    let r = u1.sqrt();
    let phi = 2.0 * std::f64::consts::PI * u2;
    let x = r * phi.cos();
    let y = r * phi.sin();
    let z = (1.0 - u1).sqrt();

    Some(tangent * x + bitangent * y + n * z)
}

// ─── Zone-local FlatScene ───────────────────────────────────────────────

/// Build a FlatScene containing only the specified polygons (all opaque).
fn build_zone_flat_scene(polygons: &[Polygon], voxel_size: f64) -> FlatScene {
    let poly_refs: Vec<&Polygon> = polygons.iter().collect();
    let voxel_grid = VoxelGrid::new(&poly_refs, voxel_size);

    let all_pts: Vec<Point> = polygons
        .iter()
        .flat_map(|p| p.vertices().iter().copied())
        .collect();
    let (bbox_min, bbox_max) = if all_pts.is_empty() {
        (Point::new(0.0, 0.0, 0.0), Point::new(1.0, 1.0, 1.0))
    } else {
        bounding_box(&all_pts)
    };

    FlatScene {
        polygons: polygons.to_vec(),
        paths: vec![String::new(); polygons.len()],
        transparent: std::collections::HashSet::new(),
        voxel_grid,
        bbox_min,
        bbox_max,
    }
}

// ─── View factor computation ────────────────────────────────────────────

/// Compute view factors for one zone via Monte Carlo ray casting.
///
/// For each surface i, shoot `rays_per_surface` cosine-weighted rays from
/// random points on i and count hits on each other surface j.
/// `F_ij = count_ij / N`.
fn compute_zone_view_factors(
    zone_polygons: &[Polygon],
    zone_polygon_uids: &[UID],
    zone_polygon_cos_tilts: &[f64],
    rays_per_surface: usize,
) -> ZoneViewFactors {
    let n = zone_polygons.len();

    let surfaces: Vec<ViewFactorSurface> = (0..n)
        .map(|i| ViewFactorSurface {
            handle: SurfaceHandle::Polygon(zone_polygon_uids[i].clone()),
            area_m2: zone_polygons[i].area(),
            cos_tilt: zone_polygon_cos_tilts[i],
        })
        .collect();

    if n == 0 {
        return ZoneViewFactors {
            surfaces,
            f_matrix: vec![],
            n: 0,
        };
    }

    // Build per-zone FlatScene for ray casting.
    let voxel_size = estimate_voxel_size(zone_polygons);
    let scene = build_zone_flat_scene(zone_polygons, voxel_size);

    // Map from polygon UID to index in zone_polygons.
    let uid_to_idx: HashMap<&UID, usize> = zone_polygon_uids
        .iter()
        .enumerate()
        .map(|(i, uid)| (uid, i))
        .collect();

    let mut f_matrix = vec![0.0_f64; n * n];
    let mut rng = rand::thread_rng();

    for i in 0..n {
        let poly_i = &zone_polygons[i];
        // Use the INWARD-facing normal for interior radiation exchange.
        // Building polygons have outward-facing normals by convention.
        let normal_i = Vector::new(-poly_i.vn.dx, -poly_i.vn.dy, -poly_i.vn.dz);
        let mut hit_count = vec![0_u64; n];
        let mut valid_rays = 0_u64;

        for _ in 0..rays_per_surface {
            let Some(origin) = sample_point_on_polygon(poly_i, &mut rng) else {
                continue;
            };
            let Some(dir) = cosine_weighted_hemisphere_dir(normal_i, &mut rng) else {
                continue;
            };

            // Offset origin slightly along inward normal to avoid self-intersection.
            let origin_offset = Point::new(
                origin.x + normal_i.dx * 1e-4,
                origin.y + normal_i.dy * 1e-4,
                origin.z + normal_i.dz * 1e-4,
            );

            if let Some((hit_idx, _dist)) = scene.find_target_surface(origin_offset, dir) {
                // Map hit polygon back to our surface index.
                let hit_uid = &scene.polygons[hit_idx].uid;
                if let Some(&j) = uid_to_idx.get(hit_uid)
                    && j != i
                {
                    hit_count[j] += 1;
                }
                valid_rays += 1;
            } else {
                valid_rays += 1; // Ray escaped — counts towards denominator.
            }
        }

        if valid_rays > 0 {
            for j in 0..n {
                f_matrix[i * n + j] = hit_count[j] as f64 / valid_rays as f64;
            }
        }
    }

    enforce_reciprocity_and_normalize(&mut f_matrix, &surfaces, n);

    ZoneViewFactors {
        surfaces,
        f_matrix,
        n,
    }
}

/// Estimate a suitable voxel size from polygon bounding box.
fn estimate_voxel_size(polygons: &[Polygon]) -> f64 {
    if polygons.is_empty() {
        return 1.0;
    }
    let all_pts: Vec<Point> = polygons
        .iter()
        .flat_map(|p| p.vertices().iter().copied())
        .collect();
    let (mn, mx) = bounding_box(&all_pts);
    let diag = ((mx.x - mn.x).powi(2) + (mx.y - mn.y).powi(2) + (mx.z - mn.z).powi(2)).sqrt();
    (diag / 10.0).max(0.1)
}

/// Enforce reciprocity (A_i F_ij = A_j F_ji) and normalize rows.
fn enforce_reciprocity_and_normalize(
    f_matrix: &mut [f64],
    surfaces: &[ViewFactorSurface],
    n: usize,
) {
    if n <= 1 {
        return;
    }

    // Step 1: Symmetrize via reciprocity.
    // Average A_i*F_ij and A_j*F_ji, then redistribute.
    for i in 0..n {
        for j in (i + 1)..n {
            let ai = surfaces[i].area_m2;
            let aj = surfaces[j].area_m2;
            if ai <= 0.0 || aj <= 0.0 {
                continue;
            }
            let af_ij = ai * f_matrix[i * n + j];
            let af_ji = aj * f_matrix[j * n + i];
            let avg = 0.5 * (af_ij + af_ji);
            f_matrix[i * n + j] = avg / ai;
            f_matrix[j * n + i] = avg / aj;
        }
    }

    // Step 2: Zero out diagonal (no self-view for flat/convex surfaces).
    for i in 0..n {
        f_matrix[i * n + i] = 0.0;
    }

    // Step 3: Normalize rows to sum to 1.0 (enclosure assumption).
    for i in 0..n {
        let row_sum: f64 = (0..n).map(|j| f_matrix[i * n + j]).sum();
        if row_sum > 1e-15 {
            let scale = 1.0 / row_sum;
            for j in 0..n {
                f_matrix[i * n + j] *= scale;
            }
        }
    }
}

// ─── Internal mass matching ─────────────────────────────────────────────

/// Match internal mass surfaces to geometric polygons in the same zone.
///
/// - Floor mass (cos_tilt <= -0.5) → match to polygon with vn.dz <= -0.5
/// - Ceiling mass (cos_tilt >= 0.5) → match to polygon with vn.dz >= 0.5
///
/// Returns mapping from mass surface index to polygon UID.
fn match_internal_mass_to_polygons(
    mass_surfaces: &[(usize, f64)], // (mass_index, cos_tilt)
    zone_polygons: &[Polygon],
    zone_polygon_uids: &[UID],
) -> HashMap<usize, UID> {
    let mut result = HashMap::new();

    for &(mass_idx, cos_tilt) in mass_surfaces {
        if cos_tilt <= -0.5 {
            // Floor: match to polygon with downward-facing normal (interior side up).
            for (i, poly) in zone_polygons.iter().enumerate() {
                if poly.vn.dz <= -0.5 {
                    result.insert(mass_idx, zone_polygon_uids[i].clone());
                    break;
                }
            }
        } else if cos_tilt >= 0.5 {
            // Ceiling: match to polygon with upward-facing normal (interior side down).
            for (i, poly) in zone_polygons.iter().enumerate() {
                if poly.vn.dz >= 0.5 {
                    result.insert(mass_idx, zone_polygon_uids[i].clone());
                    break;
                }
            }
        }
        // Vertical mass: no match → will use equal view factors.
    }

    result
}

// ─── Per-surface MRT ────────────────────────────────────────────────────

/// Compute per-surface mean radiant temperature from view factors.
///
/// `T_mrt_i = sum_j(F_ij * T_j)` for each surface.
pub fn compute_per_surface_mrt(
    vf_data: &ViewFactorData,
    surface_temps: &HashMap<SurfaceHandle, f64>,
) -> HashMap<SurfaceHandle, f64> {
    let mut result = HashMap::new();

    for zone_vf in &vf_data.zones {
        let n = zone_vf.n;
        for i in 0..n {
            let handle_i = &zone_vf.surfaces[i].handle;
            let mut mrt = 0.0;
            for j in 0..n {
                let f_ij = zone_vf.f_matrix[i * n + j];
                if f_ij <= 0.0 {
                    continue;
                }
                let handle_j = &zone_vf.surfaces[j].handle;
                let t_j = surface_temps.get(handle_j).copied().unwrap_or(20.0);
                mrt += f_ij * t_j;
            }
            result.insert(handle_i.clone(), mrt);
        }
    }

    result
}

/// Linearized radiative heat transfer coefficient [W/(m² K)].
///
/// `h_rad = 4 * eps * sigma * (T_mean + 273.15)^3`
///
/// This is **uniform** for all surfaces (unlike the failed star-network
/// which used non-uniform `h_iso - h_conv`). Uniform h_rad is key to
/// energy conservation via reciprocity.
pub fn linearized_h_rad(emissivity: f64, t_mean_c: f64) -> f64 {
    let t_k = (t_mean_c + 273.15).max(1.0);
    4.0 * emissivity * SIGMA * t_k.powi(3)
}

// ─── Top-level builder ──────────────────────────────────────────────────

/// Information about an internal mass surface needed for view factor computation.
pub struct InternalMassInfo {
    /// Index of this mass surface in the simulation's internal_mass_surfaces vec.
    pub index: usize,
    /// Zone UID that this mass belongs to.
    pub zone_uid: UID,
    /// Cosine of surface tilt.
    pub cos_tilt: f64,
    /// Exposed face area [m²].
    pub area_m2: f64,
}

/// Compute view factors for the entire building.
///
/// Builds a per-zone FlatScene containing only that zone's interior-facing
/// polygons, runs Monte Carlo ray casting to compute F_ij, and matches
/// internal mass surfaces to geometric polygons.
pub fn compute_building_view_factors(
    building: &Building,
    _index: &SurfaceIndex,
    boundaries: &SurfaceSemantics,
    mass_infos: &[InternalMassInfo],
    rays_per_surface: usize,
) -> ViewFactorData {
    let mut zones = Vec::new();
    let mut handle_to_index = HashMap::new();

    // Group mass surfaces by zone UID.
    let mut mass_by_zone: HashMap<UID, Vec<(usize, f64, f64)>> = HashMap::new(); // (idx, cos_tilt, area)
    for m in mass_infos {
        mass_by_zone
            .entry(m.zone_uid.clone())
            .or_default()
            .push((m.index, m.cos_tilt, m.area_m2));
    }

    for zone in building.zones() {
        // Collect interior-facing polygons for this zone.
        // These are exterior surfaces + inter-zone interfaces (NOT same-zone interfaces).
        let mut zone_polygons: Vec<Polygon> = Vec::new();
        let mut zone_polygon_uids: Vec<UID> = Vec::new();
        let mut zone_polygon_cos_tilts: Vec<f64> = Vec::new();

        for solid in zone.solids() {
            for wall in solid.walls() {
                for poly in wall.polygons() {
                    // Include exterior surfaces and inter-zone interfaces.
                    // Exclude same-zone interfaces (transparent internal partitions).
                    if boundaries.is_same_zone_interface(&poly.uid) {
                        continue;
                    }
                    zone_polygons.push(poly.clone());
                    zone_polygon_uids.push(poly.uid.clone());
                    zone_polygon_cos_tilts.push(poly.vn.dz);
                }
            }
        }

        if zone_polygons.is_empty() {
            continue;
        }

        // Compute geometric view factors.
        let mut zone_vf = compute_zone_view_factors(
            &zone_polygons,
            &zone_polygon_uids,
            &zone_polygon_cos_tilts,
            rays_per_surface,
        );

        let zone_idx = zones.len();

        // Register polygon handles.
        for (i, surf) in zone_vf.surfaces.iter().enumerate() {
            handle_to_index.insert(surf.handle.clone(), (zone_idx, i));
        }

        // Match internal mass surfaces to geometric polygons.
        if let Some(mass_list) = mass_by_zone.get(&zone.uid) {
            let mass_cos_tilts: Vec<(usize, f64)> =
                mass_list.iter().map(|&(idx, ct, _)| (idx, ct)).collect();
            let matches =
                match_internal_mass_to_polygons(&mass_cos_tilts, &zone_polygons, &zone_polygon_uids);

            for &(mass_idx, _cos_tilt, area_m2) in mass_list {
                let handle = SurfaceHandle::InternalMass { index: mass_idx };
                if let Some(matched_uid) = matches.get(&mass_idx) {
                    // Mass inherits the matched polygon's position in the matrix.
                    let matched_handle = SurfaceHandle::Polygon(matched_uid.clone());
                    if let Some(&(zi, si)) = handle_to_index.get(&matched_handle) {
                        handle_to_index.insert(handle, (zi, si));
                    }
                } else {
                    // Unmatched mass: add as a new surface with equal view factors.
                    let n_old = zone_vf.n;
                    let n_new = n_old + 1;
                    let new_idx = n_old;

                    // Expand matrix.
                    let mut new_matrix = vec![0.0; n_new * n_new];
                    for i in 0..n_old {
                        for j in 0..n_old {
                            new_matrix[i * n_new + j] = zone_vf.f_matrix[i * n_old + j];
                        }
                    }

                    // Equal view factors for the new surface.
                    if n_old > 0 {
                        let f_equal = 1.0 / n_old as f64;
                        for j in 0..n_old {
                            new_matrix[new_idx * n_new + j] = f_equal;
                        }
                        // Other surfaces see the new mass proportionally to its area.
                        // Redistribute: scale existing rows and add new column.
                        let total_area: f64 =
                            zone_vf.surfaces.iter().map(|s| s.area_m2).sum::<f64>() + area_m2;
                        if total_area > 0.0 {
                            let mass_fraction = area_m2 / total_area;
                            for i in 0..n_old {
                                let old_sum: f64 =
                                    (0..n_old).map(|j| new_matrix[i * n_new + j]).sum();
                                new_matrix[i * n_new + new_idx] = old_sum * mass_fraction;
                                // Re-normalize row.
                                let row_sum: f64 =
                                    (0..n_new).map(|j| new_matrix[i * n_new + j]).sum();
                                if row_sum > 1e-15 {
                                    let scale = 1.0 / row_sum;
                                    for j in 0..n_new {
                                        new_matrix[i * n_new + j] *= scale;
                                    }
                                }
                            }
                        }
                    }

                    zone_vf.surfaces.push(ViewFactorSurface {
                        handle: handle.clone(),
                        area_m2,
                        cos_tilt: _cos_tilt,
                    });
                    zone_vf.f_matrix = new_matrix;
                    zone_vf.n = n_new;

                    handle_to_index.insert(handle, (zone_idx, new_idx));
                }
            }
        }

        zones.push(zone_vf);
    }

    ViewFactorData {
        zones,
        handle_to_index,
    }
}

// ─── Tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn make_unit_cube_polygons() -> Vec<Polygon> {
        // Create 6 faces of a unit cube [0,1]^3, normals pointing OUTWARD
        // (matching the building convention — the VF code negates them internally).
        let faces = [
            // Floor (z=0, normal = -z outward)
            (
                vec![
                    Point::new(0.0, 0.0, 0.0),
                    Point::new(0.0, 1.0, 0.0),
                    Point::new(1.0, 1.0, 0.0),
                    Point::new(1.0, 0.0, 0.0),
                ],
                Vector::new(0.0, 0.0, -1.0),
                "floor",
            ),
            // Ceiling (z=1, normal = +z outward)
            (
                vec![
                    Point::new(0.0, 0.0, 1.0),
                    Point::new(1.0, 0.0, 1.0),
                    Point::new(1.0, 1.0, 1.0),
                    Point::new(0.0, 1.0, 1.0),
                ],
                Vector::new(0.0, 0.0, 1.0),
                "ceiling",
            ),
            // Front wall (y=0, normal = -y outward)
            (
                vec![
                    Point::new(0.0, 0.0, 0.0),
                    Point::new(1.0, 0.0, 0.0),
                    Point::new(1.0, 0.0, 1.0),
                    Point::new(0.0, 0.0, 1.0),
                ],
                Vector::new(0.0, -1.0, 0.0),
                "front",
            ),
            // Back wall (y=1, normal = +y outward)
            (
                vec![
                    Point::new(0.0, 1.0, 0.0),
                    Point::new(0.0, 1.0, 1.0),
                    Point::new(1.0, 1.0, 1.0),
                    Point::new(1.0, 1.0, 0.0),
                ],
                Vector::new(0.0, 1.0, 0.0),
                "back",
            ),
            // Left wall (x=0, normal = -x outward)
            (
                vec![
                    Point::new(0.0, 0.0, 0.0),
                    Point::new(0.0, 0.0, 1.0),
                    Point::new(0.0, 1.0, 1.0),
                    Point::new(0.0, 1.0, 0.0),
                ],
                Vector::new(-1.0, 0.0, 0.0),
                "left",
            ),
            // Right wall (x=1, normal = +x outward)
            (
                vec![
                    Point::new(1.0, 0.0, 0.0),
                    Point::new(1.0, 1.0, 0.0),
                    Point::new(1.0, 1.0, 1.0),
                    Point::new(1.0, 0.0, 1.0),
                ],
                Vector::new(1.0, 0.0, 0.0),
                "right",
            ),
        ];

        faces
            .into_iter()
            .map(|(pts, normal, name)| Polygon::new(name, pts, Some(normal)).unwrap())
            .collect()
    }

    #[test]
    fn test_cube_view_factors_symmetry() {
        let polygons = make_unit_cube_polygons();
        let uids: Vec<UID> = polygons.iter().map(|p| p.uid.clone()).collect();
        let cos_tilts: Vec<f64> = polygons.iter().map(|p| p.vn.dz).collect();

        let vf = compute_zone_view_factors(&polygons, &uids, &cos_tilts, 50_000);

        // For a unit cube, F_ij should be ~0.2 for all i != j.
        let n = vf.n;
        assert_eq!(n, 6);
        for i in 0..n {
            for j in 0..n {
                let f = vf.f_matrix[i * n + j];
                if i == j {
                    assert!(
                        f.abs() < 1e-12,
                        "Self view factor F[{i},{j}] = {f}, expected 0"
                    );
                } else {
                    assert!(
                        (f - 0.2).abs() < 0.03,
                        "Cube F[{i},{j}] = {f:.4}, expected ~0.2"
                    );
                }
            }
        }
    }

    #[test]
    fn test_reciprocity() {
        let polygons = make_unit_cube_polygons();
        let uids: Vec<UID> = polygons.iter().map(|p| p.uid.clone()).collect();
        let cos_tilts: Vec<f64> = polygons.iter().map(|p| p.vn.dz).collect();

        let vf = compute_zone_view_factors(&polygons, &uids, &cos_tilts, 20_000);

        let n = vf.n;
        for i in 0..n {
            for j in (i + 1)..n {
                let ai = vf.surfaces[i].area_m2;
                let aj = vf.surfaces[j].area_m2;
                let af_ij = ai * vf.f_matrix[i * n + j];
                let af_ji = aj * vf.f_matrix[j * n + i];
                assert!(
                    (af_ij - af_ji).abs() < 0.01,
                    "Reciprocity violated: A[{i}]*F[{i},{j}] = {af_ij:.4}, A[{j}]*F[{j},{i}] = {af_ji:.4}"
                );
            }
        }
    }

    #[test]
    fn test_row_sums() {
        let polygons = make_unit_cube_polygons();
        let uids: Vec<UID> = polygons.iter().map(|p| p.uid.clone()).collect();
        let cos_tilts: Vec<f64> = polygons.iter().map(|p| p.vn.dz).collect();

        let vf = compute_zone_view_factors(&polygons, &uids, &cos_tilts, 20_000);

        let n = vf.n;
        for i in 0..n {
            let row_sum: f64 = (0..n).map(|j| vf.f_matrix[i * n + j]).sum();
            assert!(
                (row_sum - 1.0).abs() < 0.01,
                "Row {i} sum = {row_sum:.4}, expected ~1.0"
            );
        }
    }

    #[test]
    fn test_energy_conservation() {
        let polygons = make_unit_cube_polygons();
        let uids: Vec<UID> = polygons.iter().map(|p| p.uid.clone()).collect();
        let cos_tilts: Vec<f64> = polygons.iter().map(|p| p.vn.dz).collect();

        let vf = compute_zone_view_factors(&polygons, &uids, &cos_tilts, 20_000);

        // Assign varying temperatures.
        let temps = [18.0, 22.0, 20.0, 21.0, 19.0, 23.0];
        let h_r = linearized_h_rad(0.9, 20.0);

        let n = vf.n;
        let mut total_exchange = 0.0;
        for i in 0..n {
            let ai = vf.surfaces[i].area_m2;
            let t_i = temps[i];
            let mut t_mrt_i = 0.0;
            for j in 0..n {
                t_mrt_i += vf.f_matrix[i * n + j] * temps[j];
            }
            total_exchange += ai * h_r * (t_i - t_mrt_i);
        }

        assert!(
            total_exchange.abs() < 1.0,
            "Energy conservation violated: sum = {total_exchange:.4} W, expected ~0"
        );
    }

    #[test]
    fn test_point_sampling() {
        let poly = Polygon::new(
            "test",
            vec![
                Point::new(0.0, 0.0, 0.0),
                Point::new(2.0, 0.0, 0.0),
                Point::new(2.0, 3.0, 0.0),
                Point::new(0.0, 3.0, 0.0),
            ],
            None,
        )
        .unwrap();

        let mut rng = rand::thread_rng();
        for _ in 0..1000 {
            let pt = sample_point_on_polygon(&poly, &mut rng).unwrap();
            assert!(
                pt.x >= -1e-10 && pt.x <= 2.0 + 1e-10,
                "Point x={} out of range",
                pt.x
            );
            assert!(
                pt.y >= -1e-10 && pt.y <= 3.0 + 1e-10,
                "Point y={} out of range",
                pt.y
            );
            assert!((pt.z).abs() < 1e-10, "Point z={} should be ~0", pt.z);
        }
    }

    #[test]
    fn test_linearized_h_rad() {
        let h = linearized_h_rad(0.9, 20.0);
        // At 20°C: 4 * 0.9 * 5.67e-8 * (293.15)^3 ≈ 5.13 W/(m²·K)
        assert!(
            (h - 5.13).abs() < 0.1,
            "h_rad = {h:.3}, expected ~5.13 W/(m²·K)"
        );
    }

    #[test]
    fn test_per_surface_mrt() {
        let polygons = make_unit_cube_polygons();
        let uids: Vec<UID> = polygons.iter().map(|p| p.uid.clone()).collect();
        let cos_tilts: Vec<f64> = polygons.iter().map(|p| p.vn.dz).collect();

        let vf = compute_zone_view_factors(&polygons, &uids, &cos_tilts, 20_000);

        // Build ViewFactorData.
        let mut handle_to_index = HashMap::new();
        for (i, surf) in vf.surfaces.iter().enumerate() {
            handle_to_index.insert(surf.handle.clone(), (0, i));
        }
        let vf_data = ViewFactorData {
            zones: vec![vf],
            handle_to_index,
        };

        // All surfaces at the same temperature → MRT = that temperature.
        let mut temps = HashMap::new();
        for uid in &uids {
            temps.insert(SurfaceHandle::Polygon(uid.clone()), 25.0);
        }

        let mrts = compute_per_surface_mrt(&vf_data, &temps);
        for uid in &uids {
            let handle = SurfaceHandle::Polygon(uid.clone());
            let mrt = mrts.get(&handle).copied().unwrap_or(0.0);
            assert!(
                (mrt - 25.0).abs() < 0.5,
                "MRT = {mrt:.2}, expected ~25.0"
            );
        }
    }
}
