pub mod absorption;
pub mod find_transparent;
pub mod propagation;
pub mod reflection;
pub mod voxel_grid;

use std::collections::HashSet;

use crate::geom::bboxes::bounding_box;
use crate::geom::ray::intersect_triangles_tolerant;
use crate::{Building, Point, Polygon, Vector};

use self::voxel_grid::VoxelGrid;

/// Barycentric tolerance for ray-triangle intersection.
const BARY_TOLERANCE: f64 = 1e-3;

/// Flattened scene representation for fast indexed access during simulation.
pub struct FlatScene {
    /// All polygons in the scene.
    pub polygons: Vec<Polygon>,
    /// Path for each polygon (zone/solid/wall/polygon).
    pub paths: Vec<String>,
    /// Indices of transparent polygons (internal interfaces).
    pub transparent: HashSet<usize>,
    /// Spatial acceleration structure.
    pub voxel_grid: VoxelGrid,
    /// Scene bounding box minimum.
    pub bbox_min: Point,
    /// Scene bounding box maximum.
    pub bbox_max: Point,
}

impl FlatScene {
    /// Creates a flat scene from a building hierarchy.
    pub fn new(building: &Building, voxel_size: f64, search_transparent: bool) -> Self {
        let mut polygons = Vec::new();
        let mut paths = Vec::new();

        for zone in building.zones() {
            for solid in zone.solids() {
                for wall in solid.walls() {
                    for poly in wall.polygons() {
                        let path =
                            format!("{}/{}/{}/{}", zone.name, solid.name, wall.name, poly.name);
                        polygons.push(poly.clone());
                        paths.push(path);
                    }
                }
            }
        }

        let transparent = if search_transparent {
            let transparent_paths = find_transparent::find_transparent_polygons(building);
            paths
                .iter()
                .enumerate()
                .filter(|(_, path)| transparent_paths.contains(*path))
                .map(|(idx, _)| idx)
                .collect()
        } else {
            HashSet::new()
        };

        let poly_refs: Vec<&Polygon> = polygons.iter().collect();
        let voxel_grid = VoxelGrid::new(&poly_refs, voxel_size);

        let all_pts: Vec<Point> = polygons
            .iter()
            .flat_map(|p| p.vertices().iter().copied())
            .collect();
        let (bbox_min, bbox_max) = bounding_box(&all_pts);

        Self {
            polygons,
            paths,
            transparent,
            voxel_grid,
            bbox_min,
            bbox_max,
        }
    }

    /// Finds the closest non-transparent polygon in the ray's direction.
    ///
    /// Uses 3D-DDA ray marching through the voxel grid to find candidate
    /// polygons along the ray path, up to the scene diagonal distance.
    ///
    /// Returns (polygon_index, distance) or None if no target found.
    pub fn find_target_surface(&self, origin: Point, direction: Vector) -> Option<(usize, f64)> {
        let dir = direction.normalize().ok()?;

        // Compute max search distance: from origin to the farthest corner of the bbox.
        // This handles rays starting both inside and far outside the scene.
        let corners = [
            self.bbox_min,
            self.bbox_max,
            Point::new(self.bbox_min.x, self.bbox_min.y, self.bbox_max.z),
            Point::new(self.bbox_min.x, self.bbox_max.y, self.bbox_min.z),
            Point::new(self.bbox_max.x, self.bbox_min.y, self.bbox_min.z),
            Point::new(self.bbox_min.x, self.bbox_max.y, self.bbox_max.z),
            Point::new(self.bbox_max.x, self.bbox_min.y, self.bbox_max.z),
            Point::new(self.bbox_max.x, self.bbox_max.y, self.bbox_min.z),
        ];
        let max_dist = corners
            .iter()
            .map(|c| {
                let d = Vector::new(c.x - origin.x, c.y - origin.y, c.z - origin.z);
                d.length()
            })
            .fold(0.0_f64, f64::max)
            + 1.0; // small margin

        let candidates = self.voxel_grid.find_along_ray(origin, dir, max_dist);

        if candidates.is_empty() {
            return None;
        }

        let mut closest: Option<(usize, f64)> = None;

        for &idx in &candidates {
            if self.transparent.contains(&idx) {
                continue;
            }

            let vertices = self.polygons[idx].vertices();
            let triangles = match self.polygons[idx].triangles() {
                Some(t) => t,
                None => continue,
            };

            if let Some(t) =
                intersect_triangles_tolerant(origin, dir, vertices, triangles, BARY_TOLERANCE)
            {
                match closest {
                    None => closest = Some((idx, t)),
                    Some((_, best_t)) if t < best_t => closest = Some((idx, t)),
                    _ => {}
                }
            }
        }

        closest
    }

    /// Finds the closest non-transparent polygon using brute-force search.
    ///
    /// Use this when the ray origin may be far from surfaces (e.g., lighting).
    /// Returns (polygon_index, distance) or None if no target found.
    pub fn find_target_surface_global(
        &self,
        origin: Point,
        direction: Vector,
    ) -> Option<(usize, f64)> {
        let dir = direction.normalize().ok()?;

        let mut closest: Option<(usize, f64)> = None;

        for (idx, polygon) in self.polygons.iter().enumerate() {
            if self.transparent.contains(&idx) {
                continue;
            }

            let vertices = polygon.vertices();
            let triangles = match polygon.triangles() {
                Some(t) => t,
                None => continue,
            };

            if let Some(t) =
                intersect_triangles_tolerant(origin, dir, vertices, triangles, BARY_TOLERANCE)
            {
                match closest {
                    None => closest = Some((idx, t)),
                    Some((_, best_t)) if t < best_t => closest = Some((idx, t)),
                    _ => {}
                }
            }
        }

        closest
    }

    /// Checks if a point is within the scene bounding box (with margin).
    pub fn is_in_bounds(&self, pos: Point, margin: f64) -> bool {
        pos.x >= self.bbox_min.x - margin
            && pos.x <= self.bbox_max.x + margin
            && pos.y >= self.bbox_min.y - margin
            && pos.y <= self.bbox_max.y + margin
            && pos.z >= self.bbox_min.z - margin
            && pos.z <= self.bbox_max.z + margin
    }
}

/// State of a single ray during simulation.
#[derive(Debug, Clone)]
pub struct RayState {
    pub position: Point,
    pub velocity: Vector,
    pub energy: f64,
}

/// A batch of rays for simulation.
pub struct RayBatch {
    pub rays: Vec<RayState>,
}

impl RayBatch {
    /// Creates a new batch of rays from a source point with random directions.
    pub fn new(source: Point, speed: f64, num_rays: usize) -> Self {
        let mut rng = rand::thread_rng();
        let rays = (0..num_rays)
            .map(|_| {
                let dir = random_unit_vector(&mut rng);
                RayState {
                    position: source,
                    velocity: dir * speed,
                    energy: 1.0,
                }
            })
            .collect();
        Self { rays }
    }

    /// Number of rays in the batch.
    pub fn len(&self) -> usize {
        self.rays.len()
    }

    /// Returns true if the batch is empty.
    pub fn is_empty(&self) -> bool {
        self.rays.is_empty()
    }
}

/// Generate a random unit vector uniformly distributed on the sphere.
fn random_unit_vector(rng: &mut impl rand::Rng) -> Vector {
    loop {
        let x: f64 = rng.gen_range(-1.0..1.0);
        let y: f64 = rng.gen_range(-1.0..1.0);
        let z: f64 = rng.gen_range(-1.0..1.0);
        let len2 = x * x + y * y + z * z;
        if len2 > 1e-6 && len2 <= 1.0 {
            let len = len2.sqrt();
            return Vector::new(x / len, y / len, z / len);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Solid, Zone};

    #[test]
    fn test_flat_scene_construction() {
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "s0").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let scene = FlatScene::new(&building, 0.5, false);
        assert!(!scene.polygons.is_empty());
        assert_eq!(scene.polygons.len(), scene.paths.len());
        assert!(scene.transparent.is_empty());
    }

    #[test]
    fn test_flat_scene_transparent() {
        let s0 = Solid::from_box(1.0, 1.0, 1.0, None, "s0").unwrap();
        let s1 = Solid::from_box(1.0, 1.0, 1.0, Some((1.0, 0.0, 0.0)), "s1").unwrap();
        let zone = Zone::new("z", vec![s0, s1]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let scene = FlatScene::new(&building, 0.5, true);
        assert!(!scene.transparent.is_empty());
    }

    #[test]
    fn test_flat_scene_find_target() {
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "s0").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let scene = FlatScene::new(&building, 0.5, false);

        // From center, should find a wall in any direction
        let origin = Point::new(1.0, 1.0, 1.0);
        let direction = Vector::new(1.0, 0.0, 0.0);
        let result = scene.find_target_surface(origin, direction);
        assert!(result.is_some());
    }

    #[test]
    fn test_find_target_surface_from_outside() {
        // A ray from far outside should find the box surface
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "s0").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let scene = FlatScene::new(&building, 0.5, false);

        // Ray from far away aimed at the box
        let origin = Point::new(1.0, 1.0, 50.0);
        let direction = Vector::new(0.0, 0.0, -1.0);
        let result = scene.find_target_surface(origin, direction);
        assert!(
            result.is_some(),
            "Should find box surface from distant origin"
        );

        if let Some((_, dist)) = result {
            // Distance from z=50 to z=2 surface should be ~48
            assert!(
                (dist - 48.0).abs() < 0.1,
                "Distance should be ~48, got {dist}"
            );
        }
    }

    #[test]
    fn test_ray_batch_creation() {
        let batch = RayBatch::new(Point::new(0.0, 0.0, 0.0), 343.0, 10);
        assert_eq!(batch.len(), 10);
        assert!(!batch.is_empty());
        for ray in &batch.rays {
            let speed = ray.velocity.length();
            assert!((speed - 343.0).abs() < 1e-6);
        }
    }

    #[test]
    fn test_is_in_bounds() {
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "s0").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let scene = FlatScene::new(&building, 0.5, false);

        assert!(scene.is_in_bounds(Point::new(1.0, 1.0, 1.0), 0.0));
        assert!(!scene.is_in_bounds(Point::new(10.0, 10.0, 10.0), 0.0));
    }
}
