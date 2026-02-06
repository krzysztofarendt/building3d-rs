use std::collections::HashSet;

use anyhow::Result;
use rand::Rng;

use crate::geom::bboxes::bounding_box;
use crate::{Building, Point, Polygon, Vector};

use super::config::SimulationConfig;
use super::find_transparent::find_transparent_polygons;
use super::voxel_grid::VoxelGrid;

/// Flattened polygon representation for fast indexed access during simulation.
struct FlatPolygons {
    polygons: Vec<Polygon>,
    #[allow(dead_code)]
    paths: Vec<String>,
    transparent: HashSet<usize>,
    absorption: Vec<f64>,
}

/// Result of a ray tracing simulation.
pub struct SimulationResult {
    /// Ray positions per time step: positions[step][ray]
    pub positions: Vec<Vec<Point>>,
    /// Ray energies per time step: energies[step][ray]
    pub energies: Vec<Vec<f64>>,
    /// Absorber hits per time step: hits[step][absorber]
    pub hits: Vec<Vec<f64>>,
    /// Configuration used for this simulation
    pub config: SimulationConfig,
}

pub struct Simulation {
    config: SimulationConfig,
    flat: FlatPolygons,
    voxel_grid: VoxelGrid,
    bbox_min: Point,
    bbox_max: Point,
}

impl Simulation {
    pub fn new(building: &Building, config: SimulationConfig) -> Result<Self> {
        // Collect all polygons with paths
        let mut flat_polygons = Vec::new();
        let mut flat_paths = Vec::new();

        for zone in building.zones() {
            for solid in zone.solids() {
                for wall in solid.walls() {
                    for poly in wall.polygons() {
                        let path =
                            format!("{}/{}/{}/{}", zone.name, solid.name, wall.name, poly.name);
                        flat_polygons.push(poly.clone());
                        flat_paths.push(path);
                    }
                }
            }
        }

        // Map absorption values from config
        let absorption: Vec<f64> = flat_paths
            .iter()
            .map(|path| {
                config
                    .absorption
                    .get(path)
                    .copied()
                    .unwrap_or(config.default_absorption)
            })
            .collect();

        // Find transparent polygons
        let transparent_paths = if config.search_transparent {
            find_transparent_polygons(building)
        } else {
            HashSet::new()
        };

        let transparent: HashSet<usize> = flat_paths
            .iter()
            .enumerate()
            .filter(|(_, path)| transparent_paths.contains(*path))
            .map(|(idx, _)| idx)
            .collect();

        // Build voxel grid
        let poly_refs: Vec<&Polygon> = flat_polygons.iter().collect();
        let voxel_grid = VoxelGrid::new(&poly_refs, config.voxel_size);

        // Compute building bounding box
        let all_pts: Vec<Point> = flat_polygons
            .iter()
            .flat_map(|p| p.vertices().iter().copied())
            .collect();
        let (bbox_min, bbox_max) = bounding_box(&all_pts);

        let flat = FlatPolygons {
            polygons: flat_polygons,
            paths: flat_paths,
            transparent,
            absorption,
        };

        Ok(Self {
            config,
            flat,
            voxel_grid,
            bbox_min,
            bbox_max,
        })
    }

    pub fn run(self) -> SimulationResult {
        let num_rays = self.config.num_rays;
        let num_steps = self.config.num_steps;
        let num_absorbers = self.config.absorbers.len();
        let dt = self.config.time_step;
        let speed = self.config.ray_speed;
        let absorber_r2 = self.config.absorber_radius * self.config.absorber_radius;

        let mut rng = rand::thread_rng();

        // Initialize ray positions and velocities
        let mut positions: Vec<Point> = vec![self.config.source; num_rays];
        let mut velocities: Vec<Vector> = (0..num_rays)
            .map(|_| random_unit_vector(&mut rng) * speed)
            .collect();
        let mut energies: Vec<f64> = vec![1.0; num_rays];

        // Output buffers
        let mut all_positions: Vec<Vec<Point>> = Vec::with_capacity(num_steps);
        let mut all_energies: Vec<Vec<f64>> = Vec::with_capacity(num_steps);
        let mut all_hits: Vec<Vec<f64>> = Vec::with_capacity(num_steps);

        let eps = 1e-10;
        let bbox_margin = 1e-6;
        let reflection_dist = speed * dt * 1.5;

        for _step in 0..num_steps {
            let mut step_hits = vec![0.0; num_absorbers];

            // Check absorbers
            for ray in 0..num_rays {
                if energies[ray] <= eps {
                    continue;
                }
                for (ai, absorber) in self.config.absorbers.iter().enumerate() {
                    let dx = positions[ray].x - absorber.x;
                    let dy = positions[ray].y - absorber.y;
                    let dz = positions[ray].z - absorber.z;
                    let dist2 = dx * dx + dy * dy + dz * dz;
                    if dist2 <= absorber_r2 {
                        step_hits[ai] += energies[ray];
                        energies[ray] = 0.0;
                        break;
                    }
                }
            }

            // Propagate rays
            for ray in 0..num_rays {
                if energies[ray] <= eps {
                    continue;
                }

                // Boundary check (with small margin like Python version)
                let pos = positions[ray];
                if pos.x < self.bbox_min.x - bbox_margin
                    || pos.x > self.bbox_max.x + bbox_margin
                    || pos.y < self.bbox_min.y - bbox_margin
                    || pos.y > self.bbox_max.y + bbox_margin
                    || pos.z < self.bbox_min.z - bbox_margin
                    || pos.z > self.bbox_max.z + bbox_margin
                {
                    energies[ray] = 0.0;
                    continue;
                }

                // Find nearby polygons via voxel grid
                let nearby = self.voxel_grid.find_nearby(pos);
                if nearby.is_empty() {
                    // No nearby surfaces, just advance
                    let dv = velocities[ray] * dt;
                    positions[ray] = pos + dv;
                    continue;
                }

                // Find target surface (closest non-transparent polygon in ray direction)
                let vel_norm = match velocities[ray].normalize() {
                    Ok(v) => v,
                    Err(_) => {
                        energies[ray] = 0.0;
                        continue;
                    }
                };

                if let Some((target_idx, target_dist)) =
                    find_target_surface(pos, vel_norm, &nearby, &self.flat)
                {
                    if target_dist <= reflection_dist {
                        // Reflect off the surface
                        let mut current_target = Some((target_idx, target_dist));
                        // Track current search position (moves to reflection point)
                        let mut search_pos = pos;

                        while let Some((tidx, _tdist)) = current_target {
                            if energies[ray] <= eps {
                                break;
                            }

                            // Apply absorption
                            energies[ray] *= 1.0 - self.flat.absorption[tidx];
                            if energies[ray] <= eps {
                                energies[ray] = 0.0;
                                break;
                            }

                            // Reflect velocity: v = v - 2*(v·n)*n
                            let normal = self.flat.polygons[tidx].vn;
                            let dot = velocities[ray].dot(&normal);
                            velocities[ray] = velocities[ray] - 2.0 * dot * normal;

                            // Find next target after reflection
                            let new_vel_norm = match velocities[ray].normalize() {
                                Ok(v) => v,
                                Err(_) => break,
                            };

                            // Look for next target from current search position
                            current_target =
                                find_target_surface(search_pos, new_vel_norm, &nearby, &self.flat);

                            // Only continue reflecting if the next target is very close
                            if let Some((_, next_dist)) = current_target {
                                if next_dist > reflection_dist {
                                    current_target = None;
                                } else {
                                    // Move search position closer to the next surface
                                    search_pos = search_pos + new_vel_norm * (next_dist * 0.5);
                                }
                            }
                        }

                        // Advance position
                        if energies[ray] > eps {
                            let dv = velocities[ray] * dt;
                            positions[ray] = pos + dv;
                        }
                    } else {
                        // Target is far away, just advance
                        let dv = velocities[ray] * dt;
                        positions[ray] = pos + dv;
                    }
                } else {
                    // No target found, advance
                    let dv = velocities[ray] * dt;
                    positions[ray] = pos + dv;
                }
            }

            // Record state
            all_positions.push(positions.clone());
            all_energies.push(energies.clone());
            all_hits.push(step_hits);
        }

        SimulationResult {
            positions: all_positions,
            energies: all_energies,
            hits: all_hits,
            config: self.config,
        }
    }
}

/// Generate a random unit vector uniformly distributed on the sphere.
fn random_unit_vector(rng: &mut impl Rng) -> Vector {
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

/// Barycentric tolerance for ray-triangle intersection.
///
/// The Python version uses atol=1e-3 in is_point_inside_projection.
/// This tolerance prevents rays from slipping through cracks between
/// adjacent triangles in STL meshes.
const BARY_TOLERANCE: f64 = 1e-3;

/// Minimum ray parameter to avoid self-intersection after reflection.
const MIN_T: f64 = 1e-6;

/// Finds the closest non-transparent polygon in the ray's direction.
///
/// Uses Moller-Trumbore ray-triangle intersection with barycentric tolerance
/// to prevent rays from escaping through cracks between adjacent triangles.
///
/// Returns (polygon_index, distance) or None if no target found.
fn find_target_surface(
    origin: Point,
    direction: Vector,
    nearby: &HashSet<usize>,
    flat: &FlatPolygons,
) -> Option<(usize, f64)> {
    let dir = direction.normalize().ok()?;

    let mut closest: Option<(usize, f64)> = None;

    for &idx in nearby {
        if flat.transparent.contains(&idx) {
            continue;
        }

        if let Some(t) =
            ray_intersect_polygon_tolerant(origin, dir, &flat.polygons[idx], BARY_TOLERANCE)
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

/// Ray-polygon intersection using Moller-Trumbore with barycentric tolerance.
///
/// Iterates through all triangles in the polygon's mesh and returns the
/// closest intersection distance, or None if no intersection.
fn ray_intersect_polygon_tolerant(
    origin: Point,
    direction: Vector,
    polygon: &Polygon,
    tolerance: f64,
) -> Option<f64> {
    let vertices = polygon.vertices();
    let triangles = polygon.triangles()?;

    let mut best_t: Option<f64> = None;

    for tri in triangles {
        let p0 = vertices[tri.0];
        let p1 = vertices[tri.1];
        let p2 = vertices[tri.2];

        if let Some(t) = moller_trumbore(origin, direction, p0, p1, p2, tolerance) {
            match best_t {
                None => best_t = Some(t),
                Some(bt) if t < bt => best_t = Some(t),
                _ => {}
            }
        }
    }

    best_t
}

/// Moller-Trumbore ray-triangle intersection with barycentric tolerance.
///
/// The tolerance parameter allows the barycentric coordinates (u, v) to
/// slightly exceed the [0, 1] range, which prevents rays from slipping
/// through cracks between adjacent triangles in triangle meshes.
fn moller_trumbore(
    origin: Point,
    direction: Vector,
    p0: Point,
    p1: Point,
    p2: Point,
    tolerance: f64,
) -> Option<f64> {
    let edge1 = p1 - p0; // Vector
    let edge2 = p2 - p0; // Vector
    let h = direction.cross(&edge2);
    let a = edge1.dot(&h);

    if a.abs() < 1e-10 {
        return None; // Ray parallel to triangle
    }

    let f = 1.0 / a;
    let s = origin - p0; // Vector (Point - Point = Vector)
    let u = f * s.dot(&h);

    if u < -tolerance || u > 1.0 + tolerance {
        return None;
    }

    let q = s.cross(&edge1);
    let v = f * direction.dot(&q);

    if v < -tolerance || u + v > 1.0 + tolerance {
        return None;
    }

    let t = f * edge2.dot(&q);

    if t > MIN_T { Some(t) } else { None }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Solid, Zone};

    #[test]
    fn test_simulation_basic() {
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "s0").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config = SimulationConfig::new();
        config.num_steps = 10;
        config.num_rays = 5;
        config.source = Point::new(1.0, 1.0, 1.0);

        let sim = Simulation::new(&building, config).unwrap();
        let result = sim.run();

        assert_eq!(result.positions.len(), 10);
        assert_eq!(result.energies.len(), 10);
        assert_eq!(result.positions[0].len(), 5);
    }

    #[test]
    fn test_simulation_with_absorber() {
        let s0 = Solid::from_box(2.0, 2.0, 2.0, None, "s0").unwrap();
        let zone = Zone::new("z", vec![s0]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config = SimulationConfig::new();
        config.num_steps = 100;
        config.num_rays = 50;
        config.source = Point::new(1.0, 1.0, 1.0);
        config.absorbers = vec![Point::new(0.5, 0.5, 0.5)];
        config.absorber_radius = 0.3;

        let sim = Simulation::new(&building, config).unwrap();
        let result = sim.run();

        // Some rays should have been absorbed
        let total_hit: f64 = result.hits.iter().map(|h| h[0]).sum();
        assert!(total_hit > 0.0, "Some rays should hit the absorber");
    }

    #[test]
    fn test_simulation_transparent() {
        // Two adjacent boxes in the same zone — rays should pass through the shared face
        let s0 = Solid::from_box(1.0, 1.0, 1.0, None, "s0").unwrap();
        let s1 = Solid::from_box(1.0, 1.0, 1.0, Some((1.0, 0.0, 0.0)), "s1").unwrap();
        let zone = Zone::new("z", vec![s0, s1]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let mut config = SimulationConfig::new();
        config.num_steps = 200;
        config.num_rays = 200;
        config.source = Point::new(0.5, 0.5, 0.5);
        // Place absorber in the second box with larger radius
        config.absorbers = vec![Point::new(1.5, 0.5, 0.5)];
        config.absorber_radius = 0.4;

        let sim = Simulation::new(&building, config).unwrap();
        let result = sim.run();

        // Some rays should reach the absorber through the transparent face
        let total_hit: f64 = result.hits.iter().map(|h| h[0]).sum();
        assert!(
            total_hit > 0.0,
            "Rays should pass through transparent face and hit absorber"
        );
    }

    #[test]
    fn test_random_unit_vector() {
        let mut rng = rand::thread_rng();
        for _ in 0..100 {
            let v = random_unit_vector(&mut rng);
            let len = v.length();
            assert!(
                (len - 1.0).abs() < 1e-10,
                "Random unit vector should have length 1, got {}",
                len
            );
        }
    }
}
