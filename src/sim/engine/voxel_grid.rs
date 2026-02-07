use std::collections::{HashMap, HashSet};

use crate::geom::bboxes::{are_bboxes_overlapping, bounding_box};
use crate::{Point, Polygon, Vector};

pub struct VoxelGrid {
    grid: HashMap<(i32, i32, i32), Vec<usize>>,
    step: f64,
}

impl VoxelGrid {
    pub fn new(polygons: &[&Polygon], step: f64) -> Self {
        let mut grid: HashMap<(i32, i32, i32), Vec<usize>> = HashMap::new();

        if polygons.is_empty() {
            return Self { grid, step };
        }

        // Compute overall bounding box from all polygon vertices
        let all_pts: Vec<Point> = polygons
            .iter()
            .flat_map(|p| p.vertices().iter().copied())
            .collect();
        let (bbox_min, bbox_max) = bounding_box(&all_pts);

        // Compute voxel grid extent
        let imin = (bbox_min.x / step).floor() as i32 - 1;
        let jmin = (bbox_min.y / step).floor() as i32 - 1;
        let kmin = (bbox_min.z / step).floor() as i32 - 1;
        let imax = (bbox_max.x / step).ceil() as i32 + 1;
        let jmax = (bbox_max.y / step).ceil() as i32 + 1;
        let kmax = (bbox_max.z / step).ceil() as i32 + 1;

        // Precompute polygon bounding boxes
        let poly_bboxes: Vec<(Point, Point)> = polygons
            .iter()
            .map(|p| bounding_box(p.vertices()))
            .collect();

        // For each voxel cell, check which polygon bounding boxes overlap
        for i in imin..=imax {
            for j in jmin..=jmax {
                for k in kmin..=kmax {
                    let vmin = Point::new(i as f64 * step, j as f64 * step, k as f64 * step);
                    let vmax = Point::new(
                        (i + 1) as f64 * step,
                        (j + 1) as f64 * step,
                        (k + 1) as f64 * step,
                    );

                    let mut cell_polys = Vec::new();
                    for (idx, (pmin, pmax)) in poly_bboxes.iter().enumerate() {
                        if are_bboxes_overlapping(vmin, vmax, *pmin, *pmax) {
                            cell_polys.push(idx);
                        }
                    }

                    if !cell_polys.is_empty() {
                        grid.insert((i, j, k), cell_polys);
                    }
                }
            }
        }

        Self { grid, step }
    }

    /// Returns polygon indices from the voxel cell containing `pos` plus 26 neighbors.
    pub fn find_nearby(&self, pos: Point) -> HashSet<usize> {
        let ci = (pos.x / self.step).floor() as i32;
        let cj = (pos.y / self.step).floor() as i32;
        let ck = (pos.z / self.step).floor() as i32;

        let mut result = HashSet::new();

        for di in -1..=1 {
            for dj in -1..=1 {
                for dk in -1..=1 {
                    if let Some(indices) = self.grid.get(&(ci + di, cj + dj, ck + dk)) {
                        result.extend(indices);
                    }
                }
            }
        }

        result
    }

    /// Returns polygon indices from all voxel cells along a ray (Amanatides-Woo 3D-DDA).
    ///
    /// Marches through cells from `origin` in `direction` up to `max_distance`,
    /// collecting polygon indices from each visited cell.
    pub fn find_along_ray(
        &self,
        origin: Point,
        direction: Vector,
        max_distance: f64,
    ) -> HashSet<usize> {
        let mut result = HashSet::new();

        // Current cell indices
        let mut ci = (origin.x / self.step).floor() as i32;
        let mut cj = (origin.y / self.step).floor() as i32;
        let mut ck = (origin.z / self.step).floor() as i32;

        // Collect polygons from the starting cell and its neighbors
        // (to handle rays starting near cell boundaries)
        for di in -1..=1 {
            for dj in -1..=1 {
                for dk in -1..=1 {
                    if let Some(indices) = self.grid.get(&(ci + di, cj + dj, ck + dk)) {
                        result.extend(indices);
                    }
                }
            }
        }

        // Step direction: +1 or -1 per axis
        let step_i: i32 = if direction.dx > 0.0 {
            1
        } else if direction.dx < 0.0 {
            -1
        } else {
            0
        };
        let step_j: i32 = if direction.dy > 0.0 {
            1
        } else if direction.dy < 0.0 {
            -1
        } else {
            0
        };
        let step_k: i32 = if direction.dz > 0.0 {
            1
        } else if direction.dz < 0.0 {
            -1
        } else {
            0
        };

        // t_max: distance along ray to the next cell boundary for each axis
        // t_delta: distance along ray to traverse one full cell for each axis
        let (mut t_max_x, t_delta_x) = if direction.dx.abs() > 1e-30 {
            let next_boundary = if step_i > 0 {
                (ci + 1) as f64 * self.step
            } else {
                ci as f64 * self.step
            };
            let t_max = (next_boundary - origin.x) / direction.dx;
            let t_delta = (self.step / direction.dx).abs();
            (t_max, t_delta)
        } else {
            (f64::INFINITY, f64::INFINITY)
        };

        let (mut t_max_y, t_delta_y) = if direction.dy.abs() > 1e-30 {
            let next_boundary = if step_j > 0 {
                (cj + 1) as f64 * self.step
            } else {
                cj as f64 * self.step
            };
            let t_max = (next_boundary - origin.y) / direction.dy;
            let t_delta = (self.step / direction.dy).abs();
            (t_max, t_delta)
        } else {
            (f64::INFINITY, f64::INFINITY)
        };

        let (mut t_max_z, t_delta_z) = if direction.dz.abs() > 1e-30 {
            let next_boundary = if step_k > 0 {
                (ck + 1) as f64 * self.step
            } else {
                ck as f64 * self.step
            };
            let t_max = (next_boundary - origin.z) / direction.dz;
            let t_delta = (self.step / direction.dz).abs();
            (t_max, t_delta)
        } else {
            (f64::INFINITY, f64::INFINITY)
        };

        // March through grid
        loop {
            // Advance the axis with the smallest t_max
            if t_max_x < t_max_y {
                if t_max_x < t_max_z {
                    ci += step_i;
                    if t_max_x > max_distance {
                        break;
                    }
                    t_max_x += t_delta_x;
                } else {
                    ck += step_k;
                    if t_max_z > max_distance {
                        break;
                    }
                    t_max_z += t_delta_z;
                }
            } else if t_max_y < t_max_z {
                cj += step_j;
                if t_max_y > max_distance {
                    break;
                }
                t_max_y += t_delta_y;
            } else {
                ck += step_k;
                if t_max_z > max_distance {
                    break;
                }
                t_max_z += t_delta_z;
            }

            // Collect polygons from this cell
            if let Some(indices) = self.grid.get(&(ci, cj, ck)) {
                result.extend(indices);
            }

            // Safety: if all t_max are infinity, we're stuck
            if t_max_x == f64::INFINITY && t_max_y == f64::INFINITY && t_max_z == f64::INFINITY {
                break;
            }
        }

        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_voxel_grid_basic() {
        // Create a simple polygon
        let pts = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
        ];
        let poly = Polygon::new("square", pts, None).unwrap();
        let polygons: Vec<&Polygon> = vec![&poly];

        let grid = VoxelGrid::new(&polygons, 0.5);

        // Point inside the polygon area should find nearby polygons
        let nearby = grid.find_nearby(Point::new(0.5, 0.5, 0.0));
        assert!(nearby.contains(&0));
    }

    #[test]
    fn test_voxel_grid_far_point() {
        let pts = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
        ];
        let poly = Polygon::new("square", pts, None).unwrap();
        let polygons: Vec<&Polygon> = vec![&poly];

        let grid = VoxelGrid::new(&polygons, 0.5);

        // Point far away should not find any polygons
        let nearby = grid.find_nearby(Point::new(100.0, 100.0, 100.0));
        assert!(nearby.is_empty());
    }

    #[test]
    fn test_find_along_ray_hits_distant_polygon() {
        let pts = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
        ];
        let poly = Polygon::new("square", pts, None).unwrap();
        let polygons: Vec<&Polygon> = vec![&poly];

        let grid = VoxelGrid::new(&polygons, 0.5);

        // Ray from far away aimed at the polygon
        let origin = Point::new(0.5, 0.5, 50.0);
        let direction = Vector::new(0.0, 0.0, -1.0);

        // find_nearby would miss this (origin too far)
        let nearby = grid.find_nearby(origin);
        assert!(nearby.is_empty(), "find_nearby should miss distant origin");

        // find_along_ray should find it
        let along_ray = grid.find_along_ray(origin, direction, 100.0);
        assert!(
            along_ray.contains(&0),
            "find_along_ray should find the polygon"
        );
    }

    #[test]
    fn test_find_along_ray_misses_off_direction() {
        let pts = vec![
            Point::new(0.0, 0.0, 0.0),
            Point::new(1.0, 0.0, 0.0),
            Point::new(1.0, 1.0, 0.0),
            Point::new(0.0, 1.0, 0.0),
        ];
        let poly = Polygon::new("square", pts, None).unwrap();
        let polygons: Vec<&Polygon> = vec![&poly];

        let grid = VoxelGrid::new(&polygons, 0.5);

        // Ray from far away aimed in the wrong direction
        let origin = Point::new(0.5, 0.5, 50.0);
        let direction = Vector::new(0.0, 0.0, 1.0); // going away

        let along_ray = grid.find_along_ray(origin, direction, 100.0);
        assert!(
            !along_ray.contains(&0),
            "find_along_ray should not find polygon in wrong direction"
        );
    }
}
