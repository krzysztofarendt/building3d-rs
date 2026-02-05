use std::collections::{HashMap, HashSet};

use crate::geom::bboxes::{are_bboxes_overlapping, bounding_box};
use crate::{Point, Polygon};

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
}
