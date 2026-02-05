use crate::Point;
use crate::Polygon;
use crate::TriangleIndex;
use crate::UID;
use crate::Vector;
use crate::Wall;
use crate::geom::IsClose;
use crate::geom::point::check::is_point_in_sequence;
use crate::geom::tetrahedron::tetrahedron_volume;
use crate::{HasMesh, Mesh};
use crate::{HasName, SortByName};
use anyhow::{Context, Result, anyhow};
use std::collections::{HashMap, HashSet};

pub mod containment;

#[derive(Debug, Clone)]
pub struct Solid {
    pub name: String,
    pub uid: UID,
    pub parent: Option<UID>,
    walls: HashMap<String, Wall>,
}

impl HasName for Solid {
    fn get_name(&self) -> &str {
        &self.name
    }
}

impl HasMesh for Solid {
    fn copy_mesh(&self) -> Mesh {
        let polygons = self.polygons();
        let mut vertices: Vec<Point> = Vec::new();
        let mut triangles: Vec<TriangleIndex> = Vec::new();
        let mut num_vertices = 0;

        for &poly in polygons.iter() {
            let mesh = poly.copy_mesh();
            vertices.extend(mesh.vertices);
            let mut tri: Vec<TriangleIndex> = poly.copy_mesh().faces.unwrap();
            tri = tri
                .into_iter()
                .map(|t| TriangleIndex(t.0 + num_vertices, t.1 + num_vertices, t.2 + num_vertices))
                .collect();
            triangles.extend(tri.into_iter());
            num_vertices += poly.mesh_ref().vertices.len();
        }

        Mesh {
            vertices,
            faces: Some(triangles),
        }
    }
}

impl Solid {
    pub fn new(name: &str, mut walls: Vec<Wall>) -> Self {
        let uid = UID::new();
        for w in walls.iter_mut() {
            w.parent = Some(uid.clone());
        }
        let parent = None;
        let walls: HashMap<String, Wall> = walls.into_iter().map(|x| (x.name.clone(), x)).collect();
        Self {
            name: name.to_string(),
            uid,
            parent,
            walls,
        }
    }

    pub fn walls(&self) -> Vec<&Wall> {
        let mut walls: Vec<&Wall> = self.walls.values().collect();
        walls.as_mut_slice().sort_by_name();

        walls
    }

    pub fn polygons(&self) -> Vec<&Polygon> {
        let walls = self.walls();
        let polygons: Vec<&Polygon> = walls.iter().flat_map(|w| w.polygons()).collect();
        // NOTE: Polygons are already sorted

        polygons
    }

    pub fn add_wall(&mut self, mut wall: Wall) -> Result<()> {
        if self.walls.contains_key(&wall.name) {
            return Err(anyhow!("Wall is already present: {}", &wall.name));
        }
        wall.parent = Some(self.uid.clone());
        self.walls.insert(wall.name.clone(), wall);

        Ok(())
    }

    /// Rotates the solid based on the rotation vector and angle (in place)
    pub fn rotate(&mut self, angle: f64, rot_vec: &Vector) {
        let mut walls: Vec<&mut Wall> = self.walls.values_mut().collect();
        for wall in walls.iter_mut() {
            wall.rotate(angle, rot_vec);
        }
    }

    /// Moves the solid by a vector (in place)
    pub fn translate(&mut self, vec: &Vector) {
        let mut walls: Vec<&mut Wall> = self.walls.values_mut().collect();
        for wall in walls.iter_mut() {
            wall.translate(vec);
        }
    }

    /// Calculates the volume of the solid.
    ///
    /// Based on: http://chenlab.ece.cornell.edu/Publication/Cha/icip01_Cha.pdf
    pub fn volume(&self) -> f64 {
        let mut total_volume = 0.;
        let polygons = self.polygons();
        let p0 = Point::new(0., 0., 0.);
        for poly in polygons.iter() {
            for tri in poly.mesh_ref().faces.as_ref().unwrap().iter() {
                let p1 = poly.mesh_ref().vertices[tri.0];
                let p2 = poly.mesh_ref().vertices[tri.1];
                let p3 = poly.mesh_ref().vertices[tri.2];
                let v = tetrahedron_volume(p0, p1, p2, p3);

                let mut pos_wrt_origin = poly.vn.dot(&(p1 - p0));
                if pos_wrt_origin.is_close(0.) {
                    pos_wrt_origin = poly.vn.dot(&(p2 - p0));
                }

                let mut sign = -1.;
                if pos_wrt_origin > 0. {
                    sign = 1.;
                }
                total_volume += sign * v;
            }
        }
        total_volume
    }

    /// Checks if a point lies inside the solid.
    ///
    /// Uses the ray casting algorithm: cast a ray from the point and count
    /// surface crossings. An odd count means the point is inside.
    pub fn is_point_inside(&self, ptest: Point) -> bool {
        containment::is_point_inside_solid(self, ptest)
    }

    /// Checks if a point lies on the boundary (surface) of the solid.
    ///
    /// A point is on the boundary if it lies on any of the solid's polygon surfaces.
    pub fn is_point_at_boundary(&self, ptest: Point) -> bool {
        containment::is_point_at_boundary(self, ptest)
    }

    /// Checks if a point lies strictly inside the solid (not on boundary).
    pub fn is_point_strictly_inside(&self, ptest: Point) -> bool {
        containment::is_point_strictly_inside(self, ptest)
    }

    /// Return a solid with given dimensions and location.
    ///
    /// `x` is the dimension along the X axis.
    /// `y` is the dimension along the Y axis.
    /// `z` is the dimension along the Z axis.
    ///
    /// Rotation is currently not supported.
    ///
    /// The corner `(min(x), min(y), min(z))` will be located at `origin`.
    ///
    /// The polygon and wall names are hardcoded:
    /// - floor/floor
    /// - wall_0/poly_0 (XZ at ymin)
    /// - wall_1/poly_1 (YZ at xmax)
    /// - wall_2/poly_2 (XZ at ymax)
    /// - wall_3/poly_3 (YZ at xmin)
    /// - ceiling/ceiling
    ///
    /// The solid will be named `name` (random if not given).
    pub fn from_box(x: f64, y: f64, z: f64, origin: Option<(f64, f64, f64)>, name: &str) -> Self {
        let origin_vec = match origin {
            Some((dx, dy, dz)) => Vector::new(dx, dy, dz),
            None => Vector::new(0., 0., 0.),
        };

        let p0 = Point::new(0., 0., 0.) + origin_vec;
        let p1 = Point::new(x, 0., 0.) + origin_vec;
        let p2 = Point::new(x, y, 0.) + origin_vec;
        let p3 = Point::new(0., y, 0.) + origin_vec;
        let p4 = Point::new(0., 0., z) + origin_vec;
        let p5 = Point::new(x, 0., z) + origin_vec;
        let p6 = Point::new(x, y, z) + origin_vec;
        let p7 = Point::new(0., y, z) + origin_vec;

        let poly_fl = Polygon::new("floor", vec![p0, p3, p2, p1], None).unwrap();
        let poly_w0 = Polygon::new("poly_0", vec![p0, p1, p5, p4], None).unwrap();
        let poly_w1 = Polygon::new("poly_1", vec![p1, p2, p6, p5], None).unwrap();
        let poly_w2 = Polygon::new("poly_2", vec![p3, p7, p6, p2], None).unwrap();
        let poly_w3 = Polygon::new("poly_3", vec![p0, p4, p7, p3], None).unwrap();
        let poly_rf = Polygon::new("ceiling", vec![p4, p5, p6, p7], None).unwrap();

        let wall_fl = Wall::new("floor", vec![poly_fl]);
        let wall_0 = Wall::new("wall_0", vec![poly_w0]);
        let wall_1 = Wall::new("wall_1", vec![poly_w1]);
        let wall_2 = Wall::new("wall_2", vec![poly_w2]);
        let wall_3 = Wall::new("wall_3", vec![poly_w3]);
        let wall_rf = Wall::new("ceiling", vec![poly_rf]);

        let mut walls: HashMap<String, Wall> = HashMap::new();
        for w in [wall_fl, wall_0, wall_1, wall_2, wall_3, wall_rf] {
            walls.insert(w.name.clone(), w);
        }

        let name: String = name.to_string();
        let uid = UID::new();
        let parent = None;

        Self {
            name,
            uid,
            parent,
            walls,
        }
    }

    pub fn from_floor_plan(fp: FloorPlan) -> Result<Self> {
        // Sanity checks
        if fp.wall_names.is_some() {
            let num_unique_names = fp.wall_names.iter().collect::<HashSet<_>>().len();
            let num_names = fp.wall_names.as_ref().unwrap().len();
            let num_points = fp.plan.len();
            assert_eq!(num_points, num_names);
            assert_eq!(num_names, num_unique_names);
        }

        // Prepare wall names
        let wall_names: Vec<String> = match fp.wall_names {
            Some(names) => names,
            None => {
                let mut wall_names: Vec<String> = Vec::new();
                for i in 0..fp.plan.len() {
                    wall_names.push(format!("wall-{i}"));
                }
                wall_names
            }
        };
        let floor_name: String = match fp.floor_name {
            Some(name) => name,
            None => String::from("floor"),
        };
        let ceil_name: String = match fp.ceiling_name {
            Some(name) => name,
            None => String::from("ceiling"),
        };

        // Set up floor and ceiling points
        let mut floor_pts: Vec<Point> = Vec::new();
        let mut ceil_pts: Vec<Point> = Vec::new();
        for (x, y) in fp.plan.iter() {
            floor_pts.push(Point::new(*x, *y, 0.));
            ceil_pts.push(Point::new(*x, *y, fp.height));
        }

        // Floor height
        let z0: f64 = 0.;

        // Make polygons and walls
        let mut walls: Vec<Wall> = Vec::new();
        for (i, w_name) in wall_names.iter().enumerate() {
            let ths = i; // This point
            let mut nxt = ths + 1; // Next point
            if nxt >= fp.plan.len() {
                nxt = 0;
            }

            let p0 = floor_pts[ths];
            let p1 = floor_pts[nxt];
            let p2 = ceil_pts[nxt];
            let p3 = ceil_pts[ths];

            let poly = Polygon::new(w_name, vec![p0, p1, p2, p3], None).context(format!(
                "Failed to create Polygon from {p0}, {p1}, {p2}, {p3}"
            ))?;
            walls.push(Wall::new(w_name, vec![poly]));
        }

        let mut floor_poly = Polygon::new(&floor_name, floor_pts, None)?;
        // Floor's normal should point downwards
        if !floor_poly.vn.is_close(&Vector::new(0., 0., -1.)) {
            floor_poly = floor_poly.flip(&floor_poly.name)?;
        }

        let mut ceil_poly = Polygon::new(&ceil_name, ceil_pts, None)?;
        // Ceiling normal should point upwards
        if !ceil_poly.vn.is_close(&Vector::new(0., 0., 1.)) {
            ceil_poly = ceil_poly.flip(&ceil_poly.name)?;
        }

        let floor = Wall::new(&floor_name, vec![floor_poly.clone()]);
        let ceil = Wall::new(&ceil_name, vec![ceil_poly]);

        // Make sure all polygon normals point outwards the zone.
        // Compare the order of wall bottom vertices to the order
        // of floor vertices - they should be opposite.
        let f_pts: &Vec<Point> = &floor_poly.mesh_ref().vertices;
        let mut to_flip: Vec<(String, String)> = Vec::new();

        for w in walls.iter() {
            let w_poly: &Polygon = w.polygons()[0]; // There is only 1 polygon
            let w_pts: &Vec<Point> = &w_poly.mesh_ref().vertices;

            // Wall bottom vertices
            let mut wall_pts: Vec<Point> = Vec::new();
            for pt in w_pts.iter() {
                if pt.z.is_close(z0) {
                    wall_pts.push(*pt);
                }
            }

            let mut floor_adjacent_pts: Vec<Point> = Vec::new();

            let mut prev_taken: Option<usize> = None;
            for (i, fpt) in f_pts.iter().enumerate() {
                if is_point_in_sequence(fpt, &wall_pts) {
                    floor_adjacent_pts.push(*fpt);
                    if prev_taken.is_none() {
                        prev_taken = Some(i);
                    } else if prev_taken.is_some_and(|x| x != i - 1) {
                        // Need to flip the vector, because the first and last points
                        // were taken and they are now in the wrong order
                        floor_adjacent_pts.reverse();
                    }
                }

                if floor_adjacent_pts.len() == 2 {
                    break;
                }
            }

            if !wall_pts[0].is_close(&floor_adjacent_pts[1])
                || !wall_pts[1].is_close(&floor_adjacent_pts[0])
            {
                // Wrong direction. Need to reverse the order of polygon vertices
                to_flip.push((w.name.clone(), w_poly.name.clone()));
            }
        }

        // Apply flipping
        for (w_name, p_name) in to_flip.iter() {
            for w in walls.iter_mut() {
                if w.name == *w_name {
                    let w_poly = w.polygons()[0];
                    let flipped_poly = w_poly.flip(p_name)?;
                    w.replace_polygon(p_name, vec![flipped_poly])?;
                }
            }
        }

        // Add floor and ceiling
        walls.push(floor);
        walls.push(ceil);

        // Make solid
        let solid_name = fp.name;
        let solid = Solid::new(&solid_name, walls);

        Ok(solid)
    }
}

pub struct FloorPlan {
    pub plan: Vec<(f64, f64)>,
    pub height: f64,
    pub name: String,
    pub wall_names: Option<Vec<String>>,
    pub floor_name: Option<String>,
    pub ceiling_name: Option<String>,
}

impl Default for FloorPlan {
    fn default() -> Self {
        let plan = vec![(0., 0.), (1., 0.), (1., 1.), (0., 1.)];
        let height = 2.;
        let name = "default_solid".to_string();
        let wall_names = None;
        let floor_name = None;
        let ceiling_name = None;

        Self {
            plan,
            height,
            name,
            wall_names,
            floor_name,
            ceiling_name,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_box() {
        let sld = Solid::from_box(1., 2., 3., None, "box");
        let expected_vol = 1. * 2. * 3.;
        assert!((sld.volume() - expected_vol).abs() < 1e-4);
    }
}
