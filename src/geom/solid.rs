use crate::geom::point::check::is_point_in_sequence;
use crate::geom::polygon::Polygon;
use crate::geom::rotation::rotate_points_around_vector;
use crate::geom::vector::Vector;
use crate::geom::wall::Wall;
use crate::geom::{IsClose, point::Point};
use crate::random_id;
use anyhow::{Result, anyhow};
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone)]
pub struct Solid {
    pub name: String,
    pub uid: String,
    pub parent: Option<String>,
    walls: HashMap<String, Wall>,
}

impl Solid {
    pub fn new(name: String, mut walls: Vec<Wall>) -> Self {
        let uid = random_id();
        for w in walls.iter_mut() {
            w.parent = Some(uid.clone());
        }
        let parent = None;
        let walls: HashMap<String, Wall> = walls.into_iter().map(|x| (x.name.clone(), x)).collect();
        Self {
            name,
            uid,
            parent,
            walls,
        }
    }

    pub fn walls(&self) -> Vec<&Wall> {
        self.walls.values().collect()
    }

    pub fn polygons(&self) -> Vec<&Polygon> {
        self.walls.values().flat_map(|w| w.polygons()).collect()
    }

    pub fn add_wall(&mut self, mut wall: Wall) -> Result<()> {
        if self.walls.contains_key(&wall.name) {
            return Err(anyhow!("Wall is already present: {}", &wall.name));
        }
        wall.parent = Some(self.uid.clone());
        self.walls.insert(wall.name.clone(), wall);

        Ok(())
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
    pub fn from_box(
        x: f64,
        y: f64,
        z: f64,
        origin: Option<(f64, f64, f64)>,
        name: Option<&str>,
    ) -> Self {
        // TODO: Add rotation

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

        let poly_fl = Polygon::new("floor".to_string(), vec![p0, p3, p2, p1], None).unwrap();
        let poly_w0 = Polygon::new("poly_0".to_string(), vec![p0, p1, p5, p4], None).unwrap();
        let poly_w1 = Polygon::new("poly_1".to_string(), vec![p1, p2, p6, p5], None).unwrap();
        let poly_w2 = Polygon::new("poly_2".to_string(), vec![p3, p7, p6, p2], None).unwrap();
        let poly_w3 = Polygon::new("poly_3".to_string(), vec![p0, p4, p7, p3], None).unwrap();
        let poly_rf = Polygon::new("ceiling".to_string(), vec![p4, p5, p6, p7], None).unwrap();

        let wall_fl = Wall::new("floor".to_string(), vec![poly_fl]);
        let wall_0 = Wall::new("wall_0".to_string(), vec![poly_w0]);
        let wall_1 = Wall::new("wall_1".to_string(), vec![poly_w1]);
        let wall_2 = Wall::new("wall_2".to_string(), vec![poly_w2]);
        let wall_3 = Wall::new("wall_3".to_string(), vec![poly_w3]);
        let wall_rf = Wall::new("ceiling".to_string(), vec![poly_rf]);

        let mut walls: HashMap<String, Wall> = HashMap::new();
        for w in [wall_fl, wall_0, wall_1, wall_2, wall_3, wall_rf] {
            walls.insert(w.name.clone(), w);
        }

        let name: String = match name {
            Some(n) => n.to_string(),
            None => random_id(),
        };

        let uid = random_id();
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

        // Define the rotation vector (it is hardcoded, floor and ceiling must be horizontal)
        let rot_vec: Vector = Vector::new(0., 0., 1.);

        // Prepare wall names
        let wall_names: Vec<String> = match fp.wall_names {
            Some(names) => {
                let mut wall_names: Vec<String> = Vec::new();
                for i in 0..names.len() {
                    wall_names.push(format!("wall-{i}"));
                }
                wall_names
            }
            None => fp.wall_names.as_ref().unwrap().clone(),
        };
        let floor_name: String = match fp.floor_name {
            Some(name) => name,
            None => String::from("floor"),
        };
        let ceil_name: String = match fp.ceiling_name {
            Some(name) => name,
            None => String::from("floor"),
        };

        // Set up floor and ceiling points
        let mut floor_pts: Vec<Point> = Vec::new();
        let mut ceil_pts: Vec<Point> = Vec::new();
        for (x, y) in fp.plan.iter() {
            floor_pts.push(Point::new(*x, *y, 0.));
            ceil_pts.push(Point::new(*x, *y, fp.height));
        }

        // Rotate
        if !fp.rot_angle.is_close(0.) {
            floor_pts = rotate_points_around_vector(&floor_pts, &rot_vec, fp.rot_angle);
            ceil_pts = rotate_points_around_vector(&ceil_pts, &rot_vec, fp.rot_angle);
        }

        // Translate
        let mut z0: f64 = 0.;

        if !fp.translate.length().is_close(0.) {
            for pt in floor_pts.iter_mut() {
                *pt = *pt + fp.translate;
            }
            for pt in ceil_pts.iter_mut() {
                *pt = *pt + fp.translate;
            }
            z0 = fp.translate.dz;
        }

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
            let p2 = ceil_pts[ths];
            let p3 = ceil_pts[nxt];

            let poly = Polygon::new(w_name.clone(), vec![p0, p1, p2, p3], None)?;
            walls.push(Wall::new(w_name.clone(), vec![poly]));
        }

        let mut floor_poly = Polygon::new(floor_name.clone(), floor_pts, None)?;
        // Floor's normal should point downwards
        if !floor_poly.vn.is_close(&Vector::new(0., 0., -1.)) {
            floor_poly = floor_poly.flip(floor_poly.name.clone())?;
        }

        let mut ceil_poly = Polygon::new(ceil_name.clone(), ceil_pts, None)?;
        // Ceiling normal should point upwards
        if !ceil_poly.vn.is_close(&Vector::new(0., 0., 1.)) {
            ceil_poly = ceil_poly.flip(ceil_poly.name.clone())?;
        }

        let floor = Wall::new(floor_name, vec![floor_poly.clone()]);
        let ceil = Wall::new(ceil_name, vec![ceil_poly]);

        // Make sure all polygon normals point outwards the zone.
        // Compare the order of wall bottom vertices to the order
        // of floor vertices - they should be opposite.
        let f_pts: &Vec<Point> = &floor_poly.pts;
        let mut to_flip: Vec<(String, String)> = Vec::new();

        for w in walls.iter() {
            let w_poly: &Polygon = w.polygons()[0]; // There is only 1 polygon
            let w_pts: &Vec<Point> = &w_poly.pts;

            // Wall bottom vertices
            let mut wall_z0_pts: Vec<Point> = Vec::new();
            for pt in w_pts.iter() {
                if pt.z.is_close(z0) {
                    wall_z0_pts.push(*pt);
                }
            }

            let mut floor_adjacent_pts: Vec<Point> = Vec::new();

            let mut prev_taken: Option<usize> = None;
            for (i, fpt) in f_pts.iter().enumerate() {
                if is_point_in_sequence(fpt, &wall_z0_pts) {
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

            if !wall_z0_pts[0].is_close(&floor_adjacent_pts[1])
                || !wall_z0_pts[1].is_close(&floor_adjacent_pts[0])
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
                    let flipped_poly = w_poly.flip(p_name.clone())?;
                    w.replace_polygon(p_name, vec![flipped_poly])?;
                }
            }
        }

        // Add floor and ceiling
        walls.push(floor);
        walls.push(ceil);

        // Make solid
        let solid_name = fp.name.unwrap_or(random_id());
        let solid = Solid::new(solid_name, walls);

        Ok(solid)
    }
}

pub struct FloorPlan {
    pub plan: Vec<(f64, f64)>,
    pub height: f64,
    pub translate: Vector, // TODO: this should be a method on a solid
    pub rot_angle: f64,    // TODO: this should be a method on a solid
    pub name: Option<String>,
    pub wall_names: Option<Vec<String>>,
    pub floor_name: Option<String>,
    pub ceiling_name: Option<String>,
}

impl Default for FloorPlan {
    fn default() -> Self {
        let plan = vec![(0., 0.), (1., 0.), (1., 1.), (0., 1.)];
        let height = 2.;
        let translate = Vector::new(0., 0., 0.);
        let rot_angle = 0.;
        let name = None;
        let wall_names = None;
        let floor_name = None;
        let ceiling_name = None;

        Self {
            plan,
            height,
            translate,
            rot_angle,
            name,
            wall_names,
            floor_name,
            ceiling_name,
        }
    }
}
