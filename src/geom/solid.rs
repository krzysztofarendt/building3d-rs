use crate::geom::point::Point;
use crate::geom::polygon::Polygon;
use crate::geom::vector::Vector;
use crate::geom::wall::Wall;
use crate::random_id;
use anyhow::{Result, anyhow};
use std::collections::HashMap;

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
}
