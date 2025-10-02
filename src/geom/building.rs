use crate::Polygon;
use crate::Solid;
use crate::Vector;
use crate::Wall;
use crate::random_id;
use std::collections::HashMap;

#[derive(Debug, Clone)]
pub struct Building {
    pub name: String,
    pub uid: String,
    pub parent: Option<String>,
    solids: HashMap<String, Solid>,
}

impl Building {
    pub fn new(name: String, mut solids: Vec<Solid>) -> Self {
        let uid = random_id();
        for s in solids.iter_mut() {
            s.parent = Some(uid.clone());
        }
        let parent = None;
        let solids: HashMap<String, Solid> =
            solids.into_iter().map(|x| (x.name.clone(), x)).collect();

        Self {
            name,
            uid,
            parent,
            solids,
        }
    }

    pub fn solids(&self) -> Vec<&Solid> {
        self.solids.values().collect()
    }

    pub fn walls(&self) -> Vec<&Wall> {
        self.solids.values().flat_map(|s| s.walls()).collect()
    }

    pub fn polygons(&self) -> Vec<&Polygon> {
        self.solids
            .values()
            .flat_map(|s| s.walls())
            .flat_map(|w| w.polygons())
            .collect()
    }

    pub fn rotate(&mut self, angle: f64, rot_vec: &Vector) {
        let mut solids: Vec<&mut Solid> = self.solids.values_mut().collect();
        for sld in solids.iter_mut() {
            sld.rotate(angle, rot_vec);
        }
    }

    pub fn translate(&mut self, vec: &Vector) {
        let mut solids: Vec<&mut Solid> = self.solids.values_mut().collect();
        for sld in solids.iter_mut() {
            sld.translate(vec);
        }
    }

    pub fn volume(&self) -> f64 {
        let mut v = 0.;
        for sld in self.solids() {
            v += sld.volume();
        }
        v
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_volume() {
        let s0 = Solid::from_box(1., 1., 1., None, None);
        let s1 = Solid::from_box(1., 2., 3., None, None);
        let bdg = Building::new("building".to_string(), vec![s0, s1]);
        let expected_vol = 1. * 2. * 3. + 1.;
        assert!((bdg.volume() - expected_vol).abs() < 1e-4);
    }
}
