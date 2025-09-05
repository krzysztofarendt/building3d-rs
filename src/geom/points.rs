use std::fmt;

#[derive(Debug)]
pub struct Point {
    x: f64,
    y: f64,
    z: f64,
}

impl Point {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
}

impl fmt::Display for Point {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let prec = f.precision().unwrap_or(2);  // Default 2 decimals
        write!(f, "Point({:.prec$}, {:.prec$}, {:.prec$})", self.x, self.y, self.z, prec=prec)
    }
}
