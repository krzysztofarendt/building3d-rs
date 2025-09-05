use crate::Point;
// use crate::Vector;

#[derive(Debug, Clone, PartialEq)]
pub struct Polygon {
    pts: Vec<Point>,
    name: String,
    uid: String,
    tri: Vec<usize>,
    vn: Vec<Point>,
}
