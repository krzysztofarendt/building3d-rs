use crate::Point;
use ndarray as nd;

pub fn points_to_array(points: &[Point]) -> nd::Array2<f64> {
    let mut arr = nd::Array2::from_elem((points.len(), 3), 0.);

    for (i, p) in points.iter().enumerate() {
        arr[[i, 0]] = p.x;
        arr[[i, 1]] = p.y;
        arr[[i, 2]] = p.z;
    }

    arr
}

pub fn array_to_points(arr: nd::Array2<f64>) -> Vec<Point> {
    let mut pts: Vec<Point> = Vec::new();
    let num_pts = arr.shape()[0];

    for i in 0..num_pts {
        let (x, y, z) = (arr[[i, 0]], arr[[i, 1]], arr[[i, 2]]);
        pts.push(Point::new(x, y, z));
    }

    pts
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_conversion() {
        let a0 = Point::new(1., 2., 3.);
        let a1 = Point::new(4., 5., 6.);
        let a2 = Point::new(7., 8., 9.);
        let arr = points_to_array(&[a0, a1, a2]);
        let v = array_to_points(arr);
        assert!(v[0] == a0);
        assert!(v[1] == a1);
        assert!(v[2] == a2);
    }
}
