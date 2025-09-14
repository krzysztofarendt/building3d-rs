use building3d::{Point, Polygon, draw_polygon};
use anyhow::Result;

fn main() -> Result<()> {
    let p0 = Point::new(0., 0., 0.);
    let p1 = Point::new(1., 0., 0.);
    let p2 = Point::new(1., 1., 0.);
    let p3 = Point::new(2., 1., 0.);
    let p4 = Point::new(2., 0., 0.);
    let p5 = Point::new(3., 0., 0.);
    let p6 = Point::new(3., 2., 0.);
    let p7 = Point::new(0., 2., 0.);
    let pts = vec![p0, p1, p2, p3, p4, p5, p6, p7];
    let poly = Polygon::new("polygon".to_string(), pts, None)?;

    println!("{:.1}", p0);
    println!("{:.1}", poly);

    // Draw the polygon in a 3D window
    draw_polygon(&poly)?;
    Ok(())
}
