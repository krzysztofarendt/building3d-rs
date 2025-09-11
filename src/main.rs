use building3d::{Point, Polygon};
use anyhow::Result;

fn main() -> Result<()> {
    let p0 = Point::new(0., 0., 0.);
    let p1 = Point::new(1., 0., 0.);
    let p2 = Point::new(1., 1., 0.);
    let p3 = Point::new(0., 1., 0.);

    let poly = Polygon::new("polygon".to_string(), vec![p0, p1, p2, p3], None)?;

    println!("{:.1}", p0);
    println!("{:.1}", poly);

    Ok(())
}
