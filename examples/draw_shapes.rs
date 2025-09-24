use anyhow::Result;
use building3d::{Point, Polygon, Wall, draw_polygons};

fn main() -> Result<()> {
    // Draw multiple polygons in a 3D window
    // U-shaped polygon
    let p0 = Point::new(0., 0., 0.);
    let p1 = Point::new(1., 0., 0.);
    let p2 = Point::new(1., 1., 0.);
    let p3 = Point::new(2., 1., 0.);
    let p4 = Point::new(2., 0., 0.);
    let p5 = Point::new(3., 0., 0.);
    let p6 = Point::new(3., 2., 0.);
    let p7 = Point::new(0., 2., 0.);
    let pts = vec![p0, p1, p2, p3, p4, p5, p6, p7];
    let u_shape = Polygon::new("u-shape".to_string(), pts, None)?;

    // A triangle offset along X
    let tri = Polygon::new(
        "triangle".to_string(),
        vec![
            Point::new(2.0, 0.0, 0.0),
            Point::new(3.0, 0.0, 0.0),
            Point::new(2.5, 1.0, 0.0),
        ],
        None,
    )?;
    // A right-angled L-shape polygon lifted in Z
    let l_shape = Polygon::new(
        "l-shape".to_string(),
        vec![
            Point::new(0.0, 2.0, 0.5),
            Point::new(1.0, 2.0, 0.5),
            Point::new(1.0, 3.0, 0.5),
            Point::new(2.0, 3.0, 0.5),
            Point::new(2.0, 4.0, 0.5),
            Point::new(0.0, 4.0, 0.5),
        ],
        None,
    )?;

    let wall = Wall::new("wall".to_string(), vec![u_shape, tri, l_shape]);
    println!("{wall:?}");

    draw_polygons(&wall.polygons())?;

    Ok(())
}
