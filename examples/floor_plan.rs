use anyhow::Result;
use building3d::FloorPlan;
use building3d::Solid;
use building3d::Vector;
use building3d::draw::simple::draw_polygons;

fn main() -> Result<()> {
    let plan: Vec<(f64, f64)> = vec![(0., 0.), (5., 0.), (5., 2.), (3., 2.), (3., 7.), (0., 5.)];
    let height: f64 = 1.8;
    let rot_angle: f64 = 1.0;
    let translation: Vector = Vector::new(1., 1., 1.);
    let name = Some("building".to_string());

    let fp = FloorPlan {
        plan,
        height,
        name,
        ..Default::default()
    };
    let mut sld = Solid::from_floor_plan(fp)?;

    sld.rotate(rot_angle, &Vector::new(0., 0., 1.));
    sld.translate(&translation);

    draw_polygons(&sld.polygons())?;

    Ok(())
}
