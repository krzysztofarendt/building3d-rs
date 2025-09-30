use building3d::Solid;
use building3d::FloorPlan;
use building3d::draw_polygons;
use anyhow::Result;


fn main() -> Result<()> {

    let plan: Vec<(f64, f64)> = vec![
        (0., 0.), (5., 0.), (5., 2.), (3., 2.), (3., 7.), (0., 5.)
    ];
    let height: f64 = 1.8;
    let name = Some("building".to_string());

    let fp = FloorPlan {
        plan,
        height,
        name,
        ..Default::default()
    };
    let sld = Solid::from_floor_plan(fp)?;
    draw_polygons(&sld.polygons())?;

    Ok(())
}
