use building3d::FloorPlan;
use building3d::Solid;
// use building3d::Vector;
use anyhow::Result;
use building3d::Building;
use building3d::draw::rerun::show;

fn main() -> Result<()> {
    let plan: Vec<(f64, f64)> = vec![(0., 0.), (5., 0.), (5., 2.), (3., 2.), (3., 7.), (0., 5.)];
    let height: f64 = 1.8;
    let name = Some("building".to_string());

    let fp = FloorPlan {
        plan,
        height,
        name,
        ..Default::default()
    };
    let sld1 = Solid::from_floor_plan(fp)?;

    let sld2 = Solid::from_box(5., 5., height, Some((-5., 0., 0.)), Some("box"));

    let bdg = Building::new("building".to_string(), vec![sld1, sld2]);

    show(&bdg)?;

    Ok(())
}
