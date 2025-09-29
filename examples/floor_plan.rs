use building3d::Solid;
use building3d::FloorPlan;
use building3d::draw_polygons;
use anyhow::Result;


fn main() -> Result<()> {

    let fp = FloorPlan {
        ..Default::default()
    };
    let sld = Solid::from_floor_plan(fp)?;
    draw_polygons(&sld.polygons())?;

    Ok(())
}
