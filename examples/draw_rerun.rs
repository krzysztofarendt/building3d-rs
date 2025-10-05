use anyhow::Result;
use building3d::Building;
use building3d::FloorPlan;
use building3d::Solid;
use building3d::draw::rerun::{draw_edges, draw_surfaces, start_session};

fn main() -> Result<()> {
    // Make model
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
    let building = Building::new("building".to_string(), vec![sld1.clone(), sld2.clone()]);

    // Draw
    let session = start_session()?;

    let rgba: (f32, f32, f32, f32) = (1., 1., 1., 0.2);
    draw_surfaces(&session, &building, rgba)?;

    let radius: f32 = 0.01;
    let rgba: (f32, f32, f32, f32) = (0., 0., 1., 0.5);
    draw_edges(&session, &sld1, radius, rgba)?;

    let radius: f32 = 0.01;
    let rgba: (f32, f32, f32, f32) = (1., 0., 0., 0.5);
    draw_edges(&session, &sld2, radius, rgba)?;

    Ok(())
}
