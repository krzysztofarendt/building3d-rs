use anyhow::Result;
use building3d::{Solid, draw_polygons};

fn main() -> Result<()> {
    let sld = Solid::from_box(1., 1., 1., None, None);
    draw_polygons(&sld.polygons())?;

    Ok(())
}
