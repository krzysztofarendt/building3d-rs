pub mod geom;
mod id;
pub mod vecutils;
mod draw;

// Prelude
pub use geom::point::Point;
pub use geom::polygon::Polygon;
pub use geom::wall::Wall;
pub use geom::vector::Vector;
pub use geom::solid::Solid;
pub use geom::solid::FloorPlan;
use id::random_id;
// Drawing utility
pub use draw::draw_polygons;
