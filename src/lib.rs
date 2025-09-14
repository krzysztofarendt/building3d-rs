pub mod geom;
mod id;
pub mod vecutils;
mod draw;

// Prelude
pub use geom::point::Point;
pub use geom::polygon::Polygon;
pub use geom::vector::Vector;
use id::random_id;
// Drawing utility
pub use draw::draw_polygon;
