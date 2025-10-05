pub mod draw;
pub mod geom;
mod name;
mod id;
pub mod vecutils;

// Prelude
pub use geom::building::Building;
pub use geom::mesh::{HasMesh, Mesh};
pub use geom::point::Point;
pub use geom::polygon::Polygon;
pub use geom::solid::FloorPlan;
pub use geom::solid::Solid;
pub use geom::triangles::TriangleIndex;
pub use geom::vector::Vector;
pub use geom::wall::Wall;
pub use name::{HasName, SortByName};
use id::random_id;
