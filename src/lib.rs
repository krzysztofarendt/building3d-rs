pub mod draw;
pub mod geom;
mod uid;
mod name;
pub mod vecutils;
pub mod world;

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
pub use uid::UID;
pub use name::{HasName, SortByName};
// pub use world::World; // TODO
