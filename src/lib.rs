//! # building3d
//!
//! A 3D building modeling library with a strict hierarchical composition model:
//!
//! **Building → Zone → Solid → Wall → Polygon → Mesh**
//!
//! Each entity carries a name, a unique identifier ([`UID`]), and an optional
//! parent reference. Geometry is built bottom-up from [`Point`]s and
//! [`Polygon`]s, while containers ([`Wall`], [`Solid`], [`Zone`], [`Building`])
//! aggregate children via `HashMap<String, Child>` for O(1) access by name.
//!
//! ## Quick start
//!
//! ```
//! use building3d::{Building, Solid};
//!
//! let s = Solid::from_box(3.0, 4.0, 5.0, None, "room").unwrap();
//! let b = Building::from_solids("house", vec![s]).unwrap();
//! assert!((b.volume() - 60.0).abs() < 1e-10);
//! ```

pub mod draw;
pub mod geom;
pub mod io;
mod name;
pub mod sim;
mod uid;
pub mod vecutils;
pub mod world;

// Prelude
pub use draw::config::RerunConfig;
pub use geom::building::Building;
pub use geom::mesh::{HasMesh, Mesh};
pub use geom::point::Point;
pub use geom::polygon::Polygon;
pub use geom::solid::FloorPlan;
pub use geom::solid::Solid;
pub use geom::triangles::TriangleIndex;
pub use geom::vector::Vector;
pub use geom::wall::Wall;
pub use geom::zone::Zone;
pub use name::{HasName, SortByName};
pub use uid::UID;
