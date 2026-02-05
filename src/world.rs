//! World module - placeholder for future global registry functionality.
//!
//! The Building struct already provides path-based access to zones, solids,
//! walls, and polygons through methods like:
//! - `building.get_zone(path)`
//! - `building.get_solid(path)`
//! - `building.get_wall(path)`
//! - `building.get_polygon(path)`
//!
//! A World struct could be added in the future if needed for:
//! - Global registry of multiple buildings
//! - Cross-building queries
//! - UID-based lookup across all objects
