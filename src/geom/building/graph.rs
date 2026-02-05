//! Building graph construction and analysis.
//!
//! This module provides functions for building adjacency graphs at different
//! levels of the hierarchy (polygon, wall, solid, zone) and for various
//! relationship types (facing, overlapping, touching).

use std::collections::HashMap;

use crate::Building;
use crate::Polygon;
use crate::Solid;
use crate::Wall;
use crate::geom::polygon::relations::{are_polygons_facing, are_polygons_touching};

/// Represents the level of granularity for graph construction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GraphLevel {
    /// Graph edges connect polygons
    Polygon,
    /// Graph edges connect walls
    Wall,
    /// Graph edges connect solids
    Solid,
    /// Graph edges connect zones
    Zone,
}

/// Parameters for graph construction.
#[derive(Debug, Clone, Copy)]
pub struct GraphParams {
    /// Level of granularity for the graph
    pub level: GraphLevel,
    /// Include edges for facing relationships (opposite normals, coplanar, overlapping)
    pub facing: bool,
    /// Include edges for touching relationships (shared edges, no overlap)
    pub touching: bool,
}

impl Default for GraphParams {
    fn default() -> Self {
        Self {
            level: GraphLevel::Solid,
            facing: true,
            touching: false,
        }
    }
}

/// Represents an edge in the building graph.
#[derive(Debug, Clone)]
pub struct GraphEdge {
    /// Path to the first entity
    pub path1: String,
    /// Path to the second entity
    pub path2: String,
    /// Type of relationship
    pub relationship: String,
}

/// Builds an adjacency graph for the building.
///
/// The graph maps entity paths to lists of adjacent entity paths.
/// The level parameter determines whether entities are polygons, walls, solids, or zones.
///
/// # Arguments
/// * `building` - The building to analyze
/// * `params` - Parameters controlling graph construction
///
/// # Returns
/// A HashMap where keys are entity paths and values are lists of adjacent entity paths.
pub fn get_graph(building: &Building, params: GraphParams) -> HashMap<String, Vec<String>> {
    let mut graph: HashMap<String, Vec<String>> = HashMap::new();

    match params.level {
        GraphLevel::Polygon => build_polygon_graph(building, params, &mut graph),
        GraphLevel::Wall => build_wall_graph(building, params, &mut graph),
        GraphLevel::Solid => build_solid_graph(building, params, &mut graph),
        GraphLevel::Zone => build_zone_graph(building, params, &mut graph),
    }

    graph
}

/// Returns all edges in the building graph.
///
/// Unlike `get_graph()` which returns adjacency lists, this returns a flat list
/// of edges with relationship information.
pub fn get_graph_edges(building: &Building, params: GraphParams) -> Vec<GraphEdge> {
    let mut edges = Vec::new();

    match params.level {
        GraphLevel::Polygon => collect_polygon_edges(building, params, &mut edges),
        GraphLevel::Wall => collect_wall_edges(building, params, &mut edges),
        GraphLevel::Solid => collect_solid_edges(building, params, &mut edges),
        GraphLevel::Zone => collect_zone_edges(building, params, &mut edges),
    }

    edges
}

/// Collects all polygons from a building with their paths.
fn collect_polygons_with_paths(building: &Building) -> Vec<(String, Polygon)> {
    let mut result = Vec::new();

    for zone in building.zones() {
        for solid in zone.solids() {
            for wall in solid.walls() {
                for poly in wall.polygons() {
                    let path = format!("{}/{}/{}/{}", zone.name, solid.name, wall.name, poly.name);
                    result.push((path, poly.clone()));
                }
            }
        }
    }

    result
}

/// Collects all walls from a building with their paths.
fn collect_walls_with_paths(building: &Building) -> Vec<(String, Wall)> {
    let mut result = Vec::new();

    for zone in building.zones() {
        for solid in zone.solids() {
            for wall in solid.walls() {
                let path = format!("{}/{}/{}", zone.name, solid.name, wall.name);
                result.push((path, wall.clone()));
            }
        }
    }

    result
}

/// Collects all solids from a building with their paths.
fn collect_solids_with_paths(building: &Building) -> Vec<(String, Solid)> {
    let mut result = Vec::new();

    for zone in building.zones() {
        for solid in zone.solids() {
            let path = format!("{}/{}", zone.name, solid.name);
            result.push((path, solid.clone()));
        }
    }

    result
}

fn build_polygon_graph(
    building: &Building,
    params: GraphParams,
    graph: &mut HashMap<String, Vec<String>>,
) {
    let polygons = collect_polygons_with_paths(building);

    for (i, (path1, poly1)) in polygons.iter().enumerate() {
        for (path2, poly2) in polygons.iter().skip(i + 1) {
            let mut connected = false;

            if params.facing && are_polygons_facing(poly1, poly2) {
                connected = true;
            }

            if params.touching && are_polygons_touching(poly1, poly2) {
                connected = true;
            }

            if connected {
                graph.entry(path1.clone()).or_default().push(path2.clone());
                graph.entry(path2.clone()).or_default().push(path1.clone());
            }
        }
    }
}

fn build_wall_graph(
    building: &Building,
    params: GraphParams,
    graph: &mut HashMap<String, Vec<String>>,
) {
    let walls = collect_walls_with_paths(building);

    for (i, (path1, wall1)) in walls.iter().enumerate() {
        for (path2, wall2) in walls.iter().skip(i + 1) {
            // Check if any polygons between the walls have the specified relationship
            let mut connected = false;
            'outer: for poly1 in wall1.polygons() {
                for poly2 in wall2.polygons() {
                    if params.facing && are_polygons_facing(poly1, poly2) {
                        connected = true;
                        break 'outer;
                    }
                    if params.touching && are_polygons_touching(poly1, poly2) {
                        connected = true;
                        break 'outer;
                    }
                }
            }

            if connected {
                graph.entry(path1.clone()).or_default().push(path2.clone());
                graph.entry(path2.clone()).or_default().push(path1.clone());
            }
        }
    }
}

fn build_solid_graph(
    building: &Building,
    params: GraphParams,
    graph: &mut HashMap<String, Vec<String>>,
) {
    let solids = collect_solids_with_paths(building);

    for (i, (path1, solid1)) in solids.iter().enumerate() {
        for (path2, solid2) in solids.iter().skip(i + 1) {
            let mut connected = false;
            'outer: for poly1 in solid1.polygons() {
                for poly2 in solid2.polygons() {
                    if params.facing && are_polygons_facing(poly1, poly2) {
                        connected = true;
                        break 'outer;
                    }
                    if params.touching && are_polygons_touching(poly1, poly2) {
                        connected = true;
                        break 'outer;
                    }
                }
            }

            if connected {
                graph.entry(path1.clone()).or_default().push(path2.clone());
                graph.entry(path2.clone()).or_default().push(path1.clone());
            }
        }
    }
}

fn build_zone_graph(
    building: &Building,
    params: GraphParams,
    graph: &mut HashMap<String, Vec<String>>,
) {
    let zones = building.zones();

    for (i, zone1) in zones.iter().enumerate() {
        for zone2 in zones.iter().skip(i + 1) {
            let mut connected = false;
            'outer: for solid1 in zone1.solids() {
                for solid2 in zone2.solids() {
                    for poly1 in solid1.polygons() {
                        for poly2 in solid2.polygons() {
                            if params.facing && are_polygons_facing(poly1, poly2) {
                                connected = true;
                                break 'outer;
                            }
                            if params.touching && are_polygons_touching(poly1, poly2) {
                                connected = true;
                                break 'outer;
                            }
                        }
                    }
                }
            }

            if connected {
                graph
                    .entry(zone1.name.clone())
                    .or_default()
                    .push(zone2.name.clone());
                graph
                    .entry(zone2.name.clone())
                    .or_default()
                    .push(zone1.name.clone());
            }
        }
    }
}

fn collect_polygon_edges(building: &Building, params: GraphParams, edges: &mut Vec<GraphEdge>) {
    let polygons = collect_polygons_with_paths(building);

    for (i, (path1, poly1)) in polygons.iter().enumerate() {
        for (path2, poly2) in polygons.iter().skip(i + 1) {
            if params.facing && are_polygons_facing(poly1, poly2) {
                edges.push(GraphEdge {
                    path1: path1.clone(),
                    path2: path2.clone(),
                    relationship: "facing".to_string(),
                });
            }

            if params.touching && are_polygons_touching(poly1, poly2) {
                edges.push(GraphEdge {
                    path1: path1.clone(),
                    path2: path2.clone(),
                    relationship: "touching".to_string(),
                });
            }
        }
    }
}

fn collect_wall_edges(building: &Building, params: GraphParams, edges: &mut Vec<GraphEdge>) {
    let walls = collect_walls_with_paths(building);

    for (i, (path1, wall1)) in walls.iter().enumerate() {
        for (path2, wall2) in walls.iter().skip(i + 1) {
            let mut found_facing = false;
            let mut found_touching = false;

            for poly1 in wall1.polygons() {
                for poly2 in wall2.polygons() {
                    if params.facing && !found_facing && are_polygons_facing(poly1, poly2) {
                        found_facing = true;
                    }
                    if params.touching && !found_touching && are_polygons_touching(poly1, poly2) {
                        found_touching = true;
                    }
                }
            }

            if found_facing {
                edges.push(GraphEdge {
                    path1: path1.clone(),
                    path2: path2.clone(),
                    relationship: "facing".to_string(),
                });
            }

            if found_touching {
                edges.push(GraphEdge {
                    path1: path1.clone(),
                    path2: path2.clone(),
                    relationship: "touching".to_string(),
                });
            }
        }
    }
}

fn collect_solid_edges(building: &Building, params: GraphParams, edges: &mut Vec<GraphEdge>) {
    let solids = collect_solids_with_paths(building);

    for (i, (path1, solid1)) in solids.iter().enumerate() {
        for (path2, solid2) in solids.iter().skip(i + 1) {
            let mut found_facing = false;
            let mut found_touching = false;

            for poly1 in solid1.polygons() {
                for poly2 in solid2.polygons() {
                    if params.facing && !found_facing && are_polygons_facing(poly1, poly2) {
                        found_facing = true;
                    }
                    if params.touching && !found_touching && are_polygons_touching(poly1, poly2) {
                        found_touching = true;
                    }
                }
            }

            if found_facing {
                edges.push(GraphEdge {
                    path1: path1.clone(),
                    path2: path2.clone(),
                    relationship: "facing".to_string(),
                });
            }

            if found_touching {
                edges.push(GraphEdge {
                    path1: path1.clone(),
                    path2: path2.clone(),
                    relationship: "touching".to_string(),
                });
            }
        }
    }
}

fn collect_zone_edges(building: &Building, params: GraphParams, edges: &mut Vec<GraphEdge>) {
    let zones = building.zones();

    for (i, zone1) in zones.iter().enumerate() {
        for zone2 in zones.iter().skip(i + 1) {
            let mut found_facing = false;
            let mut found_touching = false;

            for solid1 in zone1.solids() {
                for solid2 in zone2.solids() {
                    for poly1 in solid1.polygons() {
                        for poly2 in solid2.polygons() {
                            if params.facing && !found_facing && are_polygons_facing(poly1, poly2) {
                                found_facing = true;
                            }
                            if params.touching
                                && !found_touching
                                && are_polygons_touching(poly1, poly2)
                            {
                                found_touching = true;
                            }
                        }
                    }
                }
            }

            if found_facing {
                edges.push(GraphEdge {
                    path1: zone1.name.clone(),
                    path2: zone2.name.clone(),
                    relationship: "facing".to_string(),
                });
            }

            if found_touching {
                edges.push(GraphEdge {
                    path1: zone1.name.clone(),
                    path2: zone2.name.clone(),
                    relationship: "touching".to_string(),
                });
            }
        }
    }
}

/// Information about a stitched interface between solids.
#[derive(Debug, Clone)]
pub struct StitchInfo {
    /// Path to the first solid (zone/solid)
    pub solid1_path: String,
    /// Path to the second solid (zone/solid)
    pub solid2_path: String,
    /// Name of the polygon in solid1
    pub poly1_name: String,
    /// Name of the polygon in solid2
    pub poly2_name: String,
    /// Whether the interface is correct (matching vertices)
    pub interface_correct: bool,
}

/// Finds all interfaces between adjacent solids in the building.
///
/// Returns detailed information about each interface, including whether
/// the interface is correct (faces match exactly).
pub fn stitch_solids(building: &Building) -> Vec<StitchInfo> {
    let mut stitches = Vec::new();
    let solids = collect_solids_with_paths(building);

    for (i, (path1, solid1)) in solids.iter().enumerate() {
        for (path2, solid2) in solids.iter().skip(i + 1) {
            // Get shared polygons between these solids
            let shared = solid1.get_shared_polygons(solid2);

            for sp in shared {
                let interface_correct = solid1.has_correct_interface(solid2);
                stitches.push(StitchInfo {
                    solid1_path: path1.clone(),
                    solid2_path: path2.clone(),
                    poly1_name: sp.poly1_name,
                    poly2_name: sp.poly2_name,
                    interface_correct,
                });
            }
        }
    }

    stitches
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Zone;

    fn make_adjacent_building() -> Building {
        // Create two adjacent boxes sharing a face
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1");
        let s2 = Solid::from_box(1.0, 1.0, 1.0, Some((1.0, 0.0, 0.0)), "box2");
        Building::from_solids("building", vec![s1, s2])
    }

    fn make_building_with_zones() -> Building {
        // Create building with two zones
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1");
        let z1 = Zone::new("zone1", vec![s1]);

        let s2 = Solid::from_box(1.0, 1.0, 1.0, Some((1.0, 0.0, 0.0)), "box2");
        let z2 = Zone::new("zone2", vec![s2]);

        Building::new("building", vec![z1, z2])
    }

    #[test]
    fn test_solid_graph_adjacent() {
        let building = make_adjacent_building();
        let params = GraphParams {
            level: GraphLevel::Solid,
            facing: true,
            touching: false,
        };

        let graph = get_graph(&building, params);

        // Should have two solids in the graph
        assert_eq!(graph.len(), 2);

        // Each solid should be connected to the other
        for (_path, neighbors) in &graph {
            assert_eq!(neighbors.len(), 1);
        }
    }

    #[test]
    fn test_solid_graph_not_adjacent() {
        // Two solids with a gap
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1");
        let s2 = Solid::from_box(1.0, 1.0, 1.0, Some((3.0, 0.0, 0.0)), "box2");
        let building = Building::from_solids("building", vec![s1, s2]);

        let params = GraphParams::default();
        let graph = get_graph(&building, params);

        // Graph should be empty (no adjacencies)
        assert!(graph.is_empty());
    }

    #[test]
    fn test_zone_graph() {
        let building = make_building_with_zones();
        let params = GraphParams {
            level: GraphLevel::Zone,
            facing: true,
            touching: false,
        };

        let graph = get_graph(&building, params);

        // Two zones should be connected
        assert_eq!(graph.len(), 2);

        // Each zone connected to the other
        assert!(graph.contains_key("zone1"));
        assert!(graph.contains_key("zone2"));
    }

    #[test]
    fn test_graph_edges() {
        let building = make_adjacent_building();
        let params = GraphParams {
            level: GraphLevel::Solid,
            facing: true,
            touching: false,
        };

        let edges = get_graph_edges(&building, params);

        // Should have one edge between the two solids
        assert_eq!(edges.len(), 1);
        assert_eq!(edges[0].relationship, "facing");
    }

    #[test]
    fn test_stitch_solids() {
        let building = make_adjacent_building();
        let stitches = stitch_solids(&building);

        // Should have one stitch (one interface)
        assert_eq!(stitches.len(), 1);
        assert!(stitches[0].interface_correct);
    }

    #[test]
    fn test_stitch_solids_incorrect_interface() {
        // Two solids with different face sizes
        let s1 = Solid::from_box(1.0, 1.0, 1.0, None, "box1");
        let s2 = Solid::from_box(1.0, 2.0, 2.0, Some((1.0, 0.0, 0.0)), "box2");
        let building = Building::from_solids("building", vec![s1, s2]);

        let stitches = stitch_solids(&building);

        // Should have one stitch but incorrect interface
        assert_eq!(stitches.len(), 1);
        assert!(!stitches[0].interface_correct);
    }

    #[test]
    fn test_polygon_graph() {
        let building = make_adjacent_building();
        let params = GraphParams {
            level: GraphLevel::Polygon,
            facing: true,
            touching: false,
        };

        let graph = get_graph(&building, params);

        // Two facing polygons should be in the graph
        assert_eq!(graph.len(), 2);
    }

    #[test]
    fn test_wall_graph() {
        let building = make_adjacent_building();
        let params = GraphParams {
            level: GraphLevel::Wall,
            facing: true,
            touching: false,
        };

        let graph = get_graph(&building, params);

        // Two walls (wall_1 and wall_3) should be connected
        assert_eq!(graph.len(), 2);
    }
}
