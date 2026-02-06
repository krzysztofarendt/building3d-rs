use std::collections::HashSet;

use crate::Building;
use crate::geom::building::graph::{GraphLevel, GraphParams, get_graph};

/// Finds polygons that are transparent (internal interfaces between solids in the same zone).
///
/// Uses the building's polygon-level facing graph to find pairs of facing polygons.
/// A polygon is transparent if it faces another polygon and both belong to the same zone.
pub fn find_transparent_polygons(building: &Building) -> HashSet<String> {
    let params = GraphParams {
        level: GraphLevel::Polygon,
        facing: true,
        touching: false,
    };

    let graph = get_graph(building, params);
    let mut transparent = HashSet::new();

    for (path1, neighbors) in &graph {
        let zone1 = path1.split('/').next().unwrap_or("");

        for path2 in neighbors {
            let zone2 = path2.split('/').next().unwrap_or("");

            if zone1 == zone2 {
                transparent.insert(path1.clone());
                transparent.insert(path2.clone());
            }
        }
    }

    transparent
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Solid, Zone};

    #[test]
    fn test_transparent_same_zone() {
        // Two adjacent boxes in the same zone — their shared face should be transparent
        let s0 = Solid::from_box(1.0, 1.0, 1.0, None, "s0").unwrap();
        let s1 = Solid::from_box(1.0, 1.0, 1.0, Some((1.0, 0.0, 0.0)), "s1").unwrap();
        let zone = Zone::new("z", vec![s0, s1]).unwrap();
        let building = Building::new("b", vec![zone]).unwrap();

        let transparent = find_transparent_polygons(&building);
        assert!(
            !transparent.is_empty(),
            "Should find transparent polygons between adjacent solids in the same zone"
        );
        // All transparent paths should start with "z/"
        for path in &transparent {
            assert!(path.starts_with("z/"));
        }
    }

    #[test]
    fn test_no_transparent_different_zones() {
        // Two adjacent boxes in different zones — no transparent polygons
        let s0 = Solid::from_box(1.0, 1.0, 1.0, None, "s0").unwrap();
        let s1 = Solid::from_box(1.0, 1.0, 1.0, Some((1.0, 0.0, 0.0)), "s1").unwrap();
        let z0 = Zone::new("z0", vec![s0]).unwrap();
        let z1 = Zone::new("z1", vec![s1]).unwrap();
        let building = Building::new("b", vec![z0, z1]).unwrap();

        let transparent = find_transparent_polygons(&building);
        assert!(
            transparent.is_empty(),
            "Should not find transparent polygons between different zones"
        );
    }
}
