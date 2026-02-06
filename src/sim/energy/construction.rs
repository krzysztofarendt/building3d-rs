use crate::sim::materials::{Layer, ThermalMaterial};

/// A wall construction defined by material layers (outside to inside).
///
/// Computes steady-state thermal resistance and U-value following
/// ISO 6946 simplified method.
#[derive(Debug, Clone)]
pub struct WallConstruction {
    pub name: String,
    pub layers: Vec<Layer>,
    /// External surface resistance in m^2*K/W (default: 0.04 for walls).
    pub r_se: f64,
    /// Internal surface resistance in m^2*K/W (default: 0.13 for walls).
    pub r_si: f64,
}

impl WallConstruction {
    pub fn new(name: &str, layers: Vec<Layer>) -> Self {
        Self {
            name: name.to_string(),
            layers,
            r_se: 0.04,
            r_si: 0.13,
        }
    }

    /// Creates a construction for a floor (different surface resistances).
    pub fn floor(name: &str, layers: Vec<Layer>) -> Self {
        Self {
            name: name.to_string(),
            layers,
            r_se: 0.04,
            r_si: 0.17, // ISO 6946 for downward heat flow
        }
    }

    /// Creates a construction for a roof/ceiling.
    pub fn roof(name: &str, layers: Vec<Layer>) -> Self {
        Self {
            name: name.to_string(),
            layers,
            r_se: 0.04,
            r_si: 0.10, // ISO 6946 for upward heat flow
        }
    }

    /// Total thermal resistance in m^2*K/W (including surface resistances).
    pub fn total_resistance(&self) -> f64 {
        let r_layers: f64 = self
            .layers
            .iter()
            .map(|l| {
                if l.conductivity > 0.0 {
                    l.thickness / l.conductivity
                } else {
                    0.0
                }
            })
            .sum();
        self.r_se + r_layers + self.r_si
    }

    /// U-value in W/(m^2*K).
    pub fn u_value(&self) -> f64 {
        let r = self.total_resistance();
        if r > 0.0 { 1.0 / r } else { 0.0 }
    }

    /// Total thermal capacity per unit area in J/(m^2*K).
    pub fn thermal_capacity(&self) -> f64 {
        self.layers
            .iter()
            .map(|l| l.density * l.specific_heat * l.thickness)
            .sum()
    }

    /// Converts this construction to a ThermalMaterial.
    pub fn to_thermal_material(&self) -> ThermalMaterial {
        ThermalMaterial {
            name: self.name.clone(),
            u_value: self.u_value(),
            thermal_capacity: self.thermal_capacity(),
            layers: self.layers.clone(),
        }
    }

    /// Creates a simple single-layer construction.
    pub fn single_layer(name: &str, layer: Layer) -> Self {
        Self::new(name, vec![layer])
    }
}

/// Preset constructions for common wall types.
pub fn concrete_wall() -> WallConstruction {
    WallConstruction::new(
        "concrete_wall",
        vec![Layer {
            name: "concrete".to_string(),
            thickness: 0.20,
            conductivity: 1.4,
            density: 2300.0,
            specific_heat: 880.0,
        }],
    )
}

pub fn insulated_wall() -> WallConstruction {
    WallConstruction::new(
        "insulated_wall",
        vec![
            Layer {
                name: "plaster_ext".to_string(),
                thickness: 0.02,
                conductivity: 0.87,
                density: 1800.0,
                specific_heat: 840.0,
            },
            Layer {
                name: "insulation".to_string(),
                thickness: 0.10,
                conductivity: 0.04,
                density: 30.0,
                specific_heat: 1030.0,
            },
            Layer {
                name: "concrete".to_string(),
                thickness: 0.15,
                conductivity: 1.4,
                density: 2300.0,
                specific_heat: 880.0,
            },
            Layer {
                name: "plaster_int".to_string(),
                thickness: 0.015,
                conductivity: 0.87,
                density: 1800.0,
                specific_heat: 840.0,
            },
        ],
    )
}

pub fn double_glazing() -> WallConstruction {
    // Simplified: treating glazing as a single equivalent layer
    WallConstruction::new(
        "double_glazing",
        vec![Layer {
            name: "glass_unit".to_string(),
            thickness: 0.024,     // 4mm + 16mm gap + 4mm
            conductivity: 0.035,  // effective (gap-dominated)
            density: 2500.0,      // glass density for capacity
            specific_heat: 840.0, // glass specific heat
        }],
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_concrete_wall_u_value() {
        let wall = concrete_wall();
        let u = wall.u_value();
        // R = 0.04 + 0.20/1.4 + 0.13 = 0.04 + 0.1429 + 0.13 = 0.3129
        // U = 1/0.3129 ≈ 3.20
        assert!((u - 3.20).abs() < 0.1, "U-value should be ~3.2, got {u}");
    }

    #[test]
    fn test_insulated_wall_u_value() {
        let wall = insulated_wall();
        let u = wall.u_value();
        // Insulated wall should have much lower U-value
        assert!(u < 0.5, "Insulated wall should have U < 0.5, got {u}");
        assert!(
            u > 0.1,
            "Insulated wall U-value should be realistic, got {u}"
        );
    }

    #[test]
    fn test_thermal_capacity() {
        let wall = concrete_wall();
        let cap = wall.thermal_capacity();
        // 2300 * 880 * 0.20 = 404,800 J/(m^2*K)
        assert!(
            (cap - 404800.0).abs() < 1.0,
            "Capacity should be ~404800, got {cap}"
        );
    }

    #[test]
    fn test_to_thermal_material() {
        let wall = insulated_wall();
        let mat = wall.to_thermal_material();
        assert!((mat.u_value - wall.u_value()).abs() < 1e-10);
        assert!((mat.thermal_capacity - wall.thermal_capacity()).abs() < 1e-10);
    }

    #[test]
    fn test_surface_resistances() {
        let wall = WallConstruction::new("test", vec![]);
        let floor = WallConstruction::floor("test", vec![]);
        let roof = WallConstruction::roof("test", vec![]);
        assert!((wall.r_si - 0.13).abs() < 1e-10);
        assert!((floor.r_si - 0.17).abs() < 1e-10);
        assert!((roof.r_si - 0.10).abs() < 1e-10);
    }
}
