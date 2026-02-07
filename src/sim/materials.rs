use std::collections::HashMap;

/// Number of octave bands used for frequency-dependent simulation (125 Hz to 4 kHz).
pub const NUM_OCTAVE_BANDS: usize = 6;

/// Center frequencies of the octave bands in Hz.
pub const OCTAVE_BAND_FREQUENCIES: [f64; NUM_OCTAVE_BANDS] =
    [125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0];

/// Acoustic surface material with frequency-dependent absorption and scattering.
///
/// Absorption and scattering coefficients are specified per octave band
/// (125, 250, 500, 1000, 2000, 4000 Hz).
#[derive(Debug, Clone)]
pub struct AcousticMaterial {
    pub name: String,
    /// Absorption coefficients per octave band [0.0, 1.0].
    pub absorption: [f64; NUM_OCTAVE_BANDS],
    /// Scattering coefficients per octave band [0.0, 1.0].
    pub scattering: [f64; NUM_OCTAVE_BANDS],
}

/// Optical surface material for lighting simulation.
///
/// RGB reflectance, specular reflectance, and transmittance are specified
/// as [R, G, B] triplets in [0.0, 1.0].
#[derive(Debug, Clone)]
pub struct OpticalMaterial {
    pub name: String,
    /// Diffuse reflectance [R, G, B].
    pub diffuse_reflectance: [f64; 3],
    /// Specular reflectance [R, G, B].
    pub specular_reflectance: [f64; 3],
    /// Transmittance [R, G, B] (for glazing).
    pub transmittance: [f64; 3],
}

/// A single layer in a wall construction.
#[derive(Debug, Clone)]
pub struct Layer {
    pub name: String,
    /// Thickness in meters.
    pub thickness: f64,
    /// Thermal conductivity in W/(m*K).
    pub conductivity: f64,
    /// Density in kg/m^3.
    pub density: f64,
    /// Specific heat capacity in J/(kg*K).
    pub specific_heat: f64,
}

/// Thermal material properties for energy simulation.
#[derive(Debug, Clone)]
pub struct ThermalMaterial {
    pub name: String,
    /// Overall U-value in W/(m^2*K). Calculated from layers if provided.
    pub u_value: f64,
    /// Thermal capacity per unit area in J/(m^2*K).
    pub thermal_capacity: f64,
    /// Construction layers (outside to inside).
    pub layers: Vec<Layer>,
}

/// Combined material that can hold properties for any or all simulation domains.
#[derive(Debug, Clone)]
pub struct Material {
    pub name: String,
    pub acoustic: Option<AcousticMaterial>,
    pub optical: Option<OpticalMaterial>,
    pub thermal: Option<ThermalMaterial>,
    /// Whether this surface is glazing (for solar gain calculations).
    /// When true, the material library can identify glazing surfaces without
    /// relying on fragile name-pattern matching.
    pub is_glazing: bool,
}

/// Library of named materials with path-pattern assignment to building surfaces.
#[derive(Clone)]
pub struct MaterialLibrary {
    /// Material definitions by name.
    materials: HashMap<String, Material>,
    /// Assignments: path pattern -> material name.
    assignments: Vec<(String, String)>,
}

impl AcousticMaterial {
    pub fn new(
        name: &str,
        absorption: [f64; NUM_OCTAVE_BANDS],
        scattering: [f64; NUM_OCTAVE_BANDS],
    ) -> Self {
        Self {
            name: name.to_string(),
            absorption,
            scattering,
        }
    }

    /// Creates a material with uniform absorption across all bands.
    pub fn uniform(name: &str, absorption: f64, scattering: f64) -> Self {
        Self {
            name: name.to_string(),
            absorption: [absorption; NUM_OCTAVE_BANDS],
            scattering: [scattering; NUM_OCTAVE_BANDS],
        }
    }
}

impl OpticalMaterial {
    pub fn new(
        name: &str,
        diffuse_reflectance: [f64; 3],
        specular_reflectance: [f64; 3],
        transmittance: [f64; 3],
    ) -> Self {
        Self {
            name: name.to_string(),
            diffuse_reflectance,
            specular_reflectance,
            transmittance,
        }
    }

    /// Creates a purely diffuse material with no specular or transmittance.
    pub fn diffuse(name: &str, reflectance: [f64; 3]) -> Self {
        Self {
            name: name.to_string(),
            diffuse_reflectance: reflectance,
            specular_reflectance: [0.0; 3],
            transmittance: [0.0; 3],
        }
    }
}

impl Layer {
    pub fn new(
        name: &str,
        thickness: f64,
        conductivity: f64,
        density: f64,
        specific_heat: f64,
    ) -> Self {
        Self {
            name: name.to_string(),
            thickness,
            conductivity,
            density,
            specific_heat,
        }
    }

    /// Thermal resistance of this layer in m^2*K/W.
    pub fn resistance(&self) -> f64 {
        self.thickness / self.conductivity
    }
}

impl ThermalMaterial {
    pub fn new(name: &str, u_value: f64, thermal_capacity: f64, layers: Vec<Layer>) -> Self {
        Self {
            name: name.to_string(),
            u_value,
            thermal_capacity,
            layers,
        }
    }

    /// Creates a thermal material from layers, computing U-value and capacity.
    pub fn from_layers(name: &str, layers: Vec<Layer>) -> Self {
        let total_resistance: f64 = layers.iter().map(|l| l.resistance()).sum();
        let u_value = if total_resistance > 0.0 {
            1.0 / total_resistance
        } else {
            f64::INFINITY
        };
        let thermal_capacity: f64 = layers
            .iter()
            .map(|l| l.density * l.specific_heat * l.thickness)
            .sum();
        Self {
            name: name.to_string(),
            u_value,
            thermal_capacity,
            layers,
        }
    }
}

impl Material {
    pub fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            acoustic: None,
            optical: None,
            thermal: None,
            is_glazing: false,
        }
    }

    pub fn with_acoustic(mut self, acoustic: AcousticMaterial) -> Self {
        self.acoustic = Some(acoustic);
        self
    }

    pub fn with_optical(mut self, optical: OpticalMaterial) -> Self {
        self.optical = Some(optical);
        self
    }

    pub fn with_thermal(mut self, thermal: ThermalMaterial) -> Self {
        self.thermal = Some(thermal);
        self
    }

    /// Marks this material as glazing for solar gain calculations.
    pub fn with_glazing(mut self) -> Self {
        self.is_glazing = true;
        self
    }
}

impl MaterialLibrary {
    pub fn new() -> Self {
        Self {
            materials: HashMap::new(),
            assignments: Vec::new(),
        }
    }

    /// Adds a material to the library.
    pub fn add(&mut self, material: Material) {
        self.materials.insert(material.name.clone(), material);
    }

    /// Assigns a material to surfaces matching a path pattern (substring match).
    pub fn assign(&mut self, path_pattern: &str, material_name: &str) {
        self.assignments
            .push((path_pattern.to_string(), material_name.to_string()));
    }

    /// Looks up the material for a given full path.
    ///
    /// Searches assignments in reverse order (last assignment wins).
    pub fn lookup(&self, full_path: &str) -> Option<&Material> {
        for (pattern, material_name) in self.assignments.iter().rev() {
            if full_path.contains(pattern.as_str()) {
                return self.materials.get(material_name);
            }
        }
        None
    }

    /// Returns a reference to a material by name.
    pub fn get(&self, name: &str) -> Option<&Material> {
        self.materials.get(name)
    }

    /// Creates a library pre-populated with common building materials.
    pub fn with_presets() -> Self {
        let mut lib = Self::new();

        // Concrete
        lib.add(
            Material::new("concrete")
                .with_acoustic(AcousticMaterial::new(
                    "concrete",
                    [0.01, 0.01, 0.02, 0.02, 0.02, 0.03],
                    [0.10, 0.10, 0.10, 0.10, 0.10, 0.10],
                ))
                .with_optical(OpticalMaterial::diffuse("concrete", [0.40, 0.40, 0.40]))
                .with_thermal(ThermalMaterial::from_layers(
                    "concrete",
                    vec![Layer::new("concrete_200mm", 0.200, 1.40, 2300.0, 1000.0)],
                )),
        );

        // Glass (single pane)
        lib.add(
            Material::new("glass")
                .with_acoustic(AcousticMaterial::new(
                    "glass",
                    [0.18, 0.06, 0.04, 0.03, 0.02, 0.02],
                    [0.05, 0.05, 0.05, 0.05, 0.05, 0.05],
                ))
                .with_optical(OpticalMaterial::new(
                    "glass",
                    [0.08, 0.08, 0.08],
                    [0.10, 0.10, 0.10],
                    [0.75, 0.75, 0.75],
                ))
                .with_thermal(ThermalMaterial::from_layers(
                    "glass",
                    vec![Layer::new("glass_6mm", 0.006, 1.05, 2500.0, 840.0)],
                ))
                .with_glazing(),
        );

        // Gypsum board
        lib.add(
            Material::new("gypsum")
                .with_acoustic(AcousticMaterial::new(
                    "gypsum",
                    [0.29, 0.10, 0.05, 0.04, 0.07, 0.09],
                    [0.10, 0.10, 0.10, 0.10, 0.10, 0.10],
                ))
                .with_optical(OpticalMaterial::diffuse("gypsum", [0.75, 0.75, 0.75]))
                .with_thermal(ThermalMaterial::from_layers(
                    "gypsum",
                    vec![Layer::new("gypsum_13mm", 0.013, 0.25, 900.0, 1000.0)],
                )),
        );

        // Carpet
        lib.add(
            Material::new("carpet")
                .with_acoustic(AcousticMaterial::new(
                    "carpet",
                    [0.02, 0.06, 0.14, 0.37, 0.60, 0.65],
                    [0.40, 0.40, 0.40, 0.50, 0.50, 0.50],
                ))
                .with_optical(OpticalMaterial::diffuse("carpet", [0.30, 0.30, 0.30]))
                .with_thermal(ThermalMaterial::from_layers(
                    "carpet",
                    vec![Layer::new("carpet_10mm", 0.010, 0.06, 200.0, 1300.0)],
                )),
        );

        // Wood
        lib.add(
            Material::new("wood")
                .with_acoustic(AcousticMaterial::new(
                    "wood",
                    [0.15, 0.11, 0.10, 0.07, 0.06, 0.07],
                    [0.10, 0.10, 0.10, 0.10, 0.10, 0.10],
                ))
                .with_optical(OpticalMaterial::diffuse("wood", [0.45, 0.35, 0.25]))
                .with_thermal(ThermalMaterial::from_layers(
                    "wood",
                    vec![Layer::new("wood_20mm", 0.020, 0.13, 600.0, 1700.0)],
                )),
        );

        // Metal (steel)
        lib.add(
            Material::new("metal")
                .with_acoustic(AcousticMaterial::new(
                    "metal",
                    [0.04, 0.04, 0.03, 0.03, 0.03, 0.03],
                    [0.05, 0.05, 0.05, 0.05, 0.05, 0.05],
                ))
                .with_optical(OpticalMaterial::new(
                    "metal",
                    [0.10, 0.10, 0.10],
                    [0.70, 0.70, 0.70],
                    [0.0, 0.0, 0.0],
                ))
                .with_thermal(ThermalMaterial::from_layers(
                    "metal",
                    vec![Layer::new("steel_3mm", 0.003, 50.0, 7800.0, 450.0)],
                )),
        );

        lib
    }
}

impl Default for MaterialLibrary {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_acoustic_material_uniform() {
        let m = AcousticMaterial::uniform("test", 0.5, 0.1);
        assert_eq!(m.absorption, [0.5; NUM_OCTAVE_BANDS]);
        assert_eq!(m.scattering, [0.1; NUM_OCTAVE_BANDS]);
    }

    #[test]
    fn test_layer_resistance() {
        let layer = Layer::new("concrete", 0.2, 1.4, 2300.0, 1000.0);
        let r = layer.resistance();
        assert!((r - 0.2 / 1.4).abs() < 1e-10);
    }

    #[test]
    fn test_thermal_from_layers() {
        let layers = vec![
            Layer::new("outer", 0.1, 1.0, 2000.0, 1000.0),
            Layer::new("insulation", 0.05, 0.04, 30.0, 1400.0),
            Layer::new("inner", 0.013, 0.25, 900.0, 1000.0),
        ];
        let tm = ThermalMaterial::from_layers("wall", layers);
        // R = 0.1/1.0 + 0.05/0.04 + 0.013/0.25 = 0.1 + 1.25 + 0.052 = 1.402
        let expected_r = 0.1 / 1.0 + 0.05 / 0.04 + 0.013 / 0.25;
        assert!((tm.u_value - 1.0 / expected_r).abs() < 1e-6);
        assert!(tm.thermal_capacity > 0.0);
    }

    #[test]
    fn test_material_builder() {
        let m = Material::new("test")
            .with_acoustic(AcousticMaterial::uniform("test", 0.3, 0.1))
            .with_optical(OpticalMaterial::diffuse("test", [0.5, 0.5, 0.5]));
        assert!(m.acoustic.is_some());
        assert!(m.optical.is_some());
        assert!(m.thermal.is_none());
    }

    #[test]
    fn test_material_library_lookup() {
        let mut lib = MaterialLibrary::new();
        lib.add(
            Material::new("concrete")
                .with_acoustic(AcousticMaterial::uniform("concrete", 0.02, 0.1)),
        );
        lib.add(
            Material::new("carpet").with_acoustic(AcousticMaterial::uniform("carpet", 0.5, 0.4)),
        );

        lib.assign("floor", "carpet");
        lib.assign("wall", "concrete");

        let m = lib.lookup("zone/room/wall/north");
        assert!(m.is_some());
        assert_eq!(m.unwrap().name, "concrete");

        let m = lib.lookup("zone/room/floor/floor");
        assert!(m.is_some());
        assert_eq!(m.unwrap().name, "carpet");

        let m = lib.lookup("zone/room/ceiling/top");
        assert!(m.is_none());
    }

    #[test]
    fn test_material_library_last_wins() {
        let mut lib = MaterialLibrary::new();
        lib.add(Material::new("mat_a"));
        lib.add(Material::new("mat_b"));

        lib.assign("wall", "mat_a");
        lib.assign("wall", "mat_b");

        let m = lib.lookup("zone/room/wall/north");
        assert_eq!(m.unwrap().name, "mat_b");
    }

    #[test]
    fn test_presets() {
        let lib = MaterialLibrary::with_presets();
        assert!(lib.get("concrete").is_some());
        assert!(lib.get("glass").is_some());
        assert!(lib.get("gypsum").is_some());
        assert!(lib.get("carpet").is_some());
        assert!(lib.get("wood").is_some());
        assert!(lib.get("metal").is_some());

        let concrete = lib.get("concrete").unwrap();
        assert!(concrete.acoustic.is_some());
        assert!(concrete.optical.is_some());
        assert!(concrete.thermal.is_some());
        assert!(!concrete.is_glazing);

        let glass = lib.get("glass").unwrap();
        assert!(glass.is_glazing, "Glass preset should be marked as glazing");
    }

    #[test]
    fn test_octave_band_frequencies() {
        assert_eq!(OCTAVE_BAND_FREQUENCIES.len(), NUM_OCTAVE_BANDS);
        assert!((OCTAVE_BAND_FREQUENCIES[0] - 125.0).abs() < 1e-10);
        assert!((OCTAVE_BAND_FREQUENCIES[5] - 4000.0).abs() < 1e-10);
    }
}
