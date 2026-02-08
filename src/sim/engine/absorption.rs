use crate::sim::materials::NUM_OCTAVE_BANDS;

/// Defines how ray energy is absorbed upon surface interaction.
pub trait AbsorptionModel {
    /// Apply absorption to ray energy at a given surface.
    /// Returns the remaining energy after absorption.
    fn apply(&self, energy: f64, surface_index: usize) -> f64;
}

/// Scalar absorption: single coefficient per surface.
pub struct ScalarAbsorption {
    /// Absorption coefficient per polygon index.
    pub coefficients: Vec<f64>,
}

impl ScalarAbsorption {
    pub fn new(coefficients: Vec<f64>) -> Self {
        Self { coefficients }
    }
}

impl AbsorptionModel for ScalarAbsorption {
    fn apply(&self, energy: f64, surface_index: usize) -> f64 {
        let alpha = self.coefficients.get(surface_index).copied().unwrap_or(0.0);
        energy * (1.0 - alpha)
    }
}

/// Frequency-dependent absorption: per-band coefficients per surface.
pub struct FrequencyDependentAbsorption {
    /// Absorption coefficients per polygon index, per octave band.
    pub coefficients: Vec<[f64; NUM_OCTAVE_BANDS]>,
}

impl FrequencyDependentAbsorption {
    pub fn new(coefficients: Vec<[f64; NUM_OCTAVE_BANDS]>) -> Self {
        Self { coefficients }
    }

    /// Apply frequency-dependent absorption to a multi-band energy array.
    pub fn apply_bands(
        &self,
        energy: &[f64; NUM_OCTAVE_BANDS],
        surface_index: usize,
    ) -> [f64; NUM_OCTAVE_BANDS] {
        let default = [0.0; NUM_OCTAVE_BANDS];
        let alpha = self.coefficients.get(surface_index).unwrap_or(&default);
        let mut result = [0.0; NUM_OCTAVE_BANDS];
        for i in 0..NUM_OCTAVE_BANDS {
            result[i] = energy[i] * (1.0 - alpha[i]);
        }
        result
    }
}

/// Scalar interface for frequency-dependent absorption (uses mean across bands).
impl AbsorptionModel for FrequencyDependentAbsorption {
    fn apply(&self, energy: f64, surface_index: usize) -> f64 {
        let default = [0.0; NUM_OCTAVE_BANDS];
        let alpha = self.coefficients.get(surface_index).unwrap_or(&default);
        let mean_alpha: f64 = alpha.iter().sum::<f64>() / NUM_OCTAVE_BANDS as f64;
        energy * (1.0 - mean_alpha)
    }
}

/// Air absorption model based on ISO 9613-1.
///
/// Applies distance-dependent attenuation per frequency band.
pub struct AirAbsorption {
    /// Attenuation coefficients in dB/m per octave band.
    pub attenuation_per_meter: [f64; NUM_OCTAVE_BANDS],
}

impl AirAbsorption {
    /// Creates an air absorption model for standard conditions (20°C, 50% RH).
    pub fn standard() -> Self {
        // Pure-tone attenuation coefficients in dB/m at 20°C, 50% RH, 101.325 kPa
        // computed from ISO 9613-1:1993 for octave band center frequencies 125-4000 Hz
        Self {
            attenuation_per_meter: [0.000440, 0.001310, 0.002728, 0.004665, 0.009887, 0.029666],
        }
    }

    /// Apply air absorption over a given distance.
    /// Returns the energy multiplier per band.
    pub fn apply_distance(&self, distance: f64) -> [f64; NUM_OCTAVE_BANDS] {
        let mut factors = [0.0; NUM_OCTAVE_BANDS];
        for (i, factor) in factors.iter_mut().enumerate().take(NUM_OCTAVE_BANDS) {
            // Convert dB/m attenuation to linear factor
            let db_loss = self.attenuation_per_meter[i] * distance;
            *factor = 10.0_f64.powf(-db_loss / 10.0);
        }
        factors
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scalar_absorption() {
        let abs = ScalarAbsorption::new(vec![0.2, 0.5]);
        assert!((abs.apply(1.0, 0) - 0.8).abs() < 1e-10);
        assert!((abs.apply(1.0, 1) - 0.5).abs() < 1e-10);
        // Unknown surface: no absorption
        assert!((abs.apply(1.0, 99) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_frequency_dependent_absorption() {
        let coeffs = vec![[0.1, 0.2, 0.3, 0.4, 0.5, 0.6]];
        let abs = FrequencyDependentAbsorption::new(coeffs);

        let energy = [1.0; NUM_OCTAVE_BANDS];
        let result = abs.apply_bands(&energy, 0);
        assert!((result[0] - 0.9).abs() < 1e-10);
        assert!((result[5] - 0.4).abs() < 1e-10);
    }

    #[test]
    fn test_air_absorption_standard() {
        let air = AirAbsorption::standard();
        let factors = air.apply_distance(10.0);
        // Low frequency should have less attenuation than high frequency
        assert!(factors[0] > factors[5]);
        // All factors should be between 0 and 1
        for f in &factors {
            assert!(*f > 0.0 && *f <= 1.0);
        }
    }

    #[test]
    fn test_air_absorption_zero_distance() {
        let air = AirAbsorption::standard();
        let factors = air.apply_distance(0.0);
        for f in &factors {
            assert!((*f - 1.0).abs() < 1e-10);
        }
    }
}
