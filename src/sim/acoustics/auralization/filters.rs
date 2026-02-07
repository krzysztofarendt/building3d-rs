use crate::sim::materials::{NUM_OCTAVE_BANDS, OCTAVE_BAND_FREQUENCIES};

/// Second-order (biquad) filter coefficients in Direct Form I.
///
/// Transfer function: H(z) = (b0 + b1*z^-1 + b2*z^-2) / (1 + a1*z^-1 + a2*z^-2)
#[derive(Debug, Clone, Copy)]
pub struct BiquadCoeffs {
    pub b0: f64,
    pub b1: f64,
    pub b2: f64,
    pub a1: f64,
    pub a2: f64,
}

/// Biquad filter state for Direct Form I processing.
#[derive(Debug, Clone)]
pub struct BiquadState {
    x1: f64,
    x2: f64,
    y1: f64,
    y2: f64,
}

impl Default for BiquadState {
    fn default() -> Self {
        Self {
            x1: 0.0,
            x2: 0.0,
            y1: 0.0,
            y2: 0.0,
        }
    }
}

impl BiquadState {
    pub fn new() -> Self {
        Self::default()
    }

    /// Process a single sample through the filter.
    pub fn process(&mut self, sample: f64, coeffs: &BiquadCoeffs) -> f64 {
        let output =
            coeffs.b0 * sample + coeffs.b1 * self.x1 + coeffs.b2 * self.x2
                - coeffs.a1 * self.y1
                - coeffs.a2 * self.y2;
        self.x2 = self.x1;
        self.x1 = sample;
        self.y2 = self.y1;
        self.y1 = output;
        output
    }
}

/// Designs a 2nd-order Butterworth bandpass filter for a given octave band.
///
/// The bandwidth is one octave centered at `center_freq`.
pub fn design_octave_bandpass(center_freq: f64, sample_rate: f64) -> BiquadCoeffs {
    // One-octave bandwidth: lower = center / sqrt(2), upper = center * sqrt(2)
    let sqrt2 = std::f64::consts::SQRT_2;
    let f_low = center_freq / sqrt2;
    let f_high = center_freq * sqrt2;

    // Bandwidth in Hz
    let bw = f_high - f_low;

    // Cookbook BPF formula (constant 0 dB peak gain)
    let w0 = 2.0 * std::f64::consts::PI * center_freq / sample_rate;
    let sin_w0 = w0.sin();
    let cos_w0 = w0.cos();
    let alpha = sin_w0 / (2.0 * (center_freq / bw));

    let b0 = alpha;
    let b1 = 0.0;
    let b2 = -alpha;
    let a0 = 1.0 + alpha;
    let a1 = -2.0 * cos_w0;
    let a2 = 1.0 - alpha;

    BiquadCoeffs {
        b0: b0 / a0,
        b1: b1 / a0,
        b2: b2 / a0,
        a1: a1 / a0,
        a2: a2 / a0,
    }
}

/// Filters a signal through the 6 standard octave bands.
///
/// Returns an array of 6 filtered signals, one per band (125, 250, 500, 1000, 2000, 4000 Hz).
pub fn filter_octave_bands(signal: &[f64], sample_rate: f64) -> [Vec<f64>; NUM_OCTAVE_BANDS] {
    std::array::from_fn(|band| {
        let coeffs = design_octave_bandpass(OCTAVE_BAND_FREQUENCIES[band], sample_rate);
        let mut state = BiquadState::new();
        signal.iter().map(|&s| state.process(s, &coeffs)).collect()
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_biquad_silence() {
        let coeffs = design_octave_bandpass(1000.0, 44100.0);
        let mut state = BiquadState::new();
        for _ in 0..100 {
            let out = state.process(0.0, &coeffs);
            assert!((out).abs() < 1e-15);
        }
    }

    #[test]
    fn test_bandpass_passes_center_frequency() {
        let sample_rate = 44100.0;
        let center = 1000.0;
        let coeffs = design_octave_bandpass(center, sample_rate);
        let mut state = BiquadState::new();

        // Generate a 1000 Hz sine wave
        let n = 44100; // 1 second
        let signal: Vec<f64> = (0..n)
            .map(|i| (2.0 * std::f64::consts::PI * center * i as f64 / sample_rate).sin())
            .collect();

        let output: Vec<f64> = signal.iter().map(|&s| state.process(s, &coeffs)).collect();

        // Measure energy in the steady-state portion (skip transient)
        let skip = n / 4;
        let input_energy: f64 = signal[skip..].iter().map(|x| x * x).sum::<f64>();
        let output_energy: f64 = output[skip..].iter().map(|x| x * x).sum::<f64>();

        let ratio = output_energy / input_energy;
        // Should pass through with close to unity gain at center
        assert!(
            ratio > 0.5,
            "Center frequency should pass through, ratio={ratio}"
        );
    }

    #[test]
    fn test_bandpass_rejects_out_of_band() {
        let sample_rate = 44100.0;
        let center = 1000.0;
        let coeffs = design_octave_bandpass(center, sample_rate);
        let mut state = BiquadState::new();

        // Generate a 100 Hz sine wave (well below the 1000 Hz band)
        let n = 44100;
        let signal: Vec<f64> = (0..n)
            .map(|i| (2.0 * std::f64::consts::PI * 100.0 * i as f64 / sample_rate).sin())
            .collect();

        let output: Vec<f64> = signal.iter().map(|&s| state.process(s, &coeffs)).collect();

        let skip = n / 4;
        let input_energy: f64 = signal[skip..].iter().map(|x| x * x).sum::<f64>();
        let output_energy: f64 = output[skip..].iter().map(|x| x * x).sum::<f64>();

        let ratio = output_energy / input_energy;
        // Should strongly attenuate
        assert!(
            ratio < 0.1,
            "Out-of-band frequency should be attenuated, ratio={ratio}"
        );
    }

    #[test]
    fn test_filter_octave_bands_length() {
        let signal = vec![0.0; 1000];
        let bands = filter_octave_bands(&signal, 44100.0);
        assert_eq!(bands.len(), NUM_OCTAVE_BANDS);
        for band in &bands {
            assert_eq!(band.len(), signal.len());
        }
    }

    #[test]
    fn test_band_selectivity() {
        let sample_rate = 44100.0;
        let n = 44100;

        // Generate a 500 Hz sine
        let signal: Vec<f64> = (0..n)
            .map(|i| (2.0 * std::f64::consts::PI * 500.0 * i as f64 / sample_rate).sin())
            .collect();

        let bands = filter_octave_bands(&signal, sample_rate);

        let skip = n / 4;
        let energies: Vec<f64> = bands
            .iter()
            .map(|b| b[skip..].iter().map(|x| x * x).sum::<f64>())
            .collect();

        // Band 2 (500 Hz) should have the most energy
        let max_band = energies
            .iter()
            .enumerate()
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            .unwrap()
            .0;
        assert_eq!(max_band, 2, "500 Hz signal should peak in band 2 (500 Hz)");
    }
}
