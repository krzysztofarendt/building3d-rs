use crate::sim::materials::NUM_OCTAVE_BANDS;

use super::receiver::Receiver;

/// An impulse response extracted from a receiver's energy-time histogram.
#[derive(Debug, Clone)]
pub struct ImpulseResponse {
    /// Time resolution in seconds.
    pub time_resolution: f64,
    /// Energy per time bin per frequency band.
    pub bands: Vec<[f64; NUM_OCTAVE_BANDS]>,
    /// Broadband (summed) energy per time bin.
    pub broadband: Vec<f64>,
}

impl ImpulseResponse {
    /// Extracts an impulse response from a receiver.
    pub fn from_receiver(receiver: &Receiver) -> Self {
        let bands: Vec<[f64; NUM_OCTAVE_BANDS]> = receiver.histogram().to_vec();
        let broadband: Vec<f64> = receiver.scalar_histogram().to_vec();
        Self {
            time_resolution: receiver.time_resolution,
            bands,
            broadband,
        }
    }

    /// Creates an impulse response from raw data.
    pub fn new(
        time_resolution: f64,
        bands: Vec<[f64; NUM_OCTAVE_BANDS]>,
        broadband: Vec<f64>,
    ) -> Self {
        Self {
            time_resolution,
            bands,
            broadband,
        }
    }

    /// Number of time samples.
    pub fn len(&self) -> usize {
        self.broadband.len()
    }

    /// Whether the impulse response is empty.
    pub fn is_empty(&self) -> bool {
        self.broadband.is_empty()
    }

    /// Returns the time axis in seconds.
    pub fn time_axis(&self) -> Vec<f64> {
        (0..self.broadband.len())
            .map(|i| i as f64 * self.time_resolution)
            .collect()
    }

    /// Returns the total energy in a specific band.
    pub fn band_energy(&self, band: usize) -> f64 {
        self.bands.iter().map(|b| b[band]).sum()
    }

    /// Returns the total broadband energy.
    pub fn total_energy(&self) -> f64 {
        self.broadband.iter().sum()
    }

    /// Converts to a pressure-squared time series for a given sample rate.
    ///
    /// Each sample is the energy in the corresponding time bin,
    /// distributed across the samples that fall within that bin.
    pub fn to_time_series(&self, sample_rate: f64) -> Vec<f64> {
        let total_time = self.broadband.len() as f64 * self.time_resolution;
        let num_samples = (total_time * sample_rate).ceil() as usize;
        let mut output = vec![0.0; num_samples];

        let samples_per_bin = (self.time_resolution * sample_rate).max(1.0);

        for (bin, &energy) in self.broadband.iter().enumerate() {
            let start_sample = (bin as f64 * self.time_resolution * sample_rate) as usize;
            let end_sample =
                ((bin as f64 + 1.0) * self.time_resolution * sample_rate).ceil() as usize;
            let end_sample = end_sample.min(num_samples);

            let energy_per_sample = energy / samples_per_bin;
            for sample in &mut output[start_sample..end_sample] {
                *sample = energy_per_sample;
            }
        }

        output
    }

    /// Schroeder backward integration for a given band.
    ///
    /// Returns the energy decay curve in dB, normalized to 0 dB at time 0.
    pub fn schroeder_decay(&self, band: usize) -> Vec<f64> {
        let band_data: Vec<f64> = self.bands.iter().map(|b| b[band]).collect();
        let total: f64 = band_data.iter().sum();
        if total <= 0.0 {
            return vec![f64::NEG_INFINITY; band_data.len()];
        }

        let mut decay = vec![0.0; band_data.len()];
        let mut cumulative = total;
        for (d, &e) in decay.iter_mut().zip(band_data.iter()) {
            *d = 10.0 * (cumulative / total).log10();
            cumulative -= e;
        }

        decay
    }

    /// Schroeder backward integration for broadband.
    pub fn schroeder_decay_broadband(&self) -> Vec<f64> {
        let total: f64 = self.broadband.iter().sum();
        if total <= 0.0 {
            return vec![f64::NEG_INFINITY; self.broadband.len()];
        }

        let mut decay = vec![0.0; self.broadband.len()];
        let mut cumulative = total;
        for (d, &e) in decay.iter_mut().zip(self.broadband.iter()) {
            *d = 10.0 * (cumulative / total).log10();
            cumulative -= e;
        }

        decay
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Point;

    #[test]
    fn test_ir_from_receiver() {
        let mut receiver = Receiver::new(Point::new(0.0, 0.0, 0.0), 0.5, 0.001, 0.1);
        receiver.record_scalar_hit(0.005, 0.8);
        receiver.record_scalar_hit(0.010, 0.3);

        let ir = ImpulseResponse::from_receiver(&receiver);
        assert_eq!(ir.len(), 100);
        assert!((ir.total_energy() - 1.1).abs() < 1e-10);
    }

    #[test]
    fn test_ir_time_axis() {
        let ir = ImpulseResponse::new(0.001, vec![[0.0; NUM_OCTAVE_BANDS]; 10], vec![0.0; 10]);
        let t = ir.time_axis();
        assert_eq!(t.len(), 10);
        assert!((t[0] - 0.0).abs() < 1e-10);
        assert!((t[9] - 0.009).abs() < 1e-10);
    }

    #[test]
    fn test_schroeder_decay() {
        // Create simple exponential decay
        let mut broadband = vec![0.0; 100];
        broadband[0] = 1.0;
        broadband[10] = 0.5;
        broadband[20] = 0.25;

        let ir = ImpulseResponse::new(0.001, vec![[0.0; NUM_OCTAVE_BANDS]; 100], broadband);

        let decay = ir.schroeder_decay_broadband();
        // Decay should start at 0 dB
        assert!((decay[0] - 0.0).abs() < 1e-10);
        // Decay should be monotonically decreasing
        for i in 1..decay.len() {
            assert!(decay[i] <= decay[i - 1] + 1e-10);
        }
    }

    #[test]
    fn test_to_time_series() {
        let mut broadband = vec![0.0; 10];
        broadband[0] = 1.0;
        let ir = ImpulseResponse::new(0.001, vec![[0.0; NUM_OCTAVE_BANDS]; 10], broadband);

        let ts = ir.to_time_series(44100.0);
        assert!(!ts.is_empty());
        // Total energy in time series should approximate the IR energy
        let total: f64 = ts.iter().sum();
        assert!(total > 0.0);
    }
}
