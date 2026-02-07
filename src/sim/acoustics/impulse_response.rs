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
    /// Each sample receives the portion of each time bin's energy that overlaps
    /// the sample's time interval. This avoids rounding artifacts (e.g. bins
    /// spanning 44 vs 45 samples) and preserves total energy under resampling.
    pub fn to_time_series(&self, sample_rate: f64) -> Vec<f64> {
        if self.is_empty() || sample_rate <= 0.0 {
            return Vec::new();
        }
        if self.time_resolution <= 0.0 {
            // No meaningful bin duration -> cannot distribute energy in time.
            return vec![0.0; self.broadband.len()];
        }

        let total_time = self.broadband.len() as f64 * self.time_resolution;
        let num_samples = (total_time * sample_rate).ceil() as usize;
        let mut output = vec![0.0; num_samples];

        let sample_dt = 1.0 / sample_rate;

        for (bin, &energy) in self.broadband.iter().enumerate() {
            if energy == 0.0 {
                continue;
            }

            let bin_t0 = bin as f64 * self.time_resolution;
            let bin_t1 = (bin as f64 + 1.0) * self.time_resolution;

            let start_sample = (bin_t0 * sample_rate).floor() as usize;
            let end_sample = (bin_t1 * sample_rate).ceil() as usize;
            let start_sample = start_sample.min(num_samples);
            let end_sample = end_sample.min(num_samples);

            if start_sample >= end_sample {
                continue;
            }

            let energy_density = energy / self.time_resolution; // energy per second within the bin
            for (offset, out) in output[start_sample..end_sample].iter_mut().enumerate() {
                let sample_idx = start_sample + offset;
                let sample_t0 = sample_idx as f64 * sample_dt;
                let sample_t1 = sample_t0 + sample_dt;
                let overlap = bin_t1.min(sample_t1) - bin_t0.max(sample_t0);
                if overlap > 0.0 {
                    *out += energy_density * overlap;
                }
            }
        }

        output
    }

    /// Converts a single frequency band to a pressure-squared time series.
    ///
    /// Same algorithm as [`to_time_series`](Self::to_time_series) but uses
    /// `self.bands[i][band]` instead of `self.broadband[i]`.
    pub fn band_to_time_series(&self, band: usize, sample_rate: f64) -> Vec<f64> {
        if self.bands.is_empty() || sample_rate <= 0.0 || band >= NUM_OCTAVE_BANDS {
            return Vec::new();
        }
        if self.time_resolution <= 0.0 {
            return vec![0.0; self.bands.len()];
        }

        let band_data: Vec<f64> = self.bands.iter().map(|b| b[band]).collect();
        let total_time = band_data.len() as f64 * self.time_resolution;
        let num_samples = (total_time * sample_rate).ceil() as usize;
        let mut output = vec![0.0; num_samples];

        let sample_dt = 1.0 / sample_rate;

        for (bin, &energy) in band_data.iter().enumerate() {
            if energy == 0.0 {
                continue;
            }

            let bin_t0 = bin as f64 * self.time_resolution;
            let bin_t1 = (bin as f64 + 1.0) * self.time_resolution;

            let start_sample = (bin_t0 * sample_rate).floor() as usize;
            let end_sample = (bin_t1 * sample_rate).ceil() as usize;
            let start_sample = start_sample.min(num_samples);
            let end_sample = end_sample.min(num_samples);

            if start_sample >= end_sample {
                continue;
            }

            let energy_density = energy / self.time_resolution;
            for (sample_idx, out) in output[start_sample..end_sample].iter_mut().enumerate() {
                let sample_idx = sample_idx + start_sample;
                let sample_t0 = sample_idx as f64 * sample_dt;
                let sample_t1 = sample_t0 + sample_dt;
                let overlap = bin_t1.min(sample_t1) - bin_t0.max(sample_t0);
                if overlap > 0.0 {
                    *out += energy_density * overlap;
                }
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
    fn test_band_to_time_series() {
        let mut bands = vec![[0.0; NUM_OCTAVE_BANDS]; 10];
        bands[0][2] = 0.7; // 500 Hz band
        bands[0][4] = 0.3; // 2000 Hz band
        let broadband = vec![0.0; 10];
        let ir = ImpulseResponse::new(0.001, bands, broadband);

        let ts_500 = ir.band_to_time_series(2, 44100.0);
        let ts_2000 = ir.band_to_time_series(4, 44100.0);
        assert!(!ts_500.is_empty());
        assert!(!ts_2000.is_empty());

        let total_500: f64 = ts_500.iter().sum();
        let total_2000: f64 = ts_2000.iter().sum();
        assert!((total_500 - 0.7).abs() < 1e-10);
        assert!((total_2000 - 0.3).abs() < 1e-10);
    }

    #[test]
    fn test_band_to_time_series_with_empty_broadband() {
        let mut bands = vec![[0.0; NUM_OCTAVE_BANDS]; 10];
        bands[0][3] = 0.4;
        let ir = ImpulseResponse::new(0.001, bands, vec![]);

        let ts = ir.band_to_time_series(3, 44100.0);
        assert!(!ts.is_empty());
        let total: f64 = ts.iter().sum();
        assert!((total - 0.4).abs() < 1e-10);
    }

    #[test]
    fn test_band_to_time_series_invalid_band() {
        let ir = ImpulseResponse::new(0.001, vec![[0.0; NUM_OCTAVE_BANDS]; 10], vec![0.0; 10]);
        assert!(ir.band_to_time_series(NUM_OCTAVE_BANDS, 44100.0).is_empty());
    }

    #[test]
    fn test_to_time_series() {
        let mut broadband = vec![0.0; 10];
        broadband[0] = 1.0;
        let ir = ImpulseResponse::new(0.001, vec![[0.0; NUM_OCTAVE_BANDS]; 10], broadband);

        let ts = ir.to_time_series(44100.0);
        assert!(!ts.is_empty());
        // Total energy in time series should match the IR energy (up to fp error)
        let total_ts: f64 = ts.iter().sum();
        assert!((total_ts - ir.total_energy()).abs() < 1e-10);
    }
}
