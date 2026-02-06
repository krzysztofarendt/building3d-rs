use crate::sim::materials::{NUM_OCTAVE_BANDS, OCTAVE_BAND_FREQUENCIES};

use super::impulse_response::ImpulseResponse;

/// Reverberation time from Schroeder backward integration.
///
/// Fits a line to the decay curve between `start_db` and `end_db`,
/// then extrapolates to -60 dB.
fn rt_from_decay(decay: &[f64], time_resolution: f64, start_db: f64, end_db: f64) -> Option<f64> {
    // Find indices where decay crosses start_db and end_db
    let mut i_start = None;
    let mut i_end = None;

    for (i, &val) in decay.iter().enumerate() {
        if val <= start_db && i_start.is_none() {
            i_start = Some(i);
        }
        if val <= end_db && i_end.is_none() {
            i_end = Some(i);
        }
    }

    let i_start = i_start?;
    let i_end = i_end?;

    if i_end <= i_start {
        return None;
    }

    // Linear regression on the decay curve between start and end
    let n = (i_end - i_start + 1) as f64;
    let mut sum_x = 0.0;
    let mut sum_y = 0.0;
    let mut sum_xy = 0.0;
    let mut sum_xx = 0.0;

    for (i, &y) in decay.iter().enumerate().take(i_end + 1).skip(i_start) {
        let x = i as f64 * time_resolution;
        sum_x += x;
        sum_y += y;
        sum_xy += x * y;
        sum_xx += x * x;
    }

    let slope = (n * sum_xy - sum_x * sum_y) / (n * sum_xx - sum_x * sum_x);

    if slope >= 0.0 {
        return None; // Decay should have negative slope
    }

    // Extrapolate to -60 dB
    Some(-60.0 / slope)
}

/// T20: Reverberation time estimated from -5 dB to -25 dB range.
pub fn t20(ir: &ImpulseResponse, band: usize) -> Option<f64> {
    let decay = ir.schroeder_decay(band);
    rt_from_decay(&decay, ir.time_resolution, -5.0, -25.0)
}

/// T30: Reverberation time estimated from -5 dB to -35 dB range.
pub fn t30(ir: &ImpulseResponse, band: usize) -> Option<f64> {
    let decay = ir.schroeder_decay(band);
    rt_from_decay(&decay, ir.time_resolution, -5.0, -35.0)
}

/// RT60: Reverberation time (uses T30 by default, falls back to T20).
pub fn rt60(ir: &ImpulseResponse, band: usize) -> Option<f64> {
    t30(ir, band).or_else(|| t20(ir, band))
}

/// EDT: Early Decay Time, estimated from 0 dB to -10 dB range.
pub fn edt(ir: &ImpulseResponse, band: usize) -> Option<f64> {
    let decay = ir.schroeder_decay(band);
    rt_from_decay(&decay, ir.time_resolution, 0.0, -10.0)
}

/// C80: Clarity (80 ms) in dB.
///
/// C80 = 10 * log10(E_early / E_late) where early = 0-80ms, late = 80ms+.
pub fn c80(ir: &ImpulseResponse, band: usize) -> Option<f64> {
    let boundary_bin = (0.080 / ir.time_resolution).round() as usize;
    if boundary_bin >= ir.len() {
        return None;
    }

    let early: f64 = ir.bands[..boundary_bin].iter().map(|b| b[band]).sum();
    let late: f64 = ir.bands[boundary_bin..].iter().map(|b| b[band]).sum();

    if late <= 0.0 {
        return None;
    }

    Some(10.0 * (early / late).log10())
}

/// D50: Definition (50 ms), ratio of early to total energy.
///
/// D50 = E_early / E_total where early = 0-50ms.
pub fn d50(ir: &ImpulseResponse, band: usize) -> Option<f64> {
    let boundary_bin = (0.050 / ir.time_resolution).round() as usize;
    if boundary_bin >= ir.len() {
        return None;
    }

    let early: f64 = ir.bands[..boundary_bin].iter().map(|b| b[band]).sum();
    let total: f64 = ir.bands.iter().map(|b| b[band]).sum();

    if total <= 0.0 {
        return None;
    }

    Some(early / total)
}

/// Simplified Speech Transmission Index estimation.
///
/// Based on modulation transfer function derived from the impulse response.
/// Uses octave bands 125 Hz through 4 kHz with standard weightings.
pub fn sti(ir: &ImpulseResponse) -> Option<f64> {
    // STI band weights (simplified, based on male speech)
    let weights = [0.129, 0.143, 0.114, 0.114, 0.186, 0.171];

    // Modulation frequencies for STI (14 per band in full STI, simplified here to 1)
    // We estimate apparent SNR from T/60 relationship
    let mut weighted_snr = 0.0;
    let mut total_weight = 0.0;

    for (band, &w) in weights.iter().enumerate() {
        let t = rt60(ir, band).unwrap_or(0.0);
        if t <= 0.0 {
            continue;
        }

        // Simplified modulation transfer function
        // m(f_mod) = 1 / sqrt(1 + (2*pi*f_mod*T/13.8)^2)
        // Using f_mod = 1 Hz as representative
        let f_mod = 1.0;
        let m = 1.0 / (1.0 + (2.0 * std::f64::consts::PI * f_mod * t / 13.8).powi(2)).sqrt();

        // Convert to apparent SNR
        let snr = 10.0 * (m / (1.0 - m).max(1e-10)).log10();
        let snr_clipped = snr.clamp(-15.0, 15.0);

        weighted_snr += w * snr_clipped;
        total_weight += w;
    }

    if total_weight <= 0.0 {
        return None;
    }

    let mean_snr = weighted_snr / total_weight;
    // Convert to STI (linear mapping from [-15, 15] dB to [0, 1])
    let sti_value = (mean_snr + 15.0) / 30.0;
    Some(sti_value.clamp(0.0, 1.0))
}

/// Aggregated room acoustic metrics for all frequency bands.
pub struct RoomAcousticReport {
    /// RT60 per band (seconds). None if not computable.
    pub rt60: [Option<f64>; NUM_OCTAVE_BANDS],
    /// EDT per band (seconds).
    pub edt: [Option<f64>; NUM_OCTAVE_BANDS],
    /// C80 per band (dB).
    pub c80: [Option<f64>; NUM_OCTAVE_BANDS],
    /// D50 per band (ratio 0-1).
    pub d50: [Option<f64>; NUM_OCTAVE_BANDS],
    /// Speech Transmission Index (0-1).
    pub sti: Option<f64>,
    /// Band center frequencies.
    pub frequencies: [f64; NUM_OCTAVE_BANDS],
}

impl RoomAcousticReport {
    /// Computes all metrics from an impulse response.
    pub fn from_ir(ir: &ImpulseResponse) -> Self {
        let mut report = Self {
            rt60: [None; NUM_OCTAVE_BANDS],
            edt: [None; NUM_OCTAVE_BANDS],
            c80: [None; NUM_OCTAVE_BANDS],
            d50: [None; NUM_OCTAVE_BANDS],
            sti: sti(ir),
            frequencies: OCTAVE_BAND_FREQUENCIES,
        };

        for band in 0..NUM_OCTAVE_BANDS {
            report.rt60[band] = rt60(ir, band);
            report.edt[band] = edt(ir, band);
            report.c80[band] = c80(ir, band);
            report.d50[band] = d50(ir, band);
        }

        report
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_test_ir() -> ImpulseResponse {
        // Create a decaying impulse response
        let num_bins = 1000;
        let dt = 0.001; // 1 ms bins
        let mut bands = vec![[0.0; NUM_OCTAVE_BANDS]; num_bins];
        let mut broadband = vec![0.0; num_bins];

        // Exponential decay with different rates per band
        for (i, (bb, band)) in broadband.iter_mut().zip(bands.iter_mut()).enumerate() {
            let t = i as f64 * dt;
            for b in 0..NUM_OCTAVE_BANDS {
                // Faster decay at higher frequencies
                let decay_rate = 10.0 + b as f64 * 2.0;
                let energy = (-decay_rate * t).exp() * 0.01;
                band[b] = energy;
                *bb += energy;
            }
        }

        ImpulseResponse::new(dt, bands, broadband)
    }

    #[test]
    fn test_c80() {
        let ir = make_test_ir();
        for band in 0..NUM_OCTAVE_BANDS {
            let c = c80(&ir, band);
            assert!(c.is_some(), "C80 should be computable for band {}", band);
            // Early energy should be higher than late, so C80 > 0
            assert!(c.unwrap() > 0.0, "C80 should be positive for decaying IR");
        }
    }

    #[test]
    fn test_d50() {
        let ir = make_test_ir();
        for band in 0..NUM_OCTAVE_BANDS {
            let d = d50(&ir, band);
            assert!(d.is_some(), "D50 should be computable for band {}", band);
            let val = d.unwrap();
            assert!(val > 0.0 && val <= 1.0, "D50 should be between 0 and 1");
        }
    }

    #[test]
    fn test_rt60() {
        let ir = make_test_ir();
        // At least some bands should have computable RT60
        let mut any_computed = false;
        for band in 0..NUM_OCTAVE_BANDS {
            if let Some(rt) = rt60(&ir, band) {
                assert!(rt > 0.0, "RT60 should be positive");
                any_computed = true;
            }
        }
        assert!(
            any_computed,
            "At least one band should have computable RT60"
        );
    }

    #[test]
    fn test_edt() {
        let ir = make_test_ir();
        let mut any_computed = false;
        for band in 0..NUM_OCTAVE_BANDS {
            if let Some(e) = edt(&ir, band) {
                assert!(e > 0.0, "EDT should be positive");
                any_computed = true;
            }
        }
        assert!(any_computed, "At least one band should have computable EDT");
    }

    #[test]
    fn test_room_acoustic_report() {
        let ir = make_test_ir();
        let report = RoomAcousticReport::from_ir(&ir);
        assert_eq!(report.frequencies.len(), NUM_OCTAVE_BANDS);
    }

    #[test]
    fn test_empty_ir() {
        let ir = ImpulseResponse::new(0.001, vec![[0.0; NUM_OCTAVE_BANDS]; 100], vec![0.0; 100]);
        assert!(rt60(&ir, 0).is_none());
        assert!(c80(&ir, 0).is_none());
        assert!(d50(&ir, 0).is_none());
    }
}
