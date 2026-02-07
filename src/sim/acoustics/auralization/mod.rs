//! Auralization: convert acoustic impulse responses to audible audio.
//!
//! Provides WAV export of impulse responses and convolution with dry audio
//! for both broadband and per-band auralization.

pub mod convolve;
pub mod filters;
pub mod wav;

use std::path::Path;

use anyhow::{Context, Result};

use crate::sim::materials::NUM_OCTAVE_BANDS;

use super::impulse_response::ImpulseResponse;
use convolve::convolve;
use filters::filter_octave_bands;
use wav::{read_wav, write_wav};

fn resample_linear(signal: &[f64], from_sample_rate: u32, to_sample_rate: u32) -> Vec<f64> {
    if signal.is_empty() || from_sample_rate == 0 || to_sample_rate == 0 {
        return Vec::new();
    }
    if from_sample_rate == to_sample_rate {
        return signal.to_vec();
    }

    let out_len = ((signal.len() as f64) * (to_sample_rate as f64) / (from_sample_rate as f64))
        .round()
        .max(1.0) as usize;

    let from_sr = from_sample_rate as f64;
    let to_sr = to_sample_rate as f64;

    let mut out = Vec::with_capacity(out_len);
    for out_idx in 0..out_len {
        let src_pos = (out_idx as f64) * from_sr / to_sr;
        let src_idx = src_pos.floor() as usize;
        let frac = src_pos - src_idx as f64;

        let sample = if src_idx + 1 >= signal.len() {
            signal[signal.len() - 1]
        } else {
            (1.0 - frac) * signal[src_idx] + frac * signal[src_idx + 1]
        };
        out.push(sample);
    }
    out
}

fn normalize_to_unit_peak(samples: &[f64]) -> Vec<f64> {
    let peak = samples.iter().cloned().fold(0.0_f64, |a, b| a.max(b.abs()));
    if peak > 0.0 {
        samples.iter().map(|&s| s / peak).collect()
    } else {
        samples.to_vec()
    }
}

/// Converts pressure-squared values to pressure (sqrt), clamping negatives to 0.
fn energy_to_pressure(energy: &[f64]) -> Vec<f64> {
    energy.iter().map(|&e| e.max(0.0).sqrt()).collect()
}

/// Converts pressure-squared values to pressure and normalizes to [-1, 1].
fn energy_to_normalized_pressure(energy: &[f64]) -> Vec<f64> {
    normalize_to_unit_peak(&energy_to_pressure(energy))
}

/// Exports a broadband impulse response as a WAV file.
///
/// The IR energy values are converted to pressure (sqrt) and normalized to [-1, 1].
pub fn write_ir_wav<P: AsRef<Path>>(path: P, ir: &ImpulseResponse, sample_rate: u32) -> Result<()> {
    let time_series = ir.to_time_series(sample_rate as f64);
    let samples = energy_to_normalized_pressure(&time_series);
    write_wav(path, &samples, sample_rate).context("Failed to write IR WAV")
}

/// Convolves a broadband impulse response with dry audio and writes the result.
///
/// Reads the dry audio WAV, convolves it with the broadband IR, normalizes,
/// and writes the output WAV at the given sample rate.
pub fn auralize<P: AsRef<Path>>(
    ir: &ImpulseResponse,
    dry_audio_path: P,
    output_path: P,
    sample_rate: u32,
) -> Result<()> {
    let (dry_samples, dry_sample_rate) =
        read_wav(dry_audio_path).context("Failed to read dry audio")?;
    let dry_samples = resample_linear(&dry_samples, dry_sample_rate, sample_rate);
    let ir_series = ir.to_time_series(sample_rate as f64);
    let ir_pressure = energy_to_normalized_pressure(&ir_series);

    let convolved = convolve(&dry_samples, &ir_pressure);
    let normalized = normalize_to_unit_peak(&convolved);

    write_wav(output_path, &normalized, sample_rate).context("Failed to write auralized WAV")
}

fn auralize_per_band_from_filtered(
    ir: &ImpulseResponse,
    filtered_bands: &[Vec<f64>; NUM_OCTAVE_BANDS],
    sample_rate: f64,
) -> Vec<f64> {
    // Find maximum output length across all bands
    let mut max_len = 0;
    let mut band_results: Vec<Vec<f64>> = Vec::with_capacity(NUM_OCTAVE_BANDS);

    for (band, filtered) in filtered_bands.iter().enumerate() {
        let ir_series = ir.band_to_time_series(band, sample_rate);
        let ir_pressure = energy_to_pressure(&ir_series);
        let convolved = convolve(filtered, &ir_pressure);
        if convolved.len() > max_len {
            max_len = convolved.len();
        }
        band_results.push(convolved);
    }

    // Sum all bands
    let mut output = vec![0.0; max_len];
    for band in &band_results {
        for (i, &s) in band.iter().enumerate() {
            output[i] += s;
        }
    }

    normalize_to_unit_peak(&output)
}

/// Convolves per-band impulse responses with filtered dry audio and writes the result.
///
/// The dry audio is filtered into 6 octave bands. Each band is convolved with the
/// corresponding band's IR. The 6 convolved signals are summed and normalized.
pub fn auralize_per_band<P: AsRef<Path>>(
    ir: &ImpulseResponse,
    dry_audio_path: P,
    output_path: P,
    sample_rate: u32,
) -> Result<()> {
    let (dry_samples, dry_sample_rate) =
        read_wav(dry_audio_path).context("Failed to read dry audio")?;
    let dry_samples = resample_linear(&dry_samples, dry_sample_rate, sample_rate);
    let filtered_bands = filter_octave_bands(&dry_samples, sample_rate as f64);
    let output = auralize_per_band_from_filtered(ir, &filtered_bands, sample_rate as f64);
    write_wav(output_path, &output, sample_rate).context("Failed to write auralized WAV")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sim::materials::NUM_OCTAVE_BANDS;

    #[test]
    fn test_energy_to_normalized_pressure() {
        let energy = vec![0.0, 0.25, 1.0, 0.0];
        let pressure = energy_to_normalized_pressure(&energy);
        assert_eq!(pressure.len(), 4);
        assert!((pressure[0] - 0.0).abs() < 1e-10);
        assert!((pressure[2] - 1.0).abs() < 1e-10); // peak normalized
        assert!((pressure[1] - 0.5).abs() < 1e-10); // sqrt(0.25) / sqrt(1.0) = 0.5
    }

    #[test]
    fn test_energy_to_normalized_pressure_all_zero() {
        let energy = vec![0.0, 0.0, 0.0];
        let pressure = energy_to_normalized_pressure(&energy);
        assert!(pressure.iter().all(|&p| p == 0.0));
    }

    #[test]
    fn test_write_ir_wav() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("ir.wav");

        let mut broadband = vec![0.0; 100];
        broadband[0] = 1.0;
        broadband[10] = 0.5;
        let ir = ImpulseResponse::new(0.001, vec![[0.0; NUM_OCTAVE_BANDS]; 100], broadband);

        write_ir_wav(&path, &ir, 44100).unwrap();

        // Verify the file can be read back
        let (samples, sr) = wav::read_wav(&path).unwrap();
        assert_eq!(sr, 44100);
        assert!(!samples.is_empty());

        // Peak should be normalized to ~1.0
        let peak = samples.iter().cloned().fold(0.0_f64, |a, b| a.max(b.abs()));
        assert!(peak > 0.9, "Peak should be close to 1.0, got {peak}");
    }

    #[test]
    fn test_auralize_broadband() {
        let dir = tempfile::tempdir().unwrap();
        let dry_path = dir.path().join("dry.wav");
        let out_path = dir.path().join("wet.wav");

        // Create a short dry signal (click)
        let dry: Vec<f64> = (0..4410).map(|i| if i == 0 { 1.0 } else { 0.0 }).collect();
        wav::write_wav(&dry_path, &dry, 44100).unwrap();

        // Create an IR with some energy
        let mut broadband = vec![0.0; 50];
        broadband[0] = 1.0;
        broadband[5] = 0.3;
        let ir = ImpulseResponse::new(0.001, vec![[0.0; NUM_OCTAVE_BANDS]; 50], broadband);

        auralize(&ir, &dry_path, &out_path, 44100).unwrap();

        let (samples, sr) = wav::read_wav(&out_path).unwrap();
        assert_eq!(sr, 44100);
        assert!(!samples.is_empty());
    }

    #[test]
    fn test_auralize_per_band() {
        let dir = tempfile::tempdir().unwrap();
        let dry_path = dir.path().join("dry.wav");
        let out_path = dir.path().join("wet_bands.wav");

        // Create a short dry signal
        let dry: Vec<f64> = (0..4410)
            .map(|i| (2.0 * std::f64::consts::PI * 500.0 * i as f64 / 44100.0).sin())
            .collect();
        wav::write_wav(&dry_path, &dry, 44100).unwrap();

        // Create an IR with per-band energy
        let mut bands = vec![[0.0; NUM_OCTAVE_BANDS]; 50];
        bands[0] = [0.5; NUM_OCTAVE_BANDS];
        let broadband: Vec<f64> = bands.iter().map(|b| b.iter().sum()).collect();
        let ir = ImpulseResponse::new(0.001, bands, broadband);

        auralize_per_band(&ir, &dry_path, &out_path, 44100).unwrap();

        let (samples, sr) = wav::read_wav(&out_path).unwrap();
        assert_eq!(sr, 44100);
        assert!(!samples.is_empty());
    }

    #[test]
    fn test_auralize_resamples_dry_audio_to_requested_sample_rate() {
        let dir = tempfile::tempdir().unwrap();
        let dry_path = dir.path().join("dry_48k.wav");
        let out_path = dir.path().join("wet_44k.wav");

        // 0.1 seconds of silence with an initial click @ 48k
        let dry_sample_rate = 48_000u32;
        let dry_len = (dry_sample_rate as f64 * 0.1) as usize;
        let mut dry = vec![0.0; dry_len];
        if !dry.is_empty() {
            dry[0] = 1.0;
        }
        wav::write_wav(&dry_path, &dry, dry_sample_rate).unwrap();

        // IR = 50 ms
        let mut broadband = vec![0.0; 50];
        broadband[0] = 1.0;
        let ir = ImpulseResponse::new(0.001, vec![[0.0; NUM_OCTAVE_BANDS]; 50], broadband);

        let out_sample_rate = 44_100u32;
        auralize(&ir, &dry_path, &out_path, out_sample_rate).unwrap();

        let (samples, sr) = wav::read_wav(&out_path).unwrap();
        assert_eq!(sr, out_sample_rate);

        let expected_dry_len = ((dry_len as f64) * (out_sample_rate as f64)
            / (dry_sample_rate as f64))
            .round() as usize;
        let expected_ir_len = (0.05 * out_sample_rate as f64).ceil() as usize;
        let expected_len = expected_dry_len + expected_ir_len - 1;
        assert_eq!(samples.len(), expected_len);
    }

    #[test]
    fn test_auralize_per_band_preserves_relative_band_levels() {
        // Construct filtered bands directly (avoid filter design behavior).
        let filtered_bands: [Vec<f64>; NUM_OCTAVE_BANDS] = std::array::from_fn(|_| vec![1.0]);

        let sample_rate = 1000.0;
        let time_resolution = 1.0 / sample_rate;

        let mut bands = vec![[0.0; NUM_OCTAVE_BANDS]; 20];
        bands[0][0] = 1.0; // strong band at t=0
        bands[10][1] = 0.01; // weak band at t=10 samples

        let broadband: Vec<f64> = bands.iter().map(|b| b.iter().sum()).collect();
        let ir = ImpulseResponse::new(time_resolution, bands, broadband);

        let output = auralize_per_band_from_filtered(&ir, &filtered_bands, sample_rate);
        assert!(output.len() > 10);

        // Global normalization keeps the largest peak at 1.0; the weaker band peak remains ~0.1.
        assert!((output[0] - 1.0).abs() < 1e-12);
        assert!((output[10] - 0.1).abs() < 1e-12, "got {}", output[10]);
    }

    #[test]
    fn test_write_ir_wav_empty() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("empty_ir.wav");

        let ir = ImpulseResponse::new(0.001, vec![], vec![]);
        write_ir_wav(&path, &ir, 44100).unwrap();

        let (samples, _) = wav::read_wav(&path).unwrap();
        assert!(samples.is_empty());
    }
}
