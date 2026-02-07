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

/// Converts pressure-squared values to pressure and normalizes to [-1, 1].
fn energy_to_normalized_pressure(energy: &[f64]) -> Vec<f64> {
    let pressure: Vec<f64> = energy.iter().map(|&e| e.max(0.0).sqrt()).collect();
    let peak = pressure.iter().cloned().fold(0.0_f64, f64::max);
    if peak > 0.0 {
        pressure.iter().map(|&p| p / peak).collect()
    } else {
        pressure
    }
}

/// Exports a broadband impulse response as a WAV file.
///
/// The IR energy values are converted to pressure (sqrt) and normalized to [-1, 1].
pub fn write_ir_wav<P: AsRef<Path>>(
    path: P,
    ir: &ImpulseResponse,
    sample_rate: u32,
) -> Result<()> {
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
    let (dry_samples, _) = read_wav(dry_audio_path).context("Failed to read dry audio")?;
    let ir_series = ir.to_time_series(sample_rate as f64);
    let ir_pressure = energy_to_normalized_pressure(&ir_series);

    let convolved = convolve(&dry_samples, &ir_pressure);
    let peak = convolved
        .iter()
        .cloned()
        .fold(0.0_f64, |a, b| a.max(b.abs()));
    let normalized = if peak > 0.0 {
        convolved.iter().map(|&s| s / peak).collect()
    } else {
        convolved
    };

    write_wav(output_path, &normalized, sample_rate).context("Failed to write auralized WAV")
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
    let (dry_samples, _) = read_wav(dry_audio_path).context("Failed to read dry audio")?;
    let filtered_bands = filter_octave_bands(&dry_samples, sample_rate as f64);

    // Find maximum output length across all bands
    let mut max_len = 0;
    let mut band_results: Vec<Vec<f64>> = Vec::with_capacity(NUM_OCTAVE_BANDS);

    for (band, filtered) in filtered_bands.iter().enumerate() {
        let ir_series = ir.band_to_time_series(band, sample_rate as f64);
        let ir_pressure = energy_to_normalized_pressure(&ir_series);
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

    // Normalize
    let peak = output
        .iter()
        .cloned()
        .fold(0.0_f64, |a, b| a.max(b.abs()));
    if peak > 0.0 {
        for s in &mut output {
            *s /= peak;
        }
    }

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
        let ir = ImpulseResponse::new(
            0.001,
            vec![[0.0; NUM_OCTAVE_BANDS]; 100],
            broadband,
        );

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
        let dry: Vec<f64> = (0..4410)
            .map(|i| if i == 0 { 1.0 } else { 0.0 })
            .collect();
        wav::write_wav(&dry_path, &dry, 44100).unwrap();

        // Create an IR with some energy
        let mut broadband = vec![0.0; 50];
        broadband[0] = 1.0;
        broadband[5] = 0.3;
        let ir = ImpulseResponse::new(
            0.001,
            vec![[0.0; NUM_OCTAVE_BANDS]; 50],
            broadband,
        );

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
    fn test_write_ir_wav_empty() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("empty_ir.wav");

        let ir = ImpulseResponse::new(0.001, vec![], vec![]);
        write_ir_wav(&path, &ir, 44100).unwrap();

        let (samples, _) = wav::read_wav(&path).unwrap();
        assert!(samples.is_empty());
    }
}
