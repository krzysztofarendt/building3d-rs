use std::path::Path;

use anyhow::{Context, Result};

/// Reads a WAV file and returns mono samples normalized to [-1, 1] and the sample rate.
///
/// Stereo files are downmixed to mono by averaging channels.
pub fn read_wav<P: AsRef<Path>>(path: P) -> Result<(Vec<f64>, u32)> {
    let reader =
        hound::WavReader::open(path.as_ref()).context("Failed to open WAV file for reading")?;
    let spec = reader.spec();
    let sample_rate = spec.sample_rate;
    let channels = spec.channels as usize;

    let raw_samples: Vec<f64> = match spec.sample_format {
        hound::SampleFormat::Int => {
            let max_val = (1i64 << (spec.bits_per_sample - 1)) as f64;
            reader
                .into_samples::<i32>()
                .map(|s| s.map(|v| v as f64 / max_val))
                .collect::<Result<Vec<f64>, _>>()
                .context("Failed to read integer samples")?
        }
        hound::SampleFormat::Float => reader
            .into_samples::<f32>()
            .map(|s| s.map(|v| v as f64))
            .collect::<Result<Vec<f64>, _>>()
            .context("Failed to read float samples")?,
    };

    let mono = if channels == 1 {
        raw_samples
    } else {
        raw_samples
            .chunks(channels)
            .map(|frame| frame.iter().sum::<f64>() / channels as f64)
            .collect()
    };

    Ok((mono, sample_rate))
}

/// Writes mono samples to a 16-bit PCM WAV file.
///
/// Samples are clamped to [-1, 1] before quantization.
pub fn write_wav<P: AsRef<Path>>(path: P, samples: &[f64], sample_rate: u32) -> Result<()> {
    let spec = hound::WavSpec {
        channels: 1,
        sample_rate,
        bits_per_sample: 16,
        sample_format: hound::SampleFormat::Int,
    };
    let mut writer =
        hound::WavWriter::create(path.as_ref(), spec).context("Failed to create WAV file")?;

    for &s in samples {
        let clamped = s.clamp(-1.0, 1.0);
        let quantized = (clamped * i16::MAX as f64) as i16;
        writer
            .write_sample(quantized)
            .context("Failed to write sample")?;
    }

    writer.finalize().context("Failed to finalize WAV file")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_wav_round_trip() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.wav");

        let original = vec![0.0, 0.5, -0.5, 1.0, -1.0, 0.25];
        write_wav(&path, &original, 44100).unwrap();

        let (read_back, sr) = read_wav(&path).unwrap();
        assert_eq!(sr, 44100);
        assert_eq!(read_back.len(), original.len());

        // 16-bit quantization error: 1 / 32767 ≈ 3.05e-5
        for (a, b) in original.iter().zip(read_back.iter()) {
            assert!(
                (a - b).abs() < 1e-4,
                "Mismatch: original={a}, read_back={b}"
            );
        }
    }

    #[test]
    fn test_wav_clamping() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("clamp.wav");

        let samples = vec![2.0, -3.0];
        write_wav(&path, &samples, 44100).unwrap();

        let (read_back, _) = read_wav(&path).unwrap();
        // Values should be clamped to [-1, 1]
        assert!((read_back[0] - 1.0).abs() < 1e-4);
        assert!((read_back[1] - (-1.0)).abs() < 1e-4);
    }

    #[test]
    fn test_wav_empty() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("empty.wav");

        write_wav(&path, &[], 44100).unwrap();
        let (read_back, sr) = read_wav(&path).unwrap();
        assert_eq!(sr, 44100);
        assert!(read_back.is_empty());
    }
}
