use rustfft::FftPlanner;
use rustfft::num_complex::Complex;

/// FFT overlap-add convolution.
///
/// Returns a vector of length `signal.len() + kernel.len() - 1`.
/// Uses overlap-add with FFT for efficient O(N log N) computation.
pub fn convolve(signal: &[f64], kernel: &[f64]) -> Vec<f64> {
    if signal.is_empty() || kernel.is_empty() {
        return Vec::new();
    }

    let output_len = signal.len() + kernel.len() - 1;

    // For very short signals, use direct convolution
    if signal.len() <= 64 || kernel.len() <= 64 {
        return convolve_direct(signal, kernel);
    }

    // Block size: next power of 2 >= 2 * kernel.len()
    let fft_size = (2 * kernel.len()).next_power_of_two();
    let block_size = fft_size - kernel.len() + 1;

    let mut planner = FftPlanner::<f64>::new();
    let fft = planner.plan_fft_forward(fft_size);
    let ifft = planner.plan_fft_inverse(fft_size);

    // Pre-compute kernel FFT (zero-padded to fft_size)
    let mut kernel_fft: Vec<Complex<f64>> = kernel
        .iter()
        .map(|&x| Complex::new(x, 0.0))
        .chain(std::iter::repeat_n(
            Complex::new(0.0, 0.0),
            fft_size - kernel.len(),
        ))
        .collect();
    fft.process(&mut kernel_fft);

    let mut output = vec![0.0; output_len];
    let scale = 1.0 / fft_size as f64;

    // Process signal in blocks
    let mut pos = 0;
    while pos < signal.len() {
        let end = (pos + block_size).min(signal.len());
        let chunk_len = end - pos;

        // Zero-padded block
        let mut block: Vec<Complex<f64>> = signal[pos..end]
            .iter()
            .map(|&x| Complex::new(x, 0.0))
            .chain(std::iter::repeat_n(
                Complex::new(0.0, 0.0),
                fft_size - chunk_len,
            ))
            .collect();

        // Forward FFT
        fft.process(&mut block);

        // Pointwise multiply
        for (b, k) in block.iter_mut().zip(kernel_fft.iter()) {
            *b *= k;
        }

        // Inverse FFT
        ifft.process(&mut block);

        // Overlap-add
        let valid_len = (chunk_len + kernel.len() - 1).min(output_len - pos);
        for i in 0..valid_len {
            output[pos + i] += block[i].re * scale;
        }

        pos += block_size;
    }

    output
}

/// Direct convolution for short signals.
fn convolve_direct(signal: &[f64], kernel: &[f64]) -> Vec<f64> {
    let output_len = signal.len() + kernel.len() - 1;
    let mut output = vec![0.0; output_len];
    for (i, &s) in signal.iter().enumerate() {
        for (j, &k) in kernel.iter().enumerate() {
            output[i + j] += s * k;
        }
    }
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_convolve_unit_impulse() {
        let signal = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let kernel = vec![1.0]; // unit impulse
        let result = convolve(&signal, &kernel);
        assert_eq!(result.len(), signal.len());
        for (a, b) in signal.iter().zip(result.iter()) {
            assert!((a - b).abs() < 1e-10, "Mismatch: {a} vs {b}");
        }
    }

    #[test]
    fn test_convolve_delayed_impulse() {
        let signal = vec![1.0, 2.0, 3.0];
        let kernel = vec![0.0, 1.0]; // delay by 1
        let result = convolve(&signal, &kernel);
        assert_eq!(result.len(), 4);
        assert!((result[0] - 0.0).abs() < 1e-10);
        assert!((result[1] - 1.0).abs() < 1e-10);
        assert!((result[2] - 2.0).abs() < 1e-10);
        assert!((result[3] - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_convolve_known_result() {
        let signal = vec![1.0, 2.0, 3.0];
        let kernel = vec![1.0, 1.0];
        // Expected: [1, 3, 5, 3]
        let result = convolve(&signal, &kernel);
        let expected = [1.0, 3.0, 5.0, 3.0];
        assert_eq!(result.len(), expected.len());
        for (a, b) in result.iter().zip(expected.iter()) {
            assert!((a - b).abs() < 1e-10, "Mismatch: {a} vs {b}");
        }
    }

    #[test]
    fn test_convolve_empty() {
        assert!(convolve(&[], &[1.0]).is_empty());
        assert!(convolve(&[1.0], &[]).is_empty());
        assert!(convolve(&[], &[]).is_empty());
    }

    #[test]
    fn test_convolve_long_signal() {
        // Test with a long enough signal to trigger FFT path
        let signal: Vec<f64> = (0..1000).map(|i| (i as f64 * 0.01).sin()).collect();
        let kernel: Vec<f64> = (0..200).map(|i| if i == 0 { 1.0 } else { 0.0 }).collect();

        let result = convolve(&signal, &kernel);
        assert_eq!(result.len(), signal.len() + kernel.len() - 1);

        // Convolving with unit impulse (+ zeros) should return original signal padded
        for (i, &s) in signal.iter().enumerate() {
            assert!(
                (result[i] - s).abs() < 1e-10,
                "Mismatch at {i}: {} vs {s}",
                result[i]
            );
        }
    }

    #[test]
    fn test_convolve_fft_vs_direct() {
        // Compare FFT path with direct for a moderately sized input
        let signal: Vec<f64> = (0..500).map(|i| (i as f64 * 0.05).sin()).collect();
        let kernel: Vec<f64> = (0..100).map(|i| (-0.01 * i as f64).exp()).collect();

        let fft_result = convolve(&signal, &kernel);
        let direct_result = convolve_direct(&signal, &kernel);

        assert_eq!(fft_result.len(), direct_result.len());
        for (i, (a, b)) in fft_result.iter().zip(direct_result.iter()).enumerate() {
            assert!(
                (a - b).abs() < 1e-8,
                "Mismatch at {i}: fft={a} vs direct={b}"
            );
        }
    }
}
