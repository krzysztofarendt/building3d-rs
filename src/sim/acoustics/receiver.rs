use crate::Point;
use crate::sim::materials::NUM_OCTAVE_BANDS;

/// A spherical receiver that accumulates energy-time-frequency data.
pub struct Receiver {
    /// Position of the receiver in 3D space.
    pub position: Point,
    /// Collection radius (sphere around the position).
    pub radius: f64,
    /// Time resolution in seconds (bin width).
    pub time_resolution: f64,
    /// Number of time bins.
    pub num_bins: usize,
    /// Energy-time histogram per frequency band: histogram[bin][band].
    histogram: Vec<[f64; NUM_OCTAVE_BANDS]>,
    /// Total scalar energy histogram: scalar_histogram[bin].
    scalar_histogram: Vec<f64>,
}

impl Receiver {
    /// Creates a new receiver.
    ///
    /// - `position`: center of the collection sphere
    /// - `radius`: sphere radius in meters
    /// - `time_resolution`: time bin width in seconds
    /// - `max_time`: maximum recording time in seconds
    pub fn new(position: Point, radius: f64, time_resolution: f64, max_time: f64) -> Self {
        let num_bins = (max_time / time_resolution).ceil() as usize;
        Self {
            position,
            radius,
            time_resolution,
            num_bins,
            histogram: vec![[0.0; NUM_OCTAVE_BANDS]; num_bins],
            scalar_histogram: vec![0.0; num_bins],
        }
    }

    /// Records a ray hit at a given time with frequency-band energies.
    pub fn record_band_hit(&mut self, time: f64, energy: &[f64; NUM_OCTAVE_BANDS]) {
        let bin = (time / self.time_resolution) as usize;
        if bin < self.num_bins {
            for (b, e) in self.histogram[bin].iter_mut().zip(energy.iter()) {
                *b += e;
            }
            self.scalar_histogram[bin] += energy.iter().sum::<f64>();
        }
    }

    /// Records a scalar ray hit at a given time.
    pub fn record_scalar_hit(&mut self, time: f64, energy: f64) {
        let bin = (time / self.time_resolution) as usize;
        if bin < self.num_bins {
            self.scalar_histogram[bin] += energy;
        }
    }

    /// Returns the energy-time histogram per band.
    pub fn histogram(&self) -> &[[f64; NUM_OCTAVE_BANDS]] {
        &self.histogram
    }

    /// Returns the scalar energy-time histogram.
    pub fn scalar_histogram(&self) -> &[f64] {
        &self.scalar_histogram
    }

    /// Checks if a point is within the receiver's collection sphere.
    pub fn contains(&self, point: Point) -> bool {
        let dx = point.x - self.position.x;
        let dy = point.y - self.position.y;
        let dz = point.z - self.position.z;
        let dist2 = dx * dx + dy * dy + dz * dz;
        dist2 <= self.radius * self.radius
    }

    /// Resets all recorded data.
    pub fn reset(&mut self) {
        for bin in &mut self.histogram {
            *bin = [0.0; NUM_OCTAVE_BANDS];
        }
        for bin in &mut self.scalar_histogram {
            *bin = 0.0;
        }
    }

    /// Returns the total energy recorded across all bins and bands.
    pub fn total_energy(&self) -> f64 {
        self.scalar_histogram.iter().sum()
    }

    /// Returns a normalization factor that converts raw collected energy to
    /// energy density (W/m^2) independent of receiver size and ray count.
    ///
    /// `E_normalized = E_raw * normalization_factor`
    ///
    /// where `normalization_factor = 4*pi / (pi * r^2 * N_rays)`
    ///
    /// This accounts for:
    /// - Receiver cross-sectional area (`pi * r^2`)
    /// - Total solid angle of emission (`4*pi` steradians)
    /// - Number of rays (`N_rays`)
    pub fn normalization_factor(&self, num_rays: usize) -> f64 {
        let cross_section = std::f64::consts::PI * self.radius * self.radius;
        if cross_section <= 0.0 || num_rays == 0 {
            return 1.0;
        }
        4.0 * std::f64::consts::PI / (cross_section * num_rays as f64)
    }

    /// Returns a normalized copy of the scalar histogram (energy density, W/m^2).
    pub fn normalized_scalar_histogram(&self, num_rays: usize) -> Vec<f64> {
        let factor = self.normalization_factor(num_rays);
        self.scalar_histogram.iter().map(|&e| e * factor).collect()
    }

    /// Returns a normalized copy of the band histogram (energy density per band).
    pub fn normalized_histogram(&self, num_rays: usize) -> Vec<[f64; NUM_OCTAVE_BANDS]> {
        let factor = self.normalization_factor(num_rays);
        self.histogram
            .iter()
            .map(|bin| {
                let mut out = [0.0; NUM_OCTAVE_BANDS];
                for (i, &e) in bin.iter().enumerate() {
                    out[i] = e * factor;
                }
                out
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_receiver_creation() {
        let r = Receiver::new(Point::new(0.0, 0.0, 0.0), 0.5, 0.001, 1.0);
        assert_eq!(r.num_bins, 1000);
        assert!((r.radius - 0.5).abs() < 1e-10);
    }

    #[test]
    fn test_receiver_contains() {
        let r = Receiver::new(Point::new(1.0, 1.0, 1.0), 0.5, 0.001, 1.0);
        assert!(r.contains(Point::new(1.0, 1.0, 1.0)));
        assert!(r.contains(Point::new(1.3, 1.0, 1.0)));
        assert!(!r.contains(Point::new(2.0, 1.0, 1.0)));
    }

    #[test]
    fn test_record_scalar_hit() {
        let mut r = Receiver::new(Point::new(0.0, 0.0, 0.0), 0.5, 0.001, 1.0);
        r.record_scalar_hit(0.005, 0.8);
        r.record_scalar_hit(0.005, 0.2);
        // bin 5
        assert!((r.scalar_histogram()[5] - 1.0).abs() < 1e-10);
        assert!((r.total_energy() - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_record_band_hit() {
        let mut r = Receiver::new(Point::new(0.0, 0.0, 0.0), 0.5, 0.001, 1.0);
        let energy = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6];
        r.record_band_hit(0.010, &energy);
        // bin 10
        assert!((r.histogram()[10][0] - 0.1).abs() < 1e-10);
        assert!((r.histogram()[10][5] - 0.6).abs() < 1e-10);
        let expected_total: f64 = energy.iter().sum();
        assert!((r.total_energy() - expected_total).abs() < 1e-10);
    }

    #[test]
    fn test_receiver_reset() {
        let mut r = Receiver::new(Point::new(0.0, 0.0, 0.0), 0.5, 0.001, 1.0);
        r.record_scalar_hit(0.005, 1.0);
        assert!(r.total_energy() > 0.0);
        r.reset();
        assert!((r.total_energy() - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_out_of_range_hit_ignored() {
        let mut r = Receiver::new(Point::new(0.0, 0.0, 0.0), 0.5, 0.001, 1.0);
        r.record_scalar_hit(2.0, 1.0); // beyond max_time
        assert!((r.total_energy() - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_normalization_factor() {
        let r = Receiver::new(Point::new(0.0, 0.0, 0.0), 0.5, 0.001, 1.0);
        let factor = r.normalization_factor(1000);
        // factor = 4*pi / (pi * 0.25 * 1000) = 4 / 250 = 0.016
        let expected = 4.0 * std::f64::consts::PI / (std::f64::consts::PI * 0.25 * 1000.0);
        assert!((factor - expected).abs() < 1e-10);
    }

    #[test]
    fn test_normalized_histogram() {
        let mut r = Receiver::new(Point::new(0.0, 0.0, 0.0), 0.5, 0.001, 1.0);
        r.record_scalar_hit(0.005, 1.0);
        let normalized = r.normalized_scalar_histogram(1000);
        let factor = r.normalization_factor(1000);
        assert!((normalized[5] - 1.0 * factor).abs() < 1e-10);
    }

    #[test]
    fn test_normalization_scales_inversely_with_radius() {
        // Larger receiver should produce smaller normalization factor
        let r1 = Receiver::new(Point::new(0.0, 0.0, 0.0), 0.1, 0.001, 1.0);
        let r2 = Receiver::new(Point::new(0.0, 0.0, 0.0), 1.0, 0.001, 1.0);
        let f1 = r1.normalization_factor(1000);
        let f2 = r2.normalization_factor(1000);
        // f1/f2 = r2^2/r1^2 = 100
        assert!((f1 / f2 - 100.0).abs() < 1e-6);
    }
}
