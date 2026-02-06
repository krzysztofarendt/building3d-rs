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
}
