use crate::{Point, Vector};

/// Defines how rays propagate through space.
pub trait PropagationModel {
    /// Advance a ray position given its current position, velocity, and time step.
    /// Returns the new position.
    fn advance(&self, position: Point, velocity: Vector, dt: f64) -> Point;

    /// Returns the maximum distance a ray can travel before needing to check for reflections.
    fn reflection_distance(&self, speed: f64, dt: f64) -> f64;
}

/// Fixed time-step propagation: position += velocity * dt.
pub struct FixedTimeStep;

impl PropagationModel for FixedTimeStep {
    fn advance(&self, position: Point, velocity: Vector, dt: f64) -> Point {
        position + velocity * dt
    }

    fn reflection_distance(&self, speed: f64, dt: f64) -> f64 {
        speed * dt * 1.5
    }
}

/// Event-driven propagation: advances to the exact intersection point.
pub struct EventDriven;

impl PropagationModel for EventDriven {
    fn advance(&self, position: Point, velocity: Vector, dt: f64) -> Point {
        // Same as fixed time step for basic advance
        position + velocity * dt
    }

    fn reflection_distance(&self, _speed: f64, _dt: f64) -> f64 {
        // Event-driven checks all distances, no threshold
        f64::INFINITY
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fixed_time_step_advance() {
        let prop = FixedTimeStep;
        let pos = Point::new(0.0, 0.0, 0.0);
        let vel = Vector::new(1.0, 0.0, 0.0);
        let new_pos = prop.advance(pos, vel, 0.5);
        assert!((new_pos.x - 0.5).abs() < 1e-10);
        assert!((new_pos.y - 0.0).abs() < 1e-10);
        assert!((new_pos.z - 0.0).abs() < 1e-10);
    }

    #[test]
    fn test_fixed_time_step_reflection_distance() {
        let prop = FixedTimeStep;
        let rd = prop.reflection_distance(343.0, 2.5e-5);
        assert!((rd - 343.0 * 2.5e-5 * 1.5).abs() < 1e-10);
    }

    #[test]
    fn test_event_driven_infinite_distance() {
        let prop = EventDriven;
        let rd = prop.reflection_distance(343.0, 2.5e-5);
        assert!(rd.is_infinite());
    }
}
