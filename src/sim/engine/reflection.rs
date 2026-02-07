use crate::Vector;

/// Defines how rays reflect off surfaces.
pub trait ReflectionModel {
    /// Computes the reflected direction given incident direction and surface normal.
    fn reflect(&self, incident: Vector, normal: Vector) -> Vector;
}

/// Perfect specular (mirror) reflection.
pub struct Specular;

impl ReflectionModel for Specular {
    fn reflect(&self, incident: Vector, normal: Vector) -> Vector {
        let dot = incident.dot(&normal);
        incident - 2.0 * dot * normal
    }
}

/// Lambertian diffuse reflection (cosine-weighted hemisphere sampling via Malley's method).
pub struct Diffuse;

impl ReflectionModel for Diffuse {
    fn reflect(&self, incident: Vector, normal: Vector) -> Vector {
        use rand::Rng;
        // Flip the hemisphere so the reflected ray stays on the same side of the surface
        // as the incident ray (handles outward-facing normals for interior propagation).
        let hemisphere_normal = if incident.dot(&normal) >= 0.0 {
            normal * -1.0
        } else {
            normal
        };

        // Build orthonormal basis (tangent, bitangent) around the hemisphere normal.
        let n = hemisphere_normal;
        let arbitrary = if n.dx.abs() < 0.9 {
            Vector::new(1.0, 0.0, 0.0)
        } else {
            Vector::new(0.0, 1.0, 0.0)
        };
        let tangent = n
            .cross(&arbitrary)
            .normalize()
            .unwrap_or(Vector::new(1.0, 0.0, 0.0));
        let bitangent = n.cross(&tangent);

        // Malley's method: sample uniformly on a disk, then project onto hemisphere.
        // This produces a cosine-weighted distribution (pdf = cos(theta) / pi).
        let mut rng = rand::thread_rng();
        let u1: f64 = rng.r#gen();
        let u2: f64 = rng.r#gen();
        let r = u1.sqrt();
        let phi = 2.0 * std::f64::consts::PI * u2;
        let x = r * phi.cos();
        let y = r * phi.sin();
        let z = (1.0 - u1).sqrt(); // = sqrt(1 - r^2)

        // Transform from local to world coordinates.
        tangent * x + bitangent * y + n * z
    }
}

/// Hybrid reflection mixing specular and diffuse based on a scattering coefficient.
pub struct Hybrid {
    /// Scattering coefficient [0, 1]: 0 = pure specular, 1 = pure diffuse.
    pub scattering: f64,
}

impl Hybrid {
    pub fn new(scattering: f64) -> Self {
        Self {
            scattering: scattering.clamp(0.0, 1.0),
        }
    }
}

impl ReflectionModel for Hybrid {
    fn reflect(&self, incident: Vector, normal: Vector) -> Vector {
        use rand::Rng;
        let mut rng = rand::thread_rng();
        let r: f64 = rng.r#gen();
        if r < self.scattering {
            Diffuse.reflect(incident, normal)
        } else {
            Specular.reflect(incident, normal)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_specular_reflection() {
        let specular = Specular;
        // Ray coming straight down onto a horizontal surface (normal pointing up)
        let incident = Vector::new(0.0, 0.0, -1.0);
        let normal = Vector::new(0.0, 0.0, 1.0);
        let reflected = specular.reflect(incident, normal);
        assert!((reflected.dx - 0.0).abs() < 1e-10);
        assert!((reflected.dy - 0.0).abs() < 1e-10);
        assert!((reflected.dz - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_specular_45_degrees() {
        let specular = Specular;
        let incident = Vector::new(1.0, 0.0, -1.0);
        let normal = Vector::new(0.0, 0.0, 1.0);
        let reflected = specular.reflect(incident, normal);
        assert!((reflected.dx - 1.0).abs() < 1e-10);
        assert!((reflected.dy - 0.0).abs() < 1e-10);
        assert!((reflected.dz - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_diffuse_reflects_in_hemisphere() {
        let diffuse = Diffuse;
        let normal = Vector::new(0.0, 0.0, 1.0);
        let incident = Vector::new(0.0, 0.0, -1.0);
        for _ in 0..100 {
            let reflected = diffuse.reflect(incident, normal);
            assert!(
                reflected.dot(&normal) > 0.0,
                "Diffuse reflection should be in the same hemisphere as normal"
            );
        }
    }

    #[test]
    fn test_diffuse_cosine_weighted_distribution() {
        // Verify that Malley's method produces cosine-weighted samples:
        // mean(cos(theta)) should be 2/3 for cosine-weighted hemisphere sampling.
        let diffuse = Diffuse;
        let normal = Vector::new(0.0, 0.0, 1.0);
        let incident = Vector::new(0.0, 0.0, -1.0);
        let n = 10000;
        let mut cos_sum = 0.0;
        for _ in 0..n {
            let reflected = diffuse.reflect(incident, normal);
            cos_sum += reflected.dot(&normal);
        }
        let mean_cos = cos_sum / n as f64;
        // For cosine-weighted: E[cos(theta)] = 2/3
        assert!(
            (mean_cos - 2.0 / 3.0).abs() < 0.05,
            "Mean cos(theta) should be ~0.667 for cosine-weighted sampling, got {mean_cos}"
        );
    }

    #[test]
    fn test_diffuse_respects_incident_side() {
        let diffuse = Diffuse;
        let normal = Vector::new(0.0, 0.0, 1.0);
        // Incident coming from the same side as the normal (e.g. inside a solid with outward normal)
        let incident = Vector::new(0.0, 0.0, 1.0);
        for _ in 0..100 {
            let reflected = diffuse.reflect(incident, normal);
            assert!(
                reflected.dot(&normal) < 0.0,
                "Diffuse reflection should be in the opposite hemisphere when incident·normal > 0"
            );
        }
    }

    #[test]
    fn test_hybrid_produces_valid_directions() {
        let hybrid = Hybrid::new(0.5);
        let normal = Vector::new(0.0, 0.0, 1.0);
        let incident = Vector::new(0.0, 0.0, -1.0);
        for _ in 0..100 {
            let reflected = hybrid.reflect(incident, normal);
            let len = reflected.length();
            assert!(len > 0.0, "Reflected vector should have non-zero length");
        }
    }
}
