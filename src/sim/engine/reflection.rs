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

/// Lambertian diffuse reflection (random hemisphere direction).
pub struct Diffuse;

impl ReflectionModel for Diffuse {
    fn reflect(&self, _incident: Vector, normal: Vector) -> Vector {
        use rand::Rng;
        let mut rng = rand::thread_rng();
        loop {
            let x: f64 = rng.gen_range(-1.0..1.0);
            let y: f64 = rng.gen_range(-1.0..1.0);
            let z: f64 = rng.gen_range(-1.0..1.0);
            let len2 = x * x + y * y + z * z;
            if len2 > 1e-6 && len2 <= 1.0 {
                let len = len2.sqrt();
                let v = Vector::new(x / len, y / len, z / len);
                // Ensure the reflected direction is in the same hemisphere as the normal
                if v.dot(&normal) > 0.0 {
                    return v;
                } else {
                    return v * -1.0;
                }
            }
        }
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
