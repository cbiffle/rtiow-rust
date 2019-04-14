use rand::prelude::*;

use crate::object::HitRecord;
use crate::ray::Ray;
use crate::texture::Texture;
use crate::vec3::{reflect, refract, Vec3};

/// Material options for a rendered object.
#[derive(Clone)]
pub enum Material {
    /// An opaque material with a matte surface, where lighting is calculated
    /// using [Lambertian reflectance][lambert].
    ///
    /// [lambert]: https://en.wikipedia.org/wiki/Lambertian_reflectance
    Lambertian { albedo: Texture },
    /// A reflective material that looks like polished or frosted metal.
    Metal {
        /// The amount of light energy reflected in each color component, so
        /// `Vec3(1., 1., 1.)` is a white surface, and `Vec3(0., 0., 0.)` is
        /// totally black.
        albedo: Vec3,
        /// The amount of randomness introduced into reflected rays. A `fuzz` of
        /// 0 makes the surface look polished and mirror-smooth, while a `fuzz`
        /// of 1 produces a frosted, almost matte surface.
        fuzz: f32,
    },
    /// A transparent refractive material like glass or water.
    Dielectric {
        /// [Refractive index][ref-idx] of the material, which determines how
        /// much light is bent when traveling into or out of an object.
        ///
        /// [ref-idx]: https://en.wikipedia.org/wiki/Refractive_index
        ref_idx: f32,
    },
    /// Diffuse light.
    DiffuseLight { emission: Texture, brightness: f32 },
    /// Isotropic scattering.
    Isotropic { albedo: Texture },
}

impl Material {
    /// Performs surface scattering from a material.
    ///
    /// When light traveling along `ray` reaches a surface made out of this
    /// material (intersection described by `hit`), some of it will be absorbed,
    /// and the rest will either be reflected or refracted. If 100% of the light
    /// is absorbed, `scatter` returns `None`; otherwise, it returns a new `Ray`
    /// giving the reflected/refracted direction of the light, and a `Vec3` with
    /// the amount of energy reflected/refracted in each of red, green, and
    /// blue.
    ///
    /// (In reality, light would be *both* reflected and refracted, but we
    /// choose one or the other randomly and use over-sampling to produce a
    /// blend.)
    pub fn scatter(&self, ray: &Ray, hit: &HitRecord, rng: &mut impl Rng) -> Option<(Ray, Vec3)> {
        match self {
            Material::Lambertian { albedo } => {
                let target = hit.p + hit.normal + Vec3::in_unit_sphere(rng);
                let scattered = Ray {
                    origin: hit.p,
                    direction: target - hit.p,
                    time: ray.time,
                };
                Some((scattered, albedo(hit.p)))
            }
            Material::Metal { albedo, fuzz } => {
                let scattered = Ray {
                    origin: hit.p,
                    direction: reflect(ray.direction.into_unit(), hit.normal)
                        + *fuzz * Vec3::in_unit_sphere(rng),
                    ..*ray
                };
                if scattered.direction.dot(hit.normal) > 0. {
                    Some((scattered, *albedo))
                } else {
                    // TODO: this is in the original, but has the odd effect of
                    // making metal an emitter.
                    None
                }
            }
            Material::Dielectric { ref_idx } => {
                let (outward_normal, ni_over_nt, cosine) = if ray.direction.dot(hit.normal) > 0. {
                    (
                        -hit.normal,
                        *ref_idx,
                        *ref_idx * ray.direction.dot(hit.normal) / ray.direction.length(),
                    )
                } else {
                    (
                        hit.normal,
                        1.0 / *ref_idx,
                        -ray.direction.dot(hit.normal) / ray.direction.length(),
                    )
                };

                let direction = refract(ray.direction, outward_normal, ni_over_nt)
                    .filter(|_| rng.gen::<f32>() >= schlick(cosine, *ref_idx))
                    .unwrap_or_else(|| reflect(ray.direction, hit.normal));

                let attenuation = Vec3::from(1.);
                let ray = Ray {
                    origin: hit.p,
                    direction,
                    time: ray.time,
                };
                Some((ray, attenuation))
            }
            Material::DiffuseLight { .. } => None,
            Material::Isotropic { albedo } => Some((
                Ray {
                    origin: hit.p,
                    direction: Vec3::in_unit_sphere(rng),
                    ..*ray
                },
                albedo(hit.p),
            )),
        }
    }

    pub fn emitted(&self, p: Vec3) -> Vec3 {
        match self {
            Material::DiffuseLight {
                emission,
                brightness,
            } => *brightness * emission(p),
            _ => Vec3::default(),
        }
    }
}

impl std::fmt::Debug for Material {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        f.write_str("...")
    }
}

/// [Schlick's approximation][schlick] for computing reflection vs. refraction
/// at a material surface.
///
/// [schlick]: https://en.wikipedia.org/wiki/Schlick%27s_approximation
#[inline]
fn schlick(cos: f32, ref_idx: f32) -> f32 {
    let r0 = (1. - ref_idx) / (1. + ref_idx);
    let r0 = r0 * r0;
    r0 + (1. - r0) * f32::powf(1. - cos, 5.)
}
