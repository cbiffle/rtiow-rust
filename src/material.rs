use rand::prelude::*;

use crate::vec3::{Vec3, reflect, refract};
use crate::ray::Ray;

// TODO
use super::HitRecord;

#[derive(Clone, Debug)]
pub enum Material {
    Lambertian { albedo: Vec3 },
    Metal { albedo: Vec3, fuzz: f32 },
    Dielectric { ref_idx: f32 },
}

impl Material {
    pub fn scatter(&self, ray: &Ray, hit: &HitRecord) -> Option<(Ray, Vec3)> {
        match self {
            Material::Lambertian { albedo } => {
                let target = hit.p + hit.normal + Vec3::in_unit_sphere();
                let scattered = Ray {
                    origin: hit.p,
                    direction: target - hit.p,
                };
                Some((scattered, *albedo))
            }
            Material::Metal { albedo, fuzz } => {
                let scattered = Ray {
                    origin: hit.p,
                    direction: reflect(ray.direction.into_unit(), hit.normal)
                        + *fuzz * Vec3::in_unit_sphere(),
                };
                if scattered.direction.dot(hit.normal) > 0. {
                    Some((scattered, *albedo))
                } else {
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
                    .filter(|_| rand::thread_rng().gen::<f32>() >= schlick(cosine, *ref_idx))
                    .unwrap_or_else(|| reflect(ray.direction, hit.normal));

                let attenuation = Vec3::from(1.);
                let ray = Ray {
                    origin: hit.p,
                    direction,
                };
                Some((ray, attenuation))
            }
        }
    }
}

fn schlick(cos: f32, ref_idx: f32) -> f32 {
    let r0 = (1. - ref_idx) / (1. + ref_idx);
    let r0 = r0 * r0;
    r0 + (1. - r0) * f32::powf(1. - cos, 5.)
}
