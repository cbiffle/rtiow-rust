use crate::vec3::Vec3;
use crate::material::Material;
use crate::ray::Ray;

#[derive(Clone, Debug)]
pub struct HitRecord<'m> {
    pub t: f32,
    pub p: Vec3,
    pub normal: Vec3,
    pub material: &'m Material,
}

pub enum Object {
    Sphere {
        center: Vec3,
        radius: f32,
        material: Material,
    },
}

impl Object {
    pub fn hit<'m>(&'m self, ray: &Ray, t_range: std::ops::Range<f32>) -> Option<HitRecord<'m>> {
        let Object::Sphere {
            center,
            radius,
            material,
        } = self;

        let oc = ray.origin - *center;
        let a = ray.direction.dot(ray.direction);
        let b = oc.dot(ray.direction);
        let c = oc.dot(oc) - radius * radius;
        let discriminant = b * b - a * c;
        if discriminant > 0. {
            for &t in &[
                (-b - discriminant.sqrt()) / a,
                (-b + discriminant.sqrt()) / a,
            ] {
                if t < t_range.end && t >= t_range.start {
                    let p = ray.point_at_parameter(t);
                    return Some(HitRecord {
                        t,
                        p,
                        normal: (p - *center) / *radius,
                        material: &*material,
                    });
                }
            }
        }
        None
    }
}

pub fn hit_slice<'m>(
    slice: &'m [Object],
    ray: &Ray,
    t_range: std::ops::Range<f32>,
) -> Option<HitRecord<'m>> {
    slice.iter().fold(None, |hit, obj| {
        if let Some(rec) = obj.hit(ray, t_range.clone()) {
            let hit_t = hit.as_ref().map(|h| h.t).unwrap_or(t_range.end);
            if rec.t < hit_t {
                return Some(rec);
            }
        }
        hit
    })
}


