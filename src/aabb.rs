use crate::ray::Ray;
use crate::vec3::Vec3;

#[derive(Copy, Clone, Debug)]
pub struct Aabb {
    pub min: Vec3,
    pub max: Vec3,
}

impl Aabb {
    pub fn merge(self, other: Aabb) -> Self {
        Aabb {
            min: self.min.zip_with(other.min, f32::min),
            max: self.max.zip_with(other.max, f32::max),
        }
    }

    pub fn hit(&self, ray: &Ray, mut t_range: std::ops::Range<f32>) -> bool {
        use crate::vec3::Axis::*;

        for &axis in &[X, Y, Z] {
            let inv_d = 1. / ray.direction[axis];
            let t0 = (self.min[axis] - ray.origin[axis]) * inv_d;
            let t1 = (self.max[axis] - ray.origin[axis]) * inv_d;
            let (t0, t1) = if inv_d < 0. { (t1, t0) } else { (t0, t1) };
            t_range.start = t_range.start.max(t0);
            t_range.end = t_range.end.min(t1);
            if t_range.end <= t_range.start {
                return false;
            }
        }
        true
    }

    // experimental, totally untested, vector rewrite of the above algorithm
    pub fn hit_v(&self, ray: &Ray, mut t_range: std::ops::Range<f32>) -> bool {
        let inv_d = ray.direction.map(|x| 1. / x);
        let t0 = (self.min - ray.origin) * inv_d;
        let t1 = (self.max - ray.origin) * inv_d;

        let t0_ = inv_d.zip_with3(t0, t1, |i, a, b| if i < 0. { b } else { a });
        let t1_ = inv_d.zip_with3(t1, t0, |i, a, b| if i < 0. { b } else { a });

        t_range.start = t_range.start.max(t0_.reduce(f32::max));
        t_range.end = t_range.end.max(t1_.reduce(f32::min));

        t_range.end > t_range.start
    }
}
