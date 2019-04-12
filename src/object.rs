use std::ops::Range;

use crate::material::Material;
use crate::ray::Ray;
use crate::vec3::{
    Axis::{self, *},
    Vec3,
};

fn other_axes(axis: Axis) -> (Axis, Axis) {
    match axis {
        Z => (X, Y),
        X => (Y, Z),
        Y => (X, Z),
    }
}

/// An object in a scene.
///
/// The primary purpose of an `Object` is to interact with rays of light using
/// the `hit` method.
pub enum Object {
    /// A sphere.
    Sphere {
        /// Center point of the sphere.
        center: Vec3,
        /// Radius of the sphere.
        radius: f32,
        /// Material of the sphere.
        material: Material,
        /// Motion vector displacing sphere from `center` over time.
        motion: Vec3,
    },
    Rect {
        orthogonal_to: Axis,
        range0: Range<f32>,
        range1: Range<f32>,
        k: f32,
        material: Material,
    },
}

impl Object {
    /// Tests if `ray` intersects the object `self`, and if so, if that
    /// intersection occurs within `t_range` along the ray. (Recall that `Ray`
    /// is defined in terms of a `t` value that refers to points along the ray.)
    ///
    /// The `t_range` serves two purposes here. First, if the intersection
    /// occurs at *negative* `t`, the object is behind the photons instead of in
    /// front of them, and the intersection is an illusion. Second, while the
    /// upper end of `t_range` starts out as infinity, we adjust it down as we
    /// find objects along `ray`. Once we've found an object at position `t`, we
    /// can ignore any objects at positions greater than `t`.
    ///
    /// This function returns a `HitRecord1`, which can be turned into a full
    /// `HitRecord` by calling its `finish` method. (We don't return a full
    /// `HitRecord` here because we may find hits on multiple objects, but we
    /// only want to do all the work for the closest one.)
    #[inline]
    pub fn hit<'o>(&'o self, ray: &Ray, t_range: Range<f32>) -> Option<HitRecord1<'o>> {
        match self {
            Object::Sphere {
                center,
                radius,
                material,
                motion,
            } => {
                let center = ray.time * *motion + *center;

                let oc = ray.origin - center;
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
                            return Some(HitRecord1 {
                                t,
                                object: self,
                                material,
                            });
                        }
                    }
                }
                None
            }
            Object::Rect {
                orthogonal_to,
                range0,
                range1,
                k,
                material,
            } => {
                // The *_axis names are correct for orthogonal_to=Z. Use your
                // imagination for the other cases.
                let z_axis = *orthogonal_to;
                let (x_axis, y_axis) = other_axes(z_axis);

                let t = (k - ray.origin[z_axis]) / ray.direction[z_axis];
                if t < t_range.start || t >= t_range.end {
                    return None;
                }

                let x = ray.origin[x_axis] + t * ray.direction[x_axis];
                let y = ray.origin[y_axis] + t * ray.direction[y_axis];
                if x < range0.start || x >= range0.end || y < range1.start || y >= range1.end {
                    return None;
                }

                Some(HitRecord1 {
                    t,
                    object: self,
                    material,
                })
            }
        }
    }

    /// Computes the surface normal at a given point. If the point is not on the
    /// surface, the result will be bogus. However, this is only available
    /// within this module, and we use it carefully.
    #[inline]
    fn normal_at(&self, p: Vec3) -> Vec3 {
        match self {
            Object::Sphere { center, radius, .. } => (p - *center) / *radius,
            Object::Rect { orthogonal_to, .. } => {
                let mut normal = Vec3::default();
                normal[*orthogonal_to] = 1.;
                normal
            }
        }
    }
}

/// Initial cheap ray-object intersection record, used before we've decided
/// which object was actually hit.
#[derive(Clone)]
pub struct HitRecord1<'a> {
    /// Position along the ray, expressed in distance from the origin.
    t: f32,
    /// Object that was hit.
    object: &'a Object,
    /// Material that was hit.
    material: &'a Material,
}

impl<'o> HitRecord1<'o> {
    pub fn finish(self, ray: &Ray) -> HitRecord<'o> {
        let p = ray.point_at_parameter(self.t);
        HitRecord {
            t: self.t,
            p,
            normal: self.object.normal_at(p),
            material: self.material,
        }
    }
}

/// A description of a `Ray` hitting an `Object`. This stores information needed
/// for rendering later.
///
/// The `'m` lifetime refers to the `Material` of the `Object`, which we capture
/// by reference. Thus, a `HitRecord` cannot outlive the `Object` it refers to.
#[derive(Clone)]
pub struct HitRecord<'m> {
    /// Position along the ray, expressed in distance from the origin.
    pub t: f32,
    /// Position along the ray, as an actual point.
    pub p: Vec3,
    /// Surface normal of the object at the hit position.
    pub normal: Vec3,
    /// Material of the object at the hit position.
    pub material: &'m Material,
}

/// Runs through the `slice` of `Object`s looking for the closest hit for `ray`.
pub fn hit_slice<'m>(slice: &'m [Object], ray: &Ray) -> Option<HitRecord<'m>> {
    const NEAR: f32 = 0.001;

    let mut nearest = std::f32::MAX;
    let mut hit = None;

    for obj in slice {
        if let Some(rec) = obj.hit(ray, NEAR..nearest) {
            nearest = rec.t;
            hit = Some(rec);
        }
    }

    hit.map(|h| h.finish(ray))
}
