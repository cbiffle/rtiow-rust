use crate::material::Material;
use crate::ray::Ray;
use crate::vec3::Vec3;

/// An object in a scene.
///
/// The primary purpose of an `Object` is to interact with rays of light using
/// the `hit` method.
#[derive(Clone, Debug)]
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
    pub fn hit<'o>(&'o self, ray: &Ray, t_range: std::ops::Range<f32>) -> Option<HitRecord1<'o>> {
        // Since there's only one kind of object right now, we can simplify the
        // match:
        let Object::Sphere {
            center,
            radius,
            material,
            motion,
        } = self;

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

    /// Computes the surface normal at a given point. If the point is not on the
    /// surface, the result will be bogus. However, this is only available
    /// within this module, and we use it carefully.
    #[inline]
    fn normal_at(&self, p: Vec3) -> Vec3 {
        let Object::Sphere { center, radius, .. } = self;
        (p - *center) / *radius
    }
}

/// Initial cheap ray-object intersection record, used before we've decided
/// which object was actually hit.
#[derive(Clone, Debug)]
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
#[derive(Clone, Debug)]
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
