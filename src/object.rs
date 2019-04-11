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
    #[inline]
    pub fn hit<'m>(&'m self, ray: &Ray, t_range: std::ops::Range<f32>) -> Option<HitRecord<'m>> {
        // Since there's only one kind of object right now, we can simplify the
        // match:
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
