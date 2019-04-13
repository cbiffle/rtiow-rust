use std::ops::Range;

use crate::aabb::Aabb;
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
#[derive(Debug)]
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
    FlipNormals(Box<Object>),
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
    pub fn hit<'o>(&'o self, ray: &Ray, t_range: Range<f32>) -> Option<HitRecord<'o>> {
        match self {
            Object::Sphere {
                center,
                radius,
                material,
                motion,
            } => {
                let t_center = ray.time * *motion + *center;

                let oc = ray.origin - t_center;
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

                let p = ray.point_at_parameter(t);
                let mut normal = Vec3::default();
                normal[*orthogonal_to] = 1.;
                Some(HitRecord {
                    t,
                    p,
                    material,
                    normal,
                })
            }
            Object::FlipNormals(o) => o
                .hit(ray, t_range)
                .map(|h| HitRecord { normal: -h.normal, ..h }),
        }
    }

    pub fn bounding_box(&self, exposure: std::ops::Range<f32>) -> Aabb {
        match self {
            Object::Sphere {
                center,
                radius,
                motion,
                ..
            } => {
                let p1 = *center + exposure.start * *motion;
                let bb1 = Aabb {
                    min: p1 - Vec3::from(*radius),
                    max: p1 + Vec3::from(*radius),
                };
                let p2 = *center + exposure.end * *motion;
                let bb2 = Aabb {
                    min: p2 - Vec3::from(*radius),
                    max: p2 + Vec3::from(*radius),
                };

                bb1.merge(bb2)
            }
            Object::Rect {
                orthogonal_to,
                k,
                range0,
                range1,
                ..
            } => {
                let mut min = Vec3::default();
                let mut max = Vec3::default();
                let z_axis = *orthogonal_to;
                let (x_axis, y_axis) = other_axes(z_axis);

                // TODO: this uses epsilon fudge factors and I hate it
                min[z_axis] = *k - 0.0001;
                max[z_axis] = *k + 0.0001;
                min[x_axis] = range0.start;
                max[x_axis] = range0.end;
                min[y_axis] = range1.start;
                max[y_axis] = range1.end;

                Aabb { min, max }
            }
            Object::FlipNormals(o) => o.bounding_box(exposure),
        }
    }
}

pub fn rect_prism(p0: Vec3, p1: Vec3, material: Material) -> Vec<Object> {
    vec![
        Object::Rect {
            orthogonal_to: Z,
            range0: p0[X]..p1[X],
            range1: p0[Y]..p1[Y],
            k: p1[Z],
            material: material.clone(),
        },
        Object::FlipNormals(Box::new(Object::Rect {
            orthogonal_to: Z,
            range0: p0[X]..p1[X],
            range1: p0[Y]..p1[Y],
            k: p0[Z],
            material: material.clone(),
        })),
        Object::Rect {
            orthogonal_to: Y,
            range0: p0[X]..p1[X],
            range1: p0[Z]..p1[Z],
            k: p1[Y],
            material: material.clone(),
        },
        Object::FlipNormals(Box::new(Object::Rect {
            orthogonal_to: Y,
            range0: p0[X]..p1[X],
            range1: p0[Z]..p1[Z],
            k: p0[Y],
            material: material.clone(),
        })),
        Object::Rect {
            orthogonal_to: X,
            range0: p0[Y]..p1[Y],
            range1: p0[Z]..p1[Z],
            k: p1[X],
            material: material.clone(),
        },
        Object::FlipNormals(Box::new(Object::Rect {
            orthogonal_to: X,
            range0: p0[Y]..p1[Y],
            range1: p0[Z]..p1[Z],
            k: p0[X],
            material: material.clone(),
        })),
    ]
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
