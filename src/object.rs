use std::ops::Range;

use crate::aabb::Aabb;
use crate::material::Material;
use crate::ray::Ray;
use crate::vec3::{
    Axis::{self, *},
    Vec3,
};

pub trait StaticAxis: std::fmt::Debug + Send + Sync {
    const AXIS: Axis;
    const OTHER1: Axis;
    const OTHER2: Axis;
}

#[derive(Debug)]
pub struct StaticX;

impl StaticAxis for StaticX {
    const AXIS: Axis = Axis::X;
    const OTHER1: Axis = Axis::Y;
    const OTHER2: Axis = Axis::Z;
}

#[derive(Debug)]
pub struct StaticY;

impl StaticAxis for StaticY {
    const AXIS: Axis = Axis::Y;
    const OTHER1: Axis = Axis::X;
    const OTHER2: Axis = Axis::Z;
}

#[derive(Debug)]
pub struct StaticZ;

impl StaticAxis for StaticZ {
    const AXIS: Axis = Axis::Z;
    const OTHER1: Axis = Axis::X;
    const OTHER2: Axis = Axis::Y;
}

/// An object in a scene.
///
/// The primary purpose of an `Object` is to interact with rays of light using
/// the `hit` method.
pub trait Object: std::fmt::Debug + Sync + Send {
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
    fn hit<'o>(&'o self, ray: &Ray, t_range: Range<f32>) -> Option<HitRecord<'o>>;

    /// Computes the bounding box for the object at the given range of times.
    /// This is called during scene setup, not rendering, and so it may be
    /// expensive.
    fn bounding_box(&self, exposure: std::ops::Range<f32>) -> Aabb;
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

/// A sphere.
#[derive(Debug)]
pub struct Sphere {
    /// Radius of the sphere.
    pub radius: f32,
    /// Material of the sphere.
    pub material: Material,
}

impl Object for Sphere {
    #[inline]
    fn hit<'o>(&'o self, ray: &Ray, t_range: Range<f32>) -> Option<HitRecord<'o>> {
        let a = ray.direction.dot(ray.direction);
        let b = ray.origin.dot(ray.direction);
        let c = ray.origin.dot(ray.origin) - self.radius * self.radius;
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
                        normal: p / self.radius,
                        material: &self.material,
                    });
                }
            }
        }
        None
    }

    fn bounding_box(&self, _exposure: std::ops::Range<f32>) -> Aabb {
        Aabb {
            min: -Vec3::from(self.radius),
            max: Vec3::from(self.radius),
        }
    }
}

#[derive(Debug)]
pub struct Rect<A> {
    pub orthogonal_to: A,
    pub range0: Range<f32>,
    pub range1: Range<f32>,
    // TODO: replace with Translate?
    pub k: f32,
    pub material: Material,
}

impl<A: StaticAxis> Object for Rect<A> {
    #[inline]
    fn hit<'o>(&'o self, ray: &Ray, t_range: Range<f32>) -> Option<HitRecord<'o>> {
        // The names x and y below are correct for orthogonal_to=Z. Use your
        // imagination for the other cases.

        let t = (self.k - ray.origin[A::AXIS]) / ray.direction[A::AXIS];
        if t < t_range.start || t >= t_range.end {
            return None;
        }

        let x = ray.origin[A::OTHER1] + t * ray.direction[A::OTHER1];
        let y = ray.origin[A::OTHER2] + t * ray.direction[A::OTHER2];
        if x < self.range0.start
            || x >= self.range0.end
            || y < self.range1.start
            || y >= self.range1.end
        {
            return None;
        }

        let p = ray.point_at_parameter(t);
        let mut normal = Vec3::default();
        normal[A::AXIS] = 1.;
        Some(HitRecord {
            t,
            p,
            material: &self.material,
            normal,
        })
    }

    fn bounding_box(&self, exposure: std::ops::Range<f32>) -> Aabb {
        let mut min = Vec3::default();
        let mut max = Vec3::default();

        // TODO: this uses epsilon fudge factors and I hate it
        min[A::AXIS] = self.k - 0.0001;
        max[A::AXIS] = self.k + 0.0001;
        min[A::OTHER1] = self.range0.start;
        max[A::OTHER1] = self.range0.end;
        min[A::OTHER2] = self.range1.start;
        max[A::OTHER2] = self.range1.end;

        Aabb { min, max }
    }
}

#[derive(Debug)]
pub struct FlipNormals<T>(pub T);

impl<T: Object> Object for FlipNormals<T> {
    #[inline]
    fn hit<'o>(&'o self, ray: &Ray, t_range: Range<f32>) -> Option<HitRecord<'o>> {
        self.0.hit(ray, t_range).map(|h| HitRecord {
            normal: -h.normal,
            ..h
        })
    }

    fn bounding_box(&self, exposure: std::ops::Range<f32>) -> Aabb {
        self.0.bounding_box(exposure)
    }
}

#[derive(Debug)]
pub struct Translate<T> {
    pub offset: Vec3,
    pub object: T,
}

impl<T: Object> Object for Translate<T> {
    #[inline]
    fn hit<'o>(&'o self, ray: &Ray, t_range: Range<f32>) -> Option<HitRecord<'o>> {
        let t_ray = Ray {
            origin: ray.origin - self.offset,
            ..*ray
        };
        self.object.hit(&t_ray, t_range).map(|hit| HitRecord {
            p: hit.p + self.offset,
            ..hit
        })
    }

    fn bounding_box(&self, exposure: std::ops::Range<f32>) -> Aabb {
        let b = self.object.bounding_box(exposure);
        Aabb {
            min: b.min + self.offset,
            max: b.max + self.offset,
        }
    }
}

#[derive(Debug)]
pub struct RotateY<T> {
    pub object: T,
    pub sin_theta: f32,
    pub cos_theta: f32,
}

impl<T: Object> Object for RotateY<T> {
    #[inline]
    fn hit<'o>(&'o self, ray: &Ray, t_range: Range<f32>) -> Option<HitRecord<'o>> {
        fn rot(p: Vec3, sin_theta: f32, cos_theta: f32) -> Vec3 {
            Vec3(
                p.dot(Vec3(cos_theta, 0., sin_theta)),
                p.dot(Vec3(0., 1., 0.)),
                p.dot(Vec3(-sin_theta, 0., cos_theta)),
            )
        }

        let rot_ray = Ray {
            origin: rot(ray.origin, -self.sin_theta, self.cos_theta),
            direction: rot(ray.direction, -self.sin_theta, self.cos_theta),
            ..*ray
        };

        self.object.hit(&rot_ray, t_range).map(|hit| HitRecord {
            p: rot(hit.p, self.sin_theta, self.cos_theta),
            normal: rot(hit.normal, self.sin_theta, self.cos_theta),
            ..hit
        })
    }

    fn bounding_box(&self, exposure: std::ops::Range<f32>) -> Aabb {
        fn rot(p: Vec3, sin_theta: f32, cos_theta: f32) -> Vec3 {
            Vec3(
                p.dot(Vec3(cos_theta, 0., sin_theta)),
                p.dot(Vec3(0., 1., 0.)),
                p.dot(Vec3(-sin_theta, 0., cos_theta)),
            )
        }

        let (min, max) = self.object.bounding_box(exposure).corners().fold(
            (Vec3::from(std::f32::MAX), Vec3::from(std::f32::MIN)),
            |(min, max), c| {
                let rot_c = rot(c, self.sin_theta, self.cos_theta);
                (min.zip_with(rot_c, f32::min), max.zip_with(rot_c, f32::max))
            },
        );
        Aabb { min, max }
    }
}

#[derive(Debug)]
pub struct Composite<T>(pub T, pub Vec<T>);

impl<T: Object> Object for Composite<T> {
    fn hit<'o>(&'o self, ray: &Ray, mut t_range: Range<f32>) -> Option<HitRecord<'o>> {
        let mut hit = None;
        for object in self.1.iter().chain(std::iter::once(&self.0)) {
            if let Some(rec) = object.hit(ray, t_range.clone()) {
                t_range.end = rec.t;
                hit = Some(rec)
            }
        }
        hit
    }

    fn bounding_box(&self, exposure: std::ops::Range<f32>) -> Aabb {
        self.1.iter().fold(self.0.bounding_box(exposure.clone()), |a, b| {
            a.merge(b.bounding_box(exposure.clone()))
        })
    }
}

#[derive(Debug)]
pub struct And<T, S>(pub T, pub S);

impl<T: Object, S: Object> Object for And<T, S> {
    fn hit<'o>(&'o self, ray: &Ray, mut t_range: Range<f32>) -> Option<HitRecord<'o>> {
        let hit0 = self.0.hit(ray, t_range.clone());
        if let Some(h) = &hit0 {
            t_range.end = h.t
        }

        let hit1 = self.1.hit(ray, t_range);
        hit1.or(hit0)
    }

    fn bounding_box(&self, exposure: std::ops::Range<f32>) -> Aabb {
        self.0
            .bounding_box(exposure.clone())
            .merge(self.1.bounding_box(exposure))
    }
}

pub fn rect_prism(p0: Vec3, p1: Vec3, material: Material) -> impl Object {
    And(
        And(
            Rect {
                orthogonal_to: StaticZ,
                range0: p0[X]..p1[X],
                range1: p0[Y]..p1[Y],
                k: p1[Z],
                material: material.clone(),
            },
            And(
                Rect {
                    orthogonal_to: StaticY,
                    range0: p0[X]..p1[X],
                    range1: p0[Z]..p1[Z],
                    k: p1[Y],
                    material: material.clone(),
                },
                Rect {
                    orthogonal_to: StaticX,
                    range0: p0[Y]..p1[Y],
                    range1: p0[Z]..p1[Z],
                    k: p1[X],
                    material: material.clone(),
                },
            ),
        ),
        And(
            FlipNormals(Rect {
                orthogonal_to: StaticZ,
                range0: p0[X]..p1[X],
                range1: p0[Y]..p1[Y],
                k: p0[Z],
                material: material.clone(),
            }),
            And(
                FlipNormals(Rect {
                    orthogonal_to: StaticY,
                    range0: p0[X]..p1[X],
                    range1: p0[Z]..p1[Z],
                    k: p0[Y],
                    material: material.clone(),
                }),
                FlipNormals(Rect {
                    orthogonal_to: StaticX,
                    range0: p0[Y]..p1[Y],
                    range1: p0[Z]..p1[Z],
                    k: p0[X],
                    material: material.clone(),
                }),
            ),
        ),
    )
}

pub fn rotate_y<O: Object>(degrees: f32, object: O) -> RotateY<O> {
    let radians = degrees * std::f32::consts::PI / 180.;
    RotateY {
        object,
        sin_theta: radians.sin(),
        cos_theta: radians.cos(),
    }
}

#[derive(Debug)]
pub struct LinearMove<O> {
    pub object: O,
    pub motion: Vec3,
}

impl<O: Object> Object for LinearMove<O> {
    #[inline]
    fn hit<'o>(&'o self, ray: &Ray, t_range: Range<f32>) -> Option<HitRecord<'o>> {
        self.object.hit(&Ray { origin: ray.origin - ray.time * self.motion, ..*ray}, t_range)
    }

    fn bounding_box(&self, exposure: std::ops::Range<f32>) -> Aabb {
        let bb = self.object.bounding_box(exposure.clone());

        let bb_start = Aabb {
            min: bb.min + exposure.start * self.motion,
            max: bb.max + exposure.start * self.motion,
        };
        let bb_end = Aabb {
            min: bb.min + exposure.end * self.motion,
            max: bb.max + exposure.end * self.motion,
        };

        bb_start.merge(bb_end)
    }
}

