#[derive(Copy, Clone, Default, Debug)]
struct Vec3(f32, f32, f32);

impl Vec3 {
    pub fn dot(&self, other: &Self) -> f32 {
        self.0 * other.0 + self.1 * other.1 + self.2 * other.2
    }

    pub fn cross(&self, other: &Self) -> Self {
        Vec3(
            self.1 * other.2 - self.2 * other.1,
            -(self.0 * other.2 - self.2 * other.0),
            self.0 * other.1 - self.1 * other.0,
        )
    }

    pub fn length(&self) -> f32 {
        f32::sqrt(self.dot(self))
    }

    pub fn into_unit(self) -> Self {
        self / self.length()
    }
}

impl std::ops::Mul<Vec3> for f32 {
    type Output = Vec3;

    fn mul(self, rhs: Vec3) -> Self::Output {
        Vec3(self * rhs.0, self * rhs.1, self * rhs.2)
    }
}

impl std::ops::Div<f32> for Vec3 {
    type Output = Vec3;

    fn div(self, rhs: f32) -> Self::Output {
        Vec3(self.0 / rhs, self.1 / rhs, self.2 / rhs)
    }
}

impl std::ops::Add for Vec3 {
    type Output = Vec3;

    fn add(self, rhs: Vec3) -> Self::Output {
        Vec3(self.0 + rhs.0, self.1 + rhs.1, self.2 + rhs.2)
    }
}

impl std::ops::Sub for Vec3 {
    type Output = Vec3;

    fn sub(self, rhs: Vec3) -> Self::Output {
        Vec3(self.0 - rhs.0, self.1 - rhs.1, self.2 - rhs.2)
    }
}

#[derive(Copy, Clone, Debug)]
pub enum Channel { R, G, B }

use Channel::*;


impl ::std::ops::Index<Channel> for Vec3 {
    type Output = f32;

    fn index(&self, idx: Channel) -> &Self::Output {
        match idx {
            R => &self.0,
            G => &self.1,
            B => &self.2,
        }
    }
}

#[derive(Copy, Clone, Debug)]
pub enum Axis { X, Y, Z }

use Axis::*;


impl ::std::ops::Index<Axis> for Vec3 {
    type Output = f32;

    fn index(&self, idx: Axis) -> &Self::Output {
        match idx {
            X => &self.0,
            Y => &self.1,
            Z => &self.2,
        }
    }
}

#[derive(Copy, Clone, Debug)]
struct Ray {
    origin: Vec3,
    direction: Vec3,
}

impl Ray {
    pub fn point_at_parameter(&self, t: f32) -> Vec3 {
        self.origin + t * self.direction
    }
}

#[derive(Copy, Clone, Debug)]
pub struct HitRecord {
    t: f32,
    p: Vec3,
    normal: Vec3,
}

trait Object {
    fn hit(&self, ray: &Ray, t_range: std::ops::Range<f32>) -> Option<HitRecord>;
}

#[derive(Clone, Debug)]
pub struct Sphere {
    center: Vec3,
    radius: f32,
}

impl Object for Sphere {
    fn hit(&self, ray: &Ray, t_range: std::ops::Range<f32>) -> Option<HitRecord> {
        let oc = ray.origin - self.center;
        let a = ray.direction.dot(&ray.direction);
        let b = oc.dot(&ray.direction);
        let c = oc.dot(&oc) - self.radius * self.radius;
        let discriminant = b * b - a * c;
        if discriminant > 0. {
            for &t in &[(-b - discriminant.sqrt()) / a, (-b + discriminant.sqrt()) / a] {
                if t < t_range.end && t >= t_range.start {
                    let p = ray.point_at_parameter(t);
                    return Some(HitRecord {
                        t,
                        p,
                        normal: (p - self.center) / self.radius,
                    })
                }
            }
        }
        None
    }
}

impl<T> Object for [T] where T: std::ops::Deref<Target = dyn Object> {
    fn hit(&self, ray: &Ray, t_range: std::ops::Range<f32>) -> Option<HitRecord> {
        self.iter()
            .fold(None, |hit, obj| {
                if let Some(rec) = obj.hit(ray, t_range.clone()) {
                    let hit_t = hit.map(|h| h.t).unwrap_or(t_range.end);
                    if rec.t < hit_t {
                        return Some(rec)
                    }
                }
                hit
            })
    }
}

fn color(ray: &Ray, world: &[Box<dyn Object>]) -> Vec3 {
    if let Some(hit) = world.hit(ray, 0. .. std::f32::MAX) {
        return 0.5 * Vec3(hit.normal[X] + 1., hit.normal[Y] + 1., hit.normal[Z] + 1.)
    }

    let unit_direction = ray.direction.into_unit();
    let t = 0.5 * (unit_direction[Y] + 1.0);
    debug_assert!(t >= 0. && t <= 1., "{:?}", unit_direction);
    (1. - t) * Vec3(1., 1., 1.) + t * Vec3(0.5, 0.7, 1.0)
}

fn main() {
    const NX: usize = 200;
    const NY: usize = 100;

    println!("P3\n{} {}\n255", NX, NY);

    let lower_left_corner = Vec3(-2., -1., -1.);
    let horizontal = Vec3(4., 0., 0.);
    let vertical = Vec3(0., 2., 0.);
    let origin = Vec3::default();

    let world: [Box<dyn Object>; 2] = [
        Box::new(Sphere { center: Vec3(0., 0., -1.), radius: 0.5 }),
        Box::new(Sphere { center: Vec3(0., -100.5, -1.), radius: 100. }),
    ];

    for j in (0..NY).rev() {
        for i in 0..NX {
            let u = i as f32 / NX as f32;
            let v = j as f32 / NY as f32;
            let r = Ray {
                origin,
                direction: lower_left_corner + u*horizontal + v*vertical,
            };
            let col = color(&r, &world);

            let ir = (255.99 * col[R]) as i32;
            let ig = (255.99 * col[G]) as i32;
            let ib = (255.99 * col[B]) as i32;

            println!("{} {} {}", ir, ig, ib);
        }
    }
}
