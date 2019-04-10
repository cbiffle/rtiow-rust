//! A three-vector of floats, used as a color, coordinate, etc.

use rand::prelude::*;

#[derive(Copy, Clone, Default, Debug)]
pub struct Vec3(pub f32, pub f32, pub f32);

impl Vec3 {
    pub fn in_unit_sphere() -> Self {
        let mut rng = rand::thread_rng();
        loop {
            let v = 2. * rng.gen::<Vec3>() - Vec3::from(1.);
            if v.dot(v) < 1. {
                return v;
            }
        }
    }

    pub fn in_unit_disc() -> Self {
        let mut rng = rand::thread_rng();
        loop {
            let v = 2. * Vec3(rng.gen(), rng.gen(), 0.) - Vec3(1., 1., 0.);
            if v.dot(v) < 1. {
                return v;
            }
        }
    }

    pub fn dot(&self, other: Self) -> f32 {
        self.zip_with(other, core::ops::Mul::mul)
            .reduce(core::ops::Add::add)
    }

    pub fn cross(&self, other: &Self) -> Self {
        Vec3(
            self.1 * other.2 - self.2 * other.1,
            -(self.0 * other.2 - self.2 * other.0),
            self.0 * other.1 - self.1 * other.0,
        )
    }

    pub fn length(&self) -> f32 {
        self.dot(*self).sqrt()
    }

    pub fn into_unit(self) -> Self {
        self / self.length()
    }

    pub fn map(self, mut f: impl FnMut(f32) -> f32) -> Self {
        Vec3(f(self.0), f(self.1), f(self.2))
    }

    pub fn zip_with(self, other: Vec3, mut f: impl FnMut(f32, f32) -> f32) -> Self {
        Vec3(f(self.0, other.0), f(self.1, other.1), f(self.2, other.2))
    }

    pub fn reduce(self, f: impl Fn(f32, f32) -> f32) -> f32 {
        f(f(self.0, self.1), self.2)
    }
}

/// Broadcasts a single value to all vector lanes.
impl From<f32> for Vec3 {
    fn from(v: f32) -> Self {
        Vec3(v, v, v)
    }
}

/// Element-wise multiplication (Hadamard product). I have reservations about
/// making this available as `*`, but it sure is convenient...
impl std::ops::Mul for Vec3 {
    type Output = Vec3;

    fn mul(self, rhs: Vec3) -> Self::Output {
        self.zip_with(rhs, std::ops::Mul::mul)
    }
}

/// `scalar * vector`
impl std::ops::Mul<Vec3> for f32 {
    type Output = Vec3;

    fn mul(self, rhs: Vec3) -> Self::Output {
        Vec3::from(self) * rhs
    }
}

/// `vector / scalar`
impl std::ops::Div<f32> for Vec3 {
    type Output = Vec3;

    fn div(self, rhs: f32) -> Self::Output {
        self.map(|x| x / rhs)
    }
}

impl std::ops::Add for Vec3 {
    type Output = Vec3;

    fn add(self, rhs: Vec3) -> Self::Output {
        self.zip_with(rhs, std::ops::Add::add)
    }
}

impl std::ops::Add<Vec3> for f32 {
    type Output = Vec3;

    fn add(self, rhs: Vec3) -> Self::Output {
        rhs.map(|x| self + x)
    }
}

impl std::ops::Sub for Vec3 {
    type Output = Vec3;

    fn sub(self, rhs: Vec3) -> Self::Output {
        self.zip_with(rhs, std::ops::Sub::sub)
    }
}

impl std::ops::Neg for Vec3 {
    type Output = Vec3;

    fn neg(self) -> Self::Output {
        self.map(std::ops::Neg::neg)
    }
}

impl std::iter::Sum for Vec3 {
    fn sum<I>(iter: I) -> Self
    where
        I: Iterator<Item = Self>,
    {
        iter.fold(Vec3::default(), std::ops::Add::add)
    }
}

/// Allow `Vec3` to be produced by `Rng::gen`.
///
/// The resulting vector has each component in the half-open range `[0,1)`. Note
/// that this is *not* a unit vector.
impl rand::distributions::Distribution<Vec3> for rand::distributions::Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> Vec3 {
        Vec3(rng.gen(), rng.gen(), rng.gen())
    }
}

/// Names for vector lanes when used as a color.
#[derive(Copy, Clone, Debug)]
pub enum Channel {
    R,
    G,
    B,
}

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

/// Names for vector lanes when used as a coordinate.
#[derive(Copy, Clone, Debug)]
pub enum Axis {
    X,
    Y,
    Z,
}

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

pub fn reflect(v: Vec3, n: Vec3) -> Vec3 {
    v - 2. * v.dot(n) * n
}

pub fn refract(v: Vec3, n: Vec3, ni_over_nt: f32) -> Option<Vec3> {
    let uv = v.into_unit();
    let dt = uv.dot(n);
    let discriminant = 1.0 - ni_over_nt * ni_over_nt * (1. - dt * dt);
    if discriminant > 0. {
        Some(ni_over_nt * (uv - dt * n) - discriminant.sqrt() * n)
    } else {
        None
    }
}


