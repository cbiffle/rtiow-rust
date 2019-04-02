use crate::perlin;
use crate::vec3::Vec3;

pub type Texture = Box<dyn Fn(Vec3) -> Vec3 + Sync>;

pub fn constant(color: Vec3) -> Texture {
    Box::new(move |_| color)
}

pub fn checker(t0: Texture, t1: Texture) -> Texture {
    Box::new(move |p| {
        let s = (10. * p).map(f32::sin).reduce(std::ops::Mul::mul);
        if s < 0. {
            t1(p)
        } else {
            t0(p)
        }
    })
}

pub fn perlin(scale: f32) -> Texture {
    //Box::new(move |p| Vec3::from(perlin::turb(scale * p, 7)))
    Box::new(move |p| Vec3::from(0.5 * (1. + f32::sin(scale * p.2 + 10. * perlin::turb(p, 7)))))
}
