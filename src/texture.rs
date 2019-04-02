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
