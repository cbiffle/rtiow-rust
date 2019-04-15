#![deny(unsafe_code)]

mod aabb;
pub mod bvh;
pub mod camera;
pub mod material;
pub mod object;
mod perlin;
pub mod ray;
pub mod texture;
pub mod vec3;

use rand::prelude::*;
use rayon::prelude::*;

use crate::camera::Camera;
use crate::material::Material;
use crate::object::Object;
use crate::ray::Ray;
use crate::vec3::{Channel::*, *};

pub trait World: Send + Sync {
    fn hit_top<'a>(&'a self, ray: &Ray, rng: &mut impl Rng) -> Option<object::HitRecord<'a>>;
}

impl<'r, T: World + ?Sized> World for &'r T {
    fn hit_top<'a>(&'a self, ray: &Ray, rng: &mut impl Rng) -> Option<object::HitRecord<'a>> {
        (*self).hit_top(ray, rng)
    }
}

impl World for [Box<dyn Object>] {
    fn hit_top<'a>(&'a self, ray: &Ray, rng: &mut impl Rng) -> Option<object::HitRecord<'a>> {
        const NEAR: f32 = 0.001;

        let mut nearest = std::f32::MAX;
        let mut hit = None;

        for obj in self {
            if let Some(rec) = obj.hit(ray, NEAR..nearest, &mut || rng.gen()) {
                nearest = rec.t;
                hit = Some(rec);
            }
        }

        hit
    }
}

impl World for bvh::Bvh {
    fn hit_top<'a>(&'a self, ray: &Ray, rng: &mut impl Rng) -> Option<object::HitRecord<'a>> {
        self.hit(ray, 0.001..std::f32::MAX, &mut || rng.gen())
    }
}

/// Computes the pixel color along `ray` for the scene of objects `world`.
///
/// This is the actual ray-tracing routine.
pub fn color(world: &impl World, mut ray: Ray, rng: &mut impl Rng) -> Vec3 {
    // Accumulates contribution of each surface we reach.
    let mut accum = Vec3::default();
    // Records the cumulative (product) attenuation of each surface we've
    // visited so far.
    let mut strength = Vec3::from(1.);

    let mut bounces = 0;

    // Iterate until one of the following conditions is reached:
    // 1. The ray escapes into space (i.e. no objects are hit).
    // 2. The ray reaches a surface that does not scatter.
    // 3. The ray bounces more than 50 times.
    while let Some(hit) = world.hit_top(&ray, rng) {

        // Record this hit's contribution, attenuated by the total attenuation
        // so far.
        accum = accum + strength * hit.material.emitted(hit.p);

        // Check whether the material scatters light, generating a new ray. In
        // practice this is true for everything but the emission-only
        // DiffuseLight type.
        //
        // TODO: and also for frosted metal, which effectively makes frosted
        // metal an emitter. That can't be right.
        if let Some((new_ray, attenuation)) = hit.material.scatter(&ray, &hit, rng) {
            // Redirect flight, accumulate the new attenuation value.
            ray = new_ray;
            strength = strength * attenuation;
        } else {
            // Locally absorbed; we're done.
            return accum
        }

        if bounces == 50 { return accum }

        bounces += 1;
    }

    Vec3::default()
}

pub fn cornell_box() -> Vec<Box<dyn Object>> {
    fn diffuse_color(c: Vec3) -> Material {
        Material::Lambertian {
            albedo: texture::constant(c),
        }
    }

    let red = diffuse_color(Vec3(0.65, 0.05, 0.05));
    let white = diffuse_color(Vec3::from(0.73));
    let green = diffuse_color(Vec3(0.12, 0.45, 0.15));
    let light = Material::DiffuseLight {
        emission: texture::constant(Vec3::from(1.)),
        brightness: 15.,
    };
    vec![
        Box::new(object::Rect {
            orthogonal_to: object::StaticY,
            range0: 213. ..343.,
            range1: 227. ..332.,
            k: 554.,
            material: light,
        }),
        // floor
        Box::new(object::Rect {
            orthogonal_to: object::StaticY,
            range0: 0. ..555.,
            range1: 0. ..555.,
            k: 0.,
            material: white.clone(),
        }),
        // rear wall
        Box::new(object::FlipNormals(object::Rect {
            orthogonal_to: object::StaticZ,
            range0: 0. ..555.,
            range1: 0. ..555.,
            k: 555.,
            material: white.clone(),
        })),
        // ceiling
        Box::new(object::FlipNormals(object::Rect {
            orthogonal_to: object::StaticY,
            range0: 0. ..555.,
            range1: 0. ..555.,
            k: 555.,
            material: white.clone(),
        })),
        // right wall
        Box::new(object::Rect {
            orthogonal_to: object::StaticX,
            range0: 0. ..555.,
            range1: 0. ..555.,
            k: 0.,
            material: red,
        }),
        // left wall
        Box::new(object::FlipNormals(object::Rect {
            orthogonal_to: object::StaticX,
            range0: 0. ..555.,
            range1: 0. ..555.,
            k: 555.,
            material: green,
        })),
    ]
}

pub fn cornell_box_with_boxes() -> Vec<Box<dyn Object>> {
    fn diffuse_color(c: Vec3) -> Material {
        Material::Lambertian {
            albedo: texture::constant(c),
        }
    }

    let mut scene = cornell_box();
    let white = diffuse_color(Vec3::from(0.73));

    scene.push(
        Box::new(object::Translate {
            offset: Vec3(130., 0., 65.),
            object: object::rotate_y(
                -18.,
                object::rect_prism(Vec3(0., 0., 0.), Vec3(165., 165., 165.), white.clone()),
            ),
        })
    );
    scene.push(
        Box::new(object::Translate {
            offset: Vec3(265., 0., 295.),
            object: object::rotate_y(
                15.,
                object::rect_prism(Vec3(0., 0., 0.), Vec3(165., 330., 165.), white),
            ),
        })
    );
    scene
}

/*
pub fn simple_light() -> Vec<Object> {
    vec![
        Object::Sphere {
            center: Vec3(0., -1000., 0.),
            radius: 1000.,
            material: Material::Lambertian {
                albedo: texture::perlin(4.),
            },
            motion: Vec3::default(),
        },
        Object::Sphere {
            center: Vec3(0., 2., 0.),
            radius: 2.,
            material: Material::Lambertian {
                albedo: texture::perlin(4.),
            },
            motion: Vec3::default(),
        },
        Object::Sphere {
            center: Vec3(0., 7., 0.),
            radius: 2.,
            material: Material::DiffuseLight {
                emission: texture::constant(Vec3(1., 1., 1.)),
                brightness: 4.,
            },
            motion: Vec3::default(),
        },
        Object::Rect {
            orthogonal_to: object::StaticAxis::Z,
            range0: 3. ..5.,
            range1: 1. ..3.,
            k: -2.,
            material: Material::DiffuseLight {
                emission: texture::constant(Vec3(1., 1., 1.)),
                brightness: 4.,
            },
        },
    ]
}
*/

/*
pub fn random_scene(rng: &mut impl Rng) -> Vec<Object> {
    let mut world = vec![Object::Sphere {
        center: Vec3(0., -1000., 0.),
        radius: 1000.,
        material: Material::Lambertian {
            albedo: texture::perlin(4.),
        },
        motion: Vec3::default(),
    }];

    for a in -11..11 {
        for b in -11..11 {
            let center = Vec3(
                a as f32 + 0.9 * rng.gen::<f32>(),
                0.2,
                b as f32 + 0.9 * rng.gen::<f32>(),
            );
            if (center - Vec3(4., 0.2, 0.)).length() > 0.9 {
                let choose_mat = rng.gen::<f32>();

                let obj = if choose_mat < 0.8 {
                    Object::Sphere {
                        center,
                        radius: 0.2,
                        material: Material::Lambertian {
                            albedo: texture::constant(rng.gen::<Vec3>() * rng.gen::<Vec3>()),
                        },
                        motion: Vec3(0., rng.gen_range(0., 0.5), 0.),
                    }
                } else if choose_mat < 0.95 {
                    Object::Sphere {
                        center,
                        radius: 0.2,
                        material: Material::Metal {
                            albedo: 0.5 * (1. + rng.gen::<Vec3>()),
                            fuzz: 0.5 * rng.gen::<f32>(),
                        },
                        motion: Vec3::default(),
                    }
                } else {
                    Object::Sphere {
                        center,
                        radius: 0.2,
                        material: Material::Dielectric { ref_idx: 1.5 },
                        motion: Vec3::default(),
                    }
                };
                world.push(obj);
            }
        }
    }

    world.push(Object::Sphere {
        center: Vec3(0., 1., 0.),
        radius: 1.0,
        material: Material::Dielectric { ref_idx: 1.5 },
        motion: Vec3::default(),
    });

    world.push(Object::Sphere {
        center: Vec3(-4., 1., 0.),
        radius: 1.0,
        material: Material::Metal {
            albedo: Vec3(0.7, 0.6, 0.5),
            fuzz: 0.,
        },
        motion: Vec3::default(),
    });

    world.push(Object::Sphere {
        center: Vec3(4., 1., 0.),
        radius: 1.0,
        material: Material::DiffuseLight {
            emission: texture::perlin(10.),
            brightness: 4.,
        },
        motion: Vec3::default(),
    });

    world
}
*/

pub struct Image(Vec<Vec<Vec3>>);

impl Image {
    pub fn par_compute(nx: usize, ny: usize, f: impl Fn(usize, usize) -> Vec3 + Sync) -> Image {
        Image(
            (0..ny)
                .into_par_iter()
                .rev()
                .map(|y| (0..nx).map(|x| f(x, y)).collect())
                .collect(),
        )
    }

    pub fn compute(nx: usize, ny: usize, mut f: impl FnMut(usize, usize) -> Vec3) -> Image {
        Image(
            (0..ny)
                .rev()
                .map(|y| (0..nx).map(|x| f(x, y)).collect())
                .collect(),
        )
    }
}

pub fn print_ppm(image: Image) {
    println!("P3\n{} {}\n255", image.0[0].len(), image.0.len());
    for scanline in image.0 {
        for col in scanline {
            let col = Vec3(col.0.sqrt(), col.1.sqrt(), col.2.sqrt());

            fn to_u8(x: f32) -> i32 {
                ((255.99 * x) as i32).max(0).min(255)
            }

            let ir = to_u8(col[R]);
            let ig = to_u8(col[G]);
            let ib = to_u8(col[B]);

            println!("{} {} {}", ir, ig, ib);
        }
    }
}

pub fn par_cast(nx: usize, ny: usize, ns: usize, camera: &Camera, world: impl World) -> Image {
    Image::par_compute(nx, ny, |x, y| {
        let col: Vec3 = (0..ns)
            .map(|_| {
                let mut rng = rand::thread_rng();
                let u = (x as f32 + rng.gen::<f32>()) / nx as f32;
                let v = (y as f32 + rng.gen::<f32>()) / ny as f32;
                let r = camera.get_ray(u, v, &mut rng);
                color(&world, r, &mut rng)
            })
            .sum();
        col / ns as f32
    })
}

pub fn cast(
    nx: usize,
    ny: usize,
    ns: usize,
    camera: &Camera,
    world: impl World,
    rng: &mut impl Rng,
) -> Image {
    Image::compute(nx, ny, |x, y| {
        let col: Vec3 = (0..ns)
            .map(|_| {
                let u = (x as f32 + rng.gen::<f32>()) / nx as f32;
                let v = (y as f32 + rng.gen::<f32>()) / ny as f32;
                let r = camera.get_ray(u, v, rng);
                color(&world, r, rng)
            })
            .sum();
        col / ns as f32
    })
}
