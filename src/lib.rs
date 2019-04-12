pub mod camera;
mod material;
mod object;
mod perlin;
pub mod ray;
pub mod texture;
pub mod vec3;

use rand::prelude::*;
use rayon::prelude::*;

use crate::camera::Camera;
use crate::material::Material;
use crate::object::{hit_slice, Object};
use crate::ray::Ray;
use crate::vec3::{Channel::*, *};

/// Computes the pixel color along `ray` for the scene of objects `world`.
///
/// This is the actual ray-tracing routine.
pub fn color(world: &[Object], mut ray: Ray, rng: &mut impl Rng) -> Vec3 {
    let mut strength = Vec3::from(1.);
    let mut emitted = Vec3::default();
    let mut bounces = 0;

    while let Some(hit) = hit_slice(world, &ray) {
        if bounces < 50 {
            if let Some((new_ray, attenuation)) = hit.material.scatter(&ray, &hit, rng) {
                ray = new_ray;
                emitted = emitted + strength * hit.material.emitted(hit.p);
                strength = strength * attenuation;
                bounces += 1;
                continue;
            } else {
                return emitted + strength * hit.material.emitted(hit.p);
            }
        }
        return Vec3::default();
    }

    Vec3::default()
}

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
        Object::RectXY {
            x_range: 3. .. 5.,
            y_range: 1. .. 3.,
            k: -2.,
            material: Material::DiffuseLight {
                emission: texture::constant(Vec3(1., 1., 1.)),
                brightness: 4.,
            },
        },
    ]
}

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

pub fn par_cast(nx: usize, ny: usize, ns: usize, camera: &Camera, world: &[Object]) -> Image {
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
    world: &[Object],
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
