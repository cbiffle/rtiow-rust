mod material;
mod object;
mod ray;
mod vec3;

use rand::prelude::*;
use rayon::prelude::*;

use crate::material::Material;
use crate::object::{HitRecord, Object, hit_slice};
use crate::vec3::{*, Axis::*, Channel::*};
use crate::ray::Ray;

fn color(world: &[Object], mut ray: Ray) -> Vec3 {
    let mut strength = Vec3::from(1.);
    let mut bounces = 0;

    while let Some(hit) = hit_slice(world, &ray, 0.001..std::f32::MAX) {
        if bounces < 50 {
            if let Some((new_ray, attenuation)) = hit.material.scatter(&ray, &hit) {
                ray = new_ray;
                strength = strength * attenuation;
                bounces += 1;
                continue;
            }
        }
        return Vec3::default();
    }

    let unit_direction = ray.direction.into_unit();
    let t = 0.5 * (unit_direction[Y] + 1.0);
    let col = (1. - t) * Vec3::from(1.) + t * Vec3(0.5, 0.7, 1.0);
    strength * col
}

struct Camera {
    origin: Vec3,
    lower_left_corner: Vec3,
    horizontal: Vec3,
    vertical: Vec3,
    u: Vec3,
    v: Vec3,
    lens_radius: f32,
}

impl Camera {
    fn look(
        look_from: Vec3,
        look_at: Vec3,
        up: Vec3,
        fov: f32,
        aspect: f32,
        aperture: f32,
        focus_dist: f32,
    ) -> Self {
        let lens_radius = aperture / 2.;
        let theta = fov * std::f32::consts::PI / 180.;
        let half_height = f32::tan(theta / 2.);
        let half_width = aspect * half_height;
        let origin = look_from;
        let w = (look_from - look_at).into_unit();
        let u = up.cross(&w).into_unit();
        let v = w.cross(&u);
        let lower_left_corner =
            origin - half_width * focus_dist * u - half_height * focus_dist * v - focus_dist * w;
        let horizontal = 2. * half_width * focus_dist * u;
        let vertical = 2. * half_height * focus_dist * v;
        Camera {
            origin,
            lower_left_corner,
            horizontal,
            vertical,
            u,
            v,
            lens_radius,
        }
    }

    fn get_ray(&self, s: f32, t: f32) -> Ray {
        let rd = self.lens_radius * Vec3::in_unit_disc();
        let offset = rd[X] * self.u + rd[Y] * self.v;
        Ray {
            origin: self.origin + offset,
            direction: self.lower_left_corner + s * self.horizontal + t * self.vertical
                - self.origin
                - offset,
        }
    }
}

fn random_scene() -> Vec<Object> {
    let mut world = vec![Object::Sphere {
        center: Vec3(0., -1000., 0.),
        radius: 1000.,
        material: Material::Lambertian {
            albedo: Vec3::from(0.5),
        },
    }];

    let mut rng = rand::thread_rng();

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
                            albedo: rng.gen::<Vec3>() * rng.gen::<Vec3>(),
                        },
                    }
                } else if choose_mat < 0.95 {
                    Object::Sphere {
                        center,
                        radius: 0.2,
                        material: Material::Metal {
                            albedo: 0.5 * (1. + rng.gen::<Vec3>()),
                            fuzz: 0.5 * rng.gen::<f32>(),
                        },
                    }
                } else {
                    Object::Sphere {
                        center,
                        radius: 0.2,
                        material: Material::Dielectric { ref_idx: 1.5 },
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
    });

    world.push(Object::Sphere {
        center: Vec3(-4., 1., 0.),
        radius: 1.0,
        material: Material::Lambertian {
            albedo: Vec3(0.4, 0.2, 0.1),
        },
    });

    world.push(Object::Sphere {
        center: Vec3(4., 1., 0.),
        radius: 1.0,
        material: Material::Metal {
            albedo: Vec3(0.7, 0.6, 0.5),
            fuzz: 0.,
        },
    });

    world
}

struct Image(Vec<Vec<Vec3>>);

impl Image {
    fn compute(nx: usize, ny: usize, f: impl Fn(usize, usize) -> Vec3 + Sync) -> Image {
        Image(
            (0..ny)
                .into_par_iter()
                .rev()
                .map(|y| (0..nx).map(|x| f(x, y)).collect())
                .collect(),
        )
    }
}

fn print_ppm(image: Image) {
    println!("P3\n{} {}\n255", image.0[0].len(), image.0.len());
    for scanline in image.0 {
        for col in scanline {
            let col = Vec3(col.0.sqrt(), col.1.sqrt(), col.2.sqrt());

            let ir = (255.99 * col[R]) as i32;
            let ig = (255.99 * col[G]) as i32;
            let ib = (255.99 * col[B]) as i32;

            println!("{} {} {}", ir, ig, ib);
        }
    }
}

fn main() {
    const NX: usize = 1200;
    const NY: usize = 800;
    const NS: usize = 10;

    let world = random_scene();

    let look_from = Vec3(13., 2., 3.);
    let look_at = Vec3(0., 0., 0.);
    let dist_to_focus = 10.;
    let aperture = 0.1;

    let camera = Camera::look(
        look_from,
        look_at,
        Vec3(0., 1., 0.),
        20.,
        NX as f32 / NY as f32,
        aperture,
        dist_to_focus,
    );

    let image = Image::compute(NX, NY, |x, y| {
        let col: Vec3 = (0..NS)
            .map(|_| {
                let mut rng = rand::thread_rng();
                let u = (x as f32 + rng.gen::<f32>()) / NX as f32;
                let v = (y as f32 + rng.gen::<f32>()) / NY as f32;
                let r = camera.get_ray(u, v);
                color(&world, r)
            })
            .sum();
        col / NS as f32
    });
    print_ppm(image);
}
