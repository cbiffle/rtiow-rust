use rand::prelude::*;
use std::ops::Range;
use std::time::Instant;

use rtiow::camera::Camera;
use rtiow::object::{self, Object};
use rtiow::vec3::Vec3;
use rtiow::*;

#[allow(unused)]
fn cornell_box_scene(nx: usize, ny: usize) -> (Vec<Box<dyn Object>>, Camera, Range<f32>) {
    let look_from = Vec3(278., 278., -800.);
    let look_at = Vec3(278., 278., 0.);
    let dist_to_focus = 10.;
    let aperture = 0.0;
    let exposure = 0. ..1.;

    let camera = Camera::look(
        look_from,
        look_at,
        Vec3(0., 1., 0.),
        40.,
        nx as f32 / ny as f32,
        aperture,
        dist_to_focus,
        exposure.clone(),
    );

    (cornell_box_with_boxes(), camera, exposure)
}

#[allow(unused)]
fn motion_test(nx: usize, ny: usize) -> (Vec<Box<dyn Object>>, Camera, Range<f32>) {
    let look_from = Vec3(278., 278., -800.);
    let look_at = Vec3(278., 278., 0.);
    let dist_to_focus = 10.;
    let aperture = 0.0;
    let exposure = 0. ..1.;

    let camera = Camera::look(
        look_from,
        look_at,
        Vec3(0., 1., 0.),
        40.,
        nx as f32 / ny as f32,
        aperture,
        dist_to_focus,
        exposure.clone(),
    );

    let mut scene = cornell_box();

    scene.push(Box::new(object::Translate {
        offset: Vec3(278., 278., 278.),
        object: object::LinearMove {
            motion: Vec3(0., 100., 0.),
            object: object::Sphere {
                radius: 65.,
                material: rtiow::material::Material::Lambertian {
                    albedo: rtiow::texture::constant(Vec3::from(0.73)),
                },
            },
        },
    }));

    (scene, camera, exposure)
}

#[allow(unused)]
fn volume_test(nx: usize, ny: usize) -> (Vec<Box<dyn Object>>, Camera, Range<f32>) {
    let look_from = Vec3(278., 278., -800.);
    let look_at = Vec3(278., 278., 0.);
    let dist_to_focus = 10.;
    let aperture = 0.0;
    let exposure = 0. ..1.;

    let camera = Camera::look(
        look_from,
        look_at,
        Vec3(0., 1., 0.),
        40.,
        nx as f32 / ny as f32,
        aperture,
        dist_to_focus,
        exposure.clone(),
    );

    let mut scene = cornell_box();

    scene.push(Box::new(object::Translate {
        offset: Vec3(278., 278., 278.),
        object: object::ConstantMedium {
            boundary: object::Sphere {
                radius: 180.,
                // material does not matter here
                material: rtiow::material::Material::Lambertian {
                    albedo: rtiow::texture::constant(Vec3::from(0.73)),
                },
            },
            density: 0.01,
            material: rtiow::material::Material::Isotropic {
                albedo: rtiow::texture::constant(Vec3(0.2, 0.2, 1.0)),
            },
        },
    }));

    (scene, camera, exposure)
}

#[allow(unused)]
fn simple_light_scene(
    nx: usize,
    ny: usize,
    rng: &mut impl Rng,
) -> (Vec<Box<dyn Object>>, Camera, Range<f32>) {
    let look_from = Vec3(278., 278., -800.);
    let look_at = Vec3(278., 278., 0.);
    let dist_to_focus = 10.;
    let aperture = 0.0;
    let exposure = 0. ..1.;

    let camera = Camera::look(
        look_from,
        look_at,
        Vec3(0., 1., 0.),
        40.,
        nx as f32 / ny as f32,
        aperture,
        dist_to_focus,
        exposure.clone(),
    );

    use rtiow::material::Material;

    let mut world = cornell_box();

    const SPHERES: usize = 1000;
    for _ in 0..SPHERES {
        world.push(Box::new(object::Translate {
            object: object::Sphere {
                radius: 20.,
                material: Material::Lambertian {
                    albedo: rtiow::texture::constant(Vec3::from(0.3)),
                },
            },
            offset: 277. + 257. * rng.gen::<Vec3>(),
        }));
    }

    world.push(Box::new(object::FlipNormals(object::Sphere {
        radius: 1000.,
        material: Material::DiffuseLight {
            emission: rtiow::texture::constant(Vec3::from(0.1)),
            brightness: 1.,
        },
    })));

    (world, camera, exposure)
}

fn book_final_scene(
    nx: usize,
    ny: usize,
    rng: &mut impl Rng,
) -> (Vec<Box<dyn Object>>, Camera, Range<f32>) {
    let look_from = Vec3(478., 278., -600.);
    let look_at = Vec3(278., 278., 0.);
    let dist_to_focus = 10.;
    let aperture = 0.0;
    let exposure = 0. ..1.;

    let camera = Camera::look(
        look_from,
        look_at,
        Vec3(0., 1., 0.),
        40.,
        nx as f32 / ny as f32,
        aperture,
        dist_to_focus,
        exposure.clone(),
    );

    use rtiow::material::Material;
    use rtiow::texture;

    let ground = Material::Lambertian {
        albedo: texture::constant(Vec3(0.48, 0.83, 0.53)),
    };

    let mut world: Vec<Box<dyn Object>> = vec![];

    // Make random floor.
    for i in 0..20 {
        for j in 0..20 {
            const W: f32 = 100.;
            let c0 = Vec3(-1000. + i as f32 * W, 0., -1000. + j as f32 * W);
            let c1 = c0 + Vec3(W, 100. * (rng.gen::<f32>() + 0.01), W);
            world.push(Box::new(object::rect_prism(c0, c1, ground.clone())));
        }
    }

    // Make light.
    world.push(Box::new(object::Rect {
        orthogonal_to: object::StaticY,
        range0: 123. ..423.,
        range1: 147. ..412.,
        k: 554.,
        material: Material::DiffuseLight {
            emission: texture::constant(Vec3::from(1.)),
            brightness: 7.,
        },
    }));

    // Brown blurry sphere.
    world.push(Box::new(object::Translate {
        offset: Vec3(400., 400., 200.),
        object: object::LinearMove {
            motion: Vec3(30., 0., 0.),
            object: object::Sphere {
                radius: 50.,
                material: Material::Lambertian {
                    albedo: texture::constant(Vec3(0.7, 0.3, 0.1)),
                },
            },
        },
    }));

    let glass = Material::Dielectric { ref_idx: 1.5 };

    // Glass sphere.
    world.push(Box::new(object::Translate {
        offset: Vec3(260., 150., 45.),
        object: object::Sphere {
            radius: 50.,
            material: glass.clone(),
        },
    }));

    // Silvery sphere.
    world.push(Box::new(object::Translate {
        offset: Vec3(0., 150., 145.),
        object: object::Sphere {
            radius: 50.,
            material: Material::Metal {
                albedo: Vec3(0.8, 0.8, 0.9),
                fuzz: 10.,
            },
        },
    }));

    // Blue glass sphere.
    let boundary = object::Translate {
        offset: Vec3(360., 150., 145.),
        object: object::Sphere {
            radius: 70.,
            material: glass.clone(),
        },
    };
    world.push(Box::new(boundary.clone()));
    world.push(Box::new(object::ConstantMedium {
        boundary,
        density: 0.2,
        material: Material::Isotropic {
            albedo: texture::constant(Vec3(0.2, 0.4, 0.9)),
        },
    }));

    // Fog.
    world.push(Box::new(object::ConstantMedium {
        boundary: object::Sphere {
            radius: 5000.,
            material: glass.clone(), // doesn't matter
        },
        density: 0.0001,
        material: Material::Isotropic {
            albedo: texture::constant(Vec3::from(1.)),
        },
    }));

    // Perlin marbled sphere.
    world.push(Box::new(object::Translate {
        offset: Vec3(220., 280., 300.),
        object: object::Sphere {
            radius: 80.,
            material: Material::Lambertian {
                albedo: texture::perlin(0.05),
            },
        },
    }));

    // Cube made of random spheres.
    world.push(Box::new({
        const SPHERES: usize = 1000;
        let white = Material::Lambertian {
            albedo: texture::constant(Vec3::from(0.73)),
        };
        let spheres = (0..SPHERES).map(|_| {
            Box::new(object::Translate {
                offset: 165. * rng.gen::<Vec3>(),
                object: object::Sphere {
                    radius: 10.,
                    material: white.clone(),
                },
            }) as Box<dyn Object>
        }).collect();
        let bvh = rtiow::bvh::from_scene(spheres, exposure.clone(), rng);
        object::Translate {
            offset: Vec3(-100., 270., 395.),
            object: object::rotate_y(
                15.,
                bvh,
            ),
        }
    }));

    (world, camera, exposure)
}

const USE_BVH: bool = true;

fn main() {
    const NX: usize = 300;
    const NY: usize = 300;
    const NS: usize = 100;

    eprintln!(
        "Parallel casting {} x {} image using {}x oversampling.",
        NX, NY, NS
    );

    let mut rng = rand::rngs::SmallRng::seed_from_u64(0xDEADBEEF);

    //let (world, camera, exposure) = cornell_box_scene(NX, NY);
    //let (world, camera, exposure) = simple_light_scene(NX, NY, &mut rng);
    //let (world, camera, exposure) = volume_test(NX, NY);
    let (world, camera, exposure) = book_final_scene(NX, NY, &mut rng);

    let (image, time) = if USE_BVH {
        eprintln!("Generating bounding volume hierarchy.");
        let world = rtiow::bvh::from_scene(world, exposure, &mut rng);
        eprintln!("Done.");
        let start = Instant::now();
        (par_cast(NX, NY, NS, &camera, world), start.elapsed())
    } else {
        eprintln!("Testing every ray against every object.");
        let world: &[Box<dyn Object>] = &world;
        let start = Instant::now();
        (par_cast(NX, NY, NS, &camera, world), start.elapsed())
    };

    eprintln!("Took {:?} wall time", time);

    print_ppm(image);
}
