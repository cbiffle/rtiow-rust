use rand::prelude::*;
use std::ops::Range;
use std::time::Instant;

use rtiow::camera::Camera;
use rtiow::object::{self, Object};
use rtiow::vec3::Vec3;
use rtiow::*;

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
                motion: Vec3::default(),
            },
            offset: 277. + 257. * rng.gen::<Vec3>(),
        }));
    }

    (world, camera, exposure)
}

const USE_BVH: bool = true;

fn main() {
    const NX: usize = 300;
    const NY: usize = 300;
    const NS: usize = 500;

    eprintln!(
        "Parallel casting {} x {} image using {}x oversampling.",
        NX, NY, NS
    );

    let mut rng = rand::rngs::SmallRng::seed_from_u64(0xDEADBEEF);

    //let (world, camera, exposure) = cornell_box_scene(NX, NY);
    let (world, camera, exposure) = simple_light_scene(NX, NY, &mut rng);

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
