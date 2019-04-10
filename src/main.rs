use rand::prelude::*;

use rtiow::*;
use rtiow::vec3::Vec3;

fn main() {
    const NX: usize = 1200;
    const NY: usize = 800;
    const NS: usize = 10;

    let mut rng = rand::rngs::SmallRng::seed_from_u64(0xDEADBEEF);

    let world = random_scene(&mut rng);

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

    let image = par_cast(NX, NY, NS, &camera, &world);
    print_ppm(image);
}
