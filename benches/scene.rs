use criterion::{criterion_group, Criterion, BatchSize};
use rand::prelude::*;

use rtiow::vec3::Vec3;
use rtiow::*;

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("scene/seq/10x10x4", |b| {
        const NX: usize = 10;
        const NY: usize = 10;
        const NS: usize = 4;

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

        b.iter_batched(
            || (),
            |_| cast(NX, NY, NS, &camera, &world, &mut rng),
            BatchSize::SmallInput,
        );
    });
    c.bench_function("scene/par/10x10x4", |b| {
        const NX: usize = 10;
        const NY: usize = 10;
        const NS: usize = 4;

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

        b.iter_batched(
            || (),
            |_| par_cast(NX, NY, NS, &camera, &world),
            BatchSize::SmallInput,
        );
    });
}

criterion_group!(benches, criterion_benchmark);
criterion::criterion_main!(benches);
