# Rust One-Weekend Ray Tracer

![A demo scene rendered by this ray tracer at 1200x800 with 50x
oversampling](img/demo-scene.jpg)

*Demo scene from the book rendered at 1200x800 at 50x oversampling, in 52s.*

This is an implementation of the basic ray tracer described in Peter Shirley's
book [Ray Tracing In One Weekend][rtiow]. While that book describes an
implementation in C++, I don't believe in producing new C++ code, so this
version is in Rust.

## Trying it out

[Install Rust][install-rust], clone this repo, and type:

```shell
cargo run --release > out.ppm
```

Now view `out.ppm` with the image viewer of your choice.

If you make changes, run `cargo bench` before and after to check for performance
regressions.

## Differences from the book

The overall structure of the code is very similar to the C++, but there are some
differences.

1. The algorithm implementations use idiomatic Rust. For example,
   "out-parameters" have been eliminated, functions that may or may not return a
   result now use `Option`, and so on.

2. Images are rendered into memory before being printed to `stdout`. This makes
   concurrency (below) easier.

3. The ray propagation routine is now iterative, not recursive, which lets me
   play with higher bounce limits without blowing the stack.

4. The ray tracer will distribute rendering over available CPU cores. Because
   we're in Rust, this took [about one line of code][smp-commit] and is
   statically free of data races.

## Performance

At the time of this writing, on my Skylake Thinkpad (Intel i7-8550U, 4 cores / 8
threads), this implementation casts about 1M rays/second. The demo scene at
1200x800, 10x oversampling, takes a bit over 8 seconds. The image at the top of
this README took 52s.

[rtiow]: http://www.realtimerendering.com/raytracing/Ray%20Tracing%20in%20a%20Weekend.pdf
[smp-commit]: https://github.com/cbiffle/rtiow-rust/commit/fcebb909ae95a8b9ad83cd970be54128c1f4c629
