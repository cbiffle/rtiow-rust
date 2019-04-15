# Rust One-Week-ish Ray Tracer

This is an implementation of the basic ray tracer described in Peter Shirley's
books [Ray Tracing In One Weekend][rtiow] and [Ray Tracing: The Next
Week][rttnw]. While those books describe an implementation in C++, I don't
believe in producing new C++ code, so this version is in Rust.

More information after the two pretty pictures:

![A demo scene rendered by this ray tracer at 1200x800 with 50x
oversampling](img/demo-scene.jpg)

*Demo scene from the first book rendered at 1200x800 at 50x oversampling, in
52s.*

![Demo scene from the second book, rendered by this ray tracer at 800x800, 5000x
oversampling](img/rttnw-final.jpg)

*Demo scene from the second book, showing sub-surface scattering, volumetric
fog, motion blur, etc. 5000x oversampled, 100 minutes.*

## Trying it out

[Install Rust][install-rust], clone this repo, and type:

```shell
cargo run --release > out.ppm
```

Now view `out.ppm` with the image viewer of your choice.

If you make changes, run `cargo bench` before and after to check for performance
regressions.

## Comparison

This section has two goals:

1. To help you read the original C++ codebase and the Rust codebase for
   comparison purposes.

2. To make my larger point about why I do not produce new C/C++ code.

### Structure and implementation

The overall structure of the code is vaguely similar to the C++, but there are
some differences, and those differences are growing with time. This list is not
exhaustive.

- The algorithm implementations use idiomatic Rust. For example,
  "out-parameters" have been eliminated, functions that may or may not return a
  result now use `Option`, and so on.

- Images are rendered into memory before being printed to `stdout`. This makes
  concurrency (below) easier.

- The ray propagation routine is now iterative, not recursive, which lets me
  play with higher bounce limits without blowing the stack. (It also improves
  code generation.)

- Essentially all of the math code, including bounding box intersection tests,
  is phrased so that it gets auto-vectorized into Intel AVX instructions. (No
  actual vector intrinsics are used, so the code can be compiled for older CPUs
  or ARM.)

- The ray tracer will distribute rendering over available CPU cores. Because
  we're in Rust, this took [about one line of code][smp-commit] and is
  statically free of data races.

- Material is an enum, not a class hierarchy. C++ doesn't have Rust-style
  enums, but they're a useful way of modeling a closed set of options, and
  matching on an enum is significantly cheaper than dynamic dispatch.

- Object (`hitable` in the original) uses Rust's `trait` concept to do dynamic
  dispatch where required, but static dispatch where possible. In particular,
  transformation nodes like `Translate` and `FlipNormals` integrate with the
  object they're transforming, which not only eliminates an indirection and
  heap allocation, but allows the compiler to optimize combinations like
  `Translate+Sphere` together.

- The C++ codebase contains a lot of anti-pattern pointer usage. That's all
  gone. In particular, the data structures in this implementation can be safely
  deallocated; this was not true in the original.

- Random number generator state is explicitly passed around, so that the entire
  system can be made deterministic for benchmarking.

### Performance

To compare the performance of the Rust codebase to Shirley's C++ codebase (on
Github), I've used the following settings:

- Scene: final scene from book 2, with the texture-mapped Earth sphere removed
  (because I couldn't be bothered to implement it in Rust yet).
- Computer: Skylake Thinkpad (Intel i7-8550U, 4 cores / 8 threads).
- `rustc` 1.33.0.
- GCC 8.2.1.
- Rust built with `cargo build --release`, and restricted to a single thread by
  exporting `RAYON_NUM_THREADS=1` at runtime.
- C++ code built with `g++ -O3 -march=native main.cc -o main`. (Adding
  `-ffast-math` and/or `-fomit-frame-pointer` doesn't change things
  significantly.)

Note that the C++ code is the best possible case for GCC's optimizer
(effectively a single source file with all definitions inlined), while the Rust
code is split across many files, a library target, a binary target, and uses
upstream libraries. To level the playing field, I switched on LTO in
`Cargo.toml`.

**At the time of this writing, the performance of the two programs is identical
on the test scene when each is restricted to a single CPU.** Both programs
complete in exactly 32.43s when rendering at 300x300x100. (The Rust version uses
about 10% less RAM.)

This is *despite* the Rust code technically doing more work: all array/vector
accesses are bounds-checked, certain corners of floating-point math are checked
more rigorously than in C, every potentially null pointer is checked before use,
and all memory operations are both memory-safe and thread-safe. **Remember this
next time a C programmer insists that they need to do unsafe tricks "for
performance."**

Plus, note that the Rust program is *not* a single-threaded build: it's the
SMP-aware multithreaded build being told to use only one core. When it is
allowed use of all 4 cores (8 threads), performance improves linearly with the
number of cores made available. The scene completes in 6.92s, showing that we
benefit only a little from hyperthreading.

### Lines of code

As measured by [cloc][cloc], the Rust implementation is somewhat longer than the
C++. To do a fair comparison, I excised the C++ code responsible for image
format decoding, which I didn't implement. The results:

- C++: 1,219 LOC.
- Rust: 1,647 LOC.

While Rust is generally less boilerplatey than C++, the fact that it contained
more lines of code here didn't shock me, for three reasons:

1. I was deliberately more verbose in how objects are declared, using Rust's
   struct literal syntax with named fields. (I would have done the same in C++
   -- names are nice -- but Shirley used constructor functions.)

2. `cloc` is sensitive to formatting, and I've used `rustfmt` to enforce a
   somewhat sparse style. If I run, for example, `clang-format` on the C++
   codebase, it grows to 1,330 lines.

3. The C++ code's organization into a single `.cc` file with all definitions
   inlined made it less boilerplatey than a "real" C++ codebase using separate
   header and implementation files. (On the other hand, this means it actually
   takes longer to compile than the Rust code.)

### Compile time

As noted above, the C++ codebase is organized into a single source file with
includes. As a result, while GCC is typically faster than `rustc`, the C++ ray
tracer takes longer to compile than the Rust code (3 seconds vs 2).

### Safety / reliability

Consider the amount of work required to review each codebase for possible
memory-related errors, such as buffer overruns, dangling pointers or
use-after-free, reads of uninitialized memory, null pointer dereference, and the
like. The codebases are roughly the same size (see above).

In C++, every non-blank line of code could potentially contain such errors. (In
this codebase in particular, there are a bunch of potential use-after-frees
waiting to happen, and basically every type has a default constructor that
leaves its member variables entirely uninitialized, virtually guaranteeing reads
of uninitialized memory.)

In Rust, such errors can only occur in `unsafe` blocks, and an attribute
(pragma) at the top of the ray tracer codebase *bans them.* You don't even have
to read the code to know there is no `unsafe` in it; your review is complete.

As if to make my point, when I first checked out and built the C++ code, it
segfaulted immediately. (The required `earthmap.jpg` file is not distributed
with the code, and the program handles this error by dereferencing a null
pointer.)

[rtiow]: http://www.realtimerendering.com/raytracing/Ray%20Tracing%20in%20a%20Weekend.pdf
[rttnw]: http://www.realtimerendering.com/raytracing/Ray%20Tracing_%20The%20Next%20Week.pdf
[smp-commit]: https://github.com/cbiffle/rtiow-rust/commit/fcebb909ae95a8b9ad83cd970be54128c1f4c629
[install-rust]: https://www.rust-lang.org/tools/install
[cloc]: https://github.com/AlDanial/cloc
