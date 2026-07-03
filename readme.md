# faster? f32 math functions
Attempting to provide faster implementations of common f32 math functions with a similar level of accuracy to the standard library.
There are the beginnings of a perfectly rounded tier (cbrt_accurate) and a
faster-but-sloppier tier (cbrt_throughput ~6.7 avg ulp, cbrt_fast ~57 avg
ulp — experiments, unpolished); eventually both tiers should cover every
function.

All functions are full-range correct — negatives, denormals, zero, inf and nan
are handled (previously cbrt was NaN for negative/huge/denormal inputs, log2
broke on denormals, sin returned 0 for tiny arguments, and cbrt_accurate
overflowed to NaN above 2^127) — except exp2, where the default is the fast
unchecked version: it is only valid for x in [-126, 128) (normal, finite,
nonzero results) and returns garbage outside; exp2_checked handles the full
range (overflow to inf, denormal underflow, nan) for ~2.5 ns extra latency.
cbrt_accurate is perfectly rounded on every input tested.

All edge handling is branchless (selects / clamps, single evaluation of the
core path) so loops over arrays auto-vectorize — that's a hard requirement
here. Branchy edge handling blocks LLVM's if-conversion and silently drops
you back to scalar code, ~5-14x slower in throughput.

# precision (ulp error, dense full-domain sweep, see examples/accuracy.rs)
```
              | jodie avg | jodie max | std avg | std max
--------------|-----------|-----------|---------|--------
         cbrt |   0.085   |     1     |    0    |    0
cbrt_accurate |   0.000   |     0     |    0    |    0
         exp2 |   0.059   |     2     |  0.000  |    1
 exp2_checked |   0.059   |     2     |  0.000  |    1
         log2 |   0.006   |     2     |  0.000  |    1
 sin (|x|<1e3)|   0.027   |     2     |  0.002  |    1
 cos (|x|<1e3)|   0.082   |     2     |  0.001  |    1
```
(exp2 was previously reported as max 1: the accuracy harness silently
swept zero samples on negative ranges — measure() walks bit patterns
upward — so negative inputs were never tested. With the sweep fixed,
exp2 is 2 ulp at worst for x just below 0, from fma rounding near
result = 1.0 where the output ulp shrinks; sin/cos numbers were
confirmed unchanged by the fix.)
sin/cos use Cody-Waite reduction (4-constant fma chain mod pi); like all
non-Payne-Hanek reductions they degrade for |x| beyond ~1e6.

# benchmarks
Run on i5-1145G7, -C target-cpu=native (now set in .cargo/config.toml — without
hardware fma, mul_add becomes a libm call and everything is ~2x slower).
This machine thermally throttles up to ~2.5x: only trust ratios measured
within a single run, and check the `nop` baseline (~0.017 ns throughput when
cool now that the throughput loop vectorizes).

Three different metrics; earlier readme "latency" numbers were actually the
criterion per-call numbers (independent inputs overlap in the pipeline, so
they measure pipelined per-call time, not serial latency).

Serial latency (dependency chain, examples/quickbench.rs; lower is better)
              | jodie   | std     | improvement
--------------|---------|---------|------------
         cbrt | 11.6 ns | 21.2 ns | 1.8x
cbrt_accurate | 20.0 ns | 21.2 ns | 1.1x
          cos | 12.7 ns | 12.6 ns | 1.0x
         exp2 |  8.5 ns |  9.6 ns | 1.1x
 exp2_checked | 11.0 ns |  9.6 ns | 0.9x
         log2 |  8.9 ns | 10.0 ns | 1.1x
          sin | 10.8 ns | 12.8 ns | 1.2x

Throughput (independent array evals over [f32; 4096], auto-vectorized with
AVX2; examples/quickbench.rs; lower is better)
              | jodie    | std     | improvement
--------------|----------|---------|------------
         cbrt | 0.60 ns  | 3.89 ns | 6.5x
cbrt_accurate | 0.69 ns  | 3.89 ns | 5.6x
          cos | 0.27 ns  | 2.90 ns | 10.7x
         exp2 | 0.23 ns  | 2.27 ns | 9.9x
 exp2_checked | 0.38 ns  | 2.27 ns | 6.0x
         log2 | 0.36 ns  | 2.75 ns | 7.7x
          sin | 0.21 ns  | 2.98 ns | 13.9x

The throughput gap vs std comes almost entirely from vectorization: std's
functions have branches, so LLVM can't vectorize loops that call them.
Note the throughput harness must use fixed-size arrays ([f32; N], zip
iterators): a Vec through black_box has an opaque length, the bounds
checks survive, and the whole loop silently stays scalar.

exp2_checked pays ~2.5 ns of serial latency for branchless full-range
handling (2^k is split into two exact power-of-two factors so overflow
and denormal underflow fall out of the final multiplies, and nan
propagates through the polynomial — no selects at all). log2 folds its
denormal-rescale correction into the polynomial's exponent term, which
made it faster than the old branchy version in both metrics.

The criterion benches (benches/benches.rs, N=1) measure roughly the same
thing as the throughput table but with more per-call overhead and no
repetition-minimum, so they're noisier under thermal drift. Note their
input used to be 12345.0, which is outside exp2's domain (the old 2.1 ns
exp2 number was the wrap-around garbage path); the input is now 12.345.
Also, std sin/cos are ~30% faster for small arguments than for large
ones, so bench inputs matter for the std columns.

# tools
- `cargo run --release --example accuracy [filter]` - avg/max ulp over dense strided sweeps of each domain
- `cargo run --release --example quickbench [filter]` - latency (serial dependency chain) + throughput, min of 7 reps
- `cargo run --release --example edgecheck` - bit-exact checks of edge cases (0, -0, denormals, inf, nan, domain boundaries)
- `cargo run --release --example tune` - coordinate-descent ulp tuning of polynomial coefficients

# todo:
- do principled and thourough analysis of dependency chains and rounding errors to find optimizations
- add inverse trig functions
- add tan()
- perfectly rounded versions
- Payne-Hanek fallback for huge trig arguments
