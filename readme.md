# faster? f32 math functions
Attempting to provide faster implementations of common f32 math functions with a similar level of accuracy to the standard library.
There are also perfectly rounded variants and faster-but-sloppier variants.

All functions auto-vectorize, it's a hard requirement

# precision (see examples/accuracy.rs)

```
                        | jodie avg  | jodie max | std avg | std max
------------------------|------------|-----------|---------|--------
                   cbrt |    0.326   |     3     |    0    |    0
          cbrt_accurate |    0.000   |     1     |    0    |    0
                   exp2 |    0.030   |     1     |  0.000  |    1
           exp2_checked |    0.016   |     1     |  0.000  |    1
                   log2 |    0.003   |     3     |  0.000  |    1
        sin (|x|<1.3e7) |    0.065   |   1183    |  0.003  |    1
        cos (|x|<1.3e7) |    0.293   |   2780    |  0.002  |    1
 sin_checked (|x|<=1e6) |    0.036   |     2     |  0.000  |    1
 cos_checked (|x|<=1e6) |    0.081   |     3     |  0.000  |    1
         sinpi (|x|<1e6)|    0.336   |     2     | (no std sinpi)
         cospi (|x|<1e6)|    0.094   |     2     | (no std cospi)
```

```
                        | jodie avg  | jodie max | std avg | std max
------------------------|------------|-----------|---------|--------
                     ln |    0.126   |     3     |  0.000  |    1
                  log10 |    0.286   |     4     |  0.000  |    0
                  log1p |    0.106   |     4     |  0.000  |    0
        exp (in-domain) |    0.091   |     4     |  0.000  |    1
      expm1 (in-domain) |    0.138   |     6     |  0.000  |    0
       sinh (in-domain) |    0.081   |     5     |  0.000  |    1
       cosh (in-domain) |    0.071   |     5     |  0.000  |    0
sinh_throughput (in-domain)| 0.084   |     7     |  0.000  |    1
cosh_throughput (in-domain)| 0.061   |     4     |  0.000  |    0
       tanh (in-domain) |    0.148   |     9     |  0.000  |    0
                  asinh |    0.173   |     4     | 
                  acosh |    0.063   |     4     |  0.000  |    1
                  atanh |    0.032   |     3     | 
                   asin |    0.030   |     9     |  0.000  |    0
                   acos |    0.496   |     4     |  0.000  |    0
                   atan |    0.068   |     4     |  0.000  |    0
        tan (in-domain) |    0.331   |  2967     |  0.000  |    0
                   erf  |    0.319   |     5     | (no std erf)
         erfc (|x|<=10) |    0.311   |   109     | (no std erfc)
                  atan2 |    0.069   |     3     |  0.000  |    0
        hypot (bounded) |    0.034   |     1     |  0.000  |    0
        powf (in-domain)|    0.181   |   127     |  0.000  |    1
```

# benchmarks
Run on i5-1145G7, -C target-cpu=native (now set in .cargo/config.toml)
```
Serial latency (dependency chain, examples/quickbench.rs; lower is better)
              | jodie   | std     | improvement
--------------|---------|---------|------------
         cbrt | 12.9 ns | 22.0 ns | 1.7x
cbrt_accurate | 17.0 ns | 22.0 ns | 1.3x
          cos | 12.7 ns | 12.8 ns | 1.0x
  cos_checked | 30.3 ns | 12.8 ns | 0.4x
         exp2 |  8.5 ns | 13.3 ns | 1.6x
 exp2_checked | 13.3 ns | 13.3 ns | 1.0x
         log2 | 13.1 ns | 14.9 ns | 1.1x
          sin | 10.9 ns | 12.9 ns | 1.2x
  sin_checked | 29.2 ns | 12.9 ns | 0.4x
           ln | 15.3 ns | 14.5 ns | 0.9x
        log10 | 14.3 ns | 16.6 ns | 1.2x
        log1p | 16.5 ns | 21.1 ns | 1.3x
          exp |  9.9 ns |  9.5 ns | 1.0x
        expm1 | 11.3 ns | 13.7 ns | 1.2x
         sinh | 14.0 ns | 14.0 ns | 1.0x
         cosh | 13.8 ns | 14.0 ns | 1.0x
         tanh | 17.7 ns | 16.5 ns | 0.9x
        asinh | 43.7 ns | 86.3 ns | 2.0x
        acosh | 35.1 ns | 27.0 ns | 0.8x
        atanh | 23.2 ns |  4.2 ns | 0.2x
         asin | 19.4 ns |  4.1 ns | 0.2x
         acos | 13.0 ns |  4.0 ns | 0.3x
         atan | 14.0 ns | 19.3 ns | 1.4x
        atan2 | 14.1 ns | 23.7 ns | 1.7x
          tan | 23.5 ns | 27.3 ns | 1.2x
        sinpi |  9.0 ns |     -   |  -
        cospi | 10.7 ns |     -   |  -
          erf | 21.0 ns |     -   |  -
         erfc | 22.9 ns |     -   |  -
        hypot |  8.4 ns | 13.9 ns | 1.7x
         powf | 24.2 ns |  3.2 ns | 0.1x
    remainder | 12.7 ns |     -   |  -
```

```
Throughput (independent array evals over [f32; 4096], examples/quickbench.rs; lower is better)
              | jodie    | std     | improvement
--------------|----------|---------|------------
         cbrt | 0.37 ns  | 3.93 ns | 10.6x
cbrt_accurate | 0.67 ns  | 3.93 ns | 5.9x
          cos | 0.27 ns  | 3.41 ns | 12.6x
  cos_checked | 1.61 ns  | 3.41 ns | 2.1x
         exp2 | 0.23 ns  | 3.13 ns | 13.7x
 exp2_checked | 0.48 ns  | 3.13 ns | 6.6x
         log2 | 0.53 ns  | 4.10 ns | 7.8x
          sin | 0.22 ns  | 3.10 ns | 14.4x
  sin_checked | 1.49 ns  | 3.10 ns | 2.1x
           ln | 0.55 ns  | 3.82 ns | 7.0x
        log10 | 0.53 ns  | 5.30 ns | 10.1x
        log1p | 0.60 ns  | 6.28 ns | 10.4x
          exp | 0.31 ns  | 2.47 ns | 8.0x
        expm1 | 0.46 ns  | 5.56 ns | 12.1x
         sinh | 0.58 ns  | 5.62 ns | 9.7x
         cosh | 0.50 ns  | 5.55 ns | 11.2x
         tanh | 0.64 ns  | 5.23 ns | 8.2x
        asinh | 2.19 ns  | 41.42 ns | 18.9x
        acosh | 1.32 ns  |  6.17 ns | 4.7x
        atanh | 0.79 ns  |  4.54 ns | 5.7x
         asin | 0.46 ns  |  4.01 ns | 8.7x
         acos | 0.26 ns  |  4.01 ns | 15.6x
         atan | 0.33 ns  |  6.32 ns | 19.2x
        atan2 | 0.35 ns  |  9.97 ns | 28.2x
          tan | 0.76 ns  |  9.07 ns | 11.9x
        sinpi | 0.18 ns  |     -    |  -
        cospi | 0.25 ns  |     -    |  -
          erf | 0.70 ns  |     -    |  -
         erfc | 0.68 ns  |     -    |  -
        hypot | 0.24 ns  |  2.92 ns | 12.0x
         powf | 0.95 ns  |  0.07 ns | 0.07x
    remainder | 0.22 ns  |     -    |  -
```

```
theoretical cost from llvm-mca (-mcpu=native, 100 iterations)
                    | latency (cyc)  | throughput (cyc)
--------------------|----------------|------------------
cbrt                |          35.06 |             1.629
cbrt_accurate       |          59.06 |             3.129
exp2                |          35.00 |             0.841
exp2_checked        |          43.06 |             1.399
log2                |          34.23 |             1.556
sin                 |          46.00 |             1.151
sin_checked         |         109.02 |             5.476
cos                 |          54.00 |             1.406
cos_checked         |         113.00 |             4.537
sinpi               |          38.00 |             0.925
cospi               |          46.00 |             1.150
ln                  |          56.91 |             1.626
log10               |          56.91 |             1.626
log1p               |          61.16 |             2.276
exp                 |          42.00 |             1.327
expm1               |          74.00 |             1.779
sinh                |          58.00 |             2.523
cosh                |          57.00 |             2.074
sinh_throughput     |          62.00 |             1.943
cosh_throughput     |          61.00 |             1.616
tanh                |          94.91 |             2.567
asinh               |         120.99 |             7.716
acosh               |         127.75 |             6.588
atanh               |          79.99 |             4.331
asin                |          59.03 |             0.968
acos                |          37.11 |             0.820
atan                |          61.09 |             1.491
atan2               |          61.17 |             1.532
tan                 |          71.02 |             2.532
erf                 |          91.74 |             2.871
erfc                |          78.09 |             2.599
hypot               |          21.11 |             0.766
powf                |          98.03 |             3.898
remainder           |          33.02 |             0.647
```

# tools
- `cargo +nightly run --release --example accuracy [thorough] [filter]` - avg/max ulp against an f64
  reference. Requires nightly: the reference is computed via the `sleef` crate's SIMD functions (cheapest
  ULP bucket available per function, u35 where it exists -- still ~1e8x tighter than f32 needs), which
  depends on the unstable `portable_simd` feature -- this also means `cargo test`/
  `cargo build --tests` now need nightly, since Cargo builds all dev-dependencies together regardless of
  which target you're building. Default mode fuzzes 100M random f32 bit patterns per function (a few
  seconds); `thorough` exhaustively sweeps all 2^32 bit patterns instead (every denormal, every NaN
  payload, both signs -- a few minutes). Runs on half the machine's cores at low OS scheduling priority
  (`nice`) so it doesn't compete with foreground work while iterating; refuses to run in a debug build.
- `cargo run --release --example quickbench [filter]` - latency (serial dependency chain) + throughput, min of 7 reps
- `cargo run --release --example edgecheck` - bit-exact checks of edge cases (0, -0, denormals, inf, nan, domain boundaries)
- `cargo run --release --example tune` - coordinate-descent ulp tuning of polynomial coefficients
- `cargo run --release --example mca` - theoretical latency/throughput straight from llvm-mca's scheduler
  model for the host CPU (requires `llvm-mca` on PATH). No wall-clock timing, so no thermal-throttling
  noise, and much faster to iterate on than quickbench; see examples/mca_target.rs for the marker
  functions it analyzes and why each region is built the way it is (llvm-mca has no branch predictor,
  so branchy edge-case handling has to be routed around, not just measured through).

# todo:
- do principled and thourough analysis of dependency chains and rounding errors to find optimizations
- perfectly rounded versions
- vary both arguments in quickbench's two-argument benchmarks (atan2, hypot,
  powf, remainder currently fix one argument, which may be letting LLVM
  constant-fold std's side of a couple of comparisons -- see the benchmark
  notes above)
