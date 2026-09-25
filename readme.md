# faster? f32 math functions
Attempting to provide faster implementations of common f32 math functions with a similar level of accuracy to the standard library.
There are also perfectly rounded variants and faster-but-sloppier variants.

All functions auto-vectorize, it's a hard requirement.

# precision

```
                        | jodie avg  | jodie max | std avg | std max
------------------------|------------|-----------|---------|--------
                   cbrt |    0.281   |     3     |    0    |    0
         cbrt_unchecked |    0.282   |     3     | (bit-identical to cbrt on its domain)
          cbrt_accurate |    0.000   |     1     |    0    |    0
cbrt_accurate_unchecked |    0.000   |     1     | (bit-identical to cbrt_accurate on its domain)
                  rcbrt |    0.117   |     1     | (no std rcbrt)
                   exp2 |    0.026   |     1     |  0.000  |    1
           exp2_checked |    0.014   |     1     |  0.000  |    1
                  exp10 |    0.031   |     2     | (no std exp10)
          exp10_checked |    0.008   |     1     | (no std exp10)
                   log2 |    0.003   |     1     |  0.000  |    1
   log2_unchecked (+)   |    0.006   |     1     | (bit-identical to log2 on its domain)
       sin (|x|<2^24*pi) |    0.046   |     2     |  0.003  |    1
       cos (|x|<2^22*pi) |    0.083   |     2     |  0.002  |    1
      sin_wide (all f32)|    0.127   |     2     |  0.000  |    1
      cos_wide (all f32)|    0.150   |     2     |  0.000  |    1
        sinpi (all f32) |    0.197   |     2     | (no std sinpi)
        cospi (all f32) |    0.058   |     2     | (no std cospi)
        tanpi (all f32) |    0.033   |     2     | (no std tanpi)
```


```
                         | jodie avg  | jodie max | std avg | std max
-------------------------|------------|-----------|---------|--------
                      ln |    0.004   |     1     |  0.000  |    1
     ln_unchecked (+)    |    0.007   |     1     | (bit-identical to ln on its domain)
                   log10 |    0.003   |     1     |  0.000  |    0
   log10_unchecked (+)   |    0.007   |     1     | (bit-identical to log10 on its domain)
                   log1p |    0.025   |     2     |  0.000  |    1
                 log2p1  |    0.032   |     2     | 
                 log10p1  |    0.039   |     2     | 
          exp (in-domain) |    0.071   |     3     |  0.000  |    1
              exp_checked |    0.004   |     1     |  0.000  |    1
        expm1 (in-domain) |    0.008   |     2     |  0.000  |    0
            expm1_checked |    0.004   |     2     |  0.000  |    0
                  exp2m1  |    0.040   |     2     | (no std exp2m1)
                 exp10m1  |    0.095   |     3     |
         sinh (in-domain) |    0.061   |     3     |  0.000  |    1
         cosh (in-domain) |    0.024   |     2     |  0.000  |    0
sinh_throughput (in-domain)| 0.080    |     4     |  0.000  |    1
cosh_throughput (in-domain)| 0.049    |     3     |  0.000  |    0
        sinh_checked (all f32) | 0.032 |     3     |  0.000  |    1
        cosh_checked (all f32) | 0.013 |     2     |  0.000  |    0
         tanh (in-domain) |    0.044   |     2     |  0.000  |    0
                  sigmoid |    0.091   |     3     | (no std sigmoid)
       logsigmoid_checked |    0.039   |     3     | (no std logsigmoid)
   silu_checked (all f32)|    0.045   |     5     | (no std silu; 4.69 on a dense scan -- `silu` is 3.85, the avg/max trade in its doc)
    gelu (|x|<=10*sqrt2) |    0.131   |     6     | (no std gelu; exhaustive)
  gelu (all finite f32)  |    0.067   |     6     |
                  asinh |    0.034   |     2     | 
                  acosh |    0.003   |     3     |  0.000  |    1
                  atanh |    0.004   |     2     | 
                   asin |    0.016   |     2     |  0.000  |    0
                   acos |    0.002   |     2     |  0.000  |    0
                   atan |    0.063   |     3     |  0.000  |    0
           atan_latency |    0.052   |     3     |  0.000  |    0
      tan (|x|<2^22*pi) |    0.118   |     4     |  0.000  |    0
      tan_wide (all f32)|    0.250   |     4     |  0.000  |    0
                   erf  |    0.027   |     3     | (no std erf)
         erfc (|x|<=10) |    0.122   |     6     |
       erfinv (|x|<1)   |    0.047   |     3     | (no std erfinv; exhaustive. Almost all of the average is one number: the effective leading coefficient `1+c0` of the central poly, which is what every `|x|` under ~0.06 computes)
      erfc_inv (0<y<2)  |    0.380   |     4     | (no std erfc_inv; exhaustive. 81% of the domain's bit patterns take the far-tail branch, so its poly sets this)
                  atan2 |    0.066   |     4     |  0.000  |    0
    atan2_unchecked (+) |    0.066   |     4     | (bit-identical to atan2 on its domain)
              atan2_pos |    0.062   |     3     | (no std atan2_pos)
                  rsqrt |    0.260   |     1     | (no std rsqrt)
        powf (in-domain)|    0.001   |   1 (s)  |  0.000  |    1
    powf_unchecked (+)  |    0.001   |   1 (s)  | (bit-identical to powf on its domain)
                powf_pos|    0.001   |   1 (s)  | (bit-identical to powf for x > 0, see its doc comment)
                   fmod (|x/y|<1000, near-int excluded) | 0.000 | 0 | (matches Rust's `%`)
        fmod_unchecked (+) |    0.000   |     0     | (bit-identical to fmod on its domain)
```

# benchmarks
Run on i5-1145G7, -C target-cpu=native (now set in .cargo/config.toml)

```
Serial latency (dependency chain, examples/quickbench.rs; lower is better)
              | jodie   | std     | improvement
--------------|---------|---------|------------
         cbrt | 19.7 ns | 27.9 ns | 1.4x
cbrt_unchecked| 13.8 ns | 27.9 ns | 2.0x
cbrt_accurate | 21.8 ns | 27.9 ns | 1.3x
cbrt_accurate_unchecked| 19.7 ns | 27.9 ns | 1.4x
        rcbrt | 16.9 ns |    -    |  -
          cos | 20.1 ns | 18.8 ns | 0.9x
     cos_wide | 22.2 ns | 20.6 ns | 0.9x
         exp2 | 10.6 ns | 11.9 ns | 1.1x
 exp2_checked | 13.9 ns | 11.9 ns | 0.9x
        exp10 | 15.0 ns |    -    |  -
exp10_checked | 15.7 ns |    -    |  -
         log2 | 12.8 ns | 12.8 ns | 1.0x
log2_unchecked| 11.9 ns |    -    |  -
          sin | 19.5 ns | 16.6 ns | 0.8x
     sin_wide | 20.8 ns | 21.5 ns | 1.0x
           ln | 13.0 ns | 13.5 ns | 1.0x
 ln_unchecked | 12.2 ns |    -    |  -
        log10 | 13.2 ns | 15.2 ns | 1.1x
log10_unchecked| 11.8 ns |    -    |  -
        log1p | 16.6 ns | 18.5 ns | 1.1x
       log2p1 | 15.6 ns |    -    |  -
          exp | 12.2 ns | 11.7 ns | 1.0x
  exp_checked | 15.3 ns | 11.7 ns | 0.8x
        expm1 | 13.4 ns | 16.7 ns | 1.2x
       exp2m1 | 14.8 ns |    -    |  -
         sinh | 16.6 ns | 17.6 ns | 1.1x
         cosh | 16.3 ns | 18.1 ns | 1.1x
 sinh_checked | 18.5 ns | 17.6 ns | 1.0x
 cosh_checked | 17.6 ns | 18.1 ns | 1.0x
         tanh | 16.9 ns | 18.6 ns | 1.1x
      sigmoid | 18.0 ns |    -    |  -
        asinh | 27.9 ns | 27.4 ns | 1.0x
        acosh | 22.9 ns | 25.5 ns | 1.1x
        atanh | 23.1 ns | 24.4 ns | 1.1x
         asin | 12.6 ns | 17.3 ns | 1.4x
         acos | 12.3 ns | 18.4 ns | 1.5x
         atan | 17.4 ns | 23.5 ns | 1.4x
 atan_latency | 18.3 ns |    -    |  -
        atan2 | 21.7 ns | 28.6 ns | 1.3x
atan2_unchecked| 21.7 ns |    -    |  -
          tan | 21.6 ns | 23.5 ns | 1.1x
     tan_wide | 26.5 ns | 32.6 ns | 1.2x
        sinpi | 15.5 ns |    -    |  -
        cospi | 16.1 ns |    -    |  -
        tanpi | 18.5 ns |    -    |  -
          erf | 17.1 ns |    -    |  -
         erfc | 19.2 ns |    -    |  -
        rsqrt |  8.6 ns |    -    |  -
         powf | 34.5 ns | 24.0 ns | 0.7x
powf_unchecked| 32.3 ns | 24.0 ns | 0.7x
         fmod | 11.8 ns |    -    |  -
fmod_unchecked| 10.9 ns |    -    |  -
```

```
Throughput (independent array evals over [f32; 4096], examples/quickbench.rs; lower is better)
              | jodie    | std     | improvement
--------------|----------|---------|------------
         cbrt | 0.40 ns | 5.21 ns | 12.9x
cbrt_unchecked| 0.29 ns | 5.21 ns | 17.9x
cbrt_accurate | 0.62 ns | 5.21 ns | 8.4x
cbrt_accurate_unchecked| 0.45 ns | 5.21 ns | 11.6x
        rcbrt | 0.31 ns |    -    |  -
          cos | 0.34 ns | 4.79 ns | 14.1x
     cos_wide | 2.93 ns | 8.40 ns | 2.9x
         exp2 | 0.24 ns | 3.05 ns | 12.7x
 exp2_checked | 0.37 ns | 3.05 ns | 8.2x
        exp10 | 0.31 ns |    -    |  -
exp10_checked | 0.40 ns |    -    |  -
         log2 | 0.35 ns | 3.72 ns | 10.5x
log2_unchecked| 0.29 ns |    -    |  -
          sin | 0.35 ns | 4.55 ns | 13.0x
     sin_wide | 2.84 ns | 7.00 ns | 2.5x
           ln | 0.35 ns | 3.79 ns | 10.8x
 ln_unchecked | 0.33 ns |    -    |  -
        log10 | 0.35 ns | 5.11 ns | 14.4x
log10_unchecked| 0.29 ns |    -    |  -
        log1p | 0.43 ns | 5.96 ns | 13.7x
       log2p1 | 0.47 ns |    -    |  -
          exp | 0.32 ns | 3.10 ns | 9.6x
  exp_checked | 0.35 ns | 3.10 ns | 8.9x
        expm1 | 0.26 ns | 6.30 ns | 24.0x
       exp2m1 | 0.30 ns |    -    |  -
         sinh | 0.48 ns | 7.15 ns | 14.8x
         cosh | 0.42 ns | 7.10 ns | 17.0x
 sinh_checked | 0.49 ns | 7.15 ns | 14.5x
 cosh_checked | 0.43 ns | 7.10 ns | 16.6x
         tanh | 0.42 ns | 6.12 ns | 14.6x
      sigmoid | 0.36 ns |    -    |  -
        asinh | 1.08 ns | 8.18 ns | 7.6x
        acosh | 0.92 ns | 8.19 ns | 8.9x
        atanh | 0.68 ns | 7.75 ns | 11.4x
         asin | 0.34 ns | 7.14 ns | 20.9x
         acos | 0.27 ns | 6.96 ns | 26.3x
         atan | 0.39 ns | 7.71 ns | 20.0x
 atan_latency | 0.33 ns |    -    |  -
        atan2 | 0.56 ns | 11.85 ns | 21.4x
atan2_unchecked| 0.55 ns |    -    |  -
          tan | 0.60 ns | 7.54 ns | 12.6x
     tan_wide | 3.70 ns | 14.48 ns | 3.9x
        sinpi | 0.37 ns |    -    |  -
        cospi | 0.38 ns |    -    |  -
        tanpi | 0.79 ns |    -    |  -
          erf | 0.55 ns |    -    |  -
         erfc | 0.78 ns |    -    |  -
        rsqrt | 0.41 ns |    -    |  -
         powf | 1.73 ns | 8.83 ns | 5.1x
powf_unchecked| 1.64 ns | 8.83 ns | 5.4x
         fmod | 0.29 ns |    -    |  -
fmod_unchecked| 0.23 ns |    -    |  -
```

```
theoretical cost from llvm-mca (-mcpu=native, 100 iterations)
                    | latency (cyc)  | throughput (cyc)
--------------------|----------------|------------------
cbrt                |          35.06 |             1.629
cbrt_unchecked      |          35.06 |             0.906
cbrt_accurate       |          63.00 |             3.129
cbrt_accurate_unchecked |      63.00 |             2.067
rcbrt               |          60.95 |             1.690
exp2                |          35.00 |             0.854
exp2_checked        |          47.00 |             1.399
exp10               |          52.00 |             1.461
exp10_checked       |          55.00 |             1.657
log2                |          38.06 |             1.583
log2_unchecked      |          38.06 |             1.021
sin                 |          64.00 |             1.778
sin_wide            |          57.83 |             3.794
cos                 |          61.00 |             1.654
cos_wide            |          60.03 |             4.041
tan_wide            |          68.05 |             4.553
sinpi               |          42.02 |             1.133
cospi               |          47.00 |             1.226
tanpi               |          77.05 |             2.155
ln                  |          49.13 |             1.625
ln_unchecked        |          38.06 |             1.021
log10               |          48.16 |             1.617
log10_unchecked     |          38.06 |             1.022
log1p               |          51.25 |             1.903
log2p1              |          51.36 |             1.876
exp                 |          42.00 |             1.195
exp_checked         |          54.00 |             1.526
expm1               |          47.00 |             1.087
expm1_checked       |          55.00 |             1.320
exp2m1              |          52.00 |             1.278
exp10m1             |          55.00 |             1.342
sinh                |          54.00 |             1.824
cosh                |          53.00 |             1.695
sinh_throughput     |          62.00 |             1.689
cosh_throughput     |          61.00 |             1.498
sinh_checked        |          62.00 |             2.094
cosh_checked        |          61.00 |             2.086
tanh                |          63.00 |             1.759
sigmoid             |          61.00 |             1.222
logsigmoid          |          75.11 |             2.604
logsigmoid_checked  |          75.56 |             2.699
silu                |          65.02 |             1.407
silu_checked        |          60.97 |             1.466
gelu                |          71.88 |             3.434
asinh               |          74.48 |             4.159
acosh               |          97.88 |             3.948
atanh               |         100.83 |             2.974
asin                |          60.91 |             0.961
acos                |          44.02 |             0.961
atan                |          61.27 |             1.491
atan_latency        |          61.99 |             1.591
atan2               |          67.19 |             1.694
tan                 |          78.00 |             2.985
erf                 |          87.00 |             2.040
erfc                |          63.03 |             3.175
rsqrt               |          28.00 |             1.381
powf                |         125.81 |             6.996
powf_unchecked      |         118.00 |             6.171
fmod                |          35.13 |                 ? (*)
fmod_unchecked      |          28.00 |             0.643
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
- `cargo run --release --example edgecheck` - bit-exact checks of edge cases (0, -0, denormals, inf, nan, domain boundaries); exits nonzero if any pin fails
- `cargo run --release --example tune` - coordinate-descent ulp tuning of polynomial coefficients
- `cargo run --release --example mca` - theoretical latency/throughput straight from llvm-mca's scheduler
  model for the host CPU (requires `llvm-mca` on PATH). No wall-clock timing, so no thermal-throttling
  noise, and much faster to iterate on than quickbench; see examples/mca_target.rs for the marker
  functions it analyzes and why each region is built the way it is (llvm-mca has no branch predictor,
  so branchy edge-case handling has to be routed around, not just measured through).

Standing gates. Each exits nonzero on failure, so they can be run as a suite -- the whole set is about 20
seconds (measured: edgecheck 0.1s, worst_corpus 0.2s, special_matrix 0.4s, special_matrix2 0.4s,
saturation_pins 0.3s, eft_contract_check 2.9s, denormal_audit 3.5s, approx_bounds 12.6s). They cover
contracts the ulp sweeps above structurally cannot -- special values, clamp boundaries, denormal handling,
and documented-bound drift.

- `cargo run --release --example worst_corpus` - 92 public 1-arg functions x 90 historically-hard inputs
  (special values, every branch seam, every clamp boundary, recorded worst-x values) checked bit-identical
  against a blessed golden file, in ~0.06s. The fast counterpart to the hours-long exhaustive sweeps.
  `-- --bless` regenerates; an intentional accuracy change is *expected* to fail this, and the diff is
  meant to be eyeballed. Note a pass is not "nothing changed": see the file header for two measured cases
  that slip through (a 1-ulp coefficient nudge, and a seam move where both branches agree at the corpus point).
- `cargo run --release --example powfsearch` - adversarial worst case for `powf`/`powf_checked`, where a
  blind 2-arg fuzz is structurally weak: their error is `ln2 * |y*log2(x)| * relerr(log2)`, maximised only
  where both factors are extreme at once, so this derives `y` from `x` to pin the first and sweeps the
  octave that maximises the second exhaustively instead of sampling. Finds twice the max ulp the fuzz does.
- `cargo run --release --example clogsearch` - the same idea for `clog`'s real part, where the blind 2-arg
  fuzz is weak in *both* halves. `Re clog = ln|z|` is a near-total cancellation on the unit circle, which
  independent uniform `re`/`im` never reach -- and `accuracy.rs`'s own reference there is
  `(re*re + im*im).sqrt().ln()`, whose relative error is `2^-53/|re^2+im^2-1|`, i.e. already ~2 f32 ulp at `|v| ~ 1e-9`
  and unbounded below that, so it could not score the region even if it sampled it. This walks the manifold
  and builds `re^2+im^2` exactly in f64 instead. Found **4096 ulp** where the standing row read 3. Exits
  nonzero above a 4 ulp budget.
- `cargo run --release --example special_matrix` - +-0/+-inf/NaN in/out matrix over every public 1-arg
  function, asserting NaN propagation and quietness (the `_unchecked`/`_approx` tiers that promise nothing
  off-domain are exempted by name, so a *new* function inheriting garbage NaN behaviour still fails).
  Found the `wrap_pi(-0.0)` sign bug.
- `cargo run --release --example special_matrix2` - the same for 2-arg functions, over the full 9x9
  +-0/+-inf/NaN cross product, comparing bitwise against f64 std where a counterpart exists.
- `cargo run --release --example saturation_pins` - every input clamp and overflow threshold swept +-64 ulp
  against an f64 reference, plus far-outside pins asserting `f(+-1e30)` equals the mathematical limit.
  Generalizes the exp10_checked overflow-at-the-boundary bug.
- `cargo run --release --example denormal_audit` - which functions carry denormal outputs correctly vs
  flush early, split into normal-input-denormal-output and denormal-in-denormal-out. Reports how *early*
  each flush begins relative to the true zero, which is the metric that matters. Both of its lists are
  hand-maintained, and both were incomplete: `logsigmoid`/`silu` flush their *entire* denormal
  range (16.97 and 20.28 premature in `x`, all worse than the `sigmoid` case that prompted this
  example) and neither of them was listed. Adding a function here is not optional bookkeeping -- nothing else
  in the suite reports this. Both now have a `_checked` sibling that carries the band ("never" flushes,
  0 ulp at every integer `x` through it); the defaults keep their throughput.
- `cargo run --release --example approx_bounds` - asserts the `_approx` tier's doc-comment error bounds
  (relative/absolute, not ulp -- that tier is deliberately outside the 0.5/2 budget).
- `cargo run --release --example eft_contract_check` - the public EFT toolkit's exactness contracts against
  an f64 reference, including exhaustive `mulsign` over 8 y-values x all 2^32 x.
- `cargo run --release --example nan_payload` - record rather than a gate: whether each function preserves an
  input NaN's payload and sign, canonicalizes it, or (for the off-domain tiers) returns no NaN at all. IEEE754
  permits either of the first two; `special_matrix` asserts the part that *is* required.
- `cargo run --release --example error_profile` - diagnostic rather than a gate: per-function ulp histogram
  and per-magnitude-band avg/max, annotated by which side of the seam each band falls on. Distinguishes
  refit / seam-move / new-sub-branch, which an avg+max pair cannot.
- `cargo run --release --example codegen_check` - greps each `*_throughput` asm region for the specific
  de-vectorization signatures this crate has actually been bitten by (scalar `call`, scalar float->int
  converts, scalar divide/sqrt), and that at least one packed SIMD arithmetic instruction is present so an
  optimized-away region can't pass vacuously.

# todo:
- do principled and thourough analysis of dependency chains and rounding errors to find optimizations
- perfectly rounded versions
