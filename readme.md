# faster? f32 math functions
Attempting to provide faster implementations of common f32 math functions with a similar level of accuracy to the standard library.
There are also perfectly rounded variants and faster-but-sloppier variants.

All functions auto-vectorize, it's a hard requirement

# precision (see examples/accuracy.rs)

```
                        | jodie avg  | jodie max | std avg | std max
------------------------|------------|-----------|---------|--------
                   cbrt |    0.281   |     3     |    0    |    0
     cbrt_unchecked (+) |    0.282   |     3     | (bit-identical to cbrt on its domain)
          cbrt_accurate |    0.000   |     1     |    0    |    0
cbrt_accurate_unchecked (+) | 0.000 |     1     | (bit-identical to cbrt_accurate on its domain)
                  rcbrt |    0.418   |     5     | (no std rcbrt)
                   exp2 |    0.026   |     1     |  0.000  |    1
           exp2_checked |    0.014   |     1     |  0.000  |    1
                  exp10 |    0.031   |     2     | (no std exp10)
          exp10_checked |    0.008   |     1     | (no std exp10)
                   log2 |    0.003   |     3     |  0.000  |    1
   log2_unchecked (+)   |    0.006   |     3     | (bit-identical to log2 on its domain)
        sin (|x|<1.3e7) |    0.059   |    219    |  0.003  |    1
        cos (|x|<1.3e7) |    0.289   |   2769    |  0.002  |    1
 sin_checked (|x|<=1e6) |    0.036   |     2     |  0.000  |    1
 cos_checked (|x|<=1e6) |    0.081   |     3     |  0.000  |    1
        sinpi (all f32) |    0.197   |     2     | (no std sinpi)
        cospi (all f32) |    0.281   | 8.7e8 (near a zero of cospi -- tiny absolute error, huge ulp) | (no std cospi)
        tanpi (all f32) |    0.275   | 3.6e6 (near a pole of tanpi -- tiny denominator, huge ulp) | (no std tanpi)
        sind (|x|<4.7e7)|    0.124   |     2     | (no std sind)
        cosd (|x|<4.7e7)|    0.073   |     2     | (no std cosd)
     tand (|x|<4.7e7)|    0.177   |     3     | (no std tand)
         sinc (|x|<1e6)|    0.094   |     4     | (no std sinc)
```

```
                        | jodie avg  | jodie max | std avg | std max
------------------------|------------|-----------|---------|--------
                     ln |    0.117   |     3     |  0.000  |    1
    ln_unchecked (+)    |    0.235   |     3     | (bit-identical to ln on its domain)
                  log10 |    0.127   |     3     |  0.000  |    0
  log10_unchecked (+)   |    0.255   |     3     | (bit-identical to log10 on its domain)
                  log1p |    0.097   |     4     |  0.000  |    0
                log2p1  |    0.102   |     3     | (no std log2p1)
        exp (in-domain) |    0.075   |     3     |  0.000  |    1
            exp_checked |    0.039   |     3     |  0.000  |    1
      expm1 (in-domain) |    0.130   |     6     |  0.000  |    0
          expm1_checked |    0.068   |     6     |  0.000  |    0
exp_m1_over_x (in-domain)| 0.073   |     6     | (no std exp_m1_over_x)
                exp2m1  |    0.077   |     4     | (no std exp2m1)
       sinh (in-domain) |    0.082   |     5     |  0.000  |    1
       cosh (in-domain) |    0.059   |     5     |  0.000  |    0
sinh_throughput (in-domain)| 0.082   |     5     |  0.000  |    1
cosh_throughput (in-domain)| 0.051   |     4     |  0.000  |    0
       sinh_checked (all f32) | 0.043 |     5     |  0.000  |    1
       cosh_checked (all f32) | 0.031 |     5     |  0.000  |    0
       tanh (in-domain) |    0.146   |     6     |  0.000  |    0
                sigmoid |    0.093   |     4     | (no std sigmoid)
       softplus (|x|<80)|    0.077   |     4     | (no std softplus)
logaddexp (|a|,|b|<80)|    0.141   |  ~1e3-1e5, heavy-tailed (real, narrow cancellation -- see its own doc comment) | (no std logaddexp)
                  asinh |    0.149   |     3     | 
                  acosh |    0.060   |     4     |  0.000  |    1
                  atanh |    0.004   |     2     | 
                   asin |    0.020   |     6     |  0.000  |    0
                   acos |    0.056   |     4     |  0.000  |    0
                   atan |    0.067   |     4     |  0.000  |    0
           atan_latency |    0.052   |     3     |  0.000  |    0
        tan (in-domain) |    0.331   |  2967     |  0.000  |    0
                   erf  |    0.318   |     4     | (no std erf)
         erfc (|x|<=10) |    0.311   |   109     | (no std erfc)
        erfcx (|x|<=10) |    0.377   |   125     | (no std erfcx)
                  atan2 |    0.069   |     3     |  0.000  |    0
    atan2_unchecked (+) |    0.069   |     3     | (bit-identical to atan2 on its domain)
        hypot (bounded) |    0.034   |     1     |  0.000  |    0
hypot_unchecked (bounded, +) | 0.034 |     1     | (bit-identical to hypot on its domain)
       hypot_checked |    0.015   |     1     | (no std comparison needed, no domain restriction)
                rhypot |    0.065   |     2     | (no std rhypot)
                  rsqrt |    0.260   |     1     | (no std rsqrt)
          pown (|n|<=8) |    0.158   |    11     | (no std pown)
         pown (|n|<=64) |    0.209   |    90     | (no std pown)
   pown_small (|n|<=255)|    0.216   |   304     | (bit-identical to pown on its domain)
        powf (in-domain)|    0.181   |  >=312    |  0.000  |    1
    powf_unchecked (+)  |    0.363   |  >=312    | (bit-identical to powf on its domain)
        powf_checked (in-domain)|  0.046   |  >=203    | (no std comparison needed, same domain as powf)
powf_checked_unchecked (+)|    0.046   |  >=203    | (bit-identical to powf_checked on its domain)
        remainder (|x/y|<1000, near-tie excluded) | 0.000 | 0 | (no std comparison needed)
remainder_unchecked (+) |    0.000   |     0     | (bit-identical to remainder on its domain)
    remainder_checked (|x/y|<1e7, near-tie excluded) | 0.000 | 0 | (no std comparison needed)
       remainder_wide (|x/y|<2e14, near-tie excluded) | 0.0003 | 4 | (no std comparison needed)
      remainder_ieee (|x/y|<1000, near-tie excluded) | 0.000 | 0 | (no std comparison needed)
                   fmod (|x/y|<1000, near-int excluded) | 0.000 | 0 | (matches Rust's `%`)
        fmod_unchecked (+) |    0.000   |     0     | (bit-identical to fmod on its domain)
```

# benchmarks
Run on i5-1145G7, -C target-cpu=native (now set in .cargo/config.toml)
```
Serial latency (dependency chain, examples/quickbench.rs; lower is better)
              | jodie   | std     | improvement
--------------|---------|---------|------------
         cbrt | 12.9 ns | 22.0 ns | 1.7x
cbrt_unchecked| 10.1 ns | 21.2 ns | 2.1x
cbrt_accurate | 17.0 ns | 22.0 ns | 1.3x
cbrt_accurate_unchecked| 18.7 ns | 27.2 ns | 1.5x
        rcbrt | 15.5 ns |     -   |  -
          cos | 12.7 ns | 12.8 ns | 1.0x
  cos_checked | 30.3 ns | 12.8 ns | 0.4x
         exp2 |  8.5 ns | 13.3 ns | 1.6x
 exp2_checked | 13.3 ns | 13.3 ns | 1.0x
        exp10 | 12.1 ns |     -   |  -
exp10_checked | 17.5 ns |     -   |  -
         log2 | 13.1 ns | 14.9 ns | 1.1x
log2_unchecked |  8.5 ns |     -   |  -
          sin | 10.9 ns | 12.9 ns | 1.2x
  sin_checked | 29.2 ns | 12.9 ns | 0.4x
           ln | 15.3 ns | 14.5 ns | 0.9x
  ln_unchecked |  8.8 ns |     -   |  -
        log10 | 14.3 ns | 16.6 ns | 1.2x
log10_unchecked|  8.8 ns |     -   |  -
        log1p | 16.5 ns | 21.1 ns | 1.3x
       log2p1 | 12.8 ns |     -   |  -
          exp |  9.9 ns |  9.5 ns | 1.0x
  exp_checked | 12.0 ns |  9.5 ns | 0.8x
        expm1 | 11.3 ns | 13.7 ns | 1.2x
exp_m1_over_x | 12.8 ns |     -   |  -
       exp2m1 | 11.5 ns |     -   |  -
         sinh | 14.0 ns | 14.0 ns | 1.0x
         cosh | 13.8 ns | 14.0 ns | 1.0x
 sinh_checked | 21.5 ns | 14.0 ns | 0.7x
 cosh_checked | 20.1 ns | 14.0 ns | 0.7x
         tanh | 17.7 ns | 16.5 ns | 0.9x
      sigmoid | 15.4 ns |     -   |  -
     softplus | 24.7 ns |     -   |  -
   logaddexp | 24.5 ns |     -   |  -
        asinh | 43.7 ns | 86.3 ns | 2.0x
        acosh | 35.1 ns | 27.0 ns | 0.8x
        atanh | 23.2 ns |  4.2 ns | 0.2x
         asin | 19.4 ns |  4.1 ns | 0.2x
         acos | 13.0 ns |  4.0 ns | 0.3x
         atan | 14.0 ns | 19.3 ns | 1.4x
 atan_latency | 14.6 ns |     -   |  -
      atan2 (*) | 23.6 ns | 29.9 ns | 1.3x
atan2_unchecked | 18.1 ns |     -   |  -
          tan | 23.5 ns | 27.3 ns | 1.2x
        sinpi |  9.0 ns |     -   |  -
        cospi | 10.7 ns |     -   |  -
        tanpi | 16.0 ns |     -   |  -
         sind | 10.8 ns |     -   |  -
         cosd | 12.6 ns |     -   |  -
         tand | 17.3 ns |     -   |  -
          erf | 21.0 ns |     -   |  -
         erfc | 21.6 ns |     -   |  -
        erfcx | 13.4 ns |     -   |  -
      hypot (*) |  5.1 ns | 13.9 ns | 2.7x
hypot_unchecked |  5.0 ns |     -   |  -
  hypot_checked | 19.9 ns | 11.9 ns | 0.6x
       rhypot |  7.7 ns |     -   |  -
        rsqrt |  6.7 ns |     -   |  -
         pown | 50.3 ns |     -   |  -
   pown_small | 15.3 ns |     -   |  -
 pown_const<N>|  4.2 ns |     -   |  -
      powf (*)| 22.2 ns | 16.9 ns | 0.8x
powf_unchecked | 18.7 ns | 16.8 ns | 0.9x
 powf_checked | 34.4 ns | 16.9 ns | 0.5x
powf_checked_unchecked| 38.1 ns | 17.1 ns | 0.4x
 remainder (*)|  8.6 ns |     -   |  -
remainder_unchecked|  8.0 ns |     -   |  -
remainder_checked| 12.3 ns |     -   |  -
   remainder_ieee|  7.1 ns |     -   |  -
   remainder_wide| 45.4 ns |     -   |  -
         fmod|  7.6 ns |     -   |  -
fmod_unchecked|  6.6 ns |     -   |  -
```
(*) atan2/hypot/powf/remainder's jodie numbers jumped here vs. older
recordings of this table -- not a regression, a benchmark fix: a literal
`1.0`/`2.0`/`3.0` 2nd argument let LLVM fold away their own special-case
branches at compile time, silently hiding the real cost (see readme's own
former todo note on this). `std powf`'s number changed even more
dramatically for a related but distinct reason, confirmed via
`--emit=asm`: `x.powf(2.0)` with a literal exponent isn't even a real
`powf` call -- LLVM recognizes the constant integer exponent and replaces
the whole thing with a single `x*x` multiply, so the old "3.2 ns" was
never measuring std's actual powf cost at all. All four functions (plus
`std powf`) now use a `black_box`'d 2nd argument for an honest number.

```
Throughput (independent array evals over [f32; 4096], examples/quickbench.rs; lower is better)
              | jodie    | std     | improvement
--------------|----------|---------|------------
         cbrt | 0.37 ns  | 3.93 ns | 10.6x
cbrt_unchecked| 0.25 ns  | 4.04 ns | 16.2x
cbrt_accurate | 0.67 ns  | 3.93 ns | 5.9x
cbrt_accurate_unchecked| 0.55 ns  | 5.27 ns | 9.6x
        rcbrt | 0.39 ns  |     -    |  -
          cos | 0.27 ns  | 3.41 ns | 12.6x
  cos_checked | 1.61 ns  | 3.41 ns | 2.1x
         exp2 | 0.23 ns  | 3.13 ns | 13.7x
 exp2_checked | 0.48 ns  | 3.13 ns | 6.6x
        exp10 | 0.33 ns  |     -   |  -
exp10_checked | 0.57 ns  |     -   |  -
         log2 | 0.53 ns  | 4.10 ns | 7.8x
log2_unchecked | 0.23 ns  |     -    |  -
          sin | 0.22 ns  | 3.10 ns | 14.4x
  sin_checked | 1.49 ns  | 3.10 ns | 2.1x
           ln | 0.55 ns  | 3.82 ns | 7.0x
  ln_unchecked | 0.26 ns  |     -    |  -
        log10 | 0.53 ns  | 5.30 ns | 10.1x
log10_unchecked| 0.26 ns  |     -    |  -
        log1p | 0.60 ns  | 6.28 ns | 10.4x
       log2p1 | 0.54 ns  |     -    |  -
          exp | 0.31 ns  | 2.47 ns | 8.0x
  exp_checked | 0.39 ns  | 2.47 ns | 6.3x
        expm1 | 0.46 ns  | 5.56 ns | 12.1x
exp_m1_over_x | 0.44 ns  |     -    |  -
       exp2m1 | 0.49 ns  |     -    |  -
         sinh | 0.58 ns  | 5.62 ns | 9.7x
         cosh | 0.50 ns  | 5.55 ns | 11.2x
 sinh_checked | 0.81 ns  | 5.62 ns | 6.9x
 cosh_checked | 0.73 ns  | 5.55 ns | 7.6x
         tanh | 0.64 ns  | 5.23 ns | 8.2x
      sigmoid | 0.43 ns  |     -    |  -
     softplus | 1.33 ns  |     -    |  -
   logaddexp | 1.36 ns  |     -    |  -
        asinh | 2.19 ns  | 41.42 ns | 18.9x
        acosh | 1.32 ns  |  6.17 ns | 4.7x
        atanh | 0.79 ns  |  4.54 ns | 5.7x
         asin | 0.46 ns  |  4.01 ns | 8.7x
         acos | 0.26 ns  |  4.01 ns | 15.6x
         atan | 0.33 ns  |  6.32 ns | 19.2x
 atan_latency | 0.36 ns  |     -    |  -
  atan2 (*) | 0.65 ns  | 13.06 ns | 20.2x
atan2_unchecked | 0.48 ns  |     -    |  -
          tan | 0.76 ns  |  9.07 ns | 11.9x
        sinpi | 0.18 ns  |     -    |  -
        cospi | 0.25 ns  |     -    |  -
        tanpi | 0.65 ns  |     -    |  -
         sind | 0.24 ns  |     -    |  -
         cosd | 0.31 ns  |     -    |  -
         tand | 0.58 ns  |     -    |  -
          erf | 0.70 ns  |     -    |  -
         erfc | 0.84 ns  |     -    |  -
        erfcx | 0.78 ns  |     -    |  -
  hypot (*) | 0.17 ns  |  2.42 ns | 14.2x
hypot_unchecked | 0.17 ns  |     -    |  -
  hypot_checked | 0.40 ns  |  2.72 ns | 6.8x
       rhypot | 0.31 ns  |     -    |  -
        rsqrt | 0.31 ns  |     -    |  -
         pown | 1.24 ns  |     -    |  -
   pown_small | 0.22 ns  |     -    |  -
 pown_const<N>| 0.05 ns  |     -    |  -
      powf (*)| 1.16 ns  |  5.43 ns | 4.7x
powf_unchecked | 0.70 ns  |  5.86 ns | 8.3x
 powf_checked | 1.79 ns  |  5.43 ns | 3.0x
powf_checked_unchecked| 1.88 ns  |  6.32 ns | 3.4x
 remainder (*)| 0.16 ns  |     -    |  -
remainder_unchecked| 0.15 ns  |     -    |  -
remainder_checked| 0.29 ns  |     -    |  -
   remainder_ieee| 0.17 ns  |     -    |  -
   remainder_wide| 1.76 ns  |     -    |  -
         fmod| 0.15 ns  |     -    |  -
fmod_unchecked| 0.15 ns  |     -    |  -
```
(*) see the latency table's own footnote above -- same benchmark fix,
not a regression. `powf`'s throughput ratio flips especially hard here
(was "0.07x", now "5.7x") since the old std comparison point was a bare
`x*x` multiply, not real `powf`.

```
theoretical cost from llvm-mca (-mcpu=native, 100 iterations)
                    | latency (cyc)  | throughput (cyc)
--------------------|----------------|------------------
cbrt                |          35.06 |             1.629
cbrt_unchecked      |          35.06 |             0.906
cbrt_accurate       |          59.06 |             3.129
cbrt_accurate_unchecked |      59.06 |             2.067
rcbrt               |          46.30 |             1.666
exp2                |          35.00 |             0.841
exp2_checked        |          43.06 |             1.399
exp10               |          52.00 |             1.565
exp10_checked       |          51.06 |             1.897
log2                |          34.23 |             1.584
log2_unchecked      |          34.23 |             0.958
sin                 |          46.00 |             1.151
sin_checked         |         117.02 |             5.232
cos                 |          54.00 |             1.406
cos_checked         |         122.00 |             4.603
sinpi               |          42.02 |             1.133
cospi               |          51.00 |             1.283
tanpi               |          71.14 |             2.152
sinc                |          53.02 |             1.256
sind                |          46.00 |             1.151
cosd                |          54.00 |             1.406
tand                |          70.02 |             2.533
ln                  |          52.86 |             1.611
ln_unchecked        |          34.06 |             1.018
log10               |          56.86 |             1.635
log10_unchecked     |          38.22 |             1.113
log1p               |          47.24 |             1.857
log2p1              |          48.14 |             1.886
exp                 |          42.00 |             1.327
exp_checked         |          46.06 |             1.729
expm1               |          71.00 |             1.695
expm1_checked       |          75.06 |             1.595
exp_m1_over_x       |          83.00 |             1.798
exp2m1              |          76.06 |             1.843
sinh                |          52.00 |             1.949
cosh                |          51.00 |             1.770
sinh_throughput     |          62.00 |             1.943
cosh_throughput     |          61.00 |             1.616
sinh_checked        |          56.06 |             2.340
cosh_checked        |          55.06 |             2.051
tanh                |          82.72 |             1.859
sigmoid             |          57.06 |             1.283
softplus            |          73.14 |             2.534
logaddexp           |          73.14 |             2.534
asinh               |          81.02 |             4.169
acosh               |          89.02 |             3.569
atanh               |          69.42 |             2.903
asin                |          59.03 |             0.968
acos                |          37.11 |             0.820
atan                |          61.09 |             1.491
atan_latency        |          59.11 |             1.591
atan2               |          61.28 |             1.694
tan                 |          71.02 |             2.532
erf                 |          76.11 |             2.034
erfc                |          64.00 |             2.437
erfcx               |          39.36 |             2.278
hypot               |          21.11 |             0.766
hypot_checked       |          57.19 |             1.178
rhypot              |          32.02 |             1.389
rsqrt               |          28.00 |             1.381
pown                |         176.00 |             3.805
powf                |         102.99 |             5.651
powf_unchecked      |          79.05 |             3.095
powf_checked        |         130.33 |             9.105
powf_checked_unchecked |      129.74 |             7.234
remainder           |          34.11 |                 ? (*)
remainder_unchecked |          33.00 |             0.646
remainder_checked   |          46.24 |             1.544
remainder_ieee      |          29.11 |                 ? (*)
remainder_wide      |         165.30 |             8.186
fmod                |          29.11 |                 ? (*)
fmod_unchecked      |          28.00 |             0.643
```
(*) remainder/remainder_ieee/fmod: throughput no longer measurable via
llvm-mca after their 2026-07-09 zero/nan fix (backlog idea #85's own
follow-up) -- LLVM branch-specializes these short functions' vectorized
loop on the harness's own shared black_box'd `y`, producing multiple
physical exit paths that corrupt llvm-mca's region parser (same class
of harness limitation as `pown_small`'s own precedent; confirmed via a
standalone `--emit=asm` probe with genuinely per-lane-varying inputs
that the real functions still vectorize cleanly, no scalar fallback).
See the throughput table above for their real wall-clock numbers
instead.

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

Standing gates. Each exits nonzero on failure, so they can be run as a suite -- the whole set is about 20
seconds (measured: edgecheck 0.1s, worst_corpus 0.2s, special_matrix 0.4s, special_matrix2 0.4s,
saturation_pins 0.3s, eft_contract_check 2.9s, denormal_audit 3.5s, approx_bounds 12.6s). They cover
contracts the ulp sweeps above structurally cannot -- special values, clamp boundaries, denormal handling,
and documented-bound drift.

- `cargo run --release --example worst_corpus` - 108 public 1-arg functions x 90 historically-hard inputs
  (special values, every branch seam, every clamp boundary, recorded worst-x values) checked bit-identical
  against a blessed golden file, in ~0.06s. The fast counterpart to the hours-long exhaustive sweeps.
  `-- --bless` regenerates; an intentional accuracy change is *expected* to fail this, and the diff is
  meant to be eyeballed. Note a pass is not "nothing changed": see the file header for two measured cases
  that slip through (a 1-ulp coefficient nudge, and a seam move where both branches agree at the corpus point).
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
  each flush begins relative to the true zero, which is the metric that matters.
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
