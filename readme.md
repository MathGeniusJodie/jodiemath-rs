# faster? f32 math functions
Attempting to provide faster implementations of common f32 math functions with a similar level of accuracy to the standard library.
There are also perfectly rounded variants and faster-but-sloppier variants.

All functions auto-vectorize, it's a hard requirement

# precision (see examples/accuracy.rs)

```
                        | jodie avg  | jodie max | std avg | std max
------------------------|------------|-----------|---------|--------
                   cbrt |    0.281   |     3     |    0    |    0
         cbrt_unchecked |    0.282   |     3     | (bit-identical to cbrt on its domain)
          cbrt_accurate |    0.000   |     1     |    0    |    0
cbrt_accurate_unchecked |    0.000   |     1     | (bit-identical to cbrt_accurate on its domain)
                  rcbrt |    0.418   |     5     | (no std rcbrt)
                   exp2 |    0.026   |     1     |  0.000  |    1
           exp2_checked |    0.014   |     1     |  0.000  |    1
                  exp10 |    0.031   |     2     | (no std exp10)
          exp10_checked |    0.008   |     1     | (no std exp10)
                   log2 |    0.003   |     3     |  0.000  |    1
   log2_unchecked (+)   |    0.006   |     3     | (bit-identical to log2 on its domain)
         sin (|x|<=1e6) |    0.036   |     3     |  0.002  |    1
         cos (|x|<=1e6) |    0.078   |     3     |  0.002  |    1
 sin_checked (|x|<=1e6) |    0.036   |     2     |  0.000  |    1
 cos_checked (|x|<=1e6) |    0.081   |     3     |  0.000  |    1
        sinpi (all f32) |    0.197   |     2     | (no std sinpi)
        cospi (all f32) |    0.281   | 8.7e8 (near a zero of cospi -- tiny absolute error, huge ulp) | (no std cospi)
        tanpi (all f32) |    0.386   | 3.6e6 (near a pole of tanpi -- tiny denominator, huge ulp) | (no std tanpi)
        sind (|x|<4.7e7)|    0.124   |     2     | (no std sind)
        cosd (|x|<4.7e7)|    0.073   |     2     | (no std cosd)
        tand (|x|<4.7e7)|    0.177   |     3     | (no std tand)
          sinc (|x|<1e6)|    0.094   |     4     | (no std sinc)
```

```
                         | jodie avg  | jodie max | std avg | std max
-------------------------|------------|-----------|---------|--------
                      ln |    0.117   |     3     |  0.000  |    1
     ln_unchecked (+)    |    0.235   |     3     | (bit-identical to ln on its domain)
                   log10 |    0.127   |     3     |  0.000  |    0
   log10_unchecked (+)   |    0.255   |     3     | (bit-identical to log10 on its domain)
                   log1p |    0.097   |     4     |  0.000  |    0
                 log2p1  |    0.102   |     3     | (no std log2p1)
         exp (in-domain) |    0.071   |     3     |  0.000  |    1
             exp_checked |    0.037   |     3     |  0.000  |    1
       expm1 (in-domain) |    0.129   |     5     |  0.000  |    0
           expm1_checked |    0.067   |     5     |  0.000  |    0
exp_m1_over_x (in-domain)|    0.071   |     5     | (no std exp_m1_over_x)
                 exp2m1  |    0.077   |     4     | (no std exp2m1)
        sinh (in-domain) |    0.078   |     4     |  0.000  |    1
        cosh (in-domain) |    0.052   |     4     |  0.000  |    0
sinh_throughput (in-domain)| 0.080    |     4     |  0.000  |    1
cosh_throughput (in-domain)| 0.049    |     3     |  0.000  |    0
       sinh_checked (all f32) | 0.041 |     4     |  0.000  |    1
       cosh_checked (all f32) | 0.027 |     4     |  0.000  |    0
        tanh (in-domain) |    0.145   |     5     |  0.000  |    0
                 sigmoid |    0.091   |     3     | (no std sigmoid)
        softplus (|x|<80)|    0.075   |     3     | (no std softplus)
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
        powf (in-domain)|    0.181   |  >=312 (s)|  0.000  |    1
    powf_unchecked (+)  |    0.363   |  >=312 (s)| (bit-identical to powf on its domain)
        powf_checked (in-domain)|  0.019   |   3 (s)  | (no std comparison needed, same domain as powf)
powf_checked_unchecked (+)|    0.038   |   3 (s)  | (bit-identical to powf_checked on its domain)
        remainder (|x/y|<1000, near-tie excluded) | 0.000 | 0 | (no std comparison needed)
remainder_unchecked (+) |    0.000   |     0     | (bit-identical to remainder on its domain)
    remainder_checked (|x/y|<1e7, near-tie excluded) | 0.000 | 0 | (no std comparison needed)
       remainder_wide (|x/y|<2e14, near-tie excluded) | 0.0003 | 4 | (no std comparison needed)
      remainder_ieee (|x/y|<1000, near-tie excluded) | 0.000 | 0 | (no std comparison needed)
                   fmod (|x/y|<1000, near-int excluded) | 0.000 | 0 | (matches Rust's `%`)
        fmod_unchecked (+) |    0.000   |     0     | (bit-identical to fmod on its domain)
```
(s) from `examples/powfsearch.rs`, not from the fuzz above. powf's error is
`ln2 * |y*log2(x)| * (relative error of the log2 feeding exp2)`, so its max
lives in one corner -- `|y*log2(x)|` just inside the largest finite
exponent, `x` inside the single octave where log2's own relative error is
undiluted by the integer exponent -- which random pairs essentially never
hit. The blind fuzz reports `powf_checked` at 2 ulp; pinning `y` to `x`
and sweeping that octave exhaustively finds 4.

# benchmarks
Run on i5-1145G7, -C target-cpu=native (now set in .cargo/config.toml)

Both wall-clock tables below were re-recorded together in one sitting (min
over 28 reps), so their rows are comparable to each other but not to older
copies of this table -- std's own numbers moved too, glibc 2.43's `asinhf`
in particular being ~3x faster than whatever was last recorded here.

The re-record was forced by a benchmark fix. quickbench's dependency chain
and its throughput input array were both hardcoded to `|x|` in `[2,4)`,
which is *outside the domain* of `asin`/`acos`/`atanh`. libm takes its
domain check and returns NaN in a couple of ns, while this crate's
branchless versions compute the full result either way -- so those three
published "0.2x vs std" on latency while measuring an early `ret`. In
domain they are 1.1-1.4x *ahead* of std, and their throughput win had been
undersold by roughly 2x. Every benched function now names its own input
`Band` (examples/support/mca_common.rs), and quickbench walks the chain
untimed first and complains on stderr if the band ever leaves the domain --
which is how `log1pmx` turned up as a fourth case, out of domain via its
own negative output rather than its magnitude. The llvm-mca table further
down needed no change: its numbers came back bit-identical, since llvm-mca
reads the instruction stream and never sees a value.
```
Serial latency (dependency chain, examples/quickbench.rs; lower is better)
              | jodie   | std     | improvement
--------------|---------|---------|------------
         cbrt | 15.8 ns | 23.7 ns | 1.5x
cbrt_unchecked| 14.1 ns | 23.7 ns | 1.7x
cbrt_accurate | 22.6 ns | 23.7 ns | 1.0x
cbrt_accurate_unchecked| 21.6 ns | 23.7 ns | 1.1x
        rcbrt | 22.4 ns |    -    |  -
          cos | 18.3 ns | 16.4 ns | 0.9x
  cos_checked | 37.9 ns | 16.4 ns | 0.4x
         exp2 | 11.1 ns | 12.8 ns | 1.1x
 exp2_checked | 14.8 ns | 12.8 ns | 0.9x
        exp10 | 16.4 ns |    -    |  -
exp10_checked | 17.6 ns |    -    |  -
         log2 | 12.3 ns | 13.9 ns | 1.1x
log2_unchecked | 11.8 ns |    -    |  -
          sin | 16.3 ns | 17.1 ns | 1.0x
  sin_checked | 37.1 ns | 17.1 ns | 0.5x
           ln | 12.9 ns | 12.6 ns | 1.0x
  ln_unchecked | 11.0 ns |    -    |  -
        log10 | 12.2 ns | 14.5 ns | 1.2x
log10_unchecked| 10.8 ns |    -    |  -
        log1p | 15.4 ns | 17.1 ns | 1.1x
       log2p1 | 13.4 ns |    -    |  -
          exp | 11.2 ns | 10.5 ns | 0.9x
  exp_checked | 12.8 ns | 10.5 ns | 0.8x
        expm1 | 11.5 ns | 15.0 ns | 1.3x
exp_m1_over_x | 15.3 ns |    -    |  -
       exp2m1 | 12.4 ns |    -    |  -
         sinh | 13.5 ns | 15.7 ns | 1.2x
         cosh | 13.3 ns | 15.8 ns | 1.2x
 sinh_checked | 16.9 ns | 15.7 ns | 0.9x
 cosh_checked | 15.9 ns | 15.8 ns | 1.0x
         tanh | 15.9 ns | 17.3 ns | 1.1x
      sigmoid | 15.8 ns |    -    |  -
     softplus | 18.3 ns |    -    |  -
   logaddexp | 18.8 ns |    -    |  -
        asinh | 23.8 ns | 22.8 ns | 1.0x
        acosh | 20.2 ns | 22.1 ns | 1.1x
        atanh | 20.3 ns | 22.2 ns | 1.1x
         asin | 12.8 ns | 15.5 ns | 1.2x
         acos | 12.1 ns | 16.6 ns | 1.4x
         atan | 15.9 ns | 21.6 ns | 1.4x
 atan_latency | 16.5 ns |    -    |  -
      atan2 (*) | 20.6 ns | 26.7 ns | 1.3x
atan2_unchecked | 20.4 ns |    -    |  -
          tan | 19.3 ns | 22.8 ns | 1.2x
        sinpi | 12.8 ns |    -    |  -
        cospi | 15.5 ns |    -    |  -
        tanpi | 15.2 ns |    -    |  -
         sind | 15.8 ns |    -    |  -
         cosd | 18.2 ns |    -    |  -
         tand | 22.9 ns |    -    |  -
          erf | 16.2 ns |    -    |  -
         erfc | 17.2 ns |    -    |  -
        erfcx | 10.9 ns |    -    |  -
      hypot (*) | 6.0 ns  | 11.2 ns | 1.9x
hypot_unchecked | 6.2 ns  |    -    |  -
  hypot_checked | 18.2 ns | 11.2 ns | 0.6x
       rhypot | 8.9 ns  |    -    |  -
        rsqrt | 7.5 ns  |    -    |  -
      powf (*)(r)| 21.4 ns | 16.9 ns | 0.8x
powf_unchecked (r)| 18.8 ns | 16.9 ns | 0.9x
 powf_checked (r)| 33.4 ns | 16.9 ns | 0.5x
powf_checked_unchecked (r)| 30.3 ns | 16.9 ns | 0.6x
 remainder (*)| 10.7 ns |    -    |  -
remainder_unchecked| 9.3 ns  |    -    |  -
remainder_checked| 13.0 ns |    -    |  -
   remainder_ieee| 9.2 ns  |    -    |  -
   remainder_wide| 46.1 ns |    -    |  -
         fmod| 9.1 ns  |    -    |  -
fmod_unchecked| 7.3 ns  |    -    |  -
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

(r) the powf block in both tables (including its own `std powf` reference
column) was re-recorded in a *later* sitting than the rest, for
`powf_checked`'s atanh-form log2 -- so those rows are comparable to each
other but read ~11-15% fast against the rows around them. The unchanged
`powf` row is the conversion factor: same code, 24.0 -> 21.4 ns latency
and 1.37 -> 1.13 ns throughput, with `std powf` moving the same way
(19.3 -> 16.9, 6.45 -> 5.42).

Do not read `powf_checked`'s own cost off this table by subtracting across
sittings. Measured properly -- both binaries built once, then run
alternately, min of 6 each -- the atanh-form log2 costs it **-2.6%
latency** (34.3 -> 33.4 ns; it is *faster*) and **+6.8% throughput**
(1.76 -> 1.88 ns), and `powf_checked_unchecked` -2.6% / +8.5%. That
interleaving is what makes the numbers trustworthy on a machine that
throttles this much: the three rows whose code did *not* change (`powf`,
`powf_unchecked`, `std powf`) came back within 0.4% across the same runs,
which is the noise floor the deltas above are measured against. llvm-mca,
which needs no sitting caveat at all, agrees on both signs.

```
Throughput (independent array evals over [f32; 4096], examples/quickbench.rs; lower is better)
              | jodie    | std     | improvement
--------------|----------|---------|------------
         cbrt | 0.51 ns | 4.46 ns | 8.7x
cbrt_unchecked| 0.32 ns | 4.46 ns | 14.1x
cbrt_accurate | 0.91 ns | 4.46 ns | 4.9x
cbrt_accurate_unchecked| 0.65 ns | 4.46 ns | 6.8x
        rcbrt | 0.55 ns |    -    |  -
          cos | 0.39 ns | 3.65 ns | 9.5x
  cos_checked | 1.69 ns | 3.65 ns | 2.2x
         exp2 | 0.26 ns | 2.95 ns | 11.4x
 exp2_checked | 0.46 ns | 2.95 ns | 6.4x
        exp10 | 0.41 ns |    -    |  -
exp10_checked | 0.53 ns |    -    |  -
         log2 | 0.49 ns | 3.87 ns | 7.8x
log2_unchecked | 0.30 ns |    -    |  -
          sin | 0.33 ns | 4.15 ns | 12.8x
  sin_checked | 1.63 ns | 4.15 ns | 2.5x
           ln | 0.53 ns | 3.41 ns | 6.5x
  ln_unchecked | 0.31 ns |    -    |  -
        log10 | 0.45 ns | 4.88 ns | 10.9x
log10_unchecked| 0.31 ns |    -    |  -
        log1p | 0.54 ns | 5.19 ns | 9.6x
       log2p1 | 0.49 ns |    -    |  -
          exp | 0.34 ns | 2.44 ns | 7.3x
  exp_checked | 0.40 ns | 2.44 ns | 6.2x
        expm1 | 0.46 ns | 5.51 ns | 12.0x
exp_m1_over_x | 0.52 ns |    -    |  -
       exp2m1 | 0.51 ns |    -    |  -
         sinh | 0.49 ns | 6.21 ns | 12.6x
         cosh | 0.42 ns | 6.18 ns | 14.7x
 sinh_checked | 0.63 ns | 6.21 ns | 9.9x
 cosh_checked | 0.53 ns | 6.18 ns | 11.7x
         tanh | 0.47 ns | 5.57 ns | 11.7x
      sigmoid | 0.32 ns |    -    |  -
     softplus | 0.58 ns |    -    |  -
   logaddexp | 0.60 ns |    -    |  -
        asinh | 1.07 ns | 7.37 ns | 6.9x
        acosh | 0.98 ns | 6.25 ns | 6.4x
        atanh | 0.81 ns | 7.27 ns | 9.0x
         asin | 0.26 ns | 6.22 ns | 23.7x
         acos | 0.21 ns | 6.43 ns | 30.8x
         atan | 0.35 ns | 6.83 ns | 19.3x
 atan_latency | 0.40 ns |    -    |  -
  atan2 (*) | 0.55 ns | 11.04 ns | 20.1x
atan2_unchecked | 0.52 ns |    -    |  -
          tan | 0.60 ns | 7.21 ns | 11.9x
        sinpi | 0.35 ns |    -    |  -
        cospi | 0.40 ns |    -    |  -
        tanpi | 0.65 ns |    -    |  -
         sind | 0.32 ns |    -    |  -
         cosd | 0.44 ns |    -    |  -
         tand | 0.79 ns |    -    |  -
          erf | 0.56 ns |    -    |  -
         erfc | 0.69 ns |    -    |  -
        erfcx | 0.65 ns |    -    |  -
  hypot (*) | 0.20 ns | 2.55 ns | 12.9x
hypot_unchecked | 0.21 ns |    -    |  -
  hypot_checked | 0.38 ns | 2.55 ns | 6.7x
       rhypot | 0.36 ns |    -    |  -
        rsqrt | 0.35 ns |    -    |  -
      powf (*)(r)| 1.13 ns | 5.42 ns | 4.8x
powf_unchecked (r)| 0.70 ns | 5.42 ns | 7.7x
 powf_checked (r)| 1.88 ns | 5.42 ns | 2.9x
powf_checked_unchecked (r)| 1.42 ns | 5.42 ns | 3.8x
 remainder (*)| 0.23 ns |    -    |  -
remainder_unchecked| 0.17 ns |    -    |  -
remainder_checked| 0.34 ns |    -    |  -
   remainder_ieee| 0.20 ns |    -    |  -
   remainder_wide| 1.82 ns |    -    |  -
         fmod| 0.20 ns |    -    |  -
fmod_unchecked| 0.16 ns |    -    |  -
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
cbrt_accurate       |          63.00 |             3.129
cbrt_accurate_unchecked |      63.00 |             2.067
rcbrt               |          60.03 |             1.666
exp2                |          35.00 |             0.854
exp2_checked        |          47.00 |             1.399
exp10               |          52.00 |             1.461
exp10_checked       |          55.00 |             1.897
log2                |          34.23 |             1.584
log2_unchecked      |          34.23 |             0.958
sin                 |          48.00 |             1.151
sin_checked         |         117.02 |             5.232
cos                 |          56.00 |             1.406
cos_checked         |         122.00 |             4.603
sinpi               |          42.02 |             1.133
cospi               |          51.00 |             1.283
tanpi               |          78.88 |             2.031
sinc                |          53.02 |             1.256
sind                |          48.00 |             1.151
cosd                |          56.00 |             1.406
tand                |          70.02 |             2.533
ln                  |          44.14 |             1.611
ln_unchecked        |          34.06 |             1.018
log10               |          48.14 |             1.635
log10_unchecked     |          38.22 |             1.113
log1p               |          47.24 |             1.857
log2p1              |          48.14 |             1.886
exp                 |          42.00 |             1.230
exp_checked         |          50.00 |             1.607
expm1               |          70.00 |             1.620
expm1_checked       |          78.00 |             1.556
exp_m1_over_x       |          82.00 |             1.720
exp2m1              |          80.00 |             1.843
sinh                |          51.00 |             1.780
cosh                |          50.00 |             1.647
sinh_throughput     |          62.00 |             1.943
cosh_throughput     |          61.00 |             1.616
sinh_checked        |          59.00 |             2.274
cosh_checked        |          58.00 |             1.943
tanh                |          85.64 |             1.731
sigmoid             |          61.00 |             1.222
softplus            |          74.11 |             2.449
logaddexp           |          74.11 |             2.449
asinh               |          69.41 |             4.136
acosh               |          89.08 |             3.828
atanh               |          96.83 |             2.903
asin                |          60.99 |             0.968
acos                |          39.99 |             0.820
atan                |          61.27 |             1.491
atan_latency        |          61.99 |             1.591
atan2               |          67.19 |             1.694
tan                 |          71.02 |             2.532
erf                 |          83.98 |             2.037
erfc                |          64.00 |             2.437
erfcx               |          62.02 |             2.278
hypot               |          21.11 |             0.766
hypot_checked       |          57.19 |             1.178
rhypot              |          32.02 |             1.389
rsqrt               |          28.00 |             1.381
powf                |         104.95 |             5.651
powf_unchecked      |          79.05 |             3.095
powf_checked        |         152.89 |             9.918
powf_checked_unchecked |      128.77 |             7.255
remainder           |          34.11 |                 ? (*)
remainder_unchecked |          33.00 |             0.646
remainder_checked   |          45.17 |             1.357
remainder_ieee      |          29.11 |                 ? (*)
remainder_wide      |         171.17 |             7.688
fmod                |          29.11 |                 ? (*)
fmod_unchecked      |          28.00 |             0.643
```

(*) remainder/remainder_ieee/fmod: throughput no longer measurable via
llvm-mca after their 2026-07-09 zero/nan fix (backlog idea #85's own
follow-up) -- LLVM branch-specializes these short functions' vectorized
loop on the harness's own shared black_box'd `y`, producing multiple
physical exit paths that corrupt llvm-mca's region parser (the same
class of harness limitation `ldexp`/`frexp`/`rootn` document too;
confirmed via a standalone `--emit=asm` probe with genuinely per-lane-varying inputs
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
- `cargo run --release --example powfsearch` - adversarial worst case for `powf`/`powf_checked`, where a
  blind 2-arg fuzz is structurally weak: their error is `ln2 * |y*log2(x)| * relerr(log2)`, maximised only
  where both factors are extreme at once, so this derives `y` from `x` to pin the first and sweeps the
  octave that maximises the second exhaustively instead of sampling. Finds twice the max ulp the fuzz does.
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
