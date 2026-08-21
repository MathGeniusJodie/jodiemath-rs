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
                  rcbrt |    0.117   |     1     | (no std rcbrt)
                   exp2 |    0.026   |     1     |  0.000  |    1
           exp2_checked |    0.014   |     1     |  0.000  |    1
                  exp10 |    0.031   |     2     | (no std exp10)
          exp10_checked |    0.008   |     1     | (no std exp10)
                   log2 |    0.003   |     1     |  0.000  |    1
   log2_unchecked (+)   |    0.006   |     1     | (bit-identical to log2 on its domain)
      sin (|x|<2^24*pi) |    0.046   |     2     |  0.003  |    1
      cos (|x|<2^22*pi) |    0.083   |     2     |  0.002  |    1
    sin_fast (|x|<=1e6) |    0.036   |     2     |  0.002  |    1
    cos_fast (|x|<=1e6) |    0.078   |     3     |  0.002  |    1
  sin_fast (|x|<2^22*pi)|    0.059   |   219     |  0.003  |    1
  cos_fast (|x|<2^22*pi)|    0.288   |  2780     |  0.002  |    1
 sin_checked (|x|<=1e6) |    0.036   |     2     |  0.000  |    1
 cos_checked (|x|<=1e6) |    0.076   |     2     |  0.000  |    1
      sin_wide (all f32)|    0.127   |     2     |  0.000  |    1
      cos_wide (all f32)|    0.150   |     2     |  0.000  |    1
        sinpi (all f32) |    0.197   |     2     | (no std sinpi)
        cospi (all f32) |    0.058   |     2     | (no std cospi)
        tanpi (all f32) |    0.033   |     2     | (no std tanpi)
        sind (|x|<4.7e7)|    0.124   |     2     | (no std sind)
        cosd (|x|<4.7e7)|    0.073   |     2     | (no std cosd)
        tand (|x|<4.7e7)|    0.038   |     2     | (no std tand)
          sinc (|x|<1e6)|    0.094   |     4     | (no std sinc)
```

`sin_wide`/`cos_wide`/`tan_wide`'s rows are the only trig rows in this
table measured over the **entire** f32 range, and that is the whole point
of the tier.
`sin_checked`/`cos_checked` are accurate to ~1e13 and then degrade, and
past `2^51*pi` they return a value in `[-1,1]` with no relationship to the
answer -- `accuracy`'s own `sin_checked (all f32)` row reads **avg
314265980 / max 2130706432** on a quick fuzz (that max is `2*0x3f800000`,
the ulp distance from `+1` to `-1`, i.e. the worst a clamped output can
be). `tan_checked` is worse still, avg 406004054 / max 2324484283, and has no
row in this table at all. The `_wide` tier reduces against a window of
`1/pi` selected by `x`'s exponent instead of a fixed two-word constant,
which costs four 32-bit gathers and ~2x throughput; see `sin_wide`'s doc
comment and graveyard.md. All three `_wide` rows are exhaustive over all
2^32 patterns, not sampled.

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
                compound |    0.196   |   ~200    | 
       compound_accurate |    0.001   |     1     | 
         exp (in-domain) |    0.071   |     3     |  0.000  |    1
             exp_checked |    0.004   |     1     |  0.000  |    1
       expm1 (in-domain) |    0.008   |     2     |  0.000  |    0
           expm1_checked |    0.004   |     2     |  0.000  |    0
exp_m1_over_x (in-domain)|    0.017   |     2     | (no std exp_m1_over_x)
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
        softplus (|x|<80)|    0.075   |     3     | (no std softplus)
        softplus_checked |    0.039   |     3     | (no std softplus)
      logsigmoid_checked |    0.039   |     3     | (no std logsigmoid)
   silu_checked (all f32)|    0.045   |     5     | (no std silu; 4.69 on a dense scan -- `silu` is 3.85, the avg/max trade in its doc)
    gelu (|x|<=10*sqrt2) |    0.131   |     6     | (no std gelu; exhaustive)
  gelu (all finite f32)  |    0.067   |     6     |
   logaddexp (|a|,|b|<80)|    0.141   |  ~1e3-1e5, heavy-tailed (real, narrow cancellation -- see its own doc comment) | (no std logaddexp)
logaddexp_checked (all f32)|  0.037   |  ~1e3-1e5, heavy-tailed (same cancellation, bit-identical over |a-b|<=87; the lower avg is the wider domain, not a better answer) | (no std logaddexp)
logaddexp_accurate (all f32)| 0.000   |     0     | (no std logaddexp; the cancellation fixed -- f64 correction, ~3e-16 absolute. 10M-pair fuzz is 0/0 against the f64 reference, and on a corpus built *on* the zero curve `e^a+e^b=1` it is avg 0.27 / max 72 against an 80-digit oracle where the other two tiers are avg 1.7e7 / max 3.4e10)
                  asinh |    0.034   |     2     | 
                  acosh |    0.003   |     3     |  0.000  |    1
                  atanh |    0.004   |     2     | 
                   asin |    0.016   |     2     |  0.000  |    0
                   acos |    0.002   |     2     |  0.000  |    0
                   atan |    0.063   |     3     |  0.000  |    0
           atan_latency |    0.052   |     3     |  0.000  |    0
                 atanpi |    0.031   |     4     | (no std atanpi; `1/pi` folds in *before* the quadrant reflection, whose constant is then an exact `0.5`. The max is `atan_poly`'s own, on the `|x|<1` arm this leaves bit-identical)
      tan (|x|<2^22*pi) |    0.118   |     4     |  0.000  |    0
      tan_wide (all f32)|    0.250   |     4     |  0.000  |    0
                   erf  |    0.027   |     3     | (no std erf)
         erfc (|x|<=10) |    0.122   |     6     |
        erfcx (|x|<=20) |    0.134   |     4     | (no std erfcx; the 4 is on the *negative* arm, where the exponential's error arrives amplified by the reflection -- see erfcx's doc comment)
         erfcx (x>=20)  |    0.268   |     2     | (no std erfcx)
     norm_cdf (all f32) |    0.062   |     6     | (no std norm_cdf; exhaustive. Same two terms as `erfc`, and now in the same proportion -- mostly the erfcx polynomial)
     norm_pdf (all f32) |    0.018   |     2     | (no std norm_pdf; exhaustive. The most exponential-bound function in the family: the rest of the chain is an exactly-split square and a two-word `1/sqrt(2*pi)`)
       dawson (all f32) |    0.048   |     6     | (no std dawson; exhaustive. The 6 is the central branch at x=1.421, where the rational's own fit is already 3; the `|x|>4` tail is peeled and sits exactly on its own fit-only floor, max 3)
       erfinv (|x|<1)   |    0.047   |     3     | (no std erfinv; exhaustive. Almost all of the average is one number: the effective leading coefficient `1+c0` of the central poly, which is what every `|x|` under ~0.06 computes)
      erfc_inv (0<y<2)  |    0.380   |     4     | (no std erfc_inv; exhaustive. 81% of the domain's bit patterns take the far-tail branch, so its poly sets this)
        probit (0<p<1)  |    0.494   |     4     | (no std probit; exhaustive)
                  atan2 |    0.066   |     4     |  0.000  |    0
    atan2_unchecked (+) |    0.066   |     4     | (bit-identical to atan2 on its domain)
              atan2_pos |    0.062   |     3     | (no std atan2_pos)
                 atan2d |    0.105   |     4     | (no std atan2d)
                atan2pi |    0.070   |     3     | (no std atan2pi)
        hypot (bounded) |    0.034   |     1     |  0.000  |    0
          hypot_checked |    0.015   |     1     | (no std comparison needed, no domain restriction)
                 rhypot |    0.065   |     2     | (no std rhypot)
                  rsqrt |    0.260   |     1     | (no std rsqrt)
        powf (in-domain)|    0.001   |   1 (s)  |  0.000  |    1
    powf_unchecked (+)  |    0.001   |   1 (s)  | (bit-identical to powf on its domain)
                powf_pos|    0.001   |   1 (s)  | (bit-identical to powf for x > 0, see its doc comment)
remainder (|x/y|<1000, near-tie excluded) | 0.000 | 0 | (no std comparison needed)
remainder_unchecked (+) |    0.000   |     0     | (bit-identical to remainder on its domain)
    remainder_checked (|x/y|<1e7, near-tie excluded) | 0.000 | 0 | (no std comparison needed)
       remainder_wide (|x/y|<4e15, near-tie excluded) | 0.000 | 0 | (no std comparison needed)
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
         log2 | 14.4 ns | 14.2 ns | 1.0x
log2_unchecked| 13.4 ns |    -    |  -
          sin | 16.3 ns | 17.1 ns | 1.0x
  sin_checked | 37.1 ns | 17.1 ns | 0.5x
           ln | 12.9 ns | 12.6 ns | 1.0x
 ln_unchecked | 11.0 ns |    -    |  -
        log10 | 12.2 ns | 14.5 ns | 1.2x
log10_unchecked|10.8 ns |    -    |  -
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
        atan2 | 20.6 ns | 26.7 ns | 1.3x
atan2_unchecked |20.4 ns|    -    |  -
          tan | 19.3 ns | 22.8 ns | 1.2x
        sinpi | 12.8 ns |    -    |  -
        cospi | 15.5 ns |    -    |  -
        tanpi | 15.2 ns |    -    |  -
         sind | 15.8 ns |    -    |  -
         cosd | 18.2 ns |    -    |  -
         tand | 14.8 ns |    -    |  -
          erf | 16.2 ns |    -    |  -
      erfc (!)| 20.2 ns |    -    |  -
     erfcx (!)| 14.0 ns |    -    |  -
        hypot |  6.0 ns | 11.2 ns | 1.9x
hypot_checked | 18.2 ns | 11.2 ns | 0.6x
       rhypot |  8.9 ns |    -    |  -
        rsqrt |  7.5 ns |    -    |  -
         powf | 42.7 ns | 23.9 ns | 0.6x
powf_unchecked |39.2 ns | 23.9 ns | 0.6x
    remainder | 10.7 ns |    -    |  -
remainder_unchecked| 9.3 ns |    -    |  -
remainder_checked | 13.0 ns |    -    |  -
   remainder_ieee |  9.2 ns |    -    |  -
   remainder_wide | 46.1 ns |    -    |  -
             fmod |  9.1 ns |    -    |  -
   fmod_unchecked |  7.3 ns |    -    |  -
```

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
      log2 (r)| 0.51 ns | 3.86 ns | 7.6x
log2_unchecked| 0.34 ns |    -    |  -
          sin | 0.33 ns | 4.15 ns | 12.8x
  sin_checked | 1.63 ns | 4.15 ns | 2.5x
           ln | 0.53 ns | 3.41 ns | 6.5x
 ln_unchecked | 0.31 ns |    -    |  -
        log10 | 0.45 ns | 4.88 ns | 10.9x
log10_unchecked|0.31 ns |    -    |  -
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
        atan2 | 0.55 ns | 11.04 ns | 20.1x
atan2_unchecked|0.52 ns |    -    |  -
          tan | 0.60 ns | 7.21 ns | 11.9x
        sinpi | 0.35 ns |    -    |  -
        cospi | 0.40 ns |    -    |  -
        tanpi | 0.65 ns |    -    |  -
         sind | 0.32 ns |    -    |  -
         cosd | 0.44 ns |    -    |  -
         tand | 0.74 ns |    -    |  -
          erf | 0.56 ns |    -    |  -
         erfc | 0.94 ns |    -    |  -
        erfcx | 0.88 ns |    -    |  -
    hypot (*) | 0.20 ns | 2.55 ns | 12.9x
hypot_checked | 0.38 ns | 2.55 ns | 6.7x
       rhypot | 0.36 ns |    -    |  -
        rsqrt | 0.35 ns |    -    |  -
         powf | 2.27 ns | 7.74 ns | 3.4x
powf_unchecked| 1.89 ns | 7.74 ns | 4.1x
    remainder | 0.23 ns |    -    |  -
remainder_unchecked| 0.17 ns |    -    |  -
remainder_checked| 0.34 ns |    -    |  -
   remainder_ieee| 0.20 ns |    -    |  -
   remainder_wide| 1.82 ns |    -    |  -
          fmod | 0.20 ns |    -    |  -
fmod_unchecked | 0.16 ns |    -    |  -
```

```
theoretical cost from llvm-mca (-mcpu=native, 100 iterations)
(*L) region kept its loop; divided by the real step, not ARR_LEN -- see above
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
sin_fast            |          48.00 |             1.151
sin_checked         |          82.00 |             2.495
sin_wide            |          94.14 |             4.920
cos                 |          61.00 |             1.654
cos_fast            |          56.00 |             1.406
cos_checked         |          87.00 |             3.157
cos_wide            |          95.16 |             5.434
tan_checked         |         101.00 |             4.289
tan_wide            |         109.17 |             6.982 (*L)
sinpi               |          42.02 |             1.133
cospi               |          47.00 |             1.226
tanpi               |          77.05 |             2.155
sinc                |          53.02 |             1.256
sind                |          48.00 |             1.151
cosd                |          56.00 |             1.406
tand                |      (no row)  |             2.290
ln                  |          49.13 |             1.625
ln_unchecked        |          38.06 |             1.021
log10               |          48.16 |             1.617
log10_unchecked     |          38.06 |             1.022
log1p               |          51.25 |             1.903
log2p1              |          51.36 |             1.876
compound            |          99.69 |             3.829
compound_accurate   |         120.77 |             6.173
exp                 |          42.00 |             1.195
exp_checked         |          54.00 |             1.526
expm1               |          47.00 |             1.087
expm1_checked       |          55.00 |             1.320
exp_m1_over_x       |          65.03 |             1.279
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
softplus            |          74.11 |             2.449
softplus_checked    |          73.36 |             2.506
logsigmoid          |          75.11 |             2.604
logsigmoid_checked  |          75.56 |             2.699
silu                |          65.02 |             1.407
silu_checked        |          60.97 |             1.466
gelu                |          71.88 |             3.434
logaddexp           |          74.11 |             2.449
logaddexp_checked   |          73.35 |             2.506
logaddexp_accurate  |         121.52 |             7.616
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
erfcx               |          67.99 |             3.006
hypot               |          21.11 |             0.766
hypot_checked       |          57.19 |             1.178
rhypot              |          32.02 |             1.389
rsqrt               |          28.00 |             1.381
powf                |         125.81 |             6.996
powf_unchecked      |         118.00 |             6.171
remainder           |          40.13 |                 ? (*)
remainder_unchecked |          33.00 |             0.646
remainder_checked   |          45.17 |             1.357
remainder_ieee      |          35.13 |                 ? (*)
remainder_wide      |          58.13 |             2.266
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

- `cargo run --release --example worst_corpus` - 105 public 1-arg functions x 90 historically-hard inputs
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
  `hypot(re,im).ln()`, whose relative error is `2^-53/|re^2+im^2-1|`, i.e. already ~2 f32 ulp at `|v| ~ 1e-9`
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
  hand-maintained, and both were incomplete: `softplus`/`logsigmoid`/`silu` flush their *entire* denormal
  range (16.97, 16.97 and 20.28 premature in `x`, all worse than the `sigmoid` case that prompted this
  example) and none of them was listed. Adding a function here is not optional bookkeeping -- nothing else
  in the suite reports this. All three now have a `_checked` sibling that carries the band ("never" flushes,
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
