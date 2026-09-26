# faster? f32 math functions
Attempting to provide faster implementations of common f32 math functions with a similar level of accuracy to the standard library.
There are also perfectly rounded variants and faster-but-sloppier variants.

All functions auto-vectorize, it's a hard requirement.

# precision

Every number below is generated: `tools/readme_stats.sh` re-runs all of it and rewrites the tables in
place (`--tables` just re-renders the last run). ulp distance against an f64 reference rounded to f32
(`examples/accuracy.rs`). 1-arg rows are exhaustive over their domain (`accuracy thorough`: every f32 bit
pattern the row's filter admits, `2^32` = all of them); 2-arg rows are fuzzed, see `inputs`. A row's
domain is its `measure!` filter in `examples/accuracy.rs`. `(+)` = measured on the function's documented
precondition domain. `std` is Rust's std over the same domain, where it has the function.

<!-- BEGIN generated:accuracy -->
| function                      | avg ulp | max ulp | std avg | std max |   inputs |                                                    worst x |
|:------------------------------|--------:|--------:|--------:|--------:|---------:|-----------------------------------------------------------:|
| `cbrt`                        |  0.2813 |       3 |  0.0000 |       0 |     2^32 |                                                2.04189e-40 |
| `cbrt_accurate`               |  0.0000 |       1 |  0.0000 |       0 |     2^32 |                                              2.4120957e-38 |
| `cbrt_fast (+)`               | 57.3676 |     554 |  0.0000 |       0 | 2.13e+09 |                                               4.2535296e37 |
| `cbrt_unchecked (+)`          |  0.2824 |       3 |  0.0000 |       0 | 4.26e+09 |                                               1.303716e-38 |
| `cbrt_accurate_unchecked (+)` |  0.0000 |       1 |  0.0000 |       0 | 3.07e+09 |                                                1.42385e-17 |
| `rcbrt`                       |  0.1174 |       1 |       - |       - |     2^32 |                                                      7e-45 |
| `pow_3_2`                     |  0.1692 |       1 |       - |       - | 2.14e+09 |                                              1.6408968e-30 |
| `pow_2_3`                     |  0.1032 |       1 |       - |       - |     2^32 |                                                      4e-45 |
| `log_2`                       |  0.0028 |       1 |  0.0001 |       1 |     2^32 |                                                   1.12e-42 |
| `log_2_unchecked (+)`         |  0.0057 |       1 |  0.0001 |       1 | 2.13e+09 |                                              1.1807558e-38 |
| `exp2`                        |  0.0264 |       1 |  0.0001 |       1 | 2.25e+09 |                                               8.5991324e-8 |
| `exp2_checked`                |  0.0138 |       1 |  0.0001 |       1 |     2^32 |                                               8.5991324e-8 |
| `exp10`                       |  0.0307 |       2 |       - |       - | 2.22e+09 |                                               -3.7629665e1 |
| `exp10_checked`               |  0.0082 |       1 |       - |       - | 2.22e+09 |                                               2.5885969e-8 |
| `sin (\|x\|<2^24*pi)`         |  0.0457 |       2 |  0.0028 |       1 | 2.56e+09 |                                                1.3414527e0 |
| `cos (\|x\|<2^22*pi)`         |  0.0833 |       2 |  0.0022 |       1 | 2.53e+09 |                                               1.7422438e-4 |
| `sin_wide (all f32)`          |  0.1269 |       2 |  0.0068 |       1 |     2^32 |                                                1.3414527e0 |
| `cos_wide (all f32)`          |  0.1270 |       2 |  0.0066 |       1 |     2^32 |                                                 2.835277e0 |
| `sinpi (all f32)`             |  0.1969 |       2 |       - |       - |     2^32 |                                                9.947858e-3 |
| `cospi (all f32)`             |  0.0578 |       2 |       - |       - |     2^32 |                                               5.5447224e-5 |
| `tanpi (all f32)`             |  0.0327 |       2 |       - |       - |     2^32 |                                               3.2719934e-1 |
| `sinc_unnormalized (\|x\|<1e6)` |  0.0716 |       2 |       - |       - | 2.46e+09 |                                                 6.25657e-2 |
| `ln`                          |  0.0036 |       1 |  0.0001 |       1 |     2^32 |                                                   9.14e-43 |
| `ln_unchecked (+)`            |  0.0073 |       1 |  0.0001 |       1 | 2.13e+09 |                                              1.1772184e-38 |
| `log10`                       |  0.0034 |       1 |  0.0000 |       1 |     2^32 |                                                   5.93e-43 |
| `log10_unchecked (+)`         |  0.0069 |       1 |  0.0000 |       1 | 2.13e+09 |                                              1.1778674e-38 |
| `log1p`                       |  0.0249 |       2 |  0.0000 |       1 |     2^32 |                                               1.3032487e-1 |
| `log1pmx`                     |  0.0632 |       3 |       - |       - |     2^32 |                                                3.722855e-1 |
| `log2p1`                      |  0.0324 |       2 |       - |       - |     2^32 |                                                7.450616e-8 |
| `log10p1`                     |  0.0387 |       2 |       - |       - |     2^32 |                                                6.306496e-2 |
| `exp`                         |  0.0707 |       3 |  0.0001 |       1 | 2.24e+09 |                                                1.5595897e1 |
| `exp_narrow`                  |  0.0707 |       3 |  0.0001 |       1 | 2.24e+09 |                                                1.5595897e1 |
| `expm1`                       |  0.0081 |       2 |  0.0000 |       0 | 2.24e+09 |                                                3.621128e-1 |
| `expm1_narrow`                |  0.0081 |       2 |  0.0000 |       0 | 2.24e+09 |                                                3.621128e-1 |
| `expm1_checked`               |  0.0042 |       2 |  0.0000 |       0 |     2^32 |                                                3.621128e-1 |
| `exp_checked`                 |  0.0042 |       1 |  0.0001 |       1 |     2^32 |                                               5.9604645e-8 |
| `exp_m1_over_x_narrow`        |  0.0170 |       2 |       - |       - | 2.24e+09 |                                               3.5792866e-1 |
| `exp2m1`                      |  0.0400 |       2 |       - |       - |     2^32 |                                                5.000016e-1 |
| `exp10m1`                     |  0.0946 |       3 |       - |       - |     2^32 |                                               1.5052298e-1 |
| `sinh`                        |  0.0608 |       3 |  0.0000 |       1 | 2.24e+09 |                                               2.0394178e-1 |
| `cosh`                        |  0.0244 |       2 |  0.0000 |       0 | 2.24e+09 |                                               3.4671992e-1 |
| `sinh_throughput`             |  0.0795 |       4 |  0.0000 |       1 | 2.24e+09 |                                                7.560597e-1 |
| `cosh_throughput`             |  0.0487 |       3 |  0.0000 |       0 | 2.24e+09 |                                               9.6372455e-1 |
| `sinh_narrow`                 |  0.0608 |       3 |  0.0000 |       1 | 2.24e+09 |                                               2.0394178e-1 |
| `cosh_narrow`                 |  0.0244 |       2 |  0.0000 |       0 | 2.24e+09 |                                               3.4671992e-1 |
| `sinh_checked`                |  0.0317 |       3 |  0.0000 |       1 |     2^32 |                                               2.0394178e-1 |
| `cosh_checked`                |  0.0128 |       2 |  0.0000 |       0 |     2^32 |                                               3.4671992e-1 |
| `coshm1`                      |  0.0220 |       3 |       - |       - |     2^32 |                                                 2.085061e0 |
| `tanh`                        |  0.0440 |       2 |  0.0000 |       0 | 2.22e+09 |                                               6.2427483e-2 |
| `tanh_grad`                   |  0.1318 |       4 |       - |       - | 2.22e+09 |                                               8.9032984e-1 |
| `sigmoid`                     |  0.0909 |       3 |       - |       - | 2.24e+09 |                                              -5.0497115e-2 |
| `sigmoid_grad`                |  0.1308 |       4 |       - |       - | 2.24e+09 |                                                1.7806597e0 |
| `logsigmoid (\|x\|<80)`       |  0.0752 |       3 |       - |       - | 2.24e+09 |                                               4.3275338e-1 |
| `logsigmoid_checked`          |  0.0394 |       3 |       - |       - | 4.28e+09 |                                               4.3275338e-1 |
| `gelu (\|x\|<=10*sqrt2)`      |  0.1309 |       6 |       - |       - | 2.19e+09 |                                                -1.321733e0 |
| `gelu (all finite f32)`       |  0.0671 |       6 |       - |       - | 4.28e+09 |                                                -1.321733e0 |
| `silu`                        |  0.1082 |       4 |       - |       - | 2.24e+09 |                                              -2.6796955e-1 |
| `silu_checked`                |  0.0453 |       5 |       - |       - | 4.28e+09 |                                               -3.2044985e0 |
| `softsign`                    |  0.0412 |       1 |       - |       - | 4.28e+09 |                                                4.214685e-8 |
| `asinh`                       |  0.0340 |       2 |  0.0000 |       1 |     2^32 |                                               3.8987396e-3 |
| `acosh`                       |  0.0031 |       3 |  0.0000 |       1 |     2^32 |                                                1.0018609e0 |
| `atanh`                       |  0.0037 |       2 |  0.0369 |  363409 |     2^32 |                                               2.5002798e-1 |
| `asin`                        |  0.0158 |       2 |  0.0000 |       0 |     2^32 |                                                6.062731e-2 |
| `asind`                       |  0.0195 |       4 |       - |       - |     2^32 |                                               5.0003934e-1 |
| `asinpi`                      |  0.0057 |       2 |       - |       - |     2^32 |                                                5.000009e-1 |
| `acos`                        |  0.0022 |       2 |  0.0000 |       1 |     2^32 |                                               5.4030454e-1 |
| `acosd`                       |  0.0641 |       4 |       - |       - |     2^32 |                                                8.481021e-1 |
| `acospi`                      |  0.0438 |       3 |       - |       - |     2^32 |                                                8.796542e-3 |
| `atan`                        |  0.0627 |       3 |  0.0000 |       1 |     2^32 |                                                9.360248e-1 |
| `atan_latency`                |  0.0516 |       3 |  0.0000 |       1 |     2^32 |                                                1.0077639e0 |
| `atan_bounded`                |  0.0431 |       3 |  0.0000 |       1 | 2.13e+09 |                                                9.360248e-1 |
| `atand`                       |  0.0582 |       3 |       - |       - |     2^32 |                                               1.0878966e-3 |
| `tan (\|x\|<2^22*pi)`         |  0.0942 |       4 |  0.0000 |       0 | 2.53e+09 |                                                 6.110602e5 |
| `tan_wide`                    |  0.2307 |       3 |       - |       - |     2^32 |                                               9.9448454e-1 |
| `erf`                         |  0.0270 |       3 |       - |       - |     2^32 |                                               2.1674643e-1 |
| `erfc (\|x\|<=10)`            |  0.1217 |       6 |       - |       - | 2.19e+09 |                                                1.7139885e0 |
| `erfinv`                      |  0.0472 |       3 |       - |       - | 2.13e+09 |                                                 7.95631e-1 |
| `erfc_inv`                    |  0.3801 |       4 |       - |       - | 1.07e+09 |                                               1.5448519e-8 |
| `logit`                       |  0.0504 |       3 |       - |       - | 1.07e+09 |                                                3.456839e-1 |
| `srgb_to_linear`              |  0.0823 |       7 |       - |       - | 1.07e+09 |                                                8.461098e-2 |
| `linear_to_srgb`              |  0.1209 |       8 |       - |       - | 1.07e+09 |                                               3.4453226e-3 |
| `atan2`                       |  0.0661 |       3 |  0.0000 |       0 | 1.00e+07 |                                               3.373309e-19 |
| `atan2_latency`               |  0.0622 |       3 |  0.0000 |       0 | 1.00e+07 |                                                5.711285e18 |
| `atan2_unchecked (+)`         |  0.0660 |       3 |  0.0000 |       0 | 1.00e+07 |                                             -1.2636623e-33 |
| `atan2_pos`                   |  0.0621 |       3 |  0.0000 |       0 | 1.00e+07 |                                               2.0524904e-8 |
| `xlogy`                       |  0.1293 |       2 |       - |       - | 1.00e+07 |                                                5.118248e36 |
| `xlog1py`                     |  0.0975 |       3 |       - |       - | 1.00e+07 |                                             -4.2693671e-16 |
| `rsqrt`                       |  0.2599 |       1 |       - |       - | 2.14e+09 |                                                      4e-45 |
| `sqrt1pm1`                    |  0.2053 |       2 |       - |       - | 3.20e+09 |                                               4.1295314e-7 |
| `hypot3`                      |  0.0038 |       1 |       - |       - | 1.00e+07 |                x=1.3378555e14,y=2.4303084e-2,z=1.373468e12 |
| `rnorm3`                      |  0.0072 |       2 |       - |       - | 1.00e+07 |                 x=3.885292e6,y=6.9928085e6,z=-1.9862841e10 |
| `hypot4`                      |  0.0021 |       1 |       - |       - | 1.00e+07 | w=-5.155415e10,x=5.3304037e10,y=-9.05402e-5,z=-9.706332e13 |
| `rnorm4`                      |  0.0039 |       2 |       - |       - | 1.00e+07 |  w=-9.738718e12,x=6.147192e10,y=-3.520599e13,z=2.3639505e2 |
| `diff_of_products`            |  0.0240 |       1 |       - |       - | 6.49e+07 |                                             -1.3023296e-34 |
| `cross2`                      |  0.0241 |       1 |       - |       - | 6.49e+07 |                                                 4.53504e24 |
| `ldexp`                       |  0.0000 |       0 |       - |       - | 2.00e+07 |                                                  x=0e0,n=0 |
| `powf`                        |  0.0007 |       1 |  0.0000 |       1 | 4.99e+06 |                                              5.9030617e-32 |
| `powf_pos`                    |  0.0015 |       1 |  0.0000 |       1 | 2.50e+06 |                                               5.122882e-16 |
| `powf_unchecked (+)`          |  0.0015 |       1 |  0.0000 |       1 | 2.49e+06 |                                                 3.865901e5 |
| `fmod`                        |  0.0000 |       0 |       - |       - | 8.84e+05 |                                                        0e0 |
| `fmod_unchecked (+)`          |  0.0000 |       0 |       - |       - | 8.83e+05 |                                                        0e0 |
| `fmod_checked`                |  0.0000 |       0 |       - |       - | 5.34e+06 |                                                        0e0 |
<!-- END generated:accuracy -->

- Every `_unchecked` function is bit-identical to its checked sibling on its own domain, and `powf_pos`
  to `powf` for x > 0 (see their doc comments; `examples/unchecked_parity.rs`).
- `silu_checked` trades max for avg against `silu` -- see its doc comment.
- `erfinv`: almost all of the average is one number, the effective leading coefficient `1+c0` of the
  central poly, which is what every `|x|` under ~0.06 computes.
- `erfc_inv`: 81% of the domain's bit patterns take the far-tail branch, so its poly sets this row.
- `fmod`: `|x/y| < 1000` with near-integer quotients excluded; matches Rust's `%`.

Sub-domain bands and naive baselines (the narrow `sin`/`cos` rows past their domain are garbage by
contract; the `_wide` tier is what covers all of f32):

<!-- BEGIN generated:accuracy-bands -->
| row                      |         avg ulp |    max ulp |   inputs |
|:-------------------------|----------------:|-----------:|---------:|
| `sin \|x\|<=pi/4`        |          0.0010 |          1 | 2.12e+09 |
| `sin \|x\|<=10`          |          0.0077 |          2 | 2.19e+09 |
| `sin \|x\|<=1000`        |          0.0198 |          2 | 2.30e+09 |
| `sin \|x\|<=1e6`         |          0.0357 |          2 | 2.46e+09 |
| `sin [1e3,1e5)`          |          0.2527 |          2 | 1.10e+08 |
| `sin_wide [1e3,1e5)`     |          0.2523 |          2 | 1.10e+08 |
| `sin [1e5,1.3e7)`        |          0.2809 |          2 | 1.18e+08 |
| `sin_wide [1e5,1.3e7)`   |          0.2522 |          2 | 1.18e+08 |
| `sin [1.3e7,1e8)`        |  337896086.8740 | 2194039715 | 4.93e+07 |
| `sin_wide [1.3e7,1e8)`   |          0.2522 |          2 | 4.93e+07 |
| `sin [1e8,1e10)`         | 1277967960.5234 | 2671421137 | 1.12e+08 |
| `sin_wide [1e8,1e10)`    |          0.2522 |          2 | 1.12e+08 |
| `sin [1e10,1e13)`        | 1828371281.3905 | 3204448256 | 1.67e+08 |
| `sin_wide [1e10,1e13)`   |          0.2522 |          2 | 1.67e+08 |
| `sin [1e13,1e15)`        | 2132492643.6429 | 3204448256 | 1.11e+08 |
| `sin_wide [1e13,1e15)`   |          0.2522 |          2 | 1.11e+08 |
| `sin [1e15,3.4e38)`      | 2138952448.2860 | 3204448256 | 1.31e+09 |
| `sin_wide [1e15,3.4e38)` |          0.2522 |          2 | 1.31e+09 |
| `cos \|x\|<=pi/4`        |          0.0494 |          2 | 2.12e+09 |
| `cos \|x\|<=10`          |          0.0555 |          2 | 2.19e+09 |
| `cos \|x\|<=1000`        |          0.0650 |          2 | 2.30e+09 |
| `cos \|x\|<=1e6`         |          0.0779 |          2 | 2.46e+09 |
| `cos [1e3,1e5)`          |          0.2524 |          2 | 1.10e+08 |
| `cos_wide [1e3,1e5)`     |          0.2520 |          2 | 1.10e+08 |
| `cos [1e5,1.3e7)`        |          0.2809 |          2 | 1.18e+08 |
| `cos_wide [1e5,1.3e7)`   |          0.2522 |          2 | 1.18e+08 |
| `cos [1.3e7,1e8)`        | 1015144977.4156 | 2214573574 | 4.93e+07 |
| `cos_wide [1.3e7,1e8)`   |          0.2522 |          2 | 4.93e+07 |
| `cos [1e8,1e10)`         | 1291443513.8390 | 2711729072 | 1.12e+08 |
| `cos_wide [1e8,1e10)`    |          0.2523 |          2 | 1.12e+08 |
| `cos [1e10,1e13)`        | 1886395367.2903 | 3204448256 | 1.67e+08 |
| `cos_wide [1e10,1e13)`   |          0.2523 |          2 | 1.67e+08 |
| `cos [1e13,1e15)`        | 2137113767.9392 | 3204448256 | 1.11e+08 |
| `cos_wide [1e13,1e15)`   |          0.2522 |          2 | 1.11e+08 |
| `cos [1e15,3.4e38)`      | 2139094843.7518 | 3204448256 | 1.31e+09 |
| `cos_wide [1e15,3.4e38)` |          0.2522 |          2 | 1.31e+09 |
| `asin(x)/PI (naive)`     |          0.2402 |          4 |     2^32 |
| `acos(x)/PI (naive)`     |          0.0486 |          3 |     2^32 |
<!-- END generated:accuracy-bands -->

Non-ulp checks, identities, and the adversarial searches (`powfsearch`, `clogsearch`):

<!-- BEGIN generated:accuracy-other -->
```
normalize2               max |magnitude-1| 1.192093e-7
normalize3               max |magnitude-1| 1.788139e-7
normalize4               max |magnitude-1| 1.788139e-7
cexp                     max ulp re          4 im          4  worst (re,im)=-5.239318e1,6.117658e3 (5000000 samples)
clog                     max ulp re          2 im          3  worst (re,im)=1.3339183e-20,7.955754e-1 (5000000 samples)
clog(cexp(.))            max |re-back|  4.7684e-7 max |im-back|  3.0215e-7 (5000000 samples)
frexp                    bad reconstructions          0 / 20000000
ok   sin_wide^2+cos_wide^2=1            max |residual|    2.5059e-7 worst x=2.8352924e7 (n=19922075)
ok   tanh(x)=sinh(x)/cosh(x)            max |residual|    2.4976e-7 worst x=-4.222106e0 (n=10409293)
ok   exp(ln(x))=x                       max |residual|    7.4334e-6 worst x=1.88514e-40 (n=9962400)
ok   ln(exp(x))=x                       max |residual|    2.3842e-7 worst x=2.6388688e0 (n=10410597)
ok   sigmoid(x)+sigmoid(-x)=1           max |residual|    1.4901e-7 worst x=2.741174e-1 (n=19921952)
ok   erf(x)+erfc(x)=1                   max |residual|    3.0256e-7 worst x=-1.6626171e-3 (n=19921998)
ok   erf(-x)=-erf(x)                    max |residual|     0.0000e0 worst x=0e0 (n=19921996)
ok   atan2(sin(x),cos(x))=x [|x|<pi]    max |residual|    2.3842e-7 worst x=2.0715063e0 (n=10039744)
ok   log1p(x)=ln(1+x)                   max |residual|    5.9605e-8 worst x=7.9167265e-1 (n=10698265)

powfsearch: worst over all sweeps:  powf 1  powf_unchecked 1
clogsearch: worst 2 ulp at re=1.608668e-4 im=9.9999994e-1  (|v| = 9.333116015492823e-8)
clogsearch: PASS (budget 4)
```
<!-- END generated:accuracy-other -->

# benchmarks

<!-- BEGIN generated:machine -->
```
AMD RYZEN AI MAX+ 395 w/ Radeon 8060S, 32 threads
rustc 1.100.0-nightly (5ceaf6608 2026-09-25)
llvm-mca 21.1.8, -mcpu=znver5
RUSTFLAGS from .cargo/config.toml: ["-C", "target-cpu=native", "-C", "llvm-args=-force-vector-width=16"]
commit 607f1ec, 2026-09-26
```
<!-- END generated:machine -->

Serial latency: dependency chain, output feeds the next input (`examples/quickbench.rs`, min of 7; lower
is better; `nop` is the chain's own overhead).

<!-- BEGIN generated:bench-latency -->
| function                  | jodie ns | std ns | speedup |
|:--------------------------|---------:|-------:|--------:|
| `nop`                     |     0.00 |      - |       - |
| `cbrt`                    |    11.51 |  16.10 |    1.4x |
| `cbrt_unchecked`          |     8.86 |  16.10 |    1.8x |
| `cbrt_accurate`           |    14.79 |  16.10 |    1.1x |
| `cbrt_accurate_unchecked` |    12.84 |  16.10 |    1.3x |
| `cbrt_fast`               |     6.24 |  16.10 |    2.6x |
| `rcbrt`                   |    11.73 |      - |       - |
| `pow_3_2`                 |     5.04 |      - |       - |
| `pow_2_3`                 |    12.37 |      - |       - |
| `smoothstep`              |     3.94 |      - |       - |
| `smootherstep`            |     4.69 |      - |       - |
| `exp2`                    |     6.25 |   6.78 |    1.1x |
| `exp2_kf`                 |     5.08 |      - |       - |
| `exp2_checked`            |     7.66 |   6.78 |    0.9x |
| `exp10`                   |    10.21 |      - |       - |
| `exp10_checked`           |    10.11 |      - |       - |
| `log_2`                   |     8.34 |   8.17 |    1.0x |
| `log_2_unchecked`         |     7.85 |   8.17 |    1.0x |
| `sin`                     |    12.09 |  10.28 |    0.9x |
| `cos`                     |    12.45 |  10.30 |    0.8x |
| `sin_wide`                |    14.81 |  12.31 |    0.8x |
| `cos_wide`                |    14.25 |  10.81 |    0.8x |
| `tan_wide`                |    17.28 |  15.30 |    0.9x |
| `sinpi`                   |     7.34 |      - |       - |
| `cospi`                   |     8.36 |      - |       - |
| `tanpi`                   |     9.42 |      - |       - |
| `sinc_unnormalized`       |    17.07 |      - |       - |
| `ln`                      |     8.42 |   8.28 |    1.0x |
| `ln_unchecked`            |     7.59 |   8.28 |    1.1x |
| `log10`                   |     8.60 |   8.70 |    1.0x |
| `log10_unchecked`         |     7.81 |   8.70 |    1.1x |
| `log1p`                   |    11.76 |  10.47 |    0.9x |
| `log1pmx`                 |    13.55 |      - |       - |
| `log2p1`                  |    11.87 |      - |       - |
| `exp`                     |     8.50 |   7.40 |    0.9x |
| `exp_scaled`              |     8.55 |      - |       - |
| `exp_narrow`              |     7.77 |   7.40 |    1.0x |
| `exp_checked`             |    10.13 |   7.40 |    0.7x |
| `expm1`                   |     9.49 |   9.45 |    1.0x |
| `expm1_narrow`            |     8.91 |   9.45 |    1.1x |
| `exp_m1_over_x_narrow`    |    10.73 |      - |       - |
| `exp2m1`                  |     8.58 |      - |       - |
| `sinh`                    |    10.86 |  12.44 |    1.1x |
| `sinh_narrow`             |    10.11 |  12.44 |    1.2x |
| `cosh`                    |    10.46 |  12.49 |    1.2x |
| `cosh_narrow`             |     9.68 |  12.49 |    1.3x |
| `sinh_throughput`         |    11.95 |  12.44 |    1.0x |
| `cosh_throughput`         |    11.56 |  12.49 |    1.1x |
| `sinh_checked`            |    11.54 |  12.44 |    1.1x |
| `cosh_checked`            |    11.18 |  12.49 |    1.1x |
| `coshm1`                  |    12.09 |      - |       - |
| `tanh`                    |    11.53 |  10.94 |    0.9x |
| `tanh_grad`               |    13.86 |      - |       - |
| `sigmoid`                 |    11.11 |      - |       - |
| `sigmoid_fast`            |     5.40 |      - |       - |
| `sigmoid_grad`            |    13.40 |      - |       - |
| `logsigmoid`              |    14.03 |      - |       - |
| `asinh`                   |    17.80 |  15.42 |    0.9x |
| `acosh`                   |    15.72 |  14.82 |    0.9x |
| `atanh`                   |    15.05 |  13.65 |    0.9x |
| `asin`                    |     7.16 |   9.15 |    1.3x |
| `asind`                   |     8.75 |      - |       - |
| `asinpi`                  |     7.07 |      - |       - |
| `acos`                    |     9.25 |   9.40 |    1.0x |
| `acosd`                   |     9.77 |      - |       - |
| `acospi`                  |     8.44 |      - |       - |
| `atan`                    |    10.11 |  13.38 |    1.3x |
| `atan_latency`            |    11.12 |  13.38 |    1.2x |
| `atan_bounded`            |     7.03 |  13.38 |    1.9x |
| `atand`                   |    12.24 |      - |       - |
| `atan2`                   |    14.34 |  17.08 |    1.2x |
| `atan2_latency`           |    14.40 |  17.08 |    1.2x |
| `atan2_unchecked`         |    14.25 |  17.08 |    1.2x |
| `atan2_pos`               |    15.83 |  17.08 |    1.1x |
| `tan`                     |    14.95 |  13.63 |    0.9x |
| `erf`                     |    10.50 |      - |       - |
| `erfc`                    |    13.23 |      - |       - |
| `logit`                   |    11.20 |      - |       - |
| `xlogy`                   |     2.29 |      - |       - |
| `xlog1py`                 |     2.31 |      - |       - |
| `ldexp`                   |     0.00 |      - |       - |
| `frexp`                   |     8.08 |      - |       - |
| `erfinv`                  |    13.63 |      - |       - |
| `erfc_inv`                |    21.39 |      - |       - |
| `cabs`                    |    13.81 |      - |       - |
| `carg`                    |    14.17 |      - |       - |
| `cexp`                    |     9.74 |      - |       - |
| `clog`                    |    23.13 |      - |       - |
| `normalize2`              |     8.53 |      - |       - |
| `hypot3`                  |     5.28 |      - |       - |
| `rnorm3`                  |     7.22 |      - |       - |
| `normalize3`              |     8.91 |      - |       - |
| `hypot4`                  |     5.29 |      - |       - |
| `rnorm4`                  |     7.22 |      - |       - |
| `normalize4`              |     9.29 |      - |       - |
| `rsqrt`                   |     6.39 |      - |       - |
| `powf`                    |    22.71 |  12.89 |    0.6x |
| `srgb_to_linear`          |    17.85 |      - |       - |
| `linear_to_srgb`          |    15.88 |      - |       - |
| `powf_unchecked`          |    20.51 |  12.89 |    0.6x |
| `fmod`                    |     6.27 |      - |       - |
| `fmod_unchecked`          |     5.03 |      - |       - |
<!-- END generated:bench-latency -->

Throughput: independent evals over `[f32; 4096]` (`examples/quickbench.rs`, min of 7, ns/element; lower
is better).

<!-- BEGIN generated:bench-throughput -->
| function                  | jodie ns | std ns | speedup |
|:--------------------------|---------:|-------:|--------:|
| `nop`                     |    0.018 |      - |       - |
| `cbrt`                    |    0.137 |  2.383 |   17.4x |
| `cbrt_unchecked`          |    0.089 |  2.383 |   26.8x |
| `cbrt_accurate`           |    0.265 |  2.383 |    9.0x |
| `cbrt_accurate_unchecked` |    0.170 |  2.383 |   14.0x |
| `cbrt_fast`               |    0.061 |  2.383 |   39.1x |
| `rcbrt`                   |    0.136 |      - |       - |
| `pow_3_2`                 |    0.055 |      - |       - |
| `pow_2_3`                 |    0.154 |      - |       - |
| `smoothstep`              |    0.018 |      - |       - |
| `smootherstep`            |    0.032 |      - |       - |
| `exp2`                    |    0.049 |  1.183 |   24.1x |
| `exp2_kf`                 |    0.048 |      - |       - |
| `exp2_checked`            |    0.083 |  1.183 |   14.3x |
| `exp10`                   |    0.112 |      - |       - |
| `exp10_checked`           |    0.116 |      - |       - |
| `log_2`                   |    0.124 |  1.337 |   10.8x |
| `log_2_unchecked`         |    0.073 |  1.337 |   18.3x |
| `sin`                     |    0.126 |  1.486 |   11.8x |
| `cos`                     |    0.134 |  1.508 |   11.3x |
| `sin_wide`                |    0.302 |  2.629 |    8.7x |
| `cos_wide`                |    0.309 |  2.701 |    8.7x |
| `tan_wide`                |    0.414 |  4.964 |   12.0x |
| `sinpi`                   |    0.074 |      - |       - |
| `cospi`                   |    0.075 |      - |       - |
| `tanpi`                   |    0.163 |      - |       - |
| `sinc_unnormalized`       |    0.317 |      - |       - |
| `ln`                      |    0.126 |  1.372 |   10.9x |
| `ln_unchecked`            |    0.073 |  1.372 |   18.8x |
| `log10`                   |    0.128 |  2.000 |   15.6x |
| `log10_unchecked`         |    0.079 |  2.000 |   25.3x |
| `log1p`                   |    0.151 |  2.149 |   14.2x |
| `log1pmx`                 |    0.200 |      - |       - |
| `log2p1`                  |    0.157 |      - |       - |
| `exp`                     |    0.082 |  1.188 |   14.5x |
| `exp_scaled`              |    0.082 |      - |       - |
| `exp_narrow`              |    0.057 |  1.188 |   20.8x |
| `exp_checked`             |    0.098 |  1.188 |   12.1x |
| `expm1`                   |    0.079 |  2.250 |   28.5x |
| `expm1_narrow`            |    0.073 |  2.250 |   30.8x |
| `exp_m1_over_x_narrow`    |    0.096 |      - |       - |
| `exp2m1`                  |    0.083 |      - |       - |
| `sinh`                    |    0.151 |  4.748 |   31.4x |
| `sinh_narrow`             |    0.122 |  4.748 |   38.9x |
| `cosh`                    |    0.133 |  4.737 |   35.6x |
| `cosh_narrow`             |    0.098 |  4.737 |   48.3x |
| `sinh_throughput`         |    0.133 |  4.748 |   35.7x |
| `cosh_throughput`         |    0.107 |  4.737 |   44.3x |
| `sinh_checked`            |    0.169 |  4.748 |   28.1x |
| `cosh_checked`            |    0.145 |  4.737 |   32.7x |
| `coshm1`                  |    0.183 |      - |       - |
| `tanh`                    |    0.127 |  2.397 |   18.9x |
| `tanh_grad`               |    0.156 |      - |       - |
| `sigmoid`                 |    0.085 |      - |       - |
| `sigmoid_fast`            |    0.024 |      - |       - |
| `sigmoid_grad`            |    0.134 |      - |       - |
| `logsigmoid`              |    0.201 |      - |       - |
| `asinh`                   |    0.380 |  2.938 |    7.7x |
| `acosh`                   |    0.340 |  2.741 |    8.1x |
| `atanh`                   |    0.246 |  2.996 |   12.2x |
| `asin`                    |    0.096 |  3.054 |   31.8x |
| `asind`                   |    0.111 |      - |       - |
| `asinpi`                  |    0.104 |      - |       - |
| `acos`                    |    0.084 |  3.094 |   36.8x |
| `acosd`                   |    0.090 |      - |       - |
| `acospi`                  |    0.065 |      - |       - |
| `atan`                    |    0.110 |  3.529 |   32.1x |
| `atan_latency`            |    0.125 |  3.529 |   28.2x |
| `atan_bounded`            |    0.062 |  3.529 |   56.9x |
| `atand`                   |    0.125 |      - |       - |
| `atan2`                   |    0.158 |  5.918 |   37.5x |
| `atan2_latency`           |    0.179 |  5.918 |   33.1x |
| `atan2_unchecked`         |    0.155 |  5.918 |   38.2x |
| `atan2_pos`               |    0.173 |  5.918 |   34.2x |
| `tan`                     |    0.190 |  3.231 |   17.0x |
| `erf`                     |    0.156 |      - |       - |
| `erfc`                    |    0.258 |      - |       - |
| `logit`                   |    0.279 |      - |       - |
| `xlogy`                   |    0.017 |      - |       - |
| `xlog1py`                 |    0.017 |      - |       - |
| `ldexp`                   |    0.163 |      - |       - |
| `frexp`                   |    0.050 |      - |       - |
| `erfinv`                  |    0.369 |      - |       - |
| `erfc_inv`                |    0.581 |      - |       - |
| `cabs`                    |    0.116 |      - |       - |
| `carg`                    |    0.212 |      - |       - |
| `cexp`                    |    0.100 |      - |       - |
| `clog`                    |    1.489 |      - |       - |
| `normalize2`              |    0.099 |      - |       - |
| `hypot3`                  |    0.054 |      - |       - |
| `rnorm3`                  |    0.085 |      - |       - |
| `normalize3`              |    0.101 |      - |       - |
| `hypot4`                  |    0.054 |      - |       - |
| `rnorm4`                  |    0.085 |      - |       - |
| `normalize4`              |    0.103 |      - |       - |
| `rsqrt`                   |    0.097 |      - |       - |
| `powf`                    |    0.766 |  2.646 |    3.5x |
| `srgb_to_linear`          |    0.290 |      - |       - |
| `linear_to_srgb`          |    0.280 |      - |       - |
| `powf_unchecked`          |    0.605 |  2.646 |    4.4x |
| `fmod`                    |    0.076 |      - |       - |
| `fmod_unchecked`          |    0.031 |      - |       - |
<!-- END generated:bench-throughput -->

Theoretical cost from llvm-mca (`examples/mca.rs`, 100 iterations). Latency: branchless core, 64-deep
serial chain, cycles/call. Throughput: the real public function over a 16-element auto-vectorized block,
cycles/element; `?` = no throughput region (branchy, see `examples/mca_target.rs`).

<!-- BEGIN generated:mca -->
| function                  | latency (cyc) | throughput (cyc/elem) |
|:--------------------------|--------------:|----------------------:|
| `nop`                     |          1.00 |                 0.069 |
| `fast_round_int`          |         11.00 |                 0.170 |
| `std_round`               |         11.00 |                 0.140 |
| `cbrt`                    |         32.80 |                 1.371 |
| `cbrt_wrapped`            |         44.02 |                     ? |
| `cbrt_unchecked`          |         32.80 |                 1.067 |
| `cbrt_accurate`           |         56.00 |                 2.092 |
| `cbrt_accurate_unchecked` |         56.00 |                 1.591 |
| `cbrt_fast`               |         27.03 |                 0.592 |
| `rcbrt`                   |         61.00 |                 1.239 |
| `pow_3_2`                 |         22.00 |                 0.637 |
| `pow_2_3`                 |         49.13 |                 1.364 |
| `smoothstep`              |         29.00 |                 0.704 |
| `smootherstep`            |         33.00 |                 0.893 |
| `exp2`                    |         27.00 |                 0.642 |
| `exp2_kf`                 |         21.02 |                 0.647 |
| `exp2_checked`            |         34.02 |                 2.318 |
| `exp10`                   |         45.00 |                 1.147 |
| `exp10_checked`           |         46.02 |                 3.444 |
| `log2`                    |         36.02 |                 1.488 |
| `log2_unchecked`          |         36.02 |                 0.929 |
| `sin`                     |         58.98 |                 1.107 |
| `sin_wide`                |         53.31 |                 2.207 |
| `cos_wide`                |         51.86 |                 2.508 |
| `cos`                     |         56.98 |                 1.411 |
| `sinpi`                   |         32.03 |                 0.923 |
| `cospi`                   |         36.98 |                 0.902 |
| `tanpi`                   |         77.98 |                 2.118 |
| `sinc_unnormalized`       |         62.31 |                 2.506 |
| `ln`                      |         47.19 |                 1.574 |
| `ln_unchecked`            |         36.02 |                 0.926 |
| `log10`                   |         47.19 |                 1.552 |
| `log10_unchecked`         |         36.03 |                 0.988 |
| `log1p`                   |         45.00 |                 1.521 |
| `log1pmx`                 |         53.11 |                 2.389 |
| `log2p1`                  |         46.00 |                 1.706 |
| `log10p1`                 |         46.06 |                 1.768 |
| `exp`                     |         38.00 |                 1.071 |
| `exp_scaled`              |         38.00 |                 1.080 |
| `exp_narrow`              |         35.00 |                 0.807 |
| `exp_checked`             |         46.02 |                 3.319 |
| `expm1`                   |         43.02 |                 0.951 |
| `expm1_narrow`            |         40.02 |                 0.911 |
| `expm1_checked`           |         47.02 |                 1.140 |
| `exp_m1_over_x_narrow`    |         67.03 |                 1.351 |
| `exp2m1`                  |         49.00 |                 0.894 |
| `exp10m1`                 |         52.03 |                 1.187 |
| `sinh`                    |         48.00 |                 1.732 |
| `sinh_narrow`             |         45.00 |                 1.489 |
| `cosh`                    |         47.00 |                 1.417 |
| `cosh_narrow`             |         44.00 |                 1.161 |
| `sinh_throughput_fn`      |         56.00 |                 1.704 |
| `cosh_throughput_fn`      |         55.00 |                 1.380 |
| `sinh_checked`            |         52.00 |                 1.917 |
| `cosh_checked`            |         51.00 |                 3.694 |
| `coshm1`                  |         59.00 |                 2.186 |
| `tanh`                    |         61.02 |                 3.509 |
| `tanh_grad`               |         66.13 |                 1.656 |
| `sigmoid`                 |         53.00 |                 3.258 |
| `sigmoid_fast`            |         22.00 |                 1.754 |
| `sigmoid_grad`            |         64.11 |                 1.597 |
| `logsigmoid`              |         76.03 |                 2.108 |
| `logsigmoid_checked`      |         70.22 |                 2.166 |
| `gelu`                    |         67.50 |                 3.241 |
| `silu`                    |         56.13 |                 1.347 |
| `silu_checked`            |         58.99 |                 1.400 |
| `softsign`                |         19.02 |                 0.449 |
| `sqrt1pm1`                |         36.02 |                 1.074 |
| `asinh`                   |         72.44 |                 3.530 |
| `acosh`                   |         95.13 |                 3.326 |
| `atanh`                   |        101.03 |                 2.284 |
| `asin`                    |         61.66 |                 1.447 |
| `asind`                   |         74.92 |                 1.436 |
| `asinpi`                  |         65.86 |                 1.416 |
| `acos`                    |         43.05 |                 1.242 |
| `acosd`                   |         50.05 |                 1.312 |
| `acospi`                  |         36.00 |                 0.901 |
| `atan`                    |         48.28 |                 1.434 |
| `atan_latency`            |         49.02 |                 1.399 |
| `atan_bounded`            |         31.99 |                 0.766 |
| `atand`                   |         56.11 |                 1.502 |
| `atan2`                   |         52.14 |                 1.449 |
| `atan2_unchecked`         |         52.11 |                 1.398 |
| `atan2_latency`           |         52.05 |                 1.464 |
| `atan2_pos`               |         56.23 |                 1.475 |
| `tan`                     |         71.00 |                 5.571 |
| `tan_wide`                |         65.25 |                 2.997 |
| `erf`                     |         81.05 |                 2.114 |
| `erfc`                    |         61.84 |                 3.087 |
| `logit`                   |         60.14 |                 3.516 |
| `xlogy`                   |         60.73 |                 1.518 |
| `xlog1py`                 |         48.11 |                 1.583 |
| `erfinv`                  |         48.13 |                 3.843 |
| `erfc_inv`                |        117.81 |                 5.217 |
| `clog_re`                 |        264.89 |                 9.264 |
| `cabs`                    |         36.13 |                 0.862 |
| `carg`                    |         63.72 |                 1.852 |
| `normalize2`              |         40.02 |                 1.103 |
| `hypot3`                  |         23.13 |                 0.796 |
| `rnorm3`                  |         34.13 |                 1.140 |
| `normalize3`              |         43.02 |                 1.186 |
| `hypot4`                  |         23.13 |                 0.796 |
| `rnorm4`                  |         34.13 |                 1.140 |
| `normalize4`              |         46.02 |                 1.191 |
| `diff_of_products`        |         11.00 |                 0.236 |
| `cross2`                  |         11.00 |                 0.236 |
| `rsqrt`                   |         30.00 |                 1.072 |
| `powf`                    |        103.31 |                 5.218 |
| `powf_pos`                |         96.80 |                 5.154 |
| `srgb_to_linear`          |         55.02 |                 2.605 |
| `linear_to_srgb`          |         74.30 |                 2.431 |
| `signed_pow`              |         98.75 |                 4.907 |
| `powf_unchecked`          |         92.59 |                 4.781 |
| `fmod`                    |         26.08 |                     ? |
| `fmod_checked`            |          8.88 |                     ? |
| `fmod_unchecked`          |         22.00 |                 0.453 |
| `rem_euclid`              |         30.08 |                     ? |
| `div_euclid`              |         28.00 |                     ? |
<!-- END generated:mca -->

# tools
- `cargo run --release --example accuracy [thorough] [filter]` - avg/max ulp against an f64
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

- `cargo run --release --example worst_corpus` - every public 1-arg function (90, shared list in `examples/support/unary_fns.rs`) x 90 historically-hard inputs
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
