# faster? f32 math functions
Attempting to provide faster implementations of common f32 math functions with a similar level of accuracy to the standard library.
There are also perfectly rounded variants and faster-but-sloppier variants.

cbrt, cbrt_accurate, log2 and exp2_checked are full-range correct: negatives,
denormals, zero, inf and nan are all handled. exp2's default is the fast
unchecked version, only valid for x in [-126, 128); exp2_checked handles the
full range (overflow to inf, denormal underflow, nan) for ~2.5 ns extra
latency. cbrt_accurate is within 1 ulp on every input tested (exhaustively,
see below) -- **not** perfectly rounded on every input, despite an earlier
claim here: the full 2^32-pattern sweep found a systematic 1-ulp miss at one
specific mantissa, recurring at every octave from ~2^-126 to ~2^127.
**sin/cos still have a real domain limit, though much wider than before**:
q = round(x/pi) must be an exact integer for the reduction to land
accurately, and a single f32 only represents integers exactly up to 2^24 --
past that the old reduction silently mis-rounded q and, for large enough x,
returned outright `inf` for an ordinary finite input. Fixed by splitting q
into an exact double-float integer pair and pi into 3 error-free-transform
words (not the crate's other PI_A..D, which use trailing-zero padding valid
only for *bounded* q -- these use genuine `two_prod`/`two_sum`, exact for
any magnitude). Fully accurate out to ~1e12-1e13 (vs. the old code's ~5e7),
then degrades *gradually* -- bounded, finite output, not `inf` -- out past
1e19. This is a fundamental limit of any fixed-width reduction, not a bug:
q needs `bits(x) + ~24` bits to stay exact, and no *finite* width covers
every f32; std's sin/cos stay accurate to f32::MAX using a much heavier
reduction we didn't replicate here (verified against an arbitrary-precision
reference). An f64-based reduction was tried first (simpler, one hardware
fma) and worked, but was only fully accurate to ~1e9 and cost more
throughput than the all-f32 version below, since f64 vectors are half the
width of f32 ones on this hardware.

All functions auto-vectorize, it's a hard requirement

# precision (see examples/accuracy.rs)
Fuzz-mode (100M random f32 bit patterns/function); pass `thorough` for an
exhaustive sweep of all 2^32 patterns instead (few minutes, needs --release).
```
                  | jodie avg  | jodie max | std avg | std max
------------------|------------|-----------|---------|--------
             cbrt |    0.326   |     3     |    0    |    0
    cbrt_accurate |    0.000   |     1     |    0    |    0
             exp2 |    0.030   |     1     |  0.000  |    1
     exp2_checked |    0.016   |     1     |  0.000  |    1
             log2 |    0.003   |     3     |  0.000  |    1
     sin (|x|<1e3)|    0.021   |     2     |  0.000  |    1
     cos (|x|<1e3)|    0.066   |     2     |  0.000  |    1
    sin (all f32) | (degrades gradually above |x| ~ 1e12-1e13 -- see note above; 0.007/1 for std)
    cos (all f32) | (degrades gradually above |x| ~ 1e12-1e13 -- see note above; 0.007/1 for std)
```
sin's magnitude-bucketed avg/max ulp (2M random samples/bucket), showing how
that degradation actually looks -- low and flat for a long stretch, then a
rising tail of rare bad cycles, well before it turns into the old code's
outright `inf`:
```
        range | avg ulp | max ulp
--------------|---------|--------
     [1e6,1e7)|   0.252 |      2
     [1e9,1e10)|  0.252 |     48
    [1e11,1e12)|  0.254 |    159
    [1e12,1e13)|  0.370 | 154382
    [1e13,1e14)|  1.218 | 492480
    [1e14,1e15)| 679010  |   2.1e9
```

# benchmarks
Run on i5-1145G7, -C target-cpu=native (now set in .cargo/config.toml)
```
Serial latency (dependency chain, examples/quickbench.rs; lower is better)
              | jodie   | std     | improvement
--------------|---------|---------|------------
         cbrt | 12.9 ns | 22.0 ns | 1.7x
cbrt_accurate | 17.0 ns | 22.0 ns | 1.3x
          cos | 44.1 ns | 17.9 ns | 0.4x
         exp2 |  8.5 ns | 13.3 ns | 1.6x
 exp2_checked | 13.3 ns | 13.3 ns | 1.0x
         log2 | 13.1 ns | 14.9 ns | 1.1x
          sin | 42.5 ns | 18.1 ns | 0.4x
```
sin/cos's latency is now clearly worse than std (was faster before):
extending the accurate domain by ~5 orders of magnitude means enough extra
serial arithmetic (double-float q and pi, several compensated
add/multiply steps) that the reduction alone costs more than std's whole
implementation. Deliberate trade for going from "returns `inf` on ordinary
finite input past ~1e10" to "accurate past 1e12, degrades gracefully to
1e19" -- see the accuracy note above.
```
Throughput (independent array evals over [f32; 4096], examples/quickbench.rs; lower is better)
              | jodie    | std     | improvement
--------------|----------|---------|------------
         cbrt | 0.37 ns  | 3.93 ns | 10.6x
cbrt_accurate | 0.67 ns  | 3.93 ns | 5.9x
          cos | 2.56 ns  | 3.99 ns | 1.6x
         exp2 | 0.23 ns  | 3.13 ns | 13.7x
 exp2_checked | 0.48 ns  | 3.13 ns | 6.6x
         log2 | 0.53 ns  | 4.10 ns | 7.8x
          sin | 2.52 ns  | 4.11 ns | 1.6x
```
sin/cos's throughput went from ~13-14x std to ~1.6x -- still ahead, but
narrowly. This went through several iterations (see jodiemath-workflow
memory for the full trail): an f64-based reduction worked and was cheaper
(0.7 ns/1.6x std throughput) but only fully accurate to ~1e9; switching to
an all-f32 double-float reduction pushed accuracy out to ~1e12-1e13 at
this higher cost. Two vectorization/perf pitfalls along the way:
- a `q as i64` cast for the parity bit doesn't vectorize at all -- Rust's
  float-to-int cast is saturating, so LLVM falls back to a scalar
  convert-with-NaN/range-check per lane (~9x slower than std). Fixed by
  computing parity with plain float ops (`floor`-based "mod 2") instead.
- combining several small correction terms via plain adds first (instead
  of feeding all 7 through a sequential compensated-sum loop) cut the loop
  from 7 iterations to 4 with no accuracy cost -- verified bit-for-bit
  identical output. A *further* restructuring meant to shorten the serial
  dependency chain per element actually made both latency and throughput
  slightly worse: mca's bottleneck-analysis showed this region is
  resource-pressure-bound (16 independent vectorized elements in flight),
  not dependency-chain-bound, so shortening one element's chain doesn't
  help when the ports are already the binding constraint.
The throughput gap vs std comes almost entirely from vectorization: std's
functions have branches, so LLVM can't vectorize loops that call them.
Absolute ns swing session-to-session with CPU thermal state (the laptop
throttles up to ~2.5x mid-session) -- only trust jodie-vs-std ratios measured
in the same run. examples/mca.rs gives a thermal-noise-free second opinion in
cycles instead of ns:
```
theoretical cost from llvm-mca (-mcpu=native, 100 iterations)
                    | latency (cyc) | throughput (cyc)
--------------------|----------------|------------------
cbrt                |          35.06 |             1.629
cbrt_accurate       |          59.06 |             3.129
exp2                |          35.00 |             0.841
exp2_checked        |          43.06 |             1.399
log2                |          34.23 |             1.556
sin                 |         128.00 |             7.354
cos                 |         132.00 |             6.223
```

# tools
- `cargo run --release --example accuracy [thorough] [filter]` - avg/max ulp against an f64 reference.
  Default mode fuzzes 100M random f32 bit patterns per function (a few seconds); `thorough` exhaustively
  sweeps all 2^32 bit patterns instead (every denormal, every NaN payload, both signs -- a few minutes,
  multi-threaded, refuses to run in a debug build)
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
- add inverse trig functions
- add tan()
- perfectly rounded versions
