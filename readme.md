# faster? f32 math functions
Attempting to provide faster implementations of common f32 math functions with a similar level of accuracy to the standard library.
There are also perfectly rounded variants and faster-but-sloppier variants.

All functions are full-range correct. Negatives, denormals, zero, inf and nan
are handled, except exp2, where the default is the fast
unchecked version whicj is only valid for x in [-126, 128). exp2_checked handles the full
range (overflow to inf, denormal underflow, nan) for ~2.5 ns extra latency.
cbrt_accurate is perfectly rounded on every input tested.

All functions auto-vectorize, it's a hard requirement

# precision (see examples/accuracy.rs)
```
              | jodie avg | jodie max | std avg | std max
--------------|-----------|-----------|---------|--------
         cbrt |   0.329   |     2     |    0    |    0
cbrt_accurate |   0.000   |     0     |    0    |    0
         exp2 |   0.037   |     1     |  0.000  |    1
 exp2_checked |   0.037   |     1     |  0.000  |    1
         log2 |   0.006   |     2     |  0.000  |    1
 sin (|x|<1e3)|   0.027   |     2     |  0.002  |    1
 cos (|x|<1e3)|   0.082   |     2     |  0.001  |    1
```

# benchmarks
Run on i5-1145G7, -C target-cpu=native (now set in .cargo/config.toml)
```
Serial latency (dependency chain, examples/quickbench.rs; lower is better)
              | jodie   | std     | improvement
--------------|---------|---------|------------
         cbrt | 12.9 ns | 22.0 ns | 1.7x
cbrt_accurate | 19.5 ns | 22.0 ns | 1.1x
          cos | 18.5 ns | 18.5 ns | 1.0x
         exp2 |  8.5 ns | 13.3 ns | 1.6x
 exp2_checked | 13.3 ns | 13.3 ns | 1.0x
         log2 | 13.1 ns | 14.9 ns | 1.1x
          sin | 15.9 ns | 18.7 ns | 1.2x
```
```
Throughput (independent array evals over [f32; 4096], examples/quickbench.rs; lower is better)
              | jodie    | std     | improvement
--------------|----------|---------|------------
         cbrt | 0.37 ns  | 3.93 ns | 10.6x
cbrt_accurate | 0.69 ns  | 3.93 ns | 5.7x
          cos | 0.40 ns  | 4.13 ns | 10.4x
         exp2 | 0.23 ns  | 3.13 ns | 13.7x
 exp2_checked | 0.48 ns  | 3.13 ns | 6.6x
         log2 | 0.53 ns  | 4.10 ns | 7.8x
          sin | 0.31 ns  | 4.26 ns | 13.8x
```
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
cbrt_accurate       |          63.06 |             3.132
exp2                |          35.00 |             0.841
exp2_checked        |          43.06 |             1.399
log2                |          34.23 |             1.556
sin                 |          46.00 |             1.022
cos                 |          54.00 |             1.360
```

# tools
- `cargo run --release --example accuracy [filter]` - avg/max ulp over dense strided sweeps of each domain
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
