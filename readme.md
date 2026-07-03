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
         cbrt |   0.085   |     1     |    0    |    0
cbrt_accurate |   0.000   |     0     |    0    |    0
         exp2 |   0.059   |     2     |  0.000  |    1
 exp2_checked |   0.059   |     2     |  0.000  |    1
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
         cbrt | 11.6 ns | 21.2 ns | 1.8x
cbrt_accurate | 20.0 ns | 21.2 ns | 1.1x
          cos | 12.7 ns | 12.6 ns | 1.0x
         exp2 |  8.5 ns |  9.6 ns | 1.1x
 exp2_checked | 11.0 ns |  9.6 ns | 0.9x
         log2 |  8.9 ns | 10.0 ns | 1.1x
          sin | 10.8 ns | 12.8 ns | 1.2x
```
```
Throughput (independent array evals over [f32; 4096], examples/quickbench.rs; lower is better)
              | jodie    | std     | improvement
--------------|----------|---------|------------
         cbrt | 0.60 ns  | 3.89 ns | 6.5x
cbrt_accurate | 0.69 ns  | 3.89 ns | 5.6x
          cos | 0.27 ns  | 2.90 ns | 10.7x
         exp2 | 0.23 ns  | 2.27 ns | 9.9x
 exp2_checked | 0.38 ns  | 2.27 ns | 6.0x
         log2 | 0.36 ns  | 2.75 ns | 7.7x
          sin | 0.21 ns  | 2.98 ns | 13.9x
```
The throughput gap vs std comes almost entirely from vectorization: std's
functions have branches, so LLVM can't vectorize loops that call them.

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
