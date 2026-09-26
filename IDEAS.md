# IDEAS.md

Open ideas that nobody has implemented and measured yet. Everything that has been tried, shipped or
rejected, is in `graveyard.md`. Read the graveyard before starting anything: re-running an old negative
result is the most common way to waste a session here.

## How to screen an idea

### Is there any fit left in this polynomial?

1. **Oracle test first.** A better polynomial fit only helps if the fit is what limits the final error.
   Plug a clearly better fit (a higher degree, or an older, larger version from git) into the real
   function and measure. If the error doesn't drop, stop. Example: `ln_normal`'s old degree-9 polynomial
   fits 22x better and measures worse. The strongest version of this test feeds the chain a correctly
   rounded polynomial value, which shows the best any fit could do (for `ln_normal`: max 1 / avg 0.272
   against the shipped 3 / 0.453).
2. **Headroom ratio.** Compare the error of the shipped f32 coefficients with a fresh ulp-weighted minimax
   fit of the same degree, rounded to f32. Weight both by how much the final combine amplifies polynomial
   error. A ratio of 1 or less means the shipped coefficients already win. Always look at the absolute
   error too: `atan_poly` has a ratio of 14.2 but an absolute error of 0.228 ulp, so there's nothing to
   gain. Nine polynomials are closed on this test (table in the graveyard).
3. **How the coefficients were tuned predicts the outcome.** Tuned against the real function: nothing to
   gain. Least-squares fit: a minimax refit always trades avg for max (4 of 4 cases). A rescaled copy of
   another polynomial's coefficients, or a least-squares fit whose error sits where samples are sparse:
   usually improves both avg and max.

### Is this speed change real?

llvm-mca's throughput number is reproducible but often wrong (4 of 4 apparent regressions in one session
were artifacts). Go down this list only as far as you need to:

1. Instruction count of the `LLVM-MCA-BEGIN <name>` ... `LLVM-MCA-END` region in the newest
   `mca_target*.s` (`tools/mca_region.py <region>` finds it and prints steps 1-3).
2. Opcode histogram. Same expensive ops (`vdivps`, `vsqrtps`, fma, `vmulps`) and fewer total means mca is
   wrong.
3. uOps and Block RThroughput (`llvm-mca --bottleneck-analysis`). Fewer instructions and fewer uOps with
   flat RThroughput: take it. Fewer instructions but more uOps: a real regression. This is how `erfinv`'s
   polynomial fold was correctly rejected and `acosh`'s select merge correctly adopted.
4. Wall clock last. Run-to-run noise is about ±8%, so it can't resolve changes under 10%. For small
   deltas use `perf stat -e cycles` on a tight loop, alternating old and new builds (`tools/ab.sh` does
   this for the trig functions).

What llvm-mca gets right and wrong (measured with perf-stat cycles on the old i5-1145G7): packed f64 code
runs within 0.93-1.24x of mca's estimate, so don't avoid f64 on the theory that it's secretly slow. Gathers
are the exception: mca charged 1.5 cycles for three `vpgatherdd` that really cost about 10.5, which made
the old table-based `sin_wide` 4.5x slower than predicted.

Don't run `examples/mca` at the same time as another `--emit=asm` build in the same worktree; they write
the same file.

### Is this accuracy change real?

- **Diff the function's asm first.** If the code didn't change, the accuracy change isn't real. `cos2pi`'s
  quick-fuzz avg once moved 59% and its max 64x on identical code.
- **Quick-mode max is too optimistic, not just noisy.** It once reported `sin`'s true max of 58 as 4. Any
  max needs `accuracy thorough`. Quick-mode avg is trustworthy to about ±0.0001.
- Never tune on a sample of a domain you could test exhaustively. The tuner finds exactly the points you
  skipped.

### Where wins usually come from

- **A guard can prove a bound, so the code behind it doesn't need the general version.** Look for this
  through call chains (`srgb_to_linear` reaches `log_2` two levels down through `powf_pos`), by code shape
  rather than name (`tanh` has its own copy of a polynomial that a grep for the macro misses), and
  downstream (`compound`'s `log1p` zero-select didn't matter because `exp_checked` mapped both zeros to
  `1.0`). You can also move a clamp to create the bound: `erfc` already clamped `|x|` before squaring, and
  tightening that clamp (11 to 10.21) proved the exponent in range, so `exp_checked`'s own clamp could go.
  Same output, 11.7% faster. Any function that clamps an argument and then calls a checked function is a
  candidate.
- **Remove a rounding, not just an instruction.** Folding a polynomial's top coefficients one level lower,
  so `x^4` is never formed, saves a multiply and a rounding. It shrank 30 regions and improved accuracy
  everywhere it applied, including `tanh`'s max 6 to 5, which three refits had failed to reach.

## Open: accuracy

These improve a function without costing speed, so any win is worth taking.

- **`acosd` should use the two-word `RAD_TO_DEG`.** It still computes `acos(x) * RAD_TO_DEG_HI`, and
  `RAD_TO_DEG_HI` is 0.46 ulp low, which biases every result. The comment says the low word doesn't help
  because `acos`'s error dominates. That was true before `acos` was rewritten and isn't now. Measured over
  every f32 in `[-1, 1]`: two-word version avg 0.0527 / max 3, against the shipped 0.1291 / 4. Costs one
  fma and one multiply (not yet priced with mca). `asind` and `atand` already do this.
- **Two-word `LOG2_E` in `log_2_normal`.** The same change for `log10_normal` was measured and rejected
  (graveyard): the low word only helps after a refit, needs degree 8, and costs 2 fma and 11-12% on every
  mca metric for `log10p1` max 2 to 1. `log_2` has less to gain (smaller constant error; `log2` is already
  max 1), so expect the same cost for less. Run the oracle test first.
- **Degree 8 for `ln_normal`'s `Q`.** Measured and not shipped: avg 0.0063 against 0.0073, +3
  instructions, Block RThroughput 16 to 17 on `ln_unchecked`, no max improvement. The coefficients are in
  the graveyard if a caller ever wants the extra avg.
- **`erfc_inv`'s remaining average comes from its argument path.** In the far branch the polynomial fit
  contributes about 0.26 avg ulp and the roundings feeding it about 0.28: `v = sqrt(w)` costs 0.262, `w`'s
  own rounding 0.130. Two unmeasured options:
  - Take the square root last: `z = sqrt(w * S(v))` with `S = (1 + Q)^2` fitted directly, so the final
    rounding lands on the answer. Costs one sqrt and one multiply, and the two sqrts run in series, so the
    question is speed, not accuracy. Estimated avg ~0.35 to ~0.22.
  - A polynomial in `1/sqrt(w)` fits about 8x better but needs a divide and only fixes the fit half
    (~0.35 to ~0.29). The first option is better value.
- **`erfinv_tail_poly_m1` can't be screened with an LP in the monomial basis.** It's degree 11 in
  `t = sqrt(w) - 1` with `t` up to 3, so the Vandermonde matrix spans `3^11` and HiGHS returns a
  "minimax" (5.11 ulp) worse than the shipped coefficients (0.875). Use a Chebyshev basis or rescale `t`.
- **`erfc` (and `gelu`, which inherits `erfc`'s max 6 at `x = -1.32`) is limited by how `erfcx_pos`
  evaluates its polynomial**, not by the fit or its argument. Horner order is rejected: +19-34% latency,
  and only ~0.001 ulp better. What's untried is an order between the current Estrin and Horner. Any
  reorder needs a refit: the coefficients were tuned against the exact current order, and reordering
  alone breaks the `erfc(0)` edgecheck pins.
- **`atan2` is limited by its `y/x` division.** Measured on the (since removed) `atan2pi`, whose `x > 0`
  arm is the same code: 86% of that arm's error is the rounding of `y/x`. A compensated two-word `y/x` is
  the only lever left. Unpriced, and watch for `1/x` overflowing at legitimately small `x`.
- **Denormal intermediates in `x * f(x)` products.** When `f` decays exponentially there's a window where
  `f(x)` is denormal but `x * f(x)` isn't, and the result loses bits. `gelu` (fixed: max 9 to 6, 13.7%
  faster) and `silu` had this. To screen a function, compare where `f` underflows with where `x * f(x)`
  does; if they differ, reorder the multiply. Not yet checked: `xlogy`, `xlog1py`, `sinc_unnormalized`,
  and anything else multiplying a decaying function by a growing one. `examples/denormal_audit.rs` only
  checks outputs, so it can't see this.
- **`sin`/`cos`/`tan_wide` on `|x| >= 1` average ~0.25 ulp** (the exhaustive 0.13 is diluted by the
  `|x| < 1` half). That's the f32 rounding of the reduced argument `r` plus the polynomial. An `r_lo`
  correction term would cost about 3 ops.

## Open: speed

- **Wide trig tier.** After the two gather-free passes (graveyard, 2026-09-26):
  - AVX2 is about 2.4x slower per element than AVX-512 for half the lanes. Each funnel shift is 4 ops
    without `vpshldvd`, and about 20 constants spill to the stack. LLVM has turned every source-level
    trick tried so far back into the same code.
  - Scalar latency is still ~3 ns worse than the old table version (window select and int-to-float).
  - Measure with `examples/trig_bench.rs`, not quickbench: quickbench's band mask lets LLVM constant-fold
    the window.
- **`exp_r_poly!`: compute `e^r - 1` and fold the `+ 1` into the caller.** `exp_reduce!` already does this.
  All six call sites (`exp`, `exp_scaled`, `exp_narrow`, `sigmoid`, `exp_neg_scaled64!`, `silu_checked`)
  end in a multiply by an exact power of two, so each can absorb the `+1` as `fma(s, t, t)`. On the
  polynomial alone: max 2.392 to 2.148, avg 0.7165 to 0.6783, no extra instructions, and at the call sites
  it probably saves one add. `exp_r_poly!` is a core kernel, so this needs `jm core claim` and a
  whole-crate accuracy sweep. Don't combine it with a degree-6 bump (measured: 5-23% slower).
- **`silu_checked` still inlines the `exp_neg_scaled64!` reduction.** It matches the macro except for the
  final `* 2^-64`. Switch it over and check the asm region is byte-identical.
- **Two-for-one `sincos`.** Share the reduction and `r^2` between sin and cos when the caller wants both.
  Only worth it as a slice API; fusing scalar calls has already failed (graveyard).

### Before writing any f64 kernel

Several functions moved an f32 double-float chain to f64 (details in the graveyard). What applies to all of
them, measured on the i5-1145G7:

1. **f64 costs 2x per lane.** `vfmadd213pd zmm` has Block RThroughput 1.0, `vfmadd213ps ymm` 0.5, both
   8 lanes. f64 only wins when it makes the algorithm shorter by dropping error-free transforms. Halving
   the instruction count is a big win (`remainder_wide` -70%). Cutting 10-18% is a wash (`powf`). Know
   which case you're in before writing code. Default to `Df32`/fma for short wide steps.
2. **f64 divides are slow** (`vdivpd zmm` RThroughput 16). Use a seed plus a Newton step (~6 ops) unless
   exact division is the point.
3. **f64 constants cause register pressure.** One `powf` draft had 23 and LLVM reloaded 24
   `vbroadcastsd` inside the loop. Trim polynomials to the accuracy actually needed.
4. **Never use f64 tables or f64 gathers.** They drop the loop to 4 lanes and ran about 4.5x worse than
   predicted. `u32` tables are the only gather shape allowed.
5. **Keep f64 paths branch-free.** Scalar f64 costs about 15x the packed form (`powf`, measured both ways).
   Use selects, and split `_checked`/`_unchecked` versions rather than branching per lane.

Don't retry `cbrt_accurate` in f64: its error-free block is 9% of the region and it's already correctly
rounded.

## Open: fitting tools

- **Search the joint f32 rounding of a coefficient set.** The headroom test rounds each coefficient to f32
  independently, so it can't see this. Two ways: a mixed-integer fit (HiGHS supports it; a homegrown
  `fpminimax`, since Sollya isn't installed) or LLL lattice reduction (6-10 coefficients, small enough to
  hand-roll). Visible loss to target: `acos_poly` goes 1.49 to 1.54 when rounded, `acospi_poly` 1.61 to
  1.62. Best targets are polynomials already tuned against the real function, where rounding is the only
  thing left.
- **Sollya `fpminimax`** (not installed). Candidate: `acos_poly`.
- **Exhaustive rounding-interval LP (rlibm-style) for `log_2`, then `ln`/`log10`.** Tractable (2^23
  distinct reduced inputs; `exp2` has 2^29 and isn't). Much less attractive since `log_2` reached max 1:
  the most it can buy is correctly rounded instead of faithful.
- **Reusable ulp-weighted minimax fitting.** Weight by `1/ulp(f(x))` and fit the target the combine
  actually needs, not the textbook function. `asin_poly`'s regression came from fitting
  `acos(a)/sqrt(1-a)` when the chain wants `(fl(pi/2) - asin(a))/sqrt(1-a)`, a 4.4e-8 difference that is
  1.47 ulp at the branch edge.
- **Fit in a transformed variable.** `erfcx_pos`'s `v = 1/(2 + x)` took `erfc` from max 109 to 7 where
  every same-variable refit failed: a reciprocal maps the half-line to a finite interval, and a leading `v`
  factor gets the asymptote right automatically. Worth trying for anything with a `c/x` tail (e.g.
  `erfc_inv`'s tail).
- **Denser fit grids near output powers of two**, where the ulp size doubles. Or Sobol grids, which avoid
  aliasing against the ulp steps.
- **In-repo Remez exchange** with certified equioscillation. Candidate: `atan_poly`.
- **Tuner improvements in `tune.rs`:**
  - Scan each polynomial's final combine constant ±a few hundred ulp.
  - Search branch thresholds and coefficients together, instead of retuning thresholds against frozen
    polynomials.
  - Restart from several perturbed seeds and keep the best.
- **`exp2_q_poly` weighted by its combine's derivative.** First check where the worst `f` is: this method
  fails when the derivative vanishes exactly at the hard region.
- **Alternate coefficient sets behind a cargo feature**, weighted for a common input range (e.g. `sin` on
  `[-2pi, 2pi]`). Same code and cost, different constants.

## Open: new functions

- **`tgamma`/`lgamma`** (Lanczos or Stirling).
- **More `_approx` functions** (`sin_approx`, `tanh_approx`) for ML inference, outside the normal ulp
  budget. State both the error and the domain: in the existing tier, the unrefined bit tricks
  (`sqrt_approx` 4.5%, `rcp_approx` 6.6%) are about equally accurate everywhere, while the refined
  `cbrt_approx` (5e-6) fails completely at the smallest normal. A Newton step buys accuracy and a domain
  edge together. Add each one to `examples/approx_bounds.rs`.

## Open: tooling and harness

- **2-arg worst-case search for `atan2` and `fmod`.** Uniform sampling rarely hits the hard inputs.
  `examples/powfsearch.rs` shows how: read the hard region off the error model, then sweep it
  exhaustively. It found twice the max that a 10M-sample fuzz did.
- **`accuracy.rs` per-branch mode**: report which branch produced each worst case.
- **`denormal_audit` only covers 1-arg functions.** `atan2`, `powf`, `xlogy`, `xlog1py` and `fmod` are
  missing. Each needs a fixed second argument and a reference.
- **`saturation_pins` tests `exp2m1` at the wrong bound.** It probes `-151.0`, but `exp2m1` clamps at
  `-126.0`, so the real lower transition is never tested. Change it to `-126.0` (it should pass either way,
  since the value is `-1.0` on both sides).
- **Controlled wall-clock benchmarking** (pinned CPU frequency, alternating A/B batches). mca and wall
  clock have disagreed on direction at least once.
- **Other targets.** Every scheduling decision here was tuned on x86 (the i5-1145G7, now Zen 5).
  Division-vs-polynomial and Estrin groupings need re-measuring before claiming anything on aarch64/NEON.
