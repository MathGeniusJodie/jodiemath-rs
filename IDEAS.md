# IDEAS.md

What is still open. Everything tried -- shipped or rejected -- moved to
`graveyard.md` on 2026-07-28. **Read the graveyard before proposing
anything**: this crate has a long history of measured negative results, and
re-running one of them is the most common way to waste a session.

Nothing below has been implemented and measured. Entries carry whatever
screening they have already had.

## How to screen an idea here

The cheap checks, in the order that has actually saved the most time.

### Is there any fit left in this poly?

1. **Oracle screen first -- one run, usually dispositive.** Headroom only
   converts if the fit is the *binding* term. Score a known-better-fitting
   poly through the **real chain**: a higher degree, or the pre-shed
   ancestor in git history. If the real error does not move, stop.
   `ln_normal`'s pre-shed degree 9 is a **22x better fit and measures
   worse**. The strongest form adds an oracle row -- feed the chain a
   correctly-rounded `p` and read off the floor (for `ln_normal`, max 1 /
   avg 0.272 against a shipped 3 / 0.453: two of three ulps were never the
   polynomial's to give).
2. **Headroom ratio.** Exact-arithmetic error of the shipped f32
   coefficients, weighted by the final combine's own sensitivity
   `|d(result)/dP| / ulp(result)`, against the same-degree ulp-weighted LP
   minimax re-quantised to f32. Ratio <= 1 means the shipped coefficients
   already beat a fresh LP optimum. **Read the ratio together with the
   absolute** -- `atan_poly` is 14.2x but its absolute is 0.228 ulp, so
   there is nothing there. The full table is in the graveyard; nine polys
   are closed on it.
3. **Objective class predicts whether headroom converts.** Real-chain-tuned
   -> nothing to take. Least-squares-tuned -> a minimax refit *always*
   trades avg for max (4 of 4). A rescale of another poly's coefficients,
   or an LS fit whose error concentrates where the weight is high but
   sample density is low -> wins on both axes.

### Is this perf change real?

Escalate only as far as you need to. mca's throughput column is
reproducible but often wrong here -- 4 of 4 apparent regressions in one
session were artifacts.

1. **Instruction count** of the region between `LLVM-MCA-BEGIN <name>` and
   `LLVM-MCA-END` in the newest `target/release/examples/mca_target-*.s`.
2. **Opcode histogram.** Identical expensive ops (`vdivps`/`vsqrtps`/fma/
   `vmulps`) plus fewer total = mca is wrong.
3. **uOps + `Block RThroughput`** via `llvm-mca --bottleneck-analysis`.
   This is the one that separates artifact from real: instructions down +
   uOps down + RThroughput flat means take it; instructions down but uOps
   *up* is a real regression (that is exactly how `erfinv`'s poly fold was
   correctly rejected while `acosh`'s select merge was correctly adopted).
4. Wall-clock last, and it is unusable on this machine.

Region-trimmed llvm-mca is ~200x faster than the full harness and exact:
strip all but the `BEGIN`/`END` pairs you care about into a small `.s`.
Never run `examples/mca` concurrently with another `--emit=asm` build --
they write the same file.

### Is this accuracy change real?

**Diff the function's asm region before interpreting any delta.** If the
region did not change, the delta is not real. `cos2pi`'s quick-fuzz avg
moved +59% and its max 64x on byte-identical code; that one check killed
three of four false alarms in a single command.

**Quick-fuzz max is systematically optimistic, not merely noisy** -- it
reported `sin`'s true max of 58 as 4. Any max claim needs `thorough`. Avg
is the trustworthy quick-mode signal and its noise floor is ~+-0.0001.
Never descend on a subsample of a domain you can enumerate: the descent
finds precisely the points you skipped.

### Where do wins usually come from?

The two shapes with the best hit rate, both in the graveyard with worked
examples:

- **The guard is the licence.** A branch whose guard proves a tight bound
  does not need the general kernel. Sweep for it by **reachability** (a
  composite calling a composite calling the kernel -- `srgb` reaches
  `log_2` two levels down through `powf_pos`), by **shape** rather than
  name (`tanh` carries a *standalone copy* of a poly, invisible to a grep
  for the macro), and **downstream** (`compound`'s `log1p` zero-select is
  unobservable because `exp_checked` maps both zeros to `1.0`). Newest
  variant, and the one to try next elsewhere: **the guard can be
  *retuned* to create the licence.** `erfc`/`erfcx` already clamped `|x|`
  before squaring; moving those clamps (11 -> 10.21 and 9.41) made the
  exponent provably in-range, so `exp_checked`'s own clamp could be
  dropped entirely -- bit-identical output, -11.7% throughput on `erfc`.
  Any function that clamps an argument *and* calls a checked primitive
  downstream is a candidate.
- **Delete a rounding, not just an op.** Folding a poly's top coefficient
  group one level lower so `x^4` is never formed removes a multiply *and* a
  rounding -- it made 30 regions smaller while *improving* accuracy on
  every function it touched, including `tanh`'s max 6 -> 5 that three
  dedicated refit attempts could not reach.

- **Wide tier, after the two gather-free passes (graveyard 2026-09-26).**
  Remaining levers, priced: (a) AVX2 is ~2.4x AVX-512 per element for 2x the
  lanes -- each funnel shift is 4 ops (no `vpshldvd`) and ~20 constants live
  on the stack; LLVM canonicalises every source-level trick tried so far
  back to the same code; (b) scalar latency is still ~+3 ns over the old
  table version (window select + int->float); (c) on |x| >= 1 both tiers
  average ~0.25 ulp (the exhaustive 0.13 is diluted by the |x| < 1 half);
  that is `r`'s f32 rounding plus the poly, and an `r_lo` correction would
  cost ~3 ops. Measure with `examples/trig_bench.rs`, not
  quickbench, whose band mask constant-folds the window.
- **mca vs reality, measured 2026-08 (perf-stat cycle counts,
  i5-1145G7):** attribution correction to the entry below: the wide
  tier's ~3x understatement is **the gathers, not the f64 arithmetic**.
  Every pure-f64 hot path matches its mca region (real/simulated:
  powf 0.97 @ 6.8 cyc/elem vs std libm 24.6; compound_accurate 0.93;
  remainder_wide 1.01; logaddexp_accurate 1.00 @ 7.6; sin_checked 1.24
  @ 3.1 vs plain sin 1.6). Zero scalar-f64 instrs in any throughput
  region -- all f64 paths vectorize packed at VF=8. Only sin_wide sits
  far off (4.5x), and a gather-only kernel prices the three
  `vpgatherdd` per 8 lanes at **~10.5 of its ~19.6 cyc/elem** (mca books
  them at 1.5): the tier is gather-throughput/latency bound, which is
  why deleting f64 pipe ops moves Block RThroughput but not real cycles
  for sin/cos while tan_wide (2x poly + divide competing with the
  gathers) reliably gains. Landed from this: the chain's two mul+add
  pairs fused to two exact-product `mul_add`s in `reduce_pi_wide` --
  bit-identical rounding sites, tan_wide -3..5% across five interleaved
  old/new rounds, sin/cos parity. Screened and rejected the same
  session: a shared-scale table re-slice (uniform 2^-28 units so one
  `mm` serves all planes) is arithmetically impossible with u32 chunks
  -- deep-bit chunks are large integers whose smallness lives only in
  the per-plane power-of-two scale, so dropping that scale makes planes
  1/2 contribute 2^29/2^58x too much (caught by exhaustive sweep: avg
  7e8 ulp); selecting the sub-cut bypass residual at the f32 level plus
  an explicit NaN tail select measured +2% on sin_wide despite -12%
  instrs, and was reverted; magic-add rounding replacing reduce_pi64's
  two vrndscale peels measured dead even (±0.7%) and was reverted.
  Measuring technique: wall-clock ns drifts ±8% run-to-run on this
  machine and cannot resolve <10% deltas -- use `perf stat -e cycles`
  on a monomorphized tight loop (examples/quickbench-style, black_box
  the second argument ONCE outside the loop, not per element: a
  per-element black_box inflated powf 9x), interleaved old/new rounds
  via jm stash/unstash. Hard rules distilled from this are in CLAUDE.md
  "f64 in hot paths". Screened and rejected the same session: a shared-scale table
  re-slice (uniform 2^-28 units so one `mm` serves all planes) is
  arithmetically impossible with u32 chunks -- deep-bit chunks are large
  integers whose smallness lives only in the per-plane power-of-two
  scale, so dropping that scale makes planes 1/2 contribute 2^29/2^58x
  too much (caught by exhaustive sweep: avg 7e8 ulp); selecting the
  sub-cut bypass residual at the f32 level plus an explicit NaN tail
  select measured +2% on sin_wide despite -12% instrs, and was reverted.

## Open: fit and accuracy search

These would improve an existing function at zero perf cost, so they clear
the bar by construction if they find anything.

- **Peel the scaling constant out of `tand`, the way `tanpi` now does.**
  `tanpi` went max 5 -> 2, avg 7x better, and an instruction *cheaper*, by
  fitting the polynomial to `tan(pi*w)/w - fl(pi)` and handing the leading
  `pi*w` to a closing `fma` instead of forming `t = fl(PI*w)` first. The
  win is not the deleted rounding -- a prior screen priced that at 0.44
  ulp and correctly rejected it -- it is that the polynomial then carries
  at most 0.215 of the result instead of all of it. `tand` (max 3) ships
  as `sind(x)/cosd(x)`, and the graveyard's direct-poly attempt for it
  failed on its *reduction* (it reused `sind`'s `d` as the pole distance),
  not on this. Same shape, with `K = fl(pi/180)` and an exact `45-|e|`
  pole distance in degrees. `sind` domain, not `tanpi`'s.
- **`log_2_normal`'s two-word `LOG2_E`.** The `log10_normal` half of this
  is **measured and closed** -- see graveyard.md. Summary of what it costs
  whoever picks up the `log_2` half: the low word is *not* a drop-in,
  because `Q` is fitted against whichever leading constant the code uses
  and so already absorbs everything about the one-word value that a
  polynomial can (adding `LOG10_E_LO` alone measured **worse** than
  shipping neither). Refitted, and only at degree 8, it is worth
  `log10p1` max 2 -> 1 and avg -15/-21% -- for +2 fma and a uniform
  +11-12% on every mca rung, which is why it was rejected. `log_2` has
  strictly less to gain (`LOG2_E`'s offset is 1.33e-8 against `LOG10_E`'s
  2.33e-8, and `log2`/`log2_unchecked` are already 0.003/1 and 0.006/1),
  so expect the same cost against a smaller win. Do the oracle screen
  first: for `log10` the correctly-rounded-mantissa-term oracle capped the
  whole lever at -19% aggregate avg before any code was written.

- **`probit`'s `SQRT_2 * erfinv(...)`** is the last un-split single-word
  irrational multiply in the crate. `norm_cdf`'s `xa * FRAC_1_SQRT_2`, which
  used to be listed beside it, is **measured and closed** -- see
  graveyard.md: removing that argument's rounding *entirely* (two-word `z`,
  residual carried into `erfcx_pos`) is worth avg 0.0657 -> 0.0656 and max
  7 -> 7 for +10.9% throughput, because `erfcx`'s `d(ln erfcx)/d(ln z)` is
  only -0.73 at the binding point. `probit` is the opposite case --
  `erfinv` is steep near the ends -- so the result does not transfer, but
  price it by the callee's condition number at the binding point before
  writing code, and note it still needs a harness row first (`accuracy.rs`
  scores `probit` by a round-trip residual, not in ulp).

- **Sollya `fpminimax`**: not installed. Candidate: `acos_poly` (max 4).
  (erfc's n/d used to be listed here as "max ~100, root cause already
  traced — a tighter fit alone won't fix it". Both halves of that were
  wrong: the rational's *own* fit error was ~15 ulp, and replacing it
  with `erfcx_pos` took erfc 109 → 7. See graveyard.md.)

- **True rlibm-style exhaustive rounding-interval LP** (not a continuous
  minimax stand-in) for `log_2`: input-multiplicity confirmed tractable
  (2^23 distinct reduced values). The continuous-LP variant already tried
  found no benefit (log_2's error is rounding-chain-dominated, not
  fit-dominated) — this discrete formulation is the one remaining
  unexplored lever. `exp2`'s own reduction is infeasible regardless (2^29
  distinct values). **Much less attractive now**: the rounding chain was
  the diagnosis and it has been fixed (max 3 -> 1, see graveyard #202's
  peel), so what is left is 1 ulp against std's own 1, i.e. at most a
  correctly-rounded-vs-faithful distinction.

- **Degree 8 for `ln_normal`'s peeled `Q` is a real, measured, unshipped
  Pareto point.** Aggregate avg 0.006315 vs the shipped degree 7's
  0.007262, and 0.0598 vs 0.2161 over the `k == 0` octave, for +3
  instructions and Block RThroughput 16 -> 17 on `ln_unchecked` (21 -> 22
  `log1p`, 42 -> 44 `asinh`, 40 -> 42 `acosh`). Rejected because on the
  real public functions it is worth only 3-14% of the average and no max
  at all. Coefficients are in graveyard.md if a caller ever wants it.

- **Exponent-splitting the other `log_2`-then-`exp2` composites.**
  `rootn` was fixed (45 -> 1 max ulp, 10.1 -> 0 avg) not by the
  double-float log this entry used to propose but by never forming
  `log2(|x|)` as a single `f32` at all: split `|x| = m * 2^e` first,
  keep `e` in `i32`, and hand the integer part of `e/n` straight to
  `exp2_kf`'s exponent field. That costs no precision anywhere and is
  far cheaper than a `Df32` division. What makes it work is that the
  outer exponent is `1/n` and integer division is exact, so the split
  is exact; `srgb_to_linear`/`linear_to_srgb`'s `2.4` and `1/2.4` are
  not (`2.4*e` is not an integer), so the same trick would need a
  two-term split of `2.4*e` and is not obviously free -- and their
  13/7 max ulp is already argued down to a bounded round trip on a
  `[0,1]` domain (see graveyard.md). `powf` proper is already on
  `log2_df` and does not want this.

- **Ulp-staircase-aware LP grids**: densify fit grids near output
  power-of-2 boundaries where ulp weight steps 2x. Complements the
  ulp-weighted-fit idea above (that's weights; this is node placement).

2. **MIP quantized fit**: HiGHS supports mixed-integer — fit with the
   f32 quantization of each coefficient as integer variables (a
   homegrown fpminimax, since sollya isn't installed).
   **Explicitly NOT closed by the headroom screen, and this is the one
   place that screen has a blind spot.** The screen's `LP(f32)` column is
   *the continuous LP optimum with each coefficient then rounded to f32
   independently*. #2 and #102 search the **joint** rounding of the whole
   set, which is a different problem. The gap is visible in the screen's
   own numbers: for `acos_poly`, LP in f64 reaches 1.49 and independent
   rounding to f32 gives back 1.54; for `acospi_poly`, 1.61 -> 1.62. That
   quantization loss is exactly what a MIP or lattice search targets, and
   nothing in this repo has ever looked for it. Every *other* fit-search
   idea in this section is a better search for continuous headroom, which
   the oracle screen can rule out in one run — these two are not.

3. **1-D exhaustive ±few-hundred-ulp scan of each poly's final combine
   constant** scored on the real fuzz — the cheap slice of #1.

6. **Joint threshold+coefficient coordinate descent** in tune.rs
   (crossover as a continuous search parameter) — automates the asin
   fix-5 lesson instead of retuning thresholds against frozen polys.

10. **exp2_q_poly combine-sensitivity LP** (weight by the
    `fma(q, t1*f, t1)` combine's local derivative) — first check where
    the real worst f sits: the technique's documented failure mode is
    the sensitivity vanishing exactly at the hard region.

36. **rlibm-style discrete rounding-interval LP extended to ln/log10**
    (same 2^23 reduced-input multiplicity as the existing log_2 entry).

102. **LLL/lattice reduction over the coefficient quantization step**:
     finds good *simultaneous* f32 roundings of a whole coefficient set
     — the cheap cousin of the MIP idea (#2), and open for the same
     reason: see #2 for why the headroom screen does not cover this.
     Dimension is small (6-10 coefficients), so a hand-rolled LLL is
     feasible without `fpylll`. Best targets are the polys the screen
     scored at ratio <= 1 *because they were already real-chain tuned* —
     there the continuous fit is optimal and the only thing left to win
     **is** the quantization.

103. **Low-discrepancy (Sobol) fit/verification grids** — avoids uniform
     grids' aliasing against ulp staircases; complements the
     staircase-node-placement entry (weights vs placement vs sequence).

104. **True Remez exchange in-repo** (f64, certified equioscillation) —
     the LP is a discretized stand-in; Remez gives certificates and
     better conditioning on rationals. `atan_poly` is the remaining
     candidate; `erfc_rational` used to be the other and no longer
     exists (see graveyard.md).

106. **Caller-profile-weighted alternate coefficient sets** behind a
     cargo feature (e.g. sin weighted toward [−2π,2π]) — same shapes
     and cost, different literals.

194. **Multi-start coordinate descent in tune.rs** (N ulp-perturbed
     seeds around the LP solution, keep best) — cheap robustness
     against the single-seed local-optimum traps already documented.

- **Standing ulp-weighted minimax fit infrastructure**: weight coefficient
  fits by 1/ulp(f(x)) instead of plain relative error as a reusable,
  built-in tool rather than a one-off per-function LP script (the ad hoc
  version of this has already found real wins for exp_pos_neg/erf_poly and
  real regressions for asin_poly and for erfc's retired rational — see
  graveyard.md for when it does/doesn't transfer). Any such tool should
  take the *combine's* target function, not the mathematical one: the
  `asin_poly` regression above is now traced to fitting
  `acos(a)/sqrt(1-a)` where the chain actually wants
  `(fl(pi/2) - asin(a))/sqrt(1-a)`, a `4.371e-8` offset that is 1.47 ulp
  of the result at the branch edge. A single-word constant folded into a
  `fma`'s addend is exactly where this hides.

- **Per-function transformed-variable fit search**: fit in u=s/(s+2),
  u=s·(s+a), etc., searching over the transform family — distinct from
  centered-variable refits (already rejected, that only moved the origin).
  A nonlinear transform changes curvature matching. **Now has one shipped
  precedent, so this is no longer speculative**: `erfcx_pos`'s
  `v = 1/(2+x)` took erfc 109 → 7 and erfcx 126 → 6 where every
  same-variable refit of the old rational had failed. The transferable
  part is *why* it worked — a reciprocal variable maps the whole
  half-line into a finite interval, and an explicit leading `v` factor
  makes the function's own asymptote fall out of `P(0)` by construction
  instead of having to be fitted. Any function with a `c/x` tail
  (`erfcx`, `dawson`, `erfc_inv`'s tail, Mills-ratio shapes) is a
  candidate; a function with no asymptote to reproduce is not.

99. **tgamma** companion to the lgamma entry (Lanczos/Stirling, shares
    machinery).

182. **Compensated-Estrin generic infra** (EFT-based poly evaluation) as
     reusable machinery for future _accurate tiers — compensated-Horner
     was hand-rolled once (erfc's retired rational, max stayed flat);
     infra makes the next attempt nearly free to run. **Weaker now**:
     the crate no longer has any `_accurate` tier for it to serve
     (`erfc_accurate`/`erfcx_accurate` were retired when the base
     functions reached single digits), so this needs a caller before it
     needs machinery.

183. **Full Df32/Df32 division primitive** (div_to_f32 exists) — needed
     by future rational _accurate tiers. Same caveat as #182: no
     rational `_accurate` tier remains in the crate.

188. **_approx tier new members**: exp2_approx/log2_approx/rsqrt_approx
     now have real doc-comment error bounds (max relative/absolute
     error, see lib.rs/git log) — still open: add sin_approx/tanh_approx
     members for ML-inference users, explicitly outside the 0.5/2 budget
     (`cbrt_fast`'s tier, done properly -- `cbrt_throughput` was deleted
     as pareto-dominated, see graveyard.md). sigmoid_approx folds into
     #191's own PWL+correction design instead of a separate bit-trick.
     **Those three bounds are now verified, and the tier's remaining three
     members documented, 2026-07-27** (`examples/approx_bounds.rs`; ulp is
     the wrong metric for this tier, so it lives outside accuracy.rs and
     asserts each doc's own claimed relative/absolute figure).
     - All three documented bounds **hold**, and tightly enough to be worth
       trusting: `exp2_approx` measures 6.149% against a documented ~6.1%,
       `log2_approx` 0.086100 absolute against ~0.086, `rsqrt_approx`
       4.8419% against ~4.8%. Exits nonzero if any drifts.
     - The other three members had **no doc comment at all**; measured and
       documented now. `sqrt_approx`: ~4.5% relative, and uniquely that
       holds over the *whole* positive-normal range (halving the exponent
       field via `>> 1` cannot leave it). `rcp_approx`: ~6.6% but only over
       `[1e-30, 1e30]` — near `f32::MAX` the true reciprocal is denormal
       and it has no range handling, so relative error is unbounded
       (~1e80 at x=3.29e38). `cbrt_approx`: **~5e-6**, far tighter than any
       sibling because it actually refines (two rational steps) — but it
       degrades to **100%** error at the smallest normal, where the seed's
       exponent-field addition leaves the normal range with no rescale.
     - Generalisable observation for adding future members: the two
       *unrefined* single-bit-trick seeds are uniform-ish across the whole
       range (4.5%/6.6%), while the *refined* one is 4 orders of magnitude
       better mid-range but has a hard domain edge. So "add a Newton step"
       buys accuracy and buys a domain restriction at the same time; any
       new `sin_approx`/`tanh_approx` should state both numbers, not just
       the good one.

- **Intermediate sin/cos tier (|x|≲1e5)**: single extra correction word
  over the fast tier's 4-fma Cody-Waite, well short of checked's full
  double-float q — a third point on the speed/domain curve if any user
  workload actually sits there. **Now has measured numbers** from idea
  #47's implementation (2026-07-27, see its entry below): the two_prod
  `PI_HI`/`PI_LO`/`PI_TINY` reduction hits `sin_checked`'s exact
  `|x|<=1e6` accuracy (avg/max 0.0356/2) at 1.278 cyc/elem vs
  `sin_checked`'s 4.546 — **~3.6x cheaper for the same accuracy** over
  that range, and ~11% dearer than fast `sin`'s 1.151. (The ratio was
  4.2x when this was written; `sin_checked` has since come down from
  5.311, so re-measure both sides before quoting it again.) Correctly rejected
  as a fast-tier *replacement*; as a *new* `sin_mid`/`cos_mid` it's
  additive and zero-risk. Caveat that bounds the value: it does not move
  the ~1.3e7 cliff (same `q` failure as fast sin), so the tier's domain
  is fast sin's, not wider — the sell is 7 -> 2 max ulp inside it, not
  more range.

- **Bit-sliced two-for-one sincos**: evaluate sin and cos polynomials
  sharing y=r² registers across the same vector when the caller wants
  both — a sincos slice API where lane pairing amortizes the reduction.
  Only viable inside a slice tier (scalar fusion attempts already failed,
  see rejected section).

### The f64 lever: read this before writing any f64 kernel

Five functions have now moved an f32 double-float/EFT chain into f64
(`sin_checked`/`cos_checked`, `remainder_wide`, `powf`,
`compound_accurate` -- all in graveyard.md with numbers). Three facts about
the lever itself, none of them function-specific:

1. **f64 buys no lane throughput on this machine. It costs exactly 2x per
   lane.** `vfmadd213pd %zmm` is Block RThroughput **1.0**;
   `vfmadd213ps %ymm` is **0.5**; both process 8 lanes. So an f64 port wins
   *only* by shortening the algorithm, and the screen is **what fraction of
   the f32 version is bookkeeping rather than arithmetic**. Halve the
   instruction count and it is a rout (`remainder_wide` -70%); cut 10-18%
   and it is a wash on throughput and a small loss on latency (`powf`, whose
   Block RThroughput went *up* 17% even as its measured cycles went down).
   Know which case you are in before writing code, not after.
2. **`vdivpd %zmm` is 16.0 Block RThroughput.** A division stays the wrong
   answer in f64 too -- a seed plus one Newton step is ~6 ops.
3. **f64 constants cause real register pressure.** One `powf` draft had 23
   of them and llvm-mca showed 24 `vbroadcastsd` *inside* the unrolled loop
   body: LLVM ran out of ZMM registers and rematerialized. Trim every
   polynomial to the accuracy actually required.

**Remaining EFT sites, re-screened with measurements** (an earlier version
of this entry asserted an op-count ratio for both from eyeballing, and was
wrong about one of them -- see graveyard.md):

- `cbrt_accurate`: **fails**, measured. Its `Df32::from_mul(y,y)` / `y2*y`
  / residual block is 7 instructions of a 77-instruction region (9%), so
  2x on the other 70 swamps it. And there is no accuracy lever either --
  `cbrt_accurate` already scores avg 0.000 / max 1, i.e. correctly
  rounded. Do not re-run this one.
- `clog`: **open, and the opposite of what this entry used to say.** Its
  near-1 branch builds `v = re^2+im^2-1` from two `fma` splits and a
  `two_sum`, but the correction word `(e1+e2)+es` is itself summed *in
  f32*, and those roundings sit at the `e`-terms' own ~`2^-24` scale
  rather than at `v`'s. Probed against an exact rational reference over
  the `|z|=1` manifold: forming `v` in plain f64 (`(re*re+im*im)-1.0`,
  four ops, no EFT at all) is **16x more accurate** than the shipped
  ten-op version at `|v| ~ 1e-9` (relative 5.4e-7 vs 8.7e-6) and exact
  where the shipped one reads 1.1e-3. So this is a *simplification and an
  accuracy win at once*, not a screen failure. What still has to be
  checked is cost: `clog` has no mca region (its cost is its
  constituents -- `cabs` 1.178, `ln` 1.611, `carg` 1.694 cyc/elem), so a
  before/after needs one adding, and only the `v` block should move to
  f64 -- porting `ln`/`atan2` too is a different and much larger project.

- **`acosd` wants the two-word `RAD_TO_DEG` now, and its doc comment says
  the opposite.** Measured, exhaustive over every f32 in `[-1,1]` (these
  are in-domain figures; readme's column would be ~0.496x them):

  | `acosd` variant | avg | max |
  |---|---|---|
  | one-word `K`, old `acos` | 0.10975 | 5 |
  | one-word `K`, new `acos` (shipped) | 0.12913 | 4 |
  | two-word `K`, old `acos` | 0.10996 | 5 |
  | **two-word `K`, new `acos`** | **0.05269** | **3** |

  `acosd`'s doc comment records it as "the one member of this family where
  `RAD_TO_DEG_HI`'s two-word `fma` measures no benefit: `acos`'s own error
  dominates". Column 3 shows that was *correct when measured* -- 0.10975 ->
  0.10996, nothing. It is now stale for exactly the reason it gave: after
  the half-angle rewrite `acos` is avg 0.00435 / max 2, so its error no
  longer dominates anything, and what is left is `RAD_TO_DEG_HI` being a
  **truncated** Cody-Waite hi word sitting 5.49e-8 = 0.46 ulp low. That is
  a systematic bias on every result, and removing the noise it was hiding
  under is what exposed it. Textbook [[re-stale-check cross-function
  deps]]: re-verify, do not cite the old entry.

  Note the shipped regression this repairs. The `acos` rewrite moved
  `acosd` **max 5 -> 4 but avg 0.10975 -> 0.12913**, same mechanism -- the
  bias stopped being averaged against `acos`'s own larger, more random
  error. The two-word `fma` takes it to 0.05269 / 3, i.e. 2.1x better than
  before the rewrite rather than merely restoring it. `acosd` is its own
  domain (`acos`'s claim lists it only as "used by"), which is why this is
  a note and not a commit. Cost is one `fma` and one multiply, unpriced.
  `acosd` has no readme precision row; its current exhaustive numbers are
  avg 0.0641 / max 4 over all 2^32 patterns.

- **`exp2m1`'s `saturation_pins` bounds are stale and the gate is not
  doing its job.** `examples/saturation_pins.rs` sweeps `exp2m1` around
  `[-151.0, 128.0]`, but `exp2m1` clamps to `[-126.0, 128.0]`. The lower
  probe therefore sits 25 units *inside* the saturated region and the real
  transition is never tested; the upper one is correct. The fix is the one
  character `-151.0` -> `-126.0`, and it is expected to pass either way --
  the value is `-1.0` on both sides -- which is the point: the gate should
  be exercising the constant that is actually in the code. `exp10m1`'s
  bounds were corrected already; this one was left because it is a
  different domain.

## Open: infrastructure, build and harness

Tooling, feature flags, and portability. None of these change an existing
function's speed or accuracy directly; several unblock ideas above.

93. **2-arg importance-sampling harness** for atan2/hypot/remainder:
    structured lattices near known-hard manifolds — better worst-case
    discovery than uniform sampling. *powf is done*
    (`examples/powfsearch.rs`): deriving `y` from `x` to pin
    `|y·log2(x)|` at the largest finite exponent, then sweeping the
    `k == 0` octave exhaustively, found twice the max ulp the blind
    10M-sample fuzz reported. The pattern to copy is that the manifold
    was read off the *error model* (`ln2·|y·log2(x)|·relerr(log2)`, one
    factor per argument) rather than guessed.

95. **accuracy.rs per-branch attribution mode**: report which select
    arm produced each worst case — speeds every future refit's
    diagnosis step.

96. **Thermal-controlled wall-clock harness** (cpufreq pinning,
    alternating A/B batches): de-risks the one recorded mca-vs-
    wall-clock *direction* disagreement (sincos_checked) and the
    documented 3-9x environment swing.

178. **NEON/aarch64 re-audit**: fma is native there, but every
     mca-derived scheduling decision in this crate is Tiger-Lake-
     specific — the decided tradeoffs (division-vs-poly, Estrin
     groupings) need re-measuring before claiming portability.
179. *Closed, shipped `5db7e4b` -- see graveyard.md.* `erfinv`'s tail
     factors `1-x^2` as `n*(2-n)` on `n = 1-|x|`: 72 -> 11 max ulp, and
     one instruction and one `vdivps` cheaper.

180. *Closed, shipped `efd92b9` -- see graveyard.md, and note the degree
     was **not** the lever.* `erfinv_tail_poly` is fitted in
     `t = sqrt(w) - 1` at degree 9. A degree bump in `w` alone measures
     worse than degree 9 does; the monomial terms reach 15.6x the value
     at degree 10, so the quantisation amplifier outruns the fit.
     15.6 -> 5.9 max ulp on all three callers.

181. **`probit` can drop its `sqrt(2)` multiply -- but only if it drops
     it on *both* arms.** The tail half is built and measured (see
     graveyard.md): `sqrt(2)*sqrt(w)*Q` is `sqrt(2w)*Q` with `2*w` exact,
     a `const SQRT2: bool` on `erfc_inv_half` keeps `erfinv`/`erfc_inv`
     byte-identical, and it was worth -4.1% max / -6.6% avg. Then the
     tail poly's twelfth coefficient landed and moved `probit`'s max to
     the **central** arm (4.98 against the tail's 4.59), where the
     remaining `sqrt(2) * x*P(x^2)` lives -- so the tail fold alone now
     buys no max ulp at all, and it was not shipped. The package that
     still pays is both together: the central multiply folds into a
     scaled, re-descended copy of `erfinv_central_poly` for **zero** extra
     operations (9 duplicated constants), and the two arms together are
     4.98 -> ~4.3. Note the central arm's *other* ulp while you are in
     there, and it is the bigger one: `x = fl(1 - n)` costs 2.66 -> 3.60
     because `n < 0.5` puts `x` on a coarser grid than `n`, and undoing
     that needs a compensated `x` (~3 ops). The old "`erfinv_central_poly`
     itself has no headroom" note is **dead** -- that verdict belonged to
     the pre-peel coefficient grid; see graveyard.md.

182. **`erfc_inv`/`probit`'s remaining average is the argument path, and
     it is now the larger half.** With `erfinv_far_poly_m1` refitted, the
     far branch's own fit is ~0.26 avg ulp against ~0.28 for the `s`/`w`/`v`
     roundings that feed it, measured by oracle rows on a model of the
     chain that reproduces both exhaustive baselines to three decimals.
     Split of that 0.28: `v = fl32(sqrt(w))` costs **0.262**, `w`'s own
     rounding to f32 **0.130**, and `s = fma(-n, n, 2n)` essentially
     nothing. Two priced-but-unmeasured levers:
     - *Make the sqrt last.* `z = sqrt(w * S(v))` with `S = (1+Q)^2`
       fitted directly puts the final rounding on the answer itself, where
       it costs 0 instead of being carried in and rounded twice. The poly
       argument `v` can stay sloppy (only `dQ/dv` sees it). Cost is +1
       `vsqrtps` +1 `vmulps`, and the two sqrts serialise -- so this is a
       throughput/latency question mca has to answer, not an accuracy one.
       Ceiling is ~0.26 of the ~0.35 avg, minus the new `w*S` rounding
       (~0.12), so realistically ~0.35 -> ~0.22.
     - *A poly in `1/sqrt(w)`* fits ~8x tighter (deg 6: 0.257 idealized
       max against 1.64 for the shipped deg 7 in `sqrt(w)`), so it could
       take the fit's 0.26 to ~0.03 -- but it needs a `vdivps`, and it
       only reaches the *fit* half of the budget, i.e. ~0.35 -> ~0.29.
       The sqrt-last lever is strictly the better-value one.

183. **`erfinv_tail_poly_m1` cannot be screened by LP in the monomial
     basis.** It is degree 11 in `t = sqrt(w) - 1` with `t` up to 3, so
     `vander` spans `3^11` and HiGHS reports infeasible well above the
     true minimax value -- the "optimum" it converges to (5.11 result-ulp)
     is *worse* than the shipped coefficients already achieve (0.875), which
     is the tell. Anyone re-screening this poly needs a Chebyshev basis
     (or a rescaled `t`) before the LP means anything. The naive-basis
     front was measured anyway and is a reject; see graveyard.md.

- **`worst_corpus` covers `erfinv` but not `erfc_inv` or `probit`.** Both
  are public, both have exhaustive rows in `accuracy.rs`, and both just
  moved by a change that the corpus gate could not see (0 entries; the
  `erfinv` change moved 12 of its 90). They reach a whole branch
  `erfinv` cannot -- `erfinv_far_poly_m1` is unreachable from `erfinv`,
  so that poly currently has **no bit-exactness gate at all**. Same
  "gate whose coverage is a hand-written list" smell as the
  `denormal_audit` item below.

- **`denormal_audit`'s 2-arg coverage is still one function wide.** Both
  of its hand-maintained lists take `fn(f32) -> f32`. `logaddexp` is now
  in the (A) table via the `logaddexp(x, 0) == softplus(x)` curry, which
  is what `logaddexp_checked` was built against -- but `hypot`, `atan2`,
  `powf`, `remainder`, `xlogy`, `xlog1py` and `compound` are still simply
  absent, and no 2-arg function has a curry that exact, so each needs its
  own choice of pinned argument and its own reference. Same "gate whose
  coverage is a hand-written list" smell that produced every finding in
  that file so far.
- **`silu_checked` still writes out the `exp_neg_scaled64!` reduction
  inline.** The macro exists now (used by `softplus_checked` and
  `logaddexp_checked`) and `silu_checked`'s four lines are textually the
  same, minus the trailing `* 2^-64` it deliberately does not do. Folding
  it in is a one-line change to a domain nobody held at the time; verify
  with the full pre/post asm-region diff, which is what showed the first
  two conversions were byte-identical.

- **The `1 + r` peel is free for the core-locked `exp_r_poly!` too, and
  is screened but not taken.** `exp_reduce!` now evaluates `e^r - 1` and
  folds the `+ 1` into its reconstruction as `fma(s, t1, t1)`; the same
  transform applies verbatim to `exp_r_poly!`. There are **six** call
  sites, not the ~20 first recorded here, and reading all six confirms
  every one terminates in a multiply by an *exact power of two*, so the
  `fma(s, t, t)` absorption applies at each without an awkward caller to
  design around: `exp` and `exp_scaled` (`p * t1 * t2` off
  `exp2_field_split`), `exp_narrow` and `sigmoid` (`p *
  exp2int_field!(k)`), and `exp_neg_scaled64!` and `silu_checked` (`p *
  exp2int_field!(k + 64.0)`). Measured on the shipped f32 evaluation
  order over
  `|r| <= ln2/2`: max 2.392 -> 2.148, avg 0.7165 -> 0.6783, at **zero**
  instruction change (the `r + 1` add pays for the multiply-to-fma
  swap). **That "zero" may be conservative and should not be treated as
  confirmation**: it was measured on the isolated poly, whereas at the
  call sites the peel deletes the `l0 = r + 1.0` `vaddps` *and* turns one
  `vmulps` into a `vfmadd` -- the same -1 arithmetic instruction the
  `exp_reduce!` commit measured. Unverified (it needs the lock and an asm
  diff), so expect a possible small win rather than a wash, and do not
  read a flat instruction count as agreement. The mechanism is that
  `fl(1 + r)`'s rounding enters at
  `(1 + r)/e^r ~ 0.92` of full weight while the peeled one is attenuated
  by `|e^r - 1|/e^r <= 0.415`. Not taken here because it is a *core*
  edit -- it needs `jm core claim` and a whole-crate sweep -- and on its
  own ~0.25 ulp it probably does not move `exp`'s integer max off 3. It
  is worth doing opportunistically the next time somebody holds the core
  lock for another reason. Do **not** pair it with a degree-6 bump for
  those callers: that was measured and correctly rejected at +5-23%
  throughput, and unlike the erfc family they have nothing downstream
  amplifying them.
- **`erfc` and `norm_cdf` are now `erfcx_pos`-bound, not exponential-
  bound.** Both sit at max 6 with the Gaussian half no longer the larger
  term. `erfcx_pos`'s remaining error is its own f32 *evaluation*, not
  its fit (handing it an exactly-rounded `v` does not lower the max --
  see its comment), so the levers already closed are the ones aimed at
  `v`. What has not been tried is a different evaluation order for the
  degree-10 polynomial itself under the current, now-quieter, error
  budget. The Horner endpoint is now closed on accuracy as well:
  re-priced against the current coefficients it is worth ~0.001 ulp
  on max, not the ~0.5 the old note records -- see graveyard.md.

  **But only part of that rejection went stale, and the rest is still
  binding.** Horner's +19-34% latency on four public functions does not
  depend on how quiet the exponential half became, so the endpoint stays
  rejected on its own terms. And a **refit is mandatory, not optional**:
  the coefficients are coordinate-descent polished against the exact
  Estrin order, so reordering without refitting breaks the `erfc(0)` /
  `erfcx(0)` / `norm_cdf(0)` edgecheck pins immediately. What is
  genuinely open is an *intermediate* evaluation order, not the Horner
  endpoint. (Same shape as the two-word-constant rule: a polynomial
  tuned against one arrangement has already absorbed it, so changing
  the arrangement without refitting is not the same experiment.)

## `dawson`'s max 6: the fit is losing 3x to *coefficient quantization*, not to degree

Left open after the constant-peel landed (see graveyard.md). The chain is
no longer the binding term: at both worst inputs (`x = 1.4041231` for the
old form, `x = 1.4212224` for the shipped one) the rational evaluated in
f64 from an exact `u` already scores **3 ulp on its own**. Getting the
reported max under 6 therefore needs the fit under 3, and nothing about
the evaluation order can do it.

The numbers that make this a lead rather than a closed door are already
in the `[6/6]` entry's own table:

| | continuous minimax | after f32 quantization |
|---|---|---|
| `[6/6]` (shipped) | 1.25e-7 (2.09 ulp-eq) | **2.36** |
| `[7/6]` | 6.45e-8 (1.08) | 3.07 |
| `[7/7]` | -- | quantization-dominated |

`[7/6]`'s *continuous* optimum is 1.08 ulp-equivalent -- roughly 3x
better than what `[6/6]` ships -- and the entire gain is destroyed by
rounding the coefficients to f32, which is why the search stopped at
`[6/6]`. But that stop was on **independent** rounding of each
coefficient plus a descent that was only run on `[6/6]`. A quantization
descent over `[7/6]`/`[7/7]`'s larger coefficient space is the one thing
that was never tried, and it is searching for ~3x of headroom that
provably exists in the continuous problem. If it lands under ~2.0
quantized, `dawson`'s max should follow.

Two smaller ones from the same measurements:

- **`u = fl(x*x)` costs one ulp of max in `[2,4)`** (fit alone 2, fit
  plus `u`'s rounding 3) and ~0.02 avg overall. A double-`f32` `u` fixes
  it but the polynomial then has to consume `u_hi + u_lo` on both sides,
  which is not obviously worth one ulp on one octave -- price it only
  after the fit is under 3, since it cannot be the binding term before.
- **Folding the closing `x *` into the numerator's peel**
  (`fma(x*u, A, x) / (1 + u*B)`) is accuracy-equivalent at the same op
  count (central avg 0.5730 vs 0.5686 sampled) and was passed over on the
  documented numerator/denominator-symmetry argument, *without* measuring
  its codegen. If the symmetry rule is ever re-tested, this is a free
  candidate sitting next to it.

## `atan2pi`/`atan2`'s `x > 0` arm is division-limited

After the half-turn fold landed, the `x > 0` half-plane's residual
decomposes (20M samples of that arm) as:

```
new atan2pi total          avg 0.1177  max 3
y/x rounding alone         avg 0.1011  max 1
atanpi(q), q as given      avg 0.0289  max 3
```

**86% of it is the single `y/x` division's own rounding**, which no
rearrangement of the quadrant fold can reach -- the fold work is done.
The `x < 0` arm is at 0.0237 because adding an exact `+-1.0` puts the
result in `[0.5, 1]` while `atanpi(q)` may be tiny, damping the
quotient's relative error rather than exposing it.

A two-word `y/x` is the only remaining lever, and it would apply
identically to `atan2` (0.066) and `atan2d` (0.105), not just `atan2pi`.
Not screened: the cost of a compensated division here is unpriced, and
`compensated division reciprocal overflow` is a known trap on this exact
shape (a standalone `1.0/x` can overflow for legitimate small `x`).

## `log1p_exp_neg_f64`'s two polynomials are both over-degreed for the floor they feed

`logaddexp_accurate` is the crate's most expensive region by Block
RThroughput -- 152 instrs / 214 uOps / **94.00** BlockRT, 3.1x
`logaddexp`'s throughput and 1.6x its latency (table in graveyard.md).
Almost all of that is the two f64 polynomials in `log1p_exp_neg_f64`, and
neither was ever sized against the accuracy the chain actually needs.

Both targets have **exact rational Taylor series**, so the tail can be
bounded without any function evaluation: expand the Taylor poly in the
Chebyshev basis on the working interval, drop the tail, convert back.
50-digit `Decimal`, no fitting tools required. The tails:

| | shipped degree | Chebyshev tail at shipped | degree reaching ~5e-18 |
|---|---|---|---|
| `EXP_F64_P`, `p(f) = (e^f-1)/f` on `[-ln2/2, ln2/2]` | 12 (13 coeffs) | 1.4e-21 rel | **10** (8.79e-18 rel) |
| `ATANH_B64`, `Q(u) = sum u^k/(2k+3)` on `[0, 1/9]` | 14 (15 coeffs) | 7.5e-26 rel | **9** (4.89e-18 rel) |

The chain's own documented floor is `~3e-16` relative, so degree 10 and
degree 9 would contribute **~3% and ~1.6%** of a budget that is already
spent elsewhere. The shipped degrees are ~5 and ~13 orders of magnitude
tighter than anything downstream can observe -- `logaddexp_accurate`
narrows to f32 at the end, and already measures avg 0.000 / max 0 against
the f64 reference. This is a **pure throughput/latency lead with no
accuracy axis to trade**: the win, if it is real, is free.

Dropping 2 + 5 coefficients removes roughly 7 `fma` from a 152-instruction
region (~5%), and shortens both Estrin trees, which is the part that
should matter to the 121.5-cycle latency.

**Screened only this far**: the tail bounds above are reproduced and
correct. What has *not* been done is everything that decides it --
coefficients not generated, no accuracy sweep, no `edgecheck`, no mca. Two
specific reasons it may not convert, both of them this crate's own
repeated lesson:

- A Chebyshev *truncation* bound says nothing about the shorter
  polynomial's **evaluation rounding**, and evaluation roundings, not fit,
  are usually what sets the floor here. The margin is enormous (3e-18
  against 3e-16), so this is unlikely to bind -- but it is the thing to
  measure, not assume.
- Both polys are Estrin-grouped and the grouping is hand-written around
  the current degree (`f2`/`f4`, and `u2`/`u4`/`u8` plus the `1e-30` clamp
  that exists precisely to keep `u^4`/`u^8` normal). A degree change
  rewrites those trees, and the clamp's justification is stated in terms
  of the current tail term -- re-derive it rather than porting it.

The economization script is straightforward enough to rewrite from the
description above; it needs `Decimal` only.

## Denormal intermediates in `x * f(x)` composites: one found, the class not swept

`gelu` carried a normal-result-through-a-denormal-intermediate defect for
as long as it has existed (max 9 -> 6 and -13.7% throughput once
reassociated -- see graveyard.md). `silu`'s was found by audit. Both are
the same shape and neither was found by looking for the shape:

> any `x * f(x)` whose `f` decays exponentially has a window where `f` is
> denormal and `x * f(x)` is not, because the extra `|x|` keeps the
> product representable for `ln|x|`-worth longer.

The window is ~0.4 wide in `x` for `gelu` and ~20 for `silu`, and inside
it the result keeps only the bits the denormal had left. **The screen is
cheap and mechanical**: for each candidate, compare where `f` underflows
against where `x*f(x)` does; if they differ, the product is being formed
in the wrong order. Where a reassociation exists it costs nothing --
`gelu`'s was folding an exact `0.5*|x|` into the Gaussian factor rather
than into the finished product.

Not swept: `norm_pdf`-scaled products, `xlogy`/`xlog1py`, `sinc` near its
own decay, `dawson`'s tail, and anything else multiplying a decaying
kernel by a growing factor. `examples/denormal_audit.rs` already reports
per-function flush fractions, but it audits **outputs**, not
intermediates, so it cannot see this class -- `gelu` read as "flushes 11%"
there and the 11% was a *symptom*, not the defect.

## `gelu` is a pass-through *now*

With the denormal band and the argument repair both gone, `gelu`'s
remaining max 6 sits at `x = -1.32`, inside `erfc`'s own worst region, and
its avg 0.1309 is `erfc`'s 0.1289 plus composition rounding. The earlier
closure said this and was right about the arithmetic while being wrong
about the location; it is now true of the shipped code. What is left moves
when `erfcx_pos`'s polynomial *evaluation* does, and that entry is
unchanged.
