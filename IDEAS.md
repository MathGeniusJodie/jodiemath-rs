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
  unobservable because `exp_checked` maps both zeros to `1.0`).
- **Delete a rounding, not just an op.** Folding a poly's top coefficient
  group one level lower so `x^4` is never formed removes a multiply *and* a
  rounding -- it made 30 regions smaller while *improving* accuracy on
  every function it touched, including `tanh`'s max 6 -> 5 that three
  dedicated refit attempts could not reach.

## Open: fit and accuracy search

These would improve an existing function at zero perf cost, so they clear
the bar by construction if they find anything.

- **Sollya `fpminimax`**: not installed. Candidates: `acos_poly` (max 4),
  erfc's n/d (max ~100, root cause already traced — a tighter fit alone
  won't fix it).

- **True rlibm-style exhaustive rounding-interval LP** (not a continuous
  minimax stand-in) for `log_2`: input-multiplicity confirmed tractable
  (2^23 distinct reduced values). The continuous-LP variant already tried
  found no benefit (log_2's error is rounding-chain-dominated, not
  fit-dominated) — this discrete formulation is the one remaining
  unexplored lever. `exp2`'s own reduction is infeasible regardless (2^29
  distinct values).

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
     better conditioning on rationals (atan_poly, erfc_rational).

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
  real regressions for asin_poly/erfc — see rejected section for when it
  does/doesn't transfer).

- **Per-function transformed-variable fit search**: fit in u=s/(s+2),
  u=s·(s+a), etc., searching over the transform family — distinct from
  centered-variable refits (already rejected, that only moved the origin).
  A nonlinear transform changes curvature matching.

99. **tgamma** companion to the lgamma entry (Lanczos/Stirling, shares
    machinery).

182. **Compensated-Estrin generic infra** (EFT-based poly evaluation) as
     reusable machinery for future _accurate tiers — compensated-Horner
     was hand-rolled once (erfc, max stayed flat); infra makes the next
     attempt nearly free to run.

183. **Full Df32/Df32 division primitive** (div_to_f32 exists) — needed
     by future rational _accurate tiers.

188. **_approx tier new members**: exp2_approx/log2_approx/rsqrt_approx
     now have real doc-comment error bounds (max relative/absolute
     error, see lib.rs/git log) — still open: add sin_approx/tanh_approx
     members for ML-inference users, explicitly outside the 0.5/2 budget
     (cbrt_throughput's tier, done properly). sigmoid_approx folds into
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
  `sin_checked`'s 5.311 — **~4.2x cheaper for the same accuracy** over
  that range, and ~11% dearer than fast `sin`'s 1.151. Correctly rejected
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

## Open: infrastructure, build and harness

Tooling, feature flags, and portability. None of these change an existing
function's speed or accuracy directly; several unblock ideas above.

93. **2-arg importance-sampling harness** for powf/atan2/hypot/
    remainder: structured lattices near known-hard manifolds (e.g.
    y·log2(x) near integers) — better worst-case discovery than
    uniform sampling.

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