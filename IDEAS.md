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

## Open: fit and accuracy search

These would improve an existing function at zero perf cost, so they clear
the bar by construction if they find anything.

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

- **Peel the leading term out of `ln_normal`/`log10_normal` too**
  (graveyard #202's log_2 lever, transplanted). Both are the same
  `k*c + s*P(s)` shape with the same rounding-chain-dominated 3-ulp max
  that `log_2` had, and both currently evaluate their leading coefficient
  *inside* `P`, so for `x` near 1 all of `P`'s full-weight roundings land
  straight on the answer. `ln` is the more promising of the two by some
  margin: its leading coefficient is exactly `1.0`, so the peeled form is
  `ln(m) = s + s^2*Q(s)` with an **exact** leading term -- strictly
  better than `log_2`'s, which still rounds `s*log2(e)`. Reuse all three
  things that had to be right there: an ulp-weighted LP fit of `Q`
  (weight `s^2/ln(1+s)`, *not* a term dropped off `P`); `k`'s Cody-Waite
  combine joining last so it keeps its single rounding; and the `s^2`
  factor riding into the poly's own low group so the shape stays three
  Estrin levels deep. Budget the +1 dependency level as certain (log_2
  paid +11% latency for it at flat throughput). Big blast radius --
  `ln`/`log10`/`log1p`/`log2p1`/`asinh`/`acosh`/`atanh`/`compound`/
  `logit`/`xlogy`/`xlog1py`/`softplus`/`logaddexp`/`logsigmoid`/`clog`
  all route through these -- so it wants its own exhaustive pass per
  function, which is why it was not done alongside `log_2`.

- **`rootn` is now the family's odd one out**: 43-45 max ulp at `|n|<=3`,
  where `powf` is 3. It is still `exp2_checked(log_2(ax) / n)`, i.e. the
  single-f32 route `powf` just left, and its error is the same
  `y`-amplified one (`|log2(x)/n|` reaches ~63 at n=2). The df route
  fixes it, but not for free and not by simply multiplying by `1/n as
  f32` -- that constant's own 2^-24 error is amplified right back, so it
  needs `log2_df(ax) / (n as f32)` through `Df32`'s real division (a
  hardware `divps` plus ~5 ops). Worth doing only if a caller cares;
  `rootn` gets *more* accurate as `|n|` grows, so the bad region is
  exactly the small `n` a caller could write as `cbrt`/`sqrt` instead.
  `srgb_to_linear`/`linear_to_srgb` (13/7 max ulp) sit on the same
  fence with the same fix available and a much stronger case for leaving
  them alone -- their exponent is a constant `2.4`, so the amplification
  is bounded and their inputs are `[0,1]`.

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
  graveyard.md for when it does/doesn't transfer).

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