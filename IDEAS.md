# Ideas tried and rejected

Budget: ≤0.5 avg ulp, ≤2 max ulp. Ideas tried and reverted/not adopted, so
they aren't re-attempted without new information. Adopted changes live in
git history / readme.md, not here. Untested backlog is at the bottom.

## Cross-cutting

- **Degree-reduction probes (2026-07-07)**: lolremez rejected all 3 —
  `log_2` 9→8 (max ulp 3-5 vs 2 cap), `exp2` Q 5→4 (est. max rel error 40x
  worse), `sinf_poly` 9→7 (178x worse). Not implemented.

- **exp2/log_2 coefficient refit (2026-07-07)**: coordinate-descent tuner
  found bit-identical (zero-move) coefficients on both — already optimal
  from an earlier session.

- **exp2's Q(f) poly, max-capped LP refit (2026-07-09)**: isolated fit
  predicted a clean win (avg weighted error 0.043→0.013), but regressed
  `exp2`/`exp2_checked`/`exp10`/`exp10_checked`'s max ulp 1→2 across the
  board (exhaustive-confirmed). The poly's shipped max ulp was already 1 —
  no real margin for the LP's continuous model to protect an actual ulp
  boundary. Reverted. *Don't try this LP technique on a poly already at/near
  max ulp 1 — no margin for the isolated model's own slop.*

- **Tuner "zero-move" trap (2026-07-08)**: `tune.rs` moves coefficients by
  integer bit-steps; a coefficient seeded at `0.0` can only reach
  denormal-scale perturbations, so it's structurally stuck, not evidence of
  no headroom. Confirmed on `atan_poly`'s degree bump: zero-seed reported
  max 18 (unchanged); a real scipy least-squares seed found max 3 on the
  same grid. An arbitrary nonzero seed (`1e-3`) isn't a fix either (never
  recovered, converged to max 565). *Always seed a new coefficient with a
  real fit, not 0.0 or an arbitrary constant.*

## exp2 / exp / sin / cos / tan

- **Fused sincos / direct tan via shared reduction (2026-07-07)**: two
  attempts, both regressed catastrophically near tan's poles (avg ulp
  0.33→1.05/max 3000→32M, worse for the cofunction-identity variant) —
  the additive combine forms cancel badly near r=π/2. Reverted; needs
  multi-word reduction to fix properly.

- **sinf_poly quantized refit (2026-07-07)**: tuned against f64::sin, max
  ulp unchanged, avg barely moved — already near f32's precision floor.

- **sinf_poly LP refit (2026-07-09)**: two variants (plain minimax,
  L1-with-max-cap), both regressed `cos_checked` for real (avg ulp
  0.081→0.79 and 0.0495→0.0542) despite looking better isolated. Root
  cause: `cos_checked`'s reduction lands `r` near the poly's domain *edge*
  for small `x`; `sin_checked` lands `r` near the *center* — the two
  callers stress opposite ends of the shared poly, and a domain-uniform LP
  grid doesn't account for that. Reverted. *For a poly shared across
  callers with different domain-region emphasis, check all callers/buckets,
  not just the isolated metric.*

- **expm1 Pade degree bump 3→5 (2026-07-08)**: scipy-seeded refit found
  real avg-ulp headroom in the near-zero branch (0.138→0.135), but the
  function's actual max ulp (6) lives in expm1's other branch (`exp(x)-1`),
  untouched. Cost real throughput (+7%). Reverted.

- **tanh direct rational P(x²)/Q(x²) over [0,~9] (2026-07-08)**: needs 13
  free coefficients to converge over the full domain — far more than any
  poly in the crate. A 2-domain split needs 14 total, likely more work than
  the current expm1-based formula. Not implemented.

## cbrt family

- **Seed constant + degree-2 poly joint search (2026-07-07)**: best across
  41 seeds still ~50x over budget (max ulp 112). 3 coefficients can't
  correct this seed's error regardless of seed choice.

- **Integer-division-free seed (2026-07-07)**: cheaper codegen (3 vs 5
  instructions), but too coarse for the correction poly to compensate —
  max ulp 33, ~16x over budget.

## log_2 / ln / log10

- **Integer koff fold (2026-07-08)**: bit-exact, but mca showed zero
  measurable change — LLVM already performs this reordering. Reverted to
  avoid an f32→i32 public signature change for no benefit.

- **ln_normal/log10_normal: fuse trailing `+k*LN2_LO` into the fma
  (2026-07-08)**: "one fewer rounding" in isolation, but mca latency got
  deterministically *worse* by 1 cycle (reproducible); accuracy unchanged
  (the term was already off the critical path). Reverted.

- **log1p small-|x| dedicated branch (2026-07-08 + follow-up 2026-07-09)**:
  adds a whole extra poly eval every call (branchless convention evaluates
  every branch unconditionally). Screened as marginal (max 3 vs 4, avg
  0.068 vs 0.073, unconfirmed at full density); implemented and measured
  for real: mca throughput **+48.1%** for that marginal gain. Rejected on
  mca alone, before checking accuracy.

- **Direct minimax refits for ln/log10 (2026-07-08)**: coordinate-descended
  against ln/log10 directly instead of log_2's rescaled coefficients —
  zero-move local optimum, no headroom.

- **ln_normal's poly LP refit (2026-07-09)**: isolated fit predicted a 75%
  avg improvement (largest of the session) but real exhaustive result was
  ~0.09% — noise. Reverted. *Isolated LP predictions don't reliably predict
  real magnitude, regardless of how large the prediction is.*

## hypot / misc

- **Compensated hypot (2026-07-08)**: `e=fma(r,-r,s); r+e/(2r)` — no
  measurable accuracy improvement (already near correctly-rounded) but real
  cost: latency +109%, throughput +87%. Reverted before edgecheck.

## asin / acos / atan / atan2

- **acos_poly Horner→Estrin (2026-07-07)**: real latency win, but fma
  reassociation regressed asin max ulp 9→12, acos 4→5 (retuning made it
  worse, →6). Reverted — acos's accuracy is a protected invariant.

- **erfc's n/d rational Horner→Estrin (2026-07-07)**: small theoretical
  win, measured as a wash on speed plus a real accuracy cost (avg +2.7%).

- **erf's near-zero Padé branch refit (2026-07-07)**: max ulp unchanged,
  avg moved <0.3%. No headroom.

- **Same branch, LP numerator refit (2026-07-09)**: isolated fit predicted
  an 84% avg improvement, but real fuzz found a *regression* (0.317→0.325),
  worst-case landing right at the 0.28 branch crossover with `erf_poly`.
  Reverted. Second confirmed case (after sinf_poly) of the isolated metric
  getting the *direction* wrong — both involved a branch-crossover
  boundary. *Check the crossover neighborhood for any poly next to a
  domain split.*

- **erf's tail branch (`erf_poly`) refit (2026-07-07)**: max ulp unchanged,
  avg moved <0.3%.

- **acos_poly unconstrained joint acos+asin objective (2026-07-07)**:
  improved joint score but regressed acos's own max ulp 4→5. Rejected in
  favor of a constrained variant.

- **acos_poly degree 6→7 (2026-07-08)**: zero-seeded 8th coefficient
  converged to exactly 0.0 — useless. Re-seeded with a real scipy fit:
  found real headroom (acos avg/max 0.496/4→0.490/3) but asin was unmoved
  and mca cost was real (+7-11% both functions) — not worth it at this
  magnitude. Reverted.

- **acos_poly Df32 leading-term split, pi/2 hi+lo (2026-07-08)**: dramatic
  acos avg improvement (0.496→0.068) but acos max ulp regressed (4→5) and
  asin got worse on both axes (max 9→12). Retuning recovered asin but
  acos max still regressed (4→6), plus a real mca cost. Not adopted.

- **acos_poly joint LP, both acos+asin max ulp capped (2026-07-09)**, two
  attempts: attempt 1 let the pi/2 constant drift enough to land on a worse
  f32 value at `a=0`, corrupting near-zero calls (acos avg 0.496→1.905).
  Attempt 2 forced the constant exact, re-solved the rest — still
  regressed for real: asin max ulp 9→12, the *identical* regression
  signature as the Df32-split entry above (a real fragility at this
  crossover, not a fluke). Both reverted. *acos_poly/asin's sharing
  (subset-domain, both callers directly exposed) is much harder to model
  via a domain-uniform LP grid than cbrt/cbrt_accurate's sharing (where
  Newton's step makes the second caller insensitive to seed error).*

- **atan2 division-residual correction (2026-07-07, re-tested 2026-07-08)**:
  looked like a 10x win until a NaN-sentinel bug was fixed — evaporated to
  noise-level (atan's own poly error dominates). Real cost (latency +14%,
  throughput +43%) for zero benefit, both before and after atan's own
  accuracy later improved. Closed for good.

- **exp: k1/k2 split → k-clamp (2026-07-08)**: the "only matters
  out-of-contract" premise was wrong — k=128 is reachable from ordinary
  in-domain x, and clamping silently drops a factor of 2 for real inputs.
  Killed by a scratch-test before mca.

## Codegen & build hygiene

- **Forcing zmm-width AVX-512 (2026-07-07)**: CPU supports it, but LLVM
  deliberately avoids it (likely downclocking avoidance) — no rustc knob
  to force it independent of target-cpu tuning. Not forced; inconclusive.

## sin_checked / cos_checked internals

- **round_x_over_pi: remove dead pre_offset=0.0 add (2026-07-07)**:
  instructions confirmed gone via asm, but throughput got *worse* —
  removing an op let the scheduler pick worse elsewhere. Reverted.

- **round_x_over_pi: qh → round_ties_even (2026-07-06)**: regressed
  cos_checked's max ulp 2→6 (rare exact-half ties clash with cos's -0.5
  offset). Reverted to `f32::round`.

- **reduce_pi: rebalance 4-deep chain to depth 2 (2026-07-07)**: bit-exact
  but +3 cyc latency both functions. Reverted to the flat chain.

- **reduce_pi: downgrade e2's two_prod to plain multiply (2026-07-07)**:
  regressed sin_checked's in-domain max ulp badly (2→51,054 at |x|≤1e6).
  Reverted before mca.

- **reduce_pi: downgrade e3's two_prod to plain multiply (2026-07-07)**:
  "provably zero" premise held, expected latency win — but mca diverged by
  caller (sin_checked throughput improved, cos_checked's got worse) and the
  off-contract tail got dramatically worse. Reverted.

- **parity() via integer bit-ops (2026-07-07)**: bit-exact, latency
  unchanged, throughput worse for both — FP-port ops beat integer ops once
  scheduled. Reverted.

- **sin_checked/cos_checked clamp: move bound into poly's y.min()
  (2026-07-07)**: saves one op, but the raw residual `x` still enters
  unclamped, reintroducing the "inf for finite input" bug the clamp
  prevents. Not adopted.

- **Fast sin/cos: shorten PI_A..D chain 4→3-deep (2026-07-07)**: error
  estimate was wrong by orders of magnitude — near sin's zeros the
  single-shot rounding becomes a huge relative error (avg ulp 0.06→1.48,
  max 220→866M). Reverted before mca.

- **cos's q fold into one constant**: dies on representability (folded
  constant doesn't exist as f32); half-magnitude magic quantizes q wrong;
  folding elsewhere reintroduces rounding near cos's zeros. Not attempted.

## Missed fma contractions

- **asin: a*a-a → fma(a,a,-a) (2026-07-07)**: bit-identical, but throughput
  worse (latency unchanged). Reverted.

- **Small-poly Estrin audit, asin_small/sinh_small (2026-07-08)**: not
  bit-exact, measured backwards on every axis for both functions (sinh:
  worse latency+throughput+accuracy; asin: max ulp regressed 9→10).
  Reverted.

## Other spots

- **atan2's bothzero/hpisignx boolean simplification (2026-07-07)**:
  compiled output byte-for-byte identical — LLVM's InstCombine already
  does this.

- **erfc's final fma→mulsign+add (2026-07-07)**: intended codegen shift
  confirmed, but throughput got worse. Reverted.

- **atanh via single log1p call (2026-07-08)**: algebraically simpler, but
  real regression — max ulp 3→31303 (tail-only, near x≈-1). The two-log1p
  form passes exact literals; the one-log1p form's division rounds once,
  and log1p's derivative diverges near the singularity, amplifying that
  tiny rounding error ~5 orders of magnitude. Reverted immediately.

- **erfc: exact exponent via two_prod (2026-07-08)**: real accuracy win
  (avg/max 0.311/109→0.306/93) but real mca cost (+5-7%) for shaving an
  already-far-over-budget max ulp further. Not adopted.

- **erfc in log space, single poly over [0,10] (2026-07-08)**: lolremez
  convergence poor even at degree 15 (16 coefficients). `erf_poly` can't be
  reused directly (catastrophic error, ~297000 avg ulp — only ever fit for
  erf's own objective). Not implemented.

- **powf_checked's ~150 max ulp investigated (2026-07-08)**: backlog's
  "exp2_checked double-rounding into denormals" diagnosis doesn't survive
  measurement — exp2_checked is bit-exact even at the actual worst-case
  input. Residual lives in the log2_df/exp2_checked_df double-float chain;
  not root-caused further. No code changed.

- **ln/log10 fuse trailing fma, re-tested (2026-07-08)**: duplicate of the
  entry above, picked again without cross-referencing, initially
  mismeasured as a win (baseline taken from a stale readme number instead
  of re-measured; an mca delta wrongly dismissed as noise). Re-verified:
  real, reproducible +1-cycle regression. Reverted for real. *Grep this
  file before retrying an idea; never trust a written number over
  re-measuring; repeat-and-confirm any mca delta before calling it noise.*

- **Centered-variable refit for exp2's f (2026-07-08)**: backlog claimed
  centering shrinks coefficients — scipy showed the opposite (centered
  coefficients are *larger*). `exp2`'s R(f) is monotonic with no interior
  minimum, so centering moves away from the edge minimum the current fit
  already exploits. Premise refuted before any Rust written.

- **Tune cbrt_throughput's magic constants (2026-07-08)**: single-octave
  grid looked like a win (max 16→12) but the function's error doesn't
  repeat across octaves like cbrt_normal's — implementing it made the real
  fuzz sweep *worse* (avg 6.73→7.01). Wide 60-octave grid retry exposed
  that `tune()`'s `score()` comparison is max-first (Rust tuple ordering),
  wrecking the average to shave the worst case (avg 6.83→22.73). Not
  adopted. *Check a target's error actually repeats across a partial grid;
  `tune()`'s max-first comparison can hide a real average regression.*

- **log_2: atanh-form reduction t=(m-1)/(m+1) (2026-07-08)**: 1000x tighter
  in the underlying math, but worse for real — accuracy slightly worse
  (log_2 already deep in f32-rounding-noise territory) and latency +43%
  (the one division depends on `s` from the first step, nothing to overlap
  it against, unlike cbrt's early-starting reciprocal). Reverted.

- **erf: joint boundary+coefficients refit (2026-07-08)**: cheap proxy
  sweep (8 threshold candidates against unrefit coefficients) found a flat
  plateau around the current 0.28 boundary — no headroom, matching
  separate fixed-boundary refits. Didn't build the full joint optimizer.

- **"Resurrect the f64-reduction sin/cos" (2026-07-08)**: backlog's claim
  that working code already existed in git history didn't hold up (nothing
  in `git log --all`/reflog). Built one from scratch anyway — failed on
  accuracy by a wide margin (max ulp 637 even in the smallest bucket, vs.
  the claimed ~1e9). An f64 two-product compensation helped (637→70) but
  couldn't close the gap — f64's own PI constant has irreducible ~2^-53
  error no surrounding compensation fixes. Not implemented. *Verify "the
  code already exists" claims via git log before trusting them.*

- **cbrt: joint seed-constant + degree-3 search (2026-07-08)**: coordinate
  descent and an independent Python sweep (±2000 seed offsets) both found
  ~0.3% or less movement — shipped combination already essentially
  optimal.

- **erfc domain split [0,2]/[2,10] (2026-07-08)**: underlying math
  100-380x tighter per-domain; real implementation improved avg ulp a lot
  (0.311→0.189, ~39%) but max ulp barely moved (106→105), same worst-case
  neighborhood. Bottleneck lives downstream of the correction term
  entirely. Not adopted (bar was 109→single digits).

- **Select-tree "LUT" for exp2 (2026-07-08)**: backlog's 2-level/degree-3
  combination doesn't reach competitive accuracy per scipy (24x worse);
  even the corrected 3-level/degree-8 version regressed on real ulp (max
  2→3, avg 0.20→0.43). Not implemented.

- **Dedicated sinh/cosh kernels via reassociation (2026-07-08)**: real
  latency win for both (~11-12%) but throughput split by function (sinh
  better, cosh worse) and accuracy split the other way (cosh fine, sinh
  regressed ~12x avg ulp). Neither clears the bar. Not adopted. *When two
  functions share a reassociated intermediate, check accuracy for both —
  don't assume symmetry.*

- **Three-interval atan reduction (2026-07-08)**: accuracy case fully
  confirmed, but implementing it nearly doubled cost (latency +72%,
  throughput +99%) — the u-transform's division must resolve before the
  poly's own division can start for half the domain, two sequential
  full-latency divisions instead of one. Not adopted. *A flagged
  division-on-critical-path risk deserves an mca check before the accuracy
  work, not after.*

- **atan_latency's poly LP refit (2026-07-09)**: isolated fit predicted a
  modest 6% improvement — too weak a signal, and found nothing real (avg
  ulp 0.0516 vs 0.0517, statistically identical). Reverted. *An isolated LP
  prediction under ~10-15% has repeatedly turned out not worth the
  round-trip.*

- **Centered-variable refit for erfc's xa (2026-07-08)**: confirmed the
  exp2-centering finding transfers — centering measured ~73x worse despite
  erfc's different function shape.

- **Compensated-Horner accuracy tier for erfc (2026-07-08)**: double-float
  evaluation improved avg ulp modestly (~14%) but left max ulp flat — same
  shape as the domain-split rejection. Chased the root cause directly:
  erfc's own exponent expression has up to 87 ulp of error in plain f32
  arithmetic, before exp2_checked even runs — the same cause an
  already-rejected two_prod fix already found (not worth its mca cost). No
  code changed; closes a question three separate refits had each left open.

- **atan_poly Horner→Estrin (2026-07-08)**: measured backwards on every
  axis — fuzz accuracy worse (avg 0.068→0.072, max 3→5), mca throughput
  much worse (+156%). Not adopted. Closes the "audit every Horner poly for
  an Estrin win" line — every candidate now tried.

- **Rational (P/Q) refits for log_2/acos_poly (2026-07-08)**: log_2 ruled
  out by analogy (needs a division in the same fatal critical-path
  position as an already-rejected idea). acos_poly's fit failed to
  converge within 60s at two degrees — a reproducible dead end.

- **Rational (P/Q) refits for sinf_poly and erf_poly (2026-07-08)**:
  sinf_poly's fit converged but landed 5 orders of magnitude worse than the
  incumbent poly; erf_poly's fit timed out at 30s both degrees (same dead
  end as acos_poly). Closes the backlog entry: 4 for 4 rejected.

- **cos (fast tier) dedicated even poly in r² (2026-07-08)**: backlog
  itself predicted failure (relative error blows up near cos's zeros) — a
  quick scipy check confirmed it (unbounded relative error at the zero).
  Not implemented.

- **Simulated annealing / basin-hopping (2026-07-08)**: implemented as
  `tune_basin_hop`, tested on acos_poly. Coarse-grid candidate reported max
  ulp 2 (vs. descent's 3) — real fuzz showed it was actually a regression
  (avg 0.961 vs shipped 0.496, nearly double). Same `score()`-tuple
  max-first bias already documented for cbrt_throughput. Not adopted; kept
  as reference infra with the failure documented.

---

# Brainstorm backlog — UNTESTED

Every candidate needs: auto-vectorization check, exhaustive/fuzz accuracy
sweep, edgecheck.rs, and mca latency+throughput (both axes). FP divider is
nearly idle in every measured hot loop while FMA/mul ports are the
bottleneck, so "add a division to remove fmas" is a legitimate direction.

## Cross-cutting / methodology

- **Sollya `fpminimax`**: not installed. Candidates: `acos_poly` (max 4),
  erfc's n/d (max ~100, root cause already traced elsewhere — a tighter
  fit alone won't fix it).

- **Exhaustive/rlibm-style LP coefficient search**: tested on cbrt's
  correction poly — plain-minimax found max ulp 3→2 but avg *regressed*
  43% (minimax trades typical-case error for worst-case). A max-capped
  variant (minimize weighted-L1 subject to max weighted error never
  exceeding shipped) fixes this — use the max-capped version first on any
  future poly of this kind.

- **Batch/slice API tier** (`exp2_slice` etc.): lets the crate own
  vectorization instead of the caller's loop; natural home for a fast
  sincos too.

- **Explicit core::simd fallback tier**: only worth it if a LUT idea ever
  survives screening.

## log_2 / ln / log10

- **Shared denormal-rescale helper**: pure hygiene (log_2/ln/log10
  triplicate the tiny/xs/koff dance), no perf claim.

## sin / cos / tan

- **sincos_checked, shared reduction (2026-07-08)**: bit-exact refactor
  verified exhaustively, mca predicted ~19% throughput win — real
  wall-clock (after fixing two bench-shape artifacts) found it ~1.4%
  *slower*. Not adopted. *The one case this session where mca and real
  hardware disagreed in direction, not just magnitude — triangulate with
  more than mca on closures/CSE-sensitive changes.*

- **tan via mod-pi/2 reduction + dedicated poly (2026-07-08)**: backlog's
  premise (division near a pole amplifies error) was wrong — spot-checks
  showed the old sin(x)/cos(x) form was already 0-1 ulp at the actual
  poles. The real ~3000 max ulp comes from sin/cos's own large-x domain
  cliff, inherited regardless. New reduction measured worse everywhere and
  has a narrower safe domain. Reverted. *Check the actual behavior at a
  claimed failure point before implementing a fix.*

- **cbrt: 2/2 rational correction (2026-07-08)**: ~100x tighter in
  isolation, but both forms are already deep in f32's rounding-noise
  floor, and the extra division sits after the seed on the critical path
  with nothing to overlap (latency +28.5%, avg ulp regressed slightly).

## hypot / misc

- **remainder_checked beyond 2^24**: double-float q like sin_checked's
  reduction. Only worth it if a real use case needs it.

## Backlog round 2

### Cross-cutting / approximation theory

- **Automated evaluation-order search per poly**: fma reassociation is a
  per-poly coin flip (acos_poly's Estrin cost accuracy, atan_poly's cost
  throughput). Could enumerate Horner/Estrin groupings and score
  automatically.

- **ulp-weighted minimax fits**: weight coefficient fits by 1/ulp(f(x))
  instead of plain relative error, since ulp is a staircase near output
  power-of-two boundaries.

- **Worst-case patch lists**: compare-select against known bad bit
  patterns for functions with a tiny, stable set of failures. Checked
  erfc — not viable (~4900 bad points spread continuously, not
  concentrated). cbrt_accurate's one bad mantissa fits the bar but is an
  accepted won't-fix. Parked.

- **FTZ/DAZ feature flag**: a cargo feature assuming caller-side FTZ/DAZ
  would let every denormal branch (log family, cbrt, exp2_checked) become
  dead code.

- **f32x16/AVX-512 via #[target_feature]**: an explicit core::simd path
  sidesteps LLVM's refusal to force zmm-width; never measured.

- **Cross-check mca with uiCA and real perf counters**: mca's scheduling
  model has already produced caller-dependent surprises; uiCA is
  reportedly more accurate for this CPU.

### sin / cos

- **mod-pi/2 reduction with paired even/odd polys (2026-07-08)**: checked
  both backlog claims before writing Rust — degree only drops by 1
  coefficient (not "hard"), and real op count is ~25% *more*, not less,
  once both polys are evaluated unconditionally. Accuracy was never the
  blocker; killed by op-count reasoning alone. A hypothetical combined
  sincos might still pay off, but sincos_checked's own real-hardware result
  discourages building a new API for it.

- **Vectorized Payne-Hanek "exact" tier**: full-range correct reduction via
  a 2/pi mantissa product with a per-lane variable shift. Big job, optional
  given the graceful-degradation contract.

### atan / asin / acos

- **atan_poly denominator LP refit, numerator fixed (2026-07-09)**:
  essentially no movement (<1%, coefficients matched shipped to 6+ digits)
  — decisive enough to skip a Rust round-trip. Combined with a separate
  numerator-only refit (isolated fit predicted 11-15x, largest of the
  session, but real result was only ~0.9% with max ulp unmoved and the
  identical worst-case x before/after): atan's true worst-case residual
  isn't reachable by tuning either half; likely lives in the division
  itself or the a<1/a>=1 reciprocal-fold boundary.

- **Retune asin's 0.25 crossover after acos_poly changes (2026-07-08)**:
  true crossover sits around 0.26, but the reported max-ulp point sits
  inside asin_small's own domain regardless of threshold placement. 0.25
  already close enough; no change.

- **asin_small: one more Taylor term (2026-07-08)**: real 16% avg
  improvement, but max ulp just relocated to a tied co-bottleneck in the
  acos_poly branch, unchanged overall, plus a real perf cost. Not adopted.
  *When two branches tie at the same max-ulp value, fixing one just exposes
  the other.*

## Backlog round 3 (2026-07-09) — kitchen-sink brainstorm, all UNTESTED

Constraints as always: must autovectorize (branchless selects only), budget
≤0.5 avg / ≤2 max ulp, accuracy-for-speed trades fine inside that. FP
divider still nearly idle. Cross-checked against every rejection above —
each entry is either new or explicitly distinguished from its rejected
cousin.

### Approximation theory / fitting (cross-cutting)

1. **Round-based reduction for exp2/exp2_checked/exp10**: fit Q over
   f∈[-0.5,0.5] (k=round) instead of [0,1) (k=floor). Halving the interval
   scales minimax error ~2^-(d+1) — may allow degree 5→4 where the [0,1)
   degree probe was rejected. Also kills exp10's floor-adjust select pair
   (its reduction is natively round-based). Magic-round already proven
   faster than vroundps in exp.
2. **Gather-based real LUTs (vgatherdps)**: 8/16-entry table indexed by top
   mantissa bits — log family (per-interval rcp + log pair, shorter poly),
   exp2 (2^(i/16) exact). Distinct from the rejected select-tree "LUT"
   (blend emulation); a hardware gather is a different cost model. Screen
   with mca before any fitting.
3. **Coefficient ulp-neighborhood exhaustive search**: for each shipped
   poly, enumerate all coefficient tuples within ±k ulp (k~2-4) of the LP
   solution, scored on the real crate — rlibm-lite. Catches the
   f32-quantization effects the LP's continuous model misses (the exact
   failure mode of the exp2 LP rejection).
4. **Per-function transformed-variable fit search**: fit in u=s/(s+2),
   u=s·(s+a), etc., searching over the transform family. Distinct from
   centered-variable refits (rejected — that only moved the origin);
   a nonlinear transform changes curvature matching.
5. **Ulp-staircase-aware LP grids**: densify fit grids near output
   power-of-2 boundaries where ulp weight steps 2x. Complements the
   existing 1/ulp-weighting backlog entry (that's weights; this is node
   placement).
6. **Joint threshold+both-branch LP for every 2-branch function** (sinh
   0.5, asin 0.25, erf 0.28, expm1 0.5, atanh if split): optimize crossover
   and both coefficient sets together, each branch max-capped. The erf
   boundary sweep that found "no headroom" held coefficients fixed.
7. **Round-off budget audit per function**: enumerate every rounding on the
   critical path with a bound, attack the largest term. This is exactly how
   exp's Cody-Waite fix was found; do it systematically for the remaining
   >2-max-ulp functions (asin 9, expm1 6, tanh 6, sinh/cosh 5, erf 5).
8. **Binary-function worst-case mining**: unary functions get exhaustive
   sweeps; powf/atan2/hypot/remainder only get fuzz. Guided search
   (branch-and-bound over exponent-pair classes, or fixed y/x ratio
   classes for atan2) would find real worst cases fuzz misses.
9. **Structured-error probes**: plot per-function error vs mantissa and vs
   exponent separately; periodic structure invites a cheap structural
   correction (one select or exponent-derived fma) instead of a refit.
10. **Monotonicity/oddness harness metrics**: flag outputs non-monotonic
    where the true function is monotone, and f(-x)≠-f(x) for odd
    functions — localizes bug classes that avg/max ulp only aggregates.
11. **Differential testing vs sleef/core-math/rlibm** built locally, not
    just f64-rounded references — also catches double-rounding artifacts in
    accuracy.rs's own reference path.
12. **Standing test: every `_unchecked` bit-matches its checked sibling on
    the documented domain** (several doc comments promise this; nothing
    enforces it).
13. **Generalize the cbrt_accurate recipe** (cheap ≤1-ulp core + one Df32
    Newton step) into a template: candidates rsqrt_accurate,
    exp_accurate/ln_accurate (each is the other's Newton residual),
    sin_accurate near zeros. Paper-screen the residual budget first.
14. **ln_accurate/log2_accurate tier from existing log2_df**: log2_df
    already returns a Df32; collapsing it carefully (Cody-Waite vs LN2 as
    a Df32 constant) gives a higher-accuracy log tier nearly for free —
    the machinery exists, only the collapse is new.
15. **Integer fixed-point poly evaluation** for mantissa-only reductions
    (log's s): i32 mul-high chains free up FMA ports. Precedent warning:
    parity()'s integer version lost to FP ports — but that was 3 ops, not
    a whole poly; port-pressure math differs at scale.

### exp family

17. **exp: weave t1 into the poly like exp2_checked does (tried 2026-07-09,
    rejected)** — implemented exactly as described (Q(r) = 1 + c0·r + ... +
    c3·r^4, `p = fma(q, t1*r, t1); p*t2`). Real accuracy win, confirmed
    exhaustively: exp avg/max ulp 0.0745/3→0.0522/2, cascading to expm1
    (0.1304/6→0.1273/5), sinh_throughput (0.0723/5→0.0677/4),
    cosh_throughput (0.0507/4→0.0466/3), tanh (0.1457/6→0.1447/5). But mca
    showed a real, 3x-reproducible throughput *regression* on exp itself
    (1.327→1.393 cyc/elem, latency flat) and expm1 latency +2 cyc
    (74→76). Root cause (bottleneck-analysis + hand depth/instruction-count
    accounting): unlike exp2_checked's Q(f), which has no cheap way to fold
    in its "+1" (the un-weaved form would need an *extra* fma to
    materialize `1+f·q` before scaling), exp's original P(r) already got
    its "+1" for free from a plain, fully-parallel `r+1.0` add outside the
    poly's dependency chain. Weaving replaces that free add with an extra
    `t1*r` multiply landing on the same contended fma/mul ports as the
    poly itself — old and new both total exactly 7 instructions at the
    same critical-path depth (5), so it's a pure lateral shift onto busier
    ports, not a real reduction. Also confirmed on paper (not implemented,
    correctness-fatal): precomputing `t12=t1*t2` off the critical path to
    collapse the tail to one multiply would break the `exp(88.37628)`
    edgecheck (k=128 boundary) by prematurely overflowing `t1*t2` to inf
    before the sub-1 polynomial factor brings it back down to a finite
    result — the two-stage multiply is load-bearing for that boundary, not
    incidental. Reverted (bit-identical to prior HEAD). *A "drop one
    rounding" premise from one function's poly shape doesn't transfer to a
    sibling with a differently-derived poly — check whether the "+1" (or
    equivalent identity term) was already free before assuming the weave
    saves anything.*
18. **exp_checked tier**: exp currently has no full-range sibling (exp10
    does). Clamp + the existing k1/k2 split — trivial, closes an API gap,
    and callers like sigmoid/tanh could then drop their own ad-hoc clamps.
19. **exp10 third Cody-Waite word**: LOG10_2 reduction is 2-word; a third
    word is one fma off the critical path. Survey exp10's actual max ulp
    first to see if there's anything to collect.
20. **exp2int construction via pure-integer path**: k is already an exact
    small integer in f32 form; compare `(k+383).to_bits()<<8` against
    cvt-to-i32 + shift-23 + add on vector ports (saturating-cast problem
    doesn't apply — k is bounded). Probably a wash; cheap to check asm.
21. **exp2_checked: k1 from bit-twiddled k instead of a second
    magic-round** — k1b's fma+sub pair might be replaceable with an
    integer halving of k's already-integer bits. Screen via asm/mca only.

### log family

22. **log1p output-side correction refit**: the rejected small-|x| branch
    added a whole poly (+48% mca). Instead refit *ln's own* coefficients
    jointly with log1p's c/u correction term as part of the objective —
    zero new ops, may shave log1p's max 4.
23. **log10 exact-decade check**: verify log10(10^n) is exact for n∈[-38,38];
    if not, decide whether a fixup is worth it (probably document instead).
24. **log_2/ln/log10 shared-core macro**: the three `_normal` bodies are
    identical modulo constants (and log2_df quadruplicates it). A macro or
    generic-const core removes 3 hand-synced copies — hygiene, zero perf
    claim, reduces refit-application errors.
25. **log_2 denormal path: fold the ×2^24 rescale into the wrapping_sub
    magic** — the exponent extract could absorb a constant offset, but the
    mantissa of a denormal isn't normalized, so this needs the multiply
    anyway for the mantissa bits. Probably dead on arrival; kill it on
    paper in 5 minutes.
26. **koff-free unchecked-log fast path audit**: log_2_unchecked passes
    koff=0.0 through an add that's provably dead (k + 0.0 exact) — check
    LLVM actually deletes it (the round_x_over_pi dead-add precedent says
    removal can even *help* scheduling... or hurt; asm check only).

### sin / cos / tan (radians, half-turns, degrees)

27. **sinpi/cospi: two_prod(π, r) + derivative correction** — π·r's single
    rounding is likely these functions' dominant error; e = two_prod
    residual, add e (or e·(1-y/2) using the poly's own y) into the result.
    Two cheap ops, targets the one inexact step in an otherwise-exact
    reduction.
28. **sind/cosd: same trick for d·DEG_TO_RAD_SMALL** — 2-constant split of
    π/180 (hi with zeroed tail bits so d·HI is exact, lo folded via fma).
29. **tanpi / tand**: new functions from existing pieces (sinpi/cospi,
    sind/cosd ratios). tan's period-π means parity cancels in the ratio —
    check whether the parity xors can be skipped entirely for the pi/deg
    variants.
30. **sinpi/cospi/sind/cosd harness coverage**: confirm these four are in
    accuracy.rs's sweep at all; they were added later than the harness.
31. **cospi large-x probe**: kb=(x-0.5)+MAGIC — for x with ulp>1 the -0.5
    rounds away entirely (the exact pre_offset bug class fixed in
    round_x_over_pi). Past 2^23 every f32 is an even integer so cos(πx)=1;
    check the parity/select path actually lands there rather than by luck.
32. **Intermediate sin/cos tier (|x|≤~1e5)**: single extra correction word
    over the fast tier's 4-fma Cody-Waite, well short of checked's full
    double-float q — a third point on the speed/domain curve if any user
    workload actually sits there. Only on demand.
33. **sin fast tier: drop PI_D for a fitted 3.5-word split** — distinct
    from the rejected 4→3 chain cut (that reused the same constants and
    died near zeros); a re-*fitted* 3-word split with the third word chosen
    to minimize worst-case residual near zeros might survive. Screen in
    Python against the zeros of sin specifically before touching Rust.

### asin / acos / atan / atan2

34. **Dedicated asin-only coefficient copy of acos_poly**: every joint-fit
    attempt died protecting acos (three separate rejections). Duplicating
    the 7 literals decouples the two callers permanently — asin gets its
    own LP fit over its own domain weighting, acos keeps its protected
    values, zero runtime cost (same instruction count, different
    constants).
35. **atan: correct the 1/a fold's division rounding** — e = fma(y, a,
    -1.0) gives the reciprocal's residual; Δatan ≈ -e·y/(1+y²), and
    atan_poly already computes a denominator ≈(1+y²)-shaped value. Two-ish
    extra ops aimed exactly at the "worst case lives in the division /
    reciprocal-fold boundary" conclusion from the 2026-07-09 numerator
    refit.
36. **acos_accurate opt-in tier**: Df32 π/2 constant + two_prod'd
    sqrt(1-a)·poly product. Distinct from the rejected Df32 leading-term
    split — that changed the *shared default* and broke asin; a separate
    tier touches nothing shared.
37. **asin small branch: check 1-a rounding in [0.25,0.5)** — 1-a is only
    Sterbenz-exact for a≥0.5; quantify what the sub-0.5 rounding costs
    through sqrt+poly before deciding anything.
38. **asin via atan2(x, sqrt((1-x)(1+x)))**: different algorithm entirely
    (correctly-rounded sqrt, atan max 4). Probably slower (division inside
    atan) but it's a one-evening accuracy ceiling probe for asin's max 9.
39. **atan2 octant-symmetric exhaustive harness**: sweep all x at a few
    thousand fixed y/x ratios (and vice versa) — turns the binary-input
    problem into affordable near-exhaustive slices.
40. **atan_latency: fold FRAC_PI_2-p select into sign trickery** — the
    a<1.0 select and final mulsign might merge into one xor+select. Asm
    check; remember the mca mix() sign blind spot — throughput mode only.

### hyperbolics / sigmoid

42. **tanh via exp_pos_neg ratio**: (ep-en)/(ep+en) + a small-|x| Taylor
    branch (x - x³/3 + 2x⁵/15 - 17x⁷/315), replacing the expm1 route and
    its 2·x clamp interplay. Division is idle; may beat expm1's max 6.
44. **sigmoid: branch on sign for conditioning** — for x≫0, 1/(1+e^-x)
    is well-conditioned; for x≪0 compute e^x/(1+e^x) instead (both
    branchless-selected). Survey whether current accuracy actually needs
    it first.
45. **softplus/log1pexp(x) = ln(1+e^x)**: new function, ML-relevant;
    branchless as max(x,0) + log1p(exp(-|x|)). All pieces exist.
46. **atanh via single log1p on |x| + mulsign**: atanh is odd — compute
    0.5·log1p(2a/(1-a)) on a=|x| only, restore sign. The rejected
    single-log1p form died near x≈-1 where 2x/(1-x) cancels; on the
    positive side there is no cancellation (argument →+∞), so odd-symmetry
    sidesteps the entire failure mode. Halves the log1p count per call.
48. **sinh_accurate/cosh_accurate tier**: Df32 through the exp combine —
    only if a user asks; max 5 is comfortably documented.

### erf / erfc

49. **erfc exponent via mantissa-mask hi/lo split** (Cephes trick): xh =
    xa with low ~12 mantissa bits masked → xh² exact; exponent =
    -xh²·L - (xa-xh)(xa+xh)·L. Cheaper than the rejected two_prod fix
    (mask+sub+fma vs two_prod's mul+fma each use) aimed at the same traced
    87-ulp exponent error behind erfc's max 109.
50. **erf tail via expm1-shape**: b = -expm1_2(p)·sign form instead of
    1 - exp2(p) — removes the 1-(≈1) subtraction that's mildest exactly at
    the 0.28 crossover where erf's confirmed fragile spot sits. Needs an
    exp2m1 helper (see below).
51. **erfcx(x) = e^{x²}·erfc(x)**: new function — just the n/d rational,
    no exp at all; sidesteps the exponent-error bottleneck entirely and is
    what numerics users often actually want in the tail.
52. **erf_poly's a0 ≈ 3.4e-5**: suspiciously near zero — refit with a0
    pinned to exactly 0.0 (frees a degree of freedom for the other
    coefficients and deletes one fma if it holds). LP with max-cap,
    exhaustive-verify; cheap experiment.
53. **erfc negative-side accuracy survey**: the w=2 branch computes 2-y;
    quantify whether the x<0 half has its own error structure — every
    refit so far scored the whole domain blended.

### cbrt / sqrt / hypot / pow / remainder

54. **cbrt seed division codegen check**: confirm LLVM lowers the u32 `/3`
    to multiply-high in the vectorized loop (vpmuludq path) rather than
    anything worse; if not, hand-write the exact magic-multiply. Asm-only
    check, no accuracy risk (exact either way).
55. **rcbrt(x) = x^(-1/3)**: new function — negate-exponent-third bit
    seed + its own correction poly; division-free, useful in physics
    kernels, and 1/cbrt(x) costs an extra rounding this avoids.
56. **cbrt_accurate via Halley from a cheaper seed**: cubic convergence
    might let a cbrt_fast-grade (~5 ulp) seed reach 0.5 ulp in one Df32
    Halley step, deleting cbrt_normal's poly from the accurate tier.
    Paper-screen the error budget first.
57. **pown_small tier (|n| ≤ 255)**: 8 unrolled iterations instead of 32 —
    4x fewer squarings for the overwhelmingly common case, still
    branchless. Also **pown_const<const N: i32>**: compile-time exponent
    unrolls exactly (the "same exponent, varying base" pattern that
    motivated pown's redesign).
58. **powf: root-cause the log2_df/exp2_checked_df ~150 max ulp** (the
    2026-07-08 entry stopped short): candidate mechanisms — exp2's
    double rounding into denormals via the t2 multiply, or Df32 mul's
    dropped lo·lo term. Instrument before fixing.
59. **exp2_checked_df: second-order lo correction** — exp2(lo) ≈ 1 +
    lo·ln2 + (lo·ln2)²/2; one fma. Only matters if #58 says the
    first-order truncation is the binding term.
60. **powf special-exponent select tier**: y ∈ {1, 2, 0.5, -1} handled by
    exact ops behind one select ladder. Costs every call; likely rejected
    on mca — but powf is expensive enough that the relative cost may be
    tolerable. Screen with mca first.
61. **rhypot(x,y) = 1/hypot**: new function, divider is idle; normalizing
    2D vectors is the dominant hypot use case and this deletes the
    caller's division.
62. **hypot: fma pairing choice (tried 2026-07-09, rejected)** — implemented
    max-first pairing for hypot/hypot_unchecked/hypot_checked (compare+select
    on the signed operands). Real, consistent avg-ulp win (~15% better: hypot
    0.0338→0.0288, hypot_unchecked 0.0339→0.0288, hypot_checked 0.0149→0.0127,
    stable across repeat fuzz runs) but max ulp never moved off 1 as
    predicted, and mca showed a real, deterministic perf cost: hypot latency
    21.11→26.20 cyc (+24%), throughput 0.766→0.797 (+4%); hypot_checked
    58.22/1.278 (+1.8%/+8.5%). Tried a second variant using hardware
    `.max()`/`.min()` instead of compare+select (cheaper instructions in
    principle) — measured *worse* latency still (29.20 cyc, +38% vs baseline),
    confirming the real cost isn't the compare-vs-max/min instruction choice
    but that determining operand order at all forces a serial step (abs +
    compare/max) in front of the fma, which the naive `fma(x,x,y*y)` avoids
    entirely by squaring both raw inputs immediately in parallel — nothing
    to reorder before the multiply can start. (`.max()`/`.min()` also isn't a
    correctness-neutral swap on its own: Rust's `f32::max`/`min` follow IEEE
    maxNum/minNum and *discard* NaN — return the other operand — rather than
    propagate it, unlike the plain compare+select `if a>=b` used here, which
    happens to preserve NaN because the same condition picks both outputs
    complementarily. Would have needed an extra NaN-restoring select on top,
    even before the latency finding killed it outright.) Reverted, bit-
    identical to prior HEAD. *A pairing/reordering change that looks free in
    isolation (same op count) can still cost real latency if it inserts a
    decision in front of an op that previously had no upstream dependency at
    all — count what's on the critical path *before* the op, not just at
    it.*
63. **fmod family**: trunc-based sibling of remainder/remainder_checked —
    C-parity gap in the API, same machinery, mostly copy-paste.
64. **remainder_checked: widen past 2^24 with a 2-word q** — already in
    backlog round 1; still unclaimed.

### New API surface / tiers

65. **exp2m1 / log2p1** (C23): exp2m1 falls out of the expm1 technique on
    exp2's machinery (and unblocks idea #50); log2p1 = log1p·LOG2_E with
    the usual Cody-Waite care.
66. **exp_m1_over_x(x) = expm1(x)/x**: the well-conditioned primitive
    behind financial/ODE kernels; expm1's Pade branch is literally already
    this shape internally (numer/denom both have the x factored).
67. **sinc(x) = sin(πx)/(πx) via sinpi**: removable singularity handled by
    one select; sinpi's exact reduction makes this accurate everywhere —
    DSP users currently hand-roll it badly.
68. **lgamma (Stirling + reflection)**: big job, listed for completeness —
    the largest gap vs libm's function set that fits this crate's
    branchless style.
69. **logaddexp(a,b)**: max + log1p(exp(-|a-b|)) — two existing calls,
    branchless, ML-relevant.
70. **Slice/batch API + runtime multiversioning**: backlog round 1 has the
    slice tier; add `is_x86_feature_detected` dispatch at the slice level
    (per-call dispatch is un-inlinable, per-slice is free).

### Codegen / build / measurement

71. **Blend lowering audit under AVX-512VL**: check whether f32 selects
    lower to vblendvps or mask registers (vmovaps{k}) with target-cpu=
    native, and whether mask ops relieve port 5 pressure in the hottest
    loops.
72. **`-C llvm-args=-force-vector-interleave=N` sweep** (2/4/8) on the mca
    and wall-clock harnesses — interleave choice is LLVM's guess; several
    past "scheduler made a worse choice elsewhere" surprises suggest the
    default isn't always right.
73. **PGO (+ BOLT) probe on the bench binaries**: mostly affects branchy
    code, which this crate avoids — cheap to try once, likely a null
    result, worth knowing.
74. **codegen-units=1 + lto sweep for the bench profile**: confirm the
    harness isn't measuring artifact boundaries.
75. **Asm-grep CI test**: assert zero scalar fallbacks (no `vsqrtss`/
    `vdivss`/call instructions) in the vectorized loop bodies of every
    public function — turns the "must autovectorize" rule into a test.
76. **Continue the AVX-512 probe** (examples/scratch_avx512_probe.rs is
    sitting untracked): decide whether explicit `core::simd` f32x16 tiers
    beat the auto-vectorized ymm baseline enough to justify a feature
    flag, then delete or commit the scratch file.
77. **mca/wall-clock disagreement detector**: script that runs both on
    every candidate and flags direction disagreements automatically
    (sincos_checked precedent) instead of relying on remembering to
    triangulate.
78. **Throughput harness with N independent input streams**: current
    quickbench shape may be ILP-limited in ways that hide or exaggerate
    wins; 4 interleaved streams approximates a real vectorized caller
    better. (Also sidesteps part of the mix() sign blind spot — but only
    part; sign-only work still needs the exhaustive bit-check.)

### Accuracy micro-fixes / surveys (cheap to check, maybe nothing there)

79. **Survey sinpi/cospi/sind/cosd/exp10/sigmoid max ulp** — several newer
    functions have no recorded exhaustive numbers in readme.md; can't
    prioritize what isn't measured.
80. **exp2 poly evaluated as 1+f·Q vs direct P(f)=2^f with c0=1 pinned**:
    the current Q form multiplies by exp2int·f then adds exp2int; a direct
    P form changes the rounding structure of the last two ops. Paper
    analysis first — the max is already 1, so margin is thin (see the LP
    rejection lesson).
81. **sinf_poly's copysign(x)**: now that flip-before-poly is used in
    checked tiers, verify the copysign is still load-bearing for every
    remaining caller (it was added for the x=±0 case; sinpi/cospi/sind/
    cosd route sign differently).
82. **log1p at x=+1.0 boundary**: u=2.0, Sterbenz window edge — one-input
    check that c is still exact there.
83. **hypot_checked denormal-pair path**: pre-scale by 2^24 then exponent
    trick — sweep a grid of denormal×denormal pairs specifically (the
    Python prototype sampled broadly, maybe thinly there).
84. **erf/erfc NaN sign convention audit**: mulsign-based paths can flip
    NaN payload signs; C99 doesn't care but a cheap consistency check
    against std across all specials would close the book.
85. **atan2(±0, negative-finite) etc. full C99 special-case matrix as a
    test table** — atan2's specials were fixed piecemeal; one table test
    locks all 16+ cases.
86. **remainder's ties-away vs IEEE ties-even**: documented divergence
    from IEEE 754 remainder — add a `remainder_ieee` variant via
    round_ties_even (vroundps has the mode; may even be the same cost).
87. **powf(±1, huge y)**: log_2(1)=0, 0·y=0, exp2(0)=1 — fine; but
    powf(1+ulp, 3e38): log2≈8.5e-8, ·3e38 overflows f? No — f32 holds it.
    Check the k-clamp path saturates correctly rather than wrapping.
88. **exp10 near the decade boundaries**: k=round(x·log2_10) crossing
    integers at x = n·log10(2)-ish points — sweep a dense band around
    every such crossing (the floor-adjust select is the only branchy-ish
    logic in the function).

### Longer shots / research-flavored

89. **Bit-sliced two-for-one**: evaluate sin and cos polynomials sharing
    y=r² registers across the *same* vector when the caller wants both —
    a sincos slice API (not scalar API, which already failed) where lane
    pairing amortizes the reduction. Only viable inside a slice tier.
90. **Newton-free correctly-rounded sqrt-composites**: rsqrt/rhypot final
    accuracy via one fma-based residual step (e = fma(r,r·x,-1)-style) —
    the divider-idle finding means the extra division for the step is
    nearly free in throughput terms.
91. **Exhaustive-verified minimax over *reduced* domains** (true rlibm):
    for polys whose reduced input takes ≤2^26-ish distinct values (exp2's
    f after quantization? log's s per exponent?), solve the actual integer
    LP over rounding intervals instead of a continuous fit. Check the
    input-multiplicity math per function first — most reductions don't
    quantize enough.
92. **Domain-specific fast-math contract tiers**: a `finite-math-only`
    cargo feature gating away every inf/nan select in checked functions
    (complements the FTZ/DAZ backlog entry, which only covers denormals).
93. **Auto-generated `_unchecked` variants via macro**: every checked/
    unchecked pair is hand-maintained; a macro emitting both from one body
    with cfg'd guards removes drift risk (the doc-promise test in #12
    becomes trivial).
94. **Cross-function CSE tiers**: sigmoid+tanh (share exp), sin+cos of the
    same angle (slice-level), erf+erfc — paired entry points for callers
    that need both, amortizing the shared prefix that scalar fusion
    attempts kept failing on (the win the failed sincos sought lives at
    the API level, not inside one scalar call).
95. **Stochastic rounding harness mode**: run accuracy sweeps with the
    final fma's rounding perturbed ±1 ulp to measure how close each
    function sits to a rounding boundary — identifies which maxes are
    "one lucky rounding" vs structural, prioritizing refit targets.
96. **Interval-arithmetic self-audit build**: a cfg that swaps f32 for an
    interval type in the `_normal` cores to machine-verify the "this add
    is exact / Sterbenz applies" claims scattered through the comments —
    several past bugs (pre_offset, e3t sign) were exactly wrong claims of
    this kind.
97. **exp/exp2 denormal-output double-rounding**: the t2 multiply into the
    denormal range rounds once by design — verify against a
    correctly-rounded reference specifically on outputs in [2^-149,
    2^-126) (powf's residual max-ulp neighborhood, per #58).
98. **Karatsuba-style Df32 multiply**: doublefloat.rs's mul does 2-3
    two_prods; for the powf chain, an error-bounded cheaper mul (drop
    lo·lo, keep cross terms in one fma) might cut powf_checked's +61%
    throughput cost — audit what Df32::mul actually does first.
99. **Precision-tapered polys**: evaluate the high-order (small-magnitude)
    poly tail in a *cheaper* form (fewer fmas, plain Horner) and only the
    dominant terms carefully — the tail terms' own rounding is provably
    below the result ulp. Candidate: log family's l3/l4 tier.
100. **A cost model for "add a division"**: the divider-idle finding keeps
    paying off (cbrt rcp, sinh_throughput, log1p) — write down the actual
    reciprocal-throughput arithmetic (divider ports vs FMA ports per
    vector width) so candidates can be paper-screened instead of
    mca-round-tripped one at a time.
