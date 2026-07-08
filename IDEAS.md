# Ideas tried and rejected

Budget: ≤0.5 avg ulp, ≤2 max ulp. This file records ideas that were
actually implemented and/or measured, then reverted or not adopted — kept so
future sessions don't re-attempt them without new information. Adopted
changes live in git history / readme.md, not here. A separate untested
brainstorm backlog lives at the bottom of this file.

## Cross-cutting

- **Degree-reduction probes (2026-07-07)**: lolremez screening rejected all
  three candidates. `log_2` deg-9→8 already failed earlier (max ulp 3-5 vs a
  2 cap). `exp2`'s Q poly deg 5→4: estimated max relative error
  1.01e-8→4.07e-7 (40x worse). `sinf_poly` deg 9→7: 6.97e-9→1.24e-6 (178x
  worse). None carried through to a full implementation.

- **exp2/log_2 coefficient refit (2026-07-07)**: ran the coordinate-descent
  tuner on both polys' existing coefficients. Result was bit-identical to the
  shipped coefficients — a literal zero-move local optimum. Both had already
  been tuned in an earlier session; no headroom left.

## exp2 / exp / sin / cos / tan

- **Fused sincos / direct tan via shared reduction (2026-07-07), two
  attempts, both reverted**: sin(x)=(-1)^q·sin(r), cos(x)=(-1)^q·cos(r) means
  tan(r)=sin(r)/cos(r) needs no parity, just one shared reduction. Attempt 1:
  fit a dedicated `cosf_poly` (excellent isolated fit, ~3.6e-10) —
  catastrophic regression near tan's poles (avg ulp 0.33→1.05, max ulp
  ~3000→32M) because the additive `1+y·R(y)` form cancels badly near
  r=π/2. Attempt 2: use the cofunction identity `cos(r)=sin(π/2-|r|)` with
  `sinf_poly` for both — worse (max ulp →3.4B), because `FRAC_PI_2 -
  r.abs()` is a plain f32 subtract that itself cancels near a pole. Both
  reverted; a real fix needs multi-word-precision reduction for the
  cofunction transform, comparable in complexity to `reduce_pi` itself.

- **sinf_poly quantized refit (2026-07-07)**: tuned against f64::sin over the
  poly's domain. Max ulp unchanged (2→2), avg ulp barely moved
  (0.00248→0.00244) — already near f32's precision floor. Not applied.

- **expm1 Pade degree bump, numerator degree 3 → 5 (2026-07-08)**: added a
  new odd term (new `expm1_near0_deg5_c` in `tune.rs`, `"expm1_near0"`
  arg), seeded at 0.0 so it starts bit-identical to the shipped 5-
  coefficient form. Coordinate descent left the new coefficient at exactly
  0.0 and every other coefficient unmoved — a zero-move local optimum,
  same signature as the `exp2`/`log_2` coefficient-refit entry above and
  this session's `acos_poly8`/ln/log10 refits. No headroom found on
  `tune.rs`'s own grid (known to be coarser than `accuracy.rs`'s exhaustive
  sweep, but a *zero*-move result doesn't need the denser sweep to trust —
  there's nothing for it to reveal). Not adopted; `src/lib.rs` never
  touched.

## cbrt family

- **Seed constant joint search, degree-2 poly (2026-07-07)**:
  coordinate-descended a degree-2 correction (dropping c4) across 41 seed
  offsets around the shipped constant. Best found across all seeds: max ulp
  112 / avg 29, ~50x over budget. 3 coefficients can't correct this seed's
  error regardless of seed choice. Also moot as a speed idea — degree-2
  Horner is the same fma-depth as degree-3 here. Not adopted.

- **Integer-division-free seed, `(bits>>16)*0x5556` instead of `ax/3`
  (2026-07-07)**: codegen confirmed cheaper (3 instructions vs 5 for the
  division). But refitting the degree-3 correction poly for the new seed's
  error distribution only reached max ulp 33 / avg 5.07 — ~16x over budget.
  The seed itself is too coarse for this correction poly to compensate. Not
  adopted; would need the seed's own magic constants jointly re-derived to
  go further.

## log_2 / ln / log10

- **ln_normal/log10_normal: fuse the trailing `+ k*LN2_LO` into the
  preceding fma, `fma(k, LN2_LO, fma(p, s, k_hi))` instead of `fma(p, s,
  k_hi) + k * LN2_LO` (2026-07-08)**: the backlog framed this as "one op
  shorter and strictly one fewer rounding" (the standalone `k * LN2_LO`
  multiply plus the final add are two separate roundings; fusing them
  into one outer fma removes one). True in isolation, but measured
  outcomes on both axes were negative or flat: mca latency got
  *deterministically worse* by exactly 1 cycle on both functions
  (ln/log10 55.86→56.86 cyc, reproduced twice, not noise — same
  "removing an op frees the scheduler to make a worse choice elsewhere"
  pattern logged repeatedly elsewhere in this file), throughput was a
  wash (ln 1.714→1.709, log10 1.714→1.714 exactly unchanged). Worse,
  the predicted accuracy win didn't materialize either: a 100M-sample
  fuzz (paired via `git stash`) gave avg/max ulp 0.1168/3 → 0.1168/3
  (ln) and 0.1271/3 → 0.1270/3 (log10) — identical within sampling
  noise. Root cause: `k * LN2_LO` doesn't depend on the poly result, so
  it was already computed off the critical path in parallel with the
  fma before this change; the "extra rounding" it removes is on a term
  small enough (LN2_LO is the tiny low word of the Cody-Waite split)
  that it isn't actually contributing to the measured ulp in practice.
  Reverted; `src/lib.rs` restored to `fma(p, s, k_hi) + k * LN2_LO`.

- **log1p small-|x| (<0.25) dedicated Taylor/minimax branch, screened via
  a partially-implemented `tune.rs` scratch infra found already in the
  working tree (2026-07-08)**: the backlog's framing ("free perf-wise if
  it replaces work rather than adding a third arm") doesn't hold given
  this crate's established branchless-select convention (`asin`/`acos`'s
  own history documents this explicitly: every domain branch is
  computed *unconditionally*, then blended with a select) — `log1p`
  currently has exactly one arm (the `ln(u)+corr` formula, used for
  every x), so adding a small-x poly branch would be a strict *addition*
  of a whole extra poly evaluation to every call, not a replacement.
  Confirmed the accuracy upside is marginal anyway before spending more
  effort: `log1p` restricted to `|x|<0.25` currently measures avg/max
  ulp 0.0732/4 (accuracy.rs fuzz, ad hoc probe); the
  WIP tune.rs candidate (`log1p_small_c`, a 9-coefficient Horner fit
  directly against `x.ln_1p()`, forced c0=1.0) reached max 3/avg 0.0683
  on tune.rs's own coarser grid — a small, unconfirmed-at-full-density
  edge, not obviously worth a whole extra poly's cost on every call.
  Not implemented in `src/lib.rs`; the found WIP infra in `examples/
  tune.rs` (`log1p_small_c` + the `"log1p_small"` tuning branch) was
  reverted rather than committed, since it was never brought to a
  real verdict and this session's reasoning already closes the question
  well enough to not re-attempt without new information (e.g. an actual
  measured throughput/latency cost from implementing it, if someone
  wants to check whether the extra branch is cheap enough to be worth
  0.3-1 ulp).

- **Direct minimax refits for ln/log10, tuned against their own objective
  instead of log_2's rescaled coefficients (2026-07-08)**: coordinate-
  descended `ln_poly_c`/`log10_poly_c` (new permanent `tune.rs`
  infrastructure, `"lnlog10"` arg) against `ln`/`log10` directly, matching
  the shipped Estrin structure and Cody-Waite k-combine exactly, c[0]
  fixed at its mathematically-required exact value (same as `log_2`'s own
  tuning). Same outcome as the erf near-zero/tail refits above: essentially
  a zero-move local optimum (ln avg ulp 0.2344→0.2341, log10
  0.2551→0.2539, both max ulp unchanged) — individually rounding log_2's
  own coefficients by a fixed constant was already close enough to a
  direct fit that there's no meaningful headroom left. Not adopted.

## hypot / misc

- **Compensated hypot, `e = fma(r, -r, s); r + e/(2r)` (2026-07-08)**:
  implemented exactly as scoped (guarded at `r == 0` for the `x==y==0`
  case). Fuzz accuracy (paired via `git stash`, 10M samples both sides)
  showed **no measurable improvement at all** — 0.0337→0.0338 avg ulp,
  max ulp 1 unchanged both before and after, i.e. within pure sampling
  noise. `hypot`'s existing single-fma-then-sqrt was already close enough
  to correctly rounded (sqrt is itself correctly rounded relative to the
  once-rounded `s`, and that one rounding rarely flips the sqrt's own
  rounding direction) that there was no real residual left to recover.
  Real, large cost for zero benefit: mca latency 21.11→44.05 cyc (+109%),
  throughput 0.766→1.436 cyc/elem (+87%) — the extra fma/division/select
  chain roughly doubled hypot's cost on both axes. Reverted before even
  checking edgecheck.rs.

## asin / acos / atan / atan2

- **acos_poly Horner→Estrin restructuring (2026-07-07)**: regrouped 6-deep
  Horner into 3-deep Estrin (same coefficients). Real latency win (asin
  -15.4%, acos -21.6%) but not accuracy-neutral: fma reassociation regressed
  asin max ulp 9→12, acos 4→5. Retuning coefficients for the Estrin form
  specifically made it worse (acos max ulp →6). Reverted — `acos`'s accuracy
  was specifically protected by an earlier constrained refit, and this would
  have undone that guarantee for a latency-only win.

- **erfc's n/d rational chains Horner→Estrin (2026-07-07)**: same
  restructuring on the two degree-4 Padé chains. Smaller theoretical win (5
  coefficients only saves 1 depth level) and it showed: latency -1.3%,
  throughput +1.3% (a wash), and a real accuracy cost (avg ulp +2.7%). Not
  adopted.

- **erf's near-zero Padé branch refit (2026-07-07)**: tuned against the
  `numer/denom` formula for |x|<0.28. Max ulp unchanged (3→3), avg ulp moved
  <0.3%. No meaningful headroom; not applied.

- **erf's tail branch (erf_poly) refit (2026-07-07)**: same recipe against
  the tail formula for xa in [0.28,10]. Max ulp unchanged (4→4), avg ulp
  moved <0.3%, stable across a 10x denser grid. Not applied.

- **acos_poly refit against joint acos+asin objective, unconstrained variant
  (2026-07-07)**: an unconstrained joint metric (`max(acos ulp, asin ulp)`)
  improved the joint score but let acos's own exhaustive max ulp regress
  4→5 — a real cross-function tradeoff. Rejected in favor of a constrained
  variant (kept; not recorded here) that protected acos's metric by
  construction.

- **acos_poly degree 6 → 7, i.e. an 8th coefficient (2026-07-08)**: the
  backlog framed this as "one fma of throughput for a whole extra degree
  of freedom" (implicitly assuming Estrin, which acos_poly isn't —
  it's shipped Horner, so an 8th term would cost one full extra depth
  level too, not just one throughput fma). Turned out moot regardless:
  coordinate-descending a new `acos_poly8_c` (new permanent `tune.rs`
  infra, `"acos8"` arg) against the same joint acos+asin objective as the
  entry directly above, starting from the shipped 7 coefficients plus a
  prepended 0.0, converged with **that 8th coefficient at exactly 0.0** —
  the search found no use for the extra degree at all. The remaining 7
  coefficients it did move to match, bit-for-bit, the *already-rejected*
  unconstrained joint-objective result immediately above (same known
  4→5 acos regression). The extra degree of freedom buys nothing beyond
  what degree 6 already offers under this objective; not worth the
  latency cost of even testing in `src/lib.rs`. Not adopted.

- **atan2 division-residual correction (2026-07-07)**: added a first-order
  Taylor correction (`atan'(d)*e`) for atan2's `y/x` division rounding.
  First measurement looked like a 10x win (avg ulp 0.136→0.0134) until the
  max-ulp column showed a NaN sentinel — `d=inf` made the correction NaN.
  After guarding with `corr.is_finite()`, the "10x win" evaporated entirely
  (0.136 vs 0.1362, statistically identical) — atan's own poly-fit error
  already dominates atan2's total error. Real cost for zero benefit: latency
  +13.9%, throughput +42.9%. Reverted.

- **exp: replace the k1/k2 split with a k-clamp (2026-07-08)**: the
  framing ("split exists only because round can push k to 128, an
  out-of-contract case since exp is unchecked everywhere else") was wrong
  — k=128 is reachable from ordinary in-domain x (x*log2e in [127.5,128)
  is still inside the documented [-126,128) range), and the split isn't
  there to avoid inf for out-of-contract inputs, it's there to give the
  *correct* answer for this legitimate slice of the domain. Verified by a
  temporary exhaustive-ish (every 97th f32 bit pattern) old-vs-new
  comparison: first mismatch at x=88.376686 (comfortably inside the ~88.7
  ceiling), old=2.4071711e38 (correct, matches e^x), new=1.2035856e38 —
  exactly half, i.e. clamping k to 127 silently drops a whole factor of 2
  for real in-domain inputs near the top of the range. Never reached mca;
  killed by the fast falsification step the task description recommends
  (fuzz/scratch-test before the expensive latency/throughput cycle). Not
  adopted; reverted before touching mca.

## Codegen & build hygiene

- **Forcing zmm-width (AVX-512) codegen (2026-07-07)**: confirmed this CPU
  (Tiger Lake) has full AVX-512 available but LLVM deliberately emits ymm
  (256-bit) throughout — no `prefer-512-bit` knob exposed at the rustc level
  independent of `-C target-cpu`'s own tuning table. Client Intel parts are
  documented to downclock under sustained AVX-512, so this is very likely
  LLVM's informed choice, not an oversight. Not forced; the register-pressure
  question this was meant to investigate was never actually tested as a
  result.

## sin_checked / cos_checked internals

- **round_x_over_pi: remove sin_checked's dead `pre_offset=0.0` add
  (2026-07-07)**: monomorphized `round_x_over_pi` so sin_checked's
  instantiation skips the always-zero add. Confirmed via asm diff the
  instructions were gone. Latency unchanged (109 cyc — wasn't gating the
  critical path); throughput got *worse* (5.289→5.410 cyc/elem) despite
  fewer instructions — removing the op let the scheduler make different
  choices elsewhere that cost more than it saved. Reverted.

- **round_x_over_pi: switch qh (`p0.round()`) to round_ties_even
  (2026-07-06)**: ql's switch was safe and kept separately, but qh's
  regressed cos_checked's max ulp 2→6 in-domain (worst x≈252.9) — qh's rare
  exact-half ties interact badly with cos's -0.5 pre_offset in a way sin's
  zero offset never hits. qh reverted to `f32::round`.

- **reduce_pi: rebalance the 4-deep serial err chain to depth 2
  (2026-07-07)**: `(e1b+e2b)+(e3b-e3t)` instead of left-to-right.
  Analytically zero accuracy risk (confirmed bit-exact) but +3 cyc latency
  on both sin_checked and cos_checked — freeing part of the chain let the
  scheduler make a worse choice elsewhere. Reverted to the flat chain.

- **reduce_pi: downgrade e2's two_prod to a plain multiply (2026-07-07)**:
  unlike e3 (below), qh is unbounded, so there's no "provably zero"
  argument. Regressed sin_checked's *in-domain* (|x|≤1e6) max ulp badly:
  2→8 at |x|≤10, 2→144 at |x|≤1000, 2→51,054 at |x|≤1e6. Reverted before
  even reaching mca.

- **reduce_pi: downgrade e3's two_prod to a plain multiply (2026-07-07)**:
  the "e3 is provably zero when ql∈{-1,0,1}" premise checked out exactly,
  and gave the expected latency win (-4 cyc both functions). But: (1) mca's
  simulated scheduling diverged between callers — sin_checked throughput
  improved (-5.9%) while cos_checked's got *worse* (+12.3%), a
  caller-dependent regression not visible in the port-pressure-only
  estimate; (2) the already off-contract tail (|x|≥2^25) got dramatically
  worse (e.g. the [1e12,1e13) bucket: avg ulp 0.28→2.35M). Reverted.

- **parity() via integer bit-ops instead of mul/floor/fma (2026-07-07)**:
  read the LSB of qh/ql directly from their bit patterns instead of the
  float `parity()` formula. Verified bit-exact against the old formula on
  25M+ integer test values. Latency unchanged (parity was already off the
  critical path) but throughput got *worse* for both sin_checked (+2.5%)
  and cos_checked (+0.6%) — trading FP-port ops for integer ops didn't pay
  off once actually scheduled. Reverted.

- **sin_checked/cos_checked clamp: move the bound into `y.min()` inside the
  poly instead of the outer `r.clamp()` (2026-07-07)**: saves one op (2→1)
  by only bounding y=r². Fails: the raw residual `x` still enters
  `fma(p,x3,x)` unclamped and linearly, so `r≈-9.7e29` still overflows to
  -inf — reintroducing the exact "returns inf for finite input" bug the
  clamp exists to prevent. Not adopted.

- **Fast sin/cos: shorten the PI_A..D reduction chain from 4-deep to 3-deep
  (2026-07-07)**: computing `t=fma(q,PI_C,q*PI_D)` in parallel to the main
  chain. The error estimate ("floor roughly doubles/triples") was wrong by
  orders of magnitude once measured: avg ulp 0.0645→1.4818, max ulp
  220→866,390,494. Near sin's zeros, the correctly-reduced residual is
  tiny, so `t`'s new single-shot rounding becomes a huge *relative* error
  exactly where sin is most sensitive. Reverted before even running mca.

- **cos's `q = (kb - ROUND_MAGIC) + 0.5`, fold into one constant**: dies on
  representability (`ulp(1.5·2^23)=1`, so the folded constant doesn't exist
  as an f32). A half-magnitude magic quantizes q to halves (wrong). Folding
  +0.5 into the PI_A chain instead reintroduces a rounding near cos's
  zeros. Documented as investigated-and-probably-not; not attempted.

## Missed fma contractions

- **asin: `a*a - a` → `fma(a, a, -a)` (2026-07-07)**: same fusion pattern
  that won on hypot/asinh/acosh, bit-for-bit identical here too — but mca
  showed throughput getting *worse* (1.433→1.479 cyc/elem) with latency
  unchanged. Reverted; left as `a*a - a`.

## Other spots

- **atan2's `bothzero`/`hpisignx` boolean simplification (2026-07-07)**:
  `A || (¬A∧B) ≡ A∨B` simplifies away `bothzero`/`nonzeroy` at the source
  level. Full asm diff showed the compiled output is byte-for-byte
  identical — LLVM's InstCombine already does this simplification. No
  measurable change either way; not committed (nothing to gain, but
  recorded so it isn't retried).

- **erfc's final `fma(y, z, w)` → `mulsign(y, x) + w` (2026-07-07)**: trades
  an FMA-port multiply for a sign-xor + add, intended to relieve erfc's
  saturated FMA/mul ports. Asm confirmed the intended codegen shift
  happened, but mca showed throughput getting *worse* (2.097→2.158
  cyc/elem) instead of better. Reverted.

- **atanh via a single log1p call, `0.5*log1p(2x/(1-x))` instead of
  `0.5*(log1p(x)-log1p(-x))` (2026-07-08)**: algebraically equivalent (the
  identity checks out, and unlike the backlog's framing no `mulsign` turned
  out to be needed at all — the direct formula handles every edge, incl.
  signed zero, x=±1, and |x|>1, for free). Implemented and fuzzed: **real
  regression**, max ulp 3→31303 (avg only 0.031→0.046, so this is very much
  a tail-only blowup, worst x≈-0.999998). Root cause: the two-log1p form
  passes x and -x to log1p *completely unrounded* (they're the literal
  input, no arithmetic before the call), so whatever error exists is only
  log1p's own baseline error for that input. The one-log1p form instead
  computes u=2x/(1-x) via a division that rounds once — a tiny, ordinary
  rounding error — but log1p's derivative 1/(1+u) diverges as u→-1 (exactly
  atanh's own singularity, which u inherits), so that tiny upstream
  rounding error gets amplified by ~5 orders of magnitude before log1p even
  starts its own computation. Halving the op count traded away the
  "feed log1p an exact literal" property that was quietly doing a lot of
  work. Reverted immediately (before mca — the accuracy regression alone
  disqualifies it).

---

# Brainstorm backlog (2026-07-08) — UNTESTED

Kitchen-sink candidates, none implemented or measured. Every one must
survive the usual gauntlet before adoption: auto-vectorization check
(`--emit=asm`, no scalar fallbacks), exhaustive/fuzz accuracy.rs sweep,
edgecheck.rs, and mca latency+throughput (both axes — this file is full of
"latency won, throughput lost" reversals). Ideas already tried and rejected
above are deliberately absent. Hardware context that shapes several of
these: the FP divider is nearly idle in every measured hot loop while the
FMA/mul ports are the bottleneck, so "add a division to remove fmas" is a
legitimate direction here, unlike on most targets.

## Cross-cutting / methodology

- **Sollya `fpminimax` instead of lolremez + coordinate descent**: the
  coordinate-descent tuner has hit literal zero-move local optima on
  exp2/log_2 (see above), but fpminimax solves the coefficient-quantization
  problem *jointly* (lattice reduction over the f32 grid), which routinely
  beats round-then-tune, especially at higher degree. Candidates where max
  ulp is the open residual: acos_poly (max 4), atan_poly (max 18), erfc's
  n/d (max 109), exp's degree-5 (max 4-8 via callers). Would need
  `sollya` installed.

- **Simulated annealing / basin-hopping over coefficient space**: same
  motivation as above but no new tooling — perturb 2-3 coefficients at
  once (the coordinate-descent tuner only moves one axis at a time, so it
  can't cross diagonal valleys). Cheap to bolt onto tune.rs.

- **Exhaustive/rlibm-style correctly-rounded coefficient search for the
  smallest polys**: for a 4-coefficient poly over a bounded f32 domain, the
  set of coefficient vectors that round correctly at every domain point is
  an intersection of half-planes (linear in the coefficients) — an LP/
  interval search can find the *global* optimum rather than a local one.
  Realistic for sinf_poly (4 coeffs), cbrt's correction (4), expm1's Pade
  (5), atan_poly (4).

- **Even/odd poly sharing for ±r pairs**: `exp(x)` and `exp(-x)` share
  their reduction exactly (round is odd, so k(-x) = -k, r(-x) = -r), and a
  poly evaluated at both ±r splits into E(r²) ± r·O(r²) — one poly's worth
  of fmas produces both values. sinh/cosh/tanh currently pay for two full
  `exp` evaluations (or one exp + one division); this gets the second
  exponential for ~2 fmas + one exponent-field negation (2^-k is a bit
  trick). Biggest single-function win candidate in the file.

- **Rational (P/Q) refits of pure polys to exploit the idle divider**: a
  degree-(m/n) rational typically matches a degree-(m+n) poly's accuracy,
  so e.g. sinf_poly's 4 coeffs → 2/2 rational could cut fma count and
  Estrin depth at the cost of one division. Latency risk (division sits on
  the critical path, and unlike cbrt's rcp it can't start early), so this
  is a throughput idea, not a latency one. Candidates: log_2's degree-9
  (the deg-8 poly cut failed, but a 4/4 rational was never tried),
  acos_poly, erf_poly.

- **Batch/slice API tier (`exp2_slice(&[f32], &mut [f32])` etc.)**: the
  crate's whole perf story assumes the *caller's* loop auto-vectorizes;
  fixed-chunk slice entry points make that the crate's job instead
  (process in `[f32; 16]` chunks internally, same idiom as the bench
  harness), and open the door to explicit `core::simd` internals later
  without changing the scalar API. Also the natural place for a `sincos`
  (below).

- **Explicit `core::simd` internals as a fallback tier**: not to replace
  autovectorization, but for the few constructs autovectorization can't
  form at all (gathers for LUT variants, per-lane shuffles). Only worth it
  if a LUT idea below ever survives screening.

## log_2 / ln / log10

- **atanh-form reduction**: t = (m-1)/(m+1), log2(m) = (2/ln2)·atanh(t),
  an *odd* series in t — poly in t² needs ~5 coefficients where s = m-1
  needs 10 (t is bounded by ~0.172 vs s's ~0.414, and odd symmetry halves
  the terms). One division buys ~half the poly and one less Estrin level.
  Division-on-critical-path caveat as above; the classic tradeoff every
  libm makes the other way, but this CPU's idle divider may flip it.

- **log1p small-|x| dedicated poly branch**: log1p currently always pays a
  full ln poly + division; a Taylor/minimax branch for |x| < 0.25 (like
  asin_small) selected against the existing path could beat it on accuracy
  where log1p matters most, and is free perf-wise if it replaces work
  rather than adding a third arm. Feeds asinh/acosh/atanh accuracy too.

- **Shared denormal-rescale helper**: log_2/ln/log10 triplicate the
  tiny/xs/koff dance. Pure hygiene, no perf claim — only worth doing if
  touching these anyway.

## exp family

- **tanh: direct rational x·P(x²)/Q(x²)**: replaces expm1(2x) + division
  (which drags in the whole exp reduction+poly) with one even rational fit
  on [0, ~9.02] (tanh saturates to 1.0f32 past that; clamp handles the
  tail, same POLY_SAFE_BOUND pattern). One division, ~6-8 fmas total,
  odd symmetry via mulsign. Likely both faster *and* more accurate than
  the current route; the standard ML-workload tanh shape.

- **sinh/cosh via the shared-reduction even/odd trick** (see cross-cutting
  entry): sinh(x) = (2^k·(E+rO) − 2^-k·(E−rO))/2 needs care when the two
  scales differ hugely (for |x| > ~9 one side vanishes — which is also
  when cancellation is impossible, so it's benign), but eliminates an
  entire exp evaluation from sinh/cosh and both from tanh-via-expm1 if
  kept. sinh_throughput/cosh_throughput's 1/e division route becomes
  obsolete if this works.

## sin / cos / tan

- **`sincos` / `sincos_checked` returning both values from one
  reduction**: the rejected fused-tan attempts above failed on *poly*
  grounds (pole cancellation in the division); a plain sincos that runs
  the existing reduction once and sinf_poly twice (r and its ±pi/2
  counterpart... no — r with sin parity, and separately the cos offset
  needs its own q) — the cheap version shares only x/pi; the checked
  version shares the expensive two_prod machinery with two different
  pre_offsets, where round_x_over_pi's dominant two_prod(x, RPI_HI) and
  reduce_pi's qh-side products genuinely coincide when qh matches.
  Needs design work, but callers wanting both currently pay 2x the most
  expensive reduction in the crate.

- **tan via mod-pi/2 reduction + dedicated tan rational with reciprocal
  branch (sleef-style)**: reduce r to [-pi/4, pi/4] with q = round(x·2/pi),
  then tan(x) = tan(r) or -1/tan(r) by q's parity. This dodges the exact
  failure mode of both rejected fused-tan attempts (they evaluated near
  the pole; here the pole becomes a well-conditioned reciprocal of a
  near-zero-argument value). Costs: a new reduction constant set (pi/2
  split), a rational poly, one division, one select. Current tan max ulp
  is ~3000 near poles — this is the only idea on the list that could fix
  that.

- **cos (fast tier): dedicated even poly in r²**: q = round(x/pi) same as
  sin, cos(x) = ±cos(r) with cos poly = even, killing the -0.5/+0.5
  offset ops and the x³ multiply. Known risk (why sleef does't do it):
  near r = ±pi/2, cos's value → 0 while an absolute-error even fit stays
  O(2^-24), so *relative* error blows up near cos's zeros — exactly the
  failure the shifted-sin form avoids. Probably dies for that reason;
  listed so the reasoning is recorded rather than re-derived.

## cbrt

- **Rational correction, (1+r)^(-1/3) ≈ P(r)/Q(r) degree 1/1 or 2/1**: the
  degree-2 poly cut failed (see above) because 3 poly coefficients can't
  bend enough; a rational with 4 coefficients has a different (usually
  better) approximation class for algebraic functions like this one, and
  cbrt already spends one division (rcp) — mca showed 2 divisions/element
  still leaves the divider ~idle (cbrt_accurate measurement). Could reach
  degree-3-poly accuracy at 1 less fma depth, or beat its max ulp 2 at
  equal cost.

## asin / acos / atan

- **atan_poly degree bump (2/2 → 3/3 rational, or higher)**: max ulp 18 is
  the crate's worst in-budget-relevant offender after erfc; the 2026-07-07
  refit confirmed the *current form* is at its floor, so more degrees of
  freedom is the only way down. Each degree adds one fma per chain (numer
  and denom evaluate in parallel, so latency impact is one fma level per
  +1 degree on both). atan2 inherits any gain directly — and if atan gets
  good enough, the previously-useless division-residual correction
  (rejected above, atan error dominated) becomes worth re-testing.

## erf / erfc

- **erfc in log space**: fit log2(erfc(x)·2^(x²·log2e))'s rational part —
  i.e. fold the n/d rational *into* the exp2_checked exponent as an
  additive poly: erfc(x) = exp2(-x²·LOG2_E + R(x)). One exp2_checked call,
  no division, no separate rational, and the huge dynamic range (down to
  ~1e-45) lives where it's linear. erf's own tail already validates the
  exp2(poly) shape (max 5 there vs erfc's 109). Denormal outputs still
  round coarsely — check whether the 109 is actually *output-ulp-at-
  denormal* noise before crediting any fix.

- **erfc domain split**: if log-space fails, two rationals ([0,2] /
  [2,10]) with one select — both arms computed branchlessly, so ~2x the
  poly cost; only worth it if 109 → single digits.

- **erf: joint boundary+coefficients refit**: the 0.28 Pade/tail boundary
  was inherited, both branches were refit *separately* at fixed boundary
  (no headroom found, see above) — but moving the boundary itself while
  refitting both (the asin fix-5 lesson: measure where each branch
  actually degrades, don't trust the inherited threshold) was never done.

## hypot / misc

- **remainder_checked beyond 2^24**: double-float q (qh, ql) like
  sin_checked's reduction, with two correction candidates instead of one.
  Heavy; only worth it if a real use case needs |x/y| > 2^24.

- **powf_checked denormal-boundary fixup**: the remaining ~150 max ulp is
  exp2_checked's double-rounding into denormals; computing 2^(k+64)·p and
  multiplying by 2^-64 at the end (one extra multiply, only when the
  result is denormal-range) would make the final rounding single. Also
  fixes the same characteristic in exp2_checked itself if done there.

- **Small-poly Estrin audit (asin_small, sinh_small)**: both are 3-deep
  Horner in x²; Estrin gets them to 2-deep + one extra multiply. Latency-
  only candidates, and the asin fusion note above shows these sometimes
  measure backwards — cheap to test, low expected value.

- **exp10 / pown / rsqrt API additions**: exp10 = exp2(x·LOG2_10) with a
  Cody-Waite combine (same shape as ln's fix, in reverse); pown(x, i32)
  via repeated squaring for exact integer powers (vectorizes if the
  exponent is uniform); rsqrt = 1/sqrt (divider idle, likely just works
  at ~1 ulp). New surface, not optimization of existing functions.

## Backlog round 2 (2026-07-08) — also untested

### Cross-cutting / approximation theory

- **Automated evaluation-order search per poly**: the log shows fma
  reassociation is a per-poly coin flip (erf_poly Estrin was free,
  acos_poly's wasn't; exp2's regroup won). Enumerate all valid
  parenthesizations/groupings of each fixed coefficient set (there are
  only dozens for degree ≤ 9), score each with the exhaustive sweep + mca
  automatically, keep the Pareto set. Turns the recurring hand-tries into
  one script.

- **ulp-weighted minimax fits**: lolremez weights target relative error,
  but the actual objective is *ulp*, which is a staircase in the result's
  exponent — near an output power-of-two boundary, relative error and ulp
  disagree by up to 2x. A weight of 1/ulp(f(x)) (piecewise-constant) in
  the fit could buy back exactly the boundary cases that show up as
  max-ulp outliers.

- **Df32 leading coefficients**: for polys whose *leading* coefficient's
  own f32 rounding dominates the fit error (check: does the infinite-
  precision minimax beat the f32-quantized one by a lot?), store c0 (or
  the constant term) as a hi+lo pair applied with one extra fma. 1-op
  accuracy lever, cheaper than a whole Df32 evaluation. atan_poly's
  denominator constant and acos_poly's pi/2 term are candidates (pi/2
  is famously not well-representable in f32).

- **Centered-variable refits**: exp2's f lives in [0,1) — refit in
  g = f − 0.5 so coefficients shrink and their individual f32 roundings
  matter less (coefficient rounding error scales with coefficient
  magnitude). Same idea for erfc's xa ∈ [0,10] (hugely off-center today).
  Free at runtime (the shift folds into existing adds) but changes every
  coefficient, so it's a refit-and-sweep job.

- **Compensated-Horner accuracy tier**: run the poly with error-free
  transformations (two_prod/two_sum per step, like reduce_pi does for the
  reduction) to get a near-correctly-rounded result at ~2-3x cost.
  As `_checked`-style opt-in tiers for the max-ulp offenders (atan: 18,
  erfc: 109) rather than a default.

- **Worst-case patch lists**: functions whose exhaustive sweep leaves a
  literal handful of failing inputs could compare-select against the known
  bad bit patterns (1-2 vcmpps+blend if the misses share a mantissa or
  cluster). Fragile (any refit invalidates the list) and only sane where
  the count is tiny and stable; record which functions actually have
  concentrated misses first.

- **Unchecked tiers for the log family and others**: log_2/ln/log10 pay
  denormal-rescale + two special-case selects on every call; a
  `log_2_unchecked` (positive-normal-only contract, like exp2 vs
  exp2_checked) drops ~4 ops + 2 selects. Same argument for a
  `hypot`/`atan2` that skips inf special cases. The tier split is already
  this crate's established pattern; it just hasn't been applied uniformly.

- **FTZ/DAZ feature flag**: under a cargo feature declaring "caller runs
  with FTZ+DAZ on" (the common game/audio configuration), every denormal
  branch (log family's tiny rescale, cbrt's, exp2_checked's low clamp)
  becomes dead code. Free speed for users who already flush anyway;
  compile-time contract, no runtime cost for anyone else.

- **f32x16/AVX-512 via function-level `#[target_feature]`**: the rejected
  zmm attempt went through rustc-wide flags; an explicit `core::simd`
  f32x16 slice-API path under `#[target_feature(enable = "avx512f")]`
  sidesteps LLVM's tuning-table choice entirely. Downclocking concern
  stands, but it was never actually measured — a slice-API experiment
  would settle the register-pressure question the flags route couldn't.

- **Cross-check mca with uiCA and real counters**: llvm-mca's scheduling
  model has already produced caller-dependent surprises (e3 downgrade).
  uiCA is measurably more accurate for Tiger Lake, and `perf stat` on the
  quickbench loops validates either. Methodology, but it de-risks every
  other entry here.

### log family

- **Fuse ln/log10's final add into an fma**: `fma(p, s, k_hi) + k*LN2_LO`
  is a mul + add that won't contract on its own — `fma(k, LN2_LO,
  fma(p, s, k_hi))` is one op shorter and strictly one fewer rounding.
  Same in log10_normal. Small, likely-free, possibly measurable on
  everything downstream (asinh/acosh/atanh/log1p/powf all route through
  these).

- **Integer koff**: the tiny-branch offset is applied as a float add
  (`e as f32 + koff`); folding it into `e` as an integer subtract before
  the convert drops a float op from the k path. k is off the critical
  path, so this is a throughput-only micro-candidate.

### exp family

- **Select-tree "LUT" for exp2**: split f = f_hi + f_lo where f_hi takes 4
  (or 8) quantized values; 2^f_hi comes from a 2-level vblendvps tree over
  4-8 constants (no gather needed), and f_lo's range shrinks 4-8x so the
  poly drops from degree 5 to 2-3. Trades ~3 fma levels for ~2-3 blends +
  1 compare chain. Same trick applies to exp's e^r poly. This is the
  gather-free version of the LUT idea round 1 wrote off.

- **tanh via negative-argument expm1**: `tanh(x) = -expm1(-2|x|) /
  (expm1(-2|x|) + 2)` with mulsign. The argument is always ≤ 0, so expm1
  stays in [-1, 0] and *never overflows* — the current form inherits exp's
  unchecked-domain garbage for |x| > ~44, this one saturates to ±1
  naturally over the whole f32 range. Fixes a real (documented) domain
  hole at zero-ish cost (abs + mulsign vs. nothing).

- **Dedicated sinh/cosh kernels on the reduced argument**: instead of
  composing exp twice (or the round-1 even/odd trick), do the reduction
  once (k, r) and fit sinh(r)/cosh(r) minimax polys directly on
  [-ln2/2, ln2/2]: sinh(x) = sinh(r)·cosh_k + cosh(r)·sinh_k where
  cosh_k/sinh_k = (2^k ± 2^-k)/2 come from two exponent bit-tricks. The
  polys are even/odd so they share r² powers. Compare against round 1's
  shared-exp idea; one of the two should win.

- **exp's magic-constant round**: exp uses `.round()` (vroundps) while
  sin gets its round *and* its parity from one ROUND_MAGIC add. exp needs
  no parity, so vroundps is probably already optimal — but the magic-add
  version makes k1b available one op earlier (it *is* k+383 pre-shifted).
  Worth one mca look since exp is the hottest composite dependency
  (sinh/cosh/tanh/expm1 all inherit).

### sin / cos

- **mod-pi/2 reduction with paired even/odd polys**: reduce with
  q = round(x·2/pi) so |r| ≤ pi/4, then select sin(r)/cos(r) by quadrant.
  The poly degree drops hard on the halved range (deg-7 failed on
  ±pi/2 but has huge margin on ±pi/4; cos's even poly similar), and the
  two polys share r²/r⁴ powers, so "evaluate both + blend" costs much
  less than 2x. Costs: q's parity becomes 2 bits (two selects), and the
  reduction constants change (pi/2 splits, q twice as large so the
  fast tier's cliff moves down 2x — check that against the |x| ≤ 1e6
  contract before adopting). Applies to fast *and* checked tiers, and
  makes the round-1 sincos idea nearly free (both polys already
  evaluated).

- **Resurrect the f64-reduction sin/cos as a middle tier**: the abandoned
  f64 version (git history, ~2026-07-06) measured ~1.9/2.3 cyc/elem —
  roughly 3x faster than today's checked tier — while staying accurate to
  ~1e9 (vs. the fast tier's 1.3e7 cliff and checked's 1e13). Three-tier
  fast/medium/checked would cover the "audio/DSP phase accumulator"
  middle ground the current pair brackets awkwardly. The code already
  exists and was verified; it's a revert-and-rename.

- **Vectorized Payne-Hanek "exact" tier**: full-range correct reduction
  needs the 2/pi product against x's mantissa with the window selected by
  x's exponent — per-lane variable shifts exist (vpsrlvd, AVX2) and the
  2/pi table is small enough for a 4-8 constant select tree instead of a
  gather. Would make a `sin_exact` with no accuracy cliff anywhere in
  f32. Big job, listed for completeness (the "graceful degradation"
  contract makes it optional, not required).

### cbrt

- **Joint seed-constant + degree-3 coefficient search**: the rejected
  joint search was seed × *degree-2* (3 coeffs, hopeless). Seed × the
  shipped degree-3 was never searched — the seed constant currently in
  use was inherited, then the poly tuned around it. A basin-hopping pass
  over (seed offset, c1..c4) jointly targets avg 0.31 → lower at zero
  runtime cost.

- **Tune cbrt_throughput's magic constants**: the experimental 2-iteration
  Halley-ish form (5.5 avg ulp) has never seen the coordinate-descent
  tuner. Its 3 magic constants + 2 fma constants are all free parameters;
  even landing at ~2-3 avg would make it a legitimate documented tier
  instead of an experiment.

### atan / asin / acos

- **Pure-poly atan latency tier**: atan_poly's division sits on the
  critical path (unlike cbrt's early rcp, it can't start until the poly
  numerator resolves... actually until x2 does). An odd degree-13..17
  polynomial needs no division: worse throughput (more fmas on saturated
  ports), likely better latency (divider's ~11 cyc + dependency removed).
  Opposite tradeoff to everything else here, so it's a per-caller tier
  question, not a replacement.

- **Three-interval atan reduction**: on top of a.min(1/a), split [0,1] at
  tan(pi/8) using atan(x) = pi/8 + atan((x−t)/(1+t·x)) for the upper part
  — argument range shrinks to ±tan(pi/8) ≈ 0.414, poly degree drops
  ~2 levels. Costs one division (idle divider) + one select + one
  constant add. The standard next step every scalar libm takes; never
  tried here.

- **Retune asin's 0.25 crossover after any acos_poly change**: fix 5's
  lesson (measure both branches' actual curves, don't trust the inherited
  threshold) applies automatically after the round-1 degree-7 idea or any
  refit. Bookkeeping entry so the follow-through isn't forgotten.

### erf / erfc

- **erfc: exact exponent via two_prod**: `-(xa*xa)*LOG2_E` rounds twice
  before exp2_checked even starts, and each ulp of exponent error is
  ~0.69 ulp of result error — likely a real chunk of the 109. Compute
  x² and its fma residual (two_prod), fold `residual·LOG2_E` into
  exp2_checked's poly argument (or the round-1 log-space form's R(x)).
  Bounded cost: one fma + one add.

### hypot / atan2 / new surface

- **hypot_checked with branchless exponent rescale**: extract
  max(exp(x), exp(y)) with integer ops, scale both inputs by 2^-e via
  exponent-field subtraction (exact), sqrt, scale back by 2^(e/2)... e
  odd needs a sqrt2 factor — cleaner: scale by 2^-2⌊e/2⌋. All bit tricks
  + selects, fully vectorizable, kills the documented overflow/underflow
  tradeoff for callers who need std-grade hypot without std-grade scalar
  code.

- **sinpi/cospi (argument in half-turns)**: q = round(x) and r = x − q are
  *both exact* in plain f32 — the entire multi-word reduction apparatus
  evaporates, leaving round + subtract + poly. Full-range accurate (no
  cliff anywhere, even at f32::MAX), faster than the fast tier, and the
  poly is a trivial refit of sinf_poly onto [-1/2, 1/2] scaled by pi.
  Callers doing phase accumulation in turns (audio, FFT twiddles) get
  strictly better everything. Probably the highest value-per-effort new
  function possible in this crate. sind/cosd (degrees, mod 180 exact for
  the same reason) falls out of the same shape if wanted.

- **sigmoid/logistic**: 1/(1 + exp(-x)) — one exp tier + one division
  (idle divider), or expm1-based near 0 if the cancellation check demands
  it. Pairs with the tanh rational idea (logistic(x) =
  0.5 + 0.5·tanh(x/2), so whichever fit wins serves both).

- **exp10**: exp2(x·LOG2_10) has the same argument-rounding flaw exp2(x·
  LOG2_E) had before exp's Cody-Waite fix; do it properly from day one
  with a LOG10_2_HI/LO-style split (constants already exist for log10).

- **atan2_checked**: revisit the rejected division-residual correction
  *after* any atan_poly accuracy bump — the rejection reason was "atan's
  own error dominates", which stops being true if atan drops to ~2 ulp.
  Recorded dependency, not an independent idea.
