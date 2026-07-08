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

- **Coordinate-descent tuner methodology finding, "zero-move" isn't always
  "no headroom" (2026-07-08)**: `tune.rs`'s tuner moves each coefficient by
  integer *bit-pattern* steps (deltas of ±1..16 ulp). For a brand-new
  coefficient seeded at exactly `0.0` (bit pattern `0x0`), those steps only
  reach denormal-scale values, which have ~zero effect on the polynomial
  and essentially never score better — the search is *structurally unable*
  to leave 0, not evidence the extra degree of freedom is useless. Confirmed
  by direct comparison on `atan_poly`'s degree bump (see git history): a
  zero-seeded 3/3 rational reported "tuned max 18" (identical to the 2/2
  form, a textbook zero-move result) while a *real* least-squares Pade fit
  (scipy) of the exact same shape found max ulp 3 on the same grid — a
  genuine, large local optimum the zero-seed could never reach. An
  arbitrary nonzero seed (e.g. `1e-3`) isn't a fix either: it can start
  from a *worse* point than the existing form and the greedy per-
  coordinate search may never recover (measured: converged to max ulp 565,
  far worse than not bumping the degree at all). **When a "zero-move"
  result is reported for a genuinely new coefficient (not a re-tune of an
  existing one), don't trust it as "no headroom" without first trying a
  properly-computed nonzero starting point (scipy least_squares / a real
  Pade or minimax fit) as the tuner's seed.** This may retroactively call
  into question other zero-move entries in this file that used a bare
  `0.0` seed for a new coefficient (e.g. `acos_poly8`, `expm1_near0_deg5`)
  — not re-litigated here, but worth remembering if revisiting them.

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

- **expm1 Pade degree bump, numerator degree 3 → 5 (2026-07-08, later
  re-checked with a proper scipy seed instead of 0.0, still not
  adopted)**: originally added a new odd term (`expm1_near0_deg5_c` in
  `tune.rs`, `"expm1_near0"` arg), seeded at 0.0 — a zero-move local
  optimum, later understood (see this file's Cross-cutting tuner-
  methodology finding) to likely be the zero-seed trap rather than a
  genuine floor. Re-tested with a real
  scipy `least_squares` fit as the seed (found max abs error 3.4e-9 vs
  the shipped form's 5.3e-8, ~15x better in isolation) — this time real
  headroom *did* show up in avg ulp (fuzz + exhaustive: 0.1382→0.1345)
  but **not in max ulp** (exhaustive: 6→6 unchanged, worst case at
  x≈1.09, which is in `expm1`'s *other* branch — plain `exp(x)-1`
  for `|x|≥0.5` — entirely untouched by this change). So unlike `atan`'s
  degree bump, the extra numerator degree here improves the branch it
  targets but doesn't move the function's actual worst case at all, and
  it cost real throughput (mca 1.779→1.905 cyc/elem, +7.1% worse,
  cascading to `tanh` since it calls `expm1` internally). Reverted;
  `src/lib.rs` and `tune.rs` restored. If `expm1`'s max ulp 6 is ever
  worth chasing, the `b = exp(x)-1` branch (or `exp`'s own accuracy near
  x≈1) is where the actual headroom would need to come from, not this
  branch.

- **tanh: direct rational x·P(x²)/Q(x²) over the whole [0, ~9.02] domain
  (2026-07-08), doesn't converge at a practical degree**: the backlog
  framed this as "~6-8 fmas total," but a scipy `least_squares` fit of
  this exact shape needed **13 free coefficients** (P deg 13 / Q deg 12)
  to reach ~6.5e-8 max abs error (the scale needed for a few-ulp result)
  over the full domain — far more than any poly in this crate, the same
  wide-dynamic-range convergence problem already found for erfc's log-
  space idea. Tried splitting into two domains instead (matching erfc's
  own "domain split" fallback pattern): [0,4] converges beautifully with
  a 3/3 rational (7 coefficients, max abs err 1.7e-8), but [4,9.02] (the
  near-saturation tail) needs its own 3/3 (another 7 coefficients, max
  abs err 4.4e-7, borderline) to get there — 14 coefficients total across
  two branches plus a select, all evaluated unconditionally per this
  crate's branchless convention. That's likely *more* total work than
  the current `expm1(2x)/(expm1(2x)+2)` route, especially now that `exp`
  itself is much faster after this session's magic-round fix — the
  backlog's assumed win doesn't obviously hold once `exp`'s own cost
  dropped. Not implemented; no code changed.

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

- **Integer koff: fold the tiny-branch exponent offset into `e` as an
  integer add before the int->float convert (2026-07-08)**: changed
  `let k = e as f32 + koff;` (koff: f32, 0.0 or -24.0) to `let k = (e +
  koff) as f32;` (koff: i32, 0 or -24) across `log_2_normal`/`ln_normal`/
  `log10_normal`/`log2_df` (a public signature change, `koff: f32` ->
  `i32`, updated at the one `mca_target.rs` call site too). Verified
  bit-exact (expected -- both forms compute the same exact-integer sum,
  just in a different order) and mca showed **zero measurable change on
  any axis, to the reported decimal** (log2/ln/log10/powf all identical
  before and after) -- LLVM already performs this exact reordering
  itself regardless of the source-level int-then-convert vs. convert-
  then-add ordering, so the "micro-optimization" was already happening
  for free. Reverted rather than keep a public API signature change
  (`f32`->`i32`) for literally zero benefit.

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
  **Re-checked 2026-07-08 (later, same day) with a proper scipy seed**
  instead of 0.0, per this file's own Cross-cutting tuner-methodology
  finding (the zero-seed trap that also affected `atan_poly`'s degree
  bump) — this time real, if modest, headroom showed up: implemented
  directly in `src/lib.rs` and measured against the real crate (not just
  `tune.rs`'s grid), `acos` avg/max ulp 0.4962/4 → 0.4904/3 (genuine
  improvement, *not* the previously-seen regression) and `asin` avg
  improved slightly with max ulp unchanged at 9. But unlike `atan_poly`'s
  dramatic 18→4 cut, this is a 1-ulp acos improvement with no max-ulp
  movement at all for `asin` — for a similar-sized real mca cost (asin
  59.03/0.968 → 63.03/1.039 cyc lat/throughput, +6.8%/+7.3%; acos
  37.11/0.820 → 41.11/0.862, +10.8%/+5.1%). Judged not worth it at this
  magnitude (unlike `tanh`'s domain-hole fix or `atan`'s large cut, there
  isn't a strong enough gain to justify the cost here). Reverted;
  `src/lib.rs` untouched, `tune.rs`'s scratch scipy-seed addition also
  reverted.

- **acos_poly: Df32 leading (constant) term, pi/2 split into an exact
  hi+lo pair (2026-07-08)**: `fma(u, x, 1.5707963)` (single f32-rounded
  pi/2) replaced with `fma(u, x, PI_2_HI) + PI_2_LO` (same trick as
  `LN2_HI`/`LN2_LO`, `PI_2_HI`/`PI_2_LO` together capturing pi/2 to
  ~1e-12 instead of f32's own ~1e-7). With the *same* 6 leading
  coefficients unchanged: a real but genuinely mixed result -- `acos`
  avg ulp improved dramatically (0.4961→0.0675, ~7.3x) but max ulp got
  *worse* (4→5), and `asin` (which reuses this poly) got worse on
  *both* (avg 0.0303→0.0376, max 9→12). Retuning the 6 leading
  coefficients for the new split-constant structure (`tune.rs`'s
  `acos_poly_df_c`, coordinate-descended from the shipped values)
  recovered `asin` back to baseline exactly (avg/max statistically
  unchanged) and kept `acos`'s avg win, but `acos`'s own max ulp still
  regressed, now 4→6 (exhaustive-confirmed both fuzz and thorough sweep
  agree). This is exactly the tradeoff shape the *first* `acos_poly`
  entry in this file already tested and rejected (an unconstrained
  joint objective that improves the joint score by letting acos's own
  protected max ulp regress) — `acos`'s own accuracy has been treated as
  a protected invariant in this codebase since that entry, not something
  to trade away for a joint or single-caller average-ulp win. Also a
  real mca cost on both functions for the extra `+PI_2_LO` add: asin
  59.03/0.968→63.03/1.039 cyc lat/throughput (+6.8%/+7.3%), acos
  37.11/0.820→41.11/0.858 (+10.8%/+4.6%) — and `asin` pays that cost for
  *zero* net accuracy benefit once retuned back to baseline. Not
  adopted; reverted (`src/lib.rs` and `tune.rs`'s scratch addition both
  restored).

- **atan2 division-residual correction (2026-07-07, re-tested 2026-07-08
  after atan_poly's degree bump, same conclusion holds)**: added a
  first-order Taylor correction (`atan'(d)*e`) for atan2's `y/x` division
  rounding. First measurement looked like a 10x win (avg ulp 0.136→0.0134)
  until the max-ulp column showed a NaN sentinel — `d=inf` made the
  correction NaN. After guarding with `corr.is_finite()`, the "10x win"
  evaporated entirely (0.136 vs 0.1362, statistically identical) —
  atan's own poly-fit error already dominates atan2's total error. Real
  cost for zero benefit: latency +13.9%, throughput +42.9%. Reverted.
  Re-tested after `atan_poly` dropped from max ulp 18 to 4 (2026-07-08,
  same session as the degree bump) on the theory that "atan's own error
  dominates" might no longer hold at the lower error level — it still
  does: avg/max ulp unchanged (0.0684/3 → 0.0693/3, noise-level), and the
  cost was almost identically bad (latency +13.1%, throughput +41.8%,
  nearly the exact same percentages as the original 2026-07-07 test).
  Even atan's much-improved ~4 ulp residual still swamps a sub-ulp
  division-rounding correction. Reverted again; this dependency is now
  closed for good barring a much larger atan accuracy improvement.

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

- **Small-poly Estrin audit, asin_small/sinh_small 3-deep Horner → 2-deep
  Estrin (2026-07-08)**: `let lo = fma(c1,x2,c0); let hi = fma(c3,x2,c2);
  fma(hi, x4, lo)` instead of the nested Horner chain. Unlike the fma-
  contraction idea above, this one isn't bit-exact (fma reassociation
  changes rounding, as it does roughly half the time elsewhere in this
  file) — and it measured backwards on every axis that matters: `sinh`
  got *worse* on both mca latency (58.00→59.00 cyc) and throughput
  (2.523→2.588 cyc/elem), plus a real accuracy cost (avg ulp 0.0805→
  0.0843, `cosh` itself unaffected since it doesn't call `sinh_small`).
  `asin` was a genuine mixed result — latency improved (59.03→56.00 cyc,
  -5.1%) but throughput got worse (0.968→1.033 cyc/elem, +6.7%) *and*
  max ulp regressed (9→10, `acos` unaffected, doesn't call `asin_small`).
  Neither survives on net. Reverted; exactly the "sometimes measures
  backwards" outcome the backlog itself predicted for this idea.

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

- **erfc: exact exponent via two_prod (2026-07-08)**: computed `xa*xa`'s
  exact rounding residual (`e = fma(xa,xa,-p)`) and folded `e*LOG2_E` into
  the exponent (`fma(-e, LOG2_E, -p*LOG2_E)`) instead of dropping it, as
  the backlog proposed. Real, exhaustively-confirmed accuracy win (avg/max
  ulp 0.3106/109 → 0.3055/93), but a real mca cost too (latency
  78.09→82.14 cyc, +5.2%; throughput 2.599→2.788 cyc/elem, +7.3% worse) —
  unlike `tanh`'s domain-hole fix (a genuine NaN-for-legitimate-input bug),
  this is just trimming an already-accepted, already-far-over-any-nominal-
  budget max ulp a bit further (109→93 is still nowhere near the top-of-
  file "≤2 max" budget this function was never going to hit anyway), so
  the cost isn't clearly justified by the gain. Not adopted; reverted.

- **erfc in log space, single polynomial over the full [0,10] clamped
  domain (2026-07-08)**: fit `R(x) = log2(erfc(x)) + x²·log2(e)` via
  lolremez to replace the n/d rational + division entirely (`erfc(x) =
  exp2_checked(-x²·log2(e) + R(x))`, same shape `erf`'s own tail branch
  already uses for `erf`, though `erf_poly` itself turned out unusable
  here directly — see below). Convergence was poor: degree 10 only
  reached estimated max error 6.5e-6 (need roughly 1e-7 for a few-ulp
  result), degree 15 (16 coefficients — far more than any poly in this
  crate) got to 3.75e-7, still short, and degrees beyond that took over
  2 minutes without converging. The backlog's own fallback ("erfc domain
  split... two rationals") is likely necessary for this shape to work at
  a practical degree; not attempted (bigger scope than fits one pass).
  Separately confirmed `erf_poly` (already log2(erfc)-shaped internally
  for `erf`'s own tail branch) can't just be reused for `erfc` directly:
  swapping it in gave catastrophic error (avg ulp ~297000, max ~3.6e8) —
  `erf_poly` was fit against *erf's* accuracy objective, where erfc's
  absolute tininess for large x barely moves erf's own ulp count (erf is
  already ~1 there), so it was never actually accurate as a direct erfc
  approximation, just close enough in erf's shadow. Not adopted.

- **powf_checked's ~150 max ulp investigated (2026-07-08): the backlog's
  diagnosis ("exp2_checked's double-rounding into denormals") doesn't
  survive direct measurement — the real residual lives somewhere in the
  `log2_df`/`exp2_checked_df` double-float argument chain, not in
  `exp2_checked`'s own rounding, and isn't a simple fix.** First checked
  the stated premise directly: a domain-restricted sweep of
  `exp2_checked` alone, restricted to `x` in `[-151,-120]` (its actual
  denormal-output zone), still measured max ulp **1** — i.e.
  `exp2_checked` itself is *not* measurably broken near the denormal
  boundary, contradicting the "double-rounding into denormals" theory at
  face value. Next, found `powf_checked`'s actual worst case via a
  targeted 20M-sample search (correcting a reference-computation bug
  along the way — `x.abs().powf(y)` alone doesn't reproduce the correct
  NaN for negative-base/non-integer-exponent, needed an explicit
  integer/parity check matching the crate's own convention): the worst
  case found was `x≈0.895`, `y≈-789.4`, landing near the exponent range's
  *upper* boundary (~125.7, close to `exp2_checked`'s +128 ceiling), not
  the lower/denormal one at all. Checked `exp2_checked` directly at that
  exact bit-pattern argument (not a re-derived f64 approximation, which
  gave a misleadingly different value the first time) — it was
  bit-exact. So the ~150 max ulp residual isn't `exp2_checked`
  mis-rounding its argument; it's more likely a precision limit in how
  `log2_df(ax) * y` (the `Df32` multiply) or `exp2_checked_df`'s own
  reconstruction handles this specific regime. Not root-caused further
  this session (would need tracing through the `Df32` arithmetic
  step-by-step at this exact input, a bigger investigation than fits
  here) — the backlog's "one extra multiply, only when denormal" framing
  is not the right fix given where the actual worst case lives. No code
  changed.

- **sigmoid/logistic, implemented (2026-07-08)**: the backlog's own
  suggested identity, `0.5 + 0.5·tanh(x/2)`, was tried first and is a real
  *bug*, not just imprecise — for x≈-17.33, `tanh(x/2)` (≈-8.66) already
  correctly rounds to exactly `-1.0f32` (true value within half a ulp of
  `-1.0`, so that's tanh's own correct output), but `0.5 + 0.5·(-1.0)` then
  computes to exactly `0.0` even though the true sigmoid value there
  (≈2.98e-8) is nowhere near it or f32's underflow threshold. Same bug
  class as the rejected atanh single-log1p-fusion above: composing through
  an intermediate function that has already saturated/rounded to an
  extreme value discards precision the outer formula still needs, even
  though the identity is algebraically exact. Fixed by computing directly
  instead — `1.0/(1.0+exp((-x).clamp(-87.0, 88.0)))` — which has no
  cancellation anywhere and gracefully saturates to 0.0/1.0 over the whole
  domain. Fuzzed clean after the fix (avg/max ulp 0.0996/4) and confirmed
  exhaustively (0.0996/5 over all 2^32 patterns). mca: 65.09 cyc latency,
  2.713 cyc/elem throughput (cheaper than tanh's 94.91/2.567, as expected
  for one exp + one division vs. a rational approx).

- **ln/log10: fuse the final `+ k*LN2_LO`/`+ k*LOG10_2_LO` into an fma —
  re-tested 2026-07-08, mistakenly adopted, then corrected back to
  rejected the same day.** This is the exact same idea the entry just
  above (commit-dated earlier the same day) already tested and rejected —
  picked again later in the session without cross-referencing the
  existing entry, and this time mis-measured as a win. Two errors, both
  now corrected: (1) the "before" accuracy baseline was taken from
  `readme.md`'s numbers (ln avg 0.126, log10 avg 0.286) instead of
  re-measuring the actual current unfused code directly — a fresh
  exhaustive sweep of the genuinely-unfused code gives ln avg 0.1168/max 3
  and log10 avg 0.1270/max 3, **identical to the fused code's own
  numbers** (`readme.md`'s log10 figure had simply gone stale at some
  earlier point, unrelated to this change — there is no accuracy
  difference between the two forms at all). (2) mca's ln/log10 latency
  going from 55.86→56.86 cyc was dismissed as "run-to-run noise" without
  checking reproducibility; re-run 3x on each form just now, it's
  perfectly deterministic both ways (unfused always 55.86, fused always
  56.86) — a real, reproducible +1-cycle regression, exactly matching
  what the entry above already found. Reverted for real this time;
  `src/lib.rs` restored to `fma(p, s, k_hi) + k * LN2_LO`. **General
  lesson, two of them: (a) before re-testing an idea, grep this file for
  whether it's already been tried — a duplicate test is wasted effort at
  best and, as happened here, a chance to overwrite a correct prior
  finding with a wrong one; (b) "before" measurements must come from
  re-running the current code, never from a written number (a readme
  table, a prior IDEAS.md entry, this file's own history) — those can go
  stale, and trusting one instead of re-measuring is exactly how this
  mistake happened. Also: any mca difference, however small, needs a
  repeat-and-confirm before being called "noise" — this file has now
  logged the opposite mistake too (calling a real effect noise) alongside
  its many entries logging noise mistaken for a real effect.**

- **log_2_unchecked/ln_unchecked/log10_unchecked, implemented (2026-07-08)**:
  the log family's "checked" public functions (`log_2`, `ln`, `log10`) pay
  a denormal-rescale multiply + two post-hoc selects on every call, even
  though the underlying `*_normal` cores (already existing, previously
  `#[doc(hidden)]`) don't need any of that for positive normal finite
  input. Exposed each core directly as a new public one-argument function
  (`log_2_unchecked(x) = log_2_normal(x, 0.0)`, etc. — the `koff`
  denormal-offset parameter is always exactly 0.0 for this contract, so
  it's fully applied rather than exposed), same fast/full-safety split as
  `exp2`/`exp2_checked`, just with the safe name already taken by the
  no-suffix function, hence the `_unchecked` suffix instead of a bare
  `log_2`/`ln`/`log10` swap. Verified bit-for-bit identical to the checked
  versions over 99.2M positive-normal samples directly (not just
  aggregate-ulp-similar) — expected, since it's the literal same core
  function called with the same always-correct koff, so there is zero
  accuracy risk by construction; exhaustive sweep confirms matching max
  ulp on both sides (log_2 3→3, ln 3→3, log10 3→3; the *avg* ulp numbers
  differ because the checked functions' "everywhere" sweep average also
  includes free/exact special-cased inputs like negatives→NaN that the
  domain-restricted unchecked sweep doesn't get to include — not a real
  accuracy difference, confirmed by the direct bit-for-bit check). Real,
  substantial speed win: mca throughput log2 1.556→0.958 cyc/elem (-38%),
  ln/log10 1.714→1.145 (-33%); ln/log10 latency 55.86→38.39 cyc (-31%,
  log2's own latency number was already measuring the unchecked core per
  this file's own mca_target.rs convention, so no further latency win
  there specifically). Confirmed on real wall-clock too via quickbench
  (e.g. log10 throughput 0.383→0.269 ns, latency 10.13→8.87 ns). Backlog
  entry for the remaining half of this idea (hypot/atan2 unchecked tiers)
  left open above. (Numbers corrected 2026-07-08, same day: originally
  quoted against the ln/log10-fma-fusion entry's numbers, since that
  change was active when this one was first measured — since reverted as
  a mistaken adoption, see that entry's own correction above.)

- **Centered-variable refit for exp2's f (2026-07-08), checked via a
  10-line scipy script before writing any Rust — premise refuted before
  implementation, nothing to revert.** The backlog's claim was "refit in
  g = f − 0.5 so coefficients shrink" (coefficient rounding error scales
  with coefficient magnitude, an established lesson elsewhere in this
  file). Fit both forms with scipy `lstsq` (uncentered `R(f) = (2^f-1)/f`
  over f ∈ [0,1) vs. centered `R(g) = (2^(g+0.5)-1)/(g+0.5)` over
  g ∈ [-0.5,0.5)) and compared every coefficient directly: the centered
  fit's coefficients are *larger* across the board (leading term
  0.693→0.828, and every other coefficient likewise bigger), the opposite
  of the claimed effect. Root cause: `R(f)` is strictly monotonically
  increasing over the whole domain (checked explicitly, 0.693 at f=0 up
  to 1.0 at f=1) with no interior minimum — its smallest magnitude
  already sits exactly at the domain's own edge (f=0), which the
  *current, uncentered* fit already exploits directly (c0 = R(0) = ln2).
  "Centering" necessarily moves the evaluation point away from that edge
  minimum toward the domain's middle, which is *higher*, not lower, for a
  monotonic function. The general principle (centering shrinks
  coefficients) only holds when the function has a genuine interior
  minimum/root to center on — `log_2`'s own `s = m-1` decomposition works
  for exactly this reason (`log2(1) = 0` is a real root at the center),
  but `exp2`'s `R(f)` has no analogous root anywhere in its domain. Not
  implemented; no code changed. The erfc half of this same backlog entry
  is a different function shape (P/Q rational, not a monomial poly) and
  is **not** ruled out by this finding — left open below, split out from
  the exp2 case.

- **Tune cbrt_throughput's magic constants (2026-07-08), tested and
  rejected on two different grids — a real methodological finding about
  `tune.rs` itself, not just this one function.** Added `cbrt_throughput_c`
  and a `cbrtthroughput` tune.rs target (kept as reference infra, same
  precedent as `cbrt_shiftmul_c`). First attempt reused `cbrt_normal`'s
  own single-octave-is-representative grid — looked like a win on-grid
  (max ulp 16→12) but a direct octave-by-octave accuracy check showed
  this function's error does *not* repeat across octaves the way
  `cbrt_normal`'s does (max ulp ranges from ~7 to ~52 depending which
  octave gets sampled) — and implementing the "tuned" constants for real
  confirmed the grid was misleading: full fuzz sweep got *worse*, not
  better (avg 6.73→7.01, max 74→76). Reverted immediately, tried again
  with a properly wide ~60-octave grid instead — this exposed a second,
  more general problem: `score()` returns `(max, sum)` and `tune()`'s `if
  s < best` uses Rust's default tuple ordering, which compares `max`
  *first* and only falls back to `sum` (~avg) as a tiebreak. For
  `cbrt_normal`'s smooth, octave-periodic error surface this never
  mattered (minimizing max there also happens to minimize avg), but for
  `cbrt_throughput`'s rougher landscape the search happily wrecked the
  average to shave the worst case: max ulp did drop (68→43 on-grid) but
  avg ulp *exploded* (6.83→22.73) — a real, severe regression by this
  crate's own primary metric, hidden by a max-first comparison that never
  surfaces it. Not adopted either way; `src/lib.rs` unchanged. **General
  lesson for any future `tune.rs` use, beyond just cbrt_throughput: (1)
  before trusting a single-octave (or any partial-domain) grid as
  "representative," check directly whether the target function's error
  actually repeats across the domain the way the assumption requires —
  it's a real, checkable property, not something to inherit from a
  different function's own comment; (2) `tune()`'s max-first tuple
  comparison is a silent trap for any function whose error surface isn't
  smooth/uniform — it can trade away average accuracy for a better
  worst-case number without ever showing that tradeoff in the printed
  output, so a real end-to-end fuzz check of any "tuned" result (not just
  trusting tune.rs's own reported numbers) is not optional, it's load-
  bearing.**

- **log_2: atanh-form reduction, `t=(m-1)/(m+1)`, `log2(m)=c·t·Q(t²)`
  (2026-07-08), implemented and measured — mathematically a huge win,
  practically a clear loss on both axes.** Confirmed the backlog's core
  claim first with a quick scipy fit: a degree-4 `Q` (5 coefficients) over
  `t∈[-0.172,0.172]` (odd symmetry, half the domain of `s=m-1`'s ~0.414)
  hits max relative error 1.2e-11 — about 1000x tighter than the shipped
  degree-9/10-coefficient `s·P(s)` form's 4.1e-9, with half the
  coefficients. Implemented directly in `log_2_normal` and measured for
  real: accuracy actually got slightly *worse* (fuzz avg ulp
  0.0030→0.0053, max 3→5) — `log_2` was already deep in pure f32-rounding-
  noise territory (avg 0.003 ulp), far past where a tighter *mathematical*
  fit moves the *measured* result. Latency got substantially worse, not
  better: mca 34.23→49.05 cyc (+43%), with throughput barely moving
  (1.556→1.521, only ~2%). Root cause: the one division this form needs
  (`t = s/(m+1)`) depends on `s`, available from the very first step of
  the critical path — unlike cbrt's early-starting `rcp` (independent of
  the seed/correction chain, so its latency hides behind other work),
  there's nothing here for the division to overlap with, so its latency
  (~11 cyc per this file's other division-cost notes) plus the reduced-
  but-still-serial poly evaluation came out slower overall than the
  original's longer, division-free chain. Reverted; `src/lib.rs`
  unchanged (kept `log2_atanh_c` in `tune.rs` as reference infra, same
  precedent as `cbrt_shiftmul_c`/`cbrt_throughput_c`). **General lesson,
  same shape as the exp2-centering and cbrt_throughput-tuning entries
  above: a dramatically tighter *mathematical* approximation doesn't
  automatically translate into a better *measured* result once a function
  is already accurate enough that f32 rounding, not polynomial degree, is
  the binding constraint — and "one division buys half the poly" only
  pays off if the division can start early enough to overlap with other
  work; a division gated on the very first reduction step never gets that
  chance.**

- **erf: joint boundary+coefficients refit (2026-07-08), screened cheaply
  before investing in a full joint optimizer — no headroom, matching
  what the separate fixed-boundary refits had already found.** The
  backlog's premise: the inherited 0.28 Pade/tail crossover was never
  swept as a free parameter, only refit-at-fixed-boundary (per this
  file's own asin-crossover precedent, moving the threshold itself can
  sometimes unlock headroom a fixed-boundary refit can't reach). Added a
  temporary boundary-parameterized `erf` variant and swept 8 candidate
  thresholds (0.20 through 0.40) against the *current, unrefit*
  coefficients first, cheaply, via accuracy.rs's existing fuzz
  infrastructure (100M samples each) — before spending time building a
  full joint boundary+coefficient optimizer. Result: max ulp sits at a
  flat 5 across the whole 0.26–0.32 range (current 0.28 already
  comfortably inside it), only degrading outside that window (0.24→8,
  0.20→52, 0.35→8, 0.40→24) — a wide, flat plateau, not a narrow optimum
  the inherited value happens to miss. Combined with the backlog's own
  already-noted finding (both branches refit *separately* at 0.28 found
  no headroom), two independent pieces of evidence now agree this
  function is already near its accuracy ceiling for this degree-6-tail +
  Pade-near-zero architecture, regardless of exactly where the boundary
  sits. Didn't build the full joint optimizer (existing `tune.rs` infra
  for `erf_tail_c`/`erf_near0_c` already hardcodes the 0.28 boundary into
  each grid's construction, so a true joint search would need new
  infrastructure) — the cheap screen already answers the question with
  reasonable confidence, and the effort/expected-payoff ratio for
  building the fuller version doesn't look favorable given both signals
  point the same way. Not implemented; `src/lib.rs`/`accuracy.rs` scratch
  additions reverted, nothing shipped. **General lesson: when a
  refit-oriented idea has an inexpensive proxy check available (here,
  sweeping the free parameter alone against unrefit coefficients, using
  infrastructure that already exists), run that first — it can settle
  the question well enough to skip building a bigger joint optimizer
  entirely, the same way this file's scipy pre-checks have repeatedly
  settled centered-refit and rational-form questions before any Rust
  was written.**

- **"Resurrect the f64-reduction sin/cos as a middle tier" (2026-07-08):
  the premise didn't survive contact with either git history or a real
  implementation — rejected on two independent grounds.** First: the
  backlog's own claim ("the abandoned f64 version... the code already
  exists and was verified; it's a revert-and-rename") doesn't hold up —
  a full search of `git log --all` and the reflog found no commit, stash,
  or working-tree file anywhere containing a genuine hardware-`f64`-based
  sin/cos reduction; the only "double" reduction in this repo's history
  is `sin_checked`'s own Df32 (emulated double-float, already shipped),
  which the backlog itself distinguishes from what it's proposing. There
  was nothing to revert. Implemented one from scratch anyway to check the
  underlying idea on its own merits (`sin_mid`/`cos_mid`: cast to `f64`,
  `q = round(x*FRAC_1_PI)`, `r = (x - q*PI) as f32`, reusing `sinf_poly`
  directly) — and it fails on accuracy by a wide margin, the second,
  independent rejection: even in the smallest bucket tested (`|x|<1e6`)
  max ulp was 637 (avg a reasonable-looking 0.0356, but that average
  hides the outliers near sin's own zeros, where a tiny absolute
  reduction error becomes a huge *relative*/ulp error), degrading to
  millions of ulp by `1e10` and low billions by `1e11`-`1e12` — nowhere
  close to the claimed "~1e9". Added an f64 two-product compensation for
  `q*PI`'s own rounding (`e = q.mul_add(PI, -q*PI)`, correcting the
  *multiplication's* rounding error) before giving up: this helped (max
  ulp 637→70 in the smallest bucket) but didn't come close to closing the
  gap, because the dominant remaining error isn't from the multiply's
  rounding at all — it's that `f64`'s `PI` constant itself is only one
  ~52-bit word, off from true π by its own fixed ~2^-53 relative error,
  which no amount of *compensating the arithmetic around it* can correct.
  Reaching real accuracy would need a genuine hi+lo (Cody-Waite-style)
  split of π across *two* f64 words — at which point the design is no
  longer a simple, cheap single-f64 reduction, and starts converging on
  something not obviously cheaper than `sin_checked`'s existing Df32
  approach (which already does exactly this kind of split, just in
  native f32 arithmetic instead). Not implemented; `src/lib.rs`/
  `accuracy.rs` scratch additions reverted, nothing shipped. **General
  lesson: "the code already exists, it's a revert" is itself a claim
  worth checking (`git log --all`/reflog/stash) before assuming the
  implementation work is already done — here it wasn't, and a from-
  scratch build was the only way to find out the accuracy claim was also
  wrong. Separately: a single hardware `f64` word is not a free
  drop-in replacement for a proper double-double reduction — it has
  exactly one rounding's worth of headroom over `f32` (roughly 2^-52 vs
  2^-24, not the "infinite precision" intuition might suggest), and for
  an argument-reduction problem specifically (subtracting a huge multiple
  of an irrational constant from a huge value to get a tiny, sign-
  sensitive residual), that headroom runs out far sooner than expected.**

- **cbrt: joint seed-constant + degree-3 coefficient search (2026-07-08),
  double-checked via two independent methods, no headroom found.** The
  previously-rejected joint search only tried a *degree-2* correction
  (hopeless regardless of seed); this one keeps the shipped degree-3 (4
  coefficients) and adds the seed offset as a 5th free parameter. Added
  `cbrt_normal_joint_c` to `tune.rs` (seeded from the shipped values, not
  zero) and ran the existing coordinate descent on the same grid
  `cbrt_normal` already uses: essentially no movement (avg 0.33871→0.33856,
  ~0.04% relative, max ulp unchanged at 2). Since coordinate descent alone
  can miss a genuinely different basin the backlog's "basin-hopping"
  framing was reaching for, double-checked independently with a
  from-scratch Python sweep of the underlying continuous math (not
  tune.rs) across a wide range of seed offsets (±2000, well beyond
  coordinate descent's ±16-per-step reach) with a fresh least-squares
  refit of the poly at each: best found was 1.407e-7 max relative error
  vs. the shipped combination's own refit at 1.411e-7 — also ~0.3%,
  negligible, and only in the underlying math (before any f32-rounding
  effects that would likely wash out even that). Two independent methods
  now agree the shipped seed+poly combination is already essentially
  optimal for this architecture; not adopted. Not implemented in
  `src/lib.rs`; kept `cbrt_normal_joint_c` in `tune.rs` as reference infra
  (same precedent as `cbrt_shiftmul_c`/`cbrt_throughput_c`).

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

- **log1p small-|x| dedicated poly branch**: log1p currently always pays a
  full ln poly + division; a Taylor/minimax branch for |x| < 0.25 (like
  asin_small) selected against the existing path could beat it on accuracy
  where log1p matters most, and is free perf-wise if it replaces work
  rather than adding a third arm. Feeds asinh/acosh/atanh accuracy too.

- **Shared denormal-rescale helper**: log_2/ln/log10 triplicate the
  tiny/xs/koff dance. Pure hygiene, no perf claim — only worth doing if
  touching these anyway.

## exp family

- **tanh via the shared-reduction even/odd trick**: tanh currently reuses
  `expm1(2x)` (already a single exp evaluation, not two), so the sinh/cosh
  version of this idea (now adopted, see git history/src/lib.rs's
  `exp_pos_neg`) doesn't directly apply here.

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

- **tan via mod-pi/2 reduction + dedicated tan poly with reciprocal branch
  (2026-07-08), implemented and measured — premise was wrong, real
  regression on every axis, reverted**: the backlog framed this as fixing
  "max ulp ~3000 near poles" via `sin(x)/cos(x)` allegedly amplifying
  error when dividing by a near-zero `cos(x)`. Implemented in full (new
  `PI_HALF_A..D` Cody-Waite split, `q=round(x*2/pi)`, `r` in `[-pi/4,
  pi/4]`, a dedicated degree-6 `tanf_poly(r)` fitted with scipy, `-1/t`
  reciprocal branch by q's parity) — and direct spot-checks at exact odd
  multiples of pi/2 (`k*pi/2` for `k` up to 1001) showed the *old*
  `sin(x)/cos(x)` form was already 0-1 ulp exactly at the actual poles.
  The premise doesn't hold for this crate: `sin`/`cos` preserve *relative*
  precision even as their value shrinks toward zero (their own reduction
  + poly don't suffer cancellation), so dividing two independently-
  accurate values doesn't catastrophically amplify error the way the
  backlog assumed — the real ~3000 max ulp turned out to come entirely
  from `sin`/`cos`'s own well-known large-`|x|` degradation near their
  *domain* cliff (~1.3e7), inherited by any tan built on top, not from
  evaluating near a pole at moderate `x` at all. Measured comparison
  (exhaustive-adjacent fuzz, bucketed by `|x|`): the new reduction is
  *worse* in every practical bucket (`|x|<10`: max ulp 4→5; `|x|<1000`:
  4→5; `|x|<1e6`: 9→17; `|x|<2^22*pi/2≈6.6e6`: 136→9656) and additionally
  has a *narrower* safe domain than the old form (its own cliff sits at
  `2^22*(pi/2)≈6.6e6`, half of `sin`/`cos`'s `2^22*pi≈1.3e7`, since `q`'s
  magnitude scales with `2/pi` instead of `1/pi`). Reverted; `src/lib.rs`
  untouched. **General lesson: before implementing a fix framed around
  "operation X amplifies error near Y," check the actual current
  behavior at Y directly (a handful of spot-check values) — here it took
  under a minute and would have saved the whole implementation effort.**

- **cos (fast tier): dedicated even poly in r²**: q = round(x/pi) same as
  sin, cos(x) = ±cos(r) with cos poly = even, killing the -0.5/+0.5
  offset ops and the x³ multiply. Known risk (why sleef does't do it):
  near r = ±pi/2, cos's value → 0 while an absolute-error even fit stays
  O(2^-24), so *relative* error blows up near cos's zeros — exactly the
  failure the shifted-sin form avoids. Probably dies for that reason;
  listed so the reasoning is recorded rather than re-derived.

- **cbrt: rational correction, (1+r)^(-1/3) ≈ P(r)/Q(r), 2/2 with the same
  4 coefficients as the shipped degree-3 poly (2026-07-08)**: scipy
  `least_squares` found a 2/2 rational fitting the underlying math
  ~100x tighter than the shipped poly (max abs error 1.9e-9 vs 1.9e-7,
  same coefficient count) — looked very promising in isolation. Didn't
  survive contact with the real crate: `tune.rs`'s own coordinate-descent
  (seeded from the scipy fit, not 0.0) converged to max ulp unchanged (2)
  with avg *slightly worse* (0.339→0.359) on the tuning grid, and
  implementing it directly confirmed this on the real exhaustive/fuzz
  pipeline too (avg ulp 0.3113→0.4675, max unchanged at 3) — the
  underlying approximation is tighter, but both forms are already deep
  enough into f32's own rounding-noise floor (a handful of fma/mul/div
  ops each contributing up to 0.5 ulp) that a 100x-tighter *mathematical*
  fit doesn't move the *computed* result's ulp count; the rational's
  extra division adds one more rounding step, roughly a wash on accuracy
  or slightly worse. Also failed on speed, contrary to the backlog's
  "divider is idle" framing: mca latency 35.06→45.06 cyc (+28.5%),
  throughput 1.629→1.657 cyc/elem (+1.7%, also worse) — unlike
  `cbrt_normal`'s existing `rcp` (started immediately, independent of the
  seed chain, so its latency hides behind other work), this new division
  sits *after* the seed/r computation on the critical path, with nothing
  to hide behind. Not adopted; reverted.


## erf / erfc

- **erfc domain split**: if log-space fails, two rationals ([0,2] /
  [2,10]) with one select — both arms computed branchlessly, so ~2x the
  poly cost; only worth it if 109 → single digits (the log-space attempt
  above already found the single-poly approach doesn't converge at a
  practical degree, so this fallback is the natural next thing to try).

## hypot / misc

- **remainder_checked beyond 2^24**: double-float q (qh, ql) like
  sin_checked's reduction, with two correction candidates instead of one.
  Heavy; only worth it if a real use case needs |x/y| > 2^24.


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


- **Centered-variable refit for erfc's xa ∈ [0,10] (hugely off-center
  today)**: the same idea tried for exp2 (above, rejected) — untested here.
  erfc's correction is a P/Q rational, not a plain monomial poly, and its
  target decays across many orders of magnitude over the domain, a very
  different shape from exp2's monotonic-but-bounded R(f) — the exp2
  rejection's reasoning (no interior point beats the domain's own edge
  minimum) doesn't obviously transfer, so this is still worth checking
  independently rather than assuming it fails the same way.

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

- **Unchecked tiers for hypot/atan2 and others**: same argument as
  log_2_unchecked/ln_unchecked/log10_unchecked (implemented, see above) —
  a `hypot`/`atan2` that skips inf special cases would drop real ops the
  same way. The tier split is already this crate's established pattern;
  it just hasn't been applied uniformly everywhere it could be.

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

### exp family

- **Select-tree "LUT" for exp2**: split f = f_hi + f_lo where f_hi takes 4
  (or 8) quantized values; 2^f_hi comes from a 2-level vblendvps tree over
  4-8 constants (no gather needed), and f_lo's range shrinks 4-8x so the
  poly drops from degree 5 to 2-3. Trades ~3 fma levels for ~2-3 blends +
  1 compare chain. Same trick applies to exp's e^r poly. This is the
  gather-free version of the LUT idea round 1 wrote off.

- **Dedicated sinh/cosh kernels on the reduced argument**: instead of
  composing exp twice (or the round-1 even/odd trick), do the reduction
  once (k, r) and fit sinh(r)/cosh(r) minimax polys directly on
  [-ln2/2, ln2/2]: sinh(x) = sinh(r)·cosh_k + cosh(r)·sinh_k where
  cosh_k/sinh_k = (2^k ± 2^-k)/2 come from two exponent bit-tricks. The
  polys are even/odd so they share r² powers. Compare against round 1's
  shared-exp idea; one of the two should win.

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

- **Vectorized Payne-Hanek "exact" tier**: full-range correct reduction
  needs the 2/pi product against x's mantissa with the window selected by
  x's exponent — per-lane variable shifts exist (vpsrlvd, AVX2) and the
  2/pi table is small enough for a 4-8 constant select tree instead of a
  gather. Would make a `sin_exact` with no accuracy cliff anywhere in
  f32. Big job, listed for completeness (the "graceful degradation"
  contract makes it optional, not required).

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

### hypot / atan2 / new surface

- **hypot_checked with branchless exponent rescale**: extract
  max(exp(x), exp(y)) with integer ops, scale both inputs by 2^-e via
  exponent-field subtraction (exact), sqrt, scale back by 2^(e/2)... e
  odd needs a sqrt2 factor — cleaner: scale by 2^-2⌊e/2⌋. All bit tricks
  + selects, fully vectorizable, kills the documented overflow/underflow
  tradeoff for callers who need std-grade hypot without std-grade scalar
  code.

