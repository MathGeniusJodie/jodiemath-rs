# Ideas: faster and/or more accurate (brainstorm only, untested)

Budget: ≤0.5 avg ulp, ≤2 max ulp; accuracy→speed trades allowed up to that.
Hard constraint: everything must autovectorize (select-based, no data-dependent
branches, no scalar-only intrinsics unless the vector form exists).

## Cross-cutting

- **`round_ties_even` instead of `round`**: `f32::round` is ties-away-from-zero,
  which LLVM lowers to a multi-instruction sequence; `round_ties_even` lowers to
  a single `vroundps`. `round_x_over_pi` uses `.round()` twice on its critical
  path — any consistent rounding works there (q just needs to be *an* integer
  consistently), so this is a free latency cut. Audit `remainder` too (there the
  tie behavior is documented, so it'd be a semantic change).
- **Exploit the idle divider**: mca already showed FP divide is nearly free
  while FMA/mul ports are the bottleneck. Several rewrites below deliberately
  move work onto the divider (atanh-form log, tanh/sinh via `1/e`, tan as
  sin/cos of one shared reduction).
- **Quantized-coefficient refits**: lolremez fits real coefficients then rounds
  to f32; Sollya's `fpminimax` (or a final coordinate-descent pass like the
  existing `run_descent`, but driven by the exhaustive accuracy.rs metric)
  optimizes *in the rounded space* and typically buys 0.1–0.3 avg ulp for zero
  runtime cost. Applies to every poly in the crate.
- **Weighted refits**: minimax over the *actual* distribution of reduced
  arguments (including the reduction's own error) instead of the ideal
  interval; the avg-ulp metric rewards this even when max ulp doesn't move.
- **Degree-reduction probes — re-checked via lolremez screening, all three
  still fail by a wide margin, not adopted (2026-07-07).** The stated
  premise (re-litigate under a "looser" 0.5-avg/2-max budget) didn't
  actually change anything: `log_2` deg-9→8 was already conclusively
  tested and rejected earlier (max ulp 3-5 vs. a 2 cap, even after
  tuning — see the "Degree-reduction search" results elsewhere in this
  file, unaffected by a looser *average* budget since the failure is a
  *max* violation). The other two, screened fresh with lolremez (the
  crate's own established pre-check methodology — fit the reduced degree,
  compare estimated max relative error against the current degree's
  fitted error, only bother with a full coordinate-descent+Rust
  implementation if the gap looks borderline): `exp2`'s Q poly, degree
  5→4, `(2^x-1)/x` on `[0,1]` — estimated max relative error 1.01e-8
  (degree 5, matches the ~1e-8 already shipped) → 4.07e-7 (degree 4), a
  **40x** degradation, several ulp worth of poly-introduced error alone
  against exp2's current max-ulp-1 accuracy. `sinf_poly`, degree 9→7,
  `(sin(sqrt(y))-sqrt(y))/y^1.5` on `[0, (pi/2)^2]` — 6.97e-9 (degree 3
  in y, i.e. degree 9 overall, matches sinf_poly's own doc comment's
  ~6.1e-9) → 1.24e-6 (degree 2 in y, degree 7 overall), a **178x**
  degradation — consistent with (probably the same underlying fact as)
  this file's own already-recorded "sin/cos: degree-2 correction poly
  fails hard (82 ulp)" finding from the "Degree-reduction search"
  results. Neither gap is remotely borderline, so neither was carried
  through to a full Rust+coordinate-descent verification — the lolremez
  signal alone is decisive here, and building out a doomed candidate
  just to watch it fail the real sweep would be busywork, not rigor.
  Each dropped term would have been one fma of depth-or-width, but with
  no accuracy budget left to spend it on.
- **Small LUTs via in-register permute**: a 8/16-entry table indexed by top
  mantissa bits autovectorizes as `vpermps` (AVX2) / `vpermi2ps` (AVX-512) if
  written as a const array indexed by `usize` — LLVM does turn small
  constant-table lookups into shuffles, but verify codegen. Table+low-degree
  poly is the classic degree killer for log2/exp2. Risk: LLVM may emit gathers
  instead; needs a godbolt check before investing.
- **AVX-512 instruction targets** (only if `target_feature` allows):
  `vgetexpps`/`vgetmantps` replace the log/cbrt bit-trick exponent splits and
  handle denormals for free; `vscalefps` replaces exp2's exponent-construction
  bit trick *and* the checked variant's two-multiply split (correct overflow
  /underflow built in); `vreduceps` gives exp2's `x - floor(x)` in one op;
  `vrcp14ps`/`vrsqrt14ps` as seeds. LLVM won't synthesize these from portable
  code; would need `cfg(target_feature)` paths with intrinsics that still
  vectorize in loops (they do — they're inherently vector ops).
- **Batch/slice API**: nothing forces callers into a vectorizable loop today;
  a `fn sin_slice(&[f32], &mut [f32])`-style layer (or just doc'd examples with
  `chunks_exact`) makes the autovectorization contract explicit and testable in
  CI via `--emit=asm` grep.
- **Double-float final step as a pattern**: log_2 already gets "one rounding"
  from a final fma. The same trick (keep hi + lo until the last op, fold lo in
  with one fma) applies to ln/log10 (see below), powf, and sin near zeros —
  usually 1–2 extra fmas for a large avg-ulp cut.

## log_2 / ln / log10 / log1p

- **atanh form**: log2(m) = (2/ln2)·atanh(t), t = (m−1)/(m+1). t spans
  ~±0.171 vs s's [−0.293, 0.414], and atanh is odd, so the poly is in t² —
  degree ~4–5 in t² replaces degree 9 in s. Costs one division (idle divider!)
  and the t computation; shorter fma chain, likely both faster *and* more
  accurate (the odd form has better conditioning). Needs the same
  final-fma-with-k trick; t is not exact but computable with a df correction
  if needed.
- **Table+poly**: index top 3–4 mantissa bits, s = m·r_i − 1 with r_i ≈ 1/m_i
  from a table (plus log2(m_i) table folded into k's fma). Range shrinks
  ~8–16×, poly drops to degree 3–4. Vectorization: `vpermps` (see cross-cutting
  caveat).
- **ln/log10: fold the rescale into the poly**: `log_2(x) * LN_2` adds a full
  extra rounding on top of log_2's error (~doubles avg ulp). Instead refit the
  poly with coefficients pre-multiplied by ln2 (i.e. fit s·P(s) → ln directly)
  and handle k·ln2 as Cody–Waite: `fma(k, LN2_HI, poly_ln) + k*LN2_LO` — two
  extra fmas, roughly halves ln's error. Same for log10.
- **log_2 degree-8 probe** with quantized refit (see cross-cutting).
- **log1p, real fix (branchless) — done, tested, kept (2026-07-07).** u = 1
  + x; c = x − (u − 1) (exact by Sterbenz for the interesting range);
  result = ln(u) + c/u. This wasn't a marginal "over budget" case: a fuzz
  sweep across the whole crate (100M samples/fn) turned up log1p at
  **172,186,304 avg ulp / 864,026,618 max ulp** (worst x ≈ 6e-8, exactly
  the documented small-x collapse-to-1.0 cliff) — the naive `ln(1.0 + x)`
  was returning garbage across most of the domain log1p exists to serve,
  not a tuning shortfall. Implemented as specified, plus two edge guards
  the original idea's one-liner glossed over: x = −1 exactly (u = 0) makes
  c/u a literal 0/0, and x = +inf makes it inf/inf (via c = x−(u−1) =
  inf−inf = NaN first) — both should contribute nothing (ln(u) alone is
  already the correct −inf/+inf) but without a guard they poison the sum to
  NaN. Both collapse c/u itself to NaN, so a single `corr.is_finite()`
  check (not two separate checks on u) suppresses both at once. Verified:
  exhaustive (all 2^32 patterns) sweep now shows avg ulp 0.1061 / max ulp 4
  (worst x ≈ 0.13, unrelated to the small-x fix — this is just ln's own
  ~3-ulp fit error compounding by one step, an expected residual, not a new
  bug); edgecheck's `log1p(0)/(-1)/(-2)` and a temporary test covering
  ±inf/NaN/f32::MAX/a small-x sanity check all passed before the temp test
  was removed. mca cost, measured (not just reasoned about, after an
  initial two-guard version cost noticeably more before collapsing to the
  one-check form above): 59.98→61.14 cyc latency (+1.9%), 2.084→2.339
  cyc/elem throughput (+12.2%) — a real but modest cost, mostly one
  division (the crate's established idle-divider finding) plus the
  is_finite guard. Kept despite the nonzero perf cost: this crate's stated
  budget (≤0.5 avg / ≤2 max ulp) treats accuracy as the primary constraint
  and speed as the thing traded *within* it, and a 172-million-ulp average
  is a correctness bug, not a trade to weigh against a ~12% throughput
  cost on one function.
- **Denormal pre-scale via multiply**: current tiny-path is fine; with AVX-512
  `vgetmantps`/`vgetexpps` the whole tiny/koff dance disappears.

## exp2 / exp2_checked / exp / expm1

- **exp: skip the log2e pre-multiply's rounding**. `exp2(x * LOG2_E)` puts a
  rounding error on the *argument*, amplified by the derivative — the dominant
  error for |x| ≳ 1. Standard fix: k = round(x·log2e), r = x − k·LN2_HI −
  k·LN2_LO (Cody–Waite, exact), poly for e^r on [−ln2/2, ln2/2], scale by 2^k
  with the existing bit trick. ~2 extra fmas, big avg-ulp win for exp and
  everything built on it (sinh/cosh/tanh/erf/erfc/powf).
- **Table+poly**: f = j/16 + f', 2^(j/16) from a 16-entry `vpermps` table,
  poly degree drops from 5 to ~2–3. Cuts both latency and width.
- **Q(f) degree-4 probe** with quantized refit against the ≤0.5-avg budget.
- **exp2_checked with `vscalefps`** (AVX-512): one instruction replaces the
  k1/k2 split, both bit-trick constructions, and both multiplies, with correct
  inf/denormal semantics. The portable version stays as fallback.
- **expm1, proper reduction**: k = round(x/ln2), r reduced; expm1 =
  2^k·(P(r) + 1) − 1 where for k=0 the +1/−1 cancel exactly if P returns
  (e^r − 1) directly. Do the 2^k·P + (2^k − 1) combine with fma; fixes both
  the Padé/exp seam and the inherited exp domain limit in one shape. All
  select-based.
- **expm1's Padé division** is fine (idle divider), but the select boundary
  at 0.5 could be re-tuned after the exp accuracy fix.
- **expm1's Padé branch refit — done, tested, kept (2026-07-07), ninth use
  of this session's tuning recipe, fifth real win, and a methodology
  lesson about the tuner's own grid resolution.** The 5 coefficients
  (-2, -120, -12, 60, -120 -- an exact closed-form Padé identity to e^x,
  not an empirical fit, given the small integers and the shared -120)
  extended into `expm1_near0_c`, tuned as 5 independent free parameters
  against `f64::exp_m1` over `|x| < 0.5`. tune.rs's coordinate-descent
  grid (~1.7M points, its normal density for a 5-coefficient search)
  reported max ulp 3→2; a targeted follow-up exhaustive check (a
  temporary scalar test, 302M samples via step-7 over the full branch
  domain, not just tune.rs's sparser grid) found the *real* numbers are
  max ulp 4→3, avg ulp 0.111→0.109 — the coarse grid had missed the true
  worst input on *both* sides, undercounting by exactly 1 ulp each time
  (a wash on the reported delta's shape, but a reminder that a tuner's
  own scoring grid isn't the same guarantee as an exhaustive or
  near-exhaustive sweep, and applying a tuned result deserves its own
  verification pass even when the grid-reported numbers look clean).
  Still a genuine, if smaller-than-first-reported, improvement. Confirmed
  zero perf cost via mca (71.00 cyc / 1.441 cyc/elem, bit-for-bit
  unchanged) — same instructions, only the 5 literal constants differ.
- **exp2's and log_2's own polys — checked, zero headroom, confirming
  they're already at their coordinate-descent local optimum
  (2026-07-07).** tune.rs already had `exp2_c`/`log2_c` tuners from
  before this session (evidence in themselves that these two were tuned
  this way previously) — their `main()` init arrays were stale
  pre-tuning starting points, not the currently-shipped coefficients,
  so updated both to match `src/lib.rs` exactly before re-running. Result
  for both: the tuner's "tuned" output was *bit-identical* to "start" —
  not just a negligible move like acos_poly/erf_tail/sinf_poly, a
  literal zero-move local optimum, the strongest confirmation yet that
  "already been through dedicated tuning before" predicts no further
  headroom. (log_2's leading coefficient, exactly `LOG2_E` by
  mathematical necessity — the poly's Taylor-derivative leading term, not
  an empirical fit — was excluded from tuning via a new `tune_fixed0`
  helper, to avoid the search ever suggesting breaking that identity.)
  This session's refit scorecard, final tally: 5 real wins (asin, atan,
  erfc, cbrt, expm1) vs. 6 no-ops (acos_poly, erf_tail, erf_near0,
  sinf_poly, exp2, log_2) — every no-op case had a documented prior
  tuning history; every win was on coefficients that (as far as this
  session found) hadn't been touched since the original C port.

## sin / cos / sin_checked / cos_checked / tan

- **`round_ties_even` in `round_x_over_pi`** (see cross-cutting) — two
  `.round()` calls on the ql critical path become single instructions.
- **Fused `sincos` / direct `tan` — tried, measured, reverted (2026-07-07),
  a real numerical wall, not just a missed optimization.** The identity
  is sound: for x = r + q·π (r in [-π/2, π/2], sin's own reduction),
  sin(x) = (-1)^q·sin(r) and cos(x) = (-1)^q·cos(r), so the `(-1)^q`
  factors cancel *exactly* in tan(x) = sin(x)/cos(x) = sin(r)/cos(r) --
  no parity computation needed at all, just one shared reduction. Fit a
  dedicated `cosf_poly` via lolremez (`--degree 4 --range 0:(pi/2)^2
  "(cos(sqrt(x))-1)/x"`, mirroring sinf_poly's own "poly in x²" shape) --
  the fit itself was excellent in isolation, ~3.6e-10 relative error,
  *tighter* than sinf_poly's own 6.1e-9 (cos's even series converges
  faster at a comparable degree). Implemented `tan(x) = sinf_poly(r) /
  cosf_poly(r)` using sin's shared q/r. Exhaustive sweep: catastrophic
  regression, avg ulp 0.33→1.05, max ulp ~3000→32,434,460, worst x ≈
  252.9 (≈80.5π, right next to a tan pole). Root cause, confirmed by
  direct probe: `cosf_poly` is a `1 + y·R(y)` additive form, and near
  r = π/2 (cos's own zero -- and exactly where tan's poles put the most
  weight), `y·R(y)` must land within a hair of exactly -1 for the sum to
  be small and accurate; a 42%-relative-error result at
  `cosf_poly(π/2_f32)` despite a ~2e-8 *absolute* error shows the
  additive form cancelling exactly like every other "near-zero" bug
  fixed elsewhere this session (log1p, sinh, asinh, acosh, asin) -- just
  here the "small x" is `r` near the *domain edge* of a fresh poly, not
  near 0. This is precisely the trap the crate's *existing* `cos()`
  already engineers around: it phase-shifts its own reduction (a
  differently-decomposed q, k = round(x/π - 0.5)) specifically so it can
  reuse `sinf_poly` — already accurate near *its own* comfortable zero at
  r=0 — instead of ever needing a fresh poly evaluated near a
  cancellation point. Sharing sin's literal q/r for tan throws that
  design away and re-introduces the exact problem it was built to avoid.
  A real fix would need the same kind of Sterbenz/rationalization
  correction near cos's zero that fixed asinh/acosh/asin (e.g. compute
  `cosf_poly` via a term that stays well-conditioned as `r → π/2`, not a
  bare `1 + y·R(y)`), which is a bigger, riskier redesign than this pass
  attempted -- reverted rather than half-fixed. `cosf_poly` and the new
  `tan()` were both fully removed, not left disabled.
  **Second attempt, same session: tried the "obvious" fix, made it
  worse.** Instead of a dedicated `cosf_poly`, tried the cofunction
  identity directly — `cos(r) = sin(pi/2 - |r|)`, reusing `sinf_poly` for
  both numerator and denominator, no second poly at all. This looked like
  exactly the right fix (it's *literally* what `cos()` computes, just
  derived post-hoc from sin's `r` instead of from a dedicated reduction).
  Measured instead of trusting the symmetry: avg ulp 0.33 → **28.2**, max
  ulp ~3000 → **3,437,483,831** — worse than the first attempt, not
  better, same worst-x (≈252.9). Root cause: `FRAC_PI_2 - r.abs()` is a
  single-precision f32 subtraction, and near a pole `r` is close to
  `FRAC_PI_2` itself — the exact same catastrophic-cancellation shape
  being chased, just relocated from evaluating `cosf_poly` near its zero
  to computing the *input* to `sinf_poly` near zero. `cos()`'s own
  reduction avoids this because it computes its phase-shifted residual
  via the *same* Cody-Waite-precision PI_A..D chain used for everything
  else in this reduction family (effectively a `pi/2` that carries much
  more than f32's ~24 bits), not a single further subtraction from a
  bare f32 constant afterward. A real fix would need to synthesize that
  same multi-word precision for the cofunction transform specifically
  (splitting FRAC_PI_2 into HI/LO words à la PI_A/PI_B, Sterbenz-correct
  subtract, etc.) — comparable in complexity to `reduce_pi` itself, not a
  small delta on top of the first attempt. Reverted again. Two failed
  attempts at the same underlying idea in one session is a strong signal
  this needs the full double-float reduction machinery or nothing —
  worth remembering before a third attempt without that machinery.
- **Reduction mod π/2 instead of π**: residual lands in [−π/4, π/4]; sinf_poly
  drops ~1 term and a cos poly (even, degree 8 → 4 terms) appears. For plain
  sin you then need a select between sin-poly and cos-poly by octant — both
  get computed and blended (vectorizes, but doubles poly width), so it's only
  a win where the polys are shared (sincos/tan) or where latency dominates
  throughput. Worth an mca measurement, not obvious.
- **Residual as double-float into the poly**: reduce_pi already computes
  s3 + err; instead of collapsing, keep (rh, rl) and do
  `s = fma(rl, cos_approx(rh), sinf_poly(rh))` where cos_approx ≈ 1 − rh²/2
  (rh² is already computed inside sinf_poly — reuse y). ~2 extra fmas; cuts
  the max ulp near sin's zeros at multiples of π, which is where the
  worst cases live. Same for the fast sin via its r (the PI_A..D chain
  already produces the residual in pieces — the last fma's rounding is
  the loss).
- **Parity via bit ops in the checked versions**: `s * (1.0 - 2.0 * par)`
  is two mul/fma-port ops; parity(qh)^parity(ql) can become a sign-bit XOR
  (`f32::from_bits(s.to_bits() ^ ((pq_bits ^ pl_bits) << 31))`) like the fast
  path already does with qb's mantissa bit. parity() itself (floor-based) might
  reduce to grabbing the low mantissa bit of `q + ROUND_MAGIC` when |q| is in
  magic range — qh isn't (it can be huge), but its low bit is also recoverable
  from `qh * 0.5`'s fract; needs care, but the multiplies saved are on the
  final dependency chain.
- **Fast sin/cos: one more Cody–Waite word vs. documented cliff**: PI_A..D is
  already 4 words; the cliff is in q's exactness, not π's. A cheap middle tier
  — q from `fma(x, FRAC_1_PI, ROUND_MAGIC)` but with one `two_prod`-corrected
  x·(1/π) term — might push the fast path's valid range from ~1.3e7 to ~1e9
  for 2–3 fmas, without the full checked machinery. Worth mapping the
  speed/range Pareto point.
- **sinf_poly quantized refit — tried, no meaningful headroom found, not
  applied (2026-07-07).** Eighth use of this session's tune.rs recipe,
  extended with `sinf_poly_c` against `f64::sin` over `[-pi/2, pi/2]`
  (the poly's fitted domain). Max ulp unchanged (2→2), avg ulp barely
  moved (0.00248→0.00244) — both numbers already near f32's own precision
  floor, unlike the other "no headroom" cases (acos_poly/erf_tail/
  erf_near0) which had more room to begin with but still didn't move.
  Not surprising in hindsight: sinf_poly and the reduction it feeds have
  already been through multiple dedicated tuning/rewrite passes earlier
  in this crate's history (see jodiemath-workflow memory), unlike
  asin/atan/cbrt/erfc which hadn't. Confirmed the coarse-grid result was
  trustworthy without waiting on a denser-grid re-check this time (a 10x
  grid run stalled and was killed after ~250s, matching a pattern from
  three earlier "no headroom" cases in this file, all of which had their
  coarse-grid conclusion independently confirmed by a completed dense
  run) — treated as sufficient given the string of prior confirmations.
  Sin/cos's poly evaluation itself is also already Estrin-scheduled (see
  the "Polynomial evaluation schemes" section), another sign this
  specific poly isn't the crate's low-hanging fruit anymore. Also try
  fitting sin/π-scaled variants so the reduction constant folds in.

## cbrt family

- **Quantized/descent refit of the 4 correction coefficients** in
  cbrt_normal -- already done, and the legacy tool this entry pointed at
  removed as dead code (2026-07-07). The random-sampling `run_descent`
  this entry named was a leftover early exploratory tool: unused (its
  only caller, `#[test] fn descent2`, was already commented out) and
  fully superseded by `examples/tune.rs`'s later, proper coordinate
  descent (`cbrt_normal_c`, driven by an exhaustive grid, not random
  samples -- exactly what this entry asked for). The current shipped
  `cbrt_normal` coefficients (`-0.33333147, 0.22220612, -0.17394388,
  0.14823665`) already came from that tuner, from an earlier session
  (see readme.md/jodiemath-workflow memory). Deleted `run_descent` and
  the commented-out `descent2` test (with it, the now-unused `use
  rand::Rng`/`use rand::RngExt` imports in that module) -- confirmed
  dead via the compiler's own "never used" warning, not just unreferenced
  in this file. `rand` itself stays a dev-dependency (still actively used
  by `examples/accuracy.rs`).
- **Seed constant joint search — tried, measured, rejected (2026-07-07):
  degree-2 (dropping c4) can't come close to budget, at any seed.** Wrote
  a standalone (isolated, not touching the shipped code) joint search:
  for each of 41 seed offsets around the shipped `0x2a509a07` (every
  other value in `-40..=40`), coordinate-descended a degree-2 poly (3
  coefficients, one fewer than shipped) against `f64::cbrt` over a grid
  spanning `x in [1,2)` (representative of every octave). First attempt
  used a naive coefficient starting point (`[-0.3333, 0.2, -0.1]`) and
  got catastrophic results (max ulp 2607, avg 921) — looked like the
  search was just failing to converge, so retried from a much better
  starting point (the shipped degree-3's own first 3 coefficients,
  simply dropping c4) with 4x more coordinate-descent rounds allowed.
  Still nowhere close: best found across all 41 seeds was max ulp 112,
  avg 29 — over 50x worse than the shipped degree-3's budget (avg 0.33,
  max 2), confirming this isn't a search-quality artifact but a genuine
  expressiveness gap: 3 coefficients (degree-2) can't correct this seed's
  error to budget regardless of which seed is chosen. Also moot as a perf
  idea even if it *had* worked: degree-2 via Horner (`fma(fma(c2,r,c1),
  r,c0)`) is 2 fma's at depth 2, vs. the shipped degree-3's 3 fma's *also*
  at depth 2 (`ceil(log2(3))==ceil(log2(4))==2`, the same "doesn't cross
  a power-of-2 boundary" non-savings already documented elsewhere in this
  file) — dropping c4 would only ever have saved one throughput-only op,
  never latency. Not adopted; no changes to shipped code (tested in an
  isolated standalone copy, same pattern as the sin_checked/cos_checked
  clamp-idea rejection earlier this session).
- **AVX-512 `vgetexpps/vgetmantps/vscalefps`** kills cbrt's tiny/scale select
  dance and cbrt_accurate's three-way rescale entirely: reduce mantissa to
  [1,2), cbrt it, scale by e/3 via vscalef with the e mod 3 residue folded
  into a 3-entry table.
- **cbrt_accurate: Halley instead of Newton** from a cheaper seed — one Halley
  step from the bit seed + deg-1 correction might match Newton-from-deg-3 with
  less total width. The df32 machinery already exists.
- **Integer-division-free seed — codegen confirmed as predicted, but the
  cheaper seed doesn't hold budget under cbrt_normal's own poly degree
  (2026-07-07).** Checked `cbrt_normal`'s `ax / 3` fresh in `--emit=asm`:
  it does vectorize cleanly, exactly as this entry hoped -- the
  throughput region shows the classic constant-division-by-multiply-high
  sequence (2x `vpmuludq` + `vpshufd` + `vpermt2d` + `vpsrld`, 5
  instructions, no real `div`), not a per-lane scalar fallback.
  `cbrt_fast`'s `(bits>>16)*0x5556` alternative is genuinely cheaper in
  raw instruction count, confirmed the same way (1 `vpsrld` + 1 native
  `vpmulld` + 1 add, ~3 instructions) -- native 32-bit multiply-low needs
  none of division's widening/shuffle/merge dance. So the codegen
  question resolves in favor of trying the cheaper seed. But swapping
  `cbrt_normal`'s seed to that form (keeping its sign-safe wrapper and
  degree-3 correction poly, only the seed construction changed) and
  refitting the poly's 4 coefficients via coordinate descent for the new
  seed's error distribution doesn't recover cbrt's budget: best found
  (even after retuning, coarse grid) is max ulp 33 / avg 5.07, ~16x over
  budget (avg<=1, max<=2) -- the seed itself is just too coarse for a
  degree-3 correction to fully compensate, at least reusing `cbrt_fast`'s
  own magic constants (`0x2a4ddef1`/`0x5556`, presumably tuned for
  `cbrt_fast`'s own different downstream double-Newton-refinement
  structure, not for this poly-correction shape). Not adopted. A further
  attempt would need to jointly re-derive the seed's own magic constants
  for this specific combination (same spirit as the rejected "seed
  constant joint search" entry above, but with degree-3 kept instead of
  cut to degree-2) -- bigger scope than this pass, not attempted.
  `examples/tune.rs` gained `cbrt_shiftmul_c` / `which.contains(
  "cbrtshift")` as reference infrastructure for that future attempt, not
  deleted since it's small and inert unless invoked. No `src/lib.rs`
  changes.

## asin / acos / atan / atan2

- **asin small-x fix — done, tested, kept (2026-07-07), seventh and last of
  this session's correctness fixes.** Two bugs, not the one documented,
  found and fixed one at a time (same pattern as acosh/asinh): (169,116,394
  avg ulp / 852,038,214 max ulp overall, worst x ≈ 2.3e-8).
  1. `sqrt(1-a) - 1` (the final step of the acos-style sqrt identity) has
     the same cancellation as asinh/acosh's `sqrt(1+t) - 1`. Fixed with the
     same rationalization, `sqrt(1-a) - 1 = -a / (sqrt(1-a) + 1)`, leaving
     the fitted rational correction `a = (a²-a)/d + a` untouched. This
     alone dropped avg/max ulp to 293.15/825 — a huge improvement (~577,000x
     on avg) but still nowhere near the 0.5 budget.
  2. What's left after fix 1 isn't cancellation at all: the rational
     correction's own fit carries a small persistent *relative* bias
     (~3e-5, traceable to the C original's coarser ~1e-6-relative goal)
     that the much-larger cancellation error had been hiding. In relative
     terms this bias is roughly constant across the domain, but it
     dominates specifically where the true answer is smallest (x → 0),
     exactly mirroring why log1p/tanh/sinh needed small-x branches. Fixed
     the same way: a genuine 4-term Taylor branch (`x + x³/6 + 3x⁵/40 +
     15x⁷/336`, exact rational coefficients) below `|x| < 0.1` — a much
     smaller cutoff than sinh's 0.5, since asin's Taylor series converges
     far slower near its sqrt singularity at x=1 (checked: 2 terms alone
     leave ~63 ulp at x=0.1, not enough; 4 terms leave ~2.7e-10 relative,
     comfortable margin at that cutoff). The rational correction was left
     alone rather than refit, since it's only inaccurate in the regime the
     Taylor branch now covers.
  Sign-of-zero bonus: the old formula's zero sign was backwards from the
  rest of this crate's odd functions (`asin(+0)` came out `-0.0`) — the
  rationalized form happens to produce the mathematically conventional
  sign instead, confirmed by checking atan/sinh/asinh/tanh(-0) all agree
  with the new convention; edgecheck updated to match (this was a
  pre-existing, separately-documented quirk, not something this fix set
  out to change, but worth noting since the expected values flipped).
  Verified exhaustive: avg ulp now 0.328 (in budget), max ulp 468 (real
  progress from 852 million, but not under the stated ≤2) — concentrated
  at x → 1, a genuinely different mechanism (asin's derivative singularity
  amplifying the same rational-correction's own fit inaccuracy, this time
  at the *other* domain edge) that a small-x branch can't reach. Left open
  as a documented follow-up (see asin's doc comment) rather than bundling
  in a second, near-1 branch this round. mca cost, in line with the
  session's other fixes: 56.11→89.02 cyc latency (+58.6%), 1.433→2.124
  cyc/elem throughput (+48.2%). Kept for the same reasoning as the other
  six: the avg-ulp fix alone is a correctness fix by any reasonable bar,
  even though the max-ulp story isn't fully resolved.
- **asin's x→1 residual, closed the same session — done, tested, kept
  (2026-07-07), a direct continuation of the fix above.** Rather than a
  refit, reused acos's *already-fitted* `acos_poly` via the identity
  `asin(x) = π/2 − acos(x)`: for a ≥ 0, `acos(a) = sqrt(1-a)·acos_poly(a)`
  is exactly acos's own well-conditioned formula (a shrinking sqrt factor
  times a smooth bounded poly, no cancelling subtraction), and `π/2 -
  acos(a)` doesn't cancel either since acos(a) is small near a=1 while
  π/2 is O(1) — no new poly fit needed, gated to `a > 0.9`. Verified
  exhaustive: max ulp 468→121 (avg ulp barely moved, 0.328→0.325, already
  in budget). The worst case *relocated* rather than vanishing — to the
  small/mid branch boundary (x ≈ 0.1) — confirming the theory that both
  residuals share the same root cause (the mid branch's ~3e-5 relative
  bias, which scales with x and therefore can't be shrunk by moving
  thresholds, only by refitting the mid branch or extending the Taylor
  branch with more terms). mca: a genuine mixed result, latency actually
  *improved* on top of the accuracy gain (89.02→43.24 cyc, -51%) while
  throughput got worse (2.124→2.899 cyc/elem, +36.5%) — three branches
  computed unconditionally costs real throughput, but apparently let the
  scheduler shorten the critical path further, the same non-monotonic-
  scheduling surprise logged elsewhere in this file (asinh's fix showed
  the identical shape). Net win on 3 of 4 measured axes (avg ulp, max ulp,
  latency), real cost on the fourth (throughput) — kept.
- **asin mid-branch refit — done, tested, kept (2026-07-07), a direct
  continuation of the two entries above using real numerical tooling for
  the first time this session.** Extended `examples/tune.rs` (previously
  only wired up for exp2/log2) with `asin_mid_c`, replicating the mid
  branch's shipped formula (rational correction + the already-fixed
  rationalized sqrt step) so the coordinate-descent tuner scores it
  against `f64::asin` over a grid matching exactly where the mid branch is
  used in production (`a` in `[0.1, 0.9)`). First attempt used a much
  denser grid (10.7M points) and simply didn't finish in reasonable time;
  settled on ~3.6M points as the practical ceiling for this coordinate-
  descent implementation (evaluates the *whole* grid per coefficient-nudge
  trial). Tuner's lexicographic (max, then avg) objective found max
  120→83 (on-grid) with avg getting *worse* (25.1→29.3) — tried to find a
  Pareto-better point with a custom avg-first/max-capped variant, which
  failed cleanly and informatively: for every cap tried between the start
  and the tuned max, zero valid first moves existed, meaning the real path
  to max=83 crosses through intermediate states worse than those caps, not
  reachable by a search that must stay under a fixed ceiling at every
  step. Removed that variant (dead end for this case, not generically
  useless, just didn't help here) rather than leave unused code around.
  Applied the plain tuned result to `src/lib.rs`: exhaustive sweep
  confirmed the grid-based prediction translated directly — max ulp
  121→84 (~30%), avg ulp 0.325→0.377 (still comfortably under the 0.5
  budget). Zero perf cost (mca bit-for-bit unchanged, 43.24 cyc / 2.899
  cyc/elem — same instructions, only the 4 literal constants differ).
  Not a complete fix (the doc comment is explicit that a genuinely
  different correction shape, not just retuned coefficients, would be
  needed to go further), but real, measured, free progress on the
  specifically-flagged open residual from the entry above.
- **asin small/mid threshold widened 0.1 -> 0.3 — done, tested, kept
  (2026-07-07, later same day), and a case of a documented "can't be
  fixed this way" claim turning out to be untested and wrong.** The mid
  branch refit's own doc comment (fix 4, entry above) asserted the
  x~0.1 boundary residual "can't be shrunk further by moving thresholds
  around" since the bias "scales with x" — a plausible-sounding claim
  that was never actually checked against `asin_small`'s own error curve.
  Probed both branches independently instead of trusting it: the mid
  branch's error isn't a narrow spike right at x=0.1, it's elevated
  (~55-82 ulp) across the *whole* `[0.1, 0.3)` band, while `asin_small`
  (an exact-coefficient degree-7 Taylor series, no fitting involved)
  stays at 1-2 ulp all the way out to x<0.2 and is still comparable to
  mid's own error around x=0.3 (that's where the two curves cross).
  Coordinate-searched the threshold directly (a standalone harness
  scoring `asin_small` vs. the mid formula against `f64::asin` over
  `[0.0005, 0.9)`) rather than guessing: max ulp bottoms out around
  `0.30-0.32` (39-42), degrading sharply on either side (81 at 0.2, 68 at
  0.34, since asin_small's own truncation error grows fast for `x >
  0.35`). Picked the clean `0.3` from that flat minimum. Since *both*
  branches are already computed unconditionally in this branchless
  select regardless of which one gets chosen, moving the threshold
  constant is a pure accuracy change with zero perf implication no
  matter where it lands — confirmed via mca (43.24 cyc / 2.899 cyc/elem,
  bit-for-bit unchanged) and via `--emit=asm` reasoning (same
  instructions, only the comparison immediate differs). Exhaustive sweep
  (all 2^32 f32 bit patterns): max ulp 84 -> 41, avg ulp 0.377 -> 0.105 —
  both improved together, not a tradeoff, unlike the mid-branch refit
  above which traded avg for max. `asinh` (same accuracy.rs filter
  substring match) confirmed unaffected: 0.173/4, bit-for-bit identical
  to its own pre-change baseline. Lesson worth generalizing: a doc
  comment's own causal explanation for *why* something can't be improved
  is itself a claim, not a fact, unless it was actually measured — this
  crate has plenty of prior entries where trusting analysis over
  measurement went wrong in the "adopt something that doesn't help"
  direction, this is the first one this session caught in the opposite
  direction (a *rejection* that turned out to be wrong, leaving a real,
  free win on the table until re-checked).
- **asin's `mid` branch (a separate rational-correction formula) removed
  entirely, replaced by reusing `acos_poly` down to a < 0.25 — done,
  tested, kept (2026-07-07, immediate follow-up to the threshold entry
  above).** Investigating *why* `near1` (the `pi/2 - acos_poly(a)`
  branch, previously only used for `a > 0.9`) stayed accurate well below
  its own cutoff raised the obvious question the original 3-branch design
  never asked: `acos_poly` is `acos`'s *own* poly, fit and used across
  `acos`'s entire `[0,1]` domain in `pub fn acos` itself, not something
  tuned specifically "for near 1" — the name was just an accident of
  where it was first reused. Measured `near1_branch` in isolation across
  the full domain: garbage near x=0 (the identical `pi/2 - poly(0)`
  cancellation the `small`/`mid` split exists to avoid), but from x~0.25
  onward it's *already better than `mid` ever was anywhere in `mid`'s own
  former `[0.1, 0.9)` domain* — `mid` wasn't covering a gap `near1`
  couldn't, it was simply never tried there. Coordinate-searched a direct
  2-branch (`asin_small`/`near1`) crossover instead of guessing: flat
  minimum around `a < 0.25`. Exhaustive sweep: max ulp 41 -> 11, avg ulp
  0.105 -> 0.033 (both improved again, third accuracy win in a row on
  this function today). mca is a genuine mixed result — throughput
  improved a lot (2.899 -> 0.968 cyc/elem, -67%, `mid`'s whole 3-fma +
  2-division + sqrt chain removed from every vectorized call) but latency
  got *worse* (43.24 -> 59.03 cyc, +37%): the same non-monotonic-
  scheduling surprise logged elsewhere in this file (asinh's own fix
  showed the identical shape) — with 3 branches computed unconditionally,
  the scalar latency chain had independent work to fill cycles otherwise
  spent waiting on the sqrt; with only 2, less slack exists. Kept anyway:
  accuracy and throughput (the metric this crate's vectorization-first
  design prioritizes) both improved by a wide margin, only the secondary
  diagnostic latency number regressed, the same shape of tradeoff already
  accepted for asinh/acosh's log1p-fix side effect earlier this session.
  This also mostly fulfills the "acos/asin shared kernel" idea below —
  asin now directly reuses acos's own poly for its whole non-small
  domain, not a separate fit.
- **acos/asin shared kernel — mostly done as a side effect of the `mid`
  removal above.** Originally: both reduce to sqrt(1−a)·poly; a shared
  computation with different post-transforms would halve the code and
  enable a combined refit at f32-quantized precision. `asin` now calls
  `acos_poly` directly (`pi/2 - sqrt(1-a)*acos_poly(a)`, sign-restored),
  the exact shared-kernel shape described here — and the combined refit
  was attempted immediately after, see the entry directly below (done,
  small but real upside for asin, zero cost to acos).
- **`acos_poly` refit against a joint acos+asin objective — done, tested,
  kept (2026-07-07), immediate follow-up to the `mid`-removal entry
  above.** `acos_poly` had only ever been tuned against `acos`'s own ulp
  error (`examples/tune.rs`'s existing `acos_poly_c` tuner). But `acos(x)`
  shrinks to 0 as `x -> 1` while `asin(x)` (now built directly on the
  same poly via `pi/2 - sqrt(1-a)*acos_poly(a)`) grows to `pi/2` there —
  the *same absolute* poly error carries a *different relative* (ulp)
  weight depending on which caller's output magnitude it's measured
  against. A poly tuned purely for acos's own hardest region (`x -> 1`,
  where acos's own value is tiny, so absolute error there is heavily
  penalized in ulp terms) might be spending precision acos doesn't
  strictly need at the cost of precision asin does need elsewhere.
  Extended `acos_poly_c` in tune.rs with two variants: (1) an
  unconstrained joint objective (per-point score = `max(acos's ulp
  error, asin's ulp error)`) — found real improvement on the joint
  metric (grid max 9→6), but the resulting coefficients let acos's own
  exhaustive max ulp regress from 4 to 5, a genuine cross-function
  tradeoff. (2) A constrained variant instead: minimize asin's error
  *subject to* acos's own on-grid max ulp never exceeding its
  already-tuned best — found asin gains (grid max 9→7) with acos's own
  metric provably untouched by construction, not just empirically close.
  Applied the constrained result. Exhaustive sweep confirms both
  predictions: asin max ulp 11→9, avg 0.033→0.030; acos itself exactly
  unchanged (max ulp 4, avg 0.496, bit-for-bit identical to its pre-refit
  values) — not merely "close", genuinely untouched, since the search
  never accepted a move that would have regressed it. Zero perf cost for
  either function (same instructions, only the 7 literal constants
  differ). The rejected unconstrained variant is worth remembering as
  its own small lesson: an unconstrained joint objective can quietly
  trade one caller's accuracy for another's even when both callers share
  the exact same underlying computation — worth checking explicitly
  (constrain or report per-function breakdowns) whenever a poly serves
  more than one caller with differently-shaped error sensitivity, rather
  than trusting a single combined metric to mean "both improved."
- **acos_poly: Horner -> Estrin restructuring — tried, measured, real
  latency win, rejected anyway for costing protected accuracy
  (2026-07-07), immediately after the joint refit above.** `acos_poly`
  is a plain 6-deep sequential Horner chain (7 coefficients, no
  regrouping) — the odd one out in this crate, which elsewhere
  consistently uses Estrin for exactly this kind of poly (e.g. `log_2`'s
  own 10-coefficient poly, `ceil(log2(N))`-deep instead of `N-1`-deep).
  Regrouped into the standard Estrin split (3 pairs at depth 1, combined
  via `x²`/`x⁴` at depths 2-3 — 3 fma's deep instead of 6, same 6 fma
  total plus 2 extra plain multiplies for `x²`/`x⁴`), keeping the *exact
  same* coefficients (a pure regrouping, the same trick that was free for
  `exp2`'s Q poly and `log_2` itself). mca: a genuine, substantial
  latency win on both callers — `asin` 59.03→49.97 cyc (-15.4%,
  recovering over half of fix 6's earlier latency regression), `acos`
  37.11→29.11 cyc (-21.6%, an unrelated bonus, since acos_poly is also
  `acos`'s own core computation) — at a small throughput cost (`asin`
  +5.5%, `acos` +1.7%, the 2 extra multiplies). But unlike `exp2`/`log_2`'s
  own successful regroupings, this one is *not* accuracy-neutral: `fma`
  reassociation changes which intermediate values get rounded when, and
  the shipped coefficients were coordinate-descended specifically against
  the Horner evaluation order. Exhaustive sweep with unchanged
  coefficients: asin max ulp 9→12, acos max ulp 4→5 — both regressed.
  Tried recovering it by re-tuning the coefficients for the *Estrin* form
  specifically (reusing the constrained joint acos/asin search from the
  refit above) — made it *worse*, not better: asin max ulp 9→11 (closer,
  still regressed), but acos max ulp 9→**6** (worse than the unretuned
  Estrin's 5). Root cause of the retune failing to help: the constrained
  search's "acos must not regress" check only runs against the tuning
  grid (~12,900 points), not the full 2^32 exhaustive sweep — the earlier
  Horner-form joint refit's grid-level guarantee happened to hold
  exhaustively too, but this Estrin-form attempt's didn't, exposing that
  the earlier "genuinely untouched, not just close" claim was only ever a
  grid-level guarantee, not a proven exhaustive one. Reverted both
  `src/lib.rs`'s `acos_poly` and `examples/tune.rs`'s `acos_poly_c` back
  to Horner with the pre-Estrin coefficients (confirmed via mca: numbers
  back to bit-for-bit the same as before this attempt). Not adopted:
  `acos`'s exact accuracy was specifically protected by name in the
  immediately-preceding commit (a constrained search built *for the
  purpose* of guaranteeing it wouldn't regress), and this would undo that
  guarantee for a latency-only win on a metric this crate's own readme
  treats as secondary to throughput. A future attempt would need either a
  denser tuning grid (closer to exhaustive) or accepting the Horner
  form's accuracy as a hard constraint while searching for a different
  Estrin-compatible coefficient set — not attempted further this round.
- **erf_poly: Horner -> Estrin restructuring — done, tested, kept
  (2026-07-07), immediate retry of the acos_poly idea just above on a
  different poly with the same 6-deep-Horner shape.** `erf_poly` (erf's
  tail branch, degree-6, 7 coefficients) is the same "odd one out" shape
  acos_poly was, but with one important difference: it's used by exactly
  one caller (`erf` itself), so there's no cross-function accuracy
  guarantee to protect the way asin's use of acos_poly was — lower risk
  to try. Applied the identical Estrin regrouping (3 fma's deep instead
  of 6, same coefficients, 2 extra plain multiplies for `x²`/`x⁴`).
  Unlike the acos_poly attempt, this one turned out to be genuinely
  accuracy-neutral: exhaustive sweep avg/max ulp exactly unchanged
  (0.319/5, matching the pre-change documented baseline precisely) —
  fma reassociation doesn't uniformly cost accuracy, it has to be
  checked per poly rather than assumed either way from one data point.
  mca showed an even better result than acos_poly's own attempt: both
  latency *and* throughput improved together (102.74->91.74 cyc, -10.7%;
  3.163->2.871 cyc/elem, -9.2%) rather than trading a little throughput
  for latency. `erfc` (same accuracy.rs filter substring, doesn't even
  use `erf_poly`) confirmed unaffected as a sanity check: 0.311/109,
  bit-for-bit matching its own pre-change baseline. Kept.
- **erfc's n/d rational chains: Horner -> Estrin — tried, marginal/mixed,
  not adopted (2026-07-07), immediate third attempt at the same
  restructuring, this time on the two smaller degree-4 (5-coefficient)
  Horner chains inside `erfc`'s Padé-style tail.** Estrin for 5
  coefficients only saves 1 depth level (`ceil(log2(5))=3` vs Horner's
  4), a smaller theoretical win than acos_poly/erf_poly's 6-deep-to-3
  halving, and `n`/`d` were already independent (computable in parallel)
  even in Horner form, so there was less obvious slack to reclaim.
  Regrouped both (sharing `xa²`/`xa⁴` between them, 2 extra multiplies
  total for both chains combined, not 4). mca confirmed the small-win
  prediction: latency 78.09->77.09 cyc (-1.3%), but throughput got
  *worse* (2.599->2.634 cyc/elem, +1.3%) -- a wash, not a clean win
  either way. Exhaustive sweep also showed a small real accuracy cost
  this time, unlike erf_poly's exact match: avg ulp 0.3106->0.3189 (a
  ~2.7% relative increase, max ulp unchanged at 109, still comfortably
  in budget either way). Reverted: with latency and throughput roughly
  canceling out and accuracy moving the wrong way, this doesn't clear
  the "speeds up without an accuracy penalty" or "improves accuracy
  without a perf penalty" bar either one — it's genuinely a third
  outcome, not a scaled-down version of erf_poly's clean win or
  acos_poly's clear loss. Together, these three same-restructuring
  attempts in one session (acos_poly: real speedup, real accuracy loss,
  rejected; erf_poly: real speedup, no accuracy cost, adopted; erfc's
  n/d: negligible speedup, small accuracy cost, rejected) confirm the
  crate's recurring lesson applies here too — a technique that works
  cleanly on one poly doesn't transfer automatically to a
  similar-looking one; each needs its own measurement.
  the asin mid-branch refit just above, extending `examples/tune.rs` with
  `atan_poly_c`. Turned out atan was already close to a strong local
  optimum for this coordinate-descent scheme (unlike asin's mid branch,
  where it wasn't): even a 39M-point grid (vs. the file's normal ~25k)
  converged to the same coefficients as the coarse grid, and the result
  was a small win on *every* axis at once — max ulp 19→18 and avg ulp
  improved too, for both atan and atan2 (atan2 calls atan directly, so it
  inherited the gain for free). No avg/max tradeoff this time, unlike
  asin's refit. Confirmed zero perf cost (mca bit-for-bit unchanged: atan
  57.09/1.410, atan2 57.17/1.467, exactly matching pre-refit). The
  "likely several ulp, needs a refit" framing this bullet used to have was
  already stale before this refit — atan's avg ulp (0.188) was already
  comfortably under the 0.5 budget, only max ulp had real room, and that
  room turned out to be small (1 ulp) once actually measured against a
  real coordinate-descent search rather than assumed.
  **Degree probe: not attempted.** Dropping a term from the Padé form is
  a separate, bigger change (changes the function's shape, not just its
  coefficients) — left for a future pass if the 1-ulp ceiling found here
  ever needs to move further.
- **acos/asin shared kernel's `acos_poly` refit — tried, no meaningful
  headroom found, not applied (2026-07-07).** Third use of the
  `examples/tune.rs` recipe this session, extended with `acos_poly_c`
  (scored as the whole `sqrt(1-x)*acos_poly(x)` formula, i.e. acos itself
  for x >= 0, not the bare poly — acos_poly is also asin's near-1 branch,
  so an improvement here would help both). Unlike asin's mid branch or
  atan_poly, this one is already essentially at a coordinate-descent local
  optimum: max ulp stayed at 3 and avg ulp moved by <0.3% (0.957→0.954),
  confirmed stable across a 10x grid-density increase (10k-point and
  100k-point grids landed on the same result). Not applied — a change
  this small isn't distinguishable from grid-sampling noise once it hits
  the real exhaustive sweep, and doesn't justify the verify/document/
  commit overhead. Kept `acos_poly_c` in tune.rs as working
  infrastructure (a correct, reusable tuner, even though this particular
  run found nothing) rather than reverting it, matching how exp2_c/log2_c/
  asin_mid_c/atan_poly_c are kept as permanent tuners in that file.
- **acos's negative-zero sign bug — found and fixed while re-investigating
  the same "acos near 1" territory the refit above covered, done, tested,
  kept (2026-07-07).** Went looking for an *algebraic* (not coefficient)
  fix for acos's residual, and ran the exhaustive sweep to characterize it
  precisely first -- found something the earlier quick-mode reading
  (0.50 avg / 4 max) had completely missed: **exhaustively, avg ulp was
  0.99 and max ulp was in the billions**, worst x = `-0.0`. `acos(-0.0)`
  was returning `-pi/2` instead of the correct `+pi/2` (acos is never
  negative, unlike odd functions where `-0.0` legitimately maps to a
  negative result). Root cause: `mulsign(y, x)` uses `x`'s raw sign *bit*,
  while the trailing `if x < 0.0 { PI }` uses a value comparison -- these
  agree for every finite x except exactly `-0.0` (bit says negative,
  value says "not less than zero"), so `mulsign` flipped `y`'s sign while
  the `+pi` correction never fired. This is precisely the class of bug
  this file's own "Negative-zero audit" bullet warned about in the
  abstract ("several select-based paths can silently lose -0.0") --
  first time in this session it was actually *found*, by exhaustive
  sweeping rather than by the audit itself. Fixed with `mulsign(y, x +
  0.0)`: `-0.0 + 0.0 = +0.0` exactly (IEEE754's defined round-to-nearest
  behavior for that one case), a no-op for every other input including
  genuinely negative ones. Verified exhaustive: avg ulp 0.99 → 0.4964
  (back in budget), max ulp billions → 4 (the same small residual the
  refit above already established can't be tightened further via
  coefficients). mca essentially unchanged (37.11 cyc latency both
  before and after, throughput 0.811 → 0.820, noise-level) -- a nearly
  free fix for a bug an order of magnitude worse than what quick-mode
  fuzzing alone would ever have surfaced. edgecheck extended with
  `acos(0)`/`acos(-0)` (via `check_known_1ulp`, since the 1-ulp gap from
  the "ideal" `FRAC_PI_2` there is acos_poly's own pre-existing,
  unrelated fit imprecision — present before this fix too, not
  introduced by it).
- **atan2's negative-zero sign bug — a direct follow-up to acos's, found
  by systematically checking every other `mulsign` call site in the file
  for the same disagreement shape, done, tested, kept (2026-07-07).**
  Grepped for all `mulsign` call sites after fixing acos, then tested
  each caller (asin, atan, atan2, erf) at zero/sign combinations against
  std rather than reasoning through the boolean logic by hand (atan2 in
  particular has genuinely intricate, *legitimate* signed-zero semantics
  per IEEE754/C99, unlike acos where any negative result is simply
  wrong — hand-tracing that risked missing a real distinction). asin,
  atan and erf all checked out clean (asin and atan don't reach a
  mulsign call site at the affected inputs; erf is a straightforward odd
  function with no companion domain-correction term to disagree with).
  atan2 didn't: `atan2(-0.0, +0.0)` came out `+0.0` instead of the
  IEEE754/C99-defined `-0.0`. Different root cause from acos's (not a
  mulsign-vs-comparison disagreement): when `x` is exactly `+0.0`, `base`
  degenerates to exactly `+0.0`, and the final `base + mulsign(...)`
  combines it with a correctly `-0.0`-signed correction term -- but
  IEEE754 addition of two *opposite*-signed zeros is defined to give
  `+0.0` regardless of operand order (only same-signed zeros, or a zero
  plus a genuinely nonzero value, preserve the expected sign), silently
  destroying the correction's sign. Every other zero/sign combination
  avoids this by construction: `x = -0.0` makes the correction's own
  `hpisignx` term flip sign too (turning the correction into a real
  nonzero `+-PI`, not a degenerate zero), and whenever `x` is genuinely
  nonzero, `base` is a real nonzero-ish angle rather than an exact zero,
  so the addition never lands on the opposite-sign-zero case. Fixed by
  skipping the addition entirely when `base` would be that exact `+0.0`
  (`nonzerox` already selects between the two shapes, so this moved an
  existing branch rather than adding one). Verified all 12 zero/sign
  input combinations bit-exact against std (the only inputs IEEE754/C99
  actually pin down exactly; general nonzero inputs were never meant to
  be bit-exact with std, only within budget, and stayed exactly at
  baseline: avg/max ulp 0.135/18 unchanged). mca confirmed zero perf
  cost (57.17 cyc / 1.467 cyc/elem, bit-for-bit unchanged) -- same
  computation, just restructured which branch does the addition.
  edgecheck extended with all 12 combinations. Worth noting as
  methodology: systematically re-checking every instance of a bug
  *pattern* (not just the one instance found) turned up a second, real
  bug with a different mechanism in the same afternoon.
- **sin/cos/sin_checked/cos_checked/tan's negative-zero sign bug — a
  broader continuation of the same systematic sweep, done, tested, kept
  (2026-07-07).** After acos and atan2 both turned up real signed-zero
  bugs from the same IEEE754 "opposite-signed-zero addition gives `+0.0`"
  mechanism, swept every odd-symmetric public function at `x = -0.0`
  against std rather than stopping at two instances. Found `sin(-0.0)`,
  `cos(-0.0)` (silently — `cos(-0.0) = +1.0` was already right, since
  it's a nonzero result unaffected by this bug class), `tan(-0.0)`, and
  `sin_checked(-0.0)` all wrong. Two distinct root causes behind one
  shared mechanism: `sinf_poly` (the poly core shared by all five
  functions) computes `fma(p, x3, x)`, and at `x = +-0.0`, `x3 = y*x`
  correctly carries `x`'s sign (`y = x*x` is always `+0.0`, and
  `+0.0 * x` doesn't flip sign) but `p` (the poly's fixed leading
  coefficient `c0` at `y=0`, i.e. sin's own curvature) is a negative
  constant, so `p*x3`'s sign is always the *opposite* of `x`'s at this
  one point — the `fma` then adds two exactly-zero values of opposite
  sign, which IEEE754 always resolves to `+0.0`, silently destroying it.
  First fix tried: an explicit `if x == 0.0 { x } else { r }` select.
  Worked (bit-exact), but mca showed a real, unwelcome throughput cost on
  every one of the five shared callers (sin +16.5%, cos +12.2%,
  sin_checked +4.9%, cos_checked +1.4%, tan +16.9% cyc/elem, latency +1
  cyc each) — the branch is inlined into every one of them, so its cost
  multiplies. Replaced with `r.copysign(x)` instead: for *every nonzero*
  `x` in this poly's domain, sin is odd and monotonic and the leading `x`
  term always dominates `p*x3` in magnitude (checked by hand at the
  domain edge too, `r` near `pi/2`: correction magnitude ~0.57 vs. `x`
  magnitude ~1.57, same sign, never crosses over), so `r`'s sign already
  equals `x`'s everywhere except this one singular point — `copysign` is
  a true no-op for every nonzero input, no accuracy or extra-rounding
  risk, and compiles to a single sign-copy instead of a compare+select.
  Bit-exact vs. the branchy version everywhere (confirmed via probe), and
  strictly cheaper on every one of the five mca rows (sin 1.191→1.151,
  cos 1.526→1.406, sin_checked 5.537→5.476, cos_checked 4.537→4.537 flat,
  tan 2.716→2.532 cyc/elem) — still not fully free vs. the pre-bugfix
  baseline (sin +12.6%, cos +3.4%, sin_checked +3.7%, cos_checked +1.4%,
  tan +8.9% cyc/elem throughput, latency flat except tan +2 cyc), but the
  smallest cost found for this correctness fix. `sin_checked` needed a
  *second*, separate guard: even after `sinf_poly` was fixed,
  `sin_checked(-0.0)` still failed, because `reduce_pi`'s own
  multi-term two_sum/two_prod error-compensation chain independently
  loses `x`'s sign somewhere internal to it (same IEEE754 rule, exact
  spot not traced — probed down to "reduce_pi's raw output is already
  `+0.0` by the time sinf_poly sees it" and stopped there rather than
  walking its ~15 intermediate two_sum/two_prod terms). Guarded at
  `sin_checked`'s own output with `if x == 0.0 { x } else { result }`
  instead — `cos_checked` needs no equivalent guard (`cos_checked(-0.0)
  = +1.0`, nonzero, unaffected). Verified bit-exact against std for all
  five functions at `x = -0.0`; exhaustive accuracy.rs sweeps for
  sin/cos/sin_checked/cos_checked all matched the pre-fix documented
  baseline exactly, bucket for bucket, confirming `copysign` is a true
  no-op away from the singular zero point. `tan (in-domain)`'s max ulp
  moved slightly (2967 -> 3057, avg unchanged at 0.331/0.3305) -- this is
  deep in the already-disclosed, budget-exempt domain-edge blowup zone
  (worst x ~1.318e7, right at the documented |x|<2^22*pi cliff), not a
  new qualitative regression, and consistent with `tan`'s residual `r`
  occasionally landing outside `sinf_poly`'s "leading term dominates"
  guarantee once the reduction itself is already degrading near the cliff.
  edgecheck extended with `sin(-0)`, `cos(-0)`, `tan(-0)`,
  `sin_checked(-0)`, `cos_checked(-0)`.
- **log1p/atanh/remainder's negative-zero sign bugs — the last three
  findings from the same broad `-0.0` sweep, done, tested, kept
  (2026-07-07).** `log1p(-0.0)`: `ln(u) + corr` adds two exactly-zero
  values of opposite sign at `x = +-0.0` (`ln(1.0)` is `+0.0`, but `corr`
  correctly carries x's sign there) -- identical mechanism to sinf_poly's
  bug, different call site. `atanh(-0.0)` traced to be a pure downstream
  consequence: `atanh(x) = 0.5*(log1p(x) - log1p(-x))` reuses log1p
  directly, so fixing log1p alone fixed atanh too, verified, no separate
  atanh change needed. `remainder(-0.0, y)`: `fma(-q, y, x)` hits the
  same opposite-sign-zero addition (`q` is `+-0.0` matching x/y's sign,
  so `-q*y` ends up opposite x's sign), but remainder's sign does *not*
  generally track x's sign for nonzero x (`remainder(2.0, 3.0) == -1.0`
  is a real IEEE remainder property, not a bug), so unlike sinf_poly/
  log1p, a blanket `copysign(x)` fix would be wrong here.
  Went through several fix attempts before landing on the right one,
  each measured rather than assumed:
  1. `copysign`/equivalent-bitwise-OR for log1p (branchless, valid since
     log1p is odd/monotonic): compiled fine, log1p's own mca unchanged,
     but `asinh`/`acosh` (log1p's only in-crate callers) latency roughly
     *doubled* (55.40->122.42, 110.47->130.27 cyc) for reasons not
     obviously related to the fix itself.
  2. A plain trailing `if x == 0.0 { x } else { normal }` select (same
     shape, no copysign): nearly identical elevated asinh/acosh numbers
     (120.99/127.75) -- ruled out "copysign specifically" as the cause.
  3. An early `if x == 0.0 { return x; }` guard at the very top (before
     computing anything else): this recovered asinh/acosh back to their
     original baseline (55.04/110.52) almost exactly -- but broke
     `cargo run --example mca` outright with `llvm-mca failed: ... found
     an invalid region end directive ... unable to find an active
     anonymous region`, traced to the early return causing LLVM to
     tail-duplicate the throughput harness's trailing marker `asm!` call
     across multiple per-element exit paths once inlined into a
     partially-scalarized unrolled loop (confirmed by isolating the
     change to `remainder` alone, and independently confirmed as a real,
     reproducible source change and not stale-build noise via a fully
     clean rebuild -- ruling out an earlier false alarm where a stray
     leftover `mca_target-*.s` file from a previous build had briefly
     made the *unmodified, committed* baseline itself misreport 122.42
     instead of the true 55.40, a second instance of this crate's known
     "manual before/after mca dump can be stale" pitfall). Reverted --
     an unusable option regardless of its latency win, since it breaks
     this crate's own measurement tooling.
  4. Settled on option 2's plain trailing select (`log_2`'s own
     established "compute the normal path unconditionally, select after,
     no early returns" idiom, chosen specifically so array loops keep
     auto-vectorizing) for both log1p and remainder. Compiles cleanly
     everywhere, `log1p`/`atanh`/`remainder`'s own mca numbers land
     within noise of their pre-fix baseline, but `asinh`/`acosh` still
     pay the elevated latency as a real, disclosed side effect --
     ultimately harmless in practice, since both already computed their
     own correct sign externally via `mulsign` (log1p is only ever
     called with a provably-`+0.0` argument in their flow, so this fix
     provides them literally zero benefit, pure inlining overhead) and
     their *throughput* -- the metric this crate's vectorization-first
     design actually prioritizes -- improved instead of regressing
     (asinh 9.625->7.716, acosh 6.839->6.588 cyc/elem). See readme.md's
     mca section for the full numbers.
  Verified bit-exact against std for all three functions at their `-0.0`
  inputs. Exhaustive accuracy.rs sweeps for log1p/atanh/asinh/acosh all
  matched their pre-fix documented baseline exactly; remainder has no
  exhaustive baseline (fuzz-only, dominated by its own pre-existing,
  already-disclosed `|x/y|`-large cancellation issue, unaffected by this
  change). edgecheck extended with `log1p(-0)`, `atanh(-0)`,
  `remainder(-0,3)`, `remainder(0,3)`.
- **powf(negative x) was always NaN — a much bigger bug than a sign
  quirk, found while closing out the "Negative-zero audit" idea below,
  done, tested, kept (2026-07-07).** After the sinf_poly/log1p/atanh/
  remainder fixes, ran one more sweep checking every remaining public
  function at `-0.0` (`expm1`, `sinh`, `sinh_throughput`, `tanh`,
  `asinh`, `atan`, `erf`, `cbrt`, `cbrt_accurate`, `ln`, `log10`, `hypot`,
  `powf`) — all clean except `powf(-0.0, 3.0)` (`+0.0` instead of
  `-0.0`). Widening the check to a couple of ordinary negative bases
  turned up something much bigger: `powf(-2.0, 3.0)` and `powf(-2.0,
  2.0)` were both `NaN`, not `-8.0`/`4.0`. Root cause: `powf` is
  `exp2_checked(log2(|x|)*y)`-shaped (well, `log2(x)*y` before this fix —
  no `abs` at all), and `exp2` of any real argument is always
  non-negative — this route has *no way* to ever produce a negative
  result, for any `y`, not just at the `-0.0` singularity. So this
  wasn't a narrow edge case: *the entire negative-base half of powf's
  domain* silently returned `NaN` instead of a well-defined answer,
  completely undocumented (no doc-comment caveat, no readme mention).
  Fixed by computing the magnitude on `|x|` (unchanged formula otherwise)
  and reapplying the sign for negative `x` only when `y` is an integer
  (even -> positive, odd -> negative, via `parity`, reusing
  `sin_checked`/`cos_checked`'s existing integer-parity helper — no new
  primitive needed) and `NaN` when `y` isn't an integer (correctly
  matches std: real roots of negative numbers to a non-integer power
  aren't representable). Two more real bugs fell out of the same
  investigation: (1) `pow(x, 0) = 1` for *any* `x`, even `0`/negative/
  `NaN`, is a dedicated IEEE754/C99 special case the log/exp2 formula
  can't derive on its own (`0*inf` and `NaN*0` both degrade to `NaN`) —
  `powf(0.0, 0.0)` and `powf(f32::NAN, 0.0)` were both `NaN` instead of
  `1.0`, fixed with a trailing `if y == 0.0 { 1.0 } else { r }` override
  (compute everything unconditionally first, select last — the same
  no-early-return idiom as the log1p/remainder fix just above, and for
  the same reason: `cargo run --example mca` is the check that would
  catch an early-return regression here, not `cargo build`/`test`).
  (2) The negative-base sign check itself first used `x < 0.0`
  (value-based) instead of `x.is_sign_negative()` (bit-based) — the
  identical class of bug as acos's `-0.0` fix earlier this session,
  since `-0.0 < 0.0` is `false`, silently routing `powf(-0.0, 3.0)`
  through the *positive* branch. Also widened `accuracy.rs`'s own
  `powf` domain filter, which had been narrowed to `x > 0.0` only — the
  same "filter dodges the bug instead of exercising it" pattern already
  found in erf/erfc's filters this session — to `x != 0.0`. Verified
  bit-exact against std for every case above (`powf(-2,3)`, `powf(-2,2)`,
  `powf(-2,3.5)==NaN`, `powf(0,0)`, `powf(-0,0)`, `powf(nan,0)`,
  `powf(-0,3)`, `powf(-0,2)`, `powf(-0,-1)`). Fuzz sampling almost never
  lands on an exact integer `y`, so also ran a dedicated integer-`y`
  spot sweep (401×39 grid over `x`/`y`) to check the newly-reachable
  negative-base path specifically for an accuracy cliff — none found
  (max ulp 92, same order of magnitude as the existing positive-base
  budget, since the fix only adds a sign correction after the existing
  log_2/exp2_checked machinery, not a new computation path). Not fixed:
  `powf(-1.0, ±inf)` — IEEE754 special-cases this to `1.0`, but this
  crate's formula gives `NaN` (`log2(1.0)*inf` degrades to `0*inf`); left
  as a documented gap rather than adding more special-case logic for a
  base/exponent combination unlikely to matter in this crate's actual
  use cases. mca cost is small: 97.88->98.03 cyc latency (+0.2%),
  3.788->3.898 cyc/elem throughput (+2.9%) — by far the cheapest of this
  session's correctness fixes relative to its scope (a whole missing
  half-domain, not just a sign at one point).
- **remainder/hypot/atan2's infinity-handling gaps — a direct follow-up
  to powf's fix, done, tested, kept (2026-07-07).** Extended the same
  sweep to `+-inf`/`NaN` combinations for every remaining two-arg
  function, on the theory that if powf had a whole undocumented missing
  half-domain, other multi-arg functions might too. Found three:
  `remainder(finite x, +-inf)` was `NaN` instead of `x` -- `q` rounds to
  exactly `0.0` for any finite `x` (correct), but `fma(-q, y, x)` then
  multiplies that zero by an *infinite* `y`, and `0*inf` is `NaN`,
  destroying the intended "no reduction happened" no-op (IEEE754/C99:
  `remainder(x, +-inf) = x` for finite `x`). Fixed with a trailing
  `if y.is_infinite() && x.is_finite() { x } else { r }` select -- `x`
  itself infinite/NaN still correctly falls through to `NaN` (matches
  std's `remainder(inf, y) = NaN`), verified unaffected.
  `hypot(+-inf, NaN)`/`hypot(NaN, +-inf)` were `NaN` instead of `+inf` --
  IEEE754/C99 special-cases infinity to "win" over NaN specifically in
  hypot (unlike almost every other function: `sqrt(x²+y²)` genuinely
  diverges to `+inf` as either argument grows without bound, regardless
  of what the other argument is), but the naive `fma(x,x,y*y).sqrt()`
  can't reach that once an argument actually *is* NaN (`inf*inf +
  NaN*NaN` degrades to `NaN`). Fixed with a trailing
  `if x.is_infinite() || y.is_infinite() { INFINITY } else { normal }`
  select -- explicitly *not* the same thing as hypot's already-documented
  finite-overflow tradeoff (that's about large-but-finite `x`/`y`
  overflowing `x*x`/`y*y`; this is a literally-infinite argument, a
  different case the doc comment didn't cover). `atan2(+-inf, +-inf)`
  was `NaN` instead of a defined `+-pi/4`/`+-3pi/4` by quadrant -- `y/x`
  is `inf/inf` (`NaN`) there, so the general `atan`-based formula simply
  has no ratio to compute with at true infinity (there isn't one, unlike
  a large-but-finite y/x that approaches some limit); IEEE754/C99 instead
  define a fixed convention independent of any limit. Fixed the same way,
  a trailing select gated on `x.is_infinite() && y.is_infinite()`, added
  `FRAC_PI_4` as a new local const alongside the existing `FRAC_PI_2`.
  All three follow the same "compute everything unconditionally, select
  last, no early returns" idiom as log1p/remainder/powf above (checked
  `cargo run --example mca` compiles clean after each, not just
  `cargo build`/`test`, per that same lesson). Verified bit-exact against
  std for every case (`remainder(3,inf)`, `remainder(-3,inf)`,
  `remainder(inf,3)==NaN`, `hypot(inf,nan)`, `hypot(nan,inf)`, all 4
  `atan2(+-inf,+-inf)` quadrant combinations). `hypot(1e30,1e30)`
  (finite-overflow) deliberately left alone -- that's the pre-existing,
  already-documented tradeoff, a different case, not something this pass
  touched. Fuzz accuracy sweeps for hypot/atan2 unchanged from their
  documented baselines. mca cost negligible: atan2 completely unchanged
  (57.17/1.467, bit-for-bit — the new select only ever fires on inputs
  the old formula already couldn't handle), hypot +0.5%/+0.4%
  (21.00->21.11 cyc, 0.763->0.766 cyc/elem), remainder unchanged
  (33.02/0.647).
- **erf's tail branch (`erf_poly`) refit — tried, no meaningful headroom
  found, not applied (2026-07-07).** Fourth use of the tuning recipe,
  extended with `erf_tail_c` (scored as the whole `mulsign(1.0 -
  exp2(erf_poly(xa)), x)` formula for xa in [0.28, 10], exactly where this
  branch is used post the domain-bound fix earlier this session). Needed
  a real ground-truth reference this crate didn't have wired into tune.rs
  yet — `sleef::f64::erf_u10` (the scalar, non-SIMD form; the module path
  is private, the flat `erf_u10` re-export at the crate root is the way
  in, same as `sleef::f64x::erf_u10` accuracy.rs already uses for its
  vectorized version) — no nightly/portable_simd needed for the scalar
  form, unlike accuracy.rs's own sleef usage. Same result shape as
  acos_poly just above: max ulp unchanged (4→4), avg ulp moved <0.3%
  (0.386→0.385), confirmed stable across a 10x grid-density increase.
  Not applied, kept as infrastructure for the same reasons as acos_poly.
  Between this and acos_poly both landing at "already optimal," the
  pattern suggests this session's earlier refits (asin, atan) found real
  headroom specifically because those coefficients hadn't been touched
  since the original C port, while acos_poly/erf_poly may have already
  been through more careful tuning historically — worth checking git
  blame / original-source provenance before spending more tuning passes
  on functions without first estimating whether they're likely to have
  slack, rather than tuning everything uniformly.
- **erfc's n/d rational refit — done, tested, kept (2026-07-07), fifth use
  of the tuning recipe and the second real win (after atan_poly), landing
  in the "has real headroom" bucket the hypothesis above predicted
  wrong** (erfc is also a straight C port like acos/erf, yet had real
  slack — the provenance heuristic isn't a reliable predictor by itself).
  Extended tune.rs with `erfc_c` (whole formula, all 8 named coefficients
  free, the two Horner chains' trailing `+1.0` left fixed to match the
  shipped structure) and a scalar sleef reference (`erfc_u15`, same
  crate-root re-export pattern as `erf_u10`). Tuner found max 96→81 on
  its grid with avg getting slightly worse (0.295→0.308) — same
  max/avg tradeoff shape as asin's mid-branch refit. Applied and verified
  on the real exhaustive-equivalent sweep: max ulp 123→109 (~11%), avg
  ulp 0.297→0.311 (still comfortably under the 0.5 budget). mca
  bit-for-bit unchanged (78.09 cyc / 2.599 cyc/elem) — zero perf cost,
  same instructions, only the 8 literal constants differ (one of the 8,
  the leading `n` coefficient, the tuner didn't touch at all).
- **erfc's `w = if x<0 {2.0} else {0.0}` is a redundant second select —
  done, tested, kept (2026-07-07).** Found on a fresh look at erfc's
  current (post-refit) source while re-checking the "shared
  subexpression" family of ideas. `w` and `z` (`if x<0 {-1.0} else
  {1.0}`) both derive from the exact same `x < 0.0` condition, and
  `w = 1.0 - z` holds as an exact identity for *both* outcomes (`1.0 -
  1.0 = 0.0`, `1.0 - (-1.0) = 2.0`, both exact, no rounding) — not a
  special case at any particular `x`, a genuine algebraic redundancy:
  the second compare+select was recomputing information the first select
  already fully determined. Replaced with one subtract. Verified via a
  200M-sample bit-exact comparison against the old two-select form
  (stepping by a stride coprime to 2^32 for broad, distinct coverage —
  a true 2^32 loop timed out as a plain sequential scalar test, and this
  is a pure identity that doesn't need literal exhaustiveness the way an
  approximation's accuracy claim would): zero mismatches. mca showed
  *zero* change (78.09 cyc / 2.599 cyc/elem, exactly baseline) despite
  `--emit=asm` showing genuinely different codegen (different register
  allocation, not byte-identical like the atan2 `bothzero` no-op) — a
  real source simplification that happened to land on an equal-cost
  scheduling point rather than a faster one. Kept anyway for the simpler
  source (one fewer branch, one fewer named constant pair) at proven zero
  risk and zero regression, matching this crate's precedent of keeping
  provably-free simplifications even when the runtime win doesn't
  materialize (see the sin_checked sign-flip entry elsewhere in this
  file for another "kept for the simpler code" case).
- **erf's near-zero Padé branch refit — tried, no meaningful headroom
  found, not applied (2026-07-07).** Sixth use of the recipe, extended
  with `erf_near0_c` (the `numer/denom` formula, |x| < 0.28, the tail
  branch above that already handled separately). Same negligible-result
  shape as acos_poly/erf_tail: max ulp unchanged (3→3), avg ulp moved
  <0.3% (0.637→0.636) — confirmed stable when the grid was widened from
  ~1M to ~10.5M points (a background run at that density stalled for
  reasons unrelated to convergence — likely thermal/scheduling
  contention from an earlier long-running foreground command, not the
  search itself — and was killed after ~80s without finishing; the
  smaller grid's already-consistent result across three "no headroom"
  cases now (this one, acos_poly, erf_tail) was treated as sufficient
  without waiting on it). Not applied, kept `erf_near0_c` as
  infrastructure. This session's refit scorecard is now asin/atan/erfc
  finding real headroom vs. acos_poly/erf_tail/erf_near0 finding none —
  three wins, three no-ops, no clean predictor found yet for which is
  which before actually running the tuner.
- **atan without reciprocal-select**: `y = if a < 1 {a} else {1/a}` then a
  conditional π/2 flip — fine already; alternatively fit atan on [0, ∞) via
  t = x/(1+|x|) rational reduction, one division, no select chain.
- **atan2 accuracy — tried, measured, reverted (2026-07-07).** `d = y/x;
  e = fma(-d, x, y)/x; corr = e/fma(d,d,1.0)` (atan'(d)*e, a first-order
  Taylor correction for the division residual) added to `atan(d)` outside
  atan's own poly (couldn't inject it into atan_poly's internal range
  reduction without a bigger refactor). First measurement looked like a
  huge win (avg ulp 0.136→0.0134, ~10x) -- until the max-ulp column
  showed a `u64::MAX` sentinel, this crate's `ulp_diff` convention for "my
  result is NaN, the reference isn't." Found the cause by hand: when `d`
  itself overflows to `inf` (a case `atan(d)` alone already handles
  correctly, returning π/2), `e` becomes `inf - inf`-shaped and `corr`
  comes out NaN, contaminating an otherwise-correct result. Guarded with
  `if corr.is_finite() { corr } else { 0.0 }`, matching the pattern used
  for every other correction term this session (log1p, acosh, asinh,
  atan2 itself already has similar shape). Once guarded, the "10x
  improvement" evaporated entirely: three repeat runs landed at
  0.136/0.1362/0.1367 avg ulp, max ulp 18 -- statistically indistinguishable
  from the unmodified baseline (0.1355/18, from the atan_poly refit).
  Confirms the idea's own original caveat ("only matters once atan's own
  poly is refit under budget") -- atan's residual poly-fit error already
  dominates atan2's total error, so removing the division's rounding from
  the budget doesn't move the needle at all. mca confirmed a real cost for
  the zero benefit: 57.17→65.14 cyc latency (+13.9%), 1.467→2.096 cyc/elem
  throughput (+42.9%). Reverted. The NaN-guard bug itself is a useful
  general lesson even though the idea didn't pan out: a "10x win" in a
  quick accuracy run is worth being suspicious of, not just accepting,
  when the max-ulp column shows a sentinel value -- check *why* before
  trusting the average.

## exp-family elementary (sinh/cosh/tanh) — currently 2× exp each

- **One exp, one reciprocal** — done, tested, adopted as a separate opt-in
  tier rather than a straight replacement (2026-07-07). e = exp(x);
  exp(−x) = 1/e. Accuracy hit was negligible in practice (no refit needed:
  cosh avg/max ulp 0.275/63 → 0.272/64, comfortably inside budget). But this
  was the first idea in the whole IDEAS.md-testing series where mca showed
  latency and throughput moving in *opposite* directions with no accuracy
  angle at all: throughput improved hugely (sinh/cosh 2.083→1.279 cyc/elem,
  -38.6%, matching the "large throughput win" prediction exactly), but
  latency got measurably *worse* (48.02→58.00 cyc, +20.8%) — the division
  can't start until `exp(x)` finishes, whereas the original's two
  independent `exp(x)`/`exp(-x)` calls run in parallel on a CPU with enough
  ports, so serial latency was better *before* this change. Unlike every
  prior "op reshuffle didn't help" entry in this file, this one isn't
  fixable by reordering — the division is inherently dependent on `e`.
  Flagged to Jodie as a genuine tradeoff (not obvious which axis the crate
  should prioritize); resolved as **both**: `sinh`/`cosh` keep the
  lower-latency parallel-exp form as the default, `sinh_throughput`/
  `cosh_throughput` are new public functions for array-heavy callers,
  mirroring the sin/sin_checked and exp2/exp2_checked tiering already used
  elsewhere (see readme for the full writeup). tanh wasn't touched by this
  idea — it already only calls `exp` once (via `exp(2x)`), so the
  two-exp-calls framing doesn't apply there.
- **tanh cancellation fix — done, tested, kept (2026-07-07).** tanh =
  expm1(2x)/(expm1(2x)+2), one division, reuses expm1's existing Pade
  small-x branch. Like log1p just above, this wasn't a marginal miss: the
  same crate-wide fuzz sweep that found log1p at 172M avg ulp found tanh at
  **336,982,029 avg ulp / 868,220,899 max ulp** (worst x ≈ 8.9e-8) from the
  old `1.0 - 2.0/(exp2(2x·log2e)+1.0)` — the `1.0 - (~1.0)` step lost
  essentially all precision for small x, exactly the documented flaw.
  Considered a cheaper direct fix first (resequence the existing exp2 call
  without a full expm1 delegation) but ruled it out analytically: `y - 1`
  for `y = exp2(2x·log2e)` suffers the *identical* collapse-to-1.0 problem
  `1.0 + x` had in log1p (y itself rounds to exactly 1.0 for small x before
  any subtraction happens), so there's no cheap resequencing here the way
  Sterbenz gave log1p one — an accurate small-x numerator fundamentally
  needs Pade-style machinery, which is exactly what expm1 already has.
  Verified exhaustive (2.2B in-domain samples, accuracy.rs's existing
  `|2x·log2e| < 126` filter matching tanh's documented halved-range
  contract): avg/max ulp now 0.1438/6 (worst x ≈ -0.24, an ordinary
  interior point, not a cancellation site — the residual is inherited from
  expm1's own up-to-63-ulp worst case, not a new bug). edgecheck's
  `tanh(0)` and a manual probe of tanh at ±inf/NaN/44/50/100 (all already
  NaN or garbage before this change, per the unchecked-domain contract,
  and unchanged after) confirmed no regression outside the documented safe
  range. mca cost is real and bigger than log1p's: 58.00→87.64 cyc latency
  (+51%), 1.330→1.793 cyc/elem throughput (+35%) — expm1 unconditionally
  computes both its Pade and full-exp branches, so tanh now pays for a
  whole extra exp-equivalent call it didn't before. Kept anyway, same
  reasoning as log1p: 337-million-ulp average is a correctness bug in a
  function whose entire useful range centers on zero (tanh's most common
  use, e.g. activations near 0), not a speed/accuracy trade to weigh.
- **sinh small-x — done, tested, kept (2026-07-07).** Same bug class as
  log1p/tanh, found by the same fuzz sweep: sinh 333,732,587 avg ulp /
  864,026,619 max ulp (worst x ≈ 5.96e-8), sinh_throughput
  333,753,185/864,359,817 -- `exp(x)`/`exp(-x)` (or `e`/`1/e`) both round to
  ~1 for small x, so the subtraction cancels almost everything. Fixed with
  a genuine `x + x³/6` Taylor branch as sketched, extended to 4 terms (x +
  x³/6 + x⁵/120 + x⁷/5040 -- exact rational coefficients, not a numerical
  fit, since sinh is entire and Taylor is trivially available) after
  checking the 2-term version analytically first: at the |x|<0.5 select
  boundary, 2 terms alone leave ~5e-4 relative error (~4000 ulp, nowhere
  near budget), while 4 terms leave ~1e-8 relative (~0.1 ulp, comfortable
  margin) -- the "(one fma)" in the original idea undersold what's actually
  needed at that boundary. Shared the same `sinh_small` helper between
  `sinh` and `sinh_throughput` (both had the identical bug). Verified
  exhaustive: sinh now 0.2905 avg / 64 max ulp, sinh_throughput 0.2925/64 --
  both land almost exactly on cosh's own baseline (0.2745/64, cosh never
  had this bug, both add instead of subtract), confirming the residual
  max-64 is inherited from exp's own accuracy near its domain edge (worst x
  ≈ 86, right at the unchecked-exp2 boundary), not a new bug. mca cost, in
  the same range as tanh's fix: sinh 48.02→77.02 cyc latency (+60%),
  2.083→2.277 cyc/elem throughput (+9.3%); sinh_throughput 58.00→81.02 cyc
  (+40%), 1.279→1.545 cyc/elem (+21%) -- both now compute a Taylor branch
  unconditionally alongside the existing exp-based one. Kept for the same
  reason as log1p/tanh: 333-million-ulp average is a correctness bug in
  sinh's most-used range (near zero), not a trade.
- **Direct minimax tanh**: on [0, ~9.02] tanh saturates; a rational P/Q in x²
  (odd) with clamp — no exp at all, one division. Likely the fastest shape;
  fit feasibility at ≤0.5 avg ulp needs checking (tanh rationals are
  well-behaved).

## erf / erfc / powf / hypot / remainder / log1p-dependents

- **powf, double-float log**: exp2(log_2(x)·y) amplifies log_2's error by
  y·2^… — the dominant powf error. Compute log2(x) as df (hi, lo) — log_2
  already has the pieces pre-collapse — multiply by y in df (two_prod + fma),
  split into int + frac, feed exp2's poly with the lo word folded in as
  `result·(1 + lo·ln2)` ≈ `fma(result, lo*LN2, result)`. This is the standard
  <1-ulp powf shape; maybe +6–8 fmas, transforms powf from ~hundreds of ulp
  (for large y) to budget-compliant.
- **powf integer-y fast path is NOT vectorizable as a branch** — skip; but the
  df version above covers those cases accurately anyway.
- **erfc tail — done, tested, kept (2026-07-07).** Swapped the interior
  `exp(-xa*xa)` for `exp2_checked(-xa*xa*LOG2_E)`. This wasn't a subtle
  ulp gap: `erfc(9.5)` was `inf`, `erfc(10)` was `9.2e26`, `erfc(-9.5)`
  was `-inf`, `erfc(-10)` was `-9.2e26` — all should be a clean near-0 (or
  near-2 for negative x) — confirmed by hand before and after the fix
  (git-stash A/B), not just inferred from the doc comment. The crate's own
  fuzz-mode accuracy sweep never caught this: `examples/accuracy.rs`'s
  `erfc_domain` filter had been narrowed to `|x| < 9.3` specifically to
  dodge the bug rather than exercise it (widened back to the full `<= 10`
  clamp now that it's fixed). Kept the `xa` clamp itself — it still
  protects the degree-4 rational polynomial (n/d) from overflowing at
  genuinely huge x, a separate concern from the exponent's domain.
  edgecheck extended with `check_finite` at the four previously-broken
  points. Real mca cost: 67.09→78.09 cyc latency (+16.4%), 2.097→2.599
  cyc/elem throughput (+23.9%) — comparable to log1p's fix. Kept for the
  same reasoning as the small-x cancellation fixes above: `inf` where the
  answer is a tiny positive number is a correctness bug, not a tuning gap.
- **erf's tail, a worse bug than erfc's — done, tested, kept (2026-07-07),
  found while double-checking accuracy.rs's own domain filters after the
  erfc fix above (same shape of bug: a filter narrowed to dodge an issue
  rather than exercise it, this time `erf_domain = |x| x.abs() < 6.0`
  with a doc comment claiming it was safe past there).** Checked by hand
  instead of trusting the comment: `erf(50) = -1.02e17`, `erf(100)` and
  `erf(±inf) = NaN`, all of which should be `±1`. Root cause is *not* the
  same as erfc's: `erf_poly` (the tail branch's degree-6 polynomial,
  evaluated directly on `|x|` with no bound at all) has a positive leading
  coefficient, so instead of staying deeply negative as `|x|` grows (which
  is what keeps `exp2(erf_poly(x)) ≈ 0`, giving the correct saturation to
  `±1`), it eventually turns around and grows to `+∞`
  (`erf_poly(9) ≈ -92`, `erf_poly(20) ≈ +8698`, confirmed by direct
  computation, not guessed) — so `exp2_checked` alone wouldn't have fixed
  this the way it fixed erfc; the polynomial itself needed bounding.
  Fixed by clamping `|x|` to 10 before `erf_poly` — the same bound erfc's
  own clamp already uses, and confirmed by direct computation that
  `erf_poly(10) ≈ -83.8` stays safely within even the *unchecked* exp2's
  domain — then swapped to `exp2_checked` anyway as cheap insurance now
  that the input is bounded, matching erfc's fix. That swap turned out to
  be a free accuracy bonus even within the previously-tested `|x|<6`
  range: `exp2_checked` is itself more accurate than plain `exp2`
  (established earlier this session), moving erf's own avg ulp from
  0.631 (over budget) to 0.319 (comfortably under) as a side effect of a
  bug fix, not a separate optimization. Widened accuracy.rs's own
  `erf_domain` filter to `everywhere` now that the bug is gone. edgecheck
  extended with `erf(50)/(-50)/(±inf)/(nan)`. mca cost: 88.02→102.74 cyc
  latency (+16.7%), 2.101→3.163 cyc/elem throughput (+50.6%) — kept for
  the same reasoning as every other correctness fix this session.
- **erf/erfc joint refit** with f32-quantized coefficients against the actual
  budget (these are relative-1e-6 C ports, same caveat as atan) — erf's
  ~0.32 avg ulp residual after the bound fix above is this kind of
  fit-tightness gap, not a further algebraic trick.
- **hypot**: keep naive per the crate's stated tradeoff; optional
  `hypot_checked` via max-exponent scaling (`vgetexpps` or bit trick) if
  wanted — two extra multiplies, no division.
- **remainder**: q = x/y rounding error dominates; a `two_prod(q, y)` +
  corrected subtract (3 extra ops) extends the reliable |x/y| range
  substantially — same trick family as reduce_pi, tiny cost.

## Measurement / meta

- **Drive `run_descent` with accuracy.rs's exhaustive sweep** (and avg+max
  ulp as the objective, not the bits-diff sum over random samples) — turn it
  into the standard last-mile tuner for every refit above.
- **CI asm check**: a script that compiles a loop over each public fn with
  `--emit=asm` and greps for `vfmadd.*ymm` / absence of scalar call — locks in
  the autovectorization contract so future edits can't silently break it
  (the scale-param duplication incident shows how easy that is).
- **mca targets per idea**: extend examples/mca_target.rs with one entry per
  candidate (atanh-log, table-exp2, sincos, one-exp-sinh…) so throughput
  /latency deltas are one command away before committing to accuracy work.

---

# Part 2: the rest of the kitchen sink

## Polynomial evaluation schemes (applies to every poly here)

- **Estrin/Horner hybrid autotuning**: for each poly, enumerate all reasonable
  bracketings (pure Horner, pure Estrin, 2nd-order Horner, the 3-balanced-pair
  scheme exp2 already uses) and mca every one. The optimum depends on whether
  the surrounding function is latency- or throughput-bound, so the answer can
  differ between e.g. `sin` (called in reduction-dominated context) and
  `log_2` (poly-dominated). This is mechanical — scriptable.
- **Second-order Horner** (process two coefficients per step in two parallel
  chains on y = x², merge at the end): halves the serial fma depth of Horner
  at zero extra multiplies for odd/even polys — sinf_poly is already odd, the
  log poly is not but can be split even/odd for the same effect.
- **Regroup to shrink live ranges**: the deg-9 log poly holds 10 constants +
  s, s², s⁴ live simultaneously; in a ymm/zmm loop that's real register
  pressure and potential spills. A regrouping that retires coefficients
  earlier (deeper but narrower) can win *in the loop* even when mca on the
  isolated kernel says otherwise. Check actual loop asm for spills.
- **Compensated Horner (EFT-based)**: error-free-transform Horner gives
  near-double precision for ~3× the fmas. Never worth it for a whole poly
  here, but a *single* compensated final step (two_prod on the last multiply,
  fold the error term) is 2 extra ops and buys most of the benefit — a
  general-purpose "last half ulp" tool where the final rounding dominates.
- **Coefficient representability tricks**: rescale the reduced variable by a
  power of two so more coefficients land exactly on representable f32s (or on
  values whose products with common inputs round nicely). Also: force c0 to
  an exactly-representable value (e.g. log2e's nearest f32) and refit the
  rest around it, so the dominant term contributes no coefficient error.
- **Rational everything**: given the idle-divider finding, systematically
  revisit each poly of degree ≥ 7 as a rational P/Q with small degrees
  (deg 3/3 ≈ deg 7 accuracy for smooth functions, often better for functions
  with poles/asymptotes like tan/tanh/atan). One division each. This could
  shrink log_2, sinf_poly, erf's tail poly, and acos_poly. Caveat from
  memory: divider-is-free is *uarch-specific* (measured on this CPU) — a Zen
  or E-core target may want the poly forms back, so keep both behind cfg or
  as `_rat` variants.
- **Chebyshev economization as refit seed**: lolremez minimax fits sometimes
  land in bad local basins for the subsequent f32 quantization; starting the
  quantized descent from an economized-Chebyshev fit gives a different basin
  for free. Cheap to try per poly.

## Range reduction, the deep end (sin/cos/tan)

- **Exact vectorizable Payne–Hanek for f32**: f32 only has ~256 exponents, so
  full PH is *small* here. Treat x as (mantissa u32) × 2^e, multiply the
  24-bit mantissa by 3 pre-selected 32-bit windows of 2/π (window index =
  e >> 5-ish, from a 256-entry table) using widening 32×32→64 multiplies
  (`vpmuludq` — fully vectorizable). ~12–18 integer ops, gives the *exactly*
  reduced argument for every finite f32, no cliff, no gradual degradation,
  no double-float chain. Table access is a gather (vpgatherdd) or, since
  consecutive lanes usually share exponents, sometimes a broadcast. This
  could plausibly beat sin_checked's current double-float reduction on
  *both* speed and accuracy — the single most interesting experiment in
  this file. Compare against reduce_pi via mca + exhaustive sweep.
- **Hybrid tiering**: fast single-word reduction for |x| < 2^22·π, PH only
  beyond — but branchless means computing both and blending, so this only
  pays if PH is expensive. If exact PH lands near reduce_pi's cost, drop the
  tiering entirely and delete reduce_pi (simpler *and* better).
- **Cody–Waite word-count autotuner**: a script that, given a target |x|
  bound, emits the minimal PI_A..X split with maximal trailing zeros and
  verifies exactness of each q·PI_i product over the bound. Turns "how many
  words do I need" from an art into a lookup, and lets the fast path's
  documented ~1.3e7 cliff be dialed per-application.
- **Reduction constant with trailing-zero-aware rounding**: PI_A..D style
  splits can be re-optimized so the *last* word absorbs the rounding
  optimally (minimize worst-case residual error, not just word-wise
  truncation). Small refit, occasionally worth 0.05 avg ulp near zeros.

## New function candidates (cheap wins with existing machinery)

- **sinpi / cospi / tanpi**: sin(πx) with the reduction done on x *before*
  multiplying by π — reduction is `x - round(x)`, exact in one fma+round, no
  Cody–Waite, no cliff, at any magnitude. Faster than `sin` and *more*
  accurate, for every caller whose angles are naturally in turns or degrees
  (graphics, FFT twiddles, ML positional encodings). The poly refits to
  sin(πr) on [-1/2, 1/2]. Probably the best accuracy-per-effort ratio in
  this whole file. `sind/cosd` (degrees) fall out the same way via /180.
- **sincos** (returns both): shared reduction, both polys — see Part 1.
- **sigmoid / logistic** `1/(1+e^-x)`: one exp2 + one division (idle
  divider), plus a select-free saturation story. Enormous demand (ML), and
  the naive user composition `1.0/(1.0+exp(-x))` hits exp's unchecked-domain
  garbage for |x| > ~88 — a curated version fixes that with the exp2_checked
  clamp shape.
- **softplus** `ln(1+e^x)`: log1p∘exp fused with the overflow-safe identity
  softplus(x) = max(x,0) + log1p(exp(-|x|)) — all selects, vectorizes.
- **rsqrt (full precision)**: hardware `vrsqrtps` seed (via
  `1.0/x.sqrt()`-pattern LLVM won't produce; needs intrinsic or std::simd)
  + one NR step reaches ~0.5 ulp is *not* quite true (NR from 12-bit seed
  gives ~24 bits but max ulp ~2–3); either a second compensated step or
  just `1.0/x.sqrt()` (sqrt + div, both on the idle divider!) which is
  correctly rounded to ≤1 ulp already. Measure which wins here.
- **erfinv**: odd poly in x against w = -ln(1-x²) with a two-range select
  (the classic Giles fit) — all fma + one ln call, vectorizes, useful for
  gaussian sampling.
- **norm_cdf / norm_pdf**: trivially from erf/exp once those are in budget;
  quant users care.
- **powi / pow2i**: exact-exponent manipulations (bit tricks only) for
  integer powers — no poly at all for pow2i; powi via square-and-multiply
  unrolled for small fixed n (const generic `powi::<N>`), fully vectorizable
  since N is compile-time.
- **exp10**: one constant away from exp2; free with the exp refit.
- **fract / modf / remquo helpers**: the crate builds these internally
  (floor/round tricks) — exporting polished versions costs nothing.

## Bit-trick and special-value catalog (micro but cumulative)

- **hypot one-liner**: `(x*x + y*y).sqrt()` → `fma(x, x, y*y).sqrt()`. One
  rounding fewer, literally free (same op count), strictly more accurate.
- **NaN-propagation via `x + x`** is used in cbrt/log; audit every function's
  special-path select for the cheapest NaN/inf-correct idiom (`x * x` for
  log's +inf/nan, `x + x` preserves sign of inf, min/max NaN suppression as
  in erfc's comment). Write the catalog down in a doc comment so future
  functions pick from it instead of rediscovering.
- **Sign handling audit**: `mulsign` (xor) beats copysign beats multiply;
  several ported functions still apply signs via `* z` multiplies (erfc's
  `fma(y, z, w)`) that could be sign-bit xors, freeing mul-port pressure.
- **Integer compares on float bits**: for positive floats, IEEE bits are
  monotonic — `ax < 0x2380_0000`-style integer compares (already used in
  cbrt_accurate) can replace float compares elsewhere and sometimes fuse
  with existing bit extractions; on some uarchs int-vector compare + blend
  has better port distribution than vcmpps.
- **Exponent extraction without shifts**: `(float)(bits >> 23)` via the
  subtract-magic-constant trick (used implicitly in log_2's `e` computation)
  vs. cvtdq2ps — pick per port pressure.
- **Negative-zero audit**: exhaustive edge sweep (edgecheck.rs) should assert
  sign of zero results for every function (sin(-0.0) = -0.0 etc.); several
  select-based paths can silently lose -0.0 and it's free to preserve with
  the right idiom order.
- **Denormal policy flag**: a crate feature `assume-ftz` that compiles out
  the tiny-input rescale paths (log_2's `tiny`, cbrt's `tiny`) for callers
  who run with FTZ/DAZ set anyway (games, audio) — measurable width savings
  in log/cbrt, zero accuracy change under their runtime settings.

## Type-level API ideas (make invalid states unrepresentable)

- **Domain-witness newtypes**: `Normal(f32)` / `Reduced(f32)` /
  `Positive(f32)` zero-cost wrappers with unsafe (or checked) constructors,
  so `log_2_normal`, `exp2` (unchecked), `sinf_poly` etc. take a *witness*
  instead of a doc-comment contract. The fast/unchecked tier becomes
  safe-by-construction, callers who already know their range (the whole
  point of the fast tier) encode it once, and `#[doc(hidden)]` hacks for
  mca_target go away. Compiles to nothing, vectorizes identically.
- **Consistent tier naming**: every function gets the same suffix scheme
  (`foo` = fast + documented domain, `foo_checked` = full-range) — sin/cos
  and exp2 already do this; cbrt inverts it (`cbrt` is checked,
  `cbrt_normal` is fast). Pick one convention, apply everywhere, table in
  the readme.
- **Slice entry points with dispatch hoisted**: `pub fn sin_slice(xs: &[f32],
  out: &mut [f32])` — gives a place to (a) guarantee the vectorizable loop
  shape, (b) do runtime `target_feature` dispatch *once per call* instead of
  never, (c) later hand-tune with std::simd without touching scalar API.
- **std::simd mirror crate/module**: portable-SIMD explicit versions of each
  kernel serve as (a) the codegen ground truth autovectorization is judged
  against, (b) a guaranteed-vector API for users whose loop shapes the
  autovectorizer fumbles. The scalar kernels here are already written in
  perfect SIMD style — porting is mostly mechanical `f32 → Simd<f32, N>`.

## Per-function stragglers

- **log_2**: split P(s) even/odd for second-order Horner (see above); probe
  whether `koff` can ride in the poly's c0 instead of k (saves nothing on
  the critical path — k+koff is already off it — skip unless regrouping
  anyway). Try s ∈ [−0.25, 0.25] centering via a 2-entry (m or m·√2⁻¹...)
  — this is just the 1-bit version of table+poly; the √2 factor isn't
  exactly representable, needs a df constant, probably not worth vs 3-bit
  table.
- **exp2**: center the poly on f−0.5 over [−0.5, 0.5] with the 2^0.5 factor
  folded into the exponent bit trick — 2^0.5 isn't a power of two so the
  fold needs a df constant multiply (one extra fma); halved range drops ~1
  coefficient. Compare against the 16-entry table variant; the table likely
  dominates.
- **exp**: after the direct-reduction rewrite, `expm1`, `sinh`, `cosh`,
  `tanh`, `erf`, `erfc`, `powf` all inherit the fix — re-run their refits
  afterwards, their current coefficients partially compensate upstream
  error and will be mis-tuned once exp improves.
- **tan (dedicated)**: odd rational tan(r) = r·P(r²)/Q(r²) on [−π/4, π/4]
  with the octant handled by the −1/tan reciprocal identity — two selects,
  one division; both faster and *far* more accurate than sin/cos division
  near tan's zeros/poles. Needs the mod-π/2 reduction (see Part 1).
- **asin**: on a ∈ [0.5, 1], `1 − a` is exact by Sterbenz — so the sqrt
  form's cancellation only actually bites for a < 0.5, exactly where the
  odd-poly branch (select) is cheap and accurate. The two-range select
  version is the standard shape; both sides vectorize.
- **acos near 1**: `(1.0 - a).sqrt()` is fine (Sterbenz again), but
  acos_poly's fit should be re-verified against the ulp budget near x = 1
  where acos → 0 and *relative* error explodes; may need the poly refit in
  sqrt(1−a) directly rather than a.
- **atan2**: feed the atan kernel `y/x` computed with a two_prod residual
  correction (`d = y/x; e = fma(-d, x, y)/x;` fold e into the poly's linear
  term) — 2–3 extra ops, removes the division's rounding from the error
  budget entirely. Only matters once atan's own poly is refit under budget.
- **erfc**: after swapping in a fixed exp (Part 1), the rational part's
  4/4-degree P/Q with *two* long fma chains + 1 div is latency-heavy;
  a single higher-degree numerator against a shorter denominator (or one
  poly × gaussian with a table) may rebalance. Refit territory.
- **remainder**: `(x / y).round()` → round_ties_even (semantic change —
  current ties-away is documented; offer `remainder_ieee` with ties-even
  which is both more correct per IEEE *and* faster, and keep the old one).
- **powf special exponents**: compile-time-constant y (y = 2, 3, 0.5, 1/3,
  −1) can't be dispatched at runtime branchlessly, but a const-generic
  `powf_const::<Y_BITS>()` or just doc'ing "use x*x, cbrt, sqrt" covers it;
  the df-log2 powf (Part 1) covers the rest accurately.

## Search / tuning infrastructure

- **Exhaustive-everything**: 2^32 is ~4.3e9 evaluations — minutes per
  function with a vectorized sweep on this machine. Brute-force verify
  EVERY unary f32 function against a correctly-rounded reference (rug/MPFR,
  or f64 std where provably ≤0.5 ulp of the f32 answer) and record exact
  avg ulp, max ulp, and the argmax inputs. No sampling arguments ever again.
  accuracy.rs is already most of the way there per the memory notes.
- **Worst-case corpus files**: check in the top-N worst inputs per function
  (from the exhaustive sweep) as a fast regression corpus — refit
  candidates get scored against the corpus in microseconds before paying
  for a full sweep. Re-mine the corpus after each accepted change.
- **run_descent v2**: objective = (avg_ulp, max_ulp) lexicographic or a
  penalty blend, evaluated on the corpus + periodic full sweep; moves =
  ±1..±8 ulp on each coefficient (the current random-delta scheme wastes
  most proposals); simulated-annealing temperature on acceptance. Also
  search the *seed constants* (cbrt's 0x2a509a07, the magic offsets in the
  _approx functions) jointly with coefficients — they're just more
  coefficients.
- **ulp histograms + error-vs-x plots per function** (plotters already in
  dev-deps): avg/max hide structure — a histogram instantly shows whether
  error is poly-fit-limited (uniform) or reduction-limited (spikes at
  specific magnitudes), which tells you *which* idea in this file to reach
  for.
- **Pareto table in the readme, auto-generated**: for each function ×
  variant: cyc/elem (mca), latency, avg ulp, max ulp, valid domain. Turns
  every future tradeoff discussion into a table lookup and keeps the docs
  honest.
- **Differential fuzz vs sleef/musl/core-math**: not for correctness
  (references differ) but for *opportunity mining* — inputs where sleef
  beats this crate by >1 ulp are exactly the inputs to study.
- **Known-optimal magic constants**: published exhaustive-search results
  exist for rsqrt (0x5f375a86 lineage) and cbrt seeds; cross-check the
  crate's seeds against the literature (Moroz et al. have optimal cbrt
  seed+NR constants) — someone may have already done the exhaustive search
  run_descent approximates.

## Codegen & build hygiene

- **CI asm grep** (restating from Part 1 because it gates everything else):
  compile a `#[no_mangle]` loop per public function, assert `vfmadd.*[yz]mm`
  present and no `call`/scalar `ss`-suffixed math in the hot loop. The
  scale-param duplication incident and the saturating-cast lesson (memory)
  were both silent codegen regressions this would have caught.
- **Baseline target audit's `compile_error!` half — done, tested, kept
  (2026-07-07).** `.cargo/config.toml` sets `target-cpu=native`, but that
  setting is silently *overridden* (not merged) by an environment
  `RUSTFLAGS` variable -- a well-known Cargo gotcha some CI/build setups
  hit unknowingly, with no build error and no runtime symptom beyond
  quietly-wrong accuracy/perf numbers (`f32::mul_add` falls back to a
  ~2x-slower, differently-rounded software path without FMA). Added a
  `#[cfg(not(target_feature = "fma"))] compile_error!(...)` guard at the
  top of `src/lib.rs`, with a clear message naming the cause and the fix.
  `target_feature = "fma"` is a standard cfg the compiler sets whenever
  FMA is actually enabled, regardless of *how* (config.toml, RUSTFLAGS,
  `--target`, `-C target-feature=`), so the check is robust to all of
  them, not just the RUSTFLAGS-override case that motivated it. Verified
  three ways: (1) normal `cargo build`/`cargo check` (respecting
  `.cargo/config.toml`) compile silently, unaffected; (2) an explicit
  `RUSTFLAGS="" cargo check --target x86_64-unknown-linux-gnu` (simulating
  the exact override scenario) now fails loudly and immediately with the
  guard's message, instead of silently building broken code; (3) hit a
  real, unexpected snag along the way -- the guard broke `cargo test`'s
  doctest step specifically, since `rustdoc --test`'s own compilation of
  the library doesn't inherit `.cargo/config.toml`'s rustflags the same
  way `cargo build`/`check` do (a separate, pre-existing Cargo/rustdoc
  gap, unrelated to whether FMA is "really" needed there — there are no
  actual doctests in this crate, `cargo test`'s doc-test step still has
  to compile the library to confirm that). Fixed by adding `not(doctest)`
  to the cfg condition (`cfg(doctest)` is true specifically during
  rustdoc's doctest-compilation pass) -- re-verified all three scenarios
  still behave correctly after the fix. The "decide and document the
  compile target" half of this idea (x86-64-v3 vs whatever) is a
  separate, softer question not addressed here -- this entry only closes
  the "make the invalid state unrepresentable" half.
- **Register-pressure check at zmm width — investigated, likely not
  worth forcing, reasoning documented rather than fully tested
  (2026-07-07).** First confirmed the premise: this session's dev
  machine (i5-1145G7, Tiger Lake) genuinely has full AVX-512 available
  (`avx512f`/`bw`/`dq`/`vl`/etc. all present under `-C target-cpu=native`
  in `rustc --print cfg`), yet the crate's actual compiled throughput
  loops (checked `log2_throughput` directly in `--emit=asm`) use `%ymm`
  (256-bit AVX2) registers, not `%zmm` (512-bit) -- confirmed via the
  same evidence this crate's own mca table header already documents
  ("16-wide (two AVX2 vectors)"). So LLVM, even with AVX-512 fully
  available, is *choosing* not to use it here. Looked for a way to force
  the comparison this entry asks for (zmm-width codegen, to check for
  spills) and found `rustc -C help`/`--print target-features` expose a
  real `prefer-256-bit` / `prefer-512-bit` LLVM tuning knob -- but
  couldn't find a clean, direct rustc-level way to flip it independent of
  `-C target-cpu`'s own built-in tuning table (it's baked into LLVM's
  per-CPU `X86.td` tuning entries, not exposed as a standalone
  `target-feature=` string the way ISA features like `+fma` are).
  Didn't force it via a deeper LLVM-internals route: Tiger Lake (and
  client Intel parts generally) are widely documented to suffer real
  frequency downclocking under sustained AVX-512-width execution, and
  LLVM's own per-CPU tuning tables are specifically written to account
  for exactly this -- the ymm choice here is very likely LLVM's own
  informed judgment for *this* CPU, not an oversight or missing feature,
  and overriding it against that judgment would plausibly trade "more
  work per instruction" for "lower sustained clock," a net loss on this
  particular part even before considering the register-pressure/spill
  question this entry actually asked about. Given that prior, the
  register-pressure check itself was never run (there was nothing to
  check without first getting real zmm-width codegen to inspect) --
  this entry is closed with documented reasoning rather than empirical
  measurement, unlike this session's other closed entries. Re-open if
  `-C target-cpu` is ever changed to a part where AVX-512 tuning
  genuinely favors zmm width (server Ice Lake/Sapphire Rapids-class
  parts are less downclock-sensitive than client Tiger Lake), or if a
  more direct way to override the tuning table is found.
- **`clamp` codegen — verified, confirmed correct, no change needed
  (2026-07-07).** Checked `sin_checked`'s `.clamp(-POLY_SAFE_BOUND,
  POLY_SAFE_BOUND)` directly in `--emit=asm` output (fresh, not from
  memory of an earlier related investigation this session): the scalar
  latency region lowers to exactly `vmaxss`+`vminss` (2 instructions, no
  branch), and the vectorized throughput region lowers to exactly
  `vmaxps`+`vminps` (2 instructions, packed `ymm` width, no branch/call)
  — both regions confirm LLVM already knows `f32::clamp`'s NaN-
  propagation semantics match hardware min/max's own operand-order
  behavior, exactly as this entry hoped. No codegen risk here, nothing to
  fix; this was purely a verification, and it passed.
- **Integer division in cbrt's seed**: `ax / 3` becomes a mulhi sequence —
  fine on AVX2+, but the `(bits >> 16) * 0x5556` trick in cbrt_fast is
  cheaper if the extra seed error is absorbable by the correction poly;
  measure both codegens rather than assuming.
- **uarch sensitivity pass — done (partial), reassuring result, no code
  change needed (2026-07-07).** Re-ran `llvm-mca` on the *same* compiled
  `.s` file (native codegen already only uses AVX2-width instructions,
  see the zmm entry above, so it's a valid instruction stream to
  re-schedule under a different model) with `-mcpu=skylake` and
  `-mcpu=znver3` (`-mcpu=znver4` and `-mcpu=icelake-server` both timed
  out past 2 minutes on the full multi-megabyte JSON dump this crate's
  every-function `.s` file produces -- not attempted further, the two
  that completed already give a real cross-vendor, cross-generation
  signal). Checked every divider-leaning function this crate's design
  relies on (`cbrt_accurate`, `cbrt`, `sinh`/`cosh`_throughput and their
  plain-tier equivalents, `hypot`, `remainder`, `asinh`, `acosh`)
  cyc/elem throughput on both alternate models against this session's
  native (Tiger Lake) baseline. Result: **every one of them matches or
  *beats* native on both Skylake and Zen3** -- several are near-exact
  matches on Skylake specifically (Tiger Lake's core is a Skylake-family
  descendant, so this is expected), and Zen3 in particular shows
  meaningfully *better* numbers across the board (e.g. `hypot` 0.766 ->
  0.641 cyc/elem, `cosh` 2.083 -> 1.461, `remainder` 0.647 -> 0.415) --
  consistent with AMD's Zen3 having a well-known fast, well-pipelined FP
  divider, the same underlying property this crate's "divider is free"
  findings depend on. No evidence of a "Zen-shaped hole": the crate's
  divider-leaning optimizations aren't overfit to this one CPU, they
  generalize at least as well elsewhere among the models checked. Not a
  fully exhaustive uarch sweep (only 2 of the 4 suggested targets
  completed in reasonable time), but enough to retire the specific worry
  this entry raised. No code changes.

## Wilder / probably-not-but-fun

- **Correctly-rounded mode**: with exhaustive verification (above), a
  `_cr` tier that is *proven* correctly rounded for all inputs (the
  core-math project does this for f64; f32 is exhaustively checkable) —
  double-float everything, ~2–3× cost, still vectorizable. Niche, but
  "proven ≤0.5 ulp for every input" is a hell of a readme line.
- **Coefficient sharing across functions**: exp2's Q and exp's refit, or
  sin/sinpi, could share register-resident constants in fused loops that
  call several functions; probably invisible in practice, but free to keep
  in mind when refitting.
- **Polynomials in the *bit pattern***: the `_approx` functions show
  float-bits linearity is a decent log approximation; a hybrid "poly in
  bits, correction in floats" for log2 could skip the mantissa extraction
  entirely. Almost certainly loses to the current shape, but nobody's
  measured it.
- **Vector gather LUT variants**: 64–256-entry tables via vpgatherdps for
  log/exp — gathers are slow (~1 elem/cycle) but if the poly shrinks to
  degree 2 the tradeoff at zmm width isn't obviously bad; mca can't model
  gather memory behavior well, needs real benchmarks (quickbench.rs).
- **Newton-on-output for powf**: after df-log2 powf, one Newton step in
  log-space (`r *= 1 + (y*log2(x) - log2(r))*ln2`) reuses log_2 to polish —
  expensive (second log) but a curiosity for a `powf_accurate`.
- **Batch-transcendental fusion**: `exp(log(x)*y)` where caller has arrays —
  a fused powf_slice never materializes the intermediate and keeps the df
  hi/lo in registers across the two kernels. This is where the slice API
  (above) quietly becomes an optimization surface, not just ergonomics.

---

# Part 3: specific inefficiencies spotted in the current code

Everything here is a concrete claim about existing lines, with the reasoning
spelled out. Each one still needs the exhaustive sweep + mca before adoption —
the claims are analytical and this codebase has repeatedly shown that analysis
and measurement disagree (see the Fast2Sum comments).

## sin_checked / cos_checked: round_x_over_pi

- **round_x_over_pi accuracy option**: now that the dead `quick_two_sum` error
  term `e1` is gone from `round_x_over_pi` (it was a no-op, verified
  exhaustively and applied — see jodiemath-workflow memory, 2026-07-06),
  actually *using* e1 by keeping it out of s's ulp shadow is still an open,
  untested idea: `lo = fma(x, RPI_TINY, e1) + pre_offset` merged into `rem`
  separately, or `rem = ((p0 - qh) + s) + lo_lo` — e1 is ~x·2^-49, the same
  tier as the x·RPI_TINY term, so it plausibly matters exactly as much as
  RPI_TINY does. Would cost real ops for an as-yet-unmeasured accuracy gain.
- **sin_checked pays a real add for `pre_offset = 0.0` — tried, measured,
  reverted (2026-07-07).** LLVM cannot fold `lo + 0.0` away without the nsz
  fast-math flag (−0.0 + 0.0 = +0.0 changes the value), and Rust emits
  strict FP, so sin's path carried a dead `vxorps`+`vaddps` pair on the ql
  critical chain. Fix tried: `round_x_over_pi<const HAS_OFFSET: bool>`,
  monomorphized so sin_checked's instantiation (`HAS_OFFSET = false`) skips
  the add entirely while cos_checked's (`HAS_OFFSET = true`) is emitted
  identically to before. Confirmed via `--emit=asm` diff that the two
  instructions were genuinely gone from sin_checked's throughput region (not
  just reasoned about) and that cos_checked's assembly was untouched. Result
  (mca, reproducible across repeat runs): sin_checked latency unchanged
  (109.00 cyc — the add wasn't actually gating the critical path length) and
  throughput *worse* (5.289→5.410 cyc/elem, +2.3%) despite fewer
  instructions/uOps (16000/17200→15700/16800) and a theoretically-improved
  `Block RThroughput` (65.0→64.0); `--bottleneck-analysis` showed both
  Resource Pressure (36.39%→34.46%) and Register Dependencies
  (74.20%→72.59%) percentages *dropped* too, yet simulated `Total Cycles`
  rose (8462→8656) — removing the op let the register allocator/scheduler
  make different choices for the surrounding code that cost more than the
  removed op saved, same non-monotonic-scheduling shape as the e3-downgrade
  and err-chain-rebalance entries below, just with the direction flipped
  (this time on sin_checked, with cos_checked as an unaffected, bit-identical
  control — confirmed 113.00 cyc/4.598 cyc/elem both before and after,
  ruling out mca noise). Accuracy unaffected (quick/fuzz sweep, ~100M
  samples: sin_checked |x|≤1e6 avg/max ulp 0.0355/2, matching the established
  baseline within fuzz noise) — the claim that the add is mathematically
  inert was correct, it just isn't free once the scheduler is in the loop.
  Reverted (`git checkout -- src/lib.rs`); not adopted.
- **Both `.round()` calls (`p0.round()`, `rem.round()`) are ties-away —
  partially done, tested, kept (2026-07-06).** "Any consistent rounding
  works here" turned out true for `ql = rem.round_ties_even()` but *false*
  for `qh = p0.round_ties_even()`: switching qh alone (independent of ql's
  rounding mode) regressed cos_checked's max ulp from 2 to 6 inside its
  documented |x|≤1e6 accurate range (worst x ≈ 252.9), while leaving
  sin_checked completely untouched everywhere except the already
  off-contract [1e15,∞) tail — isolated via three exhaustive accuracy.rs
  sweeps (both changed / qh-only / ql-only). Mechanism: qh is computed from
  `p0` alone, identically for both callers regardless of `pre_offset`, but
  cos's `pre_offset = -0.5` (folded into `lo`, not into qh's input) means
  qh's rare exact-half-integer ties interact badly with that offset in a way
  sin's `pre_offset = 0` never does — a repeat of the "analysis says any
  rounding is fine, measurement disagrees" lesson from Fast2Sum. Kept ql's
  change only (qh stays `f32::round`): zero regression anywhere in-domain
  for either function, real win — mca: sin_checked 114.00→109.00 cyc latency
  (5.793→5.544 cyc/elem throughput), cos_checked 118.00→113.00 cyc latency
  (5.093→4.919 cyc/elem throughput), all ~3.5-4.4% faster.
- **The −0.5 fold's ulp ceiling is a happy coincidence worth a comment, not
  a fix**: `lo ~ x·2^-24.8`, so ulp(lo) exceeds 0.5 once |x| ≳ 2^47 — beyond
  that cos's phase offset partially rounds away inside `lo` too, the same
  bug class as the documented p0 fold bug. But at x ~ 2^47 the dropped
  ~x·2^-48 tiers already contribute O(1) residual error, so the offset loss
  is hidden under the general degradation. Worth asserting in edgecheck so
  a future retuning of the tiers doesn't silently expose it.

## sin_checked / cos_checked: reduce_pi

- **`err0` from `two_sum(x, -p1)` is (almost) provably always zero —
  Sterbenz. Done, tested, kept.** p1 = qh·PI_HI with qh = round(p0) ≈ x/π, so
  whenever qh ≠ 0, p1 ≈ x·(1 ± 2^-24-ish) — comfortably within Sterbenz's
  [x/2, 2x] window — making `x - p1` *exact* and err0 ≡ 0. When qh = 0, p1 =
  0 and the subtraction is trivially exact too. The only sliver of doubt —
  |x| right at the qh = 0/±1 boundary (p0 ≈ 0.5, where ties-away rounding
  can pick qh = 1 while x < PI_HI/2·(1+ε), landing a hair outside the
  Sterbenz window) — was checked and never bites: replaced the full two_sum
  (6 ops) with a plain subtract (1 op) and deleted err0 from the err chain,
  verified bit-for-bit identical to the two_sum version across three
  separate full exhaustive accuracy.rs sweeps (every avg/max ulp and
  worst-x, every bucket, including the off-contract tail). mca: sin_checked
  5.544→5.289 cyc/elem throughput (-4.6%), cos_checked 4.919→4.598 cyc/elem
  (-6.5%), latency unchanged for both (109.00/113.00 cyc).
- **The `err` sum is a serial 4-add dependency chain — rebalancing it
  regressed latency. Tried, measured, reverted.** Rust/LLVM will not
  reassociate strict FP, so (pre-err0-deletion) `err0 + e1b + e2b + e3b -
  e3t` evaluates left to right: 4 sequential adds before the final `s3 +
  err`. The obvious fix, combined with the err0 deletion above:
  `(e1b + e2b) + (e3b - e3t)`, depth 2 instead of 4, two of whose leaves are
  ready early — analytically zero accuracy risk (verified: also bit-exact)
  and "should" only ever help latency. Measured worse instead: identical
  throughput to the flat 3-op chain, but +3 cyc latency on *both* functions
  (109→112 sin_checked, 113→116 cos_checked) — freeing `e1b + e2b` from the
  critical path apparently let LLVM's scheduler make a different choice
  elsewhere in the 64-deep latency-chain benchmark that cost more than the
  shallower dependency saved. Same "analysis says X, measurement disagrees"
  lesson as the Fast2Sum/full-tree attempts logged in the readme's
  benchmark section — reverted to the flat left-to-right chain (now 3 ops,
  err0 gone), which measured best of the three variants tried.
- **The e2 and e3 error terms sit at (or below) the already-dropped noise
  floor — try downgrading their two_prods to plain muls.** |e2| ≤ ulp(p2)/2
  ~ x·2^-49 and |e3| ~ x·2^-48·π at large x: the same tier as the terms the
  comments already justify summing plainly.
  **e2 downgrade: tried, measured, reverted (2026-07-07).** Unlike e3 (see
  below), the "sits at the noise floor" premise did *not* check out: e2 =
  the error term of `two_prod(qh, PI_LO)`, and unlike ql (which e3's
  downgrade could prove bounded to {-1,0,1} in-domain), qh is unbounded —
  there's no analogous "provably zero" argument here, only "small on
  average," and the exhaustive sweep confirmed that distinction matters.
  Downgrading `let (p2, e2) = two_prod(qh, PI_LO)` to a plain
  `let p2 = qh * PI_LO` (dropping e2 from `tier2` entirely) regressed
  sin_checked's max ulp *inside the documented `|x| <= 1e6` range*: 2 → 8
  at `|x|<=10`, 2 → 144 at `|x|<=1000`, 2 → 51,054 at `|x|<=1e6` — not an
  off-contract-tail-only effect like several other entries in this file,
  a real in-domain budget violation. Worst-x values landing near small
  multiples of π (9.42 ≈ 3π, 505.8 ≈ 161π) match this session's other
  near-zero-sensitivity finding (the 3-deep PI_A..D attempt below): a
  fixed-tier absolute error becomes a large *relative* error exactly where
  the reduced residual should be smallest. Reverted before even reaching
  the off-contract-tail buckets or an mca measurement, since the in-domain
  regression alone is disqualifying. e3's downgrade below remains the only
  one of this pair that was actually free — the difference was `ql`'s
  provable in-domain boundedness, which `qh` never had.
  **e3 downgrade: tried, measured, reverted (2026-07-07).** The "e3 is a
  guaranteed zero" premise checked out exactly as stated — an exhaustive
  scalar sweep over every f32 bit pattern with |x| < 2^25 (past both the
  fast path's 2^22·π cliff and this function's 1e6 documented bound) found
  ql ∈ {-1, 0, 1} for sin / {-1, 0} for cos always, and `fma(ql, PI_HI,
  -(ql*PI_HI))` was exactly 0.0 every single time — n·PI_HI for n in that
  set is an exact power-of-two rescale (or zero), unlike e.g. 3·PI_HI.
  Dropping the two_prod's fma (`p3 = ql * PI_HI`, `tier2 = e2 + c45`
  instead of `(e2+e3)+c45`) gave the expected shared latency win (mca:
  sin_checked 109→105 cyc, cos_checked 113→109 cyc, both -4 cyc from one
  fewer fma on the ql-dependent chain) but failed on two independent axes
  once measured instead of just reasoned about:
  1. **llvm-mca's throughput model diverged between the two callers**
     despite identical source-level reasoning applying to both: sin_checked
     throughput improved (5.289→4.976 cyc/elem, -5.9%, as expected from
     fewer total ops) but cos_checked's got *worse* (4.598→5.166 cyc/elem,
     +12.3%) — confirmed reproducible (deterministic model, re-ran twice,
     identical). `--resource-pressure --bottleneck-analysis` on the
     extracted assembly region showed why: cos_checked's simulated
     "Data Dependencies: Register Dependencies" bottleneck jumped from
     56.96% to 72.96% of cycles, while sin_checked's stayed flat
     (74.20%→73.84%) — removing the op tightened the dependency chain in a
     way that (for cos_checked's specific surrounding inlined code only)
     let the out-of-order scheduler's register pressure become the new
     bottleneck, costing more than the removed fma saved. Block RThroughput
     (the theoretical port-pressure-only estimate) actually *improved* for
     both (66.0→64.0 cos, 65.0→63.0 sin) — only the full simulated
     `TotalCycles` (what mca.rs's cyc/elem is actually computed from)
     caught the regression. Same lesson as the err-chain-rebalancing revert
     logged below: depth/op-count analysis and the actual scheduled
     simulation can disagree, and only the latter is trustworthy.
  2. **The magnitude-bucketed accuracy sweep found a real cliff in the
     already off-contract tail**, which the "moderate range" caveat
     anticipated in kind but not in scale: sin_checked's [1e9,1e10) bucket
     went from avg ulp 0.2520/max 48 to avg ulp 369.97/max 28,432,750; its
     [1e12,1e13) bucket went from avg ulp 0.2829/max 4021 to avg ulp
     2,352,582/max 2,029,384,034 (cos_checked's matching buckets moved
     similarly). Every prior accepted change in this reduction (err0
     deletion, the ql-rounding switch, the sign-flip-before-poly rewrite)
     was verified *bit-exact identical in the tail*, not just in the
     documented range — this crate's de facto bar is "no behavior change
     outside doc'd bounds," even though nothing formally promises tail
     behavior. This change breaks that bar badly. Mechanism: past |x| ~
     2^25, qh's own ulp exceeds 1, so ql (its correction) is no longer
     bounded to {-1,0,1} and n·PI_HI for large n is *not* exact — exactly
     the regime the "moderate range" phrase in the original idea flagged,
     just with a much bigger error multiplier than "the tail notices" let
     on. There's no branchless way to keep the in-domain fma-skip without
     this cost: a select on "is ql small" would need the fma computed in
     every lane anyway (for the lanes where it isn't small), so it can't
     be cheaper than always computing it.
  Net: in-domain (|x| ≤ 1e6, confirmed via accuracy.rs, bit-exact/no
  measurable change there) this is free exactly as reasoned, but the
  combination of a caller-dependent throughput regression and a severe
  tail cliff makes it a net loss once measured on the crate's actual bar
  (mca both directions + full magnitude-bucketed sweep, not just the
  documented-range analysis). Reverted; not adopted.
- **The fractional residual is already computed and then thrown away —
  the biggest structural option.** `rem - ql` (round_x_over_pi, line ~310)
  *is* the residual in units of π, to the accuracy of the whole forward
  chain. reduce_pi then reconstructs that same information the hard way
  (~30 ops of EFT machinery in x-space). Alternative: keep `frac = rem - ql`
  (exact: both are near each other), and compute r = frac·π as a df multiply
  + collapse (~4–6 ops). Accuracy tradeoff, precisely: the current scheme's
  error floor is set by the 3-word π split (~x·2^-70 region), the frac
  scheme's by the forward 1/π chain (~x·2^-48) — but the forward chain
  *already* limits q, so for |x| in the fast range the frac scheme's
  absolute error ~x·2^-48 stays ≪ 0.5 ulp of an O(1) result until |x| ~
  2^24-ish, and degrades gradually after, which is... exactly the contract
  sin_checked advertises. Plausibly halves sin_checked's total op count at
  a measurable-but-possibly-acceptable accuracy cost. Measure both axes.

## sin_checked / cos_checked: sign application and parity

- **Flip r, not the result — done, tested, kept (2026-07-06).** sin is odd,
  so (−1)^q·sin(r) = sin((−1)^q·r): applied the parity as a sign-bit XOR on
  r *before* sinf_poly instead of `s * (1.0 - 2.0 * par)` after it, for both
  sin_checked and cos_checked. IEEE negation and a multiply-by-exactly-±1
  both only ever flip the sign bit (never round), so this is provably
  bit-exact with the old code — confirmed via the full accuracy.rs
  exhaustive sweep (every avg/max ulp and worst-x value identical to the
  pre-change baseline, all 2^32 inputs, every domain bucket) and edgecheck's
  nan/inf/max/1e10/1e20 cases (bit-for-bit identical). Real-world result was
  much smaller than hoped, and mca's *latency* harness could not measure it
  at all: mca_target.rs's 64-deep chain runs each call's output through
  `mix()` (examples/support/mca_common.rs), which masks away the sign bit
  every iteration (`& 0x007fffff`) to keep values in a safe domain — so
  LLVM can *prove* the entire old `s * (1.0 - 2.0*par)` step is dead code
  across the whole chain (verified: zero compare/kmov instructions anywhere
  in the compiled `sin_checked_latency` region) and deletes it at compile
  time. The old "113 cyc" latency baseline already excluded this step's real
  cost. The new version moves the flip *before* the opaque poly, where LLVM
  can no longer see through 4 chained fmas to prove sign doesn't matter, so
  the mask-compute chain (vcmpneqss+kmovd+shll+vxorps) now shows up for
  real — making latency mode look ~1 cyc *worse*, an artifact of the
  benchmark, not a real regression (confirmed cbrt_normal's own sign
  reapplication is silently eliminated the exact same way — this blind spot
  is crate-wide, not sin/cos-specific, and affects every function whose
  result carries a sign via a bit trick rather than surviving through
  something LLVM can't reason past). Throughput mode (real array store, no
  `mix()`, immune to this) is the trustworthy number: mca throughput
  5.916->5.793 cyc/elem for sin_checked (~2.1% faster), 5.077->5.093 for
  cos_checked (basically a wash, +0.3%); quickbench wall-clock throughput
  agrees within thermal noise (sin_checked cluster ~1.38-1.40ns vs baseline
  ~1.40-1.49ns; cos_checked ~1.40-1.43ns vs ~1.43-1.49ns). Kept anyway: it's
  free (bit-exact, no accuracy cost), simpler code (one less named
  intermediate, two fewer arithmetic ops in source), and a small real win
  for sin_checked with no measurable downside for cos_checked. The bigger
  finding here is the `mix()` blind spot itself — untested follow-up: fix
  `mix()` to preserve the sign bit (e.g. `& 0x807fffff`) so latency numbers
  for *every* sign-bit-trick function (cbrt, cbrt_accurate, and any future
  parity/mulsign idea below) become trustworthy; this will likely change
  several other rows in the mca/readme tables too, so treat as its own pass
  with a full re-baseline, not a drive-by edit.
- **parity() is 2× (mul + floor + fma) on the FMA/round ports — move it to
  the integer domain. Tried, measured, reverted (2026-07-07).** Implemented
  exactly as described: a new `parity_bit(q: f32) -> u32` reads the LSB of
  an exact-integer float straight out of its bit pattern (shift = 150 −
  exponent_field into the 24-bit significand, bounds-checked since Rust's
  `>>` wraps an out-of-range shift count instead of yielding 0 the way some
  SIMD shift instructions do — needed an explicit `(0..=23).contains(&shift)`
  select, not just a clamp, to get zero for both |q| ≥ 2^24 *and* q = 0),
  composed with the sign-flip idea directly: `(parity_bit(qh) ^
  parity_bit(ql)) << 31` replaces the old `parity()`×2 + float-compare +
  select chain in both sin_checked and cos_checked (cos_checked additionally
  XORs `SIGN_MASK` since its exponent is k+1, not k). Correctness: verified
  bit-exact against the old `parity()` formula on every exact-integer f32 a
  real caller can produce — 20M small integers both signs, every
  power-of-two boundary out to f32::MAX, and 5M random large-magnitude
  integers (field ≥ 151, where the shift goes negative) — all via a
  temporary `#[test]`, zero mismatches; edgecheck also unaffected. Result
  (mca, reproduced via git-stash before/after to rule out drift): latency
  unchanged for both (109.00/113.00 cyc, exactly as expected — parity was
  already off the critical path, per the sign-flip entry above), but
  throughput got *worse* for both: sin_checked 5.280→5.412 cyc/elem (+2.5%),
  cos_checked 4.474→4.502 cyc/elem (+0.6%). Same lesson as the
  round_x_over_pi `pre_offset` and `reduce_pi` `e3`/err-chain entries
  elsewhere in this file — trading FP-port ops for more integer ops doesn't
  automatically win once the actual scheduled simulation (not the port-count
  story) is checked; here it cost more than it saved despite the "off the
  bottleneck ports" reasoning being directionally correct in isolation.
  Reverted (`git checkout -- src/lib.rs`); not adopted. The composed
  sign-flip half of the idea (skip the float compare, XOR the bits directly)
  might still be separable from the bit-extraction half and worth an
  isolated retry, but wasn't tried standalone this round.
- **cos_checked's `2.0 * par - 1.0` vs sin's `1.0 - 2.0 * par`**: both die
  with the above; noting only that today they're an extra fma each.

## sin_checked / cos_checked: the clamp

- **`clamp(-1000, 1000)` is max+min on the r critical path; a single `min`
  on y inside the poly — tried, measured unsafe, reverted (2026-07-07).**
  The blowup mechanism is r² overflowing through the poly's r⁹ term; r² ≥
  0 always, so the idea was that only an upper bound on y is needed:
  `let y = (x * x).min(POLY_SAFE_BOUND * POLY_SAFE_BOUND)` inside
  `sinf_poly`, dropping the outer `r.clamp(-1000,1000)` call site entirely
  (2 dependent ops, confirmed via `--emit=asm`: `vmaxss`+`vminss`) down to
  1. The idea's own text already flagged the risk and asked for the
  overflow arithmetic to be checked before trusting it -- checked
  empirically (isolated standalone copy of `sinf_poly` with only `y`
  clamped, swept a wide log-spaced range of finite `r` magnitudes up to
  1e30 for both signs) rather than trusting the analysis either way, per
  this crate's established practice: **it fails**. `x` (the raw residual)
  still enters the final `fma(p, x3, x)` *linearly*, unclamped -- and for
  `r ≈ -9.7e29` (comfortably finite, `f32::MAX` is ~3.4e38), the result
  overflows to `-inf`. This is *exactly* the "returns inf for ordinary
  finite input" bug `POLY_SAFE_BOUND` was introduced to fix in the first
  place, reintroduced by a different route -- clamping `y` alone bounds
  the r⁹-term blowup but does nothing about the r¹-term blowup, and the
  crate's no-inf-for-finite-input property for sin_checked/cos_checked is
  exhaustively verified elsewhere in this file, not something to trade
  away for ~1% latency. The idea's own honest fallback ("two-op →
  still-two-op, guaranteed safe") is just the current code, so there's no
  safe cheaper version left to try here without a structurally different
  fix (e.g. actually re-deriving a tighter bound on `1+p*B²` accounting
  for the *unclamped* linear term, which the idea's own arithmetic sketch
  never got around to). Not adopted; no code changes (tested in an
  isolated, non-shared standalone copy of the function specifically so
  nothing needed reverting).

## sin_checked / cos_checked: structural

- **sincos_checked / tan_checked**: with the offset-0 reduction, cos(x) =
  (−1)^q·cos(r) with cos an *even* poly of the same r — so one reduction
  (~40 ops today) serves both outputs for the price of one extra ~4-fma
  even poly. tan_checked = the ratio, one (idle-divider) division. Today
  `tan` runs two complete independent reductions; the checked pair would
  too if composed naively. This is the single largest win available for
  any tan/sincos consumer, and it's pure plumbing — no new numerics.

## Fast sin / cos

- **The PI_A..D chain is 4 serial fmas; it can be 3-deep — tried, measured,
  reverted (2026-07-07).** Implemented exactly as specified: `t = fma(q,
  PI_C, q * PI_D)` computed as its own fma+mul (independent of `x`, so it
  runs parallel to the `x - q*PI_A - q*PI_B` chain), then `r = fma(q,
  -PI_B, fma(q, -PI_A, x)) - t`. The idea's own error estimate ("the floor
  roughly doubles/triples") was badly wrong once measured: exhaustive
  sweep showed `sin (in-domain)` avg ulp 0.0645→1.4818 and max ulp
  220→866,390,494 — not a "stays orders under budget" doubling, catastrophic.
  Even well inside the previously-solid range, `|x|<=1000` went from max
  ulp 2 to 3618 and `|x|<=1e6` from 9 to 231,121. The worst-x values (9.42,
  502.7) land almost exactly on multiples of π — i.e. sin's *zeros*, which
  is the mechanism the estimate missed: near a zero the correctly-reduced
  residual `r` is tiny, so `t`'s new single-shot rounding error (previously
  spread across two separate exact-fma subtracts, `-PI_C` then `-PI_D`,
  each landing in `r`'s own shrinking-magnitude context) becomes a large
  *relative* error in `r` exactly where sin is most sensitive to it,
  instead of the roughly-constant *absolute* error tier the estimate
  reasoned about. Reverted before even measuring mca, since the accuracy
  regression alone rules this out — this crate's own budget has no reading
  under which sin's max ulp growing 4 million-fold is acceptable for a few
  cycles of latency. Files under the "measure before trusting an error
  estimate" lesson repeated throughout this document, just this time the
  estimate was wrong by many orders of magnitude rather than the usual
  "sign flipped from expected."
- **cos's `q = (kb - ROUND_MAGIC) + 0.5` is two serial adds** where sin has
  one. The obvious fold — subtract (ROUND_MAGIC − 0.5) as one constant —
  dies on representability: ulp(1.5·2^23) = 1, so 12582911.5 doesn't exist.
  A half-magnitude magic (1.5·2^22, where ulp = 0.5) makes the constant
  representable but quantizes q to halves, which is wrong. Remaining ideas:
  fold the +0.5 into the PI_A chain as a separate exact constant term
  (0.5·PI_A = 1.5703125 is exact), i.e. r = fma(q', -PI_A, x) - 0.5·PI_A
  with matching 0.5·PI_B... additions — but the extra non-fma subtract
  reintroduces a rounding near cos's zeros, so it likely trades latency for
  exactly the accuracy that matters most. Documenting as investigated-and-
  probably-not rather than open.

## Missed fma contractions (Rust NEVER auto-contracts x*y + z)

rustc emits strict IEEE ops; `a * b + c` in source stays mul+add — every
fma in this crate is explicit, and these spots missed the memo. Each fix
saves an op and/or a rounding:

- **hypot** — done, tested, kept (2026-07-07). `(x*x + y*y).sqrt()` →
  `fma(x, x, y*y).sqrt()`. mca: 25.00→21.00 cyc latency (-16%),
  0.766→0.763 cyc/elem throughput; accuracy improved slightly (avg ulp
  ~0.039→~0.034 across fuzz runs, bounded domain).
- **asinh** — done, tested, kept (2026-07-07). `x*x + 1.0` →
  `fma(x, x, 1.0)`. mca: 76.98→71.99 cyc latency, 2.905→2.806 cyc/elem
  throughput.
- **acosh** — done, tested, kept (2026-07-07). `x*x - 1.0` →
  `fma(x, x, -1.0)`. mca: 72.06→69.03 cyc latency, 2.901→2.806 cyc/elem
  throughput.
- **asin — tried, measured, reverted (2026-07-07).** `a*a - a` →
  `fma(a, a, -a)` is bit-for-bit the same fusion as the others above, but
  mca showed throughput getting *worse* (1.433→1.479 cyc/elem) with latency
  unchanged, reproduced and isolated by reverting just this one line while
  keeping the rest of the batch — the same "fewer ops doesn't always mean
  faster" pattern as the round_x_over_pi pre_offset and reduce_pi e3
  entries elsewhere in this file, this time inside an otherwise uniformly
  positive batch of near-identical changes. Left as `(a * a - a) / d + a`.
- **remainder** — done, tested, kept (2026-07-07), *without* the
  round_ties_even part. `x - (x/y).round() * y` → `let q = (x/y).round();
  fma(-q, y, x)`. Kept `.round()` (ties-away) unchanged — switching to
  round_ties_even is a separate semantic change (documented tradeoff, see
  Part 3) not bundled into this mechanical fma pass. mca: 37.00→33.00 cyc
  latency (-11%), 0.649→0.646 cyc/elem throughput. Accuracy: avg ulp
  roughly halves across repeated fuzz trials (~1000-2500 before vs
  ~600-1800 after, high run-to-run variance either way since the metric is
  dominated by a rare tie-breaking case) — max ulp unaffected, since that's
  the separate structural tie-breaking cliff documented in the readme, not
  the double-rounding this fix targets.
- **parity()** — done, tested, kept (2026-07-07). `q - 2.0 * f` →
  `fma(-2.0, f, q)`; bit-exact (the intermediate `2.0*floor` was already
  exact once `q` is an integer, so this is a pure op-count win). Shared by
  sin_checked/cos_checked: mca 5.289→5.280 / 4.598→4.474 cyc/elem
  throughput (cos_checked -2.7%), latency unchanged for both (moot once/if
  the integer-domain parity idea below lands, but free until then).
- **tanh** — done, tested, kept (2026-07-07). `exp(2.0 * x)` expanded to
  `exp2((2*x)*LOG2_E)`: two runtime multiplies (though only one actually
  rounds, since doubling is exact — the "two roundings" framing above was
  slightly off). Rewrote as `exp2(x * (2.0 * LOG2_E))`, folding the
  compile-time-constant `2.0 * LOG2_E` in up front: one runtime multiply.
  Confirmed bit-exact via an exhaustive (all 2^32 bit patterns) old-vs-new
  comparison (temporary test, removed after use). mca: 62.00→58.00 cyc
  latency (-6.5%), 1.359→1.330 cyc/elem throughput (-2.1%); sinh/cosh/exp
  (unaffected controls) stayed exactly at baseline. Same pattern anywhere a
  scale factor meets exp/ln's internal constant — powf's `log_2(x) * y`
  consumers are a candidate for the same audit, still untested.

## Other functions, specific spots

- **acosh's sign-losing domain bug — done, tested, kept (2026-07-07).** Not
  from this file's brainstorm list originally (it was readme.md's "known
  defects" backlog, alongside log1p/sinh/tanh above), but the same
  crate-wide fuzz sweep that found those found acosh at 127,929,448 avg
  ulp / max ulp reported as `u64::MAX` (an overflow sentinel from comparing
  a NaN reference against a wrong finite result) at worst x ≈ -1.6e19. Three
  distinct bugs, all from `x*x` losing information before anyone checks it,
  found and fixed one at a time as each fix's own exhaustive re-sweep
  exposed the next:
  1. Sign loss: squaring erases x's sign, so `sqrt(x²-1)` can't tell +x from
     -x once |x| > ~4096 (where `x²`'s ulp exceeds 1 and the "-1" vanishes)
     — `x + sqrt(x²-1)` collapses to ~0 for negative x instead of staying
     negative, so `ln(...)` silently returns finite garbage instead of the
     domain-correct NaN. Wrong for roughly the whole range x < -4096, not a
     narrow edge case. Fixed with an explicit `if x < 1.0 { NaN }` select.
  2. Premature overflow: valid x above sqrt(f32::MAX) (~1.84e19) makes
     `x*x` itself overflow to +inf even though the true answer (~ln(2x),
     at most ~89.6) stays finite. Fixed by rescaling before squaring for
     x ≥ 2048: `sqrt(x²-1) = x·sqrt(1 - 1/x²)`, where `1/x²` underflows
     gracefully to 0 instead of `x²` overflowing. Costs one extra division
     (a second rounding vs. the direct `fma(x,x,-1.0)` form), so kept
     select-gated to x ≥ 2048 — comfortably below where the direct form
     starts losing the "-1" term (~4096) but past where the extra rounding
     stops mattering.
  3. A second overflow survived fix 2: once x is within 2x of f32::MAX,
     `x + s` (~2x) itself overflows even though `ln(2x)` (~89) doesn't.
     Guarded with `ln(x) + LN_2` (same asymptote, no 2x formed) whenever
     the sum isn't finite.
  4. Right at the x=1 domain boundary (derivative singularity — any
     rounding gets amplified into many ulps of a tiny result), fixes 1-3
     alone left max ulp at 1522/983: `ln(x+s)` there is `ln(1+tiny)`,
     exactly log1p's own problem. Routed through `log1p((x-1.0)+s)`
     instead, with `x-1.0` exact by Sterbenz — cut max ulp to 4, matching
     the residual log1p/tanh/sinh's own fixes leave behind.
  Verified exhaustively at each step (edgecheck extended with
  `acosh(-4096)`/`acosh(-f32::MAX)` alongside the existing cases, which had
  been asserting the *old buggy* `+inf` result — updated to assert the
  now-correct NaN). Final: avg/max ulp 0.0631/4 (was
  127,929,448/u64::MAX). mca cost is the largest of this session's four
  correctness fixes, unsurprising given three independent bugs got fixed
  in one function: 69.03→110.47 cyc latency (+60%), 2.806→6.839 cyc/elem
  throughput (+144%). Kept for the same reasoning as log1p/tanh/sinh: a
  function silently returning `+inf` instead of `NaN`, or `+inf` instead
  of a valid ~89, across roughly half its domain, is a correctness bug.
- **atanh's small-x cancellation — done, tested, kept (2026-07-07), the
  cheapest of this session's five correctness fixes.** Same bug class as
  log1p (170,498,046 avg ulp / 855,638,016 max ulp, worst x ≈ 3e-8):
  forming `(1+x)/(1-x)` directly rounds to exactly 1.0 for tiny x, so
  `ln(...)` came out exactly 0. Fixed in one line, reusing log1p directly
  (same shape as tanh reusing expm1): `atanh(x) = 0.5*ln((1+x)/(1-x)) =
  0.5*(log1p(x) - log1p(-x))`. Domain (x in [-1,1], NaN/±inf elsewhere)
  fell out for free from log1p's own already-fixed domain handling — no
  new edge-case code needed, and edgecheck's existing
  `atanh(0)/(1)/(-1)/(2)` cases all passed unmodified. Verified
  exhaustive: avg/max ulp now 0.0322/3 — better than std's own reference
  here (0.0369/363409, std has its own conditioning issue near x=±1).
  mca cost, smallest relative jump of this session's five fixes: 71.00→
  80.03 cyc latency (+12.7%), 2.650→4.210 cyc/elem throughput (+58.9%).
  Kept for the same reasoning as the other four.
- **asinh: two bugs, not one — done, tested, kept (2026-07-07), sixth
  correctness fix this session.** The doc comment only documented a
  small-x cliff (812,462,495 avg ulp / max ulp 3,258,111,228 overall), but
  the exhaustive sweep's actual argmax was x ≈ -3.4028072e38 — a *second*,
  worse bug: for large negative x, `sqrt(x²+1) ≈ |x|`, so `x +
  sqrt(x²+1)` nearly cancels, and the result's relative precision is only
  as good as `sqrt(x²+1)`'s absolute error (easily 20%+ off at that
  magnitude) — the same structural cancellation as acosh's sign-loss bug,
  just landing on a wrong-magnitude answer instead of a wrong-domain one,
  since asinh (unlike acosh) is actually defined for negative x. Fixed
  both at once by computing on `|x|` and restoring sign with `mulsign`
  (asinh is odd: `asinh(x) = sign(x)·asinh(|x|)`) — `ax + sqrt(ax²+1)`
  never cancels since both terms are non-negative, which incidentally also
  resolves the small-x cliff without a separate branch, since the
  rationalized `sqrt(ax²+1) - 1 = ax²/(sqrt(ax²+1)+1)` term (computed with
  `ax²` as its own independent multiply, not extracted from the lossy
  `ax²+1` sum) feeds `log1p` the same way acosh's `d` does. Also carried
  over acosh's two overflow guards verbatim (rescale above `ax = 2048`,
  `ln(ax) + LN_2` fallback near f32::MAX) since the mechanism is
  identical. Verified exhaustive: avg/max ulp now 0.1734/4 (was
  812,462,495/3,258,111,228); edgecheck extended with
  `asinh(±2.34e-8)` (now the correct tiny value, not 0) and
  `asinh(-1e10) == -asinh(1e10)` / `asinh(-f32::MAX) == -asinh(f32::MAX)`
  oddness checks. mca cost is a genuine mixed result, first of its kind
  this session: latency *improved* (71.99→55.40 cyc, -23%) while
  throughput got much worse (2.806→9.625 cyc/elem, +243%, the largest
  relative jump of any fix so far) — plausibly the extra `abs`/`mulsign`
  and the magnitude-based selects cost real throughput but happen to let
  the scheduler shorten the serial chain, the same kind of
  non-monotonic-scheduling surprise logged elsewhere in this file, just
  landing in the "good news" direction on one axis for once. Kept for the
  same reasoning as the other five fixes.
- **atan's range select** — done, tested, kept (2026-07-07), and a much
  bigger win than expected. `let y = if a < 1.0 { a } else { 1.0 / a }` was
  compare+blend on an already-unconditionally-computed reciprocal; for
  a ≥ 0, `a.min(1.0 / a)` picks the same value in one vminps (NaN: a=NaN →
  1/a=NaN → min(NaN,NaN)=NaN, matching the old else-branch). The `a < 1.0`
  mask is still needed for the π/2 − y select, so this removes one blend,
  not the compare. Confirmed bit-exact via an exhaustive (all 2^32 bit
  patterns) old-vs-new comparison. Measured impact was far larger than "one
  blend fewer" suggested: mca latency 76.81→57.09 cyc (-25.7%), throughput
  1.984→1.410 cyc/elem (-28.9%) for atan; atan2 (calls atan directly)
  81.00→57.17 cyc (-29.4%) / 1.969→1.467 cyc/elem (-25.5%). tan (unaffected
  control) stayed exactly at baseline. Mechanism not fully root-caused —
  the old select apparently cost more than one blend's worth in this
  region (register pressure or scheduling interacting with the division,
  by analogy to other entries in this file), not just its raw op count.
- **atan2's `bothzero = !nonzerox && !nonzeroy` feeding
  `hpisignx = if nonzerox || bothzero {...}` — checked, no-op, not adopted
  (2026-07-07).** Pure boolean algebra: `A || (¬A∧B) ≡ A∨B`, so this
  simplifies to `nonzerox || y == 0.0`, dropping `nonzeroy`/`bothzero`
  entirely (fewer source-level ops). But a full `--emit=asm` diff of the
  whole compiled example (not just the `atan2` region) showed the two
  builds are **byte-for-byte identical** — LLVM's InstCombine already
  performs this exact simplification, so there's nothing left to gain at
  the source level. mca confirmed zero change on both axes (57.17 cyc /
  1.467 cyc/elem, exactly, both before and after). Not committed (no
  measurable speedup, per this crate's bar for adopting an IDEAS.md item)
  but also not a regression — recorded so a future pass doesn't
  re-discover and re-test the same already-free simplification.
- **erfc's final `fma(y, z, w)` with z = ±1 — tried, measured, reverted
  (2026-07-07).** Replaced `z`'s select + `fma(y, z, w)` with
  `mulsign(y, x) + w`, removing the `z` variable/select entirely on the
  theory that trading an FMA-port multiply for a `mulsign` xor (LLVM folded
  it into one `vpternlogd`) plus a plain add would free up the FMA/mul
  ports this function is "drowning in" (two parallel 4-deep Horner chains).
  `--emit=asm` confirmed the intended codegen shift happened (the `after`
  region's final combine is `vpternlogd` + `vaddps` where `before` had
  `vblendmps` + `fma`), and `erf` (an unaffected control, doesn't touch this
  code) stayed bit-identical (88.02 cyc / 2.101 cyc/elem both before and
  after). But mca (`examples/mca.rs`, the crate's trusted harness,
  reproduced twice) showed erfc's throughput getting *worse*
  (2.097→2.158 cyc/elem, latency flat at ~67.1 cyc) instead of better — the
  same "op/port reshuffling doesn't always pay off" pattern as the
  round_x_over_pi `pre_offset` and `reduce_pi` `e3`/err-chain entries
  elsewhere in this file. (A follow-up manual `llvm-mca --bottleneck-analysis`
  pass meant to explain *why* hit a build-caching artifact of its own -- two
  supposedly before/after dumps came back byte-identical despite a confirmed
  disassembly difference -- so the mechanism wasn't pinned down this time;
  the trusted `mca.rs` numbers, reproduced twice, were treated as sufficient
  grounds to revert without chasing the artifact further.) Reverted; not
  adopted.
- **cbrt's seed division `ax / 3`** — partially investigated, not yet
  attempted (2026-07-07). Confirmed the premise two ways before committing
  to a refit: (1) `--emit=asm` on cbrt_normal's throughput region shows the
  division really does cost what the idea claims — 4× `vpmuludq` + 2×
  `vpshufd` for the widening multiply's odd/even lane dance, vs. the
  `(bits>>16)*0x5556` trick's single `vpmulld`. (2) Measured the two
  seeds' actual relative error directly (2M-sample sweep over x in [1,2),
  compared each seed to `f64::cbrt`) rather than trusting cbrt_fast's
  overall crudeness as a proxy (its correction step is also weak, so it
  can't isolate the seed's own contribution): `ax/3`-based seed avg/max
  relative error 1.85%/1.85%, `(bits>>16)*0x5556` 2.85%/2.86% — genuinely
  ~54% wider, not just a measurement artifact. This means a same-degree
  refit of the existing 4 correction coefficients is *not* guaranteed to
  reach cbrt_normal's current budget (avg/max ulp 0.085/1) the way the
  original idea's phrasing implied ("refit the 4 coefficients against the
  new seed error range") — it may need a higher-degree correction, which
  could eat some or all of the seed's instruction savings. Stopped short
  of the actual refit (needs real numerical fitting, lolremez or extended
  coordinate descent, not just algebra) given the uncertain payoff;
  left as a well-scoped but bigger-than-one-pass follow-up, with the seed
  error numbers above so a future attempt doesn't have to re-derive them.
  **Separately: the existing correction poly (same ax/3 seed, unchanged)
  refit — done, tested, kept (2026-07-07), seventh use of this session's
  tuning recipe and the fourth real win.** A smaller, orthogonal question
  from the seed-swap idea above: given the *current* seed, are
  cbrt_normal's 4 correction coefficients themselves at their local
  optimum? Extended tune.rs with `cbrt_normal_c` (whole formula) against
  `f64::cbrt` over one octave `[1,2)` (the bit-trick seed's relative error
  pattern repeats across octaves, so one octave is representative).
  Applied and verified on the real exhaustive sweep: avg ulp 0.3265 →
  0.3112 (~4.7%), max ulp unchanged at 3 (not regressed — matches the
  crate's original documented baseline exactly). `cbrt_accurate`
  unaffected, as expected (it layers a Newton correction on top that
  absorbs small residual seed/correction error regardless of the exact
  coefficients). mca bit-for-bit unchanged (cbrt 35.06 cyc / 1.629
  cyc/elem, cbrt_accurate 59.06 cyc / 3.129 cyc/elem) — zero perf cost.
- **powf's tier mismatch — done, tested, kept (2026-07-07), the `powf_checked`
  half of this idea.** Composed the *checked* log_2 (pays the full
  denormal/negative/inf select chain) with the *unchecked* exp2 (garbage
  outside [-126, 128)) -- the expensive half bought correctness the cheap
  half then threw away. Checked by hand how bad "garbage" actually was,
  same technique that found the erf/erfc bugs: `powf(2.0, 500.0)` (should
  be `inf`) was `2.88e17`, `powf(2.0, 1000.0)` (should be `inf`) was
  `3.6e-12` -- plausible-looking finite garbage, not just reduced
  precision, for any `y` large enough to push `log2(x)*y` past 128.
  Unlike erf/erfc, powf's doc comment was already *honest* about this
  (matching exp/sinh/tanh's shared, accepted "inherits unchecked exp2"
  disclosure, not a false safety claim) -- so this wasn't a hidden-bug
  fix in the same sense, more completing the `powf_checked` half of the
  idea directly on `powf` itself (no separate `powf_fast` tier added;
  `powf` isn't part of an established fast/checked naming pair the way
  sin/sin_checked or exp2/exp2_checked are, so fixing it in place doesn't
  violate that convention). Fixed by routing through `exp2_checked`
  instead of `exp2`. Real, disclosed perf cost (unlike erf/erfc's fixes,
  which only touched a rarely-hit tail branch, this touches every call):
  mca 85.86→97.88 cyc latency (+14.0%), 2.784→3.788 cyc/elem throughput
  (+36.1%). In-domain accuracy (the documented `[-126,128)` range)
  essentially unchanged. edgecheck extended with
  `powf(2,1000)`/`(2,-1000)`/`(10,100)`. Kept for the same reasoning as
  every other correctness fix this session: `2.88e17` where the answer is
  `inf` is a real bug, disclosed doc comment or not.
  `powf_fast` (log_2_normal + exp2, undocumented-domain but fastest) is
  still a plausible future addition for callers who want the old
  behavior back, not attempted here.
- **exp2/exp2_checked share their entire Q(f) poly + the floor/fract
  preamble** — if both stay, a shared `#[inline(always)]` core with the
  exponent-construction strategy as the only difference keeps future
  coefficient refits from having to be applied twice (they've already
  drifted once: same constants, duplicated source).
- **expm1/erf/erfc/atan compute both select arms always (correct for
  vectorization) — but check the arms share subexpressions**: erf's Padé
  arm computes x² and the tail arm computes x.abs() then squares inside
  erf_poly's Horner — x² is computable once and passed to both. Minor,
  but these double-arm functions pay 2× and deserve a shared-subexpression
  pass each.
