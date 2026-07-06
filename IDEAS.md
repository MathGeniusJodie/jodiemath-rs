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
- **Degree-reduction probes**: with the 0.5-avg/2-max budget, retry dropping
  one coefficient from each poly (log_2 deg-9→8, exp2 Q deg-5→4, sinf_poly
  deg-9→7) with a tuned refit; earlier attempts predate the budget being
  stated this loosely. Each dropped term is one fma of depth and/or width.
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
- **log1p, real fix (branchless)**: u = 1 + x; c = x − (u − 1) (exact by
  Sterbenz for the interesting range); result = ln(u) + c/u (or `fma(c, 1/u,
  ln(u))` — another division for the idle divider). Fixes the documented
  small-x cliff with selects only.
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

## sin / cos / sin_checked / cos_checked / tan

- **`round_ties_even` in `round_x_over_pi`** (see cross-cutting) — two
  `.round()` calls on the ql critical path become single instructions.
- **Fused `sincos` / direct `tan`**: `tan(x) = sin(x)/cos(x)` today runs *two*
  full range reductions. One reduction mod π/2 (octant), then evaluate both
  the sin-poly and cos-poly of the same r and select/divide: tan gets ~2×
  faster; a public `sincos` returning both helps rotation-matrix-style users.
  The division is idle-divider food.
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
- **sinf_poly quantized refit**; also try fitting sin/π-scaled variants so the
  reduction constant folds in.

## cbrt family

- **Quantized/descent refit of the 4 correction coefficients** in cbrt_normal
  (`run_descent` exists but used random sampling; drive it with the exhaustive
  sweep instead).
- **Seed constant joint search**: the magic 0x2a509a07 and the poly were tuned
  separately; a joint pass (seed ± few ulp × coefficient descent) may buy the
  max-2-ulp headroom needed to drop c4.
- **AVX-512 `vgetexpps/vgetmantps/vscalefps`** kills cbrt's tiny/scale select
  dance and cbrt_accurate's three-way rescale entirely: reduce mantissa to
  [1,2), cbrt it, scale by e/3 via vscalef with the e mod 3 residue folded
  into a 3-entry table.
- **cbrt_accurate: Halley instead of Newton** from a cheaper seed — one Halley
  step from the bit seed + deg-1 correction might match Newton-from-deg-3 with
  less total width. The df32 machinery already exists.
- **Integer-division-free seed**: `ax / 3` is an integer div; LLVM turns
  constant division into mulhi, which vectorizes on AVX2 (`vpmulhuw` tricks)
  but check codegen — the `(x.to_bits()>>16)*0x5556` variant in cbrt_fast
  exists precisely for this; measure whether the full-width div actually
  vectorizes well or needs that form in cbrt_normal too.

## asin / acos / atan / atan2

- **asin small-x fix (branchless)**: select between the existing sqrt form and
  an odd poly x + x³·P(x²) for |x| < 0.5. Kills the documented near-zero
  precision loss; both sides vectorize, cost is the blended second poly
  (or refit a single form that doesn't cancel: asin(x) = atan2-style rewrite).
- **acos/asin shared kernel**: both reduce to sqrt(1−a)·poly; a shared
  computation with different post-transforms would halve the code and enable a
  combined refit at f32-quantized precision.
- **atan poly refit + degree probe**: the Padé form costs one division (fine);
  check whether a quantized refit reaches the 0.5-avg budget (currently
  ~1e-6 *relative* goal inherited from C — likely several ulp; these ports
  probably don't meet the stated budget at all and need refits, which is
  an accuracy project on its own).
- **atan without reciprocal-select**: `y = if a < 1 {a} else {1/a}` then a
  conditional π/2 flip — fine already; alternatively fit atan on [0, ∞) via
  t = x/(1+|x|) rational reduction, one division, no select chain.
- **atan2 accuracy**: y/x division error feeds atan directly; a df32 division
  correction term (`fma` residual of the quotient) into the poly's linear term
  is cheap insurance if the refit still misses budget.

## exp-family elementary (sinh/cosh/tanh) — currently 2× exp each

- **One exp, one reciprocal**: e = exp(x); exp(−x) = 1/e (idle divider, ~0.5
  ulp extra). sinh/cosh/tanh all drop from two exp2 evaluations to one — a
  large throughput win for the price of a division and a small accuracy hit
  that a refit absorbs.
- **tanh cancellation fix**: tanh = expm1(2x)/(expm1(2x)+2) — exact-ish for
  small x, one division, reuses the fixed expm1. Kills the documented
  small-x flaw branchlessly.
- **sinh small-x**: select to x + x³/6 (one fma) below ~0.5, same pattern as
  expm1's existing select.
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
- **erfc tail**: swap the interior `exp` for the fixed direct-reduction exp
  (or exp2_checked) to close the documented |x| ≥ ~9.35 unreliability; the
  clamp then becomes unnecessary.
- **erf/erfc joint refit** with f32-quantized coefficients against the actual
  budget (these are relative-1e-6 C ports, same caveat as atan).
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
- **Baseline target audit**: decide and document the compile target
  (x86-64-v3? +fma is already assumed in the godbolt header). Without
  `-C target-feature=+fma` the entire crate's accuracy story is *wrong*
  (mul_add falls back to a libm call — slow AND differently rounded); a
  `compile_error!`/build.rs check that fma is enabled would make the
  invalid build state unrepresentable.
- **Register-pressure check at zmm width**: AVX-512 doubles the register
  file but the 10-constant log poly + reduction constants in a fused loop
  can still spill at unroll factors LLVM likes; inspect `-C target-cpu`
  variants' loop bodies, not just the kernel in isolation.
- **`clamp` codegen**: verify `f32::clamp` lowers to vmaxps+vminps in the
  POLY_SAFE_BOUND path and not a select chain (its NaN-propagation semantics
  happen to match maxps/minps operand-order behavior — confirm LLVM knows).
- **Integer division in cbrt's seed**: `ax / 3` becomes a mulhi sequence —
  fine on AVX2+, but the `(bits >> 16) * 0x5556` trick in cbrt_fast is
  cheaper if the extra seed error is absorbable by the correction poly;
  measure both codegens rather than assuming.
- **uarch sensitivity pass**: the divider-is-free result underpins several
  ideas here and was measured on *this* CPU. Re-run the mca suite with
  `-mcpu=znver4`, `-mcpu=skylake`, `-mcpu=icelake-server` (mca supports any
  scheduling model) and tag each divider-leaning idea with where it wins.
  Cheap insurance against optimizing into a Zen-shaped hole.

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
- **sin_checked pays a real add for `pre_offset = 0.0`.** LLVM cannot fold
  `lo + 0.0` away without the nsz fast-math flag (−0.0 + 0.0 = +0.0 changes
  the value), and Rust emits strict FP. So sin's path carries a dead add on
  the ql critical chain. Fix: `const PRE_OFFSET: bool` generic /
  two monomorphized variants of round_x_over_pi, so sin's version simply
  doesn't have the add.
- **Both `.round()` calls (`p0.round()`, `rem.round()`) are ties-away** —
  multi-instruction lowering, and both sit on the serial chain toward ql
  (the documented last-ready value). Any consistent rounding works here
  (q just has to be *an* integer near x/π), so `round_ties_even` = single
  `vroundps` each. This is probably several cycles off sin_checked's
  latency for free.
- **The −0.5 fold's ulp ceiling is a happy coincidence worth a comment, not
  a fix**: `lo ~ x·2^-24.8`, so ulp(lo) exceeds 0.5 once |x| ≳ 2^47 — beyond
  that cos's phase offset partially rounds away inside `lo` too, the same
  bug class as the documented p0 fold bug. But at x ~ 2^47 the dropped
  ~x·2^-48 tiers already contribute O(1) residual error, so the offset loss
  is hidden under the general degradation. Worth asserting in edgecheck so
  a future retuning of the tiers doesn't silently expose it.

## sin_checked / cos_checked: reduce_pi

- **`err0` from `two_sum(x, -p1)` is (almost) provably always zero —
  Sterbenz.** p1 = qh·PI_HI with qh = round(p0) ≈ x/π, so whenever qh ≠ 0,
  p1 ≈ x·(1 ± 2^-24-ish) — comfortably within Sterbenz's [x/2, 2x] window —
  making `x - p1` *exact* and err0 ≡ 0. When qh = 0, p1 = 0 and the
  subtraction is trivially exact too. The only sliver of doubt: |x| right at
  the qh = 0/±1 boundary (p0 ≈ 0.5, where ties-away rounding can pick qh = 1
  while x < PI_HI/2·(1+ε), landing a hair outside the Sterbenz window). So:
  replace the full two_sum (6 ops) with a plain subtract (1 op) and delete
  err0 from the err chain, then let the exhaustive sweep judge the boundary
  sliver — or keep a quick_two_sum (3 ops) as the coward's version. This is
  5 ops off the function and one term off the err reduction.
- **The `err` sum is a serial 4-add dependency chain.** Rust/LLVM will not
  reassociate strict FP, so `err0 + e1b + e2b + e3b - e3t` evaluates left to
  right: 4 sequential adds (~12–16 cycles) before the final `s3 + err`.
  These are all noise-tier terms whose mutual rounding is irrelevant —
  manually balance the tree: `(err0 + e1b) + ((e2b + e3b) - e3t)` → depth 2.
  Combined with the err0 deletion above it becomes `(e1b + e2b) + (e3b - e3t)`,
  depth 2, two of whose leaves are ready early. Several cycles of latency,
  zero accuracy risk beyond noise-tier reshuffling.
- **The e2 and e3 error terms sit at (or below) the already-dropped noise
  floor — try downgrading their two_prods to plain muls.** |e2| ≤ ulp(p2)/2
  ~ x·2^-49 and |e3| ~ x·2^-48·π at large x: the same tier as the terms the
  comments already justify summing plainly. Better: for the entire
  moderate range (|x| ≲ 1.4e7), ql ∈ {−1, 0, 1}, and ql·PI_HI is *exact* —
  e3 is a full fma spent computing a guaranteed zero. Each downgrade saves
  one fma; the sweep's magnitude-bucketed mode says whether the large-|x|
  tail notices. (The existing comment establishes the *smallest* tier was
  free to keep; that's a statement about latency, not about these two fmas'
  throughput cost.)
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

- **Flip r, not the result.** sin is odd, so
  (−1)^q·sin(r) = sin((−1)^q·r): apply the parity as a sign-bit XOR on r
  *before* sinf_poly instead of `s * (1.0 - 2.0 * par)` after it. The
  parity mask depends only on qh/ql, which are ready long before r exits
  reduce_pi — so the xor hides completely in the reduction's shadow, and
  the fma+mul currently sitting *after* the poly's last fma (~8 cycles of
  pure tail latency) vanish. Same trick works for cos_checked (its
  (−1)^(k+1)·sin(r) is likewise odd in r). The mask itself: `let flip = if
  pq == pl { 0u32 } else { SIGN_MASK }; r = f32::from_bits(r.to_bits() ^
  flip)` — cmp + blend + xor, all cheap non-FMA-port ops. Also deletes the
  `1.0 - 2.0*par` arithmetic entirely. (The fast sin/cos already xor the
  *output* sign bit — they could equally pre-xor r to shave the final
  1-cycle xor off the tail, though there it's marginal.)
- **parity() is 2× (mul + floor + fma) on the FMA/round ports — move it to
  the integer domain.** For an integer-valued float q, the parity bit is
  readable directly from the bit pattern: shift = 150 − exponent_field,
  bit = (mantissa | 1<<23) >> shift & 1, with shift clamped (for exponent ≥
  150, i.e. |q| ≥ 2^24, ulp ≥ 2 so q is even and parity = 0 — the clamp
  gives that for free). Variable per-lane shifts exist (`vpsrlvd`, AVX2).
  ~5 integer ops per parity vs 3 FP ops, but entirely off the bottleneck
  ports, and it composes with the sign-flip idea above: parity(qh) ^
  parity(ql) as an integer bit, shifted to the sign position, xor'd into r.
  The whole sign story becomes integer-domain and latency-invisible.
- **cos_checked's `2.0 * par - 1.0` vs sin's `1.0 - 2.0 * par`**: both die
  with the above; noting only that today they're an extra fma each.

## sin_checked / cos_checked: the clamp

- **`clamp(-1000, 1000)` is max+min on the r critical path; a single `min`
  on y inside the poly does the same job.** The blowup mechanism is r²
  overflowing through the poly's r⁹ term; r² ≥ 0 always, so only an upper
  bound is needed: `let y = (x * x).min(POLY_SAFE_BOUND * POLY_SAFE_BOUND)`
  in sinf_poly (or a clamped-poly variant used only by the checked pair).
  One op instead of two. NaN still propagates: f32::min(NaN, c) returns c,
  but x3 = y·x re-poisons the result through the untouched x factor, so
  sin/cos(NaN/±inf) still come out NaN. Needs edgecheck confirmation, and
  note x (the residual) still enters linearly via the final fma — with y
  clamped the result is bounded by ~|r|·(1+p·B²)… check the actual bound
  with B² = 1e6: p ~ 0.17·1e6·… — redo the overflow arithmetic before
  trusting it; if it doesn't hold, clamp still reduces to `r.min(B)` +
  `.max(-B)` where only one of the two can instead ride on |r| via abs
  tricks. (Flagging honestly: the one-op version needs its own overflow
  proof, the two-op → still-two-op fallback is guaranteed safe.)

## sin_checked / cos_checked: structural

- **sincos_checked / tan_checked**: with the offset-0 reduction, cos(x) =
  (−1)^q·cos(r) with cos an *even* poly of the same r — so one reduction
  (~40 ops today) serves both outputs for the price of one extra ~4-fma
  even poly. tan_checked = the ratio, one (idle-divider) division. Today
  `tan` runs two complete independent reductions; the checked pair would
  too if composed naively. This is the single largest win available for
  any tan/sincos consumer, and it's pure plumbing — no new numerics.

## Fast sin / cos

- **The PI_A..D chain is 4 serial fmas; it can be 3-deep.** The two
  smallest terms are additive corrections: t = fma(q, PI_C, q * PI_D)
  computes q·PI_C + q·PI_D in one fma + one mul that run *parallel* to the
  first two chain links, then r = fma(q, -PI_B, fma(q, -PI_A, x)) - t.
  Cost: t is now rounded (the exact-product property is lost for the C/D
  tier), adding ~q·2^-24·6e-7 ≈ q·4e-14 error — comparable to the split's
  own truncation floor (~q·1e-14), so the floor roughly doubles/triples but
  stays orders under the fast path's budget in its valid range. Saves ~4
  cycles of reduction latency on the crate's *default* sin/cos.
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

- **hypot**: `(x*x + y*y).sqrt()` → `fma(x, x, y*y).sqrt()`. Same op count,
  one rounding fewer, strictly more accurate. One-line change.
- **asinh**: `x*x + 1.0` → `fma(x, x, 1.0)`.
- **acosh**: `x*x - 1.0` → `fma(x, x, -1.0)`.
- **asin**: `a*a - a` → `fma(a, a, -a)`.
- **remainder**: `x - (x/y).round() * y` → `fma(-q, y, x)` with
  q = (x/y).round_ties_even(). Single rounding on the subtract — this is
  *the* accuracy-critical op of the function (massive cancellation by
  construction), so the contraction materially extends the reliable range,
  on top of the two_prod idea from Part 1.
- **parity()**: `q - 2.0 * f` → `fma(-2.0, f, q)` (moot if the integer-
  domain parity lands, but free until then).
- **tanh**: `exp(2.0 * x)` expands to exp2((2·x)·LOG2_E) — two multiplies
  and two roundings on the argument. Write `exp2(x * (2.0 * LOG2_E))`
  (constant folds at compile time): one multiply, one rounding. Same
  pattern anywhere a scale factor meets exp/ln's internal constant —
  audit powf's `log_2(x) * y` consumers similarly.

## Other functions, specific spots

- **atan's range select**: `let y = if a < 1.0 { a } else { 1.0 / a }` is
  compare+blend; for a ≥ 0, this is exactly `a.min(1.0 / a)` — one vminps
  (NaN: a = NaN → 1/a = NaN → min NaN, correct; a = 0 → min(0, inf) = 0,
  correct). The `a < 1.0` mask is still needed for the π/2 − y select, so
  this saves one blend, not the compare — small but real, and min has
  friendlier port options.
- **erfc's final `fma(y, z, w)` with z = ±1**: a full fma spent multiplying
  by a sign. `mulsign(y, x) + w` — xor + add, or with w's select folded as
  integer ops. Frees an fma-port slot in a function that's drowning in them
  (two parallel 4-deep Horner chains).
- **cbrt's seed division `ax / 3`**: u32-by-3 lowers to a magic-multiply
  needing 32×32→64 (vpmuludq + odd/even shuffle dance, ~4–5 uops
  vectorized). The `(bits >> 16) * 0x5556` trick already used in cbrt_fast
  is a single vpmulld — coarser seed, but cbrt_normal's degree-3 correction
  has ~3% seed error headroom already; check whether the coarser seed's
  extra error is absorbable (refit the 4 coefficients against the new seed
  error range). If yes: several uops off both cbrt variants.
- **powf's tier mismatch**: it composes the *checked* log_2 (pays the full
  denormal/negative/inf select chain) with the *unchecked* exp2 (returns
  garbage outside [-126, 128)). The expensive half buys correctness that
  the cheap half then throws away. Either `powf_fast` = log_2_normal +
  exp2 (document the domain, much faster) or `powf_checked` = log_2 +
  exp2_checked (actually correct at the edges); the current middle serves
  neither caller.
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
