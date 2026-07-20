# IDEAS.md

Ideas tried and empirically rejected, kept here so they aren't retried for
the same reason. Adopted ideas are not listed — see git log / lib.rs for
what shipped.

## Tried and rejected

### Cross-cutting / methodology

- **lolremez degree-reduction probes**: `log_2` 9→8 (max ulp 3-5 vs 2 cap),
  `exp2` Q 5→4 (est. max rel error 40x worse), `sinf_poly` 9→7 (178x worse).
  All rejected.
- **exp2/log_2 coordinate-descent coefficient refit**: zero-move, already
  at the local optimum.
- **exp2's Q(f) poly, max-capped LP refit**: predicted a clean win (avg
  weighted error 0.043→0.013) but regressed `exp2`/`exp2_checked`/`exp10`/
  `exp10_checked` max ulp 1→2 across the board — poly was already at max
  ulp 1, no margin for the LP's continuous model to protect a real ulp
  boundary.
- **`tune.rs` zero-seed trap**: a coefficient seeded at `0.0` can only reach
  denormal-scale perturbations under integer bit-step search, not evidence
  of no headroom. On `atan_poly`'s degree bump, zero-seed reported max 18
  (unchanged) while a real scipy least-squares seed found max 3 on the same
  grid; an arbitrary nonzero seed (1e-3) isn't a fix either (converged to
  max 565). Always seed with a real fit.
- **Plain (non-max-capped) minimax LP on cbrt's correction poly**: max ulp
  3→2 but avg regressed 43% — minimax trades typical-case error for
  worst-case. Use a max-capped weighted-L1 variant instead.
- **`-C llvm-args=-force-vector-interleave=N` global sweep**: N=4/8
  bit-identical to default (loop already fully unrolled). N=2 gave real
  but mixed results (~60% of ~70 functions improved, several regressed,
  `pown` catastrophically +165.6%). No per-function scoping on stable
  Cargo, so a single catastrophic outlier vetoes crate-wide adoption.
- **PGO (+BOLT) probe on bench binaries**: inconclusive — this machine is
  too thermally noisy for wall-clock PGO-vs-baseline comparison (same
  binary swung 6.183→11.935 ns/op run to run), and `llvm-mca`'s
  single-region static model can't evaluate whole-program PGO effects
  (inlining/layout) at all.
- **Force zmm-width AVX-512 globally (`-C target-feature=-prefer-256-bit`)**:
  real, near-universal win (60/70 functions improved, median -19.4%) but
  `exp_checked` (+88.3%) and `pown` (+84.5%) regressed severely — Tiger
  Lake only has one 512-bit-wide FMA port, and functions with a long serial
  poly chain lose the free cross-copy ILP masking the default double-unrolled
  ymm build gets. No per-function RUSTFLAGS scoping on stable Cargo, so the
  catastrophic outliers veto crate-wide adoption despite the lopsided win
  ratio.
- **Integer fixed-point poly evaluation** (mantissa reductions, freeing FMA
  ports): screened via `llvm-mca --resource-pressure` before writing any
  fitting code — false premise on this CPU. `vpmulld` is 2 uOps/RThroughput
  1.00 (worse than fma's 1 uOp/0.50); a real widening mul-high+shift+add
  term chain costs ~2.6x more port pressure than one fma per poly term,
  and lands on the *same* ports fma already uses (Tiger Lake's vector
  int/fp domains share execution ports).
- **Karatsuba-style `Df32*Df32` multiply speedup**: premise doesn't apply —
  `Df32*Df32` has zero call sites anywhere in the crate (`powf_checked`
  only ever multiplies a `Df32` by a plain `f32`); the assumed dependency
  doesn't exist. (The dead `impl Mul for Df32` was removed as a side
  effect.)
- **Precision-tapered polys, `log_2`'s full poly as plain Horner instead of
  Estrin**: real accuracy win (avg 0.0031→0.0019, max 3→2, fewer ops) but
  mca killed it — latency +55.0%, throughput +42.0%. The fully serial
  9-deep Horner chain isn't hidden even in vectorized/throughput mode; a
  narrower Horner-only-for-the-tail variant is provably no better once the
  reconstruction algebra is worked through (same op count, longer critical
  path than the existing Estrin tail).
- **Coefficient ulp-neighborhood search via `tune.rs` basin-hopping**,
  applied broadly after two real fidelity bugs were found and fixed
  (`exp2_c`'s stale A/B split, `erfc_c`'s multiply/divide reordering,
  leading to real adopted wins elsewhere): `cbrt_normal`'s own coefficient
  search looked like a clean win on a 12-octave spot-check (max ulp 3→2,
  avg +1.75% only) but reversed on the real, full `cbrt`: the true worst
  point sits at the tiny/denormal rescale boundary (`x=1.3057394e-38`)
  none of the sampled octaves came near — max stayed 3, avg got *worse*
  (0.303 vs documented 0.281, ~7.8% regression). A representative-octave
  spot-check isn't a substitute for testing the whole wrapped function.
  `erf_poly`'s own apparent "improvement" (max 5→4) turned out to be an
  artifact of `tune.rs`'s `erf_tail_c` using a plain Horner chain where the
  shipped `erf_poly` uses Estrin (different rounding structure) — once the
  probe was fixed to match, the improvement evaporated completely (grid and
  dense-domain max both stayed 5, avg noise-level ±0.02%). `log_2`,
  `sinf_poly`, `expm1_near0`, `exp`'s own poly, `asin_poly`, `ln`, `log10`,
  `erf_near0`, `atan_poly7`, `atan_pure_poly` all separately confirmed
  already at their coordinate-descent optimum (noise-level moves <0.5%
  each) — this crate's heavily-retuned polys have little coefficient
  headroom left for local search alone.
- **Branch-threshold crossover retuning**, checked directly (not just
  refit) for `sinh` (0.5), `tanh` (0.25), `expm1` (0.5), `asin` (0.25),
  `erf` (0.28): all five thresholds already sit almost exactly where the
  two branches' own independently-measured error curves cross; neither an
  earlier nor later cutoff helps in any case.
- **ulp-weighted Chebyshev LP refit for `log_2`'s poly** (real scipy,
  properly row-normalized to avoid HiGHS ill-conditioning): idealized
  metric looked dramatic (max weighted residual 8.4x tighter) but a real
  ~533M-point dense sweep through the actual Estrin chain came back
  statistically identical to shipped (avg 0.25128 vs 0.25137, LP
  infinitesimally worse) — `log_2`'s real error is dominated by the f32
  rounding accumulated through the 5-fma Estrin chain itself, not by the
  underlying continuous poly's fit quality, so a tighter continuous fit
  can't move the real number.
- **Same LP technique applied to `sigmoid`'s poly** (combine-sensitivity
  weighted): idealized gain was already an order of magnitude smaller than
  `exp_pos_neg`'s analogous win (13.7% vs ~29%) and evaporated to noise on
  real verification (max ulp unchanged at 4, avg 0.11482→0.11465, ~0.15%).
- **erf_poly/asin_poly/erfc_rational combine-sensitivity LP weighting**
  (weight the fit by the final nonlinear combine's own local derivative
  instead of the poly's raw target error): worked cleanly for `erf_poly`
  (adopted elsewhere) but the *same* technique regressed `asin_poly` (grid
  max 5→6, avg 0.887→1.390; real dense sweep avg ulp 0.05681→0.06438,
  +13.3%, max unchanged) and `erfc_rational`'s numerator (max ulp
  109→129, +18.3%, avg only -4.6%). Root cause both times: the
  linearized sensitivity term (`-sqrt(1-a)` for asin, the Gaussian factor
  for erfc) itself vanishes exactly in the region that already dominates
  the function's real worst case (`a→1` for asin, large `xa` for erfc),
  so the LP trades real accuracy there for tightness in the easier
  interior. Check whether the sensitivity term vanishes near the known
  hard region before trusting this technique on a new target.
- **Same LP technique screened on paper for `tanh`**: `exp2int`'s own
  ~76-order-of-magnitude dynamic range breaks HiGHS outright ("Model
  error") even after row/column normalization; also independently capped
  by the round-off audit finding tanh's error splits evenly between poly
  and final combine, so even a perfect poly refit only closes half the
  gap. Not pursued.
- **Worst-case compare-select patch lists**: checked on `erfc` — not
  viable, ~4900 bad points spread continuously across the domain, not
  concentrated enough for a lookup-table-style patch. `cbrt_accurate`'s
  one bad mantissa fits the bar but is an already-accepted won't-fix.

### exp / exp2 family

- **Fused sincos / direct tan via shared reduction**: two attempts, both
  regressed catastrophically near tan's poles (avg ulp 0.33→1.05, max
  3000→32M) — additive combine forms cancel badly near r=π/2. Needs
  multi-word reduction to fix properly.
- **expm1 Padé degree bump 3→5**: real avg-ulp headroom in the near-zero
  branch (0.138→0.135) but the function's actual max ulp (6) lives in the
  *other* branch (`exp(x)-1`), untouched. Cost real throughput (+7%).
- **expm1 round-off budget audit**: at the worst point, rounding (2.946
  ulp) and truncation/coefficient error (2.579 ulp) are roughly balanced —
  no single dominant term to attack, unlike `exp`'s Cody-Waite fix.
- **sinh round-off budget audit**: rounding dominates truncation (3.055 vs
  1.670 ulp) but traces to the poly's own ~1.3e-7 relative fit error
  (already near f32's precision floor) amplified ~8x by the exponent-field
  reconstruction, not a fixable single rounding step.
- **tanh direct rational P(x²)/Q(x²) over full domain**: needs 13 free
  coefficients to converge — far more than any poly in the crate; a
  2-domain split needs 14, likely more work than the current expm1-based
  formula.
- **exp2/exp2_checked/exp10/exp10_checked/exp2m1/exp2_checked_df:
  round-based `Q(f)` reduction** (centered `f∈[-0.5,0.5]`, k=round(x),
  instead of `f∈[0,1)`, k=floor(x)): degree 5→4 falsified outright (max
  ulp 13). Degree 5 itself won on `tune.rs`'s coarse grid but failed real
  verification differently on every function: `exp2`/`exp10` (unchecked)
  produce **NaN inside the documented domain** (`round(x)` can land k=128
  for x strictly inside `[-126,128)`, e.g. x=127.6); `exp2_checked` has no
  structural saving to offer, pure cost (throughput 1.399→3.007 cyc/elem,
  more than double, plus max ulp 1→2 the coarse grid didn't predict);
  `exp2m1` real max-ulp regression (3→6/7) plus real cost (+3.3%
  throughput); `exp2_checked_df` real accuracy win (`powf_checked` avg
  0.0229→0.0118, max 123→118) but real ~3.4% throughput cost;
  `exp10_checked` looked like the one clean win (avg 0.0343→0.0107, max
  1→2, latency -29%, throughput -30%) until `edgecheck.rs` caught a real
  overflow-saturation bug: `exp10_checked(inf)` came out `3.237e38`
  (finite) instead of `inf`, since the round convention lets `f<0` pull an
  overflowing product back under `f32::MAX` at the clamp boundary. All six
  reverted.
- **exp: weave t1 into the poly like exp2_checked does** (`Q(r)=1+c0·r+...`,
  `p=fma(q,t1*r,t1); p*t2`): real accuracy win confirmed exhaustively (exp
  avg/max 0.0745/3→0.0522/2, cascading to expm1/sinh/cosh/tanh) but mca
  showed a real throughput regression on exp itself (1.327→1.393
  cyc/elem, 3x-reproducible). `exp`'s original `r+1.0` "+1" term is a free,
  fully-parallel add outside the poly's dependency chain; weaving replaces
  it with an extra multiply on the same contended fma/mul ports — a
  lateral shift onto busier ports, not a real reduction. Also confirmed on
  paper: precomputing `t1*t2` off the critical path would break the
  `exp(88.37628)` edgecheck by prematurely overflowing before the sub-1
  poly factor brings it back down.
- **exp10 third Cody-Waite reduction word**: already tight (avg
  0.0343/max ulp 2, 2.2B+ samples) — minimal headroom for a third word to
  collect, not worth the extra fma.
- **exp10_reduction: floor-based adjust** (backlog idea #28,
  `a = fr.floor(); k = kr + a; f = fr - a` instead of the compare+select+
  2 add/subs): expected bit-exact by construction (`fr.floor()` is
  exactly `-1.0`/`0.0` on each half of `fr`'s `[-0.5,0.5]` range,
  matching `-adjust`), but real measurement contradicted the derivation
  on both axes -- `exp10`'s own max ulp regressed 1→2 (not bit-exact
  after all, exact mechanism not tracked down), and `exp10_checked`'s
  mca throughput blew up 2.736→4.389 cyc/elem (+60.4%, reproducible),
  worse latency for both functions too. Reverted; the compare+select
  form's own codegen is apparently better-scheduled here than a
  seemingly-cheaper single `vroundps` would suggest.
- **exp2_checked: k1 from bit-twiddled k instead of a second magic-round /
  pure-integer exp2int construction**: `k as i32` (Rust's saturating
  float-to-int cast) does not vectorize even though `k` is runtime-bounded
  by an earlier `.clamp()` — LLVM can't statically prove that bound at
  codegen time, so it falls back to scalar `cvttss2si`. De-vectorized
  `exp2_checked_throughput` and 4 downstream callers (`erf`/`erfc`/`powf`
  throughput regions). Any `as <int>` cast in a hot path needs a
  codegen_check run immediately, regardless of apparent runtime bounds.
- **exp2 poly evaluated as direct `P(f)=2^f` (c0 pinned to exactly 1.0)
  instead of `1+f·Q(f)`**: paper-screening predicted a modest 1-multiply
  combine saving; actually running it (via `tune.rs`'s `tune_fixed0`) found
  the real cost is far worse than "modest" — losing one degree of freedom
  over the *entire* `[0,1)` domain (not a small sub-octave) reached only
  max ulp 463/avg 232 after full convergence, three orders of magnitude
  worse than shipped (max 2/avg 0.203), nowhere close to exp2's
  essentially-zero headroom.
- **exp_m1_over_x: trailing `/x` → `* (1.0/x)` issued at entry** (backlog
  idea #31, divider idle so the reciprocal was meant to overlap the
  reduction+poly chain instead of sitting fully exposed after it): real
  latency win (83.00→75.00 cyc, -9.6%, reproducible) but a real mixed
  result on every other axis — throughput regressed (1.798→1.842
  cyc/elem, +2.4%) and the extra rounding cost real accuracy too (avg
  ulp 0.0729→0.0757, max 6→7). Not a clean win on either non-latency
  axis; reverted.
- **softplus/logaddexp: two_sum the final `m + corr`** (backlog idea
  #40): rejected on inspection alone, no build/measurement needed --
  this file's own `two_sum(a, b)` computes `s = a + b` as its literal
  first line, then derives `e` as a *separate* auxiliary error term;
  `s` is bit-identical to plain `a + b` by construction, always, for any
  inputs. Since `logaddexp`/`softplus` return a single `f32` and the
  idea never uses the auxiliary `e` for anything, swapping `m + corr`
  for `two_sum(m, corr).0` cannot change a single output bit -- the same
  class of proven no-op as `koff-free unchecked-log fast path` and the
  `wrapping_sub` range-compare entries below, just provable from
  `two_sum`'s own definition instead of an asm diff. The real
  cancellation source is upstream: `corr` (`log1p_unit(e)`, itself
  `ln_normal(u,0.0) + corr_inner`) carries a few-ULP-of-*its-own-value*
  approximation error inherited from `ln_normal`/the division correction
  -- when `m` and `corr` are large, opposite-signed, and nearly cancel to
  a tiny true sum (`logaddexp`'s documented ~1e4 max-ulp outlier), that
  absolute error becomes a huge *relative* error in the tiny result
  regardless of how exactly the final addition itself rounds. Fixing
  this for real would need a genuinely more precise `corr` (a
  double-float `log1p_unit`, i.e. an `_accurate` tier, not a one-line
  `two_sum` swap) -- a much bigger undertaking than the idea as stated.

### log family

- **`log_2`/`ln`/`log10` integer koff fold**: bit-exact, but mca showed
  zero measurable change — LLVM already performs this reordering.
- **`ln_normal`/`log10_normal`: fuse trailing `+k*LN2_LO` into the fma**:
  mca latency deterministically *worse* by 1 cycle (reproducible,
  re-confirmed on a second attempt after an initial stale-baseline
  false-positive); accuracy unchanged, the term was already off the
  critical path.
- **log1p small-|x| dedicated branch**: adds a whole extra poly eval every
  call (branchless convention evaluates every branch unconditionally).
  Marginal accuracy gain (max 3 vs 4, avg 0.068 vs 0.073) but mca
  throughput **+48.1%**. Rejected on mca alone.
- **Direct minimax refits for ln/log10** (vs. rescaling log_2's
  coefficients): zero-move local optimum, no headroom.
- **ln_normal's poly LP refit**: isolated fit predicted a 75% avg
  improvement (largest of the session) but real exhaustive result was
  ~0.09% — noise. Isolated LP predictions don't reliably predict real
  magnitude, regardless of how large the prediction looks.
- **log1p output-side correction refit** (joint `ln(u)+c/u` objective,
  constrained so ln never regresses): coarse-grid (~613-step) result
  looked like a modest real win (log1p avg -1.6%, ln max improved 3→2) but
  wiring into the real 100M-sample fuzz showed a genuine regression
  instead: log1p avg 0.0966→0.1085, max 4→7. The coarse grid doesn't
  sample densely enough near the real fuzz's actual worst region.
- **log1p: `corr = c * (1.0/u)` with the division issued right after
  `u`** (idea #35, the same division→reciprocal-multiply reordering
  already shown to win latency at a real accuracy/throughput cost for
  `exp_m1_over_x`): opposite result here, worse on every axis measured,
  not just the two `exp_m1_over_x` traded off. mca: latency +1.7%
  (52.19→53.09 cyc), throughput +1.6% (2.337→2.374 cyc/elem) — even
  latency regressed this time, unlike `exp_m1_over_x`'s real -9.6% win.
  Real fuzz: avg ulp 0.0966→0.0969 (noise), max unchanged at 4. Likely
  cause: `ln(u)`'s own poly evaluation is a much longer, more complex
  chain than `exp_m1_over_x`'s reduction, already giving the scheduler
  plenty of slack to hide the division regardless of source-order
  hints — so early-issuing the reciprocal has nothing extra to overlap
  with, leaving only the reciprocal-plus-multiply's extra rounding as
  pure cost. Reverted. Same mechanism, different site, opposite verdict
  — a reminder that "shown to win elsewhere" still needs its own
  measurement, not just "same technique."
- **log_2 denormal path: fold ×2^24 rescale into the wrapping_sub magic**:
  killed on paper, not just measurement — the exponent bit-trick relies on
  every input being IEEE754 *normal* (fixed relationship between bit
  pattern and value); denormals have no implicit leading 1, and the needed
  renormalizing shift is different for every denormal (1 bit at
  2^-127, 23 bits at 2^-149) — a single constant offset can't replicate
  the real multiply's data-dependent renormalization.
- **koff-free unchecked-log fast path**: confirmed no-op via `--emit=asm`
  — LLVM already inlines `koff=0.0` through, no dead add anywhere.
- **log_2: atanh-form reduction `t=(m-1)/(m+1)`**: 1000x tighter in the
  underlying continuous math, but worse for real — accuracy slightly
  worse (log_2 already deep in f32-rounding-noise territory) and latency
  +43% (the division depends on `s` from the first step with nothing to
  overlap it against, unlike cbrt's early-starting reciprocal).
- **ln_accurate/log2_accurate tier from log2_df**: no measurable accuracy
  benefit — 20M-sample fuzz gave identical avg/max ulp to plain `log_2`
  (0.0061/3 both), and a 1.63B-sample strided sweep found only 7
  bit-differing outputs (~4e-9 of the domain). `log2_df`'s extra
  double-float precision only matters once something *downstream*
  amplifies the preserved low-order bits (e.g. `powf_checked`'s multiply
  by y) — collapsing straight back to a single f32 with no such
  amplification lands on the same correctly-rounded result almost every
  time.

### sin / cos / tan / sinpi / cospi / sind / cosd / tanpi / tand

- **round_x_over_pi: remove dead pre_offset=0.0 add**: instructions
  confirmed gone via asm, but throughput got *worse* — removing an op let
  the scheduler pick worse elsewhere.
- **round_x_over_pi: qh → round_ties_even**: regressed cos_checked's max
  ulp 2→6 (exact-half ties clash with cos's -0.5 offset). Reverted to
  `f32::round`.
- **reduce_pi: rebalance 4-deep chain to depth 2**: bit-exact but +3 cyc
  latency both functions.
- **reduce_pi: downgrade e2's two_prod to plain multiply**: regressed
  sin_checked's in-domain max ulp badly (2→51,054 at |x|≤1e6).
- **reduce_pi: downgrade e3's two_prod to plain multiply**: mca diverged
  by caller (sin_checked throughput improved, cos_checked's got worse) and
  the off-contract tail got dramatically worse.
- **parity() via integer bit-ops**: bit-exact, latency unchanged,
  throughput worse for both — FP-port ops beat integer ops once scheduled.
- **sin_checked/cos_checked clamp: move bound into poly's `y.min()`**:
  saves one op, but the raw residual x still enters unclamped,
  reintroducing the "inf for finite input" bug the clamp exists to
  prevent.
- **Fast sin/cos: shorten PI_A..D chain 4→3-deep**: error estimate was
  wrong by orders of magnitude — near sin's zeros the single-shot rounding
  becomes a huge relative error (avg ulp 0.06→1.48, max 220→866M).
- **cos's q fold into one constant**: dies on representability (folded
  constant isn't representable as f32); half-magnitude magic quantizes q
  wrong; folding elsewhere reintroduces rounding near cos's zeros. Never
  implemented.
- **sinf_poly quantized refit / LP refit**: quantized refit found no
  improvement (already near f32's precision floor). Separate LP refit (two
  variants) regressed cos_checked for real (avg ulp 0.081→0.79 and
  0.0495→0.0542) despite looking better isolated — cos_checked's reduction
  lands r near the poly's domain edge for small x while sin_checked lands
  r near center, so a domain-uniform grid doesn't account for the
  caller split.
- **sinf_poly: decoupled per-caller copy for sin_checked/cos_checked**
  (paper-screened): numpy screen showed both callers actually land on the
  *exact same* worst-case r once each caller's reduction formula is
  correct (max abs err 2.169e-08, identical). Decoupled per-caller refit
  only bought ~1.25x-3x tighter continuous-math error — below the
  established "<10-15% isolated signal isn't worth the round-trip"
  threshold, especially given sinf_poly is already documented near f32's
  precision floor.
- **exp_pos_neg: decoupled per-caller e/o copy for sinh/cosh**
  (paper-screened): sinh's `|x|<0.5` restriction only selects which
  integer k the reduction picks, it doesn't narrow the reduced residual r's
  own range (`[-ln2/2,ln2/2]` for both callers) — no domain difference to
  exploit, current coefficients give the exact same max error (1.325e-07)
  either combine direction.
- **sinpi/cospi: two_prod(π, r) + derivative correction**: sinpi's avg
  ulp was bit-for-bit unchanged — `sinf_poly`'s own fit error already
  dominates its budget, not the reduction's rounding. cospi did improve
  (avg 0.0861→0.0768, ~11%) since its extra reduction step leaves more
  rounding to recover. But both cost a real throughput regression (sinpi
  +21.8%, cospi +28.8% cyc/elem) — sinpi fails on both axes, cospi's real
  gain comes at real cost.
- **sind/cosd: same two_prod trick for d·DEG_TO_RAD_SMALL**: zero
  measurable accuracy improvement on both (sind ~unchanged, cosd
  bit-for-bit identical), real throughput cost both (+27%/+27.9%).
  **Follow-up, the literal HI/LO constant-split variant** (feeds a more
  accurate reduced argument directly into the poly, a different mechanism
  than the additive two_prod correction): real accuracy win for sind (avg
  ulp 0.1237→0.0675, ~45%, confirmed exhaustively), cosd flat (noise); but
  mca cost real on both (sind +11%, cosd +13.2% throughput) — cosd gets
  zero gain for a real cost, sind gets a real gain but with a real cost.
  Neither variant ships. (A zero-benefit result from one fix *mechanism*
  is not reliable evidence against a mechanistically-different literal
  proposal targeting the same reduction step — the HI/LO split had to be
  tested for real, separately, to close this out.)
- **sin fast tier: drop PI_D for a fitted 3.5-word split**: a bare
  reduction-residual probe looked genuinely better (max abs error 1.04e-7
  vs 1.19e-7) but wiring into the real sin/cos construction and the full
  100M-sample fuzz told a different story: sin avg ulp 0.0645→0.1399, max
  412→**1,824,546**; cos avg 0.2917→0.3320, max 2780→**386,929** — a
  near-a-zero-of-sin relative-error blowup a uniformly-sampled bare-residual
  probe can't see. An entire word's worth of pi's own precision genuinely
  can't be recovered by any single replacement word across sin's whole
  domain.
- **atan2's bothzero/hpisignx boolean simplification**: compiled output
  byte-for-byte identical — LLVM's InstCombine already does this.
- **sincos_checked, shared reduction**: bit-exact refactor verified
  exhaustively, mca predicted ~19% throughput win — real wall-clock (after
  fixing two bench-shape artifacts) found it ~1.4% *slower*. The one case
  this session where mca and real hardware disagreed in direction, not
  just magnitude.
- **tan via mod-pi/2 reduction + dedicated poly**: premise (division near
  a pole amplifies error) was wrong — spot-checks showed the old
  sin(x)/cos(x) form was already 0-1 ulp at the actual poles; the real
  ~3000 max ulp comes from sin/cos's own large-x domain cliff, inherited
  regardless. New reduction measured worse everywhere with a narrower safe
  domain.
- **mod-pi/2 reduction with paired even/odd polys**: degree only drops by
  1 coefficient (not "hard" as backlog claimed), and real op count is
  ~25% *more*, not less, once both polys are evaluated unconditionally.
  Killed by op-count reasoning alone.
- **sin_checked/cos_checked: ql via magic-round, parity from the low
  mantissa bit** (backlog idea #45): the premise ("rem is small, so the
  magic add is exact") is false — `round_x_over_pi`'s `rem` is
  `(p0-qh) + lo`, and `lo` includes `x*RPI_LO`/`x*RPI_TINY`, which grow
  linearly with `x`, not a bounded double-double correction. A direct
  probe (scalar Rust, exact f32 semantics) found `|rem|` reaching
  ~9.4e30 near `f32::MAX`, ~24 orders of magnitude past the magic
  trick's ~2^22 exact-rounding range. Implemented anyway to confirm via
  real fuzz (bit-exact per-branch reasoning suggested it might still
  work by accident at that scale): catastrophic failure, not a
  precision nit — `sin_checked`/`cos_checked` over `[1e15,1e16)` came
  back avg ulp ~7.3e8/1.0e9, max ulp 2130706432 (both), and the full
  `(all f32)` sweep was avg ~3.3e8 both functions, same max. Root cause:
  once `|rem|` swamps `ROUND_MAGIC` (1.5·2^23), `rem + ROUND_MAGIC`
  rounds back to plain `rem` (the magic constant vanishes in the sum),
  so the extracted "parity bit" is just bit 0 of `rem`'s own float
  encoding at that magnitude — unrelated to true integer parity, unlike
  the existing floor-based `parity()`, which stays correct (if
  trivially "always even") at any magnitude since it's a real mod-2
  computation, not a mantissa-alignment trick. This is a genuine wrong-
  *sign* bug, the same class already fixed once for this pair (see the
  sin_checked range-invariant entry in git log) — reverted immediately,
  `git checkout -- src/lib.rs`. Any future attempt at this idea needs
  the same large-magnitude fallback idea #46 already plans for `qh`,
  applied to `ql`/`rem` too, not a magic-round-everywhere assumption.

### asin / acos / atan / atan2

- **acos: mulsign-reassociate the trailing `+select(PI,0)`**: not
  bit-identical to the original (60,142/49.6M samples differ), and a real
  accuracy regression: max ulp jumped 6→**121**. The new
  `FRAC_PI_2 - y` subtraction suffers catastrophic cancellation near x≈0
  (where acos(x)≈pi/2) — unlike a pure reassociation of already-existing
  values, this introduces a genuinely new subtraction.
  **A safe variant** (`correction = FRAC_PI_2 - mulsign(FRAC_PI_2, xn)`,
  never touching y in a subtraction) is bit-identical to shipped and
  accuracy-neutral, but mca showed a real, small *regression* (throughput
  0.820→0.834 cyc/elem, +1.7%) — a provably-safe reassociation still needs
  an actual mca measurement, "replaces select with sign-bit arithmetic"
  doesn't automatically transfer from atan_latency's own successful case.
- **acos_poly Horner→Estrin**: real latency win, but fma reassociation
  regressed asin max ulp 9→12, acos 4→5 (retuning made it worse, →6).
  acos's accuracy is a protected invariant.
- **erfc's n/d rational Horner→Estrin**: small theoretical win, measured
  as a wash on speed plus a real accuracy cost (avg +2.7%).
- **erf's near-zero Padé branch refit / tail branch (erf_poly) refit**:
  both — max ulp unchanged, avg moved <0.3%. No headroom.
- **Same near-zero branch, LP numerator refit**: isolated fit predicted an
  84% avg improvement, but real fuzz found a *regression* (0.317→0.325),
  worst case landing right at the 0.28 branch crossover with erf_poly.
  Check the crossover neighborhood for any poly next to a domain split.
- **acos_poly unconstrained joint acos+asin objective**: improved joint
  score but regressed acos's own max ulp 4→5.
- **acos_poly degree 6→7**: zero-seeded 8th coefficient converged to
  exactly 0.0 (useless). Re-seeded with a real scipy fit found real
  headroom (acos avg/max 0.496/4→0.490/3) but asin was unmoved and mca
  cost was real (+7-11% both functions) — not worth it at this magnitude.
- **acos_poly Df32 leading-term split, pi/2 hi+lo**: dramatic acos avg
  improvement (0.496→0.068) but acos max ulp regressed (4→5) and asin got
  worse on both axes (max 9→12). Retuning recovered asin but acos max
  still regressed (4→6), plus a real mca cost.
- **acos_poly joint LP, both acos+asin max ulp capped** (two attempts):
  attempt 1 let the pi/2 constant drift onto a worse f32 value, corrupting
  near-zero calls (acos avg 0.496→1.905). Attempt 2 forced the constant
  exact but still regressed for real: asin max ulp 9→12 — the identical
  fragility signature as the Df32-split entry above, not a fluke.
- **atan2 division-residual correction**: looked like a 10x win until a
  NaN-sentinel bug was fixed — evaporated to noise-level (atan's own poly
  error dominates). Real cost (latency +14%, throughput +43%) for zero
  benefit.
- **atan: correct the 1/a fold's division rounding**: zero accuracy
  benefit — avg ulp exactly unchanged (0.0675), max ulp *regressed* 3→4
  right at the a=1 fold boundary. Severe mca cost: atan throughput
  1.491→**4.196** cyc/elem (+181%), atan2 +57%.
- **atan2: dual up-front divisions** (backlog idea #19, `z=y/x` and
  `q=x/y` both computed independently instead of atan's internal serial
  `1/(y/x)`): verified correct via edgecheck (every zero/sign/inf/nan
  combination still passes) and a genuine, real accuracy improvement
  confirmed on a 99M-sample same-seed old-vs-new comparison (avg ulp
  0.068735→0.068683, max ulp unchanged at 4 both) -- the predicted
  "one less rounding" effect is real, just tiny (atan_poly's own fit
  error still dominates the budget). But the predicted *latency* win
  never materialized: mca showed latency flat (61.28→61.22 cyc, noise)
  and a real throughput *regression* (1.662→1.729 cyc/elem, +4.0%,
  reproducible) -- a second division doubles the shared divider port's
  pressure across a vectorized throughput run even though it removes a
  dependency from any single call's critical path, the same "helps
  latency, costs throughput" tradeoff already documented for
  atan_latency's own division-free design, just the opposite direction
  here (adding a division instead of removing one). Reverted -- real
  but negligible accuracy gain isn't worth a real throughput cost.
- **acos_accurate opt-in tier** (Df32 pi/2 + two-product sqrt(1-a)*poly
  combine): measured zero improvement — extra precision doesn't survive
  collapsing back to f32 without a downstream user. (While investigating,
  found and fixed a real, unrelated 1-ulp transcription bug in acos_poly's
  shared pi/2 constant — the shared default's accuracy improved for free,
  making a separate accurate tier unnecessary.)
- **asin small branch: check 1-a rounding in [0.25,0.5)**: real but modest
  — a0.25-0.5 avg ~1.4/max 6-7, a disproportionate contributor to the
  domain-wide average (0.0251) via octave-count weighting, but not the
  binding max-ulp constraint. Not actioned given effort/payoff.
  **Chased via ulp-weighted Chebyshev LP** targeting this exact band:
  idealized metric looked dramatic (max weighted residual 2.13→0.10, ~20x)
  but real construction *regressed* on every axis (full domain avg
  0.731→1.165, max 7→8; the flagged band avg 1.426→2.046, max 7→8). Root
  cause: the poly's own idealized fit-only residual at the worst point was
  only ~0.53 ulp-equivalent — nowhere near the real ~1.4 avg/max-7 — so the
  continuous fit was never the bottleneck; the LP was optimizing a residual
  that was never the actual problem.
- **asin via atan2(x, sqrt((1-x)(1+x)))**: max ulp better (9→4) but avg
  nearly doubled (0.0506→0.0953) — a real tradeoff, not a clean win. And
  mca was catastrophic regardless: latency +69.8%, throughput
  **+274%** (nearly 4x). Routing a single-branch function through a full
  binary function costs far more than "one extra call" intuition suggests.
- **atan_poly denominator LP refit (numerator fixed)**: essentially no
  movement (<1%, coefficients matched shipped to 6+ digits). A separate
  numerator-only refit predicted 11-15x but the real result was only
  ~0.9% with max ulp unmoved at the identical worst-case x.
- **asin_small: one more Taylor term**: real accuracy improvement
  confirmed after re-verification (max ulp 9→7, avg 0.0251→0.0202, ~20%
  better — a stale prior "tied co-bottleneck" finding no longer held after
  an intervening asin/acos_poly decoupling). But not adopted: real,
  unavoidable throughput cost on both evaluation orders tried (plain Horner
  extension +7.9% throughput; Estrin regroup +12.3%, worse than Horner).
  Since `asin_small` is unconditionally evaluated (branchless), there's no
  way to pay for the extra precision only where needed.
- **erf_poly's a0≈3.4e-5, naive zero-out**: confirmed genuinely
  load-bearing, not a fitting coincidence — exhaustive sweep avg ulp
  0.3166→2.3367, max ulp 5→539.
- **asin_poly combine-sensitivity LP refit**: see cross-cutting section —
  real regression (avg ulp +13.3%) from the sensitivity term vanishing
  near asin's real hard region (a→1).
- **erfc_rational numerator combine-sensitivity LP refit** (denominator
  fixed): real max-ulp *regression*, 109→129 (+18.3%), avg only -4.6% —
  see cross-cutting section for the shared vanishing-sensitivity mechanism.
- **tanh combine-sensitivity LP**: screened, not pursued — see
  cross-cutting section (HiGHS scaling failure + capped payoff per the
  round-off audit).
- **erfc's final fma→mulsign+add**: intended codegen shift confirmed, but
  throughput got worse.
- **atanh via single log1p call**: algebraically simpler, but real
  regression — max ulp 3→31303 (tail only, near x≈-1). The two-log1p form
  passes exact literals; the one-log1p form's division rounds once, and
  log1p's derivative diverges near the singularity, amplifying that tiny
  rounding error ~5 orders of magnitude.
- **erfc: exact exponent via two_prod**: real accuracy win (avg/max
  0.311/109→0.306/93) but real mca cost (+5-7%) for shaving an
  already-far-over-budget max ulp further.
- **erfc in log space, single poly over [0,10]**: lolremez convergence
  poor even at degree 15 (16 coefficients). erf_poly can't be reused
  directly (catastrophic error, ~297000 avg ulp — only ever fit for erf's
  own objective).
- **Three-interval atan reduction**: accuracy case fully confirmed, but
  implementing it nearly doubled cost (latency +72%, throughput +99%) —
  the u-transform's division must resolve before the poly's own division
  can start for half the domain, two sequential full-latency divisions
  instead of one. A flagged division-on-critical-path risk deserves an mca
  check before the accuracy work, not after.
- **atan_latency's poly LP refit**: isolated fit predicted a modest 6%
  improvement — too weak a signal, found nothing real (avg ulp 0.0516 vs
  0.0517, statistically identical).
- **atan_poly 4/4 Pade bump, properly seeded** (backlog idea #59):
  scipy `least_squares` (plain L2, multi-start, and an IRLS-style
  max-reweighted approximation to minimax, three separate attempts)
  all converged to excellent *continuous*-math fits (max abs error
  ~5e-11 to ~1e-10 over [0,1]) but none beat the existing 3/3 fit once
  actually wired through the real fma-chain and scored on a 200M-sample
  real f32 fuzz: best attempt tied the existing max ulp (3) but still
  lost on avg (0.0688 vs 0.0678, ~1.6% worse); the other two attempts
  lost on both axes (avg ~0.069, max 4). Same "isolated fit doesn't
  predict real magnitude" lesson as several entries above, one level
  further -- here the isolated fit didn't even predict the right
  *sign* of the comparison. The existing 3/3 coefficients are
  documented as coordinate-descent-tuned on top of their own
  least-squares seed; matching or beating that with a naive one-shot
  refit at a higher degree needs the same discipline (or a proper
  weighted-LP/minimax tool), not just a bigger scipy fit. Not pursued
  further this session.
- **Plain atan: port atan_latency's own adopted mulsign-reassociation**
  (backlog idea #16's atan slice): apply `mulsign` to the poly result and
  `FRAC_PI_2` individually, then select/subtract, instead of selecting
  first and applying one final `mulsign`. Bit-exact as expected (same
  identity atan_latency already verified; edgecheck + fuzz confirmed
  atan's well-documented avg ulp 0.0675 unchanged) but real mca
  regression: atan throughput 1.491→1.529 cyc/elem (+2.5%, reproducible
  across repeat runs), latency/atan2 flat. Unlike atan_latency (no
  division, poly-only), plain atan's `1.0/a` reciprocal changes the port-
  pressure picture enough that the same reassociation lands as a net
  loss here — the win doesn't transfer between the two constructions
  despite the identical algebraic shape. Reverted.

### hyperbolics / sigmoid

- **sinh_small minimax refit** (the untuned-Taylor lever that shipped for
  `asin_small`, applied to `sinh`'s `|x|<0.5` branch): refit the three
  non-leading coefficients (`1/6`, `1/120`, `1/5040`) as an equal-degree
  minimax fit over `[0,0.5]`, `c0` pinned to 1.0. Zero perf cost by
  construction (same op count, bit-identical codegen), but the premise
  that carried `asin_small` doesn't transfer. The idealized fit signal is
  already sub-0.1 ulp (max ulp-equiv 0.087→0.009): `sinh` is *entire* and
  its Taylor series converges fast, so the truncation error over `[0,0.5]`
  is negligible and the branch is rounding-chain-dominated, unlike
  `asin_small` whose `[0,0.25]` truncation error is genuinely ~4 ulp (asin's
  sqrt singularity at |x|=1). Real exhaustive `[0,0.5)` sweep: avg ulp
  0.02952→0.02951 (noise), in-branch max ulp 2→1 — but `sinh`'s *headline*
  max (4-5) lives entirely in the `|x|>=0.5` `exp_pos_neg` reconstruction
  (worst x 4.67 / -0.828), untouched, and its avg is unchanged, so the
  2→1 is an invisible sub-region move. Rejected: swapping self-documenting
  exact Taylor coefficients for opaque magic numbers buys zero movement in
  any reported metric. Closes the "untuned Taylor branches" lead —
  `asin_small` was the only branch of that pair with real headroom.
  *Sequel (shipped, see lib.rs):* the opposite lever — *dropping* a term to
  degree-2 — did land as a **speed** win. `sinh_small`'s branch max (2) sits
  well under `sinh`'s headline max (5, in the other branch), so a shorter poly
  can trade its spare accuracy for one fewer fma: throughput -5.6%/-6.1%/-7.2%
  on sinh/sinh_throughput/sinh_checked, latency flat, headline max unchanged,
  branch max 2→3, sinh avg 0.073→0.082. (This is why the same-degree refit
  above is a no-op but the degree *drop* is a win — the branch had spare
  accuracy to spend, just not spare fit quality to gain.)
- **Compensated hypot** (`e=fma(r,-r,s); r+e/(2r)`): no measurable
  accuracy improvement (already near correctly-rounded) but real cost:
  latency +109%, throughput +87%.
- **tanh via exp_pos_neg ratio** ((ep-en)/(ep+en) + small-x Taylor branch):
  real accuracy win (avg ulp 0.024 vs 0.146, ~6x tighter, same max ulp)
  but the premise ("division is idle") missed the actual cost driver —
  this route needs `exp_pos_neg` to evaluate *two* full polynomials where
  tanh's existing expm1(2x) route evaluates one. Real regression: latency
  +2.7%, throughput +25.8%.
- **atanh via single log1p on |x| + mulsign**: sidesteps the catastrophic
  failure mode of the signed-x version (max ulp only 3→4, not 3→31303),
  but still a real net accuracy loss (avg 0.0313→0.0352, max 3→4) with a
  mixed perf result (throughput -29% better, latency +2.9% worse) — not a
  clean win on either axis.
- **Dedicated sinh/cosh kernels via reassociation**: real latency win for
  both (~11-12%) but throughput split by function (sinh better, cosh
  worse) and accuracy split the other way (cosh fine, sinh regressed ~12x
  avg ulp). When two functions share a reassociated intermediate, check
  accuracy for both independently.
- **exp_pos_neg_checked_half: fold the ×0.5 into the integer reciprocal
  trick** (backlog idea #17): `t1*0.5`/`t1n*0.5` as exact exponent-field
  decrements (`t1.to_bits() - 0x0080_0000`, `0x7E80_0000 - t1.to_bits()`)
  instead of float multiplies. Bit-exact (edgecheck's known boundary
  probes plus a full exhaustive 2^32-pattern sweep of both callers, not
  just quick fuzz, given this site's documented history of exponent-
  boundary bugs — max ulp landed at exactly 5/5, matching the existing
  documented magnitude, no blowup). But mca showed the same asymmetric
  split as the entry above, same shared-callee mechanism: sinh_checked
  throughput improved (2.527→2.346 cyc/elem, -7.2%, reproducible) while
  cosh_checked got much worse (2.212→3.203, +44.8%, reproducible) — a
  `#[inline(always)]` shared body gets independently rescheduled per
  caller once inlined, and the register-pressure/port-assignment
  consequences aren't predictable from op count alone (cosh_checked's
  throughput region actually has *fewer* total vector instructions than
  sinh_checked's, which has sinh's own extra small-x branch on top and
  still came out faster). Reverted — cosh_checked's regression is too
  large to accept for sinh_checked's smaller win.
- **Newton-free correction for rsqrt** (`e=fma(r,r*x,-1)`,
  `r_new=fma(-0.5*r,e,r)`): real accuracy win (avg ulp 0.2599→0.1226, ~2x
  tighter, max unchanged at 1) at a real modest cost (latency +43.0%,
  throughput +8.8%). Doesn't clear the bar — accuracy is already near its
  practical ceiling, this is polish, not a documented-defect fix.
  **Same technique on rhypot**: a real *regression*, not a smaller version
  of rsqrt's win — avg ulp 0.065→0.175 (2.7x worse). `rhypot`'s `h2` is
  already a singly-rounded approximation of the true `x²+y²`; refining `r`
  to more precisely satisfy `r²=1/h2` just makes it a better reciprocal of
  an already-wrong target, destroying incidental cancellation the naive
  form got for free.
- **erfcx for |x|>10: asymptotic-tail fix**: the correctness gap is real
  and unbounded, not just imprecise — `erfc_rational`'s internal clamp
  freezes at `erfc_rational(10.0)` forever past x=10, and unlike `erfc`
  itself (which has a decaying exp(-x²) factor masking it), `erfcx`'s
  construction exists specifically to cancel that factor, so the frozen
  value comes back naked: relative error ~10% at x=11, ~50% at x=15, ~99%
  at x=20, ~895% by x=100, unboundedly worse beyond. A working asymptotic-
  tail fix was verified (bit-identical ≤10, max rel error 0.00015% up to
  200) but costs a real +17.7% mca throughput (an unconditional second
  division) — rejected on cost; doc comment corrected to describe the real
  "freeze, not decay" mechanism instead.
- **erfc round-off audit, [0.03,0.25] secondary bump**: root-caused to fit
  truncation error (3.96e-7) exceeding combined rounding (2.22e-7) at the
  worst point. A coordinate-descent refit biased toward this band found a
  real local grid gain (max 80→79) but *reversed* on the real, unbiased
  full [-10,10] domain sweep: shipped max 109/avg 0.321 vs tuned max
  111/avg 0.342 — worse on both axes; the biased grid was borrowing
  accuracy from elsewhere. A follow-up unbiased-grid coordinate descent
  against the real formula found essentially zero movement and a
  bit-for-bit-identical dense verification — erfc_rational is confirmed at
  a genuine local optimum, no headroom at all.
- **erfcx round-off audit fix** (Df32-precision `x²·LOG2_E` via
  `exp2_checked_df`): real, substantial improvement (worst point 121→20
  ulp, domain avg 0.20→0.15 ~25%, max 122→20 ~84%) but real modest mca
  cost (+18.0% throughput); still doesn't clear this function family's own
  established "must land in single digits" bar (20 is still double
  digits).
- **erfc exponent via mantissa-mask hi/lo split** (Cephes trick): zero
  effect — exhaustive sweep bit-identical to baseline in every field. Only
  compensates the squaring step; the very next op (`-xa²*LOG2_E`) is still
  an uncompensated multiply that dominates wherever the true worst case
  lives.
- **erf tail via expm1-shape** (`exp2m1` instead of `1-exp2_checked`):
  zero accuracy effect — literally zero differing bit patterns over 711k
  samples. erf_poly's output at the branch's worst point already lands on
  exp2m1's same direct-branch rounding that exp2_checked effectively gets
  close to; the theorized cancellation doesn't occur at the real worst
  point. Also a mixed perf regression: latency -18% but throughput +13.6%
  (paying for exp2m1's unused Padé branch on every call).
- **erfc domain split [0,2]/[2,10]**: underlying math 100-380x tighter
  per-domain; real implementation improved avg ulp ~39% (0.311→0.189) but
  max ulp barely moved (106→105) — the bottleneck lives downstream of the
  correction term entirely. Bar was 109→single digits.
- **Select-tree "LUT" for exp2**: backlog's 2-level/degree-3 combination
  doesn't reach competitive accuracy (24x worse per scipy); even a
  corrected 3-level/degree-8 version regressed real ulp (max 2→3, avg
  0.20→0.43).
- **erfc negative-side accuracy survey**: confirmed real asymmetry (x≥0
  avg 0.4815/max 109 vs x<0 avg 0.1397/max 6) is structural, not a
  fixable blended-refit artifact — ulp is a relative measure and the two
  sides sit at very different output magnitudes (near 0 vs near 2) for the
  same absolute error in the shared rational. A sign-split refit wouldn't
  help.
- **Centered-variable refit for erfc's xa**: confirmed the exp2-centering
  finding transfers — centering measured ~73x worse.
- **Compensated-Horner accuracy tier for erfc**: double-float evaluation
  improved avg ulp modestly (~14%) but left max ulp flat. Root cause:
  erfc's own exponent expression has up to 87 ulp of error in plain f32
  before exp2_checked even runs — the same cause the already-rejected
  two_prod fix targets, not worth its mca cost.
  **Re-tested inside `erfc_accurate`** (idea #64's own sequencing note,
  now that the exponent fix is real, not hypothetical): compensated
  Horner for `erfc_rational`'s numerator/denominator (Df32 multiply +
  `quick_add` at each step, `div_to_f32` at the end) still only moved
  the needle modestly even with the exponent no longer dominating —
  avg ulp 0.245→0.227 (~7.6%), max ulp 12→10 (barely). Real mca cost
  was severe this time: erfc_accurate's own throughput more than
  doubled (3.033→6.189 cyc/elem), latency +17.7% (74.00→87.11 cyc) on
  top of the exponent fix's own already-accepted cost. Rejected —
  `erfc_rational`'s degree-4 fit itself (not its evaluation rounding)
  is the real remaining bottleneck; a genuine degree bump (the rest of
  idea #64) is the more promising untried lever, not compensated
  evaluation of the existing fit.
  **The degree bump itself, tried next (idea #64's remaining half)**:
  degree 5/5 rational (one extra fma per chain), seeded from a real
  scipy/HiGHS minimax LP fit against `erfcx(xa)` (the rational's own
  target, `erfc(xa)=exp(-xa^2)*erfc_rational(xa)`; residual is linear in
  the 10 coefficients since there's no p*q cross term, so a true
  Chebyshev LP applies directly, no L2-vs-minimax tradeoff needed) and
  polished by `tune.rs`'s coordinate descent (new `erfc_c5` target, kept
  in tune.rs). Real, consistent avg-ulp win across all five
  `erfc_rational` consumers on the full exhaustive `|xa|<=10` sweep
  (erfc 0.3055→0.2486, erfc_accurate 0.2452→0.1948, erfcx 0.3768→0.3280,
  erfcx_accurate/erfcx_checked 0.3263→0.2775/0.2779, all ~13-20% better)
  but failed on both of the bar's other axes: max ulp stayed flat for
  four of the five and *regressed* on `erfc_accurate` specifically
  (13→14) — the exact function idea #64 was sequenced to target — and
  mca showed a real, broad cost: throughput +3.9% to +6.5% on all five
  (two extra fma, as expected), plus latency costs that didn't scale
  with op count the way throughput did (erfcx/erfcx_accurate +10%,
  `erfcx_checked` +46.9% (41.00→60.22 cyc) with erfc/erfc_accurate's
  own latency flat) -- the extra denominator term (where the LP put
  nearly all the freedom; the numerator's 5th-degree coefficient
  converged to ~1e-8, effectively unused) lands on each caller's
  critical path differently depending on what else that caller's own
  branch does around the shared `erfc_rational` call, the same
  per-caller-rescheduling class of surprise already documented for
  `exp_pos_neg_checked_half`'s fold. Reverted (`git checkout --
  src/lib.rs`); a real avg-only win with a real max-ulp regression on
  the targeted function and broad mca cost doesn't clear this session's
  bar. `examples/tune.rs`'s `erfc_c5` scaffold is kept for any future
  attempt (e.g. a denominator-only degree bump, or dropping the
  near-zero numerator term to claw back one of the two fmas).
- **erf joint boundary+coefficients refit**: cheap proxy sweep (8
  threshold candidates against unrefit coefficients) found a flat plateau
  around the current 0.28 boundary — no headroom, matching separate
  fixed-boundary refits. Full joint optimizer never built.
- **asinh/acosh single-sqrt restructure** (backlog idea #54): select the
  sqrt *argument* before the call instead of computing both branches'
  own unconditional `.sqrt()` and selecting the *result* -- verified
  bit-exact per-branch by construction (each branch's own fused/separate
  rounding steps untouched, confirmed exhaustively over all 2^32 inputs
  for both functions, identical avg/max ulp and worst-x). But real mca
  showed *zero* movement on either function (asinh/acosh latency and
  throughput both landed on the exact same numbers as before, to the
  hundredth), and a direct `--emit=asm` diff of the compiled region
  confirmed byte-for-byte identical machine code -- LLVM's own optimizer
  already hoists the branch above the pure `sqrt` call automatically,
  same class of already-established no-op as the `wrapping_sub`-combined-
  range-compares and `koff-free unchecked-log` findings above. Reverted,
  no reason to carry the less-obvious source form for zero benefit.

### cbrt / sqrt / hypot / powf / remainder

- **Seed constant + degree-2 poly joint search**: best across 41 seeds
  still ~50x over budget (max ulp 112). 3 coefficients can't correct this
  seed's error regardless of seed choice.
- **Integer-division-free cbrt seed**: cheaper codegen (3 vs 5
  instructions), but too coarse for the correction poly to compensate —
  max ulp 33, ~16x over budget.
- **cbrt: joint seed-constant + degree-3 search**: coordinate descent and
  an independent Python sweep (±2000 seed offsets) both found ~0.3% or
  less movement — shipped combination already essentially optimal.
- **cbrt: 2/2 rational correction**: ~100x tighter in isolation, but both
  forms are already deep in f32's rounding-noise floor, and the extra
  division sits after the seed on the critical path with nothing to
  overlap (latency +28.5%, avg ulp regressed slightly).
- **cbrt degree-4 correction poly** (idea #52, the untried midpoint
  between shipped degree-3 and the previously-rejected degree-5): seeded
  from a real scipy/HiGHS minimax LP fit of `p(r) = ((1+r)^(-1/3)-1)/r`
  against the seed's real probed r-range (`[-0.0998,0.0893]` over one
  representative octave, matching `cbrt_normal`'s own documented
  octave-repeating error), polished via a new `cbrt_normal_c5` tune.rs
  target. This is a genuine, real win, not an idealized-only signal --
  confirmed on the full exhaustive 2^32 sweep (learning from this exact
  function's own prior "clean on a spot-check, reversed at the
  denormal-rescale boundary" trap, see the coordinate-descent entry
  below): `cbrt` avg ulp 0.28→0.0884, max ulp 3→**1**, `cbrt_unchecked`
  the same (0.28→0.0887, 3→1) -- comfortably inside the documented
  avg<=1/max<=2 budget for the first time (`cbrt_accurate` unaffected,
  already 0 either way). But mca showed a real, consistent cost across
  all five affected functions: latency +4.00 cyc everywhere (cbrt
  35.06→39.06, cbrt_unchecked same delta, cbrt_accurate(_unchecked)
  59.06→63.08, rcbrt 46.30→50.30, +7-11%), throughput +2.6% to +19.6%
  (cbrt_unchecked hit hardest proportionally, from the lowest baseline).
  Reverted (`git checkout -- src/lib.rs`) -- a real cost on every
  caller, even for closing a *documented* budget overage, doesn't clear
  this session's no-penalty bar the way the degree-3-over-degree-5
  choice was itself explicitly made for speed. `examples/tune.rs`'s
  `cbrt_normal_c5` scaffold is kept (coefficients:
  `-3.3333338e-1, 2.2221813e-1, -1.7280723e-1, 1.4543797e-1,
  -1.295184e-1` for `c1..c5`) in case the budget is ever judged worth
  the cost.
- **cbrt_accurate via Halley from a cheaper seed**: accuracy parity with
  the shipped Newton-based form was fully achieved (bit-for-bit, after
  fixing two overflow-ordering bugs in the correction term and lowering the
  large-magnitude rescale threshold 2^127→2^100) — the "cheaper seed can
  reach the same precision in one step" premise is true. But mca showed a
  real regression, not a win: latency 59.06→75.16 cyc (+27.3%), throughput
  essentially a wash (+1.2%). The extra division and sign-reconstruction
  cost more than cbrt_normal's own polynomial refinement step ever did.
- **Tune cbrt_throughput's magic constants**: single-octave grid looked
  like a win (max 16→12) but the function's error doesn't repeat across
  octaves like cbrt_normal's — implementing it made the real fuzz sweep
  *worse* (avg 6.73→7.01). A wide 60-octave grid retry exposed that
  `tune()`'s max-first tuple comparison wrecked the average shaving the
  worst case (avg 6.83→22.73).
- **hypot: fma pairing choice** (max-first operand pairing via
  compare+select): real, consistent avg-ulp win (~15% better) but max ulp
  never moved off 1, and mca showed a real, deterministic cost: hypot
  latency +24%/throughput +4%, hypot_checked +1.8%/+8.5%. A `.max()`/
  `.min()` variant measured *worse* still (latency +38%) — the real cost
  isn't the instruction choice but that determining operand order at all
  forces a serial step in front of the fma, which the naive
  `fma(x,x,y*y)` avoids entirely. (`.max()`/`.min()` also isn't
  correctness-neutral: they discard NaN per IEEE maxNum/minNum semantics,
  unlike the compare+select form.)
- **powf special-exponent select tier** (`y==2.0 → x*x` single-select
  probe): real accuracy win where it hits (avg ulp 5.42/max 45 → exact
  0/0), but this crate's branchless-select style computes every branch
  unconditionally, so the cost is paid on *every* call regardless of y:
  throughput +4.3%. Narrow benefit (exactly one of a handful of values)
  for a cost paid always; callers who need exact squaring can already
  write `x*x` directly.
- **powf's ~150 (later found ≥312) max ulp, compensated two-product
  multiply fix**: the multiply contributes essentially nothing —
  100% of the error traces to `log_2(ax)`'s own single f32 precision (~24
  bits) being exponentially amplified by `y` before `exp2_checked`.
  Manually applying the two-product correction recovered only ~13% of the
  error (e.g. -312→-270 ulp at the first found worst point). The only real
  fix is a higher-precision log2 in the first place — which already exists
  as `powf_checked`'s `log2_df`/`exp2_checked_df` route, at the real extra
  cost that tier already pays.
- **ln/log10 fuse trailing fma, re-tested**: duplicate of an
  already-logged idea, initially mismeasured as a win from a stale
  baseline instead of a fresh re-measurement; re-verified as a real,
  reproducible +1-cycle regression. Always re-measure a baseline fresh,
  never trust a written number.
- **Centered-variable refit for exp2's f**: backlog claimed centering
  shrinks coefficients — scipy showed the opposite (centered coefficients
  are *larger*). exp2's R(f) is monotonic with no interior minimum, so
  centering moves away from the edge minimum the current fit already
  exploits.
- **"Resurrect the f64-reduction sin/cos"**: backlog's claim that working
  code already existed in git history didn't hold up (nothing in `git log
  --all`/reflog). Built one from scratch anyway — failed on accuracy by a
  wide margin (max ulp 637 even in the smallest bucket). An f64
  two-product compensation helped (637→70) but couldn't close the gap —
  f64's own PI constant has irreducible ~2^-53 error no surrounding
  compensation fixes. Verify "the code already exists" claims via git log
  before trusting them.
- **atan_poly Horner→Estrin**: measured backwards on every axis — fuzz
  accuracy worse (avg 0.068→0.072, max 3→5), mca throughput much worse
  (+156%).
- **Rational (P/Q) refits for log_2/acos_poly**: log_2 ruled out by
  analogy (needs a division in the same fatal critical-path position as an
  already-rejected idea). acos_poly's fit failed to converge within 60s at
  two degrees — a reproducible dead end.
- **Rational (P/Q) refits for sinf_poly and erf_poly**: sinf_poly's fit
  converged but landed 5 orders of magnitude worse than the incumbent
  poly; erf_poly's fit timed out at 30s both degrees.
- **cos (fast tier) dedicated even poly in r²**: backlog itself predicted
  failure (relative error blows up near cos's zeros) — a quick scipy
  check confirmed unbounded relative error at the zero.
- **Simulated annealing / basin-hopping (`tune_basin_hop`)**, tested on
  acos_poly: coarse-grid candidate reported max ulp 2 (vs. descent's 3) —
  real fuzz showed it was actually a regression (avg 0.961 vs shipped
  0.496, nearly double). Same max-first tuple bias as cbrt_throughput's
  own finding.
- **asin: a*a-a → fma(a,a,-a)**: bit-identical, but throughput worse
  (latency unchanged).
- **Small-poly Estrin audit, asin_small/sinh_small**: not bit-exact,
  measured backwards on every axis for both functions (sinh: worse
  latency+throughput+accuracy; asin: max ulp regressed 9→10).
- **`pown` large-|n| overflow, fix attempts**: `pown(0.997296, -32767)`
  returns `inf` instead of the true finite `3.4025991e38` — 14 compounded
  squaring-rounding steps push the intermediate just over `f32::MAX` while
  the true answer sits just under it (~10.2% of a structured
  boundary-focused sweep hit this). A `Df32`-precision-only fix doesn't
  help (buys back precision, not range — a `Df32` pair is already
  `(inf, x)` once its primary word overflows). A pure exponent-tracking
  range-extension fix also fails alone (renormalizing doesn't reduce the
  number of rounding events, just relabels the same error at a different
  threshold). The combined precision+range fix (`WideFloat`: Df32 mantissa
  + tracked exponent) works correctness-wise and was eventually made to
  auto-vectorize (needs 4-iteration loop groups, not 8, and several
  codegen fixes for saturating-cast/scalar-fallback traps) — but the real
  `mca` numbers are decisive: latency ~11x (1922 vs 176 cyc), throughput
  ~30x (115 vs 3.8 cyc/elem). Fully rejected on cost, not feasibility. A
  cheaper reciprocal-only Df32 seed (option 2) barely touches the original
  bug (6/18k cases fixed) though it does halve max ulp on the existing
  documented |n|≤8/64 ranges — also separately blocked by an mca-tooling
  limitation (llvm-mca can't parse the region once LLVM hoists n's sign
  branch outside the loop). `pown` remains unfixed for large |n|.
- **remainder_wide vs remainder_checked "bit-identical" claim**: two real,
  narrow exceptions found via denser standing-test re-runs (not fixed,
  both documented and excluded from the test domain instead): (1) denormal
  rescale precision loss when `max(|x|,|y|) > f32::MAX/4` *and* `|x|` is
  already near the denormal boundary — up to ~4 ulp, 2798/29M samples; (2)
  a genuine sign flip at exact mathematical half-integer ties in x/y
  (3/292M samples) — the `adj` middle stage's blind `.round()` double-
  corrects an already-correctly-resolved tie before `remainder_checked`'s
  own correct tie-break logic runs. Neither fixed: narrow, and a safe fix
  risks breaking `remainder_wide`'s actual job of recovering large-ratio
  quantization gaps.

### Codegen / build / measurement

- **codegen-units=1 + lto sweep for the bench profile**: bit-for-bit
  identical to the default profile everywhere (one function off by 0.001
  cyc/elem, floating-point noise in mca's own arithmetic) — confirms the
  default profile isn't losing anything to a compilation-unit boundary,
  but no benefit either.
- **Combined range compares via wrapping_sub** (backlog idea #11):
  rewrote cbrt/cbrt_accurate's `ax == 0 || ax >= EXPONENT_MASK` and
  powf_checked's `axb != 0 && axb < EXPONENT_MASK` as single
  `wrapping_sub(1)` unsigned range compares. Confirmed bit-exact
  (edgecheck + full/quick accuracy fuzz, identical avg/max ulp on both
  functions) but a genuine no-op: `--emit=asm` diff of
  `examples/mca_target.rs`'s whole compiled output before/after was
  byte-for-byte identical (not just mca-equal) — LLVM's InstCombine
  already canonicalizes this exact OR/AND-of-two-comparisons-on-one-
  variable pattern into the same range check, same class as the
  already-rejected `koff-free unchecked-log fast path` and `atan2`'s
  bothzero/hpisignx entries above. Reverted (zero benefit, no reason to
  carry the less-obvious source form).
- **Literal-transcription standing test** (idea #109, "assert every
  decimal literal against its intended bit pattern"): investigated
  rather than built. Grepped for every other hand-typed decimal literal
  matching a well-known mathematical constant (pi, pi/2, ln2, sqrt2,
  log2(e), log10(e), 1/pi, etc.) anywhere in `src/lib.rs` -- found none
  besides `acos_poly`'s own leading term, and that one specific bug is
  *already* fixed and already has a standing regression pin
  (`edgecheck.rs`'s `acos(0)`/`acos(-0)` checks against
  `std::f32::consts::FRAC_PI_2`, with the bug's history documented right
  above them). The two other "bit-identical to this literal" comments
  in the file (`LOG2_COEFFS`'s leading term, `log10`'s analogous entry)
  aren't at risk the same way -- they use the named `std::f32::consts::*`
  constant directly, not a hand-transcribed decimal, so there's nothing
  to transcribe wrong. The fully generic version of this idea (checking
  *every* decimal literal in the file) isn't really buildable: the vast
  majority of this crate's literals are minimax-fitted polynomial
  coefficients with no independently-known "intended" value to check
  against -- only the narrow "well-known named constant, hand-typed as
  decimal" sub-case has a ground truth at all, and that sub-case is
  already covered.

## Untried backlog

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
- **Batch/slice API tier** (`exp2_slice` etc.): lets the crate own
  vectorization instead of the caller's loop; natural home for a fast
  sincos too. Also enables `is_x86_feature_detected` runtime
  multiversioning at the slice level (per-call dispatch is un-inlinable,
  per-slice is free).
- **Explicit core::simd fallback tier / f32x16 AVX-512 via
  #[target_feature]**: only worth it if a LUT idea ever survives
  screening, or as a per-function-scoped way to capture the AVX-512
  zmm-width win found (and rejected as a global default) above.
- **remainder_checked beyond 2^24**: double-float q like sin_checked's
  reduction. Only worth it if a real use case needs it.
- **Automated evaluation-order search per poly**: fma reassociation is a
  per-poly coin flip (acos_poly's Estrin cost accuracy, atan_poly's cost
  throughput). Could enumerate Horner/Estrin groupings and score
  automatically.
- **Standing ulp-weighted minimax fit infrastructure**: weight coefficient
  fits by 1/ulp(f(x)) instead of plain relative error as a reusable,
  built-in tool rather than a one-off per-function LP script (the ad hoc
  version of this has already found real wins for exp_pos_neg/erf_poly and
  real regressions for asin_poly/erfc — see rejected section for when it
  does/doesn't transfer).
- **FTZ/DAZ feature flag**: a cargo feature assuming caller-side FTZ/DAZ
  would let every denormal branch (log family, cbrt, exp2_checked) become
  dead code.
- **Vectorized Payne-Hanek "exact" tier** for sin/cos: full-range correct
  reduction via a 2/pi mantissa product with a per-lane variable shift.
  Big job, optional given the graceful-degradation contract.
- **Intermediate sin/cos tier (|x|≲1e5)**: single extra correction word
  over the fast tier's 4-fma Cody-Waite, well short of checked's full
  double-float q — a third point on the speed/domain curve if any user
  workload actually sits there.
- **sinh_accurate/cosh_accurate tier**: Df32 through the exp combine —
  only if a user asks; max 5 is comfortably documented already.
- **lgamma (Stirling + reflection)**: big job, listed for completeness —
  the largest gap vs. libm's function set that fits this crate's
  branchless style.
- **Domain-specific fast-math contract tiers**: a `finite-math-only` cargo
  feature gating away every inf/nan select in checked functions
  (complements the FTZ/DAZ idea above, which only covers denormals).
- **Auto-generated `_unchecked` variants via macro**: every checked/
  unchecked pair is hand-maintained; a macro emitting both from one body
  with cfg'd guards removes drift risk.
- **Bit-sliced two-for-one sincos**: evaluate sin and cos polynomials
  sharing y=r² registers across the same vector when the caller wants
  both — a sincos slice API where lane pairing amortizes the reduction.
  Only viable inside a slice tier (scalar fusion attempts already failed,
  see rejected section).
- **tanpi/tand: fuse the two independent parity computations**: currently
  pays for both `sinpi`'s and `cospi`'s own parity computation even though
  they cancel algebraically in the ratio; a shared reduction producing
  both sin/cos values from one parity computation would be more invasive
  than the current drop-in composition but could recover real throughput.
- **Stochastic rounding harness mode**: run accuracy sweeps with the final
  fma's rounding perturbed ±1 ulp to measure how close each function sits
  to a rounding boundary — identifies which maxes are "one lucky rounding"
  vs. structural, prioritizing refit targets.
- **Interval-arithmetic self-audit build**: a cfg that swaps f32 for an
  interval type in the `_normal` cores to machine-verify "this add is
  exact / Sterbenz applies" claims scattered through the comments —
  several past bugs (pre_offset, e3 sign) were exactly wrong claims of
  this kind.
- **Per-function transformed-variable fit search**: fit in u=s/(s+2),
  u=s·(s+a), etc., searching over the transform family — distinct from
  centered-variable refits (already rejected, that only moved the origin).
  A nonlinear transform changes curvature matching.
- **Ulp-staircase-aware LP grids**: densify fit grids near output
  power-of-2 boundaries where ulp weight steps 2x. Complements the
  ulp-weighted-fit idea above (that's weights; this is node placement).

### Kitchen-sink brainstorm (2026-07-11, unscreened)

~100 brainstormed ideas, deliberately not filtered for likelihood — none
tested, screened, or measured. Every one must keep the crate's branchless
auto-vectorizing contract; items with known-risky codegen are flagged
`codegen_check`. Each was checked against the rejected list above; where
an idea revisits a rejection, the differing mechanism is stated.

#### Fitting & search techniques

1. **Real-chain refit**: build the fit Jacobian by finite-differencing
   coefficient ulp steps through the *actual compiled construction*
   scored on the real fuzz, not a continuous poly model — attacks the
   "error is rounding-chain-dominated" wall that killed log_2's LP.
   Candidates: log_2, acos_poly (max 5), sinf_poly.
2. **MIP quantized fit**: HiGHS supports mixed-integer — fit with the
   f32 quantization of each coefficient as integer variables (a
   homegrown fpminimax, since sollya isn't installed).
3. **1-D exhaustive ±few-hundred-ulp scan of each poly's final combine
   constant** scored on the real fuzz — the cheap slice of #1.
4. **Degree-shed sweep with the max-capped weighted-LP machinery** (the
   rejected degree probes used plain lolremez minimax, a weaker tool):
   ln/log10 (deg 9), atan_latency (deg 17), erf_poly — can any drop a
   term while holding max ulp?
5. **Per-caller Pade refit where domains genuinely differ**: exp2m1's
   Pade argument is y=x·ln2, |y|<0.347 — narrower than expm1/tanh's
   shared |y|<0.5 fit. (The rejected decouplings — sinf_poly,
   exp_pos_neg — had *identical* caller domains; this one doesn't.)
6. **Joint threshold+coefficient coordinate descent** in tune.rs
   (crossover as a continuous search parameter) — automates the asin
   fix-5 lesson instead of retuning thresholds against frozen polys.
7. **Seam retunes not yet done**: exp2m1's 0.5, exp_m1_over_x's 0.5,
   sinh_checked/cosh_checked's 0.5 (post-checked-half construction),
   softplus/logaddexp's 87.0 cutoff, asinh/acosh's 2048 rescale
   threshold. (The 5-function crossover audit covered sinh/tanh/expm1/
   asin/erf only.)
8. **atan_poly joint numerator+denominator nonlinear refit** (scipy
   least_squares on the true rational) — only separate num-only/
   denom-only LPs were tried; the max-4 worst point was diagnosed as
   denominator-or-division-bound.
9. **erf joint Pade+erf_poly refit with explicit crossover-region
   weighting** — the rejected num-only LP regressed exactly at the 0.28
   seam, which a joint objective would score directly.
10. **exp2_q_poly combine-sensitivity LP** (weight by the
    `fma(q, t1*f, t1)` combine's local derivative) — first check where
    the real worst f sits: the technique's documented failure mode is
    the sensitivity vanishing exactly at the hard region.

#### Codegen & micro-optimizations

12. **Exponent fields from magic-round bits via integer ops**: after any
    magic-round, k already sits in kb's low mantissa bits —
    `((kb_bits + C) << 23) & EXPONENT_MASK` replaces the `(k+383)`
    float-add/shift chain with vpaddd/vpslld on less-contended ports.
    Sites: exp10/exp10_checked, exp_r_singlefield (tanh/sigmoid),
    exp2_field_split. Verify negative-k two's-complement wrap;
    codegen_check. **Naive first attempt tried and rejected 2026-07-20**:
    collapsing `exp2_field_split`'s own `k1b = fma(k,0.5,ROUND_MAGIC) -
    (ROUND_MAGIC-383.0)` into one `fma(k, 0.5, ROUND_MAGIC-383.0)` (using
    a pre-offset magic constant to skip the de-bias subtract) is simply
    *wrong*, confirmed by a direct scalar bit-comparison probe: 9999993/
    10000001 mismatches over a wide integer-k sweep, diverging as early
    as k=272 within the documented usage range. Root cause: the magic
    constant's own raw bits (`fma(...)` before de-biasing) live at the
    magic-round's ~2^23 magnitude scale, but the exponent-field-
    extraction trick (`.to_bits() << 8 & EXPONENT_MASK`) needs its input
    at *ordinary* magnitude (a plain small float like `383.0±126`) — the
    de-bias subtract isn't redundant scaffolding, it's what moves the
    value between two different bit-trick magnitude regimes, and the two
    tricks' shift amounts (8 vs. this idea's own `<<23`) reflect that
    they aren't interchangeable representations of the same bits. This
    only falsifies the naive single-fma merge, not the idea's actual
    `<<23`-based construction (never implemented/verified) — that would
    need a from-scratch derivation of what `kb_bits + C` actually encodes
    before trusting it, not an assumption that it's "the same k, just
    with fewer steps."
13. **`core::hint::assert_unchecked` range hints after clamps** so LLVM
    can prove bounds it currently can't — directly targets the rejected
    "k1 from bit-twiddled k" (killed only because LLVM couldn't prove
    the clamp bound at codegen time).
14. **`f32::to_int_unchecked` where the range is guaranteed**: lowers to
    plain vcvttps2dq (vectorizes), unlike the saturating `as i32` that
    de-vectorized exp2_checked's rejected variant. codegen_check
    mandatory, pairs with #13.
15. **lzcnt-based denormal normalization** (`u32::leading_zeros` →
    vplzcntd, AVX-512CD): replace the compare+select 2^24 rescale in the
    log family/cbrt/hypot_checked with an exact shift-based normalize.
    codegen_check that scalar leading_zeros autovectorizes here.
16. **Sign-algebra select→xor audit**: screen erf's
    `mulsign(1 − exp2(...), x)` combine and powf_sign_combine's
    6-select tree for a bit-mask classify form. (acos's analog measured
    slightly *worse*, and plain `atan`'s port of atan_latency's own
    mulsign-reassociation also measured worse — see rejected section —
    so per-site mca is required, not assumed, for each remaining site
    here too.)
22. **Toolchain-bump re-screen list**: tag the rejections that were pure
    scheduling artifacts (pre_offset dead-add removal, ln/log10
    trailing-fma fuse +1cyc, reduce_pi depth-2 rebalance) and re-measure
    after each nightly bump — these can silently flip.

#### exp / log family

23. **exp/exp_checked floor-domain reduction**: k = floor(x·LOG2_E)
    (native vroundps-down), r ∈ [0, ln2), poly refit — guarantees
    k ≤ 127 in-domain so the k1/k2 split collapses to one exp2int
    field (~5 ops deleted). Risks: doubled poly domain may need degree
    +1 (scipy-screen first); breaks exp_pos_neg's p(−r) even/odd reuse
    (that caller keeps round); edgecheck x=88.37628. Note the rejected
    round-based experiments went the *opposite* direction (round where
    floor existed).
24. **exp t1-weave revisit with different port placement** (weave into
    t2 instead, or pre-scale p): the rejected version's accuracy win
    (max 3→2, cascading to expm1/sinh/cosh/tanh) was fully real — only
    fma/mul port contention killed it.
25. **exp_r_poly degree 5→6**: attacks expm1's direct-branch max 6;
    cost one fma at already-amortized Estrin depth.
26. **exp_pos_neg even/odd poly degree bump**: sinh/cosh's round-off
    audit traced max 5 to the poly's own ~1.3e-7 fit error — degree is
    the one lever there not yet tried.
27. **exp10: treat LOG10_2_LO as a free fitted parameter** (end-to-end
    exp10 objective) instead of the rounded mathematical residual.
29. **exp10: fold LOG2_10 into a dedicated Q(d) poly in d directly**
    (the sinpi/sind constant-folding trick) — deletes the `d*LOG2_10`
    multiply and its rounding; the floor-adjust step's units need care
    (speculative).
30. **exp2m1: refit the Pade directly in x** (absorb ln2 into the
    coefficients) — deletes the `y = x*LN_2` multiply and its rounding
    from the seam-owning branch.
32. **exp_pos_neg: return halves pre-scaled by 0.5 for plain sinh/cosh
    too** (like checked_half already does) — deletes the caller-side
    `0.5*(ep±en)` multiply.
33. **Public `sinhcosh` pair function**: exp_pos_neg already computes
    both — callers needing both pay one reduction instead of two.
34. **log1p/log2p1: two_sum the `ln(u) + corr` combine** (targets log1p
    max 4, log2p1 max 3); ship as `_accurate` twins if the mca cost is
    real.
36. **rlibm-style discrete rounding-interval LP extended to ln/log10**
    (same 2^23 reduced-input multiplicity as the existing log_2 entry).
37. **log_family_wrapper: cheaper special-case classify** — derive both
    spec selects from one bits-domain compare pair; audit-grade micro.
38. **softplus fused kernel**: one fitted poly for ln(1+2^-t) over the
    k/f-reduced domain, replacing exp (full poly) → log1p (division +
    deg-9 ln poly). Big throughput candidate.
39. **logaddexp: same fused kernel** on |a−b|.
41. **logaddexp2** (base-2 sibling, ML/audio) — near-free variant of
    whatever #39 lands on.
42. **tanh division-residual correction on e/(e+2)** (round-off audit:
    error splits ~evenly poly/combine, so this targets the combine
    half). Caution: atan's analogous division fix measured zero benefit
    at +181% throughput — mca screen before any fitting work.

#### sin / cos family

43. **sinpi/cospi dedicated poly — fold π into sinf_poly's
    coefficients** (fit sin(πr) in r over [−0.5,0.5]): deletes the π·r
    multiply and its rounding at zero added ops; leading πr term via
    `fma(r, PI_HI, …)` with the PI_LO·r residual absorbed into the fit.
    Mechanism differs from the rejected two_prod correction (which
    *added* ops for the same rounding). Targets cospi's avg 0.281.
44. **sind/cosd/tand: same fold for DEG_TO_RAD_SMALL** — deletes the
    d·c multiply and captures the accuracy the rejected HI/LO split
    found (sind avg −45%) at *negative* op cost instead of +11%.
46. **parity(qh) bit-derivation**: |p0| < 2^22 → magic bits; p0 ≥ 2^24
    → deterministically even (every f32 there is an even integer); only
    the 2^22..2^24 window needs a select. Screen vs the floor-based
    parity. **Applies equally to `ql`, not just `qh`** — see the
    rejected `ql`-via-magic-round entry below: `rem` (what `ql` rounds)
    is *not* bounded the way its name suggests, so any bit-derivation
    scheme needs the same large-magnitude fallback this idea already
    plans for `qh`, on both words.
47. **Fast-tier reduction upgrade via two_prod**: replace the bounded-q
    PI_A..D 4-fma chain with one two_prod(q, PI_HI) + a PI_LO word —
    exact at any q, could push the fast tier's ~1.3e7 cliff far out at
    similar op count. Concrete design for the backlog's "intermediate
    tier" (distinct from the rejected *word-dropping* 3-word/3.5-word
    attempts, which reduced precision; this adds none of that risk).
48. **tan_checked**: sin_checked/cos_checked composition (the
    shared-reduction fusion already measured slower wall-clock; plain
    composition mirrors tanpi/tand).
49. **sinf_poly real-chain refit** (#1's method) scoring sin_checked +
    cos_checked's actual reductions jointly — the rejected LPs used
    continuous grids that mis-weighted the caller split.
50. **sinc: attribute the max-4 error** (division vs sinpi) and screen
    an fma-residual correction on the division if division-bound.
51. **tan/tand/tanpi cross-inline CSE audit**: each pays two full
    reductions differing only in the q offset — asm-diff whether LLVM
    shares the common work; hand-share at composition level if not.

#### cbrt / sqrt / hypot

53. **rcbrt direct seed**: negated-exponent bit trick + refit of the
    same deg-3 correction shape — deletes rcbrt's trailing division
    entirely (its own doc explicitly defers this as "real fitting
    work"). Targets rcbrt's 0.418/5.
55. **hypot3/rnorm3** (3-arg vector norm): fma chain + sqrt,
    graphics/physics staple, trivially vectorizes.
56. **Slice-tier FTZ/DAZ via MXCSR**: a slice entry point can set
    FTZ/DAZ around its own loop and restore — gets the FTZ
    feature-flag idea's win without a global cargo feature.

#### asin / acos / atan

60. **atan2_latency tier**: atan_latency-based atan2 — atan2 currently
    stacks atan_poly's division on top of its own y/x division.
61. **atan_bounded tier** (|x| ≤ 1 contract): skips the 1/a division
    and fold select entirely; atan_poly alone is the whole answer
    there.

#### erf family

66. **erfinv** (Giles-style poly in w = ln(1−x²), two-poly branchless
    select, fully fma-based) — sampling/ML staple, vectorizes cleanly.
67. **norm_cdf/norm_pdf pair** (Φ via 0.5·erfc(−x/√2)) — thin
    composites for API breadth.
#### hyperbolics / activations

69. **atanh restructure**: small-|x| odd-poly branch (|x|<0.25, ~5 odd
    terms) + single-log1p-on-|x| form (`0.5·log1p(2a/(1−a))` +
    mulsign) for the rest — halves the two-log1p cost. The rejected
    single-log1p variants lacked the small branch; verify that's where
    their accuracy loss actually lived before fitting.
71. **logit(p)** (sigmoid's inverse) via cancellation-safe log1p forms.
72. **compound(x, n) = (1+x)^n** (financial kernel): n·log1p-based
    route avoiding pown-on-(1+x)'s tiny-x precision loss.

#### powf / pown / remainder

73. **powf_mid tier**: 1.5-word log2 (k + one two_prod hi/lo pair, no
    full Df32 arithmetic) through the y multiply + a single
    multiplicative correction — targets powf's y-amplified error at a
    fraction of powf_checked's +61% throughput cost.
75. **rootn (C23)**: x^(1/n) with odd-n negative handling — cbrt
    generalization reusing powf/pown pieces.
76. **Constant-base powf slice**: precompute log2_df(x) once per slice;
    per-element work drops to a Df32 multiply + exp2 (slice-tier
    candidate, fixed-base workloads).
77. **pown_16 tier** (|n| ≤ 65535): 16 iterations — pown_small's 8-iter
    precedent measured ~4-6x over pown.
82. **Const-y remainder/fmod slice variant**: precompute 1/y,
    `q = round(x * (1/y))` — trades the per-element division for a
    multiply; the changed q rounding needs an accuracy screen
    (slice-tier candidate).

#### New functions / API breadth

84. **xlogy / xlog1py** (entropy kernels, 0·log(0)=0 convention via
    select).
85. **atanpi/atan2pi/asinpi/acospi** (C23 half-turn inverses): NOT
    output-scaled (double rounding) — fold 1/π into each poly's own
    coefficients; the π/2-derived constants become *exact* (0.5, 0.25),
    so these could come out *more* accurate than the radian versions
    for free.
86. **exp2i/ldexp/frexp-style exact power-of-two utilities**
    (exp2_field_split is already the core; vcvtdq2ps-friendly).
87. **Public Df32 module** (log2_df/exp2_checked_df/two_prod etc.) for
    power users composing their own accurate kernels.
88. **Public checked pi-reduction API** (round_x_over_pi + reduce_pi)
    for user compositions (custom periodic kernels).

#### Infrastructure / harness / tiers

89. **Gather-based LUT tier screen**: a true vgatherdps table + short
    poly for exp2 (16/32-entry exact 2^(j/N) hi words) — distinct from
    the rejected compare-select tree; screen the gather's mca cost
    first; likely lives in the simd/slice tier.
90. **Same gather screen for log_2** (mantissa-segment table +
    low-degree poly).
91. **Slice tier × rejected-global-flags interaction**: per-function
    zmm-width and interleave=2 re-tests become possible via
    #[target_feature] once the slice tier exists — two documented
    median-win/outlier-veto rejections become recoverable.
92. **simd-tier division-free kernels via vrcpps/vrsqrtps + NR**
    (atan_poly's denominator, sigmoid/tanh's final division, rcbrt, an
    rsqrt_fast) — unreachable from scalar autovectorized code, natural
    once an explicit core::simd tier exists.
93. **2-arg importance-sampling harness** for powf/atan2/hypot/
    remainder: structured lattices near known-hard manifolds (e.g.
    y·log2(x) near integers) — better worst-case discovery than
    uniform sampling.
94. **Differential testing against sleef/rlibm-built references** for
    functions with no std twin (erfcx, sind, exp10, sinpi, …).
95. **accuracy.rs per-branch attribution mode**: report which select
    arm produced each worst case — speeds every future refit's
    diagnosis step.
96. **Thermal-controlled wall-clock harness** (cpufreq pinning,
    alternating A/B batches): de-risks the one recorded mca-vs-
    wall-clock *direction* disagreement (sincos_checked) and the
    documented 3-9x environment swing.
97. **Targeted LLVM flag screen**: -enable-unroll-and-jam,
    -extra-vectorizer-passes, SLP horizontal reductions — cheap sweep,
    same method as the interleave/zmm experiments.
98. **Auto-`_unchecked` macro (existing entry) — concrete new
    candidates**: sind/cosd (drop the POLY_SAFE_BOUND clamp under the
    4.7e7 contract), erf (drop the 10.0 bound), softplus/logaddexp
    (drop NaN guards), sinpi (drop the x==0 select).
99. **tgamma** companion to the lgamma entry (Lanczos/Stirling, shares
    machinery).
100. **Bessel j0/j1** (Cephes-style two-region rational + trig
     composition — big job, listed for completeness like lgamma).

#### Batch 2 (same session): fitting & search, continued

101. **Avg-first, max-capped simulated annealing** over coefficient ulp
     space — the rejected basin-hop failed on its max-first tuple
     objective (documented bias), not on the annealing mechanism itself.
102. **LLL/lattice reduction over the coefficient quantization step**:
     finds good *simultaneous* f32 roundings of a whole coefficient set
     — the cheap cousin of the MIP idea (#2).
103. **Low-discrepancy (Sobol) fit/verification grids** — avoids uniform
     grids' aliasing against ulp staircases; complements the
     staircase-node-placement entry (weights vs placement vs sequence).
104. **True Remez exchange in-repo** (f64, certified equioscillation) —
     the LP is a discretized stand-in; Remez gives certificates and
     better conditioning on rationals (atan_poly, erfc_rational).
105. **Gappa (or hand-rolled interval) certificates** for the crate's
     exactness claims (Sterbenz subtractions, exact po2 multiplies,
     k*LN2_HI) — the machine-checkable version of the interval
     self-audit entry; several past bugs were exactly wrong claims of
     this kind.
106. **Caller-profile-weighted alternate coefficient sets** behind a
     cargo feature (e.g. sin weighted toward [−2π,2π]) — same shapes
     and cost, different literals.
107. **PI_A..D bit-allocation joint refit**: keep all 4 words but
     re-search each word's trailing-zero budget jointly with downstream
     error — distinct from the rejected word-*dropping* attempts (3-word
     chain, fitted 3.5-word split), which reduced total precision.
108. **Clenshaw/Chebyshev-basis evaluation screen** for the worst-
     conditioned polys (erfc_rational's n/d, atan_latency's deg-17) —
     different rounding structure at comparable op count.
110. **±few-ulp exhaustive scan of every non-poly literal** (clamp
     bounds, seed constants, magic offsets, branch thresholds) scored on
     the real fuzz — #3's sibling for non-coefficient constants.

#### Batch 2: exp / log family

111. **expm1_checked / exp_m1_over_x_checked**: input clamp + single
     exponent field (the tanh/sigmoid pattern) — likely cheaper than the
     current k1/k2 split *and* total-domain, since the clamp caps k at
     127 by construction.
112. **exp single-field narrow tier**: contract x ≤ ~88.02 (k ≤ 127
     under round) — deletes exp2_field_split's second field for a
     slightly narrower documented domain.
113. **exp10_checked clamp consolidation**: re-derive the x clamp bound
     so the separate k clamp is provably redundant and deletable.
114. **FTZ-mode minimal exp2_checked/exp_checked** (rides the MXCSR
     slice-tier idea #56): lower clamp −151→−126 and the
     denormal-rounding half of the split's job disappears; same cascade
     deletes denormal_rescale from the log family and cbrt — scope #56
     to capture all of it.
115. **exp_pos_neg: p_neg from the reciprocal identity**
     p(r)·p(−r) = e² − (r·o)²: one division (idle divider, cbrt
     precedent) instead of the second combine — internal to
     exp_pos_neg, distinct from the rejected tanh-via-ratio route
     (which needed two whole polys). Accuracy screen first.
116. **Public exp2_kf(k, f) pre-reduced primitive** — exp10, powf, and
     user custom-base kernels skip the redundant floor/frac.
117. **exp_scaled(x, s) = e^x · 2^s** with s folded into the field
     split for free — softmax/normalization building block.
118. **logsigmoid(x) = −softplus(−x)** as a thin public function (loss
     kernels); rides the fused-softplus idea (#38) if that lands.

#### Batch 2: trig

122. **sin2pi/cos2pi/tan2pi** (full-turn arguments, DSP phase
     accumulators): literally sinpi(2x) — the doubling is an exact po2
     multiply — thin but discoverable API with the exactness documented.
123. **asind/acosd/atand/atan2d**: fold 180/π into the poly/combine
     constants (#85's mechanism) — 90.0/45.0 are exact where π/2 wasn't.
124. **π-constant hi/lo splits in the inverse-trig combines**: acos's
     `+PI`, atan2's correction, asin's big branch — distinct from the
     rejected 1/a *division*-rounding fix (different rounding source at
     the same location). **atan's own `FRAC_PI_2 − y` fold slice tried
     and rejected** (2026-07-20): `(FRAC_PI_2 - y) + FRAC_PI_2_LO` where
     `FRAC_PI_2_LO = PI_LO/2` (exact halving of this file's own PI_HI/
     PI_LO Cody-Waite split, verified `PI_HI/2` is bit-identical to
     `FRAC_PI_2` first). Predictable in hindsight: `FRAC_PI_2`'s own
     rounding error vs. true pi/2 is only ~4.37e-8, under 0.2 ulp at this
     magnitude -- already below atan_poly's own ~0.068-avg-ulp fit-error
     floor, the same "constant precision doesn't matter, poly fit error
     already dominates" pattern as several already-rejected LP-refit
     entries elsewhere in this file. Real measurement confirmed the
     prediction: atan avg ulp 0.0675→0.0673 (noise), max unchanged at 3
     (the reported worst case sits in the *other* branch, `a<1`, never
     touched by this fold at all); atan2 avg slightly worse
     (0.0681→0.0683) with max possibly 3→4 (small 10M-sample signal,
     not chased further given the cost side already killed it). mca
     showed a real, consistent cost: both atan and atan2 latency
     +4.0 cyc (61.09→65.09 / 61.28→65.28, +6.5%), throughput +3.2%/
     +11.9%. Reverted. The other three sites (acos/atan2/asin) remain
     untried.
125. **Integer-domain parity pipeline end-to-end** for
     sin_checked/cos_checked (parities as bits, XOR combine, direct
     sign mask) — composes #45/#46; deletes the float compare+select
     flip.
126. **wrap_pi(x) → [−π, π] public angle normalization** riding
     round_x_over_pi/reduce_pi (robotics/geo staple; user-level face of
     #88).
127. **sin_prereduced/cos_prereduced public** (r ∈ [−π/2, π/2]
     contract = sinf_poly + documented parity conventions) — for
     callers who already did their own reduction; the trig analog of
     #116.
128. **tanpi/tand direct poly**: tan(πr) poly on |r| ≤ 0.25 + cotangent
     reflection for the rest — every radian direct-tan attempt died on
     *inexact reduction* near poles; tanpi/tand have *exact* reductions,
     which removes precisely that documented blocker. Targets the
     sinpi/cospi division and the near-pole ulp blowup.
129. **sinpi: re-screen the x==0 guard after the poly fold (#43)** —
     the fold changes which op erases the sign of zero.
130. **sind/cosd: clamp on d instead of the radian product after #44**
     (|d| ≤ 90 exact-domain form) — cheaper/earlier bound.
131. **sinc_unnormalized** (sin(x)/x in radians, x=0 → 1) — DSP/physics
     sibling of the existing normalized sinc.

#### Batch 2: roots / hypot / geometry

133. **Half/third-power composition family**: x^(3/2) = x·sqrt(x),
     x^(2/3) = cbrt(x)² (or cbrt(x·x)) as provided functions — each 1-2
     correctly-rounded ops, always beats powf; also a docs pointer so
     users stop calling powf for these.
134. **rnorm4 / quaternion normalize** — companion to hypot3 (#55).
136. **normalize2/normalize3 slice kernels** (rhypot + scales — the
     operation users actually want hypot for).
137. **cbrt_throughput status decision**: document it as the approx-tier
     member it is (5.5 avg ulp) or drop it — currently in limbo with no
     doc comment.

#### Batch 2: erf / inverse-trig follow-ups

138. **Dawson function F(x)** — erfcx sibling (spectroscopy), shares the
     rational machinery.
139. **erfc_inv / probit (normal quantile)** alongside #66's erfinv —
     completes the sampling stack.
140. **erfc saturation-threshold tightening screen**: find where erfc
     first rounds to exactly 0/2 — if it's meaningfully inside |x|=10,
     a tighter clamp shrinks the rational's required fit domain and
     frees headroom.
142. **acos max 5 as an explicit real-chain-refit (#1) target** — its
     coordinate descent is exhausted, which is exactly the case #1
     exists for.
143. **atan2_pos: [0, 2π) variant** via a branchless +2π fold —
     geo/graphics convention ask.

#### Batch 2: hyperbolics / ML / graphics

144. **coshm1(x) = 2·sinh²(x/2)** — cancellation-safe cosh−1
     (catenary/relativity), reuses sinh_small.
145. **log1pmx(x) = log1p(x) − x** — statistics primitive
     (gamma/Poisson kernels) that needs a dedicated cancellation-safe
     form; near 0 it's −x²/2·(1−…), a small dedicated poly.
146. **srgb_to_linear / linear_to_srgb**: the piecewise-2.4-power sRGB
     transfer pair, branchless select over the linear toe — heavy
     real-world vectorization demand; powf composition first, dedicated
     poly if it earns it.
147. **smoothstep/smootherstep graphics pack** (trivial fma polys with
     exact endpoints; cheap API breadth).
148. **Softmax / logsumexp / normalize slice reductions** (max-pass +
     exp-pass + sum + scale in one fused traversal) — slice-tier
     flagship, plus a rotate2d(sincos) demo kernel.
149. **Monotonicity audit**: sigmoid/tanh/erf/atan/softplus checked for
     local non-monotonic ulp wiggles (breaks bisection for ML users);
     repair only if found and cheap.
150. **sigmoid_grad / tanh_grad fused pairs** (s·(1−s) reusing the
     already-computed e) — screen whether fusion beats the caller's own
     two ops before building.

#### Batch 2: powf / pown / remainder

154. **Const-arg slice family generalization**: const-y remainder
     (#82), const-base powf (#76), const-base log (precompute
     1/log2(b)) — one shared design decision.
155. **powf slice integer-y dispatch**: slice checks all-y-integral once
     and routes to the pown path — per-slice dispatch is free where
     per-call isn't.
156. **remainder_ieee as the documented default recommendation** —
     cheaper (native vroundps) *and* standard; `remainder`'s ties-away
     is an inherited port convention, not a design goal.
157. **remainder_checked/remainder_wide consolidation screen** after the
     ties-even fix (#80): can one tier serve both contracts, or does
     wide's 6.5x cost keep them split? Screen only.

#### Batch 2: AVX-512 simd-tier instruction ideas
(all unreachable from scalar autovectorized code — natural once the
core::simd tier exists; each replaces multi-op scalar idioms)

158. **vgetexpps/vgetmantps log core**: exponent + mantissa extraction
     in two instructions with denormals handled natively — replaces the
     whole wrapping_sub bit-trick *and* denormal_rescale in a simd-tier
     log family.
159. **vscalefps exp core**: x·2^k in one instruction with correct
     overflow/underflow/denormal semantics — replaces exp2_field_split
     and most of its clamp machinery.
160. **vreduceps/vrndscaleps**: fraction extraction (x − round-to-scale)
     in one instruction — replaces floor+subtract in the exp2-family
     reductions.
161. **vfixupimmps**: table-driven special-value patching (zero/inf/nan
     selects in one instruction) — collapses log_family_wrapper's select
     chain.
162. **vrangeps** for clamp pairs (single-instruction bounded
     magnitude).
163. **vpermi2ps in-register 32-entry LUTs** (two zmm registers hold the
     whole table, no memory gather) — revisits the LUT idea (#89)
     without vgatherdps' latency; the modern fast-table technique.

#### Batch 2: harness / verification

164. **Special-value matrix v2**: systematic ±0/±inf/NaN in/out matrix
     for every public function as a standing test — five ±0 bugs found
     ad hoc so far (acos, atan2, sinf_poly, sinpi, remainder).
165. **Saturation-boundary pins**: every clamp constant and overflow
     threshold gets an edgecheck pin at ±1 ulp around it — the
     exp10_checked overflow-at-the-boundary pattern, systematized.
166. **Denormal-output correctness audit**: which functions produce
     correctly-rounded denormal outputs vs garbage (exp2_checked
     documents its behavior; most others are unaudited).
167. **Identity-consistency fuzz**: sin²+cos²≈1, cosh²−sinh²≈1, tanh vs
     sinh/cosh, exp(ln x)≈x with documented tolerance bands — cheap
     cross-function bug detector.
168. **Worst-case corpus regression gate**: persist each function's
     known worst-x list, re-check every commit in seconds between the
     hours-long full sweeps.
169. **ULP-error histogram artifacts** per function (not just avg/max)
     — bimodal structure reveals branch-split opportunities.
170. **Worst-pocket auto-bisection**: given a fuzz argmax, exhaustively
     map the surrounding error pocket's shape and width — refit
     diagnosis tool.
171. **f16-lattice smoke gate**: all 65536 f16 values promoted to f32
     against the f64 reference for every 1-arg function — sub-second CI
     sanity check.
172. **wgpu/GPU compute sweeps** for 2-arg functions — makes the
     importance-sampling lattices (#93) orders of magnitude denser.
173. **Round-trip contract measurement**: published ulp bounds for
     exp(ln x), powf(powf(x,y),1/y), sin(asin x) pairs.
174. **Auto-generated rustdoc accuracy tables from harness output** —
     the readme's quickbench numbers already drifted once; generated
     docs can't go stale.
175. **NaN-payload/quietness propagation matrix** (which ops
     canonicalize payloads) — documentation-grade completeness.

#### Batch 2: portability / infrastructure

176. **no_std/core-only feature**: most rounding already uses magic-add
     tricks; audit the residual std surface (floor/round/trunc/sqrt) and
     gate via core intrinsics or libm fallback.
177. **C ABI export layer** (#[no_mangle] extern "C") — drop-in libm
     comparison target and FFI consumers.
178. **NEON/aarch64 re-audit**: fma is native there, but every
     mca-derived scheduling decision in this crate is Tiger-Lake-
     specific — the decided tradeoffs (division-vs-poly, Estrin
     groupings) need re-measuring before claiming portability.
179. **WASM relaxed-simd gate** (f32x4.relaxed_madd): without it the
     fma compile_error! fires — document/feature-gate the story.
180. **f64 sibling module** — the whole architecture transfers, polys
     refit at higher degree (big job, lgamma-class, listed for
     completeness).
181. **strict-ieee cargo feature**: swaps the documented convention
     divergences (remainder's ties-away, fmod's uncorrected quotient)
     for slower std-matching forms — escape hatch instead of a doc
     caveat.

#### Batch 2: misc API / numerics

182. **Compensated-Estrin generic infra** (EFT-based poly evaluation) as
     reusable machinery for future _accurate tiers — compensated-Horner
     was hand-rolled once (erfc, max stayed flat); infra makes the next
     attempt nearly free to run.
183. **Full Df32/Df32 division primitive** (div_to_f32 exists) — needed
     by future rational _accurate tiers.
184. **Public EFT toolkit**: two_prod/two_sum/quick_two_sum + mulsign
     (with the mulsign-vs-copysign semantics doc) — users keep
     reinventing these wrong.
185. **fast_round_int public** (the ROUND_MAGIC idiom with its |x|<2^22
     contract documented) — for callers building their own reductions.
186. **Complex pack**: cexp/clog/cabs/carg composites (SoA-friendly,
     mostly existing kernels).
187. **n-ary logaddexp slice reduction** (tree or max+sum-exp) — pairs
     with #148.
188. **_approx tier promotion**: exp2_approx/log2_approx/rsqrt_approx
     exist untabulated — publish error bounds, add sin/sigmoid/tanh
     members for ML-inference users, explicitly outside the 0.5/2
     budget (cbrt_throughput's tier, done properly).
189. **periodic_poly! dedup macro** once #43/#44 land (four
     near-identical folded-constant sinf_poly variants) — same
     macro-not-fn pattern as pi_reduce_and_poly!.
190. **#[doc(alias)] C-name annotations** (expf, atan2f, sincosf…) —
     zero-cost discoverability.

#### Batch 2: speculative / process

191. **PWL+correction sigmoid_fast inference tier**: hard-clamped
     piecewise-linear base + one poly correction — approx-tier member
     (#188), ML-inference latency play.
192. **Caller-side FTZ/DAZ behavior test**: callers often run with FTZ
     set globally; document and test what each denormal-handling path
     actually does under inherited MXCSR state.
193. **Rayon-parallel accuracy sweeps**: the exhaustive 2^32 runs are
     embarrassingly parallel — hours → minutes changes what's feasible
     to verify per idea.
194. **Multi-start coordinate descent in tune.rs** (N ulp-perturbed
     seeds around the LP solution, keep best) — cheap robustness
     against the single-seed local-optimum traps already documented.
195. **Per-function error-budget ledger**: reduction X + poly Y +
     combine Z ulp, from the round-off audits (several exist ad hoc for
     expm1/sinh/tanh/erfc) — makes attack selection data-driven instead
     of re-derived each session.
196. **NaN-propagating min/max helpers** (`max_nan_prop`) + audit — the
     IEEE maxNum NaN-discard trap has now bitten softplus,
     hypot_checked, and logaddexp separately.
197. **Seam continuity standing test**: value + one-sided-slope jump
     measured at every branch seam (expm1/exp2m1/sinh's 0.5, asin's
     0.25, erf's 0.28) — a seam regression detector for future refits.
198. **Fix the mca latency harness's mix() sign blindness**: mix()
     erases the sign bit each chain hop, so sign-dependent work
     (cbrt, sin_checked's flips, erfcx's branch) is silently deleted
     from every latency number — inject alternating sign into the chain
     instead. Directly repairs a documented harness defect.
199. **sigmoid one-sided evaluation**: compute e = exp(+|x|) ∈ [1,∞)
     (k ≥ 0, single field trivially, no denormal side), s = 1/(1+e) ∈
     (0, 0.5], select s or 1−s by sign — 1−s is safe here (result ≥
     0.5, absolute error bounded by ulp(0.5)), unlike the rejected
     0.5+0.5·tanh identity whose cancellation lived on the small side.
200. **Auto-tune CI loop**: a scheduled job re-runs the tune.rs
     coordinate descent (LP-seeded) on every poly and files a PR when a
     real fuzz-verified improvement appears — automates the crate's
     single most-repeated manual win pattern.
