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
  earlier nor later cutoff helps in any case. `asin`'s 0.25 later *did*
  become worth revisiting after `asin_small`'s own minimax refit shifted
  its branch's error curve (shipped, see the 0.25→0.27 crossover-shift
  entry elsewhere in this file) — this entry's finding for `asin` was
  genuinely stale by the time it was re-checked, not a false alarm.
  **Re-checked `sinh`'s 0.5 for the same reason** (2026-07-20, after
  `sinh_small`'s own later degree-3→2 drop, a different but also
  real change to that branch): unlike `asin`, this one's premise still
  holds — `sinh`/`sinh_throughput`/`sinh_checked`'s real fuzz worst-case
  points sit at x≈3.16/0.85/-7.32 respectively, all deep inside the
  `|x|>=0.5` `exp_pos_neg` branch, nowhere near the 0.5 boundary the way
  `asin`'s worst case sat right at its own crossover. No threshold sweep
  run — the worst-case-location evidence alone already rules out a
  boundary-shift fix here, the same screening step that made `asin`'s
  fix findable in the first place.
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
- **Monotonicity audit** (idea #149): built a scratch checker
  (exhaustive ulp-walk over `[-4,4]` -- properly handling that negative
  floats' *raw bit pattern* decreases toward zero, the opposite of the
  positive side -- plus a coarse 50M-point dense sample over each
  function's wider documented domain) for `sigmoid`/`tanh`/`erf`/`atan`/
  `softplus`. First pass reported alarming "violations" (e.g. sigmoid
  apparently dropping from 0.545 to ~1.2e-7 between adjacent ulps) that
  turned out to be a bug in the checker's own print statement (mislabeled
  the drop *magnitude* as the neighbor's actual value) -- direct isolated
  verification (`sigmoid(0.18189578)=0.54534906`,
  `sigmoid(0.18189579)=0.54534894`) confirmed the real neighbor value is
  nearly identical, not catastrophic. Fixed the checker and re-ran: all
  five functions do have genuine local non-monotonicity, but only at the
  1-2 ulp level (worst drops 8.9e-8 to 2.4e-7, i.e. exactly 1-2 ulp at
  that magnitude) and confined to a small fraction of the domain
  (0.015%-0.19% of ulp-walk steps). This is the ordinary, expected
  behavior of any polynomial/rational approximation to a smooth function
  evaluated in finite precision, not a coding defect -- and per the
  idea's own "repair only if found and cheap" bar, a fix would need a
  fundamentally different monotonicity-preserving construction (not a
  cheap patch), so not pursued. No code changes; the scratch checker
  (`examples/monotonicity_scratch.rs`) was deleted after use, not kept
  as permanent harness infrastructure.

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
- **exp_r_poly degree 5→6** (idea #25): unlike the `exp_pos_neg` degree
  bump above (rejected on conflicting evidence without a full
  implementation), this one's premise held up on every screen -- the
  audit above shows real, non-dominated truncation error to attack (not
  a "near precision floor" diagnosis), and a real scipy/HiGHS minimax LP
  fit found a genuinely strong ~33x idealized margin (comparable to
  cbrt's own degree-4 bump, the strongest-margin case this session).
  `exp_r_poly!` is a macro, not a function, shared inline by `exp`,
  `exp_checked`, `expm1`, `exp_m1_over_x`, and `sigmoid` -- wiring the
  bump in touches all five at once. Real quick-fuzz result: a genuine,
  substantial, consistent win on every one of them -- avg ulp improved
  9-37% (exp 0.0744→0.0473, exp_checked 0.0389→0.0246, expm1
  0.1304→0.1189, exp_m1_over_x 0.0729→0.0626, sigmoid 0.0925→0.0835) and
  max ulp improved on four of five (exp 3→2, exp_checked 3→2, expm1
  6→4, exp_m1_over_x 6→5; sigmoid unchanged at 4). But mca showed the
  same real, broad cost this session's other two degree bumps
  (erfc_rational, cbrt) already established as the pattern: latency
  +4-5 cyc (+6-9.5%) and throughput +5-23% on all five, `sigmoid`
  hit hardest (throughput +23.2%). Reverted (`git checkout --
  src/lib.rs`) -- a third same-shape confirmation that this crate's
  degree-bump lever reliably buys real accuracy at a real, unavoidable
  cost spread across every caller of the shared poly, never a clean win
  under this session's no-penalty bar (see the degree-bump-costs-real
  memory). `examples/tune.rs`'s `exp_r_c6` scaffold kept for reference.
- **sinh round-off budget audit**: rounding dominates truncation (3.055 vs
  1.670 ulp) but traces to the poly's own ~1.3e-7 relative fit error
  (already near f32's precision floor) amplified ~8x by the exponent-field
  reconstruction, not a fixable single rounding step.
  **Idea #26 ("exp_pos_neg even/odd poly degree bump: degree is the one
  lever there not yet tried") directly contradicts this audit's own
  diagnosis** — re-checked rather than assumed correct either way:
  a real scipy/HiGHS minimax LP fit of a degree-7 even/odd split (one
  more term each side) against `e^r` over the reduction's real `|r| <=
  ln2/2` range found only a ~5.9x idealized margin over an f32-emulated
  reconstruction of the shipped degree-5 form, and the odd side's own
  extra coefficient converged to exactly 0 (no freedom needed there,
  matching cbrt's earlier degree-bump finding where only one of two
  "branches" used the added freedom). A ~5.9x margin sits below the
  order-of-magnitude-plus threshold that reliably predicted a real
  signal elsewhere this session (see the LP-margin-strength memory), and
  the *pre-existing, independent* round-off audit already explicitly
  diagnosed this exact poly as rounding-dominated, not fit-dominated --
  the same diagnosis class that sank the sinpi/cospi and sind/cosd
  dedicated-poly refits outright (both real regressions despite similar
  or stronger margins, see the sin/cos family section). Not implemented
  given this converging evidence; idea #26's own premise doesn't survive
  scrutiny against the audit it claims to be extending.
- **exp10/exp2m1: fold the reduction's scaling constant directly into a
  dedicated poly/Pade** (ideas #29/#30, "the sinpi/sind constant-folding
  trick" applied to `exp2_q_poly!`/`pade_expm1_ratio!`): rejected on
  inspection, not full implementation, given the direct precedent this
  exact mechanism already set twice this session -- `sinpi`/`cospi`
  (idealized margin 2.6x) and `sind`/`cosd` (4.9x) both real-regressed
  outright when their intermediate-rounding step was folded away this
  same way. Screened idea #30 (`exp2m1`'s Pade, absorbing `LN_2` into
  the rational directly instead of rounding `y=x*LN_2` first) with the
  same idealized-margin check before committing to a verdict: only
  ~5.6x, squarely in the same weak-to-moderate band that already failed
  twice, not the order-of-magnitude+ margin that held up for cbrt/
  erfc/exp_r_poly. Given `exp2_q_poly!`/`pade_expm1_ratio!` are each
  shared by several downstream functions (broadening the blast radius
  further), not implemented -- the pattern is now strong enough that a
  third confirming failure adds little beyond what's already known,
  and a weak-margin case is exactly where this mechanism has
  consistently gone wrong.
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
  **Re-tested (idea #24, "different port placement"): both proposed
  variants re-fail, and a fresh direct probe re-confirms both original
  findings independently.** The overflow claim: reconfirmed with a
  standalone f32 probe (not the crate's own code) stepping through
  `exp`'s reduction near its domain ceiling -- at `x=88.376274`,
  `t1=t2=1.8446744e19` and `t1*t2` alone is `inf`, while the real
  interleaved `p*t1*t2` correctly lands on `2.4061984e38` (finite);
  "pre-scale p" is a real correctness bug, not just a paper concern. The
  port-contention claim: re-derived the exact op-count algebra by hand
  (counting the reduction's own ops separately) -- both the original
  direct form and the weave form cost the same 9 total ops (4 mul + 4
  fma + 1 add vs 4 mul + 5 fma), the weave just trades the free `add`
  for a 5th `fma`, matching "lateral shift onto busier ports" exactly.
  Swapping which of `t1`/`t2` gets woven in (idea #24's specific
  proposal) doesn't change this trade at all -- they play symmetric
  roles in the final combine, so whichever one is folded into the poly,
  the total op count and fma-port pressure are identical either way.
  Tested empirically anyway rather than trusting the symmetry argument
  alone: first attempt (two new standalone `#[doc(hidden)]` scratch
  functions, weave-into-t1 and weave-into-t2, coexisting in the same
  `mca_target.rs`) measured throughput 1.274/1.303 -- a real-looking
  *improvement*, contradicting the original entry. Directly replacing
  `exp`'s own body with the identical weave logic, tested in isolation
  (no coexisting sibling scratch fn), reproduced the original 1.393
  regression exactly, three separate times. The standalone-siblings test
  was a methodology artifact, not a real signal -- see
  [[jodiemath-mca-coexisting-scratch-artifact]]. Not shipped, on both
  counts, exactly as originally found.
- **exp10 third Cody-Waite reduction word**: already tight (avg
  0.0343/max ulp 2, 2.2B+ samples) — minimal headroom for a third word to
  collect, not worth the extra fma.
  **Idea #27's "LOG10_2_LO as a free fitted parameter"** is the same
  lever from a different angle (perturbing the existing second word's
  own value instead of adding a third) -- checked directly (zero added
  ops, so cheap to just try): a ±1 to ±8-ulp sweep of `LOG10_2_LO` (real
  fuzz each time, not tune.rs's grid) left `exp10`/`exp10_checked` at
  the same avg ulp (~0.0306-0.0307) and max ulp (1) regardless of which
  nearby value was used. Confirms the same "already tight" conclusion
  via a different probe; reverted (no perturbation tried moved anything).
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
  **Retested 2026-07-27 via ideas #13/#14 (`to_int_unchecked`), which
  exist specifically to rescue this rejection — the de-vectorization is
  genuinely fixable, and the idea still loses. Now rejected on cost, not
  codegen.**
  - `f32::to_int_unchecked` lowers exactly as #14 predicted: the region
    gets **2 packed `vcvttps2dq`** and *zero* scalar converts, with
    packed `vpsrad`/`vpsubd`/`vpaddd`/`vpslld` for the field math, and
    `codegen_check` passes all 148 regions. So the "LLVM can't prove the
    bound" obstacle is real but avoidable — this closes that question.
  - NaN is the subtlety worth recording: `to_int_unchecked` on NaN is
    **UB**, and `.clamp()` does *not* remove NaN (Rust's clamp propagates
    it). Handled by feeding the integer path a NaN-quashed copy via
    `k.max(-151.0).min(128.0)` (Rust's `f32::max`/`min` return the
    non-NaN operand), while NaN still reaches the result through `f`/`q`
    exactly as the shipped version already relies on. Verified:
    `exp2_checked(nan) = NaN` and all edgecheck pins pass.
  - But it is a real regression everywhere: `exp2_checked` throughput
    1.399 -> **1.584 (+13.2%)**, latency 43.06 -> 48.06; `erfc` 2.437 ->
    **2.651 (+8.8%)**; `powf` 5.651 -> 5.714 (+1.1%). (`erf`,
    `exp10_checked` and unchecked `exp2` unchanged — they don't route
    through this body.)
  - Root cause, which also predicts the outcome without measuring: the
    magic-round path is **8 ops** (`fma`+`sub`+`add`+`sub`, then
    2x`shift`+`and`) versus the integer path's **9** (`max`+`min`+`cvt`,
    then `>>1`+`sub`+2x(`+127`+`<<23`)). The float trick bakes the
    exponent bias into the magic constant so `<<8 & MASK` extracts a
    ready-biased field, while a pure-integer path must add `127` per
    field explicitly. Same lesson as idea #12's own entry: the de-bias
    arithmetic isn't redundant scaffolding, it's what makes the field
    extraction free. The idea #12 hope that moving work onto less
    contended integer ports would pay for the extra op does not
    materialize here — `vcvttps2dq` competes for the same ports as the
    FP ops it replaces.
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
- **softplus/logaddexp\_unchecked, drop the trailing NaN guard** (idea
  #98, the same "unchecked tier drops a redundant guard" mechanism that
  shipped for `sind_unchecked`/`cosd_unchecked`/`sinpi_unchecked`):
  rejected on inspection before any implementation, same style as the
  `erf_unchecked` rejection above (§asin/acos/atan/atan2) but a
  different concrete mechanism. Verified directly (a tiny standalone
  probe, not just reasoning): Rust's `f32::max`/`min` are NaN-*avoiding*
  (IEEE754 minNum/maxNum semantics) -- `f32::NAN.max(0.0)` returns
  `0.0`, `f32::NAN.min(87.0)` returns `87.0`, not NaN. Both `softplus`
  and `logaddexp` route their NaN-propagation entirely through this
  guard specifically *because* their own bodies use `.max()`/`.min()`
  for the ordinary control flow (`x.max(0.0)`, `d.min(87.0)`), not
  through comparisons that fall through to arithmetic the way most of
  this crate's other wrapper guards do. Dropping the guard doesn't risk
  a rare edge case or a merely-imprecise result -- it silently returns
  an ordinary, plausible-looking finite number for the single-NaN case
  (`softplus(NaN)` would compute as if `x` were absent entirely, and
  `logaddexp(NaN, b)` as if `a` were absent, both landing near `corr`/`b`
  respectively) instead of the NaN every other function in this crate
  guarantees. Unlike `sind`/`cosd`/`sinpi`'s own clamp/guard removals
  (provably no-ops for any in-domain input, narrowing only the *range*
  contract), this one changes behavior for an input class (NaN) with no
  natural "in domain" concept to narrow around. Not implemented.

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
- **log1p/log2p1: two_sum the `ln(u) + corr` combine** (idea #34):
  rejected on inspection, same proof as the softplus/logaddexp two_sum
  entry (idea #40, §exp/exp2 family above) — `two_sum(a,b)`'s `s` component is
  bit-identical to plain `a+b` by construction, and the idea never uses
  the auxiliary error term for anything, so swapping the final
  `ln(u)+corr` for `two_sum(ln(u), corr).0` cannot change a single output
  bit. Salvaging the idea for real would need `ln(u)` itself computed at
  more than single-f32 precision (a genuine double-float `ln`, the "ship
  as `_accurate` twins" half of the idea's own framing) so there's an
  actual low-order error term worth preserving through the combine —
  that's a much bigger undertaking than a one-line `two_sum` swap, not
  attempted here.
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
- **reduce_pi: downgrade the last remaining full two_sum** (`p3t`/`e3t`'s
  merge, `two_sum(p3, tier2)`) **to quick_two_sum**: same "diverges by
  caller" pattern as the e3 downgrade above. mca: sin_checked throughput
  improved slightly (5.311→5.226, -1.6%) but cos_checked's *regressed*
  (4.603→5.234, +13.7%), both latencies flat. In-domain accuracy looked
  clean on a quick fuzz (all four documented buckets matched baseline
  for both functions), consistent with this merge's own doc comment
  flagging `e3t` as "large, not negligible" near cos's zero crossings —
  the accuracy risk that comment warns about wasn't tripped here (the
  sign/subtraction order was untouched, only the merge's own exactness
  was relaxed), but the throughput cost lands the same way it did for
  the e3 downgrade regardless. Reverted -- a real win for one caller
  paired with a real, larger-magnitude loss for the other isn't a clean
  win. This was the crate's last full `two_sum` call; `two_sum` itself
  would need removing (dead code) if this were ever adopted.
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
- **sinpi/cospi dedicated poly, fold π directly into the fit** (idea
  #43, distinct mechanism from the two_prod entry above -- no added
  ops, a same-op-count refit against `sin(pi*r)` directly instead of
  routing through an intermediate `PI*r` value): seeded from a real
  scipy/HiGHS minimax LP fit (idealized max abs residual ~4.2e-8 vs. an
  f32-emulated reconstruction of the shipped route's own ~1.08e-7, only
  a ~2.6x idealized margin -- a weaker signal than most of this
  session's other LP fits going in), polished by tune.rs's coordinate
  descent (barely moved, confirming the LP seed was already good on
  that grid). Wired into new `sinpi_poly_raw`/`sinpi_poly` functions
  (mirroring `sinf_poly_raw`/`sinf_poly`'s copysign split) and real-
  fuzzed: a clear regression, not the hoped-for win. `sinpi` alone
  (simplest case, no other confound) got worse on *both* axes -- avg
  ulp 0.1969→0.2065, max 2→3. `cospi`'s avg also got worse
  (0.1079→0.1105); its already-huge near-x=-0.5 max-ulp outlier (a
  known "near a true zero of the function, ulp isn't meaningful"
  artifact, not a bug -- x≈-0.5 is exactly `cospi`'s zero) happened to
  shrink from ~3.3M to ~103k on this run, but that's still catastrophic
  either way and not attributable to a real precision fix. Reverted --
  the weak ~2.6x idealized margin didn't survive contact with the real
  f32 construction, the same "isolated fit doesn't predict real
  magnitude" lesson as several other entries in this file.
  `examples/tune.rs`'s `sinpi_poly_c` scaffold kept for reference.
- **sind/cosd dedicated poly, fold DEG_TO_RAD_SMALL directly into the
  fit** (idea #44, same mechanism as #43 above, applied to sind/cosd's
  own `d in [-90,90]` reduction instead of sinpi/cospi's `r in
  [-0.5,0.5]`): screened with the idealized LP margin check *before*
  full implementation, learning from #43 -- this one's margin (~4.9x,
  column-scaled LP to fix a HiGHS conditioning failure from `d^9`'s
  ~3.9e17 raw magnitude) was moderately stronger than #43's weak 2.6x,
  but still not the order-of-magnitude-plus margin that reliably
  predicted a real signal elsewhere this session (cbrt, erfc). Given the
  extra risk signal of it being the *same mechanism* that had just
  failed outright for sinpi, implemented and tested anyway rather than
  guessing from the margin alone (the setup effort was mostly shared
  with #43's already-built infrastructure). Same result: a clear
  regression on every axis. New `sind_poly_raw`/`sind_poly` functions
  (mirroring `sinf_poly_raw`/`sinf_poly`, with `POLY_SAFE_BOUND`
  converted to `d`-units for the equivalent clamp) wired into `sind`/
  `cosd` (`tand` inherits via composition): sind avg ulp
  0.1236→0.1453 (+17.6%), cosd avg ulp 0.0725→**0.1297** (+79%!), cosd
  max 2→3. Reverted. Two same-mechanism failures now (#43, #44) confirm
  this isn't just a margin-strength fluke -- folding an irrational
  scaling constant directly into a degree-9-in-the-unscaled-variable
  poly doesn't survive contact with real f32 rounding for this crate's
  trig functions, regardless of how strong the idealized fit looks.
  `examples/tune.rs`'s `sind_poly_c` scaffold kept for reference.
- **sinc: attribute the max-4 error and screen a division correction**
  (idea #50): attributed first via a direct probe at three real fuzz
  worst-case points (not the near-a-true-zero-of-sin artifact class a
  cruder wide scan hit initially, screened out separately): the
  division itself contributes a real, consistent ~2 ulp of error even
  against a hypothetically *perfect* `sinpi` value, with `sinpi`'s own
  ~1-2 ulp error adding on top -- a genuine mixed contribution, not
  purely division-bound but not negligible either, worth the screen the
  idea asked for. Implemented a compensated division (`rcp = 1.0/denom;
  q = s*rcp; r = fma(-q,denom,s); fma(r,rcp,q)`, reusing the one
  reciprocal instead of a second division) and ran it through mca
  first: real cost (latency +1.9%, throughput +12.3%). But fuzzing it
  found something worse than a cost/benefit tradeoff -- a genuine
  correctness regression: max ulp came back as the NaN-mismatch
  sentinel (`u64::MAX`) at `x` near the denormal floor. Root cause:
  forming `1.0/denom` as its own intermediate can overflow to `+-inf`
  for tiny `x` (`denom = PI*x` underflows further), and `tiny_s *
  inf` is an indeterminate `0*inf` that resolves to `NaN` -- exactly
  the failure class this crate is normally careful to avoid, introduced
  here by *not* going through the original direct `s/denom` division,
  which never forms a standalone reciprocal and has no such edge case.
  Reverted immediately (`git checkout -- src/lib.rs`) -- the compensated
  form isn't a safe drop-in replacement for a plain division when the
  input range includes values near over/underflow, regardless of its
  accuracy benefit at ordinary magnitudes.
- **Backlog ideas #129/#130 dropped as moot**: both were explicit
  follow-ups conditioned on their prerequisite ideas shipping (#129 on
  #43's sinpi/cospi poly fold, #130 on #44's sind/cosd poly fold) --
  both prerequisites were tried and rejected this session (real
  regressions, see the dedicated-poly entries above), so the follow-up
  questions ("does the fold change which op erases -0.0's sign", "can
  the clamp move earlier after the fold") no longer have anything to
  apply to. Not investigated independently of the fold.
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
- **tan/tand/tanpi cross-inline CSE audit** (idea #51, "each pays two
  full reductions differing only in the q offset — asm-diff whether
  LLVM shares the common work"): the premise-check half confirmed via
  real mca, not asm-diff -- `tan`'s throughput (2.532 cyc/elem) lands
  almost exactly at `sin`'s (1.151) plus `cos`'s (1.406) combined
  (2.557), meaning LLVM shares essentially nothing between the two
  reductions when `tan = sin(x)/cos(x)` inlines, confirming the idea's
  premise. But this isn't actually CSE-able in the first place: `sin`'s
  `qb = fma(x, FRAC_1_PI, ROUND_MAGIC)` and `cos`'s `kb = fma(x,
  FRAC_1_PI, -0.5) + ROUND_MAGIC` are two *different* fma instructions
  (different third operand) computing `x*FRAC_1_PI` fused with two
  different adds -- there's no shared *intermediate* value for any
  compiler to find, since the multiply never exists as a standalone SSA
  value in either fma. The idea's own fallback ("hand-share at
  composition level") means literally the same mechanism as the
  directly-adjacent `sincos_checked, shared reduction` entry directly
  above -- a bit-exact hand-shared reduction verified exhaustively but
  regressing wall-clock ~1.4% despite mca predicting a ~19% win, this
  session's *one* documented mca/hardware disagreement. Given `tan`/
  `tand`/`tanpi` would need the exact same class of restructuring
  (deriving one trig function's reduction from the other's, needing
  careful q/k tie-case analysis exactly like the tanpi/tand parity-fusion
  idea found fragile earlier this session) for a mechanism already shown
  to fool mca in the *closest* related case, not pursued further without
  a real quickbench wall-clock result in hand first -- mca alone
  wouldn't be trustworthy evidence either way here.
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
- **tanpi/tand: fuse the two independent parity computations**: rejected
  on inspection, no build needed -- the backlog idea's own premise
  ("`sinpi`'s and `cospi`'s own parity/sign corrections cancel
  algebraically in the ratio") is false as a blanket claim, confirmed by
  direct counterexample. `sinpi(x) = raw_sin(r)*sign_sin(q)`,
  `cospi(x) = raw_cos(rc)*sign_cos(k)` (`q=round(x)`, `k=round(x-0.5)`,
  each sign is `±1` from its own `parity()`) -- dropping *both* signs and
  dividing the raw poly outputs directly reproduces the true `tanpi`
  ratio only when `q` and `k` land on *opposite* parities (roughly half
  the domain, `x` in each unit interval's lower half), and gives the
  exact *negative* of the true ratio when `q`/`k` share the same parity
  (the other half, e.g. `x=0.1`: true `tan(0.1*pi)=0.3249`, raw ratio
  `-0.3249`) -- a real, magnitude-preserving sign flip on ~50% of the
  domain, not a rare edge case or a sign-of-zero cosmetic nit. The
  backlog idea's *actual* proposal (share `q`/`k`'s underlying rounding
  work via one computation instead of two independent
  `round_ties_even` calls, still applying the *correct*, non-constant
  sign relationship) remains structurally sound but would need careful
  case analysis of `k`'s exact relationship to `q` (including
  round-to-even tie boundaries) to get the sign right -- real
  implementation work, not the shortcut this entry was screening for,
  and the payoff (saving one hardware round instruction) looks modest
  next to the risk given this exact function pair's own documented
  history of sign bugs (`sin_checked` range-invariant entry above, the
  ±0 bugs listed in Batch 2 idea #164). Not pursued further.
- **tand direct poly (idea #128)**: shipped for `tanpi` (see lib.rs/git
  log — real fuzz-caught bug fixed along the way: the reflection's
  pole-distance must be computed by subtracting *before* scaling to
  radians, not after — the mathematically-equivalent post-scale form
  loses exactly the precision needed near the pole), but the same
  approach for `tand` regressed a real accuracy loss: reusing `sind`'s
  own `d` (already reduced away from `x`) for the pole-distance
  calculation isn't precise enough (real ~15% relative error, max ulp
  over a million near odd multiples of 90), and fixing that needs a
  whole second, `cosd`-style independent reduction. With that fix in
  place `tand` did work correctly (max ulp 3→12, still good), but real
  mca/quickbench numbers no longer showed a clean win over the simpler
  `sind(x)/cosd(x)` ratio once that extra reduction's cost is included
  — reverted, `tand` still ships as `sind(x)/cosd(x)`.

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
- **Sign-algebra select→xor audit, remaining sites** (idea #16, `erf`'s
  `mulsign(1 − exp2(...), x)` combine and `powf_sign_combine`'s 6-select
  tree): investigated rather than implemented. `erf`'s own site already
  uses `mulsign` (bit-based XOR, not a naive select) for the one sign
  operation in its combine -- there's no remaining select there to
  convert; the actual branch choice (`if xa < 0.28 {a} else {b}`) is a
  domain-branch select, not a sign-algebra one, and out of scope for
  this idea. `powf_sign_combine`'s 6-select tree is a much larger
  target, but every one of its selects encodes a distinct, individually-
  documented IEEE754/C99 special case (`y` odd/even, `x`'s sign bit vs.
  value, `y` infinite, `x==1`, `x==-1 && y` infinite, `y==0`) -- a bit-
  mask consolidation risks silently breaking one of these rare-but-
  load-bearing cases for a speed-only payoff (powf's own accuracy is
  already documented `>=312` max ulp, nowhere near a target this could
  move). Given this exact idea class (mulsign/select reassociation) is
  already 0-for-2 on closely related constructs in this file (acos's
  trailing correction, atan's port of atan_latency's own technique --
  both entries directly above/nearby), and neither remaining site offers
  a clear, low-risk target, not pursued further without a more specific
  mechanism than "audit for a bit-mask form."
- **acos_poly Horner→Estrin**: real latency win, but fma reassociation
  regressed asin max ulp 9→12, acos 4→5 (retuning made it worse, →6).
  acos's accuracy is a protected invariant.
- **asind: fold RAD_TO_DEG into asin_small/asin_poly's own coefficients**
  (idea #123, distinct mechanism from the rejected sinpi/sind folds --
  those changed the poly's *input* reduction domain shape; this instead
  rescales an *output*-side poly's coefficients by a constant, a pure
  linear operation that should be at worst neutral versus a
  post-multiply `asin(x)*RAD_TO_DEG` on paper): measured instead of
  assumed, and the idea's own "90.0/45.0 are exact" framing didn't pay
  off in practice -- real exhaustive fuzz found the fold a wash on
  average (0.3376 vs the naive composite's 0.3387) and *worse* on max
  ulp (16 vs 11). Larger-magnitude degree-space coefficients apparently
  accumulate more absolute rounding per fma step than the post-multiply
  costs, even though the transformation is linear and touches no domain
  shape. Shipped `asind`/`acosd`/`atand`/`atan2d` as the plain composite
  instead (see lib.rs/git log) -- not attempting the fold for the other
  three without new evidence it would fare differently for them.
- **erfc's n/d rational Horner→Estrin**: small theoretical win, measured
  as a wash on speed plus a real accuracy cost (avg +2.7%).
- **erf's near-zero Padé branch refit / tail branch (erf_poly) refit**:
  both — max ulp unchanged, avg moved <0.3%. No headroom. Same
  conclusion reached independently via a degree-shed screen (idea #4):
  erf_poly's own worst x is always in the Padé branch (xa<0.28), never
  erf_poly's own branch, so it isn't even erf's binding constraint --
  and dropping deg 6→5 idealizes to 6.35 ulp-equivalent, already past
  erf_poly's own current ~4-5 ulp range on the branch it does own. Not
  attempted.
- **Same near-zero branch, LP numerator refit**: isolated fit predicted an
  84% avg improvement, but real fuzz found a *regression* (0.317→0.325),
  worst case landing right at the 0.28 branch crossover with erf_poly.
  Check the crossover neighborhood for any poly next to a domain split.
- **erf_unchecked, drop the `xa_bounded` 10.0 clamp** (idea #98, the
  same "unchecked tier drops a redundant guard" mechanism that shipped
  for `sind_unchecked`/`cosd_unchecked`): rejected on inspection before
  any implementation -- `erf`'s own doc comment already documents *why*
  this specific clamp isn't a `sind`/`cosd`-style pure safety net.
  `erf_poly` is a plain degree-6 polynomial with a *positive* leading
  coefficient, so past its fitted domain it doesn't degrade gracefully
  or merely risk non-finite output -- it turns around and diverges
  (`erf_poly(9) ~ -92`, `erf_poly(20) ~ +8698`), which the doc comment
  records as a real, already-fixed bug: unclamped, `erf(50)` came out
  `-1.02e17` instead of the true ~1.0 -- catastrophically wrong sign
  *and* magnitude for a perfectly ordinary-looking input, not the
  "possibly non-finite past the documented domain" contract every other
  `_unchecked` tier in this crate carries. An `erf_unchecked` as the
  idea literally proposes would reintroduce that exact fixed bug as new
  public API. Not implemented; `sind`/`cosd`'s own clamp is a genuinely
  different, provably-no-op-in-domain case and doesn't generalize here.
- **acos_poly unconstrained joint acos+asin objective**: improved joint
  score but regressed acos's own max ulp 4→5.
- **acos_poly degree 6→7**: zero-seeded 8th coefficient converged to
  exactly 0.0 (useless). Re-seeded with a real scipy fit found real
  headroom (acos avg/max 0.496/4→0.490/3) but asin was unmoved and mca
  cost was real (+7-11% both functions) — not worth it at this magnitude.
- **acos max 5 as a real-chain-refit target** (idea #142, re-running
  tune.rs's existing degree-6 coordinate descent against its own already-
  fairly-dense ~106k-point grid, on the theory that "exhausted" was a
  stale/coarse-grid artifact): the descent did move on this grid (avg
  0.14045→0.12982, max unchanged at 4 there), so not literally a zero-
  move local optimum on the grid itself. Wired the "improved"
  coefficients into the real `acos_poly`: a first pass compared only
  100M-sample quick fuzz on each side and looked like a regression (max
  4→5) — but re-checking both sides on the full exhaustive 2^32 sweep
  (the shipped baseline's *own* true max turned out to be 5, not 4;
  quick fuzz had simply never sampled that rare point either way) showed
  the real comparison is avg ulp 0.0650→0.0649 (noise) with max
  unchanged at 5 both ways. Reverted (no reason to carry different
  literals for zero real movement) — this reconfirms rather than
  overturns the original "already exhausted" finding. Caught by this
  session's own "always get a same-precision baseline before comparing"
  discipline: a quick-fuzz-vs-quick-fuzz comparison on a function whose
  worst case is this rare very nearly produced a false regression.
- **acos_poly Df32 leading-term split, pi/2 hi+lo**: dramatic acos avg
  improvement (0.496→0.068) but acos max ulp regressed (4→5) and asin got
  worse on both axes (max 9→12). Retuning recovered asin but acos max
  still regressed (4→6), plus a real mca cost.
- **π-constant hi/lo splits in the other inverse-trig combines** (idea
  #124, a plain additive `+ FRAC_PI_2_LO`/`+ PI_LO` correction after each
  combine's existing subtraction/addition, distinct from the Df32-split
  entry above): tested all four named sites. `atan`'s `FRAC_PI_2 - y`
  fold: no real movement either way (avg 0.0675→0.0673, noise; max
  unchanged at 3 -- the real worst case lives in the untouched `a<1`
  branch) but a real mca cost (latency +4.0 cyc / +6.5% on both atan and
  atan2, throughput +3.2%/+11.9%). `asin`'s big branch (`FRAC_PI_2 -
  sqrt(1-a)*asin_poly(a)`): a real regression, avg ulp 0.0199→0.0227
  (+14%), max 6→7 -- confirms the Df32-split entry's "this specific
  combine shape is fragile near asin" finding via an unrelated
  mechanism. `acos`'s `+PI`: by far the worst result of the four --
  avg ulp 0.0650→**0.2722** (~4.2x worse) on the real fuzz, worst point
  at an ordinary interior x (-0.046), not a degenerate edge case; root
  cause not fully traced (the correction term's own arithmetic looked
  sound on paper) but the measured regression is large and unambiguous.
  `atan2`'s correction (`mulsign(FRAC_PI_2 - hpisignx, y)`): the one
  site that *looked* like a small real win on a single quick-fuzz run
  (avg 0.0683→0.0681, max 4→3) -- caught as sampling noise, not a real
  effect, only by re-running both the modified *and* unmodified versions
  three times each: baseline itself swings between max 3 and max 4 run
  to run on this function's 10M-sample (not exhaustive, 2-arg domain)
  harness, and the modified version's avg sits inside the same noise
  band as baseline's own run-to-run variance. Net: 1 real cost-only
  no-op, 2 real regressions, 1 near-miss that would have been
  misreported as a win without the repeat-run check. All four reverted;
  none of this idea's four sites are worth revisiting with this
  technique.
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
- **atan_latency degree-shed (idea #4)**: screened, not implemented --
  dropping deg 8→7 (in u=r^2) idealized to 1.41 ulp-equivalent against a
  real max ulp of 3, only a 2.1x margin (under this crate's own ~10x
  "automatically safe" bar), and the poly's already "at its LP optimum"
  per the entry above. Too thin to risk; skipped in favor of `ln`/log10
  (see lib.rs/git log), which had genuine headroom.
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
- **tanh division-residual correction on `e/(e+2)`** (idea #42, round-off
  audit found error splits ~evenly poly/combine, targeting the combine
  half): the idea's own caution ("atan's analogous division fix measured
  zero benefit at +181% throughput — mca screen before any fitting
  work") held here too. Compensated-division combine (`rcp=1/denom;
  q=e*rcp; r=fma(-q,denom,e); fma(r,rcp,q)`) mca-screened before any
  accuracy fitting: real, decisive regression -- latency 82.72→94.36 cyc
  (+14.1%), throughput 1.859→2.329 cyc/elem (+25.3%). Reverted
  (`git checkout -- src/lib.rs`) without ever measuring accuracy --
  matches this session's now-repeated finding (atan, atan2, sinc) that
  compensated-division corrections cost real throughput in this crate
  regardless of the accuracy payoff, since the extra reciprocal/residual
  ops land on the same contended ports a plain division already uses
  lightly.
- **atanh via single log1p on |x| + mulsign**: sidesteps the catastrophic
  failure mode of the signed-x version (max ulp only 3→4, not 3→31303),
  but still a real net accuracy loss (avg 0.0313→0.0352, max 3→4) with a
  mixed perf result (throughput -29% better, latency +2.9% worse) — not a
  clean win on either axis. **Superseded (shipped, see lib.rs)**: idea
  #69 found this variant's own worst case lands at `x~0.111`, inside
  `|x|<0.25` -- adding a dedicated small-x poly branch there (this
  variant alone, un-split) fixes exactly the loss this entry measured,
  landing a real win on every axis (avg/max ulp 0.0313/3→0.0037/2,
  throughput -25.5%) at a modest latency cost (+5.9%).
- **Dedicated sinh/cosh kernels via reassociation**: real latency win for
  both (~11-12%) but throughput split by function (sinh better, cosh
  worse) and accuracy split the other way (cosh fine, sinh regressed ~12x
  avg ulp). When two functions share a reassociated intermediate, check
  accuracy for both independently.
- **Public `sinh_cosh` pair function** (idea #33, share one `exp_pos_neg`
  call between both outputs instead of paying its reduction/poly/field
  split twice): implemented as a literal transcription of `sinh`'s and
  `cosh`'s own bodies calling `exp_pos_neg` once; bit-identical to calling
  both separately by construction, confirmed over a real exhaustive
  2^32-pattern sweep (0 mismatches). But the idea's own premise — that
  callers needing both today pay the reduction *twice* — is false: `mca`
  reported `sinh_cosh` and a naive `sinh(x) + cosh(x)` at exactly the same
  latency (60.00 cyc) and throughput (2.219 cyc/elem), and a direct
  `--emit=asm` diff of both mca regions showed not just matching op
  counts but an *identical* opcode sequence, 102/102 instructions in the
  same order — `exp_pos_neg` is `#[inline(always)]` and a pure function
  of `x`, so once `sinh(x)` and `cosh(x)` both inline at the same call
  site, LLVM's GVN/CSE already merges the two identical reduction+poly
  computations for free. Same no-op class as the `koff-free
  unchecked-log fast path` and `atan2's bothzero/hpisignx` entries above
  (compiler already does it) — not shipped, since a new public function
  with zero measured benefit over the existing composition is pure added
  surface area. Reverted (`git checkout --`). Only relevant if a future
  caller needs the two outputs from two *non-inlinable* call sites (e.g.
  across a real function-pointer boundary) where CSE can't reach.
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
- **exp_pos_neg: return halves pre-scaled by 0.5 for plain sinh/cosh
  too** (idea #32, the *plain-multiply* relocation `exp_pos_neg_checked_half`
  already ships -- not the rejected bit-trick variant above): lower risk
  than it first looked, since the plain-multiply form is already proven
  safe in the checked sibling. Implemented and real-tested anyway
  rather than assumed: accuracy unchanged (noise-level, as expected for
  a pure reassociation -- sinh avg 0.0822→0.0821, cosh identical, both
  max 5 unchanged). But mca showed a real tradeoff for *both* functions,
  not the checked_half entry's asymmetric-by-caller split: latency
  improved for both (sinh 56.00→54.00 cyc -3.6%, cosh 55.00→54.00 -1.8%)
  while throughput *regressed* for both (sinh 1.971→2.061 +4.6%, cosh
  1.754→1.778 +1.4%) -- a genuine latency/throughput tradeoff, not a
  clean win on the throughput axis this crate's own `exp_pos_neg` doc
  comment names as the priority. Reverted.
- **exp_pos_neg: p_neg from the reciprocal identity** (idea #115,
  `p_pos*p_neg = e^2-(r*o)^2` so `p_neg = (e^2-(r*o)^2)/p_pos`, an "idle
  divider" alternative to the plain `fma(-r,o,e)` combine): mca-screened
  first, before any accuracy work, per this session's own mca-first
  discipline -- and it failed immediately, not marginally. The identity
  needs `ro=r*o` and `e*e` as separate multiplies (neither is available
  from `p_pos`'s own fma, which never materializes `r*o`), plus the
  division itself, replacing one cheap fma with two multiplies + one fma
  + one division -- real, reproducible cost on both callers: `sinh`
  latency 56.00→70.00 cyc (+25.0%), throughput 1.971→2.286 cyc/elem
  (+16.0%); `sinh_checked` latency 58.06→72.06 cyc (+24.1%), throughput
  2.345→2.721 (+16.0%). Unlike `cbrt`'s own early-starting reciprocal
  (independent of the seed chain, genuinely idle), this division depends
  on `p_pos` already being computed -- removing the parallelism the
  current two independent one-fma combines already have (`p_pos` and
  `p_neg` need only `r`,`o`,`e`, computable simultaneously) rather than
  moving work off a contended port. Reverted (`git checkout --
  src/lib.rs`) without ever reaching the accuracy fuzz -- mca alone was
  decisive.
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
- **erfc saturation-threshold tightening screen** (idea #140): checked
  numerically (scipy bisection on the true `erfc`, not a build) where
  each side actually rounds to its saturated f32 output. Positive `x`
  (where the real bottleneck, max ulp 109, lives) doesn't hit exact
  `0.0` until `x≈10.05` — barely past the crate's own `10.0` clamp, no
  meaningful headroom to shrink from there without cutting off real
  nonzero output. Negative `x` rounds to exact `2.0` much earlier
  (`x≈3.83`), a real ~6.2-wide "wasted" slice of the shared rational's
  fit domain — but that's the *already-smaller*, non-binding side (see
  the negative-side accuracy survey entry above: x<0 avg 0.1397/max 6
  vs x>=0 avg 0.4815/max 109), so narrowing there wouldn't move the
  reported max ulp at all. No implementation attempted — the premise's
  own "if it's meaningfully inside |x|=10" condition is false for the
  side that would need to benefit.
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
- **exp10m1's `0.2` seam** (the sixth seam, missing from idea #7's own
  list): **no headroom, 2026-07-27.** Swept against the real *exhaustive*
  sweep (4.29e9 samples per point). Narrowing hurts (`0.15`: avg
  0.1279→0.1283, max **4→5**, worst x moves to 0.1505). Widening does
  nothing: `0.18`, `0.205`, `0.21`, `0.215` are all identical to the
  shipped `0.2` to 4 decimal places on avg, max, *and* worst-x
  (0.1279/4/-5.7932023e-2); `0.2171`, the largest value keeping `|v| =
  |x*LN_10| < 0.5`, is marginally worse (avg 0.1280). The shipped value
  sits in the middle of a flat plateau.
  Structural reason it differs from `exp2m1`, worth keeping: the entire
  *legal* widening window here is `0.2 -> 0.217`, only +8.5%, because
  `LN_10 ≈ 3.3x LN_2` compresses how much `x`-range fits inside the shared
  Pade's `|v| < 0.5` fit. `exp2m1`'s win came from a 30% widening
  (`0.5 -> 0.65`), a big enough slice of the domain to move a global
  average; 8.5% is not. So "an unaudited seam with a shipped analogue"
  isn't sufficient reason to expect a win — check how wide the legal
  window actually is first, since the win scales with the *fraction of
  inputs that change branch*, not with the seam's existence. Also note
  `exp10m1`'s worst case (x≈-0.0579) sits deep inside the Pade branch,
  nowhere near the seam, so no threshold could touch the max.
- **Seam retunes not yet done** (backlog idea #7, five sub-items): all
  five now checked directly against the real exhaustive sweep (plus mca
  where relevant), same method as the earlier 5-function crossover audit
  (sinh/tanh/expm1/asin/erf). `exp2m1`'s `0.5` **shipped as a real win**
  (see lib.rs/git log): 0.5→0.65, avg ulp 0.0769→0.0766, max unchanged at
  4, zero mca cost (branchless select, same op count). The other four
  found no headroom: `exp_m1_over_x`'s `0.5` -- 0.3/0.4/0.6/0.7 all
  clearly worse (avg ulp up to 0.0729→0.0850, max up to 19); a narrow
  0.52-0.55 window did find a genuinely reproducible tiny avg
  improvement (0.0729→0.0728, max unchanged at 6) but real mca showed it
  isn't free here (throughput 1.798→1.822 cyc/elem, +1.3%, reproducible
  at both 0.52 and 0.55) -- unlike `exp2m1`, "same op count" didn't mean
  zero mca cost, so not adopted (a ~0.14% avg win isn't worth a real
  throughput cost). `sinh_checked`'s `0.5` (`cosh_checked` has no such
  branch at all -- `ep+en` has no cancellation near 0, unlike
  `sinh_checked`'s `ep-en`) -- every alternative tried is clearly worse
  (0.3/0.4: avg 0.0429→0.0439/0.0431, max 5→7; 0.6/0.7: avg
  0.0429→0.0509/0.1019, max 5→**28**/**124**), confirming it already
  matches plain `sinh`'s own already-audited `0.5`. `softplus`/
  `logaddexp`'s `87.0` turned out to be a different *kind* of seam, not
  a genuine crossover -- it's a "the true correction has become
  negligible" cutoff, and the accuracy harness's own `softplus_domain`
  explicitly restricts testing to `|x|<80` because comparing ulp right
  at `87` is a documented "sub-denormal-scale difference reporting as
  millions of ulp" artifact (same class as cospi's near-zero artifact) --
  not actionable via this crate's accuracy methodology at all, no real
  measurement exists to retune against. `asinh`/`acosh`'s `2048` rescale
  threshold: swept 100 through 8000 for both functions against the real
  exhaustive sweep -- avg/max ulp identical (to 4 decimal places) and
  worst-x unchanged at every candidate tried, for both functions. The
  direct/rescaled branches are equally accurate across this entire
  range; the real worst case for both (`asinh` x≈0.0155, `acosh`
  x≈1.031) lives entirely elsewhere, nowhere near this seam. No headroom
  either way -- the exact threshold position simply doesn't matter here.
- **sigmoid one-sided evaluation** (idea #199, `e = exp(+|x|)` so `k`
  never goes negative, `s = 1/(1+e)` selected directly for `x<0` or as
  `1-s` for `x>=0`): edgecheck confirmed every special-value pin still
  passes bit-exact (including `sigmoid(-88)`, matched to the ulp, via the
  same fma-antisymmetry identity the shipped negation-folding already
  uses), so the premise's own correctness held up. But mca showed a
  real, reproducible cost on both axes before fuzzing was even reached:
  latency 57.06→61.11 cyc (+7.1%), throughput 1.283→1.468 cyc/elem
  (+14.4%). Reverted (`git checkout -- src/lib.rs`) without running the
  accuracy fuzz -- removing the negative-`k` side of the exponent-field
  construction doesn't pay for itself; the extra `abs()` plus the
  sign-selected `1-s` combine costs more than the single-field trick
  saves by never seeing negative `k`.

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
- **rcbrt direct seed** (idea #53): built `rcbrt_normal`, a from-scratch
  `x^(-1/3)` core mirroring `cbrt_normal`'s own shape exactly -- a
  negated-exponent seed (`MAGIC - ax/3` instead of `ax/3 + MAGIC`) and,
  provably, the *same* degree-3 correction polynomial: for any seed `y`
  approximating `a^p` with `y/a^p = 1+e`, the correction
  `p(r) = ((1+r)^(-p)-1)/r` only depends on `p`, not on which direction
  the seed approximates, and defining `r = si^3*a - 1` (pure multiplies,
  no division at all, unlike `cbrt_normal`'s own `1/a`-dependent `r`)
  gives exactly the same shape for `p=-1/3` as `cbrt_normal`'s `r =
  (s^3-a)/a` gives for `p=1/3`. Deletes *both* divisions the old
  `1.0/cbrt(x)` composition paid (`cbrt_normal`'s internal `1/a` and the
  outer `1/cbrt(x)`) down to one (an unconditionally-evaluated
  `1.0/(x+x)` for the zero/inf/nan special case the seed's bit trick
  can't reach on its own). Real mca win, confirmed: latency 46.30→31.50
  cyc (-32.0%), throughput 1.666→1.511 cyc/elem (-9.3%).
  But real fuzz caught a severe accuracy regression the mca-first,
  fuzz-second order didn't prevent: avg ulp 8.94, max ulp **62** (vs.
  shipped 0.418/5) at a magic constant chosen only from a `[1,2)`-octave
  proxy plus exactness pins at `x=1,8,-1,-8`. Root cause, confirmed by a
  full-domain per-octave sweep: `ax/3`'s truncating integer division has
  three distinct rounding classes depending on the exponent mod 3, and
  `x=1`/`x=8` both sit in the *same* class (`e=0` and `e=3`, both ≡0 mod
  3) -- the exactness pins and the `[1,2)` proxy calibrate only that one
  class, leaving the other two completely unchecked. One of them (`e≡1
  mod 3`) pushes the residual `r` for *every* octave in that class
  (deterministically, exactly 88.709 ulp every third octave, confirmed
  from `e=-125` to `e=127`) outside `cbrt_normal`'s poly's fitted domain
  (`|r|<=0.0998`) -- an extrapolation blowup, not a rare outlier.
  Re-searched the magic constant properly, this time scoring all three
  classes at once (3 consecutive octaves `[1,8)`, the minimum span that
  samples each class once) across the same ~620k exactness-preserving
  candidates: the *best* available candidate still only reaches max ulp
  ~23.2 across the three classes -- confirming this isn't a search
  failure but a real headroom shortfall in the *reused* poly. Reverted
  (`git checkout -- src/lib.rs`) -- the mathematical shape-reuse insight
  is sound and the division-elimination mechanism is real (confirmed via
  mca), but `cbrt_normal`'s poly was fit against its own forward seed's
  residual distribution, which apparently has enough margin across all
  three alignment classes only for *that* construction; the inverse
  seed's residual distribution needs its own dedicated minimax fit (not
  a reused poly) to cover all three classes simultaneously -- confirming
  the backlog entry's own "real fitting work" framing rather than
  finding a shortcut around it. A future attempt should fit directly
  against 3-octave-sampled data (not a single representative octave) to
  avoid this exact trap from the start.
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
- **Clenshaw/Chebyshev-basis evaluation screen** (idea #108), both named
  targets: killed on paper before implementation, for two different
  reasons per site. `atan_latency`'s own entire reason for existing is
  minimizing critical-path *depth* (it trades away atan_poly's division
  specifically for more, shallower fma parallelism -- see its own doc
  comment) -- its current Estrin split evaluates the degree-17 poly in
  ~5 dependency-chain-deep steps by sharing `r2`/`r4`/`r8` across two
  parallel Horner groups. Clenshaw's recurrence is inherently serial
  (each term needs the previous *two*, `b_k = 2r*b_{k+1} - b_{k+2} +
  a_k`), so evaluating a degree-17 series this way needs ~17 sequential
  steps -- directly destroying the ~5-deep parallelism this specific
  tier is built around, for a function whose whole point is a shallow
  critical path. `erfc_rational`'s n/d is already plain (already-serial)
  Horner, so it wouldn't lose parallelism the same way, but it's only
  degree 4/4 -- Clenshaw's real numerical-conditioning advantage over
  power-basis mainly shows up at much higher degrees, where the power
  basis's coefficient magnitudes span many orders of size; a degree-4
  fit's coefficients don't have that problem (confirmed already-fine
  ranges in this session's own `erfc_c5` degree-5 LP fit). Neither site
  looks like a good match; not implemented.
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
  **Idea #101's fix, tested (2026-07-20)**: built `tune_basin_hop_avg_first`
  -- same annealing mechanism, but every comparison uses a `better()`
  function instead of the default `(max, sum)` tuple ordering: candidates
  that would regress max beyond the un-hopped starting point are always
  rejected, and among candidates that hold the cap, lower avg wins. Fix
  verified mechanically correct on the same acos_poly coarse grid the
  original bug was found on: max held at 4 (matching the cap) with avg
  0.14447→0.12758, vs. the old max-first version's max 2→**avg 1.878**
  (the exact bias, reproduced fresh for comparison). But wiring the
  avg-first result into the real `acos_poly` and exhaustively verifying
  found only noise-level movement (avg 0.0650→0.0647, max unchanged at
  5) -- acos_poly's coordinate-descent optimum (already reconfirmed
  separately this session, see the "acos max 5 as a real-chain-refit
  target" entry above) holds up even against basin-hopping once its
  known bias is fixed; there just isn't more headroom here for *any*
  local-search variant to find. Reverted the `lib.rs` coefficient swap
  (no real benefit); kept `tune_basin_hop_avg_first` in `examples/
  tune.rs` as working, bias-free infrastructure for future basin-hop
  attempts on functions that might actually have exploitable headroom.
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
- **log_family_wrapper: cheaper special-case classify** (idea #37):
  considered, not implemented. `spec = if x==0.0 {-inf} else {NaN}`
  and the outer `if x<=0.0 {spec} else {r}` aren't redundant computation
  of the same condition -- they answer two different questions (what
  value should the special case be, vs. whether to use it at all), so
  there's no obvious single-compare consolidation the way the
  `wrapping_sub` entry above found for a genuine OR-of-two-comparisons-
  on-one-variable pattern. Declined given the accumulated evidence this
  session that (a) LLVM's InstCombine already canonicalizes the classes
  of compound-comparison pattern that genuinely are redundant (the
  `wrapping_sub` entry, `koff-free unchecked-log`, `atan2`'s bothzero/
  hpisignx, all confirmed byte-for-byte no-ops), and (b) hand-deriving a
  *new* bit-domain reformulation here carries real risk of getting it
  subtly wrong (see this session's own near-misses: the falsified naive
  attempt at idea #12, and idea #50's reciprocal-overflow NaN bug) for
  an idea explicitly labeled "audit-grade micro" -- low expected payoff
  even if executed correctly.
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
- **NaN-propagating min/max audit** (idea #196, "the IEEE maxNum
  NaN-discard trap has now bitten softplus, hypot_checked, and
  logaddexp separately"): investigated rather than built. Grepped every
  `.max(`/`.min(` call in `src/lib.rs` (7 total) and checked each one's
  actual NaN safety: `softplus`/`logaddexp`'s four sites (`x.max(0.0)`,
  `a.max(b)`, two `.min(87.0)` clamps) are all covered by each
  function's own trailing `is_nan()` override, which makes any
  intermediate NaN-discard moot regardless of what it computes; `atan`/
  `atan_latency`'s `a.min(1.0/a)` are safe by construction, not by a
  guard -- a NaN `a` makes `1.0/a` NaN too, so `min` sees NaN on *both*
  sides and correctly returns NaN (no non-NaN operand to wrongly
  prefer); `hypot_checked`'s `ax.max(ay)` is the already-fixed instance
  idea #196 itself cites, with its own doc comment explaining exactly
  why the narrower `is_zero` check downstream doesn't reintroduce the
  discard. Zero unguarded/unsafe sites found -- no 4th bug of this class
  currently exists to fix. The idea's other half (a shared
  `max_nan_prop` helper) wouldn't fix anything live, only reduce
  duplication across three different already-correct mechanisms; not
  pursued given nothing is actually broken.

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
5. **Per-caller Pade refit where domains genuinely differ**: exp2m1's
   Pade argument is y=x·ln2, |y|<0.347 — narrower than expm1/tanh's
   shared |y|<0.5 fit. (The rejected decouplings — sinf_poly,
   exp_pos_neg — had *identical* caller domains; this one doesn't.)
   **Rejected on a pre-implementation screen, 2026-07-27, on two
   independent grounds.** (a) *The premise is stale*: the `|y|<0.347`
   figure came from exp2m1's old `|x|<0.5` threshold, but that seam was
   since retuned to `|x|<0.65` (see exp2m1's own doc comment), so the
   real argument range is `|y| < 0.65*ln2 ≈ 0.4505` — a 10% narrowing
   vs. the shared 0.5, not the 30% the idea assumes. Same
   re-stale-check lesson as the asin_small/asin_poly entry below, in the
   other direction: a *proposal's* premise can go stale from a later
   change to the function it's about, not just a rejection's. (b) *The
   binding error isn't at the domain edge, so narrowing cannot help*:
   evaluating the current coefficients in exact arithmetic, max relative
   error over `|y|<0.4505` is 0.617 ulp-equivalent and it peaks at an
   **interior** point, `y≈+0.334`. That peak sits inside every candidate
   domain, so `|y|<0.4505`, the stale `|y|<0.347`, and exp10m1's own
   `|y|<0.4605` all score *identically* (0.617) — the shared 0.5 domain's
   larger 1.367 figure is reached only at its edge, which exp2m1 never
   visits. There is nothing for a narrower fit to recover. (c) Confirming
   this from the other side: exp2m1's real exhaustive max is 4 ulp at
   x≈-0.3991 (`y≈-0.277`), where the idealized fit error is only 0.158
   ulp — i.e. ~96% of the observed error there is rounding chain, the
   "near precision floor" diagnosis that has predicted refit failure
   every time in this crate (see the LP-margin entries above). Note idea
   #30 already rejected the *other* per-caller mechanism here (folding
   `LN_2` into the rational) at a 5.6x margin; a HI/LO split of `LN_2`
   instead (the sind idea #28 mechanism, which *adds* precision rather
   than removing a rounding) is the one variant still untried, but #28's
   own outcome — real accuracy win, real mca cost — is what to expect.
6. **Joint threshold+coefficient coordinate descent** in tune.rs
   (crossover as a continuous search parameter) — automates the asin
   fix-5 lesson instead of retuning thresholds against frozen polys.
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
    **Prior against this got substantially stronger 2026-07-27**: the
    closely-related `to_int_unchecked` experiment (see #14 and the "k1 from
    bit-twiddled k" rejected entry) built a fully-vectorized pure-integer
    exponent-field construction for `exp2_checked` — packed
    `vcvttps2dq`/`vpsrad`/`vpaddd`/`vpslld`, codegen_check clean — and it
    **regressed throughput 13.2%**. The measured reason applies to this
    idea too: an integer path must add the `127`/`383` exponent bias per
    field *explicitly*, whereas the magic-round bakes it into the constant
    so `<< 8 & MASK` extracts a ready-biased field for free. The hoped-for
    "less-contended ports" saving also failed to appear, because the
    float->int conversion competes for the same ports as the FP ops it
    displaces. This doesn't formally falsify the `<< 23` construction
    (still never implemented), but the port-pressure premise it rests on
    has now been measured and did not hold.
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
    the clamp bound at codegen time). **Not needed for that target:**
    #14's `to_int_unchecked` solves the codegen half outright (see below),
    and the target then loses on op count anyway. Still untried as a
    general technique on other clamped values, but note the motivating
    case is now closed.
14. **`f32::to_int_unchecked` where the range is guaranteed**: lowers to
    plain vcvttps2dq (vectorizes), unlike the saturating `as i32` that
    de-vectorized exp2_checked's rejected variant. codegen_check
    mandatory, pairs with #13. **Tested 2026-07-27 on that exact target.
    The mechanism works; the optimization it was meant to rescue does
    not.** Codegen came out exactly as claimed — 2 packed `vcvttps2dq`,
    zero scalar converts, packed `vpsrad`/`vpsubd`/`vpaddd`/`vpslld`,
    codegen_check clean across all 148 regions — but `exp2_checked`
    throughput still regressed 1.399 -> 1.584 (+13.2%), because the
    integer path needs 9 ops against the magic-round's 8 (it must add the
    `127` exponent bias per field explicitly, where the float trick bakes
    it into the magic constant). Full numbers and the NaN-UB caveat
    (`to_int_unchecked` on NaN is UB and `.clamp()` does *not* strip NaN
    — quash with `max`/`min`, which return the non-NaN operand) are in
    the "k1 from bit-twiddled k" entry above. So: keep this in the
    toolbox as *the* way to get a vectorizing float->int cast, but it
    doesn't make a pure-integer exponent-field construction competitive
    with the magic-round idiom.
15. **lzcnt-based denormal normalization** (`u32::leading_zeros` →
    vplzcntd, AVX-512CD): replace the compare+select 2^24 rescale in the
    log family/cbrt/hypot_checked with an exact shift-based normalize.
    **codegen screen done** (2026-07-20): a standalone `--emit=asm` probe
    (`u32::leading_zeros()` over a `[u32;16]` array loop) confirms this
    target really does lower `leading_zeros` to packed `vplzcntd`
    (`ymm`, two lanes of 8), not a scalar fallback -- the risk the idea
    itself flagged doesn't materialize, so this is a live option, not a
    dead end. Full implementation (an exact per-lane shift + exponent-
    field reconstruction to replace the current `denormal_rescale!`
    macro's plain "×2^24 unconditionally, if denormal" multiply) not
    attempted yet: unlike the multiply-based rescale, which is already
    just 1 compare + 1 multiply + 2 selects and needs no per-input shift
    amount (the ×2^24 constant works uniformly across the whole denormal
    range), an exact-shift version needs strictly more ops (lzcnt +
    shift + OR to reconstruct + still a compare to gate it) -- any win
    would have to come from moving work off a contended FP-multiply port
    onto less-contended integer ALU ports (the idea #12 theme), not from
    a lower op count, so it needs a real per-caller mca measurement
    before it's clear this is even a net win in principle, let alone
    worth the implementation risk (this exact code path has a documented
    bug history, see the "log_2 denormal path" rejected entry).
22. **Toolchain-bump re-screen list**: tag the rejections that were pure
    scheduling artifacts (pre_offset dead-add removal, ln/log10
    trailing-fma fuse +1cyc, reduce_pi depth-2 rebalance) and re-measure
    after each nightly bump — these can silently flip.

#### exp / log family

23. **exp/exp_checked floor-domain reduction**: superseded by the
    simpler idea #112 mechanism, which shipped instead (see lib.rs/git
    log, `exp_narrow`) -- rather than switching to floor + refitting the
    poly over a doubled/shifted `[0,ln2)` residual domain (this idea's
    own proposal, screened via scipy: a same-shape degree-5 refit there
    reaches only ~4.9 ulp-equivalent idealized error vs. the current
    poly's ~2.0 on its own domain, i.e. real headroom loss unless bumped
    to degree 6, more risk for the same destination), #112 just narrows
    `exp`'s *existing* round-based domain a hair (to the exact point
    where `k` still never reaches the split-requiring edge) with the
    *same* poly, unrefit. Not pursued further given #112 reaches the
    same single-field destination with strictly less risk.
36. **rlibm-style discrete rounding-interval LP extended to ln/log10**
    (same 2^23 reduced-input multiplicity as the existing log_2 entry).
38. **softplus fused kernel**: one fitted poly for ln(1+2^-t) over the
    k/f-reduced domain, replacing exp (full poly) → log1p (division +
    deg-9 ln poly). Big throughput candidate.
39. **logaddexp: same fused kernel** on |a−b|.
41. **logaddexp2** (base-2 sibling, ML/audio) — near-free variant of
    whatever #39 lands on.

#### sin / cos family

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
    **Implemented and measured 2026-07-27; rejected as a fast-tier
    replacement, but the numbers make it a live "intermediate tier"
    candidate — see below.** Real accuracy win, real perf cost, so it
    fails the accuracy-without-a-penalty bar:
    - Accuracy (real quick fuzz, both sides same harness run): `sin`
      `|x|<=1e6` avg/max **0.0409/7 -> 0.0356/2**, `|x|<=1000`
      0.0212/2 -> 0.0198/2, in-domain 0.0646/222 -> 0.0589/220. The
      `|x|<=1e6` numbers land *exactly* on `sin_checked`'s own
      (0.0356/2) — the fast tier becomes as accurate as the double-float
      tier over that range.
    - It does **not** move the ~1.3e7 cliff, and can't: the cliff is the
      magic-round producing a wrong `q` (`|x/pi| >= 2^22`), not an
      inexact `q*pi` product. The idea's "push the cliff far out" premise
      conflates the two error sources; in-domain max stays ~220 because
      the post-cliff region dominates it either way. Moving the cliff
      needs more pi bits (Payne-Hanek) or checked's double-float q.
    - Cost (mca): `sin` throughput 1.151 -> **1.278 (+11.0%)**, latency
      46.00 -> 50.00; `cos` throughput 1.406 -> **1.651 (+17.4%)**,
      latency 54.00 -> 58.00. "Similar op count" doesn't hold: the
      working form is **6 ops** (mul + fma for two_prod, two subs, then
      `PI_LO` *and* `PI_TINY` fmas) against the current 4 fmas.
    - `PI_TINY` is **not** droppable, which is what kills the op count.
      The 5-op two-word version (`PI_HI`+`PI_LO` only) is catastrophically
      worse than even the baseline: `|x|<=1e6` max **8758**, in-domain max
      **126119**. Reason the naive estimate misses this: `q*PI_TINY` is
      ~1.1e-9 at `q ~ 3.2e5`, which looks like ~0.01 ulp *if you assume
      `r ~ 1`* — but near sin's zeros `r` is itself ~1e-9, so that
      absolute term is the whole answer. Any future "drop a pi word" idea
      needs to be scored near the zeros, not at generic `r`.
    - Rejected *as a replacement* because sin/cos exist to be the cheap
      tier (their own doc comments point accuracy-seeking callers at
      sin_checked), and throughput is the axis this crate optimizes —
      paying 11-17% there to buy accuracy the checked tier already sells
      is the wrong direction for that function.
    - **But**: at 1.278 cyc/elem it delivers `sin_checked`'s exact
      `|x|<=1e6` accuracy (0.0356/2) for **~4.2x less** than
      `sin_checked`'s own 5.311 — which is precisely the standing
      "Intermediate sin/cos tier" backlog entry's value proposition, now
      with measured numbers instead of a guess. Cheap to revisit as a
      *new* `sin_mid`/`cos_mid` pair (additive, zero regression risk to
      existing callers); the open question is whether the narrow win
      (max 7 -> 2 over `|x|<=1e6`, same cliff, same domain as fast sin)
      justifies the API surface.
    - **Methodology warning worth reusing**: an isolated screen of just
      the reduction — scoring `r`'s error in units of `ulp(r_true)` with
      `q` held fixed — reported *no benefit at all* (CW 123 vs two_prod
      120 at 2^19..2^20, and two_prod slightly *worse* at 2^20..2^21) and
      would have killed this idea outright. That metric is an artifact:
      it explodes wherever `r_true` lands near a zero of sin, so it
      measures proximity-to-zero, not reduction quality. The real
      end-to-end fuzz disagreed. Same artifact class as the
      identity-fuzz entry — score the shipped function, not an
      intermediate, whenever cancellation is in play.
49. **sinf_poly real-chain refit** (#1's method) scoring sin_checked +
    cos_checked's actual reductions jointly — the rejected LPs used
    continuous grids that mis-weighted the caller split.

#### cbrt / sqrt / hypot

56. **Slice-tier FTZ/DAZ via MXCSR**: a slice entry point can set
    FTZ/DAZ around its own loop and restore — gets the FTZ
    feature-flag idea's win without a global cargo feature.

#### hyperbolics / activations


#### powf / pown / remainder

73. **powf_mid tier**: 1.5-word log2 (k + one two_prod hi/lo pair, no
    full Df32 arithmetic) through the y multiply + a single
    multiplicative correction — targets powf's y-amplified error at a
    fraction of powf_checked's +61% throughput cost.
76. **Constant-base powf slice**: precompute log2_df(x) once per slice;
    per-element work drops to a Df32 multiply + exp2 (slice-tier
    candidate, fixed-base workloads).
82. **Const-y remainder/fmod slice variant**: precompute 1/y,
    `q = round(x * (1/y))` — trades the per-element division for a
    multiply; the changed q rounding needs an accuracy screen
    (slice-tier candidate).

#### New functions / API breadth

87. **Public Df32 module** (log2_df/exp2_checked_df/two_prod etc.) for
    power users composing their own accurate kernels. Partially done:
    the EFT primitives themselves (two_prod/two_sum/quick_two_sum) went
    public via #184; what remains is the Df32 *type* and the df-suffixed
    kernels (log2_df/exp2_checked_df), which is the larger API decision.

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
99. **tgamma** companion to the lgamma entry (Lanczos/Stirling, shares
    machinery).
100. **Bessel j0/j1** (Cephes-style two-region rational + trig
     composition — big job, listed for completeness like lgamma).

#### Batch 2 (same session): fitting & search, continued

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
110. **±few-ulp exhaustive scan of every non-poly literal** (clamp
     bounds, seed constants, magic offsets, branch thresholds) scored on
     the real fuzz — #3's sibling for non-coefficient constants.

#### Batch 2: exp / log family

111. **expm1_checked / exp_m1_over_x_checked**: input clamp + single
     exponent field (the tanh/sigmoid pattern) — likely cheaper than the
     current k1/k2 split *and* total-domain, since the clamp caps k at
     127 by construction. **`expm1_checked` shipped 2026-07-27** (see
     lib.rs/git log); `exp_m1_over_x_checked` deliberately deferred, see
     below.
     - The idea's own premise needed one correction: "the clamp caps k at
       127 by construction" is exactly what must *not* happen. Capping k
       at 127 caps the output near `2.4e38`, which (a) returns short
       across the top third of an octave, `x` in `[88.376, 88.723)`,
       whose true results are finite and representable right up to
       `f32::MAX`, and (b) **saturates finite instead of overflowing to
       `inf`** above `ln(f32::MAX)` — the same failure `edgecheck.rs`
       caught in the rejected round-based `exp10_checked` reduction, and
       the reason `exp_checked` itself still carries the k1/k2 split
       (its clamp is deliberately set at `128/log2(e)` so `k` *reaches*
       128 and overflows on its own).
     - Fix that makes the idea work: emit the field at **`k-1`** (offset
       `382`, not `expm1_narrow`'s `383`) and fold the missing factor of
       two into an exact `p + p`. `k` can then reach 128 while `k-1`
       stays inside a single field's `[-126, 127]`, so the top of the
       range overflows naturally with no extra select. Clamp is
       `[-86.0, 128/log2(e)]`; the bottom bound only has to keep
       `k-1 >= -126` somewhere the answer is already exactly `-1` (true
       for any `x < -17.4`, since `e^x` is then under half an ulp of 1),
       and `-86.0` gives `k = -124`, clearing it with room. Verified over
       all 2^32 inputs that `k-1` lands in `[-125, 127]`.
     - Result: **bit-identical to `expm1` everywhere `expm1` is valid**
       (zero disagreements, exhaustive), so max ulp is `expm1`'s own 6 at
       the same worst point — totality at no accuracy cost. mca is a real
       split: throughput 1.695 -> **1.595** cyc/elem (-5.9%), latency
       71.00 -> **75.06** cyc (+5.7%). Fewer total ops buys the
       throughput; the clamp sits at the head of the dependency chain
       while the `exp2_field_split` work it displaces was partly parallel
       to it, so the path lengthens as the op count falls — the
       op-count-vs-path-depth distinction in mirror image (cf. the ln/
       log10 degree-shed entry, where it was throughput-win/latency-wash).
       Adopted on throughput, the axis this crate optimizes, and note the
       comparison scale: `exp -> exp_checked` pays **+30%** throughput for
       the same totality, where this is negative.
     - `exp_m1_over_x_checked` deferred, with a real obstacle rather than
       for lack of trying: `(e^x-1)/x` has to divide by something, and no
       single choice of divisor is correct at both ends. Dividing by the
       *original* `x` is right at `-inf` (`-1/-inf = +0`) and right for
       large finite `x` (`inf/1e10 = inf`), but gives `inf/inf = NaN` at
       exactly `+inf`; dividing by the *clamped* `x` fixes `+inf` but
       then `-inf` returns `-1/-86 = 0.0116` instead of `0`. Correcting
       the surviving `+inf` case needs an explicit compare+select, which
       would eat much of a saving that is only ~0.39 cyc/elem gross
       (`exp_m1_over_x` 1.798 vs `exp_m1_over_x_narrow` 1.411) before the
       clamp is even paid for.
       - `x.min(HI)` as the divisor (1 op, not a 2-op compare+select) does
         fix both ends of that particular problem — `+inf` numerator over
         a finite `HI` gives `inf`, and a `-1` numerator over the original
         very-negative `x` gives `+0`. But it doesn't rescue the idea,
         because of a **deeper obstacle found while checking it
         (2026-07-27)**: `(e^x-1)/x` is mathematically finite well past
         where `e^x` itself overflows. `e^x` overflows f32 at
         `x ≈ 88.7228`, but the *quotient* stays representable to
         `x ≈ 93.2582` (at `x=93` the true value is 2.64e38, comfortably
         in range). Any implementation that forms the numerator first
         returns `inf` across that whole ~4.5-wide band. Verified on the
         shipped function: `exp_m1_over_x(89.0)` gives `inf` where the
         true value is 5.0445088e36, and likewise at 88.73/90/92/93.
       - That is *within* the current contract — `exp_m1_over_x`'s doc says
         "garbage outside roughly `x in [-87.3, 88.7)`" — so it is a
         documented domain limit, not a bug. But it means a `_checked`
         tier here cannot be built by clamping alone: clamping saturates
         to `inf` exactly where honest finite answers exist, which is the
         premature-overflow defect class, not a fix for it.
       - Mechanism that would actually work, for whoever picks this up:
         fold the division into the exponent for the large-`x` arm, i.e.
         `e^x/x = e^(x - ln x)`, the same "reorder so the intermediate
         never overflows" lever as the `R(v)/(2*x)` -> halve-first fix in
         the dawson entry. Cost is the problem: a branchless select
         computes both arms, so this buys a full `ln` on *every* call to
         extend a band most callers never touch. Screen the mca cost
         before building it.
114. **FTZ-mode minimal exp2_checked/exp_checked** (rides the MXCSR
     slice-tier idea #56): lower clamp −151→−126 and the
     denormal-rounding half of the split's job disappears; same cascade
     deletes denormal_rescale from the log family and cbrt — scope #56
     to capture all of it.
#### Batch 2: trig

125. **Integer-domain parity pipeline end-to-end** for
     sin_checked/cos_checked (parities as bits, XOR combine, direct
     sign mask) — composes #45/#46; deletes the float compare+select
     flip.
#### Batch 2: roots / hypot / geometry

136. **normalize2/normalize3 *slice* kernels** (rhypot + scales — the
     operation users actually want hypot for). Single-call
     `normalize2`/`normalize3`/`normalize4`/`hypot4`/`rnorm4` shipped
     (see lib.rs/git log, ideas #55/#134) -- the slice-batched variant
     itself still needs the slice-tier infrastructure that doesn't
     exist yet (see the standing Batch/slice API tier entry above).
#### Batch 2: hyperbolics / ML / graphics

148. **Softmax / logsumexp / normalize slice reductions** (max-pass +
     exp-pass + sum + scale in one fused traversal) — slice-tier
     flagship, plus a rotate2d(sincos) demo kernel.
#### Batch 2: powf / pown / remainder

154. **Const-arg slice family generalization**: const-y remainder
     (#82), const-base powf (#76), const-base log (precompute
     1/log2(b)) — one shared design decision.
155. **powf slice integer-y dispatch**: slice checks all-y-integral once
     and routes to the pown path — per-slice dispatch is free where
     per-call isn't.
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
     **Shipped 2026-07-27** as `examples/special_matrix.rs`, covering all
     108 public 1-arg `f32 -> f32` functions, and it found a sixth ±0 bug
     on the first run — making it 6-for-6 on this bug class.
     - **Bug found and fixed: `wrap_pi(-0.0)` returned `+0.0`.** Same
       root cause as `sinpi`'s: inside `reduce_pi_checked` the residual is
       formed by subtracting equal signed zeros, which IEEE754 resolves to
       `+0.0`, so `r` arrives with the sign already erased and nothing
       downstream can recover it. Fixed with the identical `if x == 0.0
       { x }` guard `sinpi` already carries.
     - Isolating it mattered more than spotting it. `wrap_pi` is *not*
       bit-exactly the identity on all of `(-pi, pi]` — beyond `|x| >
       pi/2`, `q` flips to `±1` and the `r±pi` branch rounds, giving 7.79M
       one-ulp differences out of 2.16e9, all legitimate. Restricting to
       `|x| <= pi/2` (where `q = 0`, so `r == x` exactly) leaves exactly
       **3** violations in 2.14e9 inputs: `-0.0`, plus `±pi/2` itself,
       which are the same 1-ulp branch-boundary effect (f32's `pi/2`
       rounds *above* true `pi/2`, so `x/pi > 0.5` and `q` becomes 1).
       So the real signal was 1 anomaly in 2.14e9, not 7.79M.
     - Cost of the fix: throughput 4.103 -> **4.280 (+4.3%)**, latency
       92.00 -> 92.02 (flat). Accepted as a correctness fix rather than
       weighed as an optimization — `sinpi` sets the precedent of paying
       exactly this for exactly this bug. If the cost is ever unwanted,
       the established pattern is a `wrap_pi_unchecked` sibling excluding
       `-0.0` from its domain, mirroring `sinpi_unchecked` (idea #98);
       not built unprompted.
     - Three findings the matrix surfaced that are **not** bugs, confirmed
       against their own docs: `sinpi_unchecked(-0.0) = +0.0` is idea
       #98's deliberate documented tradeoff; `erfcx(+inf) = 5.61e-2` is
       the known frozen tail (`erfcx_checked(+inf) = +0` is correct, see
       the erfcx entry); and all 8 NaN-propagation failures are
       `_unchecked`/`_approx` tiers that promise nothing off-domain
       (`exp2_approx`, `ln_unchecked`, `log10_unchecked`, `log2_approx`,
       `log_2_unchecked`, `rcp_approx`, `rsqrt_approx`, `sqrt_approx`).
     - Gate design: NaN propagation and NaN quietness are checked
       automatically (nothing off-domain-specific needed), with the 8
       exemptions listed **by name** rather than matched on an
       `_unchecked`/`_approx` name pattern, so a *new* function that
       silently inherits garbage NaN behaviour still fails. The ±0/±inf
       columns are printed for review, since the correct value there is
       function-specific. Also prints the f(+0)/f(-0)-differ list (59
       functions) as a standing record of which functions are
       sign-of-zero-preserving.
     - **2-arg half also shipped** as `examples/special_matrix2.rs`: the
       full ±0/±inf/NaN cross product (9x9 = 81 combos per function) over
       the 2-arg surface, which is where the original `atan2` and
       `remainder` ±0 bugs lived. Result: **no new bugs.** All nine
       functions with an f64 std counterpart — `atan2`, `hypot`,
       `hypot_checked`, `powf`, `powf_checked`, `fmod`, `fmod_checked`,
       `div_euclid`, `rem_euclid` — match `f64` std *bitwise* on all 81
       combos each (729 comparisons, zero mismatches). Worth stating as a
       positive: the 2-arg special-value surface is in good shape, and it
       is now gated rather than assumed. The comparison is exact rather
       than tolerance-based because std implements the IEEE754/C99 Annex F
       rules and every correct result at these inputs is exactly
       representable in f32.
     - Everything notable in the print-only set (no std counterpart:
       `remainder*`, `rhypot`, `logaddexp`, `xlogy`, `xlog1py`,
       `compound`, `signed_pow`, `mulsign`) traced to a documented,
       deliberate convention, checked against each doc comment rather
       than assumed:
       - `compound`'s NaN at the `n=±0` corners (`compound(±inf, ±0)`,
         `compound(-2, ±0)`, `compound(NaN, ±0)` are all NaN where C99
         `pow(x, ±0)` would be 1.0) is called out verbatim in its own doc
         as "a real, deliberate deviation for a thin composite" — it is
         `exp_checked(n * log1p(x))`, so `0 * ±inf` is NaN by
         construction.
       - `xlogy(±0, y) = +0` and `xlog1py(±0, y) = +0` for *every* `y`
         including NaN and negative `y`: the documented `x == 0` override.
         Note it forces *positive* zero (`{ 0.0 }`, not `{ x }`), so the
         sign of a `-0.0` first argument is not preserved. Left alone
         deliberately — unlike `wrap_pi`, there is no natural sign to
         preserve here (`x*ln(y)` at `x=+0` is `-0` whenever `y<1`, so the
         "faithful" sign depends on `y`), the convention exists precisely
         to override the `0*±inf` case, and scipy likewise returns `+0`.
       - `signed_pow`'s corners are all internally consistent with
         `sign(x)*|x|^y` composed with IEEE `pow`'s own rules
         (`pow(1, NaN) = 1` giving `signed_pow(-1, NaN) = -1`,
         `pow(x, ±0) = 1` giving `signed_pow(NaN, ±0) = 1`).
       - `rhypot`'s `+0` at `(±inf, NaN)` follows IEEE `hypot(±inf, NaN)
         = +inf` (NaN deliberately not propagated), reciprocated.
     - Still open: the `_unchecked` 2-arg tiers were not asserted (they
       promise nothing off-domain, same reasoning as the 1-arg
       exemptions), and no 3-arg surface exists to sweep.
165. **Saturation-boundary pins**: every clamp constant and overflow
     threshold gets an edgecheck pin at ±1 ulp around it — the
     exp10_checked overflow-at-the-boundary pattern, systematized.
     **Shipped 2026-07-27** as `examples/saturation_pins.rs`. Result: no
     new bugs — every clamp is correctly placed.
     - Sweeps ±64 ulp (not ±1) around each bound, because a misplaced
       clamp is usually off by a few ulp, so the failure sits just inside
       or just outside the constant rather than exactly at it. 10
       functions x their bounds = ~3.4k points, each compared *bitwise*
       against an f64 reference, plus 20 direct far-outside saturation
       pins (`f(±1e30)` must equal the mathematical limit, which is the
       literal exp10_checked failure mode this generalizes).
     - Clean: `exp2_checked`, `exp10_checked`, `exp_checked`,
       `expm1_checked`, `exp2m1`, `exp10m1`, `tanh`, `sinh_checked`,
       `cosh_checked` all worst ≤1 ulp in-window, and all 20 saturation
       limits exact (`inf`/`0`/`-1`/`±1` as appropriate).
     - One hit, and it is a *documented accepted gap*, not a new bug:
       `sigmoid` shows 2098176 ulp at `x ≈ -88.72235`. Its own doc comment
       states it verbatim — "for `x` in roughly `(-104.7,-88.7)` the true
       answer is a nonzero denormal but this returns exactly `0`, slightly
       early saturation on a sliver of denormal-scale outputs". Exempted
       in the harness by *true-result-is-denormal*, deliberately not by an
       `x` range, so a new failure at normal output magnitudes still fails
       the gate (129 points skipped on that basis).
     - Useful framing the huge ulp number illustrates: at denormal output
       magnitudes ulp error is ~meaningless as a severity signal (ulp
       there is 1.4e-45, so a wrong answer of the same order as the value
       reads as ~2e6 ulp). Any future gate over near-zero outputs wants an
       absolute or relative-to-magnitude criterion, not raw ulp.
     - Overlaps idea #166 (denormal-output correctness audit): this pass
       incidentally establishes that `sigmoid` is the only clamped
       exp-family function that saturates early into the denormal range;
       the others reach their true denormal outputs.
166. **Denormal-output correctness audit**: which functions produce
     correctly-rounded denormal outputs vs garbage (exp2_checked
     documents its behavior; most others are unaudited).
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
     sanity check. Still untried as such, but the *motivation* (cheap
     coverage for every function) prompted an audit of which 1-arg
     functions accuracy.rs actually measures. **Two real gaps found and
     closed 2026-07-27**: `sqrt1pm1` and `wrap_pi` both had edgecheck
     pins but had never had an ulp sweep.
     - `sqrt1pm1`: avg 0.2053, **max 2** (exhaustive, 3.20e9 in-domain
       samples), clean. Reference is `x / (sqrt(1+x) + 1)`, the
       algebraically-equal form that does *not* repeat the cancellation
       the function exists to avoid — same rule as log1pmx's reference.
     - `wrap_pi`: avg **0.0244**, max 41, worst x 8953.539, over
       `|x| <= 1e4`. The max is the crate's usual near-a-true-zero
       artifact, not a defect: the worst inputs sit almost exactly on a
       multiple of `2*pi`, so the answer (`-2.3e-7`) is near-total
       cancellation of operands ~1e4. Away from those points it measures
       **<= 0.41 ulp**.
     - Methodology point worth keeping: the first version of this
       reference used a single-word f64 `TAU` and reported max **66**. The
       reference was itself the inaccurate party — a 1-word vs 2-word f64
       reduction disagree by **24.56 ulp** at that same x. Splitting TAU
       into `TAU_HI + TAU_LO` dropped the reported max to 41, matching an
       independent scalar probe exactly. When referencing a function that
       rides a double-float reduction, a single f64 word is not
       automatically good enough ground truth — check the reference
       against a wider one before believing a large max.
     - The `_approx` tier stays out of accuracy.rs (it is explicitly
       outside the 0.5/2 ulp budget, so ulp is the wrong metric), but idea
       #188's doc-comment bounds for it were **unverified** — now covered
       by `examples/approx_bounds.rs`, see #188.
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
     reinventing these wrong. **Shipped** (see lib.rs/git log): all four
     are now `pub` with full contract docs, and the emitted assembly is
     bit-identical to the private version (verified by a full pre/post
     `mca_target.s` diff — these are `#[inline(always)]`, so going public
     only adds standalone rlib symbols and cannot perturb inlined call
     sites). Zero perf/accuracy risk by construction; no mca or fuzz
     delta to weigh. The value was in the *contract audit* going public
     forced, which found a real doc bug: `two_prod`'s original internal
     note read "(no overflow)" as a *precondition*, and the first
     public-facing rewrite reinterpreted it as a *guarantee* ("exact for
     any a, b"). It isn't — `two_prod` is exact only for
     `2^-102 <= |a*b| <= f32::MAX`, and both ends fail hard rather than
     degrading (overflow: `p=inf, e=-inf`, so `p+e` is NaN; underflow:
     `e`'s bits reach to `2^(E-47)` vs. subnormals' `2^-149` floor, so
     `e` truncates, and below `2^-149` both flush to zero on a nonzero
     true product). Note the bound is `2^-102`, **not** the `2^-103` the
     "product must be normal" rule of thumb gives — `2^-103` still admits
     real failures (~21k in 42M in-range pairs), a trap worth remembering
     for any future EFT work. `two_sum`/`quick_two_sum`/`mulsign` needed
     no contract changes beyond an overflow caveat: `two_sum` is exact
     for any finite non-overflowing pair including throughout the
     subnormal band, and `mulsign` matches xor-of-sign-bits exhaustively.
     New standing verification: `examples/eft_contract_check.rs` (exits
     nonzero on violation; 40M mixed random pairs + exhaustive mulsign
     over 8 y-values x 2^32). Another instance of the internal->public
     promotion lesson (cf. `fast_round_int`/#185): internal callers used
     these only on bounded well-behaved values, so nothing in-tree ever
     exercised the contract the docs were about to promise.
187. **n-ary logaddexp slice reduction** (tree or max+sum-exp) — pairs
     with #148.
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
189. **periodic_poly! dedup macro** once #43/#44 land (four
     near-identical folded-constant sinf_poly variants) — same
     macro-not-fn pattern as pi_reduce_and_poly!.

#### Batch 2: speculative / process

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
198. **Fix the mca latency harness's mix() sign blindness**: mix()
     erases the sign bit each chain hop, so sign-dependent work
     (cbrt, sin_checked's flips, erfcx's branch) is silently deleted
     from every latency number — inject alternating sign into the chain
     instead. Directly repairs a documented harness defect.
199. **Loop-invariant-2nd-arg hoisting in mca_target — found and fixed
     2026-07-27.** A second, independent harness defect of the same class
     as #198, worth recording because the numbers it produced looked
     plausible enough to sit in the table unnoticed.
     - The 2-arg idiom `{ let y = black_box(2.0); move |x: f32| f(x, y) }`
       is only sound when the function's expensive work depends on `x`.
       That holds for `atan2`/`hypot`/`powf`. It fails for
       `xlogy(x,y) = x*ln(y)` and `xlog1py(x,y) = x*log1p(y)`, where the
       whole transcendental sits on `y` alone: LLVM correctly hoists it out
       of the loop, so the region measured little beyond
       `x * precomputed_constant`.
     - Two consequences, both live until now. (a) The published numbers
       were meaningless — `xlog1py` read **9.03** cyc latency while
       `log1p` alone is 56.09 — and wrong in *both* directions, since the
       bogus throughput (1.858) was also worse than the true 2.301. (b) It
       produced a **false-positive `codegen_check` failure**: the scalar
       `vdivss` that check flagged was the hoisted loop-invariant divide,
       not per-element de-vectorization. `codegen_check` had been failing
       on `xlog1py_throughput` for exactly this reason; it now passes all
       148 regions.
     - Fix: pass `x` for both operands. Region went from 12 packed / 61
       scalar to **73 / 18**, matching `log1p`'s own profile. No CSE
       hazard when one operand is a plain multiplier and the other is the
       log's argument. Corrected: `xlogy` 21.57/3.506 -> **58.94/1.612**
       (vs `ln`'s 56.86/1.614), `xlog1py` 9.03/1.858 -> **57.23/2.301**
       (vs `log1p`'s 56.09/2.335).
     - **Reusable sanity rule**: check every mca row against the cost of
       its most expensive component. A composite cheaper than its parts
       means something was hoisted or optimized away, not that it is fast.
     - Swept the rest of the invariant-arg benchmarks against that rule;
       `xlogy`/`xlog1py` were the only two affected. The others are sound
       because their heavy work still depends on `x`: `compound` is
       `exp_checked(n*log1p(x))` (104.13 latency vs `log1p` 56.09 +
       `exp_checked` 46.06 = ~102, consistent), `powf`/`powf_pos`/
       `powf_unchecked`/`signed_pow` are all `exp2(y*log2(x))` with the
       log on `x`, and `diff_of_products`/`cross2` take `x` as a live fma
       operand. `pown`'s invariant `n` is `black_box`ed so its
       bit-testing loop cannot constant-fold (that case has its own
       `pown_const` regions, asserted separately by codegen_check).
200. **Auto-tune CI loop**: a scheduled job re-runs the tune.rs
     coordinate descent (LP-seeded) on every poly and files a PR when a
     real fuzz-verified improvement appears — automates the crate's
     single most-repeated manual win pattern.
