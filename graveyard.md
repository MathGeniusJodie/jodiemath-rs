# graveyard.md

Closed questions. Nothing here should be retried without a specific reason
to think the world changed -- a toolchain bump, a new instruction, or a
premise that has since gone stale.

Read this before proposing anything: most of it is *negative* results, and
negative results with numbers are the expensive part of this repo's history.
Where an entry says "rejected", it means measured and rejected, not guessed.

Split out of IDEAS.md on 2026-07-28. IDEAS.md now holds only what is still
open. Two rules that recur often enough to state up front:

- **A rejection can reject one *implementation* rather than the idea.**
  Three separate entries here were re-derived from scratch and shipped
  (`exp10_checked`'s round-based reduction, the `ln_normal` fma fold, the
  `pre_offset` dead-add removal). If a recorded failure was a single edge
  case or a scheduling artifact, and the measured prize was large, it is
  worth re-deriving -- and re-deriving beats resurrecting, because the old
  code's actual defect is usually not in the record.
- **Disproving an idea's stated *reason* does not disprove its *transform*.**
  Idea #12 sat dead for a year because a measurement killed its rationale
  ("moves work to integer ports"); the transform itself was a real 3->2-op
  win at nine sites, for an entirely unrelated reason.

## Tried and failed

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
  one catastrophically at +165.6%). No per-function scoping on stable
  Cargo, so a single catastrophic outlier vetoes crate-wide adoption.
- **PGO (+BOLT) probe on bench binaries**: inconclusive — this machine is
  too thermally noisy for wall-clock PGO-vs-baseline comparison (same
  binary swung 6.183→11.935 ns/op run to run), and `llvm-mca`'s
  single-region static model can't evaluate whole-program PGO effects
  (inlining/layout) at all.
- **Force zmm-width AVX-512 globally (`-C target-feature=-prefer-256-bit`)**:
  real, near-universal win (60/70 functions improved, median -19.4%) but
  `exp_checked` (+88.3%) and one other function (+84.5%) regressed
  severely — Tiger Lake only has one 512-bit-wide FMA port, and functions
  with a long serial poly chain lose the free cross-copy ILP masking the
  default double-unrolled ymm build gets. No per-function RUSTFLAGS
  scoping on stable Cargo, so the
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
  **`exp10_checked`'s half was rescued and shipped 2026-07-27** — the
  other five stay rejected, and their failures above are structural
  (`exp2`/`exp10` produce NaN in-domain, `exp2_checked`/`exp2m1` are pure
  cost) so don't retry those. What made this one different is that its
  blocker was a single edge case, not a mechanism: re-derived from scratch
  rather than resurrected, `exp10_checked(inf)` comes out `+inf` at the
  *shipped* clamp bounds, i.e. the bug does not reproduce and no clamp
  move was needed. The exact mechanism of the original failure was never
  recorded, so what differs isn't knowable — most likely the combine or
  the fit, since this version also does not reproduce the max 1→2
  regression. Concretely: keep `round`'s own centered `f in [-0.5, 0.5]`
  and drop `exp10_reduction!`'s floor-adjust (a compare, a select and two
  add/subs), with a `Q(f)` refit for the centered domain
  (`exp2_q_poly_centered!`). Centering is an affine change of variable, so
  the same degree 5 buys the same accuracy — idealized max relative error
  1.073e-8 vs the `[0,1)` fit's 1.217e-8.
  - Exhaustive over all 2^32: avg ulp **0.0307 -> 0.0082** (3.8x better),
    max **1, unchanged**. Plain `exp10` untouched (0.0307/2), and no other
    row in the mca table moved.
  - mca: throughput **2.361 -> 1.897 (-19.7%)**, latency **63.56 -> 51.06
    (-19.7%)**. All 8 gates pass, including the `edgecheck` and
    `saturation_pins` entries that exist because of this exact bug;
    5 `worst_corpus` entries move sub-ulp (3 closer to the true value,
    2 further) and were re-blessed.
  - Reusable: this is the third time this session a rejection turned out
    to be rejecting one *implementation* rather than the idea. A recorded
    failure that is a single edge case, on an idea whose measured prize
    was -30% throughput, is worth re-deriving from scratch — and
    re-deriving beats resurrecting, since the old code's actual defect is
    usually not in the record.
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

- **logit: `2*atanh(2p-1)` central band** (2026-07-31, **SHIPPED**, a
  real bug fix). Exhaustive over the domain: avg ulp 0.2758 → 0.2625,
  **max 1024 → 3**, for llvm-mca latency 59.03 → 59.91 (+1.5%) and
  throughput 3.176 → 3.551 (+11.8%).
  - **The defect**: `ln(p) - log1p(-p)` subtracts two nearly-equal
    logarithms near `p = 0.5`. Both are `~-ln(2)`, each carrying its own
    rounding of a quantity ~0.693, while their difference is only
    `~4*(p-0.5)`. Error in ulp of the result is roughly
    `6e-8/ulp(result)`, so it grows smoothly as the result shrinks —
    there is no threshold below which it is safely small. This was
    documented as a "near a true zero, ulp isn't meaningful" artifact
    and was **not** one; it was found by re-auditing every user of that
    excuse after the same excuse turned out to be hiding a real bug in
    `cospi` (see the trig section).
  - **The fix**: `logit(p) = 2*atanh(2p-1)` on `|2p-1| < 0.25`, reusing
    `atanh_small` — the crate's existing poly, on its own fitted domain,
    no new fit. Nothing cancels: the result is proportional to `2p-1`,
    exact for `p >= 0.25`, and the doubling is exact. Outside the band
    the old difference form stays untouched (at the seam the result is
    ~0.51 against operands ~0.69, so ~1 ulp) and keeps denormal `p`, the
    endpoints and out-of-domain `p` correct. Seam gaps 9 / 5 ulp with
    matching one-sided slopes (4.2654 vs 4.2682 against a true
    `1/(p(1-p)) = 4.2667`).
  - **REJECTED alternative — one arm for the whole domain**,
    `mulsign(log1p(|2p-1| / min(p,1-p)), 2p-1)`. Mathematically the
    nicer object: a single `log1p`, no `ln` at all, both operands
    Sterbenz-exact across the central band, quotient non-negative over
    the whole domain, and it measured a genuinely *better* avg (0.2408
    vs 0.2625) at the same max of 3. Killed on cost: **+34% latency**
    (59.03 → 79.19) and +27% throughput (3.176 → 4.028). The reason is
    structural and worth remembering — it needs *two* divisions (the
    quotient, then `log1p`'s own `c/u` correction) and they land in
    **series**, where the shipped form's `ln` and `log1p` are
    independent and pipeline in **parallel**. Fewer total operations,
    longer critical path. It also needed a denormal guard the difference
    form gets for free: `1/p` overflows below `p ~ 2.9e-39`, returning
    `-inf` for a true `~-88`, fixed by scaling the denominator `2^24`
    and passing `koff = 24` to `ln_normal` (sound only because the
    quotient is past `2^25` there, where `log1p` has already degenerated
    to `ln`) — and that scale had to be gated on `d > 0.0`, or it pulls
    an out-of-domain negative quotient back above `-1` and returns a
    plausible finite number where `NaN` is owed (`logit(-0.1)` → -16.6).

- **`log_2`/`ln`/`log10` integer koff fold**: bit-exact, but mca showed
  zero measurable change — LLVM already performs this reordering.
- **`ln_normal`/`log10_normal`: fuse trailing `+k*LN2_LO` into the fma**:
  mca latency deterministically *worse* by 1 cycle (reproducible,
  re-confirmed on a second attempt after an initial stale-baseline
  false-positive); accuracy unchanged, the term was already off the
  critical path. **Re-screened under idea #22 on 2026-07-27 and the
  result flipped — `ln_normal`'s half now ships** (see lib.rs/git log).
  Two things had to change for that:
  - The *association* matters, and the original attempt evidently used
    the other one. `fma(p, s, fma(k, LN2_LO, k_hi))` (fold the LO word
    into the `k` term first, poly last) wins across the board on
    rustc 1.98.0-nightly; `fma(k, LN2_LO, fma(p, s, k_hi))` (poly first)
    still reproduces the documented +1-cycle latency regression on the
    same toolchain. Both are 2 fma against the old 1 fma + 1 mul + 1 add,
    so this is pure scheduling, not op count — which is exactly the
    "silently flips" class #22 predicted.
  - **`log10_normal`'s half is genuinely worse and was dropped.**
    Identical transform, opposite accuracy outcome, confirmed
    exhaustively: `log10` avg/max **0.1265/3 -> 0.1280/4** and `log10p1`
    **0.2319/3 -> 0.2347/4**, where `ln` and `log1p` come out *unchanged
    to every digit and at the same worst x*. Worth remembering that a
    Cody-Waite HI/LO reassociation is not automatically transferable
    between two functions with the same shape — `LOG10_2_LO` is 3.2x
    larger relative to its HI word than `LN2_LO` is, so it survives the
    early fold less cleanly. Shipping only the `ln` half keeps every
    perf win below (they all route through `ln_normal`) at zero
    accuracy cost.
  - Measured (mca, ln-only): `ln_unchecked` 38.22/1.113 ->
    **34.06/1.018** (-8.5%), `compound` -4.5%, `acosh` -3.9%, `log1pmx`
    -3.9%, `asinh` -3.3%, `log10p1`/`log1p` -2.0%, `xlogy` -1.6%,
    `logit` -1.2%, latency down 1-8 cycles on all 14 affected rows and
    up on none. Total instruction count in `mca_target.s` drops 380134
    -> 378845.
  - **mca reported a `probit` throughput regression of +15.1% that
    wall-clock flatly contradicts** — quickbench min-of-7, 3 reps each
    side: `probit` 1.599 -> **1.527 ns/op** (*faster*), `erfinv` 1.421
    -> 1.423 (wash). This is the second recorded mca-vs-wall-clock
    direction disagreement (after `sincos_checked`) and the reason to
    trust wall-clock here is independent of both: the change strictly
    *removes* 1289 instructions, so a real 15% slowdown would need a
    mechanism, and none is visible. quickbench confirms the wins too
    (`ln_unchecked` 0.307 -> 0.242 ns/op, -21%; `acosh` -12%;
    `log1pmx` -11%; `asinh` -7%). Feeds idea #96's case for a
    thermal-controlled wall-clock harness.
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
- **`log_family_edges!` non-negative variant for powf/powf_pos/rootn**:
  those callers pass `x.abs()` (or contract for `x >= 0`), which makes the
  domain-error arm unreachable — `x <= 0.0` can only mean `x == 0.0`, so
  the `-inf`-vs-`NaN` select collapses. The transform is real and the two
  instructions *do* disappear (`powf_throughput` 217 → 215), but mca
  throughput got 1.1% **worse** and latency stayed flat: another
  removing-an-op-reschedules-worse artifact, so not worth a new macro plus
  a call-site indirection. Worth re-screening on a toolchain bump — this
  is the exact shape that flipped for `round_x_over_pi`'s dead add.
- **log_2: atanh-form reduction `t=(m-1)/(m+1)`**: 1000x tighter in the
  underlying continuous math, but worse for real — accuracy slightly
  worse (log_2 already deep in f32-rounding-noise territory) and latency
  +43% (the division depends on `s` from the first step with nothing to
  overlap it against, unlike cbrt's early-starting reciprocal).
  **Still rejected for `log_2` itself, but the transform shipped in
  `log2_df`** (see §cbrt/powf below) — a third instance of the "rejection
  rejects one implementation" rule, with both halves of this entry going
  stale for the *other* caller. Accuracy: a single f32 output has nothing
  to spend 1000x tighter math on, but `log2_df`'s output is amplified ~7
  bits by `powf_checked`'s `y`, so there the extra bits are exactly what
  is short. Latency: the division was never load-bearing — `d = m+1`
  spans a factor of 1.41, narrow enough that a degree-3 minimax seed plus
  one Newton step beats it, and the quotient-refinement step squares
  whatever error is left anyway. Do not read this entry as "the atanh
  form costs a division."
- **ln_accurate/log2_accurate tier from log2_df**: no measurable accuracy
  benefit — 20M-sample fuzz gave identical avg/max ulp to plain `log_2`
  (0.0061/3 both), and a 1.63B-sample strided sweep found only 7
  bit-differing outputs (~4e-9 of the domain). `log2_df`'s extra
  double-float precision only matters once something *downstream*
  amplifies the preserved low-order bits (e.g. `powf`'s multiply by y) —
  collapsing straight back to a single f32 with no such amplification
  lands on the same correctly-rounded result almost every time.
  **Both halves of this went stale and the conclusion is now wrong.** It
  was measured against the *pre-atanh* `log2_df`, which was itself only
  ~2^-23 relative (see the powf entry in §cbrt/powf — that is the same
  category error that held `powf_checked` at 203 ulp), so "identical to
  log_2" said nothing about a real double-float log2. And `log_2` itself
  is max 1 now, not 3, so there is even less left to take. The re-read:
  this entry never tested what it claimed to, and the thing that actually
  fixed `log_2` was free — see idea #202's peel, which needed no extra
  precision at all, only a different place to put the leading term.
- **compound's "underflows to zero" excuse: measured false, third time
  this species has hidden a real defect.** `compound`'s doc claimed its
  ~200 max ulp came from `(x,n)` pairs whose true value is below any
  denormal. Instrumented over 20M uniform samples: of the 5440 worst
  (>50 ulp), **5435 have a perfectly normal result** and every one has
  `|n*log1p(x)|` between 30 and 86. It is the plain `y`-amplified error
  `powf_df_mag!` documents — `exp`'s error is `|n*log1p(x)| *
  relerr(log1p)`, and both `log1p`'s ulp *and* the `n*log1p(x)`
  product's own rounding go through that multiplier, which reaches ~88
  before the result overflows. Same excuse, same shape, as
  `cospi`/`atan2_pos`.
  - **`compound` is not redundant, though** — the obvious "just call
    `powf(1+x, n)`" scores **1.07e9 max / 5.0e6 avg** on the same
    filtered distribution (finite `n`, finite nonzero result) where
    `compound` scores 240/0.29. Forming `1+x` in f32 is exactly the
    catastrophe it exists to avoid; the 200 ulp is a second, unrelated
    problem.
  - **Fixed as a second tier, not in place**: `compound_accurate` =
    `exp2_checked_df(log2p1_df(x) * n)` with a new double-float
    `log2p1_df`. Max ulp **236 -> 4-5**, avg **0.195 -> 0.019**. Cost is
    real and is why `compound` stays: throughput **3.724 -> 9.239
    cyc/elem** (2.5x), latency 96.75 -> 153.56.
  - Two things had to be right in `log2p1_df` that `log2p1` gets away
    without. (1) The `c = x - (u-1)` correction needs its *own* low
    word: whenever `1+x` rounds back to exactly `1`, `log2_df(u)` is
    exactly `Df32(0,0)` and the correction **is** the whole answer, so
    `LOG2_E` has to be split hi/lo and the `c/u` division has to carry
    its residual (`fma(-eh, u, c) * rcp`). Skipping just the division's
    residual measured **44** max ulp instead of 5 — one rounding, 9x the
    error, because `n` amplifies it identically to the leading term's.
    (2) The `-e^2/2` Taylor term is load-bearing for the same reason.
  - And the degenerate-`u` override has to replace the **whole pair**,
    not just the high word the way `powf_df_mag!` does. There the low
    word is only ever read by `exp2_checked_df`'s non-finite guard; here
    the two-sum that folds the correction in feeds it back into the high
    word, so an `inf - inf` residual turned correct `+-inf` answers into
    NaN. Caught by edgecheck (`compound_accurate(-1,5)`, `(-1,-5)`,
    `(inf,1)`), invisible to both fuzz and mca.

### sin / cos / tan / sinpi / cospi / sind / cosd / tanpi / tand

- **cospi: reduce on `round(x)` instead of `round(x-0.5)`** (2026-07-30,
  **SHIPPED**, a real bug fix, not a tuning change). Exhaustive, all 2^32
  patterns: avg ulp 0.2813→0.0578, **max 868814811→2**, and it is
  *cheaper* -- llvm-mca latency 51.00→47.00 cyc, throughput
  1.283→1.226 cyc/elem (48 vs 51 instructions in the throughput region;
  `acospi` unchanged as a control row). `cos2pi`, which is
  `cospi(2*x)`, inherits it: 0.2835/868814811 → 0.0583/2.
  `sinpi`/`tanpi`/`tan2pi`/`sind`/`cosd` are untouched and measured
  bit-identical.
  - **The old form's actual defect**: `k = round(x-0.5)` then
    `r = (x-k)-0.5`. `x-k` lands in `[0,1]`, whose ulp is *coarser* than
    that of an `x` just inside a smaller binade, so the subtraction
    rounds away `x`'s low bits -- and `r` is exactly the quantity the
    answer is proportional to near a zero. At `x = -0.49999997` (the
    first f32 below `-0.5`) it returned a flat `0.0` for a true
    `9.36e-8`: a 100% relative error, which is where the 8.7e8 came
    from. The fix reduces like `sinpi`/`tanpi` already did (`r = x - k`,
    exact for every f32) and reflects: `cos(pi*r) = sin(pi*(0.5-|r|))`,
    where `0.5-|r|` is Sterbenz-exact over `|r| >= 0.25` -- i.e. exactly
    the half of the domain holding the zero. Where it *does* round
    (`|r| < 0.25`) the result is within a quarter turn of `+-1`, so the
    rounding lands where the function is flat.
  - **Why it survived so long**: it was written off four separate times
    as a "near a true zero, ulp isn't meaningful" artifact (see the
    `sinpi/cospi dedicated poly` entry below, and the cross-references
    that used to point at it from `sinc`/`softplus`/`tan_checked`/
    `wrap_pi`/`logit`/`compound`). That excuse is only valid **when the
    reduction is exact at the zero** -- `sinpi` scores max 2 ulp at its
    own zeros for precisely that reason, so `cospi` sitting at 8.7e8 at
    *its* zeros was the tell, not the proof. Before accepting a
    near-zero blowup anywhere else, check the sibling that has the same
    shape of zero.
  - Side effect, deliberate: every zero is now `+0.0` (was alternating
    `-0.0`/`+0.0`), matching IEEE 754-2019's `cosPi(n+1/2) = +0`.
    `round_ties_even` sends `k` to the even neighbour at a half-integer,
    so the parity sign is `+1` there for free. `cospi` is also now
    exactly even (`cospi(-x)` and `cospi(x)` agree bit for bit), and no
    longer needs `sinf_poly`'s copysign, since `0.5-|r|` is never
    negative. All pinned in edgecheck.rs.

- **round_x_over_pi: remove dead pre_offset=0.0 add**: instructions
  confirmed gone via asm, but throughput got *worse* — removing an op let
  the scheduler pick worse elsewhere. **Re-screened under idea #22 on
  2026-07-27 and flipped: now a real win, shipped** (see lib.rs/git log).
  Second of the three tagged scheduling artifacts to reverse on
  rustc 1.98.0-nightly.
  - Mechanism, which is worth stating because "dead add" undersells it:
    `+ 0.0` is *not* an identity LLVM is allowed to delete — it maps
    `-0.0` to `+0.0` — so sin's callers really did pay an add that only
    ever normalized a sign of zero. Turned into a `const HALF: bool`
    generic, so cos still gets its `- 0.5` and sin gets nothing.
  - mca: `tan_checked` 8.787 -> **8.235 (-6.3%)**, `sin_checked` 5.311 ->
    **5.232 (-1.5%)**, latencies flat, and no row anywhere regresses by
    more than 0.001 (noise).
  - Accuracy: unchanged everywhere. `worst_corpus` bit-identical across
    all 108 functions, every gate passes, `sin_checked`/`cos_checked`'s
    four documented buckets identical, `tan_checked`'s whole-domain max
    identical at the same worst x.
  - **One methodology note that nearly caused a false rejection**:
    quick-fuzz `wrap_pi` reported max ulp 12/2/12 against a baseline's
    1/1/2 over three repeats each, which reads like a real regression.
    It is entirely sampling noise — exhaustive sweeps of both sides give
    *identical* 0.0244 avg / 41 max at the same worst x. `wrap_pi`'s max
    lives on a near-total-cancellation artifact at multiples of `2*pi`
    that random sampling hits or misses, so its quick-fuzz max is
    meaningless as an A/B signal even with repeats. Same lesson as the
    2-arg repeat-run entries, but here repeats were *not* enough — the
    exhaustive sweep was.
- **round_x_over_pi: qh → round_ties_even**: regressed cos_checked's max
  ulp 2→6 (exact-half ties clash with cos's -0.5 offset). Reverted to
  `f32::round`.
- **reduce_pi: rebalance 4-deep chain to depth 2**: bit-exact but +3 cyc
  latency both functions. **Re-screened under idea #22 on 2026-07-27:
  does NOT flip — still rejected.** `(e1b + e2b) + (e3b - e3t)` against
  the shipped `e1b + e2b + e3b - e3t` reproduces the original +3 cycles
  exactly on rustc 1.98.0-nightly (`sin_checked` 117.02 -> 120.02,
  `cos_checked` 122.00 -> 125.00, `wrap_pi` 92.02 -> 95.02,
  `sinc_unnormalized` 128.00 -> 131.00, `tan_checked` 138.20 -> 141.20)
  and throughput is mixed-to-worse on top (`sin_checked` +0.8%,
  `tan_checked` +0.8%). Note the shipped order is already
  latency-optimal in the way that matters: `e3t` is the last-ready value
  (it carries the `ql` dependency) and the flat left-to-right chain
  consumes it *last*, whereas any depth-2 pairing has to consume it one
  level earlier. Shortening the tree does not help when the tree's
  depth was never the binding constraint.
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
  win. This was the crate's last full `two_sum` call inside `reduce_pi`;
  `two_sum` stays regardless as a documented public EFT, and
  `doublefloat.rs` keeps its own copy.
  - **Re-screened once the `POLY_SAFE_BOUND` clamp came out of
    `sin_checked`/`cos_checked`, and it flipped: shipped.** With both
    callers two ops shorter the "diverges by caller" split is gone --
    `sin_checked` 4.854 -> **4.698 (-3.2%)** and `cos_checked` 4.349 ->
    **4.079 (-6.2%)**, with instrs (-6/-6), uOps (-6/-5) and BlockRT
    (-3/-3) all corroborating, and both latencies improving. The +13.7%
    cos regression was a scheduling artifact of the *old* surrounding
    code, not a property of the downgrade. Same lesson as the two
    `Re-screened under idea #22` entries above: a rejection measured
    against a since-changed neighbourhood is stale evidence.
  - The ordering argument is a proof, not a hope: `|tier2|` is bounded
    by `~|qh|*2^-46 + |ql|*2^-23` while any nonzero `ql` forces
    `|p3| = |ql|*PI_HI >= pi`, so `|a| >= |b|` holds by many orders of
    magnitude across the whole f32 line -- and at `ql == 0` Fast2Sum is
    exact anyway (`e = b - (s - 0)` is exactly zero). Exhaustively
    bit-identical to the `two_sum` form over all 2^32 inputs for
    `sin_checked`, `cos_checked`, `reduce_pi_checked` and
    `reduce_pi_half_checked`.
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
  (0.1079→0.1105); its already-huge near-x=-0.5 max-ulp outlier (called
  a known "near a true zero of the function, ulp isn't meaningful"
  artifact, not a bug, at the time -- **that call was wrong, see the
  `cospi` reduction entry below: it was a real bug, fixed 2026-07-30**)
  happened to
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
- **atan2d's 30 max ulp was a denormal *intermediate*, not the poly, not
  the constant.** Instrumented: **every** sample over 8 ulp has a
  denormal `atan2(y,x)` result, and `atan2`'s own error at those inputs
  measures **0 ulp** — it is exactly right, it just cannot represent the
  answer, because its internal `y/x` underflowed. `180/pi > 1` then
  carries that value up into a binade with more mantissa bits than the
  denormal ever had, and the missing ones show. `atan2pi` scores 4 on
  the same plane purely because `1/pi < 1` moves the other way. 1.3% of
  20M uniform samples scored over 8 ulp; a sixth of those have a
  perfectly normal *output*, so this is not confined to denormal
  results.
  - Fixed with a second branch: a denormal `atan2` means the angle is
    deep inside `atan`'s linear regime, where the answer is just
    `y * ((180/pi) / x)`. **Max ulp 30 -> 5, avg 0.458 -> 0.175.** Not
    free, and the cost is real by the full ladder (instructions 65 ->
    72, uOps 67 -> 74, `Block RThroughput` flat): throughput 1.881 ->
    2.128 cyc/elem (+13.1%), latency 71.19 -> 78.10 (+9.7%).
  - **The association is the whole trick, and the obvious one is
    wrong.** `(y * (180/pi)) / x` scores **77** max ulp, worse than some
    runs of the original: for a denormal `y` the scaled numerator is
    *still denormal*, so it reintroduces the identical bug one step
    earlier. `y * ((180/pi) / x)` keeps every intermediate normal —
    `(180/pi)/x` cannot be denormal on this branch, since a nonzero `y`
    there forces `|x| >= 1.2e-7`. Same op count, 15x the accuracy.
  - Also: quick-fuzz max for the wrong association swung **7 / 11 / 77**
    across three consecutive runs while the right one is a flat 5 four
    times. A 2-arg max that moves by 10x run to run is itself the
    signal that a narrow input class is being hit at random, not noise
    to average away.
  - Untested but the same shape: `asind`/`acosd`/`atand` all multiply by
    `180/pi` and `asin(x) ~ x` can be denormal too. `asind` is 11.
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
  — **re-screened under idea #22 on 2026-07-27: does NOT flip.**
  Reproduces the original numbers almost exactly on rustc
  1.98.0-nightly: `atan` 1.491 -> 1.529 (+2.5%), `atan2` 1.694 -> 1.731
  (+2.2%), `carg` +1.4%, with only `atand`/`atanpi` improving (-0.8%).
  Latency flat. Reverted again. Original entry follows.
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
- **atan2_pos's "branch-cut artifact" excuse: measured false, it was a
  real bug.** `atan2_pos`'s doc comment used to claim its ~1.09e9 max ulp
  was a reference-comparison artifact at the wrap seam, the sign of a
  near-zero `atan2` differing between f32 and the f64 reference. It is
  not: instrumenting every sample over 1000 ulp showed **100% of them are
  one mechanism**, `y/x` underflowing to `-0.0` in f32 where f64 has
  plenty of exponent left. `atan(-0.0)` is `-0.0` and the `x>0`
  correction term is `-0.0`, so `r < 0.0` reads a genuinely negative
  angle as non-negative, skips the fold, and returns `~0` for a true
  angle of `~2*pi`. The sign never differs; the magnitude flushes. Not a
  rare corner either — **2.1% of 20M uniform-bit-pattern samples**
  (420587/20000000), which is why the *avg* was 2.3e7. Fixed by keying
  the fold on `y`'s sign bit: max ulp 1.09e9 → 3, avg 2.29e7 → 0.063,
  mca latency flat at 72.19 and throughput **1.993 → 1.856 cyc/elem
  (−6.9%)** — same instruction count and same uOps, the mask just moves
  off the dependency chain onto an argument. The tell-tale worth reusing:
  the reported max was *exactly* `ulp_diff(0.0, TAU)` = 1086918619, one
  constant rather than a distribution, which is what a systematic
  full-turn miss looks like and a genuine seam artifact does not.
  - **The exact-reference variant `y < 0.0 || r < 0.0`** scores the same
    3 max / 0.063 avg and additionally keeps `atan2_pos(-0.0, x>=0)` at
    `-0.0` (matching a plain `r < 0.0` f64 fold), but costs `vcmpltps` +
    `korb` per vector: latency 72.19 → 74.09 (+2.6%), throughput 1.993 →
    2.037 (+2.2%). Rejected — 9% of throughput for one input pair of
    measure 2^-33, and folding `-0.0` is arguably the better convention
    anyway since it keeps `-0.0` out of the range of a function whose
    name says non-negative.
  - The advertised range `[0, 2*pi)` was never achievable:
    `f32::consts::TAU` is *above* `2*pi` and is the correctly-rounded
    answer for the last half-ulp of the turn. Doc now says `[0, TAU]`
    closed.

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
- **exp_pos_neg: fold the 0.5 into the *poly constants* instead of a
  multiply** — the version of idea #32 that works. **Shipped
  2026-07-27** (see lib.rs/git log). The rejected attempt below relocated
  the halving to a `t1 * 0.5` **multiply**, which is why it read as a
  latency/throughput tradeoff: it moved an op, it never removed one.
  Halving `exp_pos_neg_core!`'s four fitted coefficients and its two
  leading `1.0`s costs **nothing** — the constants are compile-time, and
  scaling every operand of an fma by the same power of two scales its
  correctly-rounded result by exactly that power of two, so `p_pos`/
  `p_neg` come out bit-for-bit halved.
  - Removes the multiply from all four callers at once: `sinh`'s and
    `cosh`'s trailing `0.5 *`, and both of
    `exp_pos_neg_checked_half`'s `t1 * 0.5` / `t1n * 0.5`.
  - **Free correctness win as a side effect.** `exp_pos_neg_checked_half`
    documents (its own point 2) that the `0.5` must be applied before the
    field multiply or `p_pos * t1 * t2` overflows to `inf` while
    `sinh`/`cosh` are still finite. Halving in the poly satisfies that by
    construction, so *plain* `sinh`/`cosh` inherit the fix they never
    had: `sinh(88.7228)` was `inf`, is now **1.7014122e38**, and
    `sinh(89.4)` was `inf`, is now **3.348863e38** (both verified against
    the true values). 8 corpus entries move, all `inf` -> correct finite.
  - Accuracy otherwise untouched: `worst_corpus` shows *only* those 8
    entries, so `sinh_checked`/`cosh_checked`/`coshm1` are bit-identical
    (as the algebra predicts — `p*(t1/2)*t2` and `(p/2)*t1*t2` are the
    same real number at every step). Exhaustive `sinh` 0.0821 avg / max 5,
    `sinh_checked` 0.0429 / 5, and the whole `cosh` family's avgs
    unchanged to four digits.
  - mca: latency down 2-4 cycles on all five affected rows (`sinh`
    56.00 -> 52.00, `cosh` 55.00 -> 51.00, `sinh_checked` 58.06 -> 56.06,
    `cosh_checked` 58.06 -> 55.06, `coshm1` 70.06 -> 68.06). Throughput
    `cosh_checked` **-7.3%**, `sinh` -1.1%, `sinh_checked` -0.2%,
    `cosh` +0.9%.
  - **Applied to the narrow tier too** (`exp_pos_neg_narrow_half`, a
    standalone copy of the same poly): `sinh_narrow` 52.00/1.630 ->
    **48.00/1.524 (-6.5%)**, `cosh_narrow` 51.00/1.350 ->
    **47.00/1.274 (-5.6%)** — clean wins on *both* axes here, with
    `worst_corpus` bit-identical and avg/max ulp unchanged (0.0821/5 and
    0.0589/5).
  - **`coshm1` +61.6% is an llvm-mca artifact, and this one is provable
    rather than merely suspected** — worth recording as the cleanest
    example yet of the model diverging from the machine. The `coshm1`
    region's opcode histogram is *identical* before and after except for
    **4 fewer `vmulps` and 1 fewer `vbroadcastss`** (107 -> 102
    instructions), and it is 100% packed `ymm`/`zmm` both ways, so there
    is no de-vectorization, no new divide, no changed instruction class —
    just strictly less of exactly the same work. A shorter, otherwise
    identical, fully-vectorized instruction stream cannot take 1.6x the
    cycles. Whole-target instruction count also drops 378296 -> 377953.
    Wall-clock could not arbitrate (`sinh_checked` alone spanned
    0.919-2.051 ns/op *within one configuration*, the documented 3-9x
    environment swing), which is exactly why the asm-level check was the
    decisive evidence — reach for the opcode histogram, not the
    stopwatch, when mca reports something structurally impossible.
  - **Sharper arbitration rule, established on `acosh` 2026-07-28: use
    `Block RThroughput` and the uOp count, not the simulated cycle
    count.** `acosh`'s two non-finite selects merge into one exactly as
    `asinh`'s do. The merge is -6 instructions with **every expensive
    opcode identical** (`vdivps` 4, `vsqrtps` 2, 44 fma, 12 `vmulps`,
    14 `vaddps` — only mask/blend bookkeeping moves), yet mca prices it
    at **3.569 -> 3.828 (+7.3%)**, reproducibly, and a second independent
    formulation reports +5.9%. `--bottleneck-analysis` shows why, and
    exonerates it:

    | | old | new |
    |---|---|---|
    | Total uOps | 16700 | **16100** |
    | Block RThroughput | 40.0 | **40.0** |
    | Resource Pressure | 78.4% | 58.1% |
    | Register Dependencies | 42.2% | 48.6% |

    Fewer uOps, **identical** resource-limited throughput bound, and the
    reported slowdown is entirely mca's simulated schedule serializing on
    register dependencies. The throughput harness feeds 16 *independent*
    elements, so on real hardware there is nothing for those dependencies
    to serialize — the out-of-order window interleaves them and the
    binding constraint is the port pressure that `Block RThroughput`
    already says is unchanged. Adopted.
  - That is the same test that correctly *rejected* `erfinv`'s fold in
    the evaluation-order sweep, where uOps went the other way (19400 ->
    20300 while instructions fell 177 -> 164). So the pair is a clean
    discriminator: **instructions down + uOps down + RThroughput flat =
    take it; instructions down + uOps up = real regression.** The
    reported cycle count distinguishes neither case.
  - **The `coshm1` row is not just occasionally wrong on a delta, it is
    wrong in absolute terms — established 2026-07-28 from a *static*
    comparison, no A/B needed.** `coshm1` is literally
    `sinh_checked(x*0.5)` squared and doubled, and the two throughput
    regions' opcode histograms are the same instruction for instruction
    apart from `coshm1`'s **+4 `vmulps`, +2 `vaddps`, -2 `vmovups`,
    -1 `vbroadcastss`** (99 -> 102 instructions), both 100% packed, both
    16 elements wide, neither containing a divide or a sqrt. mca prices
    them at **2.340 vs 4.700 cyc/elem** — 2.0x apart for 3% more
    instructions, i.e. 2.65 IPC against 1.36 IPC on near-identical
    streams. So the earlier +61.6% delta was not a one-off: this
    function's whole row is mispriced, and its *level* should not be
    compared against any other row's.
  - Cheap standing check that finds rows like this without an A/B:
    divide each `*_throughput` row by its region's instruction count.
    Functions dominated by `vsqrtps`/`vdivps` (`rsqrt`, `rhypot`,
    `normalize2/3/4`, `sqrt1pm1`, `pow_3_2`) correctly come out high,
    and everything else lands in a tight band around 1.3-2.4 cyc per
    100 instructions — `coshm1` at 4.61 is the only row that is an
    outlier without a divide or a sqrt to explain it.
- **exp_pos_neg: return halves pre-scaled by 0.5 for plain sinh/cosh
  too, via a multiply** (idea #32, the *plain-multiply* relocation `exp_pos_neg_checked_half`
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
- **`erfcx_pos`'s `c0` as a single f32 word — *both* choices**: superseded,
  and the reason is worth keeping because the old entry read like a closed
  question. `c0 = P(0)` carries the `1/(x·sqrt(pi))` asymptote alone, so
  past |x| ~ 10 the result is `c0*v` and nothing else and c0's own error
  arrives as *bias*, not noise — a flat −0.63 ulp signed mean across every
  band from [10,20) out to [1e10, f32::MAX). The catch is that 1/sqrt(pi)
  is a near-tie in f32: 0.49 ulp above the nearer representable value,
  0.51 below the other. So neither single word helps. The shipped
  1-ulp-low choice gave `erfcx (x>=20)` 0.6478 avg / 4 max; pinning to
  correctly-rounded only flips the bias's sign (0.6245) *and* costs `erfc`
  a max ulp, which is what the earlier "pinning measured worse" note
  recorded. Holding `c0` in two words and peeling the high word into the
  final fma removes the bias rather than relocating it: tail 0.6478 →
  0.2684 avg, max 4 → 2, signed mean −0.63 → +0.00, with `erf`, `erfc`
  and `erfcx` all improving on their other rows too. **The recorded
  rejection was about pinning and was never evidence against splitting** —
  same shape as the "disproving the reason ≠ disproving the transform"
  rule at the top of this file. The split does require re-polishing: with
  an exact `c0` the unconstrained fit misses `erfcx_pos(0) == 1.0`, which
  costs ~1 ulp of flat bias over the whole near-zero region until c1..c4
  and c10 are each nudged an ulp to restore it.
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
  cost that tier already pays. **The diagnosis held; the last clause was
  wrong** — `log2_df` was *not* in fact a higher-precision log2, only a
  log2 with a spare low word, which is not the same thing and left
  `powf_checked` at ≥203 ulp. See the shipped entry below. Still true for
  plain `powf`, which has no way to buy relative precision cheaply and
  stays a fast tier at ~292.
- **powf_checked: atanh-form `log2_df`** — **SHIPPED**, ≥203 → 3 max ulp
  (`examples/powfsearch.rs`), avg 0.046 → 0.019, for *no* latency cost
  (mca -0.1%, interleaved wall-clock -2.6%) and +7-9% throughput. The bug was a category error worth remembering: a
  double-float carrying the *rounding error of the last two operations*
  reads like a higher-precision result but is not one. `log2_df` was
  `k + two_product(p, s)` over `log_2`'s own degree-9 poly, so it captured
  the final multiply and add exactly — while `P(s)` itself was still
  evaluated in plain f32, and for `x` near 1 (`k == 0`, exactly where
  large `|y|` makes the amplification bite) the result *is* `s*P(s)`, so
  `P`'s own evaluation rounding sat on the answer at ~2^-23 relative with
  nothing downstream able to recover it. The 1.5x that `powf_checked`
  beat `powf` by was the multiply and the collapse; the poly's rounding
  was common to both. Fixed by changing the *shape*, not the bookkeeping:
  see `log2_df`'s doc comment. Generalises — any `_df`/`_checked` tier
  whose extra precision comes from EFT bookkeeping around an f32 kernel
  is capped by that kernel, and the way to check is to score the df
  against f64 in *relative* terms over the octave where the exponent is
  zero, which is where a cancelling `k` stops hiding it.

  The first working version cost +14% latency / +24% throughput; a
  separate pass took that to zero and +7-9% **without touching a single
  coefficient**, purely by rescheduling, and the three levers generalise
  to any double-float kernel here:
  - **Feed the seed poly the un-offset variable.** The reciprocal seed
    fitted in `d = m+1` had to wait on that add; refitted in `m` (an
    affine change of variable, so *the same fit*) it starts a level
    earlier, and `d` is still computed for later use where it is off the
    path. Free.
  - **A correction term's poly argument may be sloppy even when its
    prefactor may not.** `t^3 * G(t^2)` was entirely downstream of the
    refined quotient. But `G` varies ~2% across its whole range, so it
    can take the *unrefined* `th^2` and lose nothing — and the prefactor
    gets its refinement without a squaring, via
    `t^3 = th^2 * (th + 3*tl)` (dropping `3*th*tl^2`, ~2^-44 relative).
    The whole correction now runs *beside* the refinement chain with one
    fma downstream of it instead of five ops. This was most of the win.
  - **Exact two-sums commute, so order them by operand readiness.**
    `(k + hi) + corr` instead of `k + (hi + corr)`: `hi` is ready long
    before `corr`, both orders keep the `quick_two_sum` precondition, and
    the result is identical.

  Two things that did *not* work, both measured: Estrin-splitting the
  correction poly (it had already left the critical path, so the extra
  multiply bought nothing), and the reverse — Horner-ing the *seed* poly
  to save that multiply, which made **both** metrics worse, the usual
  "removing an op let the scheduler pick worse elsewhere."
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

## Backlog entries that closed

Entries that started life in IDEAS.md's untried backlog and were resolved
negatively -- screened out on a number, or found moot once a prerequisite
was rejected.

- **Explicit core::simd fallback tier / f32x16 AVX-512 via
  #[target_feature]**: only worth it if a LUT idea ever survives
  screening, or as a per-function-scoped way to capture the AVX-512
  zmm-width win found (and rejected as a global default) above.

- **FTZ/DAZ feature flag**: a cargo feature assuming caller-side FTZ/DAZ
  would let every denormal branch (log family, cbrt, exp2_checked) become
  dead code.

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

9. **erf joint Pade+erf_poly refit with explicit crossover-region
   weighting** — the rejected num-only LP regressed exactly at the 0.28
   seam, which a joint objective would score directly.

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
    case is now closed. **Closed outright 2026-07-28: there is no
    remaining target.** The technique only pays where a range proof
    changes *lowering*, and in this crate that means a float->int cast.
    Grepped every `as i32`/`as u32`/`as usize`/`to_int_unchecked` in
    `src/lib.rs`: the only float-sourced one is `cbrt`'s
    `koff as i32`, where `koff` is a select between two constants, so
    LLVM folds it to a select between two integer constants and no
    runtime convert exists. Everything else is int->int or int->float
    (`e as f32`, `vcvtdq2ps`), neither of which a range hint affects.
    Clamps on plain FP values are not removable by range information —
    the clamp *is* the operation.

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

49. **sinf_poly real-chain refit** (#1's method) scoring sin_checked +
    cos_checked's actual reductions jointly — the rejected LPs used
    continuous grids that mis-weighted the caller split.

50b. **Single-reduction `tan`/`tan_checked`** — noticed 2026-07-27 while
    looking for throughput targets (`tan_checked` is the crate's 2nd most
    expensive function at 8.235 cyc/elem), **pre-screened on paper and
    parked as an accuracy-for-speed trade, not a free win.** Both `tan`
    and `tan_checked` are `sin(x)/cos(x)`, i.e. *two* full argument
    reductions. Only one is needed in principle: `tan` has period `pi`, so
    with a single `r = x - q*pi` the `(-1)^q` parities cancel in the
    quotient and `tan(x) = sin(r)/cos(r)` exactly — that would drop one
    whole `round_x_over_pi` + `reduce_pi` pair *and* both `parity` calls.
    - Why it is not free: `cos(r)` for `r` near `+-pi/2` is a
      cancellation. Today's form reduces the cosine around cosine's *own*
      zeros, so `cos_checked` keeps small *relative* error right at a pole
      of `tan` and the quotient's relative error stays bounded. A poly in
      `r` has only bounded *absolute* error, so its relative error goes
      like `1/cos(r)` and blows up exactly at the poles.
    - The obvious rescue does not work either: `cos(x) = +-sin(r - pi/2)`
      and that subtraction *is* exact by Sterbenz (both operands ~1.57),
      but `r` itself only carries ~1 ulp of absolute accuracy, so
      `r - pi/2` near zero is all error. The information is genuinely gone
      from a single-f32 `r` — recovering it needs a double-float `r`,
      which is the `sincos_checked` shared-reduction fusion that already
      measured *slower* on wall-clock.
    - So the live version of this is a deliberately-sloppier tier
      (`tan_fast`, poles excluded from its contract), not a change to
      either existing function. Measure the win before building the API:
      the saving is one reduction out of two, so expect roughly -40%, but
      it buys nothing the existing `tan` contract can keep.

54i. **The licence can be *downstream*, not just upstream — and that
    reopened one of the entries above.** `compound` was screened closed
    on the grounds that its `log1p` argument is unconstrained. True, and
    irrelevant: `compound` is `exp_checked(n * log1p(x))`, and
    `exp_checked` maps **both** signed zeros to exactly `1.0`, so
    `log1p`'s trailing signed-zero select cannot be observed no matter
    what reaches it. `log1p_nonzero!` there is bit-identical and
    `compound_throughput` drops 159 -> 155.
    - Add to the sweep recipe: the existing question is "what does the
      caller's guard prove about the argument?" The missing one is
      **"what does the consumer discard about the result?"** A function
      whose output is fed into something many-to-one — `exp` at zero,
      `abs`, a comparison, a saturating clamp — licenses dropping
      whatever distinction that map collapses.

89. **Gather-based LUT tier screen**: a true vgatherdps table + short
    poly for exp2 (16/32-entry exact 2^(j/N) hi words) — distinct from
    the rejected compare-select tree; screen the gather's mca cost
    first; likely lives in the simd/slice tier.
    **Screened and rejected 2026-07-28, in about a minute, by asking
    llvm-mca to price the one instruction the idea rests on** —
    `echo 'vgatherdps %ymm1, (%rax,%ymm2,4), %ymm0' | llvm-mca
    -mcpu=native`. On this target, Block RThroughput for one ymm op:
    `vfmadd231ps` **0.5**, `vmovups` 0.5, `vpermps`/`vpermi2ps` 1.0,
    `vgatherdps` **4.0**.
    - Decisive without implementing anything: one `vgatherdps` costs 4
      cycles per 8 elements, while the *entire* degree-5 Horner it would
      replace is 5 `vfmadd` = **2.5 cycles per 8 elements**. The gather
      alone is 1.6x the whole polynomial, and a LUT scheme still needs a
      degree-2/3 poly on top. Trading 3 fma for one gather is +2.5 cyc
      per 8 elements on a function (`exp2`) costing 0.841 cyc/elem in
      total — a ~37% regression.
    - Reusable: for any "replace arithmetic with a table" idea, price the
      table instruction against the arithmetic *first*. A one-line
      llvm-mca run on a single instruction is the cheapest screen in this
      file.

90. **Same gather screen for log_2** (mantissa-segment table +
    low-degree poly). **Falls with #89 on the same number**: log_2's
    degree-9 Estrin is ~10 ops, about 5 cycles per 8 elements, so even
    replacing *six* of its fma with one gather is a net loss.
    - What this does **not** kill: `vpermi2ps` at RThroughput **1.0**
      holds a 16-entry table in two ymm registers with no memory access
      at all — cheaper than 2 fma. So idea #163's in-register LUT is the
      one table technique that prices well on this CPU, and it stays
      gated on the simd tier only because it is unreachable from scalar
      autovectorized code, not because it is slow. That reorders the two
      LUT ideas' priority.

91. **Slice tier × rejected-global-flags interaction**: per-function
    zmm-width and interleave=2 re-tests become possible via
    #[target_feature] once the slice tier exists — two documented
    median-win/outlier-veto rejections become recoverable.

97. **Targeted LLVM flag screen**: -enable-unroll-and-jam,
    -extra-vectorizer-passes, SLP horizontal reductions — cheap sweep,
    same method as the interleave/zmm experiments. **Done 2026-07-27:
    none of them change this crate's codegen at all.** All five tried —
    `-enable-unroll-and-jam`, `-extra-vectorizer-passes`,
    `-slp-vectorize-hor`, `-vectorizer-maximize-bandwidth`,
    `-enable-masked-interleaved-mem-accesses` — produce **bit-identical**
    `mca_target.s` (383662 instructions each). No llvm-mca run needed:
    identical asm cannot differ in cycles. Plausible reading is that the
    hot loops here are already fully unrolled fixed-size array loops over
    straight-line FP, which gives these passes nothing to find.
    - **Screened by asm diff rather than by mca**, which is the cheap way
      to do this kind of sweep: a flag that leaves the asm untouched is
      disqualified in ~15s instead of ~5min.
    - **Two harness traps, and the first invalidated two whole attempts
      before a control caught it.** Recording them because any future flag
      or codegen sweep in this repo will hit both:
      1. **`RUSTFLAGS` (env) *replaces* `.cargo/config.toml`'s
         `build.rustflags`, it does not merge.** This repo's config sets
         `-C target-cpu=native`, and dropping it makes `src/lib.rs`'s own
         FMA `compile_error!` fire — so every flag build simply *failed*,
         left the previous `.s` in place, and looked like "identical asm,
         no effect". Always pass `RUSTFLAGS="-C target-cpu=native <flag>"`
         and do not redirect stderr away.
      2. Cargo's fingerprint does not track `--emit=asm` (mca.rs documents
         this and bumps the mtime for exactly this reason), and a changed
         `RUSTFLAGS` produces a *new* metadata hash, so several
         `mca_target-*.s` files coexist — select the newest by mtime, never
         `ls | head -1`.
    - **Always include a control flag.** `-force-vector-width=2` must
      change the asm; when it too came back "identical", that proved the
      mechanism was broken rather than the flags being inert. Once fixed,
      the control correctly showed 383662 -> 393186 instructions (narrower
      vectors need more of them). Without that control this entry would
      read "no LLVM flag helps" on the strength of builds that never ran.

125. **Integer-domain parity pipeline end-to-end** for
     sin_checked/cos_checked (parities as bits, XOR combine, direct
     sign mask) — composes #45/#46; deletes the float compare+select
     flip. **Dead 2026-07-28: both components it composes are now
     rejected.** #45 was already rejected (`rem` is unbounded, real
     wrong-sign catastrophe), and #46 is rejected above on a derivation
     — 7 integer ops against 3 FP ops, plus the magic-round shortcut is
     invalid exactly inside `sin_checked`'s accurate range. The third
     piece, the XOR combine, measures as a wash on its own. Nothing left
     to compose.

157. **remainder_checked/remainder_wide consolidation screen** after the
     ties-even fix (#80): can one tier serve both contracts, or does
     wide's 6.5x cost keep them split? Screen only. **Answered
     2026-07-28: the cost keeps them split, and the current table makes
     that unambiguous.** `remainder_checked` 1.544 cyc/elem against
     `remainder_wide` 8.186 — a **5.3x** gap, and 142 instructions
     against 47. Consolidating upward makes every `remainder_checked`
     caller pay 5.3x for a `|x/y| > 2^24` contract almost none of them
     need; consolidating downward silently drops that contract. This is
     precisely the fast/checked/wide split the crate uses everywhere
     else (`sin`/`sin_checked`, `exp2`/`exp2_checked`), so the answer is
     "keep them split" for the same reason those stay split.
     - The interesting question the screen *does* surface is a different
       one, not in this file: `remainder_wide` is the crate's 3rd most
       expensive function and its `Df32` chain has never been attacked
       on its own terms. That is a live target, unlike consolidation.

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

189. **periodic_poly! dedup macro** once #43/#44 land (four
     near-identical folded-constant sinf_poly variants) — same
     macro-not-fn pattern as pi_reduce_and_poly!.

## Tried and shipped

Kept for the mechanism, not the credit: each of these records *why* a change
worked, which is what makes the next one findable. The code itself is in
`src/lib.rs` and the git log.

- **Automated evaluation-order search per poly** — **swept 2026-07-28;
  4 of 6 sites shipped**, and the useful part is that the winning axis
  was not the one this entry names. The Horner<->Estrin axis really is
  exhausted (`acos_poly`, `atan_poly`, `erfc_rational`, `log_2`,
  `asin_small`/`sinh_small` all previously rejected). The unswept axis is
  the one `ln_normal` and `exp2_q_poly!` already use and nobody
  generalized: **fold the top coefficient group in one level lower so
  `x^4` is never formed at all.**
  - `l0 + l1*x2 + l2*x4` becomes `l0 + (l1 + l2*x2)*x2`. Algebraically
    identical, **same fma critical-path depth**, one plain multiply
    fewer — *and* one rounding fewer, since `fl(x2*x2)` no longer exists.
    So unlike a Horner/Estrin swap this is expected to help accuracy
    rather than trade it, which is exactly what happened.
  - Shipped at 4 sites: `exp_r_poly!` (20 callers),
    `exp_pos_neg_core!` + `exp_pos_neg_narrow_half` (7),
    `dawson_central_ratio` + `dawson_tail_poly`, and `tan_poly`.
    **Whole file 360429 -> 358643 instructions; 29 throughput regions
    shrink and not one grows.** `exp` 62->60, `exp_checked` 66->64,
    `expm1` 86->84, `tanpi` 89->85, `dawson` 82->78.
  - Accuracy improves across the board: `cosh`/`cosh_narrow` avg
    **-11.9%**, `cosh_checked` **-11.7%**, `exp`/`exp_narrow`/
    `exp_checked` **-5.2%**, `sinh` -4.5%, `coshm1` -3.7%, `logaddexp`
    -3.0%. Max ulp *drops* on eleven functions: `expm1` and
    `exp_m1_over_x` 6->5, `sinh`/`cosh`/`cosh_checked` 5->4, `sigmoid`
    and `softplus`/`logsigmoid` 4->3, `coshm1` 12->9, `tanpi`/`tan2pi`
    11->10. Exhaustive where it mattered: `norm_pdf` 0.0881/67 ->
    **0.0876/66**. Nothing regresses. 215 `worst_corpus` entries move,
    every one inside the changed families. All 9 gates, 27 tests pass.
  - **Two sites rejected, both mechanistically, which is what turns this
    into a rule.** `erf_poly` is degree 6 with an odd leading
    coefficient, so the fold cannot stay flat — it costs **+1
    critical-path level**, and `erf_throughput` is 68-79%
    dependency-bound, so it shows up directly (+3.2%). `erfinv` is the
    one region in the sweep that is genuinely *resource*-bound (86%
    resource pressure): dropping `u4*u4` freed registers, LLVM responded
    by folding 22 `{1to8}` broadcast-memory operands into the fmas, and
    **instructions fell 177->164 while uops rose 19400->20300** — fewer
    instructions, more uops, real regression (`probit` +21.8%,
    `erfc_inv` +25.3%).
  - **Rule: the fold wins iff it keeps critical-path depth AND does not
    raise uop count.** Instruction count alone is not sufficient, and
    `erfinv` is the counterexample that proves it.
  - **A fifth site, missed by the sweep and found by checking its
    coverage — and it is the most interesting one.** `tanh` carries a
    *standalone* rescaled copy of `exp_r_poly` (coefficients scaled by
    `2^degree` for `rh = r/2`), so a grep for the macro does not find it;
    it still formed `rh4` for a single use. Folding it:
    `tanh_throughput` **79 -> 77** (exactly -2 `vmulps`), latency region
    -64, nothing else moves, exactly 1 `worst_corpus` entry moves.
    - **Exhaustive: `tanh` avg 0.1457 -> 0.1452, max ulp 6 -> 5.**
    - That is worth noticing against idea #169, which studied `tanh`'s
      max-6 band at length and closed all three of its levers — a wider
      Pade fit (disqualified: idealized error alone would be 42 ulp), a
      dedicated third branch (measured +39.6% throughput), and a seam
      retune (the seam is already optimal). All three tried to *fit* the
      error better. This removes a rounding instead, and gets the 6 -> 5
      that none of them could, while making the function cheaper.
    - Reusable: when a poly is copied rather than shared, a
      macro/function-name grep will miss it. Sweep by *shape* — here,
      "forms a fourth power and uses it once" — not by call site. Same
      lesson as #54h's reachability sweep, one abstraction over.
  - Not opportunities, checked and ruled out: the degree-8/9 Estrin
    chains (`log_family_normal!`, `ln_normal`, `log1p_unit`,
    `erfinv_central_poly`) use `s4` **twice**, so it is genuinely shared
    — those already apply this same trick one level up, which is where
    the idea came from.
  - **Methodology, and this is the sharpest instance of it in the file:
    `cos2pi`'s avg ulp moved 0.0737 -> 0.1172 (+59%) and its max moved
    51472 -> 3294199 (64x) on BYTE-IDENTICAL assembly.** Same for
    `norm_cdf`'s max (264 -> 276, true exhaustive value 295 both sides)
    and `gelu`'s avg. The quick fuzz is unseeded, and for a
    cancellation-dominated metric a handful of near-zero samples landing
    differently moves the mean and the max by orders of magnitude. Three
    of four apparent regressions evaporated on one check: **diff the
    function's assembly before interpreting any accuracy delta.** If the
    region did not change, the delta is not real, full stop.
  - Caveat on `dawson`'s numbers, mine and anyone else's: the harness
    caps that group at 2M samples even in `thorough` mode (no sleef
    bucket — the reference is a hand-rolled Simpson's quadrature), so
    its avg/max are sampled, not exhaustive. Measured 0.8076/60 ->
    0.8070/60; treat +-0.01% avg and +-1 max there as no signal and rest
    the site on the structural argument (strictly one rounding fewer).

- **Crate-wide poly headroom screen — run 2026-07-28, and this table is
  the answer to "is there any fit left anywhere".** Metric: exact-
  arithmetic error of the shipped f32 coefficients, weighted by the final
  combine's own sensitivity `|d(result)/dP| / ulp(result)`, against the
  same-degree ulp-weighted LP minimax optimum **re-quantised to f32**.
  Ratio = headroom. Idea #1's `acos_poly` result generalised into a
  screen that costs seconds per poly.

  | poly | callers | cur | LP(f32) | ratio | outcome |
  |---|---|---|---|---|---|
  | `atan_poly` | atan, atan2, atand, atanpi | 0.228 | 0.016 | 14.2x | untried; abs error already << 1 ulp |
  | `asin_poly` | asin, asind | 1.869 | 0.191 | 9.8x | **REJECT** — headroom is fake, see below |
  | `asin_small` | asin, asind | 1.104 | 0.117 | 9.4x | REJECT (real chain) |
  | `asinpi_poly` | asinpi | 0.781 | 0.088 | 8.9x | **SHIPPED** |
  | `dawson_central_ratio` | dawson | 60.1 | 8.35 | 7.2x | REJECT (axis trade) |
  | `erfc_rational` | erfc, erfcx, norm_cdf | 15.1 | 2.70 | 5.6x | REJECT (axis trade) |
  | `cbrt_corr` | cbrt tiers, rcbrt | 2.99 | 0.789 | 3.8x | REJECT (axis trade) |
  | `erfinv_tail_poly` | erfinv, erfc_inv, probit | 28.6 | 8.44 | 3.4x | **SHIPPED** |
  | `asinpi_small` | asinpi | 1.882 | 0.676 | 2.8x | **SHIPPED** |
  | `ln_normal` | ln, log1p, asinh, acosh, atanh, erfinv, logit, clog, xlogy | 1.101 | 0.389 | 2.8x | **REJECT — fake headroom, see below** |
  | `tan_poly` | tanpi, tan2pi | 8.87 | 3.55 | 2.5x | **SHIPPED** |
  | `acos_poly` | acos, acosd | 3.185 | 1.584 | 2.0x | shipped (idea #1) |
  | `acospi_poly` | acospi | 3.001 | 1.614 | 1.9x | **SHIPPED** |
  | `log10_normal` | log10, log10p1 | 0.931 | 0.760 | 1.2x | no headroom |
  | `exp2_q_poly_centered` | exp10_checked | 0.140 | 0.133 | 1.05x | **exhausted** |
  | `erf_poly` | erf | 0.831 | 0.819 | 1.01x | **exhausted** |
  | `exp_r_poly` | exp, expm1, tanh, sigmoid, sinh, cosh... | 1.354 | 1.353 | 1.00x | **exhausted** (c0/c1 pinned to 1) |
  | `erfinv_central_poly` | erfinv, erfc_inv, probit | 0.503 | 0.503 | 1.00x | **exhausted** (floor is f32(sqrt(pi)/2)) |
  | `erf_pade` | erf | 0.871 | 0.891 | 0.98x | **exhausted** |
  | `log1pmx_Q` | log1pmx | 0.337 | 0.351 | 0.96x | **exhausted** |
  | `LOG2_COEFFS` | log_2, log2_df, log2p1 | 0.287 | 0.307 | 0.93x | **exhausted** |
  | `sinf_poly` | sin, cos, sind, cospi, tan, sinc... | 0.166 | 0.206 | 0.81x | **exhausted — closes idea #49** |
  | `log1p_unit_Q` | asinh, acosh | 0.072 | 0.360 | 0.20x | **exhausted** |

  - **Ratio < 1 means the shipped coefficients already beat a freshly
    quantised LP optimum** — the fingerprint of prior real-chain or
    coordinate-descent tuning. Nine polys are now closed on a number.
  - **The ratio predicts headroom; a second question predicts whether it
    converts — what objective produced the shipped coefficients.**
    - *real-chain / coordinate-descent tuned* (`asin_poly`, `sinf_poly`,
      `log1pmx`, `LOG2_COEFFS`): ratio <= 1, nothing to take.
    - *least-squares tuned* (`cbrt_corr`, `dawson_central_ratio`,
      `erfc_rational`): a minimax refit **always** trades avg for max —
      3 of 3 here, no exceptions. `cbrt` gets max 3->2 for avg
      0.2813 -> 0.3765 (+34%); `erfc` gets avg 0.3055 -> 0.2258 for max
      109 -> **121**; `dawson` max 61 -> 13 for avg 0.806 -> 2.993.
    - *a rescale of another poly's coefficients* (`acospi_poly`,
      `asinpi_poly` — both were `acos_poly`/`asin_poly` divided by pi),
      or an LS fit whose error concentrates where the weight is high but
      the sample density is low (`erfinv_tail_poly`, `tan_poly`):
      **wins on both axes.** All four ships are in this class.
  - **The cheapest screen of all, and it should be run FIRST — before
    the ratio, before the objective class.** The ratio measures fit
    headroom, but headroom only converts if **the fit is the binding
    term**. Test that directly: score a *known-better-fitting* poly
    through the real chain — a higher degree, or the pre-shed ancestor
    sitting in git history — and see whether the real error moves. One
    run. It settled `ln_normal` outright (degree 9 is a 22x better fit
    and measures *worse*), and the same shape of experiment is what
    exposed `asin_poly`'s fake 9.8x. The strongest version adds an
    **oracle row**: feed the chain a correctly-rounded `p` and see where
    the floor actually is. For `ln_normal` that floor is max 1 / avg
    0.272 against a shipped max 3 / avg 0.453 — i.e. two of the three
    ulps were never the polynomial's to give.
  - **`asin_poly`'s 9.8x is the table's largest headroom and is entirely
    fake** — worth knowing before anyone trusts the ratio alone. Real-
    chain simulation over all 16106127 f32 in `[0.27, 1)`: shipped max 6
    / avg 0.87462, ulp-weighted minimax max **8** / avg 1.17783.
    `asin`'s exhaustive worst case is `x = 0.27004012`, *exactly* the
    branch crossover, where `pi/2 - sqrt(1-a)*P(a)` cancels and
    `ulp(1.2975)/ulp(0.2734)` amplifies the product's own rounding 4x.
    No coefficient can move a chain floor.
  - **`atan_poly`'s 14.2x is the table's top ratio and is *already
    closed*** — the ratio is a relative measure and here the absolute is
    negligible. Its idealised weighted error is **0.228 ulp**, i.e. the
    entire fit budget is a fifth of an ulp, so recovering 14x of it
    cannot move `atan`'s real max of 4. Idea #8 already quantified this
    from the other direction with a full joint nonlinear refit of all six
    rational coefficients: ~99% of `atan`'s observed error is rounding
    chain. **Read the ratio and the absolute together** — a large ratio
    on a poly whose absolute error is already far under 1 ulp is not a
    lead, and this is the one row in the table where they disagree.
  - **`ln_normal` (2.8x) is closed, and its diagnostic is the sharpest
    in this file: a strictly better fit already exists in this repo's
    own git history, and it measures *worse*.** Real-chain scoring
    (hardware fma, bit-exact) over **all 8388608 mantissas** the
    decomposition can produce, no stride:

    | poly | idealised | real max | real avg |
    |---|---|---|---|
    | oracle: correctly-rounded `ln(1+s)/s` | 0 | **1** | **0.2723** |
    | degree 9 (pre-shed, `9345bd2^`) | ~0.05 | 3 | 0.4644 |
    | degree 8 shipped (LS) | 1.101 | 3 | **0.4531** |
    | degree 8 LP minimax -> f32 | 0.427 | 3 | 0.4670 |
    | degree 8 ulp-weighted LS | 0.939 | 3 | 0.4607 |
    | degree 8 plain LS | 1.122 | 3 | 0.4520 |

    A **22x better fit gives a worse real avg and the same max**, and six
    objectives spanning idealised 0.05 to 1.12 all land within +-3% of
    the same real avg. The oracle floor is max 1 / avg 0.272, so the
    whole gap to max 3 / avg 0.453 is the 8-fma Estrin chain's own
    rounding. Nothing a coefficient can reach.
  - The LP candidate itself, exhaustive on both sides, is the **fourth of
    four** confirmations of the least-squares-tuned rule: `acosh` max
    4 -> 3, bought with a visible avg regression on `ln`, `ln_unchecked`,
    `log1p`, `log1pmx` and `acosh` itself. A real-chain coordinate
    descent (exhaustive enumeration, +-32 ulp span, 8 free coefficients)
    does better — `ln`'s true all-exponent aggregate 0.235184 ->
    0.235132, **-0.022%** — but pushes **`asinh` max 3 -> 4** at
    x = 0.3538082, and the 100M quick fuzz reported max 3 on *both*
    sides. Only `thorough` found it. Not worth 0.022%.
  - Two corrections this produced, both worth keeping. **`softplus`,
    `logsigmoid` and `logaddexp` do not route through `ln_normal`** —
    they use `log1p_unit` (the already-closed 0.20x row), and `log2p1`
    uses `log_2_normal`; the caller list above is corrected. And their
    quick-fuzz avg drifted 0.0752 -> 0.0751 across builds with a
    *provably identical* code path, which incidentally calibrates the
    quick harness's own noise floor at **~+-0.0001 on avg**.
  - **Methodology: the stride-subsample trap applies to *search*, not
    just verification.** A stride-32 real-chain coordinate descent on
    `asin_poly` converged to a candidate that beat the shipped
    coefficients *on its own subsample* (max 6->5, avg 0.930->0.815) and
    lost on the full 16.1M sweep (max 6->**7**, avg 0.875->0.911). Never
    descend on a subsample of a domain you can enumerate — the descent
    finds precisely the points you skipped.

1. **Real-chain refit** — **shipped for `acos_poly` 2026-07-28**, the
   first refit in this file to beat the "rounding-chain-dominated" wall
   rather than be stopped by it. `acos` max ulp **5 -> 4**, avg
   **0.0650 -> 0.0555 (-14.6%)**, `acosd` avg 0.0634 -> 0.0545, for
   provably zero perf change (`acos_throughput`/`acos_latency`/
   `acosd_throughput` have byte-identical mnemonic streams; only the
   constant pool differs). `acospi` is untouched — it has its own poly.
   - **The reusable output is the screen, not the coefficients.**
     Compare a poly's *exact-arithmetic* ulp-weighted error against a
     same-degree weighted-LP optimum: seconds of scipy, and it cleanly
     separates "no headroom" from real headroom. `acos_poly`: shipped
     coefficients in exact f64 = **3.18** ulp-equivalent, chain floor
     (same f32 chain fed an ideal poly value) = **2**, measured = 5, LP
     optimum = **1.49** (1.54 quantised) — **2.1x real headroom**.
     `log_2`, run through the same screen: current **0.3092**, LP optimum
     **0.3092**, i.e. bit-for-bit at minimax already, its peak pinned at
     `s -> 0` by `LOG2_E`'s own f32 rounding. So this result does *not*
     overturn the class — it identifies which members of it were
     misdiagnosed.
   - **What was actually wrong was the objective, not the search.**
     `tune.rs`'s `acos_poly` grid walks raw bits in steps of 10000 from
     0 and is **positive-x only**. The fit error peaks as `x -> 1` (3.185
     inside `1-a < 1e-6`) but the grid's largest point is 0.9998 with
     ~1 point inside `1-a < 1.2e-3`, so **the grid's own max is 2.717 —
     it structurally cannot see the binding region** — and it never
     scores the `+PI` negative branch, which has a different ulp scale.
     The probe itself was fine (bit-identical to the shipped chain for
     `x >= 0` over 28M samples); this is not the `erf_poly` probe bug.
   - The weight that made "minimax" mean the right thing:
     `sqrt(1-x)/ulp(acos(x))` against the target `acos(x)/sqrt(1-x)`.
     The LP seed delivered ~99% of the win; the finite-difference descent
     through the real chain added ~1% (avg 0.112046 -> 0.111928). Worth
     recording that for this poly the idealized and real-chain optima
     essentially coincide *once the weight is right* — the expensive part
     of #1 was not what paid.
   - **Idea #3 (1-D +-few-hundred-ulp scan of a combine constant) is
     structurally dead here**: c[0..2] swept at +-64 ulp, stride 1, zero
     hits. The trailing constant is *pinned* by `edgecheck`
     (`acos(0) == FRAC_PI_2`), and 365 ulp of c[0] moves the result by
     1 ulp — the winning move was ~220000 ulp of c[0]. A local scan was
     never going to reach it.
   - **Max-claim trap, third independent sighting today**: a stride-8
     subsample of the full domain (266M points, 2500x denser than
     tune.rs's grid) still reported max 4 for candidates whose true max
     is 5. Avg was faithful; max was not. Only `bstep=1` settles a max.
   - Minimax redistributes, and honestly: big wins in `|x|` in
     `[1e-3,0.1]`, `[0.5,0.7]`, `[0.99,1]` (positive `[0.999,1]`: avg
     2.150 -> 1.263, max 5 -> 3), real regressions in `[0.1,0.5]` and
     `[0.7,0.9]`. Global avg and global max both improve. The 16
     `worst_corpus` entries that move (8 `acos`, 8 `acosd`) net *+3 ulp
     worse* because those curated inputs cluster in `|x|` in
     `[0.1,0.25]`, exactly the band a minimax refit raises — expected,
     not a regression signal. Blessed.
   - **Open lead, quantified by the same screen**: `acospi_poly` shows
     the identical signature — idealized **3.00** vs LP minimax **1.61**
     (1.62 quantised), 1.9x headroom, same construction, currently max 5
     / avg 0.0536.
   - Cost: ~1.1e12 scored real-fuzz evaluations (~460 exhaustive passes
     plus 401 stride-8 screens), ~2.5 h wall.

8. **atan_poly joint numerator+denominator nonlinear refit** (scipy
   least_squares on the true rational) — only separate num-only/
   denom-only LPs were tried; the max-4 worst point was diagnosed as
   denominator-or-division-bound. **Screened 2026-07-27 and closed: there
   is no fit to recover.** A joint nonlinear minimax over all six free
   coefficients (`a0..a2`, `b0..b2` of the `[3/3]` in `t = x^2`, both
   constant terms pinned to 1 as the shipped construction does), seeded
   from the shipped values and run to convergence:

   | | idealized max rel err | ulp-equivalent |
   |---|---|---|
   | shipped, as f64 literals | 1.949e-09 | 0.033 |
   | best joint refit, f64 | 1.567e-09 | **0.026** |
   | that refit rounded to f32 | 8.62e-09 | 0.145 |
   | shipped, as actually rounded to f32 | 1.73e-08 | 0.291 |

   - The whole idealized error budget is **0.03 ulp** and the best
     possible joint refit recovers 1.24x of it. Against `atan`'s real max
     of 4 ulp that is nothing: ~99% of the observed error is rounding
     chain, which is precisely the "denominator-or-division-bound"
     diagnosis the entry already carried, now quantified.
   - Note the fourth row is *worse* than the second-best row, and that is
     not a defect: the shipped coefficients were tuned by `tune.rs`
     against the **real f32 chain**, not against an idealized model, so
     they deliberately sit off the idealized optimum. Any future refit
     here has to be scored the same way to be comparable — an idealized
     0.29 -> 0.145 "improvement" would very likely be a real regression.
   - Methodology trap worth recording, since the same LP is the obvious
     tool for every rational in this crate: the natural
     **`linprog`-feasibility-plus-bisection** formulation of rational
     minimax **silently fails at this error scale**. HiGHS's default
     primal feasibility tolerance is ~1e-7, and the constraints here need
     to hold to ~1e-9 in quantities of size O(1), so it reported
     "infeasible" at tolerances the *shipped* coefficients demonstrably
     satisfy (verified directly: max constraint violation -5.3e-11, i.e.
     strictly feasible, while HiGHS returned status 2). Left uncaught it
     produces a confident, completely wrong answer — the first run of
     this screen "found" an optimum 25x *worse* than the shipped point.
     Either rescale the problem around the incumbent (substitute
     `a_i = a_i^ship + 1e-8*alpha_i` so residuals come out O(1)) or, as
     here, skip the LP and run a derivative-free nonlinear minimax
     seeded from the shipped values. Always sanity-check that the
     incumbent is reported feasible before trusting any LP verdict.

12b. **The `<< 23` half of #12, done right — shipped 2026-07-28, and the
    win is the *mask*, not the ports.** #12 framed this as moving work
    onto less-contended integer ports, and that framing is what made it
    look dead after #14 measured the port premise failing. The actual
    saving has nothing to do with ports: it is that the `<< 8` form
    **needs a mask and the `<< 23` form does not**.
    - Why: `(k + 383.0)` lands in the `2^8` binade, whose biased exponent
      field is 135 — an *odd* number — so `<< 8` shifts that field's low
      bit straight into the sign bit, and `& EXPONENT_MASK` exists purely
      to clear it. Choose a magic in the `2^23` binade instead (field 150,
      even) and the bit that lands in position 31 is one the construction
      already controls. Folding the `+127` exponent bias into the magic
      constant rather than adding it to the bits is the second half: the
      whole thing becomes `f32::from_bits((k + 12583039.0).to_bits() << 23)`
      — **2 ops against 3**, and no integer add at all, so #14's
      port-pressure finding never comes into play.
    - `12583039.0 = 1.5*2^23 + 127`. For integer `k`, `k + MAGIC` is exact
      in the ulp-1 binade `[2^23, 2^24)` and puts `k + 127` in the low 9
      bits; `<< 23` moves those into the exponent field with a zero sign
      bit and a zero mantissa. Verified bit-identical to the old form for
      **every `k` in [-127, 128]**, which includes `k = 128 -> +inf` (how
      the unchecked tier overflows) and `k = -127 -> +0.0`.
    - That range is the licence, and it was checked by sweeping every f32
      in each site's own documented domain rather than argued: `exp2`
      `[-126,128]`, `exp10` `[-127,127]`, the `exp_narrow` family
      `[-126,127]`, `exp_pos_neg_narrow_half` `[-126,126]`, and the two
      *clamped* sites `tanh` `[-126,127]` and `sigmoid` `[-126,128]`
      (total, so those two cannot leave it for any input at all).
    - 9 call sites, deduped into one `exp2int_field!` macro (macro, not a
      fn — the +32% shared-fn-boundary precedent). **14 public functions
      get smaller** with no other region moving: `exp2` **41 -> 38**
      instructions (-7.3%), `sigmoid` 55 -> 51 (-7.3%), `exp10` 63 -> 60,
      `exp_narrow` 47 -> 44, `tanh` 82 -> 79, plus `expm1_narrow`,
      `exp_m1_over_x_narrow`, `sinh_narrow`, `cosh_narrow`, `silu`,
      `softplus`, `logsigmoid`, `logaddexp`, `erf`. Whole file 363164 ->
      360429.
    - **The cost, stated plainly: out-of-domain garbage changes flavour
      and can now be negative.** Dropping the mask drops the accidental
      guarantee that the field construction was non-negative, so e.g.
      `exp2(-200)` goes from `+7.6e-6` to `-7.2e16`. Both are nonsense —
      the true value is ~6e-61 — but only one of them *looks* like an
      answer, and `powf_unchecked`'s own doc already names "bare `exp2`
      silently wraps around into plausible-looking garbage" as a hazard.
      Taken as an improvement on that axis too, and `exp2`'s doc now says
      so explicitly, because `exp2(x) >= 0` was a real (if unpromised)
      invariant that no longer holds off-domain.
    - `worst_corpus` moves **135 entries, and the pattern is the proof**:
      every one is in `exp2`/`exp10`/`exp_narrow`/`expm1_narrow`/
      `exp_m1_over_x_narrow`/`sinh_narrow`/`cosh_narrow` — the seven
      *unclamped* tiers — at an input outside that tier's domain (checked
      individually, not sampled). **Zero** entries move in `tanh`,
      `sigmoid`, `silu`, `softplus`, `logsigmoid`, `logaddexp` or `erf`,
      the total/clamped functions that share the same construction. Blessed.
      All 10 gates pass.
    - Generalizable, and it is the reason this sat unfound: #12's own
      framing named the wrong mechanism, so the measurement that killed
      the named mechanism (#14, ports) was allowed to kill the idea. When
      a backlog entry proposes *a transform* plus *a reason it should
      win*, disproving the reason does not disprove the transform.

12c. **`exp2_field_split` via the same maskless `<< 23` form — built,
    measured, and then REVERTED on a contract question, not a
    measurement.** This is the largest single perf lever found in this
    sweep and it is sitting there pre-verified; it needs one judgement
    call, not more work.
    - Construction: `half = fma(k, 0.5, EXP2INT_MAGIC)`, `lo = from_bits(
      half.to_bits() << 23)`, `hi = exp2int_field!(k - (half - MAGIC))`,
      return `(hi, lo)`. **6 ops against the shipped form's 8.**
    - Measured: **23 throughput regions shrink, none grow**, whole file
      356675 -> 353982. `exp` **60 -> 55 (-8.3%)**, `exp_checked` 64 ->
      59, `cosh`/`cosh_checked` -5, `sinh`/`sinh_checked` -5, `expm1` -5,
      `powf_checked` 251 -> 245, `norm_pdf`/`tanh_grad`/`sigmoid_grad`
      -5 each, `erfc_accurate`/`erfcx_checked` -5, `compound` -4.
    - **In-domain accuracy is untouched** — verified two ways, and both
      agree. (a) The two constructions are bit-identical for every `k` in
      `[-200, 200)`, which covers every in-domain `k` for every caller
      (`exp2_checked` clamps to `[-151, 128]`). (b) The real fuzz over
      16 function groups comes back digit-for-digit identical (`exp`
      0.0707/3, `exp_checked` 0.0369/3, `expm1` 0.1291/5, ...).
      `saturation_pins`, `denormal_audit`, `edgecheck`, `special_matrix`
      and 27 tests all pass.
    - **Why it was reverted anyway.** The magic that carries the `+127`
      bias must be `≡ 127 (mod 512)`, hence **odd**, which flips the
      round-half-to-even tie-break, so for odd `k` the two halves come
      out exchanged. Returning them swapped restores the original pairing
      exactly *within* `[-200, 200)` — but not outside it, and the
      unchecked tiers do not clamp. Concretely `cosh(2048)` goes from
      `+inf` to `1.0666397e35`. That input is far outside `cosh`'s
      documented domain and contractually garbage either way, but the old
      value happened to be the *mathematically correct* one, and the new
      one is a finite plausible-looking number. That is the opposite
      direction from #12b's own justification (there, garbage became
      *more* obviously garbage), and this crate has a consistent history
      of paying to avoid plausible-looking wrong output.
    - 55 `worst_corpus` entries move, **all out-of-domain, all on
      unchecked tiers** (`cosh`, `sinh`, `exp`, `expm1`,
      `exp_m1_over_x`, and the two `_throughput` siblings). Zero move on
      any `_checked` tier, because those clamp `k` into the range where
      the two forms agree.
    - **The variant that would be free, for whoever takes this up**:
      apply the new split only where the caller already clamps `k`, via a
      `const` generic (`exp2_field_split<const CLAMPED: bool>`, the same
      idiom `round_x_over_pi::<HALF>` already uses). That captures ~15 of
      the 23 regions — `exp_checked`, `exp2_checked`, `exp10_checked`,
      `sinh_checked`, `cosh_checked`, `powf_checked`, `erfc_accurate`,
      `erfcx_*`, `norm_pdf`, `tanh_grad`, `sigmoid_grad`, `compound` —
      with *provably zero* behavioural change anywhere, and leaves
      `exp`/`sinh`/`cosh`/`expm1` on the current form. The cost is two
      instantiations of a subtle bit-trick.
    - Methodology note on my own error, worth keeping: I verified the
      swap over `k in [-200, 200)` and called it bit-identical. It is —
      but `cosh(2048)` needs `k ~ 2954`, and the unchecked tiers reach
      there. **State the range a bit-trick was verified over, and check
      that every caller is actually inside it**, especially when the
      callers include tiers that deliberately do not clamp.

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
    **Rejected on the op count 2026-07-28, without building it.**
    `denormal_rescale!` is already **4 ops** and provably minimal:
    1 compare, 1 multiply, 2 selects, and every one of them is load-bearing
    (the `koff` select cannot become a post-hoc `r - 24.0` because `koff`
    folds into the *pre-combine* integer exponent `k`, where it is exact,
    whereas subtracting 24 from an already-rounded log result is not).
    An exact-shift normalize needs **lzcnt + a shift-count subtract +
    the variable shift + a mask + an exponent-field OR, and still a
    compare and selects to gate it against normal inputs** — 7-8 ops,
    strictly more. The only way that loses-on-op-count could still win is
    the idea #12 port-migration premise, and #14 **measured** that premise
    failing on a closely related construction (+13.2% on `exp2_checked`,
    because the int ops contend for the same ports as the FP ops they
    displace). Two independent reasons, no measurement needed.
    - Note this is *not* the same conclusion as #12b, which shipped: that
      one wins by **deleting a mask**, not by moving work to integer
      ports, and it needs no lzcnt.

22. **Toolchain-bump re-screen list**: tag the rejections that were pure
    scheduling artifacts (pre_offset dead-add removal, ln/log10
    trailing-fma fuse +1cyc, reduce_pi depth-2 rebalance) and re-measure
    after each nightly bump — these can silently flip. **First pass run
    2026-07-27 on rustc 1.98.0-nightly (f46ec5218): 1 of 3 flipped, and
    it flipped hard.** The ln/log10 trailing-fma fuse is now a win on
    every affected function — shipped for `ln_normal`, see its own entry
    above for the numbers and for why `log10_normal`'s half stayed
    rejected. So the premise holds: a rejection recorded as "the
    scheduler picked worse" has a real chance of being wrong on the next
    toolchain, and re-screening it costs one mca run.
    - **Second pass, same day: the `pre_offset` dead-add removal flipped
      too** (`tan_checked` -6.3%, `sin_checked` -1.5%; see its own entry
      above). That is 2 of 3 reversed on one toolchain, which upgrades
      this idea from "worth a look after a bump" to a standing chore.
      `reduce_pi`'s depth-2 rebalance was the one that held: it
      reproduces its original +3 cycles exactly, and the reason is
      structural rather than schedule-luck (see its entry) — a useful
      contrast, since it means "scheduling artifact" was the right tag
      for two of the three and the wrong tag for the third.
    - Lesson the flip actually turned on, which generalizes past
      toolchain bumps: the *association* of an fma fold is a separate
      degree of freedom from the fold itself, and both orders need
      measuring. Only one of the two orders wins here; the other still
      reproduces the original rejection's +1 cycle on today's toolchain,
      so the old entry was probably never wrong about what it measured —
      just about which of the two candidate expressions it measured.

54d. **`1+x` is never a denormal, so `log1p`/`log2p1`/`log10p1` don't
    need the log wrapper's denormal rescale** — the second hit from
    #54b's crate-wide sweep, and the widest one: **shipped 2026-07-27**
    (see lib.rs/git log).
    - `fl(1+x)` is exact by Sterbenz once `x <= -0.5`, so the smallest
      positive value it can take is exactly `2^-24` — about 10^30 above
      `f32::MIN_POSITIVE`. Verified exhaustively over all 2^32 patterns
      (zero positive-denormal cases, smallest positive `u` = 2^-24)
      rather than argued from the bound, since this is the licence the
      whole change rests on. The other three arms stay live and are
      shared verbatim through a new `log_family_edges!`: `u == 0` at
      `x == -1`, `u < 0` below it, `!(u < inf)` for `+inf`/NaN.
    - Implemented as `log_family_wrapper_no_denormal!`, a sibling of
      `log_family_wrapper!` over the shared edge macro — no duplicated
      edge logic, and the `_normal` cores are untouched.
    - **Bit-identical over all 2^32 inputs** for all three functions
      (checked directly against the old formulation, not inferred), so
      every downstream caller is too. All 8 gates + `codegen_check` pass.
    - mca throughput: `log1p` 2.289 -> **1.857 (-18.9%)**, `log2p1`
      2.328 -> **1.886 (-19.0%)**, and free for every caller —
      `xlog1py` 2.324 -> 1.857 (-20.1%), `probit` 5.021 -> 4.215
      (-16.1%), `compound` 4.570 -> 4.082 (-10.7%), `erfinv` 4.185 ->
      3.815 (-8.8%), `logit` 3.728 -> 3.457 (-7.3%), `erfc_inv` 4.372
      -> 4.129 (-5.6%).
    - **mca's latency column disagrees for exactly three of those
      (`erfinv` +13%, `probit` +14%, `erfc_inv` +13%) and is wrong
      there** — worth recording as a harness trap, since the arbitration
      generalizes. Those regions' instruction counts went *down*
      (`erfinv_latency` 6662 -> 6020, `probit_latency` 7046 -> 6340,
      `erfc_inv_latency` 6978 -> 6210) and the diff includes 64 `jbe`
      and 64 `jmp` — one per chain iteration. The scalar latency harness
      lowers `denormal_rescale!`'s select to a *branch*, and llvm-mca has
      no branch predictor, so the baseline's simulated trace skipped a
      path the branchless version has to count. Wall-clock settles it:
      3 interleaved A/B reps, new is equal-or-better every time (22.12
      vs 23.70, 21.22 vs 21.32, 20.93 vs 21.31 ns). Every *throughput*
      region (the vectorized one, where the select stays a blend) shrank
      too, and those numbers are the trustworthy ones.

54f. **erfinv/logit: the same `log1p` inlining, one step further** —
    fourth hit from #54b's sweep, **shipped 2026-07-27**. Both call
    `log1p` on an argument whose sign is fixed by construction
    (`erfinv`: `-x*x <= 0`; `logit`: `-p`), so `t = 1+arg` is never
    denormal by #54d's own bound and `arg == 0` only where the *other*
    branch is selected — `log1p`'s signed-zero select and `ln`'s rescale
    both go.
    - `erfinv` sheds one arm more than `logit` does: `t == 0` happens
      **exactly** when `|x| == 1` (verified exhaustively, not argued from
      `fl(x*x)`'s spacing), and the function already carries an explicit
      `x.abs() == 1.0 -> +-inf` override for the Estrin-overflow reason in
      its doc comment — so the wrapper's `-inf` arm is dead too.
      `logit` must keep it: `t == 0` at `p == 1` is what makes
      `logit(1) == +inf`.
    - **Bit-identical on all 4278190083 non-NaN patterns.** `erfinv`'s
      `t <= 0` and NaN arms then fold into a *single* `t > 0.0` select
      (false for both, and both want `NaN`), which canonicalises the NaN
      it returns instead of forwarding the input's payload — the only
      difference, on exactly the 16777213 NaN patterns. Taken, per
      Jodie's rule that a NaN is a NaN (see the `ulp_diff`/`worst_corpus`
      entries); `nan_payload.rs` now reports `erfinv` as canonicalising,
      which is simply the new truth, not a regression.
    - mca throughput, cumulative over #54d + #54f: `probit` 5.021 ->
      **3.514 (-30.0%)**, `erfinv` 4.185 -> **3.321 (-20.6%)**,
      `erfc_inv` 4.372 -> **3.514 (-19.6%)**, `logit` 3.728 ->
      **3.354 (-10.0%)**. The single-select collapse alone is worth
      `erfinv` -3.4%, `probit` -2.1%, `erfc_inv` -0.5%.
    - Latency, same cumulative span: `erfinv` 98.22 -> **43.68 (-55.5%)**,
      `probit` 108.33 -> **55.72**, `erfc_inv` 104.19 -> **52.57** — so
      #54d's apparent latency "regression" on exactly these three ends up
      not merely reversed but less than half its own starting point,
      further confirmation it was the branch-modelling artifact #54d
      describes. `logit`'s row did the same thing at the #54f step
      (56.00 -> 60.94) with its region's instruction count going *down*
      (5657 -> 5645); 3 interleaved wall-clock reps put the new code far
      ahead (12.93/14.60/15.67 vs 19.94/18.61/31.06 ns).

54g. **clog: both `log1p` calls are fully guarded** — fifth and last hit
    from #54b's sweep, **shipped 2026-07-27**, and the most complete one:
    here even `log1p`'s *own* two selects die, not just `ln`'s wrapper.
    `|mag-1| < 0.5` puts the near-1 branch's `u = 1+v` in `(0.5, 1.5)`
    and `ratio in [0,1]` puts the rescaled branch's in `[1, 2]`, so `c/u`
    can't be non-finite and the `v == 0.0` signed-zero guard has nothing
    to fix. What survives is the Sterbenz correction plus `ln_normal`.
    - Bit-identical to `log1p` over **every** f32 in both licensed ranges,
      with exactly one exception found by the check rather than reasoned
      about: `v == -0.0`, where `log1p` returns `-0.0` and the stripped
      form returns `+0.0`. Neither call site can produce it, verified
      exhaustively too — `mag - 1.0` is `+0.0` for every non-negative
      finite `mag` (IEEE `x - x` under round-to-nearest), and a square is
      never `-0.0`. Worth noting as the general shape of these: the
      licence usually holds, but *which* input breaks it is not always the
      one you would guess, so enumerate the range rather than spot-check.
    - Measured with interleaved `quickbench` A/B, not mca: `clog` is
      deliberately not mca-wired (its branching risks the multi-exit-path
      region-marker corruption `mca_target.rs` documents, and its cost
      used to be "just its already-measured constituents" — an argument
      this change is precisely what invalidates). 3 reps, new ahead every
      time: throughput **5.66/5.35/5.25 vs 5.92/5.62/5.73 ns** (~5-8%),
      latency flat (36.54/35.94/35.95 vs 37.38/35.82/35.95).

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

38. **softplus fused kernel**: one fitted poly replacing log1p's
    division + ln machinery. **Shipped 2026-07-27** — a real win on
    *every* axis, which is rare enough here to be worth stating plainly.
    - The idea as written proposed fitting `ln(1+2^-t)` over `exp`'s own
      `k/f`-reduced domain, i.e. fusing *through* the exponential. That
      is not what shipped and is not necessary: the whole win is
      available one level up, at `log1p_unit`, whose callers
      (`softplus`/`logaddexp`, via `e = exp(-something.min(87.0))` with
      `something >= 0`) already guarantee `e` in `(0, 1]`. A bounded
      domain means there is no reduction to do, so a single fitted
      polynomial covers it directly and `ln_normal`'s exponent-field
      extraction, its degree-8 poly, the `k*LN2_HI + k*LN2_LO` recombine
      **and** the Sterbenz correction's `c / u` division all disappear at
      once. `exp` is untouched.
    - Form is `e + e^2*Q(e)`, not `e*P(e)`: the leading `e` is then exact
      and only the smaller correction carries the poly's rounding. It
      also reproduces the old code's tiny-`e` behaviour for free — `e*e`
      flushes to zero below ~1e-19 and the result is literally `e`,
      matching the old `1.0 + e == 1.0` path bit for bit.
    - Degree 9 in `Q`, Estrin with the top two coefficients folded in at
      the `e^4` level (`ln_normal`'s own trick, so `e^8` is never
      formed): **2 mul + 10 fma = 12 ops** against the old form's ~22
      plus a division. Idealized max relative error 0.078
      ulp-equivalent *after* rounding the coefficients to f32 (0.221
      straight off the scipy minimax fit; a coordinate descent over f32
      ulp steps recovered the 2.8x).
    - Measured, mca: `softplus` 104.14/4.100 -> **78.14/3.006**
      (throughput **-26.7%**, latency -25.0%), `logsigmoid`
      105.10/4.224 -> **79.11/2.969** (-29.7%), `logaddexp`
      104.14/4.100 -> **78.14/3.006** (-26.7%).
    - Measured, accuracy: `softplus` avg **0.0833 -> 0.0768** (exhaustive,
      2.24e9 in-domain samples), max **4 unchanged**; `logsigmoid` the
      same. `logaddexp` avg **0.155 -> 0.141**, a real ~9% improvement.
    - **Correction to this entry's first version, and a warning worth
      keeping.** It originally also claimed `logaddexp`'s documented
      cancellation max improved ~1e4 -> ~2e3, on 3 repeat runs per side
      (baseline 7782/14422/42548, new 1507/2463/2835). **That claim does
      not hold at 8 repeats per side**: baseline spans 992-38183 and the
      new code spans 1013-15556, fully overlapping, and a separate
      full-sweep run of the new code produced **249305**. The max here is
      a heavy-tailed cancellation artifact and is *not* a usable A/B
      signal at any repeat count this harness can afford — only the avg
      is, and even that is partly tail-driven (a single 2.5e5 sample
      shifts a 2.7M-sample mean by 0.09, most of the effect being
      measured). The avg claim survives because the two 8-rep ranges are
      *disjoint* (baseline min 0.1509 > new max 0.1451); that
      disjointness, not a difference of means, is what makes it safe to
      state. Stronger version of the existing 2-arg repeat-run rule: for
      a heavy-tailed metric, require non-overlapping ranges, and if they
      overlap report "no measurable change" rather than the direction of
      the means.
    - Why accuracy *improved* rather than merely holding: the old form
      spent its precision recovering bits that the new one never loses.
      `u = 1 + e` rounds away `e`'s low bits and `c = e - (u - 1)` claws
      them back through a division; a polynomial in `e` itself has
      nothing to recover. The remaining error is `exp`'s own, which
      propagates at a factor of `e/((1+e)*ln(1+e)) <= 1` — never
      amplified.
    - Only the fitting metric needed care. Fitting `Q` to minimize its
      own *absolute* error is the wrong objective: the quantity that
      matters is `e^2*dQ/ln(1+e)`, which weights the `e->1` end 1.44x and
      the `e->0` end at essentially zero. Re-weighting the LP by that
      factor moved degree 8 from 0.541 to 0.086 ulp-equivalent for free.

39. **logaddexp: same fused kernel** on |a−b| — **shipped with #38**, no
    separate work: `logaddexp` calls the same `log1p_unit`, so it picked
    up the identical -26.7% and the accuracy improvement above.

46. **parity(qh) bit-derivation**: |p0| < 2^22 → magic bits; p0 ≥ 2^24
    → deterministically even (every f32 there is an even integer); only
    the 2^22..2^24 window needs a select. Screen vs the floor-based
    parity. **Applies equally to `ql`, not just `qh`** — see the
    rejected `ql`-via-magic-round entry below: `rem` (what `ql` rounds)
    is *not* bounded the way its name suggests, so any bit-derivation
    scheme needs the same large-magnitude fallback this idea already
    plans for `qh`, on both words.
    **Rejected 2026-07-28 on a derivation, and the premise turns out to
    be false besides.** The premise — "after a magic-round, k sits in the
    low mantissa bits" — does not hold here: `round_x_over_pi` produces
    `qh = p0.round()` and `ql = rem.round_ties_even()`, i.e. two
    `vroundps`, **not** a magic add. There are no magic bits to read, so
    the idea would have to *add* a magic round first.
    - Even granting that, the op count loses. The shipped
      `parity(q) = fma(-2.0, (q*0.5).floor(), q)` is **3 ops and
      uniformly correct at every magnitude** — for `|q| >= 2^24`,
      `q*0.5` is exact, `floor` is the identity, and it returns 0, which
      is right because every f32 there is an even integer. A
      bit-derivation must handle three magnitude regimes, and the
      general form needs `(M >> (150 - E)) & 1`: extract the exponent
      (shift + and), form the shift count (sub), rebuild the implicit
      mantissa (and + or), variable-shift, mask — **7 integer ops** to
      replace 3 FP ops, on a CPU where IDEAS.md's own "parity() via
      integer bit-ops" entry already measured FP-port ops winning.
    - The magic-round shortcut is only valid for `|q| < 2^22`, and `qh`
      reaches 2^22 at `|x| ~ 1.3e7` — which is *inside* `sin_checked`'s
      accurate range, the whole reason `sin_checked` exists. So the cheap
      path would break exactly the function it was meant to speed up.
    - The final XOR combine is also already optimal: `pq == pl` +select
      (2 ops) versus the bitwise `(pq.to_bits() ^ pl.to_bits()) << 8`
      (also 2 ops, and it does produce the right `SIGN_MASK`) — a wash,
      so there is nothing to take there either.

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

54b. **asinh/acosh large-|x| branch: drop the sqrt for a two-term
    expansion** — found and **shipped 2026-07-27** (see lib.rs/git log).
    Distinct from idea #54's rejected single-sqrt restructure, which only
    moved *which* value the sqrt was applied to; this removes the sqrt
    from that branch entirely.
    - Both functions computed their large-|x| arm as
      `ax*(1 +- 1/ax^2).sqrt()`. But that arm only runs for `ax >= 2048`,
      where `sqrt(ax^2 +- 1) = ax +- 1/(2*ax) - 1/(8*ax^3) + ...` and the
      first dropped term is `1/(8*ax^4) <= 7e-15` relative — seven orders
      of magnitude under f32's own 6e-8. So `fma(+-0.5, 1.0/ax, ax)` is
      not an approximation at this branch's own domain, it is exact to
      f32.
    - Op count: a division, an add, a sqrt and a multiply become a
      division and an fma. `acosh` also sheds a whole `x*x` (its
      `1.0/(x*x)` becomes `1.0/x`, and its *other* `x*x` has to stay
      separate anyway — see its body comment on why the cancellation
      needs its own fma).
    - mca: `asinh` 5.406 -> **4.169 (-22.9%)**, `acosh` 4.718 ->
      **3.569 (-24.4%)**, latency flat (`asinh` 81.00 -> 81.02, `acosh`
      87.72 -> 89.02).
    - Accuracy: exhaustive over all 2^32 inputs, `asinh` identical
      (0.1493 avg / max 3, same worst x) and `acosh` marginally better
      (0.0597 -> 0.0596 avg, max 4, same worst x). `worst_corpus`
      bit-identical, all 8 standing gates pass.
    - Generalizable: this is the *same* lever as idea #38 one level down —
      a branch whose guard already proves a tight bound doesn't need a
      general-purpose kernel (there, `ln`; here, `sqrt`). Worth sweeping
      the crate for other guarded branches that still call something
      general: the guard is the licence.
    - **That sweep was run to exhaustion on 2026-07-27 and is now closed.**
      Five hits shipped, all of them accuracy-free (bit-identical, or
      bit-identical off NaN): **#54c** `atanh`, **#54d** `log1p`/`log2p1`/
      `log10p1` (the widest — its licence is the *argument's own
      construction*, `1+x`, not a branch guard, so it pays out at every
      caller), **#54e** `softplus`/`logaddexp`, **#54f** `erfinv`/`logit`,
      **#54g** `clog`. Headline throughput: `probit` -30.0%, `erfinv`
      -20.6%, `erfc_inv` -19.6%, `xlog1py` -20.1%, `log2p1` -19.0%,
      `log1p` -18.9%, `softplus`/`logaddexp` -15.7%, `atanh` -12.7%,
      `logsigmoid` -10.3%, `logit` -10.0%, `compound` -10.7%.
    - Screened and **rejected**, so the sweep doesn't get re-run on them:
      - `erfc`'s `exp2_checked` clamp is load-bearing (it deliberately
        feeds the *unclamped* `xa`, see its body comment).
      - `erfcx`'s `exp2_checked(x*x*LOG2_E)` has a dead *lower* clamp only
        (`x*x >= 0`), worth exactly one `vmaxps`. Not taken: extracting an
        unclamped core would put a new shared boundary in front of
        `powf`/`erfc`/`exp2m1` as well, and this crate has a +32%
        regression precedent for exactly that (see the macro-vs-fn dedup
        entry). One op is not worth exposing four callers to it.
      - `tanh_grad`/`sigmoid_grad`/`norm_pdf` all feed `exp_checked` an
        argument that is `<= 0` by construction, so their upper clamp is
        dead — but this saves *nothing*. `.clamp(lo, hi)` is already two
        instructions, and a one-sided replacement that still propagates
        NaN (`if v < lo { lo } else { v }`) is also two; the one-instruction
        `v.max(lo)` is wrong, since `f32::max` returns the *other* operand
        for NaN and would turn `f(NaN)` into a finite value.
      - `norm_pdf` additionally cannot use `exp_narrow`: its argument must
        stay clampable to `-104.665` for the result to reach exactly `0`,
        which is outside `exp_narrow`'s `[-87.68, 88.38]` domain. Clamping
        at `-87` instead would freeze `norm_pdf` at ~`6.4e-39` forever —
        the same freeze bug `erfc`'s own doc comment records.

54h. **`srgb_to_linear`/`linear_to_srgb`: the toe's guard licenses a
    stripped `log_2`** — a sixth hit, **shipped 2026-07-28**, and the
    reason it is worth recording separately is that #54b's sweep
    declared itself closed while these two were still standing. The
    sweep had been run over sites that call `log1p`/`ln`/`exp_checked`/
    `sqrt` directly; these reach `log_2` one level down, through
    `powf_pos`, so a grep of the kernel names missed them.
    - Both are piecewise: a linear toe below the seam, a power curve
      above. The power arm's value is *discarded* for every input below
      the seam — which is exactly the range where its base could be
      zero, negative or denormal — so `log_2`'s `denormal_rescale!`, its
      `x == 0.0 -> -inf` spec and its `x <= 0.0` select all compute
      results nothing can observe. Only `!(x < inf)`'s `x*x` survives,
      because `+inf`/NaN *do* flow through the outer select.
      `powf_pos`'s own two overrides go too: `y == 0.0` is a compile-time
      constant at both call sites, and `x == 1.0` is redundant because
      `log_2(1)` is exactly `0.0`, so the general formula already
      returns exactly `1.0` there.
    - The licence is subtler than #54c-#54g's and worth stating in the
      form that generalizes: every arm is still *evaluated* for every
      input (branchless), so what makes this sound is not that the bad
      inputs never arrive, it is that the caller's own select provably
      throws the answer away when they do. `+inf`/NaN are the exception
      precisely because they *don't* get thrown away.
    - **Bit-identical to the `powf_pos` composition over all 2^32 inputs
      for both functions** (checked directly against the old formulation,
      not inferred). All 9 standing gates pass, `worst_corpus` unmoved.
    - Screened by instruction count *first*, deliberately — mca's
      throughput column is not trustworthy in this neighbourhood (see the
      `coshm1` note below, and note `srgb_to_linear` was itself priced
      +20% over `powf_pos` for +3 instructions). `srgb_to_linear_
      throughput` **163 -> 132** instructions (-19.0%),
      `linear_to_srgb_throughput` **164 -> 134** (-18.3%), latency
      regions -1154 and -550, whole-file total 363164 -> 361399. **No
      other region in the file moved by a single instruction**, so
      nothing downstream pays for it.
    - mca then agreed, in direction and roughly in size, on both axes:
      `srgb_to_linear` **4.772 -> 3.666 (-23.2%)** throughput, latency
      108.02 -> 105.02; `linear_to_srgb` **4.285 -> 3.411 (-20.4%)**,
      latency 111.92 -> **83.69 (-25.2%)**. Exactly 2 of 145 rows moved,
      matching the asm diff. Worth noting as the healthy case: when the
      instruction count and mca agree, neither needed arbitration — it is
      only when they *disagree* that the asm wins.
    - Reusable: sweep by *reachability*, not by call-site grep. A
      composite that calls a composite that calls the general kernel is
      the same lever, one indirection further out — and `powf_pos`/
      `powf`/`exp_checked` are the wrappers most likely to hide one.
    - The reachability re-sweep turned up exactly one other live site,
      and it is deliberately small: `log1pmx`'s direct arm is
      `log1p(x) - x` under an `|x| >= 0.5` guard, so `log1p`'s trailing
      `x == 0.0` signed-zero select can never reach the result.
      `log1pmx_throughput` **155 -> 151** instructions (-2.6%), latency
      region -25, nothing else moves. Taken because it costs no
      duplication — `log1p`'s body became a `log1p_nonzero!` macro that
      `log1p` itself wraps with the one select — not because -2.6% would
      justify a copied body on its own. Everything else `log1p` does is
      genuinely live on that arm: `u = 1+x` really is `0`, negative and
      `+inf` there.
    - Screened and found *closed* on the re-sweep: `rootn`, `xlog1py`,
      `compound` and `powf`/`signed_pow` all feed the general kernel an
      unconstrained argument (`|x|`, `y`, `x` respectively), so none of
      the edge arms are dead; `pow_2_3`/`rcbrt`/`pow_3_2` bottom out in
      `cbrt`/`sqrt`, which supply their own special cases for free.

54e. **softplus/logaddexp: `exp_narrow`, not `exp`** — third hit from
    #54b's sweep, **shipped 2026-07-27** (see lib.rs/git log). Both
    already clamp their own exponent argument (`-ax.min(87.0)` /
    `-d.min(87.0)`, landing in `[-87, 0]` for *every* input including
    NaN, since `min` follows IEEE `minNum`), and that clamp is exactly
    the guard `exp_narrow`'s single-exponent-field domain
    (`[-87.68311, 88.37627]`) asks for — so `exp`'s `k1`/`k2`
    `exp2_field_split` was dead weight at both call sites.
    - Bit-identical: `softplus` verified over all 2^32 inputs,
      `logaddexp` over 200M random pairs plus a 16M dense grid across
      the `|a-b|` band where the correction term is live. `t1 * t2` and
      the single field are the same exact power of two over this `k`
      range, so no double-rounding difference exists to find.
    - mca, both axes: `softplus`/`logaddexp` 3.006 -> **2.534 (-15.7%)**
      throughput, 78.14 -> 73.14 latency; `logsigmoid` (a `-softplus(-x)`
      composite) 2.969 -> **2.663 (-10.3%)**, 79.11 -> 74.11. No other
      row in the table moved. All 8 gates pass.

54c. **atanh: inline `log1p`, minus the branches its own guard makes
    unreachable** — the first hit from the crate-wide sweep idea #54b
    recommends ("the guard is the licence"), **shipped 2026-07-27** (see
    lib.rs/git log).
    - `atanh`'s big arm is `0.5*log1p(2a/(1-a))` with `a = |x|`, so the
      `log1p` argument `v` is non-negative for the whole in-domain half
      and `u = 1+v >= 1`. That kills three of the general kernel's
      branches outright: `ln`'s `denormal_rescale!` (`u` is never
      denormal), the wrapper's `x == 0.0 -> -inf` select (`u` is never
      zero), and `log1p`'s own `x == 0.0` signed-zero select (`v == 0`
      only at `a == 0`, where the small-poly arm is selected anyway).
      The two out-of-domain arms are still live and kept verbatim:
      `u <= 0` for `a > 1`, and `!(u < inf)`'s `u*u`, which turns
      `a == 1`'s `u = +inf` into `+inf` and any NaN back into NaN.
    - mca: throughput 3.325 -> **2.903 (-12.7%)**, latency 75.47 ->
      **69.42 (-8.0%)**.
    - Accuracy: exhaustive over all 2^32 inputs, avg/max ulp 0.0037/2 —
      unchanged to four decimals. `worst_corpus` bit-identical, all 8
      standing gates and `codegen_check` pass.
    - Two more sites the same sweep turned up, both smaller (only the
      denormal rescale and the `x == 0.0` select are dead, the zero arm
      stays live): `logit`'s `log1p(-p)` (`1-p >= 2^-24` for any `p < 1`)
      and `erfinv`'s `-log1p(-x*x)` (same bound). Untested.

107. **PI_A..D bit-allocation joint refit** — **shipped 2026-07-28, and
     it is the largest free accuracy win in this file**: `sin`'s
     exhaustive `|x| <= 1e6` max ulp **58 -> 3**, `cos`'s **88 -> 3**,
     avg -12.6% / -6.2%, for **provably zero perf change** (the whole
     `mca_target.s` differs in 24 lines, all of them `.long` rodata
     constants — not one instruction moves).
     - Only `PI_C`/`PI_D` change, to full f32 width:
       `PI_C = 6.278329465203569e-7` (`0x3528885a`),
       `PI_D = 1.0780605906948477e-14` (`0x284234c5`). `PI_A`/`PI_B` were
       already optimal — their 8 and 9 significant bits keep steps 1 and 2
       exactly representable with 32x and 3.4x margin, and widening either
       one regresses (it breaks that exactness, which the closed form
       predicts and the fuzz confirms).
     - **Root cause, and it is embarrassing in a useful way: these were
       Sleef's `PI_Af..PI_Df`, designed for Sleef's own
       `TRIGRANGEMAXf = 39000`, i.e. `|q| <= ~12400`.** This crate reuses
       them at `|q| <= 2^22` — 340x larger. `PI_C`'s 9-bit budget buys a
       step-3 exactness window of `|r| <= 0.031` that nothing needs, and
       pays for it with a residual `delta_4 = pi - sum(PI_i)` of
       -2.435e-18. Widening `PI_C`/`PI_D` drops `delta_4` to -1.906e-22
       (12800x), putting `q*delta_4` far under the half-ulp-of-`r` floor.
     - The mechanism to remember: a rounding *inside* the chain costs half
       an ulp **of the residual**, which is harmless even at sin's zeros
       because there the residual *is* the answer. What `delta_4` costs is
       an **absolute** error scaled by `q`, which near a zero is unbounded
       *relative* error. Those are not the same currency, and the shipped
       split was spending bits on the wrong one.
     - Does **not** move the in-domain max (`sin` 219, `cos` 2769) and
       cannot: that comes from `FRAC_1_PI`'s own 1.28e-8 error making `q`
       off by one past ~1.3e7, pushing `r` outside `sinf_poly`'s fit. All
       15 screened candidates hit the same in-domain max. Off-domain
       (`[1.32e7, 1e9]`) both splits are equally garbage.
     - **Methodology warning, and it is a sharp one: the 100M quick fuzz
       reported the baseline `|x|<=1e6` max as 4, where the exhaustive
       truth is 58.** A 14x under-report. The failing points are a thin
       set that random sampling essentially never lands on, so for this
       bucket quick-fuzz max is not merely noisy — it is systematically
       optimistic. Use `thorough` for any max claim about the fast trig
       tier.
     - Search: 15 splits (`n1 in {8,12,16,24}`, `n2 in {9..16,24}`,
       `n3 in {9,12,13,16,20,24}`) scored on the real sin/cos fuzz via a
       runtime-constant screener; monotone improving in `n3`. `8_11_24_24`
       measures 0.3% better still but cuts step-2's exactness margin from
       3.4x to 1.46x for that, and changes a third constant — declined.
     - `worst_corpus` moves 18 entries, all `sin`/`cos`/`tan`; every
       in-domain one moves **to 0 ulp** against an 80-digit reference.
       Blessed. All 10 gates pass.

110. **±few-ulp exhaustive scan of every non-poly literal** (clamp
     bounds, seed constants, magic offsets, branch thresholds) scored on
     the real fuzz — #3's sibling for non-coefficient constants.
     **First application run 2026-07-27 on `tanh`'s 0.25 seam, the
     best-motivated target available (`error_profile` reports `tanh`'s
     whole max-6 worst case sits in the band that straddles it, at
     x ~ 0.2553). No headroom: rejected.** The scan itself is cheap and
     worth reusing — score *each arm alone* over every f32 in a band, then
     evaluate every candidate seam at once via a prefix-max of the
     below-arm and a suffix-max of the above-arm. Result: the optimal
     combined max is reached for any seam in **[0.2499507, 0.3039100]**,
     i.e. the shipped 0.25 already sits (barely) inside the optimal run,
     and every seam in it scores max 6. The avg-minimising point, 0.2716,
     buys **0.1457 -> 0.1456 exhaustive avg** — nothing, and it moves 3
     `worst_corpus` entries with two of them locally *worse* (x = +-0.25
     exactly, 0.21/1.21 -> 1.79 ulp). Diagnosis matches the crate's usual
     one: both arms independently reach ~6 across the whole plateau, so
     the band's error is rounding-chain, not seam placement — fixing it
     needs a better arm (wider Pade refit or a third branch), which is not
     free.
     - **Methodology trap, and this one cost a false positive:** score the
       seam scan with **the harness's own `ulp_diff`** (integer bit
       distance to the f64 reference *rounded to f32*), not a
       true-value-relative `|got - exact| / ulp(exact)`. The first run
       used the latter and reported a clean max **6.084 -> 5.557** win for
       moving the seam to 0.27 — which is a real statement about the true
       error and still rounds to bit distance 6 on both sides, so the
       shipped metric showed no change at all. A sub-ulp difference in
       true error is invisible once both candidates land in the same
       rounding bin. Same family as the "score the shipped fn, not the
       intermediate" trap: score what the gate scores.
     - **Second application, `asin`'s 0.27 seam — also no headroom, but a
       more useful negative.** `asin` carries the crate's largest seam gap
       (6 ulp, `edgecheck`'s own seam report) and its max ulp is 6, so it
       looked like the same story as `tanh`. The optimal run is
       **[0.2608490, 0.2994175]**, 0.27 is inside it, and 0.27's band-avg
       (0.25473) is within 0.0002 of the run's avg-minimising point
       (0.2749 -> 0.25452). Nothing to move.
     - What the scan found *instead* is where `asin`'s max actually lives,
       which the seam gap had disguised: the big arm's max stays 6 for
       every candidate seam right up to 0.45 and only drops to 3 at 0.50,
       and **2875 distinct f32 in [0.27, 1] reach 6 ulp**. So this is a
       broad plateau across roughly `[0.27, 0.5)`, not a seam artifact.
       Cause is the big arm's own cancellation: `FRAC_PI_2 - sqrt(1-a)*P(a)`
       at `a ~ 0.3` subtracts 1.266 from 1.5708 to get 0.3047, a ~4.2x
       amplification, so ~1 ulp on the subtrahend is ~4 ulp out. No refit
       of `asin_poly` addresses that.
     - The only fix that would is a **third, mid-range branch** evaluating
       `asin` directly (an `a*P(a^2)` minimax on `[0.27, 0.5]` has no
       cancellation at all). Priced before building: it is a whole extra
       poly evaluated unconditionally plus a select, on a function that is
       already one of the crate's cheapest (0.968 cyc/elem) — an accuracy
       win with a real perf penalty, i.e. the wrong side of this session's
       bar. Extending `asin_small` instead is worse: it is degree 3 in
       `x^2` and already 7 ulp at 0.30, and asin's series decays by only
       `x^2 = 0.25` per term, so reaching f32 accuracy at 0.5 needs ~12
       terms against today's 4.
     - Reusable conclusion for the scan: **a large seam *gap* does not mean
       a misplaced seam.** Both functions scanned had their max at or near
       the seam and both seams were already optimal; in `asin`'s case the
       gap was pointing at a plateau that merely happens to start there.

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
       - **Screened 2026-07-28, and the cost estimate above is wrong by
         about 5x — in the idea's favour.** "A full `ln` on every call"
         would be `ln`'s 1.611 cyc/elem on top of `exp_m1_over_x`'s
         1.798, i.e. +56-89%, far outside this crate's precedent for
         totality (`exp -> exp_checked` pays +30%). But a full `ln` is
         not what the arm needs: it is only ever *selected* over
         `x` in `[88.7228, 93.2582]`, and across that whole band
         `ln(x)` ranges over `[4.4856, 4.5354]` — a span of **0.05**.
         An absolute error `d` in the exponent is a relative error `d` in
         the result, so f32 accuracy needs `ln(x)` to ~6e-8 absolute over
         a 4.5-wide interval, which a degree-4/5 minimax in `(x - 91)`
         supplies in **~4 fma**. Estimated ~+0.3 cyc/elem, ~+17%, and it
         is a *new* function so nothing existing regresses.
       - That is the "guard is the licence" lever (#54b) again, one more
         level out, and worth naming as a general move: **when a branch
         is only selected on a narrow band, it does not need the general
         kernel — it needs a kernel fitted to the band.** #38's fused
         `log1p_unit` and #54b's dropped `sqrt` are the same shape. The
         reason it was missed here is that the entry reached for `ln` by
         name and then priced `ln`.
       - Not built: `exp_m1_over_x_checked` does not exist yet, so this
         is new API surface rather than an improvement to a shipped
         function, and that is a scope call rather than a measurement.
         The screen above is the part that was missing.

136. **normalize2/normalize3 *slice* kernels** (rhypot + scales — the
     operation users actually want hypot for). Single-call
     `normalize2`/`normalize3`/`normalize4`/`hypot4`/`rnorm4` shipped
     (see lib.rs/git log, ideas #55/#134) -- the slice-batched variant
     itself still needs the slice-tier infrastructure that doesn't
     exist yet (see the standing Batch/slice API tier entry above).

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
     documents its behavior; most others are unaudited). **Shipped
     2026-07-27** as `examples/denormal_audit.rs`. Answer: **the tier is
     in good shape — only `sigmoid` flushes materially, and its gap is
     already documented.**
     - Audited in two separate classes, because they fail for different
       reasons. **(A) normal input -> denormal output**, where the
       function's own machinery has to carry a result below `MIN_POSITIVE`
       without flushing (a bit-trick exponent field structurally cannot,
       so this is where early saturation lives). **(B) denormal input ->
       denormal output**, the near-zero identity region, where the risk is
       the reverse: a reduction or rescale mangling a subnormal argument.
     - **(B) is completely clean**: `sin`, `tan`, `asin`, `atan`, `sinh`,
       `asinh`, `atanh`, `expm1`, `expm1_checked`, `log1p`, `softsign`,
       `wrap_pi` all carry denormals at *exactly* zero relative error
       (they are the identity there), and `erf`/`sinpi`/`dawson`/
       `sqrt1pm1`/`gelu`/`silu` are correct to within one representable
       step at the very bottom.
     - **(A): `exp2_checked`, `exp_checked`, `exp10_checked`, `erfc` and
       `erfc_accurate` all carry denormals with zero premature flushing.**
       Only two flush at all: `sigmoid` (94% of its denormal-output range)
       and `norm_pdf` (1%).
     - The decisive metric is **how early a flush starts relative to where
       the true result genuinely rounds to f32 zero**, not what fraction of
       samples flushed. By that measure `sigmoid` is premature by **15.6 in
       x** (flushes from -88.38 where the true zero is at -103.97 — i.e.
       the whole `(-104, -88.4)` band, exactly the gap its own doc states),
       while `norm_pdf` is premature by only **0.0125 in x**, one
       representable unit, not worth an op to fix. Everything else: never
       premature.
     - Two harness traps worth remembering, both of which I hit and fixed:
       - Classifying "denormal output" by the f64 magnitude alone
         (`|ref| < MIN_POSITIVE`) is wrong — it admits values under ~7e-46
         that round to f32 *zero*, where returning 0 is correct. That
         inflated `erfc` to a reported "39% flushed" when it is not
         premature at all. Require the correctly-rounded f32 answer to be a
         nonzero denormal.
       - Sweeping an `_unchecked` tier past its documented domain measures
         documented garbage. An early version swept unchecked `exp2` to
         -150 and reported a meaningless 5.8e76 relative error; `exp2`'s doc
         scopes it to non-denormal results, so it has no denormal outputs
         within its own domain at all.
     - Also note a reporting artifact that is *not* an error: several
       clean functions show "worst rel ~1.0" alongside **zero** flushes.
       At the smallest denormals only one or two representable values
       exist, so returning the correctly-rounded one is still ~100% off in
       relative terms. Read the flush columns, not the relative error, at
       the very bottom of the range.

168. **Worst-case corpus regression gate**: persist each function's
     known worst-x list, re-check every commit in seconds between the
     hours-long full sweeps. **Shipped 2026-07-27** as
     `examples/worst_corpus.rs` + `examples/support/worst_corpus.golden`:
     9720 entries (108 public 1-arg functions x 90 inputs) checked in
     **0.058s**.
     - Design choice that got this out of the backlog: it is
       **reference-free**. A ulp gate needs an f64 reference per function
       (100+ of them, which is why this sat unbuilt); comparing against
       *blessed output bits* instead answers "did anything move?" in
       milliseconds with no reference machinery. The sweeps still answer
       "is it correct?". `--bless` regenerates.
     - Corpus is the special values, the branch seams this crate actually
       has, every clamp boundary, the recorded worst-x values from IDEAS.md
       and the readme, and a magnitude spread across the exponent range.
     - **Validated by deliberate perturbation, and the negative results
       matter more than the positive one.** A 0.5% change to the shared
       Pade coefficient moves 61 entries and exits 1, as intended. But two
       smaller real changes slipped through: (a) a *1-ulp* change to that
       same coefficient is invisible, because it multiplies `v*v` against a
       `-120.0` term so it shifts the sum by ~1e-10 relative, far under
       f32's ~6e-8 resolution; (b) moving exp10m1's seam `0.2 -> 0.21` is
       invisible, because the two branches agree bit-for-bit at the corpus
       point `x=0.2` — a well-placed seam is *supposed* to, which is
       exactly why seam points make weak canaries.
     - So: a pass is not "nothing changed", and a corpus gate is only as
       good as its input list. Documented in the file header so it can't
       create false confidence. Worth noting the general lesson too — when
       validating any new gate, perturb at *several* magnitudes, because a
       single decisive perturbation proves only that the plumbing works.

169. **ULP-error histogram artifacts** per function (not just avg/max)
     — bimodal structure reveals branch-split opportunities. **Shipped
     2026-07-27** as `examples/error_profile.rs` (ulp histogram + a
     per-magnitude-band avg/max/worst-x breakdown, annotated with which
     side of each function's seam a band falls on). Ran it on the three
     joint-worst budgeted functions (`expm1`, `exp_m1_over_x`, `tanh`, all
     max ulp 6) and it produced a real, actionable, mechanistically
     explained diagnosis — **all three have their error concentrated in a
     narrow band at and just above their seam**, not spread over the
     domain:

     | function | straddling seam | just above | rest of domain |
     |---|---|---|---|
     | `expm1` | avg 1.116, max 5 | avg 0.877, **max 6** | ~0.35 |
     | `exp_m1_over_x` | avg 1.033, max 5 | avg 1.020, **max 6** | ~0.39 |
     | `tanh` | avg 0.977, **max 6** | avg 0.831, max 4 | ~0.40-0.50 |

     - Mechanism, verified numerically rather than guessed: the direct
       branch forms `e^x` and then subtracts 1, so it loses
       `log2(e^x/(e^x-1))` bits to cancellation, and that factor is
       *largest immediately above the seam*: **2.541** at `x=0.5` (1.35
       bits), 1.939 at 0.725, 1.582 at 1.0, 1.457 at 1.16, then 1.157 by
       `x=2` and 1.007 by `x=5`. The measured error elevation tracks that
       curve exactly — elevated while the factor is above ~1.5, back to
       baseline once it drops under ~1.3.
     - Neat confirmation that it really is this mechanism: `tanh`'s seam
       sits at `0.25` and it feeds `expm1(2x)`, so its amplification at
       the seam is the factor at `2*0.25 = 0.5` — **the same 2.541**, and
       its error band is correspondingly at `|x| in [0.21, 0.49]` rather
       than `[0.45, 1.16]`. Two different functions, one shared cause.
     - **Actionable levers this points at**, in preference order. (a)
       Widen the shared Pade's *fitted* domain (currently `|v| < 0.5`) so
       the seam can move up past the high-amplification region — blocked
       today because moving the seam alone pushes the Pade outside its fit
       (its own doc: max 3/avg 0.109 over `|v|<0.5`), and the Pade is
       shared by 5 callers so this needs the full mca sweep. (b) A
       dedicated third branch over just the elevated band — precisely the
       shape that rescued `atanh` (idea #69), where a loss confined to a
       narrow region became a win on every axis once it got its own
       branch. (c) Note what this rules *out*: a pure seam retune cannot
       fix it, which is consistent with the 5-function crossover audit
       already finding `0.5`/`0.25` optimal — the seam is in the right
       place, the *branch* is what degrades near it.
     - Reusable: "one max ulp figure" hid a 3x avg elevation confined to
       less than one octave in all three functions. Worth profiling before
       any future refit, since it distinguishes refit / seam-move /
       new-sub-branch, which the avg+max pair cannot.
     - **Lever (b) built and measured on `expm1` the same day; clean
       accuracy win, rejected on cost.** The identity used avoids a refit
       entirely: `expm1(2u) = a*(a+2)` with `a = expm1(u)`, so halving the
       argument puts it back inside the Pade's own fitted `|v| < 0.5`
       domain. Its own amplification `(2a+2)/(a+2)` is 1.124 at `x=0.5`
       against the direct arm's 2.541 — **2.26x better** — staying ahead
       until they cross at about `x=1.4`.
       - Two-branch form (doubling for all `|x| < 1.0`, seam moved 0.5 ->
         1.0): max **6 -> 5** but avg **0.1304 -> 0.1506, worse**. The
         regression is the extra `fma` rounding, now paid at every small
         `|x|` where there was no cancellation to fix in the first place.
         Not a clean win.
       - Three-branch form (plain Pade below 0.5, doubling over
         `[0.5, 1.0)`, direct above): **clean win on both axes** — avg
         0.1304 -> **0.1268**, max 6 -> **5**, worst x moving to 1.0389,
         just above the new seam exactly as the mechanism predicts. But
         throughput **1.695 -> 2.366 (+39.6%)**, because the doubling arm
         is a whole second Pade evaluation *including its own division*.
         Rejected: far outside anything an avg/max improvement of this size
         justifies. (Ignore the latency figure it reports, 71.00 -> 39.89:
         adding work cannot halve latency, and this is the documented mca
         latency-harness artifact. Throughput is the trustworthy axis.)
     - **Lever (a) screened and disqualified the same day, before
       implementation.** The plan was to refit the shared Pade over
       `|v| < 1.0` and move the seam to 1.0 with no doubling — attractive
       because it keeps today's op count exactly (one Pade, one select), so
       unlike lever (b) it would have been structurally free. A scipy
       minimax refit of the three free coefficients shows the fit simply
       cannot be had at this degree. Best achievable idealized max relative
       error, in f32-ulp-equivalents:

       | domain | shipped coeffs | optimal same-degree refit |
       |---|---|---|
       | `\|v\|<0.5` (today) | 1.37 | **0.65** |
       | `\|v\|<0.65` | 13.9 | 3.15 |
       | `\|v\|<0.8` | 59.6 | 10.98 |
       | `\|v\|<1.0` | 264.6 | **42.04** |

       At the seam-1.0 target the *idealized* error alone is 42 ulp, before
       a single rounding is added, against today's whole-function max of 6.
       Even a modest widening to 0.65 costs 3.15, which would roughly
       triple the Pade branch's own max (currently 3). The scaling is
       brutal — ~65x across `0.5 -> 1.0` — as expected for a [2/3] rational
       whose error term goes like `v^7`. Making the seam move work needs a
       *higher degree* Pade, which per this file's own repeated finding
       ("degree bump costs real": 3/3 such bumps gave real accuracy wins
       and real mca cost) would be paid by all 5 callers of the shared
       macro.
     - **Net: the #169 concentration is structural at the current op
       budget.** All three levers are now closed with evidence — (a)
       disqualified by the fit screen above, (b) measured at +39.6%
       throughput, (c) ruled out because the seam position is already
       optimal and it is the *branch* that degrades near it. `expm1`'s
       max 6 is not a tuning oversight.
     - Incidental finding from the same screen, noted but *not* pursued:
       the shipped coefficients are **not** minimax-optimal on their own
       `|v|<0.5` domain either — 1.37 ulp-equivalent against an achievable
       0.65, so ~2.1x idealized headroom. Deliberately left alone: that
       margin sits squarely in the weak-to-moderate band this crate has
       seen fail repeatedly (see the LP-margin entries), the coefficients
       were tuned against the *real chain* rather than an idealized model
       so a lower idealized error need not transfer, and #169's own
       profiling shows the Pade branch is not where the max lives anyway.

175. **NaN-payload/quietness propagation matrix** (which ops
     canonicalize payloads) — documentation-grade completeness.
     **Shipped 2026-07-27** as `examples/nan_payload.rs`. The quietness
     half is already *asserted* by `special_matrix.rs` (NaN in -> quiet
     NaN out, which IEEE754 does require); this covers the payload half,
     which it does not — a function may return an input NaN unchanged or
     mint a fresh canonical one, and both are legal. So this is a record,
     not a gate. Result over all 108 public 1-arg functions, feeding three
     tagged payloads (including a sign-negative NaN):
     - **74 preserve payload *and* sign exactly.** The common case, and
       the mechanism is simply that a payload survives any operation that
       merely propagates an operand — multiply, fma, select.
     - **22 keep the payload but vary the sign**, which is what any
       `copysign`/`mulsign`/negation in the tail does. Includes most of
       the trig family (`cos`, `cosd`, `tand`, ...) and the `cbrt` tiers.
     - **4 canonicalize to the default quiet NaN**: `acosh`, `asinh`,
       `softplus`, `logsigmoid`. Worth knowing *why*, since it predicts
       where else this will happen: all four route through a hardware
       `sqrt` or a comparison-driven select, and x86's `vsqrtps`
       canonicalizes rather than propagating its operand's payload. So
       "does my payload survive" reduces to "does the tail contain a
       sqrt/select that synthesizes a NaN from scratch".
       **Now 3, and the membership changed twice, 2026-07-28**: `erfinv`
       joined (idea #54f collapsed its two arms into one select), and
       `acosh`/`asinh` *left* — merging their `is_infinite`/`is_nan`
       selects into one `is_finite` select means the NaN they return is
       the input rather than a fresh `f32::NAN`. Current counts: **76
       preserve payload and sign, 3 canonicalize (`erfinv`, `softplus`,
       `logsigmoid`), 21 vary the sign, 8 return no NaN.** All legal, all
       still quiet, `special_matrix` unaffected — this file is a record,
       and the record simply moves when the tails do.
     - **8 return no NaN at all** — the `_unchecked`/`_approx` tiers whose
       docs promise nothing off-domain (`exp2_approx`, `ln_unchecked`,
       `log10_unchecked`, `log2_approx`, `log_2_unchecked`, `rcp_approx`,
       `rsqrt_approx`, `sqrt_approx`). Broken out as their own class rather
       than counted as "mixed": a first version lumped them in, which
       overstated the mixed count at 30 and read as though a third of the
       library had inconsistent payload behaviour.

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
       operand.

201. **Two more harness facts, both learned the expensive way
     2026-07-28** while running several experiments in parallel:
     - **Never run `examples/mca` concurrently with any other
       `--emit=asm` build of a *different* source state.** `mca.rs`
       touches `mca_target.rs`'s mtime, runs `cargo rustc --emit=asm`,
       then hands the resulting `.s` to `llvm-mca` — and a concurrent
       build writes the *same* path. llvm-mca then reads assembly that
       does not correspond to the tree mca thought it was measuring, with
       no error and no obvious tell. If you are testing several ideas at
       once, give each its own worktree (separate `target/`), or
       serialize the mca runs.
     - **Region-level llvm-mca is ~200x faster than the full harness and
       reproduces its numbers exactly.** Extract just the
       `LLVM-MCA-BEGIN`/`END` regions you care about into a small `.s`
       and run `llvm-mca` on that: seconds instead of ~45 min, verified
       against a full run. It also unlocks `--bottleneck-analysis`, which
       is what separates a *resource*-bound regression from the
       register-dependency artifacts (see the evaluation-order entry:
       `erfinv` was 86% resource-pressure-bound and its regression was
       real, while three same-shaped "regressions" at 50-90%
       register-dependency-bound were not).

202. **powf absorbs powf_checked; log_2 peels its leading term
     (2026-07-29).** `powf_checked`/`powf_checked_unchecked` are gone,
     not deprecated: `powf` *is* the double-float route now
     (`exp2_checked_df(log2_df(ax) * y)`), >=292 -> **3 max ulp**
     (`examples/powfsearch.rs`), avg 0.181 -> 0.019. Against the checked
     tier it replaces, mca latency **-18.1%** (152.89 -> 125.17) and
     throughput **-14.7%** (9.918 -> 8.459); against the inaccurate
     formula it replaces, +19.3% / +49.7%. The tiering was never a real
     speed/accuracy curve -- the fast route's error is `y`-amplified, so
     it is wrong by hundreds of ulp at ordinary inputs, not "approximate".
     Five mechanisms, each measured separately, and the last three
     generalise:

     - **Route degenerate `ax` through the formula, not a fallback.**
       `powf_checked` decided `ax == 0`/`inf`/NaN with an
       `is_safe`/`zero_or_inf`/`edge_mag` tree (~13 ops) because it needs
       the sign of `y` as well. Two selects on the log2 *high word only*
       (`-inf` for `+0`, `x*x` for `+inf`/NaN, ~5 ops) do the same job,
       because `exp2_checked_df`'s clamp already saturates each case
       correctly and the low word can be arbitrary -- every path that
       reads it ends at the guard that already exists for overflow.
     - **A `{+1,-1,NaN}` sign multiplier beats a select tree.**
       `powf_sign_combine!` was ~28 vector ops, most of them from
       `y_int || x == 0.0 || x.is_infinite()`: a compound `||` compiles to
       a hand-assembled AVX-512 mask (`korb`/`kmovd`/`cmovnel`), the same
       trap `powf_checked`'s own `is_safe` comment documents. Rebuilt at
       ~20, every compare feeding exactly one blend. Two structural
       simplifications fell out: `parity(y)` alone separates even/odd/
       non-integer, so the separate `y == y.trunc()` test is redundant;
       and pinning `mag` to `1.0` when `ax == 1.0` subsumes *both* the
       `pow(1,y)` and `pow(-1,+-inf)` overrides. `ax + ax == ax` picks out
       `+0` and `+inf` (the C99-exempt magnitudes) in one compare, and
       `y - y != 0.0` tests infinite y with no constants at all.
     - **The exact denominator was free all along.** `log2_df`'s quotient
       refinement needs `s - th*d` against the *exact* `d = m + 1`, which
       `m + 1` is not, and it was carrying a `(dh, dl)` split for it. But
       `d == s + 2` exactly, so `s - th*(s+2) == (s - 2*th) - th*s`, and
       `s - 2*th` is itself Sterbenz-exact (`s/(2*th)` is `d/2`, in
       `[0.854, 1.207]`). Two fma, against four ops for the split; `d`
       survives only as the Newton step's operand, where approximate was
       always enough. Bit-identical output, -2 ops.
     - **Guard the correction, not the corrected value.**
       `exp2_checked_df` ended in `fma(result, c, result)` +
       `is_finite`-select. Moving the correction onto `p` (the pre-`t2`
       partial product, always finite and *normal*) and the guard onto `c`
       (ready long before the poly) took **-10.3%** off the region's
       simulated throughput on its own, purely by getting two ops off the
       end of the dependency chain. It also fixed a denormal-output
       rounding (the correction used to be applied *after* the mantissa
       bits were gone) and, once the guard became `xs == v.0` instead of
       `c.is_finite()`, a real overflow bug: `powf(f32::MAX, 1.0000001)`
       clamps `128.0000153` to `128`, and a negative correction then
       pulled the answer back to `f32::MAX` where it must be `inf`.
     - **`exp2_field_split` on `exp2int_field!`'s magic**, replacing
       `(k + 383) << 8 & EXPONENT_MASK` per word: -2 ops, -3 constants,
       shared by ~15 functions. `sinh_checked` **-13.2%** throughput,
       `exp10_checked` -12.7%, `exp_checked` -8.8%, `exp_m1_over_x`
       -4.7%, down to `expm1` -1.0%, no row up. The two forms break ties
       oppositely (the magics differ in parity), so `k1` differs by 1 for
       odd `k` -- unobservable, because `t1` is an exact power of two and
       `fma(q, t1*f, t1)` is exactly `t1 * fma(q, f, 1)` with the same
       mantissa either way, and the single rounding into the denormal
       range still happens in the final `* t2`. Verified over every
       reachable `k` against 4096 mantissas, bit for bit.

     **`log_2`: max 3 -> 1 ulp** (matching std's own max), avg 0.0031 ->
     0.0028, exhaustively, with one instruction and one uOp *fewer*.
     Latency +11.2%, throughput -0.1% (`log2`) / +6.6% (`log2_unchecked`).
     The 3 was never the fit -- the degree-9 `P(s)` fit at 4.1e-9 was 16x
     tighter than it needed to be. Written `k + s*P(s)`, the answer for
     `x` near 1 (`k == 0`) simply *is* `s*P(s)`, so all three of `P`'s
     full-weight evaluation roundings landed on the result at ~2^-24 each.
     Peeling the leading `log2(e)*s` out of the polynomial demotes every
     one: what is left reaches the answer scaled by `s^2/log2(1+s) <=
     0.21`. Three things had to be right at once, and each was measured:
     - the peeled `Q` is fitted by an **ulp-weighted LP against
       `s^2/log2(1+s)`**, not by dropping a term off `P`; degree 8, since
       degree 7's 2.6e-8 fit (~0.44 ulp) lands back on max 2;
     - **`k` joins last, in its own rounding.** Threading it through the
       peeled combine (`fma(s, LOG2_E, fma(s2, q, k))`) rounds twice at
       `k`'s scale and costs **40x on the average** (0.0031 -> 0.1249)
       while still reaching max 2. Likewise `fma(s2, q, s * LOG2_E)`,
       which looks cheaper because `s*LOG2_E` starts early: it rounds the
       leading term on its own at full weight, which is the one thing the
       peel exists to prevent;
     - the `s^2` factor rides into the poly's own low group
       (`a = s2 * l0`) instead of multiplying the finished `Q`, which is
       what keeps the whole thing three Estrin levels deep. Multiplying
       afterwards costs a second level and a further +11% latency.

     Two things that did **not** work, both measured:
     - **Multiplicative correction in `log2_df`** (`hi * (1 + w)` via a
       `from_quick_fma`, replacing the `+ corr` two-sum): 2 ops cheaper,
       and **3 -> 10 max ulp**. Two separate reasons, and the second is
       the interesting one: the low word needs the same `(1 + w)` factor
       (~2^-30.7 on its own), *and* `w = u*G(u)` carries `u`'s full
       relative error where the absolute form only exposed `G`'s ~0.0176
       log-derivative to it. The absolute form's exactly-reconstructed
       `t^3 = th^2*(th + 3*tl)` prefactor is load-bearing, not
       bookkeeping. Fixing both costs more than the restructure saves.
     - **Dropping `log2_df`'s Newton step** for a degree-5 minimax seed
       (2^-20, one level shallower, one op fewer): powf latency -3.6%,
       throughput -2.2%, and **3 -> 4 max ulp** -- via the same `u`
       sensitivity above, not via `t` (whose `e^2` error is 2^-40 there,
       far past enough). Also Estrin-splitting `LOG2_ATANH_G` and
       regrouping to `(uh*gp) * t3`: a real level off the chain and -2.5%
       latency, at +2.2% throughput and +3 uOps.

     Harness note: `powf_checked`/`powf_checked_unchecked` are removed
     from `accuracy.rs`, `edgecheck.rs`, `mca.rs`/`mca_target.rs`,
     `quickbench.rs`, `unchecked_parity.rs`, `special_matrix2.rs` and
     `powfsearch.rs` (which now sweeps `powf` against `powf_unchecked`).
     `log_family_normal!` is gone with them -- it had exactly one caller
     left, `log_2_normal`, since `ln_normal`/`log10_normal` grew their own
     copies.

- **`erfc`/`erfcx` rebuilt on a reciprocal variable** — the entry that
  closes a dozen failed attempts above. Both functions shared a degree-4/4
  rational in `xa` (`erfc_rational`) that had to clamp `|xa| <= 10` to keep
  its own denominator from overflowing. Everything tried against it — a
  domain split, compensated Horner, a degree 5/5 bump, an ulp-weighted
  refit, centered variables, threshold tightening — moved avg and left max
  at ~109. **The reason none of them worked is that the rational's own fit
  error was ~15 ulp**, which no better *evaluation* of that rational can
  reach, and which nothing in the record had measured. That one number
  reframes the whole family of rejections above.
  - **What shipped**: `erfcx_pos(xa) = v*P(v)` with `v = 1/(2+xa)` and `P`
    a degree-10 minimax polynomial, shared by both functions. Two
    properties do the work. The reciprocal variable maps the whole
    half-line into `v` in `(0, 1/2]`, so there is no clamp to place and
    nothing to overflow. And the explicit leading `v` factor makes
    `erfcx(xa) ~ 1/(xa*sqrt(pi))` fall out of `P(0) = 1/sqrt(pi)` **by
    construction** rather than having to be fitted.
  - **The second half is the Gaussian factor**, worth as much as the
    polynomial. `exp(-x^2)` turns an *absolute* error in its exponent into
    a *relative* error in the result, and `fl(x*x)` alone, on a value up
    to 100, is ~6e-6 of absolute exponent error — dozens of ulp before any
    exponential runs. Keep the square exactly (`p = x*x`,
    `pe = fma(x,x,-p)`) and apply `exp(-(p+pe)) ~ exp(-p)*(1-pe)`, one fma
    folded into the `erfcx` factor. Then route through `exp_checked`, not
    `exp2_checked`: the naive `exp2(-x*x*LOG2_E)` rounds a *second* time
    forming that product, and f32's own `LOG2_E` is not precise enough to
    carry an exponent of magnitude ~144. Cody-Waite handles both.
  - **Measured, exhaustive, both sides against the same harness**:
    `erfc` 0.3055/109 -> **0.1993/7**, `erfcx` 0.3768/126 ->
    **0.2142/6**. `erf` untouched and confirmed unmoved (0.3176/4).
  - **`erfcx` is now correct on the whole real line**, which it never was:
    the old clamp froze it at `erfc_rational(10)` forever past `x = 10`,
    relative error growing without bound (~10% at 11, ~99% at 20, ~895% at
    100). A new exhaustive `accuracy.rs` row measures `x >= 20` out to
    `f32::MAX` against the asymptotic series (`exp(x^2)` overflows f64
    past x~26.6, so the composed reference cannot reach there): **4 /
    0.6478** over 1.04e9 samples.
  - **Cost is real and was accepted rather than argued away**: `erfc`
    throughput 2.437 -> 3.284 cyc/elem (**+34.8%**), latency 64.00 ->
    70.36 (+9.9%); `erfcx` 2.278 -> 2.792 (**+22.6%**), 62.02 -> 69.97
    (+12.8%). That breaks this crate's usual no-perf-penalty bar for a
    pure accuracy fix, on an explicit instruction to make these two good.
    The one favourable comparison: the retired `erfcx_checked` (the only
    thing that used to be correct past `x = 10`) cost 104.98/3.126, so the
    new `erfcx` does that job at **-33% latency and -10.7% throughput**.
  - **Three functions deleted with it**: `erfc_accurate` (0.2452/13),
    `erfcx_accurate` and `erfcx_checked` (both 0.3263/21) — the base
    functions now beat all three on both axes, so a separate tier had
    nothing left to offer. `tune.rs`'s four rational targets (`erfc_c`,
    `erfc_c5`, `erfc_lo_c`, `erfc_hi_c`) collapsed into one `erfcx_pos_c`
    matching the shipped Estrin order.
  - **A constraint that looked free and was not, worth the entry on its
    own.** `c0` *is* the asymptote — at large `xa` the result is `c0*v` and
    nothing else — so hardcoding it to the correctly-rounded `1/sqrt(pi)`
    (as `exp_r_c` does with its fixed 1.0s) looks like a strictly-better
    idea with a free degree of freedom given up. Measured: it **costs
    `erfc` a full max ulp** (0.1993/7 free vs 0.2003/8 pinned) because the
    other ten coefficients cannot re-absorb the constraint. The shipped
    `c0` is 1 ulp low and buys that ulp back; what it pays is 0.023 avg
    ulp in the far tail (`x >= 20` avg 0.6245 -> 0.6478, max 4 either
    way). Generalizes: pinning a coefficient to its mathematically exact
    value is a real constraint on the fit, not a free correctness win.
  - **What the residual is, so the next session does not chase it.**
    `erfcx_pos`'s own max is ~5.9 ulp and **the fit is ~0.5 of that**. The
    rest is forming `v`: rounding `2+xa` and then the reciprocal costs
    ~2.3 ulp no polynomial can recover, because `erfcx` has a nonzero
    slope at 0 while `v` is stationary in relative terms there.
    Compensating it needs the exact residual of `2+xa` (4-6 more ops) for
    ~1.5 ulp — measured, rejected on cost.
  - **A rejected alternative worth recording**, because it looks like the
    obvious upgrade: replacing `exp_checked(-p) * fma(-r,pe,r)` with a
    `Df32` exponent (`exp2_checked_df(Df32::from_mul(xs,xs) * -LOG2_E)`)
    made `erfc` **worse — max 23**, against 8 for the Cody-Waite form on
    the same coefficients at the time. Same lesson as `log2_df`: EFT
    bookkeeping around an f32 kernel is not higher precision, it is capped
    by that kernel's own evaluation rounding, and here plain Cody-Waite
    plus a one-fma exact-square correction beats it outright.
  - Two stale facts found and fixed along the way: `erfcx`'s negative-side
    overflow boundary is `|x| > 9.382` (`x^2 > 88.03`), not the `13.3` /
    `176.7` that `edgecheck.rs` claimed; and `worst_corpus.golden` was
    already stale at `3bad293` (42 entries across `cosh`/`exp`/`expm1`/
    `exp_m1_over_x`/`log_2`, all out-of-domain garbage from the exp
    reduction change, never re-blessed).

- **`erfc`/`erfcx` perf pass, all of it bit-identical** — a follow-up to the
  reciprocal-variable entry above, buying back cost without touching a
  single returned bit (`worst_corpus` passes *unblessed*, and the exhaustive
  sweep reproduces 0.1993/7 and 0.2142/6 exactly).
  - **The lever that worked was the guard, not the arithmetic.** Both
    functions clamp `|x|` before squaring it. Choosing that clamp so the
    *exponential's own* clamp is provably dead removes it entirely:
    `erfc` feeds `-xs^2`, so pick `ERFC_XS_CLAMP` with `xs^2` inside
    `[150*ln2, -EXP_CLAMP_LO]` -- above the low end `e^-p` rounds to
    exactly 0 (legal, the true `erfc` is already 0 by |x|~10.05), below
    the high end the reduction stays valid. `erfcx` feeds `+xs^2` and
    needs `2*e^(xs^2)` to still overflow to `+inf`, so its window is
    `[ln(f32::MAX/2), EXP_CLAMP_HI]`. **Four `const _: () = assert!`
    prove both windows at compile time** -- verified to fire in both
    directions, because a silent one-edit break here is a wrong-answer
    bug, not a slowdown.
  - `exp_checked`'s body became `exp_reduce!(x.clamp(LO, HI))`, the macro
    being the clamp-free reduction. Macro not fn, per `exp_r_poly!`'s own
    note; **verified by a full-file assembly diff showing byte-identical
    output**, so none of `exp_checked`'s 20+ other callers moved.
  - **Sign folding beat compare-and-select.** `erfc`'s
    `z = if x<0 {-1} else {1}; w = 1-z; fma(y,z,w)` became
    `w = from_bits((x.to_bits()>>1) & 0x4000_0000); mulsign(y,x) + w`.
    Bit-identical -- `fma(y,-1,2)` and `(-y)+2` are each one rounding of
    the same exact `2-y` -- and worth **-5.9% throughput** on its own,
    because it moves the work to integer ports the polynomial's fmas are
    not contending for.
  - Net: `erfc` **70.36 -> 62.28 cyc latency (-11.5%) and 3.284 -> 2.899
    cyc/elem throughput (-11.7%)**; `erfcx` 69.97 -> 66.99 (-4.3%) with
    its throughput column contested (see below). `erfc`'s latency is now
    *below* where it was before the accuracy rewrite (64.00), and its
    throughput regression against that baseline drops from +34.8% to
    +19.0%.
  - **`erfcx` is a clean worked example of mca's throughput column being
    wrong, arbitrated to the end.** The clamp removal deletes 2 `vminps`
    and adds 1 `vmovaps`; instructions 124 -> 123, uOps 143 -> 140, Block
    RThroughput 39 -> 38, latency -4.3%, fma/mul/div mix identical --
    every structural metric improves -- while simulated throughput reads
    **+3.7%**. That is the documented artifact signature ("identical
    expensive ops plus fewer total = mca is wrong"), so it was settled
    with the wall clock rather than argued: two binaries built and run
    **alternating A/B/A/B for four rounds**, and the clamp-removed
    version was faster on throughput in **4 of 4** (0.858/0.911/0.889/
    0.942 vs 0.952/0.952/0.984/0.953 ns/op) and on min latency. The
    published mca row is therefore *worse* than the previous one for a
    change that is really ~6% faster -- left as measured, flagged here.
    Alternating A/B on prebuilt binaries is what makes this machine's
    wall clock usable at all; the same session drifted from ~14% fast to
    ~15% slow on an unchanged control.
  - **Four things measured and rejected, all of them plausible:**
    - **Drop a polynomial degree.** The standing "a branch under the
      binding max can shed a term" lever does *not* apply here: the fit
      is the binding term, not headroom. Ideal (exact-arithmetic) fit
      error is 0.542 ulp at degree 10 but **>= 2.9 ulp at degree 9 for
      every one of 12 reciprocal offsets scanned** (A from 0.5 to 8), and
      through the real f32 chain max goes 4.049 -> 6.058 -> 10.400 for
      degree 10 -> 9 -> 8.
    - **Even/odd split in `v^2`** (11 ops vs Estrin's 12, depth 6 vs 5).
      Cheaper on *both* mca axes (`erfc` lat -2.0%, thr -3.9%) and it
      really is fewer operations -- but it costs **real accuracy**:
      exhaustive `erfc` 7 -> 8 and `erfcx` 6 -> 8, and still 8/8 after
      re-polishing the coefficients against the new rounding order.
      Estrin's shallower tree is worth a max ulp: same-grid standalone,
      4.670 vs 5.366. **Evaluation order is an accuracy decision here,
      not just a scheduling one.**
    - **Horner** (10 fma, no `v2`/`v4` at all): throughput -6.5%, but
      latency **+13.5%** from the 10-deep chain. Even/odd with Estrin
      halves was worse on every axis.
    - **Bitwise select for `erfcx`'s x<0 arm** (mask `g`, `mulsign` `r`,
      add): latency -5.5%, but instructions *and* uOps both **up** and
      throughput +1.9%. Unlike the `erfcx` clamp case above, here the
      structural metrics agree with the throughput column, so it is a
      real regression.
  - Stale fact corrected: `erfcx`'s doc claimed its mca latency number
    was untrustworthy because `mix()` masks the sign bit. `mix()` was
    fixed in idea #198 (`0x807f_ffff`); the note was removed.

---

## `tan_checked`'s 2.3e9 max ulp is `sin_checked`'s number, not a tan defect

Investigated 2026-08-01 as "the worst function in the crate". It is not a
`tan`-specific defect and it is not fixable inside `tan_checked`.

- **The headline numbers are the same number.** `accuracy quick` reports
  `sin_checked (all f32)` at **avg 3.1e8 / max 2.13e9** and `tan_checked`
  at avg 3.3e8 / max 2.34e9. `tan_checked` is `sin_checked/cos_checked`,
  so it simply inherits the shared `reduce_pi` double-float `q` running
  out at `|x| ~ 2^48*pi`. Ranking `tan_checked` above everything else is
  an artifact of which row got quoted, not of tan being worse.
- **Where the reported max/avg actually live: `|x| > 1e19`.** Per-band
  (200k log-spaced samples/band, both signs, f64 `tan` reference):
  `<1e9` max 3.3; `[1e9,1e12)` 535; `[1e12,1e13)` 2.0e4;
  `[1e13,1e14)` 3.6e4; `[1e14,8.9e14)` 1.0e11; `[8.9e14,1e16)` 6.6e12;
  `[1e16,1e19)` 2.2e12; `[1e19,3.4e38)` 1.1e12. Exponents 63..128 are
  ~25% of a bit-pattern-uniform fuzz, so **any change that leaves
  `|x| > 1e19` alone cannot move the reported max or avg**, however much
  it improves the bands below. That is the trap this entry exists to
  document.
- **The reduction is not "wrong", it is out of the poly's domain.** `q =
  qh + ql` is a sum of two integer-valued f32s, so it is an *exact*
  integer even when it is the *wrong* integer, and `reduce_pi` then
  computes `x - q*pi` honestly. Verified against exact-rational pi at
  `x = 8.2815029e14`: shipped `rs = -3.241658e-3` vs exact
  `-3.241631e-3`. The sin side is fine there; what fails is that `|r|`
  drifts past `pi/2`.

### Genuinely new: `cos_checked` dies a full binade before `sin_checked`

`round_x_over_pi::<true>` folds cos's `-0.5` into `lo`, then rounds
`rem = (p0 - qh) + lo` to an integer, then `reduce_pi(x, kh, kl + 0.5)`.
Once `|rem| >= 2^23` the ulp of `rem` is `>= 1`, so **all three of**
`lo - 0.5`, `rem.round_ties_even()` and `kl + 0.5` silently drop the
half-integer. That happens at `|x| ~ 2^47*pi = 4.4e14` -- *half* the
`2^48*pi = 8.85e14` limit `sin_checked`'s doc comment quotes.

Measured, same bands: `[1e14,8.9e14)` `sin_checked` max 2.3e8 / avg 9.1e3
vs `cos_checked` max **2.7e11** / avg 1.6e7, ~1000x worse. The mechanism
is visible directly: `reduce_pi_half_checked(x)` and
`reduce_pi_checked(x)` return the *identical* residual there (dumped at
`x = 6.2e14`, `2.43e15`, `3.36e15`), so `tan_checked` degenerates to
exactly `+-1.0`. This is the concrete reason behind the already-recorded
"sin_checked cliff earlier than documented" observation.

No cheap fix: `q_cos = k + 0.5` needs `2*q_cos` (an *odd* integer) exact,
which a two-f32 pair caps at the same place. Deriving the cos residual
from the sin residual afterwards is exactly what idea 50b rejected.

### Five rewrites measured; none dominates

All keep `tan(x) = tan(r)` for `r = x - q*pi` with **any** integer `q`
(period pi), so the second-stage fold costs *no* parity bookkeeping at
all -- that part works and is cheap.

1. **Re-reduce both `reduce_pi_checked`/`reduce_pi_half_checked`
   residuals mod pi.** No effect (1.04e11 -> 6.9e10). Useless because the
   cos residual has already collapsed onto the sin one; there is nothing
   left to fold back.
2. **Single f32 residual + derived cos** (`v = a - sgn*pi/2`, idea 50b's
   shape). `[1e14,1e16)` 1.0e11 -> **1.3e5**, but `[1,1e3)` 3.3 -> 1.5e5
   and `[1e19,)` 1.1e12 -> 2.5e19. Confirms 50b: a single-f32 `r` cannot
   carry the near-pole information.
3. **Single reduction, unnormalized double-float residual** (`reduce_pi`
   returning `(s3, err)` instead of `s3 + err`). *Worse* than (2) in
   places and `inf` for `|x| > 1e19`. Root cause worth remembering:
   **`err` is not a low word.** At `x = 6.465661e12` the pair is
   `rh = 1.59375, rl = -2.2952e-2` -- `rl` is ~200000x larger than
   `ulp(rh)`. `(rh - pi/2) + rl` then cancels catastrophically.
4. **Same + one `quick_two_sum(rh, rl)` renormalization.** The fix for
   (3), and the best variant found. Bit-identical to shipped
   `tan_checked` below 1e9 (max 3.3); `[1e14,8.9e14)` 1.0e11 ->
   **9.5e4**; `[8.9e14,1e16)` 6.6e12 -> **4.9e5**; `[1e16,1e19)` 2.2e12
   -> **7.6e8**. But `[1e9,1e12)` 535 -> 1.97e4 and `[1e13,1e14)` 3.6e4
   -> 9.6e5, and `[1e19,)` is a wash. So it trades a 30x near-pole
   regression over four decades for a 6-orders-of-magnitude win over four
   *other* decades, and **moves neither reported number**.
   - The residual 30x gap is `err`'s own accuracy: it is a plain 4-term
     f32 sum (`e1b + e2b + e3b - e3t`) tuned so that `s3 + err` rounds
     correctly, not to be a correct low word in its own right. A
     genuinely normalized double-float `reduce_pi` would close it -- and
     would also drop one whole `round_x_over_pi`+`reduce_pi` pair and
     both `parity` calls, so it is the one direction here that could win
     on accuracy *and* throughput. It needs the `sin_checked` domain.
5. **Hybrid** (shipped path below `2^47*pi`, variant 4 above). Best
   accuracy profile of the lot and the thresholds are principled rather
   than fitted, but ~12 extra ops on the crate's 2nd most expensive
   function to move a number nobody measures. Not built.

**Verdict: closed.** The reported max/avg cannot be improved without a
table-indexed Payne-Hanek (bits of `1/pi` selected by `x`'s exponent),
because past `~2^73` even a perfect `q` leaves `q*pi` short: the shipped
`PI_HI/PI_LO/PI_TINY` triple is only ~72 bits of pi, so the residual's
absolute error is `~|x| * 2^-73.65` regardless of how `q` is stored.
Three-word `q` alone buys nothing without four-word pi. And a table
lookup is a gather, which is the one thing this crate's auto-vectorized
scalar style cannot absorb.

## `erfc`/`exp` composites: the *argument's* rounding, not the kernel's error

`gelu(x) = x*0.5*erfc(-x/sqrt2)` scored 199 max ulp exhaustively with its
worst point at `x = -13.09`, deep in the tail, while `erfc` itself scores 7.
Its doc comment (and `norm_cdf`'s, and `norm_pdf`'s) blamed the kernel --
"inherits `erfc`'s own documented tail accuracy as-is". **That was wrong.**

The error is the composite's own *argument* rounding, amplified by the
kernel's condition number. `erfc`'s relative sensitivity to its argument,
`|z * dln(erfc)/dz| = 2z/(sqrt(pi)*erfcx(z))`, grows as `2z^2`: at
`z = 9.26` (`gelu(-13.09)`) the half-ulp already sitting in
`fl(-x*FRAC_1_SQRT_2)` comes back out as ~170 half-ulps of result. Nothing
`erfc` could do would fix that -- it is never handed the right argument.

Attribution, exhaustive over `gelu`'s negative octaves, splitting the total
into "f64 `erfc` at the *rounded* argument vs at the true argument"
(argument rounding) and "f32 `erfc` vs f64 `erfc` at the *same* rounded
argument" (kernel):

| `|x|` octave | total ulp | from argument | from `erfc` |
|---|---|---|---|
| `2^-1` |   7 |   1 | 6 |
| `2^0`  |  11 |   4 | 7 |
| `2^1`  |  22 |  16 | 6 |
| `2^2`  |  67 |  64 | 3 |
| `2^3`  | 199 | 197 | 2 |

**Fix (shipped for `gelu`):** hand `erfc` the argument as a double-`f32`.
`RSQRT2_HI + RSQRT2_LO` names `1/sqrt2` to a relative `2^-49`, one `fma`
recovers the product's exact residual `dz`, and the first-order term
`erfc(z+dz) = erfc(z)*(1 - 2*z*dz)` puts the lost bits back. Exhaustive
`gelu` **199 -> 10 max ulp, avg 0.3477 -> 0.2095**, and the worst point
moves out of the tail entirely (`x = -13.09` -> `-2.78`, i.e. what is left
is `erfc`'s own 7 plus composition rounding). llvm-mca cost is real but
small: latency 70.61 -> 74.85 (+6.0%), throughput 3.217 -> 3.470 (+7.9%),
throughput region 138 -> 153 instructions. Not worth a `gelu_fast` pareto
sibling at a 7% spread against a 20x accuracy gap.

Three details that are load-bearing, not incidental:

- **`2z` is only the *asymptotic* log-derivative, and it is wrong-signed
  nonsense for `x > 0`** (where `z < 0`, `erfc -> 2`, and the true
  sensitivity decays like `exp(-z^2)`). Applying the correction unclamped
  *adds* ~84 ulp at `x = +13`. `np = max(-x, 0)` gates it off; because
  `np*RSQRT2_HI` is then exactly `max(z, 0)`, the clamp costs one `vmaxps`
  and no extra multiply.
- **`max(-x,0)` also keeps `dz` finite at `x = +inf`**, where the
  unclamped `fma(nx, HI, -z)` is `inf - inf = NaN`.
- **Apply the correction to `erfc`'s result, not to `x*Phi(x)`.** The
  natural-looking `fma(-base, t, base)` is `NaN` at `x = +inf` even with
  `t` exactly `0`, because `base` is `inf` and `inf*0` is `NaN`. Folding it
  into the bounded `erfc` value first sidesteps the whole indeterminate
  form and needs no extra branch.

`norm_cdf` (186) and `norm_pdf` (65) are the same defect and are open --
`norm_pdf`'s is `-0.5*x*x`, whose rounding `exp` amplifies by `|x^2/2|`
(~85 at `x = 13`), fixable with `two_prod(x,x)` and a `(1 - 0.5*e)` factor.

---

## `norm_cdf`/`norm_pdf`: the composite's *argument* was the whole error

Both were thin composites whose max ulp had nothing to do with the
primitive they composed. Shipped 2026-08-01; exhaustive `norm_cdf`
**295 -> 8**, `norm_pdf` **67 -> 4**.

- **The screen, worth reusing on any composite:** price the rounding of
  the *argument* through the outer function's own log-derivative before
  touching anything else. `norm_cdf(x) = 0.5*erfc(-x/sqrt(2))` rounds
  `x/sqrt(2)` to `2^-24` relative; `erfc`'s dominant `e^(-z^2)` factor
  turns that into `2*z^2 * 2^-24` relative on the result, which at
  `x = -12.8` (`z^2 = 80`) is ~160 ulp on its own. Measured max was 188.
  `norm_pdf`'s `-0.5*x*x` is one rounding of magnitude `ulp(x^2/2)/2`
  landing *in an exponent*, i.e. that many ulp directly: ~70 at
  `|x| ~ 13`, measured 65. In both cases **no amount of work on `erfc`
  or `exp` could have moved the number** -- and `erfc` had just been
  taken to 7 max ulp, which is exactly why this looked like a mystery.
- **The fix is to move the rounding, not to compensate it.** Square `x`
  first and halve (exact), instead of halving/dividing first and then
  squaring -- the same "order of operations decides the error, not the
  formula" shape as `tanpi`'s subtract-before-scale. `erfcx` keeps the
  `x/sqrt(2)` argument, where `d(ln erfcx)/dz ~ -1/z` leaves it a plain
  `2^-24`. The residual `pe = fma(h, xs, -p)` then compensates the one
  remaining rounding.
- **Cost, and it is nearly free.** `norm_cdf` latency 71.63 -> 70.08,
  Block RThroughput 41 -> 42 (+2.4%); `norm_pdf` latency flat, RTh
  23 -> 24 (+4.3%). Note mca's *throughput column* claimed `norm_cdf`
  -8.6% while instructions (132 -> 135), uOps (+3.3%) and RThroughput
  (+2.4%) all said "slightly up" -- another instance of the column being
  the unreliable one, this time optimistic rather than pessimistic.
- **`erfcx` vs `erfcx_pos`, a real 37% trap.** Writing the first version
  against the public `erfcx` cost **+37.6%** throughput (132 -> 198
  instructions, `vpslld` 4 -> 8 = *two* `exp2_field_split`s). `erfcx`'s
  `x < 0` arm computes a whole second `exp_reduce!`, and LLVM does not
  prove it dead even though the argument is literally `x.abs() * c`.
  Calling `erfcx_pos` directly recovered all of it. Generalises: when a
  public wrapper's other arm contains an `exp`/`log`/division, check the
  opcode histogram for a doubled expensive op before accepting its cost.

### Tooling note: `jm check` false-positives on any net insertion

`cmd_check` matches **post-image** diff hunk line numbers against
function ranges generated from **master's** `lib.rs`. An edit that adds
lines therefore reports every domain that follows it in the file. This
change (+52 net lines, confined to `norm_cdf`/`norm_pdf`) reported
`dawson`, `erfc_inv` and `probit` as "NOT YOURS". Verify by checking
which lines the diff *removes* -- those are in pre-image coordinates and
are not shifted.

## `dawson`: the crate's worst avg ulp was one coefficient, 1 ulp above 1.0

`dawson` held the crate's worst *average* error (0.805 ulp, ~5x the next
function) and 60 max. Both had a single, unglamorous cause each.

**avg: `pc[0]` was `1.0000001`, which is exactly `1.0 + 2^-23`.** The
central branch is `x * P(u)/Q(u)` with `u = x^2`, and for small `|x|` the
whole rational collapses to `P(0)/Q(0) = pc[0]/qc[0]`. So `dawson(x)`
returned `x * (1 + 2^-23)` -- which rounds to 1 or 2 ulp above `x` depending
on the mantissa, mean exactly 1.5 -- for *every* `x` below about `2^-13`.
Measured per-octave avg was a flat `1.5000` from `2^-14` all the way down.
That is ~54% of the harness's uniform-over-bit-patterns samples, and
`0.54 * 1.5 = 0.81` is the entire reported 0.805.

Pinning `pc[0] = 1.0` makes every octave below `2^-14` score exactly **0**.

This is the counter-example to the "pinning an exact coefficient costs"
entry (which is about `erfc`, and still true there). The distinguishing
question is not "is the value famous" but **"does the exact value make an
entire region of the domain exact?"** Here `P/Q` is then exactly `1.0` and
`x * 1.0` is exact, so a whole half of the domain becomes error-free; the
cost is one constant term the minimax would have placed within its own
~6 ulp error band of 1 anyway. It did cost ~0.4 avg ulp in `|x|` in
[0.125,0.5], which the refit below took back.

**max: the fit minimised *absolute* error while ulp is *relative*.**
`dawsn(x)/x` falls 30x across `u in [0,16]`, so a least-squares objective
spends its budget near `u=0` and leaves the top of the range ~30x worse in
ulp. Per-octave avg was 0.4 in the middle and **10.2** over `[2,4)`, with
the max sitting at the `|x|=4` seam. A relative-error minimax (linearised
`P - g*Q` residual, LP via HiGHS, 12 iterations, `pc[0]`/`qc[0]` pinned)
drops the idealised fit error to 6.5 ulp.

**Net, both changes, zero cost:** harness avg **0.805 -> 0.154**, max
**60 -> 15**; a dense 14.7M-point scan of the central branch against
`scipy.special.dawsn` says **61 -> 16**. llvm-mca *improves*: throughput
1.969 -> 1.938, and the throughput region goes 78 -> **77** instructions,
because `pc[0] = 1.0` now shares the `1.0` broadcast `qc[0]` already needed.
So this supersedes the "REJECT (axis trade)" verdict in the poly-headroom
table above (`dawson_central_ratio`, 7.2x headroom, "max 61 -> 13 for avg
0.806 -> 2.993"). **That rejection was real but was measuring the `pc[0]`
mechanism, not the refit**: an unpinned minimax moves `P(0)/Q(0)` a few ulp
off 1, which wrecks the small-`|x|` half of the domain and shows up as
exactly the +2.2 avg that entry reports. Pin `pc[0]` and the trade vanishes.

Two things measured and *not* shipped:

- **`dawson_tail_poly` relative-minimax refit** (same treatment, `c[0]`
  pinned): fit error 10.94 -> 2.59 ulp in exact arithmetic, but a
  **regression in the real chain** -- harness avg 0.1537 -> 0.1692, max
  15 -> 16. Minimax spreads error uniformly over `v in (0,1/16]`, and the
  shipped least-squares fit is far better than uniform at *small* `v`,
  which is where nearly all the tail branch's samples are (`|x| > 4` is
  ~50% of all f32, and almost all of it has `v ~ 0`). Per-octave avg
  `2^3` 0.80 -> 2.37 and `2^4` 0.64 -> 2.33 for a seam max of 12 -> 7 that
  **does not move the reported number**, because the central branch's 15
  is the binding max. Revisit only if the central branch ever drops below
  ~10, and then fit L1-under-a-max-cap rather than pure minimax.
- **avg/max pareto sweep of the central fit** (minimise weighted L1 subject
  to a hard cap on the max, cap swept 7..60 ulp): the frontier is **flat**
  -- L1 moves only 1.101 -> 1.305 as the cap tightens from 60 to 7 ulp. At
  this degree there is no avg to buy back by loosening the max, so plain
  minimax is the right corner. Do not re-run this sweep.

The accuracy harness's `dawson` reference was also checked rather than
trusted, since it is the one hand-rolled reference in the file: Simpson's
rule at N=800 has relative error <1e-11 out to `x=2`, 1.5e-8 (**0.13 f32
ulp**) at `x=4`, 8.8e-8 (0.74 ulp) at the `x=5` handover to the asymptotic
series. So it is a real reference and the 0.805 avg was never a reference
artifact -- but treat anything under ~1 ulp in `[4,5]` as noise. (The
separate, still-valid caveat is that `dawson` is capped at 2M samples even
in `thorough` mode, so its max is sampled, not exhaustive.)

## `srgb_to_linear`'s 13 max ulp is the log2/exp2 round trip, not its own algebra

Attributed rather than assumed, over a 4.9M-point scan of the power arm
(`c > 0.04045`) against an f64 reference, splitting the error into the
rounding of the argument `b = (c+0.055)/1.055` and everything downstream:

| source | avg ulp | max ulp |
|---|---|---|
| total | 1.496 | 13 |
| rounding of `b` | 1.670 | 6 |
| `log_2_normal` -> `*2.4` -> `exp2_checked` | 1.730 | 13 |

**At the worst point (`c = 0.04841894`, total 13) the `b` term contributes
exactly 0 and the chain contributes all 13.** So the gelu-style fix -- carry
`b` as a double-`f32` and put the lost bits back with a first-order
correction -- is available and would take the *avg* down, but it cannot move
the reported max, because the max is entirely the round trip.

The mechanism is `exp2` amplifying an *absolute* argument error: with
`a = log2(b)*2.4`, `exp2(a+da) = exp2(a)*(1 + da*ln2)`, so `log_2_normal`'s
own ~1 ulp of a `log2(b)` near `-3.34` (absolute ~2.4e-7) becomes
`2.4 * 2.4e-7 * ln2 = 4e-7` relative, ~3.3 ulp, before `exp2_checked`'s own
error and the `*2.4` rounding are counted. Nothing in `srgb_to_linear`
itself is wrong.

Closing it properly needs a genuinely more accurate `log2`, and that is the
known dead end: `log2_df` is *not* higher precision (see the "Df32
bookkeeping is not precision" entry -- it is capped by its own f32 kernel's
evaluation rounding at ~2^-23 relative), and `log_2_normal` is a shared
kernel behind 9+ public functions, so moving it is a core-lock change that
would have to be justified by all of them, not by sRGB. IDEAS.md's existing
note already reaches the same verdict from the other direction (the exponent
is a constant 2.4, so the amplification is bounded and the inputs are
`[0,1]`); this is the measurement behind it.

`linear_to_srgb` (7 max ulp) is the same chain with `1/2.4`, which is why it
is ~2x better: the exponent multiplies the log's absolute error directly.

### sin/cos/tan

- **sin/cos/tan: two-word `1/pi` in the magic round** (2026-08-01,
  **SHIPPED**): `sin` 219 -> **2**, `cos` 2762 -> **2**, `tan` 3019 ->
  **4** max ulp over the whole documented `|x| < 2^22*pi` domain, avg
  0.0594/0.2902/0.3231 -> **0.0422/0.0833/0.1177**. Verified
  exhaustively (every f32 bit pattern, both signs) over
  `[1e-6, 2^22*pi)`, not just fuzzed.
  - The diagnosis in the "Fast-tier reduction upgrade via two_prod"
    entry above is **wrong and should not be reused**: it says the
    ~1.3e7 cliff is "the magic-round producing a wrong `q`
    (`|x/pi| >= 2^22`)". The magic round is fine everywhere in-domain --
    `x*FRAC_1_PI` tops out at 4194303.63, under `2^22`. What is wrong is
    the *constant*: a single f32 `1/pi` is only `2^-25` relative, so
    `x*FRAC_1_PI` sits up to `|x|*2^-25/pi` = **0.17** away from `x/pi`
    at the top of the domain, and `q` lands a whole integer off whenever
    `x/pi`'s fraction is within 0.17 of a half. That is ~25% of inputs
    up there, not a corner.
  - Being one off is *self-consistent* (parity comes from the same `q`),
    so this is not a wrong-answer cliff -- it is `sinf_poly` being
    evaluated at `|r|` up to `pi/2 + 0.17*pi = 2.1` when it is a
    degree-9 minimax fitted on `[-pi/2, pi/2]`. Pure extrapolation cost.
  - Screen that settled it in one run: replace `q` with
    `((x as f64)/PI).round() as f32` and re-measure. 222 -> 2, 2780 -> 2,
    3057 -> 4 on a dense sweep of `[1.2e7, 2^22*pi)`. Everything else in
    the reduction is already exact enough -- the four `PI_A..D` words
    leave `pi` to `1.9e-22`, so `q*(pi - sum)` is `8e-16` at `q = 2^22`,
    a thousand times under the fma roundings.
  - **`cos` was twice as bad as `sin` for a second, separate reason**:
    `fma(x, FRAC_1_PI, -0.5) + ROUND_MAGIC` rounds *twice*. The fma's
    own result is quantized to `ulp(x/pi)`, already 0.5 at the top of the
    domain, so the magic round downstream is handed an exact tie and
    ties-to-even sends it the wrong way -- measured `k` a full 0.79 off
    at `x = 1.3176051e7`. The new form never materializes `x/pi - 0.5`.
  - Shipped shape: `n = round(x*RPI_HI)` (magic, on the *exact* product),
    `f = fma(x, FRAC_1_PI, -n)` recovers that product's fraction, and
    `fc = fma(x, RPI_LO, f)` folds in the second word of `1/pi`. `sin`
    then re-rounds with `qb = fc + nb` (a second magic round that lands
    the parity bit back in `qb`'s low mantissa bit for free -- 3
    instructions cheaper and 4 latency cycles better than
    `n + fc.round_ties_even()`, which needs its own parity XOR). `cos`
    needs no second round at all: the half-odd-integer nearest `x/pi` is
    just `n + copysign(0.5, fc)`.
  - `RPI_TINY` is not needed -- it contributes at most `2e-9` over the
    domain against `fc`'s own `~6e-8` rounding.
  - Cost, mca, ladder-checked (instructions / uOps / Block RThroughput,
    because the throughput column overstated all three): `sin`
    1.151 -> 1.776 cyc/elem, instructions +20%, uOps +22%, RThroughput
    14 -> 18 (+29%); latency 48 -> 64. `cos` 1.406 -> 1.654, instructions
    +18%, uOps +18%, RThroughput 16 -> 17 (+6%); latency 56 -> 61. `tan`
    2.532 -> 3.153, instructions +8%, uOps +8%, RThroughput 30 -> 31
    (+3%); latency 71 -> 78.
  - Worth it because the accurate alternative is far more expensive:
    `sin_checked` is 5.232 cyc/elem and `cos_checked` 4.603, so this
    reaches their accuracy over `sin`'s own domain at ~3x less cost.
    `sin` now matches `sin_checked` exactly on `|x| <= 1e6`
    (0.0357/2 vs 0.0356/2).
## `rootn`: 45 max / 10.1 avg ulp was `log2(|x|)` being an `f32` at all

`rootn(x, n)` was `exp2_checked(log_2(|x|) / (n as f32))`, and its error was
neither a fit problem nor a `log_2` problem. `log2(|x|)` ranges over
`[-149, 128]`, so *storing it in an `f32`* costs up to `ulp(149)/2 = 7.6e-6`
absolute -- and `exp2` turns an absolute error in its argument straight into
a relative one, `ln(2) * 7.6e-6 = 5.3e-6`, which is **44 ulp**. Dividing by
`n` first only scales that down proportionally, which is why the function
got *more* accurate as `|n|` grew and was worst at the small `|n|` a caller
is most likely to write. The measured shape matched exactly:

| n | avg ulp | max ulp |
|---|---|---|
| -1 | 10.10 | 45 |
| 3 / -3 | 5.07 / 4.95 | 43 / 44 |
| 2 / -2 | 2.71 / 2.60 | 44 / 42 |
| 7 / -7 | 2.16 / 2.15 | 20 / 20 |

The standing proposal (IDEAS.md) was a `Df32` route -- `log2_df(ax)` divided
by `n as f32` through a real double-float division. That would have worked,
but it is the wrong shape of fix: the value that needs more bits is not the
logarithm, it is the *exponent*, and the exponent is an integer.

**Fix: never form `log2(|x|)` as a float.** Split `|x| = m * 2^e` on
`log_2_normal`'s own `[sqrt(2)/2, sqrt(2))` window, then use the exact
Euclidean split `e = q*n + rr`, `0 <= rr < |n|`:

```text
log2(|x|)/n = q + (rr + log2(m))/n
```

`q` is an integer and goes straight into `exp2_kf`'s exponent field, so `e`
never touches the float path. The fractional argument is *self-normalising*:
`rr < |n|` bounds the sum `rr + log2(m)`'s own rounding at `|n| * 2^-25`, and
the division by `n` scales it right back to `2^-25` -- independent of both
`e` and `n`. Accuracy becomes flat in `n` instead of degrading as `|n|` falls.

Three details are load-bearing:

- **The `[sqrt(2)/2, sqrt(2))` window, not `frexp`'s `[0.5, 1)`.** It is
  `log_2_normal`'s own window, so `log_2_unchecked(m)` computes `k == 0` and
  never pays its `lm + k` rounding -- `log2(m)` comes back in `[-0.5, 0.5)`
  with full *relative* accuracy, so an exact power of two costs nothing.
  Doing the split by hand (`denormal_rescale!` plus the same two bit ops
  `log_2_normal` uses) rather than `frexp` + an octave shift also drops
  `frexp`'s zero/infinite selects, which the degenerate arm redoes anyway:
  37 instructions.
- **`|n| >= 2` on the general path is what makes `exp2_kf` legal.** Its
  contract is `k` in `[-126, 128)` and a normal result, and at `|n| >= 2`
  the result exponent is at most `~149/2`, so it is always normal --
  no saturating scale, no `ldexp`. `|n| == 1` is the one case that can
  overflow or go denormal, and it is excluded rather than accommodated.
- **The `|n| == 1` closed forms and the zero/infinite/NaN arm are the same
  select.** `|x|` for `n > 0` and `1/|x|` for `n < 0` is simultaneously the
  exact `|n| == 1` identity (correctly rounded by hardware over the whole
  domain, including the overflowing and denormal results the general path
  could not produce) *and* the right magnitude for zero, infinite and NaN
  `x` at any `n` (where only `sign(n)` matters). The existing sign/parity
  logic finishes both, so `n == 1 -> x` and `n == -1 -> 1.0/x` need no
  selects of their own. `n == 0` falls into the same arm and is discarded
  by the domain-error select that was already there.

Result: max ulp **45 -> 1** (2 on a few `n`), avg **10.10 -> 0.00** at
`n = -1` and **5.07 -> 0.15** at `n = 3`, uniformly across a 24-value `n`
sweep spanning `+-1` to `i32::MIN`/`i32::MAX`. Confirmed **exhaustively**
(all 2^32 bit patterns against an f64 reference) rather than left on the
fuzz, since a 2-arg function has no `thorough` mode and its fuzz max is
optimistic:

| n | avg ulp | max ulp | finite samples |
|---|---|---|---|
| 2 | 0.13738 | 1 | 2.14e9 |
| -2 | 0.13613 | 1 | 2.14e9 |
| 3 | 0.15088 | 2 | 4.28e9 |
| -3 | 0.17085 | 1 | 4.28e9 |

(The even-`n` avgs read ~2x the fuzz's because the harness divides by
`n_samples`, not by the number of *scored* samples, and negative `x` with
even `n` is a skipped domain error. Pre-existing, and it cancels in a
before/after comparison.)

Cost, by instruction count in
a standalone `--emit=asm` probe (`rootn` cannot be `llvm-mca`'d -- its
multi-exit branching corrupts the region markers, see its doc comment):
`rootn(x, 3)` 91 -> 101, `rootn(x, -3)` 91 -> 111 (the `1.0/|x|` arm stays
live for negative constant `n`), dynamic `n` 104 -> 143 (Rust's
`div_euclid`/`rem_euclid` lower to a real `idiv` when `n` is not a
constant; for a literal `n`, which is the normal call, LLVM
strength-reduces it away). A 25-45x accuracy gain for 11-22% instructions
at constant `n`.

Note what this does *not* transfer to: the trick needs the outer exponent to
be one the integer part survives. `1/n` is exact integer division;
`srgb_to_linear`'s `2.4` is not (`2.4*e` is not an integer), so it would need
a two-term split of `2.4*e` and is a different, non-free change.

- **`asind`'s 11 max ulp is NOT the `atan2d` denormal mechanism** -- the
  lead left in this file's own `atan2d` entry ("same shape: `asind`/
  `acosd`/`atand` all multiply by `180/pi` and `asin(x) ~ x` can be
  denormal too") is measured **false**. Instrumented over all 8388608
  f32 in `[0.2, 0.4)`, the band both maxima live in: **0 of the 45835
  samples scoring >= 4 ulp have a denormal `asin(x)`** -- not a small
  fraction, zero. `asind`'s worst x is `2.74e-1`, nowhere near
  underflow.
  - What it actually is: `asind = asin(x) * (180/pi)` inherits `asin`'s
    *relative* error, and the two land in binades with different
    ulp-per-unit-value. `x/ulp(x)` ranges over `[2^23, 2^24)` within a
    binade; `asin(0.274) = 0.2775` sits at `1.09*2^23` (ulp buys the
    most) while `asind(0.274) = 15.9` sits at `1.96*2^23` (ulp buys the
    least). Measured mean amplification over those 45835 samples:
    **1.8836**, against the pure binade-geometry prediction **1.7905**
    at the worst point. That is the whole factor: `asin` 6 -> `asind`
    11, and after the fma fix below, `asin` 5 -> `asind` 9.
  - So `asind`/`acosd`/`atand` have no defect of their own and nothing
    to fix in the composite; they are exactly as good as `asin`/`acos`/
    `atan` are, read on an unluckier ruler. Any future work belongs in
    the radian function.
- **`asin`'s big branch: contract `pi/2 - sqrt(1-a)*P(a)` into one
  `fma`** -- shipped, a win on *every* axis, and the cheapest thing in
  this whole section. The product is ~4.7x larger than the difference it
  feeds just above the crossover (`1.297` vs `0.274`), so rounding it to
  f32 first costs ~2 ulp of the result; `fma(-s, P, FRAC_PI_2)` rounds
  once at the result's own magnitude. Exhaustive: `asin` max 6 -> **5**,
  avg 0.0199 -> 0.0188; `asind` max 11 -> **9**. llvm-mca, full ladder,
  every rung agreeing: instructions 57 -> 55, uOps 62 -> 60, Block
  RThroughput 14.0 -> 13.0, throughput 0.968 -> **0.900** cyc/elem
  (-7.0%), latency 60.99 -> **56.74** (-7.0%); `asind` 1.064 -> 0.981
  and 68.99 -> 64.74.
  - Rust does not contract `a - b*c` into an fma without fast-math, so
    any `k - p*q` written literally in this crate is still a `vmulps`
    plus a `vsubps` and still rounds the product. **The accuracy win and
    the instruction saving are the same edit** -- worth sweeping for the
    pattern wherever a product is subtracted from a constant.
  - Reached via `two_prod` first (`(pi/2 - q) - e`), which measured the
    *identical* max 5 / avg 0.0188 -- of course it does, both compute
    `pi/2 - s*P` with a single final rounding -- but cost 61
    instructions and 1.150 cyc/elem (**+18.8%**). Same value, opposite
    verdict on cost. If an error-free transform is being used only to
    un-round one product that is then added to something, the fma is
    strictly the better spelling.
  - This does not touch the diagnosis in entry 110 below (`asin`'s error
    is a broad plateau over `[0.27, 0.5)`, not a seam artifact) -- it
    removes one of the ~4 ulp-equivalent terms feeding that plateau. The
    remaining ones are `asin_poly`'s own evaluation error and the
    `1.0 - a` / `sqrt` roundings, both still amplified ~4x.

- **logaddexp: precision budget for the ~1e4-ulp cancellation**
  (2026-08-01, analysis only, no code -- redirected to sin/cos/tan
  before implementing; recorded so the next attempt does not re-derive
  it). `logaddexp` returns `m + corr`, `corr = log1p(exp(-d))` in
  `(0, ln2]`, and cancels when `m` is negative and close to `-corr`,
  i.e. `m` in `[-ln2, 0)`. Fuzz max swings 1549..27783 run to run purely
  from how close a sample lands to the zero curve `e^a + e^b = 1`
  (2-arg sampling noise -- repeat before trusting either number).
  - **`softplus` and `logsigmoid` do not share the defect and need no
    fix**: `softplus(x) = max(x,0) + corr` has `max(x,0)` and `corr` both
    non-negative, so nothing cancels, and `ax = |x|` is exact besides.
    Only `logaddexp`'s `d = |a-b|` is inexact and only its `m` can be
    negative.
  - Ranked error budget, all relative to `|m|`: (1) `d = fl(m-n)`'s
    rounding, `d*2^-25`, up to `2.6e-6` at `d=87`; (2) `exp_narrow`'s own
    ~`2^-24`; (3) `log1p_unit`'s ~`0.1` ulp; (4) the final `m + corr`
    add, `2^-25`. Only (1) is cheap to remove (`two_sum(m, -n)` then
    `e -= e*dl`, ~7 ops) and it buys only `(d+1.7)/1.7` -- a factor 2.6
    at `d~3.6`, 50 at `d~87`, and it does *not* set the floor.
  - After that the floor is `exp`'s own `2^-24`, so **any real fix needs
    a `2^-40`-ish `exp`** -- max ulp is `2^(24-b) * |m|/|result|`, and
    the observed cancellation depth `|m|/|result| ~ 1e4..1e5` wants
    `b >= 37`. Confirms the "double-float `log1p_unit`, an `_accurate`
    tier" verdict in the two_sum entry above, and sharpens it: the
    expensive half is `exp`, not `log1p`.
  - Three reformulations were checked on paper and all reduce to the
    same requirement, so do not re-try them: `log1p(expm1(m) + exp(n))`
    (both terms are themselves `~|m|`, identical budget);
    `log1p(e) - log1p(expm1(-m))` via `2*atanh` (algebraically the same
    quantity); and a Newton step `corr = c0 + ((1+E) - exp(c0))/exp(c0)`
    (limited by `exp(c0)`'s *relative* accuracy, which is the thing being
    fixed). A `2^-42` `exp(r)` with `|r| <= ln2/2` needs double-float
    terms through `r^5`; a 16-entry `2^(j/16)` table cuts `|r|` to
    `ln2/32` and lets the f32 tail start at `r^3`, which is the cheap
    construction if someone builds it.

- **Same `fma` contraction applied to `asinpi`'s `0.5 - sqrt(1-a)*P(a)`**
  -- second and last site of this shape in lib.rs (a full grep of
  `k - p*q` over the file finds exactly three: `asin`, `asinpi`, and
  `acos`/`acospi`, whose `sqrt(1-a)*poly` is not subtracted from
  anything and so has nothing to contract). Exhaustive: max ulp
  **7 -> 5**, avg 0.2375 -> 0.2364; llvm-mca 58 -> 56 instructions,
  64 -> 62 uOps, Block RThroughput 14 -> 13, throughput 0.974 ->
  **0.901** (-7.5%), latency 60.99 -> **56.74** (-7.0%). The
  cancellation is sharper here than in `asin` (`0.5 - 0.412 = 0.088`
  just above the crossover, 5.7x, against asin's 4.7x), which is why the
  same edit buys 2 max ulp rather than 1.
  - Screened for the mirror pattern `p*q + k` at the same time (Rust
    contracts neither): only 3 non-comment sites survive in lib.rs and
    none is on an accuracy-critical path. The crate's explicit `fma`
    helper has already absorbed essentially all of them -- these two
    subtractions were the leftovers, presumably because `k - p*q` does
    not *look* like an fma the way `p*q + k` does.
## `coshm1`: the half-angle identity's error doubling, priced and replaced

`coshm1` was `2*sinh_checked(x/2)^2` and sat at **max 10 ulp** where every
other member of the hyperbolic family is 4. Its own doc comment blamed the
identity ("squaring roughly doubles `sinh_checked`'s relative error"), which
is the shape of claim this crate has repeatedly found to be false -- so it
was checked rather than inherited. Here it is **true**: scoring the
predictor `2 * relerr(sinh_checked(x/2))` in ulp of the result, per octave,
against the real end-to-end error reproduces it to three digits at every
octave from `2^-14` to `2^6` (e.g. `[2,4)`: predictor 1.6581 avg / 9 max,
actual 1.6587 / 9). So the identity was the whole story, and the fix had to
be a different identity, not a better `sinh`.

Replaced by two branches, each avoiding the other's error mechanism:

- **`|x| >= 2`: `cosh_checked(x) - 1.0`.** The naive form's cancellation is
  what the half-angle identity existed to dodge, but there is no
  cancellation out here -- and better, **subtracting 1 from a value `>= 1`
  is exact in binary floating point** (the result's exponent falls by at
  most one, so `1` always sits on the coarser grid). So this branch is
  exactly `cosh_checked`'s own error times the conditioning factor
  `cosh/(cosh-1)`: 1.36 at `x = 2`, 1.0007 by `x = 4`. Measured per octave:
  0.81 avg / 4 max on `[2,4)`, 0.72 / 3 above -- against the identity's
  1.66 / 9 and 1.52 / 9. Below `|x| = 2` the factor takes over fast (7 max
  on `[1,2)`, 18 on `[0.5,1)`, 10^8 by `2^-12`), which is what sets the
  handover at 2 and not lower.
- **`|x| < 2`: `x*(0.5*x) + (x^2)^2 * Q(x^2)`**, `Q` a degree-3 minimax
  (HiGHS LP) of `(cosh(sqrt(u)) - 1 - u/2)/u^2` on `u` in `[0,4]`, weighted
  by `u^2/(cosh(sqrt(u))-1)` so the fit minimises the *result*'s relative
  error. Both leading terms are peeled, not just the constant: the poly
  reaches the answer scaled by `u^2`, which is 28% of it at `x = 2` and
  vanishes as `x -> 0`. Result: **max 1 ulp** over the whole range
  `[2^-11, 1)` and 2 on `[1,2)`, against the identity's 3-7.

Two details worth keeping:

- **`x*(0.5*x)`, not `0.5*(x*x)`.** They differ only for results below the
  denormal floor, and there the first is one rounding and the second is
  two: `x*x` lands on the `2^e` denormal grid, and halving a denormal grid
  value is *not* exact (the halves are multiples of `2^-150`, unrepresentable).
  This is also exactly what the old `2.0*s*s` did -- `(2*s)*s = x*(x/2)` --
  so the tiny-`x` octaves stay bit-comparable rather than regressing.
- **Degree 4 was fitted, measured, and dropped.** Its residual is half
  degree 3's (0.088 vs 0.177 ulp) and it measures *identically* end to end
  (avg 0.0287, max 4, and per-octave to four digits), because the binding
  constraint is the final `fma`'s rounding plus `u = x*x`'s, not the fit.
  The "drop a term where the branch's max sits under the function's binding
  max" lever, worth 3 instructions.

**Harness (100M fuzz): avg 0.0832 -> 0.0287, max 10 -> 4.** Confirmed
exhaustively on the shipped version (all 2^32 patterns): avg **0.0287**,
max **4**, worst `x = 2.1662018`. (The old body's own recorded exhaustive
figures were 0.0864 / 12, so the honest before/after is either fuzz-to-fuzz
or 12 -> 4 exhaustive-to-exhaustive; the fuzz's 10 is it being optimistic
about a max, as usual.) The remaining max is `cosh_checked`'s own 4, reached
just past the handover -- the floor for any formulation that goes through
`cosh`, and not something `coshm1` can fix from inside its own domain.

Cost is real and was not clawed back: `coshm1_throughput` **97 -> 108**
instructions (+11.3%), `coshm1_latency` 3468 -> 3600 (+3.8%), and no other
region in the file moved by one instruction. That is the price of evaluating
both arms of the select where the old body evaluated one `sinh_checked`.
Taken as the single shipped version rather than a `_fast`/`_accurate` split:
11% is well under the ~13% spread that justifies the existing `sinh` tiers,
and 10 max ulp in a family where everything else is 4 was an outlier, not a
pareto point. Note `coshm1`'s **mca throughput row is mispriced in absolute
terms** (documented above), so instruction and opcode counts are the measure
here, not cycles.

## `pow_2_3`: `cbrt(x)^2` doubles cbrt's error before it rounds

`pow_2_3` was `cbrt(x)*cbrt(x)` and carried the crate's **highest average
ulp, 0.763** (max 7, exhaustive over all 2^32 patterns). Unlike most
high-average cases in this file, it was *not* a bad coefficient or a bad
region: the number is arithmetic. `cbrt` is 0.281 avg / 3 max, squaring
doubles a relative error, and the square's own rounding adds ~0.25 avg on
top -- 2*0.281 + 0.25 is 0.81, which is the measured 0.763 to within the
error terms' correlation. There was nothing to fix inside `cbrt`.

The replacement fits `x^(2/3)` **directly** on cbrt's own bit-trick seed.
`s * (1+r)^(-1/3) == a^(1/3)` identically for `r = (s^3-a)/a`, and the same
`r` gives `s^2 * (1+r)^(-2/3) == a^(2/3)` for free -- so the only change
needed is which exponent the correction poly fits. Exhaustive: avg
**0.7626 -> 0.1032**, max **7 -> 1**.

Three things the `-2/3` fit needs that the `-1/3` one does not, each
measured rather than assumed:

- **Degree 4, not 3.** `(1+r)^(-2/3)`'s series coefficients are ~3x
  `(1+r)^(-1/3)`'s, so over the seed's `r` range (`[-0.0999, 0.0894]`,
  enumerated exhaustively across all three exponent-mod-3 seed alignment
  classes) a degree-3 minimax fits to 2.03 ulp against cbrt's own ~0.5.
  Degree 4 fits to 0.125 ulp with f32-rounded coefficients.
- **The leading term's rounding.** `s` is a bit pattern and exact; `s*s`
  is not, and at full weight its `e2` is worth up to a whole ulp. Dropping
  the `e2` correction entirely measures avg 0.190 / max 2, so it is worth
  exactly the one instruction it costs.
- **`e2` again, through `r`.** `d = fma(s2, s, -a)` is short of `s^3 - a`
  by `e2*s ~ 2^-24 * a`. The obvious fix, a second `fma(e2, s, d)`, sits on
  the critical path and measured **+11.5% latency** on its own. It is not
  needed: that shortfall reaches the result multiplied by `-2/3` (cbrt's
  own equivalent is `-1/3`, which is why `cbrt_normal` gets away with
  ignoring it) and `s2*s/a ~ 1`, so it is `-(2/3)*e2` to well inside its
  own last bit. Folding it into the tail addend that already carries `e2`
  -- `e2 - (2/3)*e2 = e2/3` -- is one *multiply*, off the critical path,
  and reproduces the two-fma version's avg and max exactly (0.1032 / 1 vs
  0.1039 / 1).

**Horner, not Estrin, and the reason is codegen, not arithmetic.** With the
denormal rescale folded into the kernel as a parameter (worth 4 cycles of
tail latency), the Estrin schedule triggers exactly the trap `cbrt_normal`'s
doc comment already warns about: LLVM's vectorizer duplicates the whole
kernel for the tiny and normal branches instead of computing once and
blending. `pow_2_3_throughput` went 7600 -> **12200** instructions, uOps
8400 -> 12900, Block RThroughput 18 -> **34**. Horner does not trip it.
Estrin *without* the scale parameter is fine and is a real latency/
throughput tradeoff (+5.8% lat / +18.8% RThr against Horner's +11.5% /
+12.5%), but folding the scale in dominates both.

Cost of the shipped version, `pow_2_3_throughput` region: instructions
7500 -> 7600 (+1.3%), uOps 7900 -> 8400 (+6.3%), Block RThroughput 16 -> 18
(+12.5%), mca cycles 2651 -> 2682 (+1.2%); `pow_2_3_latency` 339506 ->
352804 cycles (+3.9%) at +0.1% instructions. FP mul/fma 24 -> 26. Shipped
as the single version rather than a `_fast`/`_accurate` split: 7 max ulp
was outside the crate's 0.5 avg / 2 max budget to begin with, so the old
body was not a pareto point.
- **`tan_poly` degree 5 -> 6, and the Estrin grouping is the whole
  decision.** `tanpi`'s max 8 / avg 0.386 was the crate's remaining
  accuracy outlier. Two screens ran before any fitting:
  - **The `PI*r` argument rounding is not the mechanism**, contrary to
    the obvious suspicion (`PI` as f32 sits 2.78e-8 relative from real
    pi, and the product rounds again). Instrumented over 20M samples in
    `[-4, 4]`, splitting the error into "exact tan of the f32 argument
    the chain actually forms" versus everything downstream: of the
    4680400 samples above 3 ulp, **1030 (0.02%)** have the argument term
    contributing even half. Argument-only avg is 0.44 ulp against a
    total of 2.05. So folding `pi` into the coefficients -- which would
    delete both effects and two `vmulps` -- cannot pay for a refit, and
    is not attempted.
  - **Oracle screen says the fit really is binding**, which is the
    unusual verdict here: running the chain with a correctly-rounded
    `tan(t)/t` in place of `Q` gives max 3.66 / avg 0.61 against degree
    5's 7.85 / 2.05. Contrast `ln_normal` and `asin_poly`, where the
    same screen found the headroom fake. An ulp-weighted LP reproduces
    the shipped degree-5 coefficients almost exactly (idealised 3.53 vs
    shipped 3.55), confirming they were already optimal *for that
    degree*; degree 6 idealises 0.354 and degree 7 0.065.
  - Degree 6 converts nearly to the oracle floor: real-chain max 4.28 /
    avg 0.71 against the oracle's 3.66 / 0.61. Degree 7 cannot convert
    and is not worth its fma.
  - **The grouping decides the cost, and the two options do not
    dominate each other.** Keeping the shipped Horner-on-`u^2` spine
    (`l0 + u2*(l1 + u2*l2)`) forces `l2` to be a full degree-2 group, so
    the new term lands *on* the spine: exhaustive max **4**, avg 0.2309,
    but llvm-mca latency **88.72** (+12.5% on the shipped degree-5's
    78.88) for +5 instructions. Regrouping to two halves
    (`lo + u^4*hi`, `u4` reused) costs one more instruction again but no
    fourth dependency level: exhaustive max **5**, avg **0.2267**,
    latency **76.00** -- *below* the shipped degree-5 -- throughput
    2.031 -> 2.223 (+9.4%).
  - **Shipped: the two-group form.** 17% of tanpi's latency is too much
    for one max ulp when the two-group form's *average* is the better of
    the two anyway. Exhaustive: `tanpi` 0.3856/8 -> **0.2267/5**,
    `tan2pi` -> 0.2273/5. The 4-ulp variant is recorded here rather than
    shipped as a second tier: it is worse on avg, worse on latency, and
    better on max by one -- too thin a pareto point to carry an API.
  - Reusable: this is the first poly in this file where a **degree bump
    was the right answer**, and the reason it worked is exactly what the
    other bumps lacked -- `tan_poly` has two callers (`tanpi`, `tan2pi`)
    and no shared-kernel blast radius, so the extra fma is paid only by
    the function that needed it. Check the caller list *before* pricing
    a degree bump, not after.

## `asind`: `180.0 / f32::consts::PI` is one ulp low, and that was its whole average

`asind`'s max ulp is not its own -- it is `asin`'s relative error read on
an unluckier binade ruler (`asin`'s `[0.25,0.5)` against degrees'
`[8,16)`, ~1.8x), and that has been separately measured and recorded. Its
**average** turned out to be a different and entirely separate defect,
which the max-side analysis correctly does not cover.

`asin(x) * (180.0 / std::f32::consts::PI)` computes `180/fl(pi)`, not
`fl(180/pi)`. Those are different f32s: `0x42652ee0` against `0x42652ee1`.
The shipped one is a **-0.4606 ulp** relative error on the multiplier;
the correctly rounded one is +0.0979. A constant relative bias on the
last multiply flips a fixed fraction of results by exactly 1 ulp, which
is why it shows up as an average and never as a max.

Decomposed over 4M random in-domain samples, feeding a *correctly
rounded* `asin` so only the composite's own arithmetic is in view (in-domain
average; the harness reports ~0.496x of these because `|x| > 1` is half of
all bit patterns and scores 0):

| final step | avg | max |
|---|---|---|
| `* (180.0/PI_f32)` (shipped) | 0.6712 | 2 |
| `* fl(180/pi)` | 0.1538 | 2 |
| `fma(y, HI, y*LO)` two-word | 0.0252 | 1 |

0.6712 * 0.496 = 0.333 against the harness's measured 0.3378 -- the
constant is the entire average, with nothing left over for `asin`.

Shipped: the two-word form. Exhaustive over all 2^32 patterns, avg
**0.3378 -> 0.0222**, max 9 -> 9 (unchanged, and expected to be: it is
`asin`'s own max at `x` just above the 0.27 branch crossover, 0.27002
before and 0.27009 after). Cost `asind_throughput` 58 -> 61 instructions,
uOps 63 -> 66, Block RThroughput 14 -> 15; `asind_latency` +6.4% cycles.
The correctly-rounded single constant is *free* and gets avg 0.085, a
real option if 3 instructions ever matter more than 4x the average --
but the two-word form is what makes the contract clean: `asind`'s error
becomes exactly `asin`'s error on a different ruler, with nothing of its
own added.

**The split point is not free, for a reason that has nothing to do with
accuracy.** Taking HI as the correctly-rounded `fl(180/pi)` makes LO
negative, and then `y*LO` is `+0.0` for `y = -0.0`, so the `fma`'s
`-0.0 + 0.0` returns `+0.0` and `asind(-0.0)` loses its sign. Caught by
edgecheck, invisible to a 100M fuzz. Taking HI as the *low* neighbour
(which is the old, "wrong" constant) makes LO positive and both terms
`-0.0`. Either HI represents 180/pi to the same ~2^-49 once LO is added,
so this costs nothing.

`1/pi` and `pi/180` were checked for the same defect and do not have it:
`1/fl(pi)` and `fl(pi)/180` each land on the correctly-rounded constant.
So this is specific to 180/pi, not a general rule about deriving from
`PI`. **`acosd`, `atand` and `atan2d` all still spell it
`180.0 / std::f32::consts::PI`** and carry the same -0.46 ulp bias.
## `log1pmx`: 11 max ulp was the *reference*, and the fix was one shared poly

The recorded 11 max ulp at `x = 9.53e-7` was never `log1pmx`'s. The harness
reference switched to a rationalized form below `|v| < 1e-6` and that form
was the *leading term alone*, `-v^2/2`. Dropping `v^3/3` is `(2/3)|v|` in
relative terms, which at `v` just under the cutoff is **10.6 f32 ulp** --
computed exactly against a 60-digit series at the recorded worst `x`, and
that is the whole 11. Same class as `wrap_pi`'s own TAU-carried-as-two-words
note: the reference was the less accurate of the two things being compared.
The direct `log1p_u10(v) - v` arm it was protecting only cancels below
`|v| ~ 1e-7`, so the cutoff was two decades too conservative *and* the
replacement too short. Reference now carries the series to `v^7` and hands
over at `1e-3`, where both sides are under 1e-5 ulp.

Against a correct reference the real number was **8**, at `x = 5.01e-1` --
one ulp past the `|x| < 0.5` poly boundary, on the direct arm. Three things
were then measured, in order:

- **Reorder alone: no.** `(ln(u) - x) + corr` instead of `(ln(u) + corr) - x`
  is free and strictly better conditioned (`ln(u)` and `x` are within a
  factor of two over most of the arm, so the difference is Sterbenz-exact
  and `corr`'s rounding lands at the result's own scale). It moved the worst
  point from 5.42e-1 to 5.01e-1 and **left the max at 8**: the binding term
  was never `corr`, it was `ln_normal`'s own few ulp amplified by
  `|x / log1pmx(x)|`, which peaks at ~5.3 exactly there.
- **Widening the poly: priced and not needed.** The degree-12 fit is a tight
  minimax over exactly `[-0.5, 0.5]` -- 0.53 ulp-equivalent at the boundary,
  32 at 0.55, 291 at 0.6. Reaching `x = 1` (where the amplification is 3.26
  instead of 5.29) needs degree ~16 by the Bernstein-ellipse estimate for
  the `x = -1` singularity, i.e. +4 fma for 8 -> ~5.
- **Reusing the same poly on `ln`'s own reduced argument: yes, and it is
  cheaper.** `u = 1+x = 2^k * m` with `m` in `[1/sqrt2, sqrt2)`, so
  `w = m - 1` is exact *and* already inside the fitted `|z| < 0.5`. Then
  `ln(u) = k*ln2 + w + log1pmx(w)` rebuilds the answer from the series this
  function already evaluates, and the whole `ln` polynomial disappears. The
  large terms difference before anything small is added, and over the band
  where the amplification is worst both partial differences are
  Sterbenz-exact. **8 -> 2 max.**

Summation order inside that arm is a real axis, not bookkeeping. Three were
measured: `(t + p) + last` 0.0648 avg / 2 max, `t + (p + last)` **0.0631 /
2**, and small-terms-first `fma(k,LN2_HI,-x) + ((w + p) + last)` 0.0594 /
**4**. The last one wins the average and loses the max, because grouping
`w + p` rounds at 0.29's scale where the result is 0.095. Balanced tree
shipped. Folding `p` into the arm with `fma(h, q, t)` also measured (0.0648
/ 2): it removes `p`'s own rounding but costs one *more* rounding at the
result's scale in the far field, where nearly all the sample mass is.

The rewrite put the degree-12 chain on the critical path (`z` now depends on
`u`), so it was Estrined -- and *where* the split goes is the whole result:

| poly form | avg | max |
|---|---|---|
| Horner, 12 deep | 0.0631 | 2 |
| full Estrin | 0.0685 | 4 |
| Estrin from `C2`, `1 + C1*z` peeled | 0.0668 | 3 |
| Estrin from `C1`, `1.0` peeled | **0.0632** | **2** |

Only the last is free. Grouping the `1.0` in puts two or more roundings at
full weight; peeled, there is exactly one, and everything inside the grouped
part enters attenuated by `z`. Same lever as `log_2`'s own peeled `c0*s`.

Net, all three axes. Exhaustive (all 2^32 patterns, both sides, against the
corrected reference): max **8 -> 3**, avg 0.0541 -> 0.0632 -- the average is
the one real cost, and it is the far field, where the old form reached the
answer in one rounding at `x`'s scale and this one takes two.
`log1pmx_throughput` **151 -> 133** instructions (-12%), uOps 17100 -> 15000,
Block RThroughput 37 -> 29 (-22%); `log1pmx_latency` 4991 -> 4603 instructions
(-8%), total cycles 387804 -> 371712 (-4.1%). No pareto variant needed --
cost fell on both axes. (The 100M fuzz reports 2 for the shipped version and
8 for the old one, so fuzz-to-fuzz is 8 -> 2; the exhaustive 3 is the honest
number, fuzz being optimistic about a max as usual.)

## `rcbrt`: two divisions to none, by putting the reciprocal in the seed

`1.0 / cbrt(x)` looked like the free composition -- "division is already
correctly rounded, just compose", the same reasoning `rsqrt`/`rhypot` use --
and it does hand every special case over for nothing. But it pays *two*
divisions, not one: `cbrt_normal` needs its own `1.0/a` to form the residual
`r = (s^3 - a)/a`, and then the composition divides again. And it pays
`cbrt`'s error in the wrong units: a relative error costs up to twice as many
ulp after a reciprocal, since the same relative size can sit at either end of
a binade. Hence 5 max against `cbrt`'s own 3.

Both go away if the bit-trick seed goes *down* the exponent instead of up.
`t = from_bits(K - ax/3)` approximates `a^(-1/3)` directly, and then
`e = a*t^3 - 1` is an fma, not a division -- because `a*t^3` is already O(1),
which is exactly what `cbrt_normal` needed the reciprocal for. The correction
is the same `(1+e)^(-1/3)` shape, so the *structure* transfers verbatim; only
the coefficients don't.

- `K = 0x54a20d0e`, from the `(4/3)*(127-sigma)*2^23` derivation then swept.
  Residual range is **exactly** `[-0.100965, 0.103657]` over all normals, not
  sampled: `e` is exactly periodic in the exponent with period 3 (adding
  `3<<23` to `ax` adds exactly `1<<23` to `ax/3`, so `a` scales by 8, `t` by
  2, and `a*t^3` not at all), so three consecutive exponents times all 2^23
  mantissas is the whole domain. Worth knowing generally -- it makes any
  `ax/3` seed exhaustively checkable in 25M evaluations instead of 2^31.
- **Reusing `cbrt_normal`'s four coefficients does not work**: on this range
  their idealized relative error is 11 ulp-equivalent (it is 3.7 on
  `cbrt`'s own `[-0.0998, 0.0893]`). The ranges look similar and are not.
- A degree-3 correction refit by minimax-LP for this range measured **max 2,
  avg 0.768** -- max fixed, average nearly doubled. Its own idealized floor
  is 1.6 ulp-equivalent, and a minimax error curve equioscillates, so most
  of the domain sits near that bound. The average *is* the fit here, and one
  more coefficient collapses it: degree 4 has an idealized floor of 0.19
  ulp-equivalent, six times lower, which is well under the evaluation
  roundings and therefore invisible. Shipped least-squares at degree 4.
- Specials, which the composition used to get free: `EXPONENT_MASK - ax` is
  the reciprocal's own magnitude on exactly the inputs that need one
  (`0 -> inf`, `inf -> 0`) and wraps into the NaN encodings for a NaN input,
  so one subtract plus `x`'s sign bit covers all of them with no division.

**Exhaustive: max 5 -> 1, avg 0.418 -> 0.117.** Cost is a wash, and the
ladder disagrees with itself, so all of it: `rcbrt_throughput` 75 -> 76
instructions, uOps 7900 -> 8400 (+6%), **Block RThroughput 20 -> 15 (-25%)**,
simulated total cycles 2666 -> 2704 (+1.4%); `rcbrt_latency` 2505 -> 2817
instructions, total cycles 384206 -> 390106 (+1.5%). The `-25%` and the `+6%`
are the same fact seen twice: four `vdivps` leave (`ICXFPDivider` was 17.7%
of the throughput bottleneck and is now absent) and are replaced by cheaper
integer work that lands on port 5. Taken as one shipped version, not a
tier -- a wash on cost is not a pareto point.
## `wrap_pi`'s 41 max ulp was the harness's own reference, not `wrap_pi`

Exhaustive (`|x| <= 1e4`, every f32 bit pattern, both signs) reported avg
0.0244 / max **41** at `x = 8953.539`, and had done since the first sweep.
The comment above it in `examples/accuracy.rs` wrote that off as "the
crate's usual near-a-true-zero artifact", which is exactly the excuse the
`cospi` postmortem says must be verified. Verified: it is false in both
halves. Against 60-digit decimal arithmetic at the four worst points:

| x | wrap_pi | old f64 reference |
|---|---|---|
| 8953.5390625 | **0.029** ulp | 44.168 ulp |
| 1011.5928344726562 | **0.369** ulp | 10.007 ulp |
| 3185.574951171875 | **0.352** ulp | 6.864 ulp |
| 7231.9462890625 | **0.754** ulp | 7.485 ulp |

`wrap_pi` was already the *more* accurate of the two things being compared.

The reference's defect is not the one its own comment guarded against. It
already carried `tau` as two f64 words -- but two words fixes the
*truncation* of tau, and the error here is the *rounding of the product*:
`q * tau_hi` with `q ~ 1425` and `x ~ 1e4` rounds at `ulp_f64(1e4)/2 =
9.1e-13`, and the answer at those points is `~2e-7` formed by cancelling
operands of `~1e4`, so 9.1e-13 is ~64 ulp of it. Carrying a third word
would not have helped either; nothing short of an exact product does.

**Fix: Cody-Waite the reference, don't just lengthen it.** Three words with
their low mantissa bits cleared -- `TAU1` a multiple of `2^-30` (33
significant bits), `TAU2` a multiple of `2^-58` (27 bits), `TAU3` the
remainder -- make every `q * word` exact in f64 for the `|q| <= 1592` this
domain can produce, so `((x - q*TAU1) - q*TAU2) - q*TAU3` carries only
relative-`2^-53` roundings of a result that never cancels again. Verified
to ~2^-29 ulp against decimal. Reported: **0.0244 / 41 -> 0.0244 / 1**.

Generalizable: a multi-word constant in a *reference* is only worth its
words if the products are exact. Any reference of the form
`x - round(x/c)*c` over a wide `x` range needs the split chosen against
`max|q|`, not against how many digits of `c` are written down.

### Two `wrap_pi` improvements measured and not shipped

With a correct reference the remaining error is real but small, and both
candidates cost more than a sub-ulp is worth. Note the "avg ulp" column
counts *integer* bit distance, so it reads as a misrounding rate.

- **Two-word `pi` in the half-turn fold.** For odd `q` the answer is
  `r -+ pi`, which always lands in `[pi/2, pi]`, and `std::f32::consts::PI`
  is 8.7e-8 above `pi` -- **0.73 ulp of the result, one-signed across that
  entire branch**, which is the whole of the residual max (1.233 ulp
  fractional, at `|result| ~ 1.5711`). Two forms measured, exhaustively:
  - `quick_two_sum(-ph, r)` then `s + (e - pl)` (one final rounding):
    avg **0.0244 -> 0.0043**, fractional max **1.233 -> 0.766**. mca
    throughput 4.281 -> 4.687 cyc/elem (+9.5%, instrs 134->147, uOps
    144->157, BlockRT 51->55 all agreeing), latency 92.02 -> **115.02
    (+25%)**. The latency is not the known branch artifact -- the region
    has *fewer* branches after (128 vs 256); it is 5 genuinely serial FP
    ops on the tail.
  - Plain `(r - ph) - pl`: avg 0.0168, latency 98.02 (+6.5%). Strictly
    the worse trade -- 30% of the accuracy for 26% of the cost.
  - Not shipped, and not shipped as a second tier either: the *reported*
    max is 1 either way, both forms are inside the crate's 0.5/2 budget,
    and `0.766` vs `1.233` ulp is not an API's worth of difference. The
    floor is anyway `r`'s own rounding as an f32 (~0.27 ulp transported
    into the fold), so even a perfect fold cannot reach correct rounding
    without the fold moving *inside* `reduce_pi`, which is a core-lock
    change on sin/cos's shared kernel.

- **Reduce `x/2` mod pi instead of `x`, and double (REJECTED, broken).**
  Superficially strictly better: `round_x_over_pi`'s `q` becomes
  `round(x/2pi)`, so `2*r` is the wrapped angle outright -- no fold, no
  `pi` to represent, no parity sign to combine, and both scalings are
  exact powers of two. It deletes ops rather than adding them. It is
  also **wrong**, and the reason is worth keeping:
  - The decision boundary of `round(x/2pi)` sits exactly on `wrap_pi`'s
    own `+pi`/`-pi` discontinuity (odd multiples of `pi`). `round_x_over_pi`
    resolves a near-tie through `rem = (p0 - qh) + lo`, an f32 at
    magnitude 0.5 -- so once `|x/2pi - (k+1/2)|` drops below `2^-25`,
    `rem` rounds to exactly `+-0.5`, `round_ties_even` returns 0, and `q`
    keeps whatever ties-away gave it. Wrong side of the cut = wrong
    answer by `2*pi`. Exhaustive: avg **5.50**, max **2.16e9** (i.e.
    `+pi` returned where `-pi` was correct) at `x = 47.12389 ~ 15*pi`;
    a handful of points, but each worth ~2^31 ulp.
  - The shipped mod-`pi` structure is immune *by construction*, and this
    is the part to remember: its own rounding ties fall at odd multiples
    of `pi/2`, where `wrap_pi` is smooth, and a `q` off by one there
    flips the parity sign too, so the fold silently compensates. Reducing
    modulo the period puts the tie on the branch cut; reducing modulo
    *half* the period and folding puts it somewhere harmless. Any
    "reduce modulo the full period instead" idea for a function with a
    branch cut needs this checked first.
  - Also needs a guard `wrap_pi` did not: halving a denormal (or anything
    in the smallest normal binade) drops its last bit, so `x*0.5` is not
    exact there.

## `gelu`'s remaining max is `erfc`'s, and `erfcx_pos`'s is its *evaluation*

Two closures, measured rather than argued, so neither gets re-attacked from
the wrong end.

**`gelu` (max 9, avg 0.209) is a pass-through.** Substituting a correctly
rounded `erfc` for the real one, while keeping every other f32 operation
`gelu` performs (the `RSQRT2_HI`/`RSQRT2_LO` double-`f32` argument, the
`dz` first-order correction, the `x*0.5*` combine), leaves **max 2.86 ulp**.
`gelu = x*0.5*erfc(z)` passes `erfc`'s *relative* error straight through, and
a relative error can cost up to twice as many ulp on the other side of a
binade -- so `erfc`'s 7 plus `gelu`'s own ~3 is the whole 9. Its avg 0.209
against `erfc`'s 0.199 is the same statement. Nothing inside `gelu`'s own
domain can move this; it moves when `erfc` does.

**`erfcx_pos`'s residual is not the fit and not `v`.** Simulating its exact
f32 instruction sequence in numpy over a dense sweep of `[1e-4, 12]`:

| variant | max ulp | avg |
|---|---|---|
| shipped (plain `v`, Estrin) | 5.08 | 0.730 |
| plain `v`, Horner | 4.36 | 0.660 |
| exact-residual `v`, Estrin | 4.20 | 0.672 |
| exact-residual `v`, Horner | 4.32 | 0.597 |
| **exact `v` (f64), Estrin** | **3.75** | 0.639 |
| fit alone (exact `v`, f64 arithmetic) | 0.72 | 0.290 |

This independently reproduces the rejection already recorded in
`erfcx_pos`'s own comment -- compensating `2+xa`'s rounding is worth ~0.9
ulp on max for 4 more ops -- and adds the part that closes the question:
even with a *perfect* `v` the sequence still measures 3.75, against a fit of
0.72. Roughly three quarters of the remaining budget is the f32 evaluation
of the polynomial itself, so no cheaper or more accurate `v`, and no refit
at this degree, is the lever. Note also that Horner beats Estrin on both
axes with the plain `v` and *loses* on max once `v` is corrected -- the two
choices are not independent, and Estrin is there for depth anyway (the
division is already on the critical path ahead of it).
## `cbrt_throughput` was dominated on every axis; deleted

Not a tuning failure -- a pareto check nobody had run. llvm-mca (this
worktree, `tools/mca_region.py`) against the accuracy sweep:

| function | avg ulp | max ulp | latency | throughput | domain |
|---|---|---|---|---|---|
| `cbrt_fast` | 57.4 | 554 | **28.05** | **0.839** | positive normal |
| `cbrt_unchecked` | **0.282** | **3** | 35.06 | 0.906 | normal, both signs |
| `cbrt_throughput` | 6.74 | 74 | 42.05 | 0.917 | positive normal |
| `cbrt_accurate_unchecked` | **0.000** | **1** | 63.00 | 2.067 | rescaled range |

`cbrt_throughput` loses to `cbrt_unchecked` on **all four**: 24x the avg
ulp, 25x the max, +20% latency, +1.2% throughput -- and it accepts a
*narrower* domain (no negatives, no denormals) to do it. A function named
for throughput that is slower than the accurate one it was meant to
undercut has no pareto point to sit on. Deleted, along with the dead
`cbrt_constant` (a `&[u32]`-parameterised helper with no callers and no
doc comment; `tune.rs` carries its own copies of that shape).

The reason it lost is structural, and it is the transferable part.
`cbrt_throughput` iterates the inverse cube root, `r <- fma(r*r, (r*r)*x,
r*K)`, which is **four** dependent FP levels per step (`r*r`, `*x`, the
fma, and the final `r*r*x`). `cbrt_fast` iterates the coupled pair
`s <- s + r*(x - s^3)` spelled `fma(s*s, s*-r, fma(r, x, s))`, which is
**two** -- `s*s`, `s*-r` and `fma(r,x,s)` all issue together. Two prior
sessions tried to rescue `cbrt_throughput` by tuning its constants (see
the two entries above); the constants were never the problem.

### `cbrt_fast`: two tuned variants measured, neither taken

Both dominate the shipped constants on accuracy and both give back part of
the only edge the tier has. Constants from a pareto-accepting coordinate
descent (accept only if one axis improves and the other does not worsen --
the max-first tuple scoring that wrecked `cbrt_throughput`'s average is
recorded above), scored on a 130k-point grid spanning every positive
normal octave and verified on 2.08M.

| variant | avg | max | latency | throughput |
|---|---|---|---|---|
| shipped | 57.4 | 554 | **28.047** | **0.839** |
| + second reciprocal-seed offset | 53.7 | 464 | **28.047** | 0.910 (+8.5%) |
| + exact `bits/3` as well | 52.8 | 383 | 29.047 (+3.6%) | 0.854 (+1.8%) |

- **A second reciprocal-seed offset is the bigger lever and the cheaper
  one on paper**: the two Newton steps want different reciprocals -- step 1
  wants `1/(3*s0^2)`, whose error carries the `s0` seed's own `2*u0`
  (~6%), step 2 wants `1/(3*x^(2/3))` -- so one offset cannot centre both.
  A second `from_bits(B2 - w)` is one integer add, entirely off the
  critical path, and latency is indeed unchanged to the cycle. It still
  costs 8.5% of throughput (+3 instrs on a 34-instr vector body), which
  puts `cbrt_fast` at 0.910 -- level with `cbrt_unchecked`'s 0.906, i.e.
  it trades away half of what makes the tier exist.
- **Exact `bits/3`** (`mulhi`+`shr`, then `t<<1` for the `-2/3` seed
  instead of a second multiply) removes a real slope error: `(bits>>16) *
  0x5556` scales by 0.33334351 instead of 1/3, a drift reaching 0.25% of
  the seed at the top of the range -- which is exactly where the worst
  point sits (`x = 4.25e37`). Worth another 17% of the max, for one cycle
  of latency.
- Not shipped because the accuracy floor is structural, not constant-bound:
  with seed errors `u0`, `v`, two steps of a fixed-`r` Newton leave
  `~u0*(u0+v)*v`, and `u0 ~ 3%` is what *any* single magic constant gives.
  Even a perfect `r` bottoms out near `u0^4 ~ 1e-6` (~18 ulp). Tuning
  moves 554 -> 383; nothing in this shape reaches an ulp, so paying
  latency or throughput for a 400-vs-550-ulp approximation buys nothing a
  caller of this tier can use. Adopt one only if `cbrt_fast` ever acquires
  a caller that cares about its error at all.

### The same constant, in `atand`, `atan2d` and `acosd`

`asind`'s `180.0 / std::f32::consts::PI` (one ulp low, see above) was
spelled identically in all four degree wrappers. Taking the two-word
`fma` to each of them, measured separately rather than assumed to
transfer:

| | before | after | verdict |
|---|---|---|---|
| `atand` | 0.3682 avg / 5 max | 0.0628 / 4 | ship |
| `atan2d` | 0.1749 / 5 | 0.1072 / 4 | ship |
| `acosd` | 0.0545 / 5 | 0.0546 / 5 | **reject** |

(`atand` and `acosd` exhaustive, all 2^32 patterns, both sides; `atan2d`
is two-argument so fuzz only, 3 repeats, avg stable to four digits and
max swinging 4-5 run to run as it always does there.)

**`acosd` is the interesting one, because it looks immune and is not.**
Its average does not move because two things cancel. Most of the domain
by bit-pattern measure is tiny `|x|`, where `acos(x)` is exactly
`fl(pi/2)` and the true `acosd` is exactly **90** -- a representable
value, so a bias of a third of an ulp cannot push the result off it.
Restricting to `|x| > 0.01` and feeding a correctly-rounded `acos` so
only the composite's own arithmetic is in view, the constant matters as
much as anywhere else: avg **0.657** with the shipped constant, 0.261
with the correctly-rounded one, 0.241 two-word. What flattens it back out
over the full domain is that `acos`'s *own* error is ~0.055 and dominates
whatever the multiply does.

So `acosd`'s low average is partly luck, and the two-word form is
*structurally* the better code -- but it measures identically on both
axes for +3 instructions, so it does not ship. Worth knowing that if
`acos` itself is ever tightened, `acosd` should be re-measured rather
than assumed still fine.

**`atan2d`'s denormal branch keeps the old constant on purpose.** Its
`y * (K/x)` uses `K` as a *divisor*, so the two-word pair would need a
second division. Swapping in just the correctly-rounded single constant
does measure better there -- over 3.6M pairs that reach the branch, avg
0.0258 -> 0.0110, exactly-rounded 97.53% -> 98.90% -- but it costs a
1-ulp regression at `atan2d(1e-30, 1e10)`, a pinned edgecheck case whose
current answer is the correctly-rounded one. That branch exists to
recover *binade* bits from an underflowed `y/x`, not the last ulp of a
denormal, so the pinned case wins.

Cost of the two that shipped (Block RThroughput unchanged for both, so
this is +3 instructions each and no new port bottleneck): `atand`
55 -> 58 instrs, 57 -> 60 uOps, throughput 2487 -> 2651 cycles, latency
429304 -> 454904; `atan2d` 72 -> 75 instrs, 74 -> 78 uOps, throughput
3405 -> 3496 cycles, latency 499809 -> 493409.

Still on the single-word multiply and *not* investigated here: `atanpi`
(0.2782 avg), `atan2pi` (0.1382), `acospi` (0.0438). `fl(1/pi)` is the
correctly-rounded constant, so there is no wrong-constant bug in that
family -- but it still carries a -0.338 ulp relative bias that the same
two-word `fma` would remove.

### The same treatment for `1/pi`: `atanpi` and `atan2pi` ship, `acospi` was a false lead

Following up the note directly above. Two of the three named there took
the two-word `fma(y, HI, y*LO)`; the third turned out not to be a
single-word multiply at all.

| | before | after | verdict |
|---|---|---|---|
| `atanpi` | 0.2783 avg / 4 max | **0.0775 / 4** | ship |
| `atan2pi` | 0.1383 / 4 | **0.1137 / 4** | ship |
| `acospi` | -- | -- | **not applicable** |

**`acospi` never had a `1/pi` multiply to fix.** The note above listed it
from its avg alone; it is `sqrt(1-a) * acospi_poly(a)`, a *dedicated*
degree-6 fit of `acos(a)/(pi*sqrt(1-a))` with the `1/pi` already inside
the coefficients and a trailing term pinned to exactly 0.5. There is no
second rounding to remove. The `acosd`-style "bulk mass sits where the
true answer is representable" argument was lined up for it and never
needed.

**`asinpi` was written off in the same sentence and that was wrong** --
see the entry below. Reasoning from `acospi`'s shape to `asinpi`'s
skipped the fact that `asinpi` has a *second* branch that `acospi` has
not got.

**The split lands the right way round for free here.** `1/pi` needs no
hand-picked low neighbour the way `180/pi` did: `fl(1/pi)` already sits
**0.43 ulp below** the real value, so `std::f32::consts::FRAC_1_PI`
doubles as the HI word and the only new constant is the positive tail
`FRAC_1_PI_LO = 1.28412765e-8`. Positive LO is what keeps `-0.0` signed
through the `fma`, the same trap `RAD_TO_DEG_LO` documents; `atanpi(-0)`
is pinned in edgecheck and passes.

That 0.43 ulp on the constant is a 0.34-0.68 ulp bias on the *result*
depending on where in its binade the result lands, which is the whole of
the gap between `atan` (0.0675) and old `atanpi` (0.2783). After the fix
`atanpi` is 0.0775 against `atan`'s 0.0675 -- essentially just `atan`'s
own error plus binade shift.

**`atan2pi` gains less than `atanpi` because it is already near its
floor.** `atan2`'s result spans `[-pi, pi]`, and dividing by pi shifts
`~pi -> ~1` and `~pi/2 -> ~0.5`, each a binade step that *doubles* the
error measured in ulp. So `atan2pi`'s floor is about 2x `atan2`'s 0.0681,
i.e. ~0.136, and it now measures 0.1137. Nothing further to take here
without changing `atan2` itself.

**The rescaled-coefficient fold stays rejected, and this does not
resurrect it.** Folding `1/pi` into `atan_poly`'s coefficients measures
0.2432 avg / 4 max at 1.612 cyc/elem -- *worse on both axes* than the
two-word form's 0.0775 at 1.657, and barely better than the plain
composite it replaced while costing more than it. `atan_poly` is a Pade
rational whose numerator and denominator share the same unscaled
trailing `+1.0`, so one broadcast normally serves both; scaling only the
numerator's copy breaks the sharing and forces a second. This is the
opposite verdict from `acospi_poly`/`asinpi_poly`, which are plain Horner
polys with a single trailing constant and fold for free -- the shared
constant is the whole difference.

Cost, from `tools/mca_region.py` (Block RThroughput unchanged in both
throughput regions, so this is +3 instructions each and no new port
bottleneck): `atanpi` 55 -> 58 instrs, 57 -> 60 uOps, throughput
2487 -> 2651 cycles; `atan2pi` 65 -> 68 instrs, 67 -> 70 uOps,
throughput 3009 -> 3137. Latency is +4 cycles per unit on both
(`atanpi` 67.079 -> 71.079, `atan2pi` 71.188 -> 75.188) -- exactly one
dependent `fma`, so there is no branch artifact to arbitrate here.
`atanpi`'s instruction and cycle counts come out identical to `atand`'s
in the entry above, which is the expected cross-check: same `atan`, same
two-word tail.

### The half-turn family was *not* closed: `asinpi`'s small branch, 0.2375 -> 0.0159 avg

The entry above closed the family on `acospi`'s shape and said "`asinpi`
is the same". `asinpi` is not the same, and the difference is visible
without measuring anything: `acospi` is a single expression, `asinpi` is
**two branches**, and only the big one is `sqrt(1-a)*poly`. Below
`|x| < 0.27` it takes `asinpi_small`, an odd poly `x * P(x^2)` whose
leading coefficient is a lone `FRAC_1_PI` -- exactly the single-word
`1/pi` multiply the family was being swept for, just spelled as a
polynomial coefficient instead of a `*` at the end.

The tell was in the table the whole time and is the reusable part:
**`asinpi` 0.2365 avg against its own siblings `acospi` 0.0437 and
`atanpi` 0.0775, and against `asin` 0.0187.** A function 5x worse than
the sibling it shares a construction with, and 13x worse than the
function it is a rescale of, is a defect in what differs, not a floor.
`asin_small` has the identical Horner shape and does *not* have this
problem for one reason: its leading coefficient is exactly `1.0`.

Fix is the crate's now-standard peel: give `FRAC_1_PI_LO` the trailing
Horner slot the leading coefficient used to hold (so the tail poly is
still three fmas), scale the tail by `x` on its own, and finish with
`fma(x, FRAC_1_PI, x*t)` -- one rounding, at the result's own magnitude,
with `x*FRAC_1_PI` exact inside the fma. `x*t`'s own rounding is ~0.004
ulp of the result at the branch edge and falls away from there. No
refit: the coefficients' mathematical values are untouched, only the
leading one's representation.

Exhaustive over all 2^32 patterns: **avg 0.2375 -> 0.0159**, max 5 both
sides (the max is the *big* branch's, at `|x|` just above the crossover,
and this does not touch it). The naive `asin(x)/PI` control row on the
same run is 0.2430/7, so the folded form is now decisively better on
both axes rather than marginally.

Cost, all three counters agreeing so mca's cycle column needs no
arbitration: 56 -> 61 instrs, 62 -> 68 uOps, BlockRT 13 -> 14,
throughput 0.901 -> 0.981 cyc/elem (+8.9%), latency 56.735 -> 60.845
(+7.2%). One extra `fma` plus its broadcast, the same price `atanpi`,
`atand`, `asind` and `erf` each paid this session for the same class of
fix. Dropping `c3` to pay for it does not work: the branch is a minimax
of degree 3 in `x^2` over `[0, 0.0729]` and the degree-2 truncation is
~5.4e-7 relative, ~9 ulp.

## `sin_checked`/`cos_checked`: three op-level removals, all bit-identical

Ten percent of both functions' throughput came off without a single
output bit changing. Verified the strongest way available: an exhaustive
2^32 bit-for-bit diff of `sin_checked`, `cos_checked`, `reduce_pi_checked`
and `reduce_pi_half_checked` against a frozen copy of the previous code,
zero differences on all four.

mca (`tools/mca_region.py`), all four counters moving together:

```
region                    instrs      uOps    BlockRT   cyc/elem
sin_checked_throughput  159->142  173->153    61->56   5.232->4.698  -10.2%
cos_checked_throughput  162->145  178->159    63->58   4.603->4.079  -11.4%
sin_checked_latency                                  117.02->108.02   -7.7%
cos_checked_latency                                  122.00->113.00   -7.4%
```

**1. The `POLY_SAFE_BOUND` clamp on `r` is subsumed by the `clamp(-1,1)`
on the result.** Both functions clamped the residual to `+-1000` *and*
clamped the poly's output to `[-1,1]`. The first is redundant given the
second: `sinf_poly_raw` keeps `r`'s sign for every `r`, and `|r| > 1000`
already puts `|sinf_poly_raw(r)| > 2.6e21`, so the pre-clip and the free
run land on the same `+-1.0` afterwards. Nothing overflows to `NaN`
either -- the only `inf - inf` candidate is `fma(b, y2, a)`, and `b < 0`
requires `y < 76` while `y2 = inf` requires `y > 1.8e19`, so `b` is
positive wherever `y2` is infinite. Checked exhaustively over all 2^32
residuals for both the raw and the `copysign` path, then again end to
end. This is the biggest single piece: -6 instrs and -8 latency cycles
per function on its own.
- `POLY_SAFE_BOUND` itself stays: `sind`/`cosd` have no `clamp(-1,1)`
  downstream to fall back on, so for them it is still what keeps a
  finite input from producing an infinite output.

**2. `quick_two_sum(p3, tier2)` -- previously rejected, re-screened,
flipped.** See the `reduce_pi: downgrade the last remaining full
two_sum` entry above for the numbers; the short version is that its
+13.7% `cos_checked` regression was an artifact of the code that removal
(1) has since deleted.

**3. `quick_two_sum(a, -b)` -> a Fast2Diff spelling.** `quick_two_sum`'s
error term is `b - (s - a)`; passing `-b` makes LLVM materialise the
negation with a real `vxorps` at two of `reduce_pi`'s three merges (the
third folds into `PI_LO`'s constant, which is already negative). Written
as `s = a - b; e = (a - s) - b` the negation disappears: `a - s` is the
exact negative of `s - a`, and float addition commutes, so it is the
same expression reassociated -- bit-identical, one op cheaper. -5 instrs
per function.

**What was screened and is *not* takeable, priced here so it is not
re-screened blind:**
- **Merging `e1` and `p2` into one `fma(qh, PI_LO, e1)`** (3 ops instead
  of `two_prod`'s 4, and one fewer merge in the chain: ~5 ops total).
  Dead on arithmetic: the merged value is `~x*2^-23.4`, so its own
  rounding is `~x*2^-47.4`, which at `|x| = 1e10` is 3.4e-5 against a
  residual that has to be right to ~6e-8. Three orders of magnitude too
  coarse. The four separate `two_prod` words are all load-bearing at
  `[1e8,1e10)`, and the reduction's real error floor is the final
  `s3 + err` rounding at half an ulp of `r` -- i.e. it is already sitting
  on the floor, and *any* added term at the 6e-8 scale doubles the
  budget.
- **Clamping `r` to `+-pi/2` instead, to drop the output clamp**:
  `sinf_poly_raw` exceeds 1.0 by one ulp at 21 distinct `r` in
  `[0.5, pi/2]` (worst at `r = 1.5705949`), so the output clamp is doing
  real work -- and it is *improving* accuracy there, since the true
  `sin` at those points rounds to exactly 1.0.

### Follow-on: two words of `1/pi` are enough, and the third-word fma pays for itself

`round_x_over_pi` carried a third correction word, `fma(x, RPI_TINY, s)`,
on top of `s = e0 + x * RPI_LO`. Both went: `RPI_TINY` deleted, and the
remaining correction folded into a single `fma(x, RPI_LO, e0)`. -5 instrs
/ -6 uOps / -2 BlockRT on `sin_checked` (4.698 -> **4.546**, -3.2%) and
the same instruction savings on `cos_checked` (4.079 -> **4.037**, -1.0%);
latency flat at 108.02 / 113.00.

**Why a third word of `1/pi` buys nothing here, unlike a third word of
`pi`.** `RPI_TINY` reaches `r` only through *which integer* `ql` rounds
to -- `reduce_pi` never sees it. A perturbation of size `d` can only move
`ql` for residuals within `d` of a half-integer, and at a half-integer
either choice is self-consistent: `q -> q+-1` flips the parity *and*
shifts `r` by `-+pi`, and `(-1)^(q+1) sin(r-pi) == (-1)^q sin(r)`
exactly. So `RPI_TINY` is inert until its own contribution `|x|*1.47e-16`
exceeds 0.5, i.e. `|x| > 3.4e15` -- inside the region the two-word `q`
has already lost to whole-integer error. This is the opposite of the
`PI_HI`/`PI_LO`/`PI_TINY` split on the *output* side, where every word
lands directly in `r` and dropping one is an unbounded relative error at
sin's zeros (see the `single replacement word` entry above).

Measured, 4M scored samples per band, against the previous code:
identical avg *and* max *and* worst-x on every band below 1e13 for both
functions, and better above -- `sin_checked [1e13,1e15)` avg
13384 -> **2975**, `cos_checked [1e13,1e15)` 8.90e6 -> **8.54e6**, both
`[1e15,+)` rows down slightly too. The fused `fma`'s extra rounding
saved more than the deleted word was contributing. A dedicated
near-zero probe (every f32 within 4 ulp of a zero of sin or cos out to
1e6, 11.5M points -- the set where a small absolute slop in `r` becomes
a large *relative* error) is bit-for-bit unchanged: 0.16675 avg /
1.516 max for sin and 0.16692 / 3.460 for cos, same worst x as before.

The exhaustive diff says the same thing more sharply: over all 2^32
inputs the **smallest** input whose result changes at all is 7.09e8 for
`sin_checked` and 3.67e9 for `cos_checked`. Everything below that is
bit-identical, and above it the differences are exactly the predicted
tie-flips -- 7.3% of the whole f32 line, overwhelmingly the `>= 1e15`
patterns that make up most of it.

Downstream, all unchanged or marginally better: `wrap_pi` 0.0244 avg
(its exhaustive baseline), `sinc` 0.0938 / 3 max, `sinc_unnormalized`
0.0716 -> 0.0715, `tan_checked` 3.2835e8 -> 3.2750e8 avg with the same
~2.3e9 ceiling. edgecheck clean, 26 tests pass.

### Two parity levers priced and both rejected

The `parity(qh)` / `parity(ql)` / compare / masked-xor block is 8 of
`sin_checked`'s ~70 vector instructions per iteration -- the largest
single remaining block. Two ways to shrink it, both measured, both no.

**Round `p0` to the nearest *even* integer instead of the nearest
integer.** `parity(qh)` is then identically 0, so the whole sign fixup
collapses to "is `ql` odd" -- `(ql*0.5).floor() != ql*0.5` straight into
a mask, no `fma`, no second `parity`. `qh` costs the same 3 ops either
way (`mul`, `vroundps $8`, `h+h` against `vpternlogd`, `add`,
`vroundps $11`), and the `+-1` that `qh` gives up is absorbed by `ql`,
which was already unbounded. `PI_TINY`'s term has to move from `qh` to
`qh + ql` (+1 op) or the reduction silently drops it whenever the even
round sends `qh` to 0.

- **mca refuses it**: -9 instrs, -11 uOps and -3 BlockRT on both
  functions, but throughput *up* -- `sin_checked` 4.546 -> 4.701 (+3.4%)
  and `cos_checked` 4.037 -> **5.020 (+24.4%)** -- and latency up 3
  cycles on both (108.02 -> 111.03, 113.00 -> 116.02). Instructions and
  ports both say win, cycles say lose, on both callers and on both axes.
  Not adopted on a counter disagreement this size.
- **And the accuracy is a genuine pareto, not a wash.** On the near-zero
  probe (11.5M f32 within 4 ulp of a zero of sin or cos) `sin_checked`'s
  max goes 1.516 -> 2.199 while `cos_checked`'s goes 3.460 -> **2.199**;
  band averages are unchanged for sin and ~5% *better* for cos at
  `(pi/4,10]` through `[1.3e7,1e8)`. So the pair's worst case improves
  while sin's alone gets worse. Worth knowing if `cos_checked` ever needs
  its average pulled down, but it is not the "free speed" this was.
- Without the `PI_TINY` fix the same variant is a straight regression:
  near-zero max **5.801** for both functions, worst at x = 505.8, where
  the even round makes `qh = 162, ql = -1` and `fma(qh, PI_TINY, c5)`
  then applies pi's third word to the wrong integer.

**`qh = round_ties_even` (the recorded 2 -> 6 cos regression),
re-screened on the current code because that measurement predates three
changes to its neighbourhood.** It does not go stale: `cos_checked`'s
near-zero max is 3.460 with `f32::round` and **5.801** with
`round_ties_even`, same worst x (252.9) as the even-round failure above,
i.e. the same mechanism -- cos's `-0.5` folded into `lo` turns `qh`'s
exact-tie cases into a wrong integer. `sin_checked` is completely
unaffected, as recorded.
- Methodology: **4M random samples per band see none of this.** Every
  band row for `round_ties_even` matches the shipped code to five
  decimal places (one is marginally *better*), max and worst-x
  included. The near-zero probe finds it in under a second. When a
  suspected loss is a *relative* error at a near-zero result, sample the
  near-zeros directly -- uniform-in-band fuzzing is the wrong instrument
  no matter how many samples it gets.

### `wrap_pi`'s extreme tail: the contract violation, measured and closed

Blessing `worst_corpus` after the `RPI_TINY` change turned up a
pre-existing contract violation, unrelated to that change but visible in
the same rows: `wrap_pi(f32::MAX)` returns 1.2377339e24 on master and
-6.612423e23 after. **Both are outside `(-pi, pi]`**, which is
`wrap_pi`'s whole documented contract. It rides `reduce_pi_checked`'s
deliberately unclamped `r`, so once the two-word `q` gives up (`|x|` past
~8.85e14) `r` is a large number and `wrap_pi` passes it straight through
-- `sin_checked`/`cos_checked` are protected by their own
`clamp(-1, 1)`, `wrap_pi` has no equivalent. **Fixed with the same shape
of clamp; see lib.rs.** What the exhaustive sweep of all 2^32 bit
patterns added to the diagnosis:

- **Scale.** 1,273,675,032 of the 4.28e9 finite inputs -- 30% -- were out
  of range, not a tail. Every input from `|x| ~ 2.83e22` up violated; the
  fraction of a binade that violates is already 1.4% at `|x| ~ 8.4e14`
  and passes 50% by `|x| ~ 5.4e16`.
- **The cliff is at 6.544881e14 (`0x5814d039`), not the ~8.85e14 that
  `sin_checked`'s comment documents.** Same 25% discrepancy, same
  direction, as the `sin_checked` cliff entry: `2^48*pi` is where the
  two-word `q` is *guaranteed* gone, not where the first input loses it.
- **A second, unrelated violation the same clamp fixes**: 16 inputs
  across the whole line returned exactly `-f32::consts::PI`, which is
  `pi + 8.7e-8` in magnitude and so outside `(-pi, pi]` at the open end.
  These are not reduction failures -- they are `r - PI` for a small
  positive `r`, and no amount of accuracy in the fold repairs them,
  because the nearest f32 to `-pi` *is* out of range. That makes the
  legal result set symmetric and bounded by `WRAP_PI_MAX`, the largest
  f32 below `pi` (`0x40490fda`), which is what the clamp uses, and which
  is *also* the correctly-rounded answer at `x = f32::consts::PI` (the
  old `-PI` there was 0.73 ulp off). Any "wrap to a half-open interval"
  contract has this problem; state the bound as a representable constant.
- **Cost, and how it was paid for.** `clamp` is exactly `vmaxps` +
  `vminps` (verified in the asm; both operand orders propagate a NaN
  residual, so `wrap_pi(nan)`/`wrap_pi(+-inf)` still come out `nan` with
  no extra select). Naively that is throughput 3.756 -> 4.031 cyc/elem
  and latency 91.111 -> 99.111 (+8.8%, exactly the 2x4-cycle chain, not
  the known branch artifact -- both regions have zero branches).
  Rewriting the half-turn fold from an `r > 0.0` select to
  `r - PI.copysign(r)` gives back 3 of those 8 cycles: LLVM lowers the
  copysign to a single `vpternlogd`, replacing a compare plus a blend,
  which is both the same instruction count and a shorter dependency
  chain. Shipped: throughput **3.854 (+2.6%)**, latency **96.017
  (+5.4%)**, instrs 118->124, uOps 128->134, BlockRT 46->48.
  Bit-identical to the select over all 2^32 inputs -- the only input that
  could tell them apart is an odd `q` with `r` exactly `+0.0`, which the
  same sweep confirms the reduction never produces.
- **Accuracy: unchanged.** Exhaustive `accuracy thorough wrap_pi` gives
  0.0244 avg / max 1 at worst x 1.5707964 before and after.
- **Why it survived**: `edgecheck` did pin the range, but only out to
  `|x| = 1e9` -- three decades under the cliff -- and with a `1e-6` slop
  that also let the `-PI` boundary case through. The pins now run to
  `f32::MAX`, include both exhaustively-found worst inputs, and compare
  against `WRAP_PI_MAX` exactly.

### The whole `< 1e13` accuracy claim, verified exhaustively rather than sampled

`RPI_TINY`'s removal is the only change this session that is not
bit-identical, so every band below 1e13 was swept **exhaustively** --
every f32 bit pattern in the band, both signs -- with the pre-change code
and the shipped code scored on the same inputs. Avg, max and worst-x come
out **identical to every printed digit on all five bands**, for both
functions:

```
band            n            sin avg    sin max   cos avg    cos max
[1e3,1e5)       110,272,512   0.174304    2.1139   0.195165    3.4603
[1e5,1.3e7)     117,840,512   0.174261    2.1230   0.194887    2.1005
[1.3e7,1e8)      49,331,648   0.174225    6.2338   0.182969    6.2338
[1e8,1e10)      111,971,762   0.174241   47.6119   0.174281   67.4972
[1e10,1e13)     167,314,396   0.189633  154382.5   0.186941  154382.5
```

(worst-x identical too; the only movement anywhere is `[1e10,1e13)`'s
average, 0.189633 -> 0.189612 and 0.186941 -> 0.186918, i.e. slightly
*better*. Averages are in this harness's own ulp convention, so compare
them to each other, not to `accuracy.rs`'s columns.)

**The max column is the reusable finding.** `accuracy.rs`'s quick mode
draws ~1.3M of the 112M patterns in `[1e8,1e10)` and reported max ulp 2,
3 and 10 for `sin_checked` on three different runs of the *same* binary.
The true value is 47. Three of the five bands have a true max well above
anything sampling reports -- 6 where it says 2, 47 where it says 3, 67
where it says 21. Per-band max from quick mode is a lower bound on a
long tail, not an estimate of it: for a band this narrow the exhaustive
sweep is only ~30 seconds, so just run it rather than repeating the
fuzz. (The averages, by contrast, are stable to 4 digits.)

## `sin`: a coarser magic grid plus a second, fine round -- 4x the domain

`sin` is now documented and measured over `|x| < 2^24*pi = 5.2707178e7`,
four times the old `2^22*pi = 1.3176794e7`, with **max ulp still 2** and
every other function in the crate byte-identical. It cost one add.

**The old ceiling was the magic-round grid, not `1/pi`'s precision.**
`ROUND_MAGIC = 1.5*2^23` puts `v` on the integer grid only while
`|v| < 2^22` -- past there the sum leaves the `[2^23, 2^24)` binade and
the grid coarsens -- and `2^22*pi` is exactly where sin's documented
domain stopped. That is a suspicious coincidence and it was the real
cause: feeding the *unchanged* `pi_reduce_and_poly!` chain an oracle `q`
computed in f64 gives avg 0.318 / max 2 across `[1.318e7, 5.27e7)`, where
shipped `sin` was ~1e9 avg. `RPI_HI`/`RPI_LO` resolve `x/pi` to ~6e-8
absolute out there, three orders better than needed.

**The fix.** A magic constant's reach and its grid trade off exactly:
`1.5*2^k` holds `|v| < 2^(k-1)` on a grid of `2^(k-23)`, so the ratio is
always `2^22`. `ROUND_MAGIC_4 = 1.5*2^25` buys `|v| < 2^24` at a grid of
4. The coarse index `n` that comes back is a multiple of 4 rather than the
nearest integer, which costs the free parity bit -- the same magic that
parks `q` on the grid is what parks parity in a fixed mantissa bit -- and
leaves `|fc| = |x/pi - n|` up to ~2.7 instead of ~0.7. Both are cheap to
repair: a *second* magic round, of `fc` (an O(1) value, so `ROUND_MAGIC`
has room to spare), gives `round(fc)` with its parity in the low bit
again, and `q = n + round(fc)` is `round(x/pi)` exactly. `n` is a multiple
of 4 so `q`'s parity is `round(fc)`'s.

The one added op is kept off the critical path: `q = (n - ROUND_MAGIC) +
(fc + ROUND_MAGIC)`, where `n - ROUND_MAGIC` is exact (both are multiples
of 4) and hangs off `n`, which is ready three ops before `fc` is. So the
`q` chain is one add *wider* and the same depth.

```
region              instrs      uOps    BlockRT   cyc/elem   latency
sin_throughput      59->64    62->68    18->19   1.776->1.778
sin_latency                                                  64.00->64.00
```

Everything else is byte-identical asm, checked region by region:
`cos_throughput`, `cos_latency`, `tan_throughput`, `tan_latency`,
`sin_fast_throughput`, `cos_fast_throughput`, `sind_throughput`,
`cospi_throughput`.

**`|q| <= 2^24` is the hard ceiling for this architecture** and the new
domain sits exactly on it: `q` has to be an *exactly representable f32
integer* for `pi_reduce_and_poly!`'s residual to mean anything, and 2^24
is the last integer f32 counts by ones. Confirmed independently with an
oracle two-word `q` and an 8-fma chain: clean max 3 at `2^25`, max 262145
at `2^26`, where a second wall appears -- after dropping Cody-Waite word
`k` the residual is `~q*eps_k` quantised at that word's lowest bit, so the
step stays exact only while `log2|q| <= 24 + (L_k - T_{k+1}) - 1`, and the
slack there is the length of pi's binary zero-run at that position (never
more than 4 in the first 100 bits). Getting past 2^24 needs a two-word
`q`, which *is* `sin_checked`.

**Rejected on the way (all measured):**

- **The same trick for `cos`: broken, not merely expensive.** cos's
  `q = n + copysign(0.5, fc)` is the nearest half-odd-integer to `x/pi`
  only while `|fc| <= 1`, and a magic round only ever bounds
  `|n - x*RPI_HI|`; the `RPI_LO` correction inside `fc` then adds up to
  0.17 at `2^22*pi` and 0.34 at `2^23*pi` on top of that. On the even grid
  (`1.5*2^24`) cos went from 0.2806 avg / max 2 to **0.4523 avg / max 205
  on `[1e5,1.3e7)` -- inside the old domain** -- and 0.3530 / 11419 over
  `|x| < 2^23*pi`. This is not a domain-edge effect: it is wrong wherever
  `x*RPI_LO` is not negligible. Doing it properly (round `fc - 0.5` on the
  fine grid, then `q = (n + k) + 0.5`) costs +3 float ops, and cos's half
  caps it at `2^23*pi` anyway, since a half-odd-integer needs a spare
  mantissa bit and so is exact only to 2^23. cos is left alone.
- **The even grid (`1.5*2^24`, 2x) for `sin`**: identical cost to the
  multiple-of-4 grid (64 instrs / 68 uOps / BlockRT 19), because the
  second round is needed either way. Strictly dominated by 4x.
- **Letting `tan` stay `sin(x) / cos(x)`.** sin and cos sharing one
  `frac_x_over_pi!` is worth ~5 ops to `tan`; different coarse grids kill
  the CSE. Measured: tan 105 -> 120 instrs, 107 -> 123 uOps, BlockRT
  31 -> 36, latency 78 -> 89 -- for *byte-identical output*, since tan's
  domain is its denominator's. Fixed by giving `tan` a private
  `sin_over_cos_domain` built on cos's own grid; tan's asm is then
  byte-identical to before. (With sin on the even grid and cos on it too
  the CSE survives -- tan 105 instrs, BlockRT 31 -- but that is the broken
  cos above.)
- **`cvtps2dq` instead of the magic round** (`round_ties_even` +
  `to_int_unchecked` fuse into one 1-uop instruction; parity from
  `(qi as u32) << 31`). Reaches the same 2^24 and is *cheaper on
  instruction count* -- 59 -> 60 instrs, uOps unchanged at 62 -- but
  1.776 -> **1.906 cyc/elem (+7.3%)** and latency 64 -> 75, and
  `to_int_unchecked` is UB out of range so a shippable version owes a
  guard on top. The float-only scheme above is cheaper on both columns
  that matter and has no UB.
- Dead by construction, do not re-derive: a *single* larger magic (its
  ulp-1 window is only 2^23 wide however it is placed, and the `1.5*2^k`
  family's reach/grid ratio is fixed at 2^22); a pre-scale of `x` (it
  destroys exactly the low bits the reduction needs); and coarse-then-fine
  on the *residual* rather than on the *index* -- a coarse residual of
  magnitude ~2^k has to fit one f32 whose ulp there is `2^(k-24)`, so
  `k <= 0`, i.e. any two-stage residual scheme needs a two-word
  intermediate and is `sin_checked` again.

Accuracy, `accuracy quick` (4M scored samples per band), old -> new:

```
row                            avg ulp                    max ulp
sin |x|<=1e6                0.0357 -> 0.0357              2 -> 2
sin (in-domain)             0.0422 -> 0.0457*             2 -> 2
sin [1e5,1.3e7)             0.2810 -> 0.2804              2 -> 2
sin [1.3e7,1e8)              1.034e9 -> 1.903e8 (-81.6%)
sin [1e8,1e10)               1.242e9 -> 1.151e9
sin [1e10,1e13)              1.828e9 -> 1.827e9    non-finite 12.99% -> 12.99%
sin [1e13,1e15)              2.134e9 -> 2.132e9    non-finite 95.89% -> 95.90%
sin [1e15,3.4e38)            2.139e9 -> 2.140e9    non-finite  100% -> 100%
```

*`sin (in-domain)` is not the same input set on the two sides -- it is
`sin_domain`, which moved from `2^22*pi` to `2^24*pi`. 0.0422 -> 0.0457
is the average over four times as much domain, not a regression; the
`|x|<=1e6` and `[1e5,1.3e7)` rows, which *are* the same set, are
unchanged. Nothing past the new limit got worse, including the non-finite
percentages, which are the same reduction failing the same way one binade
of `q` later.

`worst_corpus`: 5 of 9360 entries moved, all `sin`, none inside the old
domain. Four are at `|x|` between 1e30 and 1e38, where both old and new
are meaningless (`inf` flipping sign). The fifth is the informative one:
`sin(2.6e7)` was `+1.2973863e-1` and is now `-1.2775949e-1`, against a
true `-0.12775947990596045` -- the corpus had a wrong *sign* pinned at a
point that is now in-domain. `cos` and `tan` entries did not move at all.
edgecheck clean.

**The in-domain max is exhaustive, not sampled.** Every one of the
2,559,713,206 f32 bit patterns with `|x| < 2^24*pi` was scored: **max ulp
2**, avg 0.0457, worst x 1.3414527e0. So the headline number is a true
maximum over four times the old domain, not the lower bound quick mode
gives.

**Why a wider `sin` is still worth having next to a cheap `sin_checked`.**
The f64-reduction `sin_checked` that landed the same day covers ~7e15 at
2.495 cyc/elem, so the case for a *middle rung* is weak whenever the rung
costs anything -- 4x of domain for +7.3% (the `cvtps2dq` scheme) does not
beat 5e8x for +40%. This change is not a rung: it is the same tier, four
times wider, at 1.776 -> 1.778 cyc/elem and identical latency. The pareto
after it:

```
             cyc/elem   accurate to
sin_fast       1.151    ~1e6
sin            1.778    5.27e7
sin_checked    2.495    ~7e15
```

## `remainder_wide`: the `Df32` chain was doing f64's job, in f32

The crate's 3rd most expensive function, and its double-float residual
chain had never been attacked on its own terms (flagged as a live target
by backlog #157's own follow-up). Replacing the whole thing with a
reduction in f64 wins on **every** axis at once -- there is no tradeoff
here to weigh, so no second pareto point and no new function.

```
                  instrs   uOps  BlockRT   cyc/elem   latency
  remainder_wide  134->67 141->85 38->32   7.688 -> 2.266   171.17 -> 58.13
                                           (-70.5%)         (-66.0%)
  control: remainder_checked 1.357 / 45.17 unchanged, remainder_unchecked 0.646 unchanged
```

All four counters agree in direction, so the standing "mca's throughput
column is unreliable" caveat does not bite.

### Why it is *cheaper*, not merely wider

The old body carried `Df32::from_f32(xs) - Df32::from_mul(q0, ys)`, a
`.to_f32()`, a second divide for `adj`, a second `Df32` subtraction, and
an exact-power-of-two rescale of both operands near `f32::MAX`. Every one
of those exists to work around an f32 limit that f64 simply does not
have:

* **The coarse-`q` grid is gone.** `q` is an exact integer to `2^53`, so
  there is no multi-integer quantization gap for an `adj` pass to
  recover. That deletes a divide, a round, and a whole `Df32` subtract.
* **`x - q*y` needs no error-free transform at all.** It is exactly
  representable in an **f32**: `x` is a multiple of `ulp(x)` and `q*y` of
  `ulp(y)`, so the difference is a multiple of `min(ulp(x), ulp(y))` with
  magnitude `<= 1.5|y|` -- 24 significant bits. One `fma` forms `q*y` to
  full width internally and lands on it exactly. This is the same shape
  of argument that made `reduce_pi64` cheaper than the double-f32
  reduction it replaced: pick a working precision where the intermediate
  is *exact* and the bookkeeping evaporates.
* **The rescale is gone.** It existed because `Df32::from_mul` rounds
  `q0*y` to a single f32 before pairing it with its error term, so it
  could overflow near `f32::MAX` where a hardware `fma` would not. f64's
  exponent range makes that structurally impossible -- and with it the
  denormal precision loss the rescale could inflict on a tiny `x`.

What survives is `remainder_checked`'s own single `|r0| > |y|/2`
correction, for the one thing f64 does not make exact: `fl(x/y)` carries
a relative `2^-53`, up to half an integer at `|x/y| ~ 2^52`, so `q` can
still land one integer off across a half-integer boundary.

### Accuracy: the contract widens from `2^48` to `2^53`

Simulating both algorithms exactly (f32 ops via `Fraction` rounded once;
f64 natively) against an exact rational ties-away reference, 6000 samples
per band:

| band | old inexact | old worst ulp | new inexact | new worst ulp |
|---|---|---|---|---|
| `[2^20,2^48]` | 0 | 0 | 0 | 0 |
| `[2^48,2^50]` | 0 | 0 | 0 | 0 |
| `[2^50,2^52]` | 1532 | 9.85e19 | 0 | 0 |
| `[2^52,2^53]` | 3544 | 2.68e20 | 0 | 0 |
| `[2^53,2^56]` | -- | -- | (degrades) | -- |

Confirmed on the real shipped code, not just the simulation: widening
`accuracy.rs`'s `remainder_wide` domain from `|x/y| < 2e14` to `< 4e15`
reads **0.0000 avg / 0 max** for the f64 version and **2.55e6 avg /
3.42e9 max** for the `Df32` one, against sleef's independent IEEE754
remainder. The old bound was hiding the failure, not proving its absence
-- worth remembering when a harness domain is set to "comfortably under"
a hand-derived limit.

Bit-identical to `remainder_checked` on 29.3M in-domain samples
(`unchecked_parity.rs`), `worst_corpus` bit-identical with no re-bless,
and edgecheck / special_matrix2 / saturation_pins / denormal_audit clean.

### Tooling: `codegen_check` could not see an f64-vectorized region

`has_packed_arith` required a mnemonic containing `ps\t` -- single
precision only. A function that does its whole job in f64 emits `vdivpd`
/ `vfnmadd213pd` / `vrndscalepd` on `zmm` and has no `ps` arithmetic
anywhere, so it failed as "loop may not have vectorized at all" while
being perfectly vectorized (`vcvtps2pd %ymm -> %zmm`, 8 f32 lanes per
iteration). Fixed to accept `pd\t` as well. This gap would have hit any
future f64-internal function; `reduce_pi64`'s callers never tripped it
only because they keep an f32 polynomial.
## `asinh`: deleting the rescale arm is free code, and mca's own model says no

**Rejected 2026-08-01.** `asinh`'s `ax >= 2048` arm (`sq = fma(0.5, 1/ax,
ax)`, graveyard #54b) is *removable*, not just cheapenable: the function
already carries an `ln(ax) + LN_2` arm for the `ax^2`/`2*ax` overflow
case, and `asinh(ax) - ln(2*ax) = 1/(4*ax^2) - ...` is 0.0625 ulp at
`ax = 2048` and 0.0039 ulp at `8192`, falling as `1/ax^2`. So the
overflow arm can serve the whole tail and the rescale arm, its division,
its fma, the `sq`/`sm1` selects and the `d.is_finite()` test all go --
`d.is_finite()` becomes `ax < 8192`, which is ready at cycle 1 instead of
after the `sqrt`+divide.

The static ladder says take it, unanimously:

| | instrs | uOps | BlockRT | vdivps | FPDiv press | Port0 | Port5 |
|---|---|---|---|---|---|---|---|
| shipped | 156 | 163 | 42.0 | 6 | 42.00 | 51.37 | 35.59 |
| this | 145 | 151 | 40.0 | 4 | 32.00 | 48.52 | 29.94 |

Every expensive opcode is down too (`vfmadd213ps` 30 -> 28, `vaddps`
18 -> 16). Accuracy is unchanged: quick fuzz avg **0.1493**, max **3**,
identical to shipped to four places.

**And mca's cycle column says +15.1% throughput** (6618 -> 7616,
4.136 -> 4.760 cyc/elem) **and +78.5% latency** (69.41 -> 123.88).

Both halves of that were run down, and they have different answers.

- **The latency number is an artifact, and it is the *shipped* side that
  is wrong.** In the scalar latency harness LLVM lowers `if small` to a
  real branch, laying the two arms out sequentially. llvm-mca has no
  branch predictor and executes the region as straight-line code, so the
  value that survives into `d = ax + sm1` is whichever arm is written
  *last* -- here `.LBB263_1`, the rescale arm (`vdivss`, `vfmadd132ss`,
  `vaddss`, ~22 cyc from `ax`) rather than the direct arm (`vmulss`,
  `vaddss`, `vsqrtss`, `vaddss`, `vdivss`, ~41 cyc). Removing the branch
  forces mca to simulate the real path. Confirmed three ways: the asm
  layout above, a hand count of the honest chain (~90-100 cyc, i.e.
  nearer 124 than 69), and wall clock, where the two builds' latency is
  statistically identical (shipped 42.50/37.86/33.20/36.72 ns, this
  37.18/37.41/35.61/37.71 ns) as it must be, since the in-domain chain
  differs only by two deleted selects. **`asinh`'s published mca latency
  of 69.41 has never been its in-domain latency.** `acosh` (89.08) has
  the same two-arm shape and is likely understated the same way.
- **The throughput number is a model effect, not an artifact, and it does
  not transfer.** Every region in this crate is scheduler-queue-bound in
  mca, not port-bound -- `SCHEDQ - Scheduler full` is 87-99% of cycles on
  `ln`/`exp`/`acosh`/`atanh`/`tanh`/`erfc`/`cbrt`/`asinh` alike, and every
  region simulates 1.06-1.74x its own `Block RThroughput`. mca models the
  ICX scheduler as 60 entries; the real Ice Lake/Tiger Lake unified RS is
  160. Under a 60-entry window, the shipped form wins because `1.0/ax` is
  an *independent* long-latency op that issues at cycle ~5 and frees its
  slots while the `sqrt`->divide chain runs; the new form has nothing to
  overlap and clogs the window. Under a 160-entry window that pressure
  mostly disappears.

**So the tie-break was real hardware, and it says neither.** 20 alternating
rounds of two prebuilt `quickbench` binaries under the bench lock, with
`acosh` (byte-identical in both builds) as an in-binary drift control:

```
        asinh min   asinh p25   asinh med | acosh min  acosh med
shipped   1.042       1.490       1.615   |   0.910      1.422
this      0.979       1.440       1.562   |   0.867      1.426
```

Normalised on the control's min, 1.145 vs 1.129: **-1.4%, i.e. a wash.**
The control's own min moved 4.7% between the interleaved runs, which
bounds what this machine can resolve. The divider and port savings do not
convert because the region is latency-bound on `ax^2 -> +1 -> sqrt -> +1
-> divide -> ... -> ln_normal` on real hardware too, and that chain is
untouched by the change.

Not shipped: no axis improves. Accuracy identical, real speed a wash, and
the crate's published mca number would get 15% worse for it.

Two things worth carrying forward:

- **`SCHEDQ - Scheduler full` at ~99% is the tell for this class of
  disagreement.** When it is pegged and `Block RThroughput` moves the
  *opposite* way from the cycle count, mca is ranking how well the code
  fills a 60-entry window, not how much work it does. Check it before
  spending a wall-clock A/B: it is one `--all-stats` run.
- **A branch-shaped region can make mca's latency column measure the
  wrong arm entirely**, not merely mis-time the right one. The existing
  memory covers mca over-reporting latency for branchy code; this is the
  other direction, and it is worse, because the too-*good* number is the
  one that gets published. The tell is a function whose mca latency is
  well under a hand count of its own dependency chain.

## `clog`: 4717 max ulp was `cabs`, and it did not have to be composed on

**Shipped 2026-08-01.** `clog`'s real part measured **max 4717 ulp**
(`re=1.0001993, im=3.7058e-4`) in the standing fuzz -- the crate's worst
undocumented max. Its doc comment attributed that to
"an information-theoretic floor of composing on top of plain `f32`
`cabs` ... not something reachable without a compensated/double-float
magnitude ... out of scope for this composite". **The premise was wrong:
the floor is real, but only for code that composes on `cabs`, and the
near-1 branch does not have to.**

`ln|z| = 0.5*log1p(|z|^2 - 1)` reaches the same answer without ever
rounding at magnitude 1. `re^2 + im^2 - 1` is built exactly:

```rust
let p1 = re * re;  let e1 = fma(re, re, -p1);
let p2 = im * im;  let e2 = fma(im, im, -p2);
let (s, es) = two_sum(p1, p2);
0.5 * log1p_guarded((s - 1.0) + ((e1 + e2) + es))
```

`s + es + e1 + e2` *is* `re^2 + im^2`, with no error at all, and
`s - 1.0` is Sterbenz-exact wherever the result is small enough to care
(`s` in `[0.5, 2]`). The branch guard still uses `cabs`, which is
harmless -- it only has to be right to a factor of two.

- **Accuracy: max ulp re 4717 -> 3**, stable across 3 reruns (im
  unchanged at 3). Both fuzz numbers are sampling-limited, so the honest
  statement is the manifold one below, not the ratio.
- **mca: latency 170.63 -> 158.93 (-6.9%)**, throughput 8.166 -> 8.227
  (**+0.75%**), for +26 instructions and +32 uOps in the region
  (`Block RThroughput` 65 -> 71). The latency *improves* because the
  value path no longer waits on `cabs`'s `sqrt` at all -- only the branch
  select does -- and the +13 ops per element land in slack on an already
  long chain. Compare `logit`'s own precedent (max 1024 -> 3 for +11.8%
  throughput): this one is nearly free.
- All 8 standing gates, `approx_bounds` and 26 tests pass; `worst_corpus`
  needed no re-bless.

**Two intermediate results worth keeping.**

- **Half the fix is not most of the fix.** The first version compensated
  only `re` (`fma(re, re, -1.0) + im*im`) and took max 4717 -> **1466**,
  not to single digits -- the new worst case simply moved to
  `re=-9.04e-13, im=1.0002066`, where the surviving rounding is `im*im`'s
  at magnitude 1. Symmetric problem, symmetric fix required. Ordering the
  pair by magnitude and compensating only the larger would have failed
  the same way on the `|re| ~ |im| ~ 1/sqrt(2)` diagonal, where both
  squares round at magnitude 1/2 and then cancel.
- **The fuzz's max is a sampling artifact for this function, in both
  directions.** `clog`'s error lives on the `|z| = 1` manifold, which
  random `(re, im)` f32 pairs essentially never land on; the reported max
  is just how close 5M samples got. A targeted 400k-point sweep *along*
  that manifold (pick `re`, solve `im`, then jitter `im` by +-40 ulp)
  puts the old form at **1.68e7 ulp** and the new one at **2048** at the
  same points -- and the old figure is generous to the old code, since
  the probe gave it a correctly-rounded `hypot` where the shipped `cabs`
  carries up to 1 ulp. Neither is bounded, and neither can be: as
  `|z| -> 1` the true answer goes to zero and so does its ulp. What the
  fix actually buys is moving the breakdown from `|z^2-1| ~ 2^-24` down
  to `~2^-48`.

**Method note: `clog` still must not be mca-wired permanently.** The
numbers above came from temporarily adding `clogtmp_latency`/
`clogtmp_throughput` to `mca_target.rs` and reverting it. That is safe as
an A/B (identical wiring on both sides, and the region's +26 instructions
are exactly the +13 added ops x 2 vector iterations), and `ln`/`cabs`'s
own regions came back bit-identical -- but a full 295-region diff shows
adding it **does** corrupt `remainder_wide_latency` (4056 -> 1495
instructions) and `remainder_wide_throughput` (134 -> 67), exactly the
multi-exit-path marker corruption `mca_target.rs` documents. So the
existing decision not to wire it stands; check the *whole* region list,
not two controls, before believing otherwise.

## `powf`: the `Df32` chain was doing f64's job, again -- 3 -> 1 max ulp

Same lever as `remainder_wide` above, on the crate's 2nd most expensive
function. The whole `log2` / multiply-by-`y` / `exp2` chain moves from
double-f32 (`log2_df` + `Df32 * f32` + `exp2_checked_df`) to f64, coming
back to f32 once, at the end.

```
                   instrs   uOps  BlockRT   cyc/elem            latency
  powf            233->209 274->267 65->76  8.459 -> 6.996 (-17.3%)  125.17 -> 125.81
  powf_unchecked  170->140 200->179 60->70  6.509 -> 6.171 ( -5.2%)  117.28 -> 118.00
```

Accuracy, both harnesses, baselines re-measured on the spot rather than
read off the readme:

| | old | new |
|---|---|---|
| `powfsearch` worst over all sweeps | **3** | **1** |
| `powfsearch` worst per-band avg | 0.3027 | 0.0106 |
| `accuracy quick` powf avg | 0.0193 | 0.0007 |
| `accuracy quick` powf_pos / powf_unchecked avg | 0.0387 | 0.0014 |

So `powf` is now faithfully rounded, matching std's own max of 1, and the
latency rows are a wash. Nothing is dominated, so no `powf_latency` tier.

### Why this needed care that `remainder_wide` did not

**f64 buys no lane throughput on this machine.** Measured directly with
llvm-mca on single instructions: `vfmadd213pd %zmm` has Block RThroughput
**1.0**, `vfmadd213ps %ymm` has **0.5**. Both process 8 lanes. So f64 ZMM
costs exactly 2x f32 YMM per lane, and an f64 rewrite only wins if it
*shortens the algorithm*, which is the entire question. `remainder_wide`
halved its instruction count and won 70%; `powf` cuts 10-18% and wins
5-17%. Do not expect the `remainder_wide` multiple anywhere the f32
version was not doing obvious bookkeeping busywork.

**`vdivpd %zmm` is 16.0 Block RThroughput**, so the atanh form's
`t = (m-1)/(m+1)` still cannot be a division even in f64 -- a seed plus
Newton refinement is ~6 ops against 16. Priced, not assumed.

**Constant pressure is a real cost in f64.** The first version had ~23
distinct f64 constants and llvm-mca showed **24 `vbroadcastsd`** inside
the unrolled loop body: LLVM ran out of ZMM registers to keep them live
and rematerialized. Shortening both polynomials to the accuracy actually
required (atanh tail 5 -> 4 coefficients at `2^-37.6` against `2^-31`
needed; exp2 tail 7 -> 6 at `2^-28.9` against `2^-25`) took `powf` from
211 to 207 instructions on its own.

**The seed/Newton split is a latency knob, and the obvious setting is the
wrong one.** A degree-3 seed (`2^-13`) plus *two* Newton steps and a
degree-5 seed (`2^-20`) plus *one* are within an instruction of each
other, but one dependent level apart:

| | powf cyc/elem | powf latency | powf_unchecked cyc/elem | powf_unchecked latency |
|---|---|---|---|---|
| deg-3 seed, 2 Newton | 7.251 | 131.45 (+5.0%) | **5.734** | 124.35 (+6.0%) |
| deg-5 seed, 1 Newton | **6.996** | **125.81** (+0.5%) | 6.171 | **118.00** (+0.6%) |

The two-Newton version is a genuine tradeoff (throughput win, ~5% latency
loss) and would have needed a pareto decision. One Newton removes the
tradeoff outright, and *also* wins `powf`'s own throughput row. Worth
remembering as a general shape: when an f64 port lands "faster throughput,
slower latency", check whether a shorter refinement chain buys the level
back before reaching for a second public function.

### What did not move

`compound_accurate` (255 instrs, 9.239 cyc/elem, now the crate's most
expensive function) still routes through `log2p1_df` -> `log2_df` and
`exp2_checked_df`, so both double-f32 helpers stay live and nothing is
dead. It is the obvious next target for the same treatment; the `log2p1`
shape (`log2_df` of `1+x` with the cancellation handled) is the part that
needs designing, not the `exp2` half.

### Methodology

The mca signals disagreed: instructions **down**, uOps **down**, Block
RThroughput **up** 13-17%, simulated cycles **down**. That combination is
not in the escalation ladder. Arbitrated with the documented prebuilt-
binary wall-clock A/B (12 alternating rounds under the bench lock, min per
variant), which agreed with the cycle column and against RThroughput:
throughput old 2.354 / new 2.263 ns, latency old 42.26 / new 43.79 ns for
the two-Newton build. The reason RThroughput is not binding is visible in
the numbers -- the loop measures 5.7-6.5 cyc/elem against a resource floor
of 3.75-4.25, so it is dependency- and front-end-bound, not port-bound.
Medians across the 12 rounds were useless (2.4-4.5 ns spread on the *same*
binary); only the min is a usable estimator on this machine.
## `dawson`: 15 max ulp was the fit, and 5/5 was already optimal -- the lever was a degree

2026-08-01. `dawson` was the crate's worst *ordinary* function on a fresh
`accuracy quick` sweep -- max 15, ahead of `srgb_to_linear`'s 14 and
`linear_to_srgb`/`gelu`/`asind`'s 8 -- once the reduction-limited trig rows
(`sin_checked`/`tan_checked` and friends, closed separately) and the
documented y-amplified ones (`compound`, `logaddexp`, `powf`) are set aside.

**Shipped: avg 0.1527 -> 0.0561, max 15 -> 5**, for llvm-mca throughput
1.938 -> 1.992 (+2.8%) and latency 65.89 -> 68.97 (+4.7%), Block
RThroughput unchanged at 30.00.

### Attribution first: the fit, not the chain, and not `u = fl(x*x)`

A four-row decomposition per octave (shipped f32 chain / the same rational
evaluated in f64 from an exact `u` / the same from `u = fl(x*x)` / an
oracle `x * fl(dawsn(x)/x)`) said the answer immediately, and it is worth
recording because two of the three plausible culprits were flat:

| band | shipped | fit only | + `u` rounding | oracle `x*fl(R)` |
|---|---|---|---|---|
| `[0.25,0.5)` | 6.73 avg / 13 max | 6.74 / 11 | 6.74 / 11 | 0.19 / 1 |
| `[1,2)` | 5.95 / 13 | 5.94 / 11 | 5.94 / 11 | 0.26 / 1 |
| `[2,3)` | 5.43 / 16 | 5.40 / 12 | 5.41 / 13 | 0.28 / 1 |

So `u = fl(x*x)` contributes ~0.01 avg (the amplification `u*R'/R` is
under 1 over the whole branch -- `dawsn(x)/x` decays like `1/(2u)`, so the
argument rounding cannot be amplified here the way `gelu`'s and
`norm_cdf`'s were), the f32 evaluation chain contributes ~0.02 avg and
~2 max, and everything else is the rational's own fit.

**The bit-pattern-uniform harness hides all of it.** The reported avg was
0.1527 while the *central branch's* real avg was 4-7 ulp across every
octave in `[0.06, 4]`: ~45% of f32 patterns are `|x|` small enough that
`P/Q` is exactly 1.0 (the `pc[0] = 1.0` pin), and ~49% are `|x| > 4` and
take the tail. Only ~6% of samples land where the central fit is even
visible. Same shape as `acosd`'s "the bulk mass sits where the answer is
representable" -- check the per-band profile before believing a low avg.

### `[5/5]` really was exhausted; `[6/6]` is 5.6x better for +6 instructions

Re-derived the minimax rational independently (HiGHS, linearised
`|P - g*Q| <= eps*g*Q_prev`, iterated to a fixed point) after the first
attempt in a plain `u`-monomial basis returned nonsense (`eps` driven to
0 with a real error of 4.4e-5 -- the LP satisfying its own constraints
numerically while the rational it names is garbage). **Rescaling the
variable to `t = u/16` fixes it outright**; the monomial basis over
`u in [0,16]` spans `16^5 ~ 1e6` and HiGHS cannot see through that. Worth
copying: if a rational LP here returns an `eps` far below the achieved
error, suspect the basis before the solver.

With that fixed, the shipped `[5/5]` measures 13.12 ulp-equivalent against
a freshly-computed `[5/5]` minimax optimum of **13.07** -- confirming the
existing "the frontier is flat, do not re-run this sweep" entry, and
confirming it was a statement about *that degree*:

| | continuous minimax | after f32 quantization |
|---|---|---|
| `[5/5]` (shipped) | 7.79e-7 (13.07 ulp-eq) | 13.04 |
| `[6/5]` | 1.67e-7 (2.80) | 2.82 |
| `[5/6]` | 2.53e-7 (4.24) | 4.49 |
| `[6/6]` | 1.25e-7 (2.09) | **2.36** (after a quantization descent) |
| `[7/6]` | 6.45e-8 (1.08) | 3.07 -- quantization dominates |

`[7/6]` and `[7/7]` are the point where independently rounding the
continuous optimum to f32 gives back more than the extra degree buys, so
the search stops at `[6/6]`.

**`[6/6]` beats `[6/5]` on both perf axes, and the reason is codegen, not
op count.** `[6/5]` measures throughput 1.939 / latency 71.83; `[6/6]`,
which is *three more instructions*, measures **1.906 / 69.81**. The
latency harness's asm explains it: at `[5/5]` and `[6/6]` LLVM packs the
numerator and denominator into the two halves of one SIMD register and
evaluates them together (`vfmadd231ps` in a scalar chain), and at `[6/5]`
the two sides no longer have the same shape, so it falls back to scalar
`vfmadd213ss` plus 64 `jmp`/`jae` -- a branch-shaped region llvm-mca then
misprices on top. **Keep a rational's two sides the same degree unless
there is a measured reason not to.**

The extra coefficient on each side rides into the existing `u^2` group
(`ln_normal`'s `c[8]` fold), so it is one `fma` per side and no new
multiply.

### The tail's 22 ulp sat exactly on the seam, and L1-under-a-max-cap took it

With the central branch down, the reported max moved to `|x| ~ 4.0`, i.e.
the *tail* branch -- exactly the condition the earlier tail entry names
("revisit only if the central branch ever drops below ~10, and then fit
L1-under-a-max-cap rather than pure minimax"). It is worth being precise
about the shape of the defect: the shipped degree-4 least-squares fit's
relative error over `v in [0,1/16]` is a monotone climb to a **spike at
the endpoint** -- 0.83 ulp-equiv at `x=16`, 2.4 at `x=8`, 5.2 at `x=4.5`,
and **21.9 at `x=4`**. That single endpoint was `dawson`'s whole reported
max.

The right objective needs both halves and neither alone works:

- **The weight.** A bit-pattern-uniform caller reaches every octave of
  `|x|` equally often, so a log-uniform `x` grid *is* the correct L1
  weight, and it puts nearly all the mass at `v ~ 0`. This is why the
  earlier pure-minimax refit regressed: it spread error uniformly in `v`.
- **The cap.** Without it the L1 optimum parks its error at `v = 1/16`,
  which is what least squares already did.

| | weighted L1 | max (ulp-eq) | at the seam |
|---|---|---|---|
| shipped degree 4 | 0.2967 | 21.86 | 21.86 |
| degree 4, cap 5 | 0.6196 | 5.52 | 4.99 |
| degree 4, cap 8 | 0.3992 | 8.00 | 8.00 |
| **degree 5, cap 1.35** | **0.0617** | **2.37** | **0.46** |

Degree 4 cannot get under ~5 max at any cap and pays 2x the L1 to get
there; degree 5 wins on **both** axes at once, by 4.8x and 9.2x. The LP's
own optimum puts the `v^4` coefficient at zero, so the shipped term set is
`1, v, v^2, v^3, v^5` and the top group carries `v^5` alone -- one extra
multiply, not an extra `fma`.

Non-negativity was imposed on the tail coefficients rather than taken:
the unconstrained optimum is marginally better (L1 0.0616 / max 1.42) but
wants a negative `v^4`, and all-positive is what keeps `v = +inf` -- the
discarded arm for small `|x|` -- combining to a consistently-signed `+inf`
instead of `erfinv`'s opposite-signed-infinity `NaN`. The cost of keeping
that invariant is 1.42 -> 2.37 ulp-equiv on a term that is no longer
binding.

### Not taken here

- **Pinning the tail's `c[1]` to exactly `0.5`** (its exact asymptotic
  value, and exactly representable). Tested as a first-class candidate on
  the `pc[0] = 1.0` precedent, and it does *not* pay: degree 5 with `c[1]`
  free reaches L1 0.0568 / max 2.01, pinned reaches 0.0949 / 2.13. Unlike
  `pc[0]`, `c[1] = 0.5` does not make any region *exact* -- `fma(0.5, v,
  1.0)` is already exactly `1.0` for `v < 1.2e-7` whatever `c[1]` is, so
  the pin buys nothing and costs a real degree of freedom.
- The `[7/6]`/`[7/7]` rationals, on the quantization result above.
## `srgb_to_linear`: split the integer power out, and the round trip stops mattering

**Shipped 2026-08-01, and it reopens a closed entry.** The
"`srgb_to_linear`'s 13 max ulp is the log2/exp2 round trip, not its own
algebra" entry above measured the attribution correctly and then drew the
wrong conclusion from it -- "closing it properly needs a genuinely more
accurate `log2`, and that is the known dead end". It does not. Exhaustive
over the whole `[0,1]` domain:

| | avg ulp | max ulp | mca latency | mca throughput |
|---|---|---|---|---|
| was | 0.1046 | 14 | 110.02 | 3.831 |
| now | **0.0823** | **7** | **104.30 (-5.2%)** | 4.289 (+12.0%) |

`linear_to_srgb` is untouched and re-measures bit-identically (0.1209 /
max 8 -- note 8, not the 7 the quick fuzz reports).

**Two changes, and neither works without the other.**

1. **`b^2.4` as `b*b * b^0.4`.** The entry above derives the round trip's
   contribution as `exp2` amplifying an *absolute* argument error:
   `relerr = ln2 * |y| * relerr(log2)` for `y = 2.4*log2(b)`. Read that
   with `y = log2(result)`, and it says the amplification depends only on
   how far the *result* is from 1 -- which is why re-associating as
   `(b^2)^1.2` or `sqrt(b)^4.8` buys nothing. But peeling off an *integer*
   power does, because that factor is then computed by a multiply instead
   of through the log: `b*b` is exact but for one rounding, and the log
   and `exp2` are left carrying `0.4*log2(b)` instead of `2.4*log2(b)`,
   a 6x smaller absolute error to amplify. Same lever as `rootn`'s
   exponent split (45 -> 1 max ulp), and it needs no more accurate `log2`
   -- which is the whole point, since `log_2_normal` is a core-locked
   kernel behind 9+ functions.
2. **`b = fma(c, INV_1055, OFF_1055)`** with both constants the nearest
   `f32` to the *exact* `1/1.055` and `0.055/1.055`, replacing
   `(c + 0.055) * (1.0 / 1.055)`. One rounding instead of two, one op
   cheaper, and it fixes a real constant bug: `1.0f32 / 1.055f32` divides
   by an already-rounded `1.055` and lands **0.53 ulp** above the true
   `1/1.055`.

**The trap, and it is the interesting part: change 2 alone is a
regression.** Keeping `*2.4` and only fixing `b` measures avg
**0.1046 -> 0.1339** and max **14 -> 15**, exhaustively worse on both,
despite being strictly fewer roundings on strictly better constants. The
two systematic errors were **cancelling**: `f32(2.4)` is 9.5e-8 high, so
`b^f32(2.4)` runs `2.4 * 9.5e-8 = 2.3e-7` *low* for `b < 1`, while the
0.53-ulp-high `1/1.055` ran the result `2.4 * 6.3e-8 = 1.5e-7` *high*.
Removing either one alone exposes the other. Change 1 happens to fix the
exponent constant too -- `2 + f32(0.4)` names 2.4 sixteen times more
closely than `f32(2.4)` does -- which is why the pair works.

Generalisable: **a systematic constant error that has been sitting in a
composite for a while may be load-bearing.** Before "correcting" a
constant to its exact value, check the other constants on the same path
for an opposing bias -- the crate already knows pinning a *fitted*
coefficient to its exact value costs (the `erfc` entry); this is the same
hazard for constants nobody ever fitted.

**Two things measured and not taken.**

- **`b = (c + 0.055) / 1.055`** (a real division, one correctly-rounded
  operation instead of a multiply by a rounded reciprocal): exhaustively
  *worse*, avg 0.1016 max 8 against the fma form's 0.0823/7, because
  `f32(1.055)` is itself 0.44 ulp low. The fma with two exact-value
  constants beats both, and costs a division less.
- **Compensating `b*b`** (`bb = b*b; be = fma(b,b,-bb); fma(be, p, bb*p)`):
  **zero movement**, avg 0.0922 max 7 to four places, identical to the
  plain `b*b*p` at that stage. `b^2`'s own rounding is not the binding
  term once the round trip is down; reverted, 2 ops cheaper.

**The cost is real and is on throughput only.** +2 instructions, +5 uOps,
`Block RThroughput` flat at 41, but port pressure genuinely up (Port0
44.04 -> 46.04, Port1 44.71 -> 47.03 -- exactly the two added multiplies)
and simulated cycles +12.0%. Latency *improves* 5.2%. Taken on the same
trade `logit` was (max 1024 -> 3 for +11.8% throughput): 14 ulp on a
`[0,1]` colour transfer function is outside this crate's stated budget in
a way 12% of throughput is not. If a caller ever needs the throughput
back, the pre-split form is a valid `srgb_to_linear_fast` -- it is
Pareto-optimal on exactly one axis, and the numbers to wire it up are the
table above.

## `compound_accurate`: the same f64 port, and this one is clean on every counter

Third and last application of the `Df32`-is-doing-f64's-job lever, on what
was (after `powf` landed) the crate's most expensive function. Unlike
`powf`, there is no counter that moves the wrong way.

```
                       instrs   uOps  BlockRT   cyc/elem            latency
  compound_accurate   255->159 302->211 84.5->76  9.239 -> 6.173 (-33.2%)  153.56 -> 120.77 (-21.4%)
  accuracy (quick):   avg 0.0194 -> 0.0007, max 5-6 -> 1
```

Block RThroughput *improves* here where `powf`'s got worse, and the reason
is the whole point: `powf`'s f32 version was already fairly lean, so the
f64 port only cut 10-18% of the instructions against f64's 2x per-lane
port cost. `log2p1_df` was not lean -- it carried an error-free recovery
of `1+x`'s lost bits, a separately-refined `c/u` correction with its *own*
low word, a split `LOG2_E`, and a full (not quick) two-sum -- so the port
cut 38%, which clears the 2x comfortably.

### The cancellation becomes structurally absent, not recovered

`log2p1_df` needed all that machinery for one reason: whenever `1+x`
rounds back to exactly `1`, `log2_df(u)` is exactly `Df32(0,0)` and the
`log2(1 + c/u)` correction *is* the entire answer -- so it had to carry
full precision of its own.

The atanh form deletes the problem instead of compensating it:

```text
log2(1+x) = 2*log2(e) * atanh(x/(2+x))
```

`x` appears as its own **exact factor** and the denominator is only ever
needed to *relative* accuracy, so `2+x` rounding at `2^-52` is a `2^-53`
relative error however small `x` is. There is no cancellation to recover.

The two regimes then differ in exactly one term, which is worth recording
because it is not obvious: `1+x` is exact in an f64 for every
`|x| >= 2^-29`, and `|x| < 0.25` also pins `k` to 0 -- so on that side
`m - 1` and `x` are *the same number* right up until the sum starts
rounding, at which point `x` is the one still carrying the information.
And `m + 1` is `2 + x` to a relative `2^-52` when `k` is 0, so it serves
as the denominator on both sides unchanged. The whole two-branch reduction
is therefore a **single select on the numerator**:

```rust
let s = if xd.abs() < 0.25 { xd } else { m - 1.0 };
let dd = m + 1.0;
```

### Dead code removed

With `powf` and `compound_accurate` both off the double-float route,
nothing in the crate uses it any more: `log2_df`, `log2p1_df`,
`exp2_checked_df`, `LOG2_ATANH_RCP`, `LOG2_ATANH_G`, `LOG2E_2_HI`,
`LOG2E_2_LO` and `LOG2E_LO` are all deleted -- 115 lines of `src/lib.rs`,
most of it the hand-derived `tl`-refinement and `s - 2*th` Sterbenz
argument that the f64 reduction no longer needs. `Df32` itself stays: it
is a public module and `sqrt1pm1`'s `two_prod` still uses it.

### The whole f64 lever, summarised

Three functions in one session, all the same shape -- an f32 chain doing
error-free-transform bookkeeping to reach a precision one f64 has for
free:

| | cyc/elem | latency | max ulp |
|---|---|---|---|
| `remainder_wide` | 7.688 -> 2.266 | 171.17 -> 58.13 | contract 2^48 -> 2^53 |
| `powf` | 8.459 -> 6.996 | 125.17 -> 125.81 | 3 -> 1 |
| `compound_accurate` | 9.239 -> 6.173 | 153.56 -> 120.77 | 5 -> 1 |

The screen that predicts the size of the win is **how much of the f32
version is bookkeeping rather than arithmetic**, because f64 costs exactly
2x per lane here (`vfmadd213pd %zmm` RThroughput 1.0 vs `vfmadd213ps
%ymm`'s 0.5, both 8 lanes). Halve the instruction count and it is a rout;
cut 10% and it is a wash. `sin_checked`'s earlier port and these three all
fit that rule.

### Follow-on: `dawson`'s tail needs one division, not two, and the scale is free

Same session, same claim. The entry above left `dawson` throughput +2.8%
over where it started; this takes it to **-12.2%**, i.e. the whole
function is now a win on every axis against the pre-session baseline.

`dawson`'s throughput region was **divide-bound and nothing else**: three
`vdivps` (`num/den`, `1.0/u`, `.../x`) at RThroughput 5.0 each is exactly
the 30.00 Block RThroughput mca reported, which is why the central
branch's +6 instructions in the entry above cost nothing measurable. The
tail was spending two of those three: `v = 1.0/u` for the polynomial's
argument, and `/x` for the `1/(2x)` factor.

Both come off one division. `w = 0.5/x` is the factor; `z = w*w` is
`1/(4x^2)`, i.e. the old `v` scaled by four -- and **a constant scale
folds into the polynomial's coefficients for free** (`c'_k = c_k * 4^k`),
so the branch is `dawson_tail_poly(w*w) * w` with no compensating
multiply anywhere. Re-running the quantization descent in the new
variable rather than just rescaling and re-rounding is worth 0.0005 ulp
of L1, so it is barely more than bookkeeping: L1 0.0617 -> 0.0612.

Measured (llvm-mca, `tools/mca_region.py`):

| | instrs | uOps | BlockRT | throughput | latency |
|---|---|---|---|---|---|
| two divisions | 85 | 94 | 30.00 | 1.992 | 68.97 |
| one division | 84 | 92 | **23.00** | **1.701** | **62.11** |

-14.6% throughput and -9.9% latency, and the RThroughput drop from 30 to
23 is the removed `vdivps` showing up exactly where the model says it
should. Against the pre-session baseline: throughput 1.938 -> 1.701
(-12.2%), latency 65.89 -> 62.11 (-5.7%).

**The cost, priced rather than waved through.** `R(z)*w` rounds twice
where `R/(2x)` rounded once, so this is not free in principle. In
practice it is nearly free, and the reason is worth recording as a
general screen: **the second rounding only exists where the polynomial is
not exactly 1.0.** `fma(c1, z, 1.0)` rounds to exactly `1.0` for every
`|x| > 2897` -- verified over all 1.96e9 patterns above that bound, not
argued -- so the result there is `w` itself, bit-identical to what the
two-division form produced, including where `w` goes subnormal near
`f32::MAX`. Only `|x|` in `(4, 2897]` pays, ~7% of the tail's inputs by
bit-pattern measure, at most 1 ulp each. Whole-function effect: avg
0.0561 -> 0.0585 (+4.3%), max unchanged at 5, and one `worst_corpus`
entry moved (`dawson(100)` -0.19 -> -1.19 ulp).

Not split into a second public function: 14.6% throughput against
0.0024 avg ulp and an identical max is far inside the spread `gelu_fast`
was already rejected at (7% against a 20x accuracy gap).

**Transferable:** any function whose mca Block RThroughput equals
`5.0 * (number of vdivps)` is divide-bound, and there polynomial degree
is free while a division is worth ~7 RThroughput. Look for two divisions
whose arguments are powers of the same quantity -- `1/x^2` and `1/x`,
`1/x` and `1/x^3` -- because one of them is a multiply away from the
other and any constant left over lands in the coefficients.

## Correcting my own f64 screen: one of the two sites was eyeballed, and wrong

The first version of IDEAS.md's f64-lever entry closed *both* remaining EFT
sites (`cbrt_accurate` and `clog`) as "~5 ops inside a ~60-op function, so
2x on the other 55 swamps it". That ratio was never measured for either.
Re-run properly:

- **`cbrt_accurate`: the claim holds.** Its `Df32::from_mul(y,y)` / `y2*y`
  / `(y3.0-x)+y3.1` block is **7 instructions of a 77-instruction**
  `cbrt_accurate_unchecked_throughput` region -- 9%, so the screen's
  verdict was right by luck rather than by measurement. It also has no
  accuracy lever: `cbrt_accurate` already scores avg 0.000 / max 1.
  (For the record, f64 *would* form its residual better -- one rounding of
  `yd*yd*yd` is `2^-53|x|` against the Df32 pair's `2^-48|x|`, i.e. ~2^-31.6
  vs ~2^-26.6 relative on `e` -- there is simply nothing left to spend it
  on.)
- **`clog`: the claim was wrong, and backwards.** `clog` has no mca region
  at all (`mca_target.rs` deliberately omits it), so no op-count ratio was
  ever available to eyeball. And on accuracy the shipped EFT is *not*
  exact: its correction word `(e1+e2)+es` is summed **in f32**, and those
  roundings sit at the `e`-terms' own scale (~`ulp(re^2)/2 ~ 6e-8`), not at
  `v`'s. Probed against an exact `Fraction` reference on the `|z|=1`
  manifold:

  | band | f64 `(re*re+im*im)-1.0` | shipped 10-op EFT |
  |---|---|---|
  | `\|v\| ~ 1e-9` | 5.4e-7 | 8.7e-6 |
  | `\|v\| ~ 1e-12` | 0 (exact) | 1.1e-3 |

  So four f64 ops beat ten f32 ops on the exact axis the EFT exists to
  serve. Moved to open.

### The probe that produced this nearly got it backwards too

First attempt reported the shipped EFT at **8.6e9 ulp** and f64 at 1 --
i.e. exactly inverted. The reference was
`math.log(float(R*R+I*I))/2`, and `float(R*R+I*I)` *is* the f64 method's
own intermediate, so the f64 route was being scored against itself. This is
the "harness reference can BE the bug" trap, from a file that already
documents it twice. The fix that made the comparison meaningful was to stop
scoring the end-to-end value at all and score **`v` directly against the
exact rational**, since `ln|z| ~ v/2` passes `v`'s relative error straight
through -- the same "score the thing whose error you are attributing"
discipline as scoring a shipped function rather than a reduction residual.

Transferable rule, since this is now the second time an IDEAS.md screen has
gone stale in one session: **a screen recorded as a conclusion needs the
number it was decided on written next to it.** "Fails the screen" with no
figure is indistinguishable from a guess, and the next instance cannot tell
which it was.

## `sigmoid_grad`: llvm-mca is wrong by 2.7x *in direction*, and IPC is the tell

2026-08-01. A crate-wide `tools/mca_region.py` scan over all 144
`*_throughput` regions ranked `sigmoid_grad` as the **7th most expensive
function in the crate** at 4.386 cyc/elem -- ahead of `srgb_to_linear`,
`asinh`, `acosh` and `compound`, for a three-line body. That is an
artifact. It is worth a full entry because the next person to run that
scan will see the same row and try to fix it.

### The screen said "artifact" before any hardware ran

`tanh_grad` is the same function with a different argument scale
(`exp_checked(-2|x|)` vs `exp_checked(-|x|)`, then `q/(1+q)^2` either
way), so it is an unusually clean control:

| | instrs | uOps | BlockRT | mca cyc/elem | `vdivps` | fma |
|---|---|---|---|---|---|---|
| `sigmoid_grad` | **67** | **73** | **22.0** | **4.386** | 2 | 16 |
| `tanh_grad` | 74 | 81 | 24.0 | 2.076 | 2 | 16 |

`sigmoid_grad` is smaller on *every* rung of the documented escalation
ladder -- fewer instructions, fewer uOps, lower Block RThroughput,
identical counts of the expensive ops -- and mca says it costs 2.11x as
much. That is precisely the "identical expensive ops plus fewer total
means mca is wrong" case IDEAS.md already names.

### Hardware, 4 runs, agrees with the ladder and not with mca

`jm bench` (exclusive machine lock), quickbench min-of-7, both functions
on the same `Band::Two` inputs in the same process:

| run | `tanh_grad` | `sigmoid_grad` | ratio |
|---|---|---|---|
| 1 | 1.062 ns/op | 0.840 | 0.79 |
| 2 | 1.010 | 0.903 | 0.89 |
| 3 | 1.602 | 0.902 | 0.56 |
| 4 | 0.970 | 0.861 | 0.89 |

`sigmoid_grad` is **faster in 4 of 4**, by 11-15% in the three clean
runs, against mca's claimed 2.11x slower -- and 0.85 is what the
instruction and uOp counts (0.90, 0.90) predicted. mca is out by a factor
of ~2.5 **in direction**, not just magnitude. This is the third recorded
mca-vs-hardware disagreement in this repo and the second where the sign
flips (`sincos_checked` was the first; `probit`'s was magnitude only).

### The diagnostic worth keeping: IPC inside a throughput region

`--bottleneck-analysis` names the failure outright:

```
sigmoid_grad  IPC 0.95   Data Dependencies [78.79%]  Resource Pressure [12.78%]
tanh_grad     IPC 2.23   Data Dependencies [82.57%]  Resource Pressure [39.98%]
```

A `*_throughput` region evaluates 16 **independent** elements; nothing in
it can legitimately be dependency-bound at IPC 0.95. mca has failed to
overlap its iterations and is charging one full serial critical path per
iteration -- `sigmoid_grad`'s 7017 cycles over 100 iterations is 70.2
cyc/iteration against a `sigmoid_grad_latency` of 70.079, i.e. exactly
one un-overlapped dependency chain, where `tanh_grad` gets 33 against its
own 73.

**So: compute `instrs/cycles` for any throughput region before quoting
its cyc/elem.** Ranking all 144 regions by IPC puts `sigmoid_grad` alone
among functions of its size -- every other region under 1.5 IPC is a 7-21
instruction body dominated by one `vsqrtps`/`vdivps` (`rsqrt` 0.32,
`pow_3_2` 0.66, `rhypot` 0.77, `sqrt1pm1` 0.92), where low IPC is
honest. Every region above ~50 instructions sits at 1.4 or better except
this one. That single ratio separates the honest rows from the lying one
without running anything.

No code change: the function is already the faster of the pair and there
is nothing to fix. `sigmoid_grad` and `tanh_grad` are not in readme's mca
table, so no published number needed correcting either.
## The mca *latency* column is wrong for 24 rows, and mostly too high

**Measured 2026-08-01.** Follow-up to the `asinh` entry above, which found
llvm-mca timing the wrong arm of one branch. The problem is not one
function's: **39 of the 151 latency regions in `mca_target.rs` contain a
real conditional branch, and for 24 of them the published figure lies above
*both* of that region's own arms measured in isolation.** A number outside
the interval its own two arms span is proof mca is simulating neither, and
needs no calibration to establish.

Shipped `tools/mca_arms.py` to make this checkable in one command: it
deletes the arm an in-domain input never takes and re-runs llvm-mca.

Two distinct failure modes, one measured example of each:

- **Concatenation, which overstates.** If the second arm reads a register
  the first clobbered, the two fuse into one artificial dependency chain.
  This is the common case, and all 24 outliers are in this direction.

  | fn | published | fall arm | target arm | overstated by |
  |---|---|---|---|---|
  | `exp2m1` | 80.00 | 37.00 | 48.00 | +67% |
  | `exp10m1` | 112.00 | 37.00 | 80.00 | +40% |
  | `expm1_checked` | 78.00 | 32.00 | 51.00 | +53% |
  | `expm1` | 69.00 | 32.00 | 46.00 | +50% |
  | `exp_m1_over_x` | 81.00 | 32.00 | 58.00 | +40% |
  | `erf` | 87.00 | 44.00 | 65.98 | +32% |
  | `tanh` | 85.64 | 53.00 | 62.00 | +38% |
  | `asin` | 56.74 | 26.99 | 40.99 | +38% |
  | `dawson` | 68.97 | 48.00 | 47.02 | +44% |
  | `acosh` | 89.08 | 46.08 | 77.02 | +16% |
  | `atanh` | 96.83 | 22.99 | 86.98 | +11% |

  (plus `exp10m1`, `expm1_narrow`, `exp_m1_over_x_narrow`, `asinpi`,
  `asind`, `tanpi`, `tan2pi`, `logit`, `erfc_inv`, `erfinv`, `div_euclid`,
  `srgb_to_linear`, `xlogy`, `rcbrt`.)

- **Last-writer-wins, which understates.** If both arms write the same
  register and the *cheap* one is laid out last, mca times the cheap one.
  `asinh` is the only measured case: published **69.41**, real in-domain
  chain **88.02** (+26.8%). It needs a *mixed* arm selection -- its first
  branch takes the fall-through in domain (`ax < 2048`, the sqrt arm) and
  its second takes the target (`d.is_finite()`) -- which is why a
  single-flag run brackets it (49.02/79.03) without containing it.

**Cross-checked against wall clock**, since "outside its own arms" proves
the published number wrong but not which arm is right. Calibrating
cycles-per-ns on the 92 *branchless* rows (whose mca latency has no such
problem) and applying it to the branchy ones picks the same arm the source
guard does, in every case checked: `asinh` implied 90.4 vs arm-isolated
88.02, `acosh` 74.6 vs 77.02, `tanh` 52.5 vs 53.00, `erf` 44.0 vs 44.00,
`exp2m1` 47.1 vs 48.00, `rcbrt` 63.5 vs 60.47, `ln` 45.7 vs 42.94. The
calibration is far too coarse to *set* a number (per-row ratios spread
1.85-3.42 cyc/ns) but it is decisive about which arm.

**Throughput is unaffected and remains the reference.** Those regions are
vectorized and LLVM if-converts them to masked selects, so there is no
branch to mis-simulate -- verified by the fact that all 39 branchy
*latency* regions have branchless throughput twins.

### Two corrections to entries above, both mine

- The `asinh` entry claims "~124 is honest". **It is not** -- 123.88 was
  the *modified* code's own concatenation artifact (that version still had
  2 branches per step). The honest figures are shipped **88.02** and the
  branch-deleted variant's own arm-isolated chain; the entry's conclusion
  (buried, wall clock a wash) is unaffected, but do not quote 124.
- The `srgb_to_linear` entry claims **latency -5.2%** (110.02 -> 104.30).
  Both are concatenation artifacts. Arm-isolated in-domain: **69.02 before,
  69.08 after -- flat.** The accuracy result (max 14 -> 7, avg 0.1046 ->
  0.0823, exhaustive) is unaffected.

## `srgb_to_linear_fast` was measured and is *not* worth shipping

The `srgb_to_linear` entry above offered the pre-split form as a valid
`_fast` variant, non-dominated on throughput (mca 3.831 vs 4.289, -10.7%).
**Checked on real hardware, and that advantage does not exist.** 8
alternating rounds of two prebuilt `quickbench` binaries under the bench
lock, with `atan` as an in-binary control:

```
              srgb thr min   med    | atan control min
shipped          0.898      1.020   |      0.319
pre-split        0.915      1.012   |      0.319
```

The control is identical to three places, and the shipped (split) form is
*faster* on the min. The mca +12.0% is the same 60-entry-scheduler window
effect as the `asinh` entry: `Block RThroughput` was **flat at 41** across
the change, so the resource limit never moved. Latency is flat too (24.17
vs 24.16 ns wall; 69.08 vs 69.02 cyc arm-isolated).

So the pre-split form is **dominated on every axis** -- worse avg, worse
max, no measurable throughput or latency advantage -- and under the rule
that deleted `cbrt_throughput` it does not get an API. No `_fast` shipped.
## `clog`: 4096 max ulp on the unit circle, invisible to the standing sweep

Found by re-screening my own IDEAS.md entry rather than by looking for it
(see the correction entry above -- I had closed `clog` as "fails the f64
screen" on an op-count ratio that was never measured). The re-screen said
plain f64 forms `v = re^2+im^2-1` more accurately than the shipped f32
error-free transform; checking that against the *real* function found a
defect, not just a tidier formula.

**The standing sweep cannot see this, and could not even if it sampled
there.** Two independent failures, which is why it survived the rewrite
that took `clog` from 4717 to 3 max ulp:

1. *Sampling.* `accuracy.rs` draws `re` and `im` as independent uniform f32
   bit patterns. Landing within a few ulp of `|z| = 1` essentially never
   happens.
2. *The reference.* It scores against `(re as f64).hypot(im as f64).ln()`.
   `hypot` is correctly rounded, so `|z|` carries a relative `2^-53` -- but
   `ln` of it is `~|z|-1`, which turns that into a relative `2^-53/|v|`:
   ~2 f32 ulp once `|v| ~ 1e-9`, unbounded below. **The reference is
   already wrong by more than the budget in the region being asked about.**

`examples/clogsearch.rs` walks the manifold and builds `re^2 + im^2`
exactly in f64 (both squares are exact, `two_sum` catches the one rounding
in their sum), so `v` is exact and `log1p` is the only rounding in the
whole reference. Measured on the shipped code:

| `\|v\|` band | before | after |
|---|---|---|
| `>= 1e-6` | 1 | 1 |
| `1e-7..1e-6` | 2 | 2 |
| `1e-8..1e-7` | 3 | 2 |
| `1e-9..1e-8` | **14** | 1 |
| `< 1e-9` | **4096** | **0** |

### Why the shipped EFT was not exact after all

Its own comment claimed "`s + es + e1 + e2` *is* `re^2 + im^2`, with no
error at all" -- true of the real numbers, false of the evaluation. The
correction word is spelled `((e1 + e2) + es)` and summed **in f32**, so its
roundings land at the `e`-terms' own scale (`~ulp(re^2)/2 ~ 6e-8`), giving
`v` an *absolute* error near `8e-15` -- catastrophic once `|v|` drops below
`1e-9`. An exact-as-reals decomposition still has to be *evaluated*.

### The fix: peel the `-1` off the larger component

```rust
let a = re.abs().max(im.abs()) as f64;
let b = re.abs().min(im.abs()) as f64;
let v = f64::mul_add(a, a, -1.0) + b * b;
```

Three f64 operations, no error-free transform, and `v` carries a
**relative** `2^-53` however small it gets. Which component gets the `-1`
is the whole trick and is not interchangeable:

- `mag` is in `(0.5, 1.5)` on this branch, so `a >= mag/sqrt(2) > 0.35`.
  Its exponent is at least `-2`, so `a*a` is a 48-bit number whose lowest
  bit sits at `2^-51` or above, and `fma(a, a, -1.0)` -- magnitude at most
  1.25 -- is **exact**.
- `b` has no such bound. `b*b - 1` would need bits far below `2^-53` and
  the identity fails.
- `b*b` is exact (24 bits squared is 48), and the final add is where the
  cancellation happens, so *its* rounding is `ulp(v)/2` -- relative, not
  absolute. That is the property the answer needs, since `ln|z| ~ v/2`.

Cheaper too, measured with a new `clog_re_throughput`/`clog_re_latency`
region (`mca_target.rs` had deliberately omitted `clog` on the grounds that
"its cost is just its already-measured constituents" -- no longer true once
its algebra changes):

```
                 instrs   uOps  BlockRT   cyc/elem   latency
  before           263    284    67.50     8.057     295.31
  after            258    284    61.00     7.728     281.88
```

So this is the rare case where the f64 port wins the accuracy axis *and*
every perf counter, despite f64's 2x per-lane cost -- because it deletes
ten f32 ops to add three f64 ones.

### The general lesson, and it is not about `clog`

This is the third time in one session that a *hand-picked harness bound or
sampling scheme* turned out to be load-bearing: `remainder_wide` scored
clean at `|x/y| < 2e14` and 3.42e9 ulp just past it; `clog` scores 3 under
uniform sampling and 4096 on the manifold. **A domain restriction in
`accuracy.rs` is a claim about where the function is good, and it is
evidence only if someone checked the other side.** Worth a systematic pass
over the remaining hand-picked bounds in that file -- and note that
`clog`'s case needed the *reference* replaced too, not just the sampling,
so "sample harder" is not by itself the audit.

## `norm_pdf`: the last unsplit inexact constant, and why `sinc` is not one

2026-08-01. The two-word-constant sweep that shipped `asind`, `atand`,
`atan2d`, `atanpi`, `atan2pi` and `asinpi`'s small branch was run by
*name* (`180/pi`, then `1/pi`). Re-running it by **shape** -- grep every
multiply by a named or hand-typed irrational constant, not by which
constant it is -- turns up exactly one more live site.

**Shipped: `norm_pdf`, avg 0.0320 -> 0.0269, max 4 -> 3.** `fl(1/sqrt(2*pi))`
is correctly rounded and still sits **0.4767 ulp above** the true value,
which the trailing `INV_SQRT_2PI * y` hands straight to the result as a
0.24-0.48 ulp bias with nothing downstream to cancel it. `fma(y, HI, y*LO)`
with `HI = 0.3989423`, `LO = -1.133517e-8` names the constant to a relative
`1.3e-9` ulp-equivalent. Two `worst_corpus` entries move and both are the
bias leaving: `norm_pdf(0.111)` +1.07 -> **+0.07** ulp, `norm_pdf(0.01)`
+0.68 -> **-0.32**.

The headline avg moves only 16% because ~94% of bit patterns give an
output that cannot carry the bias at all -- `|x| > 14.4` returns exactly
`0`, and tiny `|x|` returns the correctly-rounded peak. Over the ~6% that
can, removing a 0.35-ulp bias against a +-0.5-ulp rounding is worth ~0.13
ulp each, which is the 0.005 that showed up. Same "the bulk mass sits
where the answer is representable" caveat as `acosd`; the max moving 4 -> 3
is the honest signal here.

Cost: **+4 instructions, +4 uOps, Block RThroughput 24 -> 25, and latency
62.09 -> 66.09** (+4 cycles, exactly the one dependent `fma`). mca's
throughput *cycle* column went the other way (2.029 -> 1.909, -5.9%),
contradicting all three rungs above it, so the honest reading is that the
four instructions are free in a region that is not front-end bound -- not
that this made anything faster.

### `sinc` looks like the same bug and is not: the constant cancels

`sinc(x) = sinpi(x) / (f32::consts::PI * x)` has a single-word `PI` in the
denominator, `PI` is 0.47 ulp high, and `sinc`'s avg is 0.0937 -- it reads
as the identical defect. **Do not fix it.** `sinpi` computes
`sinf_poly_raw(PI * r)` with *the same constant*, so its own result is high
by the same 0.47 ulp relative, and the ratio cancels it. That is visible in
the numbers already published: `sinpi` alone measures avg 0.197 and `sinc`,
which is built on it and adds a division, measures **0.094**. Correcting
only the denominator would break the cancellation and roughly double
`sinc`'s average.

(`sinpi`'s own 0.197 is not the constant either, and both ways of removing
it are already closed in the trig section above: `two_prod(pi, r)` plus a
derivative correction left it bit-for-bit unchanged because `sinf_poly`'s
fit dominates, and folding `pi` into a dedicated fit regressed it on both
axes, 0.1969 -> 0.2065 and max 2 -> 3.)

### The rest of the sweep, for completeness

`probit`'s `SQRT_2 * erfinv(...)` is the one remaining single-word
multiply by an irrational (`fl(sqrt 2)` is 0.287 ulp low). Not done here
because `accuracy.rs` scores `probit` by a `norm_cdf(probit(p))` round-trip
residual rather than in ulp, so there is no standing measurement to show a
win against -- it needs a harness row first. `norm_cdf`'s
`xa * FRAC_1_SQRT_2` is a real 0.14-0.29 ulp bias by the same argument, but
the fix is not a two-word multiply: the constant feeds `erfcx_pos`, so it
needs `gelu`'s `RSQRT2_HI`/`RSQRT2_LO` + `erfc(z+dz) = erfc(z)*(1-2z*dz)`
treatment, which is several ops rather than one `fma`.

## `linear_to_srgb`'s max 8: the outer `1.055` is not it, twice over

**Rejected 2026-08-01.** After `srgb_to_linear` went 14 -> 7 by moving a
constant's error out of an amplifying position, the obvious next target was
`linear_to_srgb`'s own outer constant. `high = fma(1.055, p, -0.055)`
cancels near the toe -- at `l = 0.0031308`, `1.055*p = 0.0817` and the
result is `0.0267`, so anything wrong with `p` or with `1.055` is amplified
**3.06x**. And `f32(1.055)` is 4.97e-8 *low*, which is 1.2 ulp of the
result once amplified. That is a well-formed premise and it is wrong.

Two independent ways of removing that error, both measured (quick fuzz,
25M samples; shipped baseline avg 0.1209 max 8):

- **Fold `1.055` into the exponent**: `p = exp2(fma(l2, 1/2.4,
  log2(1.055)))`, `high = p - 0.055`. *Same op count* (an fma and a
  subtract replace a multiply and an fma), and the constant's error drops
  ~100x because inside the exponent it enters as `ln2 * ulp(0.0772)/2`.
  Measured **avg 0.1207, max 9** -- max worse.
- **Two-word `1.055`**: `fma(1.055, p, fma(SRGB_1055_LO, p, -0.055))`,
  naming the constant to a relative 2.7e-16, +1 fma, and keeping the exact
  `l = 1 -> 1.0` endpoint the fold puts at risk. Measured **avg 0.1228,
  max 9** -- worse on both.

**The same cancellation trap as `srgb_to_linear`'s own entry, now with two
more instances: the low `f32(1.055)` was load-bearing.** That makes three
measured cases in this crate where correcting a constant to its exact value
made a composite worse. The rule is now firm enough to state as a screen:
**in a composite whose constants have never been jointly fitted, a
systematic constant error is as likely to be cancelling as to be costing.
Measure before correcting, and correct them together or not at all.**

What is actually binding: the toe's 3.06x amplification applied to
`exp2_checked`'s *own* ~1 ulp plus the ~1.4 ulp its argument rounding
contributes (`ulp(3.46)/2 * ln2`). That is ~7 ulp with a perfect constant,
which is what is measured. Closing it needs a smaller argument, i.e. the
`srgb_to_linear` lever transplanted: `l^(5/12) = sqrt(l) * l^(-1/12)`
takes `|log2|` of the `exp2` part from 3.46 to 0.69. Not tried -- it costs
a `vsqrtps` (a real throughput item on the divider port, unlike
`srgb_to_linear`'s two extra multiplies) and `exp2_checked`'s own 1 ulp,
amplified 3x, is still a hard floor at ~3-4 max ulp, so the reachable prize
is 8 -> ~5 for a sqrt. Worth doing only if someone measures this as hot.

## `probit`: 6e4 max ulp in the tail, and the standing metric cannot see it

2026-08-01. Found while sweeping for the last single-word irrational
constant (`SQRT_2 * erfinv(fma(2,p,-1))`). The constant is not the
problem. `fma(2.0, p, -1.0)` is, and it is worth up to **60046 ulp**.

### The defect

`probit(p) = sqrt(2)*erfinv(2p-1)`. For small `p` the argument sits just
above `-1`, where `ulp` is `2^-24` -- but `p` itself carries information
down to `2^-24 * p`. Forming `2p-1` therefore **discards `log2(1/p)` bits
of the input**: at `p = 1e-4` that is ~13 bits, gone before `erfinv` is
called. `erfinv` then amplifies what is left by its own condition number,
`d(erfinv)/dx = sqrt(pi)/2 * exp(erfinv(x)^2)`, which is ~1052 there.

Measured against a Newton-refined f64 reference (seeded from the f32
answer, two steps on `x -= (Phi(x)-p)/phi(x)` with `Phi` from sleef's
`erfc_u15`; validated at `probit(0.975) = 1.959963985`,
`probit(0.025) = -1.959963985`, `probit(1e-5) = -4.264890794`, all to 9
digits against published quantiles):

| `p` band | shipped avg / max ulp | predicted by `fl(2p-1)`'s rounding alone |
|---|---|---|
| `[1e-7, 1e-5]` | 6772 / **60047** | 6767 / 55761 |
| `[1e-5, 1e-3]` | 109 / 695 | 105 / 691 |
| `[1e-3, 0.05]` | 6.22 / 36.4 | 2.64 / 18.7 |
| `[0.05, 0.25]` | 2.70 / 9.24 | 0.358 / 1.67 |
| `[0.25, 0.5]` | 0.329 / 2.62 | 0.00 / 0.00 |

The right-hand column is the first-order prediction from
`(fl(2p-1) - (2p-1)) * d(probit)/dx` and nothing else. Below `p = 1e-3` it
accounts for **essentially the whole error**; above `p = 0.25` it is
exactly zero, because `2p-1` is Sterbenz-exact there. This is the same
shape as `gelu`/`norm_cdf`'s already-recorded "the composite's *argument*
was the whole error", with a much larger amplifier.

### Why nothing caught it

`accuracy.rs` scores `probit` as a **round-trip residual**,
`max |norm_cdf(probit(p)) - p|`, which reports 1.93e-7 -- clean. It cannot
see this defect *even in principle*: `norm_cdf` is the inverse map, so it
un-does the argument rounding. `Phi(probit(p))` returns to `p` whether or
not `probit` used `p`'s low bits, because the information destroyed by
`2p-1` is exactly the information `Phi` is insensitive to at that point.
**A round-trip metric through an ill-conditioned inverse pair measures the
pair, not the function.** `erfinv` and `erfc_inv` are scored the same way
and want re-checking on the same suspicion.

### The fix, which is cheap and exact

`erfinv`'s tail branch needs `w = -ln(1 - x^2)`, and with `x = 2p-1`,

    1 - x^2 = (1-x)(1+x) = (2-2p)(2p) = 4*p*(1-p)

so `w = -ln(4*p*(1-p))` -- computable from `p` with **no cancellation at
all**, since `1-p` is exact for `p <= 0.5` (Sterbenz) and the `4` is an
exact scaling. The tail's own sign is `sign(p - 0.5)`, not `sign(2p-1)`,
which needs no subtraction either. `probit` would then reuse
`erfinv_tail_poly` on a `w` that never lost a bit, and keep the shipped
`sqrt(2)*erfinv(2p-1)` only for the central band where `2p-1` is exact
anyway. Not implemented here -- it needs the `erfinv` domain (for
`erfinv_tail_poly`) as well as `probit`'s, and a `probit` ulp row in
`accuracy.rs` to replace the round-trip that hid this.

### `SQRT_2`, the thing actually being swept for

`fl(sqrt 2)` is 0.287 ulp low, so `SQRT_2 * erfinv(...)` carries a
0.14-0.29 ulp bias exactly like the six two-word constants already shipped.
Deliberately **not** shipped: with the tail defect above unfixed, `probit`'s
own error is 3-5 orders of magnitude larger wherever anyone actually uses
it, so a two-word split here is `acosd`'s case -- structurally the better
code, no measurable effect, +3 instructions. Revisit it *after* the `w`
fix, when it would be the binding term.
## `denormal_audit`'s function list was hand-maintained, and it was missing the crate's three worst flushers

Started as the harness-bound audit the `clog` and `remainder_wide` results
called for: sweep the hand-picked domain bounds in `accuracy.rs` and find
the ones that are load-bearing. Widening nine of them (`sinc`, `sind`,
`erfc`, `erfcx`, `wrap_pi`, `cexp`'s `im`, `softplus`, `logsigmoid`,
`logaddexp`) plus the shared `hypot` magnitude window at once:

- **`hypot_checked` was a false alarm and is already right** -- it is
  passed `|_, _| true`, the full domain, not the `[1e-15, 1e18]` window
  its siblings get. Checked before assuming.
- **`sinc` (1e6 -> 1e8), `wrap_pi` (1e4 -> 1e6) and `cexp`'s `im`
  (1e4 -> 1e6) are simply conservative** -- max ulp 3, 1 and 4
  respectively outside their bounds, i.e. unchanged. No defect, and the
  bounds could be widened if anyone cares.
- **`sind`/`cosd`/`tand` (2.19e9 ulp past 4.7e7) and `hypot`/`rhypot`/
  `hypot3`/`rnorm3`/`hypot4`/`rnorm4` (5.4e8 past 1e18) are load-bearing
  *and documented*** -- the degree-reduction cliff and the
  no-rescaling-overflow tradeoff respectively, both already in the readme
  and in their doc comments. The bound is honest.
- **`softplus`/`logsigmoid` past 80 was neither.** 1.17e7 max ulp at
  `x = -87`, and the reason is a doc comment that is measurably false.

### `softplus` returns 0 where a *normal* f32 is owed

`softplus(x) = max(x,0) + log1p(exp(-|x|))`, and the correction is
selected to `0.0` once `|x| > 87` -- a threshold set by `exp_narrow`'s
`[-87.68311, 88.37627]` domain. The doc comment justified it as "once
`|x|` is far enough out that the true value is negligible at f32
precision". Measured, that is wrong: `ln(1+e^x)` does not reach zero in
f32 until `x ~ -103.97`.

| x | `softplus(x)` | true |
|---|---|---|
| -87 | 1.6458113e-38 | 1.6458115e-38 |
| **-87.3** | **0** | **1.2192433e-38** (normal) |
| -88 | 0 | 6.054601e-39 |
| -103 | 0 | 1e-45 |
| -104 | 0 | 0 |

So the flush begins while the true output is still **normal**, and covers
the entire denormal band below it -- **16.97 in `x` premature**.

### Why nothing caught it: two hand-maintained lists, not one

`denormal_audit` exists for exactly this and reported "functions that
flush some denormal output: [sigmoid, norm_pdf]". `softplus` was not in
its case-(A) list at all -- that list is **nine hand-picked functions with
nine hand-picked `x` windows**, and it had never been checked for
completeness. Adding the crate's other exponential-tailed 1-arg functions:

| function | flushed | premature by (in x) |
|---|---|---|
| `silu` | **100%** | **20.28** |
| `softplus` | **100%** | **16.97** |
| `logsigmoid` | **100%** | **16.97** |
| `sigmoid` (already known) | 94% | 15.60 |
| `gelu` | 11% | 0.134 |
| `norm_cdf` | 2% | 0.028 |

Three functions flushing their *entire* denormal range, two of them worse
than the `sigmoid` case that prompted the audit's creation in the first
place.

And there is a **second** hand-maintained list inside the same file: the
`width()` calls that produce the actionable "premature by" figure carry
their own six names. It silently omitted every function that flushes its
whole range -- precisely the set worth looking at -- because those never
produce a `last_ok_x` inside the sweep window. Both lists are extended
now.

### The fix is priced and declined, not overlooked

Restoring the tail is one substitution -- `exp_checked(-ax)` for
`exp_narrow(-ax.min(87.0))`, which also deletes the `min` and the select,
since `exp_checked` reaches `x ~ -104.7` and carries denormals:

```
                 instrs   uOps  BlockRT   cyc/elem   latency
  softplus         100    110    28.00     2.449     74.11
  + exp_checked    113    127    34.00     2.970     77.61   (+21.3% / +4.7%)
  logsigmoid       105    116    28.00     2.604
  + exp_checked    117    131    34.00     3.095     (+18.9%)
```

Declined on the same grounds `sigmoid`'s own 15.60-premature flush is
already accepted: +21% on every call to restore `-104 < x < -87`, where
the absolute value at stake is under `1.7e-38`. Recorded in `softplus`'s
doc comment with the numbers, so the next reader gets the tradeoff rather
than the old false claim. A `softplus_checked` tier is the obvious
alternative if a caller ever wants it -- the mechanism is a one-line
substitution and is written out above.

### Transferable

Both this and the `clog` result came from the same question -- *what is
just outside the bound the harness stops at?* -- and both times the answer
was hiding behind something that had been asserted rather than measured (a
doc comment's "negligible", an IDEAS.md entry's op ratio). The specific
smell worth grepping for: **a gate whose coverage is a hand-written list.**
`denormal_audit` had two of them in one file, `accuracy.rs`'s `clog` row
had a reference that could not score its own hard region, and
`codegen_check` accepted only `ps` mnemonics. None of those fail loudly;
they all just quietly report on less than they appear to.

## `probit`/`erfc_inv`: `+-inf` over ~80% of the domain, and the metric that could not see it

2026-08-01. **Fixed and landed.** The prior entry (same file, "`probit`:
6e4 max ulp in the tail") understated this by four orders of magnitude,
because it reasoned about the *rounding* of `fma(2.0, p, -1.0)` and never
asked what happens once that rounding consumes the entire argument.

### What was actually shipped

`probit(p) = sqrt(2)*erfinv(2p-1)`. For `p < 2^-25` the exact `2p-1`
rounds to exactly `-1.0`, so `erfinv(-1.0) = -inf`. Likewise
`erfc_inv(y) = erfinv(1-y)` returns `+inf` for every `y < 2^-24`. Measured
by the new direct ulp row, quick fuzz, uniform random bit patterns:

| function | avg ulp | max ulp | non-finite |
|---|---|---|---|
| `probit` (before) | 832718338 | 1053729553 | **79.53%** of samples |
| `probit` (after) | 1.5929 | 15 | 0 |
| `erfc_inv` (before) | 837462091 | 1057303482 | **79.68%** of samples |
| `erfc_inv` (after) | 1.5212 | 16 | 0 |
| `erfinv` (before / after) | 0.3855 / 0.3853 | 71 / 69 | 0 |

The true answers over that 80% are entirely ordinary: `probit(1e-45)` is
-14.12, `erfc_inv(1e-45)` is 10.02. `erfinv` is untouched by design and
its row moves only by fuzz sampling noise.

### Why nothing caught it, for eleven months

`accuracy.rs` scored all three as **round trips** --
`max |norm_cdf(probit(p)) - p|`, `max |erfc(erfc_inv(y)) - y|`,
`max |erf(erfinv(x)) - x|` -- and reported 1.93e-7, which reads clean.
A round trip through the *inverse* map cannot see this even in principle:
`norm_cdf(-inf)` is `0`, and `p` was 1e-45, so the residual is 1e-45.
The metric is structurally blind to exactly the failure it is watching
for, and it is blind *hardest* where the function is worst. `edgecheck`
missed it too -- its `probit` pins are round trips as well, and its
ordinary-value samples all sit at `p >= 0.001`.

**A round-trip residual through an ill-conditioned inverse pair measures
the pair, not the function.** The `Stats` struct already had the counter
that would have screamed (`nonfinite`, "the function blew up, not the
maths"); nothing was feeding it.

### The fix

`erfc` is odd about `y = 1`, so `n = min(y, 2-y)` lands in `(0,1]` exactly
(`2-y` is Sterbenz-exact for `y >= 1`), and `probit(p) = -sqrt(2) *
erfc_inv(2p)` with `2p` an exact scaling. Both now reduce to a shared
`erfc_inv_half(n)`, which forms

    w = -ln(1 - x^2) = -ln((1-x)(1+x)) = -ln(n*(2-n))

directly from `n`: `2n` is exact, `fma(-n, n, 2n)` is a single rounding of
the whole product, and `w` carries one `2^-25` relative error at any `n`,
denormals included. Nothing ever forms `1-n` in the tail.

That extends `w`'s reach from `[0.673, 15.94]` (all a 24-bit `x` can
encode) to `[0.673, 102.6]`, which `erfinv_tail_poly` is not fitted for --
hence a second tail poly. Fit results, LP minimax on relative error,
scored through the real f32 Estrin chain:

| variable | domain | deg 6 | deg 7 | deg 8 |
|---|---|---|---|---|
| `w` | `[16, 102.8]` | -- | -- | 1644 (ideal, deg 8 over full range) |
| `sqrt(w) - 7` | `[16, 102.8]` | 10.31 | **3.07** | 3.12 |
| `1/sqrt(w)` | `[16, 102.8]` | 0.37 | 0.70 | 0.68 |

`sqrt(w)` is the variable because `erfinv/sqrt(w) -> 1` with a `ln(w)/w`
tail no polynomial in `w` can follow over a 6.4x range. `1/sqrt(w)` fits
8x tighter again and is **not** used: it is a `vdivps` on the divider port
for a fit already 5x under this chain's binding term. `sqrt(w) - 7` is
exact for every `sqrt(w)` in `[4, 10.13]`, so the recentring is free.

The split sits at `w = 16`, i.e. exactly `-ln(1 - x_max^2) = 15.9424`,
which is where `erfinv`'s own reachable range stops. That is what leaves
`erfinv` and `erfinv_tail_poly` bit-identical: verified, `erfinv`'s mca
region is 171 instrs / 184 uOps / BlockRT 45 before and after.

Single-poly-for-the-whole-range was screened first and is dead: over
`[0.673, 102.6]` the best of five variables at degree 8 is 103 ulp
idealized (`1/sqrt(w)`), and `sqrt(w)-shift` needs degree 16 to reach 15.

### Cost

mca throughput, the trustworthy ladder:

| region | instrs | uOps | Block RThroughput |
|---|---|---|---|
| `erfc_inv_throughput` | 171 -> 209 | 186 -> 270 | 46 -> 58 |
| `probit_throughput` | 176 -> 213 | 193 -> 275 | 47 -> 60 |
| `erfinv_throughput` | 171 -> 171 | 184 -> 184 | 45 -> 45 |

+26% Block RThroughput, all three rungs agreeing in direction. No Pareto
variant was minted and none should be: the old code is not a faster point
on a tradeoff curve, it is `+-inf`.

An alternating wall-clock A/B on two prebuilt `quickbench` binaries under
`jm bench` **failed to arbitrate** and is recorded as such: across four
rounds the unchanged in-binary control (`erfc`) moved -21%, +39%, +3%
between the two binaries -- larger than the `erfc_inv` signal itself
(+7.5% to +15%, and *negative* in one round). The machine was not quiet.
This is the documented `quickbench`-absolute-numbers trap; mca stands.

### The reference, which is the reusable part

No sleef `erfinv` bucket exists. `accuracy.rs` now carries an f64 Newton
refinement against sleef's own `erf_u10`/`erfc_u15`, taking `(n, xt)` --
the erfc target and the erf target -- rather than one argument, because
each caller can supply one of the two exactly and reaches the other only
through a cancelling subtraction, and which one that is flips between the
branches. Two details are load-bearing:

- **Newton on `ln erfc`, not on `erfc`.** `erfc` is exponentially flat
  above its root and exponentially steep below it, so plain Newton crawls
  in from the left by `~1/(2z)` per step and needs dozens of iterations at
  `z = 10`. `ln erfc` is near-quadratic and lands in three.
- **The increment carries `log1p` of the *relative* residual.**
  `ln(erfc(z))` and `ln(n)` are both `~-100` near the root; differencing
  them directly cancels 1.1e-14 of absolute error into an answer that
  needs 2^-51.

Validated to 1e-15 relative against `scipy.special.erfcinv` over the whole
reachable range (`n` from 1e-45 to 1), and to 12 digits on published
quantiles. Iteration counts are not padding: 5 seed / 6 tail-Newton /
5 central-Newton, where (fp=6, tailN=2) still leaves 7726 ulp and
(fp=5, tailN=3) leaves 21.8.

### Left open

`erfinv` itself is now the worst of the three at max 69, and the cause is
*not* its poly: it forms `1-x^2` as `1 - fl(x*x)`, whose rounding is
`2^-25` absolute and therefore `1.7e-4` *relative* at the measured worst
case `x = 0.99983`. `(1-|x|)*(1+|x|)` is the same product
`erfc_inv_half` already computes. IDEAS.md #179-#181.

## `ln_normal`: the peel, and the `k` combine that was worth more than the peel

`log_2`'s leading-term peel (#202), transplanted to `ln_normal` as IDEAS.md
proposed -- and it works, but the entry's own framing was only half right.
Two independent mechanisms live in this function, they were measured
separately, and **the one nobody had proposed is the larger of the two.**

Exhaustive over all 2130706432 positive normal f32, hardware fma, scored
against `f64::ln` rounded to f32 (this is exactly `ln_unchecked`):

| chain | avg | max |
|---|---|---|
| shipped (`fma(p, s, fma(k, LN2_LO, k*LN2_HI))`, un-peeled deg-8 `P`) | 0.234319 | 3 |
| **shipped `P`, k-combine restructured only** | **0.009211** | 3 |
| peeled deg-8 `Q`, old k-combine | 0.232666 | 1 |
| peeled deg-7 `Q` + restructured combine (**shipped**) | **0.007262** | **1** |
| peeled deg-8 `Q` + restructured combine | 0.006315 | 1 |
| oracle: correctly-rounded `s^2*Q` + restructured combine | 0.005725 | 1 |

Restricted to the `k == 0` octave (all 8388608 mantissas the decomposition
can produce), which is the region `log1p`/`asinh`/`acosh`/`atanh` actually
live in:

| chain | avg | max |
|---|---|---|
| shipped | 0.453148 | 3 |
| shipped `P`, k-combine restructured | 0.453148 | 3 |
| peeled deg-7 (**shipped**) | 0.216130 | 1 |
| peeled deg-8 | 0.059757 | 1 |
| oracle | 0.042246 | 1 |

Read the two tables together and the attribution is clean: **the peel owns
the max and the near-1 octave; the k-combine owns the aggregate average.**
Neither one alone gets both. IDEAS.md predicted the first and did not
mention the second.

### The k combine: two full-weight roundings where one will do

`fma(p, s, fma(k, LN2_LO, k*LN2_HI))` was adopted for a good reason -- it
replaced `fma(p, s, k_hi) + k*LN2_LO`'s fma+mul+add with two fma, output
bit-identical. But *both* of its roundings land at the result's own scale:
the inner fma rounds `k*ln2`, which for `|k| >= 1` **is** the answer to
within `|s|`, and then the outer fma rounds the whole thing again.

    let base = fma(k, LN2_LO, s);      // off the critical path, |base| <= 0.415
    fma(k, LN2_HI, base + sq)          // k*LN2_HI is exact -> one rounding

Same three operations. `k*LN2_HI` is exact by construction (LN2_HI's low 9
mantissa bits are zero), so the closing fma is the *only* full-weight
rounding in the function; everything before it happens at `|s| <= 0.415`,
far under `ulp(result)` once `|k| >= 1`. Aggregate average **0.2343 ->
0.0073, a 32x drop, for zero operations.** For `k == 0` it changes nothing
at all, which is why it had never shown up: the near-1 region is where
everyone looks.

Three placements of `s` were measured and only one is right:

- `s` into the `k` word (**shipped**) -- 0.007262, and `base` is ready
  before the polynomial is, so the poly's critical path grows by one add
  and one fma, not two fma;
- `s + sq` first, then the two-fma k combine -- 0.006879, 14% better and
  **one dependency level deeper**. Not worth 4 cycles on every caller;
- `s` folded into the polynomial's own low group (`a = fma(s2, l0, s)`,
  one operation cheaper) -- 0.007143 but **max 1 -> 2**: it puts a second
  full-weight rounding back exactly where the peel removed one.

This supersedes the `fma(p, s, fma(k, LN2_LO, k_hi))` fold recorded in
"**`ln_normal`/`log10_normal`: fuse trailing `+k*LN2_LO` into the fma**"
above, and it is worth being precise about what it gives back. That fold
took `ln_unchecked` 38.22/1.113 -> 34.06/1.018 -- it bought both latency
*and* throughput. This change keeps the throughput (1.021, RThroughput
identical) and hands the latency back (38.06). So the fold's own ledger,
end to end, is: **throughput win kept, latency win spent on 32x of average
accuracy and a max of 3 -> 1.**

That entry also concluded that the same fold was "genuinely worse" for
`log10_normal` because `LOG10_2_LO` is 3.2x larger relative to its HI word.
**That conclusion does not carry over to the restructure here**, and anyone
citing it for `log10` should stop and re-measure: the failure it describes
is of a form that rounds `p*s + k_hi` before the LO word arrives, whereas
this one never forms anything at the result's scale until the closing fma.
Screened for `log10` (all positive normals, stride 251): shipped 0.254374
max 2 -> restructured-combine-only 0.009362 max 2 -> peeled degree 7
0.007301 **max 1**. Not shipped here; it needs the `exp10_checked` domain.

### The peel, and why degree 7

`ln(m) = s + s^2*Q(s)` with `Q(s) = (ln(1+s)/s - 1)/s`, fitted by an
ulp-weighted LP against `s^2/ln(1+s)`. This is a *better* peel than
`log_2`'s: `log_2` still has to round `s*log2(e)`, while `ln`'s leading
coefficient is exactly `1.0` and `s` is exact, so its leading term is
carried with no error whatsoever.

Degree 7, one lower than the `P` it replaces, which is what keeps the
instruction count flat. Degree 8 is a genuine Pareto point and was
rejected on measurement: aggregate 0.006315 vs 0.007262 and `k == 0`
0.0598 vs 0.2161, but +3 instructions, +4 uOps and **Block RThroughput
16 -> 17 on `ln_unchecked`**, 21 -> 22 on `log1p`, 42 -> 44 on `asinh`,
40 -> 42 on `acosh`. On the real public functions that 3.6x on the `k == 0`
average is worth only 3-14% (100M-sample quick fuzz, degree 8 vs degree 7:
`ln` 0.0031 vs 0.0036, `log1p` 0.0239 vs 0.0249, `asinh` 0.0330 vs 0.0341,
`acosh` and `logit` identical to four places, every max identical), so it buys a rounding error nobody can see for a throughput
regression on fourteen functions. Its coefficients, if a caller ever wants
that point:

    -0.499999881, 0.333333254, -0.250016093, 0.20002006, -0.166084245,
    0.141808674, -0.132478848, 0.129123241, -0.0761848763

evaluated with `w0 = fma(l2, s2, l1)`, `w1 = fma(c8, s2, l3)`,
`v = fma(w1, s4, w0)`, `sq = fma(v, s4, a)`.

### The LP was silently returning garbage until it was rescaled

Worth its own paragraph because it would have wasted the session. The
minimax LP minimises `t` subject to `w(s)*|P(s) - Q(s)| <= t`, and here
`t` is ~1e-8. **HiGHS's default primal feasibility tolerance is 1e-7**, so
every constraint is satisfied at `t = 0` and the solver returns a
degenerate vertex -- a *feasible arbitrary* coefficient set, reported with
`status == 0` and an objective of `-0.0`. The tell was that the reported
degree-8 optimum (2.6e-10) was ~500x better than the Chebyshev-ellipse
estimate for this function (~1e-7). Scaling the weights by `2^24` so the
objective is in ulp units fixes it, and the deg-8 optimum becomes 0.0685
ulp-equivalent. **Any LP in this repo whose objective is a raw relative
error needs the same scaling**; the accidentally-degenerate coefficients
still measured max 1 in the real chain, so a screen would not have caught
it either.

Sequential quantisation (fix one coefficient to f32, re-solve the LP over
the rest, try both f32 neighbours at each step -- a hand-rolled fpminimax,
which IDEAS.md #2/#102 propose in more sophisticated MIP/LLL forms) is
worth a real 13% here: degree 7 goes 0.5832 -> 0.5059 ulp-equivalent and
degree 8 goes 0.0855 -> 0.0703, against continuous optima of 0.5033 and
0.0685. Small, but free, and it is the first time the joint-rounding gap
has been measured in this repo rather than argued about.

### What it costs and what it buys

Instruction counts are **identical** in every changed `*_throughput`
region except `clog_re` (258 -> 266, register pressure in a 60-instruction
region), and Block RThroughput is identical in 13 of 14: `ln` 19, `log1p`
21, `asinh` 42, `acosh` 40, `atanh` 29, `logit` 47, `xlogy` 20, `xlog1py`
22, `erfinv` 45, `erfc_inv` 46, `probit` 47, `compound` 42 all unchanged,
`clog_re` 61 -> 62. The cycles column moves 0.3-3% on several rows with
instructions, uOps *and* RThroughput all flat, which is the documented
scheduler-window artifact.

Latency pays the peel's one extra dependency level, ~4 cycles, everywhere:
`ln_unchecked` 34.06 -> 38.06 (+11.7%), `ln` 44.14 -> 49.13 (+11.3%),
`log1p` 47.24 -> 51.25 (+8.5%), `acosh` 89.08 -> 97.88 (+9.9%), `asinh`
69.41 -> 74.49 (+7.3%), `logit` 59.91 -> 63.94 (+6.7%), `atanh` 96.83 ->
100.83 (+4.1%). Same bill `log_2`'s peel paid (+11.2%) for the same
reason, and it is a real +1 fma on the serial chain, not an mca artifact --
`ln_unchecked_latency` is branchless, so the arm caveat does not apply.

The readme's two **wall-clock** tables were deliberately not re-recorded.
They carry an explicit "recorded together in one sitting" caveat and a
measured 2.5x within-session swing on their own control function, so
appending rows measured today would break the only property they have.
`ln`/`ln_unchecked`/`log1p` are the three affected rows there and they now
understate latency by roughly the mca delta.

**No `ln_latency` variant was minted.** The pre-change code is dominated
32x on the average and 3x on the max for 10% of latency; if a caller ever
turns up that wants the old point, it is `fma(p, s, fma(k, LN2_LO,
k*LN2_HI))` over the un-peeled degree-8 `P` and it is recorded above.

### Transferable

- **A Cody-Waite combine can still round twice at full weight.** The
  split makes `k*HI` exact; it does not by itself stop the *rest* of the
  expression from being rounded at the result's scale first. Check where
  each rounding lands, not just whether the constant is split. `log10_normal`
  has the same shape (`fma(p, s, k_hi) + k*LOG10_2_LO`) and is untouched
  here.
- **`k == 0`-only and aggregate-only mechanisms are invisible to each
  other.** The k-combine fix is worth 32x on the aggregate and *exactly
  nothing* at `k == 0`; the peel is the reverse. A single-number sweep
  would have found either one and stopped. Score any reduced-argument
  function on its reduced octave *and* on the whole domain, as two rows.
- **A minimax LP whose objective is a raw relative error is inside its
  own solver's feasibility tolerance.** Scale the weights until the
  objective is order 1. The failure is silent: `status == 0`, an
  objective of `-0.0`, and a coefficient set that is merely feasible.
- **The oracle screen answers "is the fit binding", not "is the chain
  fixable".** `ln_normal`'s 2.8x headroom row was correctly closed on
  exactly that screen -- a degree-9 fit 22x better measures *worse* --
  and the function still had 32x of average and 2 ulp of max sitting in
  it. A closed headroom row means stop refitting, not stop looking.

### Verification status, stated precisely

**Exhaustive** (`accuracy thorough`, every f32 bit pattern, 4294967296
samples each): `ln` **0.117/3 -> 0.0036/1**, `ln_unchecked`
**0.235/3 -> 0.0073/1**, `log1p` **0.097/4 -> 0.0249/2**, `asinh`
**0.149/3 -> 0.0340/2**, `acosh` **0.060/4 -> 0.0031/3**, `atanh`
**unchanged at 0.0037/2**, `log1pmx` unchanged at 0.0632/3 (it has its own
polynomial). Every readme row this change touches is on this list. Plus
the standalone 2130706432-input real-chain enumeration the tables above
come from, and `worst_corpus` re-blessed -- 51 entries moved, every one
spot-checked from 1 ulp to 0.

**`acosh` is the case for having waited.** Quick fuzz reported max **2**;
exhaustive says **3**. One more entry for "quick-mode max is
systematically optimistic, not merely noisy" -- it would have gone into
the readme as a 2.

`logit` landed its exhaustive pass just after the commit: **0.0504 / max
3**, against a pre-change quick-fuzz 0.263 / max 3 -- average down 5.2x,
max unchanged. It has no readme accuracy row.

**Quick fuzz only** (100M samples; avg trustworthy to ~+-0.0001, max not):
`erfinv`/`probit`/`erfc_inv` round-trip metrics unchanged, `xlogy`/`xlog1py`/`compound`/
`clog` unchanged. None regressed on any axis. `logaddexp`'s quick max
moved 731 -> 897 between two runs of *identical* code -- that is its
documented heavy tail resampling, not this change (it routes through
`log1p_unit`, and no `logaddexp` region appears in the asm diff).
## The denormal flushers get `_checked` tiers, and the price that justified declining them was 10x too high

Follow-on to the `denormal_audit` result above, which found `softplus`,
`logsigmoid` and `silu` flushing their *entire* denormal range (16.97,
16.97 and 20.28 premature in `x`) and **priced the fix at +21.3%
throughput, then declined it**. That price was measured against the wrong
mechanism. Re-priced against the right one it is **+2.3%**, and all three
tiers now exist: `softplus_checked`, `logsigmoid_checked`, `silu_checked`.

### The mechanism: a single exponent field can carry denormals if you scale it

The obvious substitution is `exp_checked(-|x|)` for
`exp_narrow(-|x|.min(87.0))` — the correction needs `exp(-|x|)` *denormal*,
a single exponent field cannot produce one, so reach for the k1/k2 split.
That is what cost +21.3% (`instrs 100 -> 113, uOps 110 -> 127,
BlockRT 28 -> 34, 2.449 -> 2.970 cyc/elem`).

But a single field *can* produce `exp(-|x|) * 2^64`, and `2^-64` is an
exact power of two. `k = round(-|x|*log2(e))` reaches `-152` over the
extended domain, out of one field's `[-126, 127]` — `k + 64` does not.
So the whole extension is `exp_narrow`'s reduction with `64` added to `k`
and one trailing multiply, which lands the denormal with the single
correct rounding:

```
                     instrs  uOps  BlockRT  cyc/elem  latency
  softplus              100   110    28.00     2.449    74.11
  softplus + exp_checked 113   127    34.00     2.970    77.61   (+21.3% / +4.7%)
  softplus_checked      103   116    30.00     2.506    73.36   (+2.3% / -1.0%)
```

`min(105.0)` then replaces *both* `softplus`'s own `min(87.0)` and its
correction select: past `105` the trailing multiply underflows the scaled
field to exactly `0.0` by itself, which is the right answer there anyway
(`ln(1+e^x)` reaches zero at `x ~ -103.97`). Net one fewer select than the
function it extends.

Verified exhaustively: `softplus_checked` is **bit-identical to
`softplus` on all 2,237,399,042 f32 patterns with `|x| <= 87`** — a single
field is exact there, and scaling by `2^64` and back is exact while the
result is normal. It is 0 ulp at every integer `x` from `-88` to `-103`,
where `softplus` scores 4.3e6 down to 1.3.

### `silu` needed more than the substitution, and came out *faster*

`silu` has two defects, not one. `sigmoid` saturates at `x ~ -88.72` but
`silu(x) ~ x*e^x` carries an extra factor `|x| ~ 90`, so the true value
survives to `x ~ -108.6`; and for `-91.8 < x < -87.68` `sigmoid` returns a
*denormal* while `x*sigmoid(x)` is still *normal*, so the product inherits
17-22 significand bits over a band where nothing is out of range.

Both fall out of evaluating the whole quotient at the `2^64` offset —
numerator and denominator together, so `2^-64` never appears at all:

```
  x <  0:  x*e/(1+e) = (x*e2) / (2^64 + e2)
  x >= 0:    x/(1+e) = (x*2^64) / (2^64 + e2)      e2 = e^-|x| * 2^64
```

One select on the multiplier is the entire difference between the two
sides. Nothing overflows: the live domain is bounded (`|x| <= 110`), and
`|x|*e^-|x|*2^64` peaks at `6.8e18`, `|x|*2^64` at `2.0e21`.

Writing it with the explicit `2^-64` multiplies instead (`(x*e2)*P64` over
`1 + e2*P64`) is **bit-identical** — dividing by the scaled denominator is
the same value, and division by a power of two is exact — but costs 4 more
instrs and `BlockRT 17 -> 19`. Worth knowing: the scaled-denominator form
is free precision *and* free ops.

```
                     instrs  uOps  BlockRT  cyc/elem  latency
  silu                   57    60    17.00     1.407    65.02
  silu_checked           67    71    17.00     1.466    60.97
```

**Latency -6.2%** (`mca_arms.py`: the published figure equals the
in-domain arm, 13.00 is the saturated arm). Throughput is where the rungs
*disagree*: instrs +17.5% and uOps +18.3%, but **Block RThroughput is
identical at 17.00** and cyc/elem says +4.2%. Not arbitrated against
hardware, because it does not decide anything — see below.

### Why `silu_checked` is a second tier and not a replacement

It does not dominate. Over the domain the two share (`|x| <= 87`, dense
scan of every 64th bit pattern against an f64 reference):

| | avg ulp | max ulp |
|---|---|---|
| `silu` | 0.0904 | **3.85** |
| `silu_checked` | **0.0673** | 4.69 |

The usual avg-for-max trade, and structural rather than fittable: `silu`
evaluates `e^+|x|` and multiplies by a reciprocal, this evaluates `e^-|x|`
and divides. They differ on 232,916,327 of 2,237,399,042 patterns
(10.4%) by up to 7 ulp-of-result. So `silu` is the max-ulp tier and
`silu_checked` is the average/latency/tail tier — both on the frontier,
whatever the hardware would say about the throughput rung.

`softplus_checked`/`logsigmoid_checked` are a cleaner split: bit-identical
in-domain, so the *only* axis is throughput (+2.3%/+3.6%) against latency
(-1.0%/+0.6%) and the tail. `softplus` stays the default on throughput.

### Method notes

- **The quick-fuzz max for `silu_checked` printed `4` on one run and `5`
  on the next**, on byte-identical code. The dense scan's 4.69 is the
  number to quote; the readme carries both.
- **`x.max(0.0)` is not a signed-zero-safe saturation.** `silu_checked`
  needs `-0.0` on the negative tail (the true rounding, and what `silu`
  itself returns at `-1000`), and `x.max(-0.0)` gives it on both sides
  while still returning `x` for `x > 110`. `(-0.0f32).max(0.0)` is `+0.0`
  here, and Rust documents that case as *non-deterministic*, so it is not
  something to rely on either way.
- The saturation guard cannot be a clamp on the exp argument. Clamping
  freezes `e2` at a constant that `x` — unbounded — multiplies straight
  back into range (`silu(-1e30)` would return `-3.7e-26`). It has to be a
  select on the *result*, which is `sigmoid`'s own recorded trap.

### Still open, same defect, not fixed here

**`logaddexp` has the identical cutoff** (on `|a-b|` rather than `|x|`)
and nothing reports it, because `denormal_audit` covers 1-arg functions
only. `logaddexp(-88.0, 0.0)` returns `0.0`; the true value is
`6.054601e-39`, a representable denormal, and `edgecheck` already pins
`logaddexp(x,0) == softplus(x)` — the identity is the reference. It bites
only when `max(a,b)` is itself near zero, which is narrow but is exactly
the `logsumexp` normalization case. The fix is the same six lines as
`softplus_checked`; see IDEAS.md.

Note the trap while confirming it: the obvious reference
`(a.exp() + b.exp()).ln()` in **f64** returns `0.0` for `(0, -88)`,
because `1 + 6e-39` is `1.0` in f64 too. The reference collapses in
exactly the region being asked about — the `clog` result's lesson again.
Score against `softplus_checked` (or `a + ln1p(exp(b-a))`), not that.

## `tan`: one sign fold instead of two, and the `|sin|` shortcut that looks free and is wrong

`tan` was `sin_over_cos_domain(x) / cos(x)`, i.e. each half applied its
own sign combine and then the two signs were divided against each other.
Only their *difference* is observable, so the whole thing collapses into
one mask on the quotient.

The two halves share `frac_x_over_pi!` already, so they also share `n`:
the numerator reduces against `N = n + round(fc)` (the second magic
round, `qb = fc + nb`) and the denominator against `n + copysign(0.5,
fc)`. `sin(x) = (-1)^N sin(r_s)` and `cos(x) = -(-1)^n s sin(r_c)` with
`s = sign(fc)`, so

    tan(x) = (-1)^(N-n) * (-s) * sin(r_s)/sin(r_c)

and `N - n = round(fc)` is already sitting in `qb ^ nb`'s low bit. One
`vpternlogd`-shaped mask replaces two parity chains:

```
region              instrs      uOps    BlockRT   cyc/elem   latency
tan_throughput      105->101  107->104   31->30   3.153->2.985   78->78
                                                    (-5.3%)
```

`sin`, `cos`, `sin_fast`, `cos_fast` byte-identical (they are unchanged);
`tan`'s own output bit-identical over all 2^32 f32 patterns against a
literal replica of the old two-parity form. `sin_over_cos_domain` had one
caller and is gone -- its body is now `tan`'s first four lines.

**The version that is 13 instructions cheaper and wrong.** Both `(-1)^N`
factors *look* like they cancel outright, leaving

    tan(x) = sin(r_s) / |sin(r_c)|          -- no mask anywhere

which measures 105 -> 92 instrs, 107 -> 95 uOps, BlockRT 31 -> 29, 3.153
-> **2.658 cyc/elem (-15.7%)**, latency 78 -> 74, and even lets LLVM
delete the `copysign` inside the denominator's `sinf_poly` as dead. It is
wrong. The cancellation needs `sign(r_c) = -sign(r_s)`, which needs
`|r_s| <= pi/2`, which needs `N` to be the *true* `round(x/pi)` -- and at
`x = f32(pi/2) = 1.5707964` it is not. There `x/pi = 0.5000000139`, but
`fc` rounds to exactly `0.5` (the next f32 up from 0.5 is 3.3x further
away), so `qb = fc + nb` is an exact tie and round-to-even sends `N` to
`n` instead of `n+1`. `r_s` comes back as `1.5707964 > pi/2`, `r_c =
x - pi/2 = +4.37e-8` has the *same* sign as `r_s`, and the answer's sign
flips: `+2.2877334e7` against a true `-2.2877334e7`.

12 inputs over the domain do this, every one of them a pole, i.e. exactly
where the magnitude is ~2e7 and a sign flip is ~4.6e7 ulp. An exhaustive
bit-compare found them in one run; the quick fuzz did not, and would not
-- 12 in 2.5e9 is not a sampling target. The general lesson is the one
already in this file for `cospi`, one level up: **an identity derived
from "the reduction lands where it should" has to be checked at the point
where the reduction is a tie**, because that is the one place the two
halves of a fused sign argument can disagree. The shipped form assumes
nothing about `N` -- it is exact algebra for whichever `q` each half
happened to land on -- and that is what the extra four instructions buy.

**Also measured and rejected: `cos_fast`'s double-rounding fix.**
`cos_fast`'s `k = round(x/pi - 0.5)` rounds twice (`fma(x, FRAC_1_PI,
-0.5)` quantises to `ulp(x/pi)` before the magic round), which is why its
max ulp is 2780 over `|x| < 2^22*pi` against `sin_fast`'s 219, and 3
against `cos`'s 2 at `|x| <= 1e6`. Rewriting it in `cos`'s own shape --
`n = round(x*FRAC_1_PI)`, `f = fma(x, FRAC_1_PI, -n)`, `q = n +
copysign(0.5, f)` -- removes that second rounding entirely and costs the
*same four float ops in the source*. In vector form it is not free:

```
region                 instrs    uOps  BlockRT   cyc/elem   latency
cos_fast_throughput    57->64   61->68   16->16  1.406->1.526   56->57
                                                   (+8.5%)
```

+7 instructions, because `copysign` and the two-term parity are 3 more
vector ops than `!kb << 31` and they need their own broadcast constants.
That puts `cos_fast` at 1.526 against `cos`'s own 1.654 -- 8% cheaper
than the exact version for ~100x the error -- which is not a pareto point
worth having, so `cos_fast` keeps its cheap double-rounded `q`. (The
accuracy side was not measured: the cost alone disqualifies it.) The
underlying reason it cannot be done for free: the fix needs
`sign(x/pi - n)`, one fma's worth of information that the `-0.5` form
throws away, and there is no magic constant with a fractional part --
`ulp(M) = 1` forces `M` integral, so the half-odd grid is unreachable in
one round however the constant is chosen.

## `tan_checked`: the same sign fold, one level up

`tan_checked` was `sin_checked(x) / cos_checked(x)`. Both halves get
their sign from the *same* `reduce_pi64` parity `(-1)^n` (the checked
tier applies it by xoring the residual's sign bit before `sinf_poly`,
not by negating the result), and `cos_checked` adds one half-turn flip
on top. In the quotient the shared parity cancels outright, so the whole
`n + ROUND_MAGIC64` extraction is dead.

Spelled as two `reduce_pi64` calls with the masks xored once on the
quotient -- `a ^ (a ^ b)` is `b`, which instcombine takes from there:

```
region                     instrs      uOps    BlockRT   cyc/elem   latency
tan_checked_throughput   119->108   138->126   32->30   4.598->4.289  103.095->101.000
                                                          (-6.7%)      (-2.0%)
  control: sin_checked 2.495 / 82.00 and cos_checked 3.157 / 87.00 unchanged
```

Bit-identical to `sin_checked(x) / cos_checked(x)` over all 2^32 f32
patterns; all eight gates + `cargo test` clean. `clamp(-1, 1)` commutes
with the sign flip (it is odd), which is why the clamps can stay where
they are instead of moving to the quotient.

The `|sin|`-with-no-mask shortcut that is wrong for the unchecked `tan`
(entry above) is wrong here for the same reason and was not attempted:
`cos_checked`'s residual sign is not `-sign(r_s)` at a tie.

Also corrected on the way past: the readme's llvm-mca rows for
`sin_checked` (108.02 / 4.546) and `cos_checked` (113.00 / 4.037) still
carried their pre-f64-reduction numbers. Re-measured on current master:
82.00 / 2.495 and 87.00 / 3.157. `tan_checked` had no row at all and now
has one.
## `tanh`: fit the reciprocal, not the function, and the seam is free to move

`tanh` was the crate's worst-average hyperbolic: exhaustive avg 0.1452,
max 5, worst x **2.5000983e-1** -- its own seam, to four digits, exactly
where idea #169's error profile said it would be. #169 studied that band
at length, closed three levers on it, and concluded the concentration was
"structural at the current op budget". It was not. The premise none of
the three levers questioned is that the small arm has to be `expm1`'s
Pade.

**Shipped**: exhaustive avg **0.1452 -> 0.0440**, max **5 -> 2**, same
instruction count, one division instead of two.

### The mechanism as one number

The direct arm forms `e = expm1(2x)` and returns `e/(e+2)`. Compose the
`-1` cancellation gain `e^2x/(e^2x - 1)` with the divide's attenuation
`2/(e^2x + 1)` and the arm's end-to-end sensitivity to the exp chain's
*relative* error is exactly **`1/sinh(2x)`**:

| seam | 0.25 | 0.5 | 0.6 | 0.7 | 0.8 | 0.9 |
|---|---|---|---|---|---|---|
| `1/sinh(2*seam)` | **1.919** | 0.851 | 0.663 | 0.525 | **0.421** | 0.340 |

So the seam wants to sit as far out as the small arm can reach. Routing
the small arm through `expm1`'s shared Pade makes that impossible: that
Pade's argument is `2x`, so its fitted `|v| < 0.5` domain pins the seam
at **0.25** -- the single worst entry in the table. #169's lever (a)
tried to widen the Pade's own fit and correctly found it cannot be had at
that degree (42 ulp idealized at `|v| < 1`); lever (b) added a *second*
Pade and a *second* division, +39.6% throughput. Neither is needed. **A
fit in `x` has no such constraint**, and the 5-function crossover audit
that found `0.25` optimal was measuring a pair of arms that no longer
exists.

### Four small-arm forms, and why the reciprocal wins

All four are five instructions and one division; the difference is where
the roundings land. Simulated in numpy at f32 with exact fma semantics,
on a strided grid of real f32 values over `[2^-12, seam]` (7.4M points),
scored against f64 `tanh`, seam 0.7 for all four so they are comparable:

| small arm | max ulp | avg ulp |
|---|---|---|
| A `x*fma(N2,x2,1) / D` -- [5/4] rational fitted to `tanh` | 3.32 | 0.610 |
| A2 `fma(x*x2, N2, x) / D` -- peel the leading `x` | 2.46 | 0.499 |
| **A3 `x / D(x^2)`, `D ~ x*coth(x)`** (shipped) | **1.51** | **0.439** |
| C `x + x^3*c/D` -- full peel, needs its own division (7 instrs) | 1.03 | 0.252 |

A spends a full-weight rounding forming `1 + N2*x^2` and another on the
multiply by `x`. A2 folds both into one fma at the *result's* scale. A3
deletes the numerator entirely: fit `x/tanh(x) = x*coth(x) = 1 + x^2/3 -
x^4/45 + 2x^6/945 - ...` instead -- an even series, so the numerator *is*
`x`, carried exactly, and the only two roundings at result scale are
`D`'s outer fma and the division. Same lever as `ln_normal`'s peel,
pointed at a denominator. It costs nothing: `D` is degree 4 in `x^2`,
the same five instructions the [5/4] rational needed.

C is better still but cannot share the division, and **the obstruction is
general enough to write down**: the merged form
`select(small, x, 0) + select(small, n, e)/select(small, D, e+2)` returns
`+0.0` for `tanh(-0.0)`. The correction `x^3*c/D` always carries the
*opposite* sign to `x` (because `|tanh x| < |x|`), and `(-0) + (+0)` is
`+0` in IEEE round-to-nearest. All four sign conventions were enumerated
(`c`,`D` each positive/negative); the two that fix the zero invert the
sign of the correction itself. So **a peel whose final op is an addition
forces its own division.** A3 gets the same benefit for free because its
peel *is* the numerator and its final op is the divide. C was not
implemented -- +2 instructions and +1 division for one more ulp of
headroom under a max that A3 already got to 2.

### `D`'s degree is 4, and that is the whole search

Remez-LP minimax for `D`, idealized (exact arithmetic), in
ulp-equivalents of the result:

| domain | deg 3 | deg 4 | deg 5 |
|---|---|---|---|
| `[0,0.6]` | 0.470 | **0.0045** | 0.00003 |
| `[0,0.7]` | 1.536 | **0.018** | 0.0002 |
| `[0,0.8]` | 4.224 | **0.065** | 0.001 |
| `[0,0.9]` | 10.18 | **0.194** | 0.004 |

Degree 4 is an order under the evaluation rounding everywhere out to 0.9
and degree 5 buys nothing. Seam **0.8** was taken: the arm's own
simulated max is flat in the seam (1.50/1.51/1.52 at 0.6/0.7/0.8, 1.71 at
0.9), so the seam is decided entirely by the `1/sinh(2x)` column, and
past 0.9 the fit starts to move. Shipped coefficients (`x^2` ascending):
`0.33333313, -0.02221908, 0.0021010686, -0.00018176674`.

`ln_normal`'s LP trap recurred verbatim: at the raw relative scale HiGHS
returned *exact zero* residuals for degree 4 and 5 -- the answer was under
its feasibility tolerance. Rescaling the constraint rows into ulp units
(`1/5.96e-8`) made it real. Verify every LP fit on an independent dense
grid; the table above is verified at 400k points, not read off the LP.

### Measured, exhaustively (every f32 bit pattern in `tanh`'s domain, 2220710048 scored)

| | avg ulp | max | worst x |
|---|---|---|---|
| shipped before | 0.1452 | 5 | 2.5000983e-1 |
| A, seam 0.5 | 0.0588 | 4 | -5.1969117e-1 |
| A, seam 0.7 | 0.0586 | 3 | 1.4239925e-2 |
| A2, seam 0.7 | 0.0481 | 3 | 5.346194e-1 |
| **A3, seam 0.8 (shipped)** | **0.0440** | **2** | 6.2427483e-2 |

The worst-x column is the story: it walks off the old seam, then out of
the direct arm entirely. Two `worst_corpus` entries also went from 1 ulp
to *exact* -- `tanh(1e-20)` and `tanh(1e-5)` -- because below `|x| ~ 3e-4`
the whole polynomial vanishes into `D == 1.0` and the result is `x / 1.0`.

### And it is cheaper, because the select moved one step earlier

`tanh` used to divide **twice**: once inside the Pade, once in
`e/(e+2)`. Both arms are already a ratio, so the select can be taken on
the numerator and the denominator separately and the function divides
once.

| | instrs | uOps | BlockRT | cyc/elem | latency published / arms | `vdivps` |
|---|---|---|---|---|---|---|
| before | 77 | 89 | 22.00 | 1.731 | 85.64 / 53.00 / 62.00 | 4 |
| after | 77 | 88 | 21.00 | 1.759 | 63.00 / 45.00 / 62.00 | **2** |

**The cyc/elem column disagrees with every other rung and is not
believed**, per this file's own escalation ladder: instruction count is
identical, two `vdivps` (the most expensive op in the region) became two
fmas, uOps fell, and Block RThroughput fell. `mca_arms.py` on the latency
region: the small arm pays 36.02 -> 45.00 relative to form A for its two
extra dependent fmas, but is still well under the 53.00 it replaced, and
the direct arm is untouched at 62.00.

### Transferable

1. **A seam-position audit is only valid for the arms it measured.** This
   file's 5-function crossover retune found `tanh`'s `0.25` optimal, and
   it was -- for the arms that existed. Changing an arm's *shape*
   re-opens it. Same class as the `asin` 0.25 -> 0.27 entry, one level up:
   there a coefficient refit moved the crossover, here a whole arm did.
2. **The sensitivity of a composed reduction is often a closed form.**
   Here two amplification factors multiplied to `1/sinh(2x)`, which turned
   "where should the seam be" from a sweep into reading a table. Worth
   deriving before sweeping anything.
3. **Fit the reciprocal when the leading term is the argument itself.**
   `f(x) = x/D(x^2)` carries `x` exactly; `f(x) = x*N(x^2)/D(x^2)` pays
   two roundings for an `N` that a good `D` did not need. Candidates:
   any odd function whose series starts `x + O(x^3)` and that already
   divides.

## The same sign fold does *not* pay for `tand`: measured, reverted

Third application of the pattern that shipped for `tan` (-5.3%) and
`tan_checked` (-6.7%). `tand` is `sind(x) / cosd(x)`, so the shape is
identical -- but the two parities here are genuinely *independent*
(`sind`'s is `round(x/180)`'s, `cosd`'s is `round(x/180 - 0.5)`'s), so
nothing dies when they are folded. Only the bookkeeping moves: two
shifts and two xors become one `vpternlogd` + shift + xor.

```
region              instrs    uOps   BlockRT   cyc/elem   latency
tand_throughput     95->92   97->94  30->29.5  2.533->2.655  70.017->70.001
                                                  (+4.8%)
  control: sind 1.151, cosd 1.406, tand_unchecked 2.104 unchanged
```

Instructions down 3, uOps down 3, `Block RThroughput` down 0.5, and the
simulated cyc/elem **up** 4.8% with an unchanged opcode mix (2 `vdivps`,
18+4 fma, 16 `vmulps` on both sides) -- the documented signature of mca's
throughput column being wrong here, not a real regression. Reverted
anyway, and the reason is not the ambiguity: the fold has to inline
`sind`'s and `cosd`'s whole bodies into `tand`, and again into
`tand_unchecked`, i.e. four copies of a reduction that currently exists
twice, to buy at most 3 instructions out of 95 that mca will not confirm.
`tan`/`tan_checked` earned their duplication by *deleting* a shared
parity term (4 and 11 instructions); this one has no shared term to
delete.

**The transferable screen**: before folding two sign combines across a
quotient, check whether the two masks share a term. If they do
(`sin`/`cos` off one reduction, either tier), the fold deletes real work.
If they do not (`sind`/`cosd`, whose grids are offset by half a step
*before* the magic round), it only relocates three integer ops and is not
worth the copies.

## `reduce_pi_checked`'s `fc` sign extraction: `as f32` is 2 instructions cheaper and 3 cycles slower

`reduce_pi64::<true>` turns `fc`'s f64 sign into an f32-lane mask with
`(!(fc.to_bits() >> 32)) as u32 & SIGN_MASK`, which lowers to a `vpsrlq`
plus a `vpmovqd`. `!(fc as f32).to_bits() & SIGN_MASK` is the same mask
(an f64->f32 convert preserves the sign bit for every input, underflow
to `+-0.0` included) in one `vcvtpd2ps`.

```
region                     instrs      uOps    BlockRT   cyc/elem   latency
cos_checked_throughput    83->81    101->101   24->24   3.157->3.099  87.00->90.00
                                                          (-1.8%)      (+3.4%)
tan_checked_throughput   108->106   126->126   30->30   4.289->4.162  101.00->101.00
                                                          (-3.0%)
  control: sin_checked 65 / 80 / 18 / 2.495 / 82.00 unchanged (HALF=false
  never evaluates this arm)
```

Two instructions out, **uOps and Block RThroughput both flat**, so the
throughput column is the only thing claiming a win -- and it is paid for
on the latency side. The reason is that this mask is *not* off the
critical path: `sin_checked`/`cos_checked` xor it into the residual
*before* `sinf_poly`, and `fc` is only two or three ops upstream of `r`
itself, so a `vcvtpd2ps` (~4-6 cyc) where a `movq`/`shr` pair (~2-3) used
to be lands directly on the chain. `tan_checked` shows the throughput
side without the latency side precisely because the fold landed there
already moved its mask off the residual and onto the quotient.

Not taken. Recorded because the substitution looks unambiguously free on
instruction count and is not, and because it is the counterexample to
"the sign mask is bookkeeping, it cannot be on the critical path" -- in
this family it is, by construction, since the sign is applied to the
poly's *argument* rather than to its result.
## `erfcx_pos`: the reciprocal offset `2 + x` is not the conditioning lever

Follow-up to the "`erfcx_pos`'s residual is not the fit and not `v`" entry
above, which measured that even with a *perfect* `v` the f32 chain scores
3.75 against a fit of 0.72 -- i.e. three quarters of the budget is the
polynomial's own evaluation. That entry left one obvious question open,
and this closes it: the evaluation error is coefficient cancellation, so
change the variable that produces the coefficients.

`erfcx(x) = v*HI + v*P(v)` with `v = 1/(2 + x)`, so `v` spans `(0, 0.5]`
and the shipped degree-10 `P` has alternating coefficients up to 255.6
whose terms at `v = 0.5` sum to 4.26 against a true bracket of 2.0 --
a condition number of 2.1, which is exactly the size of the measured
residual. Raising the offset shrinks the interval (`a = 3` gives
`v <= 1/3`), and the endpoint condition number does improve, a lot:

| offset `a` | idealized LP fit (deg 10) | `sum abs(terms)/value` at `v = 1/a` | max abs(coeff) |
|---|---|---|---|
| 2 (shipped) | 0.121 ulp | 3.16 | 260 |
| 3 | 0.182 ulp | **1.19** | 4503 |
| 4 | 0.767 ulp | 4.67 | 6.3e5 |
| 6 | 216 ulp | 739 | 7.3e9 |

**And the real chain gets worse anyway.** Simulating the exact shipped
f32 instruction sequence (same Estrin grouping, same two-word `HI`
peel) over 3M f32 points log-spaced on `[1e-5, 30]`, each `a` with its
own LP-minimax coefficients:

| variant | max ulp | avg ulp |
|---|---|---|
| shipped, `a=2`, Estrin | 6.16 | 1.093 |
| shipped, `a=2`, Horner | 5.35 | 1.069 |
| LP refit `a=2`, Estrin | 5.73 | 1.076 |
| LP refit `a=2`, Horner | 5.07 | 1.054 |
| LP refit `a=2.5`, Estrin | 5.87 | 1.187 |
| LP refit **`a=3`**, Estrin | **7.04** | 1.382 |
| LP refit `a=3.5`, Estrin | 8.77 | 1.688 |

The endpoint condition number was the wrong proxy: it is measured at
`v = 1/a` only, and the `a=3` fit's damage is at *intermediate* `v`,
where coefficients of 1057/-4714/4126 produce terms far larger than the
sum at the endpoint. **A single-point conditioning number does not
screen a polynomial refit; simulate the chain.**

Two things this does confirm rather than close:

- **Horner beats Estrin here on both axes**, 5.35 vs 6.16 max on the
  shipped coefficients -- independently reproducing the earlier entry's
  finding with a different grid. It is not taken for the same reason as
  before: ten dependent fmas behind a division that is already on the
  critical path. It remains the one real Pareto point in this function.
- **An `a=2` refit is a ~7% max improvement (6.16 -> 5.73) at zero
  cost**, and is *also* not taken: 1.07x of idealized margin is far
  inside the weak-to-moderate band this file records failing repeatedly
  (see the LP-margin entries), and the shipped coefficients were tuned
  against the real chain rather than an idealized model.

So `erfcx_pos` stays where it is, and with it `erfc` (max 6-7), `erfcx`
(6), `norm_cdf` (8 on a quick fuzz) and `gelu` (9) -- all four are its
evaluation, not their own code. Moving them needs a *structurally*
cheaper representation of `erfcx` on `[0, inf)`, not another fit.

Incidental, and worth correcting in the record: an older note in this
file quotes `norm_cdf`'s true exhaustive max as **295**. That predates
the `erfcx_pos` rebuild; on current master the quick fuzz reports **8**
(avg 0.0996), and `norm_pdf` 3 (avg 0.0269).

## The reciprocal reorientation does **not** transfer to `dawson`, and the reason is the rule for when it does

`tanh`'s entry above ships `x / D(x^2)` in place of a rational fitted to
`tanh` itself, because the numerator is then `x`, carried exactly. Its
"transferable" note names the candidate class as "any odd function whose
series starts `x + O(x^3)` and that already divides". `dawson` is the
first member of that class checked, and it fails -- so the note needs the
sharper condition, which is this.

`dawson(x) = x*P(u)/Q(u)`, `u = x^2`, degree 6/5, both leading
coefficients pinned to `1.0`. Twelve coefficients, one division, one
multiply. If `E(u) = x/dawson(x)` were a polynomial the whole numerator
side would disappear: `dawson(x) = x / E(u)`, one polynomial and one
division, cheaper *and* with `x` exact. `E` even looks friendly --
smooth and near-linear, `E(0) = 1`, `E(16) = 30.92`, asymptotically
`2u - 1`.

Remez-LP minimax with `E(0)` pinned to `1.0`, verified on 400k points
over `|x| <= 4` (`t = u/16` normalised, without which HiGHS fails
outright above degree 6):

| degree in `u` | max relative error |
|---|---|
| 6 | 60588 ulp |
| 8 | 7736 ulp |
| 10 | 2205 ulp |
| 12 | 455 ulp |
| 14 | 68 ulp |

Against a shipped 6/5 rational whose fit is **2.8 ulp-equivalent for
twelve coefficients**. Degree 14 is still 24x worse for two more. Not
close, and not a tuning question.

**The condition, which is what to check next time.** The reorientation
works iff `x/f(x)` has no pole near the fit domain -- and its poles are
exactly `f`'s *complex* zeros. For `tanh`: `x*coth(x)` has poles at
`x = +-i*pi, +-2i*pi, ...`, distance `pi` from the origin, against a fit
domain of `|x| <= 0.8`; the ratio 0.25 is why degree 4 in `x^2` suffices
and why the coefficients decay. For `dawson`: `dawsn` has infinitely
many complex zeros and the nearest ones are close relative to a domain
that runs out to `|x| = 4`, so `E` is genuinely a rational and no
polynomial degree rescues it. The shipped rational is not a stylistic
choice.

So the screen is one line before any fitting: **locate `f`'s nearest
complex zero and compare it to the fit domain.** Real-axis behaviour says
nothing -- `dawsn` is strictly positive on `(0, 4]` and `E` is smooth,
monotone and bounded there, which is exactly what made this look like a
free win.

## `erfinv`/`erfc_inv`/`probit`: the reduction, then the poly's *variable* -- 15.6 -> 5.9 max ulp

2026-08-02, IDEAS.md #179 and #180, shipped as `5db7e4b` and `efd92b9`.

### #179, the reduction: a compensated `log1p` aimed at the wrong rounding

`erfinv`'s tail needs `w = -ln(1-x^2)` and formed it as
`log1p(-fl(x*x))` with an explicit `c/t` correction term. The correction
recovers the rounding of the *subtraction* `1 + nu`. The rounding that
dominates is `fl(x*x)`'s, which happens before `log1p` is called and which
no `log1p` can see: it is `2^-25` **absolute**, which at `x = 0.99983` is
`1.7e-4` relative to `1 - x^2 = 3.4e-4`, and the tail amplifies it by
`d(erfinv)/dw`.

`(1-|x|)(1+|x|) = n*(2-n)` on `n = 1-|x|` has no such loss and was already
`erfc_inv_half`'s reduction. Two things make it strictly cheaper as well as
strictly better here: `s = fma(-n, n, n+n)` is one rounding of the whole
product, and `s >= 2^-23` for every `|x| < 1`, so `erfinv` needs none of
`erfc_inv_half`'s `denormal_rescale!` -- and the `c/t` **division** goes
with the correction term.

Exhaustive over the positive tail band `[0.7, 1.0)`, 5033165 patterns,
against accuracy.rs's f64 Newton reference: **avg 5.0129 -> 4.9566, max
71.670 -> 11.312**. instrs 171 -> 166, uOps 186 -> 183, BlockRT 45 -> 44,
mca throughput 3.251 -> 3.212, latency 42.66 -> 41.94. Better on every
axis; no Pareto variant needed.

### #180, the poly: the degree is not the lever, the variable is

#180 proposed a degree bump on `erfinv_tail_poly`, shared by all three
functions, quoting idealized minimax 11.46 at degree 8 / 7.36 at 9 / 5.24
at 10 / 1.30 at 11. **The degree bump alone is a trap**, and the numbers
above are for the wrong range (see the next section).

`w` spans `[0.673, 16]`, a 24x range, and a monomial poly over it has
`sum|c_k w^k| / |Q|` reaching 5.7x at degree 8 and **15.6x at degree 10**.
That ratio is exactly the amplifier on each coefficient's f32
quantisation, so past some degree the fit improves and the evaluation
degrades faster. Measured through the real chain (`erfc_inv` max ulp,
every 64th pattern of its tail-poly band):

| poly | idealized (f32-quantised) | `erfinv` max | `erfc_inv` max |
|---|---|---|---|
| shipped degree 8 in `w`, Estrin | 13.65 | 11.31 | 15.56 |
| degree 9 in `w`, Estrin | 5.83 | 8.33 | 9.05 |
| degree 10 in `w`, Estrin (`w^8` top) | 3.71 | 6.71 | 13.20 |
| degree 10 in `w`, Estrin (Horner over `w^4`) | 3.71 | 7.84 | 12.32 |
| degree 10 in `w`, pure Horner | 3.71 | 5.95 | 8.85 |
| **degree 9 in `t = sqrt(w)-1`, Estrin** | **3.03** | **5.77** | **5.85** |
| degree 10 in `t`, Estrin | 2.78 | -- | (sim 4.74) |
| degree 10 in `t`, pure Horner | 2.78 | -- | (sim 4.11) |

Degree 10 in `w` idealizes 3.7x better than degree 8 and measures *worse*
than degree 9 on `erfc_inv`'s max, in two of the three evaluation orders.
Pure Horner rescues it (every intermediate stays bounded by the answer)
but costs depth 10 instead of 4: mca latency +21.3% / +9.0% / +7.8%.

In `t = sqrt(w) - 1` the span is `[-0.179, 2.993]` and the amplifier is
1.6-3.4x at *every* degree. `sqrt(w)` is free -- the `sqrt(w)*Q` combine
needs it anyway and `erfc_inv_half` was already forming `v - 7.0` beside
it -- and `v - 1` is exact for every `v` in `[0.5, 4]`, the same argument
`erfinv_far_poly`'s `sqrt(w) - 7` already makes.

| | avg before | avg after | max before | max after |
|---|---|---|---|---|
| `erfinv` (exhaustive, `[0.7,1)`) | 4.9566 | 1.5565 | 11.312 | 5.765 |
| `erfc_inv` (2.9M of its tail band) | 5.2835 | 1.4128 | 15.555 | 5.853 |
| `probit` (2.9M of its tail band) | 5.1376 | 1.3483 | 15.367 | 7.194 |

### The fit range was short by one seam, and that was 4 ulp of the 13.7

The shipped comment said the fit's domain "exactly matches what an f32
caller can ever actually reach" -- `w` up to `-ln(1-x_max^2) = 15.9424`.
True for `erfinv`. **`erfc_inv`/`probit` hand over to `erfinv_far_poly` at
exactly `ERFC_INV_W_FAR = 16.0`**, so they were extrapolating over
`(15.9424, 16]`, which is precisely where their worst case sat. The same
degree-8 coefficients idealize to **9.60 ulp over `[0.673, 15.9424]` and
13.65 over `[0.673, 16]`**. A shared poly's fit range has to be the union
of its callers' ranges, not the range of the caller it was named after.

### mca, and a clean case of its cycles column being unusable

The shipped change is +2 instructions per call on each of the three
functions and +1 Block RThroughput on each of the six regions
(+1.7-2.3%). The opcode histograms of all three *latency* regions differ
by exactly the same delta per chain element (+1 `vaddss`, +1
`vfmadd213ss`, +1 `vmovss`, -1 `vmulss`), with **zero** `jmp`/`jcc`
change -- and mca reports `erfinv_latency` **+14.1%** while reporting
`erfc_inv_latency` **-9.8%** and `probit_latency` **-8.5%**. Identical
instruction delta, opposite signs. Read Block RThroughput.

### Transferable

- **Before bumping a poly's degree, compute `max sum|c_k x^k| / |P(x)|`.**
  It is one line, it is the amplifier on the coefficients' own f32
  quantisation, and it is what decides whether a better fit converts. On a
  wide range it grows with degree faster than the fit shrinks.
- **A recentred variable is worth more than a degree**, and costs one
  subtract if you pick a centre the subtraction is exact for. Both of this
  crate's `erfinv` tail polys now use one.
- **A shared poly's fit range is the union of its callers' ranges.** The
  one that named it is not necessarily the one that reaches furthest.

## `logaddexp_checked`: the fourth correction-term flusher, and 10 stale readme mca rows

`logaddexp` carried `softplus`'s exact correction-term cutoff -- `corr = 0`
once `|a-b| > 87` -- and nothing in the crate reported it, because both of
`denormal_audit`'s lists take `fn(f32) -> f32`. Shipped as
`logaddexp_checked`, the same tier `softplus_checked`/`logsigmoid_checked`/
`silu_checked` already are, and the price came out exactly where IDEAS.md
predicted it would:

```
                     instrs  uOps  BlockRT  cyc/elem   latency
logaddexp               100   110    28.00     2.449     74.110
logaddexp_checked       103   116    30.00     2.506     73.345
                                              (+2.3%)    (-1.0%)
```

Digit-for-digit `softplus_checked`'s own +2.3% / -1.0%, which is expected
rather than a coincidence: `logaddexp(x, 0.0) == softplus(x)` identically,
so the mca rows for the two pairs measure the same instruction stream (the
house convention pins the 2nd arg of a 2-arg row to a constant). All three
rungs of the ladder move the same way, so the throughput column is not
doing anything on its own here.

Both tiers stay -- `logaddexp` on throughput, `logaddexp_checked` on
latency and on the tail -- the same Pareto split as `softplus`.

What the new tier buys, measured rather than argued:

- **2224564** values of `d = |a-b|` in `(87, 104]` where `logaddexp`
  returns exactly `0.0` and the true `log1p(exp(-d))` is a representable
  f32. It only shows when `max(a,b)` is itself near zero, which is narrow
  -- and is the log-sum-exp normalisation case exactly.
- `denormal_audit` now reports it: `logaddexp(x,0)` **FLUSHES ALL** from
  `x = -103.972`, `logaddexp_checked(x,0)` **carries denormals**. That row
  is the file's first 2-arg coverage of any kind, and it is not an
  approximation of the 2-arg behaviour: the curry is an identity.
- Bit-identity over the agreement band is **exhaustive**, not sampled.
  Both tiers build `m` and `d` identically and differ only in a correction
  that is a function of `d` alone, so all 1118699521 f32 bit patterns of
  `d` in `[0, 87]` at `m = 0` cover the whole band: **0 mismatches**.

Two things worth carrying forward:

- **The f64 reference collapses in the region being fixed.** The obvious
  `(a.exp() + b.exp()).ln()` returns `0.0` for `(-88, 0)` in *f64* too,
  because `1 + 6e-39` is `1.0` at f64 as well. `m + log1p(exp(-d))` is the
  reference that does not (accuracy.rs already used that form; the trap is
  for anyone writing a quick check by hand).
- **The 2-arg fuzz cannot see any of this.** Landing in the restored band
  needs `|a-b| > 87` *and* `max(a,b)` near zero simultaneously; four
  repeats of the thorough sweep put `logaddexp` max ulp at 1186 / 3453 /
  4213 / 20407 and `logaddexp_checked` at 934 / 1300 / 1895 / 11112 -- all
  of it the documented cancellation, sampling noise, and no signal about
  the tail either way. avg is the stable half: 0.139 vs 0.037, and the
  difference there is the wider domain (most of the finite plane has
  `|a-b| > 87`, where the answer is exactly `m` at 0 ulp), not a better
  answer anywhere.

The reduction is now a `macro_rules! exp_neg_scaled64`, shared by
`softplus_checked` and `logaddexp_checked`. Verified the way this crate
requires: a full pre/post diff of all 303 `LLVM-MCA-BEGIN` regions in
`mca_target`'s asm, normalised for label renumbering -- **0 changed
regions**, 2 new ones. `silu_checked` still has the same four lines
written out inline; it was in another domain (see IDEAS.md).

### Side finding: 10 of the readme's 79 llvm-mca rows were stale

A full `cargo run --release --example mca` diffed against the published
table, which is a ~25-minute background job and a 20-line script. The
drift is from other functions' rewrites landing without their rows:

```
                   readme lat -> real     readme thr -> real
erf                     83.98    87.00
remainder               34.11    40.13
remainder_ieee          29.11    35.13
fmod                    29.11    35.13
sinh_throughput             -        -         1.943    1.689
cosh_throughput             -        -         1.616    1.498
erfcx                       -        -         2.896    2.827
sin / asinh / erfc      sub-1% drift on one or both columns
```

`sinh_throughput` at -13.1% and `erf` at +3.6% are the two big ones. Note
which direction they point: the published `erf` latency was *better* than
the real one and the published `sinh_throughput` *worse*, so a stale table
can hide a regression and can also invent one. The two `*_checked`
sin/cos rows found stale earlier the same day were the same failure mode.
Re-measuring the table is cheap; do it before quoting any row of it.
## `log10_normal`: the peel, and a `k`-combine placement that pays for it

The `log10` half of `ln_normal`'s peel (IDEAS.md), and it lands strictly
better than `ln`'s did: `ln` bought its 3 -> 1 max for **+11.7% latency**,
this one buys the same for **zero**, because the third placement of
`k*LOG10_2_LO` tried here is one op *cheaper* and the same dependency
depth as the shape it replaces.

Exhaustive over all 2130706432 positive normal f32, hardware fma, scored
against `f64::log10` rounded to f32 (this is exactly `log10_unchecked`):

| chain | avg | max |
|---|---|---|
| shipped (`fma(p, s, k_hi) + k*LOG10_2_LO`, un-peeled deg-8 `P`) | 0.254004 | 3 |
| peeled deg-7 `Q`, `k*LO` into the poly's low group (**shipped**) | **0.006904** | **1** |

and at stride 251 over the same range, with the intermediate points:

| chain | avg | max |
|---|---|---|
| shipped | 0.254344 | 3 |
| `k`-combine restructured only, un-peeled `P` | 0.009033 | 3 |
| peeled deg 7, `fma(k, HI, fma(s, LOG10_E, sq + k*LO))` | 0.006927 | 1 |
| peeled deg 7, `fma(k, HI, fma(s, LOG10_E, fma(k, LO, sq)))` | 0.006929 | 1 |
| peeled deg 7, `s*LOG10_E` formed early into `a` | 0.008844 | 3 |
| peeled deg 7, `k*LO` into `a` (**shipped**) | **0.006913** | **1** |
| oracle: correctly-rounded mantissa term + restructured combine | 0.005619 | 1 |

`k == 0` octave, exhaustive over all 8388608 mantissas (this is the region
`log10p1` lives in): shipped 0.464848 max 3 -> peeled 0.290768 max 1. All
three peeled placements are identical there, since `k*LOG10_2_LO` is 0.

The attribution matches `ln`'s exactly: **the peel owns the max, the
k-combine owns the aggregate average.** Restructuring the combine alone is
36x better on the aggregate (0.2540 -> 0.0090) and changes *nothing* at
`k == 0`; peeling alone is what takes the max 3 -> 1.

### Where `k*LOG10_2_LO` goes, and why it is free here

`ln`'s peel puts `s` into the `k` word (`base = fma(k, LN2_LO, s)`), which
it can do because `ln`'s leading coefficient is exactly 1.0. `log10`'s is
not, so `s*LOG10_E` has to stay inside a closing fma and the tail is a
term longer. Written the obvious way that is three serial ops after the
polynomial where the old shape had two, i.e. `ln`'s +1 dependency level.

It does not have to be. `k*LOG10_2_LO` depends only on `k`, so it is ready
before the polynomial is, and it can ride into the poly's own low group:

    let a = fma(s2, l0, k * LOG10_2_LO);   // was `s2 * l0`
    ...
    fma(k, LOG10_2_HI, fma(s, LOG10_E, sq))

That is a mul + fma replacing a mul + mul + add + fma: **one op fewer than
the shipped shape**, tail back to two levels, poly depth unchanged at four.
Accuracy is if anything marginally the best of the three placements
(0.006913 vs 0.006927/0.006929) -- `a` sits at `|s^2*l0| <= 0.019` and
`|k*LOG10_2_LO| <= 6e-4`, both far under ulp(result) once `|k| >= 1`.

mca, every column, `log10`'s three regions plus `log10p1`:

| region | instrs | uOps | BlockRT | cyc/unit |
|---|---|---|---|---|
| log10_latency | 2949 -> 2885 | 3339 -> 3275 | 556.50 -> 545.83 | 48.141 -> 48.157 |
| log10_throughput | 95 -> 95 | 104 -> 103 | 20 -> 19 | 1.635 -> **1.617** |
| log10_unchecked_latency | 1866 -> 1802 | 1930 -> 1866 | 512 -> 480 | 38.219 -> **38.063** |
| log10_unchecked_throughput | 62 -> 60 | 66 -> 65 | 17 -> 16 | 1.113 -> **1.022** |
| log10p1_latency | 3279 -> 3214 | 3346 -> 3280 | 704 -> 672 | 51.986 -> **51.689** |
| log10p1_throughput | 104 -> 103 | 114 -> 112 | 23 -> 22 | 1.959 -> **1.936** |

Nothing is worse on any of the four columns except `log10_latency`'s
+0.03%, which is under the documented scheduler-window noise and has
instrs, uOps *and* BlockRT all moving the other way. So this dominates the
old function outright and no `log10_latency` variant is warranted.

Harness, `thorough` (exhaustive over all 2^32 patterns): `log10`
0.127/3 -> **0.0034/1**, `log10_unchecked` 0.255/3 -> **0.0069/1**,
`log10p1` 0.2319/3 -> **0.1458/2**. The `log10p1` baseline was re-measured
on the pre-change code in this session, not quoted -- it had no readme row.

### The fit, and the floor the peel cannot cross

`log10(m) = s*LOG10_E + s^2*Q(s)`, `Q(s) = (log10(1+s) - s*LOG10_E)/s^2`,
ulp-weighted LP against `s^2/log10(1+s)` with the weights scaled by `2^24`
(without that scaling HiGHS's 1e-7 primal tolerance returns a degenerate
vertex -- see the `ln_normal` entry), then sequentially quantised.

**This peel has a floor `ln`'s does not.** `LOG10_E` is not exact, and `Q`
cannot absorb the difference because `(log10(e) - LOG10_E)/s` is a `1/s`
term, not a polynomial one. Its weighted contribution is scale-invariant:
`w(s) * (L-A)/s -> (L-A)/L` as `s -> 0`, a fixed **0.39029 ulp-equivalent**.
Degree 7 reaches 0.5968 and degree 8 reaches **0.3903, i.e. exactly the
floor to five digits** -- degree 8 is not a Pareto point here the way it is
for `ln`, it is just the constant's own rounding, and the LP that returns
it is degenerate for the same reason. Degree 7 is shipped and is one
degree lower than the `P` it replaces, which is what pays for the peel.

Coefficients (degree 7, `c[0]` is the `s^0` term of `Q`):

    -0.21714722, 0.14476636, -0.10857988, 0.086721875,
    -0.07200416, 0.06459543, -0.06182998, 0.038040668

### Transferable

- **A peel's extra dependency level is not always intrinsic.** `ln`'s
  +11.7% latency was read as the peel's price; here it was the price of a
  particular placement of the Cody-Waite LO word. Any term that depends
  only on `k` (or only on the argument's exponent) *and is small* can be
  hoisted into the polynomial's low group, where it is off the critical
  path and costs a mul-to-fma upgrade instead of an add.
  **Checked on the two siblings, and it transfers to neither** -- do not
  re-run these. `log_2_normal` has no LO word at all (`log2(2^k) = k` is
  exact, it joins with a bare `+ k`), so there is nothing to hoist.
  `ln_normal` has one, but hoisting `k*LN2_LO` alone buys nothing: `s`
  still has to occupy the tail, leaving depth 6 and 11 ops either way,
  and `s` itself is *not* small -- folding it in is `a = fma(s2, l0, s)`,
  already measured at max 1 -> 2 in the `ln_normal` entry above. The
  lever needs a leading coefficient that is not 1.0, which is exactly
  what made `log10`'s peel look like the harder of the two.
- **A non-exact leading coefficient puts a hard floor under a peel**, and
  the floor is computable in one line before any fitting: `(f32(c) - c)/c`
  in the LP's own units. If a degree bump lands *on* that number, the
  extra degree is buying nothing and the LP producing it is degenerate.
## `dawson`'s harness: the reference was the thing that could not be checked, and the row was never exhaustive

2026-08-02. `dawson` was the last **1-arg** function in `accuracy.rs`
scored by a hand-rolled scalar sampling loop instead of `measure!`, so it
had no `thorough` mode at all -- its readme row carried the caveat
"sampled, not exhaustive". The reason was the reference: an 800-point
Simpson quadrature of `D(x) = x * integral_0^1 exp(-x^2(1-s^2)) ds`, i.e.
800 `exp` calls per sample, which capped it at 2M samples (0.05% of the
domain) and was *itself* only accurate to 1.5e-8 at `x = 4` and 8.8e-8
(0.74 f32 ulp) at its `x = 5` handover -- ~1 ulp of reference noise
across exactly the band `dawson`'s max lives in.

**Replaced by an all-positive series**, which is both far more accurate
and much cheaper:

    |x| <= 7:  D(x) = exp(-x^2) * sum_n x^(2n+1) / (n! * (2n+1))
    |x| >  7:  the double-factorial asymptotic series, 20 terms

Nothing cancels anywhere in the sum -- the `exp(-x^2)` that undoes its
growth is applied once, at the end -- so its relative error is just the
summation's. Against `scipy.special.dawsn`: **<= 6.2e-16 at every x from
1e-5 to 10** (under 0.006 f32 ulp) and **0.0** for the asymptotic arm at
x = 7, 8, 10, 20, so the handover is clean from both sides.

Three implementation notes, each a real bug on the way in:

- **`simd_min`/`simd_max` drop NaN** (IEEE minNum/maxNum), so clamping
  `|x|` for the two arms made *both* come back finite for a NaN input and
  the row scored every NaN pattern as a non-finite blow-up (0.39% of
  samples, `max ulp 18446744073709551615`). NaN needs its own explicit
  arm.
- **`/ (2n+1)` is a vector division per lane per term.** Hoisting it to
  `* (1.0/(2n+1))` -- one scalar division per term -- took the quick row
  from 77s to 18s. Its own error is a relative 1.1e-16 per step on an
  all-positive accumulation, under 1e-6 f32 ulp after all 200 terms.
- The clamps also stop the discarded arm overflowing: without
  `min(ax, 7)` the sum reaches `inf` for large `|x|` and `exp(-inf)*inf`
  is NaN.

### What the exhaustive pass then said: the old number was right, and one ulp light

| | avg ulp | max ulp | worst x |
|---|---|---|---|
| old row, 2M samples, Simpson reference | 0.059 | 5 | -- |
| **exhaustive, all 2^32 patterns, series reference** | **0.0581** | **6** | 1.4041231 |

495s for the whole domain, i.e. cheaper than the 2M-sample loop it
replaces. So the sampling was not hiding anything structural -- but the
max was understated by one, and that one is exactly the reference noise
the old comment described.

### Where the remaining 6 ulp sits, so the next attempt does not re-derive it

Four-row decomposition over `[1e-3, 4)` at stride 4 (8.4M points),
against the same series reference (this is the method the earlier dawson
entry used, re-run on the current code):

| what is being scored | avg | max |
|---|---|---|
| shipped f32 chain | 0.6595 | 5.694 |
| the same rational in f64 from `fl(x*x)` | 0.4571 | 4.322 |
| the same rational in f64 from an exact `x^2` | 0.4532 | 3.517 |
| oracle `x * fl(dawsn(x)/x)` | 0.1332 | 1.293 |

So of the 5.7: **~2.2 is the [6/6] fit**, ~0.8 is `u = fl(x*x)`, ~1.4 is
the f32 evaluation chain, and ~1.3 is the irreducible `x * fl(R)` pair of
roundings. The fit is still the largest single term but no longer
dominant, which is why the earlier entry's "the frontier is flat" holds.

### Rejected: a real-chain coordinate descent of the 12 free coefficients

The lever that has worked repeatedly here -- descend the f32 coefficients
against the *shipped chain* rather than against the idealised fit -- was
run over `pc[1..6]`/`qc[1..6]`, +-8 ulp each, 6 passes, scored
exhaustively at stride 16 over `[1e-4, 4]` (8.4M points), lexicographic
on (max, avg). It converges after two passes to **max 5.481 -> 5.103
(-6.9%) for avg 0.5536 -> 0.5753 (+3.9%)**: the same avg-for-max trade
this file has now recorded 4/4 times for a max-objective refit of an
already-tuned poly. Not shipped -- 7% of a max nobody is binding on, paid
for in the average, on a poly whose degree was already chosen against a
flat frontier.

## `asin`: the crossover was the whole defect; 5 -> 2 max ulp, and the big poly loses a term

`asin`'s error was not distributed. With the `0.27` crossover, **every one
of the 621495 f32 scoring >= 3 ulp had `a = |x|` in `[0.27, 0.4997]`**, and
the big branch measured max 2 over the entire rest of its domain:

```
a in [0.27,0.3)  max 5  avg 1.189      a in [0.5,0.7)  max 2  avg 0.505
a in [0.3,0.4)   max 5  avg 1.169      a in [0.7,0.9)  max 2  avg 0.276
a in [0.4,0.5)   max 5  avg 1.040      a in [0.9,1.0)  max 1  avg 0.238
```

Two things end at exactly `0.5` and both were in that window. The big
branch is `pi/2 - sqrt(1-a)*P(a)`, whose subtraction amplifies the
product's rounding by `|sqrt(1-a)*P(a)| / asin(a)` -- 2.3x at `a = 0.5`,
4.7x at `0.27`. And `1.0 - a` is itself inexact below `0.5`; Sterbenz makes
that subtraction exact only from `a >= 0.5` up. So the fix is not to refit
the big poly, it is to stop calling it there.

**Shipped: crossover 0.27 -> 0.5, `asin_small` degree 3 -> 5 in `x^2`,
`asin_poly` degree 6 -> 5.** Exhaustive over all 2^32 patterns:

```
              avg ulp            max ulp
asin      0.0188 -> 0.0158       5 -> 2
asind     0.0222 -> 0.0195       9 -> 4
asinpi    0.0159    unchanged    5    unchanged   (control: own polys)
```

mca, all rungs agreeing: instrs 55 -> 58, uOps 60 -> 64, Block
RThroughput 13 -> 14, throughput 0.900 -> 0.961 cyc/elem (+6.8%). Latency
*improves*, which the published column cannot show -- via
`tools/mca_arms.py`, the small arm goes 26.99 -> 34.99 and the big arm
40.99 -> **36.99**, and the big arm binds on both sides, so the published
fused 56.74 -> 60.91 is the usual two-arm concatenation artifact and not a
regression.

### The two beliefs this had to get past, and why both were wrong

- **"Extending `asin_small` needs ~12 terms to reach 0.5"** (this file,
  in the entry rejecting a mid-range third branch). That is the *Taylor*
  series' convergence, and asin's decays by only `x^2` per term. An
  ulp-weighted minimax over `[0, 0.5]` needs **degree 5**: 0.081
  ulp-equivalent idealized, against 1.49 at degree 4 and 29.2 at the
  shipped degree 3. The same entry priced the alternative (a third,
  mid-range branch) at "a whole extra poly evaluated unconditionally plus
  a select" and rejected it on cost -- correctly, but the cheaper
  restructuring was never priced.
- **"One more term in `asin_small` costs +7.9% throughput"** (this file,
  a separate entry, and it was measured). Two more terms cost +14.2% on
  their own -- but the crossover move *narrows `asin_poly`'s domain to
  `[0.5, 1)`*, and that is a real licence: degree 5 there, ulp-weighted LP
  then coordinate-descended over the f32 grid, measures **max 2 / avg
  0.324** through the real chain over every f32 in `[0.5, 1)`, against the
  degree-6 predecessor's **max 2 / avg 0.360** on the same inputs. Better
  on both axes with a coefficient removed, because that coefficient was
  paying for `[0.25, 0.5)`. Net +1 fma, not +2.

### The oracle screen, recorded because it points somewhere else

Before any of the above, the real chain was scored over all 16106127 f32
in `[0.27, 1)` with `P` replaced by a correctly-rounded oracle:

```
shipped poly                          max 5  avg 0.722
oracle P = fl(acos(a)/S)              max 4  avg 0.995
oracle P = fl((fl(pi/2) - asin(a))/S) max 2  avg 0.472
```

The gap between the two oracles is `fl(pi/2)`'s own representation error,
`4.371e-8`, which is **1.47 ulp of the result at `a = 0.27`** and 0.37 at
`a -> 1`; the fit target `acos(a)/sqrt(1-a)` does not contain it, so a
perfect fit to that target cannot beat max 4. Two consequences worth
carrying:

- Any future `asin_poly` refit must target `(fl(pi/2) - asin(a))/sqrt(1-a)`,
  not `acos(a)/sqrt(1-a)`. This is very likely why this file's earlier
  "ulp-weighted minimax refit of `asin_poly`" measured **max 8** through
  the real chain: the LP was solving for the wrong function.
- Do **not** "fix" it by adding a `PI_2_LO` correction term. Measured:
  `+ PI2_LO` after the fma takes the shipped chain from max 5 to **max 6**
  (avg 0.722 -> 1.135), and a hypothetical exactly-rounded `pi/2` inside
  the fma gives max 6 / avg 1.154. The shipped coefficients were
  real-chain coordinate-descent tuned and had already absorbed most of the
  bias; removing it explicitly breaks that compensation without replacing
  it.

The corrected-target LP also has a floor: `(fl(pi/2) - asin(a))/sqrt(1-a)`
carries a `4.371e-8/sqrt(1-a)` singularity at `a = 1` that no polynomial
tracks, and it pins the idealized minimax at **0.3658 ulp-equivalent for
every degree >= 5**. That is exactly `4.371e-8 / ulp(pi/2)`, i.e. the bias
itself, and it is why degree 6 buys nothing over degree 5 on `[0.5, 1)`.

### Still open, same shape

`asinpi` has the identical construction and the identical `0.27`
crossover, and its exhaustive worst case is at `x = 0.27000788` -- max 5,
unchanged by this work because it has its own `asinpi_small`/`asinpi_poly`
pair. It is a different domain; nobody held it at the time.
## `exp_pos_neg_core`: cosh's max was the joint poly, and one fma separates the two sinh/cosh tiers

2026-08-02. `sinh` and `cosh` sat at max 4 while `exp` itself is 3, which
is backwards: `cosh = ep + en` **adds** two same-signed halves, so it
cannot amplify their relative error at all, and `sinh`'s own worst case
was recorded at `x = 4.67`, where `coth(x) = 1.0002` -- no cancellation
either. Neither headline max was the reconstruction. Both were the shared
poly.

### Oracle screen, and it is the unusual verdict

Replacing `exp_pos_neg_core!`'s `(p_pos, p_neg)` with a correctly-rounded
`f32(e^(+-r)/2)` from f64 and changing nothing else, over `[0.5, 80)` at
stride 8 (2.9M points, against `f64::cosh`/`sinh`):

| | cosh avg / max | sinh avg / max |
|---|---|---|
| shipped (degree 5) | 0.6649 / 3.629 | 0.6821 / 3.581 |
| the same poly *shape* evaluated in f64 | 0.5810 / 2.576 | 0.5989 / 2.505 |
| oracle, correctly-rounded `e^(+-r)` | 0.1031 / 1.164 | 0.1302 / 1.665 |

So the fit is 0.58 avg / 2.6 max of it and the f32 evaluation only 0.08 /
1.05. Headroom, and real.

### The fit is one degree-6 minimax of `e^r/2` read two ways

`p_pos = e + r*o` and `p_neg = e - r*o` is not two fits: `e` is the even
part of a single polynomial approximation of `e^r` and `r*o` its odd part,
so the pair's *joint* relative error against `e^(+-r)` is exactly the
relative error of one polynomial `P(r) ~ e^r/2` over `|r| <= ln2/2`, read
at `+r` and `-r`. Written that way it is an ordinary 1-D minimax problem
with the `r^0`/`r^1` coefficients pinned to 0.5, and "degree 6" means one
more coefficient in `e` **alone** -- one fma, no new multiply, no term in
`o`.

Idealized relative error, ulp-scaled LP (HiGHS returns exact-zero
residuals without the ulp scaling here -- the same trap `ln_normal` and
`tanh` both hit), then f32-quantised and coordinate-descended:

| | degree 5 | degree 6 | degree 7 |
|---|---|---|---|
| idealized ulp | 1.76 | **0.053** | 0.004 |

The shipped degree-5 coefficients idealize to 2.03, so a *same-degree*
refit is worth only 1.14x and measures worse on cosh's max (3.63 ->
3.730) -- the usual avg-for-max shuffle. Degree 6 is where the fit stops
being the binding term: it measures 2.373 / 2.393 against the oracle's
1.164 / 1.665 floor, so degree 7 has nothing left to buy and is not
taken.

### Shipped, and what it costs

Every caller of `exp_pos_neg_core!` and its standalone copy in
`exp_pos_neg_narrow_half` pays exactly +1 fma:

| region | instrs | uOps | BlockRT | cyc/unit |
|---|---|---|---|---|
| `sinh_throughput` | 88 -> 91 | 96 -> 100 | 26 -> 27 | 1.720 -> 1.824 |
| `sinh_latency` | 3210 -> 3275 | 3211 -> 3277 | 768 -> 800 | 51.00 -> 54.00 |
| `cosh_throughput` | 70 -> 73 | 72 -> 76 | 22 -> 23 | 1.606 -> 1.695 |
| `cosh_latency` | 2631 -> 2696 | 2631 -> 2696 | 640 -> 672 | 50.00 -> 53.00 |
| `sinh_checked_throughput` | 93 -> 96 | 103 -> 107 | 28 -> 29 | 1.974 -> 2.094 |
| `cosh_checked_throughput` | 74 -> 77 | 80 -> 83 | 24 -> 25 | 1.938 -> 2.086 |
| `sinh_narrow_throughput` | 72 -> 75 | 76 -> 80 | 20 -> 21 | 1.383 -> 1.524 |
| `cosh_narrow_throughput` | 54 -> 57 | 56 -> 59 | 16 -> 17 | 1.150 -> 1.272 |
| `coshm1_throughput` | 108 -> 111 | 116 -> 120 | 33 -> 34 | 2.529 -> 2.591 |

Uniform: +3 instructions, +1 Block RThroughput, +3-6%. All three ladder
rungs agree, so this one is real and there is nothing to arbitrate.

### Why the trade is worth taking here: the speed tier already exists and does not pay it

`sinh_throughput`/`cosh_throughput` route through two independent `exp`
calls, not `exp_pos_neg_core!`, and their mca regions are **byte-identical
before and after** (79/84/24.00/1.689 and 62/64/20.00/1.498 on both
sides). So the +1 fma lands only on the accuracy tier, and it is what
finally separates the two -- before this, `sinh` was better than
`sinh_throughput` on latency alone and *tied* on accuracy, which is a
weak Pareto pair:

| | throughput | latency | avg ulp | max |
|---|---|---|---|---|
| `sinh_throughput` (unchanged) | **1.689** | 62.00 | 0.0795 | 4 |
| `sinh` before | 1.720 | 51.00 | 0.078 | 4 |
| `sinh` after | 1.824 | **54.00** | **0.061** | **3** |

Check for an existing speed tier before pricing a shared-poly degree
bump: if the fast path does not route through the poly, the blast radius
is exactly the functions that wanted the accuracy.

## `asinpi`: the same crossover fix as `asin`, and the transfer was exact

Flagged in IDEAS.md the moment `asin` landed: same two-branch shape, same
`0.27` crossover, exhaustive worst case at `x = 0.27000788`. The screen
asked for there (is the `>= 3` ulp population confined to `[0.27, 0.5)`?)
came back yes, with the same profile:

```
big [0.27,0.3) max 5 avg 1.373      big [0.5,0.7) max 2 avg 0.541
big [0.3,0.4)  max 5 avg 1.243      big [0.7,0.9) max 2 avg 0.189
big [0.4,0.5)  max 3 avg 0.676      big [0.9,1.0) max 1 avg 0.102
```

Shipped identically: crossover `0.27 -> 0.5`, `asinpi_small` degree 3 ->
5 in `x^2` (the `FRAC_1_PI_HI/LO` peel is untouched, the two new
coefficients just extend the tail Horner chain), `asinpi_poly` degree 6 ->
5 on the narrowed `[0.5, 1)`. Exhaustive over all 2^32 patterns:

```
             avg ulp             max ulp
asinpi   0.0159 -> 0.0057        5 -> 2
asin     0.0158    unchanged     2    unchanged   (control)
asind    0.0195    unchanged     4    unchanged   (control)
```

The avg is a 2.8x improvement, larger than `asin`'s own 1.19x, because
`asinpi`'s small branch already carried the two-word `1/pi` peel: this
hands that branch another 46% of the domain by measure.

mca: instrs 61 -> 63, uOps 68 -> 70, Block RThroughput 14 -> 15,
throughput 0.981 -> 1.059 cyc/elem (+8.0%). Per-arm latency, which is the
one difference from `asin`: 31.99 -> 39.99 small, 40.99 -> **36.99** big,
so the binding arm improves only 40.99 -> 39.99 and *swaps sides*. The
peel's trailing `fma(x, FRAC_1_PI, x*t)` is what makes the small arm two
fma deeper than `asin_small`'s at the same degree, and it is now the long
one. Anything further added to that branch will cost latency at full
price; `asin`'s will not until it too passes 36.99.

`asinpi_poly` at degree 5 over `[0.5, 1)`, ulp-weighted LP then
coordinate-descended over the f32 grid: max 2 / avg 0.308 through the real
chain over every f32 there, against the degree-6 predecessor's max 2 / avg
0.312. Same "the removed term was paying for `[0.25, 0.5)`" result as
`asin_poly`'s, and worth noting it reproduced without retuning anything by
hand -- the LP + descent recipe transferred as written.

`worst_corpus` moved 5 entries, all `asinpi`, each checked against f64:
two go 1 ulp -> correctly rounded, one goes correctly rounded -> 1 ulp,
and the pair-mirrored negatives. A refit moving individual corpus points
both ways while the exhaustive aggregate improves 2.8x is the expected
shape, not a warning sign.

### Where this stops

`acospi` is *not* the third instance. It is a single expression with no
crossover at all (`0.5 - sqrt(1-a)*P(a)` over the whole domain, trailing
constant exactly `0.5`), max 3 / avg 0.0438, and there is no small branch
to hand a window to. `acos` likewise. The family's remaining crossover is
`asind`'s, and that one is `asin`'s -- already moved.

## `erfcx_pos`: Horner is not the Pareto point the grid simulation said

The entry above ("the reciprocal offset `2 + x` is not the conditioning
lever") closed with Horner named as "the one real Pareto point in this
function" -- 5.35 vs Estrin's 6.16 max on a 3M-point log-spaced grid, at
a latency cost that had never been measured. Both halves are now measured
in the real chain, and **it is not a Pareto point at all**: it is worse on
accuracy, worse on latency, and it breaks an edgecheck pin.

Exhaustive, every f32 bit pattern, shipped coefficients, only the
evaluation order changed (degree-10 `P(v)` as ten dependent fma instead of
the 4-deep Estrin tree, closing `fma(v, HI, v*P)` untouched):

| function | shipped Estrin | Horner |
|---|---|---|
| `erfc` | **0.1948** / 7 | 1.3899 / 7 |
| `erfcx` | **0.2080** / 6 | 1.4028 / 6 |
| `erfcx (|x|<=20)` | **0.2078** / 6 | 1.3935 / 6 |
| `erf` (control, no `erfcx_pos`) | 0.0270 / 3 | 0.0270 / 3 |

**Seven times the average, and not one ulp of max.** The mechanism is
named in `erfcx_pos`'s own comment and is the thing the grid simulation
could not see: the coefficients are coordinate-descent polished *against
the exact Estrin evaluation order*, under a side constraint that `xa = 0`
reproduce `erfcx(0) = 1.0` bit-exactly, and `c1..c4`/`c10` each sit an ulp
off the LP's own values to hold that constraint. Reorder the evaluation
and the polish is not merely wasted, it is actively wrong -- the residual
it was cancelling is no longer there. `edgecheck` says so directly and
immediately:

    FAIL erfc(0)      got 1.0000001e0  want 1e0
    FAIL erfcx(0)     got 1.0000001e0  want 1e0
    FAIL norm_cdf(0)  got 4.9999994e-1 want 5e-1

A uniform ~1 ulp bias at the pin is exactly the ~1.4 average above.

**mca, and the reason the throughput column cannot be believed here.**
Horner is genuinely *fewer* instructions -- it drops `v2`/`v4` and two
combining fma -- and instrs, uOps and `Block RThroughput` fall on every
one of the eight regions. The cycles column still goes up on three of
four throughput rows:

| region | instrs | uOps | BlockRT | cyc/unit |
|---|---|---|---|---|
| erfc_latency | 4551 -> 4298 | 4721 -> 4423 | 1184 -> 1120 | 62.08 -> 83.33 (+34.2%) |
| erfc_throughput | 131 -> 123 | 148 -> 139 | 40 -> 38 | 2.885 -> 3.007 (+4.2%) |
| erfcx_latency | 4129 -> 3745 | 4869 -> 4962 | 1184 -> 1120 | 66.99 -> 79.99 (+19.4%) |
| erfcx_throughput | 125 -> 118 | 143 -> 135 | 39 -> 37 | 2.827 -> 2.761 (-2.3%) |
| norm_cdf_latency | 4980 -> 4649 | 5169 -> 4862 | 1280 -> 1216 | 70.22 -> 92.47 (+31.7%) |
| norm_cdf_throughput | 140 -> 132 | 160 -> 151 | 43 -> 41 | 3.219 -> 3.297 (+2.4%) |
| gelu_latency | 5773 -> 5619 | 6157 -> 6044 | 1504 -> 1440 | 75.56 -> 98.94 (+30.9%) |
| gelu_throughput | 156 -> 149 | 184 -> 173 | 50 -> 48 | 3.574 -> 3.577 (+0.1%) |

The latency column is the real one here and it is not an artifact: ten
dependent fma replacing a 4-deep tree is +6 levels by construction, these
regions are branchless, and +19-34% is far outside the documented noise.
So even *with* a fresh Horner-order polish -- which is the only way this
idea could be revived -- the trade on offer is at best a few percent of
max for a third of the latency on four public functions.

### Transferable

- **A polished poly's evaluation order is part of the coefficients.**
  Nine of eleven coefficients here are holding a residual specific to one
  instruction sequence. Any Estrin/Horner/reassociation experiment on a
  poly whose comment says "coordinate-descent polished" or "tuned against
  the real chain" must refit before it measures anything, and a grid
  simulation that reuses the shipped coefficients will report the
  reordering as free when it is not.
- **`cargo build --release` does not rebuild `examples/`.** The first pass
  of this measurement ran `./target/release/examples/accuracy` against a
  binary built before the edit and reported the change as bit-identical on
  `norm_cdf`/`gelu` -- same avg to four places, same max, same worst `x`.
  That is the tell, and it is the same tell as a byte-identical asm
  region: if a change that must move numbers moves none, check what you
  measured before you conclude anything. Use `cargo build --release
  --example accuracy`, or `cargo run`, and check the binary's mtime.

### Also: `erfc`'s readme max was 6 and is 7

Not a regression -- an exhaustive re-measure of unchanged code. Full
2^32-pattern figures for the whole `erfcx_pos` family as of this
session, for the next instance to diff against rather than re-run:
`erf` 0.0270/3, `erfc` **0.1948/7**, `erfcx` 0.2080/6, `erfcx (|x|<=20)`
0.2078/6, `erfcx (x>=20)` 0.2683/2, `norm_cdf` 0.0996/8, `norm_pdf`
0.0269/4, `gelu` 0.2029/9.

## `erfcx_pos`: the leading-term peel does not transfer

This crate's highest-hit-rate accuracy lever (`log_2`, `ln`, `log10`,
`coshm1` -- the last one peeling *two* leading terms, same as tried here)
applied to the function the previous two entries identify as the binding
term for `erfc`/`erfcx`/`norm_cdf`/`gelu`. **It buys ~1%.** Recording it
because the diagnosis in those entries -- "roughly three quarters of the
remaining budget is the f32 evaluation of the polynomial itself" -- reads
like an invitation to peel, and it is not one.

`erfcx(xa) = v*HI + v*P(v)`, so the same polynomial regrouped is
`v*HI + v*c0 + v^2*c1 + v^3*S(v)` with `S = c2..c10`: the poly's own
evaluation roundings then reach the answer scaled by `v^3` (0.125 at the
worst point) instead of `v` (0.5).

numpy simulation of the exact f32 instruction sequence (f64-emulated fma),
80000 points -- an even sweep of `v` over `(0, 0.5]` plus a log sweep of
`xa` over `[1e-5, 30]`. Run twice, because the shipped coefficients are
coordinate-descent polished *for the Estrin order* and scoring a regrouping
against them is the trap the Horner entry above documents; the second block
uses a fresh ulp-weighted LP fit (0.2345 ulp-equivalent after plain
round-to-f32, no polish) so both groupings are treated alike:

| grouping | shipped polished coeffs | fresh LP coeffs |
|---|---|---|
| shipped Estrin | 5.336 max / 0.8154 avg | 4.980 / 0.7932 |
| peel `c0` and `c1` | 5.241 / 0.8051 | 4.942 / 0.7829 |
| peel `c1` only | 5.241 / 0.8096 | 4.906 / 0.7893 |
| shipped Estrin, **exact `v`** | -- | 3.371 / 0.5681 |
| peel `c0`+`c1`, **exact `v`** | -- | 3.111 / 0.5572 |

**0.8% on max with the real `v`, and 7.7% even with a perfect one** --
against +1 operation and +1 dependency level on four public functions.
Not taken, and not worth a refit-and-repolish cycle to confirm at higher
precision.

### Why it transfers to the log family and not here

The peel only attenuates roundings that happen *inside* the polynomial. It
cannot touch the two that form `v*P` and add it to `v*HI`, and in this
function those are the big ones: both land at the result's own scale by
construction, where `ln`'s and `log10`'s equivalents land at `|s| <= 0.415`
and `|s*LOG10_E| <= 0.18`. On top of that the Estrin grouping *already*
attenuates the large partial sums -- the `t8`/`hi` branch reaches 156 and
-198 at `v = 0.5`, but enters through `v^4` twice, so it arrives at 0.15
ulp. What is left inside the poly for a peel to demote is the `lo` group
alone, and 0.26 ulp is exactly what it gets.

**The screening question is not "does the poly's evaluation dominate" but
"how much of that evaluation happens before the last full-weight
rounding".** `ln`'s poly *was* the whole mantissa term with every rounding
at full weight; `erfcx_pos`'s is a correction that a product and a sum then
round again regardless.

### Caveat on the grid

This simulation reports 4.98 where the real exhaustive `erfcx` max is 6, so
the grid misses the true worst points and the absolute numbers are not
comparable to the harness's. The *relative* comparison between groupings is
what it is for, and 0.8% is far enough inside the noise of any grid choice
that a real-chain measurement would not change the conclusion.
## `logaddexp_accurate`: the cancellation the other two tiers accept, in f64

`logaddexp`'s `~1e3-1e5` heavy-tailed max ulp was the worst number in the
crate that is not documented as by-design. It is real: `ln(e^a+e^b) =
m + log1p(exp(-d))` with `m = max(a,b)`, `d = |a-b|`, and the correction
confined to `(0, ln2]`, so whenever `m` lands in `[-ln2, 0)` the two very
nearly annihilate and the answer is set by the correction's **absolute**
error. An f32 correction carries `2^-24` of it however well it is fitted.

The prior entry ("logaddexp: precision budget for the ~1e4-ulp
cancellation", analysis only) derived that any real fix needs a `2^-40`-ish
`exp` and sketched a double-f32 construction with a `2^(j/16)` table. That
sketch is **not what shipped**, and the reason is worth recording: this
crate already runs whole chains in f64 (`log2_f64`, `exp2_f64_to_f32`,
`compound_accurate`), and f64 is both *more* accurate than double-f32
(`2^-53` against `2^-47` at this scale) and cheaper to write correctly. So
`logaddexp_accurate` is simply the whole thing in f64, narrowed to f32
once at the end.

### What it is

```
log1p_exp_neg_f64(d) = 2*atanh(s),   s = 1/(1 + 2*e^d)
```

Two things worth keeping:

- **`exp(+d)`, not `exp(-d)`.** The atanh form of `log1p` needs `E = e^-d`
  only through `s = E/(2+E) = 1/(1+2e^d)`, which is an `fma` and a
  reciprocal against a multiply, an add and a divide -- one operation
  shorter, identically conditioned, and `e^128 = 3.9e55` is nowhere near
  f64's range, so growing the exponential costs nothing.
- **One polynomial covers `d` in `[0, 128]`.** `s` runs over `(0, 1/3]`,
  `2s` is the leading term at *both* ends (`ln2` at `d = 0`, `e^-d` as `d`
  grows), and nothing cancels in between. 15 tail coefficients, `2s`
  pinned outside them.

`d` is formed in f64 too, and that is not cosmetic: `fl32(a-b)`'s rounding
enters the answer amplified by `dcorr/dd = -E/(1+E) ~ 0.5` at small `d`,
which is the *second* item in the old entry's ranked budget. Widening the
subtraction removes it instead of shrinking it.

### The clamp that is not a domain restriction

`u = (s*s).max(1e-30)`. `s` reaches `1.3e-56`, so the Estrin grouping's
`u^4`/`u^8` leave f64's normal range from `d ~ 44` -- inside the useful
domain -- and would drag a denormal assist through the whole vector. The
floor sits far below where the tail matters: `u <= 1e-30` makes `u*Q(u)` a
relative `3e-31` of the pinned `2s`, i.e. `2^-101`. Worth checking on any
Estrin poly whose argument is allowed to get genuinely small; the shape
that bites is `u2*u2*u2*u2`, not `u` itself.

### Measured

`accuracy.rs`, 9.92M random finite pairs: **avg 0.0000, max 0** -- but the
harness reference is `m + log1p_u10(exp_u10(-d))`, itself an f64 chain of
the same `~2e-16` absolute error, so that row proves the ordinary domain
is at the correct-rounding floor and *cannot* score the cancellation
region. For that, a corpus of 20000 pairs built **on** the zero curve
(`a` uniform in `(-ln2, 0)`, `b = ln(1 - e^a)` jittered +-4 ulp) was scored
against an 80-digit Python `decimal` oracle:

| | avg ulp | max ulp |
|---|---|---|
| `logaddexp_accurate` | **0.2675** | **72.2** |
| `logaddexp` | 16917850 | 34372019720 |
| `logaddexp_checked` | 16917850 | 34372019720 |

The max-72 sample has `|m| = 0.46` and a true result of `1.065e-11`; the
absolute error there is `6.3e-17`, i.e. **1.1x half an ulp of f64 at
`|m|`**. There is no f64 chain that does better, so the remaining tail is
the format's, not the formula's -- stated in the doc comment as a bound
(within an ulp while `|result| > ~5e-9`, degrading in proportion below)
rather than a claim of correct rounding.

Three of those oracle-derived pins are in `edgecheck.rs`, because the
blind 2-arg fuzz reaches the zero curve only by luck and the f64 reference
could not adjudicate them anyway.

### Price

| region | instrs | uOps | BlockRT | throughput | latency |
|---|---|---|---|---|---|
| `logaddexp` | 100 | 110 | 28.00 | 2.449 | 74.110 |
| `logaddexp_checked` | 103 | 116 | 30.00 | 2.506 | 73.345 |
| `logaddexp_accurate` | 152 | 214 | 94.00 | **7.616** | **121.517** |

3.1x throughput, 1.6x latency, all packed `pd` (`codegen_check` clean).
The `BlockRT` 28 -> 94 is mostly the half vector width, not the extra
work: 152 instructions is 1.5x, not 3.4x. All three tiers stay --
`logaddexp` on throughput, `logaddexp_checked` on the denormal tail at
`logaddexp`'s price, this on accuracy.

### Two reformulations re-confirmed dead, from the other side

The old entry rejected `log1p(expm1(a) + exp(b))` and friends on paper.
Building the f64 tier confirms the diagnosis constructively: what fixed it
was *only* the wider format, and every f32 reformulation reduces to making
two O(0.5) quantities cancel to O(1e-8) while each carries `2^-24`. Do not
try another algebraic rearrangement of this function -- the requirement is
on the working format, and nothing else.
## `acos`: screened after the `asin` crossover work, no lever, and one whole class of old rejection is now stale

Screened because `asin` and `asinpi` had just gone 5 -> 2 by finding their
error concentrated in one window. **`acos` has no such window.** Strided
scan (every 64th f32, so the avg column is not comparable to the
exhaustive figures elsewhere -- shape only):

```
+[0.0,0.1)   max 3      -[0.0,0.1)   max 2
+[0.1,0.27)  max 3      -[0.1,0.27)  max 3
+[0.27,0.5)  max 3      -[0.27,0.5)  max 2
+[0.5,0.7)   max 3      -[0.5,0.7)   max 2
+[0.7,0.9)   max 3      -[0.7,0.9)   max 1
+[0.9,0.99)  max 3      -[0.9,0.99)  max 1
+[0.99,1.0)  max 3      -[0.99,1.0)  max 1
```

Flat across the whole positive half. `acos` is a single expression --
`sqrt(1-a)*P(a)`, plus `+ PI` for negative `x` -- with no crossover to
move and no subtraction that cancels (`acos(|x|)` never exceeds `pi/2`,
so `PI - acos(|x|)` is a 2x amplification at worst, and the negative half
measures *better* than the positive one, not worse). The `asin`/`asinpi`
lever does not apply here, and nothing else in the shape suggests one.

**What is genuinely stale, and is left open rather than taken.** Six or
more `acos_poly` rejections in this file are of the form "improved acos
but regressed `asin` max ulp 9 -> 12" -- the degree 6 -> 7 bump, the Df32
pi/2 leading-term split, the joint LP with both maxes capped, the
unconstrained joint objective. **That coupling no longer exists.**
`acos_poly` has exactly one call site (`acos`), `asin` has carried its own
`asin_poly` since the decoupling, and after the crossover work that poly
is degree 5 on `[0.5, 1)` and shares nothing. So the constraint that
killed those attempts is gone and each is re-openable against `acos` and
`acosd` alone.

Not taken here, on price rather than on principle: the best-measured of
them (degree 6 -> 7) was worth **acos avg/max 0.496/4 -> 0.490/3** for
**+7-11% mca**, and that is a worse trade than the two this session did
take (5 -> 2 for +6.8% and +8.0%). Anyone re-opening it should re-measure
the mca side first, since it will now land on `acos`/`acosd` only, where
the old figure was quoted with `asin` also paying.

## readme precision-table audit: three stale rows, and the one that keeps going stale

A full quick sweep reconciled row by row against readme's precision tables
(52 rows matched by exact name; the four `sin_fast`/`cos_fast`/`sin_checked`/
`cos_checked` "hits" a looser matcher produced were artifacts of the two
different bands sharing a base name, and those rows are correct).

**Corrected:**

- `log2p1` **0.102 / 3 -> 0.092 / 2**. Verified exhaustively over all 2^32
  patterns, not from the quick run. Almost certainly went stale when
  `log_2_normal`'s peel shipped; nobody re-recorded the `p1` sibling.
- `atan2` and `atan2_unchecked` **max 3 -> 4**. Not a regression, and this
  is the interesting one: **both reached 4 in ordinary 10M-sample runs**,
  `atan2` on the second of four repeats and `atan2_unchecked` on the first.
  `atan2_pos` stayed at 3 across all four, and `atan2pi` (published as 4)
  came back 3 on one of them. Every one of those is the same phenomenon.

**The 2-arg rows cannot be recorded from one run.** `atan2`/`hypot`/`powf`/
`remainder` have no exhaustive mode, so their published max is whatever the
last sampled run happened to find, and a row recorded from a lucky run
understates the function indefinitely -- there is no later check that can
catch it, because re-running is exactly as lucky. Any future edit to these
rows should be the max over >=3 repeats, and should say so.

**Resolved after the fact, in the commit that follows this entry:**
`cos_checked (|x|<=1e6)` read 0.081 avg where the quick sweep said 0.0759 -- 50x the documented ±0.0001 quick-avg
noise floor, and its three band siblings (`sin_fast` 0.0357 vs 0.036,
`sin_checked` 0.0355 vs 0.036, `cos_fast` 0.0779 vs 0.078) all reproduce
their readme numbers to the published precision, so the quick/thorough avg
agree for this band shape and 0.0759 is probably right. Left alone because
the row was presumably recorded `thorough`, and `accuracy thorough
cos_checked` does not select it (the gate is `run("cos")`, and
`n.contains(filter)` needs the *gate* name, not the row label -- worth
knowing before scripting a filtered run). Run as `accuracy thorough cos`
it comes back **0.0760 / max 2** exhaustively, so that row was stale on
*both* columns and is now 0.076 / 2. `cos_fast (|x|<2^22*pi)`'s avg went
0.291 -> 0.288 in the same run; its 2780 max reproduced exactly, as did
`cos (in-domain)` 0.0833/2 and `cos_fast (|x|<=1e6)` 0.0779/3.
knowing before scripting a filtered run). Its max is untouched at 3: quick
mode found 2, and quick's max is a lower bound, so it proves nothing.
## `sin_wide`/`cos_wide`: the pi reduction that has no magnitude limit, and what the gather really costs

2026-08-02. **The worst number in this crate was `sin_checked`/`cos_checked`
over their own advertised domain.** `accuracy quick` reports

```
sin_checked (all f32)     avg ulp 314265980.5   max ulp 2130706432
cos_checked (all f32)     avg ulp 319625290.3   max ulp 2130706432
```

and `2130706432` is not a coincidence -- it is `2 * 0x3f800000`, the ulp
distance from `+1.0` to `-1.0`, i.e. the *worst a bounded output can be*.
Past `2^51*pi` those two return a value in `[-1, 1]` with no relationship
to `sin(x)`, and the `.clamp(-1.0, 1.0)` that keeps the row finite is
exactly what stops it being obvious.

`sin_checked`'s own doc comment said of that region "**nothing can** [fix
the accuracy] -- an f32's own ulp there is ~5e8 radians", and graveyard.md
said the same thing more precisely: "cannot be improved without a
table-indexed Payne-Hanek (bits of `1/pi` selected by `x`'s exponent) ...
and a table lookup is a gather, **which is the one thing this crate's
auto-vectorized scalar style cannot absorb**." Both halves are wrong. The
first confuses "the input's neighbours are far apart" with "the answer for
*this* input is unknowable" -- `x` is an exact dyadic rational and
`sin(x)` is a well-defined number, which is why glibc gets it right. The
second was never measured.

**Shipped: `sin_wide`/`cos_wide`.** `accuracy thorough`, exhaustive over
every f32 bit pattern:

```
sin_wide (all f32)        avg ulp 0.1269   max ulp 2
cos_wide (all f32)        avg ulp 0.1501   max ulp 2
```

### The identity that makes it a fixed-size problem

An f32 is `m * 2^E` for a **24-bit integer** `m`, so

```text
x/pi mod 2 = (m * beta(e)) mod 2,     beta(e) = (2^E / pi) mod 2
```

because `m * floor(2^E/pi)` is an integer whose parity `beta`'s own bit 0
already carries. `x/pi` -- the only unbounded quantity in the problem --
is never formed. `beta(e)` is tabulated per raw biased exponent (256 rows,
denormal row included) as three `f64` words of **29 / 53 / 53
significant** bits, which names it to a relative `2^-135` at every
magnitude; `m * beta` then needs 24 + 135 bits of which the low ~110 are
kept. Fixed budget, no dependence on `|x|`.

Three details that are the whole correctness argument:

- **29 bits in word 0, not 53.** `m` is 24 bits, so `m * a0` is exact at
  24 + 29 = 53, which makes `p0 - round(p0)` an exact peel of the integer
  part with no error-free transform. A 53-bit word 0 would need a
  `two_prod` there.
- **Significant bits, not fixed bit positions.** A fixed-position split
  (bits 0..-28, -29..-81, ...) is cheaper downstream -- `f0 + p1` becomes
  exact and the `two_sum` disappears -- but every word is zero once
  `beta` is small, i.e. for every `|x|` below ~0.125, which then needs a
  bypass branch. Significant-bit splitting is uniform over all 256 rows
  including denormals and needs no bypass at all.
- **The `2^(150-e)` scale is folded into every word**, so the chain
  multiplies `|x|` itself and never reconstructs `m`. An exact
  power-of-two pre-scale, and it removes an exponent-field build, a widen
  and a multiply from the hot path. It also makes the denormal row fall
  out for free: a denormal's absent implicit bit just makes `|x|` a
  23-bit integer instead of 24.

The parity is `n0 + n1`, two small exact integers, so `ROUND_MAGIC64`'s
`2^51` window -- the thing `reduce_pi64` actually falls off -- cannot
apply. Verified against exact rationals (900-bit Machin `pi`, `Fraction`
arithmetic): **0 parity mismatches and worst relative error `2^-53.00` on
the reduced fraction**, over 120000 random `(e, m)` pairs spanning every
exponent from the denormal row to `e = 254`. `2^-53` is the final
`fc + ec` add, i.e. one rounding of the answer itself -- nothing else in
the reduction is inexact.

### The gather: 3.2x throughput, 11% latency

| region | `_checked` | `_wide` | delta |
|---|---|---|---|
| `sin` throughput | 2.495 | **8.045** | +222% |
| `cos` throughput | 3.157 | **9.299** | +195% |
| `sin` latency | 82.00 | **91.05** | +11.0% |
| `cos` latency | 87.00 | **99.06** | +13.9% |

Latency is nearly free -- three gathers plus a wider chain add 9 cycles to
one that is already 82 deep, because the gathers issue early and the table
words are not on the critical path. (An earlier draft of this entry
claimed **+1.3%**, from measuring a clamp-less `sin_wide` against a
clamped `sin_checked`: the `.clamp` is two `vminps`/`vmaxps` at the very
end of a serial chain and is worth 8 cycles on its own. Compare tiers with
the same tail.)

Throughput is not free, and **the gathers themselves are not the
reason**, which is the transferable part:

```
llvm-mca Block RThroughput, one instruction:
  vgatherdpd (%rax,%xmm2,8), %ymm0 {%k1}   2.0   (4 doubles)
  vgatherdpd (%rax,%ymm2,8), %zmm0 {%k1}   4.0   (8 doubles)
  vpgatherdd/vpgatherqd                    4.0
  vpermi2q                                 1.0
  vfmadd231pd %zmm                         1.0
```

At 12 gathers per 16 elements that is 24 of the region's 66 Block
RThroughput. The other 42 is that **LLVM drops the vectorization factor
from 8 to 4 the moment an `f64` gather appears** -- `sin_checked`
vectorizes f32 in `ymm` and f64 in `zmm`, `sin_wide` falls back to `xmm`
/`ymm`, so *every arithmetic op in the whole function* costs twice as much
per element. Three measurements isolate it:

| variant | VF | cyc/elem |
|---|---|---|
| `sin_checked` (baseline) | 8 | 2.495 |
| `sin_wide`, table index replaced by a constant (no gather) | 8 | **3.286** |
| `sin_wide`, `[u32; 256]` probe tables (32-bit gather) | 8 | **5.022** |
| `sin_wide` as shipped, `[f64; 256]` tables | **4** | 7.540 |

(the three probe rows are clamp-less, so they are comparable to each other
and to `sin_wide`'s own 7.540 at that time, not to the shipped 8.045.)

So the reduction's *arithmetic* is only +32% over the baseline, a 32-bit
gather keeps VF = 8, and the shipped 3.0x is mostly a codegen cliff rather
than the algorithm. **Reducing the gather count does not move it** -- two
`f64` gathers instead of three measured 7.539, identical to three.

Left open with a measured target rather than taken, because it is a real
rewrite: `[u32; 256]` tables of 29-bit chunks at *fixed* bit positions
(the thing the shipped version deliberately avoids) would keep VF = 8, but
they need the `|x| < 0.25` bypass back, a `2^(150-e)` scale rebuild, and
four gathers instead of three to reach the same precision. The probe above
puts the ceiling at ~5.5 cyc/elem, i.e. **8.05 -> ~5.5**, against a
baseline of 2.495. Anyone taking it should re-measure the probe first.

### Kept as a separate tier, not a replacement

Neither dominates: `sin_checked` is 3x cheaper in throughput and correct
for every `|x|` a bounded caller will ever produce; `sin_wide` is correct
everywhere. `tan_wide` is *not* shipped -- `tan_checked` calls the
reduction twice, so it would pay the gathers twice, and no measurement
was made.

### The `.clamp(-1, 1)` survives an exact reduction, and the reason is not the reduction

Worth recording because the argument for dropping it is completely
convincing and completely wrong. `sin_checked`/`cos_checked` clamp because
their reduction can hand `sinf_poly` a residual of ~1000, and a degree-9
poly evaluated there is ~2.6e21 -- a real `|sin(x)| <= 1` violation for
ordinary finite input. `sin_wide`'s reduction is exact, so `|r| <= pi/2`
by construction and that mechanism cannot fire. The clamp was dropped on
exactly that reasoning.

An exhaustive `|result| <= 1` scan over all 2^32 patterns then found
**2726588 violations, worst `1.0000001`**: `sinf_poly` is a minimax fit of
`sin` near its own maximum and has no reason to stay under it. Every
pinned `edgecheck` value passed, `special_matrix` passed, and the ulp rows
cannot see it (1 ulp above 1.0 is a 1-ulp error, not an outlier). Two
transferable points: **a proof that one mechanism cannot fire is not a
proof that the invariant holds**, and the cheap check that catches this
class -- an exhaustive scan of the *invariant* with no reference function
to compute -- runs in a couple of minutes over the whole domain, which is
far cheaper than an exhaustive accuracy sweep and finds a different kind
of bug.

The clamp is also free accuracy: the true `|sin|` never exceeds 1, so the
correctly-rounded f32 answer never does either, and clamping can only move
a result towards it.

### Side finding, not fixed: `sin`/`cos`/`sin_fast`/`cos_fast` also exceed 1

The same exhaustive `|result| <= 1` scan, run across every trig tier
*in each one's own documented domain* (`|x| < 2^22*pi`):

| | patterns with `|result| > 1` | worst |
|---|---|---|
| `sin` | 660 | 1.0000001 |
| `cos` | 2720382 | 1.0000001 |
| `sin_fast` | 670 | 1.0000001 |
| `cos_fast` | 2720361 | 1.0000001 |
| `sin_checked`/`cos_checked`/`sin_wide`/`cos_wide` | 0 | 1.0 |

Same `sinf_poly` overshoot; the four `_checked`/`_wide` tiers are clean
only because they clamp. Left alone deliberately -- those four are the
fast tiers, `1.0000001` is a 1-ulp error and already inside their
published max, and a clamp is two ops on functions costing 1.15-1.78
cyc/elem in total. Recorded because "does this function respect
`|sin| <= 1`" is a different question from its ulp row, and nothing in the
harness asked it before. (Outside their domains all four return `inf`,
which is documented and is what `accuracy`'s NON-FINITE annotation
counts.)

## `erfcx_pos`: the reciprocal is an identity away from exact, and the record's own diagnosis was the thing to check

Three entries above (the reciprocal-offset scan, the Horner reorder, the
leading-term peel) all close with the same sentence: `erfc`/`erfcx`/
`norm_cdf`/`gelu` are `erfcx_pos`'s evaluation and need "a structurally
cheaper representation of `erfcx` on `[0, inf)`". That handover was
followed and it is wrong in a useful way -- there was a cheap structural
fix left, it is two instructions, and having taken it the diagnosis is now
genuinely different: **none of these four functions is limited by
`erfcx_pos` any more.**

### What shipped, and why it is not the rejected compensation

`v = 1/(2 + xa)` satisfies

    v == 0.5 - (xa/2)*v

identically (put the right side over `2*(2+xa)` and the `xa` cancels).
Substituting the *computed* `v0 = fl(1/fl(2+xa))` back into it gives a `v`
whose relative error is `(xa/2)` times `v0`'s, plus the one rounding the
`fma` itself makes. `xa/2` is exact, and it is the **unrounded** `xa` that
enters -- which is the whole point, because `fl(2+xa)` throws `xa`'s low
bits away outright and `erfcx` has a nonzero slope at 0 while `v` is
stationary in relative terms there, so that rounding arrives amplified by
`2*|d(ln erfcx)/d(ln v)| = 2.257` at `xa = 0`.

    let v0 = 1.0 / (2.0 + xa);
    let v = if xa <= 2.0 { fma(-0.5 * xa, v0, 0.5) } else { v0 };

The record already rejected "compensating `v`" at "4-6 more ops, correct
ordering included, for ~1.5 ulp". That rejection was of an **EFT**: recover
the exact residual of `2+xa` by Fast2Sum and Newton it back in. The
identity is a different mechanism, costs two arithmetic instructions
instead of four to six, and lands further -- it reaches the exact-`v`
floor, not part of the way.

- **The attenuation is `xa/2`, so above `xa = 2` the identity amplifies
  `v0`'s error instead** -- hence the select, which is the only reason it
  is not unconditional. Threshold scanned 0.75..2.5: everything in
  `[1.0, 2.5]` gives the same max, `0.75` is worse; `2.0` is shipped as the
  identity's own crossover. **Clamping the multiplier (`min(xa, 2.0)`)
  instead of selecting is not the same function** -- it breaks the identity
  and returns 9.5e8 ulp above the clamp. Measured, so it is not
  re-proposed.
- `erfcx_pos(0.0) == 1.0` is preserved bit-exactly and for free: at
  `xa = 0` the fma is `fma(-0.0, 0.5, 0.5) = 0.5`, the same `v` as before,
  so the pin the shipped coefficients are holding is untouched and no
  refit was needed.
- numpy simulation of the exact f32 instruction sequence, 6M distinct f32
  points over `[1e-7, 30]`, shipped coefficients: **max 6.223 -> 4.187,
  avg 0.9997 -> 0.6897**. The exact-`v` floor for this evaluation order
  and these coefficients is **4.187 / 0.6650** -- i.e. the identity
  recovers *all* of the available max and 93% of the available avg. That
  is the sense in which `erfcx_pos` is now purely its own polynomial
  evaluation.

### Two `fma` folds found while pricing it, both free

Both were taken and both are worth more than they look, because they pay
for the identity's select:

- **`erfc`/`norm_cdf`'s sign fold, one level up.** `y = e * t;
  mulsign(y, x) + w` became `fma(mulsign(e, x), t, w)`. Carrying the sign
  on the Gaussian factor rather than on the finished product lets the
  closing multiply and the reflection's add fuse: the `x >= 0` arm is
  bit-identical (`fma(e, t, +0.0)` is exactly `e*t`, confirmed -- the
  positive-side scan is unchanged to four decimals) and the `x < 0` arm
  rounds `2 - e*t` **once instead of twice**. `erfc`'s negative side
  2.87 -> 2.61 max, 0.216 -> 0.193 avg, for **one fewer instruction**.
- **`erfcx`'s negative arm.** `g = exp*(2+2pe); g - r` became
  `fma(exp, 2+2pe, -r)`: one rounding and one instruction cheaper. The
  `inf` argument the old comment gives for keeping `2+2*pe` as the
  multiplicand still holds (`inf*(positive) - finite = +inf`), and
  `erfcx(-inf) = +inf` is still pinned.

### Measured

Exhaustive, every f32 bit pattern, before/after on the same harness:

| row | before | after |
|---|---|---|
| `erf` (control, no `erfcx_pos`) | 0.0270 / 3 | 0.0270 / 3 |
| `erfc` | 0.1948 / 7 | **0.1289** / 7 |
| `erfcx` | 0.2080 / 6 | **0.1439** / 6 |
| `erfcx (\|x\|<=20)` | 0.2078 / 6 | **0.1442** / 6 |
| `erfcx (x>=20)` | 0.2683 / 2 | 0.2683 / 2 (bit-identical, `xa > 2`) |
| `norm_cdf` | 0.0996 / 8 | **0.0657 / 7** |
| `norm_pdf` | 0.0269 / 4 | 0.0269 / 4 (unchanged -- no `erfcx_pos`) |
| `gelu` | 0.2029 / 9 | **0.1433** / 9 |

Range-restricted probe scans, both binaries, same ranges, for where the
change actually bites: `erfc` over `|x|` in `[1e-7, 0.03125]` goes
6.540 / 1.2522 -> **4.720 / 0.7165** on the positive side and
3.513 / 0.6007 -> **2.678 / 0.3532** on the negative; over
`[0.03125, 10]` it is 7.001 / 1.0556 -> 6.790 / 0.8891 and
3.641 / 0.2923 -> 2.611 / 0.1925.

llvm-mca, all eight regions, against the same baseline the entries above
quote (62.079/2.885 and 66.986/2.827 reproduced to the digit):

| region | instrs | uOps | BlockRT | cyc/unit |
|---|---|---|---|---|
| erfc_latency | 4551 -> 4766 | 4721 -> 4942 | 1184 -> 1216 | 62.08 -> **61.42 (-1.1%)** |
| erfc_throughput | 131 -> 136 | 148 -> 157 | 40 -> 41 | 2.885 -> 3.069 (+6.4%) |
| erfcx_latency | 4129 -> 4258 | 4869 -> 5125 | 1184 -> 1216 | 66.99 -> **64.03 (-4.4%)** |
| erfcx_throughput | 125 -> 130 | 143 -> 148 | 39 -> 40 | 2.827 -> 2.903 (+2.7%) |
| norm_cdf_latency | 4980 -> 5210 | 5169 -> 5468 | 1280 -> 1312 | 70.22 -> **68.99 (-1.8%)** |
| norm_cdf_throughput | 140 -> 143 | 160 -> 172 | 43 -> 44 | 3.219 -> 3.436 (+6.7%) |
| gelu_latency | 5773 -> 6034 | 6157 -> 6419 | 1504 -> 1536 | 75.56 -> **74.02 (-2.0%)** |
| gelu_throughput | 156 -> 161 | 184 -> 193 | 50 -> 51 | 3.574 -> 3.792 (+6.1%) |

Full asm region diff: **299 of 307 regions byte-identical**, the 8 changed
ones exactly these four functions -- `exp_checked` and its other 20+
callers (`log1p`, `sigmoid_grad`, `tanh_grad`, ...) did not move.

Latency down on all four, throughput up 2.7-6.7% (`BlockRT` says +2.5%,
which is the usual direction of this crate's mca throughput artifact, so
the real number is probably nearer the structural one -- not arbitrated
with a wall-clock A/B, the machine was at load 17-22 all session).

### The diagnosis this replaces, which is the part worth keeping

**`erfc`'s max did not move and cannot be moved from inside `erfc`.** Split
its ulp error at each point into the `erfcx_pos` factor and everything the
Gaussian path contributes (probe: score `erfc(x)` against
`exp(-x^2)_f64 * erfcx(x)_f32`, i.e. take the *actual* f32 `erfcx` factor
as given), exhaustively over `x` in `[1, 4]` where the worst case lives:

    total          max 6.79
    erfcx factor   max 4.19
    gaussian+tail  max 4.19

**Neither half dominates.** The Gaussian half is `exp`'s own accuracy
(`exp` is 3 max ulp, and `exp_r_poly` is core-locked and recorded
exhausted), so roughly half of `erfc`'s remaining max is not reachable
from inside this function at all, and the other half is the polynomial
evaluation the three entries above already measured as immovable.

**`erfcx`'s worst case has moved to the negative arm and is also `exp`.**
Same split on `erfcx(-a) = 2e^(a^2) - erfcx(a)`, exhaustive over
`a` in `[0.25, 1]`:

    total          max 5.90 at x = -0.514
    erfcx_pos part max 1.74
    gaussian part  max 5.23

which is `exp_reduce!(0.264)`'s own ~2.6 ulp amplified by
`2e^(x^2)/erfcx(x) = 1.30`. The positive arm is now 4.03 max. So the
standing "`erfcx`'s worst case is `erfcx_pos` near zero" note in its doc
comment was true and is not any more.

### Also screened and not taken

- **A two-branch split** -- the specific proposal this session inherited:
  `erfcx(x) = fma(x, Q(x), 1.0)` below a seam (argument exact, no `2+xa`
  rounding at all) and `u = 1/x` above it. Idealised LP-minimax relative
  fit error, both branches, so the degree bill is on the record:

  | seam T | small branch `1 + x*Q(x)` | tail branch `u*(HI + R(u))` |
  |---|---|---|
  | 1.0 | deg 7: 0.37 ulp | deg 10: 2.9, deg 11: ~0 |
  | 1.25 | deg 8: 0.23 | deg 9: 2.6, deg 10: 0.25 |
  | 1.5 | deg 9: 0.13 | deg 9: 0.39, deg 10: ~0 |
  | 2.0 | deg 10: 0.31 | deg 8: 0.29 |
  | 2.5 | deg 11: 0.44 | deg 7: 0.60 |

  The cheapest balanced pair is ~18 degrees against the shipped 10, and
  **both branches are evaluated on every lane** (the crate is branchless
  by requirement), so it is ~+12 instructions -- roughly double the
  polynomial work. Against that: the shipped chain restricted to `xa >= 1`
  already scores 3.54 max and `xa >= 0.5` scores 4.79, so the *most* a
  split can buy on the positive arm is ~6.2 -> ~3.5, and the sections
  above show the public functions would not see it because `exp` binds
  first. Not taken.
- **Reciprocal offsets below 2.** The offset entry above scanned `a` = 2,
  2.5, 3, 3.5 through the real chain and concluded raising it loses.
  Lowering it also loses, and for the opposite reason: the `a+xa` rounding
  is amplified by `1.128*a`, so `a = 1` halves that term -- but `v` then
  spans `(0, 1]` instead of `(0, 1/2]` and the degree-10 fit degrades from
  0.003 to 5.1 ulp-equivalent. Real chain, fresh LP fit at each offset,
  5M f32 points: `a` = 0.5 / 0.75 / 1.0 / 1.25 / 1.5 / 1.75 / 2.0 / 2.5 /
  3.0 gives max 575 / 73.2 / 23.6 / 8.30 / 7.18 / 5.69 / 6.05 / 5.80 /
  7.10. Nothing outside the noise band the offset entry already
  documented. **The offset lever is closed in both directions.**

### Transferable

- **An algebraic identity on the reciprocal beats an EFT on its
  denominator.** `1/(a+x) == 1/a - (x/a)*(1/a+x)` re-substitutes the
  *unrounded* `x`, so it does not need the residual of `a+x` at all. Two
  instructions, attenuation `x/a`, valid wherever `|x| < a`.
- **The one other site of that shape in the crate is `sigmoid`, and it is
  not worth it** -- screened here so the next instance does not spend a
  session on it. `w = 1.0/(1.0 + e)` obeys `w == 1 - e*w`, attenuation
  `e`, so it helps exactly where `e < 1`, i.e. `x > 0`. But `sigmoid`'s
  error does not live there: exhaustive per-octave, max is **1.71 for
  `x > 0`** against **3.35 for `x < 0`**, and scoring the shipped result
  against an exact reciprocal of the *same* f32 `e` puts at most 1.5 ulp
  of that on the reciprocal in either direction. The rest is `exp`, on
  the arm the identity cannot reach.
- **A rejection is of a *mechanism*, not of a goal.** "Compensating `v`
  was measured and rejected on cost" reads like the goal is closed. It
  closed one route to it.
- **Split a composite's error before concluding the primitive is the
  binding term.** Three entries in a row named `erfcx_pos` as the limit
  for four public functions; scoring the composite against
  "exact-Gaussian, actual-f32-`erfcx`" takes ten lines and shows the two
  halves are equal, which reframes what is left to do.

## `dawson`: the concentrated-window screen says no, and the `u = fl(x*x)` term is not worth its ops

Screened after the `erfcx_pos` identity entry above, because that entry's
lever (kill the argument's own rounding with one fma) looks like it should
transfer to `u = x^2` -- the decomposition in `dawson`'s harness entry
prices `u = fl(x*x)` at ~0.8 of the 5.7 max. Both checks come out
negative; recording so the next instance skips them.

**The error is not concentrated, so the `asin`/`atanh` crossover lever does
not apply.** Exhaustive at stride 8 over `[1e-3, 4]` (12.6M points),
current code, series reference:

| binade | max ulp | #(err >= 4.5) |
|---|---|---|
| 2^-10 .. 2^-5 | 1.51 -> 3.55 | 0 |
| 2^-4 | 3.78 | 0 |
| 2^-3 | 4.44 | 0 |
| 2^-2 | 4.83 | 4 |
| 2^-1 | 4.72 | 4 |
| 2^0 | 5.50 | 355 |
| 2^1 | 5.54 | 238 |

It climbs monotonically with `|x|` and the worst set is spread across the
whole of `[1, 4]` -- not a pocket a seam move or a dedicated sub-branch
can excise. (`asin`'s entire >=3-ulp set sat inside `[0.27, 0.4997]`;
that is what a concentration looks like.)

**Compensating `u` is structurally expensive here, unlike `erfcx_pos`'s
`v`.** The identity that worked for `1/(2+x)` re-substitutes the exact
argument into a *closed form* of the variable. `u = x^2` has no such
identity: the exact residual is one fma (`ue = fma(x,x,-u)`), but applying
it to `x*P(u)/Q(u)` needs `d(ln P/Q)/du`, i.e. `P'` and `Q'` -- two more
polynomial evaluations -- because the correction has to enter *through*
the rational, not through its argument. The ceiling is also low: the
harness entry's own f64 rows say an exact `x^2` takes 4.322 -> 3.517, so
the whole term is worth ~0.8 ulp of a 5.7, for 5+ ops.

**The rule this sharpens:** an argument-rounding compensation is cheap
only when the variable has an algebraic identity that re-admits the
unrounded input (`1/(a+x) == 1/a - (x/a)*(1/(a+x))`). When the argument
enters an opaque approximant, compensating it costs a derivative, and the
EFT that recovers the residual is the cheap half of the job.

## Method: an out-of-tree probe crate is what makes "where does the error live" cheap

Both `erfcx_pos` entries above turn on measurements `accuracy.rs` cannot
express -- error restricted to one binade, error restricted to one *sign*,
and error split between a composite's factors -- and on A/B-ing two builds
of the same function. Adding rows to `accuracy.rs` for that is shared-file
churn, and adding scratch files to `examples/` has already put one by
accident into a commit. The alternative used here costs about a minute to
set up and is worth writing down:

    <scratch>/probe/Cargo.toml     jodiemath-rs = { path = "<worktree>" }
                                   sleef = "0.3.3"        # same f64 refs
    <scratch>/probe/.cargo/config.toml
                                   [build] rustflags = ["-C","target-cpu=native"]
    <scratch>/probe/rust-toolchain.toml   copied from the worktree

`src/bin/*.rs` then gives one binary per question, and `cargo build
--release` rebuilds only the probe (~2s) after a `src/lib.rs` edit.

Three things it bought that are hard otherwise:

- **A/B on two builds.** `git show HEAD:src/lib.rs > src/lib.rs`, build the
  probe, copy the binary aside, restore, build again. Two binaries scoring
  the *same* ranges with the *same* reference -- the cleanest before/after
  this crate has, and far cheaper than two `thorough` sweeps under load.
  (`cargo build --example accuracy` cannot do this at all while a sweep is
  running: writing a running executable is `ETXTBSY`.)
- **Per-binade and per-sign tables.** `erfcx`'s worst case turning out to
  be on the negative arm, and `erfc`'s error being flat across `[1,4]` and
  falling below `2^-3`, are both invisible in a single avg/max row.
- **Factor decomposition.** Score the composite against
  "exact-in-f64 factor A times the *actual* f32 factor B" to get B's share,
  and the complement for A's. Ten lines, and it is what shows `erfc`'s two
  halves are 4.19 and 4.19.

The `.cargo/config.toml` is not optional: without `-C target-cpu=native`
the crate's own `compile_error!` on missing FMA fires, and per the
RUSTFLAGS entry elsewhere in this file a failed build leaves the previous
artefact in place, so the probe would silently measure stale code.
## `tan_wide`, and the throughput column's own branch artifact

2026-08-02, immediately after `sin_wide`/`cos_wide`. With those landed,
the worst row in the crate was `tan_checked`, for exactly the same reason
(it calls `reduce_pi64` twice) -- `accuracy quick` reads **avg 406062669 /
max 2310525355** over all f32.

`tan_wide` is `tan_checked`'s body verbatim on `reduce_pi_wide`.
`accuracy thorough`, exhaustive over every f32 bit pattern (102s):

```
tan_wide      avg ulp 0.2495   max ulp 4   worst x 1.3138148
```

against `tan_checked`'s quick-fuzz **avg 406004054 / max 2324484283**.
A separate dense stride-4093 scan confirms the change is *only* the
reduction: on the bands where `tan_checked` is still accurate the two
agree to four digits (3.3356 / 3.3356 on `[1e-3, 1e3]`, 3.3136 / 3.3136
on `[1e3, 1e6]`), and diverge only above `1e15` (3.3259 vs 2.199e12).

Worth noting because the existing `tan_checked` row carries the opposite
expectation: its comment warns the number "will look alarming" because
`tan` has a pole every `pi` and at large `|x|` the poles sit closer than
the local float spacing, so any correct implementation shows unbounded
relative error near them. **With an exact reduction that does not
happen** -- max 4, and the worst input is `1.31`, nowhere near a pole.
The unbounded-ulp-near-a-pole argument was describing a real effect of
the *broken* reduction landing on the wrong side of a pole, not an
intrinsic property of `tan`; once the pole locations are right, the
reference and the function are near-pole together and the relative error
stays bounded. This is the same shape as this file's own "near-zero
excuse needs an exact reduction" rule, at infinity instead of zero.

**GVN does share the two reductions.** `tan_wide` calls
`reduce_pi_wide` twice, `HALF = false` and `HALF = true`, and the emitted
region contains **three** `vgatherdpd`, not six: the table lookups, the
three products, the integer peel and the two-word residual are common
subexpressions and only the half-turn shift and the final `* pi` are
duplicated. Same phenomenon as the already-recorded "inline CSE defeats
pair-fn ideas", used deliberately this time. It does *not* make the wide
tier relatively cheaper here than for `sin` -- `tan_checked`'s two
`reduce_pi64` calls share the same way, so both land at ~3.2x. Sharing is
why the ratio is not *worse*, which is a different claim, and the draft of
this entry got it wrong until the corrected numbers below arrived.

### The measurement trap this turned up, which is the more transferable half

`tan_wide_throughput` first measured **3.479 cyc/elem against
`tan_checked`'s 4.289** -- i.e. the wide tier appearing *cheaper* than the
one it fixes. That is this file's own documented tell ("a composite
cheaper than its parts means something got hoisted"), and the cause is
new:

`throughput_fn!` wraps a 16-element loop, and both `mca.rs` and
`tools/mca_region.py` divide `TotalCycles` by `ARR_LEN = 16` to get
cycles per element. That is only right when LLVM **fully unrolls** the
loop. When the body gets big enough that it keeps the loop, llvm-mca
simulates one *iteration*, and the `/16` is wrong by exactly the unroll
factor -- silently, and always in the flattering direction. `tan_wide`'s
region ends `addq $4, %rax / cmpq $16, %rax / jne`, so it covers 4
elements: **the true figure is 13.915, not 3.479.**

readme.md said of this column: "The **throughput** column has no such
problem (those regions are vectorized and if-converted to masked selects,
so there is no branch to mis-simulate)." That is true of *data-dependent*
branches and false of a retained loop, which is a different mechanism.

Audited all 153 throughput regions for a backward branch to one of their
own labels. Exactly two pre-existing cases, and one is a relief:

| region | real step | published | true |
|---|---|---|---|
| `clog_re_throughput` | 8 | 7.744 | **15.489** |
| `powf_throughput` | -- (forward `jne`/`jp` only, no loop) | 6.996 | 6.996 |

`clog_re` is in neither `mca.rs`'s `order` array nor the readme, so
nothing published was ever wrong. `powf` looked like a hit on a first,
sloppier audit that matched any jump to an in-region label without
checking direction -- **check the direction**, its `jne`/`jp` are forward
branches into a shared tail.

Both tools now read the loop's own induction step (`addq $N, %r..`, LLVM's
own unroll factor, not a guess) and divide by that, printing
`(loop kept, /4 not /16)` so the row cannot go quietly wrong again.

### Cost, corrected

| region | `tan_checked` | `tan_wide` | delta |
|---|---|---|---|
| throughput | 4.289 | **13.915** | +224% |
| latency | 101.00 | **116.06** | +14.9% |

Same shape as `sin_wide`'s 3.2x and for the same reason -- the f64 gather
drops the vectorization factor from 8 to 4. `tan_wide`'s region is `xmm`
/`ymm` where `tan_checked`'s is `ymm`/`zmm`.

## `sin`/`cos` exceed 1 by an ulp too, and the clamp is +14% -- rejected, documented instead

2026-08-02, follow-on from `sin_wide`'s own clamp finding. An exhaustive
`|result| <= 1` scan over all 2^32 patterns, run across every trig tier in
this crate *inside* `|x| < 2^22*pi` (`cos`'s documented domain, the
narrower of the two):

| | patterns with `|result| > 1` | worst |
|---|---|---|
| `sin` | 660 | 1.0000001 |
| `cos` | **2720382** | 1.0000001 |
| `sin_fast` | 670 | 1.0000001 |
| `cos_fast` | 2720361 | 1.0000001 |
| `sin_checked` / `cos_checked` / `sin_wide` / `cos_wide` | 0 | 1.0 |

Always exactly one ulp over, never more, and always where the true value
is under 1 by less than an ulp -- so it is *inside* the published 2-ulp
row and invisible to every accuracy number the crate keeps. The four
clean tiers are clean only because they clamp.

**Priced and rejected.** `.clamp(-1.0, 1.0)` on `sin` and `cos`, measured
one at a time by replacing the real function bodies:

| region | before | after | delta |
|---|---|---|---|
| `sin` throughput | 1.778 | 2.033 | **+14.3%** |
| `sin` latency | 64.00 | 72.00 | **+12.5%** |
| `cos` throughput | 1.654 | 1.883 | **+13.8%** |
| `cos` latency | 61.00 | 69.00 | **+13.1%** |

`tan` does not move (2.985 / 78.00 either way -- it goes through
`pi_reduce_and_poly!` directly, not through `sin`/`cos`), and it would not
want the bound anyway.

14% on the crate's cheapest accurate trig tier, to move an error that is
already inside its own published max, is the wrong trade. **Documented in
`sin`'s and `cos`'s doc comments instead**, with the measured cost, so the
next reader does not have to re-derive either half: a caller who needs
`|sin| <= 1` -- feeding `acos`, a `sqrt(1 - s*s)`, or a range assertion --
has four functions that guarantee it.

Two things worth reusing. **An exhaustive scan of an *invariant* costs
nothing** -- there is no reference function to evaluate, so all 2^32
patterns take a couple of minutes against the hours an exhaustive accuracy
sweep needs, and it finds a class of defect ulp rows structurally cannot
see. And **"is the invariant violated" is a separate question from "is the
answer accurate"**: 1 ulp over 1.0 is simultaneously a correct answer and
a broken postcondition.

## `erfinv_tail_poly`: the twelfth coefficient, and where the three callers' error actually lives

The `sqrt(w)-1` variable change landed earlier the same day took
`erfinv`/`erfc_inv`/`probit` from ~15 to ~6 max ulp and left the tail
polynomial at "2.8 idealized ulp, saturating". Re-screening it against a
300k-point dense grid with the same ulp weight and a wider coordinate
descent says it was not saturating -- it was one coefficient short of a
cliff:

| coefficients | 10 (shipped) | 11 | 12 | 13 |
|---|---|---|---|---|
| idealized ulp | 3.32 | 2.38 | **0.78** | 0.84 |
| amplifier `max sum|c_k t^k| / |Q|` | 3.42 | 2.96 | 6.18 | 11.76 |

10 -> 11 buys 28%, 11 -> 12 buys **3.1x**, and 13 measures worse after f32
quantisation. A degree sweep that flattens and then falls off a cliff is
not the usual shape, which is presumably why the first sweep stopped: two
consecutive small steps read as saturation.

Measured end to end (scipy `ndtri`/`erfcinv`/`erfinv` as the reference,
f32 grid strided by 128, so the same points before and after):

| | avg before | avg after | max before | max after |
|---|---|---|---|---|
| `erfinv` | 0.3973 | 0.3925 | 5.0727 | **3.6714** |
| `erfinv`, tail arm only | 1.5990 | **0.5895** | 5.0727 | 3.6714 |
| `erfc_inv` | 0.8680 | 0.7199 | 5.8861 | **4.5009** |
| `probit` | 0.9444 | 0.8003 | 7.0831 | **4.9831** |

llvm-mca, +2 `fma` per call landing on three functions:

| region | instrs | uOps | BlockRT | cyc/unit |
|---|---|---|---|---|
| `erfinv_throughput` | 165 -> 173 | 193 -> 206 | 45 -> 47 | 3.930 -> 4.314 |
| `erfc_inv_throughput` | 214 -> 219 | 275 -> 283 | 59 -> 61 | 6.166 -> 6.473 |
| `probit_throughput` | 221 -> 223 | 283 -> 289 | 61 -> 63 | 7.167 |
| `erfinv_latency` | 5506 -> 5698 | | 1408 -> 1472 | 47.870 -> 46.973 |
| `erfc_inv_latency` | 7173 -> 7365 | | 1824 -> 1888 | 110.845 -> 111.861 |

+2 Block RThroughput on each (+3.3 to +4.4%), as the "degree bump costs
real" rule predicts. Taken because 25-30% of the max on three functions is
worth 4%.

### Where the rest of it is, measured rather than guessed

Same probe, after the bump, split by branch. This is the map the next
attempt should start from:

| | avg | max |
|---|---|---|
| `erfinv` central (gets `x` exactly) | 0.392 | 2.66 |
| `erfc_inv` central (`x = fl(1-n)`) | 0.715 | 3.60 |
| `probit` central (`+ sqrt(2)*`) | 0.737 | **4.98** |
| `erfinv` tail | 0.590 | 3.67 |
| `erfc_inv` tail | 0.720 | 4.50 |
| `probit` tail | 0.801 | 4.59 |

So `probit`'s binding term is no longer the tail at all -- it is the
central arm, and the two steps that get it there are both *chain*, not
fit: forming `x = 1.0 - n` costs 2.66 -> 3.60 (one bit, since `n < 0.5`
puts `x` on a coarser grid than `n`), and the outer `sqrt(2)*` costs
3.60 -> 4.98.

### `erfinv_central_poly` has no headroom at all -- do not refit it

Screened the same way, over `u = x^2` in `[0, 0.49]`, ulp-weighted, LP +
descent:

| coefficients | 9 (shipped) | 9 refit | 11 | 12 |
|---|---|---|---|---|
| idealized ulp | **0.464** | 0.504 | 0.504 | 0.677 |

The shipped nine are already *better* than anything the same machinery
produces at any degree, and the unquantised LP optimum at 11 coefficients
is 0.0017 ulp -- i.e. the fit has been irrelevant here for a long time and
f32 coefficient quantisation is the entire floor. The amplifier is 1.006,
so there is nothing to win by changing variable either. The central arm's
2.66-4.98 is chain rounding, full stop.

## `probit`: folding `sqrt(2)` into the tail's own sqrt (IDEAS #181) -- real, then superseded

`probit(p) = -sqrt(2)*erfc_inv_half(2p)`, and on the tail arm
`sqrt(2)*sqrt(w)*Q` is `sqrt(2w)*Q` with **`2*w` exact**, so the scale can
ride the sqrt the combine already performs: one rounding instead of two
plus `fl(sqrt 2)`'s 0.287-ulp-low bias. The polys still want the unscaled
`v = sqrt(w)`, which becomes `t = fma(vc, 1/sqrt2, -1)` -- the same
instruction count as `v - 1`, with the constant's error landing on a
polynomial argument whose relative sensitivity `|v Q'/Q|` is at most 0.064
(0.052 far arm), so 0.7 ulp of it comes out as 0.045.

Implemented as a `const SQRT2: bool` generic on `erfc_inv_half`, with
`erfinv`/`erfc_inv` instantiating `false` -- their four mca regions verified
**byte-identical** to before, so the tiering costs the other two callers
nothing.

Measured **on the pre-degree-bump code**: `probit` avg 0.9444 -> 0.8825
(-6.6%), max 7.0831 -> 6.7963 (-4.1%), for `probit_throughput` 221 -> 225
instrs, 283 -> 286 uOps, Block RThroughput 61 -> 62 (+1.6%), and
`probit_latency` -3.7%. A real but thin trade.

**Then the tail poly's twelfth coefficient landed and it stopped paying.**
`probit`'s max moved to the central arm (4.98 against the tail's 4.59), so
folding the scale out of the *tail* now buys zero max ulp and only ~5% of
the average, for the same +4 instructions. Not shipped.

It becomes worth revisiting only as a package with the central arm --
`sqrt(2) * x*P(x^2)` folded into a scaled, re-descended copy of
`erfinv_central_poly` (0 extra ops there, 9 duplicated constants), which
together would be 4.98 -> ~4.3. On its own it is dominated.

Third instance of the "[re-stale-check cross-function deps]" pattern, and
the first where a *win* went stale rather than a rejection: a change worth
taking against one baseline was worth nothing against the next one landed
30 minutes later. Measure the lever against the code you are actually
shipping on, not the code the idea was written against.

### The scaled central polynomial, already fitted

So the next attempt does not have to redo it. `sqrt(2)*erfinv(x)/x` over
`u = x^2` in `[0, 0.49]`, ulp-weighted against `probit`'s own grid, LP
start from `sqrt(2) *` the shipped coefficients then descended over f32
quantisation -- **0.1672 idealized ulp**, against 0.4643 for the unscaled
poly on `erfinv`'s grid, i.e. the scaling is free (naive `sqrt(2)*c`
rounded to f32 without a descent is already 0.1923, so even that would
do):

    1.2533141, 0.3281154, 0.18047298, 0.12070113, 0.109420195,
    -0.026161268, 0.37856013, -0.50125736, 0.4776894

Caveat to check before believing the ~4.3: `probit`'s result is
`sqrt(2)` times `erfc_inv`'s, so the same *relative* error lands on a
different ulp grid -- 1.414x the ulp count where the binade does not
change and 0.707x where it does. The 3.60 -> ~4.3 estimate assumes the
max stays on a non-crossing sample, which is where it is today but is not
guaranteed after the multiply comes out.

### A strided estimate under-reports the max, on all three, every time

The before/after table above is a *strided* f32 sweep (every 128th bit
pattern) against scipy. The harness's exhaustive sweep, run afterwards on
the same code, is worse on all three:

| | strided-scipy max | exhaustive max |
|---|---|---|
| `erfinv` | 3.67 | **4** |
| `erfc_inv` | 4.50 | **5** |
| `probit` | 4.98 | **5** |

That is the expected failure mode and not a contradiction -- a stride of
128 walks past the worst bit pattern -- but it is worth writing down with
numbers, because the gap is 9-11% and in the same direction every time. A
strided scan with an independent oracle is a fine instrument for a *delta*
(same points both sides, and the deltas here held up) and is not one for
an absolute max. Quote the exhaustive row.

Exhaustive before/after exists for `erfinv` only (0.3692 avg / 6 max ->
0.3642 / 4); the earlier baseline run was killed after that first row, and
`erfc_inv`/`probit` have exhaustive numbers on the new code only.

## `exp_reduce!` range-halving: the last untried lever on the erfc family, priced and rejected

`erfc` (max 7), `erfcx` (6), `norm_cdf` (7) and `gelu` (9) are all
**`exp`-bound** -- the attribution entries above establish that, and a
fresh quick sweep on current master (2b51ddb) reproduces every published
number, so the diagnosis is not stale: `erfc` 0.1289 / 6q, `erfcx` 0.1439 /
6q, `norm_cdf` 0.0657 / 6q, `norm_pdf` 0.0269 / 3q, `erfcx (x>=20)` 0.2683
/ 2. Every `erfcx_pos`-side lever is closed (offset in both directions,
Horner, degree drop, even/odd split, leading-term peel, EFT on `v`,
two-branch split), so the only remaining move is to make `exp` itself
better for these four callers.

`exp_r_poly!` is core-locked and its degree bump is already recorded as a
real win at a real cost. But **`exp_reduce!` is not core-locked** -- its
only five call sites are `exp_checked`, `erfc`, `erfcx`, `norm_cdf` and
`norm_pdf`, i.e. exactly one domain -- so its *reduction* can be changed
without a core claim even though its polynomial cannot. That leaves one
classical lever nobody has priced here: **halve the Cody-Waite range.**

Reduce against `ln2/2` instead of `ln2`, so `|r| <= ln2/4`. The degree-5
minimax truncation error scales as `r^6`, so it drops **64x**:

| | `|r| <= ln2/2` | `|r| <= ln2/4` |
|---|---|---|
| deg 5 minimax, `e^r` | 7.5e-8 = **1.26 ulp** | 1.2e-9 = **0.02 ulp** |
| deg 4 minimax, `e^r` | 1.2e-6 = 20 ulp | 8.1e-8 = 1.37 ulp |

The 1.26 matches the headroom table's own 1.354 for `exp_r_poly`, which is
the check that the estimate is calibrated. So the accuracy case is real
and is *stronger* than the degree-6 bump (which buys ~40x, not 64x).

**It is the reconstruction that kills it.** `2^(k/2)` for integer `k` is
`2^floor(k/2)` times `1` or `sqrt(2)`, and `exp2_field_split` only builds
integer powers. Getting the odd half back costs, on top of the shipped
chain: `vroundps` (floor of `k/2`), one fma (the parity residual), a
compare and a blend (select `1.0` vs `sqrt(2)`), and **two** ops for the
scale itself -- `fma(p, S_HI, p*S_LO)` -- because a single-word `sqrt(2)`
carries 0.44 ulp of relative error, which would hand back a third of what
the range-halving just won. That is **~6 instructions**.

Priced the way idea #89's gather screen prices things, against the thing
being improved: the entire degree-5 `exp_r_poly!` is 5 `vfmadd231ps` =
**2.5 cycles Block RThroughput per 8 elements**. Six added ops are ~3
cycles. **The reconstruction costs more than the whole polynomial it is
making more accurate**, on all five callers including `exp_checked` and
therefore its own 20+ downstream users. For comparison the degree-6 bump
is +1 fma (+0.5 cyc) and was rejected at +5-23% throughput. Not
implemented.

Row 3 of the table is the other half of the verdict: **deg 4 at half range
is a dead heat with deg 5 at full range** (1.37 vs 1.26 ulp), so the
"spend the saved degree on the reconstruction" variant buys nothing at all
-- it trades one fma for six ops at equal accuracy.

### Two smaller reduction-side ideas, also priced and not taken

- **Kill the second Cody-Waite rounding.** `r` is rounded twice (`fma(-k,
  LN2_HI, x)` then `fma(-k, LN2_LO, r)`), each up to `ulp(0.347)/2 =
  1.5e-8`, i.e. 0.25 ulp each in the result. Carrying the low part
  multiplicatively instead (`rl = -k*LN2_LO`, then `fma(p, rl, p)` after
  the poly) removes the second one for +1 op and +1 dependency level.
  Worth 0.25 ulp of `exp`'s 3 -- below the noise of what the callers see,
  and the added level lands on the critical path of all five.
- **Move `erfc`'s `pe` into the exponent** instead of applying it as
  `fma(-r, pe, r)` on the `erfcx` factor. Op-count wash (erfc loses one
  fma, `exp_reduce!` gains one add) and accuracy wash: the linearisation
  `e^-pe ~ 1-pe` already errs by `pe^2/2 < 1e-11`, and `r1 - pe` inside
  the reduction introduces a fresh `1.5e-8` rounding of its own. Nothing
  on either axis.

**Transferable:** `exp_reduce!`'s five call sites are one whole `jm`
domain, so a future instance holding `exp_checked` can change the
*reduction* freely -- only the poly needs the core lock. That is a real
degree of freedom nobody had noticed; it just does not happen to contain a
win, because every way of buying accuracy in a reduction costs more
reconstruction than an f32 polynomial is worth.
## `dawson`'s tail was paying full-weight roundings for a correction worth 3% -- peeled, shipped

**Shipped 2026-08-02.** The leading-term peel, applied to `dawson_tail_poly`.
IDEAS.md lists the peel as one of the two shapes with the best hit rate here,
and this is the cleanest instance of it the file has: the branch went from
"fit error plus a chain of roundings" to "fit error", exactly, and got one
multiply cheaper doing it.

### The screen that found it

A four-row decomposition over `|x| > 4` (stride 16, ~1.3M points, the
`accuracy.rs` series reference in scalar f64), separating the *structure*
from the *coefficients*:

| row | what it is | avg | max |
|---|---|---|---|
| shipped chain | `R(z)*w` as written | 0.4637 | **5** |
| fit-only | shipped f32 coefficients evaluated exactly in f64, one final rounding | 0.2842 | **3** |
| oracle | correctly-rounded `R`, i.e. `fl32(d/w)*w` | 0.2944 | **1** |

So of the tail's 5 ulp, **2 were the evaluation chain and 2 were the fit**,
and only 1 was structural. That 2-ulp chain is what the peel removes -- all
of it, not some of it.

### Why it is total rather than partial

`R(z) = 1 + z*T(z)` with `z = w^2 <= 1/64` on the whole branch, so `z*T` is
never more than **0.033** of the result. Peeling the `1` out of the
polynomial and into the final `fma`'s addend puts every rounding inside `T`
at the scale of the correction instead of the scale of the answer: they
arrive attenuated 30x or more, and the final `fma` carries the only
full-weight rounding left. The measured peeled row is `0.284184 / 3` against
a fit-only floor of `0.284190 / 3` -- identical to five decimal places. There
is nothing left in this branch that is not its coefficients.

### It is also one op cheaper

`R(z)` then `* w` is 4 fma + 4 mul. `fma(w*z, T(z), w)` is 4 fma + 3 mul --
the trailing multiply by `w` is absorbed into the `fma` that adds the peeled
term, and `T` is one degree shorter than `R`.

| | before | after |
|---|---|---|
| `dawson_throughput` instrs | 84 | 84 |
| region `vmulps` | 14 | **12** |
| region fma (213+231) | 32 | 32 |
| uOps | 92 | 93 |
| Block RThroughput | 23.00 | **22.00** |
| cyc/elem | 1.701 | 1.645 |
| `dawson_latency` | 62.111 | 62.002 |

Instruction count is *flat* despite two fewer multiplies: LLVM spent the two
slots on `vmovaps` (6 -> 8), which is also where the +1 uOp comes from. Read
the histogram, not the total -- the arithmetic really did drop by one
multiply per element, and `Block RThroughput` is the rung that shows it.

### Accuracy

| | avg | max |
|---|---|---|
| tail only, `|x| > 4` | 0.4637 -> **0.2842** | 5 -> **3** |
| whole function, exhaustive all 2^32 | 0.0581 -> **0.0502** | 6 (unchanged) |

The headline max does not move because it never lived here: it is the
central branch at `x = 1.4041231`. Six `worst_corpus` entries move, all
`|x| > 4`, all by exactly 1 ulp and all *toward* the true value
(`dawson(100)`: `5.0002495e-3` -> `5.00025e-3` against a true
`5.00025004e-3`).

Bit-identical for `|x| > 2897`, and for the same reason as before the change:
the result is `w` there once the correction falls under half an ulp of `w`,
and `z*T < 2^-24` is the same threshold that used to round `R` to exactly
`1.0`. The `z = +inf` discarded arm still combines to a consistently-signed
infinity rather than a `NaN` -- `T`'s coefficients are all positive and
`w*z` only carries `w`'s sign in.

## The same peel does **not** transfer to `dawson`'s central branch -- rejected, measured

**Rejected 2026-08-02**, measured in the same run as the tail peel above, so
the cost of knowing was zero.

`x*P(u)/Q(u)` invites the identical treatment: `P/Q = 1 + u*S(u)/Q(u)` with
`S = (P-Q)/u`, giving `fma(x, corr, x)`. It is **free** in op count -- the
final `x * ratio` multiply becomes the peel's own `fma`, and `N = P-Q` is the
same degree as `P`. It is also *better conditioned* than `P`: at `u = 16`
every term of `N` is negative (sum `-631.4`, condition 1.00) where `P`'s
terms alternate (sum `21.1`, largest term `13.0`, condition 1.29).

None of that matters, because the attenuation runs the wrong way. The tail's
correction is 3% of its result; the central branch's is
`|ratio-1|/|ratio|`, which is 0.18 at `x = 0.5` but **30x** at `x = 4`,
where `ratio = 0.0323`. Measured as an oracle row (exactly-computed
correction, so this is the *floor* of the idea, not an implementation):

| octave | oracle peel avg/max | shipped avg/max |
|---|---|---|
| `2^-10 .. 2^-7` | 0.0000 / **0** | 0.43 / 2 |
| `2^-4` | 0.0017 / 1 | 0.60 / 4 |
| `2^-2` | 0.0285 / 1 | 1.33 / 5 |
| `2^-1` | 0.0931 / 1 | 1.01 / 4 |
| `2^0` | 0.679 / 2 | 1.09 / 5 |
| `2^1` | 2.850 / **8** | 1.00 / 5 |
| `2^2` (tail's range) | 11.90 / **32** | 0.81 / 5 |

Break-even is where `ratio = 0.5`, i.e. `x ~ 1.06`. Below it the peel is
better than the shipped chain *and* better than the shipped chain's own
oracle; above it, it loses fast.

**What this leaves open, and what it closes.** Closed: a single peeled form
for the whole central branch. Open, and priced: a *third* branch,
`fma(x, u*S(u)/Q(u), x)` for `|x| <= 1`, sharing `Q` with the existing
rational. Its oracle floor over `2^-13 .. 2^0` is max 1 against a shipped
2-5, so it is worth roughly 0.014 of the exhaustive avg (0.050 -> ~0.036)
-- but **zero** of the headline max, which is at `x = 1.404` on the other
side of the split, and it costs a select plus a second numerator on every
call. Recorded rather than taken: an avg-only win that adds a branch to a
7-op hot path is the wrong side of this crate's Pareto rule.

### Transferable

The peel's gain is not a property of the polynomial, it is
`1 / (1 + |peeled term| / |result|)` -- so **screen it by evaluating that one
ratio at the far end of the branch's domain before writing any code.** Under
~0.1 (the tail: 0.033) it removes the entire evaluation chain; over ~1 it is
an amplifier. That single number would have predicted both results here, and
it also explains the `erfcx_pos` peel's ~1%: same lever, ratio in between.

## `sinc` is *more* accurate than a correctly-rounded numerator and denominator would make it -- the shared rounding is the mechanism, and improving either operand alone is a regression

**Rejected 2026-08-02**, all four candidate directions, on one measurement.
`sinc`'s max 4 had been attributed before at three fuzz worst-case points
(the division ~2 ulp, `sinpi` ~1-2, "a genuine mixed contribution"). A
systematic per-octave attribution says something the three-point probe could
not see, and it inverts the conclusion.

### The measurement

Five rows over `|x| in [1e-3, 1e6]`, stride 64, against an f64
`(-1)^k sin(pi*r)/(pi*x)` reference reduced on the exact `r = x - round(x)`
so it stays relatively accurate at `sinc`'s own zeros:

| row | what it is | avg | max |
|---|---|---|---|
| shipped | `sinpi(x) / (PI*x)` | **0.394** | 3 |
| exact numerator | correctly-rounded `sin(pi*x)`, shipped denominator | 0.472 | 3 |
| exact denominator | shipped `sinpi`, correctly-rounded `pi*x` | 0.510 | 3 |
| both exact | correctly-rounded operands, one division | 0.368 | 2 |
| full oracle | correctly-rounded `sinc` | 0.000 | 0 |

**Improving either operand on its own makes the function worse** -- by 20%
and 29% on avg respectively. Exhaustively over the no-reduction region
`|x| in [1e-3, 0.5)` (every bit pattern, not sampled) the effect is much
larger and the shipped form beats *all three* alternatives including the
correctly-rounded one:

| row | avg | max |
|---|---|---|
| shipped | **0.3585** | **2** |
| exact numerator | 0.5959 | 3 |
| exact denominator | 0.6739 | 3 |
| both exact | 0.4636 | 2 |

### Why

For `|x| <= 0.5`, `sinpi` does not reduce (`k = round(x) = 0`, `r = x`), so
the argument it hands its own polynomial is literally the expression
`PI * x` -- the same one `sinc` writes as its denominator, and the compiler
CSEs them into a single multiply. The quotient is therefore `sin(t)/t`
evaluated *self-consistently* at `t = fl(PI*x)`, i.e. it computes
`sinc(t/PI)` exactly where it should compute `sinc(x)`. The error that
introduces is `d(ln sinc)/d(ln t) * relerr(t)`, and **that log-derivative
vanishes at the origin** (it is `-t^2/3 + O(t^4)`, about `-3e-6` at
`x = 1e-3`). So the entire argument error -- both `PI`'s own `2.78e-8`
representation bias and the multiply's rounding -- costs nothing at all.
Two independently correct roundings do not cancel and each spend a half ulp.

### What this closes

- **A two-word `pi` in the denominator** (`fma(PI_LO, x, PI*x)`). This was
  the obvious next move and it is the "exact denominator" row: worse, on
  both axes. It also does not even do what it claims -- the fma adds the
  `PI_LO*x` term but does not remove `fl(PI*x)`'s own rounding, so it trades
  a `2.78e-8` bias for a fresh `2^-24`.
- **Any accuracy work inside `sinpi` aimed at `sinc`** -- the "exact
  numerator" row. It would help `sinpi` and hurt `sinc`.
- **Double-wording both sides.** The "both exact" row is the ceiling for
  that whole family and it is *below the shipped function* over the
  no-reduction region, and worth 0.026 avg / 1 max over the full range --
  for a compensated division, which is independently closed above on mca
  cost (+12.3% throughput) and on a real `0*inf -> NaN` at the denormal
  floor.
- Above `|x| ~ 1` the sharing is gone (`sinpi` reduces, so its polynomial
  argument is `PI*r`, not `PI*x`) and the two errors simply add: shipped
  0.42-0.44 avg / max 3 against a "both exact" floor of 0.34 / max 2 in
  every octave. **That is where `sinc`'s max lives** (worst sampled point
  `x = 1.4538369`), and it is ~1 ulp deep against a floor that needs both
  operands exact. Nothing cheap reaches it.

Documented in `sinc`'s doc comment as a load-bearing invariant rather than
left implicit: the property is invisible in the source -- it lives in the
fact that two *textually different* expressions (`sinpi`'s internal
reduction and `sinc`'s denominator) evaluate to the same rounded value --
and a refactor to `sinpi(x) * (1.0/(PI*x))`, or a `sinpi` that reduced
unconditionally, would silently give back 0.24 avg ulp with nothing failing.

### Transferable

**Before trying to make one operand of a ratio more accurate, check whether
it shares a rounding with the other.** Where it does, correlated error can
cancel to *below* what correctly-rounded operands would give, and the usual
instinct -- improve the sloppier-looking side -- is a regression. The tell
is cheap and it is the one this file already recommends for a different
reason: score an oracle row per operand *separately*, not just a combined
one. Here the combined oracle looks like 2 ulp of headroom while the two
single-operand oracles are both *negative*, and only running all three shows
which it is.

## `tanpi`: peel `pi` out of the polynomial, not into it -- max 5 -> 2, avg 7x, and an instruction cheaper

`tan_poly` was a degree-6 minimax of `tan(t)/t` in `t = pi*r`, with
`tan_core` forming `t = fl(PI*r)` and `aphi = fl(PI*(0.5-|r|))` first. It
is now the *remainder* after the leading term:

    B(u) = tan(pi*w)/w - fl(pi),  u = w*w,  w in [0, 0.25]
    tan(pi*w) == fma(w, PI, w*B(u))

Same degree, same Estrin grouping, **one instruction fewer**, and the two
`PI*` multiplies are gone -- `tan_core` now takes `r` and `s = 0.5-|r|`
directly, both exact, and the seam is `|r| <= 0.25` exactly instead of a
comparison against `fl(pi/4)`.

| | avg ulp | max ulp |
|---|---|---|
| `tanpi` before | 0.2267 | 5 |
| `tanpi` after | **0.0327** | **2** |

llvm-mca, `tools/mca_region.py`:

| region | instrs | uOps | BlockRT | cyc/elem |
|---|---|---|---|---|
| `tanpi_throughput` | 94 -> **93** | 100 -> **98** | 26 -> 26 | 2.223 -> **2.155 (-3.1%)** |
| `tan2pi_throughput` | 96 -> **95** | 102 -> **101** | 27 -> 27 | 2.405 -> **2.342 (-2.6%)** |
| `tanpi_latency` | 3140 -> 3074 | 3530 -> 3654 | 832 -> 832 | 76.00 -> 77.05 |

Full-file asm region diff: **4 of 313 regions changed**, exactly
`tanpi`/`tan2pi` in both modes. Nothing else in the crate moves --
`tan_poly` and `tan_core` have no other caller.

The throughput rows are a clean take by the ladder (instrs down, uOps
down, RThroughput flat). **The latency row's +1.4% is the usual
artifact and was arbitrated rather than believed:** jmp/jcc is 256 on
both sides, so it is not arm concatenation this time; the uOp rise is
*load* uOps. The scalar latency harness had `PI` in a register
(`vbroadcast` 1 -> 0) feeding 640 `vmulss`; now 512 `vmulss` + 896
`vfmadd`, and 128 of those fmas take `PI` as a RIP-relative memory
operand -- +1 uOp each, 2 per element, off the dependency chain. The
critical path is 7 levels deep either way (`r -> t -> t^2 -> poly ->
mul` becomes `r -> r^2 -> poly -> mul -> fma`), which is the structural
reason to disbelieve the +1.4%.

### Why this was sitting there, and the rule it sharpens

The graveyard already contained a screen of "fold `pi` into `tan_poly`'s
coefficients", and that screen **rejected it on a correct measurement of
the wrong term**: instrumented over 20M samples, the argument-rounding
contribution was avg 0.44 ulp of a 2.05 total, only 0.02% of the
above-3-ulp samples had it contributing even half, so "folding `pi` into
the coefficients -- which would delete both effects and two `vmulps` --
cannot pay for a refit, and is not attempted."

Every number in that screen is right. It priced **deleting the `PI*r`
rounding**, which is indeed worth ~0.4 ulp. What it did not price is the
thing the peel actually collects:

**B's own error reaches the result attenuated by `w*B/tan(pi*w)`** --
zero at `w=0`, at most **0.215** at `w=0.25`. A `t*Q(t^2)` form has no
attenuation at all, because there the polynomial carries the *entire*
value. So the peeled polynomial has to be right to ~4.7x fewer bits for
the same result, and that slack is where the 5 -> 2 comes from. The
argument rounding is a rounding error; the attenuation is a change of
*sensitivity*, and only the first one shows up in an
error-attribution instrument.

This is also **not** the `pi`-into-the-coefficients fold twice rejected
for `sinpi` (idea #43) and `sind` (idea #44). Those keep the shape
`w*S(u)` with `S(0) = fl(pi)`, and `fl(S(0)*w)` is bit-for-bit the same
operation as `fl(PI*w)` -- the leading rounding is reproduced exactly,
so the transform buys nothing and the refit's risk is all there is. That
is why both regressed. Peeling `PI` into a closing `fma` is a different
transform with the same name: `fma(w, PI, ...)` forms `pi*w` **exactly**
and rounds once, at the end.

Oracle screen, before any fitting (probe over `r` in `[0,0.5)`, ideal
polynomial in each form, f64 reference):

| | shipped form, ideal Q | peeled form, ideal B |
|---|---|---|
| direct arm, `r` in [0, 0.25] | 2.810 | **0.962** |
| reflected arm, `r` in (0.25, 0.5) | 3.626 | **1.561** |

The 3.626 reproduces the 3.7 oracle floor this file already recorded for
the shipped chain, which is the check that the probe is calibrated. The
peel moves the *floor*, which is what no coefficient work could do.

### Fitting notes

Weighted minimax LP (scipy/HiGHS) of `B` against `u`, weight
`w/tan(pi*w)` -- the factor turning an absolute error in B into a
relative error of the result, and **the same expression for both arms**,
so one fit serves the direct `w*B` and the reflected `1/(w*B)` alike.

- **The LP must be column-scaled** (`z = u/0.0625`). In raw `u` at degree
  6, HiGHS returns a coefficient pinned to its bound and a residual of
  1.61e-6 -- 27x too large -- and reports the same number for degree 7,
  which is the tell. Same failure the `sind` entry hit on `d^9`.
- Degree 5 quantises to **4.78** ulp-equivalent, degree 6 to **1.115**,
  degree 7 no better than 6. Degree 6 idealises below the LP's own
  resolution.
- `c[0]` is pinned to `pi - fl(pi)` exactly. Free here, unlike
  `erfcx_pos`'s c0: the *unpinned* LP converges to -8.742276e-08 against
  the exact -8.742278e-08, so the pin costs nothing and buys an exactly
  reproduced asymptote as `w -> 0`.

### What this retires

The **degree-6-on-the-Horner-spine variant** recorded above (max 4, avg
0.2309, latency 88.72 vs 76.00, "too thin a pareto point to carry an
API") is now strictly dominated and needs no tier: the peeled form is max
**2**, avg **0.0327**, and *cheaper* than the form that variant was
losing to. Deleted from consideration rather than shipped as
`tanpi_accurate`.

### Transferable

**A rejection that instruments "how much does this rounding contribute"
has only closed the rounding, not the restructure.** The question that
finds these is not "what does this step's error cost" but "**what
fraction of the result does the polynomial carry**" -- if it carries all
of it, peeling the leading term off is a sensitivity change worth several
bits, independent of any rounding it also deletes. `log_2` (3 -> 1 ulp),
`erfcx_pos`'s two-word `c0`, and now `tanpi` are three instances.

Screening shape that made it cheap: an out-of-tree probe comparing the
two forms' *oracle floors* (ideal polynomial in each) settled it in one
run, before a single coefficient was fitted.

### Exhaustive confirmation and a second, independent one

Harness, all 2^32 bit patterns: `tanpi` **avg 0.0327 / max 2**, worst
`x = 3.2719934e-1` (was 0.2267 / 5). `tan2pi` quick-fuzz 0.2273 / 5 ->
**0.0324 / 2**; `sin2pi`/`cos2pi` unmoved at 2, as the asm diff requires.

Cross-checked without the harness, because it was worth knowing whether
the number depended on the reference: **every distinct internal state
`tanpi` can reach occurs for `x` in `[0, 0.5]`** -- there `round_ties_even(x)`
is 0, so `r == x` -- and larger `|x|` only ever produce `r` on coarser
subgrids of that same set. Sweeping all 1056964607 f32 in `[0, 0.5]`
against an f64 `tan(pi*x)` gives **max 2.1198 fractional ulp at
x = 4.1912597e-1** (in the reflected arm, whose oracle floor is 1.561),
avg 0.2552 over that range. The harness reports integer bit-distance and
says 2; the two agree.

That reachability argument is reusable for any function whose reduction
is `x - round(x)`: the exhaustive answer lives in one half-turn, and a
single-threaded probe over it costs a couple of minutes with no bench
lock, which matters when three instances are queued behind one.
## `log2p1`/`log10p1`: a two-word scaling constant, avg -65% and -73% -- shipped

**Shipped 2026-08-02.** Both functions compute `log1p`'s Sterbenz-exact
correction and then scale it: `corr = (c/u) * LOG2_E`. That one multiply
carried the whole gap between them and `log1p`.

### Why the constant is not diluted here

The usual reason a scaling constant's error does not matter is that it only
ever multiplies an already-small correction. That reasoning fails on exactly
the region that dominates by sample count: for `|x| < 2^-24`, `1+x` rounds to
exactly `1.0`, so `c = x`, `c/u = x` exactly, the log kernel returns an exact
`0.0`, and the final `0.0 + corr` is exact. The answer **is** `fl(x*log2(e))`
and nothing else. A one-word constant's fixed relative offset therefore
arrives undiluted, as bias rather than noise, over ~40% of all bit patterns
(exponents `-126..-25`, both signs). `log1p` has no such term -- its `corr`
is `c/u` with no scaling, so it returns `x` *exactly* there, which is the
whole of its 3-4x lead.

Sizes: `LOG2_E` is 1.334976e-8 low relative, `LOG10_E` 2.326313e-8 high,
i.e. 0.11-0.22 and 0.20-0.39 ulp of the result. The crate had already priced
the same constant independently -- `log10_normal`'s peel comment calls
`LOG10_E`'s inexactness "a fixed 0.39 ulp-equivalent, which is the floor this
fit sits on".

### The fix, and the shape that matters

    let cu = c / u;
    let corr = fma(cu, LOG2_E, cu * LOG2_E_LO);

Big product **inside** the fma, low word as the addend -- the same idiom as
`erfcx_pos`'s two-word `1/sqrt(pi)`. That is one rounding for the whole
scaled correction and no bias. The naive-looking transposition
`fma(cu, LOG2_E_LO, cu * LOG2_E)` is *worse than the original*: it rounds the
big product on its own first and then rounds again, trading a 0.13-ulp bias
for a fresh half ulp. Low words: `LOG2_E_LO = 0x32a57060`,
`LOG10_E_LO = 0xb22d91af` (note the sign -- `LOG10_E` is high, so its low
word is negative, which is why `cu = +-inf` now yields `NaN` rather than
`inf`; the existing `is_finite` guard maps both to `0.0`, so the edges are
unchanged and edgecheck's `log2p1(3)==2.0` / `log10p1(99)==2.0` exactness
pins still hold).

### Measured

Exhaustive over all 2^32, same binary either side, run unlocked:

| | avg before | avg after | change | max |
|---|---|---|---|---|
| `log2p1` | 0.0919 | **0.0324** | **-65%** | 2 -> 2 |
| `log10p1` | 0.1458 | **0.0387** | **-73%** | 2 -> 2 |

against `log1p`'s own 0.025, which is the bar these two now sit next to
rather than at 3-6x. `log10p1`'s worst-`x` moves out of the tiny region
entirely (5.3175995e-7 -> 6.306496e-2), which is the tell that the region is
now correctly rounded. Two of 9630 `worst_corpus` entries move, both
`log10p1` near 1e-10, each 1 ulp toward the true value.

**The max does not move**, and would not: it lives in the
`log_2_normal`-dominated mid-range, not in the region this fixes. This is an
avg-only win.

### Cost -- real, and it is the price of the trade

+1 arithmetic op each. Every rung agrees, so this is not an mca artifact:

| | instrs | uOps | BlockRT | cyc/elem | latency |
|---|---|---|---|---|---|
| `log2p1` | 100 -> 103 | 110 -> 112 | 22 -> 23 | 1.876 -> 1.934 (+3.1%) | +1.2% |
| `log10p1` | 103 -> 106 | 112 -> 115 | 22 -> 23 | 1.936 -> 1.977 (+2.1%) | +0.6% |

Taken as a single function rather than split into tiers: a 2.8-3.8x avg
accuracy gain for 2-3% throughput is the trade this crate makes, and a 3%
perf delta is explicitly too small to justify a second tier.

### Also fixed: a stale doc number

`log2p1`'s doc comment claimed "avg 0.102, max 3"; the readme said 0.092 /
max 2. The before-run settles it at **0.0919 / max 2** -- the readme was
right and the doc comment was stale by a max.

### Transferable

**A constant that "only scales a small correction" is not safe until you
check whether some region makes the correction the entire answer.** The
argument that protects it here (`corr` is dwarfed by the `k`-dominated log
term) is true everywhere except the region holding 40% of the inputs, where
the log term is identically zero. Look for the input set that annihilates the
*other* term. The same question is now open on `log10_normal`'s own
`fma(s, LOG10_E, sq)` -- see IDEAS.md, where the identical two-word split
would remove the 0.39-ulp floor its comment documents.
## `expm1`: fit `e^r - 1`, not `e^r`, and the seam disappears with the Pade

Idea #169 profiled `expm1`/`exp_m1_over_x`/`tanh` and found all three had
their error concentrated in a narrow band at and just above their seam,
traced it to the direct arm forming `e^x` and subtracting 1 (losing
`log2(e^x/(e^x-1))` bits, a factor peaking at **2.541** just above
`x = 0.5`), closed three levers on it -- (a) widen the Pade's fit,
(b) a doubling identity as a third branch, (c) retune the seam -- and
concluded the concentration was "structural at the current op budget".
`tanh` escaped that verdict by questioning a premise none of the three
levers touched. So does this one, and it is the same *shape* of premise:
**the direct arm does not have to compute `e^r` at all.**

**Shipped**: exhaustive avg **0.1291 -> 0.0166**, max **5 -> 2**, and
throughput **1.604 -> 1.087 cyc/elem (-32.2%)**. Not a tradeoff on
throughput; see the latency note at the bottom for the one axis that moves
the other way.

### The mechanism

The combine wants `2^k * e^r - 1`. Carry `E = e^r - 1` instead of
`p = e^r` and it becomes `2^k*E + (2^k - 1)`, one fma, with `2^k - 1`
exact for every `k <= 24` and under a quarter ulp of the answer above
that. The error analysis inverts:

| | error carried into the combine | amplification at `x = 0.5` |
|---|---|---|
| `p = e^r` | absolute, ~`ulp(1)` regardless of `r` | `e^x/(e^x-1)` = **2.541** |
| `E = e^r - 1` | absolute, ~`ulp(E)`, i.e. proportional to `E` | `2^k E/(2^k(1+E)-1)` < **1** |

The second row *de*-amplifies. That is the whole result: the polynomial's
own evaluation roundings shrink with `|E|` exactly where the combine's
sensitivity grows, and the two cancel instead of compounding.

### What that deletes

The Cody-Waite reduction already lands every input on `|r| <= ln2/2`, and
the new form is relatively accurate over that whole range, so **there is
no near-zero branch left to have a seam**. Gone: the Pade approximant,
its division, the `|x| < 0.5` select, and `exp2_field_split` (one field at
`k-1` suffices, `expm1_checked`'s own trick). At `k = 0` -- every
`|x| < 0.3466`, most of the old Pade branch -- the addend is exactly 0 and
the answer is the polynomial itself, unrounded.

`expm1_r_poly!` is `r + r^2*P(r)`, `P` degree 4, fitted by ulp-weighted LP
against `(e^r-1-r)/r^2`. Degree 4 is the natural stop: the degree-5
coefficient converges to exactly 0, the same signal `cbrt` and
`exp_pos_neg` gave at their own degree bumps. Same instruction count as
`exp_r_poly!`, same Estrin fold, `r^4` never formed.

### Numbers

Exhaustive over `exp_domain` (all 2.24e9 in-range f32, not sampled):
avg 0.004919, **max 1.5642** at `x = 3.7117928e-1`. `expm1`,
`expm1_narrow` and `expm1_checked` agree bit-for-bit over every pattern in
their shared domains.

| region | instrs | uOps | BlockRT | cyc/elem |
|---|---|---|---|---|
| `expm1_throughput` | 79 -> 56 | 86 -> 58 | 24 -> 15 | 1.604 -> **1.087** |
| `expm1_checked_throughput` | 78 -> 61 | 87 -> 64 | 22 -> 17 | 1.556 -> **1.320** |
| `expm1_narrow_throughput` | 67 -> 50 | 73 -> 52 | 19 -> 14 | 1.278 -> **0.964** |

All three rungs of the ladder move the same way, so this is not the
throughput-column artifact. The size of the throughput win is itself a
consequence of the branch deletion rather than of the op count alone: the
vectorized loop is if-converted, so the old code paid for the Pade's
`vdivps` on **every** lane including the ones that selected the other arm.

### Two things that had to be paid for, and one that did not

- **The `k-1` field makes the pre-doubling value denormal.** `b` is
  `expm1(x)/2`, so every result under `2^-125` is *formed* as a denormal
  and the doubling cannot restore the bit rounding took. Caught by an
  exhaustive bit-compare against `expm1_narrow` (which uses the field at
  `k` and so does not halve): **16777216 disagreements, exactly 2^24, the
  whole denormal input range**. `denormal_audit` would have caught it too.
  Fixed by an `|x| < 2^-125` select, which also carries `expm1(-0.0)`
  (at `k = 0` the addend is `+0.0` and `-0.0 + 0.0` is `+0.0`) -- both
  regions are exactly where `expm1(x) == x`, so one select covers both.
- **Latency, honestly.** The old `expm1_latency` published 69.00 but is
  one of the readme's flagged branch-shaped rows; its arms measure 32.00
  (Pade) and 46.00 (direct). The new region is **branchless** (0 jcc
  against the old region's 128) at 47.00. So `|x| >= 0.5` is +1 cycle and
  `|x| < 0.5` goes 32 -> 47 in a latency-bound scalar chain. Not split
  into a tier: throughput is the axis this crate optimizes and accuracy
  improves in both regions.
- **`expm1_narrow` did *not* become redundant**, which was the first
  guess. `expm1` now uses a single field too, but at `k-1`; `expm1_narrow`
  puts it at `k` and so writes `fma(e, t, t - 1.0)` with no doubling --
  genuinely one instruction cheaper, at the cost of the top of the domain
  (`k <= 127`, i.e. `x <= 88.37627`, versus `expm1`'s `88.72`). Same
  Pareto split it always had, different reason for it.

### Transferable

Any `_m1`/`m1_`-shaped function whose direct arm reconstructs `f(x)` and
then subtracts the leading term is a candidate, and the tell is a doc
comment explaining a *near-zero branch* as the fix for cancellation --
that branch is treating the symptom. `exp_m1_over_x` (max 5),
`exp2m1` (max 4) and `exp10m1` all still carry the Pade and the seam.

## `tand`: the peel transfers, but the win is the *signed* pole distance

`tand` was `sind(x)/cosd(x)`. It is now its own reduction plus `tanpi`'s
peel: `q = round(x/180)`, `d = x - q*180` in `[-90,90]`, then
`tan(d deg) == fma(d, DEG_TO_RAD_SMALL, d*B(d*d))` with
`B(u) = tan(w*pi/180)/w - fl(pi/180)`, degree 6, `c[0]` pinned exact.

| | avg ulp | max ulp |
|---|---|---|
| before | 0.1773 | 3 |
| after | **0.0383** | **2** |

`tand_unchecked` tracks it (0.0384), so the readme's "(+) bit-identical
on its domain" still holds -- it is the same `tand_core`, minus only the
out-of-range guard, which provably cannot fire in range.

### What unlocked it: a signed `s`, which the earlier attempt lacked

This file already rejected idea #128's direct-poly `tand` on cost:
"reusing `sind`'s own `d` for the pole-distance calculation loses
precision `d` never kept ... fixing that needed a whole second,
`cosd`-style independent reduction -- at which point real mca/quickbench
numbers no longer showed a clean win."

The diagnosis was right and the fix was overkill. `tanpi`'s `r = x -
round(x)` is exact, so `|r| <= 0.5` strictly; `tand`'s `q` comes from a
*rounded* `x*INV_180`, so near a pole `|d|` can land just past 90 (up to
~95.6 over the exact-reduction range). A **magnitude** pole distance
`90 - |d|` then goes negative and a closing `mulsign(.., d)` hands the
reflected arm the wrong sign -- which is what needs either a period fold
(~6 ops) or a second reduction (much worse).

Taking `s = mulsign(90, d) - d` instead makes `tan(d) = 1/tan(s)` hold
**with sign for `s` of either sign**, so the overshoot corrects itself
and the closing `mulsign` disappears too. Measured identical accuracy to
the period-fold version (max 1.87 vs 1.87 on a 1.0e8-point strided sweep)
at ~7 fewer ops. Sterbenz-exact wherever the reflected arm uses it
(`45 <= |d| <= 180`).

### Cost: the two rungs of the ladder disagree, and it is not arbitrated

`tand_throughput` (branchless on both sides, jmp/jcc 0):

| | instrs | uOps | BlockRT | cyc/elem |
|---|---|---|---|---|
| before | 95 | 97 | 30 | 2.533 |
| after | **101** | **107** | **26** | **2.290** |

Instruction and uOp counts say **+6% / +10%**; Block RThroughput and the
cycle count say **-13% / -9.6%**. `vdiv` is 2 on both sides and `vmul`
16 on both; the delta is +8 fma against the `sind`/`cosd` pair's integer
parity ops and two `.clamp()`s, which sit on different ports.

**Arbitrated: the structural rung was right, the counting rungs were
wrong.** Alternating wall-clock A/B on two prebuilt `quickbench` binaries
(same tree, only `lib.rs` swapped between this commit and its parent),
run under `jm bench`:

| band | pairs | latency | throughput |
|---|---|---|---|
| `Band::Two`, quiet box | 6 | new 6/6, min 14.76 vs 22.99 (**-36%**) | new 5/6, min 0.742 vs 0.762 (**-2.6%**), median ratio 0.947 |
| full period `[-180,180)`, quiet box | 8 | new 8/8, min 41.13 vs 51.62 (**-20%**) | new 8/8, min 1.774 vs 1.960 (**-9.5%**) |
| `Band::Two`, loaded box | 12 | new 11/12, min 26.45 vs 39.65 (**-33%**) | inconclusive, medians cross |

So the rewrite is **cheaper on both axes**, not neutral: throughput -2.6%
to -9.5%, bracketing mca's -9.6% cyc/elem. When the counting rungs and
Block RThroughput disagree *here*, believe Block RThroughput -- the +8
fma land on ports the removed parity/clamp integer ops were never
contending for, which is precisely what an instruction count cannot see.

The middle row is the load-bearing one. `quickbench`'s three `Band`s all
have `|x| <= 4` **degrees**, so every stock row takes `tand_core`'s direct
arm and never the cotangent reflection -- which flatters branchy code. A
full-period band mispredicts the `|d| <= 45` branch on ~half its inputs
and the rewrite still wins 8/8 on both axes, because the old form paid
*two* polynomials and a divide unconditionally where the new one pays one
polynomial, and a divide only on the reflected arm.

Method note, since the box moved under the measurement again: two runs an
hour apart differ by **1.7x in absolute ns** (load 11.45 on 8 cores, from
another worktree's accuracy sweep). Absolute ns is publishable only when
the *old* binary reproduces the readme's existing figure first -- old
measured 22.99 against a published 22.9, and 0.762 against a published
0.79, in the same run that produced the new 14.76/0.742. That calibration
is what licensed updating the two wall-clock rows; without it, don't.

**`tand_latency` is unusable and its readme row is withdrawn rather than
updated.** mca read 70.02 -> 37.349, i.e. *halved* while instructions
nearly doubled (1806 -> 3519). Cause: **jmp/jcc 0 -> 448.** The three
selects (`|d| > 128`, `|d| <= 45`, `s == 0`) lower to branches in the
scalar latency harness, and llvm-mca has no branch predictor -- the same
artifact this file records for `tanpi`, but firing *optimistically* here
instead of pessimistically.

The row stays withdrawn, but **the reason recorded first was wrong and is
corrected here.** It read: the new critical path is ~2-4 levels longer
(reduce -> guard -> poly -> mul -> fma -> divide in series, against two
parallel polys then one divide), so the true latency is likely slightly
*worse*. The A/B above measures it **-20% to -36%**. That chain is only
the *reflected* arm's; the direct arm has no divide at all, and the old
form's divide was unconditional. Counting depth through the longest arm
overstates a branchy function's latency the same way counting
instructions overstated its throughput.

So mca's 70.02 -> 37.35 (-47%) was directionally right and overstated by
~1.3x, not the "2x fantasy" it looked like. It is still withheld, for a
reason that survives the arbitration: **this column models branch-free
codegen for every other function in it**, and `tand_latency` is the one
region carrying 448 `jcc` with only 1 divide per element -- mca is
scoring a perfectly-predicted single arm. A number that is not comparable
to its own column's neighbours does not belong in the column, even once
its direction checks out.

Full-file asm diff: 4 of 313 regions, exactly `tand`/`tand_unchecked`.

### Transferable

- **A magnitude pole distance and a signed one are not the same
  reduction.** Whenever a reflection `f(d) = 1/f(s)` is fed a `d` from an
  *inexact* quotient, `s` must carry the sign or the overshoot needs a
  fold. The signed form costs one `mulsign` and saves the output one.
- **The peel's payoff scales with how bad the leading constant was.**
  `fl(pi)` is 0.47 ulp from `pi` and `tanpi`'s average improved 7x;
  `fl(pi/180)` is only 0.13 ulp from `pi/180` and `tand`'s max still
  improved but by way of the polynomial's attenuation, not the constant.
  Check the constant's own error first -- it predicts which half of the
  win is available.
- **A wall-clock A/B settles the ladder, and it can overturn the
  *structural* rung too, not just the counting ones.** Both readings here
  were wrong in the same direction: instrs/uOps overstated throughput
  cost, and critical-path depth overstated latency cost. Both errors have
  one cause -- reasoning about a *branchy* function as if every arm
  executed. Price a branchy candidate in a band that actually
  mispredicts; the stock harness bands may reach only one arm.
- **`c[0]`'s sign is load bearing for signed zero.** `fma(d, K, d*B(0))`
  returns `-0.0` for `d = -0.0` only because `B(0) > 0` here, so both
  addends are `-0.0`. `tan_poly`'s `c[0]` is negative, which is why
  `tanpi` needs its own `x == 0.0` pin. An ulp sweep is blind to this
  (zero-vs-zero scores 0); edgecheck is not.

## `exp_m1_over_x`: the removable singularity is `fma(r, P, 1.0)`

Same transform as the `expm1` entry above, and it pays twice here because
`(e^x - 1)/x` needed the Pade for a *second* reason: the division at
`x = 0`.

**Shipped**: exhaustive avg **0.0709 -> 0.0170**, max **5 -> 2**,
throughput **1.639 -> 1.279 cyc/elem (-22.0%)**.

The polynomial is `e^r - 1 = r + r^2*P(r)`, so the two things this
function wants are both one fma off `P`:

    (e^r - 1)/r  =  fma(r, P, 1.0)      <- r cancelled algebraically
     e^r - 1     =  fma(r2, P, r)

and the first is exactly `1.0` at `r = 0`, matching the true limit with no
`x == 0` select. Whenever `k == 0` -- every `|x| < 0.3466` -- the
Cody-Waite reduction leaves `r == x`, so that quotient *is* the answer and
no division happens at all. Only `|x| >= 0.3466` divides, where dividing
by `x` is unconditionally safe. So the select is on `k`, which the
reduction already computed, and there is no fitted seam on either side of
it: both arms are the same polynomial.

That is the general shape worth remembering: **a removable singularity
`f(x)/x` is free whenever the numerator's own polynomial carries `x` as an
explicit factor.** The Pade was doing this job (`x*N/D` divided by `x` is
`N/D`) but had to buy a division to do it.

| region | instrs | uOps | BlockRT | cyc/elem |
|---|---|---|---|---|
| `exp_m1_over_x_throughput` | 80 -> 61 | 88 -> 64 | 23 -> 16 | 1.639 -> **1.279** |
| `exp_m1_over_x_narrow_throughput` | 68 -> 60 | 75 -> 63 | 20 -> 15 | 1.347 -> **1.276** |

Latency, honestly, and the same shape as `expm1`'s: the published 81.00
was the branch artifact, arms 32.00 (Pade) / 58.00 (direct) then against
42.03 / 60.03 now. The near-zero arm is ~10 cycles longer because it runs
the full reduction where the Pade started from `x` directly. Throughput is
the axis this crate optimizes.

`expm1_r_poly!` split into `expm1_p_poly!` (returns `P`) plus a one-line
wrapper to make `P` reachable; all three `expm1` mca regions verified
**byte-identical** across that refactor, so it is purely structural.

## `log10_normal`: the two-word `LOG10_E`, which works and still does not pay

IDEAS.md filed this under *"would improve an existing function at zero perf
cost, so they clear the bar by construction"*. It does improve it. It is not
zero perf cost, and that is the whole result.

`log10_normal`'s own comment prices its floor: `LOG10_E` is 2.33e-8 high,
`Q` cannot absorb `(log10(e) - LOG10_E)/s` because that is a `1/s` term, and
it costs a fixed **0.39 ulp-equivalent**. `LOG10_E_LO` (already in `lib.rs`,
added for `log2p1`/`log10p1`) removes it: the closing
`fma(s, LOG10_E, sq)` becomes two-word and the singularity drops to
~2.6e-16/s.

### The low word alone is *worse than shipping neither*

This is the part worth remembering. Dropping `LOG10_E_LO` into the existing
chain without touching `Q` measures **worse** than the shipped code:
0.006941 -> 0.007195 avg over all positive normals at stride 251.

`Q` is fitted against `(log10(1+s) - s*A)/s^2` for whichever leading constant
`A` the code actually uses. A `Q` fitted for the one-word `A` therefore
already carries whatever compensation for that constant's error a polynomial
*can* express -- only the `1/s` part is beyond it. Adding the low word on top
double-counts everything the fit already absorbed. **A two-word constant
split is not a drop-in anywhere the polynomial beside it was fitted against
the one-word value.** (`log2p1`/`log10p1`'s split *was* a drop-in, because
there the constant scales a correction term with no co-fitted polynomial.)

### Refit, and only then does a degree bump pay

Same ulp-weighted LP as the shipped fit (weight `s^2/log10(1+s)`, scaled by
`2^24`, sequential f32 quantisation) reproduces the shipped degree-7
coefficients to 8 digits and its 0.5968 objective exactly, so the setup is
verified. Against the two-word constant:

| leading constant | deg | LP objective | note |
|---|---|---|---|
| one-word | 7 | 0.59679 | shipped |
| one-word | 8 | 0.39189 | stalls on the 0.39 floor -- the constant, not the fit |
| two-word | 7 | 0.50683 | singularity gone, degree still binding |
| two-word | 8 | **0.07434** | 5.3x past the one-word floor |

The floor is genuinely irreducible with one word, and the reason is worth
stating: the weighted singularity error is `|(log10(e)-A)/s| * s/log10(e)`,
i.e. **scale-invariant** -- it collapses to `A`'s own relative offset no
matter what `s` does. `A = fl(log10(e))` already minimises that. There is no
better single f32.

### Measured through the real chain

Standalone probe reproducing `log10_normal` exactly (it matches the
graveyard's own prior rows: shipped `k==0` 0.290768, stride-251 oracle
0.005622, both to 4 digits):

| chain | `k==0` octave, exhaustive | all normals, stride 37 |
|---|---|---|
| shipped (one-word, deg 7) | 0.290768 | 0.006899 |
| one-word, deg 8 refit | 0.184718 | 0.006389 |
| two-word, deg 7 refit | 0.213922 | 0.006546 |
| **two-word, deg 8 refit** | **0.057549** | **0.005843** |
| oracle (correctly-rounded mantissa term) | 0.000000 | 0.005601 |

Two-word deg 8 takes 80% of the `k==0` octave's error and **81% of the
entire oracle headroom** on the aggregate. Note in passing that the
comment's "degree 8 buys nothing" was a statement about the *fit objective*
and is not true of the measured chain -- one-word degree 8 is a real -36%
on the octave. It is just dominated.

Harness, `thorough` (exhaustive over all 2^32), two-word deg 8:

| | before | after |
|---|---|---|
| `log10` | 0.0034 / 1 | 0.0029 / 1 |
| `log10_unchecked` | 0.0069 / 1 | 0.0058 / 1 |
| `log10p1` | 0.0387 / **2** | 0.0304 / **1** |

`log10p1` becoming faithfully rounded is the only max that moves anywhere,
and it is real -- its worst points sit at `x ~ 0.06-0.07`, i.e. `u = 1+x` in
`log10_normal`'s `k == 0` octave at small `s`, exactly where this fit wins.

### Why it was not taken

+2 fma per call, and llvm-mca agrees on every rung (the two extra
`vbroadcastss` are the new `c[8]` and `LOG10_E_LO`):

| region | instrs | uOps | BlockRT | cyc/unit |
|---|---|---|---|---|
| log10_latency | 2885 -> 3013 | 3275 -> 3531 | 545.83 -> 588.50 | 48.157 -> 53.126 (+10.3%) |
| log10_throughput | 95 -> 101 | 103 -> 111 | 19 -> 21 | 1.617 -> 1.806 (+11.7%) |
| log10_unchecked_latency | 1802 -> 1932 | 1866 -> 1996 | 480 -> 544 | 38.063 -> 42.407 (+11.4%) |
| log10_unchecked_throughput | 60 -> 66 | 65 -> 72 | 16 -> 18 | 1.022 -> 1.149 (+12.4%) |
| log10p1_latency | 3279 -> 3408 | 3346 -> 3476 | 704 -> 768 | 52.017 -> 55.673 (+7.0%) |
| log10p1_throughput | 106 -> 111 | 115 -> 122 | 23 -> 25 | 1.977 -> 2.149 (+8.7%) |

So the whole win is one max ulp on `log10p1`, and `log10`/`log10_unchecked`
-- already max 1, i.e. already faithful, against a `std log10` that is also
max 1 -- pay 11-12% for an average they cannot spend. Set against the
precedent that filed this idea: `log2p1`/`log10p1`'s own two-word split
bought avg -65%/-73% for +2-3%, a ratio ~20x better than this one's -15%
for +11%.

**The intermediate is dominated, not a tier.** Two-word degree 7 (+1 op)
was measured end-to-end specifically to see whether the cheaper half reaches
the same max: `log10p1` 0.0387/2 -> 0.0317/**2**, throughput 1.977 -> 2.067.
It buys no max anywhere and costs 4.6-6.2%, so it is strictly worse than
both neighbours. No `log10_accurate` tier is warranted either: a second
public entry point whose only distinguishable property is one ulp on a
*third* function's max, both tiers being faithfully rounded, is API debt.

Ready to paste if `log10p1`'s faithfulness ever becomes the priority --
`Q` degree 8, refitted for the two-word constant:

    -0.21714719, 0.14476478, -0.108580895, 0.08686869, -0.07212288,
    0.06155919, -0.0575642, 0.056282505, -0.033275615

with `a = fma(s2, l0, fma(s, LOG10_E_LO, k * LOG10_2_LO))` and
`hi = fma(c[8], s2, l3)` feeding `u = fma(hi, s2, l2)`. Both low-word terms
ride the poly's low group, off the critical path, so the tail stays at two
fma and the poly depth stays at 4 -- the +2 is instruction count, not
dependency depth.

### Transferable

- **A two-word constant split is a drop-in only where no polynomial was
  co-fitted against the one-word value.** Otherwise the fit has already
  absorbed what it can and the split double-counts. Check for a co-fitted
  neighbour before assuming the `log2p1` result transfers.
- **A "fixed floor" in a fit objective is not a floor on the measured
  function.** Degree 8 against the one-word constant "reaches exactly the
  floor and buys nothing" in objective terms, and is still -36% on the
  octave that matters.
- **Price the cheaper half explicitly.** Here it was not a Pareto point but
  a dominated one, and only an end-to-end run showed that -- the fit
  objective ordering (0.507 vs 0.597) suggested it should have helped.
### Peeling a leading term (`fma(x, P-1, x)`) silently loses signed zero

**This trap is general, not an `erfinv` detail.** Any peel of the form
`x * P(u)` -> `fma(x, P_m1(u), x)` over a domain that includes `-0.0`
returns `+0.0` where the original returned `-0.0`, and the accuracy sweep
cannot see it.

The mechanism is forced, not a slip. A peel is only worth doing when
`P ~ 1`, which is exactly when `c0 - 1` is *negative*. So for `x = -0.0`:

    -0.0 * (P-1)  =  +0.0        (negative times negative zero)
    +0.0 + (-0.0) =  +0.0        (IEEE-754 round-to-nearest)

whereas `-0.0 * P` is `-0.0` and carried the sign for free. Both the
product's sign flip and the `+0 + -0` tie-break are required behaviour, so
no ordering of the fma recovers it.

The ulp sweep is structurally blind here: `+0.0` and `-0.0` compare equal,
so the difference scores 0 across an exhaustive 2^32 run. `examples/
edgecheck.rs`'s `check("erfinv(-0)", erfinv(-0.0), -0.0)` pin is what
caught it, on a change that was otherwise fuzz-clean, mca-clean and
exhaustively-verified-better. Keep that pin, and add the equivalent one
before peeling any other odd function.

The fix costs nothing where the function is odd and already signs its
other arm: build **both** arms on `|x|` and apply one `mulsign` to the
merged select, rather than letting the central arm carry a signed `x`
through the poly. In `erfinv` that is the same single `mulsign` the tail
arm alone used to pay -- `erfinv_throughput` uOps and `Block RThroughput`
both flat -- and it is bit-identical to the old form for every non-zero
`x`, since `fma` is sign-symmetric (`fma(x,p,x) == -fma(ax,p,ax)` for
`x < 0`).
## `exp2m1`: never rescale the argument, and round the reduction

Third application of the `expm1` transform, and the one where two
*separate* premises had to go.

**Shipped**: exhaustive avg **0.0769 -> 0.0400**, max **4 -> 2**,
throughput **1.843 -> 1.278 cyc/elem (-30.7%)**.

### Premise 1: it reached the approximant through `y = x*LN_2`

The old near-zero arm shared `expm1`'s Pade by substituting
`2^x - 1 = e^{x ln2} - 1`. That multiply's own rounding *was* the max:
the worst case sat at `x ~ -0.3991`, deep inside the Pade branch, which
is exactly why idea #7's seam-retune audit found the threshold could not
touch it -- it was auditing the wrong knob. Fitting `2^f - 1` in `f`
directly, with the leading `f*ln2` peeled into the closing fma, removes
the rescale entirely. Same lever as `tanpi`'s `pi` peel.

### Premise 2: the reduction has to floor, because `exp2`'s does

It does not, and for a *minus one* it must not. `k = floor(x)`,
`f in [0,1)` sends every small negative `x` to `f ~ 1`, where `F = 2^f - 1
~ 1` and `2^-1*(1+F) - 1` cancels catastrophically -- at `x = -1e-8`,
`f` rounds to `1.0` outright and the answer comes back **0**. That, not
the near-zero conditioning `expm1` has, is the real reason this function
carried a Pade branch. `k = round(x)`, `f in [-0.5, 0.5]` is continuous
through 0 and exact (`f` is a multiple of `ulp(x)` and smaller than it).

The measured cost of rounding rather than flooring is a bounded
amplification `2^k F/(2^k(1+F)-1)` reaching **1.414** at `|f| = 0.5`,
against floor's exactly 1.0 at `k = 0`. Worth noting because the
graveyard already rejected round-based `Q(f)` for `exp2`/`exp2_checked`/
`exp10`/`exp10_checked` -- correctly, since for those there is no
cancellation to fix and the 1.414 is pure loss. For `exp2m1` it buys a
whole Pade and a division. **The old entry's verdict does not transfer to
the `m1` members of that list, and this is why.**

### Numbers

| region | instrs | uOps | BlockRT | cyc/elem |
|---|---|---|---|---|
| `exp2m1_throughput` | 90 -> 59 | 103 -> 62 | 27 -> 17 | 1.843 -> **1.278** |

Latency: the old 80.00 was the branch artifact (arms 37.00/48.00); both
arms are 52.00 now, so the real move is 48.00 -> 52.00 on the arm that
carries essentially all inputs.

Two shapes were fitted and scored through the real f32 chain before
picking:

| shape | degree | wide max ulp |
|---|---|---|
| `F = f*Q(f)`, `Q` fitted on `[-0.5,0.5]` | 5 | 2.81 |
| `F = fma(f2, P(f), f*ln2)` (leading peeled) | 4 | **2.20** |

The peel wins on accuracy *and* is a degree lower. The clamp also
tightened, `[-151, 128)` -> `[-126, 128]`, so one exponent field at `k-1`
covers the range: `2^x - 1` is exactly `-1` for every `x < -25`, so the
low end had nothing to lose. "The guard is the licence", again.

### The sign-of-zero trap, hit for real

`fma(f2, P, f*ln2)` returns `+0.0` at `x = -0.0`: `f2` is `+0.0`, the
addend is `-0.0`, and `+0.0 + -0.0` is `+0.0`. The ulp sweep cannot see
this (zero-vs-zero scores 0); `edgecheck` caught it. The fix costs
nothing because the denormal arm this function needs anyway (the `k-1`
field halves the intermediate, same as `expm1`) can return the peeled
`f*ln2` term, which is already an operand of that fma and *does* carry
the sign. One select, three jobs.

## `norm_cdf`'s un-split `FRAC_1_SQRT_2`: measured at 5.9% of the binding max

The last un-split single-word irrational multiply in the crate, left open
above with the note that fixing it "needs `gelu`'s `RSQRT2_HI`/`RSQRT2_LO`
+ `erfc(z+dz) = erfc(z)*(1-2z*dz)` treatment, which is several ops rather
than one `fma`". It does. It is also worth almost nothing, and the reason
is a sensitivity that the earlier note read off the wrong function.

### The attribution, which had never been done for `norm_cdf`

`norm_cdf`'s doc comment claimed its 7 splits like `erfc`'s -- `exp`'s own
error and `erfcx_pos`'s evaluation, neither dominating -- but the standing
measurements only ever split `erfc` and `erfcx`. Done properly at
`norm_cdf`'s own worst point (`x = -1.8841501`, exhaustive max 7):

| term | ulp of the result | share |
|---|---|---|
| `erfcx_pos`'s polynomial evaluation | **2.931** | 42% |
| `exp_reduce!(-p)` | **<= 3.0** (pure factor) | 43% |
| `x/sqrt(2)` argument rounding | **0.415** | 5.9% |

So the doc's claim was right, and is now a number rather than an analogy.

**Why the argument term is so small, and why `gelu`'s reasoning does not
transfer.** `gelu` carries a two-word argument because it feeds `erfc`,
whose relative sensitivity `|d(ln erfc)/d(ln z)|` grows as `2z^2` -- ~170
half-ulps at `z = 9.2`. `norm_cdf` feeds `erfcx_pos`, and `erfcx` is
*flat*: `d(ln erfcx)/d(ln z)` is `-1` asymptotically and only **-0.729** at
the point that actually binds. That is the whole reason `norm_cdf` is
written against `erfcx` instead of `erfc` in the first place, and it is
exactly what makes the argument fix pointless here. **A two-word-argument
fix inherits its value from the callee's condition number, not from the
constant's own offset** -- `FRAC_1_SQRT_2`'s 1.711e-8 relative offset is
identical in both functions, and it is worth 170 half-ulps in one and 0.4
in the other.

Also worth stating: the constant's *bias* is not even the dominant part of
the argument term. `fl(xa * FRAC_1_SQRT_2)` carries the constant's fixed
-1.711e-8 relative offset **plus** the product's own rounding, up to
5.96e-8. Screened against exact `erfcx` and exact `exp` over 600k samples
per band, correctly rounding `z` moves that path's max from ~1.69 to ~1.40
ulp and its avg from 0.336 to 0.312 -- the bias is real and one-signed
(signed mean +0.14 to +0.20 ulp) but it sits inside a larger symmetric
rounding it cannot remove.

### The upper bound, built and measured

Not the cheap version -- the *ceiling*: `z` as a two-word pair with its
residual carried into `erfcx_pos`, i.e. the argument rounding removed
outright. `v = 1/(2+z)` in both of `erfcx_pos`'s branches, so
`dv/dz = -v^2` uniformly and one fma puts `dz` back:

    let zp = xa * RSQRT2_HI;
    let dz = fma(xa, RSQRT2_LO, fma(xa, RSQRT2_HI, -zp));
    let v  = fma(-dz, v0 * v0, vb);      // inside erfcx_pos

Gated on a `const CORRECT: bool` so `erfc`/`erfcx`/`norm_pdf` compile
identically -- confirmed, their mca regions came back byte-identical
(`erfc_throughput` 136/157/41, `erfcx_throughput` 130/148/40,
`norm_pdf_throughput` 76/83/25, all unchanged). Only `norm_cdf` pays:

| region | instrs | uOps | BlockRT | cyc/unit |
|---|---|---|---|---|
| norm_cdf_latency | 5210 -> 5571 | 5468 -> 5817 | 1312 -> 1440 | 68.985 -> 72.970 (+5.8%) |
| norm_cdf_throughput | 143 -> 152 | 172 -> 182 | 44 -> 48 | 3.436 -> 3.811 (+10.9%) |

Exhaustive over all 2^32: avg **0.0657 -> 0.0656**, max **7 -> 7**, and
the worst point does not even move (`x = -1.8841501` either way). That is
the entire return on removing the argument rounding *completely* -- 0.15%
of the average and none of the max, for +10.9% throughput. The cheap half
(a two-word constant with no residual carried, +1 fma) can only be a
fraction of that, so it is closed too. The 0.415 ulp the attribution
predicted is precisely what a 7 absorbs without moving.

**A `-0.0`-style edge trap on the way, worth recording**: the first build
returned non-finite for `x = +-inf`. `zp = inf * RSQRT2_HI` is `inf`, so
`fma(xa, RSQRT2_HI, -zp)` is `inf - inf = NaN` and the residual poisons the
result. `gelu` already guards this (its `np = max(-x, 0)` clamp is
documented as keeping `dz` finite at `x = +inf`); any *new* two-word
argument split needs the same guard, and the ulp sweep does catch it --
it reported `NON-FINITE RESULT in 2 of 4294967296 samples` rather than a
plausible-looking number.

### Transferable

- **Value a two-word argument split by the callee's condition number.**
  Same constant, same offset, 400x difference in what it is worth between
  `erfc` and `erfcx`. Compute `d(ln f)/d(ln z)` *at the point that
  actually binds* before writing any code.
- **A constant's one-signed bias is not automatically the dominant part of
  an argument's error.** Here the co-located product rounding is 3.5x
  larger and symmetric, so removing the bias alone recovers a minority of
  a minority.
- The remaining single-word irrational, `probit`'s `SQRT_2 * erfinv(...)`,
  is unaffected by this result: `erfinv` is steep near the ends, so its
  condition number is the opposite case and it still wants a harness row
  first.

### atan_poly's numerator peel: the "x + O(x^3) that already divides" class pays out again

`§tanh`'s small-arm table named a candidate class and left it open: *any
odd function whose series starts `x + O(x^3)` and that already divides*.
Only `dawson` had ever been checked against it (and failed). `atan_poly`
is a verbatim member -- `x * N(x^2) / D(x^2)` -- and had never been tried.
It pays:

    numer = fma(fma(fma(a2,x2,a1),x2,a0), x2, 1.0) * x     // before
    numer = fma(x*x2, fma(fma(a2,x2,a1),x2,a0), x)         // after

Same polynomial, same four operations. The unpeeled form spends one
full-weight rounding forming a value near 1 and a second on the multiply
by `x`; the peeled form folds both into one rounding at the result's own
scale. This is `tanh`'s row A -> row A2 (3.32 -> 2.46 max there), on a
different function.

Exhaustive (2^32) where available, 10M-sampled for the 2-arg rows:

| | avg before | avg after | max before | max after |
|---|---|---|---|---|
| `atan` | 0.0675 | **0.0627** | 4 | **3** |
| `atan_bounded` | 0.0526 | **0.0431** | 3 | 3 |
| `atand` | 0.0628 | **0.0582** | 4 | **3** |
| `atanpi` | 0.0775 | **0.0733** | 4 | 4 |
| `atan2` | 0.0683 | **0.0660** | 3 | 3 |
| `atan2_pos` | 0.0632 | **0.0621** | 3 | 3 |
| `atan2d` | 0.1072 | **0.1051** | 4 | 4 |
| `atan2pi` | 0.1139 | **0.1118** | 4 | **3** |

No regression on any axis of any caller. mca: instructions, uOps and
`Block RThroughput` are **byte-identical** on every region (`atan`
52/54/20, `atan2` 62/64/20, `atan_bounded` 37/39/10) while latency drops
3 cycles across the board -- `atan` 61.27 -> 58.27, `atan2` 67.19 ->
64.19, `atan_bounded` 37.00 -> 34.00. Same ops, one level shallower: the
numerator resolves a step earlier into the division that owns the
critical path. This is the mirror of the usual "removing an op need not
shorten the chain" -- here the chain shortens while the op count does
not move at all.

Two controls confirm the change is where it is claimed to be:
`atan_latency` (its own direct degree-17 fit, no `atan_poly`) and
`atanh` (never routes through `atan_poly`) are both bit-flat, 0.0516/3
and 0.0037/2 on each side.

Two notes for whoever reads this next:

- **The old localization of `atan`'s max was wrong.** This file said "the
  real worst case lives in the untouched `a<1` branch", on the strength
  of a quick-fuzz max of 3. Exhaustively the max is 4 and it sits at
  `x = 1.0220603`, i.e. *inside* the reflected `a >= 1` branch, just past
  the fold. Any future work on the `FRAC_PI_2 - y` fold -- including the
  pi hi/lo split rejected under idea #124 partly on that localization --
  should be re-screened against the exhaustive worst case, not the
  sampled one. (After the peel the worst case moves to `x = 0.9360248`,
  back inside `a < 1`.)
- **The 2-arg maxima above are sampled and were repeated 3x per side.**
  `atan2pi` is the only one whose max moved reproducibly (base 4/4/4,
  peel 3/3/3, no overlap). `atan2_unchecked` and `atan2_pos` each threw a
  single 4 on the *base* side across three runs while the peel side threw
  none -- suggestive, not established, so their readme maxima are left
  alone. Only the avg columns are published for the 2-arg rows.
## `exp10m1`: fit `10^d - 1` against the reduction's own residue, and the last Padé caller goes away

The fourth and last application of the `m1` treatment (`expm1`,
`exp_m1_over_x`, `exp2m1` preceded it), and the one that retires
`pade_expm1_ratio!` entirely.

The old form was two arms: a near-zero Padé fed `y = x*LN_10`, and a
direct arm that ran `exp10_checked`'s full reduction, built `2^k` through
`exp2_field_split`, evaluated `exp2_q_poly!` and *then* subtracted one.
The seam sat at `|x| < 0.2` rather than `exp2m1`'s `0.5` purely because
`LN_10` is 6x `LN_2` and the shared Padé is only fitted on `|v| <= 0.5`.

Carrying `D = 10^d - 1` through the combine instead removes the arm, the
division and the seam together: `10^x - 1 = 2^k*D + (2^k - 1)`, whose
sensitivity to `D` stays under 1.414 everywhere.

### Fit in `d`, not in `f` -- and this is where it differs from `exp2m1`

`exp10_checked`'s reduction already produces `x = k*log10(2) + d`, i.e.
`10^x = 2^k * 10^d`. So `10^d - 1` can be fitted against `d` **directly**,
and the `f = d*LOG2_10` multiply that the `2^f` form needs never happens.
That deletes a full-weight rounding *and* `fl(LOG2_10)`'s own **0.296 ulp**
representation error, against `fl(LN_10)`'s 0.134. `exp10_reduction!`'s
floor-adjust goes with it (a compare, a select and two adds): a centred
residue is what removes the near-zero branch, so nothing needs `[0,1)`.

Because `|d*ln10| <= ln2/2` is exactly `expm1`'s own `|r|` bound, the
approximant is `expm1_p_poly!` in rescaled coordinates -- same degree 4,
same Estrin fold.

### Degree 4 stops for a *different* reason than `expm1_p_poly!` does

`expm1_p_poly!` stops at degree 4 because the ulp-weighted LP's degree-5
coefficient converges to 0. Here degree 5 reaches the **same** LP margin,
0.2316 ulp. An oracle screen with an *exact* `ln10` peel separates the two
causes:

| peel constant | deg 3 | deg 4 | deg 5 |
|---|---|---|---|
| `fl(LN_10)` (shipped) | 6.301 | **0.2316** | 0.2316 |
| exact `ln10` (oracle) | 6.238 | 0.1830 | 0.0032 |

So at degree 5 the *entire* margin is the constant, and at degree 4 the
constant is worth only 0.049 of it. A two-word `LN_10` therefore buys
nothing without also going to degree 5 -- together ~0.23 ulp for +2 ops,
against an error budget whose dominant terms are the `dl` rounding and the
closing fma at ~0.5 ulp each. **Rounding-chain-dominated, not fit-limited**,
so it was not taken. Recorded because "two-word constant" pattern-matching
is exactly what `log10_normal` was rejected for the same week.

### Numbers (exhaustive, all 2^32)

| | avg ulp | max ulp |
|---|---|---|
| old (Padé + `exp2_q_poly`) | 0.1279 | 4 |
| new | **0.0946** | **3** |

| region | instrs | uOps | BlockRT | cyc/elem |
|---|---|---|---|---|
| `exp10m1_throughput` | 106 -> 64 | 116 -> 69 | 31 -> 18 | 2.629 -> **1.342** |

Better on max, avg and throughput at once, so no Pareto tier is warranted.
The residual 1 ulp against `exp2m1`'s max 2 is the peel constant and
nothing else: `fl(LN_10)` is 0.134 ulp off `ln(10)` where `fl(LN_2)` is
0.032 off `ln(2)`, and the two functions are otherwise the same object.

### The latency row is branch-shaped and *not* an artifact

`tools/mca_arms.py` reports 55.00 as-published against arms of 55.00 and
55.00. The near-zero arm returns the peeled `d*LN_10`, which is already an
operand of the other arm's closing fma, so both arms are the same
dependency chain and there is nothing for mca to mis-fuse. The old 112.00
*was* the artifact (arms 37.00/80.00), so readme's worst-offender list
loses an entry. Worth stating as a positive rule: a branch in a latency
region is a reason to *run the tool*, not a reason to assume the number is
wrong.

Against the arms it replaces: `|x| >= 0.2` gains 25 cycles, and the old
Padé arm's 37.00 becomes 55.00 -- the same near-zero latency cost the
other three `m1` rewrites paid, documented rather than split into a tier,
since throughput is the crate's stated axis and accuracy improves in both
regions.

### Clamp

`[-45.154503, 38.53184]` -> `[-37.0, 38.53184]`. The top is unchanged --
it is the exact float where `k` first reaches 128, so overflow still
saturates through the `k-1` field. The bottom tightens so one field at
`k-1` covers the range (`k >= -125`); `10^x - 1` is exactly `-1` for every
`x < -7.53`, so there was nothing to lose. "The guard is the licence",
again, and the same move `exp2m1` made.

### Signed zero, again

Same trap as `exp2m1` and it would have shipped the same bug: at `x = -0.0`
the residue `d` is `-0.0`, `d2` is `+0.0`, the poly's constant term is
positive, so `fma(+0.0, P, -0.0)` rounds to `+0.0`. The ulp sweep scores
zero-vs-zero as 0 and cannot see it. The denormal arm this function needs
anyway returns the peeled `d*LN_10`, which carries the sign for free --
one select, three jobs. `edgecheck`'s `exp10m1(-0)` pin is what proves it.

## `acos`: the single-fma combine is free, and an oracle screen says the max 4 is *not* in the combine

`asin`'s doc comment records that its big branch is "one `fma`, not a
multiply and a subtract", worth ~2 ulp of the result and an instruction
cheaper. `acos`'s negative arm still did the rejected thing --
`PI - fl(sqrt(1-a)*P(a))`, rounding the product on its own first. So the
transform transferred; the *reason* did not.

### What the oracle screen says (exhaustive, every f32 in `[-1,1]`)

Scoring the shipped combine with **both factors exact in f64**:

| | avg | max |
|---|---|---|
| exact factors, shipped combine | 0.0750 | **1** |
| exact everything | 0 | 0 |

So the combine contributes at most 1 ulp, and `acos`'s max 4 lives
entirely in `acos_poly` and the `sqrt`. Predicted consequence: the fma
cannot buy accuracy. Measured consequence, same sweep:

| variant | avg+ | max+ | avg- | max- | avg all | max all |
|---|---|---|---|---|---|---|
| shipped | 0.1048 | 3 | 0.1191 | 4 | 0.1119 | 4 |
| single-fma combine | 0.1048 | 3 | 0.1187 | 4 | **0.1118** | 4 |

Dead even, exactly as the screen predicted. **It shipped anyway, because
it is cheaper**: 47 -> 43 instrs, 50 -> 45 uOps, Block RThroughput flat at
12.00, 0.820 -> 0.771 cyc/elem, latency 39.99 -> 34.99 (both figures
honest -- `mca_arms.py` gives both arms equal to the published number on
both sides). It also deletes the `mulsign` + `x + 0.0` signed-zero dance:
with one value-based `x < 0.0` compare choosing both the sign of `s` and
the addend, there is no bit-based sign left to disagree at `-0.0`.

Worth stating as a rule: **an oracle screen that says "no accuracy here"
is not a reason to drop the transform.** Price it on the other axis.

### Rejected: a two-word `PI` on the negative arm

`fl(PI)` sits 0.367 ulp above `pi`, so a low word looks like free avg.
It is a **large regression**: avg- 0.1191 -> 0.9352, ~8x worse.

The mechanism is worth recording because it is not obvious.
`fl32(pi) == 2 * fl32(pi/2)` *exactly* -- same mantissa, exponent one
higher. `acos_poly`'s constant term is `fl32(pi/2)`, so near `x = 0` the
negative arm computes `fl(PI) - fl(pi/2) = fl(pi/2)` with no error at all,
and the whole poly was coordinate-descended against that consistency.
Adding a low word to `PI` alone breaks the pairing that the poly was
tuned to. A two-word constant is only free when nothing downstream has
already been fitted against the one-word value.

### Measured but not taken: a dedicated small branch (the real lever)

`acos(x) = pi/2 - asin(x)` below `|x| < 0.5`, reusing `asin_small`'s odd
poly, removes the `sqrt` and its inexact `1.0 - a` from the region where
Sterbenz does not apply -- the same move that took `asin` from max 5 to 2.

| variant | avg+ | max+ | avg- | max- | avg all | max all |
|---|---|---|---|---|---|---|
| single-fma combine | 0.1048 | 3 | 0.1187 | 4 | 0.1118 | 4 |
| + small branch | 0.0861 | 3 | 0.0717 | **2** | **0.0789** | **3** |

avg **-29%** and max **4 -> 3**, with the negative arm reaching max 2. The
residual max 3 then sits at `x = 0.5404`, just above the crossover, i.e.
in the big branch -- so the crossover would want its own tune afterwards,
as `asin`'s did.

Not taken here only because it is ~9 instructions on a 43-instruction
region (a whole second poly, evaluated unconditionally for the branchless
select) against a crate whose stated axis is throughput, and the free
combine win above was ready to land on its own. This is a **priced,
ready-to-implement lead**, not a rejection: the numbers above are the real
exhaustive sweep, and adding a two-word `pi/2` to the small arm on top
measured *identical* (0.0789/3), so that part is not needed.
## `cargo build --example mca_target` does not re-emit asm -- a stale `.s` serves old numbers

Cost time twice on 2026-08-02, in two different worktrees, so it is worth
its own entry. `examples/mca.rs` regenerates the assembly with

```sh
touch examples/mca_target.rs
cargo rustc --release --example mca_target -- --emit=asm -C debuginfo=0
```

and **both halves are load bearing.** Cargo's fingerprint does not record
the `--emit=asm` passed as a raw rustc arg, so if `mca_target.rs` itself
has not changed, `cargo rustc` treats the request as a no-op and leaves
the previous `.s` in place; the `touch` is what forces it. A plain
`cargo build --release --example mca_target` never emits assembly at all,
no matter what changed -- it rebuilds the *library*, reports a normal
recompile, and takes long enough to look like it worked.

The failure is silent and reads as a *result*: `mca_region.py` happily
scores the stale file and reports that a change costs exactly nothing.
Here that produced "candidate is byte-identical to baseline on all five
rungs" for an edit that demonstrably changed the function's output.

Two tells, both cheap:

- `ls -lt target/release/examples/mca_target-*.s` -- if the timestamp
  predates the edit, every number from it is the old code's. A single
  `.s` whose mtime is older than your last build is the whole diagnosis.
- The crate's standing rule in the other direction: a change that *must*
  move numbers moving none. Byte-identical asm means "not real" only once
  you have confirmed the asm was actually regenerated.

Both traps have the same root -- the `.s` is a build *side effect* that
nothing in the normal build graph depends on -- so neither `cargo build`
nor a fresh `cargo run --example mca` of the wrong shape will fix it.

## `atanpi`: fold `1/pi` *before* the quadrant reflection, so its constant is an exact `0.5`

`atanpi` was `atan(x)` followed by the two-word `1/pi` multiply, and the
entry above called what remained "essentially just `atan`'s own error plus
binade shift". That was true of the `|x| < 1` arm and wrong about the
other one, which carries a bias `atan` cannot avoid and `atanpi` can.

`atan` folds `|x| >= 1` with `FRAC_PI_2 - t`. **`fl(pi/2)` sits 0.367 ulp
*above* `pi/2`**, and that is an absolute offset, so it survives the later
scaling intact: `(fl(pi/2) - pi/2)/pi = 1.391e-8`, against an ulp of
`2^-25` for a result in `[0.25, 0.5)` -- which is exactly where every
`|x| >= 1` input lands. A fixed **+0.47 ulp**, on half the domain, that no
amount of work on `atan_poly` can reach.

Scaling *before* the fold moves the reflection into half-turns, where its
constant is an exactly representable `0.5`, and replaces both that bias
and the `FRAC_PI_2 - t` subtraction's own rounding with one exact
constant. `0.5 - h` for `h <= 0.25` then rounds once, at the result's own
magnitude.

Exhaustive over every positive f32 pattern (the function is odd, so the
negative half mirrors bit for bit), scored with `accuracy.rs`'s own
integer `ulp_diff`:

| branch | before | after |
|---|---|---|
| `\|x\| < 1` | 0.0561 avg / max 4 | **bit-identical** |
| `\|x\| >= 1` | 0.0909 avg / max 2 | **0.0059 / max 2** |
| whole domain | 0.0733 / 4 | **0.0308 / 4** |

All rows are exhaustive, scored with `accuracy.rs`'s own integer
`ulp_diff` on its scale: mirror the negative half of an exactly odd
function, divide by all `2^32` patterns so the ~0.39% that are NaN dilute
the mean as the harness lets them. **Both formulas were measured in one
binary, each defined locally, so neither depends on what `src/lib.rs`
contained at the time** -- and the old one reproduces the published
0.0733 / max 4 in that same run, which establishes the scale rather than
assuming it. Both also report the *same* worst input, `x = 9.3035436e-1`:
direct evidence the max is untouched and sits on the arm this leaves
bit-identical.

Measuring it that way was not fussiness. A first attempt at this A/B
scored the two sides against *different* `atan_poly`s, because another
worktree landed a peel of that polynomial midway through the sweep; the
baseline moved 0.0775 -> 0.0733 under it, and publishing against the
stale figure would have credited this change with that one's average.
A shared dependency going stale mid-measurement is the live form of this
file's standing rule about re-verifying a baseline that belongs to
somebody else's function.

**The max does not move, and this does not claim it.** `atanpi`'s max 4
lives on the `|x| < 1` arm, which this leaves bit-identical -- it is
`atan_poly`'s own, at `x ~ 0.930`. The `|x| >= 1` arm was already max 2
and stays there. What moves is the average, 2.4x.

Cost, `mca_region.py`: 58 -> 60 instrs, 60 -> 62 uOps, **Block
RThroughput unchanged at 20**, throughput 1.643 -> 1.724 cyc/elem
(+4.9%), latency 68.079 -> 66.189 (-2.8%). Blast radius 2 of 313 regions,
exactly `atanpi`'s own. All 7 edgecheck pins pass, `-0.0` and `+-inf`
included: both arms are computed from `|x|` and the sign is reapplied
last, so signed zero rides out on the closing `mulsign` rather than
depending on a coefficient's sign the way `tand`/`tanpi` do.

### The +2 instructions are a domain boundary, deliberately

The reduction is `atan`'s own `min(a, 1/a)`, and the value handed to
[`atan_bounded`] is already non-negative, so that function's internal
`abs` and `mulsign` are dead work. LLVM folds the `abs` only if the value
is *visibly* sign-cleared -- `.abs()` at the call site is what makes it
visible, and is worth 2 of the 4 instructions it otherwise costs (62 ->
60). A bit-mask spelling of the same thing measures identically. The last
`mulsign` does not fold.

Calling the private `atan_poly` directly removes the other two and lands
on **58 instrs / 60 uOps, dead level with the baseline** -- but it makes
`atanpi` share private code with `atan`, which merges `atan`'s ownership
domain into `asinpi`'s (`jm_domains.py` joins public functions through
shared *private* units; a call to a public one is a contract, not a
merge). Measured and rejected on that ground rather than on cost: `atan`
was claimed by another worktree and mid-edit in `atan_poly` at the time.
`atan_bounded`'s documented contract is precisely this caller -- "callers
who already know their input is bounded, e.g. already reduced via some
other identity" -- so the public route is also the designed one.

### Transferable

- **A composite's floor is its callee's error only where it shares the
  callee's branches.** The `|x| >= 1` arm of `atan` computes something
  `atanpi` does not want: an angle in radians whose reflection constant
  is irrational. Rescaling a composite's *output* can never undo a
  constant its *callee* already folded in. Look for a reflection,
  quadrant fold or offset inside the callee whose constant is exact in
  the caller's units -- `pi/2 -> 0.5`, `pi -> 1`, `pi/4 -> 0.25`.
- **"Essentially just the callee's error" is a claim about one branch
  until it is measured per branch.** Splitting the sweep by the callee's
  own branch condition was the entire diagnosis here and cost one probe.

## `atan2pi`: the same fold, and the "~2x floor" this file recorded was wrong

Immediately after `atanpi`. That entry moved `atan`'s `|x| >= 1` quadrant
fold into half-turns, where its constant is an exact `0.5` instead of
`fl(pi/2)`. `atan2pi` had the identical defect twice over and was closed
in this file on a floor estimate that turns out not to be a floor.

**The recorded claim, quoted:** "`atan2`'s result spans `[-pi, pi]`, and
dividing by pi shifts `~pi -> ~1` and `~pi/2 -> ~0.5`, each a binade step
that *doubles* the error measured in ulp. So `atan2pi`'s floor is about
2x `atan2`'s 0.0681, i.e. ~0.136, and it now measures 0.1137. Nothing
further to take here without changing `atan2` itself."

Two things are wrong with it. The arithmetic first: dividing by pi is not
a flat 2x in ulp terms, it is **0.64x-1.27x depending on where in its
binade the result lands** -- `pi -> 1` crosses a binade boundary and
*halves* the relative ulp, `3.0 -> 0.9549` does not and tightens it by
1.27x. So the estimate was high. The larger error is the second one: it
priced only the *rescaling*, and never asked whether the constants being
rescaled had to be inexact at all. They did not.

`atan2` folds quadrants with `FRAC_PI_2 - mulsign(FRAC_PI_2, x)`, giving
`0`, `fl(pi)` or `fl(pi/2)`, and its both-infinite convention is
`fl(pi/4)`/`3*fl(pi/4)`. In half-turns every one of those is exactly
representable -- `0`, `1`, `0.5`, `0.25`, `0.75` -- and
`0.5 - mulsign(0.5, x)` is exact for both signs of `x`. `fl(pi)` sits
0.367 ulp above `pi` and the `x < 0` fold adds it as an **absolute**
offset, so it survives the later `1/pi` scaling as a fixed relative bias
across the whole `x < 0` half-plane. Rescaling afterwards cannot remove
it; folding first means it is never introduced.

Measured with an out-of-tree probe, both formulas compiled into one
binary against one reference, so neither side depends on what `src/lib.rs`
held at the time. 40M pairs, uniform random bit patterns for both args:

```
branch      old avg  old max    new avg  new max    share
x > 0        0.1583        3     0.1177        3   49.61%
x < 0        0.0673        2     0.0237        1   49.61%
x == 0       0.0000        0     0.0000        0    0.00%
nonfinite    0.0000        0     0.0000        0    0.78%
WHOLE        0.1119        3     0.0702        3
```

The old side reproduces this file's published 0.1119 in that same run,
which is what makes the new column readable. `accuracy.rs` then confirms
it end to end: **0.112 -> 0.0702 avg**, three consecutive runs giving
0.0702 / 0.0701 / 0.0702, max 3 both sides. The whole-plane figure agrees
with the probe to four digits across two different references (sleef
`atan2_u35` vs std's f64 `atan2`) and two different sample sets. Sibling
rows in the same runs are controls and do not move: `atan2` 0.066,
`atan2_pos` 0.062, `atan2d` 0.105, each matching its published value.

**Both arms improve, for two different reasons.** `x < 0` is the `fl(pi)`
bias, and it is the bigger relative win (2.8x, and the only max that
moves, 2 -> 1). `x > 0` never touches `fl(pi)` at all -- it improves
because the core is now `atanpi` rather than `atan`, so the `|y/x| >= 1`
reflection inside it carries an exact `0.5` instead of `fl(pi/2)`. That
is `atanpi`'s own win arriving through the call.

**The max does not move and this does not claim it**: 3 on both sides.

**What the `x > 0` arm is now limited by**, decomposed on 20M samples of
that arm, scoring `atanpi` against `atan` of the *already-rounded*
quotient to separate the two:

```
new atan2pi total          avg 0.1177  max 3
y/x rounding alone         avg 0.1011  max 1
atanpi(q), q as given      avg 0.0289  max 3
old atan2pi total          avg 0.1582  max 3
```

86% of what remains on that arm is the single `y/x` division's own
rounding, which no rearrangement of the quadrant fold can reach. The
`x < 0` arm is far cleaner (0.0237) because adding an exact `+-1.0`
puts the result in `[0.5, 1]` while `atanpi(q)` may be tiny, so the
quotient's relative error is damped rather than exposed. A two-word
division is the only lever left and is not obviously worth it.

Cost, from `tools/mca_region.py`, all counting rungs agreeing so mca's
cycle column needs no arbitration: **68 -> 70 instrs, 70 -> 72 uOps,
Block RThroughput unchanged at 20.00**, throughput 1.867 -> 1.907
cyc/unit (+2.1%), latency 72.188 -> 72.111 (flat). Two instructions, the
same price `atanpi` and `atand` each paid for this class of fix. Blast
radius 2 of 313 asm regions, exactly `atan2pi`'s own.

**`atan2pi` now sits at 0.0702 against `atan2`'s own 0.0659** -- so the
real floor was near `atan2`'s number all along, not 2x it. The gap that
is left is the binade attenuation the old entry was reasoning about,
which is real but is a factor under 1.3, not 2.

`atan2pi` is now the third function carrying `atan2`'s full
`nonzerox || bothzero` skeleton -- `atan2` and `atan2_latency` are the
others, `atan2_unchecked` carries a deliberately reduced one and
`atan2_pos` is a plain composite. So it owns its own copy of `atan2`'s
two documented edge-case fixes: the
`nonzerox` select that keeps `atan2pi(-0.0, +0.0)` at `-0.0`, and the
trailing `y.is_nan()` override. Neither is reachable by the ulp sweep --
the two zeros compare equal, and a NaN degrading to a finite `+-0.5` is
just one sample. `edgecheck` gains 14 `atan2pi` pins covering both, on
top of the 4 it had; an edit to `atan2pi` alone would otherwise not be
caught by `atan2`'s pins at all. All 18 pass.

Reusable, and it is the third time today the same shape has paid:

- **A floor estimate is not a measurement.** This one was arrived at by
  multiplying another function's number by a factor that was itself
  wrong, and it closed the idea for a whole session.
- **Price the *rescaling*, then ask whether the thing being rescaled had
  to be inexact.** The old entry did the first and stopped. Every
  quadrant constant here was exact in the target units.
- **Splitting by branch is what makes it visible.** `x > 0` and `x < 0`
  are half the plane each, improve for unrelated reasons, and the
  whole-plane average alone would have read as one undifferentiated 37%.
## `exp_reduce!`: degree 6 and a peeled `1 + r`, at one arithmetic instruction *fewer* -- shipped

`exp_r_poly!` is core-locked and shared by ~20 public functions.
`exp_reduce!` is not: its five call sites -- `exp_checked`, `erfc`,
`erfcx`, `norm_cdf`, `norm_pdf` -- are exactly one `jm` domain. Those
five are also the crate's most exponential-bound functions. So the two
macros can carry different polynomials, and the cost of a better one
lands only on the callers that wanted it.

A degree-6 `exp_r_poly` had been measured before and rejected at +5-23%
throughput. That rejection priced the bump against `exp`/`expm1`/
`sigmoid`, which are at 2-3 ulp with nothing downstream amplifying them.
It never priced it against the erfc family, and it never combined it with
a peel -- which is what makes it free.

**Attribution first.** Scoring the shipped f32 evaluation order over the
whole reduction range `|r| <= ln2/2`, against `e^r`:

| poly | max | avg |
|---|---|---|
| shipped degree 5, `l0 = r + 1` | 2.392 | 0.7165 |
| shipped degree 5, peeled | 2.148 | 0.6783 |
| degree 6, peeled | **0.830** | **0.2654** |

Degree 5 is genuinely spent: a full f32 minimax refit of the shipped
shape moves its max 2.392 -> 2.352, i.e. 1.7%. The fit, not the
evaluation, is the binding term at degree 5, and it is ~2.0
ulp-equivalent against degree 6's ~0.15.

**The peel is what pays for the degree.** Evaluating `e^r - 1` and
folding the `+ 1` into the reconstruction as `fma(s, t1, t1)` is exact:
`t1` is a power of two out of `exp2_field_split`, so that fma *is*
`t1 * fl(1 + s)`, the same single rounding the multiply it replaces
already paid. It removes `fl(1 + r)`'s rounding, which entered at
`(1 + r)/e^r ~ 0.92` of full weight, and replaces it with one attenuated
by `|e^r - 1|/e^r <= 0.415`. Opcode histogram, per call:

    +1 vfmadd   (the degree-6 term)
    -1 vaddps   (`l0 = r + 1`, gone)
    -1 vmulps   (`p*t1`, absorbed into the fma)

Net **one arithmetic instruction fewer**, and `vbroadcastss` does not
move -- the new coefficient costs no broadcast. Total instruction count
is +1 (one `vmovaps` of register pressure) and `Block RThroughput` is
**unchanged in all five regions** (41/40/44/25/20 both sides).

**Exhaustive result** (every f32 bit pattern; `erfc`/`erfcx` over their
documented reference range):

| | max before | max after | avg before | avg after |
|---|---|---|---|---|
| `exp_checked` | 3 | **1** | 0.0370 | 0.0042 |
| `norm_pdf` | 4 | **2** | 0.0270 | 0.0178 |
| `erfcx` | 6 | **4** | 0.1439 | 0.1333 |
| `erfc` | 7 | **6** | 0.1289 | 0.1217 |
| `norm_cdf` | 7 | **6** | 0.0662 | 0.0622 |
| `erf` (control) | 3 | 3 | 0.0270 | 0.0270 |

`exp_checked` is faithfully rounded. `erf` and `erfcx (x >= 20)` are
bit-identical, which is the control: neither routes through
`exp_reduce!`. The cost is one level of fma latency, +2 to +4 cycles
(`exp_checked` 50.00 -> 54.00).

### The rejected grouping, and why the tradeoff is real

Degree 6 can be arranged two ways, and they are genuinely
Pareto-incomparable -- one arithmetic op against one dependency level:

- **Folded** (ships): top coefficient enters at the `r^2` level, `r^4`
  never formed. Serialises the top pair behind the bottom one.
  `Block RThroughput` unchanged; +1 fma level of latency.
- **Distributed** (rejected): `[r + r^2*(c0 + c1 r)] + r^4*[(c2 + c3 r)
  + c4 r^2]`, halves independent, meeting in the closing fma. Latency
  *exactly* neutral -- `exp_checked_latency` 50.001 both sides, which is
  what a depth-3 `s` predicts. But forming `r^4` raises `Block
  RThroughput` by one in **every** region (41->42, 40->41, 44->45,
  25->26, 20->21), with instrs and uOps moving the same way. Rejected:
  throughput is the metric this crate prioritizes, and the folded form's
  apparent throughput cost is *not* confirmed by the bottleneck metric
  while this one's is.

Degree 4 is the reason there is a choice at all: a depth-2 `m` can only
reach degree 3 (`A + r^2*B` with `A`, `B` degree <= 1), so degree 6 needs
either `r^4` or one more level. There is no third option.

The two are worth the **same accuracy downstream** -- 0.830 vs 1.032 ulp
in the isolated poly, but `erfc` 5.252 and `erfcx` 2.980 in simulation
for *both*. That is the lesson: the 0.2 ulp between them sits far under
`erfcx_pos`'s own ~2.9, so it never reaches the output. Price a poly
variant at the point that binds, not in isolation.

### Transferable

- **A shared kernel's rejection does not bind a private one.** The
  degree-6 result was correct and its cost was real; it was measured
  against the wrong caller set. Check whether the accuracy-wanting
  callers are a *subset* reachable behind a narrower lock before
  accepting a shared-poly cost verdict.
- **A leading-term peel can pay for a degree bump.** The peel hands back
  an add and converts a multiply to an fma; that is the whole budget for
  one more coefficient. Neither change alone is as good as the pair --
  degree 6 alone had been priced at +5-23%, and the peel alone is worth
  only ~0.25 ulp.
- **`fma(s, t1, t1)` where `t1` is a power of two is exactly
  `t1 * fl(1 + s)`.** Reconstruction multiplies by an exponent field are
  free places to put a `+ 1`.
- The simulation predicted `erfcx`'s worst point at `x = -0.514` where
  the exhaustive sweep found `-0.505`, and predicted the post-change
  worst point migrating to near zero, where it landed at `3.9e-4`. A
  numpy model of the exact f32 evaluation order is worth building before
  a multi-hour exhaustive sweep.
