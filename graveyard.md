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
needed. `asinpi` is the same shape and equally not applicable, so this
closes the whole half-turn family.

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
