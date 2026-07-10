# Ideas tried and rejected

Budget: ≤0.5 avg ulp, ≤2 max ulp. Ideas tried and reverted/not adopted, so
they aren't re-attempted without new information. Adopted changes live in
git history / readme.md, not here. Untested backlog is at the bottom.

## Cross-cutting

- **Degree-reduction probes (2026-07-07)**: lolremez rejected all 3 —
  `log_2` 9→8 (max ulp 3-5 vs 2 cap), `exp2` Q 5→4 (est. max rel error 40x
  worse), `sinf_poly` 9→7 (178x worse). Not implemented.

- **accuracy.rs coverage audit (2026-07-10, resolved -- no gap)**: cross-
  referenced every `pub fn` in `src/lib.rs` against `examples/accuracy.rs`
  to check for a silently-unmeasured public function (the same class of
  gap the `edgecheck.rs`/`codegen_check.rs` hardening passes closed for
  other tools). 12 functions have zero references there
  (`cbrt_accurate_normal`, `cbrt_approx`, `cbrt_constant`, `exp2_approx`,
  `ln_normal`, `log10_normal`, `log2_approx`, `log_2_normal`, `pown_const`,
  `rcp_approx`, `rsqrt_approx`, `sqrt_approx`) -- all confirmed benign:
  `*_normal` are shared internal cores whose public wrappers (`ln`,
  `log10`, `log_2`, `log_2_unchecked`, etc.) are already measured;
  `pown_const` is a compile-time-specialized copy of `pown`'s own already-
  measured algorithm, no new numerics to check; the `*_approx`/
  `cbrt_constant` family are deliberately-rough exploratory functions kept
  only for the `*_approx_plot` test suite/PNG generation (already
  documented elsewhere, e.g. `rsqrt`'s own doc comment, as "not a
  candidate replacement," never meant to carry the crate's real accuracy
  guarantee). No code changes.

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

- **sinf_poly: decoupled per-caller copy for sin_checked/cos_checked
  (paper-screened 2026-07-09, not implemented)**: tried the asin_poly
  playbook (full decoupling instead of joint reweighting) on the theory
  that the LP rejection above was a joint-fit constraint, not a real
  no-headroom result. numpy screen replicating each caller's *actual*
  reduction (first attempt used the wrong formula for cos_checked --
  `q = round(x/pi-0.5)` directly instead of the real `k=round(x/pi-0.5);
  q=k+0.5`, which looked like a catastrophic 3.3e-3 error at x≈0 before
  the bug was found and fixed) shows both callers land on the *exact same*
  worst-case r once the reduction is correct: max abs err 2.169e-08,
  identical for both, since both ultimately draw `r` from the same
  distribution regardless of caller. Decoupled per-caller least-squares
  refit only bought ~1.25x (sin_checked) to ~3x (cos_checked) tighter
  continuous-math error, avg ~2x -- real, but at the "isolated LP
  predictions under ~10-15% has repeatedly turned out not worth the
  round-trip" magnitude already established elsewhere in this file
  (atan_latency's own LP refit), compounded by sinf_poly already being
  documented near f32's precision floor (2026-07-07 entry above). Not
  implemented; the earlier LP rejection's "opposite ends of the domain"
  framing was directionally right about *why the joint reweighting
  failed* but doesn't imply large decoupled headroom the way asin/
  acos_poly's *domain-restriction* (asin only needing `[0.25,1)` vs
  acos_poly's `[0,1)`) did -- sinf_poly's callers don't restrict `r`'s
  own range at all, they just weight it differently, a structurally
  smaller opportunity.

- **exp_pos_neg: decoupled per-caller e/o copy for sinh/cosh
  (paper-screened 2026-07-09, not implemented)**: same playbook, ruled
  out even faster. sinh's `|x|<0.5` restriction is on the *input* `x`,
  which only selects which integer `k` the Cody-Waite reduction picks --
  it does not narrow the reduced residual `r`'s own range, which stays
  `[-ln2/2, ln2/2]` identically for sinh and cosh regardless of caller.
  Confirmed numerically: current `e`/`o` coefficients give the *exact
  same* max error (1.325e-07) evaluated as either `e+r*o` (cosh/sinh's
  `ep`) or `e-r*o` (`en`) across the shared `r` range -- no domain
  difference exists to exploit at all, unlike acos_poly/asin's genuine
  restricted-vs-full domain split. Not implemented; ruled out on paper
  before writing any Rust or running any tuner.

- **expm1 Pade degree bump 3→5 (2026-07-08)**: scipy-seeded refit found
  real avg-ulp headroom in the near-zero branch (0.138→0.135), but the
  function's actual max ulp (6) lives in expm1's other branch (`exp(x)-1`),
  untouched. Cost real throughput (+7%). Reverted.

- **expm1 direct branch round-off budget audit (idea #7, 2026-07-10,
  screened -- no actionable single term found)**: root-caused the direct
  branch's own max-ulp-6 worst point (`x=0.9652361`) the same way as
  asin_small above -- computed the identical formula in f64 to split
  rounding from truncation/coefficient error. Unlike asin (7:1 truncation-
  dominated), this one is roughly balanced: rounding 2.946 ulp, truncation/
  coefficient error 2.579 ulp, neither wildly dominant. Drilled into the
  rounding half specifically (per-step f32-vs-f64 comparison): `k`/`r`
  (the Cody-Waite reduction) are already near-exact (0.000/-0.016 ulp);
  essentially all 2.946 ulp of rounding traces to the polynomial
  evaluation of `p` itself (1.473 ulp in `p`'s own scale, doubled to 2.946
  in the final result purely because this worst point's `t2=2` exponent
  scaling amplifies it 2x -- not a separate error source). This is
  compounding rounding across a fixed, already-minimal-length fma chain
  (3 fma's + 2 multiplies + 1 add) evaluating a degree-8-ish poly in f32,
  not a single avoidable extra rounding step the way `exp`'s own
  Cody-Waite fix (this technique's proof point) removed a genuinely
  redundant single-word reduction -- no restructuring found that would
  reduce this without either adding ops (failing "no perf penalty", same
  shape as asin_small's own rejection just above) or fundamentally more
  precision (a Df32-style accurate tier, out of scope for a "no perf
  penalty" fix to the default function). Not implemented; no code
  changes. *The round-off budget technique doesn't always find a clean
  single dominant term to attack -- sometimes the honest answer is
  "many small roundings, roughly tied with truncation, nothing cheap to
  cut."*

- **sinh round-off budget audit (idea #7, 2026-07-10, screened -- no
  actionable single term found)**: same technique, applied to sinh's own
  max-ulp-5 worst point (`x=-1.8935952`, direct `exp_pos_neg`-based
  branch). Rounding dominates over truncation here (3.055 ulp vs 1.670
  ulp) -- but drilling in found the rounding isn't concentrated in one
  avoidable step either: `p_neg` (the `exp(-x)` polynomial combine)
  carries a modest 0.643 ulp of its own rounding error, which then gets
  amplified ~8x by the necessary exponent-field reconstruction
  (`t1n*t2n=2^3` at this `x`) purely because that's how relative error in
  a value multiplied by a large power-of-two scale factor translates to
  absolute-ulp terms at the *final* result's (much smaller) magnitude --
  not a separate, fixable error source, just inherent amplification
  through the reconstruction every `exp_pos_neg` caller already relies on.
  No single step stood out as cheaply improvable. Not implemented; no code
  changes. *A second confirmation (after expm1) that this technique
  doesn't always surface a clean target -- when the dominant contributor
  turns out to be "a normal amount of poly rounding, amplified by an
  unavoidable exponent scale," there's nothing left to cut without adding
  real precision (cost) somewhere.*

- **Range-invariant sweep for every function with a known mathematical
  output bound (2026-07-10, resolved -- confirms the sin_checked/
  cos_checked fix above was isolated, not a wider pattern)**: generalized
  that fix's own discovery method -- checked `tanh` (`[-1,1]`), `sigmoid`
  (`[0,1]`), `erf` (`[-1,1]`), `atan`/`asin` (`[-pi/2,pi/2]`), `acos`
  (`[0,pi]`), `atan2` (`[-pi,pi]`), `hypot` (`>=0`), plus `sinpi`/`cospi`/
  `sind`/`cosd`/`sin_checked`/`cos_checked` themselves, all at
  `{0,-0,+-1,+-1e6,...,+-1e37,+-f32::MAX,+-inf}`. Every one of them stayed
  correctly bounded except the already-known/already-fixed cases: plain
  `sin`/`cos` (unchecked, expected garbage outside their documented
  domain) and `sind`/`cosd` (whose own doc comment already explicitly
  disclaims correctness -- only finiteness -- past their `~4.7e7` limit,
  and `2.6e21` is technically still finite). No new bugs found; confirms
  the double-float-reduction class of bug is specific to `sin_checked`/
  `cos_checked`'s own `two_prod`/`two_sum`-based reduction (grep confirms
  no other function in the crate uses `two_prod`/`two_sum` at all, and no
  other function uses `POLY_SAFE_BOUND`), not a pattern requiring a
  broader sweep of the rest of the crate.

- **tanh direct rational P(x²)/Q(x²) over [0,~9] (2026-07-08)**: needs 13
  free coefficients to converge over the full domain — far more than any
  poly in the crate. A 2-domain split needs 14 total, likely more work than
  the current expm1-based formula. Not implemented.

- **exp2/exp2_checked/exp10/exp10_checked/exp2m1/exp2_checked_df:
  round-based Q(f) reduction (2026-07-09/10, rejected everywhere it was
  tried)**: backlog round 3 #1's core premise -- refit `Q(f)=(2^f-1)/f`
  over the centered `f∈[-0.5,0.5]` (k=round(x)) instead of `f∈[0,1)`
  (k=floor(x)) -- did land a real, measurable win on `tune.rs`'s own coarse
  grid (degree 5: max ulp 2→1, avg 0.20307→0.04845) and the degree-4 half
  of the idea was cleanly falsified first (max ulp 13, not competitive, no
  Rust round-trip needed). But every one of the 6 functions sharing this
  poly failed for a *different* reason once actually wired in and verified
  against the real 100M-sample fuzz, `edgecheck.rs`, and `mca` -- not one
  survived:
  - **exp2, exp10 (unchecked, single-exponent-field, no k1/k2 split)**: a
    genuine new correctness bug, not an accuracy tradeoff. `round(x)` can
    land `k=128` for `x` strictly inside the documented `[-126,128)`
    domain (e.g. `x=127.6`), which this construction can't represent --
    produces `NaN` inside the promised-safe range. Caught immediately by
    the real fuzz (`ulp_diff`'s NaN-mismatch sentinel), not by `tune.rs`'s
    own coarse grid (which only ever samples floor-reachable `x`).
  - **exp2_checked**: no structural saving available (never had a
    floor-adjust to remove), so the swap was pure cost -- `x.floor()`
    (one `vroundps`, an otherwise-idle port on this CPU) became an
    explicit `(xs+M)-M` add/sub pair that instead contends with the
    poly's own already-saturated fma/add ports. Confirmed via assembly
    diff (64→66 total instructions, `vroundps` 2→0, `vaddps` 4→8, fma
    count unchanged at 14) and `mca` (throughput 1.399→**3.007** cyc/elem,
    more than double, latency flat). Also a real max-ulp regression on the
    full fuzz (1→2) that the coarse grid didn't predict (grid showed
    max 1→1).
  - **exp2m1**: same coefficient/reduction swap, same port-contention
    mechanism, smaller but still real cost (latency +1 cyc, throughput
    +3.3%) and a real max-ulp regression (3→6/7, avg roughly flat).
  - **exp2_checked_df**: the one case with a genuine accuracy win that
    *did* survive real verification -- unlike `exp2_checked` itself,
    `powf_checked`'s overall error is dominated by `log2_df`/the
    downstream `y`-multiply, not this poly's own fit quality, so the
    tighter centered fit helped rather than getting swamped (`powf_checked`
    avg ulp 0.0229→0.0118, max 123→118; `powf_checked_unchecked` avg
    0.0459→0.0243, max 96→89). But `mca` showed the same inherited
    port-contention cost, diluted but still real: `powf_checked`
    9.105→9.418 cyc/elem (+3.4%), `powf_checked_unchecked`
    7.234→7.475 (+3.3%). Real accuracy gain *with* a real perf penalty --
    doesn't clear this loop's bar.
  - **exp10_checked**: looked like the one clean win -- real structural
    saving (deletes the floor-adjust select + 2 adds, since `exp10`'s own
    reduction is natively round-based already), confirmed by both fuzz
    (avg ulp 0.0343→0.0107, max 1→2, still inside the crate's ≤2 budget)
    and `mca` (latency 72.06→**51.06** cyc, -29%; throughput
    2.736→**1.903**, -30%). But `edgecheck.rs` caught a real correctness
    bug the fuzz never sampled: at the overflow-saturation boundary (`k`
    clamped down to exactly `128`, reached via `exp10_checked(inf)` after
    its own `x.clamp(-1000,1000)`), `t1*t2=2^128` sits right at
    `f32::MAX`. The old floor convention guarantees `f>=0` there, so the
    poly's `2^f>=1` factor only ever pushes the product *up* into
    overflow; the round convention allows `f<0`, letting `2^f<1` pull the
    product just *under* `f32::MAX` instead -- `exp10_checked(inf)` came
    out `3.237e38` (finite) instead of `inf`, a real contract violation.
  All six reverted, bit-identical to prior HEAD (verified via diff: only
  doc-comment additions remain, zero non-comment lines changed). *Two
  compounding lessons, both already documented separately in this file but
  now confirmed together on one idea: (1) always verify a `tune.rs`
  coarse-grid result against the real 100M-sample fuzz -- the grid's
  domain coverage and sampling density both differ from the real one in
  ways that can hide either a regression (exp2_checked/exp2m1's max-ulp
  surprises) or a bug (exp2/exp10's NaN, which the grid's own floor-only
  reachable-x set structurally couldn't reach). (2) `mca` and a wide fuzz
  sweep both passing is *still* not sufficient -- `edgecheck.rs`'s
  dedicated special-value pins (here, `x=inf`) exercise boundary
  conditions neither a random fuzz nor a static instruction-scheduling
  model samples at all; run all three before trusting any reduction-scheme
  change, not just two of them.*

- **sinpi/cospi/sind/cosd/tanpi/tand special-case matrix (2026-07-10,
  resolved -- mostly clean, one cosmetic non-issue found and left as-is)**:
  built the usual matrix (`{0,-0,±1,±0.5,±2,±90,±180,±inf,NaN}`) against
  hand-computed closed forms. Clean except two related, explicable, and
  ultimately non-actionable findings:
  1. **`cosd`'s NaN sign/payload is inverted relative to every sibling in
     this family** (`cosd(inf)`/`cosd(NaN)` both come out positive-NaN
     while `sinpi`/`cospi`/`tanpi`/`sind`/`tand` all come out negative-NaN
     for the same inputs). Root cause: `cos`/`cosd` compute their sign via
     `parity = !kb.to_bits() << 31` (bitwise `NOT` of the reduction
     variable's raw bits, extracting `(k+1)`'s parity from `k`'s own bit
     pattern without materializing `k+1` as a float) where `sin`/`sind`
     use the un-negated `qb.to_bits() << 31`. This trick's derivation only
     holds when `kb` represents a genuine finite `k+0.5`; for non-finite
     `kb` (`x=inf`/`NaN`), the `!` just flips whatever bit happens to sit
     in that position, mechanically producing a different NaN sign than
     the un-negated sibling functions. Per idea #84's already-established
     precedent (erf/erfc NaN sign audit): NaN sign/payload is
     implementation-defined under IEEE754/C99, not a spec violation, and
     this crate has already decided not to spend a real `.abs()`-style op
     on every NaN-producing path just to canonicalize it. Same call here
     -- left as is, not a new class of bug.
  2. **`sinpi`/`sind`/`tanpi`/`tand` don't preserve odd-function sign-of-
     -zero at nonzero integers**: `sinpi(1.0) == sinpi(-1.0) == -0.0`
     (same sign, not flipped) rather than the `-0.0`/`+0.0` pair pure odd
     symmetry would require. Root cause: `r = x - q` where
     `q=x.round_ties_even()` equals `x` itself exactly at any integer,
     and IEEE754 defines `x-x` as *always* `+0.0` regardless of `x`'s own
     sign -- the same "opposite/same-signed-zero operation erases sign"
     mechanism this crate has hit and fixed multiple times before
     (`sinf_poly`'s own `-0` fix, `atan2(-0,+0)`, `sinpi`'s existing
     `x==0.0` guard), just recurring at every nonzero integer instead of
     only at `x=0`. Confirmed by hand-tracing both `sinpi(1.0)` and
     `sinpi(-1.0)` through the actual reduction. **Not fixed**, unlike
     those prior instances, for two reasons specific to this case: (a)
     `sinpi` is a crate-specific function with no C99/std convention
     dictating which signed zero is "correct" at a nonzero integer
     crossing (unlike `sin(-0)=-0`, which *is* standard) -- there's no
     wrong answer to converge to, only an arbitrary choice; (b) confirmed
     via `accuracy.rs`'s own `ulp_diff`/`ord()` that `+0.0` and `-0.0`
     already sort identically (`ord(+0.0)==ord(-0.0)==0`), so this is
     invisible to every accuracy metric this crate tracks, and a real fix
     (detecting exact-integer `x` and selecting a sign from `parity(q)`)
     would cost a real branchless select for a benefit no measurement
     here can see. Left as is. No code changes from this audit.

## cbrt family

- **Seed constant + degree-2 poly joint search (2026-07-07)**: best across
  41 seeds still ~50x over budget (max ulp 112). 3 coefficients can't
  correct this seed's error regardless of seed choice.

- **Integer-division-free seed (2026-07-07)**: cheaper codegen (3 vs 5
  instructions), but too coarse for the correction poly to compensate —
  max ulp 33, ~16x over budget.

- **Range-invariant audit: cbrt_accurate's Df32 Newton step + powf_checked's
  Df32 machinery (2026-07-10, resolved -- both clean)**: continued the
  same technique that found the sin_checked/cos_checked bug, applied to
  the crate's other Df32-based functions. `cbrt_accurate(x)^3 ~= x`
  checked across its entire documented safe range (`2^-56` to `2^127`,
  35.7M in-domain samples) plus explicit edge probes right at `2^-56`,
  `2^100` (the rescale threshold idea #56 tuned), and `2^127`/`f32::MAX`:
  zero violations, every relative error within `~2e-7` of true (matching
  its own near-perfectly-rounded budget). `powf_checked`/
  `powf_checked_unchecked`'s own `x>0 => result>0` invariant (using the
  same `log2_df`/`exp2_checked_df` Df32 pair idea #58 fixed a real bug in
  previously) checked across 50M/24.8M samples: zero violations. Neither
  function shares the specific mechanism that broke `sin_checked`/
  `cos_checked` (a *reduction* whose integer quotient silently loses
  precision at extreme magnitude) -- `cbrt_accurate`'s Newton step and
  `powf_checked`'s log2/exp2 combine don't have an analogous "coarse
  integer count" step to lose precision in. No code changes.

## log_2 / ln / log10

- **log_2/ln/log10/log1p/log2p1 special-case matrix (2026-07-10,
  resolved -- clean)**: checked `{0,-0,±1,±2,±0.5,f32::MIN_POSITIVE and
  its negation,±1e-40 (genuine denormal),±inf,NaN}` against expected
  IEEE754/C99 behavior. All clean: `log(±0)=-inf`, `log(negative)=NaN`
  (including tiny negative denormals), `log(inf)=inf`, `log(-inf)=NaN`,
  `log(NaN)=NaN` (both NaN signs preserved), denormal inputs give
  sensible large-negative finite results (not garbage) for all three of
  `log_2`/`ln`/`log10`. `log1p`/`log2p1` likewise clean at their own
  special points (`-1`→`-inf`, `<-1`→`NaN` via the `u<0` path, `-inf`→
  `NaN`, `inf`→`inf`, and the already-documented `-0.0`→`-0.0` sign
  preservation still holds). No bugs found, no code changes -- unlike the
  newer sinpi/cospi/sind/cosd/tanpi/tand family (a real bug found here as
  recently as idea #31, 2026-07-09), the log family has evidently already
  had this class of edge case shaken out over this crate's longer
  history. *Matches idea #85's own "diminishing returns" signal for this
  technique -- two families in a row (this one and sinpi/cospi's own
  matrix) have now come back clean; worth trying a genuinely different
  angle on the next idea rather than matrix-auditing a third family
  on the assumption this technique still has an easy hit somewhere.*

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

- **Range-invariant audit: hypot family + remainder_wide (2026-07-10,
  resolved -- both already-documented, no new bug)**: applied the same
  "check a fundamental invariant at extreme inputs" technique that found
  the sin_checked/cos_checked range bug to two more double-float/extended-
  range functions. `hypot`/`hypot_unchecked` return exactly `0` for paired
  denormal inputs (violating `hypot(x,y) >= max(|x|,|y|)`) -- but `hypot`'s
  own doc comment already explicitly says "naive sqrt(x^2+y^2), no
  anti-overflow rescaling... trades the overflow/underflow edge cases for
  vectorizability," an already-accepted tradeoff (confirmed `hypot_checked`
  itself has zero such violations, matching idea #83's own prior
  "denormal-pair path, no bug" finding for that specific tier).
  `remainder_wide` showed dramatic-looking violations of `|remainder|<=|y|/2`
  (e.g. `remainder_wide(1e37,1e10)=-1.33e14`) at first -- but every single
  test case had `|x/y|` far beyond its own documented `~2^48` limit
  (its own doc comment: "degrading past that where a single correction
  pass is no longer enough... the same kind of 'harder, separate limit'...
  just much further out", i.e. already anticipated, not new). Re-tested
  strictly within the documented `|x/y| <= 2^48` domain (50M random
  samples, 33.2M in-domain): zero violations. Both functions behave exactly
  as already documented; no code changes. *Same technique, same session,
  a useful negative result this time -- confirms the sin_checked/
  cos_checked bug wasn't symptomatic of a wider pattern in this crate's
  other double-float/extended-range functions.*

## asin / acos / atan / atan2

- **acos: mulsign-reassociate the trailing `+select(PI,0)` (2026-07-10,
  new idea this session, tried and rejected)**: inspired by the
  `atan_latency` fix's `mulsign(a,x)-mulsign(b,x) == mulsign(a-b,x)`
  reassociation (a pure identity on *already-computed* values, no new
  rounding) -- tried the analogous `acos(x) = FRAC_PI_2 -
  mulsign(FRAC_PI_2 - y, x)` in place of the shipped `mulsign(y,x) +
  if x<0.0 {PI} else {0.0}`, hoping to trade a compare+select for pure
  sign-bit arithmetic the same cheap way. Verified specials match (0,
  -0, ±1, NaN, ±inf all bit-identical) but a 100M-sample fuzz found it's
  *not* bit-identical to the original (60,142/49.6M samples differ) and,
  worse, a real accuracy regression: avg ulp barely moved (0.1362->0.1371)
  but max ulp jumped 6->**121**. Root cause: unlike `atan_latency`'s
  trick, which only reassociates two values that already existed with no
  new subtraction between them, this introduces a genuinely *new*
  subtraction (`FRAC_PI_2 - y`) that suffers catastrophic cancellation
  exactly where `y` approaches `FRAC_PI_2` (i.e. `x` near `0`, where
  `acos(x)` itself is near `pi/2`) -- the same "algebraically-exact
  identity reintroduces cancellation" bug class this file's `sigmoid`/
  `atanh` entries already document, now confirmed a third time. Reverted,
  no lib.rs changes (mca not even checked -- the accuracy regression alone
  is decisive). *A mulsign-reassociation trick is only free when it
  reassociates values that already exist untouched; introducing any *new*
  subtraction as part of the reassociation reopens the door to
  cancellation and needs the same full accuracy verification as any other
  algorithm change, not just a bit-identity spot-check.*

  **Follow-up, same day: a genuinely safe variant, still not adopted
  (real but tiny perf regression, zero accuracy benefit)**. Root-caused
  *why* the first attempt failed (it subtracted the continuously-varying
  `y` from a constant) and why `atan2`'s own precedent is safe (it only
  ever subtracts between the two *discrete* values `+FRAC_PI_2` and
  `-FRAC_PI_2`, both exact) -- constructed an acos analog that keeps the
  same discreteness: `correction = FRAC_PI_2 - mulsign(FRAC_PI_2, xn)`
  (exactly `0.0` for `x>=0`, exactly `PI` for `x<0` -- confirmed
  `2*FRAC_PI_2` bit-matches the independently-rounded `PI` constant
  exactly, so this is exact either way), then `mulsign(y, xn) +
  correction`, never touching `y` in the subtraction at all. This *is*
  bit-identical to the shipped form (verified: 99.2M-sample fuzz, zero
  mismatches, every special value matches) and the real 100M-sample
  accuracy.rs sweep confirms unchanged avg/max ulp (0.0676/5, matching
  baseline). But mca showed no speedup -- a real, if small, *regression*:
  throughput 0.820->0.834 cyc/elem (+1.7%, confirmed against a fresh
  git-stash baseline matching readme.md exactly), latency unchanged
  (37.11 both ways). Since this doesn't speed up the function and the
  accuracy is bit-identical (not improved), it clears neither of this
  loop's two bars. Reverted, bit-identical to prior HEAD (confirmed via
  `git diff`). *Even a reassociation that's provably safe on accuracy
  (bit-identical, unlike the first attempt) still needs an actual mca
  measurement before adopting -- "replaces a compare+select with sign-bit
  arithmetic" was the exact reasoning that worked for `atan_latency`, but
  doesn't automatically transfer to every structurally-similar-looking
  select; LLVM's own instruction selection for `acos`'s specific
  surrounding code apparently already handles the original select as
  cheaply or more cheaply than the reassociated form.*

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

- **Dead-code + coefficient-convergence housekeeping sweep (2026-07-10,
  resolved -- crate is clean)**: three checks, all clean. (1) Every
  private `fn` and module-level `const` in `src/lib.rs` has >=1 real call
  site beyond its own definition -- no repeat of idea #98's dead
  `Df32*Df32` multiply find. (2) `cargo clippy` reports zero dead-code/
  unused warnings. (3) Ran every `examples/tune.rs` filter
  (`acos8`/`asinacos`/`asinpoly`/`atan`/`cbrtthroughput`/`erf`/`expm1`/
  `exp_r`/`lnlog10`/`log1pjoint`/`log2atanh`/`sinf`/`exp2lut`, ~20 tuning
  targets total) checking "start" (shipped coefficients) vs "tuned"
  (coordinate descent's own local optimum): every polynomial currently
  shipped is already at (or within noise of) its own coordinate-descent
  optimum on tune.rs's grid -- no free accuracy sitting unclaimed in any
  existing fit. The handful of runs that *did* show a big apparent swing
  (`acos_poly_bh`'s basin-hop: max 3->2 but avg 0.95->1.89; `cbrt_throughput`'s
  own tuner: max 68->43 but avg 6.8->22.7, converging on an absurd
  `-3.88e12`-scale coefficient) are the exact "max-first tuple-ordering
  lets avg blow up" / grid-overfitting traps this file's own
  `tune_basin_hop` doc comment already warns about, not real headroom --
  and `asin_mid`'s own "improvement" is tuning a branch removed from the
  shipped `asin` entirely (fix 6, 2026-07-07), irrelevant to current code.
  `log1p_joint`'s own result exactly reproduces idea #22's already-known
  numbers (0.29664->0.29179), already established not to survive real
  verification. No code changes.

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

- **sinf_poly copysign audit (idea #81) → sin_checked/cos_checked
  `[-1,1]` output invariant fix (2026-07-10, adopted -- real correctness
  bug found, not what the audit set out to look for)**: idea #81 asked
  whether `sinf_poly`'s internal `.copysign(x)` is still load-bearing for
  every caller now that `sin_checked` flips `r`'s sign before the poly and
  `sinpi` has its own `x==0.0` override. Checked empirically (a temporary
  copysign-free build, tested at each of 8 callers' own actual
  zero-crossing, both signs -- not just `x=+-0.0`, which is the wrong
  probe point for the cos-family functions, whose own singular point is
  at their *own* zero-crossing, e.g. `pi/2` for `cos`, not `0`): `sin`,
  `sind`, `cospi`, `cosd` all still genuinely need it (`sin`/`sind` lose
  the odd `x=-0.0` sign; `cospi`/`cosd` lose even-function sign-of-zero
  *symmetry* at their own crossings -- `cospi(0.5)` and `cospi(-0.5)` came
  out with *opposite* signs of zero without it, impossible for a genuinely
  even function). `sin_checked` and `sinpi`'s own separate guards make it
  provably redundant for them specifically -- split into a copysign-free
  `sinf_poly_raw` core plus the existing `sinf_poly` wrapper, confirmed
  bit-identical via a 100M-sample fuzz across every documented-accurate
  bucket for both functions.
  While verifying this at scale (running the full `examples/accuracy.rs`
  sweep rather than trusting the narrow spot-check), found something the
  audit wasn't looking for: `sin_checked`/`cos_checked`'s `[1e15,1e16)`
  bucket showed avg ulp in the *hundreds of millions*, not "gradually
  degraded." Root-caused (not just re-measured): `round_x_over_pi`'s
  double-float `(qh,ql)` representation of `q=round(x/pi)` only resolves
  `q` to ~48 bits total (two f32 mantissas) -- genuinely exact for the
  *transforms* used to build it (`two_prod`/`two_sum` lose no precision
  converting a sum/product to a hi/lo pair), but that doesn't give the
  pair *unlimited* resolution. Once the true `q` itself needs more than
  ~48 bits to pin down (`|x| > 2^48*pi ~ 8.85e14`), `qh+ql` comes out off
  by a few whole integers (confirmed directly: off by 1 at `x=1e15`, by 14
  at `x~1e16`) -- exactly the "relocatable cliff, not a slope" this whole
  double-word scheme was built to eliminate, just relocated from a single
  f32's `~2^24` ceiling to `~2^48` instead of actually removed (an
  incorrect claim, "exact at any magnitude," in this function's own prior
  doc comment -- corrected). Each unit of `q` error shifts the reduced
  residual `r` by a whole multiple of pi, and `POLY_SAFE_BOUND` only
  clamps `sinf_poly`'s *input* (to `+-1000`) -- it does nothing to the
  *output*, and a degree-9 poly evaluated at `|r|=1000` (dominated by its
  own leading `r^9` term) returns values up to `~2.6e21`. Confirmed via a
  direct probe: `sin_checked`/`cos_checked` already silently return values
  like `1.07e9` or `2.6e21` for ordinary (if extreme) finite input on
  *unmodified* master -- a real, pre-existing, previously-undocumented
  violation of the fundamental `|sin(x)| <= 1` invariant, not a regression
  from this session's own changes (verified against a fresh git-stash
  baseline before concluding anything). This is a substantially worse
  failure mode than reduced accuracy -- any caller relying on the basic
  sine/cosine range guarantee (e.g. `sqrt(1 - sin_checked(x)^2)`) could
  silently misbehave. `sind`/`cosd` share the identical mechanism (same
  `POLY_SAFE_BOUND` pattern) but their own doc comment already explicitly
  disclaims correctness past their own `~4.7e7` exactness limit, promising
  only finiteness (which `2.6e21` technically satisfies) -- not the same
  class of broken promise as `sin_checked`/`cos_checked`'s "full-range
  gradual degradation," so left unchanged.
  Fixed with a final `.clamp(-1.0, 1.0)` on both `sin_checked` and
  `cos_checked`'s return value -- doesn't repair the underlying accuracy
  for `|x|` this extreme (a real fix needs a wider-than-double-float `q`,
  a substantially bigger undertaking, left open), but restores the one
  invariant every caller can rely on regardless of `x`'s magnitude.
  Verified a true no-op everywhere the functions were already accurate
  (100M-sample fuzz: every bucket from `|x|<=pi/4` through `[1e12,1e13)`
  unchanged vs. baseline); the `[1e15,inf)` tail's own *ulp* numbers don't
  visibly improve in `accuracy.rs`'s own report (this crate's `ulp_diff`
  measures distance via a global sign+magnitude ordering, so two arbitrary
  values inside `[-1,1]` can still read as billions of "ulp" apart --
  that metric was never going to show this fix), but the actual returned
  values there are now confirmed bounded to exactly `[-1,1]` up to
  `f32::MAX` (previously unbounded up to `2.6e21`). Real, non-zero mca
  cost, accepted under this crate's own established "real perf cost to
  fix a wrong-for-legitimate-input defect" precedent (same shape as
  sin/cos's inf-for-large-x fix and sinh/cosh's domain-hole fix): combined
  with the `sinf_poly_raw` saving above, `sin_checked` latency
  109.02->117.02 cyc (+7.3%), throughput 5.476->5.542 (+1.2%);
  `cos_checked` (clamp only, no offsetting saving) 113.00->122.00 (+8.0%),
  4.537->5.672 (+25.0% -- confirmed via assembly diff to be a real
  port-contention/scheduling shift, not de-vectorization: total
  instruction count only grew 166->173, `codegen_check` clean). readme.md
  mca table updated to match. *A narrow, well-scoped audit (idea #81) that
  ran its verification at real scale (the actual 100M-sample fuzz, not
  just the zero-crossing spot-check it set out to confirm) surfaced a
  much bigger, previously-invisible bug purely as a side effect of
  actually running the full test -- this is exactly why "verify with the
  real harness, not a narrow spot-check" keeps paying off in this crate,
  even for changes that look self-contained going in.*

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

- **asin_small: one more Taylor term (2026-07-08, re-tested 2026-07-10 as
  part of a round-off budget audit -- this entry's own original claim was
  stale, corrected below)**: original 2026-07-08 note said this "relocated
  to a tied co-bottleneck in the acos_poly branch, unchanged overall" --
  re-verified fresh on today's code (fix 8, 2026-07-09, decoupled
  `asin_poly` from `acos_poly` the day *after* this entry was written,
  which this entry's own claim never accounted for) per this crate's own
  "verify a stale prior finding before trusting it" discipline. First,
  root-caused via a hand round-off-budget trace at the crate's own
  documented worst point (x=0.24595731): computed the exact same formula
  in f64 to isolate rounding from truncation. Result: truncation error
  dominates completely (7.03 of the 8.50 total ulp; the f32-vs-f64-same-
  formula rounding gap is only 1.47 ulp) -- the degree-7 Taylor series
  itself, not any rounding step, is the real budget-eater right at
  `asin_small`'s own domain edge (0.25). This means idea #7's "attack the
  largest rounding term" framing doesn't actually apply here (there's no
  single dominant rounding step to fix); this is a truncation problem, and
  the obvious fix (this entry's own proposal, one more Taylor term) is the
  right lever after all. Implemented and measured fresh: exhaustive sweep
  now shows a **real, non-tied improvement**, not the stuck result the
  stale note described -- max ulp 9->7, avg 0.0251->0.0202 (~20% better),
  new worst point exactly `x=0.25010535` (confirmed independently via a
  dedicated exhaustive scan of the `big` branch's own `[0.25,1.0)` domain,
  giving max 7 there in isolation, matching -- this is `asin`'s own
  pre-existing, independent worst point in the *other* branch, now
  unmasked rather than "relocated"). So fix 8's `asin_poly` decoupling
  really did break the tie this entry originally found -- the situation
  changed, and last time's negative result no longer holds.
  **Still not adopted**, but now for a cleaner, better-quantified reason:
  a real, unavoidable throughput cost. Tried two evaluation orders: plain
  Horner extension (`fma` chain now 4 deep instead of 3) measured
  latency 59.03->63.03 cyc (+6.8%), throughput 0.968->1.044 cyc/elem
  (+7.9%); an Estrin-style regroup (matching `exp2`'s own "3 balanced
  pairs" pattern, trading one extra multiply for a shorter critical path)
  measured latency 59.03->60.02 (+1.7%, much better) but throughput
  0.968->1.087 (+12.3%, *worse* than plain Horner) -- confirming this
  crate's own repeated finding that fewer/shorter critical path wins on
  latency can lose on throughput once the extra op competes for port
  bandwidth across many independent vectorized lanes, and that throughput
  (this crate's stated priority metric) is the one that matters more.
  Since this is a branchless design (`asin_small` is unconditionally
  evaluated on every call regardless of which branch's result is
  selected), there's no way to pay for the extra precision only where it's
  needed -- the cost is unavoidable if the term is added at all. A real
  accuracy gain with a real throughput cost on both tries; doesn't clear
  this loop's bar. Reverted both variants, bit-identical to prior HEAD
  (confirmed via `git diff`). *When re-testing a prior "no effect" finding
  after an intervening change to a function it depends on (here,
  `asin_poly`'s decoupling), don't just trust the old conclusion --
  re-verify on the current code, since the entire premise (a tied
  co-bottleneck) can silently stop being true.*

## Backlog round 3 (2026-07-09) — kitchen-sink brainstorm, all UNTESTED

Constraints as always: must autovectorize (branchless selects only), budget
≤0.5 avg / ≤2 max ulp, accuracy-for-speed trades fine inside that. FP
divider still nearly idle. Cross-checked against every rejection above —
each entry is either new or explicitly distinguished from its rejected
cousin.

### Approximation theory / fitting (cross-cutting)

1. ~~**Round-based reduction for exp2/exp2_checked/exp10**~~ (tried
   2026-07-09/10, rejected on all 6 functions sharing this poly --
   see the "exp2 / exp / sin / cos / tan" section near the top of this
   file for the full writeup). Degree 5→4 falsified outright (max ulp
   13, not competitive); degree 5 itself won on `tune.rs`'s coarse grid
   but failed real verification differently on every function: NaN
   inside the documented domain (exp2/exp10 unchecked), a real >2x
   throughput regression from floor→round port contention
   (exp2_checked), a real max-ulp regression with a real-if-smaller perf
   cost (exp2m1), a real accuracy win *with* a real perf cost
   (exp2_checked_df), and a genuine overflow-saturation correctness bug
   at `x=inf` only `edgecheck.rs` caught, not the fuzz or mca
   (exp10_checked, which had otherwise looked like the one clean win).
   "Magic-round already proven faster than vroundps in exp" doesn't
   transfer here: `exp`'s own magic-round is fused into a multiply the
   reduction needs anyway, while exp2/exp2_checked's would-be swap adds
   a bare, unfused add/sub pair with nothing to share it with.
2. ~~**Gather-based real LUTs (vgatherdps)**~~ (screened 2026-07-10,
   rejected before any fitting -- see the "Codegen / build / measurement"
   section near the top of this file for the full mca writeup). A
   standalone `vgatherdps`-based 16-entry LUT probe, given a generously
   *shorter* residual poly (degree 2, 2 coefficients) than the shipped
   direct fit (degree 5, 6 coefficients) to make the comparison as
   favorable as possible to the gather approach, still cost ~41% more
   throughput overall (0.647->0.913 cyc/elem) than the direct polynomial.
   Root cause: `vgatherdps` itself measures at 5 uOps / 22 cyc latency /
   4.00 cyc RThroughput per 8-wide gather on this CPU -- as expensive as
   roughly 8 FMAs. Two gathers per 16-wide throughput iteration (8.0 cyc)
   alone very nearly matches the *entire* direct-poly baseline's whole
   10.35 cyc/iteration budget. "A hardware gather is a different cost
   model" turned out to be right, just not in the hoped-for direction --
   confirms the backlog's own "screen with mca before any fitting" caution
   was exactly the right call here.
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
   `sinh`'s own `0.5` threshold checked directly (2026-07-10, no code
   changes -- confirms no headroom): measured `sinh_small` (exact Taylor)
   and the direct `exp_pos_neg`-based branch *independently* across
   `[0,3)` in half-wide buckets. The small branch wins clearly through
   `[0.25,0.5)` (avg 0.36 vs the direct branch's 1.44) and the direct
   branch only just starts winning at `[0.5,0.75)` (avg 0.89 vs 1.05) --
   the crossover is already almost exactly where the shipped threshold
   sits, same conclusion this file already reached for `asin` (idea #37,
   "0.25 already close enough") and `erf` (idea just above, "plateau...
   no headroom"). `sinh`/`asin`/`erf` are now all confirmed with no
   threshold-placement headroom; `expm1`/`atanh` not independently
   re-checked this way but `expm1`'s own round-off budget audit (idea #7)
   already found its dominant error split roughly evenly between
   rounding and truncation with no single cheap lever, a related
   (if not identical) conclusion.
7. **Round-off budget audit per function**: enumerate every rounding on the
   critical path with a bound, attack the largest term. This is exactly how
   exp's Cody-Waite fix was found; do it systematically for the remaining
   >2-max-ulp functions (asin 9, expm1 6, tanh 6, sinh/cosh 5, erf 5).
8. **Binary-function worst-case mining**: unary functions get exhaustive
   sweeps; powf/atan2/hypot/remainder only get fuzz. Guided search
   (branch-and-bound over exponent-pair classes, or fixed y/x ratio
   classes for atan2) would find real worst cases fuzz misses.

   ~~**atan2** (idea #39, resolved 2026-07-10)~~: structured octant/ratio
   sweep found nothing worse than the documented max ulp 3 -- see idea
   #39's own entry.

   **powf (tested 2026-07-10, real find: true max ulp is at least ~2.5x
   worse than documented, not fixed, doc corrected)**: unlike atan2, this
   one found something real. Random 10M-sample fuzz documents `powf`'s
   max ulp as 127; a structured sweep -- fixing `x` near specific values
   (both near 1.0 at ulp granularity, and log-spaced across the full
   range) and solving for `y` so the domain-constraining product
   `y*log2(x)` lands at chosen targets, rather than hoping uniform random
   sampling stumbles onto a bad combination -- found a first worse point
   almost immediately (207 ulp at x=1.0013627, y=63095.734, product
   ~123.96), verified two independent ways (`(x as f64).powf(y as f64)`
   and a manual `(y*ln(x)).exp()`, both agreeing to 10+ significant
   digits against each other and disagreeing with `powf`'s own output by
   exactly 207 ulp -- not a probe artifact this time, unlike idea #39's
   own first-pass bug).

   Local refinement around that point, then broader targeted sweeps
   specifically biasing `y` so the product approaches the domain's own
   upper boundary (128) from below, found progressively worse points
   across several rounds: 207 -> 229 -> 260 -> 277 -> 312 ulp (last found
   at x=1.1969347, y=485.83997, product ~125.9999). Each refinement round
   found something a little worse without fully converging -- **312 is a
   confirmed real lower bound on the true worst case, not a proven exact
   supremum**; a more exhaustive search would likely find something
   somewhat worse still. `powf_unchecked` is confirmed bit-identical to
   `powf` at every one of these points (as its own doc comment promises),
   so it inherits the same true max. `powf_checked` is measurably better
   at the same three spot-checked points (192/82/107 ulp vs. `powf`'s
   312/209/207) but still far from clean -- consistent with its own
   documented "substantially more accurate for large `|y|`" framing
   (better, not perfect).

   All of this concentrates in the region where `y*log2(x)` approaches
   the domain's own upper edge (128) rather than being uniformly spread
   across the whole valid range.

   **Root-caused properly in a follow-up pass, correcting an initial
   misattribution**: first guessed this was "a harder quantitative
   confirmation of" this file's older "residual lives in the
   log2_df/exp2_checked_df double-float chain" finding (see "Other
   spots") -- wrong. Checked `powf`'s actual source: it computes
   `exp2_checked(log_2(ax) * y)`, a single f32 `log_2` call and a single
   f32 multiply -- it never touches `log2_df`/`exp2_checked_df` at all;
   that Df32 machinery is `powf_checked`-only. Tested the real mechanism
   directly instead of assuming: computed `log_2(ax) * y` two ways at
   each of the three found worst points -- once as the plain single
   multiply `powf` actually does, once as a Dekker/fma-style compensated
   two-product (`p = l*y`, `e = l.mul_add(y, -p)`, giving `p+e` as the
   multiply's own correctly-rounded result). The two-product form's
   *residual* error (after fully compensating the multiply itself)
   matched `(log_2(ax)'s own rounding error) * y` almost to the digit at
   all three points (ratio 1.000 in all three cases, computed
   independently). **The multiply contributes essentially nothing --
   100% of the error traces to `log_2(ax)` itself only ever being
   accurate to a single f32's ~24 bits, and that fixed absolute error
   getting scaled up by whatever `y` happens to be before `exp2_checked`
   exponentially amplifies it.** This means there's no cheap compensated-
   multiply fix available (confirmed by testing one: manually applying
   the two-product correction through `exp2_checked`'s own `fma(result,
   e*LN_2, result)` trick recovered only ~13% of the error, e.g. -312 ->
   -270 ulp at the first point -- consistent with the multiply being a
   minor contributor). The *only* real fix is a higher-precision `log2`
   in the first place -- exactly what `log2_df` (routed through
   `exp2_checked_df`) already provides, at the real extra cost
   `powf_checked` already pays for exactly this reason. So this doesn't
   open a new, cheap lever after all: it's a complete, decisive
   confirmation that `powf`'s fast/accurate split is already drawn in the
   only sensible place, just with the *actual* size of the accuracy gap
   now measured for the first time (>=312 ulp, not the previously-assumed
   127) rather than a vague "residual lives somewhere in there" shrug.
   Updated readme.md's `powf`/`powf_unchecked` max ulp from `127`/`135`
   to `>=312` (avg ulp columns left alone -- these bad points are sparse
   enough that 10M-sample fuzzing's *average* almost certainly isn't
   materially affected, only its claimed *max* was wrong). *A "worse than
   documented" find via structured search doesn't need a fix to be worth
   committing -- correcting a wrong documented bound to the true, verified
   one is itself the deliverable, same as this session's several mca/
   accuracy staleness fixes elsewhere.*

   **remainder (checked 2026-07-10, confirms existing documented behavior,
   nothing new found)**: extended the same structured-search technique to
   the last binary function idea #39 left open. Targeted `q = round(x/y)`
   landing on large magnitudes (2^10 through 2^30) and on exact
   half-integer ties across many `y` scales. Found a real invariant
   violation almost immediately -- `remainder`'s own fundamental guarantee
   `|result| <= |y|/2` failing by more than 3x at `x=5.330646e37,
   y=7.9432826e29` (`remainder` returned `1.2840237e30` against a true
   value of `-3.046328e29`, `|y|/2 = 3.972e29`). But checking `remainder`'s
   *own doc comment* before treating this as new: it already explicitly
   states "`round(x/y)*y`'s absolute error scales with ulp(x), which
   swamps the true remainder... once `|x/y|` is large" and points to
   `remainder_checked`/`remainder_wide` as the designed fix. Tested both:
   `remainder_checked` *also* violates the invariant at this same point
   (`4.897e29`, still over `|y|/2`, though less badly than plain
   `remainder`) -- but `remainder_checked`'s own doc comment *also*
   already says so explicitly ("Confirmed by fuzzing... 0 max ulp for
   `|x/y|` up to `1e7`, degrading only past `2^24`... a separate, harder
   limit this correction can't reach past" -- my point's ratio is
   `~6.7e7`, past both thresholds). `remainder_wide` handles it correctly
   (`-3.0463283e29`, matching the true value to the last couple of ulp).
   So: not a new bug, just a structured re-derivation of an already-fully-
   documented, already-solved-by-tiering limitation (same
   fast/near-tie-safe/wide-safe three-tier shape as `powf`/`powf_checked`
   above). Given `remainder_checked`'s doc comment makes a specific,
   falsifiable claim ("0 max ulp for `|x/y|` up to `1e7`"), tested *that*
   directly instead: a dense structured sweep (many `y` scales, `x`
   chosen to land on or extremely near half-integer ties every ~5 units
   of `q` up to `1e7`, 2.814 billion `(x,y)` pairs total) found **zero**
   counterexamples -- the documented claim holds up under targeted search,
   not just the random fuzz that originally established it. No code
   changed; this closes idea #8's `remainder` extension with high
   confidence rather than leaving it as an unverified claim. *Before
   treating a structured-search find as new, check the target function's
   own doc comment for an existing, already-quantified acknowledgment --
   `powf`'s worse-than-documented number was genuinely new because no
   such acknowledgment existed; `remainder`'s wasn't, because it did.*
9. **Structured-error probes**: plot per-function error vs mantissa and vs
   exponent separately; periodic structure invites a cheap structural
   correction (one select or exponent-derived fma) instead of a refit.
   Tried once (2026-07-10) on `erfc` (the crate's own worst-behaved
   function, max ulp ~109) via a 50M-sample error-by-exponent-bucket
   sweep (`x`'s own `floor(log2(|x|))`, libm's `erfc` as a structural
   reference): found smooth, monotonic growth from `avg ulp ~0.6` at
   small `|x|` up to `avg ulp ~8.2`/`max ulp 109` right at the `x=10`
   clamp boundary, plus a much smaller secondary bump around `|x|` in
   `[0.03,0.25]` (avg ulp ~3.3-4.5, max ulp only up to 10, not the
   dominant contributor). No genuinely new *periodic* structure found --
   the dominant large-`x` growth is exactly the mechanism idea #53
   already root-caused (`erfc(x)` itself shrinks toward 0 as `x` grows,
   so the same absolute error reads as ever-larger ulp), not a missed
   correction term. Still open for a function this technique hasn't been
   tried on yet, or for chasing the smaller `[0.03,0.25]` bump
   specifically if someone wants the average-only, modest payoff.
10. ~~**Monotonicity/oddness harness metrics**~~ (tried 2026-07-10,
    resolved -- two real deviations found, both explained by already-
    accepted tradeoffs/noise, no new actionable fix): built a standalone
    sweep checking ~15 monotone functions for monotonicity violations and
    ~17 odd functions for exact `f(-x)==-f(x)` across large samples.
    Mostly clean (`log_2`/`ln`/`log10`/`exp2`/`exp2_checked`/`exp_checked`/
    `atan`/`sigmoid`/`erf`/`asin`/`cbrt`/`log1p`/`asinh`/`acosh` all
    monotonic; `sin`/`sin_checked`/`atan`/`asin`/`sinh`/`sinh_checked`/
    `erf`/`asinh`/`atanh`/`cbrt`/`rcbrt` all exactly odd). Four flagged
    deviations, all already explained:
    - `sinpi`/`sind`/`tanpi`/`tand` "oddness violations" at `x=0` are
      exactly the already-known, already-decided-not-to-fix sign-of-zero
      non-issue (confirmed: the reported "diff" is bit-for-bit `0.0`,
      i.e. a `+0.0`/`-0.0` bit-pattern mismatch with no value difference
      -- see this file's earlier sinpi/cospi/sind/cosd/tanpi/tand
      special-case-matrix entry).
    - `tan`'s oddness mismatch (9997/20M samples, max abs diff ~9.3e-10)
      and `tanh`'s (5.57M/20M samples -- ~28%!, max abs diff ~2.4e-7,
      ~2 ulp) both trace to the same structural cause: neither is built
      from a directly-negated-argument construction the way `sin`'s own
      `x - x^3*p(x^2)` form is (`tan=sin/cos` as an independent ratio;
      `tanh` computes `expm1(2x)`, itself *not* an odd function, and only
      becomes odd through the nonlinear `e/(e+2)` combine) -- so there's
      no natural bit-exact symmetry to preserve, and small (already
      budgeted, ~1-2 ulp) rounding noise in each independently-computed
      branch shows up as an oddness mismatch at this granularity. Fixing
      this for `tanh` specifically would need the same abs-then-mulsign
      restructuring (`mulsign(tanh_of_abs(x), x)`) this crate's own doc
      comment says was *already tried* (for a different reason, an
      overflow/domain-hole fix) and rejected on a mixed latency/throughput/
      avg-ulp tradeoff -- re-litigating that same rejected tradeoff for
      oddness alone isn't a new argument.
    - `tanh`'s reported monotonicity violation (~x=8.66, y drops by
      exactly 1 ulp) is a genuine oscillation between `1.0` and one ulp
      below it as `x` approaches saturation -- traced directly (a dense
      per-ulp scan from `x=8` to `9`): the *true* `tanh(x)` in this range
      sits so close to `1.0` that whether the correctly-rounded f32 result
      lands on `1.0` or one ulp below is genuinely sensitive to sub-ulp
      variation in the true value, and the function's own existing,
      already-accepted ~1-2 ulp noise floor is enough to flip which side
      of that boundary the *computed* result lands on for adjacent `x`
      values -- a natural consequence of existing budgeted noise meeting a
      saturation boundary, not a new, larger defect (`tanh`'s own max ulp
      budget is already 6-8, far exceeding the single-ulp oscillation
      here). No code changes; the harness itself (not committed, a
      standalone probe) is cheap to reconstruct if a future session wants
      to re-run it after some other change.
11. **Differential testing vs sleef/core-math/rlibm** built locally, not
    just f64-rounded references — also catches double-rounding artifacts in
    accuracy.rs's own reference path.
12. ~~**Standing test: every `_unchecked` bit-matches its checked sibling on
    the documented domain**~~ (implemented 2026-07-10, `examples/
    unchecked_parity.rs`): fuzzed all 10 checked/unchecked pairs
    (`log_2`/`ln`/`log10`, `cbrt`/`cbrt_accurate`, `atan2`, `fmod`,
    `remainder`, `powf`/`powf_checked`) at 50M samples each against the
    exact domain closures `accuracy.rs` already uses for its own ulp
    reporting -- all 10 came back bit-identical, confirming the doc
    comments' promises are genuinely true today. Exits nonzero on any
    mismatch so it can be run as a real regression gate, not just a
    one-off check. `edgecheck.rs`'s existing 2-3 pinned spot-values per
    pair weren't a systematic sweep; this is. No bugs found this run, but
    the standing test itself (not just the audit) is the deliverable --
    it locks the invariant in for future coefficient/logic changes to
    either half of any pair.
13. **Generalize the cbrt_accurate recipe** (cheap ≤1-ulp core + one Df32
    Newton step) into a template: candidates rsqrt_accurate,
    exp_accurate/ln_accurate (each is the other's Newton residual),
    sin_accurate near zeros. Paper-screen the residual budget first.
14. **ln_accurate/log2_accurate tier from existing log2_df (tried 2026-07-09,
    rejected)**: implemented `log2_accurate(x) = log2_df(x).to_f32()` plus
    `log_2`'s own domain selects, right after idea #58's own log2_df
    precision fix (so this got the *best-case* version of log2_df, not
    the pre-fix one). **No measurable accuracy benefit**: a 20M-sample
    fuzz gave identical avg/max ulp to plain `log_2` (0.0061/3 both), and
    a 1.63-billion-sample strided sweep across the entire positive-normal
    domain found only 7 bit-differing outputs total (~4e-9 of the
    domain). Root cause: `log2_df`'s extra double-float precision only
    matters once something *downstream* (like `powf_checked`'s multiply
    by `y`) amplifies the low-order bits it preserves -- collapsing
    straight back to a single f32 with no such amplification, `log_2`'s
    own single-rounding `fma(p,s,k)` already lands on the same
    correctly-rounded result almost every time, since f32's own 24-bit
    output resolution can't distinguish the extra precision in the first
    place. Unlike `cbrt_accurate` (which gains real accuracy from Newton
    iteration's quadratic convergence on the seed error, a genuinely
    different mechanism), collapsing a double-float log doesn't transfer
    that same win. Reverted, no lib.rs changes survived.
15. ~~**Integer fixed-point poly evaluation**~~ for mantissa-only reductions
    (log's s): i32 mul-high chains free up FMA ports. Precedent warning:
    parity()'s integer version lost to FP ports — but that was 3 ops, not
    a whole poly; port-pressure math differs at scale.

    **Screened via llvm-mca's own resource-pressure model on this crate's
    actual target (2026-07-10, rejected before writing a single line of
    the real fixed-point poly)**: rather than build a whole fixed-point
    log_2 mantissa-poly evaluator (a real, multi-hour numerical-fitting
    effort -- picking a Q-format, re-deriving/rescaling all 10
    coefficients, handling the f32-mantissa-to-fixed-point and
    fixed-point-back-to-f32 conversions at both ends) only to maybe learn
    the premise doesn't pay off, checked the idea's *entire* claim --
    "i32 mul-high chains free up FMA ports" -- directly against
    `llvm-mca -mcpu=native --resource-pressure` first, the same
    "screen with mca before any fitting" discipline idea #2's gather-LUT
    probe already validated as the right order of operations. Built a
    minimal hand-written `.s` file (no Rust, no crate code at all) with
    isolated `# LLVM-MCA-BEGIN/END` regions for the actual instructions
    involved, run through the exact `llvm-mca` binary/CPU target
    (`tigerlake`, confirmed via `llvm-mca -mcpu=native --version`) this
    crate's own `examples/mca.rs` already uses.

    The premise is false on this CPU, decisively: `vfmadd213ps` (ymm) --
    this crate's actual poly-evaluation instruction -- is 1 uOp,
    RThroughput 0.50, scheduled on `ICXPort0`/`ICXPort1` (alternating,
    i.e. 2 FMAs/cycle combined across both ports). `vpmulld` (ymm, plain
    32x32->32 truncated multiply, the most obvious "integer multiply"
    choice) is *worse*, not better: 2 uOps, RThroughput 1.00, latency 10
    (vs FMA's 4) -- and both of its uOps land on `ICXPort0`+`ICXPort1`
    *simultaneously* per instruction, i.e. double the port pressure of a
    single FMA for one multiply, before any addition/combine step is even
    considered. The idea's own more-precise reading (real fixed-point
    poly evaluation needs a *widening* mul-high, not a truncated `vpmulld`)
    fares better in isolation -- `vpmuldq` (32x32->64 signed widening) is
    1 uOp/RThroughput 0.50, matching FMA's own port cost exactly -- but a
    real mul-high *term* still needs a shift to extract the high half
    (`vpsrlq`, 1 uOp/0.50 RThroughput) and an add to fold in the next
    coefficient (`vpaddd`, 1 uOp/0.33 RThroughput) to replicate what a
    single `fma(c,s,acc)` does atomically in one op. Measured the full
    3-instruction chain back-to-back against two `vfmadd213ps`s directly
    in the same probe: the mul-high chain totals 3 uOps/~1.33 combined
    RThroughput per term vs FMA's 1 uOp/0.50 -- **~2.6x more port
    pressure per equivalent poly term, not less**, and every one of those
    extra uOps (`vpmuldq`, `vpsrlq`, `vpaddd`) is scheduled on the *exact
    same* `ICXPort0`/`ICXPort1` (plus `ICXPort5` for the add) that FMA
    already uses -- Tiger Lake's vector int and vector fp domains share
    execution ports here, they aren't separate resources the way the
    idea's "frees up FMA ports" framing assumes. This also doesn't yet
    count the real, additional cost a full implementation would still
    need on top: extracting `s` into fixed-point form from the f32
    mantissa bits, and converting the final fixed-point accumulator back
    into a normal f32 result -- both non-free steps this probe didn't
    even have to model to already lose. Consistent with, and now
    concretely quantified beyond, this file's older `parity()`
    integer-bit-ops precedent ("FP-port ops beat integer ops once
    scheduled") and idea #98's Karatsuba audit (assumed infrastructure
    cost that turned out not to apply where hoped) -- a third instance of
    "the port/instruction-count story sounds right until you check the
    actual scheduling model for *this* CPU." No lib.rs change was ever
    made (nothing to revert); the standalone `.s` probe isn't part of the
    crate. *A cross-cutting perf idea phrased at the instruction-class
    level ("integer ops free up FP ports") is itself falsifiable with
    `llvm-mca --resource-pressure` alone, with zero Rust code -- always
    check whether the target CPU's actual port model agrees before
    spending real implementation effort on a coefficient refit or
    numerical-format redesign the premise doesn't survive.*

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
19. **exp10 third Cody-Waite word (surveyed 2026-07-09, not pursued)**:
    exhaustive sweep shows exp10/exp10_checked already at avg ulp
    0.0343/max ulp 2 (2.2B+ samples each) -- already tight (max 2 is
    close to the practical floor for a non-perfectly-rounded function).
    Minimal headroom for a third reduction word to collect; not worth the
    extra fma given how little room there is to improve.
20. **exp2int construction via pure-integer path (tried 2026-07-09 via
    #21, rejected)**: this idea's own "saturating-cast problem doesn't
    apply -- k is bounded" claim was checked directly and found false --
    see #21.
21. **exp2_checked: k1 from bit-twiddled k instead of a second
    magic-round (tried 2026-07-09, rejected)**: implemented `k as i32`
    + `>>1` integer halving in place of the fma+sub magic-round, exactly
    as described. Accuracy was (as expected) unaffected -- integer
    halving of an already-exact integer is trivially exact either way.
    But `codegen_check` immediately caught a real, serious regression:
    `k as i32` (Rust's saturating float-to-int cast) does *not* vectorize
    to a simple packed convert instruction, even though `k` is runtime-
    bounded by the function's own earlier `.clamp(-151.0,128.0)` -- LLVM
    can't statically prove that bound from a `.clamp()` call at
    codegen time, so it still emits the full saturating-cast fallback
    (`cvttss2si` scalar extract+convert+check), the *exact* known failure
    class this crate's own `codegen_check` comment already documents.
    De-vectorized `exp2_checked_throughput` and, via their own calls into
    `exp2_checked`, `erf_throughput`/`erfc_throughput`/
    `powf_throughput`/`powf_unchecked_throughput` too -- 5 regions
    silently broke the crate's own hard "must auto-vectorize" requirement.
    Reverted immediately, bit-identical to prior HEAD. Also directly
    refutes idea #20's own premise ("saturating-cast problem doesn't
    apply — k is bounded") -- it does apply, runtime boundedness isn't
    visible to the codegen path that decides cast lowering. *A value
    being runtime-bounded by an earlier `.clamp()` doesn't make an `as`
    numeric cast free of the saturating-cast vectorization trap -- LLVM's
    codegen decision for `as i32` doesn't see through arbitrary earlier
    control flow, only genuinely bounded integer types or explicit
    unchecked-cast intrinsics would; always run codegen_check immediately
    after introducing any `as <int>` cast in a public function's hot
    path, not just when the value's own type looks unbounded.*

### log family

22. **log1p output-side correction refit (tried 2026-07-09, rejected --
    coarse-grid win didn't survive real-sweep verification)**: built the
    requested joint objective in `examples/tune.rs` (new
    `tune_joint_fixed0` + `log1p_via_ln_c`, matching `log1p`'s exact
    shipped construction, `ln(u)+c/u` routed through `ln_poly_c` so a
    shared coefficient set affects both) -- minimize log1p's error
    subject to ln's own max never regressing past its shipped baseline,
    same "constrained search" shape as acos_poly's own joint asin refit
    (fix 7). The coarse grid (~613-step walk over log1p's own bit
    patterns) reported a modest, real-looking win: log1p's own max ulp
    unchanged (4), avg improved ~1.6% (0.29664->0.29179), and `ln`'s own
    max even improved as a bonus (3->2). But wiring the resulting
    coefficients into `src/lib.rs`'s actual `ln_normal` and checking
    against the *real* 100M-sample `accuracy.rs` fuzz (not tune.rs's own
    coarse grid) showed a genuine regression instead: log1p avg ulp
    0.0966->0.1085 (worse) and max ulp 4->7 (worse) -- the coarse grid's
    ~613-step spacing apparently doesn't sample densely enough near
    whatever region the real random-fuzz's own worst case actually lives
    in. Reverted the coefficient change in `src/lib.rs`, bit-identical to
    prior HEAD; kept the tuning infrastructure itself in `tune.rs` (the
    *technique* -- constrained joint refit sharing one coefficient set
    across two call sites -- is reusable for a future attempt with a
    denser or differently-distributed grid, even though this specific
    run's result didn't survive verification). Commit `9b420a6`
    (tune.rs only, no src/lib.rs change kept). *This is exactly the
    "always verify a tune.rs coarse-grid result against the real
    100M-sample fuzz before trusting it" lesson this file's own
    `tune_basin_hop` doc comment already documents for acos_poly, now
    confirmed a second time for a completely different function pair --
    a coarse grid reporting "no regression, modest improvement" is not
    sufficient evidence on its own, regardless of how principled the
    constrained-search setup looks.*
23. **log10 exact-decade check (resolved 2026-07-09, no bug)**: checked --
    `log10(10^n)` is exact (bit-for-bit `== n as f32`) for every
    `n in [-38,38]`.
24. **log_2/ln/log10 shared-core macro**: the three `_normal` bodies are
    identical modulo constants (and log2_df quadruplicates it). A macro or
    generic-const core removes 3 hand-synced copies — hygiene, zero perf
    claim, reduces refit-application errors.
25. ~~**log_2 denormal path: fold the ×2^24 rescale into the wrapping_sub
    magic**~~ (killed on paper 2026-07-10, confirmed dead on arrival as
    suspected -- not just "probably", provably so): the idea's own hope
    was that `log_2_normal`'s exponent-extraction bit trick
    (`(x.to_bits() as i32).wrapping_sub(0x3f3504f3) >> 23`) could absorb
    the needed ×2^24 denormal rescale as a *constant* adjustment folded
    into that same subtraction, avoiding the real `x * 16777216.0`
    multiply `log_2`'s caller currently does before calling in. This
    can't work for a structural reason, not just a fiddly implementation
    one: the bit trick relies on every input being an IEEE754 *normal*
    number, where the stored mantissa bits mean `1.mantissa × 2^(e-127)`
    (an implicit leading 1) -- a fixed, uniform relationship between bit
    pattern and value that a constant integer offset can shift correctly.
    Denormals have no implicit leading 1: the stored mantissa bits mean
    `0.mantissa × 2^-126`, so recovering a normalized `(e, m)` pair
    requires knowing *where the mantissa's own leading 1 bit sits* --
    and that position is different for every denormal. Concretely:
    `f32::MIN_POSITIVE/2 = 2^-127` (mantissa `0.1000...0`, leading bit at
    position 22) needs a 1-bit renormalizing shift, while the smallest
    denormal `2^-149` (mantissa `0.0...01`, leading bit at position 0)
    needs a 23-bit shift -- two denormals, two completely different
    required corrections, not a shared constant. A single `wrapping_sub`
    offset can only ever apply one fixed correction to every input alike,
    so it's fundamentally the wrong tool regardless of which constant is
    chosen. The real `x * 16777216.0` multiply isn't a convenience
    method for adding 24 to an exponent -- it's doing genuine,
    data-dependent renormalization work (IEEE754 multiply hardware
    correctly shifts each denormal's mantissa by exactly however much
    *that* value needs, then produces a valid implicit-leading-1
    representation), which is exactly the class of computation a fixed
    bit-level offset cannot replicate. The only bit-level alternative
    (counting each lane's mantissa leading-zeros and shifting by a
    *variable*, data-dependent amount) is its own real op, likely no
    cheaper than the single multiply it would replace, and this crate
    has no existing lzcnt-style idiom to build on. No code changed --
    this was a pure paper analysis, no benchmarking needed since the
    premise fails on a representation argument, not a performance one.
26. **koff-free unchecked-log fast path audit (checked 2026-07-09, no-op)**:
    `--emit=asm` on `log2_unchecked_throughput`'s region confirms LLVM
    already inlines `koff=0.0` through and feeds the converted exponent
    directly into the final `fma(p,s,k)` -- no separate `vaddps` for the
    dead add anywhere in the region. Duplicate of the already-logged
    "Integer koff fold (2026-07-08)" finding (same conclusion, different
    entry point into the same question). No code changed.

### sin / cos / tan (radians, half-turns, degrees)

27. **sinpi/cospi: two_prod(π, r) + derivative correction (tried 2026-07-09,
    rejected)**: implemented exactly as described -- `(p,e) =
    two_prod(PI, r)`, `sinf_poly(p) + e*(1-p²/2)` via a single fma. The
    premise's own "π·r's single rounding is likely the dominant error"
    turned out to be only half right: `sinpi`'s own average ulp was
    *bit-for-bit unchanged* (0.1969 both ways -- the correction does
    shift ~0.23% of individual outputs by 1 ulp in a random-sample probe,
    but not enough to move the aggregate at all), meaning `sinf_poly`'s
    own ~4-term minimax fit error already dominates `sinpi`'s budget, not
    the reduction's rounding step -- "two cheap ops" was never going to
    help there regardless of cost. `cospi` did see a real, if modest,
    improvement (avg ulp 0.0861->0.0768, ~11%; near-a-zero max ulp also
    dropped from 823550->411775, though that's an already-documented
    "huge relative error at a true zero, harmless in absolute terms"
    artifact either way) -- `cospi`'s own extra reduction step
    (`k=round(x-0.5)`, `r=(x-k)-0.5`, one more subtraction than `sinpi`'s
    plain `r=x-q`) apparently does leave more rounding on the table for
    this correction to recover. But "two cheap ops" wasn't actually
    cheap: mca showed a real, substantial throughput regression on both
    -- `sinpi` 1.149->1.400 cyc/elem (+21.8%) for *zero* accuracy gain,
    `cospi` 1.283->1.653 cyc/elem (+28.8%) for the 11% avg-ulp gain.
    Neither clears this loop's own bar (speedup, or accuracy gain *without*
    a perf penalty) -- `sinpi_tp` fails on both axes, `cospi_tp`'s real
    accuracy gain comes at a real, non-trivial cost. Reverted, bit-identical
    to prior HEAD. *A "single inexact step in an otherwise-exact reduction"
    is only the dominant error term if everything downstream of it is
    exact or near-exact -- here `sinf_poly`'s own several-ulp-of-headroom
    polynomial fit swallowed the correction whole for `sinpi`, and only
    `cospi`'s extra reduction rounding gave the fix any real room to work;
    measure both functions sharing an idea separately, don't assume a
    shared reduction trick pays off identically for both.*
28. **sind/cosd: same trick for d·DEG_TO_RAD_SMALL (two_prod variant tried
    2026-07-09, rejected -- literal HI/LO split still untested)**: tested
    idea #27's `two_prod`-based correction (recovers `d*DEG_TO_RAD_SMALL`'s
    own rounding error, not literally this idea's proposed HI/LO constant
    split, but attacking the same reduction step) on `sind`/`cosd`
    separately per this session's own "test siblings independently, don't
    assume shared payoff" lesson (from idea #27's own sinpi/cospi
    finding). Result was more clear-cut than sinpi/cospi's own mixed
    outcome: **zero** measurable accuracy improvement on *both* functions
    (`sind` avg ulp 0.1237 unchanged -- actually 0.1247, marginally worse
    within noise; `cosd` avg ulp 0.0725 bit-for-bit identical), while mca
    showed a real cost on both: `sind` throughput 1.151->1.461 cyc/elem
    (+27%), `cosd` 1.406->1.798 (+27.9%). Reverted, bit-identical to
    prior HEAD. Consistent with idea #27's own conclusion that
    `sinf_poly`'s polynomial fit (not the reduction's own rounding)
    dominates the error budget here -- if recovering the *multiplication's*
    rounding error doesn't help at all, the *constant's* own quantization
    (idea #28's actual literal proposal, a HI/LO split of DEG_TO_RAD_SMALL)
    likely wouldn't either, though that specific variant remains formally
    untested if someone wants full closure. *When a backlog idea's
    literal proposal (constant split) is close in spirit but not
    identical to an already-tested one (residual correction via
    two_prod), testing the already-implemented variant first is a cheap,
    informative proxy -- a zero-benefit result there is a strong (if not
    airtight) signal the untested literal variant would fare similarly,
    since both attack the same reduction step and the dominant error
    source (the poly fit) is unaffected by either.*
29. **tanpi / tand (implemented 2026-07-09)**: added as plain
    `sinpi(x)/cospi(x)` and `sind(x)/cosd(x)` ratios -- correct by
    construction, not an approximation: `tan(pi*x)` has period 1 in
    half-turns (unlike `sin`/`cos` individually, which flip sign every
    integer), so whichever integers `sinpi`/`cospi`'s own reductions
    resolve to, their respective sign corrections cancel exactly in the
    division. Poles (`cospi(x)==0` at half-integers) are handled for free
    by IEEE754 division giving the correctly-signed `+-inf`; both inherit
    their sin/cos siblings' existing `NaN`-at-infinity convention with no
    new special-casing. Verified: fuzz accuracy clean (`tanpi` avg ulp
    0.275, `tand` avg ulp 0.177 -- `tanpi`'s max ulp is a huge but
    harmless near-pole artifact, same class as `cospi`'s own documented
    near-zero blowup, confirmed by direct inspection: the denominator is
    a genuinely tiny nonzero value there, not a bug), `codegen_check`
    clean (70 regions, up from 68), 21 new edgecheck pins, mca (`tanpi`
    66.97/2.521 cyc lat/throughput, `tand` 70.02/2.533 -- reasonable,
    expected cost for composing two existing calls plus a division, no
    surprises). Commits `883104a` (code+harness), `bcf3eee` (readme
    sync). **Not pursued**: the idea's own secondary suggestion (skip the
    parity xors entirely, fusing the two reductions into one that never
    computes an unused sign correction) -- this implementation still
    pays for `sinpi`'s and `cospi`'s own independent parity computations,
    which happen to cancel *algebraically* in the ratio but aren't
    actually eliminated from the generated code. Left open as a real,
    separate follow-up for whoever wants to chase the extra throughput
    (would need a shared reduction that produces both sin/cos values
    from one parity computation, more invasive than this drop-in
    composition).
30. **sinpi/cospi/sind/cosd harness coverage (resolved 2026-07-09)**:
    sinpi/cospi were in the sweep, but domain-restricted to |x|<1e6, well
    below 2^22 -- widened to the full f32 range, which is what caught
    idea #31's real bug below. sind/cosd's own |x|<1e6 restriction stays
    (their own, unrelated, still-accurate ~4.7e7 limit).
31. **cospi large-x probe (confirmed real bug, fixed 2026-07-09)**: the
    suspicion was right -- not "by luck", a genuine bug. sinpi/cospi's
    magic-round trick (x + 1.5*2^23) is only exact for |x|<=2^22, unlike
    every other magic-round use in this crate (which only ever rounds an
    already-small reduced value, not raw unbounded input); for
    2^22 < |x| < 2^24 it silently gave *wrong* results (e.g.
    cospi(2^22+1) ~ -0.0033 instead of -1.0), contradicting the doc
    comment's "exact out to f32::MAX" claim. Fixed with `x.round()`
    (full-range correct) + the existing `parity()` helper instead of the
    magic-constant's bit-trick sign extraction; cospi restructured via
    `cos(pi*x)=sin(pi*(x+0.5))` without ever forming `x+0.5` or `k+1` as
    single floats (same lossy-for-large-x trap). Also fixed an identical
    bug in accuracy.rs's own `cospi_ref` (same trap at f64's ~2^53 limit)
    and a `sinpi(-0.0)` sign regression the fix introduced (`x.round()`'s
    `-0.0` case loses its sign in the following subtraction -- same
    "opposite-signed-zero" mechanism as sinf_poly's own -0.0 fix,
    guarded the same way). Real mca cost (~25-45% both functions),
    accepted -- correctness for documented in-range inputs isn't
    optional. See git history (commit 5f13f7d) for full details.
32. **Intermediate sin/cos tier (|x|≤~1e5)**: single extra correction word
    over the fast tier's 4-fma Cody-Waite, well short of checked's full
    double-float q — a third point on the speed/domain curve if any user
    workload actually sits there. Only on demand.
33. ~~**sin fast tier: drop PI_D for a fitted 3.5-word split**~~ (tried
    2026-07-10, rejected -- also fails near zeros, just less catastrophically
    than the naive 4→3 cut): no Python/sollya available, so screened
    numerically in Rust instead. First pass looked genuinely promising: a
    freshly-derived third word (`PI_C_NEW = round_to_f32(pi - PI_A - PI_B)`
    computed at f64 precision, *not* a reuse of the existing `PI_C`, which
    carries trailing zero bits chosen for the old 4-word chain's own needs)
    measured *better* raw reduction accuracy than the shipped 4-word
    version across a 40M-sample random fuzz of the raw residual alone (max
    abs error 1.04e-7 vs 1.19e-7). But wiring it into the *real* sin/cos
    construction (full poly + parity combine, not just the bare residual)
    and comparing against the actual 100M-sample everywhere-domain fuzz
    told a completely different story: sin avg ulp 0.0645->0.1399 (worse),
    max 412->**1,824,546**; cos avg 0.2917->0.3320 (worse), max
    2780->**386,929**. Root-caused the worst point (`x=-7.5558225e6`,
    `q=-2405093`): this is a genuine near-a-zero-of-sin case (the true
    reduced residual is `~2.37e-7`, tiny), and while the shipped 4-word
    reduction still resolves it to within `~2.4e-8` absolute error, the
    3-word version's absolute error there is enough to be off by *~11%
    relative* to the (already tiny) true residual -- the same "huge
    relative/ulp error at a true zero, small in absolute terms" shape this
    crate has repeatedly documented elsewhere (cospi's own artifact), but
    quantitatively far worse here (both avg *and* max ulp regress, not just
    an already-tolerated near-zero spike). Confirms the earlier "4→3
    rejected" finding's premise was right even with a properly re-fit
    (not just truncated) third word: an entire word's worth of pi's own
    precision genuinely can't be recovered by any single replacement word
    across sin's whole documented domain. Reverted, no lib.rs changes.
    *A narrow probe (bare reduction-residual error, uniformly sampled) can
    look like a clean win while completely missing a context-dependent
    failure mode (relative-error blowup specifically near the target
    function's own zeros) that only shows up once wired into the real
    end-to-end construction and measured against the real, full-domain
    fuzz -- exactly the discipline this crate's own accuracy.rs was built
    to enforce, and exactly what a quick standalone probe skips by
    default.*

### asin / acos / atan / atan2

35. **atan: correct the 1/a fold's division rounding (tried 2026-07-09,
    rejected)** — implemented `e = fma(recip, a, -1.0)`, `corr = -e·y/(1+y²)`
    folded into the existing `FRAC_PI_2 - y` else-branch (needed an
    `is_finite` guard on `corr` first: `e` is a `0*inf` NaN at `a=inf`,
    which broke `atan2`'s extreme-ratio inputs badly before the guard was
    added — caught by accuracy.rs's own u64-wraparound garbage output,
    not reasoned about in advance). Once correctness was fixed: **zero
    accuracy benefit** — exhaustive sweep showed atan's avg ulp exactly
    unchanged (0.0675) and max ulp *regressed* 3→4 (right at the a=1 fold
    boundary, x≈1.022), confirming the worst case doesn't actually live
    where the backlog guessed. Also a severe mca cost: atan throughput
    1.491→**4.196** cyc/elem (+181%), atan2 1.532→2.411 (+57%) — the extra
    division landed far worse than "two-ish ops" suggested. Reverted,
    bit-identical to prior HEAD. *A root-caused "worst case lives here"
    conclusion from an unrelated refit (the numerator-only LP attempt)
    doesn't mean a *targeted* fix at that location will actually move the
    number — measure, don't just aim.*
36. **acos_accurate opt-in tier (tried 2026-07-09, superseded by a bigger
    win)**: implemented the proposed Df32 pi/2 + two-product
    sqrt(1-a)*poly(a) combine as a standalone tier -- measured *zero*
    improvement (avg/max ulp identical to plain acos, same "extra
    precision doesn't survive collapsing back to f32 without a
    downstream user" lesson idea #14's own rejection established).
    While investigating why, found the real story: `acos_poly`'s own
    leading constant (meant to be exactly pi/2) was parsed from a
    literal one ulp short of the correctly-rounded value, a real,
    previously-uninvestigated bug in the *shared default* itself (not
    the "protected coefficient" idea #36's own text worried about --
    this wasn't a refit choice, just a transcription slip). Fixed
    directly: avg ulp 0.4962->0.0676 (-86%), zero perf cost. See
    `acos_poly`'s own doc comment for the full story. Commit `0862c58`.
    No separate accurate tier needed -- the shared default's own
    accuracy already improved dramatically for free.
37. **asin small branch: check 1-a rounding in [0.25,0.5) (quantified
    2026-07-09, real but modest, not actioned)**: probed the "big"
    branch's own accuracy specifically in `a in [0.25, 0.5)` (both a
    dense linear sweep and a random-bit-pattern sweep, cross-checked
    against each other) -- avg ulp ~1.4, max ulp 6-7 in that band alone,
    noticeably higher than asin's own domain-wide average (0.0251, per
    fix 8) but *not* a new max-ulp record (the domain-wide worst case,
    max 9, sits at x~0.246, just inside `asin_small`'s own domain, per
    fix 8's own note -- confirmed still true). So this specific concern
    doesn't set the binding constraint on `asin`'s worst case, but it is
    a real, quantifiable, disproportionate contributor to the *average*:
    `asin`'s own accuracy.rs sweep uses the unrestricted `everywhere`
    domain, and random-bit-pattern sampling gives each octave inside
    `(0,1]` roughly equal weight regardless of its linear width -- with
    ~149 such octaves down through the denormal range each contributing
    comparably, a single octave averaging 1.4 (vs. a "well-behaved"
    octave's typical <0.5) is a disproportionately large single
    contributor to the overall 0.0251 average, even though no individual
    octave dominates outright. Fixing this would need another numerical
    refit (via `examples/tune.rs`) of either `asin_poly` or the mid-range
    construction specifically targeting this octave -- real work, and
    the expected payoff is average-only (this session's own refit
    attempts elsewhere have repeatedly found average-only gains modest
    relative to the effort, e.g. idea #52's own "effort/payoff ratio"
    call). Not pursued this session; left open with the concrete numbers
    above for whoever wants to chase the average further.
38. ~~**asin via atan2(x, sqrt((1-x)(1+x)))**~~ (tried 2026-07-10,
    rejected -- decisive perf regression, mixed accuracy result too):
    implemented literally as `atan2(x, ((1.0-x)*(1.0+x)).sqrt())` and
    fuzz-compared (50M in-domain samples) against shipped `asin`. The
    accuracy ceiling probe's own premise was half right: max ulp really
    is much better (9->4), but avg ulp is nearly *double* (0.0506->0.0953)
    -- a real max/avg tradeoff, not a clean win, before even considering
    speed. And speed was never close: mca showed a catastrophic
    regression, not just "probably slower" -- latency 59.03->100.22 cyc
    (+69.8%), throughput 0.968->**3.622** cyc/elem (+274%, nearly 4x).
    Decisive enough that no further exhaustive verification was needed.
    Reverted (probe only, no lib.rs changes). *"Probably slower" turned
    out to undersell it by an order of magnitude -- routing a single-
    branch function through a full binary function (`atan2`, itself a
    division plus a poly plus several selects) costs far more than
    intuition from "it's just one extra call" suggests; always mca a
    cross-function reformulation like this before trusting the "probably"
    in a backlog note.*
39. ~~**atan2 octant-symmetric exhaustive harness**~~ (resolved 2026-07-10,
    no worse case found than the existing fuzz-mode max): built exactly
    as described -- a standalone scratch probe (`examples/
    atan2_worstcase_probe.rs`, not committed) sweeping 1294 fixed y/x
    ratio classes (log-spaced from `1e-30` to `1e30` in both signs, plus
    the 8 octant-boundary angles `tan(k*pi/8)` and close neighbors) times
    601 log-spaced magnitudes in both signs, both "fix the ratio, vary x"
    and "fix the ratio, vary y" directions -- 3,105,968 total `(y,x)`
    pairs, dwarfing the *combination coverage* atan2's own 10M-sample
    uniform-random fuzz gets structurally (random pairs rarely land
    exactly on an octant boundary or an extreme ratio; this sweep targets
    exactly those). Result: max ulp 3, identical to readme.md's
    documented fuzz-mode number -- the structured sweep didn't uncover
    anything the random fuzz was missing.

    Caught a real bug in the *probe itself* before trusting this result:
    the first pass computed the f64 reference from the un-rounded sweep
    value (e.g. `y.atan2(x)` from the original f64 `x`), while
    `atan2(yf, xf)` received the *rounded* `f32` inputs -- for ratios
    near the sweep's extreme end, casting the f64 magnitude to f32
    silently overflows to `+-inf` even though the original f64 value was
    finite, so the two sides were being compared against genuinely
    different inputs. This produced a spurious "3 billion ulp" max before
    the fix (comparing `atan2(finite, inf)` against a reference computed
    from `atan2(finite, finite-but-huge)`). Fixed by computing the
    reference from `yf as f64`/`xf as f64` -- the actual post-rounding
    f32 inputs, not the pre-cast sweep values -- after which the result
    dropped to the real, unremarkable max ulp 3. *When ulp-testing against
    a higher-precision reference, the reference must be computed from the
    exact value the function under test actually receives, not from
    whatever higher-precision value was used to construct it -- a cast
    that overflows/rounds differently than expected will silently compare
    against the wrong input otherwise.* hypot's own max ulp is already 1
    (bounded) domain-wide, leaving little room for this technique to find
    anything there. Applied to `powf` next -- see idea #8's own entry
    below for a real, substantial worse-case find. `remainder` also
    checked (2026-07-10, see idea #8's own entry) -- confirms existing
    documented behavior rather than finding anything new.
40. **atan_latency: fold FRAC_PI_2-p select into sign trickery (adopted
    2026-07-09)**: implemented as described -- apply `mulsign` to `p` and
    `FRAC_PI_2` individually first (`mulsign(a,x) - mulsign(b,x) ==
    mulsign(a-b,x)` exactly: `mulsign` is a sign-bit XOR, and IEEE754
    subtraction commutes exactly with negating both operands, so this is
    a pure reassociation, not a new approximation), then select between
    `sp` and `hpisignx - sp` instead of selecting on the unsigned value
    and applying one final `mulsign`. Verified bit-identical against the
    prior formulation over ~20M random bit patterns plus every special
    value (0, -0, ±1, ±inf, NaN) before trusting it. Paper reasoning
    beforehand (counting ops in each formulation) actually suggested this
    *wouldn't* help -- the new form needs one extra `mulsign` (2 instead
    of 1) -- so this was measured anyway per this session's own
    "measure, don't just reason" precedent, and the paper reasoning
    turned out to be an incomplete predictor of the actual generated
    code: mca showed a small, real, reproducible throughput win
    (1.611->1.591 cyc/elem, -1.2%, confirmed via repeat runs) with
    latency unchanged (59.09->59.11, within noise). A modest win, not a
    dramatic one -- adopted because it's free (zero accuracy cost, no
    latency cost), not because the throughput gain alone justified the
    effort. Commits `997dec1` (code), `7bcafe7` (readme sync). *Even when
    a hand-counted op tally says a reassociation shouldn't help (or
    should hurt), LLVM's actual instruction selection/scheduling can
    still land differently than the naive op-count model predicts --
    worth a quick mca check before dismissing an idea on paper alone,
    especially when verifying bit-identical accuracy is cheap.*

### hyperbolics / sigmoid

42. **tanh via exp_pos_neg ratio (tried 2026-07-09, rejected)**: implemented
    exactly as described -- `(ep-en)/(ep+en)` via `exp_pos_neg`, plus a
    small-`|x|` Taylor branch to avoid the cancellation `sinh_small`/
    expm1's Pade branch also exist for. Two real problems, one fixed, one
    fatal to the idea's own premise. First, the literal formula as
    described is unsafe at moderate-to-large `|x|`: reusing a
    wide-saturating exp helper (or even `exp_checked`'s own `88.7228`
    overflow bound directly) lets `ep`/`en` reach literal `inf`, and
    `inf/inf = NaN` -- not the correct `+-1.0`. Fixed with a much tighter
    `+-40.0` clamp (`tanh` already saturates to exactly `1.0f32` by
    `x~11`, so this has huge margin on both sides without ever reaching
    the field-split's own overflow edge). Second, the suggested 4-term
    Taylor branch (`x - x^3/3 + 2x^5/15 - 17x^7/315`) is nowhere near
    accurate enough out to the idea's own implied `|x|<0.5` threshold
    (matching `sinh_small`'s own convention) -- `tanh`'s Taylor series has
    a much smaller analytic radius of convergence than `sinh`/`cosh`'s
    (a real pole at `+-i*pi/2` vs. entire/no poles), so at `x=0.5` the
    first dropped term alone is already ~800 ulp of error. An initial
    fuzz run at the `0.5` threshold confirmed this directly: avg ulp
    1.05, max ulp 1302. Swept the threshold empirically (hand-deriving
    the exact truncation bound seemed more effort than just measuring)
    and found `0.2` is the sweet spot -- avg ulp 0.0239, max ulp 6,
    matching plain `tanh`'s own max-ulp floor, with the 4 terms as given
    (no need for more). So the accuracy side is a genuine, substantial
    win (avg ulp 0.024 vs `tanh`'s existing 0.146, ~6x tighter on
    average, same max). But the idea's own core premise -- "division is
    idle, may beat expm1's max 6" -- didn't survive mca measurement: this
    formula needs `exp_pos_neg` to evaluate **two** full polynomials
    (`p_pos` and `p_neg`, sharing only the reduction) where `tanh`'s
    existing expm1(2x) route evaluates **one**, and both routes already
    have exactly one division at the end (`e/(e+2.0)` vs `(ep-en)/
    (ep+en)`) -- so the "idle divider" framing was never the actual cost
    driver, the extra polynomial evaluation is. Measured: latency
    86.73->89.08 cyc (+2.7%), throughput 2.039->2.565 cyc/elem (+25.8%) --
    a real, not marginal, regression on the axis this idea was supposed
    to win on. Since the loop's own bar requires *either* a speedup *or*
    an accuracy gain *without* a perf penalty, and this is a real accuracy
    gain *with* a real perf penalty (the same three-way-tradeoff shape as
    idea #46's atanh/log1p rejection), reverted -- bit-identical to prior
    HEAD. *A plausible-sounding resource-idleness argument ("division is
    idle") is only as good as identifying the actual bottleneck it's
    supposed to route around -- here the real cost was a doubled
    polynomial-evaluation count that had nothing to do with the divider at
    all, and only showed up once actually measured with mca rather than
    reasoned about on paper.*
46. **atanh via single log1p on |x| + mulsign (tried 2026-07-09, rejected)**:
    implemented exactly as described — this time it does *not* die
    catastrophically (max ulp only 3→4, not 3→31303 like the earlier
    signed-x rejection), confirming the odd-symmetry reasoning was right:
    `2a/(1-a)` for a=|x| only ever approaches log1p's benign `u→+inf` edge,
    never the derivative-diverging `u→-1` edge that killed the earlier
    attempt. But it's still a real, exhaustive-confirmed net loss on
    accuracy (avg 0.0313→0.0352, max 3→4) with a genuinely mixed perf
    result: mca throughput improved a lot (4.649→3.289 cyc/elem, -29%,
    halving the log1p count really did help) but latency got *worse*
    (74.06→76.22, +2.9%) — three-way tradeoff (better throughput, worse
    latency, worse accuracy), not a clean win on either axis. Reverted,
    bit-identical to prior HEAD. *Sidestepping a catastrophic failure mode
    doesn't mean the replacement is free of smaller, real cost — still
    check the actual before/after numbers, not just whether the disaster
    case recurred.*
48. **sinh_accurate/cosh_accurate tier**: Df32 through the exp combine —
    only if a user asks; max 5 is comfortably documented.

### erf / erfc

49. **erfc exponent via mantissa-mask hi/lo split (tried 2026-07-09,
    rejected)** (Cephes trick): implemented as described -- `xh = xa` with
    low 12 mantissa bits masked (`xh²` exact), `xa² = fma(xa-xh, xa+xh,
    xh²)` (difference-of-squares identity, one rounding instead of
    `xa*xa`'s one rounding but on a much smaller correction term).
    **Zero effect**: exhaustive sweep bit-identical to baseline in every
    field (avg 0.3106, max 109, worst x=9.00551, confirmed via git stash
    on the unmodified code). Root cause: this only compensates the
    *squaring* step; the very next op, `-xa2 * LOG2_E`, is still a single
    uncompensated multiply, and that rounding apparently dominates
    wherever the true worst case actually lives -- the already-rejected
    two_prod fix's real (if modest) 109→93 improvement must come from
    compensating *that* multiply too, not just the square. Reverted
    (5 extra ops for literally zero measured benefit, not even a
    borderline case). *A "compensate step N" idea only helps if step N is
    actually the dominant error source at the specific worst-case point --
    verify with the real worst-x, don't assume from "this step also
    rounds."*
50. **erf tail via expm1-shape (tried 2026-07-09, rejected)**: implemented
    `b = mulsign(-exp2m1(erf_poly(...)), x)` in place of
    `mulsign(1.0 - exp2_checked(erf_poly(...)), x)`, now that idea #65
    supplies the needed `exp2m1` helper. **Zero accuracy effect**: exhaustive
    sweep bit-identical to baseline (avg 0.3166, max 5, worst x=0.28000325,
    confirmed via git stash) -- a direct point-by-point comparison over
    711k samples in [0.28,10] found *literally* zero differing bit patterns,
    not just equal aggregates. Root cause: erf_poly's output at the branch's
    own worst point (x≈0.28) is already ~-0.53 in magnitude, past
    exp2m1's own |x|<0.5 Pade/direct split -- so exp2m1 lands on the
    *same* direct branch (`fma(p,t2,-1.0)`, one fused rounding) that
    `exp2_checked` effectively already gets close to; the theorized
    "1-(≈1) cancellation" doesn't actually occur at the real worst point.
    Also a real, deterministic **mixed perf regression**, not a clean
    win: mca latency 85.97->70.42 cyc (-18%, genuinely better) but
    throughput 2.788->3.167 cyc/elem (+13.6%, worse) -- paying for
    exp2m1's own unused Pade branch (an extra division) on every
    vectorized call for zero accuracy return. Reverted, bit-identical to
    prior HEAD. *A "helper X now exists, idea Y needs it" dependency
    being unblocked doesn't mean the original idea's own error-source
    theory was right -- check where the swapped formula's branch
    selection actually lands relative to its own internal thresholds
    before trusting a "removes a cancellation" story.*
51. **erfcx(x) = e^{x²}·erfc(x) (implemented 2026-07-09)**: the idea's own
    "just the n/d rational, no exp at all" premise held up exactly for
    `x>=0` -- factored `erfc`'s own n/d rational into a shared
    `erfc_rational(xa)` helper (verified by direct code inspection to be
    a pure move, bit-identical to `erfc`'s prior body: same clamp, same
    fma sequence, same exponent computation, just relocated), and for
    `x>=0` `erfcx` collapses to exactly that rational alone -- the
    `e^{x²}`/`e^{-x²}` factors cancel algebraically, so this really does
    sidestep the exponent computation entirely, not just numerically
    approximate the cancellation. For `x<0`, needed one real
    `exp2_checked` call after all (via `erfc`'s own reflection identity,
    `erfcx(x) = 2·exp(x²) - erfcx(-x)`) -- `erfcx` genuinely diverges to
    `+inf` there, so this isn't avoidable, just correctly saturating
    instead of wrapping to garbage. Verified against known reference
    values (`erfcx(1)~0.4276`, `erfcx(5)~0.1107`, `erfcx(10)~0.05614`,
    all matched) and a fresh f64 reference (`exp(x²)·erfc_u15(x)`, safe
    over this `|x|<=10` domain since `x²<=100` is nowhere near f64's own
    ~709 exponent overflow). Fuzz accuracy: avg ulp 0.377, max ulp 125
    (same order as `erfc`'s own 0.311/109-ish budget, one more rounding
    step on the negative side accounts for the difference). 7 new
    edgecheck pins, codegen_check clean (71 regions, up from 70). mca:
    throughput 2.278 cyc/elem (a real, modest win vs `erfc`'s 2.530,
    from skipping the `z`/`w` sign-handling erfc pays) -- but the
    *latency* number (39.36 vs erfc's 64.03) is a `mca mix()` sign
    blind spot artifact, not trustworthy: the harness's latency chain
    always folds its value into `[2,4)` (mask away the sign bit
    entirely), so `erfcx`'s `x<0` branch (with the real `exp2_checked`
    call) is never exercised there, same pitfall this crate has hit
    before for other sign-dependent functions -- documented directly in
    `erfcx`'s own doc comment so this doesn't get miscited later.
    Incidentally also found `erfc`'s own mca latency number in readme.md
    was stale (78.09, now re-measured at 64.03 with this refactor in
    place) -- confirmed via direct code-diff that the refactor itself
    is behavior-preserving, so this is a re-measurement catching drift,
    not a regression from the change. Commits `994dc98` (code+harness),
    `945b9ea` (readme sync). **General lesson: whenever a new function's
    correctness depends on a runtime sign branch, check whether either
    benchmark harness's own input-generation scheme (mca's `mix()`,
    quickbench's own latency chain) can actually reach both signs before
    trusting its latency number -- if the harness folds the chain into a
    fixed-sign range (as `mix()` does here), the branch with the real
    cost can go completely unmeasured in latency mode while still
    showing up correctly in throughput mode's real mixed-sign array.**
52. **erf_poly's a0 ≈ 3.4e-5 (screened 2026-07-09, naive substitution
    fails hard; full refit not attempted)**: the cheap first check --
    naively zero `a0` without refitting anything else -- confirms it's
    genuinely load-bearing, not just numerically small by fitting
    coincidence: exhaustive sweep avg ulp 0.3166->2.3367, max ulp 5->539.
    Reverted immediately. The idea's actual proposal (a full LP refit
    with `a0` *constrained* to 0, letting the other 6 coefficients
    compensate) is real numerical-fitting work -- a new scipy/linprog
    setup against an accurate erf reference, not a quick edit -- and the
    speculative payoff if it works is modest (at most one `fma`->`mul`
    substitution in `b0`, which this session's own `log2_df` and
    exp2_checked pure-integer screens both suggest often isn't a real
    speedup even when it "removes an op"). Not pursued further given the
    effort/payoff ratio; left open for a session that wants to invest in
    the full LP setup.
53. **erfc negative-side accuracy survey (resolved 2026-07-09, structural,
    not actionable)**: split the exhaustive sweep by sign (temporary
    accuracy.rs domain split, not kept) -- confirmed a real asymmetry:
    x>=0 avg ulp 0.4815/max 109 (worst x=9.00551, the already-known case),
    x<0 avg ulp 0.1397/max 6, comfortably inside this crate's usual
    budget. Root cause is structural, not a fixable blended-refit
    artifact: for x<0, erfc(x)=2-(the same n/d rational), and the output
    magnitude there is near 2 (large), so a given absolute error in the
    shared rational corresponds to far fewer ulp than the *same* absolute
    error does for large positive x, where erfc(x) itself is tiny
    (approaching underflow) -- ulp is a *relative* measure, and the two
    sides sit at very different output magnitudes for the same input
    magnitude. The max-109 worst case is already the one idea #49's own
    mask-based/two_prod investigations targeted (and partially improved,
    109->93, at real extra cost) -- this survey doesn't open a new,
    cheaper avenue, just confirms *where* the existing hard case lives and
    *why* a differential (sign-split) refit wouldn't help (the underlying
    rational is shared; the asymmetry is in how ulp itself scales with
    output magnitude, not in the fit quality per side).

### cbrt / sqrt / hypot / pow / remainder

54. **cbrt seed division codegen check (resolved 2026-07-09, no change
    needed)**: grepped `cbrt_throughput`'s own asm region for div/mul --
    LLVM already lowers `ax / 3` (u32, compile-time-constant divisor) to
    `vpmuludq` (multiply-high strength reduction), no `idiv`/`vpdivd`
    anywhere in the region. The `vdivps` instructions also present there
    are the real, unrelated fp division in cbrt_normal's own seed/Newton
    reciprocal (already known cheap on this CPU, see the divider-idle
    finding elsewhere in this file). Nothing to fix.
56. **cbrt_accurate via Halley from a cheaper seed (tried 2026-07-09,
    rejected -- but reached exact accuracy parity along the way)**:
    implemented as described -- `cbrt_fast` seed + a Halley step
    (`y_new = y - y*e/(2y^3+x)`, `e=y^3-x`, derived from the standard
    `y_{n+1}=y_n-2ff'/(2f'^2-ff'')` form) via `Df32`, instead of
    `cbrt_normal`'s own polynomial-refined seed + a Newton step. First
    correction to the idea's own premise: `cbrt_fast` measured (not
    assumed) at avg 57.3/max 554 ulp on its own positive-normal domain,
    not the "~5 ulp" the idea guessed -- an order of magnitude rougher.
    Continuous math still comfortably supported Halley closing that gap
    in one step (57 ulp relative error cubed is still ~1e9x tighter than
    f32 needs), so implemented and measured anyway rather than rejecting
    on the wrong number alone. Two real bugs found and fixed while
    getting there, both instances of the same "Df32 upgrade creates a
    new overflow path a plain fma never would" lesson `remainder_wide`
    already established this session: (1) the naive `2*y^3+x`
    denominator, computed at full magnitude, can itself approach/exceed
    f32::MAX for ordinary large x (not just extreme values) -- reformed
    to `3*x+2*e` (algebraically exact), which helped but didn't fully
    fix it; (2) the real culprit was forming `y*e` as an intermediate
    *before* dividing by the denominator -- for large x with a rough
    seed, `y*e` alone can overflow f32 even though the final ratio
    `e/den` is a small, ordinary, bounded correction (e.g. x=1e37:
    y~2.15e12, e~1.67e32, y*e~3.6e44 > f32::MAX, while e/den~5.6e-6 is
    unremarkable) -- fixed by computing `e/den` first, multiplying by
    `y` after. With both fixed, a third, narrower issue remained:
    ~0.18% of a random sweep (all concentrated in the top ~2 bits of the
    exponent range, x approaching f32::MAX) still showed real error (up
    to max ulp 422) -- traced to `cbrt_accurate`'s own large-magnitude
    rescale threshold (`2^127`) being tuned for the original Newton
    formulation, not this one; the `e` residual for this rougher seed
    needs more headroom before it starts losing relative precision.
    Lowering the "big" rescale trigger from `2^127` to `2^100` (generous
    margin past where errors actually started) fixed it completely:
    exhaustive-adjacent fuzz sweep came back avg ulp 0.0000, max ulp 1 --
    bit-for-bit matching `cbrt_accurate`'s own existing precision, plus
    every special value (0, -0, +-inf, NaN, +-f32::MAX,
    +-f32::MIN_POSITIVE) verified matching exactly. So the idea's core
    premise -- a much cheaper seed really can reach `cbrt_accurate`'s own
    precision in one Halley step -- is *true*, fully validated. But the
    implied payoff ("deleting cbrt_normal's poly... " suggesting a
    speedup) did not materialize: mca showed a real, substantial
    *regression*, not a win -- latency 59.06->75.16 cyc (+27.3%),
    throughput 3.129->3.165 cyc/elem (+1.2%, essentially a wash). The
    extra division and sign-reconstruction bit trick this formulation
    needs cost more than `cbrt_normal`'s own polynomial refinement step
    ever did. Reverted, bit-identical to prior HEAD (no code changes
    kept). *A "might let X reach Y" premise can turn out completely
    correct on the accuracy axis (worth confirming, not dismissing on a
    wrong assumed ulp count) while still failing on the axis that
    motivated trying it in the first place -- accuracy and cost are
    genuinely separate questions, and cubic convergence closing an error
    budget says nothing about whether the arithmetic needed to get there
    is actually cheaper than what it replaces.*
58. **powf: root-cause the log2_df/exp2_checked_df ~150 max ulp (fixed
    2026-07-09)**: traced a concrete worst case against a Decimal-
    precision Python reference at every intermediate step -- neither
    candidate mechanism was it. The real cause: `log2_df`'s own
    `Df32::from_add(k, p*s)` combines `k` with an *already-rounded*
    single-multiply `p*s`, so `p*s`'s own rounding error never enters
    either Df32 word -- harmless when `k` dominates, but for `x` near 1
    (`k=0`) the whole double-float result was silently no more accurate
    than a plain f32 multiply. Fixed with a real two-product
    (`Df32::from_mul(p,s)`). See `log2_df`'s own doc comment for the full
    numbers (avg ulp -13%, max ulp -19%, >100ulp count -55%, real mca
    cost accepted under the opt-in-tier precedent). Commit `994d1ff`.
59. **exp2_checked_df: second-order lo correction (moot, per #58's
    finding)**: the binding term was never exp2_checked_df's first-order
    lo truncation -- it was log2_df's own dropped p*s rounding, now
    fixed. No longer applicable as originally framed.
60. **powf special-exponent select tier (screened 2026-07-09, rejected)**:
    implemented just the `y==2.0 -> x*x` case as a single-select probe
    (matching "screen with mca first"). Real accuracy win where it hits:
    powf(x,2.0) avg ulp 5.42/max 45 -> exact 0/0. But mca showed a real,
    non-negligible cost paid on *every* call regardless of `y` -- since
    this crate's branchless-select style computes every branch
    unconditionally, llvm-mca gave identical numbers whether the
    benchmark's own `y` was `2.0` (the "hit" case) or `2.5` (a "miss"),
    confirming the extra select's cost isn't hidden by the fast path ever
    being cheaper to reach: throughput 5.098->5.319 cyc/elem (+4.3%),
    latency 103.05->101.08 (-1.9%, a wash). Extrapolating to the full
    4-case ladder the idea originally proposed (`y` in `{1,2,0.5,-1}`)
    would compound this further for benefit that only manifests when `y`
    lands on one of exactly 4 values -- narrow relative to powf's whole
    input space, and callers who specifically need exact squaring/sqrt/
    reciprocal already have a zero-cost workaround (write `x*x`,
    `x.sqrt()`, `1.0/x` directly, or use `pown_const::<2>` for a
    compile-time-known integer exponent). Reverted (single-case probe,
    no lib.rs changes survived).
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

### New API surface / tiers

68. **lgamma (Stirling + reflection)**: big job, listed for completeness —
    the largest gap vs libm's function set that fits this crate's
    branchless style.
70. **Slice/batch API + runtime multiversioning**: backlog round 1 has the
    slice tier; add `is_x86_feature_detected` dispatch at the slice level
    (per-call dispatch is un-inlinable, per-slice is free).

### Codegen / build / measurement

71. **Blend lowering audit under AVX-512VL (resolved 2026-07-09)**: grepped
    the entire compiled `mca_target` assembly -- **zero** `vblendvps`
    instructions anywhere, vs. 6915 AVX-512 mask-register-predicated
    instructions (`{%k0}`-style). With `target-cpu=native` on this
    machine (which supports AVX-512VL), LLVM exclusively lowers every
    f32 select/blend in this crate to mask-register form, never the
    older AVX2 `vblendvps`. Whether this *relieves port 5 pressure*
    can't be answered by static analysis alone (would need real
    perf-counter profiling, not `llvm-mca`'s own scheduler model) -- but
    the sheer volume of mask usage isn't uniformly "free": this
    session's own `is_finite()`/`is_nan()` sweep (idea #90's follow-up)
    already found a real subset of that mask usage is the *expensive*
    kind (per-lane scalar extraction + int comparisons + manual
    `kshiftlb`/`korb` mask reassembly, not a single vectorized compare
    producing a mask directly), fixed in `powf_checked` and confirmed
    absent elsewhere via the same audit. No further action from this
    entry specifically -- the "which lowering form" question is answered
    (mask registers, universally), and the "is it all cheap" question is
    already covered by the more targeted `vextractps`/`kshiftlb`-signature
    audit from idea #90's own follow-up, not by counting mask
    instructions in the aggregate.
72. **`-C llvm-args=-force-vector-interleave=N` sweep (tried 2026-07-09,
    rejected -- real but genuinely mixed, not adoptable as a global
    default)**: swept N=2/4/8 via `RUSTFLAGS` over the *entire* mca
    suite (~70 functions), comparing against a freshly-captured default
    baseline from the same session (same machine state, no thermal
    confound since mca is a static model). N=4 and N=8 came back
    **bit-for-bit identical** to the default in every single region --
    LLVM's own default choice for this harness's fixed-size (`ARR_LEN=16`,
    two 8-wide AVX2 vectors) loop already effectively *is* 4 (or the
    array is small enough to fully unroll regardless, making the
    interleave knob moot at that point). N=2 was the only setting that
    actually changed anything, and it's a real, substantial, genuinely
    *mixed* result across the whole suite, not a clean win: roughly 60%
    of functions improved (several by double digits --
    `powf_checked` -14.0%, `powf_checked_unchecked` -16.7%,
    `remainder_wide` -17.0%, `softplus`/`logaddexp` -15.9%,
    `cbrt_accurate` -15.1%, `sin_checked` -19.0%, `tanh` -13.9%), but a
    real minority regressed (`atan2` +19.5%, `acosh` +10.0%, `sinc`
    +9.8%, `cospi` +8.7%, `erf` +8.6%), and one function regressed
    **catastrophically**: `pown` 3.805->10.103 cyc/elem, **+165.6%**.
    Since `-C llvm-args` is a whole-crate `RUSTFLAGS` setting with no
    per-function scoping available in stable Cargo, there's no way to
    keep the broad wins while avoiding `pown`'s collapse short of a much
    more invasive build setup (per-module compilation units with
    different flags, or an opt-in Cargo feature/profile) -- out of scope
    for adopting this as the crate's actual default. Not adopted; the
    idea's own premise ("the default isn't always right") is confirmed
    true for roughly a third of the functions in this crate, but "isn't
    always right" cuts in a direction (helps most, devastates one) that
    makes a uniform crate-wide flag change a net loss for `pown`'s own
    users specifically, not a free lunch. *A compiler-flag-level
    experiment's "aggregate" result (most functions better) can hide a
    severe per-function regression that a per-function code change would
    never get away with shipping -- when a change can't be scoped to
    just the functions that benefit, a single catastrophic outlier
    (here, +165.6% on `pown`) should veto adoption even if the majority
    of the suite improves, unless the crate is prepared to accept that
    tradeoff explicitly (e.g. via an opt-in build profile, not the
    default).*
73. **PGO (+ BOLT) probe on the bench binaries (tried 2026-07-09,
    inconclusive with the tools/environment available -- matches the
    idea's own predicted-null framing, but for a different reason than
    expected)**: ran a real, full PGO pipeline on `quickbench` (build
    with `-C profile-generate`, execute to collect `.profraw`, merge
    with `llvm-profdata`, rebuild with `-C profile-use`) rather than just
    reasoning about it. Two real obstacles surfaced, both worth knowing
    for future PGO attempts here: (1) wall-clock (`quickbench`) is too
    thermally noisy on this machine to trust a PGO-vs-baseline
    comparison at all -- the *same* PGO binary run twice back-to-back
    showed `powf_checked` throughput swing from 6.183 to 11.935 ns/op
    (a ~2x difference between two runs of identical code), already
    documented elsewhere in this readme as this CPU's own known
    thermal-throttling behavior (~2.5x mid-session), which completely
    swamps whatever real PGO effect might exist. (2) This crate's
    primary, trusted measurement tool (`mca.rs`'s `llvm-mca` static
    analysis) can't evaluate PGO at all even in principle -- PGO's real
    lever is profile-guided *inlining and branch/block layout decisions
    across the whole compiled program*, not the isolated single-function
    assembly region `mca.rs` extracts and feeds to `llvm-mca`
    one region at a time. Disassembling the whole PGO vs. non-PGO
    `quickbench` binary did show a real, deterministic difference (PGO
    binary ~5.6% fewer total disassembled lines, consistent with
    different inlining choices across the many monomorphized bench
    closures), confirming PGO *did* change something -- just not
    something either of this crate's own tools can currently attribute
    to specific functions or trust as a genuine speed verdict. Cleaned
    up all scratch PGO artifacts (`.profraw`/`.profdata`/temp binaries),
    no code or config changes kept. *The idea's own "likely a null
    result" prediction held up, but the actual reason is more interesting
    than "branchless code has nothing for PGO to grab" -- it's that this
    crate's own toolchain (a static single-region analyzer plus a
    thermally-noisy wall-clock harness) genuinely can't produce a
    trustworthy verdict on a whole-program optimization technique like
    PGO at all, regardless of whether PGO itself would help. A future
    attempt would need either a quieter benchmark environment (fixed
    CPU frequency, isolated core) or a different measurement approach
    entirely (e.g. `perf stat` cycle counts averaged over many runs)
    before PGO could be fairly judged here.*
74. **codegen-units=1 + lto sweep for the bench profile (resolved
    2026-07-09, no artifact-boundary issue)**: rebuilt the entire mca
    suite (~70 functions) with `CARGO_PROFILE_RELEASE_CODEGEN_UNITS=1`
    and `CARGO_PROFILE_RELEASE_LTO=fat` (env override, no Cargo.toml
    change) and diffed the full output against a freshly-captured
    default-profile baseline from the same session -- bit-for-bit
    identical everywhere except one function off by 0.001 cyc/elem
    (`powf_checked` 9.105 vs 9.104, plainly floating-point rounding noise
    in the mca cycle-count arithmetic itself, not a real codegen
    difference). Confirms the existing default profile (codegen-units=16,
    no LTO) isn't splitting anything across a compilation-unit boundary
    in a way that costs real inlining/optimization -- this crate's
    existing mca numbers are already trustworthy as measured, no LTO
    needed to get an honest picture. No code/config change.
75. ~~**Asm-grep CI test**~~ (mostly already existed; gap closed
    2026-07-10): `examples/codegen_check.rs` already asserted zero `call`
    instructions and zero `cvttsd2si`/`cvttss2si` (the saturating-cast
    de-vectorization class), plus confirmed at least one packed op is
    present per region -- but that last check only proves "some packed
    arithmetic exists," not "no scalar sqrt/div coexists alongside it"
    (a partial de-vectorization of just one sub-computation could hide
    behind an otherwise-packed region and pass silently). Checked
    directly: zero regions currently have this (verified via a manual
    assembly scan scoped to genuine `_throughput`-suffixed regions only,
    distinguishing them from same-named `_latency` regions where scalar
    ops are expected and fine), but nothing was actually asserting it.
    Added an explicit `vsqrtss`/`vdivss`/`vsqrtsd`/`vdivsd` check
    alongside the existing two; all 71 regions still pass. Zero-risk,
    zero perf/accuracy effect (test-tooling only) -- closes the gap this
    backlog entry asked for.
76. ~~**Continue the AVX-512 probe**~~ (resolved 2026-07-10 -- real,
    large, near-universal win found, but not adoptable as a blanket
    default): `examples/scratch_avx512_probe.rs` (the untracked scratch
    file this idea pointed at) didn't actually test anything -- no
    timing, no comparison, just confirmed a `#[target_feature(enable =
    "avx512f")]` wrapper around plain `exp2` compiles and runs. Deleted
    it in favor of a real measurement, and pursued a more fundamental
    question than the idea's own literal ask (explicit `core::simd`
    f32x16 tiers): **is the existing auto-vectorizer already using the
    widest vector width this hardware supports?** `lscpu`/`rustc --print
    target-features` confirm this machine (an 11th Gen Intel Core
    i5-1145G7, "Tiger Lake") fully supports `avx512f`/`avx512dq`/
    `avx512cd`/`avx512bw`/`avx512vl`, and `.cargo/config.toml` already
    builds everything with `-C target-cpu=native`. But grepping the
    compiled `mca_target.s` showed every `*_throughput` region uses only
    256-bit `ymm` registers (e.g. `exp2_throughput`: 92 `ymm`, 0 `zmm`)
    -- LLVM's x86 backend defaults to a `prefer-256-bit` target feature
    on this CPU (a real, if `rustc`-flagged-unstable, target feature;
    confirmed via `rustc --print target-features`), overriding the
    hardware's own 512-bit capability. This is a genuinely different
    knob from idea #72's already-rejected `-force-vector-interleave`
    sweep -- that one tested unroll factor at a *fixed* vector width,
    this tests the width itself.

    Rebuilt the whole mca suite with
    `RUSTFLAGS="-C target-cpu=native -C target-feature=-prefer-256-bit"`
    -- confirmed via grep this actually flips codegen to `zmm` throughout
    (7015 `zmm` register uses in the full `.s` file, up from 0), and
    diffed every throughput number against the default build. Result:
    a **dramatic, near-universal win** -- of 70 measurable throughput
    rows, 60 improved by more than 1% (median -19.4%, mean -15.2%,
    several past -30%: `powf_checked_unchecked` -42.3%, `powf` -38.0%,
    `exp10_checked` -35.7%, `remainder_wide` -35.2%, `powf_checked`
    -33.8%, `cos_checked`/`logaddexp`/`softplus` all past -32%), 5 were
    flat (within +-1%: `rhypot`, `fmod_unchecked`, `rsqrt`, `expm1`,
    and one borderline), and only 2 regressed meaningfully:
    `remainder_unchecked` (+5.9%) and `asin` (+6.0%) are minor, but
    `exp_checked` (1.729->3.255 cyc/elem, **+88.3%**) and `pown`
    (3.805->7.022, **+84.5%**) are severe. (Latency numbers are
    completely unaffected either way, as expected -- `mca`'s latency
    chain measures the scalar `*_normal` core, never vectorized.)

    Root-caused the two severe regressions with `llvm-mca
    -resource-pressure` on the isolated regions rather than guessing:
    under the default (`ymm`) build, `exp_checked_throughput`'s
    FMA/mul/add work splits evenly across two ports (`ICXPort0` 21.17,
    `ICXPort1` 21.17 pressure/iteration -- `llvm-mca`'s scheduling model
    reuses Ice Lake's "ICX" resource names for Tiger Lake, its nearest
    documented relative). Under the `zmm` build, the *same* arithmetic
    collapses almost entirely onto `ICXPort0` alone (19.05, vs. `ICXPort1`
    at just 4.92) -- confirming this CPU has two independent 256-bit FMA
    units (one per port) but only *one* of them widens to handle a full
    512-bit FMA; the other sits nearly idle once every op is 512-bit
    wide. So going 256-bit->512-bit does *not* double this CPU's peak FMA
    throughput the way it might on hardware with two genuine 512-bit FMA
    units -- it only removes overhead (fewer, denser instructions: AVX-512's
    embedded per-element broadcast, e.g. `vfmadd231ps
    .LCPI23_4(%rip){1to16}, %zmm2, %zmm0`, folds what used to be a
    separate `vbroadcastss` into the arithmetic op itself, and the whole
    16-element array is now one pass instead of two manually-unrolled
    8-wide copies). Functions dominated by *overhead* (broadcasts, extra
    movs, loop duplication) net-win big from that reduction; functions
    whose bottleneck was already raw FMA-port throughput itself (`pown`'s
    long compile-time-unrolled integer-power multiply ladder;
    `exp_checked`'s unusually dense poly+range-reduction chain) lose the
    second port for zero compensating benefit and net-regress hard.
    Confirms and fully explains (not just observes) idea #72's own
    finding that `pown` is uniquely vulnerable to vector-width/layout
    changes on this hardware -- same function, same underlying
    single-512-bit-FMA-port constraint, different trigger.

    **Not adopted.** Even though this is a *far* stronger case than idea
    #72's own mixed 60/40 split (here: 60 real wins, only 2 severe
    losses, both now understood), the same practical blocker applies:
    `-C target-feature` is a whole-crate `RUSTFLAGS` setting with no
    per-function scoping on stable Cargo, and this crate's own precedent
    (idea #72) already established that a single catastrophic per-function
    outlier vetoes adopting a build-wide flag as the default, since
    `pown`'s (or here, also `exp_checked`'s) users would eat an 85-88%
    regression with no opt-out. Additionally: `llvm-mca`'s static cycle
    model has no way to represent Tiger Lake's real, separate AVX-512
    frequency/license downclocking behavior under *sustained* wide-vector
    load (a genuinely different, additional risk from the port-contention
    story above, which the model *does* capture correctly) -- so even the
    60 modeled wins could be smaller in practice than shown here, the
    same category of unverifiable-with-this-crate's-tools caveat idea
    #73's PGO probe already flagged for a different technique. No
    `.cargo/config.toml` or `src/lib.rs` change made; all scratch
    RUSTFLAGS-built binaries/assembly cleaned up. Left open as a concrete,
    well-quantified opportunity for whoever wants to build the
    infrastructure this would actually need: either a real opt-in Cargo
    build profile/feature (accepting `pown`/`exp_checked` regress for
    users who don't call them), or waiting for a stable per-function
    scoping mechanism finer than whole-crate `RUSTFLAGS`. *A whole-crate
    codegen-flag experiment needs the same "one bad outlier vetoes the
    default" discipline as any other crate-wide change, no matter how
    lopsided the win/loss ratio looks in aggregate -- but a decisive
    majority-win result is still worth fully root-causing (not just
    measuring) when the mechanism might explain a previously-mysterious
    related finding, as it did here for idea #72's own `pown` outlier.*
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

79. ~~**Survey sinpi/cospi/sind/cosd/exp10/sigmoid max ulp**~~ (resolved
    2026-07-10): ran a genuine `thorough` (exhaustive, every f32 bit
    pattern in each function's own documented domain) sweep for all six,
    rather than trusting that readme.md's existing numbers were already
    exhaustive rather than leftover fuzz-mode (100M sample) estimates.
    Five of six matched their documented numbers exactly: `sinpi`
    (0.1969/2 vs. documented 0.197/2), `cospi` (0.2813/868814811 vs.
    0.281/8.7e8, the already-known near-a-zero artifact), `sind`
    (0.1237/2 vs. 0.124/2), `cosd` (0.0725/2 vs. 0.073/2), `exp10`
    (0.0343/2 vs. 0.034/2) -- confirming the existing numbers really were
    already the true exhaustive worst case, not stale fuzz estimates,
    despite never having been explicitly re-verified as such.

    `sigmoid` was the exception: exhaustive measured 0.0925/4, *tighter*
    than the documented 0.100/5 -- a discrepancy in the safe direction
    (the doc overstated the error), but still wrong and worth tracing.
    Since a random-sample fuzz can never find a worse max than the true
    exhaustive sweep (fuzz here literally uses `rand::rng().random::<u32>()`,
    a raw uniform subset of the same bit-pattern space `exhaustive()`
    enumerates completely -- see `examples/accuracy.rs`'s `fuzz` fn), the
    documented "5" had to predate some real change. Root-caused via `git
    log -S"pub fn sigmoid"`, which initially (and misleadingly) showed no
    hits -- that pickaxe search only catches the function being added or
    removed, not its body being edited while the signature stays put.
    Broadening to any commit touching `sigmoid` found it immediately:
    commit `2632e14` ("sigmoid: single exponent-field construction
    instead of exp's k1/k2 split", 2026-07-09) already measured and
    reported this *exact* number (0.0925/4, bit-identical to the original,
    not a new approximation) as part of a real, already-shipped speed win
    (mca throughput 2.713->1.354 cyc/elem, halved) -- but that commit's
    diff touched only `src/lib.rs`, never `readme.md`. The mca table's own
    sigmoid row *was* correctly resynced at the time (still reads
    61.09/1.354, matching), so this was a partial doc-sync miss: one
    table updated, the sibling accuracy table forgotten, the same failure
    shape as this session's earlier erf/log1p/atanh/sinpi/tanpi/sinc mca
    staleness findings, just on the accuracy side instead of mca. Fixed
    readme.md's sigmoid accuracy row to `0.093 | 4`. *A pickaxe search for
    a function's own signature line is not sufficient to rule out "this
    function's body changed" -- it only proves the function wasn't added
    or removed. When auditing whether a specific number could be stale,
    search for any commit touching the function by name (or grep the
    commit log directly), not just the declaration line.*
80. **exp2 poly evaluated as 1+f·Q vs direct P(f)=2^f with c0=1 pinned
    (paper-screened 2026-07-09, not implemented -- real modest op savings
    identified, but needs a fresh fit, not a mechanical rewrite)**: traced
    through the actual op count for both forms. Current (`Q(f)=(2^f-1)/f`
    fit, reconstructed as `exp2int*(1+f*Q(f))`): the final combine is
    `fma(q, exp2int*f, exp2int)`, needing a separate `exp2int*f` multiply
    *before* the final fma (2 ops for the combine). A direct `P(f)=2^f`
    fit (same degree-5 shape, `c0` pinned to exactly `1.0` instead of a
    fitted constant near `ln(2)`) would combine as a single `p*exp2int`
    multiply (1 op) -- a real, if modest, saving of one multiply out of
    ~8 total ops in the poly+combine chain. But this is *not* a
    mechanical rewrite of the existing coefficients: algebraically
    expanding `P(f) = 1 + f*Q(f)` from `Q`'s own already-fitted
    coefficients raises the degree by one (a degree-5 `Q` gives a
    degree-6 `P`, more terms, not fewer) -- a same-degree direct `P(f)`
    fit needs its own independent minimax fit (via `tune.rs` or
    `lolremez`), the same real numerical-fitting effort as this session's
    other coefficient refits. Given the idea's own caution held up on
    inspection (`exp2`/`exp2_checked` are already at max ulp 1, avg
    ~0.03 -- essentially no room to spare, so a fresh fit risks a real
    regression for a payoff capped at one multiply), and no existing
    fitted coefficients could be reused directly, not undertaken this
    iteration -- the effort (a full fresh fit) vs. payoff (one multiply)
    ratio is comparable to idea #52's own "not pursued" call. Left open
    for a session that wants to invest in the fresh fit specifically.
81. ~~**sinf_poly's copysign(x)**~~ (resolved 2026-07-10 -- see the
    "sinf_poly copysign audit" entry in the "sin_checked / cos_checked
    internals" section near the top of this file for the full writeup,
    including the much bigger bug this audit's own real-scale
    verification surfaced as a side effect): `sin_checked`/`sinpi` are
    provably redundant with it (split into `sinf_poly_raw`); `sin`/`sind`/
    `cospi`/`cosd` still genuinely need it.
82. **log1p at x=+1.0 boundary (resolved 2026-07-09, no bug)**: checked --
    `log1p(1.0)` is bit-exact (`c` is exact at `u=2.0`, right at Sterbenz's
    inclusive boundary), and a dense sweep either side of `x=1.0` shows
    max ulp 1, matching log1p's own documented budget.
83. **hypot_checked denormal-pair path (resolved 2026-07-09, no bug)**:
    20M-sample fuzz of denormal x denormal pairs (both signs) against an
    f64 reference: max ulp 1. The Python prototype's broad sampling wasn't
    thin here after all.
84. **erf/erfc NaN sign convention audit (resolved 2026-07-09, no action)**:
    checked directly -- several mulsign-based paths (erf, atan2, atan,
    asin) do preserve a negative-NaN input's sign bit through to the
    output, while std canonicalizes to a positive NaN. Confirmed this is
    NaN sign/payload propagation, which IEEE754/C99 leaves
    implementation-defined (not a spec violation, matching the idea's own
    "C99 doesn't care" framing) -- fixing it to mimic std's canonicalization
    would need a real `.abs()`-style op added to every mulsign-based
    NaN-producing path, for zero standards-compliance benefit. Left as is.
85. **atan2 full C99 special-case matrix (fixed 2026-07-09)**: built the
    matrix (all zero/inf/nan/sign combinations against std) -- found a
    real bug, not just a coverage gap: `atan2(NaN, 0.0)`/
    `atan2(NaN, -0.0)` returned `+-FRAC_PI_2` instead of `NaN` (the
    `x==0` branch bypasses `atan(y/x)` and falls to a bare `mulsign`,
    which only reads `y`'s sign bit, not its NaN-ness). Every other NaN
    combination already worked. Fixed with a trailing override; 7 new
    edgecheck pins lock the matrix down going forward. Commit `5732abc`.
    Immediate follow-up (same day): applied the same matrix technique to
    `hypot`/`hypot_checked`/`rhypot` (clean, no bugs) and to
    `remainder`/`remainder_checked`/`remainder_ieee`/`remainder_wide`/
    `fmod` -- found the *same* NaN-discarding shape recurring: all five
    share `if x==0.0 {x} else {normal}`, and `normal` already correctly
    evaluates to NaN whenever `y` is `0` or NaN, but the guard fired
    unconditionally anyway, silently overriding it. `remainder(0,0)`,
    `fmod(0,0)`, etc. all returned `0`/`-0` instead of NaN, directly
    contradicting `fmod`'s own doc comment (claims to match Rust's `%`
    exactly, and `0.0f32 % 0.0f32` is NaN). Fixed uniformly (`&&
    !normal.is_nan()` added to the guard) across all five. Commit
    `bc81831`. Fourth wave (same day): built the full matrix for
    `powf`/`powf_checked` too (242 combinations) -- found `powf_checked`
    dropping NaN `y` the same way `atan2` dropped NaN `y` (a `y > 0.0`
    comparison silently false for NaN), plus a real, separate design
    gap: `x == -0.0`/`x == -inf` were treated identically to a
    genuinely negative *finite* `x` (NaN for non-integer `y`), but C99
    exempts them -- only odd-integer `y` preserves their sign, any
    other `y` gives the unsigned magnitude. A third rule: `y` infinite
    always gives an unsigned result regardless of `x`'s sign. Fixed all
    three in both functions; matrix now matches std exactly (was 15
    mismatches). Commit `47bc7c1`.
    Fifth wave (same day): applied the matrix to `sinh`/`cosh` -- found
    `sinh`/`cosh(+-inf)` both return `NaN` instead of `+-inf`/`inf`, and
    for large finite `|x|` (e.g. `1000`) they wrap around to a
    wrong-sign or garbage finite value, because `exp_pos_neg` (their
    shared `exp(x)`/`exp(-x)` helper) has no domain clamp at all, so the
    magic-round exponent-field split silently wraps outside its own safe
    range instead of saturating. Added `sinh_checked`/`cosh_checked`
    (full-range siblings, mirroring `exp_checked`'s existing pattern)
    built on a new `exp_pos_neg_checked_half` helper. Fixing this
    surfaced a second, larger bug along the way: `sinh`/`cosh` apply
    their `0.5` scale factor *after* computing the full `exp(x)`, so the
    intermediate overflows to `inf` right around `x~88.7` even though
    the true (halved) answer is still finite up to `x~89.4` -- invisible
    to the existing accuracy.rs sweep because plain `sinh`/`cosh` use
    their own `sinh_domain` closure that deliberately excludes this
    exact window. Fixed by pushing the `0.5` into the exact-power-of-two
    field split itself (`t1 * 0.5`, exact for any power-of-two float)
    before the final multiply, rather than scaling the final result.
    Also found, via a standalone exhaustive probe (hand-deriving the
    bound got too intricate to trust): the field split's own safe `k`
    range is `[-254, 254]`, not the `[-128, 128]` initially assumed by
    naively applying `exp_checked`'s own asymmetric overflow bound
    symmetrically to both `+k` and `-k` -- `exp_checked`'s `88.7228`
    bound is where `exp(x)` *alone* overflows, not where the split's bit
    trick mechanically breaks. Landed on a `170.0` clamp (`k` up to
    ~245.3), comfortably inside the proven-safe window. Verified:
    fuzz-mode accuracy clean (`sinh_checked` avg ulp 0.0421 max 5,
    `cosh_checked` avg ulp 0.0373 max 5, matching `sinh`/`cosh`'s own
    baseline), `codegen_check` clean (no scalar fallback), 16 new
    edgecheck pins (`+-inf`, `+-1000`, the `89.415`/`89.416` boundary).
    mca: `sinh_checked` 58.06/2.527 cyc (latency/throughput) vs `sinh`'s
    56.00/2.089; `cosh_checked` 58.06/2.212 vs `cosh`'s 55.00/1.754 --
    modest overhead for full-range correctness, same tradeoff shape as
    `exp_checked` vs `exp`. Commit `d3b99ac`.
    Sixth wave (2026-07-09, later the same day): applied the matrix to
    `softplus`/`logaddexp` -- clean, no bug found. Full 9x9 grid over
    `{0, -0, +-1, +-inf, NaN, +-100}` for `logaddexp` (81 combinations)
    plus the 1-argument set for `softplus`, checked by hand against known
    closed forms (`logaddexp(0,0)=ln(2)`, `logaddexp(100,100)=100+ln(2)`,
    `logaddexp(1,-1)=ln(e+1/e)`, etc.) and against the "infinity beats
    everything except NaN" pattern expected for a log-sum-exp: every
    `+-inf`/finite combination, every `+-inf`/`+-inf` combination
    (including the `inf-(-inf)=inf` and `(-inf)-inf=NaN`-then-`min`-
    discards-NaN internal path for `logaddexp(+-inf,+-inf)`'s own `d`
    computation) all resolved to the mathematically correct answer
    despite an internal NaN intermediate in two of the four `+-inf,+-inf`
    cases -- `f32::min`'s NaN-discarding semantics (the same mechanism
    `softplus`'s own doc comment already flags as a fixed hazard
    elsewhere) happens to route around itself harmlessly here, verified
    directly rather than assumed. Any input `NaN` correctly propagates
    to `NaN` unconditionally (by this function's own design choice, not
    a C99-mandated "infinity wins" exemption the way `hypot`/`atan2` have
    -- `logaddexp` isn't a standard function, so there's no external
    convention being violated either way). No code change.
    Seventh wave (2026-07-09, later still): `asinh`/`acosh`/`atanh`
    against `x in {0, -0, +-1, +-inf, NaN, +-0.5, +-2}` -- all clean,
    matching an f64-computed reference exactly except one pair
    (`asinh(+-0.5)`) off by exactly 1 ulp, comfortably inside `asinh`'s
    own documented ~0.15-avg/3-max budget, not a bug. Domain-boundary
    cases all correct: `acosh(x<1)` (including `-1`, `-0`, `0`) is `NaN`,
    `acosh(1)=0`, `atanh(+-1)=+-inf`, `atanh(|x|>1)=NaN`,
    `atanh(+-inf)=NaN`. No code change.
    Eighth wave (2026-07-09, later still): `erf`/`erfc`/`erfcx` and
    `asin`/`acos`/`atan`, same value set plus `+-10`/`+-0.28`. All clean
    -- `erf`/`erfc`/`erfcx` consistent with behavior already verified
    when `erfcx` was implemented earlier this session (`erfc(+-inf)`'s
    own clamped-rational bound, `erfcx`'s extrapolation-past-fit-domain
    caveat), and `asin`/`acos`/`atan` came back **bit-exact** against an
    f64 reference at every tested point, no even-1-ulp misses -- expected
    given both families already went through extensive dedicated
    refitting earlier this session (asin's own 8-fix history, atan's
    poly refits). No code change. Given three consecutive clean waves
    now across most of the remaining function families, the special-
    case-matrix technique has reached diminishing returns for this
    session -- worth switching to a different technique for the next
    idea rather than continuing exhaustive sweeps of already-hardened
    functions.
88. **exp10 near the decade boundaries (resolved 2026-07-09, no bug found)**:
    densely fuzzed (12M samples) right around every point where
    kr=round(x·log2_10) crosses an integer (where the floor-adjust select
    flips), plus every crossing point exactly -- clean, max ulp 1 both
    ways, matching exp10_checked's own documented budget. The
    floor-adjust logic is correct at its own boundaries; no sinpi/cospi-
    style bug here.

### Longer shots / research-flavored

89. **Bit-sliced two-for-one**: evaluate sin and cos polynomials sharing
    y=r² registers across the *same* vector when the caller wants both —
    a sincos slice API (not scalar API, which already failed) where lane
    pairing amortizes the reduction. Only viable inside a slice tier.
90. **Newton-free correctly-rounded sqrt-composites (tried 2026-07-09 for
    rsqrt, rejected -- real accuracy win, real (if modest) cost, plus a
    genuinely new codegen pitfall found)**: implemented `rsqrt`'s
    residual step exactly as described -- `e = fma(r, r*x, -1)`,
    correction `r_new = fma(-0.5*r, e, r)`. Two real bugs found before
    it worked at all: (1) computing `r*r` first (the naive reading of
    "r·r·x") overflows f32 for any x small enough that `r=1/sqrt(x)`
    itself exceeds ~1.8e19 (e.g. a denormal `x=1.95e-43` gives
    `r~2.27e21`, `r*r~5e42 > f32::MAX`), even though the true residual is
    tiny -- fixed by computing `r*x` first (~=`sqrt(x)`, always
    well-behaved) and multiplying by `r` after, same "which
    multiplication order avoids a needless overflow" lesson as idea
    #56's `cbrt_accurate` Halley investigation, tried immediately before
    this one in the same session. (2) The correction degrades to a real
    `0*inf=NaN` indeterminate form at `x==0`/`x==inf` (r is `+-inf`/`0`
    there) -- both cases `rsqrt`'s own bare `1.0/x.sqrt()` already gets
    exactly right for free via plain IEEE754 semantics, so this is a
    real regression, not a pre-existing gap. Fixing it with the obvious
    guard, `if x > 0.0 && x.is_finite() { corrected } else { r }`,
    fixed correctness but caused a *catastrophic* throughput regression
    (1.381->7.076 cyc/elem, +412%) -- inspecting the actual generated
    assembly found why: `x.is_finite()` (or this specific compound
    condition) compiles to a fully scalar per-lane sequence (extract
    each of the 8 lanes with `vextractps`, run ~10 scalar integer
    test/sub/cmp/set instructions per lane to compute "positive and
    finite," then hand-assemble an AVX-512 mask register bit-by-bit via
    a chain of `kmovd`/`kshiftlb`/`kshiftrb`/`korb`/`kandb`) instead of a
    single vectorized compare -- a genuinely new de-vectorization
    pitfall this crate's `codegen_check` doesn't currently catch at all
    (it only greps for a scalar `call` or `cvttsd2si`/`cvttss2si`, not
    this "scalar bit-tests reassembled into a mask" pattern). Rewriting
    the guard with this crate's own established bit-trick idiom instead
    (`ax = x.to_bits() & !SIGN_MASK; if ax==0 || ax>=EXPONENT_MASK {r}
    else {corrected}`, the same shape `cbrt_accurate`'s own zero/inf/nan
    guard already uses) fixed the vectorization completely: throughput
    back down to 1.502 (only +8.8% over plain `rsqrt`'s 1.381), all
    special cases verified correct again. With everything fixed:
    real accuracy win (avg ulp 0.2599->0.1226, ~2x tighter, max ulp
    unchanged at 1) at a real, modest cost (latency 28.00->40.03 cyc,
    +43.0%; throughput +8.8%). Didn't clear this loop's own bar (speedup,
    or accuracy gain *without* a perf penalty) -- both axes show a real
    cost, even though it's far short of the `is_finite()` version's
    catastrophic one. Also: `rsqrt`'s accuracy is already extremely
    tight (max 1 ulp is close to the practical ceiling for a
    two-composed-hardware-ops function), so the accuracy gain here is a
    nice-to-have polish, not fixing a documented defect -- unlike the
    "real perf cost to fix a genuine correctness gap" cases this crate
    has accepted elsewhere (sin/cos's inf-for-large-x fix, tanh's
    domain-hole fix), there's no defect being fixed, just squeezing an
    already-excellent number further. Reverted, bit-identical to prior
    HEAD. Not tried for `rhypot` this round (same technique, likely a
    similar shape of result -- left for a future session if the
    tradeoff calculus differs there). *Any time a domain guard is added
    to protect a numerically-motivated correction (division-by-zero,
    inf*0, etc.), check the *generated assembly* for the guard itself,
    not just its correctness -- `.is_finite()` (or compound conditions
    built from it) can silently de-vectorize into dozens of scalar
    per-lane instructions even when no `call` or saturating-cast pattern
    is present, a failure mode this crate's own `codegen_check` doesn't
    currently detect; prefer this crate's own established bit-trick
    idiom (`x.to_bits() & !SIGN_MASK` compared against `0`/`EXPONENT_MASK`)
    for exactly this class of zero/inf/nan guard, which is already known
    to vectorize cleanly everywhere else it's used.*

    **Follow-up (same day): swept the crate for this exact pattern in
    already-shipped functions.** `powf_checked`'s own `ax > 0.0 &&
    ax.is_finite()` guard (a compound `&&` of two conditions, the same
    shape that broke here) turned out to have the *identical* live
    regression -- confirmed via the actual assembly (`vextractps` +
    scalar int test/cmp/set chains + `kshiftlb`/`kshiftrb`/`korb`
    mask reassembly, all present) before fixing it with the same
    bit-trick rewrite: bit-identical accuracy (verified via fuzz +
    edgecheck, all special cases unchanged), real mca win, throughput
    9.293->9.105 cyc/elem (-2.0%), latency 130.67->130.33 (-0.3%).
    Commits `3d9c67d` (code), `69164a5` (readme sync). Also checked
    `remainder`/`remainder_ieee`/`remainder_checked`/`remainder_wide`/
    `fmod`, which all share an *identical-looking* `y.is_infinite() &&
    x.is_finite()` compound guard -- a naive broad grep for
    `vextractps|kmovd|kshiftlb|kshiftrb|korb|kandb` suggested
    `remainder_checked`/`remainder_wide` had the same issue (4 matches
    each) while the other three didn't (0 matches), but checking
    specifically for `vextractps`/`kshiftlb` (the actual catastrophic
    signature, not just any `kmovd`) showed **zero** in all five --
    the same textual guard pattern does *not* uniformly trigger the bad
    codegen, apparently depending on surrounding context LLVM sees
    (register pressure, nearby ops) rather than the condition's own
    shape alone. Applied the bit-trick rewrite to
    `remainder_checked`/`remainder_wide` anyway on the (wrong) assumption
    the broad grep was meaningful, then measured mca before committing
    (per this session's own discipline) -- numbers came back **bit-for-
    bit identical** to the pre-rewrite baseline (46.13/1.287 and
    179.20/8.876, exactly matching readme.md's existing values),
    confirming these two never had the expensive pattern at all; the
    broad grep's "4 matches" were benign, ordinary AVX-512 mask usage
    unrelated to this bug. Reverted that unnecessary rewrite (kept only
    the two confirmed-beneficial fixes above) rather than leave in
    complexity with no measured payoff. *A grep for `kmovd`/`kshiftlb`-
    family mnemonics alone is not a reliable detector for this
    de-vectorization class -- always confirm the specific
    `vextractps`-plus-mask-reassembly *chain* is present (and ideally
    confirm via an actual mca before/after, not just instruction-count
    grepping) before spending effort "fixing" a guard that was already
    vectorizing fine.*
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
97. **exp/exp2 denormal-output double-rounding (resolved 2026-07-09, no
    bug found)**: swept `x in [-149,-126)` (the exponent range whose
    output lands in `[2^-149,2^-126)`, i.e. the denormal zone) against a
    correctly-rounded f64 reference, both a dense linear sweep and a
    random-bit-pattern sweep -- `exp2_checked` (the relevant full-range
    variant `powf_checked`'s own residual max-ulp neighborhood per #58
    actually uses) came back clean: avg ulp 0.0217, max ulp 1, no
    evidence of a double-rounding defect from the `t2` multiply. Plain
    `exp2` (unchecked) does produce garbage there (avg ulp in the
    billions), but that's expected, already-documented pre-existing
    behavior, not a new bug -- its own doc comment explicitly states
    "valid for x in [-126,128)... outside that range the exponent
    construction wraps around and the result is garbage," and denormal
    outputs fall well outside that stated domain by design (use
    `exp2_checked` there, which this sweep confirms does the right
    thing). No code change.
98. **Karatsuba-style Df32 multiply (audited 2026-07-09, premise doesn't
    apply -- found and removed genuinely dead code instead)**: did the
    requested audit first. The full `Df32*Df32` multiply (`impl Mul for
    Df32`, the "2-3 two_prods" this idea targets) turned out to have
    **zero call sites anywhere in the crate** -- `powf_checked`'s actual
    chain (`log2_df(ax) * y`, `exp2_checked_df`'s own `v.1*LN_2`
    correction) only ever multiplies a `Df32` by a plain `f32` (the
    already-cheap 3-op `Mul<f32> for Df32` impl), never `Df32` by
    `Df32`. Confirmed via exhaustive grep across `src/` and `examples/`
    (only comment mentions, no actual type-level usage) --
    `doublefloat.rs` is a private module with no external consumers
    either, so this wasn't reachable any other way. So the idea's own
    premise (cheapen the Df32*Df32 multiply to help powf_checked) can't
    apply: that code path doesn't exist in the powf chain at all, and
    powf_checked's own real cost lives elsewhere (`log2_df`'s own poly,
    `exp2_checked_df`'s field split). Removed the now-confirmed-dead
    `impl Mul for Df32` entirely (12 lines) rather than leave it as
    unreachable complexity -- verified behavior-unchanged via
    `cargo test`, `accuracy.rs` (powf/cbrt/remainder families, all
    unchanged), and `edgecheck.rs` (0 failures). Commit `3c375f3`.
    *Auditing a target function's actual call graph before optimizing a
    piece of shared infrastructure it's assumed to use is worth doing
    first every time -- here it revealed the assumed dependency doesn't
    exist at all, redirecting the useful outcome from "cheapen a hot
    multiply" to "delete an unreachable one."*
99. **Precision-tapered polys (tried 2026-07-09 for log_2's full poly,
    rejected -- real accuracy win, dramatic perf cost)**: tested the
    simplest version of this idea first -- replace `log_2_normal`'s
    entire degree-9 Estrin-balanced tree (2 muls + 9 fma's) with a plain
    Horner chain (9 fma's, no `s2`/`s4` precompute at all), rather than
    the idea's own narrower "only the l3/l4 tail" scoping. Accuracy was a
    genuine, if modest, surprise win: avg ulp 0.0031->0.0019, max ulp
    3->2 (both improved, not a tradeoff) -- fewer total operations
    apparently gives the rounding fewer chances to compound, at least for
    this specific polynomial. But mca killed it decisively: latency
    34.23->53.05 cyc (+55.0%), throughput 1.584->2.249 cyc/elem (+42.0%)
    -- a dramatic regression on *both* axes despite genuinely fewer total
    ops (9 fma's vs. 11 ops). The fully serial 9-deep dependency chain
    (each fma must wait for the previous) costs far more than the 2 extra
    multiplies the balanced Estrin tree pays for its shorter critical
    path -- true even in throughput/vectorized mode, where a single
    call's own serial depth might be expected to matter less (many
    independent lanes should be able to fill the pipeline) but evidently
    still doesn't fully hide a 9-deep chain here. Reverted, bit-identical
    to prior HEAD. Didn't test the idea's own narrower proposal (Horner
    only for the small-magnitude l3/l4 tail, keeping Estrin for the
    dominant l0/l1 terms) empirically, but worked the algebra by hand
    afterward and it closes the question anyway: the tail's contribution
    (`l2*s4+l3*s6+l4*s8`, degree >=4) written as a Horner-evaluated
    sub-poly needs 5 serial fma's plus 1 more to combine with `r0` -- 6
    fma's total, *fully* serially dependent. The current Estrin tail
    (`l2`,`l3`,`l4` computed in parallel, then `r1`,`r2`,`p` combining
    them) is *also* exactly 6 fma's, but only 4 deep (3 independent
    leaves, then 3 dependent combine stages) instead of 6 -- same op
    count, strictly shorter critical path. So the narrower literal
    variant isn't an open question after all: it can only match or lose
    to the existing Estrin tail, never win, since there's no actual op
    reduction available once the reconstruction algebra is worked
    through (the earlier "fewer fma's" framing only held for the *full*
    poly, where Estrin's own `s2`/`s4` precompute -- 2 extra multiplies
    -- is the only place real ops are spent beyond what Horner needs;
    restricting to just the tail removes that comparison's own basis).
    Not pursued further; this fully closes idea #99 rather than leaving
    it open. *Total operation count is not a reliable proxy for vectorized
    throughput when the operations being removed also happen to shorten
    the critical dependency path -- a "fewer ops" restructuring that
    lengthens the serial chain can lose badly even in a throughput-
    oriented loop with many independent lanes, echoing this crate's own
    repeated finding elsewhere (asin's branch-count changes, the
    Newton/Halley correction-step rejections) that op-count alone rarely
    predicts mca's actual verdict.*
100. ~~**A cost model for "add a division"**~~ (resolved 2026-07-10):
    measured directly with a minimal hand-written asm probe fed straight
    to `llvm-mca -resource-pressure` (four independent instructions per
    region, both `ymm`/8-wide and `zmm`/16-wide, `-mcpu=native`) rather
    than reasoning from first principles. Results:

    | op | width | latency | RThroughput | cyc/elem |
    |----|-------|---------|-------------|----------|
    | vfmadd231ps/vmulps | ymm (8) | 4 | 0.50 | 0.0625 |
    | vfmadd231ps/vmulps | zmm (16) | 4 | 1.00 | 0.0625 |
    | vdivps | ymm (8) | 11 | 5.00 | 0.625 |
    | vdivps | zmm (16) | 18 | 10.00 | 0.625 |
    | vsqrtps | ymm (8) | 12 | 6.00 | 0.75 |
    | vsqrtps | zmm (16) | 20 | 12.00 | 0.75 |

    **Headline number: one packed division costs ~10x a single FMA/mul
    per element on this CPU (0.625 vs 0.0625 cyc/elem); one packed sqrt
    costs ~12x.** Paper-screening rule: replacing a division with N
    extra FMAs is a throughput win whenever N is meaningfully below 10
    (comfortable at N<=3-4, the shape every actual divider-idle win this
    crate has shipped -- cbrt's Newton reciprocal, sinh_throughput,
    log1p -- already has), a wash around N~8-10, and a loss past that.

    **A genuinely interesting second finding, width-invariance:** the
    divider/sqrt's own per-element cost is *identical* at ymm and zmm
    width (0.625 and 0.75 cyc/elem respectively, exactly, both widths) --
    confirmed via the resource table: a 512-bit `vdivps`/`vsqrtps`
    decomposes into *two* sequential 256-bit sub-uOps against the same
    `ICXFPDivider` resource (RThroughput exactly doubles alongside the
    doubled element count), unlike the FMA story below. **This means
    division-heavy code doesn't pay idea #76's zmm penalty at all** --
    going wider is a clean win there (fewer broadcast/overhead
    instructions, same per-element divider cost) with none of the
    downside that hit `exp_checked`/`pown`.

    **Refines idea #76's own FMA finding, doesn't contradict it:** re-ran
    the same probe for `vfmadd231ps` specifically (4 fully independent,
    mutually non-dependent FMA instructions per region, no shared
    registers) and found per-element *throughput* is actually equal at
    both widths too in this idealized case (0.0625 cyc/elem, matching
    ymm's dual-port 2-per-cycle rate against zmm's single-port 1-per-cycle
    rate exactly, since a zmm op does double the work per instruction).
    Confirmed via resource pressure that ymm's 4 independent FMAs split
    2.00/2.00 across `ICXPort0`/`ICXPort1`, while zmm's put all 4.00 on
    `ICXPort0` alone (`ICXPort1` completely idle) -- so idea #76's
    original "only one port handles 512-bit FMA" statement holds exactly.
    The missing piece idea #76 didn't fully spell out: in this idealized
    *independent* case, losing the second port costs nothing because a
    single zmm instruction already does 2x the ymm instruction's work, so
    "half the issue rate, double the work per issue" nets out even. Real
    functions regress specifically when their poly evaluation has a
    *long serial dependency chain* (Horner-style `p = fma(p, r, c)`
    repeated many times) -- the *default* ymm build's own double-unrolled
    codegen (two textually-separate copies of the whole chain, one per
    8-element half) gets *free* cross-copy parallelism from the
    scheduler interleaving copy A's and copy B's independent chains
    across both ports, fully masking each copy's own serial latency. At
    zmm width there is only *one* copy (all 16 elements share one
    register), so that free masking disappears entirely -- the single
    chain is now serially bottlenecked *and* stuck on one port, with
    nothing independent left to interleave. This is exactly why
    `exp_checked`'s dense poly and `pown`'s long squaring ladder regressed
    hard while broadcast/overhead-bound functions (most of the other 60)
    didn't: those functions' *chains* were short/shallow to begin with,
    so they were never relying on the free cross-copy masking in the
    first place, only paying for redundant broadcast/loop overhead that
    zmm genuinely eliminates.

101. **`pown`: real intermediate-overflow correctness bug for large `|n|`,
    found via structured search 2026-07-10, fix attempted and doesn't
    work as hoped, not fixed this round**: `pown` has no documented
    accuracy number for large `|n|` at all (readme.md only covers
    `|n|<=8` and `|n|<=64`, both already high -- 11 and 90 ulp
    respectively, clearly growing with range). Extended this session's
    structured worst-case mining to `pown`'s own (f32, i32) domain,
    targeting `n` values with many consecutive set bits (`2^k-1`, which
    force every one of `pown`'s 32 squaring iterations to actually
    contribute to `result`, maximizing the number of compounded rounding
    steps) crossed with `x` near 1.0 (so `x^n` can still land on a large
    but finite, representable value even for huge `n`). Found a real,
    reproducible **correctness bug**, not just an ulp gap: `pown(0.997296,
    -32767)` returns `inf`, but the true value (`3.4025991e38`) is
    finite and comfortably under `f32::MAX` (`3.4028235e38`). Traced
    exactly where: `pown` does exponentiation-by-squaring, 32 always-
    executed iterations, multiplying an accumulator by `base` whenever
    the corresponding bit of `n` is set and unconditionally squaring
    `base` every iteration. At iteration 14 (the last set bit for
    `n=32767`, whose binary form is 15 one-bits), the exact product of
    the accumulated `result` (`1.8443513e19`) and `base`
    (`1.849352e19`) is `~3.4109e38` -- genuinely over `f32::MAX` --
    while the *true* mathematically exact answer is `3.4026e38`,
    comfortably under it. **The gap between these two numbers is 14
    squaring steps' worth of compounded relative rounding error (each up
    to half a ulp, compounding multiplicatively)** landing the computed
    intermediate just over the overflow threshold when the true answer
    sits just under it -- not a wild miscalculation, a boundary case
    where accumulated imprecision is exactly the size of the gap between
    the true answer and the overflow cliff. Quantified how often this
    matters: of ~5.86M structured `(x, n)` pairs with a finite true
    answer (large all-ones `n` crossed with `x` near 1), ~10.2% came back
    as `pown` incorrectly returning `inf`/`0` -- a real, substantial rate
    *within this deliberately boundary-focused sweep* (not a claim that
    10% of all `pown` calls are broken; ordinary calls with small `|n|`
    or `x` far from 1 aren't near this cliff at all).

    **Tried the obvious fix and it doesn't work, for an instructive
    reason**: prototyped tracking `base`/`result` as `Df32` (double-float)
    through the squaring loop, with a hand-rolled `Df32*Df32` multiply
    (error-free two-product of the high words, cross terms in plain
    f32, matching this crate's own established double-float idiom).
    Result: **zero** of the 599,974 broken cases got fixed, and the
    specific traced example got *worse* (`NaN` instead of `inf`). Root
    cause of the prototype's own failure: `Df32` buys back *precision*
    (a second word to hold what a single f32 rounds away) but not
    *range* -- its primary (`.0`) word is still an ordinary f32 with the
    same `f32::MAX` ceiling, and `to_f32()` collapses back to
    `self.0 + self.1`, so if the primary term itself crosses the
    overflow boundary during the squaring chain, the extra precision in
    `.1` is irrelevant; the pair is already `(inf, something)` and stays
    `inf` once collapsed. **The real fix needs range extension (an
    explicit, separately-tracked exponent/scale, applied via integer
    arithmetic and reconstructed at the end via a bit-level scale, the
    same shape `remainder_wide`'s own `big`/`scale` rescue and
    `cbrt_accurate`'s large-magnitude rescale already use elsewhere in
    this crate), not precision extension** -- a genuinely different, more
    invasive fix than the one that worked for e.g. `log2_df`'s own
    precision bug (idea #58). Not implemented this round: this crate's
    own hard auto-vectorization requirement rules out a data-dependent
    "detect overflow, fall back to a safer path" branch (any fallback
    must be paid unconditionally on every call, in every lane, or it
    doesn't vectorize), so a real fix means redesigning the whole
    squaring loop around a uniformly-applied rescale, real engineering
    work beyond this iteration's scope -- left open as a well-defined,
    scoped follow-up, not a vague "someone should look at this."
    Currently `pown` has no `_checked`/`_wide` sibling tier at all (unlike
    `sin`/`exp2`/`powf`/`remainder`, which all offer a slower-but-safer
    variant for exactly this class of extreme input) -- that gap is
    itself worth closing, whether by fixing `pown` directly (if the
    rescale can be made cheap enough to justify changing the default) or
    by adding a `pown_wide` opt-in tier (matching this crate's own
    established default/checked split precedent). No code changed this
    round; this is a real, verified, quantified, root-caused bug report,
    not a fix. *Df32/double-float compensation fixes rounding-error bugs
    (too little precision) but is structurally the wrong tool for
    overflow bugs (too little range) -- know which class of problem
    you're looking at before reaching for double-float as the default
    hammer.*

    **Follow-up same day: tried the "obviously correct" complementary fix
    (exponent-tracking rescale, range extension instead of precision
    extension) and it *also* fails, for a genuinely instructive reason.**
    Prototyped decomposing `base`/`result` into `(mantissa in [1,2),
    exponent: i64)` pairs (frexp/ldexp-style), renormalizing every
    squaring/multiply step so the mantissa never leaves a safe range and
    the exponent -- plain integer arithmetic -- absorbs all the
    magnitude, only converting back to a single f32 (with a final
    saturating check) at the very end. This is exactly `remainder_wide`'s
    own rescue-scale idiom, generalized. Result: the specific traced case
    (`pown(0.997296, -32767)`) *still* returns `inf`, and the broader
    sweep got *worse* in one respect -- the existing documented `|n|<=64`
    max ulp regressed from `90` to `4,194,272`. Traced why: renormalizing
    the mantissa to stay in `[1,2)` doesn't reduce the *number of
    rounding events* in the mantissa chain at all -- it's the exact same
    32 iterations of f32 multiply/square, just re-expressed. The traced
    case's own numbers make this concrete: the rescue computation lands
    on `mantissa=1.0023601, exponent=128`, i.e. `~1.0024 * 2^128 ~=
    3.4108e38` -- almost exactly the *same* ~0.3%-too-high value the
    original algorithm computed (`3.4109e38`), just relabeled with an
    exponent one higher than the true answer's real exponent (`127`, not
    `128`) *because* that same ~0.3% of accumulated rounding error pushed
    the renormalized mantissa just over `2.0`, incrementing the tracked
    exponent by exactly the amount needed to reproduce the original
    overflow at a different threshold (`e>127` in the new scheme instead
    of `f32::MAX` in the old one). Same underlying number, same error,
    same failure, new representation. **This means the bug has two
    distinct, complementary failure ingredients -- insufficient range (
    which exponent-tracking alone *would* fix, for an `x` extreme enough
    that intermediate magnitudes need more exponent bits than f32 has,
    regardless of rounding) and insufficient precision (which `Df32`
    alone *would* fix, for compounding rounding error that stays within
    range) -- and this specific found case is purely the second kind, so
    neither fix alone touches it.** A real fix needs *both* combined: a
    mantissa tracked with `Df32`-level precision *and* an exponent
    tracked as a wide integer, in the same structure simultaneously --
    genuinely more engineering than either attempt above, closer to
    building a small custom extended-range double-float type than reusing
    an existing crate primitive. Confirms this is correctly scoped as
    real follow-up work, not something to force through this session --
    and the `|n|<=64` regression is a useful warning of its own: a
    plausible-looking rescale scheme can silently make an *already-good*
    range measurably worse if the renormalization itself isn't verified
    against the existing documented accuracy before considering it an
    improvement. *A bug hunt sometimes needs to identify that it's
    "solve problem A AND problem B simultaneously," not "A or B" -- two
    fixes that each correctly solve half the problem can each
    individually look like they "don't work at all" if judged only
    against a symptom that happens to require both halves.*

    **Final chapter, same day: implemented the combined fix for real
    (`pown_wide`, `Df32`-precision mantissa via a new `WideFloat` type +
    `Mul for Df32`), fully verified correct, then discovered a second,
    unrelated, genuine blocker -- LLVM will not auto-vectorize this
    computation, no matter how it's restructured. Not shipped; fully
    reverted; documented here in detail since the investigation itself
    is the reusable asset.**

    The real implementation (mantissa tracked as `Df32`, renormalized to
    `[1,2)`; exponent tracked as `f32`, not `i64` -- an `i64` field
    doesn't share this crate's 16-lanes-of-f32 vectorization width, found
    and fixed early via `codegen_check`) passed every verification this
    session's other fixes were held to: the traced bug case
    (`pown(0.997296, -32767)`) now returns the true value bit-exactly
    (`0` ulp); zero of the found overflow cases remained broken across
    the full structured sweep; the existing documented `|n|<=8`/`|n|<=64`
    ranges *improved* (11->8, 90->64 max ulp) rather than regressed;
    `cargo test` passed; a real intermediate-exponent-construction
    subtlety (`f32::powi(-128)` wrongly returns `0.0` -- std's own
    `powi` hits the *exact same bug class* this tier exists to fix,
    computing the positive exponent first and overflowing before
    reciprocating) was found and worked around with this crate's own
    `ROUND_MAGIC` bit-construction idiom (borrowed from `exp2_checked`,
    shift amount `23` confirmed by exhaustive check over `e` in
    `[-126,127]`, not derived by inspection alone). First `mca` numbers
    looked genuinely good: latency 132.28 cyc (*better* than plain
    `pown`'s 176.00), throughput 5.536 cyc/elem (+45% over `pown`'s
    3.805) -- a real cost for a real fix, in line with this crate's own
    `remainder`/`remainder_wide` precedent.

    Then `codegen_check` failed: `pown_wide_throughput` had **zero**
    packed arithmetic and a genuine scalar loop (`jmp`/`cmpq $16`/`je`
    over an explicit element index, not a masked blend) -- this crate's
    hard auto-vectorization requirement, violated. Root-caused via five
    separate isolation probes rather than guessing:
    - Removing all special-case handling (0/inf/nan, sign, `n==0`) and
      keeping just the bare 32-iteration loop: still fully scalar.
    - Removing the wide-exponent tracking entirely (pure `Df32`
      arithmetic only, no `WideFloat`): still fully scalar, identical
      instruction profile.
    - Replacing `Df32` (a struct) with raw `(f32,f32)` tuples and inlined
      `two_sum`/`two_prod` (no struct type, no method calls at all):
      still fully scalar, ruling out "struct-of-2-f32 confuses the
      vectorizer" as the cause.
    - Confirmed `pown` itself *also* compiles its `for i in 0..32`
      loop to a real runtime loop with real branches (58 `jne`/`testl`
      instructions) -- branches-around-a-vectorized-body is fine and
      already how this crate's own shipped `pown` works (LLVM hoists the
      scalar, per-call-shared `bit_set` check outside the per-lane
      arithmetic, branching around a *vectorized* multiply rather than
      blending it per lane); so "has branches" was never the
      discriminator between working and non-working.
    - Iteration count, bisected precisely: 8 iterations compiles to a
      fully vectorized loop (330 packed arithmetic instructions found);
      9 and 10 also vectorize (374, 418 packed instructions); **11
      iterations of the exact same per-iteration computation collapses
      to zero packed instructions** -- a sharp cliff, not a gradual
      falloff, confirmed by testing every iteration count from 8 to 32.

    **This means the fix's own correctness requirement (need up to 31-32
    iterations to cover `i32::MIN`/`MAX`) and this crate's own hard
    vectorization requirement are in direct, apparently irreconcilable
    tension for this specific per-iteration computation's complexity** --
    10 iterations (the vectorizing ceiling found) only covers
    `|n|<=1023`, nowhere near enough (the originally-traced bug needs 15
    iterations just for `n=32767`; the worst structured-search cases
    needed up to 30). This reads like a genuine LLVM loop-vectorizer
    cost-model threshold (total unrolled instruction count for a single
    tightly-dependent recurrence, not iteration count per se, given
    8/9/10 all vectorize with proportionally more instructions each) --
    not a workaround-able quirk of this crate's own code style, since
    three independently-restructured implementations (struct-based,
    tuple-based, exponent-free) all hit the *identical* wall at the same
    iteration count with the same per-iteration op cost.

    Reverted completely: `WideFloat`, `pown_wide`, the new `Mul for Df32`
    impl, and every harness wiring change (`mca_target.rs`, `mca.rs`'s
    `order` array, `accuracy.rs`'s sweep entries) -- confirmed via
    `git diff --numstat` that every changed file was a pure addition (0
    deletions), so `git checkout` cleanly restored the pre-existing
    state with no manual reconstruction risk. Left open as real,
    precisely-scoped future work, in decreasing order of promise:
    1. Split the computation into two (or more) independently-vectorized
       sub-loops of <=10 iterations each, combined *outside* whatever
       triggers the vectorizer's cost-model bailout (untested whether
       the combining step itself reintroduces the same wall).
    2. A lower-complexity per-iteration recurrence that still fixes the
       original bug's root cause (dominant error from the single f32
       reciprocal, amplified by `|n|`) without needing the full two_sum/
       two_prod machinery every iteration -- e.g. correcting only the
       initial reciprocal to `Df32` precision and leaving the squaring
       chain as plain `f32` might already shrink the amplified error
       enough in practice, worth checking against the same structured
       sweep before assuming the full per-iteration `Df32` treatment is
       necessary.
    3. Accept a documented, `pown_small`-style restricted-range safer
       tier (e.g. `|n|<=1023`, matching the confirmed 10-iteration
       vectorizing ceiling) instead of a full-range one -- narrower than
       hoped, but still a real, shippable, vectorizing improvement over
       nothing, if a future session decides that's worth doing.
    *A fix being fully verified correct doesn't mean it's finished --
    always confirm the target still meets every one of the crate's own
    hard requirements (here, auto-vectorization) before considering
    something ready to ship, and when a requirement conflicts with a
    fix's own inherent complexity, isolate exactly where the conflict
    lives (iteration count vs. total instruction count vs. specific
    operations) before concluding it's unfixable — a precisely
    characterized dead end is still a valuable, reusable result, even
    unshipped.*
