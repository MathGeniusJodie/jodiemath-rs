# faster? f32 math functions
Attempting to provide faster implementations of common f32 math functions with a similar level of accuracy to the standard library.
There are also perfectly rounded variants and faster-but-sloppier variants.

cbrt, cbrt_accurate, log2 and exp2_checked are full-range correct: negatives,
denormals, zero, inf and nan are all handled. exp2's default is the fast
unchecked version, only valid for x in [-126, 128); exp2_checked handles the
full range (overflow to inf, denormal underflow, nan) for ~2.5 ns extra
latency. cbrt_accurate is within 1 ulp on every input tested (exhaustively,
see below) -- **not** perfectly rounded on every input, despite an earlier
claim here: the full 2^32-pattern sweep found a systematic 1-ulp miss at one
specific mantissa, recurring at every octave from ~2^-126 to ~2^127.
sin/cos's default, like exp2's, is the fast unchecked version: a single-f32
Cody-Waite reduction (`fma`-based magic-constant rounding), accurate only
while q = round(x/pi) is exact, i.e. |x| < 2^22 * pi (~1.3e7) -- ~2.7x lower
serial latency and ~7x higher vectorized throughput than std (see
benchmarks below). **sin_checked/cos_checked trade that speed for gradual
degradation, not a cliff -- sub-ulp (avg) out to ~1e13, decaying smoothly
from there, and always finite for every finite input** -- for ~18 ns extra
serial latency (~1.3 ns extra per element in a vectorized loop).
q = round(x/pi) must be an exact integer for the reduction to land in
[-pi/2, pi/2]; a single f32 q (sin/cos's default, above) is a binary
either-or (exactly right, or off by a whole integer once |x| crosses q's
exact-integer ceiling), which shifts the residual by a whole multiple of pi
and puts the degree-9 polynomial (fit only for [-pi/2, pi/2]) hopelessly
outside its domain -- a relocatable *cliff*, not a slope, no matter how q
is rounded (tried both the classic magic-constant trick and a native
hardware `.round()`; the latter only moves the cliff from ~1.3e7 to ~2.6e7,
same shape). sin_checked/cos_checked fix this by giving q a second, small
f32 word (qh, ql) -- genuine double-float precision via
`two_prod`/`two_sum` error-free transforms (exact at any magnitude, unlike
the default's PI_A..D trick, which needs bounded q).
A cheap clamp on the *reduced residual* (before it reaches the polynomial)
still guarantees sin_checked/cos_checked's output can never overflow to
`inf`, however bad the reduction gets. Verified
exhaustively over all 2^32 bit patterns (examples/accuracy.rs,
examples/edgecheck.rs): avg ulp stays flat (~0.25-0.35) from small x out
through ~1e13, then climbs gradually (not a jump) through ~1e19 before
plateauing at essentially-uncorrelated-but-finite output for anything
larger. **Real bug found and fixed along the way**: folding cos_checked's
-0.5 phase-shift offset directly into the two_prod's dominant term breaks
once that term's own ulp exceeds 1 (|x| ~ 1.68e7) -- adding a fixed 0.5 to
an already-coarse-ulp float just rounds it away, silently corrupting
cos_checked's residual by a whole pi. This is the same class of bug as
everything else in this reduction (a small quantity lost against a big-ulp
value) and likely also affected an earlier, since-removed version of this
same double-float design that was never checked with a cos-specific
magnitude-bucketed sweep (only sin's was; the two aren't symmetric here).
Fixed by folding the offset into the small correction word instead, where
it survives regardless of the dominant term's own precision. A middle
ground that keeps *only* the two dominant tiers (dropping the third,
smallest one) was tried and measured to cost the same as keeping all three
(~130-140 cyc either way) -- there isn't a cheaper partial version once you
need genuine multi-word q precision at all, so the version here just keeps
all three tiers for the best accuracy at no extra cost.
sin_checked/cos_checked (then still named sin/cos) were the crate's default
for a while during this work; the cheap single-word version was reinstated
as the plain sin/cos afterwards, once it was clear the accurate version's
latency cost (below) wasn't worth paying unconditionally for every caller,
matching the exp2/exp2_checked split above.

Straight-ported from jodiemath's C library (github.com/MathGeniusJodie/jodiemath):
ln, log10, log1p, exp, expm1, sinh, cosh, tanh, asinh, acosh, atanh, asin, acos,
atan, atan2, tan, erf, erfc, hypot, powf, remainder. These reuse this crate's
own log_2/exp2/sin/cos internally rather than re-deriving them from scratch,
so they inherit those functions' own tradeoffs: exp, expm1, sinh, cosh and
tanh route through the fast *unchecked* exp2, so they share its
`[-126, 128)` domain limit (see exp's doc comment); tan shares sin/cos's
`|x| < 2^22*pi` domain. powf and erf/erfc's tail branch used to share that
same unchecked-exp2 limit but are now fixed (routed through exp2_checked
instead, see below). A few originally inherited real accuracy defects
from the C original's naive formulas -- confirmed by hand-tracing the
floating-point ops and by the exhaustive/fuzz sweep below, not assumed,
and identical in the C source, so these weren't translation bugs. Several
have since been fixed (each still a straight port of the *shape* of the
original algorithm, just with the specific cancellation/overflow bug
patched -- see each function's doc comment for the exact mechanism and
IDEAS.md for the before/after measurements):
- **fixed**: log1p, tanh and sinh/sinh_throughput all lost essentially all
  precision for small `|x|` (e.g. the old `sinh` returned exactly `0` well
  away from zero) from a `ln(1+tiny)` or "subtract two near-1 values"
  cancellation pattern; each now has a small-x branch (log1p: a Sterbenz
  correction term; tanh: reuses expm1; sinh: a direct Taylor branch).
- **fixed**: acosh returned `+inf` instead of `NaN` for large-magnitude
  negative x (squaring erased the sign before the domain check ever saw
  it, wrong for roughly the whole range x < -4096) and separately
  returned `+inf` instead of a small finite value for valid x above
  sqrt(f32::MAX) (`x*x` overflowing prematurely); both fixed with an
  explicit domain check and a rescaled/log1p-routed formula that avoids
  ever squaring or summing something that overflows.
- **fixed**: atanh lost essentially all precision for small `|x|` from
  forming `(1+x)/(1-x)` directly (rounds to exactly 1.0 for tiny x); now
  `0.5*(log1p(x) - log1p(-x))`, reusing log1p's fix directly (same shape
  as tanh reusing expm1).
- **fixed**: asinh had two bugs, not just the documented small-x cliff --
  the exhaustive sweep's actual worst case was large-*negative*-x
  cancellation (`x + sqrt(x^2+1)` nearly cancels for very negative x,
  losing relative precision), worse than the small-x one. Sidestepped
  entirely by computing on `|x|` (asinh is odd) with a rationalized
  `sqrt(x^2+1)-1` term routed through log1p, then restoring the sign with
  `mulsign`; also gained the same overflow guards as acosh for huge `|x|`.
- **fixed** (mostly): asin had the same near-zero cancellation as the
  others (rationalized the same way), plus two further, separate issues
  found one at a time as each fix's own re-sweep exposed the next: the
  rational correction's own fit has a small persistent relative bias that
  cancellation had been masking (fixed with a small-x Taylor branch below
  `|x| < 0.1`), and that same bias gets amplified by the sqrt singularity
  in asin's derivative right at x -> 1 (fixed by reusing acos's own
  well-conditioned formula, `asin(x) = pi/2 - acos(x)`, above `|x| > 0.9`).
  avg ulp comfortably under budget; max ulp down to 121 (from 852 million)
  but not fully closed -- the same relative bias then showed up at the
  small/mid branch boundary (x ~ 0.1) instead. Refit the mid branch's 4
  coefficients with examples/tune.rs's coordinate-descent tuner
  (zero perf cost, same instructions): max ulp 121 -> 84, avg ulp 0.325 ->
  0.377 (still comfortably under budget) -- not a complete fix, see asin's
  doc comment for why and what one would need.
- **fixed**: erfc's own `|x|<=10` clamp didn't fully protect its internal
  exp2 call: for `|x| >= ~9.35` the exponent it computes fell outside the
  *unchecked* exp2's domain, returning `inf` or huge finite garbage
  instead of a clean near-0 (or near-2 for negative x) -- not a subtle
  ulp issue, e.g. `erfc(9.5)` used to be `inf`. Fixed by routing through
  `exp2_checked` instead of `exp` for the Gaussian factor; its wider
  `[-151, 128)` domain comfortably covers the whole clamped range.
- **fixed**: erf's tail branch evaluated `erf_poly` (a plain unbounded
  degree-6 polynomial) directly on `|x|` with no bound at all -- worse
  than erfc's gap, since it wasn't just an exp2-domain issue: the
  polynomial's own positive leading coefficient makes it turn around and
  grow to `+inf` for large `|x|` instead of staying deeply negative
  (`erf_poly(9) ~ -92`, `erf_poly(20) ~ +8698`), so `erf(50)` was
  `-1.02e17`, `erf(100)` and `erf(+-inf)` were `NaN`, all of which should
  saturate to `+-1`. Fixed by bounding `|x|` to 10 before `erf_poly` (same
  bound erfc's clamp already uses) -- also swapped to `exp2_checked` as
  cheap extra insurance now that the input is bounded, which as a bonus
  turned out to be more accurate than plain exp2 even within the
  previously-tested range (avg ulp 0.631 -> 0.319, moving erf from
  *over* budget to comfortably under it).
- **fixed**: powf composed the *checked* log_2 with the *unchecked* exp2 --
  the expensive half bought correctness the cheap half then threw away.
  Unlike erf/erfc, the doc comment here was already honest about the
  limitation (not a false safety claim), but the practical effect was the
  same: `powf(2.0, 500.0)` (should be `inf`) was `2.88e17`,
  `powf(2.0, 1000.0)` was `3.6e-12` -- plausible finite garbage for any
  `y` large enough to push `log2(x)*y` past 128. Fixed by routing through
  `exp2_checked`; real, disclosed perf cost since this touches every
  call, not just an edge branch (+14% latency, +36% throughput).
- **fixed**: `acos(-0.0)` returned `-pi/2` instead of the correct `+pi/2`
  (acos is never negative -- unlike odd functions like sin/asinh, `-0.0`
  has no legitimate negative result here). Root cause: `mulsign` (bit-based
  sign) and `x < 0.0` (value-based comparison) disagree on exactly this
  one input, whose sign *bit* is set but whose *value* equals `+0.0`.
  An exhaustive-sweep-only find -- fuzz sampling essentially never lands
  on this one bit pattern (avg ulp exhaustively was 0.99 with max ulp in
  the billions, invisible in quick-mode's 0.50/4 reading). Fixed by
  normalizing `-0.0` to `+0.0` before `mulsign` sees it (`x + 0.0`, exact
  and a no-op for every other input). Essentially free (mca unchanged,
  37.11 cyc latency both before and after).
- **fixed**: `atan2(-0.0, +0.0)` returned `+0.0` instead of the
  IEEE754/C99-defined `-0.0` -- found by systematically checking every
  other `mulsign` call site for the same bug shape after fixing acos.
  Different mechanism than acos's: when `x` is exactly `+0.0`, the
  formula's `base` term degenerates to exactly `+0.0`, and combining it
  with the correctly-signed `-0.0` correction via `base + mulsign(...)`
  hits IEEE754's rule that adding two *opposite*-signed zeros always
  gives `+0.0`, silently destroying the correction's sign. Fixed by
  skipping that addition when `base` would be the degenerate `+0.0`.
  Verified all 12 zero/sign combinations bit-exact against std; zero
  perf cost (mca bit-for-bit unchanged).
- **fixed**: `sin(-0.0)`, `cos(-0.0)` (silently, since `cos(-0.0)=+1.0` was
  already the correct nonzero answer), `tan(-0.0)`, and `sin_checked(-0.0)`
  all returned `+0.0` instead of the IEEE754/C99-defined `-0.0` --
  found by the same systematic sweep as atan2's fix, checking every odd
  function at `-0.0`. Two distinct root causes, same underlying IEEE754
  rule (adding two opposite-signed exact zeros always gives `+0.0`):
  `sinf_poly` (shared by all of sin/cos/sin_checked/cos_checked/tan)
  computes `fma(p, x3, x)`, and at `x = +-0.0`, `x3` correctly carries
  `x`'s sign but `p` (the poly's fixed leading coefficient at `y=0`, sin's
  own curvature) is a negative constant, so `p*x3` always ends up the
  *opposite* sign to `x` at this one point, and the `fma` silently loses
  it. Fixed with `r.copysign(x)` -- free for every nonzero `x` (sin is odd
  and monotonic on this poly's domain, so the leading `x` term always
  dominates `p*x3` in magnitude there, meaning `r`'s sign already equals
  `x`'s), only changes the singular zero case; confirmed cheaper than an
  `x == 0.0` branch/select tried first (that cost real throughput --
  sin/cos +12-17%, sin_checked/cos_checked/tan +1-9% -- since it's inlined
  into every caller; copysign costs far less). `sin_checked` needed a
  *second*, separate fix: `reduce_pi`'s own multi-term two_sum/two_prod
  error-compensation chain independently loses `x`'s sign somewhere
  internally (same IEEE754 rule, exact spot not traced), well before
  `sinf_poly` is even reached, so `sinf_poly`'s own fix can't see the
  original sign to restore. Guarded at `sin_checked`'s own output with an
  explicit `x == 0.0` select instead; `cos_checked` needs no such guard
  (`cos_checked(-0.0) = +1.0`, a nonzero result, unaffected). Verified
  bit-exact against std for all 5 functions at `x = -0.0`.
- **fixed**: `log1p(-0.0)`, `atanh(-0.0)`, and `remainder(-0.0, y)` all lost
  their sign the same way, found by the same broader `-0.0` sweep.
  `log1p`'s `ln(u) + corr` adds two exactly-zero values of opposite sign
  at `x = +-0.0` (`ln(1.0)` is `+0.0`, but `corr` correctly carries x's
  sign there). Fixed with a trailing `if x == 0.0 { x } else { normal }`
  select (computes the normal path unconditionally first, `log_2`'s own
  established "select, not early return" idiom) -- log1p is odd and
  monotonic through the origin, so `normal`'s sign already matches x's
  for every nonzero x, making this a no-op everywhere except the singular
  zero point. `atanh(x) = 0.5*(log1p(x) - log1p(-x))` reuses log1p
  directly, so its own `-0.0` bug fell out fixed for free, no separate
  change needed. `remainder`'s `fma(-q, y, x)` hits the identical
  IEEE754 mechanism, but unlike log1p, remainder's sign does *not*
  generally track x's sign for nonzero x (`remainder(2.0, 3.0) == -1.0`
  is correct, not a bug), so a blanket copysign fix would be wrong here
  -- fixed with the same trailing-select shape instead. A plain early
  `return x;` guard was tried first for both and gave better latency, but
  broke `llvm-mca`'s region markers when inlined into a vectorized loop
  (`cargo run --example mca` failed outright with a "found an invalid
  region end directive" error) -- reverted in favor of the select form,
  which compiles cleanly everywhere. Real, disclosed side effect: since
  `asinh`/`acosh` already computed their own correct sign externally via
  `mulsign` (this fix provides them zero actual benefit, log1p is always
  called with a provably-`+0.0` argument in their flow), the extra
  in-lined select still costs their *latency* substantially (asinh
  55.40->120.99 cyc, +118%; acosh 110.47->127.75 cyc, +16%) while their
  *throughput* -- the metric this crate's vectorization-first design
  actually prioritizes -- improved instead (asinh 9.625->7.716,
  acosh 6.839->6.588 cyc/elem). `log1p`/`atanh`/`remainder` themselves are
  each within noise of their own pre-fix baseline. Verified bit-exact
  against std for all three functions at their `-0.0` inputs; exhaustive
  accuracy.rs sweeps for log1p/atanh/asinh/acosh all matched the pre-fix
  documented baseline exactly.
- **still open**: remainder's `x - round(x/y)*y` loses precision to
  cancellation once `|x/y|` is large, since `round(x/y)*y`'s absolute
  error scales with `ulp(x)`, which can exceed the true remainder's own
  magnitude (at most `|y|/2`).
- **fixed**: `powf(x, y)` for negative `x` was *always* `NaN`, even for
  well-defined, common cases like `powf(-2.0, 3.0)` (should be `-8.0`).
  Found by the same `-0.0` sweep, widened to a couple of general negative
  bases once `powf(-0.0, 3.0)` turned up wrong. Root cause: the whole
  function is `exp2_checked(log2(|x|)*y)`-shaped, and `exp2` of any real
  argument is always non-negative -- this route has no way to ever
  produce a negative result, regardless of `y`. Fixed by computing the
  magnitude on `|x|` as before, then reapplying the correct sign for
  negative `x` when `y` is an integer (even `y` -> positive, odd `y` ->
  negative, reusing `parity`, the same integer-parity helper
  `sin_checked`/`cos_checked` already use) and `NaN` when `y` isn't an
  integer (matches std -- e.g. `(-8.0)^(1/3)` is `NaN` in f32 too, real
  cube roots of negative numbers aren't reachable through this branch).
  Two more bugs surfaced from the same investigation and were fixed in
  the same pass: `pow(x, 0)` must be `1` for *any* `x` (even `0`,
  negative, or `NaN`) per IEEE754/C99's dedicated special case, not
  derivable from the log/exp2 formula (`0*inf`/`NaN*0` both degrade to
  `NaN`) -- `powf(0.0, 0.0)` and `powf(f32::NAN, 0.0)` were both `NaN`
  instead of `1.0`. And the negative-base sign check itself first used
  `x < 0.0` (value-based), which -- same pitfall as acos's earlier fix --
  disagrees with the bit-based sign exactly at `x = -0.0`
  (`-0.0 < 0.0` is `false`), silently routing `powf(-0.0, 3.0)` through
  the wrong branch; fixed with `x.is_sign_negative()` instead. Not fixed:
  `powf(-1.0, ±inf)` (IEEE754 special-cases this to `1.0`; this crate's
  formula gives `NaN` there, `0*inf` inside `log2(1.0)*inf`) -- a narrow,
  rarely-relied-on corner left as a known gap rather than adding more
  special-case logic for it. accuracy.rs's own `powf` domain filter had
  also been narrowed to `x > 0.0` only, dodging the whole negative-base
  case instead of exercising it (the same pattern erf/erfc's filters
  had) -- widened to `x != 0.0`. Verified bit-exact against std at every
  case above; a dedicated integer-`y` spot sweep (fuzz sampling almost
  never lands on an exact integer `y`) found no accuracy cliff in the
  newly-reachable negative-base path (max ulp 92, same budget as the
  existing positive-base path, since it reuses the identical
  log_2/exp2_checked machinery with only a sign correction at the end).
  mca cost is small: 97.88->98.03 cyc latency (+0.2%), 3.788->3.898
  cyc/elem throughput (+2.9%).
- **fixed**: three infinity-handling gaps in the two-arg functions, found
  by extending the same sweep to `+-inf`/`NaN` combinations after the
  powf fix above. `remainder(finite x, +-inf)` was `NaN` instead of `x`
  (IEEE754/C99: `q` rounds to exactly `0.0` for any finite `x`, but
  multiplying that zero by an *infinite* `y` gave `NaN` via `0*inf`,
  instead of the intended "no reduction happened" no-op) -- fixed with a
  trailing `if y.is_infinite() && x.is_finite() { x } else { r }` select
  (`x`/`y` themselves infinite/NaN still correctly fall through to `NaN`,
  matching std). `hypot(+-inf, NaN)` was `NaN` instead of `+inf` --
  IEEE754/C99 special-cases infinity to "win" over NaN here (unlike
  almost every other function), which the naive `x*x+y*y` formula can't
  reach on its own once either argument actually is NaN -- fixed with a
  trailing `if x.is_infinite() || y.is_infinite() { INFINITY } else { normal }`
  select; distinct from hypot's already-documented finite-overflow
  tradeoff, unaffected either way. `atan2(+-inf, +-inf)` was `NaN`
  instead of a defined `+-pi/4`/`+-3pi/4` by quadrant -- `y/x` is
  `inf/inf` (`NaN`) there, so the general atan-based formula has no ratio
  to work with at true infinity; IEEE754/C99 define a fixed convention
  instead, added as a trailing select the same way. All three follow the
  same "compute everything unconditionally, select last, no early
  returns" idiom as the log1p/remainder/powf fixes above (an early return
  here would risk the same `llvm-mca` region-marker corruption already
  found once this session). mca cost is negligible: atan2 unchanged
  (57.17/1.467), hypot +0.5%/+0.4% (21.00->21.11 cyc, 0.763->0.766
  cyc/elem), remainder unchanged (33.02/0.647).
- **improved**: `asin`'s worst-case accuracy, in three rounds. First, the
  small/mid Taylor threshold widened `0.1 -> 0.3` after re-checking (and
  disproving) a prior doc comment's own untested claim that the residual
  "can't be shrunk by moving thresholds" -- max ulp 84 -> 41, zero perf
  cost. Then, investigating *why* the acos-derived `near1` branch stayed
  accurate well below its `a > 0.9` cutoff led to removing the `mid`
  branch (a whole separate rational-correction formula) entirely: it
  turns out `acos_poly` is `acos`'s own full-`[0,1]`-domain poly, not
  something limited to "near 1" -- the name was just where it happened to
  get reused first. A 2-branch `asin_small`/`acos_poly` structure beats
  the old 3-branch one everywhere past x~0.25: max ulp 41 -> 11, avg ulp
  0.105 -> 0.033. mca is a genuine mixed result: throughput improved a
  lot (2.899 -> 0.968 cyc/elem, -67%, less total vectorized work) but
  latency got *worse* (43.24 -> 59.03 cyc, +37%) -- removing the `mid`
  branch also removed independent work the scalar chain could previously
  use to fill cycles otherwise spent waiting on the sqrt. Kept: both
  accuracy and throughput (this crate's prioritized metric) improved by a
  wide margin, only the secondary latency number regressed. Third,
  `acos_poly`'s coefficients refit against a *joint* objective (`acos`'s
  own error plus asin's use of the same poly, since the two callers'
  outputs scale oppositely near x=1, giving the same absolute poly error
  a different ulp weight for each) with a constraint that acos's own
  on-grid max ulp may never regress -- max ulp 11 -> 9, avg 0.033 ->
  0.030 for asin, `acos` itself exactly unchanged (max ulp 4, avg 0.496,
  bit-for-bit identical). Zero perf cost throughout (same instructions,
  only literal constants changed in every round) -- see asin's own doc
  comment (fixes 5-7) for the full investigation.

**sinh_throughput/cosh_throughput** are a second tier for sinh/cosh, added
after finding that computing `exp(-x)` as `1.0 / exp(x)` (instead of a
second independent `exp` evaluation) is a genuine tradeoff, not a strict
win: on this CPU the FP divider is close to idle even when the FMA/mul
ports are saturated, so a vectorized loop sees a large throughput gain
(one whole `exp` poly evaluation replaced by one division), but a single
serial call is *slower*, since the division can't start until `exp(x)` is
done -- unlike the default's two independent `exp` calls, which a
CPU with enough ports can run in parallel. Plain `sinh`/`cosh` keep the
lower-latency, parallel-exp form as the default (matching every other
tier split in this crate, where the plain name is the "normal" choice);
`sinh_throughput`/`cosh_throughput` are opt-in for callers who mainly loop
over arrays. Accuracy is essentially unaffected (one extra rounding from
the division; see the precision table below).

All functions auto-vectorize, it's a hard requirement

# precision (see examples/accuracy.rs)
Fuzz-mode (100M random f32 bit patterns/function); pass `thorough` for an
exhaustive sweep of all 2^32 patterns instead (few minutes, needs --release).
```
                        | jodie avg  | jodie max | std avg | std max
------------------------|------------|-----------|---------|--------
                   cbrt |    0.326   |     3     |    0    |    0
          cbrt_accurate |    0.000   |     1     |    0    |    0
                   exp2 |    0.030   |     1     |  0.000  |    1
           exp2_checked |    0.016   |     1     |  0.000  |    1
                   log2 |    0.003   |     3     |  0.000  |    1
        sin (|x|<1.3e7) |    0.065   |   1183    |  0.003  |    1
        cos (|x|<1.3e7) |    0.293   |   2780    |  0.002  |    1
 sin_checked (|x|<=1e6) |    0.036   |     2     |  0.000  |    1
 cos_checked (|x|<=1e6) |    0.081   |     3     |  0.000  |    1
  sin_checked (all f32) | (degrades gradually past |x| ~ 1e13 -- see note above; 0.007/1 for std)
  cos_checked (all f32) | (degrades gradually past |x| ~ 1e13 -- see note above; 0.007/1 for std)
```
The sin/cos/sin_checked/cos_checked rows are exhaustive (all 2^32 bit
patterns, examples/accuracy.rs `thorough` mode); the rest of the table is
the default 100M-sample fuzz mode. sin/cos's row covers its whole documented
domain (|x| < 2^22*pi, ~1.3e7 -- see the overview above); its max ulp
climbing into the thousands even in-domain, well before the cliff at the
domain edge, is expected -- the single-word Cody-Waite reduction's own
rounding error grows as |x| grows, independent of whether q is still
exactly rounded (std's error stays ~1 ulp across the same domain, for
comparison). Magnitude-bucketed
avg/max ulp for sin_checked/cos_checked (also now exhaustive, not just
fuzzed) shows the actual shape of the degradation -- flat and low for a very
long stretch, then a real but gradual climb, never a sudden jump to garbage
and never `inf` (confirmed by an exhaustive all-2^32-pattern check that no
finite input produces a non-finite output). These buckets are a permanent
part of `examples/accuracy.rs`, so the shape claim stays checkable after
future changes instead of just asserted:
```
        range | sin_checked avg | sin_checked max | cos_checked avg | cos_checked max
--------------|-----------------|------------------|-----------------|------------------
     |x|<=1e6 |           0.036 |                2 |           0.081 |                3
   [1e7,1e8)  |           0.252 |                6 |           0.262 |                6
  [1e9,1e10)  |           0.252 |               48 |           0.252 |               67
 [1e12,1e13)  |           0.292 |           154382 |           0.286 |           154382
 [1e15,1e16)  |           3.2e8 |             2.4e9|           9.7e8 |             2.4e9
```

Straight-ported functions (100M-sample fuzz mode / exhaustive where noted;
see the overview above for each function's real domain, and its doc comment
for known inherited defects and any fixes since). log1p, sinh/
sinh_throughput, tanh, acosh, atanh, asinh and asin all had the same class
of small-x cancellation defect but are fixed now (see above; asin's max
ulp has a separate, not-fully-closed residual near x=1, see its doc
comment); their rows are exhaustive (all 2^32 f32 bit patterns), not fuzz.
```
                        | jodie avg  | jodie max | std avg | std max
------------------------|------------|-----------|---------|--------
                     ln |    0.126   |     3     |  0.000  |    1
                  log10 |    0.286   |     4     |  0.000  |    0
                  log1p |    0.106   |     4     |  0.000  |    0
      exp (in-domain)   |    0.289   |    64     |  0.000  |    1
    expm1 (in-domain)   |    0.240   |    63     |  0.000  |    0
     sinh (in-domain)   |    0.291   |    64     |  0.000  |    0
     cosh (in-domain)   |    0.274   |    63     |  0.000  |    0
 sinh_throughput (in-domain) | 0.293 |    64     |  0.000  |    0
    cosh_throughput (in-domain) |    0.272   |    64     |  0.000  |    0
     tanh (in-domain)   |    0.144   |     6     |  0.000  |    0
                  asinh |    0.173   |     4     | (std also imperfect at extreme |x|)
                  acosh |    0.063   |     4     |  0.000  |    1
                  atanh |    0.032   |     3     |  0.037  | 363409
                   asin |    0.030   |     9     |  0.000  |    0
                   acos |    0.496   |     4     |  0.000  |    0
                   atan |    0.186   |    18     |  0.000  |    0
       tan (in-domain)  |    0.331   |  2967     |  0.000  |    0
                            erf  |    0.319   |     5     | (no std erf)
                   erfc (|x|<=10)|    0.311   |   109     | (no std erfc)
                  atan2 |    0.136   |    18     |  0.000  |    0
        hypot (bounded) |    0.034   |     1     |  0.000  |    0
        powf (in-domain)|    0.181   |   127     |  0.000  |    1
     remainder (|x/y|<1000) |    ~1000**  | ~2e9** | (no std remainder)
```
`**` remainder's max ulp stays large even inside the `|x/y|<1000` bound: a
handful of inputs land close enough to an exact half-integer quotient that
f32 rounding flips which integer `round(x/y)` picks vs. the f64 reference,
jumping the result by a whole `y` -- an inherent tie-breaking sensitivity of
any round()-based remainder, not specific to this formula. Fusing the final
`x - q*y` into one `fma(-q, y, x)` (single rounding instead of two) roughly
halves the *average* ulp across repeated fuzz runs (avg fluctuates run to
run, ~1000-2500 before vs ~600-1800 after over several trials, since the
10M-sample fuzz rarely hits the pathological tie-break inputs that dominate
max) but doesn't move the max-ulp tie-breaking cliff itself, which is a
separate, structural property of round()-based remainder.

# benchmarks
Run on i5-1145G7, -C target-cpu=native (now set in .cargo/config.toml)
```
Serial latency (dependency chain, examples/quickbench.rs; lower is better)
              | jodie   | std     | improvement
--------------|---------|---------|------------
         cbrt | 12.9 ns | 22.0 ns | 1.7x
cbrt_accurate | 17.0 ns | 22.0 ns | 1.3x
          cos | 12.7 ns | 12.8 ns | 1.0x
  cos_checked | 30.3 ns | 12.8 ns | 0.4x
         exp2 |  8.5 ns | 13.3 ns | 1.6x
 exp2_checked | 13.3 ns | 13.3 ns | 1.0x
         log2 | 13.1 ns | 14.9 ns | 1.1x
          sin | 10.9 ns | 12.9 ns | 1.2x
  sin_checked | 29.2 ns | 12.9 ns | 0.4x
           ln | 15.3 ns | 14.5 ns | 0.9x
        log10 | 14.3 ns | 16.6 ns | 1.2x
        log1p | 16.5 ns | 21.1 ns | 1.3x
          exp | 13.1 ns | 13.6 ns | 1.0x
        expm1 | 14.6 ns | 19.1 ns | 1.3x
         sinh | 16.9 ns | 20.1 ns | 1.2x
         cosh | 27.7 ns | 37.5 ns | 1.4x
         tanh | 35.9 ns | 33.8 ns | 0.9x
        asinh | 43.7 ns | 86.3 ns | 2.0x
        acosh | 35.1 ns | 27.0 ns | 0.8x
        atanh | 23.2 ns |  4.2 ns | 0.2x
         asin | 19.4 ns |  4.1 ns | 0.2x
         acos | 13.0 ns |  4.0 ns | 0.3x
         atan | 16.5 ns | 27.4 ns | 1.7x
        atan2 | 16.6 ns | 33.8 ns | 2.0x
          tan | 23.5 ns | 27.3 ns | 1.2x
          erf | 21.0 ns |     -   |  -
         erfc | 22.9 ns |     -   |  -
        hypot |  8.4 ns | 13.9 ns | 1.7x
         powf | 24.2 ns |  3.2 ns | 0.1x
    remainder | 12.7 ns |     -   |  -
```
atanh/asin/acos/powf's std comparisons here are suspiciously fast (4, 4, 4
and 3.2 ns -- close to or under this crate's own fastest functions) and are
likely partly LLVM constant-folding artifacts from quickbench's methodology
(sinh/asin/acos/atan2/powf are benched via a fixed second argument, e.g.
`|x| x.powf(2.0)`, which can let LLVM specialize the std call at compile
time in a way a real call site with a variable exponent wouldn't get) rather
than std's true per-call cost -- take these specific ratios with a grain of
salt pending a version of quickbench that varies both arguments.
sin_checked/cos_checked's latency is worse than std (2026-07-06, later same
day: reinstated a double-float reduction -- see the accuracy note above for
why the cheap single-word version moved behind these _checked names instead
of being std's default). Deliberate trade, same shape as the very first
version of this reduction: gradual, predictable degradation instead of a
cliff costs real latency. (Improved slightly from 32.9/31.7 ns by the
same-day critical-path shortening described below -- still 0.4x at this
rounding.) sin/cos (the crate's default, restored to the pre-double-float
single-word version afterwards) are back to ~1.0-1.2x std, matching the
exp2/exp2_checked fast-default/checked-full-range split.
```
Throughput (independent array evals over [f32; 4096], examples/quickbench.rs; lower is better)
              | jodie    | std     | improvement
--------------|----------|---------|------------
         cbrt | 0.37 ns  | 3.93 ns | 10.6x
cbrt_accurate | 0.67 ns  | 3.93 ns | 5.9x
          cos | 0.27 ns  | 3.41 ns | 12.6x
  cos_checked | 1.61 ns  | 3.41 ns | 2.1x
         exp2 | 0.23 ns  | 3.13 ns | 13.7x
 exp2_checked | 0.48 ns  | 3.13 ns | 6.6x
         log2 | 0.53 ns  | 4.10 ns | 7.8x
          sin | 0.22 ns  | 3.10 ns | 14.4x
  sin_checked | 1.49 ns  | 3.10 ns | 2.1x
           ln | 0.55 ns  | 3.82 ns | 7.0x
        log10 | 0.53 ns  | 5.30 ns | 10.1x
        log1p | 0.60 ns  | 6.28 ns | 10.4x
          exp | 0.32 ns  | 3.09 ns | 9.7x
        expm1 | 0.49 ns  | 6.46 ns | 13.2x
         sinh | 0.69 ns  | 14.84 ns | 21.5x
         cosh | 1.69 ns  | 22.54 ns | 13.3x
         tanh | 1.04 ns  | 18.06 ns | 17.3x
        asinh | 2.19 ns  | 41.42 ns | 18.9x
        acosh | 1.32 ns  |  6.17 ns | 4.7x
        atanh | 0.79 ns  |  4.54 ns | 5.7x
         asin | 0.46 ns  |  4.01 ns | 8.7x
         acos | 0.26 ns  |  4.01 ns | 15.6x
         atan | 0.62 ns  |  8.28 ns | 13.3x
        atan2 | 0.62 ns  | 13.30 ns | 21.4x
          tan | 0.76 ns  |  9.07 ns | 11.9x
          erf | 0.70 ns  |     -    |  -
         erfc | 0.68 ns  |     -    |  -
        hypot | 0.24 ns  |  2.92 ns | 12.0x
         powf | 0.95 ns  |  0.07 ns | 0.07x
    remainder | 0.22 ns  |     -    |  -
```
powf's std throughput row here (0.07 ns, faster than a single cycle) is
almost certainly the same fixed-second-argument constant-folding artifact
noted in the latency table above (LLVM likely reduces `x.powf(2.0)` to
`x*x` at compile time) rather than a genuine per-call cost -- the other
ratios are consistent with the rest of this crate's usual 5-20x vectorized
throughput advantage.
sin_checked/cos_checked's throughput is ~2.1x std (was ~9-11x with the
cheap single-word-plus-clamp version, ~13-14x with the original
single-word-but-inf-prone version -- i.e. today's sin/cos default,
restored unchanged, now measuring ~12.6-14.4x). This is the direct cost of genuine
gradual degradation: a double-float q needs several `two_prod`/`two_sum`
error-free transforms (each 2-6 ops) to stay accurate well past a single
f32's exact-integer range, instead of one cheap magic-constant rounding
trick. Tried to find a cheaper partial version (only the two biggest
cross-term tiers, dropping the smallest) -- measured to cost the *same* as
keeping all three (~130-140 cyc either way, see the mca table below), so
there's no meaningfully cheaper middle ground once multi-word q is needed
at all; kept all three tiers since the third one is free once you're
already paying for the other two. Two vectorization/perf pitfalls to
remember if this reduction is ever revisited (see jodiemath-workflow memory
for the full trail, including the cos pre_offset bug described above):
- a `q as i64` cast for the parity bit doesn't vectorize at all -- Rust's
  float-to-int cast is saturating, so LLVM falls back to a scalar
  convert-with-NaN/range-check per lane (~9x slower than std). Parity here
  uses a `floor`-based "mod 2" instead (`q - 2*(q*0.5).floor()`), the same
  instruction family as the `.round()` already used elsewhere, so it stays
  vectorized (confirmed via `--emit=asm`: 0 scalar convert instructions).
- combining several small correction terms via plain adds first (instead
  of feeding all of them through a sequential compensated-sum loop) cuts
  the summation loop from 7 iterations to 4 with no accuracy cost --
  verified bit-for-bit identical output over the full accuracy.rs sweep.
Two further micro-optimizations landed later the same day, found by
re-running `llvm-mca --bottleneck-analysis` on the *current* code (it now
skews more dependency-chain-bound than the ~39%-resource-pressure region an
earlier attempt at this measured, so shortening the critical path pays off
here where it didn't before -- always re-check, don't assume an old
bottleneck-analysis finding still applies after the surrounding code
changes):
- `reduce_pi`'s first correction merge (`two_sum(s, -e1)`, combining the
  residual right after subtracting the dominant term with that
  subtraction's own two_prod rounding-error) uses the cheaper 3-op
  `quick_two_sum` (Fast2Sum) instead of full 6-op `two_sum`. **Correction to
  an initial claim here**: a first check (only near exact multiples of pi,
  a narrow and unrepresentatively well-conditioned slice) seemed to show the
  `|s|>=|e1|` ordering Fast2Sum needs for *exactness* holds everywhere; a
  broader recheck (uniform random bit patterns, not just near-exact
  multiples of pi) found this is false -- real violations starting around
  `|x| ~ 1e3` and exceeding 80% of samples by `|x| ~ 1e8`+. Kept anyway,
  because Fast2Sum's failure mode when misordered is *bounded* (`e` off by
  up to ~1 ulp of `s`, not unbounded), and re-verified against the metric
  that actually matters -- the full accuracy.rs exhaustive sweep plus a
  magnitude-bucketed sweep against std out to f32::MAX -- shows no
  measurable difference from the full-`two_sum` version at any magnitude.
  Lesson: checking a theoretical invariant (an ordering assumption) is not
  the same as checking the thing that matters (final ulp error); it can be
  technically false while still being practically harmless. The other three
  terms in that loop (`p2`, `p3`, `tier2`) violate this ordering far more
  severely (`p3` in particular is usually exactly 0 but occasionally ~pi,
  comparable to the whole residual, whenever `ql != 0`) and were left on
  full `two_sum`, untested whether quick_two_sum would be harmless there too.
  Tried going further and pre-combining more of the 4 correction terms into
  a shallower tree (a single merge instead of 4 sequential ones, or even
  just pairing off `p3`+`tier2` in parallel with the dominant subtraction)
  -- both made latency *and* throughput worse (e.g. the full-tree version:
  sin 140→149 cyc latency, 7.609→8.545 cyc/elem), despite being exact and
  despite the region being dependency-bound: more total register-live-range
  pressure in the 16-wide vectorized loop apparently outweighs the shorter
  chain. Reverted both; kept only the one verified-safe single swap.
- the final parity combine (`parity(parity(qh) + parity(ql))`) doesn't need
  its 3rd `floor`-based `parity()` call: `parity(qh)` and `parity(ql)` are
  each exactly 0.0 or 1.0, so their sum mod 2 is just whether they differ,
  i.e. `if pq == pl { 0.0 } else { 1.0 }` -- a compare+select instead of a
  4-op floor chain, still branchless/vectorized (confirmed via `--emit=asm`:
  0 scalar convert instructions).
Both changes verified against the full accuracy.rs exhaustive sweep (`sin
|x|<=1e6` 0.036→0.036 avg / 2→2 max ulp, `cos` 0.077→0.077 avg / 2→2 max
ulp -- unchanged) and a magnitude-bucketed sweep confirming the same
gradual (non-cliff) shape past 1e6. Net: sin 140.00→136.00 cyc latency,
7.609→7.109 cyc/elem throughput; cos 144.00→140.00 cyc latency,
6.657→5.850 cyc/elem throughput.

**Further pass, same reduction, later session: revisits the "pairing off
p3+tier2" idea rejected just above, this time getting a real win from it --
the earlier attempt most likely had the same sign bug this session found
and fixed, though its code wasn't kept to confirm directly.** `round_x_over_pi`'s
`two_sum(e0, x*RPI_LO)` (comparable-magnitude operands, same
ordering-not-guaranteed-but-bounded situation as the swap above) moved to
`quick_two_sum` -- free, no measurable accuracy change. Then, in
`reduce_pi`: `p3` and `tier2` are the only two of the four correction terms
that depend on `ql` (round_x_over_pi's last-ready output, needing its whole
chain), while `p1`/`p2`/the dominant `x - p1` subtraction only need `qh`,
ready much earlier -- so combining `p3+tier2` into one value via `two_sum`
runs parallel to the `qh`-only work instead of stacking as two more
sequential merges after it, and the merge of that combined value into the
running residual was further downgraded to `quick_two_sum`. **Real bug hit
while building this**: `two_sum(p3, tier2)` guarantees `p3t + e3t == p3 +
tier2` exactly, so subtracting `(p3+tier2)` from the residual means
subtracting *both* `p3t` and `e3t` -- the first attempt added `e3t` into the
running low-order correction instead of subtracting it, which is a sign
error, not a rounding one. It passed small-`x` spot checks but broke `cos`
badly on the full sweep (`|x|<=1e6` avg 2.51 ulp, max 37M ulp, worst case
right at `x ≈ 2.5π` and `3.5π` -- `cos`'s zero crossings, exactly where `p3`
and `tier2` partially cancel, making `e3t` unexpectedly large instead of the
negligible correction it usually is). Fixed by flipping the sign; re-ran the
full exhaustive sweep before trusting it further. This is very likely the
same class of mistake behind the "made things worse" result recorded just
above for the identical-sounding idea -- a sign error wouldn't necessarily
look like a correctness bug in an mca run (mca has no notion of numerical
correctness, only codegen), so a broken version could plausibly still
compile to *legitimately worse* code for unrelated reasons and get reverted
without the sign bug itself ever being noticed. Not provable without the old
code, which wasn't kept, but worth recording as a caution: an idea "already
tried and found not to help" is only as trustworthy as the correctness of
the attempt, and llvm-mca's cycle counts can't tell you whether the code
being measured was actually right.
`llvm-mca --bottleneck-analysis` on `cos_throughput` (this session, before
these changes) showed a roughly even split -- 53.5% resource pressure
(`ICXPort0`/`ICXPort1`, the fma/mul ports) vs. 61.7% register/data
dependencies -- slightly skewed towards the dependency chain, consistent
with a chain-shortening restructure (rather than an op-count cut) paying
off here. Verified against the full exhaustive accuracy.rs sweep: `sin
|x|<=1e6` unchanged (0.036 avg / 2 max ulp); `cos |x|<=1e6` 0.077→0.081 avg
/ 2→3 max ulp -- a small, bounded cost from the two additional
`quick_two_sum` downgrades, still >10x under the 1-ulp-average budget. A
magnitude-bucketed sweep (see the precision table above, now a permanent
part of examples/accuracy.rs) confirms the same gradual, non-cliff shape
survives out past 1e15. Net: sin 136.00→124.00 cyc latency (-8.8%),
7.109→6.421 cyc/elem throughput (-9.7%); cos 140.00→128.00 cyc latency
(-8.6%), 5.850→5.354 cyc/elem throughput (-8.5%).
The throughput gap vs std comes almost entirely from vectorization: std's
functions have branches, so LLVM can't vectorize loops that call them.
Absolute ns swing session-to-session with CPU thermal state (the laptop
throttles up to ~2.5x mid-session) -- only trust jodie-vs-std ratios measured
in the same run. examples/mca.rs gives a thermal-noise-free second opinion in
cycles instead of ns:
```
theoretical cost from llvm-mca (-mcpu=native, 100 iterations)
                    | latency (cyc) | throughput (cyc)
--------------------|----------------|------------------
cbrt                |          35.06 |             1.629
cbrt_accurate       |          59.06 |             3.129
exp2                |          35.00 |             0.841
exp2_checked        |          43.06 |             1.399
log2                |          34.23 |             1.556
sin                 |          46.00 |             1.151
sin_checked         |         109.02 |             5.476
cos                 |          54.00 |             1.406
cos_checked         |         113.00 |             4.537
ln                  |          56.91 |             1.626
log10               |          56.91 |             1.626
log1p               |          61.16 |             2.276
exp                 |          39.00 |             0.974
expm1               |          71.00 |             1.441
sinh                |          77.02 |             2.277
cosh                |          48.02 |             2.083
sinh_throughput     |          81.02 |             1.545
cosh_throughput     |          58.00 |             1.279
tanh                |          87.64 |             1.793
asinh               |         120.99 |             7.716
acosh               |         127.75 |             6.588
atanh               |          79.99 |             4.331
asin                |          59.03 |             0.968
acos                |          37.11 |             0.820
atan                |          57.09 |             1.410
atan2               |          57.17 |             1.467
tan                 |          71.02 |             2.532
erf                 |         102.74 |             3.163
erfc                |          78.09 |             2.599
hypot               |          21.11 |             0.766
powf                |          98.03 |             3.898
remainder           |          33.02 |             0.647
```
sin/cos are the restored single-word Cody-Waite version (identical codegen
to before this session's double-float work, confirmed by these numbers
matching exactly); sin_checked/cos_checked are today's name for the
double-float version this whole section was benchmarking.
`sinf_poly`'s `-0.0` sign fix (see the known-defects list above) added a
`copysign` to every one of these five rows -- the real, small, disclosed
cost of that correctness fix (sin/cos throughput +12-16%, tan +9%,
sin_checked/cos_checked +1-4%, latency unchanged except tan +2 cyc); a
branch/select tried first cost meaningfully more on every row and was
rejected in favor of `copysign` once measured.
`log1p`'s own `-0.0` sign fix (see the known-defects list above) is why
`asinh`/`acosh`'s *latency* jumped so much here (55.40->120.99,
110.47->127.75) despite neither function's own `-0.0` handling changing
at all: they already computed their correct sign externally via
`mulsign`, so log1p's fix is pure dead weight for them once inlined, but
unlike sinf_poly's case, no branchless (copysign) form was valid for
`remainder`'s matching fix and no compiling form avoided the cost for
log1p either -- see that entry for the full tradeoff (their *throughput*
improved instead).

**Further pass: `reduce_pi`'s `err0` term deleted via Sterbenz's lemma, not
just downgraded.** `two_sum(x, -p1)` was the one remaining full 6-op merge
left computing a real error term; `p1 = qh*PI_HI` tracks `x` closely enough
(whenever `qh != 0`, comfortably inside Sterbenz's `[x/2, 2x]` exactness
window; trivially when `qh == 0`, since `p1` is then exactly `0`) that the
subtraction is already exact -- `err0` is provably (almost) always zero, not
just noise-tier like the terms downgraded to `quick_two_sum` above. Replaced
with a plain subtract; `err0` drops out of the residual sum entirely.
Verified over three separate full exhaustive accuracy.rs sweeps (every
avg/max ulp and worst-x bit-for-bit identical to the two_sum version, in
every bucket including the off-contract tail) that the one flagged sliver of
doubt -- `|x|` right at the `qh = 0/±1` boundary, where ties-away rounding
could in principle land `p1` a hair outside the Sterbenz window -- never
actually bites. A companion idea from the same brainstorm (rebalancing the
now-4-term `err` sum from a left-to-right chain into a depth-2 tree,
`(e1b + e2b) + (e3b - e3t)`, on the reasoning that shallower must be faster)
was tried alongside it and measured *worse*: identical throughput but +3 cyc
latency both functions (109->112 sin_checked, 113->116 cos_checked) --
apparently freeing `e1b + e2b` from the critical path let LLVM's scheduler
make a different choice elsewhere in the 64-deep latency chain that cost
more than the shorter dependency depth saved. Same lesson as the
bottleneck-analysis note above: a restructuring that looks strictly better
by op-count or depth still needs the measurement. Dropped the rebalance,
kept the flat chain with `err0` gone: net win with no latency cost. mca:
sin_checked 5.544->5.289 cyc/elem throughput (-4.6%), cos_checked
4.919->4.598 cyc/elem throughput (-6.5%), latency unchanged for both
(109.00/113.00 cyc).

**Missed fma contractions swept crate-wide.** Rust never auto-contracts
`a*b + c` into a single `fma` -- every fma in this crate is explicit, and a
handful of spots had missed the memo: `parity()`'s `q - 2.0*floor(...)` (the
intermediate `2.0*floor` is already exact once `q` is an integer, so this is
bit-exact, purely an op-count win, shared by both sin_checked and
cos_checked), `asinh`/`acosh`'s `x*x +/- 1.0`, `hypot`'s `x*x + y*y`, and
`remainder`'s `x - q*y`. All fused (`fma(x, x, 1.0)` etc.), each one fewer
instruction and one fewer rounding on the same critical path. mca:
sin_checked 5.289->5.280 cyc/elem, cos_checked 4.598->4.474 cyc/elem (-2.7%,
from parity() alone), asinh 76.98->71.99 cyc latency / 2.905->2.806 cyc/elem
throughput, acosh 72.06->69.03 cyc / 2.901->2.806 cyc/elem, hypot
25.00->21.00 cyc (-16%) / 0.766->0.763 cyc/elem, remainder 37.00->33.00 cyc
(-11%) / 0.649->0.646 cyc/elem -- real, if individually small, wins across
the board with zero accuracy cost (fma is provably at least as accurate as
two separate roundings; remainder's avg ulp visibly improved too, see the
accuracy table note above). One candidate from the same list, `asin`'s
`a*a - a`, was tried and reverted: bit-for-bit the same fusion, but mca
showed throughput getting *worse* (1.433->1.479 cyc/elem) with latency flat,
reproduced and isolated by reverting that single line while keeping the
rest -- yet another instance of this crate's recurring "fewer ops doesn't
always mean faster" lesson, this time in an otherwise uniformly-positive
batch of near-identical changes.

**tanh's `exp(2.0 * x)` folded its scale into a single multiply.** `exp(y)`
is `exp2(y * LOG2_E)`, so `exp(2.0 * x)` was two runtime multiplies
(`2.0 * x`, then `* LOG2_E`); `2.0 * LOG2_E` is a compile-time constant, so
computing `exp2(x * (2.0 * LOG2_E))` directly drops to one. Bit-exact,
confirmed by an exhaustive (all 2^32 f32 bit patterns) old-vs-new comparison
(temporary test, removed after use, same one-off-verification pattern used
elsewhere in this file) -- `2.0 * x` never rounds (doubling a float is exact
barring overflow), so both forms round the same real-valued product exactly
once. mca: 62.00->58.00 cyc latency (-6.5%), 1.359->1.330 cyc/elem
throughput (-2.1%); sinh/cosh/exp (unaffected controls, don't share this
code path) stayed exactly at baseline.

**atan's range-reciprocal select replaced with a single `min`.** `let y = if
a < 1.0 { a } else { 1.0 / a }` computed `1.0 / a` unconditionally either way
(this crate's branchless/vectorized style evaluates both arms), so the
`if`/`else` was purely a compare+blend choosing between two already-computed
values -- and since `a = x.abs() >= 0`, `a.min(1.0 / a)` picks the exact same
value the select did (`a` itself below 1, the reciprocal at/above 1) in one
`vminps`. The second select (`if a < 1.0 { y } else { FRAC_PI_2 - y }`)
still needs the `a < 1.0` compare, so this removes one blend, not the
compare. Bit-exact, confirmed by an exhaustive (all 2^32 bit patterns)
old-vs-new comparison (temporary test, removed after use; NaN: `a=NaN ->
1/a=NaN -> min(NaN,NaN)=NaN`, matching the old else-branch). Much bigger win
than "save one blend" suggested: mca latency 76.81->57.09 cyc (-25.7%),
throughput 1.984->1.410 cyc/elem (-28.9%) for `atan`; `atan2` (calls `atan`
directly, same reduction applies) 81.00->57.17 cyc (-29.4%) latency,
1.969->1.467 cyc/elem (-25.5%) throughput. `tan` (unaffected control, doesn't
call `atan`) stayed exactly at baseline. The magnitude suggests the old
select was interacting badly with something beyond its own instruction
count (e.g. register pressure or scheduling around the division), not just
costing one op in isolation -- consistent with this crate's recurring
finding that op-count reasoning and measured cost don't always match, this
time in the surprising direction (a small-looking change, a large real win).

# tools
- `cargo +nightly run --release --example accuracy [thorough] [filter]` - avg/max ulp against an f64
  reference. Requires nightly: the reference is computed via the `sleef` crate's SIMD functions (cheapest
  ULP bucket available per function, u35 where it exists -- still ~1e8x tighter than f32 needs), which
  depends on the unstable `portable_simd` feature -- this also means `cargo test`/
  `cargo build --tests` now need nightly, since Cargo builds all dev-dependencies together regardless of
  which target you're building. Default mode fuzzes 100M random f32 bit patterns per function (a few
  seconds); `thorough` exhaustively sweeps all 2^32 bit patterns instead (every denormal, every NaN
  payload, both signs -- a few minutes). Runs on half the machine's cores at low OS scheduling priority
  (`nice`) so it doesn't compete with foreground work while iterating; refuses to run in a debug build.
- `cargo run --release --example quickbench [filter]` - latency (serial dependency chain) + throughput, min of 7 reps
- `cargo run --release --example edgecheck` - bit-exact checks of edge cases (0, -0, denormals, inf, nan, domain boundaries)
- `cargo run --release --example tune` - coordinate-descent ulp tuning of polynomial coefficients
- `cargo run --release --example mca` - theoretical latency/throughput straight from llvm-mca's scheduler
  model for the host CPU (requires `llvm-mca` on PATH). No wall-clock timing, so no thermal-throttling
  noise, and much faster to iterate on than quickbench; see examples/mca_target.rs for the marker
  functions it analyzes and why each region is built the way it is (llvm-mca has no branch predictor,
  so branchy edge-case handling has to be routed around, not just measured through).

# todo:
- do principled and thourough analysis of dependency chains and rounding errors to find optimizations
- perfectly rounded versions
- asin's max-ulp residual (9, down from 84 via three rounds of fixes --
  widening the small-branch threshold, removing the `mid` branch entirely
  in favor of reusing acos's own full-domain poly, then a joint acos/asin
  refit of that poly's coefficients, see asin's doc comment fixes 5-7;
  closing the rest of the way needs a genuinely different correction
  shape, not just a threshold or coefficients)
- asin's latency regressed 43.24->59.03 cyc as a side effect of fix 6
  above (throughput improved a lot, 2.899->0.968 cyc/elem, and accuracy
  improved a lot too; only latency got worse) -- worth a closer look if
  serial-latency-bound callers turn out to care, see the mca notes above
- fix (or at least give a "_checked" full-range companion to) the remaining
  inherited accuracy defects in the newly-ported functions: remainder's
  tie-breaking cliff, and the exp-family's (exp/expm1/sinh/cosh/tanh/powf)
  dependence on the fast unchecked exp2 for their *main* computation, not
  just an edge-case tail (erf/erfc's own tail-specific gaps already fixed,
  see above) -- see the overview above and each function's doc comment
- vary both arguments in quickbench's two-argument benchmarks (atan2, hypot,
  powf, remainder currently fix one argument, which may be letting LLVM
  constant-fold std's side of a couple of comparisons -- see the benchmark
  notes above)
