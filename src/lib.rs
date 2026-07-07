// godbolt flags -C opt-level=3 -C target_feature=+fma

mod doublefloat;
use doublefloat::Df32;

const SIGN_MASK: u32 = 0x80000000;
const EXPONENT_MASK: u32 = 0x7f800000;

#[inline(always)]
fn fma(a: f32, b: f32, c: f32) -> f32 {
    a.mul_add(b, c)
}

#[inline(always)]
pub fn log_2(x: f32) -> f32 {
    // edge handling is done with selects (no early returns) so loops over
    // arrays of log_2 calls can auto-vectorize: scale denormals up before
    // the single normal-path evaluation, then patch specials afterwards.
    let tiny = x < f32::MIN_POSITIVE; // denormal, zero, negative; false for nan
    let xs = if tiny { x * 16777216.0 } else { x };
    let koff = if tiny { -24.0 } else { 0.0 };
    // koff is folded into the exponent term inside log_2_normal so the
    // correction stays off the serial critical path (k + koff is exact)
    let r = log_2_normal(xs, koff);
    // -inf for +-0, nan for x < 0 (includes -inf); the select input only
    // depends on x, so it resolves in parallel with the poly evaluation
    let spec = if x == 0.0 { f32::NEG_INFINITY } else { f32::NAN };
    let r = if x <= 0.0 { spec } else { r };
    // +inf and nan: x*x is inf/nan respectively (false for -inf: -inf < inf)
    if !(x < f32::INFINITY) {
        x * x
    } else {
        r
    }
}

/// Core of log_2 for positive normal finite x only: no handling for zero,
/// negative, denormal, inf, or nan input (those are the caller's job, see
/// log_2). Called directly with an out-of-domain x, this returns a
/// plausible-looking but wrong finite value rather than NaN/-inf.
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn log_2_normal(x: f32, koff: f32) -> f32 {
    // decompose x = 2^k * m with m in [sqrt(2)/2, sqrt(2)), so s = m - 1
    // is exact (Sterbenz) and centered on 0: log2 stays relatively
    // accurate near x = 1. log2(m) = s * P(s), degree-9 minimax P fitted
    // with lolremez (rel. error 4.1e-9).
    let e = (x.to_bits() as i32).wrapping_sub(0x3f3504f3) >> 23; // exponent if m in [√2/2, √2)
    let m = f32::from_bits((x.to_bits() as i32).wrapping_sub(e << 23) as u32);
    let k = e as f32 + koff; // both integers: exact, and off the poly's critical path
    let s = m - 1.0;
    let c: [f32; 10] = [
        std::f32::consts::LOG2_E, // bit-identical to this literal; not a coincidence
        -0.72134733,
        0.4808985,
        -0.36069715,
        0.288568,
        -0.23961738,
        0.20460059,
        -0.19106273,
        0.18617496,
        -0.10994955,
    ];
    let s2 = s * s;
    let s4 = s2 * s2;
    let l0 = fma(c[1], s, c[0]);
    let l1 = fma(c[3], s, c[2]);
    let l2 = fma(c[5], s, c[4]);
    let l3 = fma(c[7], s, c[6]);
    let l4 = fma(c[9], s, c[8]);
    let r0 = fma(l1, s2, l0);
    let r1 = fma(l3, s2, l2);
    let r2 = fma(l4, s4, r1);
    let p = fma(r2, s4, r0);
    // k + s * P(s) in a single rounding
    fma(p, s, k)
}

/// exp2 without domain checks: valid for x in [-126, 128), i.e. normal
/// (non-denormal, finite, nonzero) results only. Outside that range the
/// exponent construction wraps around and the result is garbage (including
/// for nan). Use exp2_checked for full-range handling; this version is
/// ~2.7 ns faster in serial latency.
#[inline(always)]
#[allow(clippy::approx_constant)] // g0's constant term is a fitted minimax
// coefficient near ln(2), not ln(2) itself (bit pattern deliberately differs)
pub fn exp2(x: f32) -> f32 {
    // exp2(floor(x)) * exp2(fract(x)) == exp2(x). exp2int must come from
    // the same floor(x) as f: computing it from x + 383 double-counts the
    // integer part when x + 383 rounds up across an integer (e.g.
    // x = 4.9999999).
    let k = x.floor();
    let f = x - k;
    let exp2int = f32::from_bits(((k + 383_f32).to_bits() << 8) & EXPONENT_MASK);
    // Q(f) = (2^f - 1)/f, degree 5, grouped into 3 balanced pairs (g0, g1,
    // g2) instead of two degree-2 Horner halves: same 6 coefficients (so
    // same accuracy target) and the same 4-deep fma critical path, but the
    // combine only ever needs f^2 (never exp2int*f^4), so it's 2 fewer
    // plain multiplies per call than the old A/B split.
    let f2 = f * f;
    let g0 = fma(2.4022985e-1, f, 6.93147e-1);
    let g1 = fma(9.678826e-3, f, 5.548333e-2);
    let g2 = fma(2.1702237e-4, f, 1.2439679e-3);
    let h = fma(g2, f2, g1);
    let q = fma(h, f2, g0);
    fma(q, exp2int * f, exp2int)
}

#[inline(always)]
#[allow(clippy::approx_constant)] // g0's constant term is a fitted minimax
// coefficient near ln(2), not ln(2) itself (bit pattern deliberately differs)
pub fn exp2_checked(x: f32) -> f32 {
    // fully branchless (auto-vectorizes): exp2(x) = P(f) * 2^k1 * 2^k2 with
    // k1 + k2 = k = floor(x). Splitting k keeps both power-of-two factors
    // representable over the whole clamped range, so overflow to inf and
    // (correctly rounded) denormal underflow fall out of the two multiplies
    // — no pre-offset, no rescale. Both multiplies are exact power-of-two
    // scalings except the final rounding into the denormal range, so the
    // result rounds exactly once. nan propagates through P(f), so there are
    // no fixup selects at all.
    // k must come from the same floor(x) as f: computing the exponent from
    // x + 383 double-counts the integer part when x + 383 rounds up across
    // an integer (e.g. x = 4.9999999).
    let xs = x.clamp(-151.0, 128.0);
    let k = xs.floor();
    let f = xs - k;
    // any split k = k1 + k2 with both halves in valid exponent range works,
    // so k1 = round(xs/2) via the magic constant is fine (and cheap: it
    // runs in parallel with the floor). k1 in [-76, 64], k2 in [-77, 65].
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k1b = fma(xs, 0.5, ROUND_MAGIC) - (ROUND_MAGIC - 383.0); // k1 + 383
    let k2b = (k + 766.0) - k1b; // k2 + 383, exact: all integers
    let t1 = f32::from_bits((k1b.to_bits() << 8) & EXPONENT_MASK);
    let t2 = f32::from_bits((k2b.to_bits() << 8) & EXPONENT_MASK);
    // same 3-balanced-pair Q(f) as exp2 (see there for the derivation): same
    // coefficients/critical-path depth as the old A/B split, 2 fewer plain
    // multiplies (never needs t1*f^4, only t1*f)
    let f2 = f * f;
    let g0 = fma(2.4022985e-1, f, 6.93147e-1);
    let g1 = fma(9.678826e-3, f, 5.548333e-2);
    let g2 = fma(2.1702237e-4, f, 1.2439679e-3);
    let h = fma(g2, f2, g1);
    let q = fma(h, f2, g0);
    // weave t1 into the fma chain (t1*f is exact: both factors normal) so
    // only one multiply (by t2) remains after the polynomial
    let p = fma(q, t1 * f, t1);
    p * t2
}
// sin(x) ~= x + x^3*p(x^2) on [-pi/2, pi/2], degree-9 minimax (relative
// error ~6.1e-9), fitted with lolremez. Estrin evaluation, 2 fma chains.
#[inline(always)]
fn sinf_poly(x: f32) -> f32 {
    let c0 = -0.16666660f32;
    let c1 = 8.3330662e-3f32;
    let c2 = -1.9809603e-4f32;
    let c3 = 2.6057806e-6f32;
    let y = x * x;
    let y2 = y * y;
    let x3 = y * x;
    let a = fma(c1, y, c0);
    let b = fma(c3, y, c2);
    let p = fma(b, y2, a);
    fma(p, x3, x)
}

// pi split into pieces with trailing zero bits so q*PI_A and q*PI_B are
// exact for moderate |q|, keeping the reduced argument accurate in a
// relative sense near the zeros of sin (Cody-Waite with fma).
const PI_A: f32 = 3.140625;
const PI_B: f32 = 0.0009670257568359375;
const PI_C: f32 = 6.277114152908325e-7;
const PI_D: f32 = 1.2154201256553421e-10;
const FRAC_1_PI: f32 = std::f32::consts::FRAC_1_PI;

// 1.5 * 2^23; adding this to |v| < 2^22 rounds v to the nearest integer
// in the low mantissa bits (round-to-nearest-even).
const ROUND_MAGIC: f32 = 12582912.0;

/// sin(x) via single-f32 range reduction: q = round(x/pi) (the round-via-fma
/// magic-constant trick) is exact only while |x| stays under ~1.3e7 (2^22 *
/// pi). Past that, q can land a whole integer off, shifting the residual by
/// a whole multiple of pi and pushing it outside sinf_poly's fitted domain
/// [-pi/2, pi/2] -- including returning inf for some large finite x, since
/// nothing here clamps the residual. Use sin_checked for full-range gradual
/// degradation instead of this cliff; this version is much faster.
#[inline(always)]
pub fn sin(x: f32) -> f32 {
    let qb = fma(x, FRAC_1_PI, ROUND_MAGIC);
    let q = qb - ROUND_MAGIC;
    let r = fma(q, -PI_A, x);
    let r = fma(q, -PI_B, r);
    let r = fma(q, -PI_C, r);
    let r = fma(q, -PI_D, r);
    let s = sinf_poly(r);
    // sin(x) = (-1)^q * sin(r); parity of q is the lowest mantissa bit of qb
    let parity = qb.to_bits() << 31;
    f32::from_bits(s.to_bits() ^ parity)
}
/// Same domain limits as sin (see its doc comment); use cos_checked for
/// full-range gradual degradation.
#[inline(always)]
pub fn cos(x: f32) -> f32 {
    // k = round(x/pi - 0.5), q = k + 0.5, r = x - q*pi in [-pi/2, pi/2]
    let kb = fma(x, FRAC_1_PI, -0.5) + ROUND_MAGIC;
    let q = (kb - ROUND_MAGIC) + 0.5;
    let r = fma(q, -PI_A, x);
    let r = fma(q, -PI_B, r);
    let r = fma(q, -PI_C, r);
    let r = fma(q, -PI_D, r);
    let s = sinf_poly(r);
    // cos(x) = (-1)^(k+1) * sin(r)
    let parity = !kb.to_bits() << 31;
    f32::from_bits(s.to_bits() ^ parity)
}

// q = round(x/pi) must be an *exact* integer for x - q*pi to land accurately
// in [-pi/2, pi/2]. A single-f32 q (tried first, and again after a
// native-round detour -- see jodiemath-workflow memory, 2026-07-06) is a
// binary either-or: either q is exactly the correctly-rounded integer, or
// (once |x| crosses q's exact-integer ceiling) it's off by a whole integer,
// which shifts the residual by a whole multiple of pi and puts sinf_poly
// (fit only for [-pi/2, pi/2]) hopelessly outside its domain -- a relocatable
// *cliff*, not a slope, no matter how q is rounded. This version instead
// gives q a second, small f32 word (qh, ql) -- genuine double-float
// precision via `two_prod`/`two_sum` error-free transforms (exact at any
// magnitude, unlike the crate's other PI_A..D trick, which needs bounded q)
// -- specifically, the dominant cross term and the two next-biggest
// (~x*2^-24) get real two_prod treatment, and the smallest tier (~x*2^-48,
// plus the reciprocal-of-pi's own 3rd correction word) is folded in with
// plain multiplies/adds (its own rounding error there is already far below
// 1 ulp of the O(1) result). A version that dropped that smallest tier
// entirely was tried and measured (via examples/mca.rs) to cost the *same*
// as keeping it (~130-140 cyc either way) -- no meaningful savings once
// genuine multi-word q precision is needed at all, so all three tiers are
// kept here for the best accuracy at no extra cost. See POLY_SAFE_BOUND for
// why the output stays finite even once this gradual degradation is severe.
//
// Bug found while building this (2026-07-06, same day): folding pre_offset
// (cos's -0.5 phase shift) directly into the two_prod's dominant term p0 is
// broken once |x| is large enough that p0's own ulp exceeds 1 (|x| ~
// 1.68e7) -- adding a fixed 0.5 to an already-coarse-ulp float rounds it
// away to nothing, silently dropping cos's phase shift and reducing to the
// wrong (off-by-a-whole-pi) residual. This is the exact same class of bug
// as everything else in this reduction (a small quantity lost against a
// big-ulp value), just one level removed from where the earlier version of
// this bug search was looking. Fixed by folding pre_offset into the *low*
// correction term (comparable magnitude to the other small corrections)
// instead, which stays precise regardless of p0's own ulp. sin (pre_offset
// = 0) was never affected -- adding exactly 0 can't round away -- which is
// why this only showed up in cos's accuracy curve, not sin's.
const RPI_HI: f32 = 0.31830987334251404;
const RPI_LO: f32 = 1.2841276486597053e-8;
const RPI_TINY: f32 = 1.4685477398157775e-16;
const PI_HI: f32 = 3.1415927410125732;
const PI_LO: f32 = -8.742277657347586e-8;
const PI_TINY: f32 = -3.4302490200117637e-15;

#[inline(always)]
fn two_sum(a: f32, b: f32) -> (f32, f32) {
    let s = a + b;
    let v = s - a;
    let e = (a - (s - v)) + (b - v);
    (s, e)
}
// Cheaper 3-op form (a.k.a. Fast2Sum): exact iff |a| >= |b|; when violated,
// `e` is off by up to ~1 ulp of `s` instead of being the exact correction
// (still finite, just not exact -- unlike two_sum, which is exact for any
// a, b, always). Used exactly once below. An initial check (a throwaway
// scratch sweep testing only "near exact multiples of pi") wrongly seemed
// to show this call site always satisfies |a|>=|b| -- a broader, later
// check (uniform random bit patterns, not just near-exact multiples of pi)
// found real violations starting around |x| ~ 1e3 and reaching >80% of
// samples by |x| ~ 1e8+. Kept anyway: re-verified with the *actual*
// accuracy metric that matters (examples/accuracy.rs's exhaustive sweep
// plus a magnitude-bucketed sweep against std out to f32::MAX) shows no
// measurable difference from the full-two_sum version at any magnitude --
// Fast2Sum's bounded (not unbounded) error here is small enough, relative
// to everything else already inexact in this reduction, to be lost in the
// noise. Lesson: "the ordering assumption holds" is a claim about one
// intermediate value, not about the thing that actually matters (the final
// ulp error) -- always re-verify against the real accuracy sweep, since a
// theoretical invariant can be technically false while still being
// practically harmless (or vice versa).
#[inline(always)]
fn quick_two_sum(a: f32, b: f32) -> (f32, f32) {
    let s = a + b;
    let e = b - (s - a);
    (s, e)
}
// error-free product: p+e == a*b exactly, for any a, b (no overflow).
#[inline(always)]
fn two_prod(a: f32, b: f32) -> (f32, f32) {
    let p = a * b;
    let e = fma(a, b, -p);
    (p, e)
}

/// round(x/pi + pre_offset), split into a double-float integer pair
/// (qh, ql). pre_offset is 0 for sin, -0.5 for cos.
#[inline(always)]
fn round_x_over_pi(x: f32, pre_offset: f32) -> (f32, f32) {
    let (p0, e0) = two_prod(x, RPI_HI);
    // quick_two_sum's error term is provably dead: its only use is `s + e1`
    // immediately after, and for a Fast2Sum pair s+e1 == a+b exactly (as
    // reals), so `fl(s + e1)` is just the correctly-rounded a+b again --
    // i.e. `s` itself. True even when the |a|>=|b| assumption is violated
    // (then e1 isn't the exact correction, but s+e1 is still noise-tier
    // relative to s, not a real correction -- verified exhaustively, see
    // IDEAS.md). So the whole quick_two_sum call collapses to a plain add.
    let s = e0 + x * RPI_LO;
    // pre_offset folded in here (not into p0 -- see the bug note above)
    let lo = fma(x, RPI_TINY, s) + pre_offset;
    // ql: ties-to-even instead of f32::round's ties-away-from-zero -- q only
    // needs to be *an* integer within 0.5 of the true residual, so any
    // consistent nearest-rounding rule keeps the qh/ql double-rounding
    // correct in principle (IDEAS.md's claim), and round_ties_even lowers to
    // a single vroundps instead of round's multi-instruction sequence, on
    // ql, the last-ready value out of this function (see reduce_pi).
    //
    // qh stays f32::round (ties-away), NOT round_ties_even, despite the same
    // "any consistent rule" argument applying to it too in principle --
    // measured exception, not applied uniformly: qh's tie-break is shared
    // unmodified between sin (pre_offset=0) and cos (pre_offset=-0.5), and
    // switching qh alone to ties-even regressed cos_checked's max ulp from 2
    // to 6 within its |x|<=1e6 documented-accurate range (worst x ~252.9),
    // while leaving sin_checked untouched -- cos's -0.5 folded into `lo`
    // (not into qh's own input p0) makes qh's rare exact-tie cases interact
    // badly with the offset in a way sin's never does. Isolated via the
    // exhaustive accuracy.rs sweep (both changed vs. each alone): ql alone
    // reproduces zero regression anywhere in-domain for either function
    // (only the already off-contract [1e15,∞) tail moves at all, and only
    // in its already-garbage average, never its max ulp or worst x).
    let qh = p0.round();
    let rem = (p0 - qh) + lo;
    let ql = rem.round_ties_even();
    (qh, ql)
}

/// x - (qh+ql)*pi, accurate well beyond a single f32's exact-integer range.
#[inline(always)]
fn reduce_pi(x: f32, qh: f32, ql: f32) -> f32 {
    let (p1, e1) = two_prod(qh, PI_HI);
    let (p2, e2) = two_prod(qh, PI_LO);
    let (p3, e3) = two_prod(ql, PI_HI);
    // smallest tier (~x*2^-48): a plain multiply/add here is fine, this is
    // exactly the tier the once-adopted full version also just summed
    // in plainly rather than two_sum'ing (see jodiemath-workflow memory)
    let c5 = ql * PI_LO;
    let c45 = fma(qh, PI_TINY, c5);
    let tier2 = (e2 + e3) + c45;
    // p3 and tier2 both depend on ql, the last-ready value out of
    // round_x_over_pi (it needs the whole reduction chain, while qh is
    // ready much earlier from a single multiply+round) -- so p1/p2, and the
    // x/-p1 subtraction below, are all ready well before p3/tier2 are.
    // Combining p3+tier2 into one value here runs fully parallel with the
    // qh-only chain instead of stacking as two more sequential merges after
    // it, shortening the ql-dependent tail by one two_sum's worth of
    // latency (verified via examples/mca.rs). NB: two_sum guarantees
    // p3t + e3t == p3 + tier2 exactly, so subtracting (p3+tier2) means
    // subtracting *both* p3t and e3t -- e3t must be subtracted from err
    // below, not added; got this backwards on the first attempt and it
    // broke cos badly near its zero crossings (where p3 and tier2 can
    // nearly cancel, making e3t large instead of negligible).
    let (p3t, e3t) = two_sum(p3, tier2);
    // x - p1 is exact, not just noise-tier -- by Sterbenz's lemma, not
    // assumption: p1 = qh*PI_HI with qh = round(x/pi) (approximately, via
    // RPI_HI), so whenever qh != 0, p1 sits within a factor of
    // ~(1 +- 2^-24) of x, comfortably inside [x/2, 2x]; when qh == 0, p1 ==
    // 0 and the subtraction is trivially exact too. So the two_sum's error
    // term is provably (almost) always zero and the full 6-op two_sum
    // collapses to a plain subtract (see IDEAS.md). The only sliver of
    // doubt -- |x| right at the qh = 0/+-1 boundary, where ties-away
    // rounding could in principle land p1 a hair outside the Sterbenz
    // window -- was confirmed clean by the exhaustive sweep (examples/
    // accuracy.rs thorough): sin_checked/cos_checked bit-exact identical to
    // the two_sum version in every bucket, including worst-x, at every
    // magnitude tested (in-domain and the off-contract tail alike).
    let s0 = x - p1;
    // e1's merge already used quick_two_sum before this session; p3t's
    // merge (the new one, below) is downgraded the same way -- both violate
    // the |a|>=|b| Fast2Sum ordering assumption somewhere in the domain
    // (verified, not assumed), but empirically cost nothing beyond what's
    // already inside budget (examples/accuracy.rs's exhaustive sweep:
    // unchanged avg ulp for sin, cos |x|<=1e6 avg 0.077->0.081, still ~12x
    // under the 1-ulp budget). p2's merge stays on full two_sum: it's not on
    // the ql-dependent critical path (p2 only needs qh, ready early) so
    // downgrading it saves no latency, only risks accuracy for nothing.
    let (s1, e1b) = quick_two_sum(s0, -e1);
    let (s2, e2b) = two_sum(s1, -p2);
    let (s3, e3b) = quick_two_sum(s2, -p3t);
    // err0 is gone (folded into the exact subtract above) -- one term
    // shorter than before, still left-to-right. Rebalancing this into a
    // depth-2 tree ((e1b + e2b) + (e3b - e3t)) was tried and measured
    // (mca): identical throughput to the flat form here but +3 cyc *worse*
    // latency (109->112 sin_checked, 113->116 cos_checked) -- apparently
    // moving e1b+e2b off the critical path let LLVM's scheduler make a
    // choice elsewhere in the 64-deep chain that cost more than it saved.
    // Same lesson as the Fast2Sum comment above: a "should be strictly
    // better" restructuring needs the measurement, not just the depth
    // count. Kept the flat form: it's simpler and it's what measured best.
    let err = e1b + e2b + e3b - e3t;
    s3 + err
}

// parity of an exact-integer float q via floor-based "mod 2" (q*0.5 and its
// floor stay exact once q is an integer), not `q as i64`: Rust's
// float-to-int cast is saturating, which LLVM can't lower to a single
// vector instruction (confirmed via --emit=asm in an earlier session --
// see jodiemath-workflow memory).
#[inline(always)]
fn parity(q: f32) -> f32 {
    fma(-2.0, (q * 0.5).floor(), q)
}

// Bound for the reduced residual right before it enters sinf_poly. Once the
// dropped smallest-precision-tier terms above start to matter (|x| beyond
// the gradual-degradation range), the residual is no longer close to
// [-pi/2, pi/2] and can grow large -- squaring that inside sinf_poly is
// where the original (pre-2026-07) single-word code's "returns inf for
// ordinary finite input" bug came from. sinf_poly's dominant term for
// large |r| is ~c3*r^9 (c3 ~ 2.6e-6), which overflows f32 around |r| ~ 8e4;
// 1000 leaves a large safety margin while still being far outside
// [-pi/2, pi/2], so a legitimately-reduced residual is never clipped.
// `.clamp` on a NaN residual (x itself nan or +-inf) returns nan unchanged,
// so sin/cos(nan/inf) still correctly come out nan with no extra selects.
const POLY_SAFE_BOUND: f32 = 1000.0;

#[inline(always)]
pub fn sin_checked(x: f32) -> f32 {
    let (qh, ql) = round_x_over_pi(x, 0.0);
    let r = reduce_pi(x, qh, ql).clamp(-POLY_SAFE_BOUND, POLY_SAFE_BOUND);
    // sin(x) = (-1)^q * sin(r); q = qh+ql, so parity(q) = (parity(qh) +
    // parity(ql)) mod 2. parity(qh) and parity(ql) are each exactly 0.0 or
    // 1.0, so their sum mod 2 is just whether they differ (XOR), cheaper
    // than a 3rd floor-based parity() call on the sum.
    //
    // sin is odd, so (-1)^q * sin(r) == sin((-1)^q * r): flip r's sign bit
    // *before* sinf_poly instead of negating its result after. Bit-exact
    // with the old `s * (1.0 - 2.0 * par)` (both IEEE negation and a
    // multiply by exactly +-1 only ever flip the sign bit, never round),
    // but the flip mask depends solely on qh/ql -- ready long before r
    // exits reduce_pi -- so it hides entirely in the reduction's shadow
    // instead of costing a real fma+mul on sinf_poly's tail.
    let pq = parity(qh);
    let pl = parity(ql);
    let flip = if pq == pl { 0 } else { SIGN_MASK };
    let r = f32::from_bits(r.to_bits() ^ flip);
    sinf_poly(r)
}
#[inline(always)]
pub fn cos_checked(x: f32) -> f32 {
    // k = round(x/pi - 0.5), q = k + 0.5, r = x - q*pi in [-pi/2, pi/2]
    let (kh, kl) = round_x_over_pi(x, -0.5);
    // q = k + 0.5; fold the 0.5 into the small word kl, not the (possibly
    // huge) kh word, for the same reason pre_offset itself is folded into
    // the low correction term above -- kl stays small enough that + 0.5
    // is always exact
    let r = reduce_pi(x, kh, kl + 0.5).clamp(-POLY_SAFE_BOUND, POLY_SAFE_BOUND);
    // cos(x) = (-1)^(k+1) * sin(r); k = kh+kl. Same sign-flip-before-the-
    // poly trick as sin_checked above (also odd in r), inverted since the
    // exponent is k+1 instead of k.
    let pk = parity(kh);
    let pl = parity(kl);
    let flip = if pk == pl { SIGN_MASK } else { 0 };
    let r = f32::from_bits(r.to_bits() ^ flip);
    sinf_poly(r)
}

/// Core of cbrt for normal finite x: bit-trick seed (~3% error), then a
/// single degree-3 correction. d = s^3 - x is exact-ish via fma at any
/// scale, and x/s^3 == 1/(1+r) exactly for r = d/x, so
/// cbrt(x) = s * (1+r)^(-1/3), approximated by a minimax poly in r (fitted
/// with lolremez on the seed error range [-0.0999, 0.0894], then
/// coordinate-descent tuned). Degree3 (4 coeffs) instead of degree5 trades
/// avg/max ulp (0.085/1 -> 0.33/2, still inside the 1 avg / 2 max budget)
/// for one fewer fma and one less critical-path depth, worth it: throughput
/// 2.41 -> 1.00 cyc/elem, latency 39 -> 35 cyc (examples/mca.rs).
/// All intermediates are O(1) or O(x): no overflow/underflow anywhere.
/// Callers needing a rescaled result (an exact power of two, or 1.0) must
/// multiply the *return value*, not thread a scale parameter through: an
/// in-function scale param that touches 2 downstream ops (the old code
/// scaled `ss`, which then feeds both `sr` and the final fma) gives LLVM's
/// vectorizer's per-branch constant-folding heuristic enough incentive to
/// fully duplicate this entire function for tiny vs. normal inputs instead
/// of computing once and blending -- confirmed via llvm-mca disassembly
/// (2x the fma/mul count and 2x the divisions in cbrt's throughput region
/// vs. cbrt_accurate_normal, which only touches its scale param in a single
/// final op and doesn't duplicate). A single post-multiply by an exact
/// power of two rounds identically to pre-scaling (same argument as
/// cbrt_accurate_normal's scale), so this is a pure codegen fix: throughput
/// 2.247 -> 1.629 cyc/elem, latency unchanged, no accuracy change
/// (examples/mca.rs).
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn cbrt_normal(x: f32) -> f32 {
    let ax = x.to_bits() & !SIGN_MASK;
    let a = f32::from_bits(ax);
    let rcp = 1.0 / a; // independent of the seed chain, starts immediately
    let s = f32::from_bits(ax / 3 + 0x2a509a07u32);
    let s2 = s * s;
    let d = fma(s2, s, -a);
    let r = d * rcp;
    let c1 = -0.33333164f32;
    let c2 = 0.22220786f32;
    let c3 = -0.17394418f32;
    let c4 = 0.1482371f32;
    let r2 = r * r;
    let a1 = fma(c2, r, c1);
    let b1 = fma(c4, r, c3);
    let p = fma(b1, r2, a1);
    let ss = f32::from_bits(s.to_bits() | (x.to_bits() & SIGN_MASK));
    let sr = ss * r;
    fma(sr, p, ss)
}

#[inline(always)]
pub fn cbrt(x: f32) -> f32 {
    let ax = x.to_bits() & !SIGN_MASK;
    let tiny = ax < 0x0080_0000; // denormal or zero: rescale by 2^24 = (2^8)^3
    let xs = if tiny { x * 16777216.0 } else { x };
    let scale = if tiny { 0.00390625 } else { 1.0 };
    // scale multiplies cbrt_normal's *return value* (see cbrt_normal's doc
    // comment for why threading it through as a parameter instead makes
    // LLVM duplicate the whole function per branch)
    let r = cbrt_normal(xs) * scale;
    // +-0, +-inf, nan propagate (also kills the rcp=inf NaN for x == +-0)
    if ax == 0 || ax >= EXPONENT_MASK { x + x } else { r }
}

/// cbrt to within ~0.5 ulp: cbrt_normal (<= 1 ulp), then one Newton step
/// carried out in double-f32 arithmetic. Only valid for x already rescaled
/// into cbrt_accurate's safe range (roughly 2^-56 to 2^127): outside it the
/// double-f32 residual denormalizes/misrounds, or the Newton step's cube
/// can overflow to inf, silently breaking the ~0.5 ulp guarantee (see
/// cbrt_accurate for the rescale).
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn cbrt_accurate_normal(x: f32, scale: f32) -> f32 {
    let ax = x.to_bits() & !SIGN_MASK;
    let a = f32::from_bits(ax);
    // Depends only on x (not y or e), so this whole expression runs fully
    // parallel with cbrt_normal's seed chain and the double-float cube
    // below -- folding scale in here too (rather than into a separate
    // multiply after y or e are ready) means the only work left on the
    // critical path after e is a single fma, and the only work left after y
    // is a single multiply (den_recip), instead of two of each. scale is an
    // exact power of two, so this rounds identically to scaling afterwards
    // (no intermediate hits the denormal range).
    let neg_rcp3_scale = -((1.0 / a) * (1.0 / 3.0)) * scale;
    let y = cbrt_normal(x);
    let y2 = Df32::from_mul(y, y);
    let y3 = y2 * y;
    // e = y^3 - x, exact-ish: |e| ~ ulp(x)
    let e = (y3.0 - x) + y3.1;
    // 1/(3y^2) without a second hardware division: y^3 ~ x (cbrt_normal's
    // seed is within ~2 ulp) gives 1/y^2 = y/y^3 ~ y/x = |y|/a, so
    // |y|*rcp3 approximates 1/(3y^2). Newton's quadratic convergence only
    // needs den to a handful of accurate bits, not a full division -- this
    // ~2-ulp-relative substitution error is far inside the rounding budget.
    // Bit-exact (0 ulp) against the old division over the full
    // examples/accuracy.rs sweep and examples/edgecheck.rs. Also a real
    // (small) throughput win here since this CPU's FP divider is nearly
    // idle while FMA/mul ports are the bottleneck (examples/mca.rs).
    let neg_den_recip_scale = y.abs() * neg_rcp3_scale;
    fma(e, neg_den_recip_scale, y * scale)
}

#[inline(always)]
pub fn cbrt_accurate(x: f32) -> f32 {
    let ax = x.to_bits() & !SIGN_MASK;
    // below 2^-56 the double-f32 residual would denormalize and misround;
    // rescale by 2^126 = (2^42)^3 (also covers denormals). Above 2^127 the
    // Newton step's y^3 can overflow to inf (NaN out), so rescale down too.
    const SCALE_UP: f32 = f32::from_bits(0x7e80_0000); // 2^126
    const SCALE_UP_OUT: f32 = f32::from_bits(0x2a80_0000); // 2^-42
    const SCALE_DN: f32 = f32::from_bits(0x0080_0000); // 2^-126
    const SCALE_DN_OUT: f32 = f32::from_bits(0x5480_0000); // 2^42
    let small = ax < 0x2380_0000;
    let big = ax >= 0x7f00_0000; // 2^127; inf/nan land here too, fixed up below
    let xs = if small { x * SCALE_UP } else if big { x * SCALE_DN } else { x };
    let scale = if small { SCALE_UP_OUT } else if big { SCALE_DN_OUT } else { 1.0 };
    let r = cbrt_accurate_normal(xs, scale);
    // +-0, +-inf, nan propagate (also kills the rcp=inf NaN for x == +-0)
    if ax == 0 || ax >= EXPONENT_MASK { x + x } else { r }
}

// higher throughput cbrt experiment, 5.5 ulp average error
pub fn cbrt_throughput(x: f32) -> f32 {
    //let r = f32::from_bits(0xd461ff81u32.wrapping_sub((x.to_bits()>>16)*0x5556u32));
    let r = f32::from_bits(0xd461ff81u32.wrapping_sub(x.to_bits() / 3));
    let r = fma(r * r, (r * r) * x, r * f32::from_bits(0x3fb6e3d7));
    let r = fma(r * r, (r * r) * x, r * f32::from_bits(0x3fe09c2a));
    r * r * x
}


pub fn cbrt_approx(x: f32) -> f32 {
	let y = f32::from_bits(0x2a509849u32 + (x.to_bits() / 3));
	let y = (x + 2.*(y*y)*y) / (3.*(y*y));
    (2.*x*y + (y*y)*(y*y))/(x + 2.*(y*y)*y)
    //(x + 2.*(y*y)*y) / (3.*(y*y))
}
pub fn sqrt_approx(x: f32) -> f32 {
    f32::from_bits(0x1FBD22DF + (x.to_bits() >> 1))
}
pub fn rcp_approx(x: f32) -> f32 {
    f32::from_bits(0x7EEF370B - x.to_bits())
}
pub fn exp2_approx(x: f32) -> f32 {
    -f32::from_bits((x + 383.).to_bits() << 8)
}
pub fn log2_approx(x: f32) -> f32 {
    f32::from_bits((x).to_bits() >> 8 | 256_f32.to_bits()) - 383.
}
pub fn rsqrt_approx(x: f32) -> f32 {
    f32::from_bits(0x5F33E79F - (x.to_bits() >> 1))
}

// 50 average ulp error 32 cycle latency 5.5 cycle rthroughput
pub fn cbrt_fast(x: f32) -> f32 {
    let s = f32::from_bits(0x2a4ddef1u32.wrapping_add((x.to_bits()>>16)*0x5556u32));
    let r = f32::from_bits(0x68ff2381u32.wrapping_sub((x.to_bits()>>16)*0xaaacu32));
    let s = fma(s * s, s * -r, fma(r, x, s));
    fma(s * s, s * -r, fma(r, x, s))
}
pub fn cbrt_constant(x: f32, c: &[u32]) -> f32 {
	let y = f32::from_bits(c[0] + (x.to_bits() / 3));
	let y = (x + 2.*(y*y)*y) / (3.*(y*y));
    y

    //h4
    //s - ((s3 - x)*(Df32(s3x.0*16.,s3x.1*16.).quick_add_df(10.*s6) + x2)).div_to_f32(
    //    (15.*s6 + 15.*x2 + 51.*s3x)*s2
    //)

    /*
    let s = f32::from_bits((x.to_bits() / 3).wrapping_add(
        c0
        //0x2a53c472
    )
    );
    let rs = f32::from_bits(
        c1
        //0x543846f5
         .wrapping_sub( x.to_bits()/ 3));
    let s = fma(fma(s*s,-s,x),rs*rs,s);

    let s32 = Df32::from_mul(s,s) * (s*2.);
    s * 2f32 + ((s32*-1.5)*s).div_to_f32(s32.quick_add(x))*/
    //s
}

// x * sign(y): an xor of sign bits, not the same as copysign (which
// replaces x's sign outright -- mulsign(-2,-3) == 2, copysign(-2,-3) == -2).
// Ported from jodiemath's mulsign.
#[inline(always)]
fn mulsign(x: f32, y: f32) -> f32 {
    f32::from_bits(x.to_bits() ^ (y.to_bits() & SIGN_MASK))
}

const LN_2: f32 = std::f32::consts::LN_2;
const LOG10_2: f32 = std::f32::consts::LOG10_2;
const LOG2_E: f32 = std::f32::consts::LOG2_E;
const FRAC_PI_2: f32 = std::f32::consts::FRAC_PI_2;

/// Straight port of jodiemath's logf: log2(x) rescaled by ln(2). Same domain
/// behavior as log_2 (its edge handling covers zero/negative/denormal/inf/nan).
#[inline(always)]
pub fn ln(x: f32) -> f32 {
    log_2(x) * LN_2
}

/// Straight port of jodiemath's log10f: log2(x) rescaled by log10(2).
#[inline(always)]
pub fn log10(x: f32) -> f32 {
    log_2(x) * LOG10_2
}

/// ln(1+x), accurate for small |x| (unlike the naive `ln(1.0 + x)`, which
/// loses x's low bits forming 1.0+x -- the exact case log1p exists to
/// handle -- and rounds to exactly 1.0, hence exactly 0, for |x| below
/// ~6e-8, half of f32's ulp(1.0)). u = 1+x still rounds away those bits,
/// but c = x - (u - 1) recovers the *exact* rounding error (u - 1 is exact
/// by Sterbenz whenever u is within a factor of 2 of 1, i.e. x in roughly
/// [-0.5, 1] -- comfortably covering the whole small-x range this matters
/// for), and d(ln)/du = 1/u folds it back in as one division (idle
/// divider) + one add. u == 0 (x == -1 exactly, log1p's other domain edge)
/// and u == inf (x == inf) both make the correction degenerate to a
/// literal 0/0 or inf-inf-over-inf NaN even though it should contribute
/// nothing there (ln(u) alone is already the correct -inf/+inf) -- both
/// edges collapse c/u itself to NaN, so one `is_finite` check on the
/// already-computed correction (not a separate check on u) suppresses both
/// at once instead of letting it poison the result.
#[inline(always)]
pub fn log1p(x: f32) -> f32 {
    let u = 1.0 + x;
    let c = x - (u - 1.0);
    let corr = c / u;
    let corr = if corr.is_finite() { corr } else { 0.0 };
    ln(u) + corr
}

/// Straight port of jodiemath's expf: exp2(x * log2(e)). Inherits exp2's
/// unchecked domain (see exp2's doc comment): only accurate while
/// x*log2(e) stays in [-126, 128), i.e. roughly x in [-87.3, 88.7) --
/// outside that, exp2's bit-trick construction produces garbage rather
/// than a clamped/overflowed value. Every other function below that's
/// built on exp (expm1, sinh, cosh, tanh, powf, erf, erfc) inherits the
/// same limit; this is a straight port of the C original, which has the
/// identical gap (its own expf also calls the unchecked exp2f).
#[inline(always)]
pub fn exp(x: f32) -> f32 {
    exp2(x * LOG2_E)
}

/// Straight port of jodiemath's expm1f: a Pade approximant near 0 (where
/// exp(x)-1 loses precision to cancellation), exp(x)-1 directly elsewhere.
/// See exp's doc comment for the inherited unchecked-exp2 domain limit.
#[inline(always)]
pub fn expm1(x: f32) -> f32 {
    let a = x * fma(-2.0, x * x, -120.0) / fma(x, fma(x, x - 12.0, 60.0), -120.0);
    let b = exp(x) - 1.0;
    if x.abs() < 0.5 { a } else { b }
}

// sinh(x) = x + x^3/6 + x^5/120 + x^7/5040 + O(x^9), the odd Taylor series
// (exact rational coefficients, not a numerical fit -- sinh is entire, so
// this converges everywhere, and truncation error at the |x|<0.5 select
// boundary below is dominated by the next (dropped) term, x^9/362880 ~
// 5.4e-9 at x=0.5, ~0.1 ulp of sinh(0.5) -- comfortable margin under
// budget for all four kept terms). Same role as expm1's Pade "a" branch:
// a cheap, cancellation-free small-x numerator.
#[inline(always)]
fn sinh_small(x: f32) -> f32 {
    let x2 = x * x;
    let c0 = 1.0f32;
    let c1 = 1.0 / 6.0f32;
    let c2 = 1.0 / 120.0f32;
    let c3 = 1.0 / 5040.0f32;
    let p = fma(fma(fma(c3, x2, c2), x2, c1), x2, c0);
    x * p
}

/// sinh(x) = 0.5*(exp(x) - exp(-x)) directly, except for |x| < 0.5 where
/// exp(x) and exp(-x) are both ~1 and the subtraction cancels almost all
/// precision (the same class of bug log1p/tanh had, see IDEAS.md) --
/// there, use the Taylor form above instead, same branchless-select
/// pattern as expm1's Pade/exp split. See exp's doc comment for the
/// inherited unchecked-exp2 domain limit (only relevant on the `b` side,
/// unconditionally evaluated but only selected for |x| >= 0.5). cosh below
/// doesn't need this: it adds instead of subtracting, so it never cancels.
#[inline(always)]
pub fn sinh(x: f32) -> f32 {
    let a = sinh_small(x);
    let b = 0.5 * (exp(x) - exp(-x));
    if x.abs() < 0.5 { a } else { b }
}

/// Straight port of jodiemath's coshf. See exp's doc comment for the
/// inherited unchecked-exp2 domain limit.
#[inline(always)]
pub fn cosh(x: f32) -> f32 {
    0.5 * (exp(x) + exp(-x))
}

/// Throughput-tier sinh: computes `exp(-x)` as `1.0 / exp(x)` instead of a
/// second full exp evaluation, trading one whole poly evaluation for one
/// division. On this CPU the FP divider is close to idle even when the
/// FMA/mul ports are saturated (same finding behind cbrt_accurate's
/// reciprocal reuse), so a vectorized loop over many elements sees a large
/// throughput win -- but a single serial call now waits on `exp(x)` before
/// the division can even start, where the two independent `exp` calls in
/// `sinh` could previously run in parallel, so *latency* is worse here, not
/// better. Same accuracy as `sinh` for practical purposes (one extra
/// rounding from the division; measured negligible impact, see readme).
/// Use `sinh` for a value on its own or a serial dependency chain, this for
/// a large array/SIMD loop. Mirrors `cosh_throughput` below. Same small-x
/// cancellation fix as `sinh` above (`e` and `1/e` are both ~1 for small
/// x), same select boundary and Taylor branch.
#[inline(always)]
pub fn sinh_throughput(x: f32) -> f32 {
    let a = sinh_small(x);
    let e = exp(x);
    let b = 0.5 * (e - 1.0 / e);
    if x.abs() < 0.5 { a } else { b }
}

/// Throughput-tier cosh: see `sinh_throughput`'s doc comment for the
/// latency/throughput tradeoff this shares (same `1.0/exp(x)` reuse, same
/// reasoning, same caveat).
#[inline(always)]
pub fn cosh_throughput(x: f32) -> f32 {
    let e = exp(x);
    0.5 * (e + 1.0 / e)
}

/// tanh(x) = (e^2x - 1) / (e^2x + 1) = expm1(2x) / (expm1(2x) + 2), reusing
/// expm1's already-correct small-x handling (its own Pade branch below
/// |x|<0.5) instead of computing exp2(2x) and cancelling `1.0 - (~1.0)`
/// directly, which lost essentially all precision for small x (a fuzz
/// sweep found this the same 300-million-ulp-average class of bug as
/// log1p's, before that fix -- see IDEAS.md). Still inherits exp's
/// unchecked-domain limit via expm1 (halved range, same as before: expm1
/// sees 2x).
#[inline(always)]
pub fn tanh(x: f32) -> f32 {
    let e = expm1(2.0 * x);
    e / (e + 2.0)
}

/// Straight port of jodiemath's asinhf: ln(x + sqrt(x^2+1)), inherited as-is
/// from the C original including its known flaw -- for small |x| (below
/// ~6e-8, half of f32's ulp(1.0)), x + sqrt(x^2+1) rounds down to exactly
/// 1.0, so ln(...) returns exactly 0 instead of the (tiny but nonzero)
/// correct answer. A real fix needs a log1p-based small-x branch (like a
/// proper libm asinh), which is beyond a straight port; this just documents
/// the cliff so it isn't mistaken for a translation bug.
#[inline(always)]
pub fn asinh(x: f32) -> f32 {
    ln(x + fma(x, x, 1.0).sqrt())
}

/// Straight port of jodiemath's acoshf: ln(x + sqrt(x^2-1)), inherited as-is
/// including its known flaw -- squaring x erases its sign before the
/// domain check, so for large-magnitude *negative* x (where the true
/// answer is NaN, acosh's domain is x >= 1) this returns +inf instead:
/// x*x overflows to +inf the same for either sign, and x + inf is +inf
/// regardless of x's sign, so the NaN that a correctly-signed negative
/// sqrt argument would otherwise produce never happens.
#[inline(always)]
pub fn acosh(x: f32) -> f32 {
    ln(x + fma(x, x, -1.0).sqrt())
}

/// Straight port of jodiemath's atanhf: 0.5*ln((1+x)/(1-x)), inherited as-is
/// including its known flaw -- same near-zero cancellation as asinh above
/// (for tiny |x|, (1+x)/(1-x) rounds to exactly 1.0, so ln(...) is exactly
/// 0 instead of the correct tiny nonzero answer). A real fix needs a
/// log1p-based small-x branch, beyond a straight port.
#[inline(always)]
pub fn atanh(x: f32) -> f32 {
    0.5 * ln((1.0 + x) / (1.0 - x))
}

// degree-6 minimax poly (Estrin via fma), fitted for acos's sqrt(1-|x|)
// factor. Ported from jodiemath's acosf_poly.
#[inline(always)]
fn acos_poly(x: f32) -> f32 {
    let u = 2.2960134e-3f32;
    let u = fma(u, x, -1.1146357e-2);
    let u = fma(u, x, 2.6900099e-2);
    let u = fma(u, x, -4.8802612e-2);
    let u = fma(u, x, 8.875567e-2);
    let u = fma(u, x, -2.1458527e-1);
    fma(u, x, 1.5707962)
}

/// Straight port of jodiemath's acosf.
#[inline(always)]
pub fn acos(x: f32) -> f32 {
    const PI: f32 = std::f32::consts::PI;
    let a = x.abs();
    let y = (1.0 - a).sqrt() * acos_poly(a);
    mulsign(y, x) + if x < 0.0 { PI } else { 0.0 }
}

/// Straight port of jodiemath's asinf: a rational correction folded into
/// the acos-style sqrt identity, inherited as-is including its known flaw --
/// the same near-zero cancellation as asinh/atanh (sqrt(1-a)-1 loses
/// precision as a -> 0), so for small |x| this returns a coarser answer
/// than the ~1e-6-relative-error goal the rest of jodiemath aims for. A
/// real fix needs a small-x series branch, beyond a straight port.
#[inline(always)]
pub fn asin(x: f32) -> f32 {
    let a = x.abs();
    let d = fma(-0.0392588, a, 0.179323);
    let d = fma(-a, d, 1.75866);
    let d = fma(-a, d, -3.66063);
    let a = (a * a - a) / d + a;
    mulsign((1.0 - a).sqrt() - 1.0, x) * (-FRAC_PI_2)
}

// Pade-style rational approximation of atan on [0,1]. Ported from
// jodiemath's atanf_poly.
#[inline(always)]
fn atan_poly(x: f32) -> f32 {
    let a = f32::from_bits(0x3d267031);
    let b = f32::from_bits(0x3f28513c);
    let c = f32::from_bits(0x3e2f725f);
    let d = f32::from_bits(0x3f7da425);
    let x2 = x * x;
    (fma(fma(a, x2, b), x2, 1.0) * x) / fma(fma(x2, c, d), x2, 1.0)
}

/// Straight port of jodiemath's atanf: reciprocates |x| > 1 into range
/// (atan(x) = pi/2 - atan(1/x)) before the poly, matching atan_poly's fit.
#[inline(always)]
pub fn atan(x: f32) -> f32 {
    let a = x.abs();
    // a >= 0, so min(a, 1/a) picks whichever branch the old a<1.0 select
    // did (a itself below 1, the reciprocal at/above 1) in one vminps
    // instead of a compare+blend; the reciprocal was already computed
    // unconditionally either way (both branches evaluate in the
    // branchless/vectorized style this crate uses). NaN: a=NaN -> 1/a=NaN
    // -> min(NaN, NaN) = NaN, matching the old else-branch's 1.0/NaN.
    let y = a.min(1.0 / a);
    let y = atan_poly(y);
    let y = if a < 1.0 { y } else { FRAC_PI_2 - y };
    mulsign(y, x)
}

/// Straight port of jodiemath's atan2f.
#[inline(always)]
pub fn atan2(y: f32, x: f32) -> f32 {
    let nonzerox = x != 0.0;
    let nonzeroy = y != 0.0;
    let bothzero = !nonzerox && !nonzeroy;
    let hpisignx = if nonzerox || bothzero { mulsign(FRAC_PI_2, x) } else { 0.0 };
    let base = if nonzerox { atan(y / x) } else { 0.0 };
    base + mulsign(FRAC_PI_2 - hpisignx, y)
}

/// Straight port of jodiemath's tanf: sin(x)/cos(x), same domain limits as
/// this crate's sin/cos (see their doc comments).
#[inline(always)]
pub fn tan(x: f32) -> f32 {
    sin(x) / cos(x)
}

// degree-6 minimax poly feeding erf's exp2-based tail (|x| >= 0.28). Ported
// from jodiemath's erff_poly.
#[inline(always)]
fn erf_poly(x: f32) -> f32 {
    let u = 3.118769e-4f32;
    let u = fma(u, x, -4.67225e-3);
    let u = fma(u, x, 3.3162573e-2);
    let u = fma(u, x, -1.5214339e-1);
    let u = fma(u, x, -9.1684705e-1);
    let u = fma(u, x, -1.6282598);
    fma(u, x, 3.1332566e-5)
}

/// Straight port of jodiemath's erff: a Pade approximant near 0 (where the
/// tail form loses precision to cancellation), the exp2-based tail elsewhere.
/// The tail branch's exp2 call inherits exp2's unchecked domain (see exp's
/// doc comment) once erf_poly(|x|)'s degree-6 growth pushes its argument
/// out of [-126, 128) -- only relevant once erf has long since saturated to
/// +-1 at f32 precision (|x| beyond roughly 4), so it doesn't affect any
/// input where the answer isn't already indistinguishable from +-1.
#[inline(always)]
pub fn erf(x: f32) -> f32 {
    let x2 = x * x;
    let numer = x * fma(f32::from_bits(0x3f174f6e), x2, f32::from_bits(0x3f906ebb));
    let denom = fma(fma(f32::from_bits(0x3e3e2be3), x2, f32::from_bits(0x3f5b6db7)), x2, 1.0);
    let a = numer / denom;
    let b = mulsign(1.0 - exp2(erf_poly(x.abs())), x);
    if x.abs() < 0.28 { a } else { b }
}

/// Straight port of jodiemath's erfcf: a rational*gaussian tail, clamped to
/// |x| <= 10 before evaluation (matching the C original). That clamp is
/// meant to keep exp(-xa*xa) in exp's safe range, but doesn't fully manage
/// it: for |x| >= ~9.35, xa*xa >= ~87.4, and exp(-87.4) needs
/// exp2(-87.4*log2(e)) =~ exp2(-126.1) -- already at, or just past, exp2's
/// unchecked [-126, 128) domain (see exp's doc comment), so the result is
/// unreliable (not a clean 0) for that whole tail rather than just far out
/// past f32's underflow point. Inherited as-is from the C original (same
/// clamp, same gap); a real fix would need exp2_checked here instead.
#[inline(always)]
pub fn erfc(x: f32) -> f32 {
    let z = if x < 0.0 { -1.0 } else { 1.0 };
    let w = if x < 0.0 { 2.0 } else { 0.0 };
    // NaN-preserving clamp: f32::min suppresses NaN (returns the other
    // operand), unlike C's `x>10.f?10.f:x` ternary (false for NaN, so it
    // takes the x branch, keeping NaN). This if/else matches the ternary.
    let xa = x.abs();
    let xa = if xa > 10.0 { 10.0 } else { xa };
    let n = fma(f32::from_bits(0x35c42f59), xa, f32::from_bits(0x3daf42dc));
    let n = fma(n, xa, f32::from_bits(0x3ee32e3c));
    let n = fma(n, xa, f32::from_bits(0x3f7a7520));
    let n = fma(n, xa, 1.0);
    let d = fma(f32::from_bits(0x3e1b69eb), xa, f32::from_bits(0x3f48fde0));
    let d = fma(d, xa, f32::from_bits(0x3fe918cc));
    let d = fma(d, xa, f32::from_bits(0x4006d465));
    let d = fma(d, xa, 1.0);
    let y = exp(-(xa * xa)) * n / d;
    fma(y, z, w)
}

/// Straight port of jodiemath's hypotf: naive sqrt(x^2+y^2), no anti-overflow
/// rescaling (unlike std's hypot) -- trades the overflow/underflow edge cases
/// for vectorizability, same tradeoff this crate makes for cbrt/sin/cos vs.
/// their std counterparts.
#[inline(always)]
pub fn hypot(x: f32, y: f32) -> f32 {
    fma(x, x, y * y).sqrt()
}

/// Straight port of jodiemath's powf: exp2(log2(x) * y). Inherits exp2's
/// unchecked domain (see exp's doc comment): only accurate while
/// log2(x)*y stays in [-126, 128).
#[inline(always)]
pub fn powf(x: f32, y: f32) -> f32 {
    exp2(log_2(x) * y)
}

/// Straight port of jodiemath's remainderf: x - round(x/y)*y (ties away from
/// zero, via f32::round -- not IEEE 754 remainder's ties-to-even).
/// round(x/y)*y's absolute error scales with ulp(x), which swamps the true
/// remainder (at most |y|/2) once |x/y| is large -- inherited from the C
/// original's identical formula, only reliable while |x/y| stays moderate.
#[inline(always)]
pub fn remainder(x: f32, y: f32) -> f32 {
    let q = (x / y).round();
    fma(-q, y, x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;
    use rand::RngExt;

    fn run_descent(
        f: impl Fn(f32, &[u32]) -> f32,
        reference: impl Fn(f32) -> f32,
        initial_consts: &[u32],
    ) {
        let mut consts: Vec<u32> = initial_consts.to_vec();
        let iters = 10_000;

        let mut best_err: u64 = 0;
        let mut steps: u64 = 0;
        for _ in 1..iters {
            let x: f32 = rand::rng().random::<f32>().abs();
            let result = f(x, &consts);
            best_err += reference(x).to_bits().abs_diff(result.to_bits()) as u64;
        }
        let mut tryagain = false;
        println!(
            "optimizing! starting error: {}",
            best_err as f64 / iters as f64
        );
        let mut deltas: Vec<u32> = consts
            .iter()
            .map(|&_| {
                let n: u32 = rand::rng().random();
                f32::from_bits(n) as i32 as u32
            })
            .collect();
        loop {
            let mut new_err: u64 = 0;
            let mut old_err: u64 = 0;
            if !tryagain {
                deltas = consts
                    .iter()
                    .map(|&_| {
                        let n: u32 = rand::rng().random();
                        //return f32::from_bits(n) as i32 as u32;
                        if steps < 10 {
                            n
                        } else {
                            if steps & 1 == 0 {
                                n
                            } else {
                                f32::from_bits(n) as i32 as u32
                            }
                        }
                    })
                    .collect();
            }
            steps += 1;
            let new_consts: Vec<u32> = consts
                .iter()
                .zip(deltas.iter())
                .map(|(&c, &d)| c.wrapping_add(d))
                .collect();
            if new_consts.iter().zip(consts.iter()).all(|(&x, &y)| x == y) {
                continue;
            }
            for _ in 1..iters {
                let x: f32 = rand::rng().random::<f32>().abs();
                let ref_val = reference(x);
                let old_result = f(x, &consts);
                let new_result = f(x, &new_consts);
                new_err += ref_val.to_bits().abs_diff(new_result.to_bits()) as u64;
                old_err += ref_val.to_bits().abs_diff(old_result.to_bits()) as u64;
            }
            if new_err <= old_err {
                consts = new_consts.clone();
                if new_err < best_err {
                    let mut new_err_nomul: u64 = 0;
                    // no mul check
                    for _ in 1..iters {
                        let x: f32 = rand::rng().random::<f32>().abs();
                        let new_result = f(x, &new_consts);
                        new_err_nomul +=
                            reference(x).to_bits().abs_diff(new_result.to_bits()) as u64;
                    }
                    if new_err_nomul < best_err {
                        best_err = new_err_nomul;
                        let const_strs: Vec<String> =
                            consts.iter().map(|c| format!("0x{:x}u32", c)).collect();
                        println!(
                            "new best consts {} with error {} nomul: {} step:{}",
                            const_strs.join(", "),
                            new_err as f64 / iters as f64,
                            new_err_nomul as f64 / iters as f64,
                            steps
                        );
                    }
                    tryagain = true;
                } else {
                    tryagain = false;
                }
            } else {
                tryagain = false;
            }
            if new_err == best_err {
                best_err = new_err;
                consts = new_consts;
            }
        }
    }
    /*
    #[test]
    fn descent2() {
        run_descent(
            |x, consts| cbrt_constant(x, consts),
            |x| (x as f64).cbrt() as f32,
            &[0x2a5063f7],
        );
    }*/

    #[test]
    fn it_works() {
        assert_eq!(log_2(1.0), 0.0);
        assert_eq!(log_2(2.0), 1.0);
        assert_eq!(log_2(4.0), 2.0);
        assert_eq!(log_2(8.0), 3.0);
    }

    fn ulp_error(
        range: std::ops::Range<i32>,
        scale: f64,
        f: impl Fn(f32) -> f32,
        reference: impl Fn(f64) -> f64,
    ) -> f32 {
        let count = (range.end - range.start) as f64;
        let err: u64 = range
            .map(|x| {
                let xf = x as f64 * scale;
                let ref_val = reference(xf) as f32;
                ref_val.to_bits().abs_diff(f(xf as f32).to_bits()) as u64
            })
            .sum();
        err as f32 / count as f32
    }

    #[test]
    fn cbrt_precision() {
        println!(
            "jodie cbrt error: {}",
            ulp_error(1..10000, 1.0, cbrt, |x| x.cbrt())
        );
        println!(
            "std   cbrt error: {}",
            ulp_error(1..10000, 1.0, |x| x.cbrt(), |x| x.cbrt())
        );
    }
    #[test]
    fn cbrt_accurate_precision() {
        println!(
            "jodie cbrt accurate error: {}",
            ulp_error(1..10000, 1.0, cbrt_accurate, |x| x.cbrt())
        );
        println!(
            "std   cbrt error: {}",
            ulp_error(1..10000, 1.0, |x| x.cbrt(), |x| x.cbrt())
        );
    }
    #[test]
    fn cbrt_throughput_precision() {
        println!(
            "jodie cbrt throughput error: {}",
            ulp_error(1..10000, 1.0, cbrt_throughput, |x| x.cbrt())
        );
        println!(
            "std   cbrt error: {}",
            ulp_error(1..10000, 1.0, |x| x.cbrt(), |x| x.cbrt())
        );
    }
    #[test]
    fn exp2_precision() {
        println!(
            "jodie exp2 error: {}",
            ulp_error(-100..100, 0.01, exp2, |x| x.exp2())
        );
        println!(
            "std   exp2 error: {}",
            ulp_error(-100..100, 0.01, |x| x.exp2(), |x| x.exp2())
        );
    }
    #[test]
    fn log2_precision() {
        println!(
            "jodie log2 error: {}",
            ulp_error(2..1000, 0.1, log_2, |x| x.log2())
        );
        println!(
            "std   log2 error: {}",
            ulp_error(2..1000, 0.1, |x| x.log2(), |x| x.log2())
        );
    }
    #[test]
    fn sin_precision() {
        println!(
            "jodie sin error: {}",
            ulp_error(-100..100, 0.01, sin, |x| x.sin())
        );
        println!(
            "std   sin error: {}",
            ulp_error(-100..100, 0.01, |x| x.sin(), |x| x.sin())
        );
    }
    #[test]
    fn cos_precision() {
        println!(
            "jodie cos error: {}",
            ulp_error(-100..100, 0.01, cos, |x| x.cos())
        );
        println!(
            "std   cos error: {}",
            ulp_error(-100..100, 0.01, |x| x.cos(), |x| x.cos())
        );
    }

    fn plot_approx(
        path: &str,
        x_start: f32,
        x_end: f32,
        approx: impl Fn(f32) -> f32,
        truth: impl Fn(f32) -> f32,
    ) {
        use plotters::prelude::*;
        let xs: Vec<f32> = (0..1000)
            .map(|i| x_start + (x_end - x_start) * i as f32 / 999.0)
            .collect();
        let all_y: Vec<f32> = xs.iter().flat_map(|&x| [approx(x), truth(x)]).collect();
        let y_min = all_y
            .iter()
            .cloned()
            .filter(|y| y.is_finite())
            .fold(f32::INFINITY, f32::min);
        let y_max = all_y
            .iter()
            .cloned()
            .filter(|y| y.is_finite())
            .fold(f32::NEG_INFINITY, f32::max);
        let root = BitMapBackend::new(path, (480, 480)).into_drawing_area();
        root.fill(&WHITE).unwrap();
        let mut chart = ChartBuilder::on(&root)
            .margin(5)
            .x_label_area_size(30)
            .y_label_area_size(30)
            .build_cartesian_2d(x_start..x_end, y_min..y_max)
            .unwrap();
        chart.configure_mesh().draw().unwrap();
        chart
            .draw_series(LineSeries::new(xs.iter().map(|&x| (x, approx(x))), &BLACK))
            .unwrap()
            .label("integer approx")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLACK));
        chart
            .draw_series(LineSeries::new(xs.iter().map(|&x| (x, truth(x))), &RED))
            .unwrap()
            .label("ground truth")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], RED));
        chart.configure_series_labels().draw().unwrap();
        root.present().expect("Unable to write result to file");
        std::process::Command::new("pngquant")
            .args(["--force", "--ext", ".png", "16", "--", path])
            .status()
            .unwrap();
    }

    fn plot_error(path: &str, x_start: f32, x_end: f32, f: impl Fn(f32) -> f32) {
        use plotters::prelude::*;
        let samples: Vec<(f32, f32)> = (0..1000)
            .map(|i| x_start + (x_end - x_start) * i as f32 / 999.0)
            .map(|x| (x, f(x)))
            .filter(|&(_, y)| y.is_finite())
            .collect();
        let y_min = samples
            .iter()
            .map(|&(_, y)| y)
            .fold(f32::INFINITY, f32::min);
        let y_max = samples
            .iter()
            .map(|&(_, y)| y)
            .fold(f32::NEG_INFINITY, f32::max);
        let root = BitMapBackend::new(path, (480, 480)).into_drawing_area();
        root.fill(&WHITE).unwrap();
        let mut chart = ChartBuilder::on(&root)
            .margin(5)
            .x_label_area_size(30)
            .y_label_area_size(50)
            .build_cartesian_2d(x_start..x_end, y_min..y_max)
            .unwrap();
        chart
            .configure_mesh()
            .y_label_formatter(&|y| format!("{:.2e}", y))
            .draw()
            .unwrap();
        chart.draw_series(LineSeries::new(samples, &BLACK)).unwrap();
        root.present().expect("Unable to write result to file");
        std::process::Command::new("pngquant")
            .args(["--force", "--ext", ".png", "16", "--", path])
            .status()
            .unwrap();
    }

    #[test]
    fn cbrt_approx_plot() {
        plot_approx("cbrt_approx.png", 1., 128., cbrt_approx, |x| x.cbrt());
    }
    #[test]
    fn sqrt_approx_plot() {
        plot_approx("sqrt_approx.png", 1., 128., sqrt_approx, |x| x.sqrt());
    }
    #[test]
    fn rcp_approx_plot() {
        plot_approx("rcp_approx.png", 1., 10., rcp_approx, |x| 1.0 / x);
    }
    #[test]
    fn exp2_approx_plot() {
        plot_approx("exp2_approx.png", 0., 10., exp2_approx, |x| x.exp2());
    }
    #[test]
    fn log2_approx_plot() {
        plot_approx("log2_approx.png", 1., 128., log2_approx, |x| x.log2());
    }
    #[test]
    fn sin_plot() {
        plot_approx("sin.png", -20., 20., sin, |x| x.sin());
    }
    #[test]
    fn cos_plot() {
        plot_approx("cos.png", -20., 20., cos, |x| x.cos());
    }
    #[test]
    fn rsqrt_approx_plot() {
        plot_approx("rsqrt_approx.png", 1., 128., rsqrt_approx, |x| {
            1.0 / x.sqrt()
        });
    }

    #[test]
    fn log_2_error() {
        plot_error("log_2_error.png", 1., 128., |x| {
            log_2(x) / (x as f64).log2() as f32 - 1.0
        });
    }
    #[test]
    fn exp2_error() {
        plot_error("exp2_error.png", 0., 10., |x| {
            exp2(x) / (x as f64).exp2() as f32 - 1.0
        });
    }
    #[test]
    fn sin_error() {
        plot_error("sin_error.png", -20., 20., |x| {
            sin(x) / (x as f64).sin() as f32 - 1.0
        });
    }
    #[test]
    fn cos_error() {
        plot_error("cos_error.png", -20., 20., |x| {
            cos(x) / (x as f64).cos() as f32 - 1.0
        });
    }
    #[test]
    fn cbrt_error() {
        plot_error("cbrt_error.png", 1., 128., |x| {
            cbrt(x) / (x as f64).cbrt() as f32 - 1.0
        });
    }
    #[test]
    fn cbrt_accurate_error() {
        plot_error("cbrt_accurate_error.png", 1., 128., |x| {
            cbrt_accurate(x) / (x as f64).cbrt() as f32 - 1.0
        });
    }
    #[test]
    fn cbrt_approx_error() {
        plot_error("cbrt_approx_error.png", 1., 128., |x| {
            cbrt_approx(x) / (x as f64).cbrt() as f32 - 1.0
        });
    }
}
