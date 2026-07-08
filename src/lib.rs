// godbolt flags -C opt-level=3 -C target_feature=+fma

// This crate's entire accuracy and performance story depends on `fma`
// (see the function just below) compiling to a single hardware
// instruction. Without FMA, `f32::mul_add` falls back to a ~2x-slower
// libm call that also rounds *differently* (two roundings instead of
// one) -- every ulp figure in this crate's doc comments, readme.md, and
// examples/accuracy.rs assumes the single-rounding hardware form.
// `.cargo/config.toml` sets `target-cpu=native` for exactly this reason,
// but that setting is silently overridden (not merged) by a `RUSTFLAGS`
// environment variable, a well-known Cargo gotcha -- some CI/build
// setups export `RUSTFLAGS` directly, which would disable FMA with no
// build error and no runtime symptom beyond quietly-wrong accuracy
// numbers. Fail loudly at compile time instead: `target_feature = "fma"`
// is set by the compiler whenever FMA is actually enabled, regardless of
// how that happened (config.toml, RUSTFLAGS, --target, etc.), so this
// check is robust to all of them.
#[cfg(all(not(target_feature = "fma"), not(doctest)))]
compile_error!(
    "jodiemath-rs requires hardware FMA (target-feature=+fma or target-cpu=native) -- \
     without it, f32::mul_add falls back to a slower, differently-rounded software path \
     and every accuracy/perf figure in this crate's docs is invalid. Build with \
     `RUSTFLAGS=\"-C target-cpu=native\"` or ensure .cargo/config.toml's rustflags \
     aren't being overridden by an environment RUSTFLAGS variable."
);

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
//
// At x = +-0.0, x3 = y*x always carries x's own sign (y = x*x is always
// +0.0, and +0.0 * x doesn't flip sign), but p (the poly's leading
// coefficient c0 at y=0) is a fixed negative constant -- sin's own
// curvature -- so p*x3 always ends up with the *opposite* sign to x at
// this one point. `fma(p, x3, x)` then adds two exactly-zero values of
// opposite sign, which IEEE754 defines to give +0.0 regardless of
// operand order or x's own sign, silently losing it (the same mechanism
// behind the atan2(-0.0,+0.0) bug fixed earlier). A branchy `x == 0.0`
// select fixed it but cost real throughput (measured, ~12-17% worse on
// sin/cos/tan) since it's inlined into every caller. `r.copysign(x)`
// fixes the same bug for free: for every *nonzero* x in this poly's
// domain, sin is odd and monotonic so the leading `x` term always
// dominates the correction term `p*x3` in magnitude, meaning r's sign
// already equals x's sign there -- copysign is a true no-op for all of
// them and only changes the singular x=+-0.0 case. Bit-exact vs. the
// branchy version everywhere, cheaper (single sign-copy instruction, no
// compare/select) on every shared caller (sin, cos, sin_checked,
// cos_checked, tan).
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
    let r = fma(p, x3, x);
    r.copysign(x)
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
    let result = sinf_poly(r);
    // reduce_pi's own multi-term error-compensation chain loses x's sign
    // at x = +-0.0 (an opposite-signed-zero addition somewhere inside it,
    // the same IEEE754 mechanism as sinf_poly's own -0.0 fix and the
    // atan2(-0.0,+0.0) bug), well before sinf_poly ever sees it -- guard
    // here rather than trace through reduce_pi's whole two_sum/two_prod
    // chain to find the exact spot.
    if x == 0.0 { x } else { result }
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
/// (examples/mca.rs). The 4 correction coefficients were refit again
/// (2026-07-07) with examples/tune.rs against the *existing* seed (the
/// seed itself -- swapping to a cheaper bit-trick -- is a separate,
/// bigger question, see IDEAS.md): avg ulp 0.3265 -> 0.3112 (~4.7%), max
/// ulp unchanged at 3, zero perf cost (same instructions, only the 4
/// literal constants differ).
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
    let c1 = -0.33333147f32;
    let c2 = 0.22220612f32;
    let c3 = -0.17394388f32;
    let c4 = 0.14823665f32;
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
const LOG2_E: f32 = std::f32::consts::LOG2_E;
const FRAC_PI_2: f32 = std::f32::consts::FRAC_PI_2;
const FRAC_PI_4: f32 = std::f32::consts::FRAC_PI_4;

// Cody-Waite split of ln(2): LN2_HI keeps its low 9 mantissa bits zeroed, so
// k*LN2_HI (k an exact small integer, this crate's log_2_normal decomposition
// never produces |k| past a couple hundred) is *exact* -- no rounding at all,
// confirmed by brute force for k in [-300, 300]. LN2_LO is the f32-rounded
// residual (LN2 - LN2_HI as f64, then rounded). This is the same trick
// reduce_pi already uses for pi/2's own hi/lo split.
const LN2_HI: f32 = 0.693145751953125;
const LN2_LO: f32 = 1.428606765330187e-6;
const LOG10_2_HI: f32 = 0.301025390625;
const LOG10_2_LO: f32 = 4.605039066518657e-6;

/// ln(x) via a poly fitted directly for ln, not log_2's poly rescaled after
/// the fact. `log_2(x) * LN_2` (the naive approach, and jodiemath's own
/// original formula) rounds *twice*: once inside log_2 to produce its own
/// f32 result, then again multiplying that already-rounded value by LN_2 --
/// and that second rounding applies to the *whole* result (dominated by the
/// integer exponent term k, not the small poly correction), so it costs
/// nearly a full ulp of avoidable error. Same domain behavior as log_2 (its
/// edge handling covers zero/negative/denormal/inf/nan).
#[inline(always)]
pub fn ln(x: f32) -> f32 {
    let tiny = x < f32::MIN_POSITIVE;
    let xs = if tiny { x * 16777216.0 } else { x };
    let koff = if tiny { -24.0 } else { 0.0 };
    let r = ln_normal(xs, koff);
    let spec = if x == 0.0 { f32::NEG_INFINITY } else { f32::NAN };
    let r = if x <= 0.0 { spec } else { r };
    if !(x < f32::INFINITY) {
        x * x
    } else {
        r
    }
}

/// Core of ln for positive normal finite x only -- see log_2_normal, same
/// contract. Reuses log_2_normal's exact decomposition (s = m - 1) and
/// poly *shape*, but with coefficients fitted for ln directly (log_2's own
/// c[i] * LN_2, each individually rounded to f32) and a Cody-Waite combine
/// with k instead of a single fma: `k*LN2_HI` is exact (see LN2_HI's own
/// comment), so `fma(poly, s, k*LN2_HI)` folds the poly correction in with
/// only one rounding, then `+ k*LN2_LO` adds back the tiny residual LN2_HI
/// dropped -- that final add's own rounding now only affects a small
/// correction term instead of the whole (k-dominated) result, unlike the
/// naive `log_2(x) * LN_2` where the second rounding scales everything.
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn ln_normal(x: f32, koff: f32) -> f32 {
    let e = (x.to_bits() as i32).wrapping_sub(0x3f3504f3) >> 23;
    let m = f32::from_bits((x.to_bits() as i32).wrapping_sub(e << 23) as u32);
    let k = e as f32 + koff;
    let s = m - 1.0;
    let c: [f32; 10] = [
        1.0,
        -0.49999988,
        0.33333343,
        -0.25001621,
        0.20002009,
        -0.16609012,
        0.14181833,
        -0.13243459,
        0.12904665,
        -0.07621122,
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
    let k_hi = k * LN2_HI; // exact, see LN2_HI's comment
    fma(p, s, k_hi) + k * LN2_LO
}

/// log10(x), same Cody-Waite-combine approach as ln (see ln's own doc
/// comment for why this avoids the naive `log_2(x) * LOG10_2`'s double
/// rounding).
#[inline(always)]
pub fn log10(x: f32) -> f32 {
    let tiny = x < f32::MIN_POSITIVE;
    let xs = if tiny { x * 16777216.0 } else { x };
    let koff = if tiny { -24.0 } else { 0.0 };
    let r = log10_normal(xs, koff);
    let spec = if x == 0.0 { f32::NEG_INFINITY } else { f32::NAN };
    let r = if x <= 0.0 { spec } else { r };
    if !(x < f32::INFINITY) {
        x * x
    } else {
        r
    }
}

/// Core of log10 for positive normal finite x only -- see ln_normal, same
/// approach with coefficients fitted for log10 (log_2's c[i] * LOG10_2).
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn log10_normal(x: f32, koff: f32) -> f32 {
    let e = (x.to_bits() as i32).wrapping_sub(0x3f3504f3) >> 23;
    let m = f32::from_bits((x.to_bits() as i32).wrapping_sub(e << 23) as u32);
    let k = e as f32 + koff;
    let s = m - 1.0;
    let c: [f32; 10] = [
        0.4342945,
        -0.2171472,
        0.14476489,
        -0.10858066,
        0.08686763,
        -0.07213202,
        0.06159092,
        -0.05751561,
        0.05604425,
        -0.03309811,
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
    let k_hi = k * LOG10_2_HI; // exact, see LN2_HI's comment (same trick)
    fma(p, s, k_hi) + k * LOG10_2_LO
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
///
/// At x = +-0.0, `ln(u) + corr` adds two exactly-zero values of opposite
/// sign (`ln(1.0)` is `+0.0`, but `corr` correctly carries x's sign
/// there), which IEEE754 always resolves to `+0.0` -- the same mechanism
/// as sinf_poly's own `-0.0` bug. Fixed with a trailing `if x == 0.0 { x }
/// else { normal }` select (compute the normal path unconditionally
/// first, `log_2`'s own established "no early returns" idiom so array
/// loops keep auto-vectorizing): log1p is odd and monotonic through the
/// origin, so for every nonzero x in its domain `normal`'s sign already
/// equals x's (this is exactly what the correction above exists to
/// guarantee even for tiny x, where `corr` alone carries the whole
/// answer), making this select a no-op everywhere except the singular
/// zero point. A branchless `copysign(x)` was tried first and also
/// worked, but cost real latency once inlined into asinh/acosh (see
/// readme.md); this select form measured no better for them, so the
/// simpler, more idiomatic form was kept instead.
#[inline(always)]
pub fn log1p(x: f32) -> f32 {
    let u = 1.0 + x;
    let c = x - (u - 1.0);
    let corr = c / u;
    let corr = if corr.is_finite() { corr } else { 0.0 };
    let normal = ln(u) + corr;
    if x == 0.0 { x } else { normal }
}

/// exp(x) via a proper Cody-Waite reduction instead of `exp2(x * LOG2_E)`.
/// The naive form rounds `x * LOG2_E` *once* before ever calling exp2 --
/// that rounding lands on the *argument*, and since exp2's derivative
/// scales with exp2 itself, a relative error `delta` in the argument
/// becomes roughly `delta * ln2` of *relative* error in the result,
/// growing with `|x|` (worst near exp's own domain ceiling, ~88.7) --
/// dozens of ulp, confirmed the dominant error source for exp and
/// everything built on it. Fixed the standard way: `k = round(x*log2e)`,
/// `r = x - k*ln2` done as an exact Cody-Waite reduction (`LN2_HI`/`LN2_LO`,
/// the same split `ln`/`log10` already use -- `k*LN2_HI` is exact for this
/// domain's k, and `x - k*LN2_HI` is exact by Sterbenz since `k*ln2` tracks
/// `x` closely), then a dedicated degree-5 minimax poly for `e^r` directly
/// on `[-ln2/2, ln2/2]` (lolremez estimated max error 7.6e-8, ~1.3 ulp),
/// scaled by `2^k`. Inherits exp2's unchecked domain (see exp2's doc
/// comment): only accurate while `x*log2(e)` stays in `[-126, 128)`, i.e.
/// roughly `x` in `[-87.3, 88.7)` -- outside that, the exponent
/// construction produces garbage rather than a clamped/overflowed value,
/// same as before. `expm1`, `sinh`, `cosh`, `sinh_throughput`, and
/// `cosh_throughput` call this directly and inherit both this fix and
/// that same domain limit (measured: max ulp 63->4 / 64->8, avg ulp
/// roughly halved, across all of them). `tanh` (via `expm1`) is only
/// indirectly affected; `powf`/`erf`/`erfc` route through `exp2_checked`
/// directly and don't call this function at all, so are unaffected
/// either way.
///
/// Scaling by `2^k` needs exp2_checked's own k1/k2 split (not exp2's
/// simpler single-multiply trick), even though this function is otherwise
/// unchecked: `round` (unlike `floor`) can push `k` one integer past
/// where a *single* exponent-field construction stays valid -- e.g.
/// `x=88.37628` gives `x*log2e=127.50002`, which floor (what exp2 itself
/// uses) keeps at `k=127` (safely representable) but round pushes to
/// `k=128`, an exponent field value reserved for inf/NaN that a single
/// bit-trick multiply can't represent at all (found via a real
/// `exp(88.37628)=inf` edgecheck-style failure, not reasoned about in
/// advance) -- splitting into two representable halves sidesteps this by
/// construction, the same way exp2_checked already does for its own,
/// wider checked range.
#[inline(always)]
pub fn exp(x: f32) -> f32 {
    let k = fma(x, LOG2_E, 0.0).round();
    let r = fma(-k, LN2_HI, x);
    let r = fma(-k, LN2_LO, r);
    // c0 and c1 both forced exactly 1.0 (were 1.0 and 1.0000000647031426):
    // exp(r) = 1 + r + r^2*P(r) for tiny r, so a c1 off from 1.0 by even
    // ~6e-8 relative is a systematic bias right where exp(x) is most
    // commonly called (x near 0). l0 = r + 1.0 needs no fma since both its
    // coefficients are now exactly 1.0. c2..c5 coordinate-descent refit
    // for this constraint (examples/tune.rs's exp_r_c/"exp_r").
    let c: [f32; 4] = [4.9999008e-1, 1.6666375e-1, 4.1917525e-2, 8.3811125e-3];
    let r2 = r * r;
    let r4 = r2 * r2;
    let l0 = r + 1.0;
    let l1 = fma(c[1], r, c[0]);
    let l2 = fma(c[3], r, c[2]);
    let r0 = fma(l1, r2, l0);
    let p = fma(l2, r4, r0);
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k1b = fma(k, 0.5, ROUND_MAGIC) - (ROUND_MAGIC - 383.0);
    let k2b = (k + 766.0) - k1b;
    let t1 = f32::from_bits((k1b.to_bits() << 8) & EXPONENT_MASK);
    let t2 = f32::from_bits((k2b.to_bits() << 8) & EXPONENT_MASK);
    p * t1 * t2
}

/// A Pade approximant near 0 (where exp(x)-1 loses precision to
/// cancellation), exp(x)-1 directly elsewhere. See exp's doc comment for
/// the inherited unchecked-exp2 domain limit. The 5 coefficients (an
/// exact closed-form Pade identity to e^x, not an empirical fit -- small
/// integers, -120 shared between numerator and denominator) were refit
/// (2026-07-07) with examples/tune.rs against f64::exp_m1 over |x| < 0.5,
/// treating them as 5 independent free parameters. tune.rs's own
/// coarse tuning grid (~1.7M points) underestimated both the before and
/// after max ulp (reported 3->2); a targeted exhaustive scalar check of
/// just this branch (302M samples, step-7 over the full |x|<0.5 range)
/// found the real numbers are max ulp 4 -> 3, avg ulp 0.111 -> 0.109 --
/// still a genuine improvement, just smaller than the coarse grid
/// suggested. Zero perf cost (same instructions, only the 5 literal
/// constants differ).
#[inline(always)]
pub fn expm1(x: f32) -> f32 {
    let a = x * fma(-1.9999927, x * x, -120.0) / fma(x, fma(x, x - 12.000030, 59.999996), -120.0);
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

/// asinh(x) = ln(x + sqrt(x^2+1)), fixed for two bugs in the straight-ported
/// form (the doc comment used to describe only the first; the exhaustive
/// sweep that found it also turned up the second, worse one):
///
/// 1. Small-x cliff: for |x| below ~6e-8, `x*x` is already too small to
///    survive being summed into `x^2+1` (it rounds to exactly 1.0 before
///    sqrt even runs), so `x + sqrt(x^2+1)` collapses to exactly 1.0 and
///    `ln(...)` returns exactly 0 instead of the correct tiny nonzero
///    value. Fixed the standard way: `sqrt(x^2+1) - 1 = x^2 / (sqrt(x^2+1)
///    + 1)` (rationalized, no cancellation -- `x^2` is computed as its own
///    multiply here, independent of the lossy `x^2+1` sum, so it keeps
///    full precision), then `asinh(x) = log1p(x + (sqrt(x^2+1) - 1))`.
/// 2. Large-negative-x cancellation (the actual worst case: the exhaustive
///    sweep's argmax was x ~ -3.4e38, not a small-x input): for x very
///    negative, `sqrt(x^2+1) ~ |x|`, so `x + sqrt(x^2+1)` nearly cancels
///    to a small value whose *relative* precision is only as good as
///    `sqrt(x^2+1)`'s absolute error -- easily 20%+ off for large |x|.
///    Sidestepped entirely by computing on `|x|` (asinh is odd, so
///    `asinh(x) = sign(x) * asinh(|x|)`, restored with `mulsign`), where
///    `ax + sqrt(ax^2+1)` never cancels (both terms are non-negative).
/// Also mirrors acosh's overflow guards: `ax^2` overflowing prematurely
/// for huge `ax` (rescaled sqrt above `ax = 2048`, matching acosh's
/// threshold) and the final sum overflowing near f32::MAX (`ln(ax) +
/// LN_2` fallback, same asymptote as acosh's).
#[inline(always)]
pub fn asinh(x: f32) -> f32 {
    let ax = x.abs();
    let small = ax < 2048.0;
    // ax*ax used to be computed 3 times independently (once fused inside
    // direct_sq's fma, once for inv_ax2, once again for sm1's small-branch
    // numerator) -- shared here, same "both branches always computed
    // unconditionally, so shareable" reasoning as erf's own x2 fix.
    let ax2 = ax * ax;
    let direct_sq = (ax2 + 1.0).sqrt();
    let inv_ax2 = 1.0 / ax2;
    let rescaled_sq = ax * (1.0 + inv_ax2).sqrt();
    let sq = if small { direct_sq } else { rescaled_sq };
    let sm1 = if small { ax2 / (sq + 1.0) } else { sq - 1.0 };
    let d = ax + sm1;
    let r = if d.is_finite() { log1p(d) } else { ln(ax) + LN_2 };
    mulsign(r, x)
}

/// ln(x + sqrt(x^2-1)), domain x >= 1 (NaN elsewhere). Two bugs in the
/// straight-ported `x*x - 1.0` form, both from the same root cause
/// (`x*x` losing information well before it looks "wrong"):
///
/// 1. Sign loss: squaring erases x's sign, so sqrt(x^2-1) is the same
///    magnitude for +x and -x. Once |x| is large enough that ulp(x^2)
///    exceeds 1 (roughly |x| > 4096, far short of actual overflow) the
///    "-1" term vanishes entirely and sqrt(x^2-1) rounds to exactly |x|,
///    so `x + sqrt(x^2-1)` collapses to ~0 for negative x instead of
///    staying reliably negative -- ln of that silently returns finite
///    garbage (or +inf, once x^2 overflows) instead of the correct NaN.
///    Not a narrow edge case: wrong for roughly the whole range x < -4096.
///    Fixed with an explicit domain select (cheap next to the sqrt+ln
///    chain -- unlike the small-x cancellation fixes elsewhere in this
///    file, this one trades no accuracy or speed, it was just a missing
///    check).
/// 2. Premature overflow: for valid x above sqrt(f32::MAX) (~1.84e19),
///    `x*x` overflows to +inf even though the true answer (~ln(2x), at
///    most ~89.6 for any finite f32) stays comfortably finite -- ln(inf)
///    then wrongly returns +inf. Fixed by rescaling before squaring for
///    large x: sqrt(x^2-1) = x*sqrt(1 - 1/x^2); `1/x^2` underflows
///    gracefully to 0 for huge x (giving the correct sqrt(1-0)=1
///    asymptote) instead of `x^2` overflowing. This rescaled form costs
///    an extra rounding (the division) that the direct `fma(x,x,-1.0)`
///    doesn't pay, so it's only used above `x = 2048` -- comfortably
///    below where the direct form starts losing the "-1" term (~4096, see
///    point 1) but far enough into "smooth, ~ln(2x)" territory that the
///    switchover itself isn't a precision cliff; below that, the exact
///    single-rounding `fma(x,x,-1.0)` form stays in use, most importantly
///    right at the domain boundary x = 1 where acosh's derivative blows up
///    and every extra rounding gets amplified.
/// 3. A third, smaller-range overflow survives fix 2: once `x` itself is
///    within a factor of 2 of f32::MAX, `s` (now ~x exactly, per fix 2's
///    own asymptote) makes `x + s` ~2x overflow even though ln(2x) (~89)
///    is nowhere near overflowing. Guarded with `ln(x) + LN_2` (the same
///    asymptote, computed without ever forming 2x) whenever the sum isn't
///    finite.
/// 4. Right at the domain boundary x = 1 (where acosh's derivative blows
///    up, so any rounding gets amplified into a lot of ulps of a tiny
///    result), `ln(x + s)` computes `ln(1 + tiny)` -- exactly log1p's own
///    reason to exist. `d = (x - 1.0) + s` is `x + s - 1` computed with
///    `x - 1.0` exact (Sterbenz, x near 1) instead of forming `x + s`
///    (rounding tiny `s` against `x`'s magnitude) and subtracting 1 from
///    that afterward; `log1p(d)` then reuses log1p's own Sterbenz
///    correction on top. Cut max ulp right at the boundary from 1522 to 4
///    (exhaustive sweep), matching the small residual log1p/tanh/sinh's
///    own fixes elsewhere in this file leave behind.
#[inline(always)]
pub fn acosh(x: f32) -> f32 {
    // NOT the same "shared x2" opportunity asinh's own fix below found:
    // `direct` needs x*x - 1.0 computed as a *single* rounding (the fma)
    // specifically because x is near 1 at acosh's own domain boundary,
    // where x*x - 1.0 is a catastrophic-cancellation subtraction -- a
    // rounded-then-reused x2 loses exactly the precision that cancellation
    // needs (tried, measured: max ulp 3 -> 700). `direct` and `inv_x2`
    // each need their own x*x in a different rounding context, so there's
    // no real redundant computation to remove here after all.
    let direct = fma(x, x, -1.0).sqrt();
    let inv_x2 = 1.0 / (x * x);
    let rescaled = x * fma(-inv_x2, 1.0, 1.0).sqrt();
    let s = if x < 2048.0 { direct } else { rescaled };
    let d = (x - 1.0) + s;
    let r = if d.is_finite() { log1p(d) } else { ln(x) + LN_2 };
    if x < 1.0 { f32::NAN } else { r }
}

/// atanh(x) = 0.5*ln((1+x)/(1-x)) = 0.5*(log1p(x) - log1p(-x)), reusing
/// log1p's already-correct small-x handling instead of forming
/// (1+x)/(1-x) directly, which rounds to exactly 1.0 for tiny |x| (so
/// ln(...) came out exactly 0 instead of the correct tiny nonzero answer)
/// -- same fix shape as tanh reusing expm1 above. Domain x in [-1, 1]
/// falls out for free: log1p(-x) is -inf at x=1 and log1p(x) is -inf at
/// x=-1 (log1p's own domain edge), giving +-inf; |x|>1 makes one of the
/// two arguments < -1, where log1p is already NaN.
#[inline(always)]
pub fn atanh(x: f32) -> f32 {
    0.5 * (log1p(x) - log1p(-x))
}

// degree-6 minimax poly (Estrin via fma), fitted for acos's sqrt(1-|x|)
// factor. Ported from jodiemath's acosf_poly. Coefficients refit
// (2026-07-07) jointly against both acos's and asin's use of this same
// poly (asin reuses it directly, see asin's own doc comment fix 7 for
// the full story) -- asin's accuracy improved, acos's own left exactly
// unchanged by construction (the refit was constrained to guarantee it).
#[inline(always)]
fn acos_poly(x: f32) -> f32 {
    let u = 2.2960256e-3f32;
    let u = fma(u, x, -1.1146317e-2);
    let u = fma(u, x, 2.6900213e-2);
    let u = fma(u, x, -4.8802543e-2);
    let u = fma(u, x, 8.8755615e-2);
    let u = fma(u, x, -2.1458544e-1);
    fma(u, x, 1.5707963)
}

/// acos(x), domain x in [-1,1] (result always in [0,pi], never negative --
/// unlike sin/asinh/etc., acos isn't an odd function, so x=-0.0 has no
/// legitimate negative result the way it does for those). `mulsign`
/// (bit-based sign) and `x < 0.0` (value-based comparison) disagree on
/// exactly one input: `-0.0`, whose sign *bit* is set but whose *value*
/// equals `+0.0`. `mulsign` flipped `y`'s sign there (per the bit), while
/// `x < 0.0` correctly saw "not negative" and skipped the `+pi`
/// correction -- so `acos(-0.0)` came out `-pi/2` instead of the correct
/// `+pi/2`, an exhaustive-sweep-only find (fuzz sampling essentially
/// never lands on this one exact bit pattern; avg ulp exhaustively was
/// 0.99 with max ulp in the billions, not the quick-mode reading of
/// 0.50/4). Fixed by normalizing `-0.0` to `+0.0` before `mulsign` sees
/// it: `x + 0.0` is `-0.0 + 0.0 = +0.0` exactly (IEEE754's defined
/// round-to-nearest behavior for that one case) and a no-op for every
/// other `x`, genuinely negative or not.
#[inline(always)]
pub fn acos(x: f32) -> f32 {
    const PI: f32 = std::f32::consts::PI;
    let a = x.abs();
    let y = (1.0 - a).sqrt() * acos_poly(a);
    mulsign(y, x + 0.0) + if x < 0.0 { PI } else { 0.0 }
}

// asin(x) = x + x^3/6 + 3x^5/40 + 15x^7/336 + O(x^9), the odd Taylor series
// (exact rational coefficients). Unlike sinh's Taylor series, this one
// converges slowly as |x| approaches 1 (asin has a sqrt singularity
// there), so it's only used below |x| < 0.25 (see asin below), where 4
// terms already leave truncation error orders of magnitude under budget:
// the next (dropped) term, 105x^9/3456, is ~1.2e-7 at x=0.25 relative to
// asin(0.25) ~ 0.2527, i.e. ~4.6e-7 relative -- a few ulp, comparable to
// (not swamped by) the other branch's own residual there, see fix 6.
#[inline(always)]
fn asin_small(x: f32) -> f32 {
    let x2 = x * x;
    let c0 = 1.0f32;
    let c1 = 1.0 / 6.0f32;
    let c2 = 3.0 / 40.0f32;
    let c3 = 15.0 / 336.0f32;
    let p = fma(fma(fma(c3, x2, c2), x2, c1), x2, c0);
    x * p
}

/// A small-x Taylor branch plus the acos-style sqrt identity (`asin(x) =
/// pi/2 - acos(x)`, reusing acos's own poly). The straight-ported form
/// used a different, rational-correction formula instead of the acos
/// identity for everything above the small-x cutoff -- multiple accuracy
/// problems in that form, found and fixed one at a time as each fix's own
/// exhaustive re-sweep exposed the next (same pattern as acosh/asinh
/// above), until fix 6 replaced it outright:
/// 1. The final step, `sqrt(1-a) - 1`, has the same near-zero cancellation
///    as asinh/acosh's `sqrt(1+t) - 1` (loses precision as `a -> 0`, i.e.
///    as `x -> 0`, exactly where asin needs to be most accurate). Fixed
///    the same way as those: rationalized to `sqrt(1-a) - 1 = -a /
///    (sqrt(1-a) + 1)`, no cancellation. The rational correction
///    `a = (a^2-a)/d + a` is left exactly as fitted -- only the final
///    cancelling subtraction is rewritten.
/// 2. Fix 1 alone still left avg/max ulp at 293/825 (down from
///    169,116,394/852,038,214, but nowhere near budget): the rational
///    correction's own fit has a small persistent *relative* bias
///    (~3e-5, inherited from the C original's coarser ~1e-6-relative
///    goal) that cancellation had been masking -- with the cancellation
///    gone, that bias is what's left, concentrated at small x where it
///    dominates the (otherwise tiny) true answer. Fixed by adding a
///    genuine Taylor branch below |x| < 0.1 (see asin_small), the same
///    shape as sinh's fix -- no refit of the rational correction needed,
///    since it's only inaccurate in the regime the Taylor branch now
///    covers instead.
/// After fixes 1-2, avg ulp was 0.328 (in budget) but max ulp still 468,
/// concentrated right at the x=1 boundary: asin's derivative has a sqrt
/// singularity there, so the rational correction's own ~1e-6-relative fit
/// gets amplified into real ulp error, the same mechanism as the small-x
/// cliff just manifesting at the *other* extreme the fit wasn't tuned for.
/// 3. Closed with a third branch above `a > 0.9`, reusing acos_poly (no
///    new fit): `asin(x) = pi/2 - acos(x)`, and for a >= 0, `acos(a) =
///    sqrt(1-a) * acos_poly(a)` is exactly acos's own well-conditioned
///    formula (a shrinking sqrt factor times a smooth bounded poly, no
///    subtraction between comparable-magnitude values) -- computing
///    `pi/2 - acos(a)` doesn't cancel either, since acos(a) is small near
///    a=1 while pi/2 is O(1). Sign restored the same way as the other two
///    branches, via `mulsign`. Result: max ulp 468 -> 121 (avg ulp barely
///    moved, 0.328 -> 0.325, since it was already in budget). The
///    remaining worst case moved to the *other* branch boundary (x ~ 0.1,
///    the small/mid transition) -- same root cause (the mid branch's own
///    ~3e-5 relative bias), now the largest surviving residual since
///    that bias scales with x and can't be shrunk further by moving
///    thresholds around; a full fix would need the mid branch refit or
///    widened *and* extended (asin's Taylor series converges too slowly
///    approaching x=1 to just push the small-x cutoff out).
/// 4. The mid branch's 4 coefficients, refit (2026-07-07) with
///    examples/tune.rs's coordinate-descent tuner against a grid over
///    a in [0.1, 0.9) (exactly the range the mid branch is used for in
///    the shipped code). Zero perf cost -- same instructions, only the
///    literal constants changed. Result: max ulp 121 -> 84 (a real ~30%
///    cut); avg ulp moved 0.325 -> 0.377, still comfortably under budget.
///    The tuner's own lexicographic (max, then avg) objective explains
///    the tradeoff -- a separate avg-first, max-capped search found only
///    a negligible avg improvement (25.115 -> 25.092 on the tuning grid)
///    with max unchanged, confirming the two objectives pull in different
///    directions around this point rather than one dominating the other;
///    kept the max-first result since max was the specifically open
///    residual from fix 3. This step's own comment previously claimed the
///    x~0.1 boundary residual "can't be shrunk further by moving
///    thresholds around" -- untested at the time, and wrong: see fix 5.
/// 5. The small/mid threshold moved from 0.1 to 0.3 (2026-07-07, later
///    same day). Re-checked fix 4's own untested claim by actually
///    probing both branches independently instead of trusting it: the
///    mid branch's error isn't a narrow spike right at its own domain
///    edge, it's elevated (up to ~80 ulp) across the whole [0.1, 0.3)
///    band and only settles below ~40 past that -- while `asin_small`
///    (an exact-coefficient degree-7 Taylor series) stays at 1-2 ulp all
///    the way out to x<0.2 and is still fine at 0.3 (its own error only
///    climbs past there, crossing mid's around x~0.3-0.35). Since both
///    branches are already computed unconditionally in this branchless
///    select, moving the threshold is a pure accuracy change, zero perf
///    cost regardless of where it lands (confirmed: same instructions,
///    only the comparison constant differs). Result (exhaustive, all
///    2^32 f32 bit patterns): max ulp 84 -> 41, avg ulp 0.377 -> 0.105
///    (both improved together, not a tradeoff). Root cause of fix 4's
///    wrong claim: it reasoned from "the bias scales with x" without
///    separately measuring where `asin_small` itself stopped being
///    trustworthy -- the two curves were never actually compared until
///    now.
/// 6. The `mid` branch (the rational correction from fixes 1-2) removed
///    entirely (2026-07-07, later still). Once fix 5's investigation
///    showed `acos_poly`'s `sqrt(1-a)*acos_poly(a)` formula (the `near1`
///    branch) staying accurate from x~0.25 onward, it raised an obvious
///    question fix 3 never actually asked: `acos_poly` is `acos`'s *own*
///    poly, fit and used across `acos`'s *entire* `[0,1]` domain, not
///    something specifically tuned "for near 1" -- the `near1` name was
///    just a historical accident of where it was *first* reused here (a
///    > 0.9), not a real limitation. Measured how far it holds up in
///    isolation: garbage near x=0 (catastrophic cancellation in
///    `pi/2 - acos_poly(0)`, the same mechanism the `mid`/`small` split
///    exists to avoid), but from x >= 0.25 it's already better than `mid`
///    ever was anywhere in `mid`'s own former domain -- so `mid` wasn't
///    filling a gap `near1` couldn't cover, it was just never tried
///    there. A direct 2-branch (`asin_small` / `near1` only) coordinate
///    search over the crossover point found a flat minimum around
///    `a < 0.25`. Collapsing to 2 branches removes `mid`'s entire
///    computation (3 fma's, 2 divisions, a sqrt) from every call, not
///    just changes a threshold. Exhaustive sweep: max ulp 41 -> 11, avg
///    ulp 0.105 -> 0.033 (both improved again). mca is a genuine mixed
///    result, not the "faster on both axes" outcome the op-count cut
///    suggested before actually measuring: throughput improved a lot
///    (2.899 -> 0.968 cyc/elem, -67%, fewer total ops in the vectorized
///    loop), but latency got *worse* (43.24 -> 59.03 cyc, +37%) -- the
///    same non-monotonic-scheduling surprise logged elsewhere in this
///    file (asinh's own fix showed the identical shape): with 3 branches
///    computed unconditionally, the scalar chain had independent work to
///    fill cycles that would otherwise sit idle waiting on the sqrt;
///    with only 2, less of that slack exists. Kept anyway -- accuracy
///    improved substantially and throughput (the metric this crate's
///    vectorization-first design actually prioritizes, per the readme)
///    improved even more, both by a wide margin; only the secondary,
///    diagnostic latency number regressed, the same shape of tradeoff
///    already accepted for asinh/acosh's log1p-fix side effect earlier
///    this session.
/// 7. `acos_poly`'s coefficients refit (2026-07-07, immediate follow-up
///    to fix 6) against a *joint* objective instead of `acos`'s error
///    alone: `acos(x)` shrinks to 0 as `x -> 1` while `asin(x)` grows to
///    `pi/2` there, so the same absolute poly error carries a different
///    *relative* (ulp) weight depending on which caller's output it's
///    measured against -- a poly tuned purely for acos's own hardest
///    region may leave accuracy on the table specifically for asin's use
///    of the same coefficients. `examples/tune.rs`'s `acos_poly_c` tuner
///    extended with a constrained search: minimize asin's error subject
///    to acos's own on-grid max ulp never exceeding its already-tuned
///    best (a stricter, cross-function-safe version of a first
///    unconstrained joint-max attempt, which found a similar asin
///    improvement but let acos's own max regress by 1 ulp -- rejected in
///    favor of this one once the safe version was confirmed to find
///    asin gains too). Exhaustive sweep: asin max ulp 11 -> 9, avg 0.033
///    -> 0.030; acos itself exactly unchanged (max ulp 4, avg 0.496,
///    bit-for-bit identical to its pre-refit values). Zero perf cost for
///    either function (same instructions, only the 7 literal constants
///    differ).
#[inline(always)]
pub fn asin(x: f32) -> f32 {
    let a = x.abs();
    let small = asin_small(x);
    let big = mulsign(FRAC_PI_2 - (1.0 - a).sqrt() * acos_poly(a), x);
    if a < 0.25 { small } else { big }
}

// Pade-style rational approximation of atan on [0,1]. Ported from
// jodiemath's atanf_poly; coefficients refit (2026-07-07) with
// examples/tune.rs's coordinate-descent tuner against f64::atan over
// [0,1] (exactly atan's own domain for this poly, via a.min(1/a)). Zero
// perf cost (same instructions, only the 4 literal constants changed):
// max ulp 19 -> 18 and avg ulp improved too for both atan and atan2
// (which calls atan directly) -- unlike asin's mid-branch refit, this one
// didn't trade one metric for the other, both moved the same direction.
// A denser tuning grid (39M points vs. this file's 25k) converged to the
// same coefficients, suggesting this is a genuine local optimum for this
// coordinate-descent scheme, not an artifact of grid resolution.
#[inline(always)]
fn atan_poly(x: f32) -> f32 {
    let a = f32::from_bits(0x3d26709d);
    let b = f32::from_bits(0x3f285132);
    let c = f32::from_bits(0x3e2f7213);
    let d = f32::from_bits(0x3f7da42d);
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

/// atan2(y, x). `atan2(-0.0, +0.0)` used to come out `+0.0` instead of
/// IEEE754/C99's defined `-0.0`: when `x` is exactly `+0.0`, `base`
/// degenerates to exactly `+0.0`, and the final `base + mulsign(...)`
/// combines it with the (correctly `-0.0`-signed) correction term --
/// but IEEE754 addition of two *opposite*-signed zeros is defined to
/// give `+0.0` regardless of operand order (only same-signed zeros, or a
/// zero plus a genuine nonzero value, preserve the expected sign), so
/// the correction's sign silently vanished. Every other zero/sign
/// combination avoids this: `x = -0.0` makes `hpisignx` flip sign too
/// (so the correction becomes a real nonzero `+-PI`, not a degenerate
/// zero), and whenever `x` is genuinely nonzero, `base` is a real
/// nonzero-ish angle, not an exact zero, so the addition never hits the
/// opposite-sign-zero case. Fixed by skipping the addition entirely when
/// `base` would be that exact `+0.0` -- `nonzerox` already selects
/// between the two shapes, so no new branch, just moved.
#[inline(always)]
pub fn atan2(y: f32, x: f32) -> f32 {
    let nonzerox = x != 0.0;
    let nonzeroy = y != 0.0;
    let bothzero = !nonzerox && !nonzeroy;
    let hpisignx = if nonzerox || bothzero { mulsign(FRAC_PI_2, x) } else { 0.0 };
    let correction = mulsign(FRAC_PI_2 - hpisignx, y);
    let r = if nonzerox { atan(y / x) + correction } else { correction };
    // atan2(+-inf, +-inf): y/x is inf/inf, which is NaN, so the general
    // formula above can't produce an answer here at all. IEEE754/C99
    // define a canonical result by quadrant regardless (+-pi/4 or
    // +-3pi/4) -- not derived from any real ratio, since there isn't one
    // at true infinity, just a fixed convention.
    let bothinf = x.is_infinite() && y.is_infinite();
    let inf_result = mulsign(if x.is_sign_negative() { 3.0 * FRAC_PI_4 } else { FRAC_PI_4 }, y);
    if bothinf { inf_result } else { r }
}

/// Straight port of jodiemath's tanf: sin(x)/cos(x), same domain limits as
/// this crate's sin/cos (see their doc comments).
#[inline(always)]
pub fn tan(x: f32) -> f32 {
    sin(x) / cos(x)
}

// degree-6 minimax poly feeding erf's exp2-based tail (|x| >= 0.28). Ported
// from jodiemath's erff_poly. Regrouped from a 6-deep Horner chain to
// Estrin (2026-07-07, same restructuring tried on acos_poly immediately
// before this) -- 3 fma's deep instead of 6, same 6 fma's total plus 2
// extra plain multiplies for x^2/x^4, same coefficients unchanged. Unlike
// acos_poly's attempt, this one turned out to be accuracy-neutral on the
// exhaustive sweep (avg/max ulp exactly unchanged, 0.319/5) -- fma
// reassociation doesn't always cost accuracy, it has to be checked per
// poly, not assumed either way. mca: 102.74->91.74 cyc latency (-10.7%),
// 3.163->2.871 cyc/elem throughput (-9.2%) -- both axes improved
// together here, unlike acos_poly's case (small throughput cost there).
#[inline(always)]
fn erf_poly(x: f32, x2: f32) -> f32 {
    let a6 = 3.118769e-4f32;
    let a5 = -4.67225e-3f32;
    let a4 = 3.3162573e-2f32;
    let a3 = -1.5214339e-1f32;
    let a2 = -9.1684705e-1f32;
    let a1 = -1.6282598f32;
    let a0 = 3.1332566e-5f32;
    let x4 = x2 * x2;
    let b0 = fma(a1, x, a0);
    let b1 = fma(a3, x, a2);
    let b2 = fma(a5, x, a4);
    let c0 = fma(b1, x2, b0);
    let c1 = fma(a6, x2, b2);
    fma(c1, x4, c0)
}

/// A Pade approximant near 0 (where the tail form loses precision to
/// cancellation), the exp2-based tail elsewhere. The tail branch used to
/// evaluate `erf_poly` directly on `|x|` with no bound, on the theory
/// (stated here previously) that it's "only relevant once erf has long
/// since saturated to +-1... so it doesn't affect any input where the
/// answer isn't already indistinguishable from +-1" -- checked by hand
/// and false: `erf_poly` is a plain degree-6 polynomial (unbounded, not
/// fitted to stay well-behaved outside where it was tuned), and its
/// leading term's positive coefficient means it eventually turns around
/// and grows to +inf for large |x| instead of staying deeply negative
/// (erf_poly(9) ~ -92, erf_poly(20) ~ +8698) -- so `exp2(erf_poly(|x|))`
/// stopped being ~0 and started being huge, giving `erf(50) = -1.02e17`,
/// `erf(100) = NaN`, `erf(+-inf) = NaN` instead of the correct +-1.
/// Fixed by clamping `|x|` to 10 before `erf_poly` (same bound erfc's own
/// clamp already uses, chosen the same way: comfortably past where erf
/// has actually saturated -- confirmed erf_poly stays safely negative,
/// -83.8, at that bound) -- `exp2` alone wasn't the bug here (its
/// [-126,128) domain covers erf_poly(10) fine); swapped to `exp2_checked`
/// anyway as cheap extra insurance now that the input is guaranteed
/// bounded, matching erfc's fix.
#[inline(always)]
pub fn erf(x: f32) -> f32 {
    let xa = x.abs();
    let xa_bounded = if xa > 10.0 { 10.0 } else { xa };
    // Shared between both branches: the Pade arm's own x2 (x*x, unclamped)
    // and erf_poly's internal x2 (xa_bounded*xa_bounded) only differ once
    // |x| > 10, but the Pade arm's result (`a`) is only ever *selected*
    // when `xa < 0.28`, comfortably inside the clamp -- so reusing the
    // already-clamped x2 here changes nothing observable, just removes a
    // redundant multiply (and bounds the discarded arm's x2 to <= 100
    // instead of letting it run up toward overflow for huge |x|, a minor
    // side benefit, not the point of the change).
    let x2 = xa_bounded * xa_bounded;
    let numer = x * fma(f32::from_bits(0x3f174f6e), x2, f32::from_bits(0x3f906ebb));
    let denom = fma(fma(f32::from_bits(0x3e3e2be3), x2, f32::from_bits(0x3f5b6db7)), x2, 1.0);
    let a = numer / denom;
    let b = mulsign(1.0 - exp2_checked(erf_poly(xa_bounded, x2)), x);
    if xa < 0.28 { a } else { b }
}

/// A rational*gaussian tail, clamped to |x| <= 10 before evaluation
/// (matching the C original) -- the polynomial (n/d, degree 4 in xa)
/// still needs that bound to avoid overflowing for huge x, but the
/// Gaussian factor no longer does: it used to be `exp(-xa*xa)`, which
/// routes through the *unchecked* exp2 ([-126, 128) domain), and for
/// |x| >= ~9.35 (xa*xa >= ~87.4) the exponent `-xa*xa*log2(e)` is already
/// past -126, so the result was unreliable garbage for that whole tail
/// instead of a clean 0. Fixed by calling `exp2_checked` directly on the
/// same exponent instead of going through `exp` -- its wider [-151, 128)
/// domain comfortably covers the full clamped range (xa in [0,10] means
/// the exponent never goes below -100*log2(e) =~ -144.3, still inside
/// exp2_checked's bound), and it's already the crate's existing
/// correctly-rounded full-range primitive, no new code needed.
/// n/d's 8 coefficients refit (2026-07-07) with examples/tune.rs against
/// sleef's erfc reference over the full clamped domain. Zero perf cost
/// (same instructions): max ulp 123 -> 109 (~11%); avg ulp moved slightly
/// worse (0.297 -> 0.311, still comfortably under budget), same kind of
/// max/avg tradeoff as asin's mid-branch refit.
#[inline(always)]
pub fn erfc(x: f32) -> f32 {
    let z = if x < 0.0 { -1.0 } else { 1.0 };
    // w = 1.0 - z exactly, for both branches (0 = 1-1, 2 = 1-(-1)) --
    // one subtract instead of a second compare+select on the same
    // condition z already resolved.
    let w = 1.0 - z;
    // NaN-preserving clamp: f32::min suppresses NaN (returns the other
    // operand), unlike C's `x>10.f?10.f:x` ternary (false for NaN, so it
    // takes the x branch, keeping NaN). This if/else matches the ternary.
    let xa = x.abs();
    let xa = if xa > 10.0 { 10.0 } else { xa };
    let n = fma(f32::from_bits(0x35c42f59), xa, f32::from_bits(0x3daf42cd));
    let n = fma(n, xa, f32::from_bits(0x3ee32e39));
    let n = fma(n, xa, f32::from_bits(0x3f7a7525));
    let n = fma(n, xa, 1.0);
    let d = fma(f32::from_bits(0x3e1b69eb), xa, f32::from_bits(0x3f48fdde));
    let d = fma(d, xa, f32::from_bits(0x3fe918da));
    let d = fma(d, xa, f32::from_bits(0x4006d464));
    let d = fma(d, xa, 1.0);
    let y = exp2_checked(-(xa * xa) * LOG2_E) * n / d;
    fma(y, z, w)
}

/// Straight port of jodiemath's hypotf: naive sqrt(x^2+y^2), no anti-overflow
/// rescaling (unlike std's hypot) -- trades the overflow/underflow edge cases
/// for vectorizability, same tradeoff this crate makes for cbrt/sin/cos vs.
/// their std counterparts.
#[inline(always)]
pub fn hypot(x: f32, y: f32) -> f32 {
    let normal = fma(x, x, y * y).sqrt();
    // hypot(+-inf, anything) and hypot(anything, +-inf) = +inf, even when
    // the other argument is NaN -- IEEE754/C99 special-cases infinity to
    // "win" over NaN here (unlike almost every other function). The naive
    // formula can't reach this on its own: once either argument actually
    // is NaN, `inf*inf + NaN*NaN` degrades to NaN instead. Distinct from
    // this function's already-documented overflow tradeoff above (that's
    // about *finite* x/y large enough to overflow x*x/y*y; this is about
    // a literally-infinite argument, unaffected by that tradeoff either way).
    if x.is_infinite() || y.is_infinite() { f32::INFINITY } else { normal }
}

/// log2(x) as a double-float (Df32) instead of a collapsed f32, for
/// positive finite x only (same domain log_2_normal assumes -- callers
/// must guard zero/negative/inf/nan themselves). Reuses log_2_normal's
/// exact decomposition and poly (`s = m - 1`, `P(s)`) but keeps the
/// integer exponent `k` and the poly correction `s*P(s)` apart via
/// `Df32::from_add` (a proper two-sum, not a naive pair) instead of
/// collapsing them with a single `fma(p, s, k)` -- `s*P(s)` alone already
/// has full relative f32 precision, but adding it to `k` in a single
/// rounding (as log_2_normal does) throws away exactly the low bits a
/// later multiply-by-y would otherwise be able to use. Denormal input is
/// handled the same way log_2's own wrapper does (scale up, offset k).
#[inline(always)]
fn log2_df(x: f32) -> Df32 {
    let tiny = x < f32::MIN_POSITIVE;
    let xs = if tiny { x * 16777216.0 } else { x };
    let koff = if tiny { -24.0 } else { 0.0 };
    let e = (xs.to_bits() as i32).wrapping_sub(0x3f3504f3) >> 23;
    let m = f32::from_bits((xs.to_bits() as i32).wrapping_sub(e << 23) as u32);
    let k = e as f32 + koff;
    let s = m - 1.0;
    let c: [f32; 10] = [
        LOG2_E,
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
    Df32::from_add(k, p * s)
}

/// exp2 of a double-float argument, reusing exp2_checked's own clamp/
/// k1-k2-split/poly machinery on the hi component, then folding the lo
/// component in as a multiplicative correction:
/// `exp2(hi + lo) = exp2(hi) * exp2(lo) ~= exp2(hi) * (1 + lo*ln2)`
/// (first-order Taylor, valid since `lo` is always tiny relative to 1 by
/// Df32's own invariant) -- `fma(result, lo*LN_2, result)`. This is the
/// step that actually captures the precision `log2_df` preserved: without
/// it, the lo component would just be silently dropped and this would be
/// no more accurate than the plain `exp2_checked(v.to_f32())`.
#[inline(always)]
fn exp2_checked_df(v: Df32) -> f32 {
    // clamp *before* floor/subtract (matching exp2_checked exactly): if v.0
    // itself is already +-inf (y large enough that y*log2(x) overflows in
    // the Df32 multiply), flooring first and clamping after would compute
    // `inf - floor(inf)` = `inf - inf` = NaN instead of correctly
    // saturating.
    let xs = v.0.clamp(-151.0, 128.0);
    let k = xs.floor();
    let f = xs - k;
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k1b = fma(k, 0.5, ROUND_MAGIC) - (ROUND_MAGIC - 383.0);
    let k2b = (k + 766.0) - k1b;
    let t1 = f32::from_bits((k1b.to_bits() << 8) & EXPONENT_MASK);
    let t2 = f32::from_bits((k2b.to_bits() << 8) & EXPONENT_MASK);
    let f2 = f * f;
    let g0 = fma(2.4022985e-1, f, 6.93147e-1);
    let g1 = fma(9.678826e-3, f, 5.548333e-2);
    let g2 = fma(2.1702237e-4, f, 1.2439679e-3);
    let h = fma(g2, f2, g1);
    let q = fma(h, f2, g0);
    let p = fma(q, t1 * f, t1);
    let result = p * t2;
    // when result saturates to 0 or +-inf (xs clamped away from its real
    // value), `fma(result, v.1*LN_2, result)` can hit an inf*0 or 0*finite
    // -> still-fine-looking-but-actually-NaN indeterminate form (e.g.
    // result=inf, v.1=0.0: inf*0.0 is NaN by itself, and NaN+inf is NaN,
    // contaminating an otherwise-correct saturated result) -- same
    // "correction term goes non-finite" class of bug as log1p/atanh/
    // acosh/atan2's own corr.is_finite() guards. The correction is
    // meaningless in the saturated regime anyway (there's no precision
    // left to refine), so skip it there.
    let corr = fma(result, v.1 * LN_2, result);
    if corr.is_finite() { corr } else { result }
}

/// exp2(log2(x) * y). Used to route through the *unchecked* exp2 for its
/// exponent -- correctly documented as inaccurate outside
/// `log2(x)*y in [-126, 128)`, but "inaccurate" undersold it: outside
/// that range the unchecked bit-trick construction wraps around instead
/// of overflowing/underflowing, so e.g. `powf(2.0, 500.0)` (should be
/// `inf`) came out `2.88e17`, `powf(2.0, 1000.0)` (should be `inf`) came
/// out `3.6e-12` -- plausible-looking finite garbage, not just reduced
/// precision. Also a real tier mismatch: `log_2` is already the
/// full-range *checked* primitive (handles zero/negative/denormal/inf
/// cleanly), so powf was already paying that cost without getting the
/// matching benefit on the exp2 side. Fixed by routing through
/// `exp2_checked` instead, so the whole function is consistently
/// checked. Real perf cost (unlike erf/erfc's fixes, which only touched
/// a rarely-hit edge branch, this touches every call): see mca numbers
/// in the readme/IDEAS.md. See [`powf_checked`] for a variant with
/// substantially better accuracy for large `|y|`, at extra cost.
#[inline(always)]
pub fn powf(x: f32, y: f32) -> f32 {
    let ax = x.abs();
    let mag = exp2_checked(log_2(ax) * y);
    // For negative x, `exp2(log2(|x|)*y)` alone can't ever be negative
    // (exp2 of any real argument is positive), so this route always gave
    // NaN for x < 0.0 -- even for a well-defined case like `(-2.0)^3.0 =
    // -8.0`. A real result only exists there when y is an integer: even y
    // -> +mag, odd y -> -mag (reusing `parity`, the same integer-parity
    // helper sin_checked/cos_checked already use), non-integer y -> NaN
    // (correctly matches std, e.g. `(-8.0)^(1/3)` is NaN in f32 too --
    // real cube roots of negative numbers aren't picked by this branch).
    let y_int = y == y.trunc();
    let y_odd = y_int && parity(y) != 0.0;
    let neg_signed = if y_odd { -mag } else { mag };
    let neg_result = if y_int { neg_signed } else { f32::NAN };
    // `x.is_sign_negative()` (bit-based), not `x < 0.0` (value-based): the
    // latter disagrees with the former exactly at x = -0.0 (same class of
    // bug as acos's earlier `-0.0` fix this session), which would silently
    // route `(-0.0)^3.0` through the wrong (positive) branch instead of
    // the correctly-signed `-0.0`.
    let r = if x.is_sign_negative() { neg_result } else { mag };
    // pow(x, 0) = 1 for *any* x -- even 0, negative, or NaN -- a
    // dedicated IEEE754/C99 special case, not derivable from the log/exp2
    // formula (0*inf and NaN*0 both degrade to NaN above). Override last.
    if y == 0.0 { 1.0 } else { r }
}

/// Higher-accuracy variant of [`powf`]: `exp2(log_2(x)*y)` amplifies
/// log_2's own rounding error by `y` -- for `|y|` large that swamps the
/// result (hundreds of ulp), since `log_2(x)` is collapsed to a single f32
/// *before* the multiply, throwing away exactly the low bits that `y`'s
/// multiplication would otherwise be able to use. Fixed by keeping
/// `log2(x)` as a double-float (Df32) through the multiply by `y` and the
/// exp2 reconstruction, only collapsing to a single f32 at the very end
/// (see [`log2_df`]/[`exp2_checked_df`]). Confirmed by fuzzing (25M
/// samples, apples-to-apples against the plain formula on the same
/// inputs): avg ulp 0.36 -> 0.05 (-87%), max ulp 170 -> 154 -- both
/// improve, though the remaining ~150 max ulp is a separate, pre-existing
/// `exp2_checked` characteristic near the denormal/underflow boundary
/// (present, and about equally large, in *both* the plain and precise
/// formula there), not something this fix introduces or was meant to
/// close. Real mca cost though (double-float bookkeeping plus the poly
/// evaluation isn't free): latency 98.03->105.48 cyc (+7.6%), throughput
/// 3.898->6.285 cyc/elem (+61.2%) -- kept as an opt-in tier rather than
/// the default, matching sin/sin_checked and exp2/exp2_checked.
#[inline(always)]
pub fn powf_checked(x: f32, y: f32) -> f32 {
    let ax = x.abs();
    // The df path is only valid for ax strictly positive and finite (same
    // domain log_2_normal's own decomposition assumes); ax == 0 or
    // non-finite needs a fallback, but *not* a second full
    // log_2+exp2_checked computation -- that would roughly double this
    // function's cost just to cover a few degenerate inputs. The only
    // three possible magnitudes there are cheap direct selects: ax == 0
    // -> 0 (y > 0) or +inf (y < 0); ax == +inf -> +inf (y > 0) or 0
    // (y < 0); ax == NaN (from x == NaN) -> NaN. (y == 0 is overridden
    // separately below regardless of any of this.)
    let is_safe = ax > 0.0 && ax.is_finite();
    let mag_precise = exp2_checked_df(log2_df(ax) * y);
    // (ax == 0) == (y > 0) picks out exactly the two "goes to zero" cases
    // (ax==0,y>0 and ax==+inf,y<0) vs. the two "goes to infinity" cases --
    // cheaper than a 4-way branch and avoids inf/inf-is-NaN traps a
    // division-based shortcut would hit for the ax==+inf,y<0 case.
    let zero_or_inf = if (ax == 0.0) == (y > 0.0) { 0.0 } else { f32::INFINITY };
    let edge_mag = if ax.is_nan() { f32::NAN } else { zero_or_inf };
    let mag = if is_safe { mag_precise } else { edge_mag };
    // Same negative-x/y-parity/y==0 handling as powf -- see its own doc
    // comment for the reasoning.
    let y_int = y == y.trunc();
    let y_odd = y_int && parity(y) != 0.0;
    let neg_signed = if y_odd { -mag } else { mag };
    let neg_result = if y_int { neg_signed } else { f32::NAN };
    let r = if x.is_sign_negative() { neg_result } else { mag };
    if y == 0.0 { 1.0 } else { r }
}

/// Straight port of jodiemath's remainderf: x - round(x/y)*y (ties away from
/// zero, via f32::round -- not IEEE 754 remainder's ties-to-even).
/// round(x/y)*y's absolute error scales with ulp(x), which swamps the true
/// remainder (at most |y|/2) once |x/y| is large -- inherited from the C
/// original's identical formula, only reliable while |x/y| stays moderate.
///
/// Even within that "moderate" range this formula has a real, low-probability
/// failure mode: `x/y`'s own single-rounding division error can occasionally
/// land `q = (x/y).round()` on the wrong side of a true half-integer tie
/// (whenever the exact mathematical `x/y` happens to fall within about half a
/// division-ulp of `N+0.5`), producing a result with the *wrong sign* and a
/// similar magnitude to the correct answer -- confirmed by fuzzing (found
/// concrete cases with `|x/y|` as small as ~100). The failure probability
/// scales with `|x/y|` (roughly `|x/y| * 2^-24`) but isn't zero at any ratio.
/// See [`remainder_checked`] for a variant that detects and self-corrects
/// this at extra cost.
///
/// At x = +-0.0 (nonzero finite y), `fma(-q, y, x)` adds two exactly-zero
/// values of opposite sign (`q` is `+-0.0` matching x/y's sign, so `-q*y`
/// ends up the *opposite* sign to x), which IEEE754 always resolves to
/// `+0.0` -- same mechanism as sinf_poly's own `-0.0` bug, but unlike
/// sinf_poly/log1p, remainder's result sign does *not* generally track
/// x's sign for nonzero x (e.g. remainder(2.0, 3.0) == -1.0, an IEEE
/// remainder property, not a bug), so `copysign(x)` isn't a valid
/// blanket fix here -- guarded with an explicit `x == 0.0` select
/// instead, correct for both this singularity and the ordinary
/// remainder(+0.0, y) case (already correctly `+0.0`, so the guard is a
/// no-op there).
#[inline(always)]
pub fn remainder(x: f32, y: f32) -> f32 {
    let q = (x / y).round();
    let normal = fma(-q, y, x);
    let r = if x == 0.0 { x } else { normal };
    // remainder(finite x, +-inf) = x (IEEE754/C99 special case): q rounds
    // to exactly 0.0 for any finite x, but `fma(-q, y, x)` then multiplies
    // that zero by an *infinite* y, giving NaN (0*inf is NaN) instead of
    // the intended "no reduction happened, answer is just x" no-op. x
    // itself infinite/nan still correctly falls through to `normal`
    // (matches std's remainder(inf, ...) = NaN) since `x.is_finite()`
    // excludes it here.
    if y.is_infinite() && x.is_finite() { x } else { r }
}

/// Self-correcting variant of [`remainder`]: detects when `x/y`'s own
/// division rounding pushed `q` to the wrong side of a half-integer tie
/// (see [`remainder`]'s doc comment for the failure mode -- a rare but
/// real sign-flip bug, not just a large-ratio accuracy gap) and nudges `q`
/// by one integer to correct it. A correctly-rounded `q` always leaves
/// `|r0| <= |y|/2`, so that inequality failing is a direct signal to move
/// toward whichever side shrinks the residual and recompute -- both
/// branches are computed unconditionally and selected, matching this
/// crate's branchless style. Confirmed by fuzzing against an f64
/// reference: 0 max ulp for `|x/y|` up to `1e7` (including every concrete
/// sign-flip case [`remainder`] can hit), degrading only past `2^24`
/// where `q` itself stops being an exactly-representable f32 integer -- a
/// separate, harder limit this correction can't reach past. Costs a
/// second `fma` plus the correction's compare/select on every call (mca:
/// +42% latency, +60% throughput vs plain `remainder`), so kept as an
/// opt-in tier for callers who need the reliability guarantee, matching
/// sin/sin_checked and exp2/exp2_checked.
#[inline(always)]
pub fn remainder_checked(x: f32, y: f32) -> f32 {
    let q0 = (x / y).round();
    let r0 = fma(-q0, y, x);
    let adj = if (r0 > 0.0) == (y > 0.0) { 1.0 } else { -1.0 };
    let r1 = fma(-(q0 + adj), y, x);
    let normal = if r0.abs() > y.abs() * 0.5 { r1 } else { r0 };
    let r = if x == 0.0 { x } else { normal };
    if y.is_infinite() && x.is_finite() { x } else { r }
}

#[cfg(test)]
mod tests {
    use super::*;

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
