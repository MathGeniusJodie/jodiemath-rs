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
        1.442695,
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
// q = round(x/pi) must be an *exact* integer for x - q*pi to land accurately
// in [-pi/2, pi/2], and a single f32 only holds exact integers up to 2^24 --
// past that the old single-word reduction silently mis-rounded q, and for
// large enough x the polynomial input overflowed to literal inf (see
// readme's sin/cos accuracy note; this was a real, previously undocumented
// bug found by examples/accuracy.rs's exhaustive/fuzz sweep).
//
// Fix: split q into an exact double-float integer pair (qh, ql), and pi
// into 3 words (PI_W0..2, unrelated to the crate's other PI_A..D -- those
// use trailing-zero padding to make q*PI_A exact for *bounded* q, which
// breaks down once qh itself is large; these are genuine double-float
// words, each computed via `two_prod`/`two_sum` error-free transforms that
// stay exact for *any* magnitude). Not every cross term needs the same
// treatment: qh*PI_W0 (dominant) and the two next-biggest cross terms
// (qh*PI_W1, ql*PI_W0, both ~x*2^-24) need a real two-product to avoid
// losing precision, but the smaller ones (qh*PI_W2, ql*PI_W1, ~x*2^-48)
// only need to be *present* in the sum, not extra-precise -- a plain
// multiply's own rounding error at their magnitude is already far below
// 1 ulp of the final O(1) result. Getting this split wrong was the actual
// struggle here (see jodiemath-workflow memory, 2026-07-06): an earlier
// attempt used two-product only for the dominant term and treated every
// other cross term as "small enough to ignore imprecision in", which
// silently dropped terms that were negligible in *rounding error* but not
// in *value*, corrupting the reduced angle by ~1e-3. This version was
// checked against an arbitrary-precision reference (wolframscript) before
// being trusted. Result: fully accurate out to ~1e12-1e13 (vs. the old
// code's ~5e7), gracefully degrading (not exploding) out past 1e19.
const PI_W0: f32 = 3.1415927410125732;
const PI_W1: f32 = -8.742277657347586e-8;
const PI_W2: f32 = -3.4302490200117637e-15;
const RPI_W0: f32 = 0.31830987334251404;
const RPI_W1: f32 = 1.2841276486597053e-8;
const RPI_W2: f32 = 1.4685477398157775e-16;

#[inline(always)]
fn two_sum(a: f32, b: f32) -> (f32, f32) {
    let s = a + b;
    let v = s - a;
    let e = (a - (s - v)) + (b - v);
    (s, e)
}
// error-free product: p+e == a*b exactly, for any a, b (no overflow) --
// unlike the crate's other PI_A..D exact-multiply trick, this doesn't need
// q to stay under some bound.
#[inline(always)]
fn two_prod(a: f32, b: f32) -> (f32, f32) {
    let p = a * b;
    let e = fma(a, b, -p);
    (p, e)
}

/// round(x/pi + pre_offset), split into an exact double-float integer pair
/// (qh, ql). pre_offset is 0 for sin, -0.5 for cos (both exact additions).
#[inline(always)]
fn round_x_over_pi(x: f32, pre_offset: f32) -> (f32, f32) {
    let (p, e0) = two_prod(x, RPI_W0);
    let (s, e1) = two_sum(e0, x * RPI_W1);
    let v_lo = fma(x, RPI_W2, s + e1);
    let p = p + pre_offset;
    let qh = p.round();
    let rem = (p - qh) + v_lo;
    let ql = rem.round();
    (qh, ql)
}

/// x - (qh+ql)*pi, accurate to near f32 ulp even when qh/ql are far beyond
/// a single f32's exact-integer range (see the module doc comment above).
#[inline(always)]
fn reduce_pi(x: f32, qh: f32, ql: f32) -> f32 {
    let (p1, e1) = two_prod(qh, PI_W0);
    let (p2, e2) = two_prod(qh, PI_W1);
    let (p3, e3) = two_prod(ql, PI_W0);
    let c5 = ql * PI_W1;
    let c45 = fma(qh, PI_W2, c5); // = qh*PI_W2 + ql*PI_W1, one op cheaper
    // e2, e3, c45 are all comparably tiny (~x*2^-48): combining them via
    // plain adds first is safe (their combination with *each other*
    // doesn't need compensation, only their interaction with the much
    // bigger p1/s does) and cuts the two_sum loop from 7 iterations to 4 --
    // verified bit-for-bit identical against the 7-term version over the
    // full accuracy.rs sweep before adopting. (A further restructuring,
    // pre-combining e1/p2/p3 independent of the x-p1 subtraction to
    // shorten the per-element critical path, was tried and *worsened*
    // both latency and throughput -- this is a throughput-bound
    // 16-wide-vectorized region (examples/mca.rs bottleneck-analysis:
    // ~39% resource pressure), so shortening one element's dependency
    // chain doesn't help when port pressure across all the parallel
    // elements is already the binding constraint.)
    let tier2 = (e2 + e3) + c45;
    let (mut s, mut err) = two_sum(x, -p1);
    for t in [e1, p2, p3, tier2] {
        let (s2, e2) = two_sum(s, -t);
        s = s2;
        err += e2;
    }
    s + err
}

// parity of an exact-integer float q via floor-based "mod 2" (q*0.5 and
// its floor stay exact since q is already an integer), not `q as i64`: a
// cast looked simpler but Rust's float-to-int cast is *saturating* (clamps
// out-of-range/NaN instead of wrapping), which LLVM can't lower to a
// single vector instruction -- confirmed via llvm-mca's --emit=asm output:
// it fell back to extracting every lane and doing a scalar vcvttsd2si plus
// a NaN/range compare-and-cmov per element, which alone was most of this
// function's throughput cost in an earlier (f64-reduction) version of this
// fix. Staying in float land (floor, same instruction family as the
// .round() above) avoids that fallback entirely.
#[inline(always)]
fn parity(q: f32) -> f32 {
    q - 2.0 * (q * 0.5).floor()
}

#[inline(always)]
pub fn sin(x: f32) -> f32 {
    let (qh, ql) = round_x_over_pi(x, 0.0);
    let r = reduce_pi(x, qh, ql);
    let s = sinf_poly(r);
    // sin(x) = (-1)^q * sin(r); q = qh+ql, so parity(q) = (parity(qh) +
    // parity(ql)) mod 2
    let par = parity(parity(qh) + parity(ql));
    s * (1.0 - 2.0 * par)
}
#[inline(always)]
pub fn cos(x: f32) -> f32 {
    // k = round(x/pi - 0.5), q = k + 0.5, r = x - q*pi in [-pi/2, pi/2]
    let (kh, kl) = round_x_over_pi(x, -0.5);
    // q = k + 0.5; fold the 0.5 into the small word kl, not the (possibly
    // huge) kh word, since kh + 0.5 silently rounds away once kh's own ulp
    // exceeds 1 -- kl stays small enough that + 0.5 is always exact
    let r = reduce_pi(x, kh, kl + 0.5);
    let s = sinf_poly(r);
    // cos(x) = (-1)^(k+1) * sin(r); k = kh+kl
    let par = parity(parity(kh) + parity(kl));
    s * (2.0 * par - 1.0)
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
