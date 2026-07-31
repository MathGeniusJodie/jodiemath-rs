// Coordinate-descent ULP tuner for polynomial coefficients.
// Scalar (non-SIMD, no portable_simd/nightly needed) 1.0-ulp reference,
// same crate accuracy.rs uses for its vectorized ground truth.
//
// This file copies fitted polynomial coefficients from src/lib.rs (and
// perturbs/tunes them), so some literals are expected to land close to
// well-known constants without actually being them -- not worth per-site
// `#[allow]`s the way the shipped library documents each one individually.
#![allow(clippy::approx_constant)]
use sleef::f64::erf_u10 as erf_ref;
use sleef::f64::erfc_u15 as erfc_ref;

#[inline(always)]
fn fma(a: f32, b: f32, c: f32) -> f32 {
    a.mul_add(b, c)
}

// ULP distance along the monotonic ordering of f32 bit patterns. NaN vs
// NaN is always exactly 0, whatever the two bit patterns are: payload,
// sign and quiet bits carry no numeric meaning, so a difference there is
// not an accuracy difference. Exactly one side being NaN stays maximal.
fn ulp_diff(a: f32, b: f32) -> u64 {
    fn ord(x: f32) -> i64 {
        let b = x.to_bits();
        if b & 0x8000_0000 != 0 { -((b & 0x7fff_ffff) as i64) } else { b as i64 }
    }
    if a.is_nan() || b.is_nan() {
        return if a.is_nan() == b.is_nan() { 0 } else { u64::MAX };
    }
    (ord(a) - ord(b)).unsigned_abs()
}

const EXPONENT_MASK: u32 = 0x7f800000;

#[inline(always)]
fn exp2_c(x: f32, c: &[f32]) -> f32 {
    // Must mirror src/lib.rs's real exp2/exp2_checked exactly: 3 balanced
    // Estrin-in-f2 pairs (g0/g1/g2), not the older "two degree-2 Horner
    // halves times exp2int*f^4/exp2int*f" A/B split this function used to
    // compute (same target polynomial in exact arithmetic, but a
    // different fma grouping -- and so different real f32 rounding).
    // Found and fixed 2026-07-10 after the same class of bug (Horner vs
    // Estrin) was caught in erf_tail_c; see IDEAS.md.
    let k = x.floor();
    let f = x - k;
    let exp2int = f32::from_bits(((k + 383_f32).to_bits() << 8) & EXPONENT_MASK);
    // c indices match the crate's own established convention (see
    // exp2_direct_c's doc comment): c[0..2) is g2 (f^4,f^5 coefficients),
    // c[2..4) is g1 (f^2,f^3), c[4..6) is g0 (f^0,f^1).
    let f2 = f * f;
    let g2 = fma(c[0], f, c[1]);
    let g1 = fma(c[2], f, c[3]);
    let g0 = fma(c[4], f, c[5]);
    let h = fma(g2, f2, g1);
    let q = fma(h, f2, g0);
    fma(q, exp2int * f, exp2int)
}

// idea #80 (2026-07-10): direct P(f)=2^f fit, same 3-balanced-pair
// Estrin-in-f2 shape and same 6 stored coefficients as the shipped Q(f)=
// (2^f-1)/f form above, but c[0] pinned to exactly 1.0 (P(0)=2^0=1) so the
// final combine is a single `p*exp2int` multiply instead of exp2_c's own
// `exp2int*f` (extra multiply) + final fma. This is genuinely NOT the same
// fit as exp2_c algebraically re-expressed: exp2_c's Q(f) is a free degree-5
// fit (6 free coefficients), and P(f)=1+f*Q(f) inherits degree 6 with all 6
// of Q's coefficients still free -- but this direct P(f) construction is
// only degree 5 with c[0] forced, leaving just 5 free coefficients (c[1]
// through c[5]), one fewer degree of freedom than exp2_c's own implicit
// degree-6 fit. Tuned with tune_fixed0 (never perturbs c[0]) to check
// whether that one fewer degree of freedom still clears exp2's own
// essentially-zero headroom (already max ulp 1, avg ~0.03) -- see the
// tune("exp2_direct", ...) call in main() and IDEAS.md's own writeup of
// the result.
#[inline(always)]
fn exp2_direct_c(x: f32, c: &[f32]) -> f32 {
    let k = x.floor();
    let f = x - k;
    let exp2int = f32::from_bits(((k + 383_f32).to_bits() << 8) & EXPONENT_MASK);
    let f2 = f * f;
    let g0 = fma(c[1], f, c[0]); // c[0] pinned to 1.0 exactly
    let g1 = fma(c[3], f, c[2]);
    let g2 = fma(c[5], f, c[4]);
    let h = fma(g2, f2, g1);
    let p = fma(h, f2, g0);
    p * exp2int
}

// IDEAS.md backlog round 3 #1, tried 2026-07-09/10 and rejected on every
// one of the 6 functions sharing this poly (see IDEAS.md's "exp2 / exp /
// sin / cos / tan" section for the full writeup): k=round(x) (f in
// [-0.5,0.5]) instead of k=floor(x) (f in [0,1)) -- same 3-balanced-pair
// Estrin-in-f2 structure and degree (6 coefficients, degree 5) as the
// shipped exp2_c, just a centered/halved domain. This coarse grid genuinely
// did show a real improvement (max ulp 2->1, avg 0.20307->0.04845, see
// below) -- the *rejection* came from real verification downstream (real
// fuzz, mca, edgecheck all disagreeing with this grid in different ways
// for different functions), not from this tuning step itself being wrong.
// Kept as reference infra, not wired into src/lib.rs. c indices match
// exp2_c's own g0/g1/g2 pairing: c[0..2) is g2 (f^4,f^5), c[2..4) is g1
// (f^2,f^3), c[4..6) is g0 (f^0,f^1).
const ROUND_MAGIC_TUNE: f32 = 12582912.0; // 1.5 * 2^23
#[inline(always)]
fn exp2_round_c(x: f32, c: &[f32]) -> f32 {
    let k = (x + ROUND_MAGIC_TUNE) - ROUND_MAGIC_TUNE;
    let f = x - k;
    let exp2int = f32::from_bits(((k + 383_f32).to_bits() << 8) & EXPONENT_MASK);
    let f2 = f * f;
    let g0 = fma(c[4], f, c[5]);
    let g1 = fma(c[2], f, c[3]);
    let g2 = fma(c[0], f, c[1]);
    let h = fma(g2, f2, g1);
    let q = fma(h, f2, g0);
    fma(q, exp2int * f, exp2int)
}

// Same idea, but degree 4 (5 coefficients: f^0..f^4) -- tests the
// backlog's "halving the interval may allow degree 5->4" claim directly.
// Grouping: g0 = c[3]+c[4]*f (f^0,f^1), g1 = c[1]+c[2]*f (f^2,f^3), plus a
// lone f^4 term c[0] (odd coefficient count, no partner) folded in via
// Horner-in-f2: q = g0 + f2*(g1 + f2*c[0]).
#[inline(always)]
fn exp2_round4_c(x: f32, c: &[f32]) -> f32 {
    let k = (x + ROUND_MAGIC_TUNE) - ROUND_MAGIC_TUNE;
    let f = x - k;
    let exp2int = f32::from_bits(((k + 383_f32).to_bits() << 8) & EXPONENT_MASK);
    let f2 = f * f;
    let g0 = fma(c[4], f, c[3]);
    let g1 = fma(c[2], f, c[1]);
    let h = fma(f2, c[0], g1);
    let q = fma(h, f2, g0);
    fma(q, exp2int * f, exp2int)
}

// IDEAS.md's "Select-tree LUT for exp2" idea, tried and rejected
// (2026-07-08): f=x-floor(x) in [0,1) split into an 8-way quantized
// f_hi (a 3-level blend over 2^(i/8), i=0..7, exact f32 constants) plus
// a residual f_lo in [0,1/8) fit directly (not the (2^f-1)/f trick,
// since f_hi already carries the "-1" baseline). c[0] is fixed at 1.0
// (2^0 at f_lo=0), tuned with tune_fixed0. A scipy feasibility check
// found k=4 (the backlog's own framing) needs degree 3 to only reach
// 2.94e-7 max relerr, ~24x worse than the shipped degree-5 form's
// 1.22e-8 -- not competitive; k=8 (one more blend level than "2-level
// vblendvps" suggested) with degree 3 reaches 1.84e-8, close. But real
// tuning against the actual ulp objective (coordinate descent, seeded
// from the scipy fit, not zero) landed at max ulp 3 / avg 0.43 on this
// file's own exp2 grid, vs. the shipped form's max 2 / avg 0.20 on the
// identical grid -- a real, measured regression even at k=8, despite
// the scipy math suggesting near-parity. Not chased further (no
// implementation bug found in a quick review of the blend-tree/f_lo
// boundary consistency); disqualified on accuracy alone before even
// checking mca for the hoped-for latency win. Not implemented in
// src/lib.rs; kept as reference infra.
const EXP2_LUT8: [f32; 8] = [
    1.0, 1.0905077, 1.1892071, 1.2968396, 1.4142135, 1.5422108, 1.6817929, 1.8340081,
];
#[inline(always)]
fn exp2_lut8_c(x: f32, c: &[f32]) -> f32 {
    let k = x.floor();
    let f = x - k;
    let exp2int = f32::from_bits(((k + 383_f32).to_bits() << 8) & EXPONENT_MASK);
    // 3-level blend tree selecting among the 8 precomputed 2^(i/8)
    // constants by comparing f against 1/8-spaced thresholds -- mirrors
    // this crate's established "compute both arms unconditionally, then
    // select" branchless idiom, just nested three deep instead of one.
    let b0 = if f < 4.0 / 8.0 {
        if f < 2.0 / 8.0 {
            if f < 1.0 / 8.0 { EXP2_LUT8[0] } else { EXP2_LUT8[1] }
        } else if f < 3.0 / 8.0 {
            EXP2_LUT8[2]
        } else {
            EXP2_LUT8[3]
        }
    } else if f < 6.0 / 8.0 {
        if f < 5.0 / 8.0 { EXP2_LUT8[4] } else { EXP2_LUT8[5] }
    } else if f < 7.0 / 8.0 {
        EXP2_LUT8[6]
    } else {
        EXP2_LUT8[7]
    };
    let idx = (f * 8.0).floor();
    let f_lo = f - idx * (1.0 / 8.0);
    let f_lo2 = f_lo * f_lo;
    let poly = fma(fma(c[3], f_lo, c[2]), f_lo2, fma(c[1], f_lo, c[0]));
    exp2int * b0 * poly
}

// Same Cody-Waite splits as src/lib.rs's private LN2_HI/LN2_LO and
// LOG10_2_HI/LOG10_2_LO (not pub, so redefined here verbatim).
const LN2_HI: f32 = 0.693145751953125;
const LN2_LO: f32 = 1.428606765330187e-6;
const LOG10_2_HI: f32 = 0.301025390625;
const LOG10_2_LO: f32 = 4.605039066518657e-6;

// ln_normal/log10_normal's exact shape (see src/lib.rs), refitting the
// poly directly against ln/log10 instead of log_2's own coefficients
// individually rescaled by LN_2/LOG10_2. c[0] is each function's
// mathematically-required exact leading term (1.0 for ln, LOG10_2 for
// log10 -- log10(1+s)'s derivative at s=0), so both are tuned with
// tune_fixed0.
#[inline(always)]
fn ln_poly_c(x: f32, c: &[f32]) -> f32 {
    let e = (x.to_bits() as i32).wrapping_sub(0x3f3504f3) >> 23;
    let m = f32::from_bits((x.to_bits() as i32).wrapping_sub(e << 23) as u32);
    let k = e as f32;
    let s = m - 1.0;
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
    let k_hi = k * LN2_HI;
    fma(p, s, k_hi) + k * LN2_LO
}

// backlog idea #22: log1p(x) = ln(u) + c/u (u=1+x, c=x-(u-1), the exact
// Sterbenz-recovered rounding error), matching src/lib.rs's own log1p
// exactly -- ln(u) routed through ln_poly_c so the SAME candidate
// coefficients affect both ln's own accuracy and log1p's, letting a
// joint objective refit ln's coefficients specifically for log1p's own
// benefit without changing ln's own op count (log1p already reuses ln
// unmodified in the shipped code, so any coefficient set that helps
// here is a free win, zero new ops).
#[inline(always)]
fn log1p_via_ln_c(x: f32, c: &[f32]) -> f32 {
    let u = 1.0 + x;
    let corr_num = x - (u - 1.0);
    let corr = corr_num / u;
    let corr = if corr.is_finite() { corr } else { 0.0 };
    ln_poly_c(u, c) + corr
}

#[inline(always)]
fn log10_poly_c(x: f32, c: &[f32]) -> f32 {
    let e = (x.to_bits() as i32).wrapping_sub(0x3f3504f3) >> 23;
    let m = f32::from_bits((x.to_bits() as i32).wrapping_sub(e << 23) as u32);
    let k = e as f32;
    let s = m - 1.0;
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
    let k_hi = k * LOG10_2_HI;
    fma(p, s, k_hi) + k * LOG10_2_LO
}

#[inline(always)]
fn log2_c(x: f32, c: &[f32]) -> f32 {
    let e = (x.to_bits() as i32).wrapping_sub(0x3f3504f3) >> 23;
    let m = f32::from_bits((x.to_bits() as i32).wrapping_sub(e << 23) as u32);
    let k = e as f32;
    let s = m - 1.0;
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
    fma(p, s, k)
}

// IDEAS.md's "atanh-form reduction" idea, tried and rejected (2026-07-08):
// t = (m-1)/(m+1) is an *odd* series in t (half the domain of s = m-1,
// ~0.172 vs ~0.414, thanks to odd symmetry), so log2(m) = c_const*t*Q(t^2)
// needs far fewer coefficients than log_2's shipped s*P(s) for the same
// *mathematical* accuracy (scipy: degree-4 Q, 5 coefficients, hits max
// relative error 1.2e-11 vs the shipped degree-9/10-coefficient form's
// 4.1e-9 -- a real ~1000x tighter fit with half the coefficients).
// Implemented directly in src/lib.rs and measured for real: accuracy
// actually regressed slightly (fuzz avg 0.003->0.0053, max 3->5 --
// log_2 was already deep in pure-rounding-noise territory, far past the
// point where a tighter *mathematical* fit moves the *measured* ulp), and
// latency got much *worse*, not better: mca 34.23->49.05 cyc (+43%),
// with throughput barely moving (1.556->1.521, ~2%). The one division
// this form needs (t = s/(m+1)) depends on `s`, available from the very
// first step of the critical path, so unlike cbrt's early-starting rcp
// there's no independent work for it to hide behind -- its latency (this
// file's own atanh/tanh entries put a division around ~11 cyc) plus the
// reduced-but-still-serial poly evaluation came out *slower* overall than
// the original's longer but division-free chain. Reverted; src/lib.rs
// unchanged. c[0] is Q(0) = atanh'(0)/1 = 1.0 exactly (mathematically
// required, like log_2_normal's own c[0]), tuned with tune_fixed0.
// c_const = 2*log2(e) is likewise exact, not a free parameter.
const LOG2_ATANH_CONST: f32 = 2.0 * std::f32::consts::LOG2_E;
#[inline(always)]
fn log2_atanh_c(x: f32, c: &[f32]) -> f32 {
    let e = (x.to_bits() as i32).wrapping_sub(0x3f3504f3) >> 23;
    let m = f32::from_bits((x.to_bits() as i32).wrapping_sub(e << 23) as u32);
    let k = e as f32;
    let s = m - 1.0;
    let t = s / (m + 1.0);
    let u = t * t;
    let u2 = u * u;
    let r0 = fma(c[1], u, c[0]);
    let r1 = fma(c[3], u, c[2]);
    let q = fma(c[4], u2 * u2, fma(r1, u2, r0));
    fma(t * q, LOG2_ATANH_CONST, k)
}

fn score(
    f: &dyn Fn(f32, &[f32]) -> f32,
    reference: &dyn Fn(f64) -> f64,
    grid: &[f32],
    c: &[f32],
) -> (u64, u64) {
    let mut sum = 0u64;
    let mut max = 0u64;
    for &x in grid {
        let r = reference(x as f64) as f32;
        let d = ulp_diff(f(x, c), r);
        sum += d;
        max = max.max(d);
    }
    (max, sum)
}

#[inline(always)]
fn mulsign_c(x: f32, y: f32) -> f32 {
    f32::from_bits(x.to_bits() ^ (y.to_bits() & 0x8000_0000))
}

// asin's mid-branch rational correction (see src/lib.rs's asin doc
// comment): a2 = (a^2-a)/d + a where d is this 3-coefficient poly (plus
// the leading constant folded in as c[3]), then the already-fixed
// rationalized sqrt step. Only the correction's 4 coefficients are being
// tuned here -- the rationalization and mulsign/branch structure are
// fixed, matching what's actually shipped in src/lib.rs.
#[inline(always)]
fn asin_mid_c(x: f32, c: &[f32]) -> f32 {
    let a = x.abs();
    let d = fma(-c[0], a, c[1]);
    let d = fma(-a, d, c[2]);
    let d = fma(-a, d, c[3]);
    let a2 = (a * a - a) / d + a;
    let sq = (1.0 - a2).sqrt();
    let sm1 = -a2 / (sq + 1.0);
    mulsign_c(sm1, x) * (-std::f32::consts::FRAC_PI_2)
}

// atan_poly (see src/lib.rs): a Pade form approximating atan(x) directly
// for x in [0,1] -- atan() calls this on a.min(1/a), always in that
// range, so that's exactly the grid to tune against.
#[inline(always)]
fn atan_poly_c(x: f32, c: &[f32]) -> f32 {
    let x2 = x * x;
    (fma(fma(c[0], x2, c[1]), x2, 1.0) * x) / fma(fma(x2, c[2], c[3]), x2, 1.0)
}

// IDEAS.md's "Three-interval atan reduction" idea, tried and rejected
// (2026-07-08): on top of a.min(1/a) (giving r in [0,1]), split further
// at t=tan(pi/8) -- atan(r) = poly(r) for r<t directly, or pi/8 + poly(u)
// for r>=t via u=(r-t)/(1+t*r) (both r and u land in [0,t], a much
// smaller domain than atan_poly's full [0,1], letting a degree-2/2
// rational match or beat the shipped degree-3/3's accuracy with 2 fewer
// coefficients -- scipy check: 8.7e-10 max relerr over [-t,t] vs the
// shipped 3/3's 1.87e-8 over [0,1], confirmed for real ulp too: this
// tuned to max 3/avg 0.275, matching/beating shipped's max 3/avg 0.286).
// The accuracy case fully held up -- but implemented for real, mca
// showed a severe regression, not a modest cost: latency 61.09->104.94
// cyc (+72%), throughput 1.491->2.974 (+99%, nearly doubled). Root
// cause: the u-transform's own division sits *before* the (still-a-
// rational) poly's own division for every r>=t (about half the domain),
// two sequential full-latency divisions instead of one, with nothing to
// overlap them against -- the same risk pattern that sank the log_2
// atanh-form idea, here worse than expected even having been flagged in
// advance. Not adopted; kept as reference infra given how clean the
// accuracy result was on its own.
const ATAN_T: f32 = 0.41421356237309503; // tan(pi/8)
#[inline(always)]
fn atan_three_c(x: f32, c: &[f32]) -> f32 {
    let r = x; // caller already reduces to a.min(1/a) in [0,1]
    let below = r < ATAN_T;
    let u = (r - ATAN_T) / fma(ATAN_T, r, 1.0);
    let arg = if below { r } else { u };
    let arg2 = arg * arg;
    let numer = fma(fma(c[0], arg2, c[1]), arg2, 1.0) * arg;
    let denom = fma(fma(arg2, c[2], c[3]), arg2, 1.0);
    let p = numer / denom;
    if below { p } else { std::f32::consts::FRAC_PI_8 + p }
}

// Degree bump on atan_poly: 3/3 instead of 2/2 (one more term each in
// numer/denom). c = [a2, a1, a0, b2, b1, b0]; a2/b2 start at 0.0 so this
// starts bit-identical to the shipped 2/2 form.
#[inline(always)]
fn atan_poly7_c(x: f32, c: &[f32]) -> f32 {
    let x2 = x * x;
    let numer = fma(fma(fma(c[0], x2, c[1]), x2, c[2]), x2, 1.0) * x;
    let denom = fma(fma(fma(c[3], x2, c[4]), x2, c[5]), x2, 1.0);
    numer / denom
}

// IDEAS.md's "Pure-poly atan latency tier" idea: a division-free odd
// degree-17 poly (c[0..8), d0=1.0 fixed, tuned with tune_fixed0 via a
// leading dummy 1.0 slot... actually c has no d0 slot, tune_fixed0
// isn't used here, c[0..8) map directly to d1..d8) instead of
// atan_poly's 3/3 rational -- removes the division atan_poly's own
// critical path pays, trading it for more fma depth. Split into a
// low half (d0..d3) and high half (d4..d8) so both halves' Horner
// chains can evaluate in parallel before one final combine, instead of
// one long degree-8 Horner chain.
#[inline(always)]
fn atan_pure_poly_c(x: f32, c: &[f32]) -> f32 {
    let x2 = x * x;
    let x4 = x2 * x2;
    let lo = fma(fma(fma(c[2], x2, c[1]), x2, c[0]), x2, 1.0);
    let hi = fma(fma(fma(fma(c[7], x2, c[6]), x2, c[5]), x2, c[4]), x2, c[3]);
    fma(hi, x4 * x4, lo) * x
}

// acos_poly (see src/lib.rs), scored as the *whole* acos(x) formula for
// x >= 0 (sqrt(1-x)*acos_poly(x)) rather than the bare poly value -- the
// sqrt factor's own rounding interacts with the poly, so tuning the poly
// in isolation could mistune it relative to what actually ships.
// acos_poly is also reused by asin's near-1 branch, so an improvement
// here benefits both callers.
#[inline(always)]
fn acos_poly_c(x: f32, c: &[f32]) -> f32 {
    let u = fma(c[0], x, c[1]);
    let u = fma(u, x, c[2]);
    let u = fma(u, x, c[3]);
    let u = fma(u, x, c[4]);
    let u = fma(u, x, c[5]);
    let poly = fma(u, x, c[6]);
    (1.0 - x).sqrt() * poly
}

// Degree-7 (8-coefficient) version of acos_poly_c, for the "acos_poly
// degree 6 -> 7" idea: one more Horner term/fma than shipped, screening
// whether the extra degree of freedom buys real accuracy headroom.
#[inline(always)]
fn acos_poly8_c(x: f32, c: &[f32]) -> f32 {
    let u = fma(c[0], x, c[1]);
    let u = fma(u, x, c[2]);
    let u = fma(u, x, c[3]);
    let u = fma(u, x, c[4]);
    let u = fma(u, x, c[5]);
    let u = fma(u, x, c[6]);
    let poly = fma(u, x, c[7]);
    (1.0 - x).sqrt() * poly
}

// Dedicated asin-only copy of acos_poly's shape (backlog idea #34):
// every attempt to jointly fit acos_poly for both callers died protecting
// acos's own accuracy (three separate rejections logged in IDEAS.md), so
// this decouples them entirely -- same 7-coefficient Estrin/Horner shape,
// independently tuned, scored as asin's *own* actual combine (FRAC_PI_2 -
// sqrt(1-x)*poly, sign restored elsewhere) rather than acos's. Only
// scored where asin actually uses this branch (x >= 0.25); acos_poly
// itself is untouched by this, so acos's protected accuracy can't regress
// no matter what this converges to.
#[inline(always)]
fn asin_poly_c(x: f32, c: &[f32]) -> f32 {
    let u = fma(c[0], x, c[1]);
    let u = fma(u, x, c[2]);
    let u = fma(u, x, c[3]);
    let u = fma(u, x, c[4]);
    let u = fma(u, x, c[5]);
    let poly = fma(u, x, c[6]);
    std::f32::consts::FRAC_PI_2 - (1.0 - x).sqrt() * poly
}

// erf's tail branch (see src/lib.rs's erf): scored as the whole
// mulsign(1.0 - exp2(erf_poly(xa)), x) formula, xa in [0.28, 10] (exactly
// where this branch is used in the shipped code; the Pade near-zero
// branch below 0.28 is untouched). Uses f32::exp2 (std) as a stand-in for
// the shipped exp2_checked -- within this bounded range erf_poly never
// leaves exp2's safe domain, so the checked/unchecked distinction doesn't
// matter for tuning purposes.
#[inline(always)]
fn erf_tail_c(x: f32, c: &[f32]) -> f32 {
    let xa = x.abs().min(10.0);
    let u = fma(c[0], xa, c[1]);
    let u = fma(u, xa, c[2]);
    let u = fma(u, xa, c[3]);
    let u = fma(u, xa, c[4]);
    let u = fma(u, xa, c[5]);
    let poly = fma(u, xa, c[6]);
    mulsign_c(1.0 - poly.exp2(), x)
}

// erf's near-zero Pade branch (see src/lib.rs's erf), used for |x| < 0.28
// -- that's exactly the grid to tune against; the tail branch above 0.28
// is untouched.
#[inline(always)]
fn erf_near0_c(x: f32, c: &[f32]) -> f32 {
    let x2 = x * x;
    let numer = x * fma(c[0], x2, c[1]);
    let denom = fma(fma(c[2], x2, c[3]), x2, 1.0);
    numer / denom
}

// cbrt_normal (see src/lib.rs): the degree-3 correction poly's 4
// coefficients, against the *existing* ax/3-based seed (a separate,
// bigger question -- swapping to the cheaper (bits>>16)*0x5556 seed --
// is deferred, see IDEAS.md's cbrt seed entry).
#[inline(always)]
fn cbrt_normal_c(x: f32, c: &[f32]) -> f32 {
    const SIGN_MASK: u32 = 0x8000_0000;
    let ax = x.to_bits() & !SIGN_MASK;
    let a = f32::from_bits(ax);
    let rcp = 1.0 / a;
    let s = f32::from_bits(ax / 3 + 0x2a509a07u32);
    let s2 = s * s;
    let d = fma(s2, s, -a);
    let r = d * rcp;
    let r2 = r * r;
    let a1 = fma(c[1], r, c[0]);
    let b1 = fma(c[3], r, c[2]);
    let p = fma(b1, r2, a1);
    let ss = f32::from_bits(s.to_bits() | (x.to_bits() & SIGN_MASK));
    let sr = ss * r;
    fma(sr, p, ss)
}

// cbrt_normal_c with the seed offset (shipped: 0x2a509a07) as a 5th
// tunable parameter (c[0], its bits used directly as the offset added to
// ax/3) instead of a fixed literal -- IDEAS.md's "joint seed-constant +
// degree-3 coefficient search" idea. The previously-rejected joint search
// only tried this against a *degree-2* correction (3 coeffs, hopeless
// regardless of seed); this keeps the full shipped degree-3.
#[inline(always)]
fn cbrt_normal_joint_c(x: f32, c: &[f32]) -> f32 {
    const SIGN_MASK: u32 = 0x8000_0000;
    let ax = x.to_bits() & !SIGN_MASK;
    let a = f32::from_bits(ax);
    let rcp = 1.0 / a;
    let s = f32::from_bits(ax / 3 + c[0].to_bits());
    let s2 = s * s;
    let d = fma(s2, s, -a);
    let r = d * rcp;
    let r2 = r * r;
    let a1 = fma(c[2], r, c[1]);
    let b1 = fma(c[4], r, c[3]);
    let p = fma(b1, r2, a1);
    let ss = f32::from_bits(s.to_bits() | (x.to_bits() & SIGN_MASK));
    let sr = ss * r;
    fma(sr, p, ss)
}

// cbrt_normal_c, but with cbrt_fast's shift-multiply seed
// ((ax>>16)*0x5556 + 0x2a4ddef1, cbrt_fast's own constants, reused as-is)
// instead of the division-based ax/3 + 0x2a509a07 -- checking whether the
// existing degree-3 correction poly, refit for this cheaper-but-cruder
// seed's different error distribution, can still meet cbrt's budget.
// Tried and rejected (2026-07-07): even after retuning, best found is
// max ulp 33 / avg 5.07 on a coarse grid, ~16x over budget (avg<=1,
// max<=2) -- see IDEAS.md's cbrt seed entry. Kept as reference
// infrastructure (`which.contains("cbrtshift")`), not deleted, since it's
// small and inert unless invoked; a future attempt would need freshly-
// derived seed constants for this specific poly-degree combination
// rather than reusing cbrt_fast's own (tuned for its different
// downstream Newton-refinement structure, not a degree-3 correction).
#[inline(always)]
fn cbrt_shiftmul_c(x: f32, c: &[f32]) -> f32 {
    const SIGN_MASK: u32 = 0x8000_0000;
    let ax = x.to_bits() & !SIGN_MASK;
    let a = f32::from_bits(ax);
    let rcp = 1.0 / a;
    let s = f32::from_bits((ax >> 16) * 0x5556u32 + 0x2a4ddef1u32);
    let s2 = s * s;
    let d = fma(s2, s, -a);
    let r = d * rcp;
    let r2 = r * r;
    let a1 = fma(c[1], r, c[0]);
    let b1 = fma(c[3], r, c[2]);
    let p = fma(b1, r2, a1);
    let ss = f32::from_bits(s.to_bits() | (x.to_bits() & SIGN_MASK));
    let sr = ss * r;
    fma(sr, p, ss)
}

// cbrt_throughput (see src/lib.rs): experimental positive-only 2-iteration
// Halley-style inverse-cbrt refinement with a magic bit-trick seed.
// c[0]'s bits are the seed subtrahend (shipped: 0xd461ff81), c[1]/c[2] are
// the two Halley iterations' own f32 multiplier constants.
// Tried and rejected (2026-07-08), IDEAS.md's "tune cbrt_throughput's
// magic constants" idea: unlike cbrt_normal, this function's error does
// NOT repeat across octaves (checked directly: max ulp ranges ~7 to ~52
// depending which octave is sampled), so a single-octave grid tune looked
// slightly better on-grid (max 16->12) but was real-world *worse* on the
// full fuzz sweep (avg 6.73->7.01, max 74->76) -- not representative. A
// proper ~60-octave grid tune found an even worse tradeoff: max ulp did
// drop (68->43 on-grid) but avg ulp exploded (6.83->22.73), because
// `score()`'s tuple ordering `(max, sum)` optimizes max first and only
// uses avg as a tiebreak -- fine for cbrt_normal's smooth error surface,
// but for this rougher one it happily wrecks the average to shave the
// worst case. Not adopted either way; src/lib.rs unchanged. Kept as
// reference infrastructure (`which.contains("cbrtthroughput")`), same as
// cbrt_shiftmul_c below.
#[inline(always)]
fn cbrt_throughput_c(x: f32, c: &[f32]) -> f32 {
    let r = f32::from_bits(c[0].to_bits().wrapping_sub(x.to_bits() / 3));
    let r = fma(r * r, (r * r) * x, r * c[1]);
    let r = fma(r * r, (r * r) * x, r * c[2]);
    r * r * x
}

// idea #52: cbrt_normal_c's degree-4 correction poly (5 coefficients
// instead of shipped's 4) -- one more fma than shipped, one fewer than
// the previously-tried (and rejected on cost) degree-5.
#[inline(always)]
fn cbrt_normal_c5(x: f32, c: &[f32]) -> f32 {
    const SIGN_MASK: u32 = 0x8000_0000;
    let ax = x.to_bits() & !SIGN_MASK;
    let a = f32::from_bits(ax);
    let rcp = 1.0 / a;
    let s = f32::from_bits(ax / 3 + 0x2a509a07u32);
    let s2 = s * s;
    let d = fma(s2, s, -a);
    let r = d * rcp;
    let r2 = r * r;
    let a1 = fma(c[1], r, c[0]);
    let b1 = fma(fma(c[4], r, c[3]), r, c[2]);
    let p = fma(b1, r2, a1);
    let ss = f32::from_bits(s.to_bits() | (x.to_bits() & SIGN_MASK));
    let sr = ss * r;
    fma(sr, p, ss)
}

// sinf_poly (see src/lib.rs): sin(x) ~= x + x^3*P(x^2) on [-pi/2, pi/2],
// the shared poly behind sin/cos/sin_checked/cos_checked.
#[inline(always)]
fn sinf_poly_c(x: f32, c: &[f32]) -> f32 {
    let y = x * x;
    let y2 = y * y;
    let x3 = y * x;
    let a = fma(c[1], y, c[0]);
    let b = fma(c[3], y, c[2]);
    let p = fma(b, y2, a);
    fma(p, x3, x)
}

// idea #43: dedicated sinpi/cospi poly, fit directly against sin(pi*r)
// in r (r in [-0.5,0.5]) instead of routing r through an intermediate
// `PI*r` value (which only keeps PI_HI's precision, discarding PI_LO*r)
// before calling the generic sinf_poly. Same op count as sinf_poly_raw
// applied to a pre-scaled `PI*r`: one multiply for the leading term
// either way, just moved from the caller (computing `PI*r` before the
// call) to here (computing `r*PI_HI` as this function's own leading
// term) -- zero added ops, only different coefficient values, refit to
// approximate sin(pi*r) directly instead of sin(t) generically.
#[inline(always)]
fn sinpi_poly_c(r: f32, c: &[f32]) -> f32 {
    const PI_HI: f32 = 3.1415927410125732;
    let y = r * r;
    let y2 = y * y;
    let r3 = y * r;
    let a = fma(c[1], y, c[0]);
    let b = fma(c[3], y, c[2]);
    let p = fma(b, y2, a);
    fma(p, r3, r * PI_HI)
}

// idea #44: same fold mechanism as #43 above (sinpi_poly_c), applied to
// sind/cosd's own reduction: fit sin(d*pi/180) directly in d (d in
// [-90,90]) instead of routing d through the intermediate `d *
// DEG_TO_RAD_SMALL` rounding step before calling the generic sinf_poly.
#[inline(always)]
fn sind_poly_c(d: f32, c: &[f32]) -> f32 {
    const DEG_TO_RAD_SMALL: f32 = std::f32::consts::PI / 180.0;
    let y = d * d;
    let y2 = y * y;
    let d3 = y * d;
    let a = fma(c[1], y, c[0]);
    let b = fma(c[3], y, c[2]);
    let p = fma(b, y2, a);
    fma(p, d3, d * DEG_TO_RAD_SMALL)
}

// expm1's Pade near-zero branch (see src/lib.rs's expm1), |x| < 0.5. The
// shipped constants (-2, -120, -12, 60, -120) look like an exact closed-
// form Pade approximant to e^x rather than an empirical lolremez fit
// (small integers, -120 shared between numerator and denominator) --
// tuning treats them as 5 independent free parameters regardless, since
// coordinate descent doesn't care about the identity's elegance.
#[inline(always)]
fn expm1_near0_c(x: f32, c: &[f32]) -> f32 {
    let numer = x * fma(c[0], x * x, c[1]);
    let denom = fma(x, fma(x, x + c[2], c[3]), c[4]);
    numer / denom
}

// expm1_near0_c with the numerator bumped from degree 3 to degree 5 (one
// more odd term, x*(c0*x^4+c1*x^2+c2) instead of x*(c0*x^2+c1)) -- for the
// "expm1 Pade degree bump" idea, screening whether the extra numerator
// degree buys real headroom over the shipped 5-coefficient form's max ulp
// 3 (in the |x|<0.5 branch). Denominator left at its shipped degree 3.
#[inline(always)]
fn expm1_near0_deg5_c(x: f32, c: &[f32]) -> f32 {
    let x2 = x * x;
    let numer = x * fma(fma(c[0], x2, c[1]), x2, c[2]);
    let denom = fma(x, fma(x, x + c[3], c[4]), c[5]);
    numer / denom
}

// erfcx_pos's degree-10 poly in the reciprocal variable (see src/lib.rs's
// erfcx_pos): `erfcx(xa) = v*P(v)` with `v = 1/(2+xa)`. The evaluation
// order below is the shipped Estrin grouping exactly, which is the point
// -- descent against a different grouping tunes for roundings the shipped
// code never makes.
//
// **All 11 coefficients are free, including `c0` -- deliberately, and
// against the obvious argument.** `c0` is `P(0)`, and at large `xa` the
// result *is* `c0*v` and nothing else, so `c0` alone carries the
// `erfcx(xa) ~ 1/(xa*sqrt(pi))` asymptote and the tempting move is to
// hardcode it to the correctly-rounded `1/sqrt(pi)` (0x3f106ebb) the way
// `exp_r_c` hardcodes its own fixed 1.0s. That was tried and **measured
// worse**: pinning it costs `erfc` a full max ulp (exhaustive 0.1993/7
// free vs 0.2003/8 pinned), because the other ten coefficients cannot
// re-absorb the constraint. The shipped `c0` is 1 ulp *below*
// correctly-rounded and buys that ulp back; what it costs is 0.023 avg
// ulp in the far tail (exhaustive `x >= 20`, 0.6245 -> 0.6478, max 4
// either way). Do not "fix" it without re-running both sweeps.
//
// This replaced a family of four targets (`erfc_c`/`erfc_c5`/
// `erfc_lo_c`/`erfc_hi_c`) for the degree-4/4 rational in `xa` that
// `erfcx_pos` retired; see graveyard.md for what they measured.
#[inline(always)]
fn erfcx_pos_c(xa: f32, c: &[f32]) -> f32 {
    let v = 1.0 / (2.0 + xa);
    let v2 = v * v;
    let v4 = v2 * v2;
    let p01 = fma(c[1], v, c[0]);
    let p23 = fma(c[3], v, c[2]);
    let p45 = fma(c[5], v, c[4]);
    let p67 = fma(c[7], v, c[6]);
    let t9 = fma(c[10], v, c[9]);
    let t8 = fma(t9, v, c[8]);
    let lo = fma(p23, v2, p01);
    let hi = fma(p67, v2, p45);
    v * fma(fma(t8, v4, hi), v4, lo)
}

// exp's e^r poly (see src/lib.rs's exp) with c0 AND c1 both forced to
// exactly 1.0 (hardcoded, not tuned) instead of just c0 -- only c2..c5 (4
// values, named c[0..4] here) are free. Matches the shipped Estrin
// evaluation order exactly (l0 = 1+r needs no fma since both its
// coefficients are exactly 1.0).
#[inline(always)]
fn exp_r_c(r: f32, c: &[f32]) -> f32 {
    let r2 = r * r;
    let r4 = r2 * r2;
    let l0 = r + 1.0;
    let l1 = fma(c[1], r, c[0]);
    let l2 = fma(c[3], r, c[2]);
    let r0 = fma(l1, r2, l0);
    fma(l2, r4, r0)
}

// idea #25: exp_r_c's degree 5->6 bump (one more coefficient, one more
// fma at the end of the Estrin chain -- see exp_r_c's own doc comment
// for the fixed-c0/c1 convention this preserves).
#[inline(always)]
fn exp_r_c6(r: f32, c: &[f32]) -> f32 {
    let r2 = r * r;
    let r4 = r2 * r2;
    let r6 = r4 * r2;
    let l0 = r + 1.0;
    let l1 = fma(c[1], r, c[0]);
    let l2 = fma(c[3], r, c[2]);
    let r0 = fma(l1, r2, l0);
    let r1 = fma(l2, r4, r0);
    fma(c[4], r6, r1)
}

// Scratch: sinh/cosh's shared exp(r)/exp(-r) even/odd split (see
// src/lib.rs's exp_pos_neg). e(u)=1+c0*u+c2*u^2, o(u)=1+c1*u+c3*u^2,
// p(r)=e+r*o. Same grid as exp_r (symmetric in r), so tuning this one
// function against exp(r) over +-r already covers both e+r*o and e-r*o.
#[inline(always)]
fn exp_r_pair_c(r: f32, c: &[f32]) -> f32 {
    let r2 = r * r;
    let r4 = r2 * r2;
    let e = fma(c[2], r4, fma(c[0], r2, 1.0));
    let o = fma(c[3], r4, fma(c[1], r2, 1.0));
    fma(r, o, e)
}

fn tune(
    name: &str,
    f: &dyn Fn(f32, &[f32]) -> f32,
    reference: &dyn Fn(f64) -> f64,
    grid: &[f32],
    init: &[f32],
) {
    let mut c: Vec<f32> = init.to_vec();
    let mut best = score(f, reference, grid, &c);
    println!("{name}: start max {} avg {:.5}", best.0, best.1 as f64 / grid.len() as f64);
    let mut improved = true;
    while improved {
        improved = false;
        for i in 0..c.len() {
            for delta in [1i32, -1, 2, -2, 4, -4, 8, -8, 16, -16] {
                let mut trial = c.clone();
                trial[i] = f32::from_bits((trial[i].to_bits() as i32 + delta) as u32);
                let s = score(f, reference, grid, &trial);
                if s < best {
                    best = s;
                    c = trial;
                    improved = true;
                }
            }
        }
    }
    println!(
        "{name}: tuned max {} avg {:.5}  coeffs: {:?}",
        best.0,
        best.1 as f64 / grid.len() as f64,
        c.iter().map(|v| format!("{v:e}")).collect::<Vec<_>>()
    );
}

// IDEAS.md's "Simulated annealing / basin-hopping over coefficient
// space" idea, tried on acos_poly and rejected (2026-07-08): tune()'s
// coordinate descent only moves one axis at a time, so it can't cross a
// "diagonal valley" (a direction where improvement requires two or more
// coefficients to move together, each individually making the score
// worse). Runs the usual single-axis descent to a local optimum first,
// then repeatedly perturbs 2-3 random coefficients simultaneously (a
// wider jump than descent's own +-16 step) and re-descends from there,
// keeping the result only if it beats the current best.
//
// Inherits the exact same max-first-tuple-comparison trap already found
// for cbrt_throughput's own tuning (`score()` returns `(max, sum)`,
// compared via Rust's default tuple ordering, max first): confirmed
// bitten by it directly -- on acos_poly with a coarse 1M-step grid (for
// speed; the full 10000-step grid makes even 50 restarts too slow to be
// practical, each restart re-running a full descent), 50 restarts found
// a candidate reporting max 2 vs. the single-axis descent's max 3 on
// that same coarse grid -- but on the *real* accuracy.rs fuzz (100M
// samples, the actual measure that matters), that candidate came out at
// max ulp 4 (same as shipped, not an improvement) and avg ulp 0.9611 vs.
// shipped's 0.4961 -- avg nearly *doubled*. Any result from this
// function needs the same real-fuzz verification `tune()`'s own results
// do; don't trust the coarse-grid "improvement" number by itself.
fn tune_basin_hop(
    name: &str,
    f: &dyn Fn(f32, &[f32]) -> f32,
    reference: &dyn Fn(f64) -> f64,
    grid: &[f32],
    init: &[f32],
    restarts: usize,
) {
    fn descend(
        f: &dyn Fn(f32, &[f32]) -> f32,
        reference: &dyn Fn(f64) -> f64,
        grid: &[f32],
        start: &[f32],
    ) -> (Vec<f32>, (u64, u64)) {
        let mut c = start.to_vec();
        let mut best = score(f, reference, grid, &c);
        let mut improved = true;
        while improved {
            improved = false;
            for i in 0..c.len() {
                for delta in [1i32, -1, 2, -2, 4, -4, 8, -8, 16, -16] {
                    let mut trial = c.clone();
                    trial[i] = f32::from_bits((trial[i].to_bits() as i32 + delta) as u32);
                    let s = score(f, reference, grid, &trial);
                    if s < best {
                        best = s;
                        c = trial;
                        improved = true;
                    }
                }
            }
        }
        (c, best)
    }

    let (mut c, mut best) = descend(f, reference, grid, init);
    println!(
        "{name}: single-axis descent max {} avg {:.5}",
        best.0,
        best.1 as f64 / grid.len() as f64
    );

    use rand::RngExt;
    let mut rng = rand::rng();
    for _ in 0..restarts {
        let mut trial_start = c.clone();
        let n_perturb = 2 + (rng.random::<u8>() % 2) as usize; // 2 or 3
        for _ in 0..n_perturb {
            let idx = (rng.random::<u32>() as usize) % trial_start.len();
            let delta = (rng.random::<i32>() % 64) - 32; // wider than descend's own +-16
            trial_start[idx] = f32::from_bits((trial_start[idx].to_bits() as i32 + delta) as u32);
        }
        let (c_trial, s_trial) = descend(f, reference, grid, &trial_start);
        if s_trial < best {
            best = s_trial;
            c = c_trial;
        }
    }
    println!(
        "{name}: basin-hopped ({restarts} restarts, max-first) max {} avg {:.5}  coeffs: {:?}",
        best.0,
        best.1 as f64 / grid.len() as f64,
        c.iter().map(|v| format!("{v:e}")).collect::<Vec<_>>()
    );
}

// idea #101: avg-first, max-capped basin-hopping -- the mechanism above
// is unchanged (single-axis descent, then random 2-3-coefficient
// perturbation + re-descend), but every comparison here uses `better()`
// instead of the default (max, sum) tuple ordering: a candidate whose
// max ulp exceeds `cap` (the un-hopped single-axis descent's own max,
// i.e. basin-hopping is never allowed to regress the max at all) is
// always worse than one that meets the cap; among candidates that meet
// the cap, lower sum (avg) wins. This directly targets the documented
// bias (`tune_basin_hop`'s own doc comment: "happily wrecks the average
// to shave the worst case") from the other direction -- here the worst
// case can only ever tie or improve, never trade away, so a real avg
// win (if found) can't be an illusion from the tuple ordering.
fn tune_basin_hop_avg_first(
    name: &str,
    f: &dyn Fn(f32, &[f32]) -> f32,
    reference: &dyn Fn(f64) -> f64,
    grid: &[f32],
    init: &[f32],
    restarts: usize,
) {
    fn better(candidate: (u64, u64), current: (u64, u64), cap: u64) -> bool {
        let cand_ok = candidate.0 <= cap;
        let cur_ok = current.0 <= cap;
        match (cand_ok, cur_ok) {
            (true, false) => true,
            (false, true) => false,
            _ => candidate.1 < current.1,
        }
    }

    fn descend_capped(
        f: &dyn Fn(f32, &[f32]) -> f32,
        reference: &dyn Fn(f64) -> f64,
        grid: &[f32],
        start: &[f32],
        cap: u64,
    ) -> (Vec<f32>, (u64, u64)) {
        let mut c = start.to_vec();
        let mut best = score(f, reference, grid, &c);
        let mut improved = true;
        while improved {
            improved = false;
            for i in 0..c.len() {
                for delta in [1i32, -1, 2, -2, 4, -4, 8, -8, 16, -16] {
                    let mut trial = c.clone();
                    trial[i] = f32::from_bits((trial[i].to_bits() as i32 + delta) as u32);
                    let s = score(f, reference, grid, &trial);
                    if better(s, best, cap) {
                        best = s;
                        c = trial;
                        improved = true;
                    }
                }
            }
        }
        (c, best)
    }

    let mut c: Vec<f32> = init.to_vec();
    let start = score(f, reference, grid, &c);
    let cap = start.0; // basin-hopping may never regress max beyond the un-hopped starting point
    let mut best = start;
    println!("{name}: start max {} avg {:.5}", best.0, best.1 as f64 / grid.len() as f64);

    use rand::RngExt;
    let mut rng = rand::rng();
    for _ in 0..restarts {
        let mut trial_start = c.clone();
        let n_perturb = 2 + (rng.random::<u8>() % 2) as usize;
        for _ in 0..n_perturb {
            let idx = (rng.random::<u32>() as usize) % trial_start.len();
            let delta = (rng.random::<i32>() % 64) - 32;
            trial_start[idx] = f32::from_bits((trial_start[idx].to_bits() as i32 + delta) as u32);
        }
        let (c_trial, s_trial) = descend_capped(f, reference, grid, &trial_start, cap);
        if better(s_trial, best, cap) {
            best = s_trial;
            c = c_trial;
        }
    }
    println!(
        "{name}: basin-hopped ({restarts} restarts, avg-first max-capped at {cap}) max {} avg {:.5}  coeffs: {:?}",
        best.0,
        best.1 as f64 / grid.len() as f64,
        c.iter().map(|v| format!("{v:e}")).collect::<Vec<_>>()
    );
}

// Like tune(), but never perturbs coefficient index 0 -- for cases where
// that coefficient is a mathematically-required exact value (e.g.
// log_2's leading term, exactly log2(e)), not a free empirical parameter.
fn tune_fixed0(
    name: &str,
    f: &dyn Fn(f32, &[f32]) -> f32,
    reference: &dyn Fn(f64) -> f64,
    grid: &[f32],
    init: &[f32],
) {
    let mut c: Vec<f32> = init.to_vec();
    let mut best = score(f, reference, grid, &c);
    println!("{name}: start max {} avg {:.5}", best.0, best.1 as f64 / grid.len() as f64);
    let mut improved = true;
    while improved {
        improved = false;
        for i in 1..c.len() {
            for delta in [1i32, -1, 2, -2, 4, -4, 8, -8, 16, -16] {
                let mut trial = c.clone();
                trial[i] = f32::from_bits((trial[i].to_bits() as i32 + delta) as u32);
                let s = score(f, reference, grid, &trial);
                if s < best {
                    best = s;
                    c = trial;
                    improved = true;
                }
            }
        }
    }
    println!(
        "{name}: tuned max {} avg {:.5}  coeffs: {:?}",
        best.0,
        best.1 as f64 / grid.len() as f64,
        c.iter().map(|v| format!("{v:e}")).collect::<Vec<_>>()
    );
}

// backlog idea #22: minimize a secondary function's error (e.g. log1p)
// via a shared coefficient set, subject to the primary function's own
// (e.g. ln's) max ulp never regressing past its already-shipped
// baseline -- same "constrained search" shape as acos_poly's own joint
// asin refit (fix 7), just generalized to take two arbitrary (f,
// reference, grid) triples instead of hardcoding acos/asin. Coefficient
// index 0 is never perturbed (mathematically-required exact leading
// term, same convention as tune_fixed0).
#[allow(clippy::too_many_arguments)]
fn tune_joint_fixed0(
    name: &str,
    primary_f: &dyn Fn(f32, &[f32]) -> f32,
    primary_ref: &dyn Fn(f64) -> f64,
    primary_grid: &[f32],
    secondary_f: &dyn Fn(f32, &[f32]) -> f32,
    secondary_ref: &dyn Fn(f64) -> f64,
    secondary_grid: &[f32],
    init: &[f32],
) {
    let mut c: Vec<f32> = init.to_vec();
    let primary_cap = score(primary_f, primary_ref, primary_grid, &c).0;
    let mut best = score(secondary_f, secondary_ref, secondary_grid, &c);
    println!(
        "{name}: start secondary max {} avg {:.5} (primary cap max {})",
        best.0,
        best.1 as f64 / secondary_grid.len() as f64,
        primary_cap
    );
    let mut improved = true;
    while improved {
        improved = false;
        for i in 1..c.len() {
            for delta in [1i32, -1, 2, -2, 4, -4, 8, -8, 16, -16] {
                let mut trial = c.clone();
                trial[i] = f32::from_bits((trial[i].to_bits() as i32 + delta) as u32);
                let primary_s = score(primary_f, primary_ref, primary_grid, &trial);
                if primary_s.0 > primary_cap {
                    continue; // would regress the primary function's own max
                }
                let s = score(secondary_f, secondary_ref, secondary_grid, &trial);
                if s < best {
                    best = s;
                    c = trial;
                    improved = true;
                }
            }
        }
    }
    let final_primary = score(primary_f, primary_ref, primary_grid, &c);
    println!(
        "{name}: tuned secondary max {} avg {:.5}  primary max {} avg {:.5}  coeffs: {:?}",
        best.0,
        best.1 as f64 / secondary_grid.len() as f64,
        final_primary.0,
        final_primary.1 as f64 / primary_grid.len() as f64,
        c.iter().map(|v| format!("{v:e}")).collect::<Vec<_>>()
    );
}

fn main() {
    let which = std::env::args().nth(1).unwrap_or_default();
    if which.contains("exp2") || which.is_empty() {
        // grid over (-126, 128)
        let mut grid = vec![];
        let mut b = 1e-6f32.to_bits();
        while b <= 126.0f32.to_bits() {
            grid.push(f32::from_bits(b));
            grid.push(-f32::from_bits(b));
            b += 997;
        }
        // current shipped coefficients (src/lib.rs's exp2), not the
        // pre-tuning starting point this init array used to be -- check
        // for headroom from where the crate actually is now.
        let init = [2.1702237e-4, 1.2439679e-3, 9.678826e-3, 5.548333e-2, 2.4022985e-1, 6.93147e-1];
        tune("exp2", &exp2_c, &|x| x.exp2(), &grid, &init);
    }
    if which.contains("exp2direct") {
        // same domain/grid as "exp2" above.
        let mut grid = vec![];
        let mut b = 1e-6f32.to_bits();
        while b <= 126.0f32.to_bits() {
            grid.push(f32::from_bits(b));
            grid.push(-f32::from_bits(b));
            b += 997;
        }
        // Taylor seed for P(f)=2^f=e^(f ln2), c[0]=1.0 exact (pinned, never
        // perturbed by tune_fixed0), c[1..6) = ln2^k/k! for k=1..5.
        let l = std::f64::consts::LN_2;
        let init = [
            1.0,
            l as f32,
            (l * l / 2.0) as f32,
            (l * l * l / 6.0) as f32,
            (l * l * l * l / 24.0) as f32,
            (l * l * l * l * l / 120.0) as f32,
        ];
        tune_fixed0("exp2_direct", &exp2_direct_c, &|x| x.exp2(), &grid, &init);
    }
    if which.contains("exp2round") {
        // same domain/grid as "exp2" above, just k=round(x) internally.
        let mut grid = vec![];
        let mut b = 1e-6f32.to_bits();
        while b <= 126.0f32.to_bits() {
            grid.push(f32::from_bits(b));
            grid.push(-f32::from_bits(b));
            b += 997;
        }
        // Taylor-series seed for Q(f)=(2^f-1)/f = sum_{k>=1} f^(k-1)*ln2^k/k!
        // (exact regardless of domain -- coordinate descent refines it into
        // a minimax fit for the centered [-0.5,0.5] domain specifically).
        let q0 = std::f64::consts::LN_2;
        let q1 = q0 * q0 / 2.0;
        let q2 = q0 * q0 * q0 / 6.0;
        let q3 = q2 * q0 / 4.0;
        let q4 = q3 * q0 / 5.0;
        let q5 = q4 * q0 / 6.0;
        let init5 = [q5 as f32, q4 as f32, q3 as f32, q2 as f32, q1 as f32, q0 as f32];
        tune("exp2_round (degree 5)", &exp2_round_c, &|x| x.exp2(), &grid, &init5);
        let init4 = [q4 as f32, q2 as f32, q3 as f32, q0 as f32, q1 as f32];
        tune("exp2_round4 (degree 4)", &exp2_round4_c, &|x| x.exp2(), &grid, &init4);
    }
    if which.contains("exp2lut") {
        let mut grid = vec![];
        let mut b = 1e-6f32.to_bits();
        while b <= 126.0f32.to_bits() {
            grid.push(f32::from_bits(b));
            grid.push(-f32::from_bits(b));
            b += 997;
        }
        // scipy-derived seed (least_squares fit of 2^f_lo directly over
        // f_lo in [0,1/8), not zero-seeded).
        let init = [1.0, 0.69315195, 0.2400347, 0.05795604];
        tune_fixed0("exp2_lut8", &exp2_lut8_c, &|x| x.exp2(), &grid, &init);
    }
    if which.contains("log2") || which.is_empty() {
        let mut grid = vec![];
        let mut b = 0x0080_0000u32;
        while b < 0x7f80_0000 {
            grid.push(f32::from_bits(b));
            b += 1499;
        }
        // current shipped coefficients (src/lib.rs's log_2_normal); c[0]
        // is std::f32::consts::LOG2_E exactly, not a coincidence (it's
        // the poly's mathematically-required leading term, log2(1+s)'s
        // derivative at s=0), so it's excluded from tuning below.
        let init = [
            std::f32::consts::LOG2_E, -0.72134733, 0.4808985, -0.36069715, 0.288568,
            -0.23961738, 0.20460059, -0.19106273, 0.18617496, -0.10994955,
        ];
        tune_fixed0("log2", &log2_c, &|x| x.log2(), &grid, &init);
    }
    if which.contains("log2atanh") {
        let mut grid = vec![];
        let mut b = 0x0080_0000u32;
        while b < 0x7f80_0000 {
            grid.push(f32::from_bits(b));
            b += 1499;
        }
        // scipy-derived seed (least_squares minimax fit of Q(u) against
        // atanh(t)/t, u=t^2, t=(m-1)/(m+1) over m in [sqrt2/2, sqrt2)) --
        // not zero-seeded, see this file's own zero-seed-trap lesson.
        // c[0] = 1.0 exactly (Q(0) = atanh'(0) = 1), excluded from tuning.
        let init = [1.0, 0.33333333, 0.20000162, 0.14269426, 0.11772234];
        tune_fixed0("log2_atanh", &log2_atanh_c, &|x| x.log2(), &grid, &init);
    }
    if which.contains("lnlog10") {
        let mut grid = vec![];
        let mut b = 0x0080_0000u32;
        while b < 0x7f80_0000 {
            grid.push(f32::from_bits(b));
            b += 1499;
        }
        // current shipped coefficients (src/lib.rs's ln_normal/
        // log10_normal): log_2's own c[i] * LN_2/LOG10_2, individually
        // rounded to f32 -- never refit directly against ln/log10 as
        // their own objective until now.
        let init = [
            1.0, -0.49999988, 0.33333343, -0.25001621, 0.20002009, -0.16609012, 0.14181833,
            -0.13243459, 0.12904665, -0.07621122,
        ];
        tune_fixed0("ln", &ln_poly_c, &|x| x.ln(), &grid, &init);
        let init = [
            0.4342945, -0.2171472, 0.14476489, -0.10858066, 0.08686763, -0.07213202, 0.06159092,
            -0.05751561, 0.05604425, -0.03309811,
        ];
        tune_fixed0("log10", &log10_poly_c, &|x| x.log10(), &grid, &init);
    }
    if which.contains("log1pjoint") {
        // backlog idea #22: found on 2026-07-09 that the coarse grid
        // here reports a modest, real-looking win (secondary/log1p max
        // ulp unchanged, avg -1.6%; primary/ln max even improves 3->2 as
        // a bonus) -- but wiring the resulting coefficients into
        // src/lib.rs's actual `ln_normal` and re-checking against the
        // real 100M-sample `accuracy.rs` fuzz (not this file's own
        // coarse grid) showed a genuine regression instead: log1p avg
        // ulp 0.0966->0.1085 (worse) and max ulp 4->7 (worse). Same
        // "always verify a tune.rs result against the real sweep before
        // trusting it" lesson this file's own `tune_basin_hop` doc
        // comment already documents for acos_poly -- the coarse grid
        // here (built from a ~613-step walk over log1p's own bit
        // patterns) apparently doesn't sample densely enough near
        // whatever region the real 100M-sample fuzz's own worst case
        // lives in. Reverted the coefficient change in src/lib.rs; kept
        // this tuning infrastructure (`log1p_via_ln_c`,
        // `tune_joint_fixed0`) since the *technique* (constrained joint
        // refit sharing one coefficient set across two call sites) is
        // reusable even though this specific run's result didn't survive
        // real verification -- a denser/differently-distributed grid
        // might do better for a future attempt.
        //
        // ln's own grid (the constraint: never let ln's own max regress
        // past its shipped baseline), same as the "lnlog10" block above.
        let mut ln_grid = vec![];
        let mut b = 0x0080_0000u32;
        while b < 0x7f80_0000 {
            ln_grid.push(f32::from_bits(b));
            b += 1499;
        }
        // log1p's own grid: x in (-1, large), built the same "iterate
        // raw bits, coarse step" way but over log1p's actual domain
        // (u=1+x staying positive) instead of ln's own x>0 domain --
        // log1p's whole point is precision for x near 0, so this needs
        // its own domain, not ln's.
        let mut log1p_grid = vec![];
        // x in (-1, -1e-4): bits just below 0xbf800000 (-1.0) up toward
        // -0.0001, walked by decreasing magnitude (increasing raw bits
        // of -x, since -x in (1e-4, 1) has the same bit-pattern-order
        // property positive floats do).
        let mut nb = (1e-4f32).to_bits();
        while nb < 0x3f800000 {
            // magnitude in (1e-4, 1.0)
            log1p_grid.push(-f32::from_bits(nb));
            nb += 613;
        }
        let mut b = (1e-4f32).to_bits();
        while b < 0x4c000000 {
            // x in (1e-4, 3.36e7)
            log1p_grid.push(f32::from_bits(b));
            b += 613;
        }
        // current shipped ln coefficients (src/lib.rs's log1p reuses
        // `ln` unmodified) -- same seed as the "lnlog10" block's own ln
        // init, since log1p_via_ln_c calls ln_poly_c directly.
        let init = [
            1.0, -0.49999988, 0.33333343, -0.25001621, 0.20002009, -0.16609012, 0.14181833,
            -0.13243459, 0.12904665, -0.07621122,
        ];
        tune_joint_fixed0(
            "log1p_joint",
            &ln_poly_c,
            &|x| x.ln(),
            &ln_grid,
            &log1p_via_ln_c,
            &|x| x.ln_1p(),
            &log1p_grid,
            &init,
        );
    }
    if which.contains("asin") {
        // asin's mid branch is only ever evaluated for a = |x| in
        // [0.1, 0.9) in the shipped code (see src/lib.rs's asin) -- tune
        // against exactly that range, not the whole [0,1] domain, so the
        // objective matches what's actually on the hot path.
        let mut grid = vec![];
        let mut b = 0.1f32.to_bits();
        while b < 0.9f32.to_bits() {
            grid.push(f32::from_bits(b));
            grid.push(-f32::from_bits(b));
            b += 37;
        }
        let init = [0.0392588, 0.179323, 1.75866, -3.66063];
        tune("asin_mid", &asin_mid_c, &|x| x.asin(), &grid, &init);
    }
    if which.contains("atan") {
        // atan_poly is only ever called on a.min(1/a), i.e. x in [0,1].
        let mut grid = vec![];
        let mut b = 1e-7f32.to_bits();
        while b <= 1.0f32.to_bits() {
            grid.push(f32::from_bits(b));
            b += 37;
        }
        let init = [0.040634338, 0.65748954, 0.17133473, 0.9907859];
        tune("atan_poly", &atan_poly_c, &|x| x.atan(), &grid, &init);
        // Zero-seeded: gets trapped exactly at the 2/2 form's own optimum.
        // 0.0's bit pattern is 0x0, so the tuner's integer-bit-step search
        // (deltas of 1..16 ulp) only reaches denormal-scale values from
        // there, which have ~zero effect on the polynomial and never score
        // better -- a "zero-move local optimum" that's really just the
        // search being unable to leave 0, not evidence of no headroom (an
        // arbitrary nonzero seed like 1e-3, tried and discarded, did
        // worse still: it starts from a point worse than the 2/2 form and
        // the greedy per-coordinate search never recovered, converging to
        // max ulp 565).
        let init7 = [0.0, 0.040634338, 0.65748954, 0.0, 0.17133473, 0.9907859];
        tune("atan_poly7 (3/3)", &atan_poly7_c, &|x| x.atan(), &grid, &init7);
        // What actually works: a real least-squares Pade fit (scipy, not
        // a guess) of the 3/3 shape against atan(x) directly over [0,1] --
        // found max abs error 1.3e-9 vs the 2/2 form's 7.8e-7 (~600x).
        // Fed as the coordinate-descent starting point instead of 0.0.
        let init7d = [0.00883003, 0.28497791, 1.12717105, 0.05016619, 0.57181574, 1.46050425];
        tune("atan_poly7 (3/3, scipy seed)", &atan_poly7_c, &|x| x.atan(), &grid, &init7d);

        // Three-interval reduction: scipy-derived degree-2/2 seed (fit
        // over u in [-tan(pi/8), tan(pi/8)] against atan(u) directly, not
        // zero-seeded).
        let init3 = [0.0616303, 0.75416583, 0.22413336, 1.08749911];
        tune("atan_three", &atan_three_c, &|x| x.atan(), &grid, &init3);

        // Pure-poly latency tier: scipy-derived degree-17 seed (least_squares
        // odd-poly fit against atan(x) directly over [0,1], not zero-seeded).
        let init_pure = [
            -0.33333168,
            0.19994266,
            -0.14216056,
            0.10689226,
            -0.07608681,
            0.04395558,
            -0.01687014,
            0.0030569030,
        ];
        tune("atan_pure_poly", &atan_pure_poly_c, &|x| x.atan(), &grid, &init_pure);
    }
    if which.contains("acos") {
        // acos/asin's near-1 branch both evaluate this for x = |input| in
        // [0,1] (acos directly; asin only above a > 0.9, but tuning the
        // whole [0,1] range keeps the poly consistent for both callers).
        let mut grid = vec![];
        let mut b = 0.0f32.to_bits();
        while b < 1.0f32.to_bits() {
            grid.push(f32::from_bits(b));
            b += 10000;
        }
        // c[6] = 1.5707964, matching FRAC_PI_2 exactly (idea #36's fix,
        // 2026-07-10 -- this init array previously held the stale
        // pre-fix 1.5707963, an 8-significant-digit transcription that
        // rounds 1 ulp short of the true nearest f32).
        let init = [2.2960256e-3, -1.1146317e-2, 2.6900213e-2, -4.8802543e-2, 8.8755615e-2, -2.1458544e-1, 1.5707964];
        tune("acos_poly", &acos_poly_c, &|x| x.acos(), &grid, &init);
        // Coarser grid (100x fewer points) specifically for the basin-hop
        // phase, so each of many restarts' full re-descent stays fast --
        // the single-axis result above already confirms the full grid's
        // own local optimum, this just needs a representative proxy to
        // screen many candidate basins quickly.
        let mut coarse_grid = vec![];
        let mut b = 0.0f32.to_bits();
        while b < 1.0f32.to_bits() {
            coarse_grid.push(f32::from_bits(b));
            b += 1_000_000;
        }
        tune_basin_hop("acos_poly_bh", &acos_poly_c, &|x| x.acos(), &coarse_grid, &init, 50);
        // idea #101: same coarse grid and restart count, but avg-first
        // max-capped instead of the default max-first tuple ordering --
        // re-checks whether the documented bias (not the annealing
        // mechanism itself) was really the reason the plain basin-hop
        // result reversed on the real fuzz.
        tune_basin_hop_avg_first("acos_poly_bh_avgfirst", &acos_poly_c, &|x| x.acos(), &coarse_grid, &init, 50);
    }
    if which.contains("asinpoly") {
        // Decoupled asin-only fit (backlog idea #34) -- grid restricted to
        // asin's own actual domain [0.25, 1.0) for this branch, seeded
        // from a real scipy least-squares fit (not 0.0 or an arbitrary
        // constant -- this file's own "zero-move trap" precedent) of the
        // same 7-coefficient shape against acos(x)/sqrt(1-x) weighted by
        // 1/target over exactly this domain.
        let mut grid = vec![];
        let mut b = 0.25f32.to_bits();
        while b < 1.0f32.to_bits() {
            grid.push(f32::from_bits(b));
            b += 4000;
        }
        // Exact current shipped coefficients (2026-07-10 -- the previous
        // seed here was a slightly-imprecise transcription, off by a few
        // ULP in c[0]; not a real bug like this session's other stale
        // seeds, but coordinate descent was spending its first move just
        // recovering the true value instead of starting there).
        let seed = [1.3137129e-3, -7.664533e-3, 2.2003133e-2, -4.5330178e-2, 8.745548e-2, -2.1434246e-1, 1.5707785];
        tune("asin_poly", &asin_poly_c, &|x| x.asin(), &grid, &seed);
    }
    if which.contains("asinacos") {
        // asin now calls acos_poly_c directly for a >= 0.25 (pi/2 -
        // acos_poly_c(a), sign-restored) instead of a separate fit -- but
        // acos_poly_c was only ever tuned against *acos*'s own ulp error.
        // acos(x) shrinks to 0 as x->1 while asin(x) grows to pi/2 there,
        // so the same absolute poly error gets a very different *relative*
        // (ulp) weight depending on which caller's output magnitude it's
        // measured against -- a poly tuned purely for acos's hardest
        // region (x->1, where acos's own value is tiny) may be leaving
        // accuracy on the table specifically for asin's use. Score jointly
        // (worst of the two callers' ulp error at each point) instead.
        let mut grid = vec![];
        let mut b = 0.0f32.to_bits();
        while b < 1.0f32.to_bits() {
            grid.push(f32::from_bits(b));
            b += 10000;
        }
        // c[6] = 1.5707964, matching FRAC_PI_2 exactly (idea #36's fix,
        // 2026-07-10 -- this init array previously held the stale
        // pre-fix 1.5707963, an 8-significant-digit transcription that
        // rounds 1 ulp short of the true nearest f32).
        let init = [2.2960256e-3, -1.1146317e-2, 2.6900213e-2, -4.8802543e-2, 8.8755615e-2, -2.1458544e-1, 1.5707964];
        let joint_score = |c: &[f32]| -> (u64, u64) {
            let mut sum = 0u64;
            let mut max = 0u64;
            for &x in &grid {
                let got = acos_poly_c(x, c);
                let d_acos = ulp_diff(got, x.acos() as f32);
                let d_asin = if x >= 0.25 {
                    let asin_got = std::f32::consts::FRAC_PI_2 - got;
                    ulp_diff(asin_got, x.asin() as f32)
                } else {
                    0
                };
                let d = d_acos.max(d_asin);
                sum += d;
                max = max.max(d);
            }
            (max, sum)
        };
        let mut c: Vec<f32> = init.to_vec();
        let mut best = joint_score(&c);
        println!("acos_asin_joint: start max {} avg {:.5}", best.0, best.1 as f64 / grid.len() as f64);
        let mut improved = true;
        while improved {
            improved = false;
            for i in 0..c.len() {
                for delta in [1i32, -1, 2, -2, 4, -4, 8, -8, 16, -16] {
                    let mut trial = c.clone();
                    trial[i] = f32::from_bits((trial[i].to_bits() as i32 + delta) as u32);
                    let s = joint_score(&trial);
                    if s < best {
                        best = s;
                        c = trial;
                        improved = true;
                    }
                }
            }
        }
        println!(
            "acos_asin_joint: tuned max {} avg {:.5}  coeffs: {:?}",
            best.0,
            best.1 as f64 / grid.len() as f64,
            c.iter().map(|v| format!("{v:e}")).collect::<Vec<_>>()
        );
        // Constrained variant: minimize asin's own error while requiring
        // acos's own max ulp never exceeds its current shipped best (4 on
        // this grid) -- checks whether asin can improve *without* costing
        // acos anything, rather than accepting a cross-function tradeoff.
        let acos_cap = {
            let mut m = 0u64;
            for &x in &grid {
                let got = acos_poly_c(x, &init);
                m = m.max(ulp_diff(got, x.acos() as f32));
            }
            m
        };
        let asin_only_score = |c: &[f32]| -> Option<(u64, u64)> {
            let mut sum = 0u64;
            let mut max = 0u64;
            let mut acos_max = 0u64;
            for &x in &grid {
                let got = acos_poly_c(x, c);
                acos_max = acos_max.max(ulp_diff(got, x.acos() as f32));
                if x >= 0.25 {
                    let asin_got = std::f32::consts::FRAC_PI_2 - got;
                    let d = ulp_diff(asin_got, x.asin() as f32);
                    sum += d;
                    max = max.max(d);
                }
            }
            if acos_max > acos_cap { None } else { Some((max, sum)) }
        };
        let mut c: Vec<f32> = init.to_vec();
        let mut best = asin_only_score(&c).expect("init must satisfy its own cap");
        println!("acos_capped_asin: start max {} avg {:.5} (acos cap {})", best.0, best.1 as f64 / grid.len() as f64, acos_cap);
        let mut improved = true;
        while improved {
            improved = false;
            for i in 0..c.len() {
                for delta in [1i32, -1, 2, -2, 4, -4, 8, -8, 16, -16] {
                    let mut trial = c.clone();
                    trial[i] = f32::from_bits((trial[i].to_bits() as i32 + delta) as u32);
                    if let Some(s) = asin_only_score(&trial) {
                        if s < best {
                            best = s;
                            c = trial;
                            improved = true;
                        }
                    }
                }
            }
        }
        println!(
            "acos_capped_asin: tuned max {} avg {:.5}  coeffs: {:?}",
            best.0,
            best.1 as f64 / grid.len() as f64,
            c.iter().map(|v| format!("{v:e}")).collect::<Vec<_>>()
        );
    }
    if which.contains("acos8") {
        // Same joint (worst-of-acos/asin) scoring as "asinacos" above, but
        // for the degree-7 (8-coefficient) acos_poly8_c -- one extra
        // coefficient prepended (0.0) as the starting point, so coordinate
        // descent starts from a value bit-identical to the shipped degree-6
        // poly and only has to find where the new degree helps.
        let mut grid = vec![];
        let mut b = 0.0f32.to_bits();
        while b < 1.0f32.to_bits() {
            grid.push(f32::from_bits(b));
            b += 10000;
        }
        let init = [
            0.0,
            2.2960256e-3,
            -1.1146317e-2,
            2.6900213e-2,
            -4.8802543e-2,
            8.8755615e-2,
            -2.1458544e-1,
            1.5707964,
        ];
        let joint_score = |c: &[f32]| -> (u64, u64) {
            let mut sum = 0u64;
            let mut max = 0u64;
            for &x in &grid {
                let got = acos_poly8_c(x, c);
                let d_acos = ulp_diff(got, x.acos() as f32);
                let d_asin = if x >= 0.25 {
                    let asin_got = std::f32::consts::FRAC_PI_2 - got;
                    ulp_diff(asin_got, x.asin() as f32)
                } else {
                    0
                };
                let d = d_acos.max(d_asin);
                sum += d;
                max = max.max(d);
            }
            (max, sum)
        };
        let mut c: Vec<f32> = init.to_vec();
        let mut best = joint_score(&c);
        println!("acos8_asin_joint: start max {} avg {:.5}", best.0, best.1 as f64 / grid.len() as f64);
        let mut improved = true;
        while improved {
            improved = false;
            for i in 0..c.len() {
                for delta in [1i32, -1, 2, -2, 4, -4, 8, -8, 16, -16] {
                    let mut trial = c.clone();
                    trial[i] = f32::from_bits((trial[i].to_bits() as i32 + delta) as u32);
                    let s = joint_score(&trial);
                    if s < best {
                        best = s;
                        c = trial;
                        improved = true;
                    }
                }
            }
        }
        println!(
            "acos8_asin_joint: tuned max {} avg {:.5}  coeffs: {:?}",
            best.0,
            best.1 as f64 / grid.len() as f64,
            c.iter().map(|v| format!("{v:e}")).collect::<Vec<_>>()
        );
    }
    if which.contains("erf") {
        // erf's tail branch is only ever used for xa in [0.28, 10] (see
        // erf's doc comment).
        let mut grid = vec![];
        let mut b = 0.28f32.to_bits();
        while b < 10.0f32.to_bits() {
            grid.push(f32::from_bits(b));
            grid.push(-f32::from_bits(b));
            b += 3000;
        }
        let init = [3.118769e-4, -4.67225e-3, 3.3162573e-2, -1.5214339e-1, -9.1684705e-1, -1.6282598, 3.1332566e-5];
        tune("erf_tail", &erf_tail_c, &erf_ref, &grid, &init);

        // erf's near-zero Pade branch, |x| < 0.28.
        let mut grid = vec![];
        let mut b = 0f32.to_bits();
        while b < 0.28f32.to_bits() {
            grid.push(f32::from_bits(b));
            grid.push(-f32::from_bits(b));
            b += 2000;
        }
        let init = [0.5910557508468628, 1.128379225730896, 0.18571428954601288, 0.8571428656578064];
        tune("erf_near0", &erf_near0_c, &erf_ref, &grid, &init);
    }
    if which == "erfcx" {
        // The fit variable is `v = 1/(2+xa)`, not `xa`, so the grid is
        // built by sweeping `v` evenly over its whole range (0, 1/2] and
        // mapping back -- an even sweep of `xa` would spend almost every
        // sample in a sliver of `v` and leave the tail (small `v`, where
        // the asymptote lives) unrepresented. `v` runs from 1/2 (xa = 0)
        // down toward 0 (xa -> inf); the smallest `v` here corresponds to
        // xa ~ 1e5, far past where the fit has anything left to resolve.
        let mut grid = vec![];
        for i in 1..=20000u32 {
            let v = 0.5 * (i as f32) / 20000.0;
            grid.push(1.0 / v - 2.0);
        }
        // erfcx(xa) = exp(xa^2)*erfc(xa) in f64 while that is safe, and
        // the standard asymptotic series past it (`exp(xa^2)` overflows
        // f64 at xa ~ 26.6). The series' first dropped term is
        // 945/(2*xa^2)^5, ~2.9e-12 relative at xa = 20 -- five orders
        // under f32 resolution, so the seam is invisible at this scale.
        // Same reference accuracy.rs's own `erfcx (x>=20)` row uses.
        let erfcx_ref = |xa: f64| {
            if xa < 20.0 {
                (xa * xa).exp() * erfc_ref(xa)
            } else {
                let t = 0.5 / (xa * xa);
                let p = 1.0 - t * (1.0 - t * (3.0 - t * (15.0 - t * 105.0)));
                p / (xa * std::f64::consts::PI.sqrt())
            }
        };
        // Seed: the shipped coefficients (a Chebyshev-basis relative-error
        // LP over v, already coordinate-descent polished here). Descending
        // from here is a no-op unless the grid or evaluation order
        // changed -- which is exactly what makes it a useful regression
        // check on a future reshape.
        //
        // The one property descent cannot see and cannot be given as an
        // objective: the shipped coefficients also satisfy
        // `erfcx_pos(0.0) == 1.0` bit-exactly (pinned in edgecheck.rs, and
        // the reason `erfc(0)` and `erfcx(0)` are exact). Re-check any
        // candidate this prints against that pin before shipping.
        let init = [
            0.56418955, 1.1283774, 1.9749641, 2.807048, 2.975676, -3.7488432, 17.02367,
            -117.490135, 255.59447, -243.95302, 90.238014,
        ];
        tune("erfcx_pos", &erfcx_pos_c, &erfcx_ref, &grid, &init);
    }
    if which.contains("cbrt") {
        // one octave [1,2) is representative: the bit-trick seed's
        // relative error pattern repeats across octaves (see
        // cbrt_normal's doc comment -- fitted against the seed's error
        // range, not a specific magnitude range).
        let mut grid = vec![];
        let mut b = 1.0f32.to_bits();
        while b < 2.0f32.to_bits() {
            grid.push(f32::from_bits(b));
            grid.push(-f32::from_bits(b));
            b += 5;
        }
        // Current shipped coefficients (2026-07-10 -- this init array
        // previously held a stale pre-2026-07-09-refit value, the same
        // class of staleness already found for exp_pos_neg/erf_tail_c/
        // acos_poly's own seeds).
        let init = [-0.3333314061164856, 0.22221335768699646, -0.1739402711391449, 0.14720453321933746];
        tune("cbrt_normal", &cbrt_normal_c, &|x| x.cbrt(), &grid, &init);

        // seed offset (c[0], bits used directly) jointly with the same
        // degree-3 poly, seeded from the shipped values (not zero -- see
        // this file's own zero-seed-trap lesson).
        let mut grid2 = vec![];
        let mut b = 1.0f32.to_bits();
        while b < 2.0f32.to_bits() {
            grid2.push(f32::from_bits(b));
            grid2.push(-f32::from_bits(b));
            b += 5;
        }
        let init_joint = [
            f32::from_bits(0x2a509a07),
            -0.33333147,
            0.22220612,
            -0.17394388,
            0.14823665,
        ];
        tune("cbrt_normal_joint", &cbrt_normal_joint_c, &|x| x.cbrt(), &grid2, &init_joint);

        // idea #52: degree-4 correction poly (5 coeffs), seeded from a
        // real scipy/HiGHS minimax LP fit of p(r) = ((1+r)^(-1/3)-1)/r
        // against the seed's real r-range (probed directly: [-0.0998,
        // 0.0893] over one octave), not zero-seeded.
        let mut grid5 = vec![];
        let mut b = 1.0f32.to_bits();
        while b < 2.0f32.to_bits() {
            grid5.push(f32::from_bits(b));
            grid5.push(-f32::from_bits(b));
            b += 5;
        }
        let init5 = [-0.33333338, 0.22221859, -0.17280712, 0.14543792, -0.12951847];
        tune("cbrt_normal_c5", &cbrt_normal_c5, &|x| x.cbrt(), &grid5, &init5);
    }
    if which.contains("cbrtshift") {
        // coarser grid for a fast first-pass screen of the shift-multiply
        // seed variant -- refine with the dense grid only if promising.
        let mut grid = vec![];
        let mut b = 1.0f32.to_bits();
        while b < 2.0f32.to_bits() {
            grid.push(f32::from_bits(b));
            grid.push(-f32::from_bits(b));
            b += 500;
        }
        let init = [-0.33333147, 0.22220612, -0.17394388, 0.14823665];
        tune("cbrt_shiftmul", &cbrt_shiftmul_c, &|x| x.cbrt(), &grid, &init);
    }
    if which.contains("cbrtthroughput") {
        // Unlike cbrt_normal, cbrt_throughput's error does NOT repeat
        // across octaves (checked directly -- max ulp ranges from ~7 to
        // ~52 depending which octave is sampled), so a single-octave grid
        // isn't representative here. Sample log-uniformly across ~60
        // octaves instead (positive only, same as above).
        let mut grid = vec![];
        let mut e = -30i32;
        while e < 30 {
            let mut b = (2.0f32.powi(e)).to_bits();
            let bmax = (2.0f32.powi(e + 1)).to_bits();
            while b < bmax {
                grid.push(f32::from_bits(b));
                b += (bmax - (2.0f32.powi(e)).to_bits()) / 20;
            }
            e += 1;
        }
        let init = [
            f32::from_bits(0xd461ff81),
            f32::from_bits(0x3fb6e3d7),
            f32::from_bits(0x3fe09c2a),
        ];
        tune("cbrt_throughput", &cbrt_throughput_c, &|x| x.cbrt(), &grid, &init);
    }
    if which.contains("sinf") {
        // sinf_poly's fitted domain, [-pi/2, pi/2].
        let mut grid = vec![];
        let mut b = 0f32.to_bits();
        let hpi = std::f32::consts::FRAC_PI_2.to_bits();
        while b < hpi {
            grid.push(f32::from_bits(b));
            grid.push(-f32::from_bits(b));
            b += 300;
        }
        let init = [-0.16666660, 8.3330662e-3, -1.9809603e-4, 2.6057806e-6];
        tune("sinf_poly", &sinf_poly_c, &|x| x.sin(), &grid, &init);
    }
    if which.contains("sinpipoly") {
        // idea #43: sinpi/cospi's own reduction always lands r in
        // [-0.5, 0.5] (see sinpi's doc comment) -- seeded from a real
        // scipy/HiGHS minimax LP fit of sin(pi*r) directly in r (not a
        // continuous refit of the generic sinf_poly), not zero-seeded.
        let mut grid = vec![];
        let mut b = 0f32.to_bits();
        let half = 0.5f32.to_bits();
        while b < half {
            grid.push(f32::from_bits(b));
            grid.push(-f32::from_bits(b));
            b += 150;
        }
        let init = [-5.167724609375, 2.550321102142334, -0.5996360778808594, 0.0800275132060051];
        tune("sinpi_poly", &sinpi_poly_c, &|r| (std::f64::consts::PI * r).sin(), &grid, &init);
    }
    if which.contains("sindpoly") {
        // idea #44: sind/cosd's own reduction always lands d in [-90,90]
        // (see sind's doc comment) -- seeded from a real scipy/HiGHS
        // minimax LP fit of sin(d*pi/180) directly in d, column-scaled
        // for HiGHS conditioning (d^9 up to ~3.9e17 unscaled), not
        // zero-seeded.
        let mut grid = vec![];
        let mut b = 0f32.to_bits();
        let ninety = 90.0f32.to_bits();
        while b < ninety {
            grid.push(f32::from_bits(b));
            grid.push(-f32::from_bits(b));
            b += 27000;
        }
        let init = [-8.860945968081069e-7, 1.3494975102668061e-11, -9.762676095181961e-17, 3.8616159998020477e-22];
        tune("sind_poly", &sind_poly_c, &|d| (std::f64::consts::PI / 180.0 * d).sin(), &grid, &init);
    }
    if which.contains("expm1") {
        // expm1's Pade branch domain, |x| < 0.5.
        let mut grid = vec![];
        let mut b = 0f32.to_bits();
        while b < 0.5f32.to_bits() {
            grid.push(f32::from_bits(b));
            grid.push(-f32::from_bits(b));
            b += 300;
        }
        let init = [-1.9999927, -120.0, -12.000030, 59.999996, -120.0];
        tune("expm1_near0", &expm1_near0_c, &|x| x.exp_m1(), &grid, &init);
        // degree-5-numerator bump: c[0] (the new x^5 term) starts at 0.0
        // so this starts bit-identical to the shipped degree-3 form.
        let init = [0.0, -1.9999927, -120.0, -12.000030, 59.999996, -120.0];
        tune("expm1_near0_deg5", &expm1_near0_deg5_c, &|x| x.exp_m1(), &grid, &init);
    }
    if which.contains("exp_r") {
        // exp's reduced-argument domain, r in [-ln2/2, ln2/2].
        let bound = std::f32::consts::LN_2 / 2.0;
        let mut grid = vec![];
        let mut b = 0f32.to_bits();
        while b < bound.to_bits() {
            grid.push(f32::from_bits(b));
            grid.push(-f32::from_bits(b));
            b += 42000;
        }
        // Current shipped exp coefficients (2026-07-10 -- this init array
        // previously held the pre-2026-07-09-LP-refit coordinate-descent
        // values exp's own doc comment says were superseded).
        let init = [4.9999300e-1, 1.6667245e-1, 4.1883811e-2, 8.3009899e-3];
        tune("exp_r (c0=c1=1 forced)", &exp_r_c, &|x| x.exp(), &grid, &init);
        // idea #25: degree 5->6 bump, seeded from a real scipy/HiGHS
        // minimax LP fit (idealized margin ~33x over an f32-emulated
        // reconstruction of the shipped degree-5 form -- much stronger
        // than the exp_pos_neg degree-bump screen that killed idea #26),
        // not zero-seeded.
        let init6 = [
            0.5000000596046448,
            0.16666492819786072,
            0.041665248572826385,
            0.00837147980928421,
            0.0013992299791425467,
        ];
        tune("exp_r6 (c0=c1=1 forced)", &exp_r_c6, &|x| x.exp(), &grid, &init6);
        // current shipped exp_pos_neg coefficients (src/lib.rs), not the
        // pre-retuning starting point above -- check for headroom from
        // where the crate actually is now (idea #7's sinh/cosh round-off
        // audit found the poly's own fit error dominates the real max-ulp-5
        // ceiling, so worth checking whether coordinate descent can do
        // better than what's shipped before assuming it's already optimal).
        let shipped = [4.999897e-1, 1.6666329e-1, 4.1917525e-2, 8.3811125e-3];
        tune("exp_r_pair (even/odd split)", &exp_r_pair_c, &|x| x.exp(), &grid, &init);
        tune("exp_r_pair (from shipped)", &exp_r_pair_c, &|x| x.exp(), &grid, &shipped);
        tune_basin_hop(
            "exp_r_pair (basin hop, from shipped)",
            &exp_r_pair_c,
            &|x| x.exp(),
            &grid,
            &shipped,
            200,
        );
    }
}
