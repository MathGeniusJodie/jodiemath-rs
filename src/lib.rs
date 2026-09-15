// godbolt flags -C opt-level=3 -C target_feature=+fma

// This crate's entire accuracy and performance story depends on `fma` (see the
// function just below) compiling to a single hardware instruction. Without FMA,
// `f32::mul_add` falls back to a ~2x-slower libm call that also rounds
// *differently* (two roundings instead of one) -- every ulp figure in this
// crate's doc comments, readme.md, and examples/accuracy.rs assumes the
// single-rounding hardware form.
#[cfg(all(not(target_feature = "fma"), not(doctest)))]
compile_error!(
    "jodiemath-rs requires hardware FMA (target-feature=+fma or target-cpu=native) -- \
     without it, f32::mul_add falls back to a slower, differently-rounded software path \
     and every accuracy/perf figure in this crate's docs is invalid. Build with \
     `RUSTFLAGS=\"-C target-cpu=native\"` or ensure .cargo/config.toml's rustflags \
     aren't being overridden by an environment RUSTFLAGS variable."
);

mod doublefloat;
mod pitable;
use doublefloat::Df32;

const SIGN_MASK: u32 = 0x80000000;
const EXPONENT_MASK: u32 = 0x7f800000;

#[inline(always)]
fn fma(a: f32, b: f32, c: f32) -> f32 {
    a.mul_add(b, c)
}

// `2^k` as a bare exponent field, for an integer-valued `k` in `[-127, 128]` --
// the range each caller's own domain contract already guarantees. `k = 128`
// deliberately yields `+inf` (field 255, mantissa 0), which is how the
// unchecked exp2 family overflows, and `k = -127` yields `+0.0`.
const EXP2INT_MAGIC: f32 = 12583039.0; // 1.5 * 2^23 + 127
macro_rules! exp2int_field {
    ($k:expr) => {
        f32::from_bits(($k + EXP2INT_MAGIC).to_bits() << 23)
    };
}

/// Round to the nearest integer (ties-to-even) for `|x| <= 2^22`.
/// Uses the magic constant `1.5 * 2^23`. Preserves the sign of zero for `x` in `[-0.5, 0)`.
#[inline(always)]
pub fn fast_round_int(x: f32) -> f32 {
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    (fma(x, 1.0, ROUND_MAGIC) - ROUND_MAGIC).copysign(x)
}

macro_rules! denormal_rescale {
    ($x:expr) => {{
        let tiny = $x < f32::MIN_POSITIVE;
        let xs = if tiny { $x * 16777216.0 } else { $x };
        let koff = if tiny { -24.0 } else { 0.0 };
        (xs, koff)
    }};
}

// exp(r) for tiny r, shared by exp/exp_checked/expm1/exp_m1_over_x/tanh/
// sigmoid (c0/c1 pinned to exactly 1.0 -- see exp's own body comment). A macro
// is pure textual substitution with no boundary at all -- verified equivalent
// via a full pre/post assembly diff.
macro_rules! exp_r_poly {
    ($r:expr) => {{
        let c: [f32; 4] = [4.9999300e-1, 1.6667245e-1, 4.1883811e-2, 8.3009899e-3];
        let r2 = $r * $r;
        let l0 = $r + 1.0;
        let l1 = fma(c[1], $r, c[0]);
        let l2 = fma(c[3], $r, c[2]);
        let m = fma(l2, r2, l1);
        fma(m, r2, l0)
    }};
}

// The `k`/`r` Cody-Waite reduction, e^r poly and 2^k reconstruction that
// `exp_checked` is, minus the input clamp -- the caller supplies an argument
// already inside `EXP_CLAMP_LO..=EXP_CLAMP_HI`. Split out for callers that can
// *prove* one side of that clamp is unreachable and so should not pay for it:
// `erfc`'s argument is `-(xs*xs)`, a negated square, so the upper bound is dead
// by construction (see its own comment).
macro_rules! exp_reduce {
    ($x:expr) => {{
        let x = $x;
        const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
        let k = fma(x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
        let r = fma(-k, LN2_HI, x);
        let r = fma(-k, LN2_LO, r);
        // e^r - 1 = r + r^2*(c0 + c1 r + c2 r^2 + c3 r^3 + c4 r^4), with the
        // top coefficient folded in at the `r^2` level rather than its own
        // `r^4` one, so `r^4` is never formed (`exp_r_poly!`'s trick). That
        // fold is why the extra degree costs no arithmetic; it pays for it by
        // serialising the top pair behind the bottom one, which is one more
        // level of fma latency.
        let c: [f32; 5] = [
            0.50000006,
            0.16666451,
            0.041665636,
            0.0083748708,
            0.0013946877,
        ];
        let r2 = r * r;
        let l1 = fma(c[1], r, c[0]);
        let l2 = fma(c[3], r, c[2]);
        let l3 = fma(c[4], r2, l2);
        let m = fma(l3, r2, l1);
        let s = fma(m, r2, r);
        let (t1, t2) = exp2_field_split(k);
        fma(s, t1, t1) * t2
    }};
}

/// `exp_checked`'s clamp bounds: `exp2_checked`'s own `k` boundary (`[-151,
/// 128)`) converted into `x`'s units.
const EXP_CLAMP_LO: f32 = -104.66522426455174;
const EXP_CLAMP_HI: f32 = 88.72283911167308;

// Q(f) = (2^f - 1)/f, shared by exp2/exp2_checked/exp10/exp10_checked/ exp2m1.
// Returns `q`; each caller does its own final combine (exp2/exp10's
// single-field `fma(q, exp2int*f, exp2int)` vs.
macro_rules! exp2_q_poly {
    ($f:expr) => {{
        let f2 = $f * $f;
        let g0 = fma(2.4022985e-1, $f, 6.93147e-1);
        let g1 = fma(9.678817e-3, $f, 5.548333e-2);
        let g2 = fma(2.1702255e-4, $f, 1.2439643e-3);
        let h = fma(g2, f2, g1);
        fma(h, f2, g0)
    }};
}

// `exp2_q_poly!`'s centered sibling: the same `Q(f) = (2^f - 1)/f` and the same
// 3-balanced-pair Estrin shape, refit for `f in [-0.5, 0.5]` instead of `[0,
// 1)`. Only `exp10_checked` uses it, to keep `round`'s own residual rather than
// paying a floor-adjust to convert it.
macro_rules! exp2_q_poly_centered {
    ($f:expr) => {{
        let f2 = $f * $f;
        let g0 = fma(2.402265e-1, $f, 6.9314719e-1);
        let g1 = fma(9.6182375e-3, $f, 5.5503574e-2);
        let g2 = fma(1.5403504e-4, $f, 1.3390731e-3);
        let h = fma(g2, f2, g1);
        fma(h, f2, g0)
    }};
}

// Shared by exp_pos_neg/exp_pos_neg_checked_half (sinh/cosh's unchecked/checked
// exp(x)/exp(-x) core): one Cody-Waite reduction, an even/odd-split poly, and
// the t1n/t2n reciprocal construction. Only each caller's optional input clamp
// differs, so that stays at the call site.
macro_rules! exp_pos_neg_core {
    ($x:expr) => {{
        const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
        let k = fma($x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
        let r = fma(-k, LN2_HI, $x);
        let r = fma(-k, LN2_LO, r);
        // The fitted coefficients, each exactly halved (see this macro's own
        // comment). Written as `* 0.5` rather than as pre-divided decimal
        // literals so the fit's own values stay legible and the halving cannot
        // be mistranscribed.
        let c: [f32; 5] = [
            0.49999994 * 0.5,
            0.16666521 * 0.5,
            0.041668329 * 0.5,
            8.3687045e-3 * 0.5,
            1.3814511e-3 * 0.5,
        ];
        // Each half is Horner in `r^2`, not a leading `c*r^4` term: `r^4`
        // is never formed, one plain multiply cheaper across both halves
        // at the same fma critical-path depth. Same fold as exp_r_poly!.
        let r2 = r * r;
        let e = fma(fma(fma(c[4], r2, c[2]), r2, c[0]), r2, 0.5);
        let o = fma(fma(c[3], r2, c[1]), r2, 0.5);
        let p_pos = fma(r, o, e);
        let p_neg = fma(-r, o, e);
        let (t1, t2) = exp2_field_split(k);
        // t1n = 1/t1, t2n = 1/t2: both are exact power-of-two fields, and for a
        // power-of-two float with bit pattern b = (127+e)<<23 the reciprocal
        // 2^-e has bit pattern (127-e)<<23 = 0x7F000000 - b. Exactly what
        // exp2_field_split(-k) would produce (round-half-to- even is
        // antisymmetric under negation), but built with two integer subtracts
        // instead of a second magic-round chain.
        let t1n = f32::from_bits(0x7F00_0000u32.wrapping_sub(t1.to_bits()));
        let t2n = f32::from_bits(0x7F00_0000u32.wrapping_sub(t2.to_bits()));
        (p_pos, p_neg, t1, t2, t1n, t2n)
    }};
}

// Shared by log_2/ln/log10: the denormal-rescale + special-case-select wrapper
// around each caller's own `_normal` fn. Edge handling uses selects (no early
// returns) so array loops auto-vectorize.
macro_rules! log_family_edges {
    ($x:expr, $r:expr) => {{
        let r = $r;
        let spec = if $x == 0.0 {
            f32::NEG_INFINITY
        } else {
            f32::NAN
        };
        let r = if $x <= 0.0 { spec } else { r };
        if !($x < f32::INFINITY) {
            $x * $x
        } else {
            r
        }
    }};
}

macro_rules! log_family_wrapper {
    ($x:expr, $normal:ident) => {
        log_family_edges!($x, {
            let (xs, koff) = denormal_rescale!($x);
            $normal(xs, koff)
        })
    };
}

// `log_family_wrapper!` for callers whose argument provably can't be a positive
// denormal, so the rescale's compare/multiply/two selects are dead code: `u =
// 1.0 + x` is such an argument for *every* f32 `x`, since `1+x` is exact by
// Sterbenz once `x <= -0.5` (making the smallest positive `u` exactly `2^-24`,
// ~10^30 above `f32::MIN_POSITIVE`) and is `>= 0.5` otherwise -- verified
// exhaustively over all 2^32 patterns, not argued from the bound alone. The
// zero/negative/inf/NaN arms are all still reachable at those call sites (`u ==
// 0` at `x == -1`, `u < 0` below it) and are shared verbatim.
macro_rules! log_family_wrapper_no_denormal {
    ($x:expr, $normal:ident) => {
        log_family_edges!($x, $normal($x, 0.0))
    };
}

// `log_family_wrapper!` for callers whose own outer select *discards* this
// value entirely unless the argument is positive and normal, so only
// `log_family_edges!`'s inf/NaN arm survives -- the zero/negative selects and
// the denormal rescale all compute results nothing can observe. Note what the
// licence is and isn't: every arm here is still *evaluated* for every input
// (branchless), so this is sound only when the caller's select provably drops
// the result for zero/negative/denormal arguments.
macro_rules! log_family_wrapper_discarded_unless_normal {
    ($x:expr, $normal:ident) => {
        if !($x < f32::INFINITY) {
            $x * $x
        } else {
            $normal($x, 0.0)
        }
    };
}

#[doc(alias = "log2f")]
#[doc(alias = "log2")]
#[inline(always)]
#[allow(clippy::neg_cmp_op_on_partial_ord)]
pub fn log_2(x: f32) -> f32 {
    log_family_wrapper!(x, log_2_normal)
}

/// Core of log_2 for positive normal finite x only: no handling for zero,
/// negative, denormal, inf, or nan input (those are the caller's job, see
/// log_2). Called directly with an out-of-domain x, this returns a
/// plausible-looking but wrong finite value rather than NaN/-inf.
const LOG2_Q_COEFFS: [f32; 9] = [
    -0.7213475,
    0.48089963,
    -0.36067435,
    0.28850868,
    -0.24009936,
    0.2058956,
    -0.18871288,
    0.17711402,
    -0.10358754,
];

#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn log_2_normal(x: f32, koff: f32) -> f32 {
    // decompose x = 2^k * m with m in [sqrt(2)/2, sqrt(2)), so s = m - 1
    // is exact (Sterbenz) and centered on 0: log2 stays relatively
    // accurate near x = 1.
    let e = (x.to_bits() as i32).wrapping_sub(0x3f3504f3) >> 23;
    let m = f32::from_bits((x.to_bits() as i32).wrapping_sub(e << 23) as u32);
    let k = e as f32 + koff;
    let s = m - 1.0;
    // log2(m) = log2(e)*s + s^2*Q(s), with the leading term kept *out* of the
    // polynomial rather than evaluated as its constant coefficient.
    let c = LOG2_Q_COEFFS;
    let s2 = s * s;
    let s4 = s2 * s2;
    let l0 = fma(c[1], s, c[0]);
    let l1 = fma(c[3], s, c[2]);
    let l2 = fma(c[5], s, c[4]);
    let l3 = fma(c[7], s, c[6]);
    let a = s2 * l0;
    let w0 = fma(l2, s2, l1);
    let w1 = fma(c[8], s2, l3);
    let v = fma(w1, s4, w0);
    let sq = fma(v, s4, a);
    // `k` joins *last*, in its own single rounding. Threading it through the
    // peeled combine instead (`fma(s, LOG2_E, fma(s2, q, k))`) rounds twice at
    // `k`'s scale, and for |k| >= 1 that is the whole error budget: it still
    // reaches max 2 but costs 40x on the average (0.0031 -> 0.1249, measured
    // exhaustively).
    let lm = fma(s, std::f32::consts::LOG2_E, sq);
    lm + k
}

/// log_2 without domain checks: valid for positive normal finite x only (no
/// handling for zero, negative, denormal, inf, or nan -- those give a
/// plausible-looking but wrong finite value instead of NaN/-inf). Mirrors
/// exp2/exp2_checked's fast/full-safety split; drops the denormal-rescale
/// multiply and both post-hoc selects log_2 pays on every call.
#[inline(always)]
pub fn log_2_unchecked(x: f32) -> f32 {
    log_2_normal(x, 0.0)
}

/// exp2 without domain checks: valid for x in [-126, 128), i.e. normal
/// (non-denormal, finite, nonzero) results only.
#[doc(alias = "exp2f")]
#[inline(always)]
#[allow(clippy::approx_constant)] // g0's constant term is a fitted minimax
pub fn exp2(x: f32) -> f32 {
    // exp2(floor(x)) * exp2(fract(x)) == exp2(x). exp2int must come from the
    // same floor(x) as f: computing it from x + 383 double-counts the integer
    // part when x + 383 rounds up across an integer (e.g.
    let k = x.floor();
    let f = x - k;
    let exp2int = exp2int_field!(k);
    // Q(f) = (2^f - 1)/f, degree 5, grouped into 3 balanced pairs (g0, g1, g2)
    // instead of two degree-2 Horner halves: same 6 coefficients and the same
    // 4-deep fma critical path, but the combine only ever needs f^2 (never
    // exp2int*f^4), so 2 fewer plain multiplies per call. Avg ulp 0.069, max 1
    // (dense sweep of the whole unchecked domain).
    let q = exp2_q_poly!(f);
    fma(q, exp2int * f, exp2int)
}

/// `2^(k+f)` for an already-integer-valued `k` and `f` in `[0,1)`: [`exp2`]'s
/// own combine step, exposed directly for user custom-base kernels (and this
/// crate's own `exp10`, though not refactored to call through here -- its
/// existing, separately verified body is untouched) that already have their own
/// `k`/`f` and want to skip re-deriving this exact exponent-field-plus-poly
/// combine. No domain check: same `[-126,128)` contract as `exp2` itself
/// (garbage out for `k` outside that range), and `f` outside `[0,1)` is simply
/// a different (still well-defined) input to the same polynomial, not a checked
/// contract.
#[inline(always)]
#[allow(clippy::approx_constant)]
pub fn exp2_kf(k: f32, f: f32) -> f32 {
    let exp2int = exp2int_field!(k);
    let q = exp2_q_poly!(f);
    fma(q, exp2int * f, exp2int)
}

#[inline(always)]
#[allow(clippy::approx_constant)] // g0's constant term is a fitted minimax
pub fn exp2_checked(x: f32) -> f32 {
    // fully branchless (auto-vectorizes): exp2(x) = P(f) * 2^k1 * 2^k2 with k1
    // + k2 = k = floor(x). Splitting k keeps both power-of-two factors
    // representable over the whole clamped range, so overflow to inf and
    // (correctly rounded) denormal underflow fall out of the two multiplies —
    // no pre-offset, no rescale.
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
    let q = exp2_q_poly!(f);
    // weave t1 into the fma chain (t1*f is exact: both factors normal) so
    // only one multiply (by t2) remains after the polynomial
    let p = fma(q, t1 * f, t1);
    p * t2
}

/// Computes `x * 2^n` (`ldexp`/`scalbn`).
/// Decomposes `x` via [`frexp`], adds `n` to the exponent, and reconstructs the float.
/// Overflows to `+/-inf`, underflows to `+/-0.0`.
#[inline(always)]
pub fn ldexp(x: f32, n: i32) -> f32 {
    let (mantissa, e) = frexp(x);
    let target_exp_wide = e as i64 + n as i64;
    let overflow = target_exp_wide > 128;
    let underflow = target_exp_wide < -151;
    let target_exp = target_exp_wide.clamp(-151, 128) as f32;
    let (t1, t2) = exp2_field_split(target_exp);
    let reconstructed = mantissa * t1 * t2;
    let saturated = if overflow {
        f32::INFINITY.copysign(x)
    } else if underflow {
        0.0f32.copysign(x)
    } else {
        reconstructed
    };
    if x == 0.0 || !x.is_finite() {
        x
    } else {
        saturated
    }
}

/// Decomposes `x` into `(mantissa, exponent)` such that `x == mantissa * 2^exponent`,
/// with `mantissa` in `[0.5, 1)`. If `x` is zero or non-finite, returns `(x, 0)`.
#[inline(always)]
pub fn frexp(x: f32) -> (f32, i32) {
    let ax = x.abs();
    let (xs, koff) = denormal_rescale!(ax);
    let bits = xs.to_bits();
    let raw_exp = (bits >> 23) & 0xFF;
    let mantissa_bits = (bits & 0x807FFFFF) | (126u32 << 23);
    let mantissa = f32::from_bits(mantissa_bits).copysign(x);
    let exponent = raw_exp as i32 - 126 + koff as i32;
    let is_special = x == 0.0 || !x.is_finite();
    (
        if is_special { x } else { mantissa },
        if is_special { 0 } else { exponent },
    )
}

/// 10^x. Naively rounding `x*LOG2_10` once before `exp2_checked` even starts
/// loses precision that grows with `|x|` (the same flaw `exp`'s own doc comment
/// describes for `exp2(x*LOG2_E)`).
macro_rules! exp10_reduction {
    ($x:expr) => {{
        const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
        let kb = fma($x, std::f32::consts::LOG2_10, ROUND_MAGIC);
        let kr = kb - ROUND_MAGIC; // round(x*log2(10)), coarse multiply is fine
        let d = fma(-kr, LOG10_2_HI, $x);
        let d = fma(-kr, LOG10_2_LO, d);
        let fr = d * std::f32::consts::LOG2_10; // small, precise correction in log2 units, in [-0.5, 0.5]
        // floor-adjust (kr, fr) from round's [-0.5,0.5] convention to
        // exp2_checked's own floor-based [0,1) convention -- both ops exact or
        // near-exact since they only ever combine values of comparable
        // magnitude (unlike the rejected single-combine above).
        let adjust = if fr < 0.0 { 1.0 } else { 0.0 };
        let k = kr - adjust;
        let f = fr + adjust;
        (k, f)
    }};
}

#[inline(always)]
#[allow(clippy::approx_constant)] // g0's constant term is a fitted minimax
pub fn exp10_checked(x: f32) -> f32 {
    // Clamped before the reduction starts (matching exp2_checked's own
    // early-clamp pattern) so +-inf can't poison `d = x - kr*LOG10_2` with an
    // inf-inf NaN -- NaN itself passes through unaffected (f32::clamp preserves
    // NaN in the receiver). See this function's own doc comment for why these
    // exact bounds make the old separate `k` clamp redundant (removed).
    let x = x.clamp(-45.154503, 38.53184);
    // Unlike `exp10`, this keeps `round`'s own centered `f in [-0.5, 0.5]`
    // instead of paying `exp10_reduction!`'s floor-adjust (a compare, a select
    // and two add/subs) to reach `exp2_q_poly!`'s `[0,1)` convention. That
    // needs its own Q refit -- see `exp2_q_poly_centered!` -- and is only safe
    // here, not in `exp10`: `round` can put `k` at 128, which the k1/k2 split
    // represents fine but a single exponent field cannot (the reason `exp10`'s
    // own doc comment gives for keeping the adjust).
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let kb = fma(x, std::f32::consts::LOG2_10, ROUND_MAGIC);
    let k = kb - ROUND_MAGIC; // round(x*log2(10))
    let d = fma(-k, LOG10_2_HI, x);
    let d = fma(-k, LOG10_2_LO, d);
    let f = d * std::f32::consts::LOG2_10;
    let (t1, t2) = exp2_field_split(k);
    let q = exp2_q_poly_centered!(f);
    let p = fma(q, t1 * f, t1);
    p * t2
}

/// Same reduction as [`exp10_checked`], but a single exponent-field
/// construction (no k1/k2 split) instead of two -- faster, narrower- domain
/// tier, same pairing as `exp2`/`exp2_checked`. Valid while `k` (see
/// `exp10_checked`'s own doc comment) stays in `[-126,128)`.
#[inline(always)]
#[allow(clippy::approx_constant)] // g0's constant term is a fitted minimax
pub fn exp10(x: f32) -> f32 {
    let (k, f) = exp10_reduction!(x);
    let exp2int = exp2int_field!(k);
    let q = exp2_q_poly!(f);
    fma(q, exp2int * f, exp2int)
}

// sin(x) ~= x + x^3*p(x^2) on [-pi/2, pi/2], degree-9 minimax (relative error
// ~6.1e-9), fitted with lolremez. Estrin evaluation, 2 fma chains.
#[inline(always)]
fn sinf_poly_raw(x: f32) -> f32 {
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

/// `sinf_poly_raw` plus the copysign fixup -- see `sinf_poly_raw`'s own doc
/// comment for the full story. Used by every caller except
/// `sin_checked`/`sinpi`, which call `sinf_poly_raw` directly since their own
/// separate zero-handling makes this copysign provably redundant for them.
#[inline(always)]
fn sinf_poly(x: f32) -> f32 {
    sinf_poly_raw(x).copysign(x)
}

// pi split four ways for Cody-Waite reduction with fma. PI_A and PI_B carry
// only 8 and 9 significant bits: the trailing zeros keep the first two `x -
// q*PI_x` steps exactly representable for every |q| the magic-round `q` is
// defined over, so neither of them rounds at all.
const PI_A: f32 = 3.140625;
const PI_B: f32 = 0.0009670257568359375;
const PI_C: f32 = 6.278329465203569e-7;
const PI_D: f32 = 1.0780605906948477e-14;
const FRAC_1_PI: f32 = std::f32::consts::FRAC_1_PI;

// 1.5 * 2^23; adding this to |v| < 2^22 rounds v to the nearest integer
// in the low mantissa bits (round-to-nearest-even).
const ROUND_MAGIC: f32 = 12582912.0;

// The same trick two binades up: adding this to |v| < 2^24 lands v on the
// multiples of 4 instead of the integers. Four times the reach for an index
// whose parity is then known a priori (even) rather than read out of the sum --
// which is the trade `sin` wants, since it has a cheaper place to get the one
// parity bit back from.
const ROUND_MAGIC_4: f32 = 50331648.0;

/// sin(x) via single-f32 range reduction, accurate for `|x| < 2^24 * pi`
/// (~5.27e7). `q = round(x/pi)` has to be an exactly-representable f32 integer
/// for `pi_reduce_and_poly!`'s residual to mean anything, and 2^24 is the
/// largest integer f32 counts by ones; past the limit q lands whole integers
/// off, shifting the residual by whole multiples of pi and pushing it outside
/// sinf_poly's fitted domain [-pi/2, pi/2] -- including returning inf for some
/// large finite x, since nothing here clamps the residual.
macro_rules! pi_reduce_and_poly {
    ($x:expr, $q:expr) => {{
        let r = fma($q, -PI_A, $x);
        let r = fma($q, -PI_B, r);
        let r = fma($q, -PI_C, r);
        let r = fma($q, -PI_D, r);
        sinf_poly(r)
    }};
}

// An index `n` on the grid `$magic` implies, plus the leftover fraction `fc =
// x/pi - n` refined to two words of 1/pi.
macro_rules! frac_x_over_pi {
    ($x:expr, $magic:expr) => {{
        let nb = fma($x, FRAC_1_PI, $magic);
        let n = nb - $magic;
        let f = fma($x, FRAC_1_PI, -n);
        (nb, n, fma($x, RPI_LO, f))
    }};
}

/// Computes `sin(x)` (radians) via single-f32 Cody-Waite range reduction.
/// Accurate for `|x| < 2^24 * pi` (~5.27e7).
#[doc(alias = "sinf")]
#[inline(always)]
pub fn sin(x: f32) -> f32 {
    let (_, n, fc) = frac_x_over_pi!(x, ROUND_MAGIC_4);
    // A second, *fine* magic round, of `fc` this time -- an O(1) value, so the
    // plain integer grid has room to spare: `qb` is `ROUND_MAGIC + round(fc)`,
    // and `q = n + round(fc)` is round(x/pi) exactly, an f32-exact integer for
    // every `|q| <= 2^24` (which is exactly the documented domain). `n -
    // ROUND_MAGIC` is exact (both are multiples of 4) and hangs off `n`, not
    // off `fc`, so the whole `q` chain is no deeper than the single round it
    // replaces -- one add wider, but the same latency.
    let nm = n - ROUND_MAGIC;
    let qb = fc + ROUND_MAGIC;
    let q = nm + qb;
    let s = pi_reduce_and_poly!(x, q);
    // sin(x) = (-1)^q * sin(r); n is a multiple of 4, so q's parity is
    // round(fc)'s, sitting in qb's lowest mantissa bit.
    let parity = qb.to_bits() << 31;
    f32::from_bits(s.to_bits() ^ parity)
}
/// Computes `cos(x)` (radians) via single-f32 Cody-Waite range reduction.
/// Accurate for `|x| < 2^22 * pi` (~1.32e7).
#[doc(alias = "cosf")]
#[inline(always)]
pub fn cos(x: f32) -> f32 {
    // q = the half-odd-integer nearest x/pi, so r = x - q*pi is in [-pi/2,
    // pi/2]. Given `n = round(x/pi)` and `fc = x/pi - n`, that is just `n +
    // copysign(0.5, fc)` -- no second rounding.
    let (nb, n, fc) = frac_x_over_pi!(x, ROUND_MAGIC);
    let q = n + 0.5f32.copysign(fc);
    let s = pi_reduce_and_poly!(x, q);
    // cos(x) = (-1)^n * cos(r +- pi/2) = -+(-1)^n * sin(r), so the
    // half-turn contributes one more sign flip, taken when `fc >= 0`
    // (i.e. when `fc`'s sign bit is clear).
    let parity = (nb.to_bits() << 31) ^ (!fc.to_bits() & SIGN_MASK);
    f32::from_bits(s.to_bits() ^ parity)
}

/// Fast `sin(x)` with single-word `1/pi` reduction. Accurate for `|x| <= 1e6`.
#[inline(always)]
pub fn sin_fast(x: f32) -> f32 {
    let qb = fma(x, FRAC_1_PI, ROUND_MAGIC);
    let q = qb - ROUND_MAGIC;
    let s = pi_reduce_and_poly!(x, q);
    let parity = qb.to_bits() << 31;
    f32::from_bits(s.to_bits() ^ parity)
}

/// Fast `cos(x)` with single-word `1/pi` reduction. Accurate for `|x| <= 1e6`.
#[inline(always)]
pub fn cos_fast(x: f32) -> f32 {
    // k = round(x/pi - 0.5), q = k + 0.5, r = x - q*pi in [-pi/2, pi/2]
    let kb = fma(x, FRAC_1_PI, -0.5) + ROUND_MAGIC;
    let q = (kb - ROUND_MAGIC) + 0.5;
    let s = pi_reduce_and_poly!(x, q);
    // cos(x) = (-1)^(k+1) * sin(r)
    let parity = !kb.to_bits() << 31;
    f32::from_bits(s.to_bits() ^ parity)
}

/// Computes `sin(pi * x)`, argument in half-turns. Exact at integers and total over all finite f32.
#[inline(always)]
pub fn sinpi(x: f32) -> f32 {
    let q = x.round_ties_even();
    let r = x - q;
    // At x=-0.0: q is -0.0 too, so r=(-0.0)-(-0.0), which IEEE754 always
    // resolves to +0.0 -- the same "opposite-signed-zero op erases sign"
    // mechanism as sinf_poly_raw's own -0.0 note, one level further out (r's
    // lost sign means a copysign inside the poly would have nothing left to
    // copy, which is why this guard, not sinf_poly's copysign, is what makes
    // sinpi(-0.0) correct). Same select idiom as log1p/log_2's x==0.0 case:
    // compute the normal path unconditionally, select x itself only at the
    // singular zero point.
    let normal = sinf_poly_raw(std::f32::consts::PI * r) * fma(-2.0, parity(q), 1.0);
    if x == 0.0 {
        x
    } else {
        normal
    }
}

/// Computes `cos(pi * x)`, argument in half-turns. Exact at half-integers and total over all finite f32.
#[inline(always)]
pub fn cospi(x: f32) -> f32 {
    let k = x.round_ties_even();
    let r = x - k;
    let s = sinf_poly_raw(std::f32::consts::PI * (0.5 - r.abs()));
    s * fma(-2.0, parity(k), 1.0)
}

/// Normalized sinc function: `sin(pi * x) / (pi * x)`, with `sinc(0) = 1.0`.
#[inline(always)]
pub fn sinc(x: f32) -> f32 {
    let normal = sinpi(x) / (std::f32::consts::PI * x);
    if x == 0.0 {
        1.0
    } else {
        normal
    }
}

/// Unnormalized sinc function: `sin(x) / x` in radians, with `sinc(0) = 1.0`.
#[inline(always)]
pub fn sinc_unnormalized(x: f32) -> f32 {
    let normal = sin_checked(x) / x;
    if x == 0.0 {
        1.0
    } else {
        normal
    }
}

// The *remainder* of tan(pi*w) after its leading term: with `u = w*w` and `w`
// in [0, 0.25], this is `B(u) = tan(pi*w)/w - fl(pi)`, so that `tan(pi*w) ==
// fma(w, PI, w*B(u))`. Degree 6, two-group Estrin (`lo + u^4*hi`) so no group
// is deeper than the `u^4` it multiplies.
#[inline(always)]
fn tan_poly(u: f32) -> f32 {
    let c: [f32; 7] = [
        -8.742278e-8,
        10.335385,
        40.82169,
        160.9828,
        741.58649,
        701.91418,
        28496.229,
    ];
    let u2 = u * u;
    let u4 = u2 * u2;
    let l0 = fma(c[1], u, c[0]);
    let l1 = fma(c[3], u, c[2]);
    let l2 = fma(c[6], u2, fma(c[5], u, c[4]));
    fma(u4, l2, fma(u2, l1, l0))
}

// tan(pi*r), `r` the exact half-turn-fraction reduction `tanpi`'s own `r`
// already is, and `s = 0.5-|r|` the exact distance to the pole.
#[inline(always)]
fn tan_core(r: f32, s: f32) -> f32 {
    let ar = r.abs();
    let direct = fma(r, std::f32::consts::PI, r * tan_poly(r * r));
    let reflected = mulsign(1.0 / fma(s, std::f32::consts::PI, s * tan_poly(s * s)), r);
    let normal = if ar <= 0.25 { direct } else { reflected };
    if s == 0.0 {
        f32::NEG_INFINITY
    } else {
        normal
    }
}

/// Computes `tan(pi * x)`, argument in half-turns. Total over all finite f32.
#[inline(always)]
pub fn tanpi(x: f32) -> f32 {
    let q = x.round_ties_even();
    let r = x - q;
    let ar = r.abs();
    let normal = tan_core(r, 0.5 - ar);
    if x == 0.0 {
        x
    } else {
        normal
    }
}

/// Computes `sin(2 * pi * x)`, argument in full turns.
#[inline(always)]
pub fn sin2pi(x: f32) -> f32 {
    sinpi(2.0 * x)
}

/// Computes `cos(2 * pi * x)`, argument in full turns.
#[inline(always)]
pub fn cos2pi(x: f32) -> f32 {
    cospi(2.0 * x)
}

/// Computes `tan(2 * pi * x)`, argument in full turns.
#[inline(always)]
pub fn tan2pi(x: f32) -> f32 {
    tanpi(2.0 * x)
}

// 1/180: precomputed reciprocal for the magic-round trick, same idiom as
// sin's own FRAC_1_PI.
const INV_180: f32 = 1.0 / 180.0;
// pi/180, applied only to the *small* (|d|<=90) reduced residual, never
// to the original x -- same "small correction multiplied by an
// irrational constant is fine" reasoning as exp10's own reduction.
const DEG_TO_RAD_SMALL: f32 = std::f32::consts::PI / 180.0;

/// Computes `sin(x)` for `x` in degrees. Accurate for `|x| < 4.7e7`.
#[inline(always)]
pub fn sind(x: f32) -> f32 {
    let qb = fma(x, INV_180, ROUND_MAGIC);
    let q = qb - ROUND_MAGIC;
    let d = fma(-q, 180.0, x);
    let s = sinf_poly((d * DEG_TO_RAD_SMALL).clamp(-POLY_SAFE_BOUND, POLY_SAFE_BOUND));
    let parity = qb.to_bits() << 31;
    f32::from_bits(s.to_bits() ^ parity)
}

/// Computes `cos(x)` for `x` in degrees. Accurate for `|x| < 4.7e7`.
#[inline(always)]
pub fn cosd(x: f32) -> f32 {
    let kb = fma(x, INV_180, -0.5) + ROUND_MAGIC;
    let q = (kb - ROUND_MAGIC) + 0.5;
    let d = fma(-q, 180.0, x);
    let s = sinf_poly((d * DEG_TO_RAD_SMALL).clamp(-POLY_SAFE_BOUND, POLY_SAFE_BOUND));
    let parity = !kb.to_bits() << 31;
    f32::from_bits(s.to_bits() ^ parity)
}

// tan(w degrees) after its leading term, peeled exactly the way `tan_poly`
// peels `tanpi`'s: with `u = w*w` and `w` in [0, 45], this is `B(u) =
// tan(w*pi/180)/w - fl(pi/180)`, so that `tan(w deg) == fma(w,
// DEG_TO_RAD_SMALL, w*B(u))`. Degree 6, two-group Estrin, `c[0]` pinned to
// `pi/180 - fl(pi/180)` exactly.
#[inline(always)]
fn tand_poly(u: f32) -> f32 {
    let c: [f32; 7] = [
        1.3519960e-10,
        1.7721820e-6,
        2.1602940e-10,
        2.6339457e-14,
        3.6819033e-18,
        1.3777568e-22,
        1.3176453e-25,
    ];
    let u2 = u * u;
    let u4 = u2 * u2;
    let l0 = fma(c[1], u, c[0]);
    let l1 = fma(c[3], u, c[2]);
    let l2 = fma(c[6], u2, fma(c[5], u, c[4]));
    fma(u4, l2, fma(u2, l1, l0))
}

// tan of an already-reduced `d` in degrees: `tand_poly` directly for `|d| <=
// 45`, else the cotangent reflection `tan(d) = 1/tan(s)`.
#[inline(always)]
fn tand_core(d: f32) -> f32 {
    let ad = d.abs();
    let s = mulsign(90.0, d) - d;
    let direct = fma(d, DEG_TO_RAD_SMALL, d * tand_poly(d * d));
    let reflected = 1.0 / fma(s, DEG_TO_RAD_SMALL, s * tand_poly(s * s));
    let normal = if ad <= 45.0 { direct } else { reflected };
    if s == 0.0 {
        f32::NEG_INFINITY
    } else {
        normal
    }
}

/// Computes `tan(x)` for `x` in degrees. Accurate for `|x| < 4.7e7`.
#[inline(always)]
pub fn tand(x: f32) -> f32 {
    let qb = fma(x, INV_180, ROUND_MAGIC);
    let q = qb - ROUND_MAGIC;
    let d = fma(-q, 180.0, x);
    tand_core(if d.abs() > 128.0 { 0.0 } else { d })
}

/// `sind` without the safety clamp: valid for `|x| <= 4.7e7`.
#[inline(always)]
pub fn sind_unchecked(x: f32) -> f32 {
    let qb = fma(x, INV_180, ROUND_MAGIC);
    let q = qb - ROUND_MAGIC;
    let d = fma(-q, 180.0, x);
    let s = sinf_poly(d * DEG_TO_RAD_SMALL);
    let parity = qb.to_bits() << 31;
    f32::from_bits(s.to_bits() ^ parity)
}

/// `cosd` without the safety clamp: valid for `|x| <= 4.7e7`.
#[inline(always)]
pub fn cosd_unchecked(x: f32) -> f32 {
    let kb = fma(x, INV_180, -0.5) + ROUND_MAGIC;
    let q = (kb - ROUND_MAGIC) + 0.5;
    let d = fma(-q, 180.0, x);
    let s = sinf_poly(d * DEG_TO_RAD_SMALL);
    let parity = !kb.to_bits() << 31;
    f32::from_bits(s.to_bits() ^ parity)
}

/// `tand` without the safety clamp: valid for `|x| <= 4.7e7`.
#[inline(always)]
pub fn tand_unchecked(x: f32) -> f32 {
    let qb = fma(x, INV_180, ROUND_MAGIC);
    let q = qb - ROUND_MAGIC;
    tand_core(fma(-q, 180.0, x))
}

// q = round(x/pi) must be an *exact* integer for x - q*pi to land accurately in
// [-pi/2, pi/2]. A single-f32 q is a binary either-or: either exactly the
// correctly-rounded integer, or (once |x| crosses q's exact-integer ceiling)
// off by a whole integer, shifting the residual by a whole multiple of pi and
// putting sinf_poly hopelessly outside its fitted domain -- a relocatable
// *cliff*, not a slope, no matter how q is rounded.
const RPI_LO: f32 = 1.2841276486597053e-8;
// 1/pi and pi, each split so the leading word has <= 26 significant bits.
const IPI64_HI: f64 = 0.31830988079309464;
const IPI64_LO: f64 = 5.390696036528002e-09;
const PI64_HI: f64 = 3.1415926218032837;
const PI64_LO: f64 = 3.178650954705639e-08;
// 1.5 * 2^52: the f32 `ROUND_MAGIC` trick one exponent range up. Adding it to
// an exact integer |q| < 2^51 changes no bits of q but parks q's parity in bit
// 0 of the f64, where a 63-bit shift turns it straight into a sign mask.
const ROUND_MAGIC64: f64 = 6755399441055744.0;

/// Error-free transformation (Knuth two-sum): returns `(s, e)` such that `s + e == a + b` exactly.
#[inline(always)]
pub fn two_sum(a: f32, b: f32) -> (f32, f32) {
    let s = a + b;
    let v = s - a;
    let e = (a - (s - v)) + (b - v);
    (s, e)
}

/// Error-free transformation (Fast2Sum): returns `(s, e)` such that `s + e == a + b` exactly.
/// Requires `|a| >= |b|`.
#[inline(always)]
pub fn quick_two_sum(a: f32, b: f32) -> (f32, f32) {
    let s = a + b;
    let e = b - (s - a);
    (s, e)
}

/// Error-free product: returns `(p, e)` such that `p + e == a * b` exactly.
#[inline(always)]
pub fn two_prod(a: f32, b: f32) -> (f32, f32) {
    let p = a * b;
    let e = fma(a, b, -p);
    (p, e)
}

/// Reduces `x` modulo `pi` in f64. Returns `(r, sign)` where `r` is in `[-pi/2, pi/2]`.
#[inline(always)]
fn reduce_pi64<const HALF: bool>(x: f32) -> (f32, u32) {
    let xd = x as f64;
    // Both products are exact: x carries 24 bits, each word at most 26.
    let t = xd * IPI64_HI;
    let tl = xd * IPI64_LO;
    // `t - nh` is exact (both are multiples of ulp(t)), so `fr` is `x/pi - nh`
    // to a relative 2^-53 with no error-free transform -- this is why the f64
    // version is *cheaper* than the double-f32 one it replaced, which needed a
    // `two_prod` at exactly this step.
    let nh = t.round_ties_even();
    let fr = (t - nh) + tl;
    let d = fr.round_ties_even();
    // n = round(x/pi), exact for |n| < 2^53. Formed as `nh + d` rather than
    // rounding `t + tl` directly: that add would quantize at ulp(t), which is
    // already past 1 for the magnitudes this tier exists to serve.
    let n = nh + d;
    // x/pi - n, exact by Sterbenz, |fc| <= 0.5.
    let fc = fr - d;
    // The half-odd-integer nearest x/pi is `n + copysign(0.5, fc)` --
    // no second rounding, the same construction the unchecked `cos`
    // uses and for the same reason.
    let q = if HALF { n + 0.5f64.copysign(fc) } else { n };
    let r = f64::mul_add(-q, PI64_HI, xd);
    let r = f64::mul_add(-q, PI64_LO, r);
    // parity of n, straight out of bit 0 of `n + ROUND_MAGIC64`.
    let par = ((n + ROUND_MAGIC64).to_bits() as u32) << 31;
    // sin: (-1)^n. cos: (-1)^(k+1) for k = round(x/pi - 0.5), which is
    // `n` when fc >= 0 and `n - 1` when fc < 0 -- so the half-turn adds
    // one more flip exactly when fc's sign bit is clear.
    let sgn = if HALF {
        par ^ ((!(fc.to_bits() >> 32)) as u32 & SIGN_MASK)
    } else {
        par
    };
    (r as f32, sgn)
}

/// 1/pi in f64, the small-exponent bypass source of `reduce_pi_wide`'s
/// fraction: for `e < CUT_PI_WIDE` the fraction of x/pi is just |x|/pi -- no
/// wrap past an integer has happened yet -- and one f64 multiply names it to a
/// flat relative 2^-53, which no truncated chunk table could.
const INV_PI_F64: f64 = 0.318309886183790671537767526745028724_f64;

/// Highest raw biased exponent `reduce_pi_wide` serves from its chunk tables;
/// below this, `m*beta < 1` never wraps past an integer, so the smallest true
/// residual at exponent `e` is `beta(e)` itself and the table's absolute 2^-86
/// truncation would be too coarse relative to it. Above (and at) the cut,
/// wrapping decorrelates residuals from beta's own magnitude and the truncation
/// bound is what matters.
const CUT_PI_WIDE: usize = 96;

/// Payne-Hanek range reduction modulo `pi` using 29-bit chunk tables.
/// Total over all finite f32 with no magnitude limit. Returns `(r, sign)`.
#[inline(always)]
fn reduce_pi_wide<const HALF: bool>(x: f32) -> (f32, u32) {
    let b = x.to_bits();
    let e = ((b >> 23) & 0xff) as usize;
    let sgnx = b & SIGN_MASK;
    // Rebuild the 24-bit mantissa integer m: clear sign AND exponent field,
    // then write biased exponent 150, i.e. value m*2^(150-150) = m.
    let exp_field = if e == 255 { 0x7f800000 } else { 0x4b00_0000 };
    let m = f32::from_bits((b & 0x007f_ffff) | exp_field) as f64;
    // One flat table, one base-pointer load: planes 1/2 ride in the
    // gathers' displacements (see pitable.rs).
    let w0 = pitable::REDUCE_PI_W[0][e & 0xff] as f64;
    // mm scales m once; mm1/mm2 differ from it only by exact powers
    // of two, matching where each chunk's bits live in beta (units
    // 2^-28, 2^-57, 2^-86; see pitable.rs for the split).
    let mm = m * 2.0f64.powi(-28);
    // Exact (24 + 29 = 53 bits), so the integer part peels off exactly.
    let p0 = mm * w0;
    let n0 = p0.round_ties_even();
    let f0 = p0 - n0;
    // The two remaining lookups are interleaved with the chain rather than
    // issued up front: all six gathers of an unrolled iteration sharing one
    // issue window made llvm-mca's scheduler queue-stall (~5% on sin_wide's
    // simulated cycles) even though every resource column improved. Same
    // machine code shape either way -- this only spreads the long-latency
    // gathers apart.
    let w1 = pitable::REDUCE_PI_W[1][e & 0xff] as f64;
    let w2 = pitable::REDUCE_PI_W[2][e & 0xff] as f64;
    // Two fmas, not two mul-then-add pairs: p1/p2's scaled products are exact
    // (24 + 29 = 53 bits), so `mul_add(acc, w, acc2)` computes acc + p_exact
    // and rounds ONCE -- at exactly the site `(acc + p)` rounds here.
    // Bit-identical results, two fewer 512-bit f64 pipe ops, and one add
    // shorter on the critical path out of the gathers.
    let s3 = (mm * 2.0f64.powi(-58)).mul_add(w2, (mm * 2.0f64.powi(-29)).mul_add(w1, f0));
    let n1 = s3.round_ties_even();
    // Exact: |s3| <= 0.5625 and |n1| <= 1, so Sterbenz applies.
    let fc = s3 - n1;
    let half = 0.5f64.copysign(fc);
    let tt = if HALF { fc - half } else { fc };
    let r_chain = (tt * std::f64::consts::PI) as f32;
    let par = ((n0 + n1 + ROUND_MAGIC64).to_bits() as u32) << 31;
    let sgn_chain = if HALF {
        par ^ ((!(fc.to_bits() >> 32)) as u32 & SIGN_MASK) ^ sgnx
    } else {
        par
    };

    let (r, sgn) = if e < CUT_PI_WIDE {
        if HALF {
            (-std::f32::consts::FRAC_PI_2, SIGN_MASK ^ sgnx)
        } else {
            (f32::from_bits(b & !SIGN_MASK), 0)
        }
    } else {
        (r_chain, sgn_chain)
    };
    (f32::from_bits(r.to_bits() ^ sgnx), sgn)
}

// parity of an exact-integer float q via floor-based "mod 2" (q*0.5 and its
// floor stay exact once q is an integer), not `q as i64`: Rust's float-to-int
// cast is saturating, which LLVM can't lower to a single vector instruction.
#[inline(always)]
fn parity(q: f32) -> f32 {
    fma(-2.0, (q * 0.5).floor(), q)
}

// Bound for the reduced residual right before it enters sinf_poly, used by
// sind/cosd (whose own reduction has no `|result| <= 1` clamp downstream to
// fall back on). Once the reduction's precision runs out (|x| beyond the
// gradual-degradation range), the residual can grow large -- squaring that
// inside sinf_poly is where an earlier "returns inf for ordinary finite input"
// bug came from.
const POLY_SAFE_BOUND: f32 = 1000.0;

// ---------------------------------------------------------------------------
// Gather-free Payne-Hanek: register-permute window extraction (x8 prototype).

/// The whole `pitable` collapsed into one 48-byte constant (top 2 words zero).
/// Layout: define C = the 256-bit string whose bit j is 1/pi's bit at weight
/// 2^(60-j) (so C covers weights 2^60 down to 2^-195 -- exactly the span the
/// chain's exponent range 96..=254 needs), then store D with D[x] = C[297 - x]
/// (bit-reversed, zero-padded to 384 bits).
#[repr(align(64))]
struct WinAlign([u64; 6]);
static REDUCE_PI_WIN: WinAlign = WinAlign([
    0x39041c0000000000,
    0x4ddc0db6295993c4,
    0x41529fc2757d1f53,
    0x00000a2f9836e4e4,
    0,
    0,
]);

/// Scalar reference for the x8 prototype's differential tests.
#[doc(hidden)]
pub fn reduce_pi_wide_ref<const HALF: bool>(x: f32) -> (f32, u32) {
    reduce_pi_wide::<HALF>(x)
}

/// Raw x8 reduction for differential debugging.
#[doc(hidden)]
pub unsafe fn reduce_pi_wide_x8_pub<const HALF: bool>(
    x: std::arch::x86_64::__m256,
) -> (std::arch::x86_64::__m256, std::arch::x86_64::__m256i) {
    reduce_pi_wide_x8::<HALF>(x)
}

/// Gather-free x8 reduction: same contract as [`reduce_pi_wide`] applied
/// lane-wise to 8 packed f32 inputs. Returns the reduced residuals (packed
/// 8xf32, already xor'd with each lane's own sign, `-0.0` carried) and the sign
/// masks (packed 8xu32 in `SIGN_MASK` position) for `sinf_poly`-style callers.
#[doc(hidden)]
#[inline]
#[target_feature(enable = "avx512f,avx512dq,avx512bw,avx512vbmi,avx512vbmi2")]
unsafe fn reduce_pi_wide_x8<const HALF: bool>(
    x: std::arch::x86_64::__m256,
) -> (std::arch::x86_64::__m256, std::arch::x86_64::__m256i) {
    use std::arch::x86_64::*;

    let xi = _mm256_castps_si256(x);
    let sgnx = _mm256_and_si256(xi, _mm256_set1_epi32(SIGN_MASK as i32));
    let e = _mm256_and_si256(_mm256_srli_epi32::<23>(xi), _mm256_set1_epi32(0xff));

    // Window position: S = 301 - e, byte offset B = S>>3, bit rho = S&7.
    let sv = _mm256_sub_epi32(_mm256_set1_epi32(301), e);
    let bd = _mm256_srli_epi32::<3>(sv);
    let rho = _mm256_and_si256(sv, _mm256_set1_epi32(7));

    // Byte index vectors for the two vpermi2b: byte p of the destination qword
    // lane l wants D byte B_l + (p&7), resp. B_l + 8 + (p&7).
    let b64 = _mm512_cvtepu32_epi64(bd);
    let brep = _mm512_mullo_epi64(b64, _mm512_set1_epi64(0x0101_0101_0101_0101));
    let idx0 = _mm512_add_epi8(brep, _mm512_set1_epi64(0x0706_0504_0302_0100));
    let idx1 = _mm512_add_epi8(brep, _mm512_set1_epi64(0x0f0e_0d0c_0b0a_0908));

    let win = _mm512_load_si512(REDUCE_PI_WIN.0.as_ptr().cast());
    let zero = _mm512_setzero_si512();
    let q0 = _mm512_permutex2var_epi8(win, idx0, zero);
    let q1 = _mm512_permutex2var_epi8(win, idx1, zero);

    let cnt = _mm512_cvtepu32_epi64(rho);
    let xx = _mm512_shrdv_epi64(q0, q1, cnt); // D bits [S, S+63]
    let x2 = _mm512_shrdv_epi64(q1, zero, cnt); // D bits [S+64, S+127]
    let m29 = _mm512_set1_epi64((1 << 29) - 1);
    let w2i = _mm512_and_si512(xx, m29);
    let w1i = _mm512_and_si512(_mm512_srli_epi64::<29>(xx), m29);
    let w0i = _mm512_and_si512(_mm512_shrdi_epi64::<58>(xx, x2), m29);
    let w0 = _mm512_cvtepu64_pd(w0i);
    let w1 = _mm512_cvtepu64_pd(w1i);
    let w2 = _mm512_cvtepu64_pd(w2i);

    // From here the chain is reduce_pi_wide's, lane-wide.
    let mant = _mm256_and_si256(xi, _mm256_set1_epi32(0x007f_ffff));
    let isnn = _mm256_cmpeq_epi32_mask(e, _mm256_set1_epi32(255));
    let expn = _mm256_or_si256(mant, _mm256_set1_epi32(0x4b00_0000));
    let m_dwords = _mm256_mask_blend_epi32(
        isnn,
        expn,
        _mm256_or_si256(mant, _mm256_set1_epi32(0x7f80_0000)),
    );
    // via f32 BITS so inf/NaN payloads poison exactly like the scalar
    // rebuild -- this is a bitcast, not an int->float conversion: the dword
    // already holds the desired f32 pattern.
    let m = _mm512_cvtps_pd(_mm256_castsi256_ps(m_dwords));

    let mm = _mm512_mul_pd(m, _mm512_set1_pd(2.0f64.powi(-28)));
    // NB: these scale mm (not m), i.e. net m*2^-57 / m*2^-86 like the scalar
    let mm1 = _mm512_mul_pd(mm, _mm512_set1_pd(2.0f64.powi(-29)));
    let mm2 = _mm512_mul_pd(mm, _mm512_set1_pd(2.0f64.powi(-58)));
    let p0 = _mm512_mul_pd(mm, w0);
    let n0 = _mm512_roundscale_pd::<{ _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC }>(p0);
    let f0 = _mm512_sub_pd(p0, n0);
    let s3 = _mm512_fmadd_pd(mm2, w2, _mm512_fmadd_pd(mm1, w1, f0));
    let n1 = _mm512_roundscale_pd::<{ _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC }>(s3);
    let fc_chain = _mm512_sub_pd(s3, n1);

    let xf = _mm512_cvtps_pd(_mm256_and_ps(
        x,
        _mm256_set1_ps(f32::from_bits(0x7fff_ffff)),
    ));
    let byp = _mm512_mul_pd(xf, _mm512_set1_pd(INV_PI_F64));
    let small = _mm256_movemask_ps(_mm256_castsi256_ps(_mm256_cmpgt_epi32(
        _mm256_set1_epi32(CUT_PI_WIDE as i32),
        e,
    )));
    let fc = _mm512_mask_blend_pd(small as u8, fc_chain, byp);

    let tt = if HALF {
        // copysign(0.5, fc): the AND already produces the sign-bit mask
        let sb = _mm512_and_si512(_mm512_castpd_si512(fc), _mm512_set1_epi64(1 << 63));
        _mm512_sub_pd(
            fc,
            _mm512_or_pd(_mm512_set1_pd(0.5), _mm512_castsi512_pd(sb)),
        )
    } else {
        fc
    };
    let r = _mm512_cvtpd_ps(_mm512_mul_pd(tt, _mm512_set1_pd(std::f64::consts::PI)));
    let par = _mm512_add_pd(_mm512_add_pd(n0, n1), _mm512_set1_pd(ROUND_MAGIC64));
    let mut sgn = _mm512_srli_epi64::<32>(_mm512_slli_epi64::<63>(_mm512_and_si512(
        _mm512_castpd_si512(par),
        _mm512_set1_epi64(1),
    )));
    if HALF {
        // scalar: (!(fc.to_bits() >> 32)) as u32 & SIGN_MASK -- i.e. the
        // complement of fc's sign bit, which is just fc_high XOR SIGN_MASK.
        let fchi = _mm512_and_si512(
            _mm512_srli_epi64::<32>(_mm512_castpd_si512(fc)),
            _mm512_set1_epi64(0x8000_0000),
        );
        sgn = _mm512_xor_si512(sgn, _mm512_xor_si512(fchi, _mm512_set1_epi64(0x8000_0000)));
    }
    let r_signed = {
        let rb = _mm256_castps_si256(r);
        _mm256_castsi256_ps(_mm256_xor_si256(rb, sgnx))
    };
    // scalar contract: sgnx is always folded into r; for HALF it is also
    // folded into sgn -- the two cancel through the caller's `r ^ flip`,
    // which is what keeps cos even in x.
    let mut sgn = _mm512_cvtepi64_epi32(sgn);
    if HALF {
        sgn = _mm256_xor_si256(sgn, sgnx);
    }
    (r_signed, sgn)
}

/// `sinf_poly` on 8 lanes, plus the `|result| <= 1` clamp (NaN-transparent,
/// same reasoning as the scalar `sin_wide`/`cos_wide` tails).
#[doc(hidden)]
#[inline]
#[target_feature(enable = "avx512f,avx512dq,avx512bw,avx512vbmi,avx512vbmi2")]
unsafe fn sinf_poly_x8(r: std::arch::x86_64::__m256) -> std::arch::x86_64::__m256 {
    use std::arch::x86_64::*;
    let y = _mm256_mul_ps(r, r);
    let y2 = _mm256_mul_ps(y, y);
    let x3 = _mm256_mul_ps(y, r);
    let a = _mm256_fmadd_ps(
        _mm256_set1_ps(8.333_066_2e-3),
        y,
        _mm256_set1_ps(-0.166_666_60),
    );
    let b = _mm256_fmadd_ps(
        _mm256_set1_ps(2.605_780_6e-6),
        y,
        _mm256_set1_ps(-1.980_960_3e-4),
    );
    let p = _mm256_fmadd_ps(b, y2, a);
    let s = _mm256_fmadd_ps(p, x3, r);
    // clamp like `.clamp(-1.0, 1.0)`: compare+blend, not min/max, so NaN
    // lanes survive unchanged (minps/maxps would swallow them)
    let klt = _mm256_cmp_ps_mask::<{ _CMP_LT_OQ }>(s, _mm256_set1_ps(-1.0));
    let s = _mm256_mask_blend_ps(klt, s, _mm256_set1_ps(-1.0));
    let kgt = _mm256_cmp_ps_mask::<{ _CMP_GT_OQ }>(s, _mm256_set1_ps(1.0));
    _mm256_mask_blend_ps(kgt, s, _mm256_set1_ps(1.0))
}

/// `sin_wide` on 8 lanes (residual sign flip folded in).
#[doc(hidden)]
#[inline]
#[target_feature(enable = "avx512f,avx512dq,avx512bw,avx512vbmi,avx512vbmi2")]
unsafe fn sin_wide_lanes_x8(x: std::arch::x86_64::__m256) -> std::arch::x86_64::__m256 {
    use std::arch::x86_64::*;
    let (r, flip) = reduce_pi_wide_x8::<false>(x);
    let r = _mm256_castsi256_ps(_mm256_xor_si256(_mm256_castps_si256(r), flip));
    sinf_poly_x8(r)
}

/// `cos_wide` on 8 lanes.
#[doc(hidden)]
#[inline]
#[target_feature(enable = "avx512f,avx512dq,avx512bw,avx512vbmi,avx512vbmi2")]
unsafe fn cos_wide_lanes_x8(x: std::arch::x86_64::__m256) -> std::arch::x86_64::__m256 {
    use std::arch::x86_64::*;
    let (r, flip) = reduce_pi_wide_x8::<true>(x);
    let r = _mm256_castsi256_ps(_mm256_xor_si256(_mm256_castps_si256(r), flip));
    sinf_poly_x8(r)
}

/// Slice driver for the x8 tier: `out[i] = sin_wide(xs[i])` (bit-identical to
/// the autovectorized scalar path on every input, see the differential sweep in
/// examples/wide_x8.rs). Tail elements use the scalar path.
#[doc(hidden)]
#[inline]
pub unsafe fn sin_wide_x8_slice(xs: &[f32], out: &mut [f32]) {
    let n = xs.len().min(out.len());
    let mut i = 0;
    while i + 8 <= n {
        let xv = std::arch::x86_64::_mm256_loadu_ps(xs.as_ptr().add(i));
        let s = sin_wide_lanes_x8(xv);
        std::arch::x86_64::_mm256_storeu_ps(out.as_mut_ptr().add(i), s);
        i += 8;
    }
    while i < n {
        out[i] = sin_wide(xs[i]);
        i += 1;
    }
}

/// Region-marked throughput driver for llvm-mca: must live in-crate because
/// LLVM refuses to inline `#[target_feature]` functions across crates, and the
/// markers only capture what's inside them.
#[doc(hidden)]
#[inline]
#[target_feature(enable = "avx512f,avx512dq,avx512bw,avx512vbmi,avx512vbmi2")]
pub unsafe fn thr_sin_wide_x8_region(input: &[f32; 16], output: &mut [f32; 16]) {
    use std::arch::x86_64::*;
    unsafe { core::arch::asm!(concat!("# LLVM-MCA-BEGIN ", "sin_wide_x8_throughput")) };
    let mut i = 0;
    while i + 8 <= input.len() {
        let xv = _mm256_loadu_ps(input.as_ptr().add(i));
        let s = sin_wide_lanes_x8(xv);
        _mm256_storeu_ps(output.as_mut_ptr().add(i), s);
        i += 8;
    }
    unsafe { core::arch::asm!("# LLVM-MCA-END") };
}

/// [`sin_wide_x8_slice`] for the cos grid.
#[doc(hidden)]
#[inline]
pub unsafe fn cos_wide_x8_slice(xs: &[f32], out: &mut [f32]) {
    let n = xs.len().min(out.len());
    let mut i = 0;
    while i + 8 <= n {
        let xv = std::arch::x86_64::_mm256_loadu_ps(xs.as_ptr().add(i));
        let s = cos_wide_lanes_x8(xv);
        std::arch::x86_64::_mm256_storeu_ps(out.as_mut_ptr().add(i), s);
        i += 8;
    }
    while i < n {
        out[i] = cos_wide(xs[i]);
        i += 1;
    }
}

#[inline(always)]
pub fn sin_checked(x: f32) -> f32 {
    // sin is odd, so (-1)^q * sin(r) == sin((-1)^q * r): flip r's sign bit
    // *before* sinf_poly instead of negating its result after. Bit-exact with a
    // `s * (1.0 - 2.0 * par)` tail (both IEEE negation and a multiply by
    // exactly +-1 only ever flip the sign bit, never round), but the mask
    // depends solely on q -- ready long before r exits the reduction -- so it
    // hides in the reduction's shadow instead of costing a real fma+mul on
    // sinf_poly's tail.
    let (r, flip) = reduce_pi64::<false>(x);
    let r = f32::from_bits(r.to_bits() ^ flip);
    // `sinf_poly`, not the copysign-free `sinf_poly_raw`: the reduction does
    // carry `-0.0` through intact (that is what the positive pi words buy), but
    // `sinf_poly_raw` then loses it on its own -- its last step is `fma(p, x3,
    // x)` with `p < 0` and `x3 = -0.0`, so the product is `+0.0` and `+0.0 +
    // -0.0` is `+0.0`. Cheaper and more local than the `if x == 0.0 { x }`
    // guard this replaced, which sat at the end of the function blaming the
    // reduction for it.
    sinf_poly(r).clamp(-1.0, 1.0)
}
#[inline(always)]
pub fn cos_checked(x: f32) -> f32 {
    // q = the half-odd-integer nearest x/pi, so r = x - q*pi is in [-pi/2,
    // pi/2] and cos(x) = +-sin(r). Same sign-flip-before-the-poly trick as
    // sin_checked above (sinf_poly is odd in r too); the half-turn's extra flip
    // is already folded into the mask.
    let (r, flip) = reduce_pi64::<true>(x);
    let r = f32::from_bits(r.to_bits() ^ flip);
    // See sin_checked's own comment for why this clamp is needed and why it is
    // the only one needed: `q` stops being an exact integer past |x| ~ 2^53*pi,
    // and without this clamp cos_checked could silently return values like
    // 2.6e21 for legitimate finite input, violating `|cos(x)| <= 1`.
    sinf_poly(r).clamp(-1.0, 1.0)
}

/// Computes `sin(x)` with no magnitude limit across all finite f32 (2 max ulp).
#[inline(always)]
pub fn sin_wide(x: f32) -> f32 {
    let (r, flip) = reduce_pi_wide::<false>(x);
    let r = f32::from_bits(r.to_bits() ^ flip);
    // `sinf_poly`, not `sinf_poly_raw`, for the same `-0.0` reason
    // `sin_checked` gives.
    sinf_poly(r).clamp(-1.0, 1.0)
}

/// Computes `cos(x)` with no magnitude limit across all finite f32 (2 max ulp).
#[inline(always)]
pub fn cos_wide(x: f32) -> f32 {
    let (r, flip) = reduce_pi_wide::<true>(x);
    let r = f32::from_bits(r.to_bits() ^ flip);
    // See `sin_wide` for why this clamp survives an exact reduction.
    sinf_poly(r).clamp(-1.0, 1.0)
}

/// Reduces `x` modulo `pi` in f64. Returns `(r, sign)`.
#[inline(always)]
pub fn reduce_pi_checked(x: f32) -> (f32, f32) {
    let (r, flip) = reduce_pi64::<false>(x);
    // the mask is already in the sign-bit position, so `+-1.0` is one xor
    // away -- no compare, no select.
    (r, f32::from_bits(1.0f32.to_bits() ^ flip))
}

/// Reduces `x` modulo `pi`, offset by half a turn. Returns `(r, sign)`.
#[inline(always)]
pub fn reduce_pi_half_checked(x: f32) -> (f32, f32) {
    let (r, flip) = reduce_pi64::<true>(x);
    (r, f32::from_bits(1.0f32.to_bits() ^ flip))
}

/// Largest magnitude [`wrap_pi`] can return: the largest `f32` whose *exact*
/// value is below `pi`, one ulp under `f32::consts::PI`.
pub const WRAP_PI_MAX: f32 = f32::from_bits(0x40490fda);

/// Wraps `x` (radians) into `(-pi, pi]`.
#[inline(always)]
pub fn wrap_pi(x: f32) -> f32 {
    let (r, sign) = reduce_pi_checked(x);
    let normal = if sign > 0.0 {
        r
    } else {
        r - std::f32::consts::PI.copysign(r)
    };
    // See WRAP_PI_MAX: this is what makes the documented range true, and it
    // costs two instructions with no branch, so it vectorizes with everything
    // above it. It also absorbs the one in-range case the `r - copysign(PI, r)`
    // step gets wrong on its own: `PI` is `pi + 8.7e-8`, so a small positive
    // `r` lands on exactly `-PI`, which is *outside* `(-pi, pi]` however it is
    // rounded (16 inputs over the whole f32 line).
    let normal = normal.clamp(-WRAP_PI_MAX, WRAP_PI_MAX);
    // x=-0.0 needs the same guard sinpi uses, for the same reason: inside
    // `reduce_pi_checked` the residual is formed by subtracting equal signed
    // zeros, which IEEE754 resolves to +0.0, so `r` arrives with the sign
    // already erased and there is nothing left downstream to recover it from.
    if x == 0.0 {
        x
    } else {
        normal
    }
}

/// Computes `sin(r)` for `r` already reduced to `[-pi/2, pi/2]`.
#[inline(always)]
pub fn sin_prereduced(r: f32) -> f32 {
    sinf_poly(r)
}

/// Computes `cos(r)` for `r` already reduced to `[-pi/2, pi/2]`.
#[inline(always)]
pub fn cos_prereduced(r: f32) -> f32 {
    sinf_poly(r)
}

/// Computes `tan(x)` with f64 range reduction (`sin_checked / cos_checked`).
#[inline(always)]
pub fn tan_checked(x: f32) -> f32 {
    let (rs, flip_s) = reduce_pi64::<false>(x);
    let (rc, flip_c) = reduce_pi64::<true>(x);
    let num = sinf_poly(rs).clamp(-1.0, 1.0);
    let den = sinf_poly(rc).clamp(-1.0, 1.0);
    f32::from_bits((num / den).to_bits() ^ (flip_s ^ flip_c))
}

/// Computes `tan(x)` with no magnitude limit across all finite f32.
#[inline(always)]
pub fn tan_wide(x: f32) -> f32 {
    let (rs, flip_s) = reduce_pi_wide::<false>(x);
    let (rc, flip_c) = reduce_pi_wide::<true>(x);
    let num = sinf_poly(rs).clamp(-1.0, 1.0);
    let den = sinf_poly(rc).clamp(-1.0, 1.0);
    f32::from_bits((num / den).to_bits() ^ (flip_s ^ flip_c))
}

/// Core of cbrt for normal finite x: bit-trick seed (~3% error), then a single
/// degree-3 correction. d = s^3 - x is exact-ish via fma at any scale, and
/// x/s^3 == 1/(1+r) exactly for r = d/x, so cbrt(x) = s * (1+r)^(-1/3),
/// approximated by a minimax poly in r.
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
    let c1 = -0.3333314061164856f32;
    let c2 = 0.22221335768699646f32;
    let c3 = -0.1739402711391449f32;
    let c4 = 0.14720453321933746f32;
    let r2 = r * r;
    let a1 = fma(c2, r, c1);
    let b1 = fma(c4, r, c3);
    let p = fma(b1, r2, a1);
    let ss = f32::from_bits(s.to_bits() | (x.to_bits() & SIGN_MASK));
    let sr = ss * r;
    fma(sr, p, ss)
}

#[doc(alias = "cbrtf")]
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
    if ax == 0 || ax >= EXPONENT_MASK {
        x + x
    } else {
        r
    }
}

/// `cbrt` without domain checks: valid for normal finite `x`.
#[inline(always)]
pub fn cbrt_unchecked(x: f32) -> f32 {
    cbrt_normal(x)
}

/// cbrt to within ~0.5 ulp: cbrt_normal (<= 1 ulp), then one Newton step
/// carried out in double-f32 arithmetic. Only valid for x already rescaled into
/// cbrt_accurate's safe range (roughly 2^-56 to 2^127): outside it the
/// double-f32 residual denormalizes/misrounds, or the Newton step's cube can
/// overflow to inf, silently breaking the ~0.5 ulp guarantee (see cbrt_accurate
/// for the rescale).
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn cbrt_accurate_normal(x: f32, scale: f32) -> f32 {
    let ax = x.to_bits() & !SIGN_MASK;
    let a = f32::from_bits(ax);
    // Depends only on x (not y or e), so this whole expression runs fully
    // parallel with cbrt_normal's seed chain and the double-float cube below --
    // folding scale in here too (rather than into a separate multiply after y
    // or e are ready) means the only work left on the critical path after e is
    // a single fma, and the only work left after y is a single multiply
    // (den_recip), instead of two of each. scale is an exact power of two, so
    // this rounds identically to scaling afterwards (no intermediate hits the
    // denormal range).
    let neg_rcp3_scale = -((1.0 / a) * (1.0 / 3.0)) * scale;
    let y = cbrt_normal(x);
    let y2 = Df32::from_mul(y, y);
    let y3 = y2 * y;
    // e = y^3 - x, exact-ish: |e| ~ ulp(x)
    let e = (y3.0 - x) + y3.1;
    // 1/(3y^2) without a second hardware division: y^3 ~ x (cbrt_normal is
    // within ~2 ulp) gives 1/y^2 = y/y^3 ~ y/x = |y|/a, so |y|*rcp3
    // approximates 1/(3y^2). Newton's quadratic convergence only needs den to a
    // handful of accurate bits -- this substitution is bit-exact against a real
    // division over the full accuracy sweep, and a small throughput win (the FP
    // divider is nearly idle while FMA/mul ports are the bottleneck).
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
    let xs = if small {
        x * SCALE_UP
    } else if big {
        x * SCALE_DN
    } else {
        x
    };
    let scale = if small {
        SCALE_UP_OUT
    } else if big {
        SCALE_DN_OUT
    } else {
        1.0
    };
    let r = cbrt_accurate_normal(xs, scale);
    // +-0, +-inf, nan propagate (also kills the rcp=inf NaN for x == +-0)
    if ax == 0 || ax >= EXPONENT_MASK {
        x + x
    } else {
        r
    }
}

/// `cbrt_accurate` without domain checks: valid for normal finite `x`.
#[inline(always)]
pub fn cbrt_accurate_unchecked(x: f32) -> f32 {
    cbrt_accurate_normal(x, 1.0)
}

/// Core of rcbrt for normal finite x, the mirror of [`cbrt_normal`]: bit-trick
/// seed, then a single degree-3 correction on the same `(1+r)^(-1/3)` shape.
/// The seed goes *down* the exponent (`K - ax/3`, not `ax/3 + K`), so `t ~
/// a^(-1/3)` directly and `e = a*t^3 - 1` is formed by an fma with no division
/// anywhere -- which is the whole point, since `1.0 / cbrt(x)` pays two
/// (`cbrt_normal`'s own `1.0/a` reciprocal for its residual, plus the final
/// reciprocal).
#[inline(always)]
fn rcbrt_normal(x: f32) -> f32 {
    let ax = x.to_bits() & !SIGN_MASK;
    let a = f32::from_bits(ax);
    let t = f32::from_bits(0x54a20d0eu32.wrapping_sub(ax / 3));
    let t2 = t * t;
    let e = fma(a * t2, t, -1.0);
    let d1 = -0.3333333134651184f32;
    let d2 = 0.22221790254116058f32;
    let d3 = -0.17284449934959412f32;
    let d4 = 0.14551934599876404f32;
    let d5 = -0.12521736323833466f32;
    let e2 = e * e;
    let a1 = fma(d2, e, d1);
    let b1 = fma(d4, e, d3);
    let p = fma(fma(d5, e2, b1), e2, a1);
    // rcbrt is odd (so is cbrt); the seed was built from |x|, so put the
    // sign back on before the correction rather than after it
    let ts = f32::from_bits(t.to_bits() | (x.to_bits() & SIGN_MASK));
    fma(ts * e, p, ts)
}

/// Computes `x^(-1/3)` (reciprocal cube root).
#[inline(always)]
pub fn rcbrt(x: f32) -> f32 {
    let ax = x.to_bits() & !SIGN_MASK;
    let tiny = ax < 0x0080_0000; // denormal or zero
    let xs = if tiny { x * 16777216.0 } else { x };
    let scale = if tiny { 256.0 } else { 1.0 };
    let r = rcbrt_normal(xs) * scale;
    let spec = f32::from_bits(EXPONENT_MASK.wrapping_sub(ax) | (x.to_bits() & SIGN_MASK));
    if ax == 0 || ax >= EXPONENT_MASK {
        spec
    } else {
        r
    }
}

/// Bit-trick approximation of `cbrt(x)` with two rational refinement steps.
pub fn cbrt_approx(x: f32) -> f32 {
    let y = f32::from_bits(0x2a509849u32 + (x.to_bits() / 3));
    let y = (x + 2. * (y * y) * y) / (3. * (y * y));
    (2. * x * y + (y * y) * (y * y)) / (x + 2. * (y * y) * y)
}
/// Single-bit-trick seed for `sqrt(x)`.
pub fn sqrt_approx(x: f32) -> f32 {
    f32::from_bits(0x1FBD22DF + (x.to_bits() >> 1))
}
/// Single-bit-trick seed for `1/x`.
pub fn rcp_approx(x: f32) -> f32 {
    f32::from_bits(0x7EEF370B - x.to_bits())
}
/// Single-bit-trick seed for `2^x`.
pub fn exp2_approx(x: f32) -> f32 {
    -f32::from_bits((x + 383.).to_bits() << 8)
}

/// Single-bit-trick seed for `log2(x)`.
pub fn log2_approx(x: f32) -> f32 {
    f32::from_bits((x).to_bits() >> 8 | 256_f32.to_bits()) - 383.
}

/// Quake-style bit-trick seed for `1/sqrt(x)`.
pub fn rsqrt_approx(x: f32) -> f32 {
    f32::from_bits(0x5F33E79F - (x.to_bits() >> 1))
}

/// Latency-optimal `cbrt` approximation using two bit-trick seeds.
#[inline(always)]
pub fn cbrt_fast(x: f32) -> f32 {
    let s = f32::from_bits(0x2a4d_def1u32.wrapping_add((x.to_bits() >> 16) * 0x5556u32));
    let r = f32::from_bits(0x68ff_2381u32.wrapping_sub((x.to_bits() >> 16) * 0xaaac));
    let s = fma(s * s, s * -r, fma(r, x, s));
    fma(s * s, s * -r, fma(r, x, s))
}

/// Computes `x * sign(y)` via an XOR of sign bits (preserves magnitude of `x`).
#[inline(always)]
pub fn mulsign(x: f32, y: f32) -> f32 {
    f32::from_bits(x.to_bits() ^ (y.to_bits() & SIGN_MASK))
}

const LN_2: f32 = std::f32::consts::LN_2;
const LOG2_E: f32 = std::f32::consts::LOG2_E;

// Low words of two-word `log2(e)` and `log10(e)`, for the one place each is not
// merely scaling an already-small correction: `log2p1`/`log10p1`'s answer for
// `|x| < 2^-24` *is* `x * log2(e)` (resp. `log10(e)`) and nothing else, because
// `1+x` is then exactly `1.0` and the log kernel contributes an exact zero.
const LOG2_E_LO: f32 = f32::from_bits(0x32a5_7060);
const LOG10_E_LO: f32 = f32::from_bits(0xb22d_91af);
const FRAC_PI_2: f32 = std::f32::consts::FRAC_PI_2;
const FRAC_PI_4: f32 = std::f32::consts::FRAC_PI_4;

// Cody-Waite split of ln(2): LN2_HI keeps its low 9 mantissa bits zeroed, so
// k*LN2_HI (k an exact small integer, this crate's log_2_normal decomposition
// never produces |k| past a couple hundred) is *exact* -- no rounding at all,
// confirmed by brute force for k in [-300, 300]. LN2_LO is the f32-rounded
// residual (LN2 - LN2_HI as f64, then rounded).
const LN2_HI: f32 = 0.693145751953125;
const LN2_LO: f32 = 1.428606765330187e-6;
const LOG10_2_HI: f32 = 0.301025390625;
const LOG10_2_LO: f32 = 4.605039066518657e-6;

/// Computes the natural logarithm of `x`.
#[doc(alias = "logf")]
#[doc(alias = "log")]
#[inline(always)]
#[allow(clippy::neg_cmp_op_on_partial_ord)]
pub fn ln(x: f32) -> f32 {
    log_family_wrapper!(x, ln_normal)
}

/// Core of ln for positive normal finite x only -- see log_2_normal, same
/// contract, same decomposition (`s = m - 1`, exact by Sterbenz) and the same
/// *peeled* poly shape: `ln(m) = s + s^2*Q(s)`, with the leading term kept out
/// of the polynomial rather than evaluated as its constant coefficient. `ln`'s
/// peel is the strictly better of the two, because its leading coefficient is
/// exactly `1.0` and `s` is exact -- so unlike `log_2`, which still has to
/// round `s*log2(e)`, the leading term here is carried with no error at all.
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn ln_normal(x: f32, koff: f32) -> f32 {
    let e = (x.to_bits() as i32).wrapping_sub(0x3f3504f3) >> 23;
    let m = f32::from_bits((x.to_bits() as i32).wrapping_sub(e << 23) as u32);
    let k = e as f32 + koff;
    let s = m - 1.0;
    // `Q(s) = (ln(1+s)/s - 1)/s`, degree 7, fitted by an ulp-weighted LP
    // against `s^2/ln(1+s)` -- the weight that makes the fit minimise the
    // *result*'s relative error, since this poly only ever reaches the answer
    // scaled by `s^2` -- and not by dropping a term off the older `P(s) =
    // ln(1+s)/s`. Written that way the answer for x near 1 (k == 0) simply
    // *was* `s*P(s)`, so all three of P's full-weight evaluation roundings
    // landed straight on the result at ~2^-24 each; peeling demotes every one
    // of them by `s^2/ln(1+s) <= 0.5`.
    let c: [f32; 8] = [
        -0.4999999,
        0.33333948,
        -0.2500179,
        0.1996218,
        -0.16569935,
        0.14916396,
        -0.1431012,
        0.08741673,
    ];
    let s2 = s * s;
    let s4 = s2 * s2;
    let l0 = fma(c[1], s, c[0]);
    let l1 = fma(c[3], s, c[2]);
    let l2 = fma(c[5], s, c[4]);
    let l3 = fma(c[7], s, c[6]);
    // The `s^2` factor rides into the poly's own low group (`a = s2 * l0`)
    // instead of multiplying the finished `Q`, the same way `log_2_normal` does
    // it: same op count, but it keeps the whole thing three Estrin levels deep.
    // Folding `s` in here as well (`a = fma(s2, l0, s)`) saves a further op and
    // is not taken -- it puts a second full-weight rounding back on the result
    // and costs max 1 -> 2.
    let a = s2 * l0;
    let u = fma(l3, s2, l2);
    let w = fma(u, s2, l1);
    let sq = fma(w, s4, a);
    // `s` joins the `k` word rather than `sq`, so the only thing left on the
    // polynomial's critical path is one add and one fma. `base` is exact for k
    // == 0 and rounds at `|s| <= 0.415` otherwise, far under ulp(result) once
    // |k| >= 1; `fma(k, LN2_HI, .)` is then the single full-weight rounding in
    // the whole function.
    let base = fma(k, LN2_LO, s);
    fma(k, LN2_HI, base + sq)
}

/// `ln` without domain checks: valid for positive normal finite `x`.
#[inline(always)]
pub fn ln_unchecked(x: f32) -> f32 {
    ln_normal(x, 0.0)
}

/// Computes the base-10 logarithm of `x`.
#[doc(alias = "log10f")]
#[inline(always)]
#[allow(clippy::neg_cmp_op_on_partial_ord)]
pub fn log10(x: f32) -> f32 {
    log_family_wrapper!(x, log10_normal)
}

/// Core of log10 for positive normal finite x only -- see ln_normal, same
/// contract, same decomposition (`s = m - 1`, exact by Sterbenz) and the same
/// *peeled* poly shape: `log10(m) = s*LOG10_E + s^2*Q(s)`, with the leading
/// term kept out of the polynomial rather than evaluated as its constant
/// coefficient, so the polynomial's own evaluation roundings all reach the
/// answer attenuated by `s^2/log10(1+s) <= 0.35` instead of landing on it at
/// full weight.
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn log10_normal(x: f32, koff: f32) -> f32 {
    let e = (x.to_bits() as i32).wrapping_sub(0x3f3504f3) >> 23;
    let m = f32::from_bits((x.to_bits() as i32).wrapping_sub(e << 23) as u32);
    let k = e as f32 + koff;
    let s = m - 1.0;
    // `Q(s) = (log10(1+s) - s*LOG10_E)/s^2`, degree 7, fitted by an
    // ulp-weighted LP against `s^2/log10(1+s)` -- the weight that makes the fit
    // minimise the *result*'s relative error, since this poly only ever reaches
    // the answer scaled by `s^2`. One degree lower than the un-peeled `P` it
    // replaces, which is what pays for the peel.
    let c: [f32; 8] = [
        -0.21714722,
        0.14476636,
        -0.10857988,
        0.086721875,
        -0.07200416,
        0.06459543,
        -0.06182998,
        0.038040668,
    ];
    let s2 = s * s;
    let s4 = s2 * s2;
    let l0 = fma(c[1], s, c[0]);
    let l1 = fma(c[3], s, c[2]);
    let l2 = fma(c[5], s, c[4]);
    let l3 = fma(c[7], s, c[6]);
    // `k*LOG10_2_LO` rides into the poly's own low group. It depends only on
    // `k`, so it is ready before the polynomial is, and joining it here rather
    // than at the end keeps the tail at two fma -- everything before the
    // closing `fma(k, LOG10_2_HI, .)` then happens at `|s*LOG10_E| <= 0.18`,
    // far under ulp(result) once `|k| >= 1`, leaving that fma as the only
    // full-weight rounding in the function.
    let a = fma(s2, l0, k * LOG10_2_LO);
    let u = fma(l3, s2, l2);
    let w = fma(u, s2, l1);
    let sq = fma(w, s4, a);
    // `s*LOG10_E` must stay *inside* the closing fma. Forming it early --
    // folding it into `a` as `fma(s2, l0, s*LOG10_E)` -- puts a second
    // full-weight rounding back exactly where the peel removed one.
    fma(k, LOG10_2_HI, fma(s, std::f32::consts::LOG10_E, sq))
}

/// `log10` without domain checks: valid for positive normal finite `x`.
#[inline(always)]
pub fn log10_unchecked(x: f32) -> f32 {
    log10_normal(x, 0.0)
}

/// ln(1+x), accurate for small |x| (unlike the naive `ln(1.0 + x)`, which loses
/// x's low bits forming 1.0+x -- the exact case log1p exists to handle -- and
/// rounds to exactly 1.0, hence exactly 0, for |x| below ~6e-8, half of f32's
/// ulp(1.0)). u = 1+x still rounds away those bits, but c = x - (u - 1)
/// recovers the *exact* rounding error (u - 1 is exact by Sterbenz whenever u
/// is within a factor of 2 of 1, i.e.
macro_rules! log1p_nonzero {
    ($x:expr) => {{
        let u = 1.0 + $x;
        let c = $x - (u - 1.0);
        let corr = c / u;
        let corr = if corr.is_finite() { corr } else { 0.0 };
        log_family_wrapper_no_denormal!(u, ln_normal) + corr
    }};
}

#[doc(alias = "log1pf")]
#[allow(clippy::neg_cmp_op_on_partial_ord)]
#[inline(always)]
pub fn log1p(x: f32) -> f32 {
    let normal = log1p_nonzero!(x);
    if x == 0.0 {
        x
    } else {
        normal
    }
}

/// Computes `ln(1 + x) - x`, accurate for small `|x|`.
#[inline(always)]
#[allow(clippy::neg_cmp_op_on_partial_ord)]
pub fn log1pmx(x: f32) -> f32 {
    let u = 1.0 + x;
    let c = x - (u - 1.0);
    let corr = c / u;
    let corr = if corr.is_finite() { corr } else { 0.0 };
    // ln_normal's own reduction, done here rather than by calling it: u = 2^k *
    // m with m in [1/sqrt2, sqrt2), so w = m - 1 is exact (Sterbenz) and lands
    // in [-0.293, 0.415] -- inside the same |z| < 0.5 the poly below is fitted
    // over. That is what lets one poly serve both arms.
    let e = (u.to_bits() as i32).wrapping_sub(0x3f3504f3) >> 23;
    let m = f32::from_bits((u.to_bits() as i32).wrapping_sub(e << 23) as u32);
    let k = e as f32;
    let w = m - 1.0;
    // |x| < 0.5 evaluates the series at x itself; past that at w, whose
    // own log1pmx the identity below turns back into log1pmx(x).
    let z = if x.abs() < 0.5 { x } else { w };
    let z2 = z * z;
    const C1: f32 = -0.6666664970675296;
    const C2: f32 = 0.4999995348856975;
    const C3: f32 = -0.40001991139756876;
    const C4: f32 = 0.33337041934430095;
    const C5: f32 = -0.28505743489396185;
    const C6: f32 = 0.24902498252380087;
    const C7: f32 = -0.23131634089321032;
    const C8: f32 = 0.21157552643272212;
    const C9: f32 = -0.12265182957195057;
    const C10: f32 = 0.09945676037713247;
    const C11: f32 = -0.32486571239903606;
    const C12: f32 = 0.32005194082375105;
    // `z` is only available after the reduction below has produced `w`, so this
    // chain sits on the critical path and a 12-deep serial Horner would
    // dominate it. Estrin instead -- five fma levels rather than twelve, for
    // two extra multiplies -- but only *below* the leading term: `C1..C12` are
    // grouped, then the `1.0` is added by a single trailing fma.
    let z4 = z2 * z2;
    let z8 = z4 * z4;
    let a0 = fma(C2, z, C1);
    let a1 = fma(C4, z, C3);
    let a2 = fma(C6, z, C5);
    let a3 = fma(C8, z, C7);
    let a4 = fma(C10, z, C9);
    let a5 = fma(C12, z, C11);
    let b0 = fma(a1, z2, a0);
    let b1 = fma(a3, z2, a2);
    let b2 = fma(a5, z2, a4);
    let c0 = fma(b1, z4, b0);
    let g = fma(b2, z8, c0);
    let q = fma(g, z, 1.0);
    let p = -0.5 * z2 * q;
    // ln(u) = k*ln2 + ln(m) = k*ln2 + w + log1pmx(w), so log1pmx(x) = ln(u) - x
    // + c/u = (k*LN2_HI - x + w) + p + (k*LN2_LO + c/u) with p the same series
    // value the |x| < 0.5 arm returns directly. Every large term is now
    // differenced *before* anything small is added: k*LN2_HI is exact, w is
    // exact, and over the whole band where |x/log1pmx(x)| is big (it peaks at
    // ~5.3 just past |x| = 0.5) both partial differences are Sterbenz-exact
    // too, so nothing is rounded at ln(u)'s scale and then amplified.
    let big = log_family_edges!(u, {
        let t = fma(k, LN2_HI, -x) + w;
        t + (p + fma(k, LN2_LO, corr))
    });
    let normal = if x.abs() < 0.5 { p } else { big };
    if x == f32::INFINITY {
        f32::NEG_INFINITY
    } else {
        normal
    }
}

/// Computes `log2(1 + x)` (C23 `log2p1`).
#[allow(clippy::neg_cmp_op_on_partial_ord)]
#[inline(always)]
pub fn log2p1(x: f32) -> f32 {
    let u = 1.0 + x;
    let c = x - (u - 1.0);
    // Two-word `log2(e)`: the big product rides inside the fma, so the whole
    // scaled correction carries a single rounding and no constant bias. See
    // `LOG2_E_LO` -- for `|x| < 2^-24` this product is the entire answer, and
    // one word leaves a fixed offset in it.
    let cu = c / u;
    let corr = fma(cu, LOG2_E, cu * LOG2_E_LO);
    let corr = if corr.is_finite() { corr } else { 0.0 };
    let normal = log_family_wrapper_no_denormal!(u, log_2_normal) + corr;
    if x == 0.0 {
        x
    } else {
        normal
    }
}

/// Computes `log10(1 + x)` (C23 `log10p1`).
#[allow(clippy::neg_cmp_op_on_partial_ord)]
#[inline(always)]
pub fn log10p1(x: f32) -> f32 {
    let u = 1.0 + x;
    let c = x - (u - 1.0);
    // Two-word `log10(e)`, exactly as `log2p1` above -- and it matters
    // more here: `LOG10_E`'s single-word offset is 2.33e-8 relative
    // against `LOG2_E`'s 1.33e-8. See `LOG10_E_LO`.
    let cu = c / u;
    let corr = fma(cu, std::f32::consts::LOG10_E, cu * LOG10_E_LO);
    let corr = if corr.is_finite() { corr } else { 0.0 };
    let normal = log_family_wrapper_no_denormal!(u, log10_normal) + corr;
    if x == 0.0 {
        x
    } else {
        normal
    }
}

/// Computes `e^x` via Cody-Waite range reduction.
#[doc(alias = "expf")]
#[inline(always)]
pub fn exp(x: f32) -> f32 {
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let r = fma(-k, LN2_HI, x);
    let r = fma(-k, LN2_LO, r);
    // exp_r_poly!'s c0 and c1 are both pinned to exactly 1.0: exp(r) = 1 + r +
    // r^2*P(r) for tiny r, so a c1 off from 1.0 by even ~6e-8 relative is a
    // systematic bias right where exp(x) is most commonly called (x near 0),
    // and l0 = r + 1.0 needs no fma. c2..c5 are an ulp-weighted Chebyshev LP
    // fit.
    let p = exp_r_poly!(r);
    let (t1, t2) = exp2_field_split(k);
    p * t1 * t2
}

/// Computes `e^x * 2^s`.
#[inline(always)]
pub fn exp_scaled(x: f32, s: i32) -> f32 {
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let r = fma(-k, LN2_HI, x);
    let r = fma(-k, LN2_LO, r);
    let p = exp_r_poly!(r);
    let (t1, t2) = exp2_field_split(k + s as f32);
    p * t1 * t2
}

/// Computes `e^x` via a single exponent field. Valid for `x` in `[-87.3, 88.7]`.
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn exp_narrow(x: f32) -> f32 {
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let r = fma(-k, LN2_HI, x);
    let r = fma(-k, LN2_LO, r);
    let p = exp_r_poly!(r);
    let exp2int = exp2int_field!(k);
    p * exp2int
}

/// Computes `e^x` across the full f32 domain with saturation.
#[inline(always)]
pub fn exp_checked(x: f32) -> f32 {
    exp_reduce!(x.clamp(EXP_CLAMP_LO, EXP_CLAMP_HI))
}

// `e^r - 1` on the Cody-Waite reduction's own `|r| <= ln2/2`, as `r +
// r^2*P(r)`. Fitting `e^r - 1` instead of `e^r` is what removes expm1's
// cancellation at the source: an absolute error carried on `e^r` survives the
// combine's `-1` and comes out amplified by `e^x/(e^x - 1)`, a factor that
// peaks at 2.541 just above `x = 0.5`.
macro_rules! expm1_p_poly {
    ($r:expr, $r2:expr) => {{
        let c: [f32; 5] = [0.5, 1.6666504e-1, 4.1666778e-2, 8.3707254e-3, 1.3916677e-3];
        let l1 = fma(c[1], $r, c[0]);
        let l2 = fma(c[3], $r, c[2]);
        let m = fma(c[4], $r2, l2);
        fma(m, $r2, l1)
    }};
}

macro_rules! expm1_r_poly {
    ($r:expr) => {{
        let r = $r;
        let r2 = r * r;
        fma(expm1_p_poly!(r, r2), r2, r)
    }};
}

// `expm1`'s exponent field, emitted at `k-1`: `exp2int_field!`'s magic with the
// `+127` bias one lower. See `expm1`'s own doc comment for why the field has to
// sit at `k-1` rather than `k`.
const EXPM1_HALF_MAGIC: f32 = 12583038.0;

// Below this, `expm1(x)` is `x` to the last bit, and the `k-1` field's
// halved intermediate would be denormal. See `expm1`'s doc comment.
// 2^-125
const EXPM1_LINEAR: f32 = 2.0 * f32::MIN_POSITIVE;

/// Computes `e^x - 1`, avoiding cancellation near zero.
#[doc(alias = "expm1f")]
#[inline(always)]
pub fn expm1(x: f32) -> f32 {
    // Deliberately a standalone copy of exp's reduction (not routed
    // through the public `exp` fn -- a shared-fn attempt regressed an
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let r = fma(-k, LN2_HI, x);
    let r = fma(-k, LN2_LO, r);
    let e = expm1_r_poly!(r);
    let t = f32::from_bits((k + EXPM1_HALF_MAGIC).to_bits() << 23);
    let b = fma(e, t, t - 0.5);
    let b = b + b;
    if x.abs() < EXPM1_LINEAR {
        x
    } else {
        b
    }
}

/// `expm1` via a single exponent field. Valid for `x` in `[-87.3, 88.7]`.
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn expm1_narrow(x: f32) -> f32 {
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let r = fma(-k, LN2_HI, x);
    let r = fma(-k, LN2_LO, r);
    let e = expm1_r_poly!(r);
    let t = exp2int_field!(k);
    let b = fma(e, t, t - 1.0);
    f32::from_bits(b.to_bits() | (x.to_bits() & SIGN_MASK))
}

/// `expm1` across the full f32 domain with saturation.
#[inline(always)]
pub fn expm1_checked(x: f32) -> f32 {
    let xc = x.clamp(-86.0, 88.72283911167308);
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(xc, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let r = fma(-k, LN2_HI, xc);
    let r = fma(-k, LN2_LO, r);
    let e = expm1_r_poly!(r);
    let t = f32::from_bits((k + EXPM1_HALF_MAGIC).to_bits() << 23);
    let b = fma(e, t, t - 0.5);
    let b = b + b;
    if x.abs() < EXPM1_LINEAR {
        x
    } else {
        b
    }
}

/// Computes `(e^x - 1) / x`, avoiding cancellation near zero.
#[inline(always)]
pub fn exp_m1_over_x(x: f32) -> f32 {
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let r = fma(-k, LN2_HI, x);
    let r = fma(-k, LN2_LO, r);
    let r2 = r * r;
    let p = expm1_p_poly!(r, r2);
    let q = fma(r, p, 1.0);
    let e = fma(r2, p, r);
    let t = f32::from_bits((k + EXPM1_HALF_MAGIC).to_bits() << 23);
    let b = fma(e, t, t - 0.5);
    let b = b + b;
    if k == 0.0 {
        q
    } else {
        b / x
    }
}

/// `exp_m1_over_x` via a single exponent field. Valid for `x` in `[-87.3, 88.7]`.
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn exp_m1_over_x_narrow(x: f32) -> f32 {
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let r = fma(-k, LN2_HI, x);
    let r = fma(-k, LN2_LO, r);
    let r2 = r * r;
    let p = expm1_p_poly!(r, r2);
    let q = fma(r, p, 1.0);
    let e = fma(r2, p, r);
    let t = exp2int_field!(k);
    let b = fma(e, t, t - 1.0);
    if k == 0.0 {
        q
    } else {
        b / x
    }
}

// `2^f - 1` for the round-based reduction's own `|f| <= 0.5`, as `f*ln2 +
// f^2*P(f)`. The leading `f*ln2` is peeled out into the closing fma rather than
// left as a polynomial coefficient, so `x` is never rescaled into `e`'s units
// at all -- the same lever that took `tanpi` from max 5 to 2, and it matters
// more here because the old form's worst case sat exactly on the rounding of `y
// = x*LN_2`.
macro_rules! exp2m1_f_poly {
    ($f:expr, $f2:expr) => {{
        let c: [f32; 5] = [
            2.402265e-1,
            5.55035e-2,
            9.618533e-3,
            1.3395752e-3,
            1.526698e-4,
        ];
        let l1 = fma(c[1], $f, c[0]);
        let l2 = fma(c[3], $f, c[2]);
        let m = fma(c[4], $f2, l2);
        fma(m, $f2, l1)
    }};
}

// 2^-124. Below this the `k-1` field's halved intermediate is denormal, exactly
// as in `expm1`.
const EXP2M1_LINEAR: f32 = 4.0 * f32::MIN_POSITIVE;

/// Computes `2^x - 1`, avoiding cancellation near zero.
#[inline(always)]
pub fn exp2m1(x: f32) -> f32 {
    let xs = x.clamp(-126.0, 128.0);
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = (xs + ROUND_MAGIC) - ROUND_MAGIC;
    let f = xs - k;
    let f2 = f * f;
    let fl = f * LN_2;
    let big = fma(f2, exp2m1_f_poly!(f, f2), fl);
    let t = f32::from_bits((k + EXPM1_HALF_MAGIC).to_bits() << 23);
    let b = fma(big, t, t - 0.5);
    let b = b + b;
    if x.abs() < EXP2M1_LINEAR {
        fl
    } else {
        b
    }
}

// `(10^d - 1 - d*LN_10) / d^2` over `|d| <= 0.5*log10(2)`, degree 4,
// ulp-weighted minimax LP. Same Estrin fold and same degree as `expm1_p_poly!`,
// and for the same reason: the reduction lands `d` where `|d*ln10| <= ln2/2`,
// exactly `expm1`'s own `|r|` bound, so the two approximants are the same
// object in rescaled coordinates.
macro_rules! exp10m1_d_poly {
    ($d:expr, $d2:expr) => {{
        let c: [f32; 5] = [2.650949, 2.0346525, 1.1712452, 0.5420898, 0.20779254];
        let l1 = fma(c[1], $d, c[0]);
        let l2 = fma(c[3], $d, c[2]);
        let m = fma(c[4], $d2, l2);
        fma(m, $d2, l1)
    }};
}

// 2^-125. Below this the `k-1` field's halved intermediate `b` is denormal (`b
// ~ x*ln10/2`, so this is only reachable for denormal `x`) and the doubling
// cannot put back the bit rounding took.
const EXP10M1_LINEAR: f32 = 2.0 * f32::MIN_POSITIVE;

/// Computes `10^x - 1`, avoiding cancellation near zero.
#[inline(always)]
pub fn exp10m1(x: f32) -> f32 {
    let xs = x.clamp(-37.0, 38.53184);
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(xs, std::f32::consts::LOG2_10, ROUND_MAGIC) - ROUND_MAGIC;
    let d = fma(-k, LOG10_2_HI, xs);
    let d = fma(-k, LOG10_2_LO, d);
    let d2 = d * d;
    let dl = d * std::f32::consts::LN_10;
    let big = fma(d2, exp10m1_d_poly!(d, d2), dl);
    let t = f32::from_bits((k + EXPM1_HALF_MAGIC).to_bits() << 23);
    let b = fma(big, t, t - 0.5);
    let b = b + b;
    if x.abs() < EXP10M1_LINEAR {
        dl
    } else {
        b
    }
}

// exp2_checked's k1/k2 exponent-field split, factored out for
// exp/expm1/exp_pos_neg and friends. Any k1 + k2 == k works as long as both
// halves stay inside the exponent field, so k1 = round(k/2).
#[inline(always)]
fn exp2_field_split(k: f32) -> (f32, f32) {
    let a = fma(k, 0.5, EXP2INT_MAGIC);
    let k1 = a - EXP2INT_MAGIC;
    let b = (k - k1) + EXP2INT_MAGIC;
    let t1 = f32::from_bits(a.to_bits() << 23);
    let t2 = f32::from_bits(b.to_bits() << 23);
    (t1, t2)
}

// Shared exp(x)/exp(-x) for sinh/cosh: the Cody-Waite reduction only needs to
// happen once, since -x's reduction is exactly (-k, -r). exp's e^r poly splits
// into even/odd parts in r^2 (p(r) = e + r*o), so p(-r) = e - r*o reuses e/o at
// the cost of one more fma instead of a whole second poly; only the final
// exponent-field scaling (2^k vs 2^-k) is genuinely duplicated -- cheap
// bit-trick work, not fma-port pressure.
#[inline(always)]
fn exp_pos_neg_half(x: f32) -> (f32, f32) {
    let (p_pos, p_neg, t1, t2, t1n, t2n) = exp_pos_neg_core!(x);
    (p_pos * t1 * t2, p_neg * t1n * t2n)
}

// exp_pos_neg_half, single-exponent-field tier: both `+k` and `-k` must fit a
// single field's own valid range simultaneously here (unlike exp_narrow's
// one-sided k), which needs a domain a hair tighter than exp_narrow's own --
// see sinh_narrow/cosh_narrow's own doc comment for the exact (symmetric)
// boundary. Standalone copy of exp_pos_neg_core!'s poly (not routed through the
// macro, which hardcodes the split) -- same standalone-copy precedent as
// expm1/exp_checked's own reductions.
#[inline(always)]
fn exp_pos_neg_narrow_half(x: f32) -> (f32, f32) {
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let r = fma(-k, LN2_HI, x);
    let r = fma(-k, LN2_LO, r);
    // Pre-halved, exactly as in `exp_pos_neg_core!` -- see its comment for
    // why this is free and exact, and `sinh`'s for the premature-overflow
    // gap it closes.
    let c: [f32; 5] = [
        0.49999994 * 0.5,
        0.16666521 * 0.5,
        0.041668329 * 0.5,
        8.3687045e-3 * 0.5,
        1.3814511e-3 * 0.5,
    ];
    let r2 = r * r;
    let e = fma(fma(fma(c[4], r2, c[2]), r2, c[0]), r2, 0.5);
    let o = fma(fma(c[3], r2, c[1]), r2, 0.5);
    let p_pos = fma(r, o, e);
    let p_neg = fma(-r, o, e);
    let t = exp2int_field!(k);
    // reciprocal via bit-subtraction, same trick exp_pos_neg_core! uses
    // for t1n/t2n: exact for any power-of-two field.
    let tn = f32::from_bits(0x7F00_0000u32.wrapping_sub(t.to_bits()));
    (p_pos * t, p_neg * tn)
}

// sinh(x) ~ x * P(x^2), a two-fma odd approximation on |x| < 0.5. The leading
// coefficient is pinned to exactly 1.0 (tiny x returns x, its correctly-rounded
// sinh).
#[inline(always)]
fn sinh_small(x: f32) -> f32 {
    let x2 = x * x;
    let c0 = 1.0f32;
    let c1 = 0.1666623055934906f32;
    let c2 = 0.00839646439999342f32;
    let p = fma(fma(c2, x2, c1), x2, c0);
    x * p
}

/// Computes the hyperbolic sine `sinh(x)`.
#[doc(alias = "sinhf")]
#[inline(always)]
pub fn sinh(x: f32) -> f32 {
    let a = sinh_small(x);
    let (ep, en) = exp_pos_neg_half(x);
    let b = ep - en;
    if x.abs() < 0.5 {
        a
    } else {
        b
    }
}

/// Computes the hyperbolic cosine `cosh(x)`.
#[doc(alias = "coshf")]
#[inline(always)]
pub fn cosh(x: f32) -> f32 {
    let (ep, en) = exp_pos_neg_half(x);
    ep + en
}

/// `sinh` via a single exponent field. Valid for `|x| <= 88.7`.
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn sinh_narrow(x: f32) -> f32 {
    let a = sinh_small(x);
    let (ep, en) = exp_pos_neg_narrow_half(x);
    let b = ep - en;
    if x.abs() < 0.5 {
        a
    } else {
        b
    }
}

/// `cosh` via a single exponent field. Valid for `|x| <= 88.7`.
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn cosh_narrow(x: f32) -> f32 {
    let (ep, en) = exp_pos_neg_narrow_half(x);
    ep + en
}

/// `exp_pos_neg_half`, but with `x` clamped first so `k = round(x*log2e)` never
/// leaves the safe range for *both* `exp2_field_split(k)` and its
/// reciprocal-based negation, AND returning `0.5*exp(x)`/`0.5*exp(-x)` directly
/// instead of the raw pair:
#[inline(always)]
fn exp_pos_neg_checked_half(x: f32) -> (f32, f32) {
    let x = x.clamp(-170.0, 170.0);
    let (p_pos, p_neg, t1, t2, t1n, t2n) = exp_pos_neg_core!(x);
    (p_pos * t1 * t2, p_neg * t1n * t2n)
}

/// `sinh` across the full f32 domain with saturation.
#[inline(always)]
pub fn sinh_checked(x: f32) -> f32 {
    let a = sinh_small(x);
    let (ep, en) = exp_pos_neg_checked_half(x);
    let b = ep - en;
    if x.abs() < 0.5 {
        a
    } else {
        b
    }
}

/// `cosh` across the full f32 domain with saturation.
#[inline(always)]
pub fn cosh_checked(x: f32) -> f32 {
    let (ep, en) = exp_pos_neg_checked_half(x);
    ep + en
}

// coshm1's own poly: `Q(u) = (cosh(sqrt(u)) - 1 - u/2) / u^2` on `u = x^2` in
// `[0, 4]`, i.e. the series with *both* leading terms peeled off, not just the
// constant.
const COSHM1_Q_COEFFS: [f32; 4] = [0.04166663, 0.0013889787, 2.4732362e-5, 2.963389e-7];

/// Computes `cosh(x) - 1`, avoiding cancellation near zero.
#[inline(always)]
pub fn coshm1(x: f32) -> f32 {
    let u = x * x;
    let u2 = u * u;
    let c = COSHM1_Q_COEFFS;
    let a = fma(c[1], u, c[0]);
    let b = fma(c[3], u, c[2]);
    let q = fma(b, u2, a);
    // x*(0.5*x), not 0.5*(x*x): the halved factor keeps the leading term
    // a single rounding even where x*x itself would land denormal.
    let small = fma(u2, q, x * (0.5 * x));
    let big = cosh_checked(x) - 1.0;
    if x.abs() < 2.0 {
        small
    } else {
        big
    }
}

/// Throughput-optimized `sinh(x)`.
#[inline(always)]
pub fn sinh_throughput(x: f32) -> f32 {
    let a = sinh_small(x);
    let e = exp(x);
    let b = 0.5 * (e - 1.0 / e);
    if x.abs() < 0.5 {
        a
    } else {
        b
    }
}

/// Throughput-optimized `cosh(x)`.
#[inline(always)]
pub fn cosh_throughput(x: f32) -> f32 {
    let e = exp(x);
    0.5 * (e + 1.0 / e)
}

/// Computes the hyperbolic tangent `tanh(x)`.
#[doc(alias = "tanhf")]
#[inline(always)]
pub fn tanh(x: f32) -> f32 {
    let xc = x.clamp(-43.5, 44.0);

    // Small arm's denominator, x*coth(x) over [0, 0.8]. The numerator is
    // `xc` itself, so it needs no instruction at all.
    const COTH1: f32 = 0.33333313;
    const COTH2: f32 = -0.02221908;
    const COTH3: f32 = 0.0021010686;
    const COTH4: f32 = -0.00018176674;
    let x2 = xc * xc;
    let dp = fma(fma(fma(COTH4, x2, COTH3), x2, COTH2), x2, COTH1);
    let ds = fma(dp, x2, 1.0);

    // Direct arm: a standalone copy of expm1 (not a call through the public
    // `expm1` fn, same shared-helper scheduling risk as everywhere else) but
    // with a single exponent-field construction instead of exp's k1/k2 split --
    // the clamp bound guarantees `k = round(x*2*log2e)` stays in [-126, 127],
    // comfortably short of the k=128 edge case the split exists for. Same
    // fma(p, exp2int, -1.0) tail fusion as expm1.
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    const LOG2_E_X2: f32 = 2.0 * LOG2_E;
    let k = fma(xc, LOG2_E_X2, ROUND_MAGIC) - ROUND_MAGIC;
    const LN2_HI_HALF: f32 = LN2_HI / 2.0;
    const LN2_LO_HALF: f32 = LN2_LO / 2.0;
    let rh = fma(-k, LN2_HI_HALF, xc);
    let rh = fma(-k, LN2_LO_HALF, rh);
    // exp_r_poly!'s c[0..3], each rescaled by 2^(degree) for `rh = r/2`.
    const D0: f32 = 4.0 * 4.9999300e-1;
    const D1: f32 = 8.0 * 1.6667245e-1;
    const D2: f32 = 16.0 * 4.1883811e-2;
    const D3: f32 = 32.0 * 8.3009899e-3;
    let rh2 = rh * rh;
    let l0 = fma(2.0, rh, 1.0);
    let l1 = fma(D1, rh, D0);
    let l2 = fma(D3, rh, D2);
    let m = fma(l2, rh2, l1);
    let p = fma(m, rh2, l0);
    let exp2int = exp2int_field!(k);
    let b = fma(p, exp2int, -1.0);

    // Two selects, one division -- not one select and two divisions.
    // Both arms are already a ratio, so the branch can be taken a step
    // earlier and the divider visited once.
    let small = xc.abs() < 0.8;
    let n = if small { xc } else { b };
    let d = if small { ds } else { b + 2.0 };
    n / d
}

/// Computes the derivative of `tanh`: `1 - tanh(x)^2`.
#[inline(always)]
pub fn tanh_grad(x: f32) -> f32 {
    let q = exp_checked(-2.0 * x.abs());
    4.0 * q / ((1.0 + q) * (1.0 + q))
}

/// Logistic sigmoid: `1 / (1 + exp(-x))`.
#[inline(always)]
pub fn sigmoid(x: f32) -> f32 {
    // Standalone copy of exp's reduction (not routed through the public `exp`
    // fn, same pattern as expm1); the poly itself is shared via `exp_r_poly!`.
    // Single exponent-field construction, NOT exp2_checked's k1/k2 split
    // (tried: it works but doubles throughput cost, an unjustified price here).
    let xc = x.clamp(-88.722839111673, 87.0);
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(xc, -LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let t1 = fma(k, LN2_HI, xc);
    let t2 = fma(k, LN2_LO, t1);
    let r = -t2;
    let p = exp_r_poly!(r);
    let exp2int = exp2int_field!(k);
    let e = p * exp2int;
    1.0 / (1.0 + e)
}

/// Fast piecewise approximation of sigmoid.
#[inline(always)]
pub fn sigmoid_fast(x: f32) -> f32 {
    let xc = x.clamp(-3.288051, 3.288051);
    let poly = fma(-0.006715598, xc * xc, 0.2178126);
    fma(poly, xc, 0.5).clamp(0.0, 1.0)
}

/// Derivative of sigmoid: `sigmoid(x) * (1 - sigmoid(x))`.
#[inline(always)]
pub fn sigmoid_grad(x: f32) -> f32 {
    let e = exp_checked(-x.abs());
    e / ((1.0 + e) * (1.0 + e))
}

/// `log1p(e)` specialized for callers whose own domain already guarantees `e`
/// in `(0, 1]`)`, `something >= 0`, so `e` never leaves that range) -- `u =
/// 1+e` then always lands in `(1, 2]`, safely away from every special case
/// `log1p`'s own wrapper exists for (never zero, negative, denormal, inf, or
/// nan), so this skips straight to `ln_normal` + the Sterbenz correction, no
/// wrapper, no `corr.is_finite()` guard, no `x==0.0` select. Same domain-bypass
/// mechanism as `asinh`/`acosh`'s own `ln_normal` calls, distinct from the
/// already-rejected "log1p small-|x| branch" (which added a poly to *every*
/// general `log1p` call regardless of caller) -- this is a separate callee only
/// reachable from callers whose domain already proves the skipped checks
/// unreachable.
#[inline(always)]
fn log1p_unit(e: f32) -> f32 {
    let c: [f32; 10] = [
        -0.499999881,
        0.333326906,
        -0.249885798,
        0.198979303,
        -0.161293283,
        0.124671057,
        -0.0830737948,
        0.041981101,
        -0.0136313466,
        0.00207291939,
    ];
    let e2 = e * e;
    let e4 = e2 * e2;
    let l0 = fma(c[1], e, c[0]);
    let l1 = fma(c[3], e, c[2]);
    let l2 = fma(c[5], e, c[4]);
    let l3 = fma(c[7], e, c[6]);
    let l4 = fma(c[9], e, c[8]);
    let r0 = fma(l1, e2, l0);
    let r1 = fma(l3, e2, l2);
    let r1b = fma(l4, e4, r1);
    let q = fma(r1b, e4, r0);
    fma(e2, q, e)
}

/// Softplus: `ln(1 + exp(x))`.
#[inline(always)]
pub fn softplus(x: f32) -> f32 {
    let ax = x.abs();
    // `exp_narrow`, not `exp`: the `min(87.0)` above is already the guard its
    // single-exponent-field domain (`x` in `[-87.68311, 88.37627]`) asks for,
    // so the k1/k2 split is dead weight here -- `-ax.min(87.0)` lands in `[-87,
    // 0]` for every input, NaN included (`min` follows IEEE `minNum` and
    // returns `87.0`, and the trailing `is_nan` below restores the NaN).
    // Bit-identical, since `t1 * t2` and the single field are the same exact
    // power of two over this k range.
    let e = exp_narrow(-ax.min(87.0));
    let corr = if ax > 87.0 { 0.0 } else { log1p_unit(e) };
    let normal = x.max(0.0) + corr;
    if x.is_nan() {
        f32::NAN
    } else {
        normal
    }
}

// `exp(-a) * 2^64` for `a` in `[0, 105]`, the reduction the `_checked` tiers of
// the softplus family need: their correction term has to survive down into the
// denormals, which a single 2^k exponent field cannot represent -- but `exp(-a)
// * 2^64` is normal over the whole range, so the caller lands the denormal
// itself with one exact `2^-64` multiply and a single rounding. Biasing the
// field by `+64` is what keeps `k`, which runs to `-152` at `a = 105`, inside
// the field's own `[-126, 127]`.
macro_rules! exp_neg_scaled64 {
    ($a:expr) => {{
        const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
        let a = $a;
        let k = fma(a, -LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
        let t1 = fma(k, LN2_HI, a);
        let t2 = fma(k, LN2_LO, t1);
        let p = exp_r_poly!(-t2);
        p * exp2int_field!(k + 64.0)
    }};
}

/// Softplus across the full f32 domain.
#[inline(always)]
pub fn softplus_checked(x: f32) -> f32 {
    const P64: f32 = 5.421010862427522e-20; // 2^-64, exact
    let e = exp_neg_scaled64!(x.abs().min(105.0)) * P64;
    let normal = x.max(0.0) + log1p_unit(e);
    if x.is_nan() {
        f32::NAN
    } else {
        normal
    }
}

/// Log-sigmoid: `ln(sigmoid(x)) = -softplus(-x)`.
#[inline(always)]
pub fn logsigmoid(x: f32) -> f32 {
    -softplus(-x)
}

/// Log-sigmoid across the full f32 domain.
#[inline(always)]
pub fn logsigmoid_checked(x: f32) -> f32 {
    -softplus_checked(-x)
}

/// Numerically stable `ln(exp(a) + exp(b))`.
#[inline(always)]
pub fn logaddexp(a: f32, b: f32) -> f32 {
    let m = a.max(b);
    let d = (a - b).abs();
    // `exp_narrow` for the same reason as `softplus`'s -- see its comment.
    let e = exp_narrow(-d.min(87.0));
    let corr = if d > 87.0 { 0.0 } else { log1p_unit(e) };
    let normal = m + corr;
    if a.is_nan() || b.is_nan() {
        f32::NAN
    } else {
        normal
    }
}

/// Numerically stable `ln(exp(a) + exp(b))` across the full f32 domain.
#[inline(always)]
pub fn logaddexp_checked(a: f32, b: f32) -> f32 {
    const P64: f32 = 5.421010862427522e-20; // 2^-64, exact
    let m = a.max(b);
    // `min` is IEEE `minNum` and returns `105.0` for a NaN `d` -- which is also
    // what a NaN `a`/`b` produces via `a - b`, and what `a = b = inf` produces
    // out of `inf - inf`. Every one of those wants a dead correction term, and
    // the trailing `is_nan` restores the NaN cases.
    let d = (a - b).abs().min(105.0);
    let e = exp_neg_scaled64!(d) * P64;
    let normal = m + log1p_unit(e);
    if a.is_nan() || b.is_nan() {
        f32::NAN
    } else {
        normal
    }
}

// `ln 2` as a two-word f64, split so that `n * LN2_HI64` is *exact* for every
// `|n| <= 2^20`: the low 21 mantissa bits of `LN2_HI64` are zero, so the
// product needs 32 + 21 = 53 bits at most. `LN2_LO64` is the f64 nearest the
// remainder, leaving a residual of `1.2e-26` -- times the `|n| <= 185` this
// file's f64 exp reduction can reach, still `2e-24`, i.e.
const LN2_HI64: f64 = 0.6931471803691238;
const LN2_LO64: f64 = 1.9082149292705877e-10;

// `(e^f - 1)/f` on `|f| <= ln2/2`, Taylor rather than minimax: this is an
// accurate tier, and the *rounding* of an f64 evaluation (`2^-53` times the sum
// of absolute terms, `e^|f| = 1.41`) already sits an order of magnitude above
// the degree-13 truncation of `4e-18`, so a minimax refit would buy nothing
// that the evaluation does not immediately spend.
const EXP_F64_P: [f64; 13] = [
    1.0,
    0.5,
    0.16666666666666666,
    0.041666666666666664,
    0.008333333333333333,
    0.001388888888888889,
    0.0001984126984126984,
    2.48015873015873e-05,
    2.7557319223985893e-06,
    2.755731922398589e-07,
    2.505210838544172e-08,
    2.08767569878681e-09,
    1.6059043836821613e-10,
];

// `(atanh(s)/s - 1)/u` in `u = s^2` over `u` in `[0, 1/9]`, i.e. `1/3`, `1/5`,
// ...
const ATANH_B64: [f64; 15] = [
    0.3333333333333333,
    0.2,
    0.14285714285714285,
    0.1111111111111111,
    0.09090909090909091,
    0.07692307692307693,
    0.06666666666666667,
    0.058823529411764705,
    0.05263157894736842,
    0.047619047619047616,
    0.043478260869565216,
    0.04,
    0.037037037037037035,
    0.034482758620689655,
    0.03225806451612903,
];

// `log1p(exp(-d))` in f64 for `d` in `[0, 128]`, the correction term of the
// `logaddexp` family computed to an *absolute* `~2^-52` rather than the
// `~2^-24` an f32 chain can reach. See [`logaddexp_accurate`] for why absolute
// is the metric that matters and f32 cannot supply it.
#[inline(always)]
fn log1p_exp_neg_f64(d: f64) -> f64 {
    // e^d = 2^n * e^f, n = round(d*log2e) in [0, 185], |f| <= ln2/2.
    // `n * LN2_HI64` is exact and within a factor of two of `d`, so `t` is
    // exact by Sterbenz and `f` carries a single rounding.
    let nm = f64::mul_add(d, std::f64::consts::LOG2_E, ROUND_MAGIC64);
    let n = nm - ROUND_MAGIC64;
    let t = f64::mul_add(-n, LN2_HI64, d);
    let f = f64::mul_add(-n, LN2_LO64, t);
    let c = EXP_F64_P;
    let f2 = f * f;
    let f4 = f2 * f2;
    let e0 = f64::mul_add(c[1], f, c[0]);
    let e1 = f64::mul_add(c[3], f, c[2]);
    let e2 = f64::mul_add(c[5], f, c[4]);
    let e3 = f64::mul_add(c[7], f, c[6]);
    let e4 = f64::mul_add(c[9], f, c[8]);
    let e5 = f64::mul_add(c[11], f, c[10]);
    let r0 = f64::mul_add(e1, f2, e0);
    let r1 = f64::mul_add(e3, f2, e2);
    let r2 = f64::mul_add(f64::mul_add(c[12], f2, e5), f2, e4);
    let p = f64::mul_add(f64::mul_add(r2, f4, r1), f4, r0);
    // Same exponent-field reconstruction `exp2_f64_to_f32` documents:
    // `nm`'s low 52 bits already hold `n + 2^51`, and `2^51` is a multiple
    // of 4096, so the 12 bits the shift keeps are exactly `n + 1023`.
    let scale = f64::from_bits(nm.to_bits().wrapping_add(1023) << 52);
    let x = f64::mul_add(f, p, 1.0) * scale;
    let s = 1.0 / f64::mul_add(2.0, x, 1.0);
    // `s` runs down to `1/(1+2*e^128) = 1.3e-56`, so `u^4` and `u^8` -- which
    // the Estrin grouping below forms -- would leave f64's normal range from `d
    // ~ 44` onward, well inside the useful domain, and drag a denormal assist
    // through the whole vector when they did. The floor is far below where the
    // tail term matters: `u <= 1e-30` makes `u*Q(u)` a relative `3e-31` of the
    // pinned `2s`, i.e.
    let u = (s * s).max(1e-30);
    let b = ATANH_B64;
    let u2 = u * u;
    let u4 = u2 * u2;
    let u8 = u4 * u4;
    let a0 = f64::mul_add(b[1], u, b[0]);
    let a1 = f64::mul_add(b[3], u, b[2]);
    let a2 = f64::mul_add(b[5], u, b[4]);
    let a3 = f64::mul_add(b[7], u, b[6]);
    let a4 = f64::mul_add(b[9], u, b[8]);
    let a5 = f64::mul_add(b[11], u, b[10]);
    let a6 = f64::mul_add(b[13], u, b[12]);
    let q0 = f64::mul_add(a1, u2, a0);
    let q1 = f64::mul_add(a3, u2, a2);
    let q2 = f64::mul_add(a5, u2, a4);
    let q3 = f64::mul_add(b[14], u2, a6);
    let q = f64::mul_add(f64::mul_add(q3, u4, q2), u8, f64::mul_add(q1, u4, q0));
    let s2 = s + s;
    f64::mul_add(s2 * u, q, s2)
}

/// Accurate `ln(exp(a) + exp(b))` using f64 intermediate correction.
#[inline(always)]
pub fn logaddexp_accurate(a: f32, b: f32) -> f32 {
    let ad = a as f64;
    let bd = b as f64;
    let m = ad.max(bd);
    let d = (ad - bd).abs().min(128.0);
    let normal = (m + log1p_exp_neg_f64(d)) as f32;
    if a.is_nan() || b.is_nan() {
        f32::NAN
    } else {
        normal
    }
}

/// Gaussian Error Linear Unit (GELU): `x * Phi(x)`.
#[inline(always)]
pub fn gelu(x: f32) -> f32 {
    let xa = x.abs();
    // `erfcx_pos`, not `erfcx`: the argument is an absolute value, so
    // erfcx's own x<0 arm is dead and LLVM does not prove that -- same
    // note as `norm_cdf`'s, which shares this factor exactly.
    let r = erfcx_pos(xa * std::f32::consts::FRAC_1_SQRT_2);
    // `NORM_CDF_XS_CLAMP` does double duty here, on the same two conditions it
    // is asserted for: `x^2/2` stays inside `exp_reduce!`'s range, and `e^-p`
    // is already exactly `0.0` at the clamp -- so a clamped input returns its
    // addend untouched, which is the exactly-right answer (`x` for `x > 0`,
    // since `x*Phi(-x)` is then ~2e-46 against a half-ulp of ~5e-7; `-0.0` for
    // `x < 0`, since `0.5*|x|*Phi(-|x|)` is ~1e-46 against the 7.0e-46 that
    // would round up to the smallest denormal).
    let xs = if xa > NORM_CDF_XS_CLAMP {
        NORM_CDF_XS_CLAMP
    } else {
        xa
    };
    let h = 0.5 * xs;
    let p = h * xs;
    let pe = fma(h, xs, -p);
    let e = exp_reduce!(-p);
    // `x` for `x >= 0`, `-0.0` for `x < 0` -- this is `0.5*x*w` with `w` the
    // reflection addend `erfc`/`norm_cdf` build the same way, except that the
    // sign bit is kept rather than masked off so `gelu(-0.0)` comes out of the
    // closing `fma` as `-0.0` and not `+0.0`.
    let s = (x.to_bits() as i32 >> 31) as u32;
    let addend = f32::from_bits(x.to_bits() & (!s | 0x8000_0000));
    // The two `mulsign`s this would otherwise need -- one to sign `0.5*x`
    // and one for the reflection, off `x` and `-x` -- always disagree, so
    // they collapse into the constant negation on `h * e`.
    fma(-(h * e), fma(-r, pe, r), addend)
}

/// SiLU / Swish activation: `x * sigmoid(x)`.
#[inline(always)]
pub fn silu(x: f32) -> f32 {
    let normal = x * sigmoid(x);
    if x == f32::NEG_INFINITY {
        0.0
    } else {
        normal
    }
}

/// SiLU across the full f32 domain.
#[inline(always)]
pub fn silu_checked(x: f32) -> f32 {
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    const TWO64: f32 = 18446744073709551616.0; // 2^64, exact
    let ax = x.abs();
    let k = fma(ax, -LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let t1 = fma(k, LN2_HI, ax);
    let t2 = fma(k, LN2_LO, t1);
    let r = -t2;
    let p = exp_r_poly!(r);
    let e2 = p * exp2int_field!(k + 64.0); // exp(-|x|) * 2^64
    let num = x * if x < 0.0 { e2 } else { TWO64 };
    let normal = num / (TWO64 + e2);
    if ax > 110.0 {
        x.max(-0.0)
    } else {
        normal
    }
}

/// Softsign: `x / (1 + |x|)`.
#[inline(always)]
pub fn softsign(x: f32) -> f32 {
    let normal = x / (1.0 + x.abs());
    if x.is_infinite() {
        x.signum()
    } else {
        normal
    }
}

/// Computes `sqrt(1 + x) - 1`, avoiding catastrophic cancellation near zero.
/// Domain: `x >= -1.0`.
#[inline(always)]
pub fn sqrt1pm1(x: f32) -> f32 {
    let normal = x / ((1.0 + x).sqrt() + 1.0);
    if x.is_infinite() {
        x
    } else {
        normal
    }
}

/// Computes `x^(3/2) = x * sqrt(x)` for `x >= 0`.
#[inline(always)]
pub fn pow_3_2(x: f32) -> f32 {
    x * x.sqrt()
}

/// Core of `x^(2/3)` for `a` positive and normal: `cbrt_normal`'s own bit-trick
/// seed `s` and residual `r = (s^3-a)/a`, but fitting `(1+r)^(-2/3)` instead of
/// `(1+r)^(-1/3)` -- `s^2 * (1+r)^(-2/3) == a^(2/3)` identically, exactly as `s
/// * (1+r)^(-1/3) == a^(1/3)`, so the two-thirds power is a *direct* fit rather
/// than a cube root squared. `scale`/`scale3` carry the caller's denormal
/// rescale (`scale3` is `scale/3`, see below); folding them in here rather than
/// multiplying the return value keeps them off the tail of the dependency
/// chain, which measures a real latency win.
#[inline(always)]
fn pow_2_3_normal(a: f32, scale: f32, scale3: f32) -> f32 {
    let ax = a.to_bits();
    let rcp = 1.0 / a; // independent of the seed chain, starts immediately
    let s = f32::from_bits(ax / 3 + 0x2a509a07u32);
    let s2u = s * s;
    let e2 = fma(s, s, -s2u) * scale3;
    let d = fma(s2u, s, -a);
    let r = d * rcp;
    let s2 = s2u * scale; // exact: scale is a power of two
    let c1 = -0.6666668057441711f32;
    let c2 = 0.5555411577224731f32;
    let c3 = -0.49370095133781433f32;
    let c4 = 0.45774292945861816f32;
    let c5 = -0.440396785736084f32;
    let p = fma(fma(fma(fma(c5, r, c4), r, c3), r, c2), r, c1);
    // s2*r is off the poly's dependency chain, so the tail after p is one
    // fma and one add; e2 rides in as the fma's addend for free.
    s2 + fma(s2 * r, p, e2)
}

/// Computes `x^(2/3)` for `x >= 0`.
#[inline(always)]
pub fn pow_2_3(x: f32) -> f32 {
    // denormal (or zero) rescale: x by 2^24 = (2^8)^3, so the result comes
    // back 2^16 too big. scale3 is scale/3, the kernel's tail weight.
    const OUT: f32 = 1.52587890625e-5; // 2^-16
    const OUT3: f32 = OUT * (1.0 / 3.0);
    let ax = x.to_bits() & !SIGN_MASK;
    let a = f32::from_bits(ax);
    let tiny = ax < 0x0080_0000;
    let ascaled = if tiny { a * 16777216.0 } else { a };
    let scale = if tiny { OUT } else { 1.0 };
    let scale3 = if tiny { OUT3 } else { 1.0 / 3.0 };
    let r = pow_2_3_normal(ascaled, scale, scale3);
    // +-0, +-inf, nan propagate (also kills the rcp=inf NaN for x == +-0)
    if ax == 0 || ax >= EXPONENT_MASK {
        a + a
    } else {
        r
    }
}

/// Hermite smoothstep on `[edge0, edge1]`.
#[inline(always)]
pub fn smoothstep(edge0: f32, edge1: f32, x: f32) -> f32 {
    let t = ((x - edge0) / (edge1 - edge0)).clamp(0.0, 1.0);
    t * t * fma(-2.0, t, 3.0)
}

/// Perlin smootherstep on `[edge0, edge1]`.
#[inline(always)]
pub fn smootherstep(edge0: f32, edge1: f32, x: f32) -> f32 {
    let t = ((x - edge0) / (edge1 - edge0)).clamp(0.0, 1.0);
    let p = fma(t, fma(t, 6.0, -15.0), 10.0);
    t * t * t * p
}

/// Computes the inverse hyperbolic sine `asinh(x)`.
#[doc(alias = "asinhf")]
#[inline(always)]
pub fn asinh(x: f32) -> f32 {
    let ax = x.abs();
    let small = ax < 2048.0;
    // ax*ax shared across direct_sq, inv_ax2, and sm1's numerator -- both
    // branches are always computed unconditionally, so it's shareable
    // (unlike acosh, see its body comment).
    let ax2 = ax * ax;
    let direct_sq = (ax2 + 1.0).sqrt();
    // sqrt(ax^2+1) = ax + 1/(2*ax) - 1/(8*ax^3) + ..., and this branch only
    // runs for ax >= 2048, where the first dropped term is 1/(8*ax^4) <= 7e-15
    // relative -- seven orders of magnitude under f32's own 6e-8, so the
    // two-term form is exact here. Costs one division and one fma instead of a
    // division, an add, a sqrt and a multiply.
    let inv_ax = 1.0 / ax;
    let rescaled_sq = fma(0.5, inv_ax, ax);
    let sq = if small { direct_sq } else { rescaled_sq };
    let sm1 = if small { ax2 / (sq + 1.0) } else { sq - 1.0 };
    let d = ax + sm1;
    // `log1p_finite(d)` and the overflow fallback `ln(ax) + LN_2` each pay
    // their own full `ln`-family poly + wrapper, unconditionally (branchless),
    // for every call -- but `ln(2*ax) = ln(ax) + ln(2)` is exactly
    // `ln_normal`'s own `koff` hook (it adds directly into the pre-combine
    // exponent field `k`, not as a post-hoc add onto an already-rounded
    // `ln(ax)`), so both branches reduce to one shared `ln_normal` call on a
    // selected (argument, koff) pair. `u = 1+d` is always `>= 1` (finite-d
    // branch, `d = ax+sm1 >= 0`) and `ax` is always a genuine positive value
    // here too, so neither ever needs `log_family_wrapper!`'s
    // zero/negative/denormal handling -- only its inf/nan handling, which the
    // trailing overrides below restore (the raw `_normal` core doesn't
    // propagate either, see its own doc comment: `ax` is `NaN`/`+inf` exactly
    // when `x` is, since `d`'s own non-finiteness routes here).
    let finite_d = d.is_finite();
    let u = 1.0 + d;
    let c = d - (u - 1.0);
    let corr = c / u;
    let arg = if finite_d { u } else { ax };
    let koff = if finite_d { 0.0 } else { 1.0 };
    let shared = ln_normal(arg, koff);
    let combined = if finite_d { shared + corr } else { shared };
    // One select restores both non-finite cases: `ax` is already `+inf`
    // for `x = +-inf` and `NaN` for `x = NaN`, and the trailing `mulsign`
    // puts the infinity's sign back.
    let combined = if x.is_finite() { combined } else { ax };
    mulsign(combined, x)
}

/// Computes the inverse hyperbolic cosine `acosh(x)` for `x >= 1`.
#[doc(alias = "acoshf")]
#[inline(always)]
pub fn acosh(x: f32) -> f32 {
    // NOT the same "shared x2" opportunity as asinh: `direct` needs x*x - 1.0
    // computed as a *single* rounding (the fma) because x is near 1 at acosh's
    // domain boundary, where it's a catastrophic- cancellation subtraction -- a
    // rounded-then-reused x2 loses exactly the precision that cancellation
    // needs (measured: max ulp 3 -> 700). `direct` and `inv_x2` each need their
    // own x*x in a different rounding context.
    let direct = fma(x, x, -1.0).sqrt();
    // sqrt(x^2-1) = x - 1/(2*x) - 1/(8*x^3) - ..., same two-term expansion
    // (and same 7e-15 bound at this branch's own x >= 2048) as asinh's --
    // see its body comment.
    let inv_x = 1.0 / x;
    let rescaled = fma(-0.5, inv_x, x);
    let s = if x < 2048.0 { direct } else { rescaled };
    let d = (x - 1.0) + s;
    // Same shared-ln_normal merge as asinh (see its own doc comment for the
    // full mechanism): `log1p_finite(d)` and the `ln(x) + LN_2` overflow
    // fallback each paid a full ln-family poly + wrapper, unconditionally,
    // every call. `x` itself (not `ax`: acosh's domain is `x >= 1`, no sign to
    // strip) and `u = 1+d` are both always positive here, so only inf/nan need
    // restoring after the raw `ln_normal` core -- `x < 1.0` (false for NaN)
    // already runs last and independently supplies the out-of-domain NaN, so
    // the trailing overrides only need to cover the in-domain `x >= 1` inf/nan
    // cases.
    let finite_d = d.is_finite();
    let u = 1.0 + d;
    let c = d - (u - 1.0);
    let corr = c / u;
    let arg = if finite_d { u } else { x };
    let koff = if finite_d { 0.0 } else { 1.0 };
    let shared = ln_normal(arg, koff);
    let combined = if finite_d { shared + corr } else { shared };
    // One select restores both non-finite cases, exactly as in `asinh`:
    // `x` is already `+inf` for `x = +inf` and `NaN` for `x = NaN`, and
    // `-inf` is discarded by the domain check below.
    let combined = if x.is_finite() { combined } else { x };
    if x < 1.0 {
        f32::NAN
    } else {
        combined
    }
}

// atanh(x) ~ x*(1 + x^2/3 + x^4/5 + x^6/7 + ...), a degree-7 minimax refit of
// the odd Taylor series over |x| < 0.25: the leading coefficient is pinned to
// exactly 1.0 (same convention as asin_small/sinh_small), and only 3
// non-leading terms are needed to stay near f32 precision over this narrow
// domain -- a fitted LP found the next two odd terms (x^8, x^10) converge to
// exactly 0, so this is cheaper than the idea's own "~5 odd terms" guess. Same
// role as asin_small/sinh_small: a cheap, cancellation-free small-x numerator
// for the branch below.
#[inline(always)]
fn atanh_small(x: f32) -> f32 {
    let x2 = x * x;
    let c0 = 0.3333338f32;
    let c1 = 0.19981473f32;
    let c2 = 0.15260197f32;
    let p = fma(fma(c2, x2, c1), x2, c0);
    fma(x * x2, p, x)
}

/// Computes the inverse hyperbolic tangent `atanh(x)` for `|x| < 1`.
#[doc(alias = "atanhf")]
#[inline(always)]
#[allow(clippy::neg_cmp_op_on_partial_ord)]
pub fn atanh(x: f32) -> f32 {
    let a = x.abs();
    let small = atanh_small(x);
    let v = 2.0 * a / (1.0 - a);
    // log1p(v), minus the branches this call site can't reach.
    let u = 1.0 + v;
    let c = v - (u - 1.0);
    let corr = c / u;
    let corr = if corr.is_finite() { corr } else { 0.0 };
    let l = ln_normal(u, 0.0) + corr;
    let l = if u <= 0.0 { f32::NAN } else { l };
    let l = if !(u < f32::INFINITY) { u * u } else { l };
    let big = mulsign(0.5 * l, x);
    if a < 0.25 {
        small
    } else {
        big
    }
}

// `asin(sqrt(t))/sqrt(t)` on `t` in `[0, 1/4]`, degree 5, ulp-weighted minimax
// (LP). One poly serves *both* of acos's branches -- see [`acos`] for the
// identity that makes the two arms ask the same question.
#[inline(always)]
fn acos_poly(t: f32) -> f32 {
    let u = 4.2285666e-2f32;
    let u = fma(u, t, 2.4075409e-2);
    let u = fma(u, t, 4.5502156e-2);
    let u = fma(u, t, 7.494872e-2);
    let u = fma(u, t, 1.6666777e-1);
    fma(u, t, 1.0)
}

/// Dedicated asin.
#[inline(always)]
fn asin_poly(x: f32) -> f32 {
    let u = -2.0342327e-3f32;
    let u = fma(u, x, 1.20692495e-2);
    let u = fma(u, x, -3.609505e-2);
    let u = fma(u, x, 8.2684554e-2);
    let u = fma(u, x, -2.1304381e-1);
    fma(u, x, 1.5706329)
}

/// Computes `acos(x)` in radians for `x` in `[-1, 1]`. Result is in `[0, pi]`.
#[doc(alias = "acosf")]
#[inline(always)]
pub fn acos(x: f32) -> f32 {
    const PI: f32 = std::f32::consts::PI;
    // (pi/2 - fl32(pi/2)) / fl32(pi/2). The same ratio recovers pi's own
    // low word, because fl32(pi) is exactly 2*fl32(pi/2) -- same mantissa,
    // exponent one higher -- so one multiply serves all three addends.
    const LO_RATIO: f32 = -2.7827534e-8;
    let na = f32::from_bits(x.to_bits() | SIGN_MASK); // -|x|
    let small = na > -0.5;
    // t = x^2 below the crossover, (1-|x|)/2 above it. The second is
    // exact: |x| >= 1/2 makes 1 - |x| Sterbenz-exact and halving is free.
    let t = if small { x * x } else { fma(na, 0.5, 0.5) };
    let y = t.sqrt();
    let m = mulsign(if small { na } else { y + y }, x);
    let c = if small {
        FRAC_PI_2
    } else if x < 0.0 {
        PI
    } else {
        0.0
    };
    fma(m, acos_poly(t), c * LO_RATIO) + c
}

/// Computes `acos(x)` in degrees for `x` in `[-1, 1]`.
#[inline(always)]
pub fn acosd(x: f32) -> f32 {
    acos(x) * RAD_TO_DEG_HI
}

// acos(x)/pi as its own degree-6 minimax poly. Fitted as an ulp-weighted
// minimax (LP) of `acos(a)/(pi*sqrt(1-a))` against `a`, weighted by the
// `sqrt(1-a)*P(a)` combine's own sensitivity `sqrt(1-a)/ulp(acospi)` for
// whichever of the two halves (`y` and `1-y`) binds harder, then
// coordinate-descended over the f32 quantisation.
#[inline(always)]
fn acospi_poly(x: f32) -> f32 {
    let u = 7.5414003e-4f32;
    let u = fma(u, x, -3.6262998e-3);
    let u = fma(u, x, 8.662372e-3);
    let u = fma(u, x, -1.55939115e-2);
    let u = fma(u, x, 2.8268332e-2);
    let u = fma(u, x, -6.830643e-2);
    fma(u, x, 5e-1)
}

/// Computes `acos(x) / pi` in half-turns for `x` in `[-1, 1]`.
#[inline(always)]
pub fn acospi(x: f32) -> f32 {
    let a = x.abs();
    let y = (1.0 - a).sqrt() * acospi_poly(a);
    mulsign(y, x + 0.0) + if x < 0.0 { 1.0 } else { 0.0 }
}

// Odd approximation asin(x) ~ x * P(x^2), degree 5 in x^2, on |x| < 0.5. The
// leading coefficient is pinned to exactly 1.0 so tiny x returns x (its
// correctly-rounded asin).
#[inline(always)]
fn asin_small(x: f32) -> f32 {
    let x2 = x * x;
    let c0 = 1.0f32;
    let c1 = 0.16666752f32;
    let c2 = 0.074952975f32;
    let c3 = 0.04547038f32;
    let c4 = 0.02417949f32;
    let c5 = 0.042166352f32;
    let p = fma(
        fma(fma(fma(fma(c5, x2, c4), x2, c3), x2, c2), x2, c1),
        x2,
        c0,
    );
    x * p
}

/// Computes `asin(x)` in radians for `x` in `[-1, 1]`.
#[doc(alias = "asinf")]
#[inline(always)]
pub fn asin(x: f32) -> f32 {
    let a = x.abs();
    let small = asin_small(x);
    let big = mulsign(fma(-(1.0 - a).sqrt(), asin_poly(a), FRAC_PI_2), x);
    if a < 0.5 {
        small
    } else {
        big
    }
}

// 180/pi as a double-f32, for the radians-to-degrees composites: HI+LO
// represents it to ~2^-49 relative, so `fma(y, HI, y*LO)` rounds once at the
// result's own magnitude where a single-word `y * K` rounds twice and carries
// whatever bias the f32 `K` has.
const RAD_TO_DEG_HI: f32 = 57.2957763671875;
const RAD_TO_DEG_LO: f32 = 3.1458948e-6;

// The same double-f32 treatment for 1/pi, for the half-turn composites. Here
// [`FRAC_1_PI`] is *already* the low neighbour -- the correctly rounded
// `fl(1/pi)` sits 0.43 ulp *below* the real 1/pi -- so it also serves as the HI
// word and only this positive tail has to be named.
const FRAC_1_PI_LO: f32 = 1.28412765e-8;

/// Computes `asin(x)` in degrees for `x` in `[-1, 1]`.
#[inline(always)]
pub fn asind(x: f32) -> f32 {
    let y = asin(x);
    fma(y, RAD_TO_DEG_HI, y * RAD_TO_DEG_LO)
}

// asin(x)/pi, seeded by rescaling each coefficient by 1/pi and then refit in
// half-turn space directly: NOT assumed safe just because 's own RAD_TO_DEG
// fold measured worse for asind -- verified separately since it's a different
// constant, not just a relabeling. Plausible reason for the opposite verdict:
// `1/pi` (~0.318) keeps rescaled coefficients in a similar-or-smaller magnitude
// range, while `180/pi` (~57.3) inflates them roughly 57x, apparently
// accumulating more absolute rounding per fma step than one final multiply
// costs.
#[inline(always)]
fn asinpi_small(x: f32) -> f32 {
    let x2 = x * x;
    let c1 = 0.053051922f32;
    let c2 = 0.023858273f32;
    let c3 = 0.014473671f32;
    let c4 = 0.0076965746f32;
    let c5 = 0.013421959f32;
    // Degree 5 in `x^2`, covering `[0, 0.5]` -- see [`asinpi`] for why the
    // crossover sits there and `asin_small` for the same sizing argument (a
    // minimax needs degree 5 to reach 0.5, not the ~12 terms the Taylor series
    // would).
    let t = fma(
        fma(fma(fma(fma(c5, x2, c4), x2, c3), x2, c2), x2, c1),
        x2,
        FRAC_1_PI_LO,
    );
    fma(x, FRAC_1_PI, x * t)
}

// asin's `0.5 - sqrt(1-a)*P(a)` branch in half-turns, degree 5 over a in [0.5,
// 1). Fitted as an ulp-weighted minimax (LP) of `acos(a)/(pi*sqrt(1-a))` --
// weight `sqrt(1-a)/ulp(asinpi)`, the combine's own sensitivity -- and
// coordinate-descended over the f32 quantisation.
#[inline(always)]
fn asinpi_poly(x: f32) -> f32 {
    let u = -0.0006497958f32;
    let u = fma(u, x, 0.0038506954);
    let u = fma(u, x, -0.011503078);
    let u = fma(u, x, 0.026329458);
    let u = fma(u, x, -0.06781764);
    fma(u, x, 0.4999485)
}

/// Computes `asin(x) / pi` in half-turns for `x` in `[-1, 1]`.
#[inline(always)]
pub fn asinpi(x: f32) -> f32 {
    let a = x.abs();
    let small = asinpi_small(x);
    let big = mulsign(fma(-(1.0 - a).sqrt(), asinpi_poly(a), 0.5), x);
    if a < 0.5 {
        small
    } else {
        big
    }
}

// 3/3 Pade-style rational approximation of atan on [0,1], seeded from a
// least-squares fit and coordinate-descent tuned. Current: atan avg/max ulp
// 0.063/3 (exhaustive).
#[inline(always)]
fn atan_poly(x: f32) -> f32 {
    let a2 = 0.008830042167832291;
    let a1 = 0.2849778513254418;
    let a0 = 1.1271711055988247;
    let b2 = 5.0166193e-2;
    let b1 = 5.718157e-1;
    let b0 = 1.4605043e0;
    let x2 = x * x;
    let numer = fma(x * x2, fma(fma(a2, x2, a1), x2, a0), x);
    let denom = fma(fma(fma(b2, x2, b1), x2, b0), x2, 1.0);
    numer / denom
}

/// Computes `atan(x)` in radians.
#[doc(alias = "atanf")]
#[inline(always)]
pub fn atan(x: f32) -> f32 {
    let a = x.abs();
    // a >= 0, so min(a, 1/a) picks whichever branch the old a<1.0 select did (a
    // itself below 1, the reciprocal at/above 1) in one vminps instead of a
    // compare+blend; the reciprocal was already computed unconditionally either
    // way (both branches evaluate in the branchless/vectorized style this crate
    // uses). NaN: a=NaN -> 1/a=NaN -> min(NaN, NaN) = NaN, matching the old
    // else-branch's 1.0/NaN.
    let y = a.min(1.0 / a);
    let y = atan_poly(y);
    let y = if a < 1.0 { y } else { FRAC_PI_2 - y };
    mulsign(y, x)
}

/// Computes `atan(x)` in degrees.
#[inline(always)]
pub fn atand(x: f32) -> f32 {
    let y = atan(x);
    fma(y, RAD_TO_DEG_HI, y * RAD_TO_DEG_LO)
}

/// Computes `atan(x) / pi` in half-turns.
#[inline(always)]
pub fn atanpi(x: f32) -> f32 {
    let a = x.abs();
    // `a.min(1.0/a)` is already non-negative, so `atan_bounded`'s own sign
    // handling is dead work -- but LLVM only folds it away if the value is
    // *visibly* sign-cleared, and `.abs()` is what makes it visible. Worth two
    // instructions in the emitted region; a bit-mask spelling of the same thing
    // measures identically, and reaching past `atan_bounded` to the private
    // `atan_poly` would save two more at the cost of merging `atan`'s ownership
    // domain into this one.
    let t = atan_bounded(a.min(1.0 / a).abs());
    let h = fma(t, FRAC_1_PI, t * FRAC_1_PI_LO);
    // Exactly `atan`'s `a < 1.0` branch, in half-turns. `0.5 - h` for
    // `h <= 0.25` rounds once, at the result's own magnitude; the sign is
    // reapplied last because both arms are computed from `|x|`.
    mulsign(if a < 1.0 { h } else { 0.5 - h }, x)
}

/// `atan(x)` for `|x| <= 1`.
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn atan_bounded(x: f32) -> f32 {
    mulsign(atan_poly(x.abs()), x)
}

/// Latency-optimized division-free `atan(x)`.
#[inline(always)]
pub fn atan_latency(x: f32) -> f32 {
    let a = x.abs();
    let r = a.min(1.0 / a);
    let r2 = r * r;
    let r4 = r2 * r2;
    let c0 = -0.33333167;
    let c1 = 0.19994265;
    let c2 = -0.14216055;
    let c3 = 0.10689225;
    let c4 = -0.07608681;
    let c5 = 0.04395558;
    let c6 = -0.01687014;
    let c7 = 0.0030569038;
    let lo = fma(fma(fma(c2, r2, c1), r2, c0), r2, 1.0);
    let hi = fma(fma(fma(fma(c7, r2, c6), r2, c5), r2, c4), r2, c3);
    let p = fma(hi, r4 * r4, lo) * r;
    let sp = mulsign(p, x);
    let hpisignx = mulsign(FRAC_PI_2, x);
    if a < 1.0 {
        sp
    } else {
        hpisignx - sp
    }
}

/// Computes the four-quadrant arctangent `atan2(y, x)` in radians.
#[doc(alias = "atan2f")]
#[inline(always)]
pub fn atan2(y: f32, x: f32) -> f32 {
    let nonzerox = x != 0.0;
    let nonzeroy = y != 0.0;
    let bothzero = !nonzerox && !nonzeroy;
    let hpisignx = if nonzerox || bothzero {
        mulsign(FRAC_PI_2, x)
    } else {
        0.0
    };
    let correction = mulsign(FRAC_PI_2 - hpisignx, y);
    let r = if nonzerox {
        atan(y / x) + correction
    } else {
        correction
    };
    let r = if y.is_nan() { f32::NAN } else { r };
    // atan2(+-inf, +-inf): y/x is inf/inf, which is NaN, so the general formula
    // above can't produce an answer here at all. IEEE754/C99 define a canonical
    // result by quadrant regardless (+-pi/4 or +-3pi/4) -- not derived from any
    // real ratio, since there isn't one at true infinity, just a fixed
    // convention.
    let bothinf = x.is_infinite() && y.is_infinite();
    let inf_result = mulsign(
        if x.is_sign_negative() {
            3.0 * FRAC_PI_4
        } else {
            FRAC_PI_4
        },
        y,
    );
    if bothinf {
        inf_result
    } else {
        r
    }
}

/// Latency-optimized `atan2(y, x)`.
#[inline(always)]
pub fn atan2_latency(y: f32, x: f32) -> f32 {
    let nonzerox = x != 0.0;
    let nonzeroy = y != 0.0;
    let bothzero = !nonzerox && !nonzeroy;
    let hpisignx = if nonzerox || bothzero {
        mulsign(FRAC_PI_2, x)
    } else {
        0.0
    };
    let correction = mulsign(FRAC_PI_2 - hpisignx, y);
    let r = if nonzerox {
        atan_latency(y / x) + correction
    } else {
        correction
    };
    let r = if y.is_nan() { f32::NAN } else { r };
    let bothinf = x.is_infinite() && y.is_infinite();
    let inf_result = mulsign(
        if x.is_sign_negative() {
            3.0 * FRAC_PI_4
        } else {
            FRAC_PI_4
        },
        y,
    );
    if bothinf {
        inf_result
    } else {
        r
    }
}

/// Computes `atan2(y, x)` in degrees in `[-180, 180]`.
#[inline(always)]
pub fn atan2d(y: f32, x: f32) -> f32 {
    let r = atan2(y, x);
    let normal = fma(r, RAD_TO_DEG_HI, r * RAD_TO_DEG_LO);
    // Single-word, and specifically the HI word rather than the
    // correctly-rounded `fl(180/pi)`: the two-word form would need a second
    // division, and swapping in the correctly-rounded single constant -- which
    // does measure better on this branch's average -- costs a 1-ulp regression
    // at a pinned denormal edge case that is currently correctly rounded.
    let tiny = y * (RAD_TO_DEG_HI / x);
    if r.abs() < f32::MIN_POSITIVE && x != 0.0 {
        tiny
    } else {
        normal
    }
}

/// Computes `atan2(y, x) / pi` in half-turns in `[-1, 1]`.
#[inline(always)]
pub fn atan2pi(y: f32, x: f32) -> f32 {
    let nonzerox = x != 0.0;
    let nonzeroy = y != 0.0;
    let bothzero = !nonzerox && !nonzeroy;
    // Exactly `atan2`'s shape with `FRAC_PI_2` replaced by `0.5`. Every
    // value this can take -- `0.5 - 0.5`, `0.5 - -0.5`, `0.5 - 0.0` -- is
    // exact, so the quadrant fold contributes no rounding of its own.
    let hsignx = if nonzerox || bothzero {
        mulsign(0.5, x)
    } else {
        0.0
    };
    let correction = mulsign(0.5 - hsignx, y);
    let r = if nonzerox {
        atanpi(y / x) + correction
    } else {
        correction
    };
    let r = if y.is_nan() { f32::NAN } else { r };
    let bothinf = x.is_infinite() && y.is_infinite();
    let inf_result = mulsign(if x.is_sign_negative() { 0.75 } else { 0.25 }, y);
    if bothinf {
        inf_result
    } else {
        r
    }
}

/// `atan2` without special zero/infinite case handling.
#[inline(always)]
pub fn atan2_unchecked(y: f32, x: f32) -> f32 {
    let hpisignx = mulsign(FRAC_PI_2, x);
    let correction = mulsign(FRAC_PI_2 - hpisignx, y);
    atan(y / x) + correction
}

/// `atan2` folded into `[0, 2*pi)` (single positive turn).
#[inline(always)]
pub fn atan2_pos(y: f32, x: f32) -> f32 {
    let r = atan2(y, x);
    if y.is_sign_negative() {
        r + std::f32::consts::TAU
    } else {
        r
    }
}

/// Computes `tan(x)` in radians for `|x| < 2^22 * pi`.
#[doc(alias = "tanf")]
#[inline(always)]
pub fn tan(x: f32) -> f32 {
    // `nb`, `n`, `fc` and both Cody-Waite chains are `sin`'s and `cos`'s own,
    // inlined here only so the shared reduction is visibly shared: sin's wider
    // multiple-of-4 grid is deliberately not used, since a second, non-shared
    // reduction costs tan ~15 instructions for byte-identical output (mca: 105
    // -> 120 instrs, BlockRT 31 -> 36).
    let (nb, n, fc) = frac_x_over_pi!(x, ROUND_MAGIC);
    // `fc + nb` is a second magic round, of `x/pi` this time instead of
    // `x*FRAC_1_PI`: `nb` is already `ROUND_MAGIC + n` on the integer grid, so
    // adding the (corrected, |fc| <= 0.67) fraction back rounds to `ROUND_MAGIC
    // + N`. `sin` cannot reuse this shape -- its coarse grid is too wide for
    // `nb` to sit on the integer grid at all.
    let qb = fc + nb;
    let num = pi_reduce_and_poly!(x, qb - ROUND_MAGIC);
    let den = pi_reduce_and_poly!(x, n + 0.5f32.copysign(fc));
    let sign = ((qb.to_bits() ^ nb.to_bits()) << 31) ^ (!fc.to_bits() & SIGN_MASK);
    f32::from_bits((num / den).to_bits() ^ sign)
}

// degree-6 minimax poly feeding erf's exp2-based tail (|x| >= 0.28), Estrin (3
// fma's deep instead of Horner's 6, accuracy-neutral here -- fma reassociation
// has to be checked per poly, not assumed either way). Coefficients are an
// ulp-weighted Chebyshev LP fit, weighted by the linearized sensitivity of
// erf's final `1 - 2^poly` combine.
#[inline(always)]
fn erf_poly(x: f32, x2: f32) -> f32 {
    let a6 = 2.8388531e-4f32;
    let a5 = -4.4954885e-3f32;
    let a4 = 3.2736249e-2f32;
    let a3 = -1.5164591e-1f32;
    let a2 = -9.1713983e-1f32;
    let a1 = -1.6281782f32;
    let a0 = 2.2989703e-5f32;
    let x4 = x2 * x2;
    let b0 = fma(a1, x, a0);
    let b1 = fma(a3, x, a2);
    let b2 = fma(a5, x, a4);
    let c0 = fma(b1, x2, b0);
    let c1 = fma(a6, x2, b2);
    fma(c1, x4, c0)
}

/// Error function `erf(x)`.
#[doc(alias = "erff")]
#[inline(always)]
pub fn erf(x: f32) -> f32 {
    let xa = x.abs();
    let xa_bounded = if xa > 10.0 { 10.0 } else { xa };
    // Shared between both branches: the Pade arm's own x2 (x*x, unclamped) and
    // erf_poly's internal x2 (xa_bounded*xa_bounded) only differ once |x| > 10,
    // but the Pade arm's result (`a`) is only ever *selected* when `xa < 0.28`,
    // comfortably inside the clamp -- so reusing the already-clamped x2 here
    // changes nothing observable, just removes a redundant multiply (and bounds
    // the discarded arm's x2 to <= 100 instead of letting it run up toward
    // overflow for huge |x|, a minor side benefit, not the point of the
    // change).
    let x2 = xa_bounded * xa_bounded;
    // The Pade numerator's constant term is `erf'(0) = 2/sqrt(pi)`, pinned by
    // the approximant rather than fitted -- and it lands almost exactly on an
    // f32 tie, 0.49 ulp above the nearer representable value. A single-word
    // constant therefore carries that 0.49 ulp as *bias*, not noise, into every
    // result small enough that `A*x2` has vanished beneath it: `x *
    // fl(2/sqrt(pi))` is systematically one ulp high over roughly half the f32
    // domain by bit pattern.
    let numer = fma(
        x,
        f32::from_bits(0x3f906ebb),
        x * fma(f32::from_bits(0x3f174f6e), x2, f32::from_bits(0xb37bd649)),
    );
    let denom = fma(
        fma(f32::from_bits(0x3e3e2be3), x2, f32::from_bits(0x3f5b6db7)),
        x2,
        1.0,
    );
    let a = numer / denom;
    let b = mulsign(1.0 - exp2(erf_poly(xa_bounded, x2)), x);
    if xa < 0.28 {
        a
    } else {
        b
    }
}

// The `|x|` clamps `erfc` and `erfcx` feed their `x*x` through. Both are set by
// one rule: make the exponent they hand `exp_reduce!` provably inside its valid
// range, so neither pays for a clamp `exp_checked` would have applied.
const ERFC_XS_CLAMP: f32 = 10.21;
const _: () = assert!((ERFC_XS_CLAMP as f64) * (ERFC_XS_CLAMP as f64) <= -(EXP_CLAMP_LO as f64));
// 150*ln2: below `e^-p` is at most half the smallest denormal, so it
// rounds to 0 and every clamped input returns the exactly-right answer.
const _: () = assert!((ERFC_XS_CLAMP as f64) * (ERFC_XS_CLAMP as f64) > 103.97207708399179);

// `erfcx`'s negative arm needs `2*e^(xs^2)` to still overflow to +inf
// (its true value out there), so `xs^2` must stay at or below
// `EXP_CLAMP_HI` while landing above `ln(f32::MAX/2)`.
const ERFCX_XS_CLAMP: f32 = 9.41;
const _: () = assert!((ERFCX_XS_CLAMP as f64) * (ERFCX_XS_CLAMP as f64) <= EXP_CLAMP_HI as f64);
// ln(f32::MAX/2) = 88.0296...: above this, `2*e^(xs^2)` leaves f32.
const _: () = assert!((ERFCX_XS_CLAMP as f64) * (ERFCX_XS_CLAMP as f64) > 88.02969187150839);

// erfcx(xa) for xa >= 0, shared by `erfc` and `erfcx` -- the scaled
// complementary error function, i.e. erfc(xa)*exp(xa^2).
#[inline(always)]
fn erfcx_pos(xa: f32) -> f32 {
    let v0 = 1.0 / (2.0 + xa);
    let v = if xa <= 2.0 {
        fma(-0.5 * xa, v0, 0.5)
    } else {
        v0
    };
    // c[0] is the *low* word of a two-word `1/sqrt(pi)`; the high word is
    // peeled out of the polynomial and applied in the final fma below.
    let c: [f32; 11] = [
        f32::from_bits(0xb2fbd649),
        1.1283773,
        1.974964,
        2.8070478,
        2.9756768,
        -3.7488432,
        17.02367,
        -117.490135,
        255.59447,
        -243.95302,
        90.238,
    ];
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
    // `erfcx(x) -> 1/(x*sqrt(pi))` as x grows, so out in the tail this whole
    // polynomial has collapsed to its constant term and the result is `v *
    // 1/sqrt(pi)` and nothing else. 1/sqrt(pi) sits 0.49 ulp above the nearer
    // f32 (a near-tie: the other side is 0.51 ulp below), so *either*
    // single-word choice leaves ~0.5 ulp of the constant as pure bias -- ~0.63
    // ulp of it in the result, in one direction, for every x past ~10.
    fma(
        v,
        f32::from_bits(0x3f106ebb),
        v * fma(fma(t8, v4, hi), v4, lo),
    )
}

/// Complementary error function `erfc(x) = 1 - erf(x)`.
#[doc(alias = "erfcf")]
#[inline(always)]
pub fn erfc(x: f32) -> f32 {
    // The x<0 reflection with no compare and no select. `w` shifts x's sign bit
    // straight into the exponent field (0x4000_0000 is 2.0f), and `mulsign`
    // applies the same bit to the Gaussian factor, so the two arms are `e*t +
    // 0` and `-e*t + 2`, in integer ops on ports the polynomial's fmas are not
    // contending for.
    let w = f32::from_bits((x.to_bits() >> 1) & 0x4000_0000);
    let xa = x.abs();
    let xs = if xa > ERFC_XS_CLAMP {
        ERFC_XS_CLAMP
    } else {
        xa
    };
    let p = xs * xs;
    let pe = fma(xs, xs, -p);
    let r = erfcx_pos(xa);
    // No clamp on the exponent at all: `ERFC_XS_CLAMP` is chosen so `-p` cannot
    // leave `exp_reduce!`'s valid range, which the two const assertions beside
    // it check. NaN reaches here as NaN (`NaN > c` is false), and comes back
    // out through `r`.
    let e = exp_reduce!(-p);
    fma(mulsign(e, x), fma(-r, pe, r), w)
}

// erfinv's central branch: erfinv(x) = x*P(x^2) for |x| <= 0.7, degree 8 in
// u=x^2 (Estrin-grouped, same idiom as erf_poly). Coefficients are a real
// least-squares fit (scipy) of erfinv(x)/x against u over [0,0.7], not a
// transcription of any published algorithm's constants.
#[inline(always)]
fn erfinv_central_poly_m1(u: f32) -> f32 {
    let c: [f32; 9] = [
        -0.11377308,
        2.3201263e-1,
        1.2761366e-1,
        8.534858e-2,
        7.737176e-2,
        -1.8498806e-2,
        2.6768243e-1,
        -3.5444248e-1,
        3.377774e-1,
    ];
    let u2 = u * u;
    let u4 = u2 * u2;
    let l0 = fma(c[1], u, c[0]);
    let l1 = fma(c[3], u, c[2]);
    let l2 = fma(c[5], u, c[4]);
    let l3 = fma(c[7], u, c[6]);
    let r0 = fma(l1, u2, l0);
    let r1 = fma(l3, u2, l2);
    fma(c[8], u4 * u4, fma(r1, u4, r0))
}

// erfinv's tail branch: `erfinv(x) = sign(x)*sqrt(w)*Q`, `w = -ln(1-x^2)`, for
// `|x| > 0.7`. An ulp-weighted minimax (LP) fit of `erfinv(x)/sqrt(w)`,
// weighted by the `sqrt(w)*Q` combine's own sensitivity `sqrt(w)/ulp(erfinv)`
// and then coordinate-descended over the f32 quantisation -- a plain
// least-squares fit left ~3.4x more idealized error than the same degree can
// reach.
#[inline(always)]
fn erfinv_tail_poly_m1(t: f32) -> f32 {
    let c: [f32; 12] = [
        -0.103669584,
        0.019397417,
        0.007896575,
        -0.0021255405,
        -0.0007704421,
        -1.969348e-5,
        -2.4305053e-5,
        0.00033242995,
        -0.00024495937,
        7.880177e-5,
        -1.2436317e-5,
        7.908929e-7,
    ];
    let t2 = t * t;
    let t4 = t2 * t2;
    let l0 = fma(c[1], t, c[0]);
    let l1 = fma(c[3], t, c[2]);
    let l2 = fma(c[5], t, c[4]);
    let l3 = fma(c[7], t, c[6]);
    let l4 = fma(c[9], t, c[8]);
    let l5 = fma(c[11], t, c[10]);
    let r0 = fma(l1, t2, l0);
    let r1 = fma(l3, t2, l2);
    let r2 = fma(l5, t2, l4);
    fma(fma(r2, t4, r1), t4, r0)
}

// erfinv's *far*-tail branch, a poly in `t = sqrt(w) - 7`. Reachable only from
// `erfc_inv`/`probit`, which build `w = -ln(1-x^2)` out of their own argument
// instead of out of `x`, and so reach `w` up to ~102.6 -- far past the `w <=
// 15.9424` that any 24-bit `x` can encode, which is all `erfinv_tail_poly_m1`
// is fitted for.
#[inline(always)]
fn erfinv_far_poly_m1(t: f32) -> f32 {
    let c: [f32; 8] = [
        -0.018711485,
        0.0039010558,
        -0.00063050824,
        9.089283e-5,
        -1.2126031e-5,
        1.5352017e-6,
        -1.7609956e-7,
        1.2450847e-8,
    ];
    let t2 = t * t;
    let t4 = t2 * t2;
    let l0 = fma(c[1], t, c[0]);
    let l1 = fma(c[3], t, c[2]);
    let l2 = fma(c[5], t, c[4]);
    let l3 = fma(c[7], t, c[6]);
    let r0 = fma(l1, t2, l0);
    let r1 = fma(l3, t2, l2);
    fma(r1, t4, r0)
}

// Where `erfc_inv_half`'s tail hands over from `erfinv_tail_poly_m1` to
// `erfinv_far_poly_m1`: the largest `w` a 24-bit `erfinv` argument can produce
// is `-ln(1 - x_max^2) = 15.9424`, so the two polys split exactly where
// `erfinv`'s own reachable range stops and `erfc_inv`/`probit`'s extra reach
// begins. Keeping the split there is what leaves `erfinv_tail_poly_m1` -- and
// `erfinv` itself -- untouched by this.
const ERFC_INV_W_FAR: f32 = 16.0;

// `|erfc_inv(n)|` for `n` in `(0, 1]`, the half both erfc_inv and probit reduce
// to, and the reason neither is written as `erfinv(1-y)` any more.
#[inline(always)]
fn erfc_inv_half(n: f32) -> f32 {
    let x = 1.0 - n;
    let central = fma(x, erfinv_central_poly_m1(x * x), x);
    let s = fma(-n, n, n + n);
    // `denormal_rescale!` and nothing else from `log_family_wrapper!`: `s`
    // reaches down to `2 * f32::MIN_POSITIVE_SUBNORMAL` for the smallest `n`,
    // so the rescale is live, but the zero/negative/inf/NaN arms are all dead
    // here -- the `n > 0.0` select below already owns every input that could
    // reach them ("the guard is the licence").
    let (ss, koff) = denormal_rescale!(s);
    let w = -ln_normal(ss, koff);
    let v = w.sqrt();
    let q = if w > ERFC_INV_W_FAR {
        erfinv_far_poly_m1(v - 7.0)
    } else {
        erfinv_tail_poly_m1(v - 1.0)
    };
    let mag = if x <= 0.7 { central } else { fma(v, q, v) };
    // `n == 0` is the pole and `n < 0` (with NaN) is a domain error. Both have
    // to be pinned rather than left to fall out: `s == 0` puts `ln_normal` off
    // its own positive-normal contract, and it returns a large finite number
    // there instead of the `+inf` the pole needs.
    let edge = if n == 0.0 { f32::INFINITY } else { f32::NAN };
    if n > 0.0 {
        mag
    } else {
        edge
    }
}

/// Inverse error function for `x` in `(-1, 1)`.
#[inline(always)]
pub fn erfinv(x: f32) -> f32 {
    let ax = x.abs();
    // `1 - x^2` factored, never subtracted: `(1-|x|)(1+|x|) = n*(2-n)`, the
    // same reduction `erfc_inv_half` uses and for the same reason. `n = 1-|x|`
    // is Sterbenz-exact for every `|x|` the tail is selected for (`|x| > 0.7`),
    // `n+n` is an exact scaling, and `fma` makes `2n - n^2` a single rounding
    // of the whole product -- so `s` carries one `2^-25` *relative* error at
    // any `n`.
    let n = 1.0 - ax;
    let s = fma(-n, n, n + n);
    let w = -if s > 0.0 { ln_normal(s, 0.0) } else { f32::NAN };
    // Both arms are built on `ax` and take `x`'s sign once, on the merged
    // select, rather than the central arm carrying a signed `x` through the
    // poly. The peeled `fma(x, P-1, x)` form cannot reproduce `erfinv(-0.0)`:
    // `P-1` is negative, so `-0.0 * (P-1)` is `+0.0` and `+0.0 + -0.0` is
    // `+0.0` under round-to-nearest -- losing the sign of zero that plain `x *
    // P` carried for free.
    let central = fma(ax, erfinv_central_poly_m1(x * x), ax);
    let v = w.sqrt();
    let tail = fma(v, erfinv_tail_poly_m1(v - 1.0), v);
    let normal = mulsign(if ax <= 0.7 { central } else { tail }, x);
    if ax == 1.0 {
        f32::INFINITY.copysign(x)
    } else {
        normal
    }
}

// `norm_cdf`'s counterpart to `ERFC_XS_CLAMP`, in `x`'s units rather than
// `x/sqrt(2)`'s. Same two conditions, and they are asserted the same way:
// `x^2/2` must stay inside `exp_reduce!`'s valid range, and must already be far
// enough below zero that `e^-x^2/2` rounds to exactly 0 -- the true `norm_cdf`
// has reached exactly 0.0f32 by `x ~ -14.2`.
const NORM_CDF_XS_CLAMP: f32 = 14.44;
const _: () = assert!(
    (NORM_CDF_XS_CLAMP as f64) * (NORM_CDF_XS_CLAMP as f64) * 0.5 <= -(EXP_CLAMP_LO as f64)
);
const _: () =
    assert!((NORM_CDF_XS_CLAMP as f64) * (NORM_CDF_XS_CLAMP as f64) * 0.5 > 103.97207708399179);

/// Standard normal cumulative distribution function `Phi(x) = 0.5 * erfc(-x / sqrt(2))`.
#[inline(always)]
pub fn norm_cdf(x: f32) -> f32 {
    let xa = x.abs();
    // `erfcx_pos`, not `erfcx`: the argument is an absolute value, so erfcx's
    // own x<0 arm is dead, and LLVM does not prove that -- taking the public
    // wrapper left a whole second `exp_reduce!` (two `exp2_field_split`s in the
    // asm) in the region. x/sqrt(2) is formed here and nowhere else, and erfcx
    // is well-conditioned in it.
    let r = erfcx_pos(xa * std::f32::consts::FRAC_1_SQRT_2);
    let xs = if xa > NORM_CDF_XS_CLAMP {
        NORM_CDF_XS_CLAMP
    } else {
        xa
    };
    // p + pe == xs^2/2 exactly: the halving is exact, so this is just
    // `two_prod`'s error term on a single multiply. The clamp above is
    // what keeps `pe` from becoming `inf*inf - inf == NaN` at x = +-inf.
    let h = 0.5 * xs;
    let p = h * xs;
    let pe = fma(h, xs, -p);
    let e = exp_reduce!(-p);
    // `-x` as a sign carrier only (Φ(x) uses erfc(-x/sqrt(2))), so the
    // negation is a bit flip and stays exact at +-0.0: x = +0.0 takes the
    // `1 - y/2` arm and x = -0.0 the `y/2` arm, and both are exactly 0.5.
    let nx = -x;
    let w = f32::from_bits((nx.to_bits() >> 1) & 0x4000_0000);
    0.5 * fma(mulsign(e, nx), fma(-r, pe, r), w)
}

/// Inverse complementary error function for `x` in `(0, 2)`.
#[inline(always)]
pub fn erfc_inv(y: f32) -> f32 {
    let n = if y < 1.0 { y } else { 2.0 - y };
    mulsign(erfc_inv_half(n), 1.0 - y)
}

/// Probit (standard normal quantile function) for `p` in `(0, 1)`.
#[inline(always)]
pub fn probit(p: f32) -> f32 {
    let m = if p < 0.5 { p } else { 1.0 - p };
    mulsign(std::f32::consts::SQRT_2 * erfc_inv_half(m + m), p - 0.5)
}

/// Standard normal probability density function `phi(x) = exp(-x^2 / 2) / sqrt(2*pi)`.
#[inline(always)]
pub fn norm_pdf(x: f32) -> f32 {
    // `1/sqrt(2*pi)` as a double-`f32` pair, same shape as
    // `FRAC_1_PI`/`RPI_LO`. The single-word constant is correctly rounded and
    // still sits 0.48 ulp above the true value, which the final multiply hands
    // straight to the result as a 0.24-0.48 ulp bias -- there is nothing else
    // in the chain to cancel it, unlike `sinc`, where the same constant appears
    // on both sides of a ratio.
    const INV_SQRT_2PI_HI: f32 = 0.3989423;
    const INV_SQRT_2PI_LO: f32 = -1.133517e-8;
    let xa = x.abs();
    let xs = if xa > NORM_CDF_XS_CLAMP {
        NORM_CDF_XS_CLAMP
    } else {
        xa
    };
    let h = 0.5 * xs;
    let p = h * xs;
    let pe = fma(h, xs, -p);
    let e = exp_reduce!(-p);
    let y = fma(-pe, e, e);
    fma(y, INV_SQRT_2PI_HI, y * INV_SQRT_2PI_LO)
}

// dawson's central branch: `x*P(u)/Q(u)`, `u=x^2`, degree 6/5, Estrin-grouped
// with each side's top group folded in at the `u^2` level so `u^4` is never
// formed (exp_r_poly!'s own fold, applied to the numerator and the denominator
// alike).
#[inline(always)]
fn dawson_central_ratio(u: f32) -> f32 {
    // `ac`/`bc` are the fitted numerator and denominator with their pinned
    // `1.0` constant terms *removed*, so `P = 1 + u*A` and `Q = 1 + u*B`. Same
    // polynomial and same coefficient values as an unpeeled `pc`/`qc` pair;
    // only where the `1.0` enters the evaluation changes.
    let ac: [f32; 6] = [
        -0.085751414,
        0.037434783,
        -0.0004054072,
        0.00019858626,
        5.6392253e-8,
        2.8313497e-8,
    ];
    let bc: [f32; 6] = [
        0.5809171,
        0.15803601,
        0.026255792,
        0.002856104,
        0.00021274923,
        5.002549e-6,
    ];
    let u2 = u * u;
    // Each side's top coefficient (degree 6 in `u` overall) rides into the
    // existing `u^2` group rather than needing a `u^6` of its own -- the
    // same trick `ln_normal`'s c[8] uses, one fma and no new multiply.
    let al0 = fma(ac[1], u, ac[0]);
    let al1 = fma(ac[3], u, ac[2]);
    let al2 = fma(ac[5], u, ac[4]);
    let ar0 = fma(al2, u2, al1);
    let a = fma(ar0, u2, al0);
    let num = fma(u, a, 1.0);
    let bl0 = fma(bc[1], u, bc[0]);
    let bl1 = fma(bc[3], u, bc[2]);
    let bl2 = fma(bc[5], u, bc[4]);
    let br0 = fma(bl2, u2, bl1);
    let b = fma(br0, u2, bl0);
    let den = fma(u, b, 1.0);
    num / den
}

// dawson's tail branch: `w * R(z)` with the leading term peeled, i.e. `w +
// w*z*T(z)`, where `w = 1/(2x)`, `z = w^2 = 1/(4x^2)` and `R(z) = 1 + z*T(z)`.
#[inline(always)]
fn dawson_tail(w: f32) -> f32 {
    let c: [f32; 4] = [2.0000212, 11.969649, 131.78029, 116698.56];
    let z = w * w;
    let z2 = z * z;
    let t0 = fma(c[1], z, c[0]);
    let t1 = fma(c[3], z2, c[2]);
    let t = fma(t1, z2, t0);
    fma(w * z, t, w)
}

/// Dawson's integral `F(x) = exp(-x^2) * integral_0^x exp(t^2) dt`.
#[doc(alias = "dawsn")]
#[inline(always)]
pub fn dawson(x: f32) -> f32 {
    let u = x * x;
    let central = x * dawson_central_ratio(u);
    let w = 0.5 / x;
    let tail = dawson_tail(w);
    if x.abs() <= 4.0 {
        central
    } else {
        tail
    }
}

/// Logit function: `ln(p / (1 - p))` for `p` in `(0, 1)`.
#[inline(always)]
pub fn logit(p: f32) -> f32 {
    // `2p-1`, exact for every `p >= 0.25` (`2p` is an exact scaling,
    // Sterbenz covers the subtraction), which is what lets the central
    // arm keep full relative accuracy at the zero.
    let a = fma(2.0, p, -1.0);
    // logit(p) = 2*atanh(2p-1), reusing `atanh`'s own small-argument polynomial
    // on its own `|x| < 0.25` domain. No subtraction of logarithms, so nothing
    // cancels: the result is proportional to `a`, and `a` is exact.
    let central = 2.0 * atanh_small(a);
    // Outside that band the difference form has nothing left to cancel (at the
    // seam the result is already ~0.51 against operands ~0.69) and it is what
    // keeps `p` denormal, `0`, `1` and out-of-domain correct, so it stays
    // exactly as it was.
    let np = -p;
    let t = 1.0 + np;
    let c = np - (t - 1.0);
    // No `corr.is_finite()` guard: the only `p` making it non-finite are `p ==
    // 1` (`t == 0`, which takes the `spec` arm) and `p == -inf` (where `ln(p)`
    // is already `NaN`, so the subtraction below is `NaN` either way).
    let corr = c / t;
    let spec = if t == 0.0 {
        f32::NEG_INFINITY
    } else {
        f32::NAN
    };
    let l = if t > 0.0 {
        ln_normal(t, 0.0) + corr
    } else {
        spec
    };
    let outer = ln(p) - l;
    if a.abs() < 0.25 {
        central
    } else {
        outer
    }
}

/// Computes `x * ln(y)`, with `0 * ln(y) = 0`.
#[inline(always)]
pub fn xlogy(x: f32, y: f32) -> f32 {
    let normal = x * ln(y);
    if x == 0.0 {
        0.0
    } else {
        normal
    }
}

/// Computes `x * ln(1 + y)`, with `0 * ln(1 + y) = 0`.
#[inline(always)]
pub fn xlog1py(x: f32, y: f32) -> f32 {
    let normal = x * log1p(y);
    if x == 0.0 {
        0.0
    } else {
        normal
    }
}

/// Computes `(1 + x)^n`.
#[inline(always)]
pub fn compound(x: f32, n: f32) -> f32 {
    // `log1p` minus its trailing signed-zero select: that select only changes
    // `log1p(-0.0)` from `+0.0` to `-0.0`, and `exp_checked` maps both zeros to
    // exactly `1.0`, so the sign never reaches the result.
    exp_checked(n * log1p_nonzero!(x))
}

// log2(1+x) in f64, `compound_accurate`'s exponent, to the same ~31 bits
// `log2_f64` has to reach and for the same reason: the result is about to be
// multiplied by an unrestricted `n`.
#[inline(always)]
fn log2p1_f64(x: f32) -> f64 {
    let xd = x as f64;
    let d = 1.0 + xd;
    // Same decomposition `log2_f64` does, one format wider: m in
    // [2^-0.5, 2^0.5), k exact. No denormal rescale is needed at any
    // width here -- `d` is either 0 (handled below) or at least 2^-24.
    let bits = d.to_bits() as i64;
    let ki = (bits - 0x3fe6a09e667f3bcd) >> 52;
    let m = f64::from_bits((bits - (ki << 52)) as u64);
    let k = ki as f64;
    let s = if xd.abs() < 0.25 { xd } else { m - 1.0 };
    let dd = m + 1.0;
    let rc = LOG2_ATANH_RCP64;
    let md2 = m * m;
    let e0 = f64::mul_add(rc[1], m, rc[0]);
    let e1 = f64::mul_add(rc[3], m, rc[2]);
    let e2 = f64::mul_add(rc[5], m, rc[4]);
    let r = f64::mul_add(e2, md2 * md2, f64::mul_add(e1, md2, e0));
    let r = r * f64::mul_add(-dd, r, 2.0);
    let t = s * r;
    let u = t * t;
    let a = LOG2_ATANH_A64;
    let u2 = u * u;
    let l0 = f64::mul_add(a[1], u, a[0]);
    let l1 = f64::mul_add(a[3], u, a[2]);
    let l2 = f64::mul_add(l1, u2, l0);
    let q = f64::mul_add(LOG2E_2_F64 * t, f64::mul_add(l2, u, 1.0), k);
    // Degenerate `d` routes around the bit-level decomposition, exactly as
    // `powf_f64_mag!` does: `d == 0` is `x == -1`, exactly `0^n`; `d < 0` is `x
    // < -1`, where the real power does not exist; `d` non-finite is `x`
    // non-finite. `exp2_f64_to_f32`'s own clamp turns each into the right
    // saturation.
    let deg = if d == 0.0 { f64::NEG_INFINITY } else { d };
    let deg = if d < 0.0 { f64::NAN } else { deg };
    if d > 0.0 && d < f64::INFINITY {
        q
    } else {
        deg
    }
}

/// Accurate `(1 + x)^n` using f64 intermediate computation.
#[inline(always)]
pub fn compound_accurate(x: f32, n: f32) -> f32 {
    exp2_f64_to_f32(log2p1_f64(x) * n as f64)
}

/// Scaled complementary error function `erfcx(x) = exp(x^2) * erfc(x)`.
#[inline(always)]
pub fn erfcx(x: f32) -> f32 {
    let xa = x.abs();
    let r = erfcx_pos(xa);
    let xs = if xa > ERFCX_XS_CLAMP {
        ERFCX_XS_CLAMP
    } else {
        xa
    };
    let p = xs * xs;
    let pe = fma(xs, xs, -p);
    let g = exp_reduce!(p);
    if x >= 0.0 {
        r
    } else {
        fma(g, fma(pe, 2.0, 2.0), -r)
    }
}

/// Computes `1 / sqrt(x)`.
#[inline(always)]
pub fn rsqrt(x: f32) -> f32 {
    1.0 / x.sqrt()
}

/// Computes `sqrt(x^2 + y^2)`.
#[doc(alias = "hypotf")]
#[inline(always)]
pub fn hypot(x: f32, y: f32) -> f32 {
    let normal = fma(x, x, y * y).sqrt();
    // hypot(+-inf, anything) and hypot(anything, +-inf) = +inf, even when the
    // other argument is NaN -- IEEE754/C99 special-cases infinity to "win" over
    // NaN here (unlike almost every other function). The naive formula can't
    // reach this on its own: once either argument actually is NaN, `inf*inf +
    // NaN*NaN` degrades to NaN instead.
    if x.is_infinite() || y.is_infinite() {
        f32::INFINITY
    } else {
        normal
    }
}

/// `hypot` with anti-overflow/underflow scaling.
#[inline(always)]
pub fn hypot_checked(x: f32, y: f32) -> f32 {
    let ax = x.abs();
    let ay = y.abs();
    let m = ax.max(ay);
    let is_zero = ax == 0.0 && ay == 0.0;
    let m_safe = if is_zero { 1.0 } else { m };
    let tiny = m_safe < f32::MIN_POSITIVE;
    let pre = if tiny { 16777216.0 } else { 1.0 }; // 2^24
    let post = if tiny { 1.0 / 16777216.0 } else { 1.0 };
    let ax = ax * pre;
    let ay = ay * pre;
    let m_safe = m_safe * pre;
    let e = ((m_safe.to_bits() >> 23) as i32) - 127;
    let es = 2 * (e >> 1);
    let scale = f32::from_bits(((127 - es) as u32) << 23);
    let descale = f32::from_bits(((127 + es) as u32) << 23);
    let xs = ax * scale;
    let ys = ay * scale;
    let normal = fma(xs, xs, ys * ys).sqrt() * descale * post;
    let normal = if is_zero { 0.0 } else { normal };
    if x.is_infinite() || y.is_infinite() {
        f32::INFINITY
    } else {
        normal
    }
}

/// Computes `1 / hypot(x, y)`.
#[inline(always)]
pub fn rhypot(x: f32, y: f32) -> f32 {
    let normal = 1.0 / fma(x, x, y * y).sqrt();
    if x.is_infinite() || y.is_infinite() {
        0.0
    } else {
        normal
    }
}

/// Computes `sqrt(x^2 + y^2 + z^2)`.
#[inline(always)]
pub fn hypot3(x: f32, y: f32, z: f32) -> f32 {
    let normal = fma(x, x, fma(y, y, z * z)).sqrt();
    if x.is_infinite() || y.is_infinite() || z.is_infinite() {
        f32::INFINITY
    } else {
        normal
    }
}

/// Computes `1 / sqrt(x^2 + y^2 + z^2)`.
#[inline(always)]
pub fn rnorm3(x: f32, y: f32, z: f32) -> f32 {
    let normal = 1.0 / fma(x, x, fma(y, y, z * z)).sqrt();
    if x.is_infinite() || y.is_infinite() || z.is_infinite() {
        0.0
    } else {
        normal
    }
}

/// Normalizes a 2D vector `(x, y)`.
#[inline(always)]
pub fn normalize2(x: f32, y: f32) -> (f32, f32) {
    let r = rhypot(x, y);
    (x * r, y * r)
}

/// Normalizes a 3D vector `(x, y, z)`.
#[inline(always)]
pub fn normalize3(x: f32, y: f32, z: f32) -> (f32, f32, f32) {
    let r = rnorm3(x, y, z);
    (x * r, y * r, z * r)
}

/// Computes `sqrt(w^2 + x^2 + y^2 + z^2)`.
#[inline(always)]
pub fn hypot4(w: f32, x: f32, y: f32, z: f32) -> f32 {
    let normal = fma(w, w, fma(x, x, fma(y, y, z * z))).sqrt();
    let any_inf = w.is_infinite() || x.is_infinite() || y.is_infinite() || z.is_infinite();
    if any_inf {
        f32::INFINITY
    } else {
        normal
    }
}

/// Computes `1 / sqrt(w^2 + x^2 + y^2 + z^2)`.
#[inline(always)]
pub fn rnorm4(w: f32, x: f32, y: f32, z: f32) -> f32 {
    let normal = 1.0 / fma(w, w, fma(x, x, fma(y, y, z * z))).sqrt();
    let any_inf = w.is_infinite() || x.is_infinite() || y.is_infinite() || z.is_infinite();
    if any_inf {
        0.0
    } else {
        normal
    }
}

/// Normalizes a 4D quaternion/vector `(w, x, y, z)`.
#[inline(always)]
pub fn normalize4(w: f32, x: f32, y: f32, z: f32) -> (f32, f32, f32, f32) {
    let r = rnorm4(w, x, y, z);
    (w * r, x * r, y * r, z * r)
}

/// Computes `a*b - c*d` using Kahan's compensated algorithm.
#[inline(always)]
pub fn diff_of_products(a: f32, b: f32, c: f32, d: f32) -> f32 {
    let w = c * d;
    let e = fma(-c, d, w);
    let f = fma(a, b, -w);
    f + e
}

/// 2D cross product (`ax*by - ay*bx`) via Kahan's algorithm.
#[inline(always)]
pub fn cross2(ax: f32, ay: f32, bx: f32, by: f32) -> f32 {
    diff_of_products(ax, by, ay, bx)
}

/// Complex modulus `|re + im * i|`.
#[doc(alias = "cabsf")]
#[inline(always)]
pub fn cabs(re: f32, im: f32) -> f32 {
    hypot_checked(re, im)
}

/// Complex argument (principal value in `(-pi, pi]`).
#[doc(alias = "cargf")]
#[inline(always)]
pub fn carg(re: f32, im: f32) -> f32 {
    atan2(im, re)
}

/// Complex exponential `e^(re + im * i) = e^re * (cos(im) + i * sin(im))`.
#[inline(always)]
pub fn cexp(re: f32, im: f32) -> (f32, f32) {
    let m = exp(re);
    (m * cos(im), m * sin(im))
}

/// Complex natural log `ln(re+im*i) = ln(|re+im*i|) + i*arg(re+im*i)`, returned
/// as `(re, im)`.
#[inline(always)]
pub fn clog(re: f32, im: f32) -> (f32, f32) {
    // `log1p(v)` for a `v` that is known to keep `1+v` positive, normal and
    // finite -- see this function's doc comment. Bit-identical to `log1p` over
    // every f32 in both call sites' licensed ranges (`|v| < 0.5` and `[0, 1]`,
    // checked exhaustively) with exactly one exception, `v == -0.0`, where
    // `log1p`'s signed-zero select returns `-0.0` and this returns `+0.0`.
    #[inline(always)]
    fn log1p_guarded(v: f32) -> f32 {
        let u = 1.0 + v;
        let c = v - (u - 1.0);
        ln_normal(u, 0.0) + c / u
    }
    let mag = cabs(re, im);
    let log_mag = if mag.is_finite() {
        if (mag - 1.0).abs() < 0.5 {
            // `re^2 + im^2 - 1` to full *relative* precision however hard it
            // cancels, in three f64 operations and with no error-free transform
            // at all.
            let are = re.abs();
            let aim = im.abs();
            let a = are.max(aim) as f64;
            let b = are.min(aim) as f64;
            let v = f64::mul_add(a, a, -1.0) + b * b;
            0.5 * log1p_guarded(v as f32)
        } else {
            ln(mag)
        }
    } else if re.is_finite() && im.is_finite() {
        let are = re.abs();
        let aim = im.abs();
        let (mx, mn) = if are > aim { (are, aim) } else { (aim, are) };
        let ratio = mn / mx;
        ln(mx) + 0.5 * log1p_guarded(ratio * ratio)
    } else {
        ln(mag)
    };
    (log_mag, carg(re, im))
}

/// `2*log2(e)`, the atanh form's leading coefficient (see `log2_f64`).
const LOG2E_2_F64: f64 = 2.8853900817779268;

/// Minimax seed for `1/(m+1)` over `m` in `[2^-0.5, 2^0.5]`, accurate to
/// ~`2^-20` -- only a seed, squared by the single Newton step that follows it,
/// so it does not need to be better. Fitted in `m` rather than in `d = m + 1`
/// (the same fit either way, an affine change of variable) so the seed does not
/// have to wait on the `m + 1` add.
const LOG2_ATANH_RCP64: [f64; 6] = [
    0.9836614733399011,
    -0.8856304007709652,
    0.6435262038352056,
    -0.3284692378061808,
    0.10056909996623936,
    -0.013656925651166552,
];

/// `(atanh(t)/t - 1)/u` in `u = t^2` over `u` in `[0, (3-2*sqrt(2))^2]`, i.e.
/// the atanh series past its own leading term, with the leading `1` pinned
/// rather than fitted: that makes `log2(1)` come out exactly `0`, which
/// `powf_unchecked` (which has no `x == 1` override) relies on.
const LOG2_ATANH_A64: [f64; 4] = [
    0.33333332824327616,
    0.20000167265984317,
    0.14268673572031404,
    0.117907343543545,
];

/// `(2^f - 1)/f` over `f` in `[-0.5, 0.5]`, leading `1` pinned so `2^0` is
/// exactly `1`. Idealized relative error `2^-28.9`, against the `2^-25` needed
/// for an f32 result.
const EXP2_F64_E: [f64; 6] = [
    0.6931472028549269,
    0.24022647913384074,
    0.05550332471225973,
    0.009618437395496837,
    0.0013398874430087457,
    0.0001535334944368378,
];

/// `log2(x)` in f64, for positive finite `x` (denormals included; callers must
/// guard zero/negative/inf/nan themselves). The `powf` family's log half.
#[inline(always)]
fn log2_f64(x: f32) -> f64 {
    let (xs, koff) = denormal_rescale!(x);
    // Same decomposition log_family_normal! does, spelled out (like
    // ln_normal/log10_normal's own copies): m in [2^-0.5, 2^0.5), k exact.
    let e = (xs.to_bits() as i32).wrapping_sub(0x3f3504f3) >> 23;
    let m = f32::from_bits((xs.to_bits() as i32).wrapping_sub(e << 23) as u32);
    let k = (e as f32 + koff) as f64;
    let md = m as f64;
    // Both exact: `md - 1` by Sterbenz, `md + 1` because 24 bits plus a
    // leading one still fits. So the whole reduction below is exact
    // except for the reciprocal itself.
    let s = md - 1.0;
    let d = md + 1.0;
    let rc = LOG2_ATANH_RCP64;
    let md2 = md * md;
    // Estrin, not Horner: this sits at the head of the chain everything
    // else waits on, so a level of depth is worth an extra multiply.
    let e0 = f64::mul_add(rc[1], md, rc[0]);
    let e1 = f64::mul_add(rc[3], md, rc[2]);
    let e2 = f64::mul_add(rc[5], md, rc[4]);
    let r = f64::mul_add(e2, md2 * md2, f64::mul_add(e1, md2, e0));
    let r = r * f64::mul_add(-d, r, 2.0);
    let t = s * r;
    let u = t * t;
    let a = LOG2_ATANH_A64;
    let u2 = u * u;
    let l0 = f64::mul_add(a[1], u, a[0]);
    let l1 = f64::mul_add(a[3], u, a[2]);
    let l2 = f64::mul_add(l1, u2, l0);
    // `1 + u*A(u)` folded into the `t` multiply, so the pinned leading
    // term stays exact and no separate `+1` rounding happens.
    f64::mul_add(LOG2E_2_F64 * t, f64::mul_add(l2, u, 1.0), k)
}

/// `2^v` for an f64 `v`, narrowed to f32. The `powf` family's exp half.
#[inline(always)]
fn exp2_f64_to_f32(v: f64) -> f32 {
    let vc = v.clamp(-200.0, 200.0);
    // ROUND_MAGIC64 does double duty: `nm - MAGIC` is round-ties-even of `v`,
    // and `nm`'s own low bits already *are* that integer, so the scale's
    // exponent field costs an integer add and a shift rather than a
    // float-to-int cast (which is saturating in Rust and does not vectorize --
    // the signature `codegen_check` watches for).
    let nm = vc + ROUND_MAGIC64;
    let n = nm - ROUND_MAGIC64;
    // Exact: |f| <= 0.5 and both operands share an exponent range.
    let f = vc - n;
    let c = EXP2_F64_E;
    let f2 = f * f;
    let l0 = f64::mul_add(c[1], f, c[0]);
    let l1 = f64::mul_add(c[3], f, c[2]);
    let l2 = f64::mul_add(c[5], f, c[4]);
    let r0 = f64::mul_add(l1, f2, l0);
    let r1 = f64::mul_add(l2, f2 * f2, r0);
    let p = f64::mul_add(r1, f, 1.0);
    // (n + 1023) << 52. `nm`'s low 52 bits hold `n + 2^51`, and 2^51 is a
    // multiple of 4096, so the 12 bits the shift keeps are exactly
    // `n + 1023` -- in range for every `n` the clamp above allows.
    let scale = f64::from_bits(nm.to_bits().wrapping_add(1023) << 52);
    (p * scale) as f32
}

/// `exp2(log2(ax) * y)`, the magnitude half of the whole `powf` family.
macro_rules! powf_f64_mag {
    ($ax:expr, $y:expr) => {{
        let ax = $ax;
        let l = log2_f64(ax);
        let l = if ax == 0.0 { f64::NEG_INFINITY } else { l };
        let l = if !(ax < f32::INFINITY) {
            (ax * ax) as f64
        } else {
            l
        };
        exp2_f64_to_f32(l * ($y as f64))
    }};
}

// Shared by powf/rootn: the negative-base/y-parity/y==0/x==+-1 special-case
// combine, given each caller's own already-computed `mag`.
macro_rules! powf_sign_combine {
    ($x:expr, $ax:expr, $y:expr, $mag:expr) => {{
        let x = $x;
        let y = $y;
        let ax = $ax;
        let mag = if ax == 1.0 { 1.0 } else { $mag };
        // 0 for an even integer y, 1 for an odd one, neither otherwise --
        // and forced to the even case for any x that isn't sign-negative,
        // where no sign flip and no domain error can apply.
        let par = fma(-2.0, (y * 0.5).floor(), y);
        let par = if x.is_sign_negative() { par } else { 0.0 };
        // what a negative base with a *non*-integer y gives: a domain error,
        // except where |x| alone decides the answer -- `ax + ax == ax` picks
        // out exactly `+0` and `+inf` (the two C99-exempt boundary magnitudes)
        // in one compare, and infinite y is the third, independent exemption.
        // Infinite y specifically, not `y - y != 0.0`: that spelling is two
        // constants cheaper and looked free, because a NaN y makes `mag` NaN
        // anyway -- except at `ax == 1`, where `mag` is *pinned* to 1 just
        // above.
        let spec = if ax + ax == ax { 1.0 } else { f32::NAN };
        let spec = if y.abs() == f32::INFINITY { 1.0 } else { spec };
        let sm = if par == 0.0 { 1.0 } else { spec };
        let sm = if par == 1.0 { -1.0 } else { sm };
        let r = mag * sm;
        if y == 0.0 {
            1.0
        } else {
            r
        }
    }};
}

/// Computes `x^y` (C99 `pow` semantics).
#[doc(alias = "pow")]
#[inline(always)]
#[allow(clippy::neg_cmp_op_on_partial_ord)] // `!(ax < inf)` catches NaN too
pub fn powf(x: f32, y: f32) -> f32 {
    let ax = x.abs();
    let mag = powf_f64_mag!(ax, y);
    powf_sign_combine!(x, ax, y, mag)
}

/// `powf` for `x > 0` (or `x == +0.0`).
#[inline(always)]
#[allow(clippy::neg_cmp_op_on_partial_ord)] // `!(x < inf)` catches NaN too
pub fn powf_pos(x: f32, y: f32) -> f32 {
    let mag = powf_f64_mag!(x, y);
    let r = if x == 1.0 { 1.0 } else { mag };
    if y == 0.0 {
        1.0
    } else {
        r
    }
}

/// Signed power: `copysign(|x|^y, x)`.
#[inline(always)]
pub fn signed_pow(x: f32, y: f32) -> f32 {
    mulsign(powf_pos(x.abs(), y), x)
}

/// `powf` without domain or sign checks.
#[inline(always)]
pub fn powf_unchecked(x: f32, y: f32) -> f32 {
    exp2_f64_to_f32(log2_f64(x) * y as f64)
}

/// Computes `x^(1/n)` for integer `n` (C23 `rootn`). Correctly handles negative bases when `n` is odd.
#[inline(always)]
pub fn rootn(x: f32, n: i32) -> f32 {
    let ax = x.abs();
    // |x| = m * 2^e with m in [sqrt(2)/2, sqrt(2)) -- log_2_normal's own window
    // and its own bit trick, so log_2_unchecked(m) sees k == 0 and returns
    // log2(m) in [-0.5, 0.5) with no exponent to add back and no `lm + k`
    // rounding. Reaching for frexp instead costs its [0.5, 1) -> window shift
    // and its zero/infinite selects, which the `degenerate` arm below re-does
    // anyway.
    let (xs, koff) = denormal_rescale!(ax);
    let bits = xs.to_bits() as i32;
    let ew = bits.wrapping_sub(0x3f35_04f3) >> 23;
    let m = f32::from_bits(bits.wrapping_sub(ew << 23) as u32);
    let e = ew + koff as i32;
    // |n| >= 2 for the general path; n in {-1, 0, 1} takes the closed-form
    // arm below, so any stand-in works and 2 keeps both the integer
    // division and exp2_kf's exponent range well defined.
    let small = n.unsigned_abs() < 2;
    let nz = if small { 2 } else { n };
    let q = e.div_euclid(nz);
    let rr = e.rem_euclid(nz);
    let t = (rr as f32 + log_2_unchecked(m)) / (nz as f32);
    // t is in (-1, 1], so this only ever moves one octave; exp2_kf wants
    // the fraction separately anyway, and folding q into the same integer
    // is what keeps e's magnitude off the float path.
    let tk = t.floor();
    let mag_normal = exp2_kf(q as f32 + tk, t - tk);
    // `|x|` for `n > 0` and `1/|x|` for `n < 0` is the whole answer for three
    // separate reasons at once: it is the exact `|n| == 1` identity, and it is
    // also the right magnitude for zero, infinite and NaN `x` (where the
    // mantissa split reads nothing meaningful and only `n`'s sign matters). `n
    // == 0` lands here too and is discarded by the domain-error select at the
    // end.
    let degenerate = ax == 0.0 || !ax.is_finite();
    let mag = if degenerate || small {
        if n > 0 {
            ax
        } else {
            1.0 / ax
        }
    } else {
        mag_normal
    };
    let n_odd = n % 2 != 0;
    let signed = if n_odd { mulsign(mag, x) } else { mag };
    let neg_even_domain_error = x < 0.0 && !n_odd;
    let r = if neg_even_domain_error {
        f32::NAN
    } else {
        signed
    };
    if n == 0 {
        f32::NAN
    } else {
        r
    }
}

// `1/1.055` and `0.055/1.055` as the nearest `f32` to each *exact* value, so `b
// = (c+0.055)/1.055` is one fma with one rounding instead of an add and a
// multiply with two. Written out rather than left to the compiler: `1.0f32 /
// 1.055f32` divides by an already-rounded `1.055` and lands 0.53 ulp above the
// true `1/1.055`, which `^2.4` turns into a systematic 1.3 ulp of
// `srgb_to_linear`.
const SRGB_INV_1055: f32 = 0.9478672742843628;
const SRGB_OFF_1055: f32 = 0.05213269963860512;

/// Converts an sRGB color component in `[0, 1]` to linear (IEC 61966-2-1).
#[inline(always)]
#[allow(clippy::neg_cmp_op_on_partial_ord)]
pub fn srgb_to_linear(c: f32) -> f32 {
    let low = c * (1.0 / 12.92);
    let b = fma(c, SRGB_INV_1055, SRGB_OFF_1055);
    let p = exp2_checked(log_family_wrapper_discarded_unless_normal!(b, log_2_normal) * 0.4);
    let high = b * b * p;
    if c <= 0.04045 {
        low
    } else {
        high
    }
}

/// Converts a linear color component in `[0, 1]` to sRGB.
#[inline(always)]
#[allow(clippy::neg_cmp_op_on_partial_ord)]
pub fn linear_to_srgb(l: f32) -> f32 {
    let low = l * 12.92;
    let p =
        exp2_checked(log_family_wrapper_discarded_unless_normal!(l, log_2_normal) * (1.0 / 2.4));
    let high = fma(1.055, p, -0.055);
    if l <= 0.0031308 {
        low
    } else {
        high
    }
}

/// Straight port of jodiemath's remainderf: x - round(x/y)*y (ties away from
/// zero, via f32::round -- not IEEE 754 remainder's ties-to-even).
/// round(x/y)*y's absolute error scales with ulp(x), which swamps the true
/// remainder (at most |y|/2) once |x/y| is large -- inherited from the C
/// original's identical formula, only reliable while |x/y| stays moderate.
macro_rules! remainder_style_combine {
    ($x:expr, $y:expr, $q:expr) => {{
        let normal = fma(-$q, $y, $x);
        // IEEE754 defines the sign of an exact-zero remainder/fmod result to
        // match x's sign (also fmod's own C99 "same sign as x" contract) -- but
        // when x is a nonzero exact multiple of y, `-q*y` exactly cancels x,
        // and IEEE754 subtraction of two equal-magnitude values always gives
        // *positive* zero regardless of the "intended" sign, silently dropping
        // it (same class of bug as sinf_poly's own -0.0 note, one level further
        // out: here the cancellation is real, not a q that happens to be zero).
        // Confirmed against libm (Python's math.remainder/math.fmod, a real C
        // library, not a hand-derived assumption) for every exact-multiple case
        // checked.
        let normal = if normal == 0.0 {
            normal.copysign($x)
        } else {
            normal
        };
        let r = if $x == 0.0 && !normal.is_nan() {
            $x
        } else {
            normal
        };
        if $y.is_infinite() && $x.is_finite() {
            $x
        } else {
            r
        }
    }};
}

#[doc(alias = "remainderf")]
#[inline(always)]
pub fn remainder(x: f32, y: f32) -> f32 {
    let q = (x / y).round();
    remainder_style_combine!(x, y, q)
}

/// IEEE 754 floating-point remainder (ties to even).
#[inline(always)]
pub fn remainder_ieee(x: f32, y: f32) -> f32 {
    let q = (x / y).round_ties_even();
    remainder_style_combine!(x, y, q)
}

/// `remainder` without domain checks: valid for `x != 0` and finite `y`.
#[inline(always)]
pub fn remainder_unchecked(x: f32, y: f32) -> f32 {
    let q = (x / y).round();
    fma(-q, y, x)
}

/// Self-correcting `remainder` for `|x/y| <= 2^24`.
#[inline(always)]
pub fn remainder_checked(x: f32, y: f32) -> f32 {
    let q0 = (x / y).round();
    let r0 = fma(-q0, y, x);
    // The correction always moves `r0` toward zero, so it is `r0 -
    // copysign(|y|, r0)` -- no `+-1` multiplier to select and no fma. (The two
    // forms differ only at `r0 == +-0.0`, where the guard below keeps `r0`
    // anyway.) `ay` is shared with that guard.
    let ay = y.abs();
    let r1 = r0 - ay.copysign(r0);
    let normal = if r0.abs() > ay * 0.5 { r1 } else { r0 };
    // Same exact-cancellation sign bug as remainder_style_combine! (see its own
    // comment): a nonzero x that's an exact multiple of y exactly cancels to
    // +0.0 regardless of x's sign, silently dropping it.
    let normal = if normal == 0.0 {
        normal.copysign(x)
    } else {
        normal
    };
    let r = if x == 0.0 && !normal.is_nan() {
        x
    } else {
        normal
    };
    if y.is_infinite() && x.is_finite() {
        x
    } else {
        r
    }
}

/// `remainder` using f64 for `|x/y| <= 2^53`.
#[inline(always)]
pub fn remainder_wide(x: f32, y: f32) -> f32 {
    let xd = x as f64;
    let yd = y as f64;
    // Ties away from zero, matching `remainder`/`remainder_checked`'s
    // documented divergence from IEEE754 (`remainder_ieee` is the ties-even
    // variant). At an exact half-integer `x/y` this lands `r0` on exactly
    // `+-|y|/2`, which the strict `>` below leaves alone.
    let q = (xd / yd).round();
    // Exact, and that is the whole point: `x` is a multiple of `ulp(x)` and
    // `q*y` a multiple of `ulp(y)`, so `x - q*y` is a multiple of `min(ulp(x),
    // ulp(y))` with magnitude `<= 1.5|y|` -- 24 significant bits, representable
    // in an *f32*, never mind an f64. The fma forms `q*y` to full width
    // internally, so no error-free transform is needed at this step at all.
    let r0 = f64::mul_add(-q, yd, xd);
    // The correction always moves `r0` toward zero, so it is `r0 -
    // copysign(|y|, r0)` -- no `+-1` multiplier to select and no fma. (The two
    // forms differ only at `r0 == +-0.0`, where the guard below keeps `r0`
    // anyway.) `ay` is shared with that guard.
    let ay = yd.abs();
    let r1 = r0 - ay.copysign(r0);
    let normal = if r0.abs() > ay * 0.5 { r1 } else { r0 };
    // Exact by the same argument as `r0`, so the narrowing rounds
    // nothing: a true IEEE remainder is always representable in its
    // operands' own format.
    let normal = normal as f32;
    // Same exact-cancellation sign bug as remainder_style_combine! (see its own
    // comment): a nonzero x that's an exact multiple of y exactly cancels to
    // +0.0 regardless of x's sign, silently dropping it.
    let normal = if normal == 0.0 {
        normal.copysign(x)
    } else {
        normal
    };
    let r = if x == 0.0 && !normal.is_nan() {
        x
    } else {
        normal
    };
    if y.is_infinite() && x.is_finite() {
        x
    } else {
        r
    }
}

/// Truncated floating-point remainder `x - trunc(x/y) * y`.
#[doc(alias = "fmodf")]
#[inline(always)]
pub fn fmod(x: f32, y: f32) -> f32 {
    let q = (x / y).trunc();
    remainder_style_combine!(x, y, q)
}

/// Self-correcting `fmod` for `|x/y| <= 2^24`.
#[inline(always)]
pub fn fmod_checked(x: f32, y: f32) -> f32 {
    let q0 = (x / y).trunc();
    let r0 = fma(-q0, y, x);
    let wrong_sign = r0 != 0.0 && (r0 > 0.0) != (x > 0.0);
    // The correction's sign depends on *both* which failure mode fired
    // (overshoot/`wrong_sign` needs q0 nudged toward zero, undershoot needs it
    // nudged away) *and* whether x/y have the same sign -- verified against all
    // 4 sign combinations crossed with both failure modes by hand (8 cases)
    // before trusting this, not derived by inspection alone (an earlier
    // y-sign-only version passed 4 of those 8 and failed the rest).
    let same_sign = (x > 0.0) == (y > 0.0);
    let adj = if wrong_sign != same_sign { 1.0 } else { -1.0 };
    let r1 = fma(-adj, y, r0);
    let needs_fix = wrong_sign || r0.abs() >= y.abs();
    let normal = if needs_fix { r1 } else { r0 };
    // Same exact-cancellation sign bug as remainder_style_combine! (see its own
    // comment): a nonzero x that's an exact multiple of y exactly cancels to
    // +0.0 regardless of x's sign, silently dropping it.
    let normal = if normal == 0.0 {
        normal.copysign(x)
    } else {
        normal
    };
    let r = if x == 0.0 && !normal.is_nan() {
        x
    } else {
        normal
    };
    if y.is_infinite() && x.is_finite() {
        x
    } else {
        r
    }
}

/// `fmod` without domain checks.
#[inline(always)]
pub fn fmod_unchecked(x: f32, y: f32) -> f32 {
    let q = (x / y).trunc();
    fma(-q, y, x)
}

/// Computes the Euclidean remainder `x - y * floor(x / y)` (`0 <= result < |y|`).
#[inline(always)]
pub fn rem_euclid(x: f32, y: f32) -> f32 {
    let r = fmod(x, y);
    if r < 0.0 {
        r + y.abs()
    } else {
        r
    }
}

/// Computes the Euclidean quotient `floor(x / y)`.
#[inline(always)]
pub fn div_euclid(x: f32, y: f32) -> f32 {
    let q = (x / y).trunc();
    let r = fmod(x, y);
    if r < 0.0 {
        if y > 0.0 {
            q - 1.0
        } else {
            q + 1.0
        }
    } else {
        q
    }
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
