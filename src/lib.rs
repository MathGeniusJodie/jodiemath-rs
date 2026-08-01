// godbolt flags -C opt-level=3 -C target_feature=+fma

// This crate's entire accuracy and performance story depends on `fma`
// (see the function just below) compiling to a single hardware
// instruction. Without FMA, `f32::mul_add` falls back to a ~2x-slower
// libm call that also rounds *differently* (two roundings instead of
// one) -- every ulp figure in this crate's doc comments, readme.md, and
// examples/accuracy.rs assumes the single-rounding hardware form.
// `.cargo/config.toml` sets `target-cpu=native` for exactly this reason,
// but an environment `RUSTFLAGS` silently overrides (not merges with)
// that setting, which would disable FMA with no build error and no
// runtime symptom beyond quietly-wrong accuracy numbers. Fail loudly at
// compile time instead: `target_feature = "fma"` is set by the compiler
// whenever FMA is actually enabled, regardless of how that happened.
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

// `2^k` as a bare exponent field, for an integer-valued `k` in
// `[-127, 128]` -- the range each caller's own domain contract already
// guarantees. `k = 128` deliberately yields `+inf` (field 255, mantissa
// 0), which is how the unchecked exp2 family overflows, and `k = -127`
// yields `+0.0`.
//
// The magic constant is `1.5*2^23 + 127`, so `k + MAGIC` lands in the
// ulp-1 binade `[2^23, 2^24)` with `k + 127` sitting in the low 9 bits
// -- a single `<< 23` then moves those into the exponent field with a
// zero sign bit and a zero mantissa, needing no mask. That is the whole
// saving over the older `(k + 383.0) << 8 & EXPONENT_MASK` form: `+383`
// lands in the `2^8` binade, whose exponent field (135) is odd, so the
// shift carries a 1 into the sign bit that a mask then has to clear.
// Folding the `+127` into the magic constant rather than adding it to
// the bits is what keeps this at two ops.
//
// Macro, not a fn -- see exp_r_poly!.
const EXP2INT_MAGIC: f32 = 12583039.0; // 1.5 * 2^23 + 127
macro_rules! exp2int_field {
    ($k:expr) => {
        f32::from_bits(($k + EXP2INT_MAGIC).to_bits() << 23)
    };
}

/// Round to the nearest integer (ties-to-even), for `|x| <= 2^22`
/// (backlog idea #185): the magic-constant idiom this crate uses
/// throughout its own reductions (`exp`'s `k`, `sinpi`'s `q`'s cheaper
/// sibling, etc.), exposed standalone for callers building their own
/// periodic reductions. Adding `1.5*2^23` forces the sum to round to an
/// integer at that magnitude (ties resolved by the hardware's default
/// round-to-nearest-even, unlike `f32::round`'s ties-away-from-zero,
/// which needs extra emulation instructions), and subtracting the same
/// constant back off reveals that integer as a plain `f32` -- one `fma`
/// and one subtract, versus `.round()`'s multi-instruction ties-away
/// path. Exact only up to `|x| <= 2^22`: past that the *combined*
/// magnitude `x + 1.5*2^23` needs more than 23 mantissa bits to keep
/// distinguishing integers one apart, so the addition itself starts
/// rounding to a coarser grid -- verified directly, not assumed: at
/// `x=5_000_001.0` (already an exact integer, no rounding even needed)
/// this returns `5_000_000.0`, off by one, not `x` unchanged. No
/// guarantee of any kind past the documented bound, same convention as
/// this crate's other `_unchecked`-style contracts.
///
/// The trailing `.copysign(x)` fixes a real sign-of-zero gap the bare
/// idiom has on its own: for any negative `x` that rounds to zero
/// (`x` in `[-0.5, 0)`), `x + 1.5*2^23` rounds to *exactly*
/// `1.5*2^23` -- subtracting the same constant back off is then
/// `1.5*2^23 - 1.5*2^23`, which IEEE754 always resolves to `+0.0`
/// regardless of `x`'s own sign, not the `-0.0` correctly-rounded
/// output needs. Caught by a real exhaustive sweep before shipping
/// (`round_ties_even`-vs-bare-idiom mismatches: 0 in magnitude, but
/// ~1.057 billion in sign, every one confined to `|x| <= 0.5`) --
/// this crate's ~16 *internal* call sites never needed this fix
/// (their own reductions only ever consume the rounded integer's
/// *value*, e.g. as an exponent, where `0` and `-0` are
/// interchangeable), which is exactly why the gap went unnoticed
/// until this function's contract had to stand on its own for a
/// general-purpose caller. `copysign` is a no-op for any nonzero
/// result (reapplying the sign the subtraction already got right), so
/// this only ever changes the zero case.
#[inline(always)]
pub fn fast_round_int(x: f32) -> f32 {
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    (fma(x, 1.0, ROUND_MAGIC) - ROUND_MAGIC).copysign(x)
}

// Shared denormal handling for log_2/ln/log10/log2_df: scale a denormal
// input up by 2^24 before the single normal-path evaluation, tracking
// the compensating exponent offset to fold back in afterwards. Macro,
// not a fn -- see exp_r_poly! for why that distinction matters here.
macro_rules! denormal_rescale {
    ($x:expr) => {{
        let tiny = $x < f32::MIN_POSITIVE;
        let xs = if tiny { $x * 16777216.0 } else { $x };
        let koff = if tiny { -24.0 } else { 0.0 };
        (xs, koff)
    }};
}

// exp(r) for tiny r, shared by exp/exp_checked/expm1/exp_m1_over_x/tanh/
// sigmoid (c0/c1 pinned to exactly 1.0 -- see exp's own body comment).
// Deliberately a macro, not a fn: a fn-based sharing attempt caused a
// real, reproduced +32% mca regression on an unrelated caller purely
// from the new function-call boundary's scheduling side effects (see
// graveyard.md §hyperbolics). A macro is pure textual substitution with no
// boundary at all -- verified equivalent via a full pre/post assembly
// diff. The same reasoning applies to every other shared-body macro in
// this file. Each caller keeps its own reduction (`k`/`r`) and exponent
// reconstruction/final combine around this, since those differ.
//
// The top coefficient pair folds in at the `r^2` level rather than its
// own `r^4` one (`ln_normal`'s trick, and `exp2_q_poly!`'s): `r^4` is
// never formed, so this is one plain multiply cheaper than the balanced
// 4-way split at the same fma critical-path depth, and one rounding
// fewer on the shared power.
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
// `exp_checked` is, minus the input clamp -- the caller supplies an
// argument already inside `EXP_CLAMP_LO..=EXP_CLAMP_HI`. Split out for
// callers that can *prove* one side of that clamp is unreachable and so
// should not pay for it: `erfc`'s argument is `-(xs*xs)`, a negated
// square, so the upper bound is dead by construction (see its own
// comment). Macro, not a fn -- see exp_r_poly! for why a new shared-fn
// boundary is the thing to avoid here; verified via a full pre/post
// assembly diff that `exp_checked` and its other callers are unchanged.
macro_rules! exp_reduce {
    ($x:expr) => {{
        let x = $x;
        const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
        let k = fma(x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
        let r = fma(-k, LN2_HI, x);
        let r = fma(-k, LN2_LO, r);
        let p = exp_r_poly!(r);
        let (t1, t2) = exp2_field_split(k);
        p * t1 * t2
    }};
}

/// `exp_checked`'s clamp bounds: `exp2_checked`'s own `k` boundary
/// (`[-151, 128)`) converted into `x`'s units.
const EXP_CLAMP_LO: f32 = -104.66522426455174;
const EXP_CLAMP_HI: f32 = 88.72283911167308;

// Q(f) = (2^f - 1)/f, shared by exp2/exp2_checked/exp10/exp10_checked/
// exp2m1/exp2_checked_df. Macro, not a fn -- see exp_r_poly!. Returns
// `q`; each caller does its own final combine (exp2/exp10's single-field
// `fma(q, exp2int*f, exp2int)` vs. the others' k1/k2-split
// `p = fma(q, t1*f, t1); p*t2`).
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

// `exp2_q_poly!`'s centered sibling: the same `Q(f) = (2^f - 1)/f` and
// the same 3-balanced-pair Estrin shape, refit for `f in [-0.5, 0.5]`
// instead of `[0, 1)`. Only `exp10_checked` uses it, to keep `round`'s
// own residual rather than paying a floor-adjust to convert it. This is
// a genuine per-caller domain difference, not the kind of shared-poly
// decoupling that got rejected elsewhere in IDEAS.md (there the callers'
// domains were identical): centering is an affine change of variable, so
// the same degree buys the same accuracy -- the idealized max relative
// error is 1.073e-8 here against the `[0,1)` fit's 1.217e-8, i.e. 0.179
// vs 0.203 ulp-equivalent.
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

// Shared by exp_pos_neg/exp_pos_neg_checked_half (sinh/cosh's
// unchecked/checked exp(x)/exp(-x) core): one Cody-Waite reduction,
// an even/odd-split poly, and the t1n/t2n reciprocal construction.
// Only each caller's optional input clamp differs, so that stays at the
// call site. Returns `(p_pos, p_neg, t1, t2, t1n, t2n)` where p_pos/p_neg
// are `0.5*e^(+-r)`, not `e^(+-r)`: every caller is a sinh/cosh that wants
// the halves, and folding the 0.5 into the *poly constants* is free where
// a `t1 * 0.5` multiply is not. Exact, not an approximation -- halving a
// float is exact, and scaling every operand of an fma by the same power of
// two scales its correctly-rounded result by exactly that power of two, so
// `p_pos`/`p_neg` are bit-for-bit half of what the unhalved poly gives.
// Halving here rather than at the end also keeps the checked caller's
// premature-overflow fix (see exp_pos_neg_checked_half's own point 2):
// the product never has to represent the unhalved `e^x` at any stage.
// Macro, not a fn -- see exp_r_poly!.
macro_rules! exp_pos_neg_core {
    ($x:expr) => {{
        const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
        let k = fma($x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
        let r = fma(-k, LN2_HI, $x);
        let r = fma(-k, LN2_LO, r);
        // The fitted coefficients, each exactly halved (see this macro's
        // own comment). Written as `* 0.5` rather than as pre-divided
        // decimal literals so the fit's own values stay legible and the
        // halving cannot be mistranscribed.
        let c: [f32; 4] = [
            4.99993e-1 * 0.5,
            1.6667245e-1 * 0.5,
            4.188372e-2 * 0.5,
            8.300987e-3 * 0.5,
        ];
        // Each half is Horner in `r^2`, not a leading `c*r^4` term: `r^4`
        // is never formed, one plain multiply cheaper across both halves
        // at the same fma critical-path depth. Same fold as exp_r_poly!.
        let r2 = r * r;
        let e = fma(fma(c[2], r2, c[0]), r2, 0.5);
        let o = fma(fma(c[3], r2, c[1]), r2, 0.5);
        let p_pos = fma(r, o, e);
        let p_neg = fma(-r, o, e);
        let (t1, t2) = exp2_field_split(k);
        // t1n = 1/t1, t2n = 1/t2: both are exact power-of-two fields, and
        // for a power-of-two float with bit pattern b = (127+e)<<23 the
        // reciprocal 2^-e has bit pattern (127-e)<<23 = 0x7F000000 - b.
        // Exactly what exp2_field_split(-k) would produce (round-half-to-
        // even is antisymmetric under negation), but built with two
        // integer subtracts instead of a second magic-round chain.
        let t1n = f32::from_bits(0x7F00_0000u32.wrapping_sub(t1.to_bits()));
        let t2n = f32::from_bits(0x7F00_0000u32.wrapping_sub(t2.to_bits()));
        (p_pos, p_neg, t1, t2, t1n, t2n)
    }};
}

// Pade approximant for e^v-1 near v=0 (an exact closed-form identity,
// not an empirical fit), shared by expm1/exp_m1_over_x/exp2m1/tanh.
// Two arms, not one: a macro invocation parses as one atomic expression,
// so writing `v * pade_expm1_ratio!(v)` against a ratio-only macro would
// silently reassociate `(v*N)/D` into `v*(N/D)` -- mathematically equal
// but not bit-identical. The `mul` arm keeps `v * fma(N)` and the final
// `/ fma(D)` inside one expansion so Rust precedence reproduces the
// original left-to-right grouping exactly.
macro_rules! pade_expm1_ratio {
    ($v:expr) => {
        fma(-1.9999927, $v * $v, -120.0) / fma($v, fma($v, $v - 12.000030, 59.999996), -120.0)
    };
    ($v:expr, mul) => {
        $v * fma(-1.9999927, $v * $v, -120.0) / fma($v, fma($v, $v - 12.000030, 59.999996), -120.0)
    };
}

// Shared by log_2/ln/log10: the denormal-rescale + special-case-select
// wrapper around each caller's own `_normal` fn. Edge handling uses
// selects (no early returns) so array loops auto-vectorize. `spec` is
// `-inf` for `+-0`, `NaN` for `x < 0` (includes `-inf`); its select
// input only depends on `x`, so it resolves in parallel with the poly.
// `!(x < f32::INFINITY)` exploits NaN's always-false comparisons to
// catch both +inf and NaN in one cheap fcmp (`x*x` is `inf`/`nan` there,
// respectively; false for `-inf`) -- hence each caller's
// `#[allow(clippy::neg_cmp_op_on_partial_ord)]`. Macro, not a fn -- see
// exp_r_poly!.
macro_rules! log_family_edges {
    ($x:expr, $r:expr) => {{
        let r = $r;
        let spec = if $x == 0.0 { f32::NEG_INFINITY } else { f32::NAN };
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

// `log_family_wrapper!` for callers whose argument provably can't be a
// positive denormal, so the rescale's compare/multiply/two selects are
// dead code: `u = 1.0 + x` is such an argument for *every* f32 `x`, since
// `1+x` is exact by Sterbenz once `x <= -0.5` (making the smallest
// positive `u` exactly `2^-24`, ~10^30 above `f32::MIN_POSITIVE`) and is
// `>= 0.5` otherwise -- verified exhaustively over all 2^32 patterns, not
// argued from the bound alone. The zero/negative/inf/NaN arms are all
// still reachable at those call sites (`u == 0` at `x == -1`, `u < 0`
// below it) and are shared verbatim. Same "the guard is the licence"
// lever as `asinh`/`acosh`/`atanh`, one level down: here the licence
// comes from the argument's own construction rather than a branch guard.
macro_rules! log_family_wrapper_no_denormal {
    ($x:expr, $normal:ident) => {
        log_family_edges!($x, $normal($x, 0.0))
    };
}

// `log_family_wrapper!` for callers whose own outer select *discards* this
// value entirely unless the argument is positive and normal, so only
// `log_family_edges!`'s inf/NaN arm survives -- the zero/negative selects
// and the denormal rescale all compute results nothing can observe. Note
// what the licence is and isn't: every arm here is still *evaluated* for
// every input (branchless), so this is sound only when the caller's select
// provably drops the result for zero/negative/denormal arguments. `+inf`
// and NaN are the exception that keeps one arm alive -- they flow through
// the outer select, so `x*x` still has to turn them into `+inf`/NaN.
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
// log_2_normal's own poly, and only its own: the double-float `log2_df`
// needs ~31 relative bits, which this `s = m - 1` shape cannot reach at
// any degree (its f32 evaluation rounding lands straight on the result
// for x near 1), so that one uses an atanh-form reduction instead.
//
// This is the *peeled* poly `Q(s) = (log2(1+s)/s - log2(e)) / s`, degree 8,
// not the degree-9 `P(s) = log2(1+s)/s` the combine used to evaluate whole.
// Fitted by an ulp-weighted LP against `s^2/log2(1+s)` -- the weight that
// makes the fit minimise the *result*'s relative error, since this poly
// only ever reaches the answer scaled by `s^2` -- and not by dropping a
// term off `P`. One degree lower than `P` because `Q` reaches the answer
// diluted, so the fit no longer has to carry the whole budget; but not two
// (degree 7's fit is 2.6e-8, ~0.44 ulp, which lands back on max 2).
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
    // log2(m) = log2(e)*s + s^2*Q(s), with the leading term kept *out* of
    // the polynomial rather than evaluated as its constant coefficient.
    //
    // That is where this function's accuracy comes from. Written as the old
    // `k + s*P(s)`, the answer for x near 1 (k == 0) simply *is* `s*P(s)`,
    // so all three of P's own full-weight evaluation roundings -- the
    // leading coefficient pair, and both Estrin joins above it -- landed
    // directly on the result at ~2^-24 each, and that, not the fit, was the
    // 3-ulp max. Peeling `log2(e)*s` out demotes every one of them: what is
    // left is scaled by `s^2/log2(1+s) <= 0.21` before it reaches the
    // answer, so the whole polynomial's rounding is worth about a fifth of
    // an ulp and only the final fma still rounds at full weight.
    // Exhaustively: max 3 -> 1 ulp (matching std's own max), avg 0.0031 ->
    // 0.0028, one instruction and one uOp *fewer*.
    //
    // The `s^2` factor rides into the poly's own low group (`a = s2 * l0`)
    // instead of multiplying the finished `Q`. Same op count, but it keeps
    // the whole thing three levels deep the way the un-peeled degree-9 form
    // was: multiplying afterwards puts `s2*q` on the end of the chain and
    // costs a second level, which measures as a further +11% latency.
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
    // `k` joins *last*, in its own single rounding. Threading it through
    // the peeled combine instead (`fma(s, LOG2_E, fma(s2, q, k))`) rounds
    // twice at `k`'s scale, and for |k| >= 1 that is the whole error
    // budget: it still reaches max 2 but costs 40x on the average
    // (0.0031 -> 0.1249, measured exhaustively). Built this way, log2(m)
    // is finished to ~2^-25 absolute *before* `k` is anywhere near it, so
    // every octave but the k == 0 one is back to a single rounding.
    // Likewise `fma(s2, q, s * LOG2_E)`, which looks cheaper because
    // `s*LOG2_E` starts early: it rounds that leading term on its own at
    // full weight, which is exactly what the peel exists to avoid, and
    // measures max 2.
    let lm = fma(s, std::f32::consts::LOG2_E, sq);
    lm + k
}

/// log_2 without domain checks: valid for positive normal finite x only
/// (no handling for zero, negative, denormal, inf, or nan -- those give a
/// plausible-looking but wrong finite value instead of NaN/-inf). Mirrors
/// exp2/exp2_checked's fast/full-safety split; drops the denormal-rescale
/// multiply and both post-hoc selects log_2 pays on every call.
#[inline(always)]
pub fn log_2_unchecked(x: f32) -> f32 {
    log_2_normal(x, 0.0)
}

/// exp2 without domain checks: valid for x in [-126, 128), i.e. normal
/// (non-denormal, finite, nonzero) results only. Outside that range the
/// exponent construction wraps around and the result is garbage (including
/// for nan) -- and note the garbage is not even sign-correct: the wrapped
/// exponent field can land in the sign bit, so an out-of-domain call can
/// return a *negative* value from a function that is mathematically
/// positive everywhere. That is deliberate (it makes a domain violation
/// unmistakable rather than plausible-looking) but it means `exp2(x) >= 0`
/// is not an invariant you can lean on here. Use exp2_checked for
/// full-range handling; this version is ~2.7 ns faster in serial latency.
#[doc(alias = "exp2f")]
#[inline(always)]
#[allow(clippy::approx_constant)] // g0's constant term is a fitted minimax
// coefficient near ln(2), not ln(2) itself (bit pattern deliberately differs)
pub fn exp2(x: f32) -> f32 {
    // exp2(floor(x)) * exp2(fract(x)) == exp2(x). exp2int must come from
    // the same floor(x) as f: computing it from x + 383 double-counts the
    // integer part when x + 383 rounds up across an integer (e.g.
    // x = 4.9999999). A k=round(x) reduction (tighter Q(f) fit) was tried
    // and rejected here and on every other exp2_q_poly! caller, each for
    // its own reason -- here, round can land k=128 inside the promised
    // [-126,128) domain, which this single-field construction can't
    // represent (NaN for a legit input). See graveyard.md §exp/exp2.
    let k = x.floor();
    let f = x - k;
    let exp2int = exp2int_field!(k);
    // Q(f) = (2^f - 1)/f, degree 5, grouped into 3 balanced pairs (g0, g1,
    // g2) instead of two degree-2 Horner halves: same 6 coefficients and
    // the same 4-deep fma critical path, but the combine only ever needs
    // f^2 (never exp2int*f^4), so 2 fewer plain multiplies per call.
    // Avg ulp 0.069, max 1 (dense sweep of the whole unchecked domain).
    let q = exp2_q_poly!(f);
    fma(q, exp2int * f, exp2int)
}

/// `2^(k+f)` for an already-integer-valued `k` and `f` in `[0,1)`
/// (backlog idea #116): [`exp2`]'s own combine step, exposed directly
/// for user custom-base kernels (and this crate's own `exp10`, though
/// not refactored to call through here -- its existing, separately
/// verified body is untouched) that already have their own `k`/`f` and
/// want to skip re-deriving this exact exponent-field-plus-poly
/// combine. No domain check: same `[-126,128)` contract as `exp2`
/// itself (garbage out for `k` outside that range), and `f` outside
/// `[0,1)` is simply a different (still well-defined) input to the same
/// polynomial, not a checked contract.
#[inline(always)]
#[allow(clippy::approx_constant)]
pub fn exp2_kf(k: f32, f: f32) -> f32 {
    let exp2int = exp2int_field!(k);
    let q = exp2_q_poly!(f);
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
    // an integer (e.g. x = 4.9999999). k=round(x) was tried and rejected
    // here too (perf and max-ulp regression, no structural savings) --
    // see graveyard.md §exp/exp2.
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

/// x*2^n (backlog idea #86), C's `ldexp`/`scalbn`. A first version tried
/// `x * exp2_checked(n as f32)`, reusing `exp2_checked` wholesale --
/// killed by a real fuzz run, not caught on paper: `exp2_checked`'s own
/// clamp (`.clamp(-151.0, 128.0)`) is only safe *in isolation*, where the
/// clamped exponent alone determines overflow. Here `x`'s own magnitude
/// can compensate for an `n` past that clamp (e.g. `x=7.27e-9, n=148`:
/// the true product `~2.59e36` is finite, but clamping `n` to `128`
/// first computes `x * 2^128`, which overflows to `inf` on its own even
/// though the real answer wouldn't) -- the same "premature overflow
/// before a compensating factor gets a chance to apply" class of bug
/// `exp`'s own t1-weave revisit (idea #24) found. Fixed by decomposing
/// `x` via [`frexp`] first (mantissa always `< 1`, matching the same
/// safety property `exp2_checked`'s own `p*t1*t2` combine relies on for
/// its own poly correction) and combining the *exponents* before any
/// clamping, so `x`'s own headroom is available to the combined value,
/// not just to `n` alone.
///
/// A second real bug, also only found by fuzzing: clamping a combined
/// exponent that's *grossly* out of range (not just past the boundary by
/// a little) down to `128`/`-151` and reconstructing anyway silently
/// computes `mantissa * 2^128` (a large but finite value) instead of the
/// true answer, which at that magnitude gap (checked: 68 past the
/// boundary in one real failing case) overflows regardless of `mantissa`
/// -- no mantissa in `[0.5,1)` can pull a target exponent of `196` back
/// under `f32::MAX`. Fixed by deciding overflow/underflow from the
/// *unclamped* combined exponent directly (`> 128` always overflows,
/// `< -151` always underflows to `0` -- verified against
/// `exp2_checked`'s own already-confirmed bit-exact boundary at exactly
/// those values, since `mantissa < 1` can only shrink the result
/// relative to that bare-power-of-two reference, never grow it) rather
/// than inferring it from whatever the clamped reconstruction happens to
/// produce. Only the genuinely in-range combined exponents reach the
/// `exp2_field_split` reconstruction below, where interleaving
/// `mantissa` between `t1` and `t2` (mirroring `exp2_checked`'s own
/// combine) is safe -- confirmed against a real exhaustive/fuzz sweep,
/// including the exact boundary (`ldexp(f32::MAX, 0)` round-trips
/// bit-exactly).
///
/// Not wired into `examples/mca.rs`/`mca_target.rs`: this function's own
/// multi-exit-path branching (the same documented harness limitation
/// `frexp`/`rootn` have) corrupts llvm-mca's inline-asm region markers
/// for the *whole* assembly file, not just this region -- confirmed by
/// removing it restores every other function's mca numbers. Use
/// quickbench for this one.
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
    if x == 0.0 || !x.is_finite() { x } else { saturated }
}

/// Decompose `x` into `(mantissa, exponent)` with `mantissa` in
/// `[0.5, 1)` such that `x == mantissa * 2^exponent` (backlog idea #86),
/// C's `frexp`. Denormal `x` is rescaled by `2^24` first (same
/// `denormal_rescale!` idiom `log_2`/`ln`/`log10` already use), tracked
/// via `koff` and folded back into the reported exponent at the end.
/// The mantissa itself is built by a fixed bit substitution, not a
/// relative adjustment of `x`'s own exponent field: replacing the
/// biased exponent with the constant `126` (representing `2^-1`)
/// directly encodes `1.mantissa_bits * 2^-1`, which is always in
/// `[0.5, 1)` regardless of `x`'s original exponent -- a *relative*
/// `-1` adjustment of the original field instead would wrap into a
/// denormal encoding right at the smallest normal exponent (where the
/// implicit leading-1 assumption breaks), a boundary bug avoided
/// entirely by never depending on the original field's value for
/// anything but the *reported* exponent. `x == 0`/non-finite `x` are
/// overridden to `(x, 0)`: zero has no `[0.5,1)` decomposition at all,
/// and `+-inf`/`NaN`'s own bit patterns would otherwise feed the same
/// fixed-substitution trick and produce a finite, wrong mantissa.
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
    (if is_special { x } else { mantissa }, if is_special { 0 } else { exponent })
}

/// 10^x. Naively rounding `x*LOG2_10` once before `exp2_checked` even
/// starts loses precision that grows with `|x|` (the same flaw `exp`'s
/// own doc comment describes for `exp2(x*LOG2_E)`). Fixed the same way
/// as `exp`: `k = round(x*LOG2_10)` (only needs to land on the right
/// *integer*, a coarse multiply is fine for that), then reduce `x`
/// itself (not `x*LOG2_10`) via a Cody-Waite split of `LOG10_2` (already
/// defined for `log10`'s own fix): `d = x - k*LOG10_2_HI - k*LOG10_2_LO`
/// stays small and precisely known in `x`'s own units (mirrors `exp`'s
/// `r = x - k*LN2_HI - k*LN2_LO`).
///
/// Unlike `exp`, this can't just hand `k + d*LOG2_10` to `exp2_checked`
/// as a single combined argument -- that recombination itself
/// reintroduces the exact bug being fixed: adding the *small* correction
/// `d*LOG2_10` to the *large* integer `k` (up to ~127) forces the sum to
/// round to `k`'s own coarse ulp (e.g. ulp(75) ~ 9e-6), silently
/// discarding the precision the careful reduction just earned (measured:
/// max ulp 45, traced to exactly this recombination step). Fixed by never
/// forming that combined value at all: `exp2_checked`'s own internal
/// split ("any split k=k1+k2 ... works") is reproduced here directly
/// against this function's precisely-known integer `k` and fractional `f`
/// (floor-adjusted from the round-based reduction into `exp2_checked`'s
/// `[0,1)` convention).
///
/// Dropping the floor-adjust (feeding `kr`/`fr` straight through with a
/// round-domain Q(f) fit) was tried and rejected: at the overflow-
/// saturation boundary the round convention allows `f < 0`, so `2^f < 1`
/// can pull `t1*t2 = 2^128`'s product back *under* `f32::MAX`, giving
/// `exp10_checked(inf)` a finite result -- a contract violation only
/// edgecheck.rs's special-value pins caught. See graveyard.md §exp/exp2.
///
/// The leading `x` clamp bound is `[-45.154503, 38.53184]` (backlog idea
/// #113), not an arbitrarily-wide safety net: these are the exact bit-
/// level boundary floats (found by stepping one ulp at a time through
/// the real reduction above) where `k` first reaches `-151`/`128`
/// respectively -- clamping `x` *to* one of these literals forces `k` to
/// land exactly on that same boundary value by construction, making the
/// separate trailing `k.clamp(-151.0, 128.0)` provably redundant over
/// the *entire* domain (every in-range `x` was already producing `k` in
/// `[-151,128]` on its own; every out-of-range `x` now clamps straight
/// to a literal that reproduces the exact `k` the old wider clamp plus
/// trailing `k`-clamp used to saturate to). Verified bit-identical
/// across the full exhaustive sweep before adopting -- this is a pure
/// codegen win, not a behavior change.
// Shared by exp10/exp10_checked: the round-based reduction (`kb`/`kr`/`d`/
// `fr`/floor-adjust to `(k, f)`) is identical between the two -- only
// exp10_checked's own leading `x` clamp and trailing `k` clamp (needed
// since it feeds the k1/k2-split combine, unlike exp10's single-field
// one) differ, both left at the call site. Macro, not a fn -- see
// exp_r_poly!.
macro_rules! exp10_reduction {
    ($x:expr) => {{
        const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
        let kb = fma($x, std::f32::consts::LOG2_10, ROUND_MAGIC);
        let kr = kb - ROUND_MAGIC; // round(x*log2(10)), coarse multiply is fine
        let d = fma(-kr, LOG10_2_HI, $x);
        let d = fma(-kr, LOG10_2_LO, d);
        let fr = d * std::f32::consts::LOG2_10; // small, precise correction in log2 units, in [-0.5, 0.5]
        // floor-adjust (kr, fr) from round's [-0.5,0.5] convention to
        // exp2_checked's own floor-based [0,1) convention -- both ops exact
        // or near-exact since they only ever combine values of comparable
        // magnitude (unlike the rejected single-combine above).
        let adjust = if fr < 0.0 { 1.0 } else { 0.0 };
        let k = kr - adjust;
        let f = fr + adjust;
        (k, f)
    }};
}

#[inline(always)]
#[allow(clippy::approx_constant)] // g0's constant term is a fitted minimax
// coefficient near ln(2), not ln(2) itself (bit pattern deliberately differs)
pub fn exp10_checked(x: f32) -> f32 {
    // Clamped before the reduction starts (matching exp2_checked's own
    // early-clamp pattern) so +-inf can't poison `d = x - kr*LOG10_2`
    // with an inf-inf NaN -- NaN itself passes through unaffected
    // (f32::clamp preserves NaN in the receiver). See this function's own
    // doc comment for why these exact bounds make the old separate `k`
    // clamp redundant (removed).
    let x = x.clamp(-45.154503, 38.53184);
    // Unlike `exp10`, this keeps `round`'s own centered `f in [-0.5, 0.5]`
    // instead of paying `exp10_reduction!`'s floor-adjust (a compare, a
    // select and two add/subs) to reach `exp2_q_poly!`'s `[0,1)`
    // convention. That needs its own Q refit -- see `exp2_q_poly_centered!`
    // -- and is only safe here, not in `exp10`: `round` can put `k` at
    // 128, which the k1/k2 split represents fine but a single exponent
    // field cannot (the reason `exp10`'s own doc comment gives for
    // keeping the adjust).
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
/// construction (no k1/k2 split) instead of two -- faster, narrower-
/// domain tier, same pairing as `exp2`/`exp2_checked`. Valid while `k`
/// (see `exp10_checked`'s own doc comment) stays in `[-126,128)`.
/// Skipping the floor-adjust was tried and rejected here for the same
/// reason as plain `exp2` (round can land k=128 at the domain edge,
/// which a single field can't represent) -- see graveyard.md §exp/exp2.
#[inline(always)]
#[allow(clippy::approx_constant)] // g0's constant term is a fitted minimax
// coefficient near ln(2), not ln(2) itself (bit pattern deliberately differs)
pub fn exp10(x: f32) -> f32 {
    let (k, f) = exp10_reduction!(x);
    let exp2int = exp2int_field!(k);
    let q = exp2_q_poly!(f);
    fma(q, exp2int * f, exp2int)
}

// sin(x) ~= x + x^3*p(x^2) on [-pi/2, pi/2], degree-9 minimax (relative
// error ~6.1e-9), fitted with lolremez. Estrin evaluation, 2 fma chains.
//
// At x = +-0.0, x3 = y*x carries x's own sign, but p (the poly's leading
// coefficient at y=0) is a fixed negative constant, so p*x3 has the
// *opposite* sign to x at this one point. `fma(p, x3, x)` then adds two
// exactly-zero values of opposite sign, which IEEE754 defines to give
// +0.0 regardless of operand order, silently losing x's sign (the same
// mechanism behind the atan2(-0.0,+0.0) bug). `r.copysign(x)` (in
// sinf_poly below) fixes this for free: for every *nonzero* x in this
// poly's domain the leading `x` term dominates `p*x3`, so r's sign
// already equals x's -- copysign only changes the singular x=+-0.0 case,
// and is cheaper than an `x == 0.0` select (which measured ~12-17% worse
// throughput on sin/cos/tan). Verified per-caller: sin, sind, cosd all
// genuinely need the copysign (cosd would lose even-function
// sign-of-zero symmetry at its own crossings); sin_checked and sinpi's
// own separate zero guards, and cospi's never-negative `0.5-|r|`
// argument, make it redundant for them, so they call this raw version
// directly.
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

/// `sinf_poly_raw` plus the copysign fixup -- see `sinf_poly_raw`'s own
/// doc comment for the full story. Used by every caller except
/// `sin_checked`/`sinpi`, which call `sinf_poly_raw` directly since their
/// own separate zero-handling makes this copysign provably redundant for
/// them.
#[inline(always)]
fn sinf_poly(x: f32) -> f32 {
    sinf_poly_raw(x).copysign(x)
}

// pi split four ways for Cody-Waite reduction with fma. PI_A and PI_B
// carry only 8 and 9 significant bits: the trailing zeros keep the first
// two `x - q*PI_x` steps exactly representable for every |q| the
// magic-round `q` is defined over, so neither of them rounds at all.
// PI_C and PI_D then take the full f32 width, because from there on a
// step's rounding is only half an ulp *of the residual itself* (harmless
// even at sin's zeros, where the residual is what the answer is), while
// what the four words together fail to capture of pi is an absolute
// error scaled by q -- an unbounded *relative* one near those zeros, and
// the term that actually binds. So the tail words spend their bits on
// pi, not on trailing zeros.
const PI_A: f32 = 3.140625;
const PI_B: f32 = 0.0009670257568359375;
const PI_C: f32 = 6.278329465203569e-7;
const PI_D: f32 = 1.0780605906948477e-14;
const FRAC_1_PI: f32 = std::f32::consts::FRAC_1_PI;

// 1.5 * 2^23; adding this to |v| < 2^22 rounds v to the nearest integer
// in the low mantissa bits (round-to-nearest-even).
const ROUND_MAGIC: f32 = 12582912.0;

// The same trick two binades up: adding this to |v| < 2^24 lands v on
// the multiples of 4 instead of the integers. Four times the reach for an
// index whose parity is then known a priori (even) rather than read out
// of the sum -- which is the trade `sin` wants, since it has a cheaper
// place to get the one parity bit back from. See `frac_x_over_pi!`.
const ROUND_MAGIC_4: f32 = 50331648.0;

/// sin(x) via single-f32 range reduction, accurate for `|x| < 2^24 * pi`
/// (~5.27e7). `q = round(x/pi)` has to be an exactly-representable f32
/// integer for `pi_reduce_and_poly!`'s residual to mean anything, and 2^24
/// is the largest integer f32 counts by ones; past the limit q lands whole
/// integers off, shifting the residual by whole multiples of pi and
/// pushing it outside sinf_poly's fitted domain [-pi/2, pi/2] -- including
/// returning inf for some large finite x, since nothing here clamps the
/// residual. Use sin_checked for full-range gradual degradation instead of
/// this cliff; this version is much faster.
///
/// Inside that domain the residual stays in [-pi/2, pi/2] to within
/// ~7e-7 and max ulp is 2 all the way to the top -- see the
/// `frac_x_over_pi!` macro below, which is what the two extra fmas over
/// a naive single-word magic round buy.
// Shared by sin/cos: the Cody-Waite pi-split reduction (given each
// caller's own `q` -- sin's plain `round(x/pi)` vs cos's phase-shifted
// half-odd-integer) plus the `sinf_poly` call. Macro, not a fn --
// see exp_r_poly!.
macro_rules! pi_reduce_and_poly {
    ($x:expr, $q:expr) => {{
        let r = fma($q, -PI_A, $x);
        let r = fma($q, -PI_B, r);
        let r = fma($q, -PI_C, r);
        let r = fma($q, -PI_D, r);
        sinf_poly(r)
    }};
}

// An index `n` on the grid `$magic` implies, plus the leftover fraction
// `fc = x/pi - n` refined to two words of 1/pi.
//
// The refinement is what the two extra fmas over a naive single-word
// magic round buy. A single f32 `1/pi` is only good to `2^-25`
// *relative*, i.e. `|x|*2^-25/pi` absolute, which reaches 0.68 at the top
// of sin's domain: an index built from it alone lands a whole integer off
// whenever `x/pi`'s fraction sits within 0.68 of a half, and the residual
// leaves `sinf_poly`'s fitted `[-pi/2, pi/2]` by that same slice of pi.
// The answer stays self-consistent (parity comes from the same `q`), so
// that is not a wrong-branch cliff -- it is the poly being extrapolated,
// and it costs hundreds of ulp near the domain edge while `|x| <= 1e6`
// barely notices.
//
// `nb` is the magic round of the *exact* `x*FRAC_1_PI` product (that is
// what the fma buys), `f` recovers that product's fraction about `n`, one
// more fma folds in `RPI_LO` -- so `fc = x/pi - n` to ~`2.4e-7` absolute
// rather than `|x|*2^-25`. A third word of 1/pi would be dead weight
// here: it contributes at most `8e-9` over the domain, against `fc`'s own
// rounding. Neither error reaches the answer anyway -- both only decide
// *which* `q` the caller lands on, and only when `x/pi` is that close to
// the midpoint between two candidates, where either choice is
// self-consistent; `pi_reduce_and_poly!` then derives the residual from
// that `q` and `x` alone.
//
// `$magic` picks the grid `n` lands on. `ROUND_MAGIC` gives the plain
// integers, and with them `|n - x*FRAC_1_PI| <= 0.5`, which is what `cos`
// needs: its `n + copysign(0.5, fc)` is the nearest half-odd-integer to
// `x/pi` only while `|fc| <= 1`, and `fc` already carries `|x*RPI_LO|`
// (0.17 at the top of cos's domain) on top of that half. `ROUND_MAGIC_4`
// gives the multiples of 4 and four times the reach, at an `|fc|` up to
// ~2.7 -- fine for `sin`, whose own second magic round absorbs any `fc`,
// and which wants the reach because its `q` is a whole integer and so
// stays exact twice as far out as cos's half-odd one.
macro_rules! frac_x_over_pi {
    ($x:expr, $magic:expr) => {{
        let nb = fma($x, FRAC_1_PI, $magic);
        let n = nb - $magic;
        let f = fma($x, FRAC_1_PI, -n);
        (nb, n, fma($x, RPI_LO, f))
    }};
}

#[doc(alias = "sinf")]
#[inline(always)]
pub fn sin(x: f32) -> f32 {
    let (_, n, fc) = frac_x_over_pi!(x, ROUND_MAGIC_4);
    // A second, *fine* magic round, of `fc` this time -- an O(1) value, so
    // the plain integer grid has room to spare: `qb` is `ROUND_MAGIC +
    // round(fc)`, and `q = n + round(fc)` is round(x/pi) exactly, an
    // f32-exact integer for every `|q| <= 2^24` (which is exactly the
    // documented domain). `n - ROUND_MAGIC` is exact (both are multiples
    // of 4) and hangs off `n`, not off `fc`, so the whole `q` chain is no
    // deeper than the single round it replaces -- one add wider, but the
    // same latency.
    let nm = n - ROUND_MAGIC;
    let qb = fc + ROUND_MAGIC;
    let q = nm + qb;
    let s = pi_reduce_and_poly!(x, q);
    // sin(x) = (-1)^q * sin(r); n is a multiple of 4, so q's parity is
    // round(fc)'s, sitting in qb's lowest mantissa bit.
    let parity = qb.to_bits() << 31;
    f32::from_bits(s.to_bits() ^ parity)
}
/// cos(x). `q` here is the nearest *half-odd-integer* to `x/pi`, which
/// costs a mantissa bit that sin's whole-integer `q` does not, so cos is
/// documented over the narrower `|x| < 2^22 * pi` (~1.32e7) -- see
/// graveyard.md for what widening it would cost. Use cos_checked for
/// full-range gradual degradation.
#[doc(alias = "cosf")]
#[inline(always)]
pub fn cos(x: f32) -> f32 {
    // q = the half-odd-integer nearest x/pi, so r = x - q*pi is in
    // [-pi/2, pi/2]. Given `n = round(x/pi)` and `fc = x/pi - n`, that is
    // just `n + copysign(0.5, fc)` -- no second rounding. (The older
    // `round(x/pi - 0.5) + 0.5` needed one, and rounded twice getting
    // there: `fma(x, FRAC_1_PI, -0.5)` quantizes to `ulp(x/pi)`, already
    // 0.5 at the top of the domain, so the magic round downstream saw a
    // tie and cos left the poly's domain roughly twice as often as sin.)
    let (nb, n, fc) = frac_x_over_pi!(x, ROUND_MAGIC);
    let q = n + 0.5f32.copysign(fc);
    let s = pi_reduce_and_poly!(x, q);
    // cos(x) = (-1)^n * cos(r +- pi/2) = -+(-1)^n * sin(r), so the
    // half-turn contributes one more sign flip, taken when `fc >= 0`
    // (i.e. when `fc`'s sign bit is clear).
    let parity = (nb.to_bits() << 31) ^ (!fc.to_bits() & SIGN_MASK);
    f32::from_bits(s.to_bits() ^ parity)
}

/// [`sin`] with `q` from a *single* word of `1/pi` and a single magic
/// round -- the whole `frac_x_over_pi!` refinement dropped, and sin's
/// second round with it: two fmas and three adds cheaper (mca: 1.778 ->
/// 1.151 cyc/elem throughput, 64 -> 48 latency).
///
/// The price is paid twice over. That one round puts `q` on the plain
/// integer grid, so this is documented over `|x| < 2^22*pi` (~1.32e7),
/// two binades narrower than `sin`'s own domain. And inside even that,
/// `q` lands one integer off once `|x|*2^-25/pi` grows comparable to the
/// distance from `x/pi`'s fraction to a half, and `sinf_poly` is then
/// evaluated outside its fitted `[-pi/2, pi/2]`. That second error scales
/// with `|x|`, so it is a *narrower accurate range*, not a uniformly
/// looser one: at `|x| <= 1e6` it is indistinguishable from `sin` (0.0358
/// avg, 2 max, both), and only past there does it separate -- 0.0592 avg
/// and 219 max over its own `|x| < 2^22*pi` domain. Prefer this whenever
/// the argument is known to stay under ~1e6.
#[inline(always)]
pub fn sin_fast(x: f32) -> f32 {
    let qb = fma(x, FRAC_1_PI, ROUND_MAGIC);
    let q = qb - ROUND_MAGIC;
    let s = pi_reduce_and_poly!(x, q);
    let parity = qb.to_bits() << 31;
    f32::from_bits(s.to_bits() ^ parity)
}

/// [`cos`]'s counterpart to [`sin_fast`], same tradeoff and the same
/// single-word `q` (mca: 1.654 -> 1.406 cyc/elem, 61 -> 56 latency).
/// Degrades faster than `sin_fast` does: `k = round(x/pi - 0.5)` rounds
/// twice (see `cos`), which costs another half-integer of `q` error on
/// top of the single-word one, so max ulp reaches 2780 over the domain
/// against `sin_fast`'s 219. At `|x| <= 1e6` it is 3 against `cos`'s 2,
/// on an identical 0.0779 average.
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

/// sin(pi*x), argument in half-turns instead of radians. Unlike `sin`'s
/// own reduction (which needs a multi-constant Cody-Waite pi split
/// because pi itself isn't exactly representable), `sinpi`'s reduction
/// is *exact*: q = round(x) and r = x - q are both plain f32 operations
/// with no rounding error to correct for (q is an exact integer by
/// construction, and r = x - q is exact whenever |x| and |q| are within
/// a factor of 2 -- the same Sterbenz argument this crate already
/// relies on elsewhere, e.g. `exp`'s `x - k*LN2_HI`). `pi*r` then lands
/// exactly in `[-pi/2, pi/2]`, `sinf_poly`'s own fitted domain, so this
/// reuses that poly directly with no new fit needed. No accuracy cliff
/// anywhere in f32 (unlike `sin`'s ~1.3e7 or even `sin_checked`'s
/// ~1e13): past `|x| ~ 2^24`, every representable f32 is already an
/// *even* integer (ulp >= 2 there, so odd integers aren't even
/// representable), so `r` becomes exactly 0 and the result is exactly
/// `+0.0` everywhere out to `f32::MAX` (correct, since `sin(pi*integer)
/// == 0`, and parity is deterministically even, not just "untracked").
///
/// Uses `x.round_ties_even()` (native, full-range-correct), NOT the
/// magic-constant `x + 1.5*2^23` trick used elsewhere in this crate:
/// that trick is only exact for `|x| <= 2^22`, and `sinpi`/`cospi` round
/// the *raw*, unbounded input directly -- unlike every other magic-round
/// use here, which only ever rounds an already-reduced small value. (A
/// real bug lived here for exactly that reason: wrong for
/// `2^22 < |x| < 2^24`.) `round_ties_even` rather than `.round()` since
/// `q`'s tie-breaking rule doesn't affect correctness (`parity(q)`'s
/// sign correction self-compensates for whichever nearby integer `q`
/// lands on), and ties-even lowers to a single native `vroundps` while
/// `.round()`'s ties-away needs extra emulation instructions -- measured
/// ~10-14% faster, the same finding `remainder_ieee` made for its own
/// `q`.
#[inline(always)]
pub fn sinpi(x: f32) -> f32 {
    let q = x.round_ties_even();
    let r = x - q;
    // At x=-0.0: q is -0.0 too, so r=(-0.0)-(-0.0), which IEEE754 always
    // resolves to +0.0 -- the same "opposite-signed-zero op erases sign"
    // mechanism as sinf_poly_raw's own -0.0 note, one level further out
    // (r's lost sign means a copysign inside the poly would have nothing
    // left to copy, which is why this guard, not sinf_poly's copysign,
    // is what makes sinpi(-0.0) correct). Same select idiom as
    // log1p/log_2's x==0.0 case: compute the normal path unconditionally,
    // select x itself only at the singular zero point.
    let normal = sinf_poly_raw(std::f32::consts::PI * r) * fma(-2.0, parity(q), 1.0);
    if x == 0.0 { x } else { normal }
}

/// sinpi without the `x==0.0` guard (backlog idea #98): bit-identical to
/// [`sinpi`] everywhere except `x=-0.0`, where dropping the guard lets
/// the `q=r=-0.0` cancellation described in [`sinpi`]'s own doc comment
/// go uncorrected, returning `+0.0` instead of the correctly-signed
/// `-0.0`. Unlike the rejected `erf_unchecked` candidate, this is a pure
/// sign-of-zero nit at exactly one input (magnitude always correct, and
/// every other input -- including every other special value -- is
/// unaffected), not a diverging/wrong-magnitude failure mode; narrowing
/// the domain to exclude `-0.0` mirrors other `_unchecked` cores that
/// already exclude zero outright (e.g. `cbrt_unchecked`).
#[inline(always)]
pub fn sinpi_unchecked(x: f32) -> f32 {
    let q = x.round_ties_even();
    let r = x - q;
    sinf_poly_raw(std::f32::consts::PI * r) * fma(-2.0, parity(q), 1.0)
}

/// cos(pi*x), argument in half-turns -- see `sinpi`'s doc comment for why
/// this reduction is exact and shares `sinf_poly` directly, same
/// full-range-accurate (no cliff) guarantee, and the same magic-round
/// range-limit bug this version fixes.
///
/// Reduces on `k = round(x)`, exactly like `sinpi`/`tanpi` (*not* on
/// `round(x-0.5)`, where a real precision bug lived until 2026-07-30),
/// so `r = x - k` inherits sinpi's own exactness argument verbatim and
/// `cos(pi*x) = (-1)^k * cos(pi*r)` with `|r| <= 0.5`. The
/// quarter-turn reflection `cos(pi*r) = sin(pi*(0.5-|r|))` then hands
/// `sinf_poly` the distance to the *zero* rather than to the peak, and
/// `0.5-|r|` is Sterbenz-exact over `|r| >= 0.25` -- precisely the half
/// of the domain containing cospi's own zero -- so full relative
/// accuracy survives right up to the crossing. Subtracting inside the
/// small `[-0.5,0.5]` domain *before* scaling by pi (never after) is the
/// same ordering `tan_core` relies on near its own poles, and for the
/// same reason: `r`'s ulp there is finer than a pi-scaled value's.
///
/// Where `0.5-|r|` does round (`|r| < 0.25`) the result is within a
/// quarter turn of `+-1`, so the induced error stays under half an ulp
/// of a near-1 value -- the rounding lands where the function is flat,
/// not where it crosses.
///
/// `sinf_poly_raw` rather than `sinf_poly`: `0.5-|r|` is never negative
/// and never `-0.0`, so the copysign fixup that other callers need for
/// their own `+-0.0` arguments has nothing left to correct here. At an
/// exact half-integer `x` the poly argument *is* `+0.0`, and ties-even
/// rounding always sends `k` to the even neighbour there, so the sign
/// multiplier is `+1` and every one of cospi's zeros comes out `+0.0` --
/// matching IEEE 754-2019's `cosPi(n+1/2) = +0` for every `n`, with no
/// special case spent on it.
#[inline(always)]
pub fn cospi(x: f32) -> f32 {
    let k = x.round_ties_even();
    let r = x - k;
    let s = sinf_poly_raw(std::f32::consts::PI * (0.5 - r.abs()));
    s * fma(-2.0, parity(k), 1.0)
}

/// The normalized sinc function, `sin(pi*x)/(pi*x)` (DSP convention),
/// with the removable singularity at `x=0` handled directly (`sinc(0) =
/// 1`, the limiting value everywhere else already converges to). Built
/// directly on `sinpi`'s own exact, full-range reduction (see its doc
/// comment), so this is accurate across sinpi's *entire* domain -- not
/// just near zero, which is the part DSP users hand-rolling this
/// (`sin(pi*x)/(pi*x)` plus a manual near-zero branch) typically get
/// right, if anything. No cancellation risk in the division: for small
/// `x`, `sinpi(x)` is already close to `pi*x` (`sin(t) ~ t` near 0), so
/// the ratio stays well-conditioned throughout, including right up to
/// `x=0` itself. `sinc` is even (`sin(-pi*x)/(-pi*x) = sin(pi*x)/(pi*x)`
/// algebraically), which falls out for free here with no extra sign
/// handling needed.
#[inline(always)]
pub fn sinc(x: f32) -> f32 {
    let normal = sinpi(x) / (std::f32::consts::PI * x);
    if x == 0.0 { 1.0 } else { normal }
}

/// The unnormalized sinc function, `sin(x)/x` in radians (backlog idea
/// #131), the DSP/physics sibling of [`sinc`]'s own half-turn
/// convention. Built on `sin_checked` (not the unchecked `sin`) for the
/// same full-range reasoning `sinc` uses `sinpi`'s own exact reduction
/// for -- accurate across the whole domain, not just where the fast
/// tier's magic-round trick stays exact. Same removable-singularity and
/// no-cancellation-risk reasoning as `sinc` applies here too (`sin(t) ~
/// t` near 0, so the ratio stays well-conditioned all the way to `x=0`).
#[inline(always)]
pub fn sinc_unnormalized(x: f32) -> f32 {
    let normal = sin_checked(x) / x;
    if x == 0.0 { 1.0 } else { normal }
}

// tan(t) = t * Q(t^2), `t = pi*e` with e in [-0.25, 0.25], degree 6,
// two-group Estrin (`lo + u^4*hi`) so no group is deeper than the `u^4`
// it multiplies, leading coefficient pinned to exactly 1.0. An
// ulp-weighted minimax (LP) fit of tan(t)/t against t^2 -- weighted by
// `|d tanpi / d Q| / ulp(tanpi)` for whichever of `tan_core`'s two
// branches (direct `t*Q`, reflected `1/(t*Q)`) binds harder at each point
// -- then coordinate-descended over the f32 quantisation against the real
// chain. Not a transcription of any published algorithm's constants.
//
// Degree 6, not 5, because this poly is `tanpi`'s binding term and
// nothing else's. An oracle screen (the chain run with a correctly-rounded
// `tan(t)/t` in place of `Q`) puts the chain's own floor at max 3.7 ulp
// against degree 5's real 7.8, so the fit was carrying the error rather
// than the roundings; degree 6 lands near that floor and degree 7
// idealises 5x better again without converting.
//
// The grouping is load-bearing, not cosmetic. Horner-on-`u^2` at the top
// (`l0 + u2*(l1 + u2*l2)`) needs `l2` a full degree-2 group, so the extra
// term lands *on* the Estrin spine: same instruction count as this form
// minus one, but a fourth dependency level, which llvm-mca prices at
// tanpi latency 88.72 against this form's 76.00 (the shipped degree-5's
// own 78.88 for reference). It buys max ulp 5 -> 4. Not taken: 17% of
// tanpi's latency is a steep price for one ulp when the average is
// already better here (0.2267 against 0.2309).
#[inline(always)]
fn tan_poly(u: f32) -> f32 {
    let c: [f32; 7] = [
        1.0, 0.3333313, 0.13338836, 0.053408977, 0.024446711, 0.0030841262, 0.009410244,
    ];
    let u2 = u * u;
    let u4 = u2 * u2;
    let l0 = fma(c[1], u, c[0]);
    let l1 = fma(c[3], u, c[2]);
    let l2 = fma(c[6], u2, fma(c[5], u, c[4]));
    fma(u4, l2, fma(u2, l1, l0))
}

// tan(pi*e) (idea #128, tanpi only -- the same direct-poly idea applied
// to tand regressed a real fuzz-found precision bug and was reverted,
// see IDEAS.md), `e` the exact half-turn-fraction reduction `tanpi`'s
// own `r` already is. `theta=PI*e`, `aphi=PI*(0.5-|e|)`: the caller
// computes both, subtracting *before* scaling by `PI`, not after --
// `0.5-|e|` in the small, bounded `[-0.5,0.5]` domain is more precise
// than the mathematically-equivalent `FRAC_PI_2-|theta|` computed after
// scaling, since `theta`'s own ulp (fixed by its larger, ~pi/2-ish
// magnitude) is coarser than `e`'s, and near the pole that coarseness
// swamps the tiny quantity being computed. Both forms are Sterbenz-exact
// *given their own inputs*, so this isn't about avoiding rounding in the
// subtraction itself, only about which domain has finer-grained ulps to
// subtract in to begin with. Found by fuzzing (max ulp 8 vs 356266+
// for the exact same reflection formula, differing only in this
// subtract-then-scale vs scale-then-subtract order).
//
// `tan_poly` directly for `|e| <= 0.25`; past that, the cotangent
// reflection `tan(pi*e) = 1/tan(aphi)` (sign matching `e`/`theta`).
//
// `aphi == 0.0` exactly (`x` a genuine half-integer, `tan`'s true pole)
// needs an explicit override, found by fuzzing, not assumed: the
// reflected formula is a real `1.0/0.0` there, but `mulsign(.., theta)`
// ties its sign to `theta`'s own sign, which alternates with which
// half-integer `x` happens to be (tracks `theta`'s sign, not tan's
// actual pole-crossing direction) -- while the previous
// `sinpi(x)/cospi(x)` form (like this crate's own f64 `tanpi_ref` test
// reference, verified exhaustively over 1000 half-integers) always
// lands on the *same* sign, `-inf`, regardless of which half-integer:
// `sinpi`'s numerator sign and `cospi`'s zero sign there both alternated
// in lockstep, canceling to a constant ratio sign that `theta`'s sign
// alone doesn't reconstruct. (That lockstep was a property of the
// *pre-2026-07-30* `cospi`, whose zeros alternated `-0.0`/`+0.0`;
// today's `cospi` gives `+0.0` at every half-integer, so the same ratio
// would now alternate. `tanpi` never took the ratio route and is
// unaffected either way -- this override is what fixes the sign, and
// `tanpi_ref` keeps the older f64 shape deliberately, see its own note.)
#[inline(always)]
fn tan_core(theta: f32, aphi: f32) -> f32 {
    let ae = theta.abs();
    let direct = theta * tan_poly(theta * theta);
    let reflected = mulsign(1.0 / (aphi * tan_poly(aphi * aphi)), theta);
    let normal = if ae <= std::f32::consts::FRAC_PI_4 { direct } else { reflected };
    if aphi == 0.0 { f32::NEG_INFINITY } else { normal }
}

/// tan(pi*x), argument in half-turns (idea #128): `tan` has period 1 in
/// half-turns (unlike `sin`/`cos` individually, which flip sign every
/// integer -- see `sinpi`'s own doc comment for why *its* reduction needs
/// a parity correction that this one doesn't), so `x`'s exact `q=round(x)`/
/// `r=x-q` reduction (identical to `sinpi`'s own) already lands directly
/// on `tan(pi*x) = tan(pi*r)`, no sign combine needed. `tan_core` handles
/// the direct-poly/cotangent-reflection split and the pole itself (`r` a
/// half-integer distance from `x`, i.e. `x` itself a half-integer):
/// `aphi=0` there, so `1.0/(aphi*tan_poly(aphi*aphi))` is a real
/// `1.0/0.0`, giving the correctly-signed `+-inf` for free, the same
/// "division by a true zero, never a `0/0`" free case the previous
/// `sinpi(x)/cospi(x)` form also had.
#[inline(always)]
pub fn tanpi(x: f32) -> f32 {
    let q = x.round_ties_even();
    let r = x - q;
    let ar = r.abs();
    let normal = tan_core(std::f32::consts::PI * r, std::f32::consts::PI * (0.5 - ar));
    if x == 0.0 { x } else { normal }
}

/// sin(2*pi*x), full-turn argument (backlog idea #122): DSP phase
/// accumulators naturally count turns in `[0,1)`, not half-turns, so
/// this is the direct convention match for `sinpi`. `2.0*x` is an exact
/// power-of-two multiply (no rounding at all) for any finite `x` up to
/// `f32::MAX/2` -- past that it overflows to `+-inf`, and `sinpi`'s own
/// domain is total but not *that* total (`sinpi(inf)=NaN`), so this
/// silently gives up (`NaN`) for `|x|` in the top single octave of the
/// f32 range instead of a meaningful answer. A real, narrow gap, not a
/// bug: no caller with a genuine turn-count that large has a
/// meaningful "which fraction of a turn" answer left in f32 precision
/// anyway (`sinpi`'s own reduction already saturates the mantissa long
/// before that point).
#[inline(always)]
pub fn sin2pi(x: f32) -> f32 {
    sinpi(2.0 * x)
}

/// cos(2*pi*x) -- see [`sin2pi`]'s own doc comment for the convention
/// and domain caveat. `2.0*x` is exact, so this inherits `cospi`'s own
/// accuracy unchanged, including at `x = 0.25 + k/2` where `2x` lands
/// exactly on one of `cospi`'s true zeros: those stay full-relative-
/// accuracy (max 2 ulp, exhaustively), not a near-zero blowup.
#[inline(always)]
pub fn cos2pi(x: f32) -> f32 {
    cospi(2.0 * x)
}

/// tan(2*pi*x) -- see [`sin2pi`]'s own doc comment for the convention
/// and domain caveat. `x = 0.25 + k/2` puts `2x` on a genuine pole of
/// `tan`, where `tanpi` returns an exact `+-inf` rather than a large
/// finite value (see [`tanpi`] and `tan_core`).
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

/// sin(x*pi/180), argument in degrees. `180.0` has 18 trailing zero
/// mantissa bits (only needs 6 significant bits, since 180 = 4*45), so
/// unlike `pi` it needs no multi-word Cody-Waite split at all -- a
/// single constant `180.0` keeps `q*180.0` exact for `q` up to ~2^18,
/// far beyond any realistic input. `q = round(x/180)` (the coarse
/// `x*INV_180` multiply only needs to land on the right *integer*, same
/// reasoning as `exp10`'s own `k`), `d = x - q*180.0` stays small and
/// exact (Sterbenz), and `d*DEG_TO_RAD_SMALL` only multiplies the
/// *small* residual by an irrational constant (not the original,
/// possibly large, `x`) -- avoiding exactly the "recombine a small
/// correction with something large" pitfall `exp10`'s own doc comment
/// describes, since here the multiply happens *before* any combination
/// with `q`, not after. Reuses `sinf_poly` directly (its domain is
/// `[-pi/2,pi/2]`, and `d` in `[-90,90]` scaled by `DEG_TO_RAD_SMALL`
/// lands exactly there), same as `sinpi`.
///
/// Unlike `sinpi` (exact for the *entire* f32 range), this reduction is
/// only exact while `q*180.0` stays representable -- true for `|q|` up
/// to ~2^18 (`|x|` up to ~4.7e7), comfortably past any realistic input.
/// Past that, `d` stops being small and `d*DEG_TO_RAD_SMALL` could land
/// far outside `sinf_poly`'s fitted domain, so the radian residual is
/// clamped to `+-POLY_SAFE_BOUND` before `sinf_poly` sees it (same guard
/// as `sin_checked`/`cos_checked`). That guarantees always-finite output
/// for finite input out to `f32::MAX` -- *not* numerical correctness
/// past the exact boundary, same as any fast tier in this crate past its
/// documented range.
#[inline(always)]
pub fn sind(x: f32) -> f32 {
    let qb = fma(x, INV_180, ROUND_MAGIC);
    let q = qb - ROUND_MAGIC;
    let d = fma(-q, 180.0, x);
    let s = sinf_poly((d * DEG_TO_RAD_SMALL).clamp(-POLY_SAFE_BOUND, POLY_SAFE_BOUND));
    let parity = qb.to_bits() << 31;
    f32::from_bits(s.to_bits() ^ parity)
}

/// cos(x*pi/180), argument in degrees -- see `sind`'s doc comment for
/// the reduction (including its exactness limit and safety clamp); same
/// `-0.5`/`+0.5` quadrant-offset idiom `cos`/`cospi` already use.
#[inline(always)]
pub fn cosd(x: f32) -> f32 {
    let kb = fma(x, INV_180, -0.5) + ROUND_MAGIC;
    let q = (kb - ROUND_MAGIC) + 0.5;
    let d = fma(-q, 180.0, x);
    let s = sinf_poly((d * DEG_TO_RAD_SMALL).clamp(-POLY_SAFE_BOUND, POLY_SAFE_BOUND));
    let parity = !kb.to_bits() << 31;
    f32::from_bits(s.to_bits() ^ parity)
}

/// tan(x*pi/180), argument in degrees -- same `sind(x)/cosd(x)` ratio
/// construction as `tanpi`'s own original doc comment describes
/// (period-180 cancellation, poles handled for free by IEEE754
/// division). See `sind`'s own doc comment for the shared reduction's
/// exactness limit (`|x|` up to ~4.7e7) and safety clamp.
///
/// idea #128's direct-poly approach (see [`tanpi`], which uses it)
/// was tried here too and reverted: reusing `sind`'s own `d` for the
/// pole-distance calculation loses precision `d` never kept (it's
/// already reduced away from `x`, unlike `tanpi`'s `r`), and fixing
/// that needed a whole second, `cosd`-style independent reduction --
/// at which point real mca/quickbench numbers no longer showed a clean
/// win over this simpler ratio form. See IDEAS.md.
#[inline(always)]
pub fn tand(x: f32) -> f32 {
    sind(x) / cosd(x)
}

/// sind without the safety clamp (backlog idea #98): valid while `sind`'s
/// own reduction stays exact, `|x|` up to ~4.7e7 (see [`sind`]'s doc
/// comment) -- within that range `|d| <= 90` always (round-to-nearest-180
/// residual), so `d*DEG_TO_RAD_SMALL` never exceeds `pi/2`, nowhere near
/// `POLY_SAFE_BOUND` (1000.0): the clamp is provably a no-op in-domain,
/// only doing real work past the documented exactness limit. Past that,
/// unlike `sind`, no guarantee of even a finite result. Mirrors this
/// crate's other `_unchecked` cores (`cbrt_unchecked`, etc.); see [`sind`]
/// for the full-domain-safe version.
#[inline(always)]
pub fn sind_unchecked(x: f32) -> f32 {
    let qb = fma(x, INV_180, ROUND_MAGIC);
    let q = qb - ROUND_MAGIC;
    let d = fma(-q, 180.0, x);
    let s = sinf_poly(d * DEG_TO_RAD_SMALL);
    let parity = qb.to_bits() << 31;
    f32::from_bits(s.to_bits() ^ parity)
}

/// cosd without the safety clamp -- see [`sind_unchecked`] for the
/// rationale (same no-op-in-domain clamp removal) and [`cosd`] for the
/// full-domain-safe version.
#[inline(always)]
pub fn cosd_unchecked(x: f32) -> f32 {
    let kb = fma(x, INV_180, -0.5) + ROUND_MAGIC;
    let q = (kb - ROUND_MAGIC) + 0.5;
    let d = fma(-q, 180.0, x);
    let s = sinf_poly(d * DEG_TO_RAD_SMALL);
    let parity = !kb.to_bits() << 31;
    f32::from_bits(s.to_bits() ^ parity)
}

/// tand without the safety clamp -- same `sind_unchecked(x)/
/// cosd_unchecked(x)` ratio construction as [`tand`] itself, so it
/// inherits both this domain (`|x|` up to ~4.7e7) and its safety
/// tradeoff for free.
#[inline(always)]
pub fn tand_unchecked(x: f32) -> f32 {
    sind_unchecked(x) / cosd_unchecked(x)
}

// q = round(x/pi) must be an *exact* integer for x - q*pi to land
// accurately in [-pi/2, pi/2]. A single-f32 q is a binary either-or:
// either exactly the correctly-rounded integer, or (once |x| crosses q's
// exact-integer ceiling) off by a whole integer, shifting the residual by
// a whole multiple of pi and putting sinf_poly hopelessly outside its
// fitted domain -- a relocatable *cliff*, not a slope, no matter how q is
// rounded. And the reduction must then carry more bits of pi than a
// working f32 can hold. Both are f64's natural size, so the `_checked`
// tier does its reduction there and returns to f32 only for the poly.
//
// Why f64 rather than a wider f32 scheme -- three separate walls, all
// measured, none of them movable by spending more f32 ops:
//
//   * the *grid*: `ROUND_MAGIC + q` only has ulp 1 for |q| < 2^22, which
//     is exactly where `sin`/`cos`'s documented domain ends. Feeding the
//     unchecked chain an oracle q extends it to |q| < 2^24 and no
//     further -- past that an integer is not representable in an f32.
//   * the *Cody-Waite chain*: after removing word k the residual is
//     ~q*eps_k but is quantized at word k's lowest bit, so the step is
//     exact only while log2|q| <= 24 + (that word's trailing zero run in
//     pi's binary expansion) - 1. Pi has no zero run longer than 4 in
//     its first 100 bits, so *no* f32 splitting of pi survives past
//     |q| ~ 2^26; the shipped PI_A..PI_D bind at 2^25.6. Measured with
//     an exact two-word q: clean at 2^25, max 262145 ulp at 2^26.
//   * the *residual*: any coarse-then-fine two-stage scheme has to store
//     an intermediate of magnitude ~2^k in one f32, whose ulp there is
//     2^(k-24) against the 2^-24 the answer needs. So k <= 0.
//
// Everything past |q| ~ 2^25 therefore needs a wider-than-f32 residual.
// The double-f32 error-free-transform version this replaced was one such
// (and cost ~2x); f64 is the other, and is both cheaper and good three
// decades further out.
//
// The splitting rule below is what makes it cheap: x carries 24 bits, so
// giving each 1/pi word at most 26 significant bits makes every
// `x * IPI64_*` product *exact* in f64 with no error-free transform at
// all. On the pi side, PI64_HI's lowest bit at 2^-23 leaves the step-1
// residual (~q*2^-24.9, quantized at 2^-23) exact for log2|q| <= 53.9 --
// so the same analysis that caps f32 at 2^25 caps this at 2^53. What
// actually binds first is the parity extraction at 2^51; see
// `reduce_pi64`'s own comment.
//
// BOTH pi words must be POSITIVE -- pi truncated *downward*, not rounded
// to nearest. This is invisible from the arithmetic and is not a
// stylistic choice: with a negative low word, `-q * PI64_LO` at q = +0.0
// is `+0.0`, and `-0.0 + 0.0` is `+0.0`, so sin(-0.0) silently returns
// +0.0. Fuzzing and llvm-mca both missed that; edgecheck caught it.
const RPI_LO: f32 = 1.2841276486597053e-8;
// 1/pi and pi, each split so the leading word has <= 26 significant bits.
const IPI64_HI: f64 = 0.31830988079309464;
const IPI64_LO: f64 = 5.390696036528002e-09;
const PI64_HI: f64 = 3.1415926218032837;
const PI64_LO: f64 = 3.178650954705639e-08;
// 1.5 * 2^52: the f32 `ROUND_MAGIC` trick one exponent range up. Adding
// it to an exact integer |q| < 2^51 changes no bits of q but parks q's
// parity in bit 0 of the f64, where a 63-bit shift turns it straight
// into a sign mask.
const ROUND_MAGIC64: f64 = 6755399441055744.0;

/// Error-free transformation (Knuth's `two_sum`, backlog idea #184):
/// returns `(s, e)` such that `s = fl(a+b)` (the ordinary rounded sum)
/// and `s + e == a + b` *exactly*, for any finite `a`, `b` whose sum
/// does not overflow -- no ordering assumption on their magnitudes,
/// unlike [`quick_two_sum`]. Exactness holds right down through the
/// subnormal range (verified over the subnormal band; sums there are
/// multiples of `2^-149`, so nothing is lost). If `a + b` *does*
/// overflow, `s` is infinite and `e` is `NaN` -- the transform has no
/// meaningful result to report, so check the range first if your inputs
/// can reach it.
///
/// `e` recovers the rounding error the plain `+` silently dropped -- the
/// building block a compensated (double-float) reduction is made of,
/// exposed for users composing their own. This crate no longer builds
/// one itself: its own wide-range pi reduction used to, and now reduces
/// in f64 instead, which is both cheaper and good three decades further
/// out (see `reduce_pi64`). That is a statement about f32 pi reduction
/// specifically, not about the technique -- reach for these whenever the
/// next-wider float is unavailable or too slow, and prefer them to
/// reinventing the transform (often incorrectly:
/// the naive `e = (a+b) - a - b` loses exactly the precision this
/// exists to keep, since `(a+b)` has already rounded away the part
/// being recovered).
#[inline(always)]
pub fn two_sum(a: f32, b: f32) -> (f32, f32) {
    let s = a + b;
    let v = s - a;
    let e = (a - (s - v)) + (b - v);
    (s, e)
}

/// Cheaper 3-op form of [`two_sum`] (a.k.a. Fast2Sum, backlog idea
/// #184): `s + e == a + b` exactly *only if* `|a| >= |b|` -- violate
/// that and `e` is off by up to ~1 ulp of `s` instead of exact (still
/// bounded, unlike a plain `+` alone, just not the unconditional
/// guarantee `two_sum` gives for any `a`, `b`). Violating the ordering
/// deliberately can still be right, once the bounded error is verified
/// to be lost in the noise of everything else already inexact at that
/// site -- but that is a per-call-site judgment call, not something this
/// function itself checks, so callers should confirm `|a| >= |b|` (or
/// that bounded slop is acceptable) rather than assume it. Same overflow
/// caveat as [`two_sum`]: if `a + b` overflows, `s` is infinite and `e`
/// is infinite too (of the opposite sign), not a usable correction.
#[inline(always)]
pub fn quick_two_sum(a: f32, b: f32) -> (f32, f32) {
    let s = a + b;
    let e = b - (s - a);
    (s, e)
}

/// Error-free product (backlog idea #184): returns `(p, e)` with
/// `p = fl(a*b)` and `p + e == a * b` exactly, **provided the exact
/// product satisfies `2^-102 <= |a*b| <= f32::MAX`**. The multiplicative
/// counterpart to [`two_sum`]/[`quick_two_sum`], same purpose: recovers
/// the rounding a plain `a * b` alone drops, using a single `fma` rather
/// than a Dekker-style splitting chain.
///
/// Unlike [`two_sum`], the range precondition here is real and worth
/// checking -- it is not a formality, and both ends fail loudly rather
/// than degrading:
///
/// - **Overflow** (`|a*b| > f32::MAX`): `p` is infinite and `e` is
///   infinite of the opposite sign, so `p + e` is `NaN`.
/// - **Underflow** (`|a*b| < 2^-102`): `e` would need the 24 bits sitting
///   immediately below `p`, reaching down to `2^(E-47)` for `2^E <= |p|`.
///   Subnormals bottom out at `2^-149`, so those bits are simply not
///   representable and `e` is silently truncated; deep enough
///   (`|a*b| < 2^-149`) both `p` and `e` flush to zero while the true
///   product is nonzero.
///
/// The bound is `2^-102`, not the `2^-103` that the "product must be
/// normal" rule of thumb suggests: exactness was verified over ~41M
/// random in-range pairs against an `f64` reference with zero
/// violations, while `2^-103` still admits real failures.
#[inline(always)]
pub fn two_prod(a: f32, b: f32) -> (f32, f32) {
    let p = a * b;
    let e = fma(a, b, -p);
    (p, e)
}

/// `x - q*pi` reduced in f64, plus the sign the caller owes the result.
///
/// `q` is `round(x/pi)` for `HALF = false` (sin's grid) and the nearest
/// half-odd-integer for `HALF = true` (cos's), so `r` lands in
/// `[-pi/2, pi/2]` either way. Returns `(r, sgn)` with `sgn` already a
/// bare sign *mask* (0 or `SIGN_MASK`) rather than a `+-1.0` multiplier:
/// every caller either xors it into `r` before an odd polynomial or xors
/// it into a `1.0` to publish, and both are one instruction from the
/// mask while the multiplier form costs a select plus a multiply.
///
/// `HALF` is a const generic rather than an `f32` pre-offset because
/// `+ 0.0` is *not* a no-op LLVM may delete -- it turns `-0.0` into
/// `+0.0` -- so sin's callers would pay a real add for nothing.
///
/// Accurate to 2 ulp through `|x| < 1e13`, then degrading gradually:
/// max 906 at `[1e13,1e14)`, 1200 at `[1e14,1e15)`, ~55000 by `7e15`.
/// The cliff is at **`2^51*pi` = 7.07e15** and it is `ROUND_MAGIC64`'s
/// window, not the reduction's -- `n` stays an exact integer to `2^53`,
/// but `n + ROUND_MAGIC64` only has ulp 1 (so only carries the parity in
/// bit 0) while `|n| < 2^51`. A sign-dependent magic would buy one more
/// octave for two more ops, and is not worth it: an f32's own ulp at
/// 7e15 is already ~5e8 radians, 8e7 whole periods between neighbours.
/// Nothing in here bounds the *output* past that -- see `sin_checked`'s
/// `.clamp(-1,1)` for what keeps `|sin| <= 1` at every magnitude.
///
/// The cost profile assumes AVX-512 (this crate builds `-C
/// target-cpu=native`): the f32 lanes vectorize 256 bits wide, so the
/// f64 half widens to 512-bit ZMM at the *same instruction count*.
/// Without AVX-512 the f64 half needs two vectors per f32 vector and
/// roughly doubles.
#[inline(always)]
fn reduce_pi64<const HALF: bool>(x: f32) -> (f32, u32) {
    let xd = x as f64;
    // Both products are exact: x carries 24 bits, each word at most 26.
    let t = xd * IPI64_HI;
    let tl = xd * IPI64_LO;
    // `t - nh` is exact (both are multiples of ulp(t)), so `fr` is
    // `x/pi - nh` to a relative 2^-53 with no error-free transform --
    // this is why the f64 version is *cheaper* than the double-f32 one
    // it replaced, which needed a `two_prod` at exactly this step.
    let nh = t.round_ties_even();
    let fr = (t - nh) + tl;
    let d = fr.round_ties_even();
    // n = round(x/pi), exact for |n| < 2^53. Formed as `nh + d` rather
    // than rounding `t + tl` directly: that add would quantize at
    // ulp(t), which is already past 1 for the magnitudes this tier
    // exists to serve.
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
    let sgn = if HALF { par ^ ((!(fc.to_bits() >> 32)) as u32 & SIGN_MASK) } else { par };
    (r as f32, sgn)
}

// parity of an exact-integer float q via floor-based "mod 2" (q*0.5 and
// its floor stay exact once q is an integer), not `q as i64`: Rust's
// float-to-int cast is saturating, which LLVM can't lower to a single
// vector instruction.
#[inline(always)]
fn parity(q: f32) -> f32 {
    fma(-2.0, (q * 0.5).floor(), q)
}

// Bound for the reduced residual right before it enters sinf_poly, used by
// sind/cosd (whose own reduction has no `|result| <= 1` clamp downstream to
// fall back on). Once the reduction's precision runs out (|x| beyond the
// gradual-degradation range), the residual can grow large -- squaring that
// inside sinf_poly is where an earlier "returns inf for ordinary finite
// input" bug came from. sinf_poly's dominant term for large |r| is
// ~c3*r^9 (c3 ~ 2.6e-6), which overflows f32 around |r| ~ 8e4; 1000 leaves
// a large safety margin while still being far outside [-pi/2, pi/2], so a
// legitimately-reduced residual is never clipped.
// `.clamp` on a NaN residual (x itself nan or +-inf) returns nan unchanged,
// so sind/cosd(nan/inf) still correctly come out nan with no extra selects.
const POLY_SAFE_BOUND: f32 = 1000.0;

#[inline(always)]
pub fn sin_checked(x: f32) -> f32 {
    // sin is odd, so (-1)^q * sin(r) == sin((-1)^q * r): flip r's sign
    // bit *before* sinf_poly instead of negating its result after.
    // Bit-exact with a `s * (1.0 - 2.0 * par)` tail (both IEEE negation
    // and a multiply by exactly +-1 only ever flip the sign bit, never
    // round), but the mask depends solely on q -- ready long before r
    // exits the reduction -- so it hides in the reduction's shadow
    // instead of costing a real fma+mul on sinf_poly's tail.
    let (r, flip) = reduce_pi64::<false>(x);
    let r = f32::from_bits(r.to_bits() ^ flip);
    // `sinf_poly`, not the copysign-free `sinf_poly_raw`: the reduction
    // does carry `-0.0` through intact (that is what the positive pi
    // words buy), but `sinf_poly_raw` then loses it on its own -- its
    // last step is `fma(p, x3, x)` with `p < 0` and `x3 = -0.0`, so the
    // product is `+0.0` and `+0.0 + -0.0` is `+0.0`. Cheaper and more
    // local than the `if x == 0.0 { x }` guard this replaced, which sat
    // at the end of the function blaming the reduction for it.
    //
    // `.clamp(-1.0, 1.0)`: past `2^51*pi` (~7.07e15) the parity comes out
    // of the wrong bit and `q` eventually stops being an exact integer at
    // all, shifting r by multiples of pi and putting it wildly outside
    // sinf_poly's fitted domain (a degree-9 poly at |r| = 1000 is
    // ~2.6e21). Without this clamp sin_checked/cos_checked could silently
    // return values that large for legitimate (if extreme) finite input,
    // a genuine `|sin(x)| <= 1` violation. The clamp doesn't fix accuracy
    // that far out -- nothing can, an f32's own ulp there is ~5e8 radians
    // -- but it restores the one invariant every caller can rely on at
    // any magnitude.
    //
    // It also subsumes a `POLY_SAFE_BOUND` clamp on `r` itself: for every
    // f32 residual, `sinf_poly_raw` keeps the sign of `r` and only ever
    // grows past 1 in magnitude as |r| grows, so pre-clipping |r| to 1000
    // and letting |r| run free give bit-identical results after this clamp
    // (verified exhaustively over all 2^32 residuals). A NaN residual
    // (x itself nan or +-inf) survives both clamps unchanged, so
    // sin/cos(nan/inf) still come out nan with no extra selects.
    sinf_poly(r).clamp(-1.0, 1.0)
}
#[inline(always)]
pub fn cos_checked(x: f32) -> f32 {
    // q = the half-odd-integer nearest x/pi, so r = x - q*pi is in
    // [-pi/2, pi/2] and cos(x) = +-sin(r). Same sign-flip-before-the-poly
    // trick as sin_checked above (sinf_poly is odd in r too); the
    // half-turn's extra flip is already folded into the mask.
    let (r, flip) = reduce_pi64::<true>(x);
    let r = f32::from_bits(r.to_bits() ^ flip);
    // See sin_checked's own comment for why this clamp is needed and why
    // it is the only one needed: `q` stops being an exact integer past
    // |x| ~ 2^53*pi, and without this clamp cos_checked could silently
    // return values like 2.6e21 for legitimate finite input, violating
    // `|cos(x)| <= 1`.
    sinf_poly(r).clamp(-1.0, 1.0)
}

/// Reduce `x` modulo `pi`, accurate well beyond a single f32's
/// exact-integer range (backlog idea #88): the same f64 `reduce_pi64`
/// machinery `sin_checked` builds on, exposed for power users composing
/// their own periodic kernels (the "public checked pi-reduction API"
/// companion to a public `Df32` module, idea #87). Returns `(r, sign)`:
/// `r` is `x - q*pi` for the nearest integer `q`, unclamped (see
/// `sin_checked`'s own comment for why *it* clamps its result to
/// `[-1, 1]` -- that's specific to `sinf_poly`'s fitted domain, not a
/// general reduction contract, so it's on the caller here); `sign` is `+-1.0`
/// such that `sin(x) = sign * f(r)` for any odd `f` approximating `sin`
/// on `[-pi/2, pi/2]` -- a plain multiplier rather than a "which way
/// does this bool go" convention deliberately, found necessary by
/// fuzzing: an earlier version returned a bare parity bool and got its
/// own doc comment's if/else direction backwards for the cos variant
/// below (verified by reconstructing `cos_checked` from it and finding
/// a real sign mismatch, not just a doc typo -- the bool's *value* was
/// right, only which case meant "negate" was swapped). Only the
/// integer-`q` (sin) convention is exposed -- see
/// [`reduce_pi_half_checked`] for the half-odd-integer (cos) one;
/// `reduce_pi64`'s `HALF` flag has only those two settings and no
/// general pre-offset contract to expose.
#[inline(always)]
pub fn reduce_pi_checked(x: f32) -> (f32, f32) {
    let (r, flip) = reduce_pi64::<false>(x);
    // the mask is already in the sign-bit position, so `+-1.0` is one xor
    // away -- no compare, no select.
    (r, f32::from_bits(1.0f32.to_bits() ^ flip))
}

/// Reduce `x` modulo `pi`, offset by half a turn (backlog idea #88):
/// the `cos_checked`-style companion to [`reduce_pi_checked`], `q` here
/// the nearest integer to `x/pi - 0.5` (so `r = x - (q+0.5)*pi`, i.e.
/// `x` reduced around cosine's own zero-crossing grid instead of sine's
/// -- see `cos_checked`'s own doc comment for why this needs its own
/// `pre_offset=-0.5` reduction rather than a `+ pi/2` shift applied
/// after the fact, same "fold small corrections in before the reduction
/// loses the precision to represent them" reasoning throughout this
/// crate's own `_checked` tier). Same `(r, sign)` shape and the exact
/// same multiplier convention as [`reduce_pi_checked`] -- `cos(x) =
/// sign * f(r)` -- so callers needing both never have to remember two
/// different sign conventions, even though the underlying parity
/// check's sense really is inverted between sin's `q` and cos's `q+1`
/// exponent (folded in here, not left for the caller to get backwards).
#[inline(always)]
pub fn reduce_pi_half_checked(x: f32) -> (f32, f32) {
    let (r, flip) = reduce_pi64::<true>(x);
    (r, f32::from_bits(1.0f32.to_bits() ^ flip))
}

/// Largest magnitude [`wrap_pi`] can return: the largest `f32` whose
/// *exact* value is below `pi`, one ulp under `f32::consts::PI`.
///
/// `f32::consts::PI` itself is `pi + 8.7e-8`, so it is not in `(-pi,
/// pi]` and cannot be a legal `wrap_pi` result at either end -- which
/// makes the legal set symmetric, `|wrap_pi(x)| <= WRAP_PI_MAX`, and
/// this a single two-sided clamp rather than an asymmetric one. The
/// cost is that a true wrapped angle in `(WRAP_PI_MAX, pi]` comes back
/// as `WRAP_PI_MAX` (under 1 ulp low) instead of rounding up out of
/// range; exposed so range assertions can be written against the same
/// bound the implementation enforces.
pub const WRAP_PI_MAX: f32 = f32::from_bits(0x40490fda);

/// Wrap `x` (radians) to `(-pi, pi]` (backlog idea #126): the public
/// angle-normalization primitive robotics/geometry callers keep
/// reinventing, riding [`reduce_pi_checked`]'s own accurate-well-beyond-
/// a-single-f32 reduction rather than a naive `x - TAU*round(x/TAU)`
/// (which would need its own wide-range double-float treatment to avoid
/// `reduce_pi_checked`'s only *because it already exists* here).
/// `reduce_pi_checked` reduces mod `pi`, giving `r` in `[-pi/2,pi/2]`
/// and `sign=(-1)^q` for `q=round(x/pi)` -- if `q` is even, `x` and `r`
/// already sit in the same `2*pi` branch, so `r` alone is the answer;
/// if `q` is odd, `x = q*pi + r` sits a half-turn away, so the true
/// wrapped angle is `r +- pi` -- whichever keeps the result in
/// `(-pi,pi]`, i.e. `r - copysign(pi, r)`. That is written as a
/// copysign rather than an `r > 0.0` select because it is one
/// `vpternlogd` instead of a compare plus a blend, and the two are
/// bit-identical over the whole f32 line (the only input that could
/// tell them apart is an odd `q` with `r` exactly `+0.0`, which the
/// reduction never produces).
///
/// The `(-pi,pi]` range is a hard guarantee for every finite `x`, not
/// just where the reduction is accurate: past `|x| ~ 6.5e14` the
/// double-float `q` starts coming out off by whole integers, so `r`
/// leaves `[-pi/2,pi/2]` and the wrapped value degrades to noise (from
/// `|x| ~ 2.8e22` up, *every* input is in that regime). The result is
/// clamped to [`WRAP_PI_MAX`] so callers relying on the range still
/// get an angle rather than something like `1e24`; accuracy out there
/// is gone either way, and a caller who needs to know should range-check
/// its own input. `wrap_pi(nan)`/`wrap_pi(+-inf)` are `nan`, which the
/// clamp passes through unchanged.
#[inline(always)]
pub fn wrap_pi(x: f32) -> f32 {
    let (r, sign) = reduce_pi_checked(x);
    let normal = if sign > 0.0 {
        r
    } else {
        r - std::f32::consts::PI.copysign(r)
    };
    // See WRAP_PI_MAX: this is what makes the documented range true, and
    // it costs two instructions with no branch, so it vectorizes with
    // everything above it. It also absorbs the one in-range case the
    // `r - copysign(PI, r)` step gets wrong on its own: `PI` is
    // `pi + 8.7e-8`, so a small positive `r` lands on exactly `-PI`, which
    // is *outside* `(-pi, pi]` however it is rounded (16 inputs over the
    // whole f32 line). Clamping those to WRAP_PI_MAX is also the more
    // accurate answer, not just the in-range one.
    let normal = normal.clamp(-WRAP_PI_MAX, WRAP_PI_MAX);
    // x=-0.0 needs the same guard sinpi uses, for the same reason: inside
    // `reduce_pi_checked` the residual is formed by subtracting equal
    // signed zeros, which IEEE754 resolves to +0.0, so `r` arrives with
    // the sign already erased and there is nothing left downstream to
    // recover it from. Found by the special-value matrix (idea #164): over
    // |x| <= pi/2, where this function is otherwise *exactly* the
    // identity, -0.0 was the single sign anomaly in 2.14e9 inputs (the
    // only other two, at +-pi/2 itself, are legitimate 1-ulp boundary
    // effects where q flips to +-1 and the r+-pi branch rounds).
    if x == 0.0 { x } else { normal }
}

/// sin(r) for r already reduced to `[-pi/2, pi/2]` (backlog idea #127,
/// the trig analog of the still-backlog "public exp2_kf pre-reduced
/// primitive", idea #116): callers who already have their own accurate
/// reduction (e.g. via
/// [`reduce_pi_checked`]/[`reduce_pi_half_checked`], or their own) skip
/// re-deriving this crate's own reduction and its sign combine, since
/// [`sin`]/[`cos`]/[`sin_checked`]/[`cos_checked`] all ultimately
/// evaluate this exact same polynomial on their own respective `r` --
/// the two functions only differ in *reduction* and the parity-based
/// sign combine after, never in this step (verified directly against
/// `cos`'s own source: it computes `r` via its own `k=round(x/pi-0.5)`
/// convention, then calls the identical poly before its own separate
/// sign flip). No domain check: garbage in, garbage out for `|r| >
/// pi/2`, same contract as this crate's other `_unchecked`/prereduced
/// primitives.
#[inline(always)]
pub fn sin_prereduced(r: f32) -> f32 {
    sinf_poly(r)
}

/// cos-side companion to [`sin_prereduced`] (backlog idea #127) --
/// bit-identical to it. Kept as its own named function anyway: a caller
/// who reduced `x` via cosine's own `k=round(x/pi-0.5)` convention (not
/// sine's `round(x/pi)`) is thinking in terms of "the cosine step", and
/// a function named `sin_prereduced` is easy to miss even though it's
/// exactly what's needed -- the same discoverability reasoning as
/// [`cabs`]/[`carg`] being named aliases for [`hypot_checked`]/[`atan2`]
/// rather than expecting callers to know the underlying identity.
#[inline(always)]
pub fn cos_prereduced(r: f32) -> f32 {
    sinf_poly(r)
}

/// tan(x), full-range gradual degradation -- `sin_checked(x) /
/// cos_checked(x)`, mirroring `tanpi`/`tand`'s own plain-composition
/// pattern (period cancellation, poles handled for free by IEEE754
/// division). A shared-reduction fusion was tried for the closely
/// related `sincos_checked` and, despite a promising mca prediction
/// (~19% throughput), measured *slower* on real wall-clock -- this
/// session's one documented mca/hardware disagreement (see IDEAS.md) --
/// so plain composition is used here directly rather than re-attempting
/// that fusion. Like any full-range tan, ulp near a true pole (every
/// `pi`, and increasingly close together relative to float spacing as
/// `|x|` grows) is unbounded by nature, not a defect -- `cos_checked`
/// correctly rounding to a tiny value there is exactly what should
/// happen, the ratio just amplifies that tiny value's own relative
/// error the same way any division does near a zero denominator.
#[inline(always)]
pub fn tan_checked(x: f32) -> f32 {
    sin_checked(x) / cos_checked(x)
}

/// Core of cbrt for normal finite x: bit-trick seed (~3% error), then a
/// single degree-3 correction. d = s^3 - x is exact-ish via fma at any
/// scale, and x/s^3 == 1/(1+r) exactly for r = d/x, so
/// cbrt(x) = s * (1+r)^(-1/3), approximated by a minimax poly in r.
/// Degree 3 (4 coeffs) instead of degree 5 trades some accuracy (still
/// inside the 1 avg / 2 max ulp budget) for one fewer fma and one less
/// critical-path depth -- a large measured throughput win. Current:
/// avg ulp 0.28, max 3 (exhaustive). All intermediates are O(1) or O(x):
/// no overflow/underflow anywhere.
///
/// Callers needing a rescaled result (an exact power of two, or 1.0) must
/// multiply the *return value*, not thread a scale parameter through: an
/// in-function scale param touching 2 downstream ops gives LLVM's
/// vectorizer enough incentive to fully duplicate this entire function
/// for tiny vs. normal inputs instead of computing once and blending
/// (confirmed via llvm-mca disassembly -- 2x the fma/mul/div counts).
/// A single post-multiply by an exact power of two rounds identically to
/// pre-scaling, so keeping the scale outside is a pure codegen fix.
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
    if ax == 0 || ax >= EXPONENT_MASK { x + x } else { r }
}

/// cbrt without domain checks: valid for `x` normal (not denormal or
/// zero) and finite (not inf/nan) -- both signs are fine, `cbrt_normal`
/// already reapplies `x`'s own sign bit internally. Skips the denormal-
/// rescale select pair and the final zero/inf/nan propagation select
/// `cbrt`'s own doc comment describes paying on every call. Mirrors this
/// crate's other `_unchecked` cores (`log_2_unchecked`, etc.); see
/// [`cbrt`] for the full-domain-safe version.
#[inline(always)]
pub fn cbrt_unchecked(x: f32) -> f32 {
    cbrt_normal(x)
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
    // 1/(3y^2) without a second hardware division: y^3 ~ x (cbrt_normal
    // is within ~2 ulp) gives 1/y^2 = y/y^3 ~ y/x = |y|/a, so |y|*rcp3
    // approximates 1/(3y^2). Newton's quadratic convergence only needs
    // den to a handful of accurate bits -- this substitution is bit-exact
    // against a real division over the full accuracy sweep, and a small
    // throughput win (the FP divider is nearly idle while FMA/mul ports
    // are the bottleneck).
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

/// cbrt_accurate without the small/big-domain rescale selects or the
/// final zero/inf/nan-propagation select: valid for `x` already inside
/// cbrt_accurate's own safe range (roughly `2^-56` to `2^127`, see its
/// doc comment for why the rescale exists at all -- outside that range
/// `scale=1.0` is no longer correct and the double-float residual
/// misrounds or the Newton step's cube can overflow). Mirrors `cbrt`/
/// `cbrt_unchecked`'s own split, applied one tier up.
#[inline(always)]
pub fn cbrt_accurate_unchecked(x: f32) -> f32 {
    cbrt_accurate_normal(x, 1.0)
}

/// Core of rcbrt for normal finite x, the mirror of [`cbrt_normal`]:
/// bit-trick seed, then a single degree-3 correction on the same
/// `(1+r)^(-1/3)` shape. The seed goes *down* the exponent
/// (`K - ax/3`, not `ax/3 + K`), so `t ~ a^(-1/3)` directly and
/// `e = a*t^3 - 1` is formed by an fma with no division anywhere -- which
/// is the whole point, since `1.0 / cbrt(x)` pays two (`cbrt_normal`'s own
/// `1.0/a` reciprocal for its residual, plus the final reciprocal).
/// `a*t^3` is O(1) at every scale, so nothing here can overflow or
/// underflow. The correction poly is `rcbrt`'s own, not `cbrt_normal`'s:
/// the two seeds' residual ranges differ enough (`[-0.101, 0.104]` here,
/// exact over all normals) that reusing those coefficients would cost
/// several ulp.
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

/// x^(-1/3). Not `1.0 / cbrt(x)`: a dedicated inverse-cbrt seed (see
/// [`rcbrt_normal`]) reaches the answer with no division at all and one
/// fewer rounding, since the reciprocal is folded into the seed's
/// exponent rather than applied to a finished `cbrt`.
///
/// Composing did make every special case fall out for free; here they are
/// restored by one bit trick instead. `EXPONENT_MASK - ax` *is* the
/// magnitude of the reciprocal on exactly the inputs that need one:
/// `0 -> inf` and `inf -> 0`, and it wraps into the NaN encodings for a
/// NaN input, so a single subtract plus `x`'s own sign bit covers
/// `rcbrt(+-0) = +-inf`, `rcbrt(+-inf) = +-0` and `rcbrt(nan) = nan`.
/// Denormals rescale by `2^24` on the way in and `2^8` on the way out
/// (`(2^24)^(-1/3) = 2^-8`), the same shape as `cbrt`'s own rescale --
/// and, as there, applied to the *return value* rather than threaded
/// through as a parameter, so LLVM doesn't duplicate the kernel per
/// branch.
#[inline(always)]
pub fn rcbrt(x: f32) -> f32 {
    let ax = x.to_bits() & !SIGN_MASK;
    let tiny = ax < 0x0080_0000; // denormal or zero
    let xs = if tiny { x * 16777216.0 } else { x };
    let scale = if tiny { 256.0 } else { 1.0 };
    let r = rcbrt_normal(xs) * scale;
    let spec = f32::from_bits(EXPONENT_MASK.wrapping_sub(ax) | (x.to_bits() & SIGN_MASK));
    if ax == 0 || ax >= EXPONENT_MASK { spec } else { r }
}

/// Bit-trick `cbrt` seed plus two rational refinement steps (backlog idea
/// #188's remaining members): max relative error **~5e-6** for `x` in
/// `[1e-30, 1e30]` -- far tighter than its `_approx` siblings, since unlike
/// them it does refine. Outside that band it degrades completely rather than
/// gracefully: relative error reaches **100%** at the smallest normal
/// (`1.175e-38`), because the seed's `+ 0x2a509849` carries the exponent
/// field out of the normal range with no rescale to bring it back, and
/// `cbrt_approx(-0.0)` is `NaN`. Deliberately outside the crate's 0.5/2 ulp
/// budget -- use [`cbrt`] or [`cbrt_accurate`] for real work; this is kept
/// for the `_approx_plot`/`_error` test suite. Bounds measured by
/// `examples/approx_bounds.rs`.
pub fn cbrt_approx(x: f32) -> f32 {
	let y = f32::from_bits(0x2a509849u32 + (x.to_bits() / 3));
	let y = (x + 2.*(y*y)*y) / (3.*(y*y));
    (2.*x*y + (y*y)*(y*y))/(x + 2.*(y*y)*y)
}
/// Classic single-bit-trick `sqrt` seed (backlog idea #188), no Newton
/// refinement: max relative error **~4.5%** (4.484%), and unlike its
/// siblings that figure holds over the *whole* positive-normal range, not
/// just a mid-range band -- halving the exponent field via `>> 1` can't
/// leave it. Deliberately outside the crate's 0.5/2 ulp budget; use
/// `f32::sqrt` (a single hardware instruction) for real work. Bound measured
/// by `examples/approx_bounds.rs`.
pub fn sqrt_approx(x: f32) -> f32 {
    f32::from_bits(0x1FBD22DF + (x.to_bits() >> 1))
}
/// Classic single-bit-trick reciprocal seed (backlog idea #188), no Newton
/// refinement: max relative error **~6.6%** (6.557%) for `x` in
/// `[1e-30, 1e30]`. The band matters here -- near `f32::MAX` the true
/// reciprocal is denormal and this has no range handling, so relative error
/// grows without bound (~1e80 at `x = 3.29e38`), and it does not saturate
/// sanely at the specials either (`rcp_approx(+0)` is `1.59e38`, not `inf`;
/// `rcp_approx(+inf)` is negative). Deliberately outside the crate's 0.5/2
/// ulp budget; use a plain `1.0 / x` for real work. Bounds measured by
/// `examples/approx_bounds.rs`.
pub fn rcp_approx(x: f32) -> f32 {
    f32::from_bits(0x7EEF370B - x.to_bits())
}
/// Classic single-bit-trick `2^x` seed (backlog idea #188), no Newton
/// refinement: max relative error ~6.1% over `|x| < 120`. Deliberately
/// rough (this crate's real `exp2`/`exp2_checked` exist for accurate
/// use) -- kept for the `_approx_plot`/`_error` test suite, demonstrating
/// the exponent-field bit-manipulation idiom those real functions build
/// on, not a candidate tier of its own.
pub fn exp2_approx(x: f32) -> f32 {
    -f32::from_bits((x + 383.).to_bits() << 8)
}

/// Classic single-bit-trick `log2(x)` seed (backlog idea #188), no
/// Newton refinement: max absolute error ~0.086 over positive normal
/// `x` (relative error is only meaningful away from `log2(x)=0` at
/// `x=1`, same "ulp isn't meaningful near a true zero" caveat this
/// crate's other log-family functions document). Deliberately rough --
/// this crate's real `log_2`/`ln`/`log2_unchecked` exist for accurate
/// use -- kept for the `_approx_plot`/`_error` test suite.
pub fn log2_approx(x: f32) -> f32 {
    f32::from_bits((x).to_bits() >> 8 | 256_f32.to_bits()) - 383.
}

/// Quake-style `1/sqrt(x)` bit-trick seed (backlog idea #188), no Newton
/// refinement: max relative error ~4.8% over positive normal `x` (the
/// famous version adds one Newton iteration to reach ~0.2% -- this one
/// deliberately doesn't, see [`rsqrt`]'s own doc comment). Kept for the
/// `_approx_plot`/`_error` test suite, not a candidate replacement.
pub fn rsqrt_approx(x: f32) -> f32 {
    f32::from_bits(0x5F33E79F - (x.to_bits() >> 1))
}

/// Latency-optimal `cbrt` approximation: two bit-trick seeds -- one for
/// `x^(1/3)`, one for `x^(-2/3)/3` -- and two coupled Newton steps
/// `s <- s + r*(x - s^3)`, spelled `fma(s*s, s*-r, fma(r, x, s))` so each
/// step is only two dependent FP levels deep. `r` is never refined, which
/// is the whole design: it is what keeps a step at two levels instead of
/// the four an inverse-cbrt iteration needs, and it is also the accuracy
/// floor (the residual after two steps is `~u0*(u0+v)*v` for seed errors
/// `u0`, `v` of a few percent each, so no amount of constant tuning gets
/// this near an ulp -- `cbrt_normal` is 0.28 avg / 3 max for 25% more
/// latency, and is what to reach for unless the dependency chain is the
/// binding constraint).
///
/// Positive, finite, normal `x` only -- no zero/negative/denormal/inf/nan
/// handling, unlike [`cbrt`]/[`cbrt_unchecked`]. Approx tier, deliberately
/// outside the crate's 0.5/2 ulp budget: avg ulp 57.4, max ulp 554
/// (positive-normal domain), llvm-mca latency 28.05 cyc and throughput
/// 0.839 cyc/elem -- against `cbrt_unchecked`'s 35.06 / 0.906.
///
/// Both seeds scale the exponent field by a shift-multiply rather than an
/// exact `bits/3`, and both Newton steps share one reciprocal seed. Two
/// tuned variants that improve on this -- a second reciprocal-seed offset,
/// and the exact `bits/3` -- are measured in graveyard.md and not taken:
/// each buys 16-31% of the max ulp but gives back part of the latency or
/// throughput edge that is the only reason this tier exists.
#[inline(always)]
pub fn cbrt_fast(x: f32) -> f32 {
    let s = f32::from_bits(0x2a4d_def1u32.wrapping_add((x.to_bits() >> 16) * 0x5556u32));
    let r = f32::from_bits(0x68ff_2381u32.wrapping_sub((x.to_bits() >> 16) * 0xaaac));
    let s = fma(s * s, s * -r, fma(r, x, s));
    fma(s * s, s * -r, fma(r, x, s))
}

/// `x * sign(y)` (backlog idea #184): an xor of sign bits, genuinely
/// different from [`f32::copysign`], not just a naming variant --
/// `copysign` *replaces* `x`'s sign outright with `y`'s, while
/// `mulsign` *combines* them (flips `x`'s sign iff `y` is negative,
/// leaves it alone iff `y` is positive). They agree only when `x >= 0`;
/// for negative `x` they diverge: `mulsign(-2.0, -3.0) == 2.0` (two
/// negatives flip back to positive) but `copysign(-2.0, -3.0) == -2.0`
/// (unconditionally takes `-3.0`'s sign). Ported from jodiemath's
/// `mulsign`; used throughout this crate wherever a sign needs to be
/// *conditionally toggled* based on some other value's sign (parity
/// flips, odd-function reconstructions) rather than *pinned* to it.
#[inline(always)]
pub fn mulsign(x: f32, y: f32) -> f32 {
    f32::from_bits(x.to_bits() ^ (y.to_bits() & SIGN_MASK))
}

const LN_2: f32 = std::f32::consts::LN_2;
const LOG2_E: f32 = std::f32::consts::LOG2_E;
// log2(e) - LOG2_E, i.e. what the f32 constant is short by. Only
// `log2p1_df` needs it: everywhere else log2(e) multiplies a value whose
// own error already swamps 2^-24 of it, but there the product *is* the
// answer whenever `1+x` rounds back to `1`.
const LOG2E_LO: f32 = 1.9259630e-8;
const FRAC_PI_2: f32 = std::f32::consts::FRAC_PI_2;
const FRAC_PI_4: f32 = std::f32::consts::FRAC_PI_4;

// Cody-Waite split of ln(2): LN2_HI keeps its low 9 mantissa bits zeroed, so
// k*LN2_HI (k an exact small integer, this crate's log_2_normal decomposition
// never produces |k| past a couple hundred) is *exact* -- no rounding at all,
// confirmed by brute force for k in [-300, 300]. LN2_LO is the f32-rounded
// residual (LN2 - LN2_HI as f64, then rounded). This is the same trick
// PI_A..PI_D and PI64_HI/PI64_LO use for pi's own splits.
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
#[doc(alias = "logf")]
#[doc(alias = "log")]
#[inline(always)]
#[allow(clippy::neg_cmp_op_on_partial_ord)]
pub fn ln(x: f32) -> f32 {
    log_family_wrapper!(x, ln_normal)
}

/// Core of ln for positive normal finite x only -- see log_2_normal, same
/// contract. Reuses log_2_normal's exact decomposition (s = m - 1) and
/// poly *shape*, but with coefficients fitted for ln directly (log_2's own
/// c[i] * LN_2, each individually rounded to f32) and a Cody-Waite combine
/// with k instead of a single fma: `k*LN2_HI` is exact (see LN2_HI's own
/// comment), and `k*LN2_LO` adds back the tiny residual LN2_HI dropped, so
/// the correction's own rounding only ever affects a small term instead of
/// the whole (k-dominated) result, unlike the naive `log_2(x) * LN_2` where
/// the second rounding scales everything.
///
/// Both correction terms are folded in as `fma(p, s, fma(k, LN2_LO, k_hi))`
/// -- two fma against the older `fma(p, s, k_hi) + k*LN2_LO`'s fma + mul +
/// add, output bit-identical (verified exhaustively for `ln` and `log1p`).
/// The association is load-bearing and was measured both ways: folding the
/// LO word into the `k` term *first* is what wins, while doing the poly
/// first reproduces a +1-cycle latency regression. `log10_normal`
/// deliberately does *not* mirror this -- same transform, but there it
/// costs real accuracy; see IDEAS.md.
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn ln_normal(x: f32, koff: f32) -> f32 {
    // Degree 8 (backlog idea #4), not `log_family_normal!`'s shared
    // degree 9 -- a real least-squares refit (scipy), not the old
    // coefficients with the last one dropped. Real max ulp (exhaustive)
    // stayed at 3 (same as degree 9): the degree-9 poly's own idealized
    // fit error was ~5e-9 (~0.05 ulp) against a real 3-ulp max that's
    // entirely rounding-chain-dominated (the reduction's own fma combine),
    // so degree 9 had spare accuracy well beyond what it needed --
    // confirmed before implementing, not assumed, by checking the
    // degree-8 refit's own idealized error (~2.4e-8, ~0.3 ulp-equivalent)
    // stays far under both 1 ulp and that 3-ulp floor. `log10` shares
    // this same reduction shape and same "already at its degree-9
    // optimum" diagnosis, but was not independently re-verified at
    // degree 8 -- do that before assuming this transfers.
    let e = (x.to_bits() as i32).wrapping_sub(0x3f3504f3) >> 23;
    let m = f32::from_bits((x.to_bits() as i32).wrapping_sub(e << 23) as u32);
    let k = e as f32 + koff;
    let s = m - 1.0;
    let c: [f32; 9] = [
        1.0,
        -0.49999994,
        0.33334,
        -0.25001356,
        0.19962999,
        -0.16583247,
        0.14908722,
        -0.14196032,
        0.08632387,
    ];
    let s2 = s * s;
    let s4 = s2 * s2;
    let l0 = fma(c[1], s, c[0]);
    let l1 = fma(c[3], s, c[2]);
    let l2 = fma(c[5], s, c[4]);
    let l3 = fma(c[7], s, c[6]);
    let r0 = fma(l1, s2, l0);
    let r1 = fma(l3, s2, l2);
    // c[8] (the odd 9th coefficient, degree 8) folds into r1 at the s4
    // level instead of needing its own s8 = s4*s4 level: r1b*s4 ==
    // (c4+c5*s+c6*s2+c7*s3+c8*s4)*s4, correctly placing c8 at s8 with
    // one more fma but no new multiply, keeping the same op count
    // `log_family_normal!`'s degree-9 form uses for s2/s4 alone.
    let r1b = fma(c[8], s4, r1);
    let p = fma(r1b, s4, r0);
    let k_hi = k * LN2_HI; // exact, see LN2_HI's comment
    fma(p, s, fma(k, LN2_LO, k_hi))
}

/// ln without domain checks: valid for positive normal finite x only, see
/// log_2_unchecked for the general rationale (same fast/full-safety split,
/// same reason for the `_unchecked` suffix instead of `ln`/`ln_checked`).
#[inline(always)]
pub fn ln_unchecked(x: f32) -> f32 {
    ln_normal(x, 0.0)
}

/// log10(x), same Cody-Waite-combine approach as ln (see ln's own doc
/// comment for why this avoids the naive `log_2(x) * LOG10_2`'s double
/// rounding).
#[doc(alias = "log10f")]
#[inline(always)]
#[allow(clippy::neg_cmp_op_on_partial_ord)]
pub fn log10(x: f32) -> f32 {
    log_family_wrapper!(x, log10_normal)
}

/// Core of log10 for positive normal finite x only -- see ln_normal, same
/// approach (own degree-8 refit, own doc comment explains the shared
/// Estrin-folding trick) with coefficients fitted directly for log10
/// (not log_2's own coefficients times a constant, and not ln_normal's
/// degree-8 refit times a constant either -- least-squares was run
/// fresh against log10(1+s)/s). Real max ulp (exhaustive) confirmed
/// unchanged at degree 8, same as ln_normal's own verification.
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn log10_normal(x: f32, koff: f32) -> f32 {
    let e = (x.to_bits() as i32).wrapping_sub(0x3f3504f3) >> 23;
    let m = f32::from_bits((x.to_bits() as i32).wrapping_sub(e << 23) as u32);
    let k = e as f32 + koff;
    let s = m - 1.0;
    let c: [f32; 9] = [
        0.43429446,
        -0.21714722,
        0.14476772,
        -0.10857951,
        0.086698204,
        -0.07202013,
        0.06474776,
        -0.06165259,
        0.037489995,
    ];
    let s2 = s * s;
    let s4 = s2 * s2;
    let l0 = fma(c[1], s, c[0]);
    let l1 = fma(c[3], s, c[2]);
    let l2 = fma(c[5], s, c[4]);
    let l3 = fma(c[7], s, c[6]);
    let r0 = fma(l1, s2, l0);
    let r1 = fma(l3, s2, l2);
    let r1b = fma(c[8], s4, r1);
    let p = fma(r1b, s4, r0);
    let k_hi = k * LOG10_2_HI; // exact, see LN2_HI's comment (same trick)
    fma(p, s, k_hi) + k * LOG10_2_LO
}

/// log10 without domain checks: valid for positive normal finite x only,
/// see log_2_unchecked for the general rationale.
#[inline(always)]
pub fn log10_unchecked(x: f32) -> f32 {
    log10_normal(x, 0.0)
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
/// as sinf_poly's own `-0.0` bug. Fixed with a trailing `if x == 0.0
/// { x } else { normal }` select (compute the normal path
/// unconditionally, the "no early returns" idiom, so array loops keep
/// auto-vectorizing): log1p is odd and monotonic through the origin, so
/// for every nonzero x `normal`'s sign already equals x's, making the
/// select a no-op everywhere except the singular zero point.
// `log1p` without its trailing `x == 0.0` signed-zero select -- correct
// for every other input, and the whole function for callers whose own
// select already discards this value at zero. Macro, not a fn -- see
// exp_r_poly!.
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
    if x == 0.0 { x } else { normal }
}

/// log1p(x) - x (backlog idea #145), the gamma/Poisson-kernel primitive
/// where the naive form cancels for small `x` (`log1p(x) ~ x` there --
/// same cancellation class `log1p`/`expm1` themselves exist to avoid,
/// one level further out, and *worse*: measured directly, `log1p(x)-x`'s
/// own real f32 cancellation error stays in the tens-of-ulp range all
/// the way out to `|x|~0.8-1.0`, not just a narrow band near zero the
/// way `sqrt1pm1`'s analogous correction is). `log(1+x) = x - x^2/2 +
/// x^3/3 - ...`, so `log1pmx(x) = -x^2/2 * Q(x)` with `Q(0)=1`; `Q`'s
/// own Taylor series converges too slowly for a pure truncation to
/// reach f32 precision at any useful radius (measured: even a degree-5
/// pure-Taylor `Q` is ~4.6 ulp-equivalent at just `|x|<0.1`), so this
/// uses a real minimax refit over `|x|<0.5` (degree 12) instead.
///
/// That one series then covers the *whole* domain, not just `|x| < 0.5`.
/// Past there the argument handed to it is `w = m - 1`, the residual of
/// `ln`'s own `u = 1+x = 2^k * m` reduction (`m` in `[1/sqrt2, sqrt2)`,
/// so `w` is exact and already inside the fitted range), and
/// `ln(u) = k*ln2 + w + log1pmx(w)` rebuilds the answer from it. So no
/// separate `ln` polynomial is evaluated here at all, and no arm ever
/// rounds `ln(1+x)` into a single f32 and *then* subtracts `x` -- that
/// is the form that hands `ln`'s own few ulp to a cancellation whose
/// amplification, `|x / log1pmx(x)|`, peaks at ~5.3 just past
/// `|x| = 0.5`. `x=+inf` is the one input the reduction alone
/// mishandles (`inf-inf` is indeterminate; the true limit is `-inf`,
/// `log` growing arbitrarily slower than `x`), corrected with a
/// trailing override, same mechanism as `sqrt1pm1`'s own `x=+inf` fix.
#[inline(always)]
#[allow(clippy::neg_cmp_op_on_partial_ord)]
pub fn log1pmx(x: f32) -> f32 {
    let u = 1.0 + x;
    let c = x - (u - 1.0);
    let corr = c / u;
    let corr = if corr.is_finite() { corr } else { 0.0 };
    // ln_normal's own reduction, done here rather than by calling it:
    // u = 2^k * m with m in [1/sqrt2, sqrt2), so w = m - 1 is exact
    // (Sterbenz) and lands in [-0.293, 0.415] -- inside the same |z| < 0.5
    // the poly below is fitted over. That is what lets one poly serve both
    // arms.
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
    // `z` is only available after the reduction below has produced `w`, so
    // this chain sits on the critical path and a 12-deep serial Horner
    // would dominate it. Estrin instead -- five fma levels rather than
    // twelve, for two extra multiplies -- but only *below* the leading
    // term: `C1..C12` are grouped, then the `1.0` is added by a single
    // trailing fma. Grouping the `1.0` in as well costs real accuracy,
    // because that is the one step whose rounding lands at full weight;
    // every rounding inside the grouped part enters attenuated by `z`.
    // `|z| < 0.5` throughout, so no grouped partial sum can overflow the
    // way an Estrin split can at an infinite argument.
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
    // ln(u) = k*ln2 + ln(m) = k*ln2 + w + log1pmx(w), so
    //   log1pmx(x) = ln(u) - x + c/u = (k*LN2_HI - x + w) + p + (k*LN2_LO + c/u)
    // with p the same series value the |x| < 0.5 arm returns directly.
    // Every large term is now differenced *before* anything small is added:
    // k*LN2_HI is exact, w is exact, and over the whole band where
    // |x/log1pmx(x)| is big (it peaks at ~5.3 just past |x| = 0.5) both
    // partial differences are Sterbenz-exact too, so nothing is rounded at
    // ln(u)'s scale and then amplified. Evaluating ln(u) as a single f32
    // and subtracting x afterwards -- the direct `log1p(x) - x` -- instead
    // hands that cancellation ln's own few ulp, magnified fourfold.
    // Reusing the series for ln(m) is also why no separate `ln` poly is
    // evaluated here at all.
    let big = log_family_edges!(u, {
        let t = fma(k, LN2_HI, -x) + w;
        t + (p + fma(k, LN2_LO, corr))
    });
    let normal = if x.abs() < 0.5 { p } else { big };
    if x == f32::INFINITY { f32::NEG_INFINITY } else { normal }
}

/// log2(1+x) (C23 `log2p1`). Same `u = 1+x` / Sterbenz-exact correction
/// `c = x - (u-1)` trick as `log1p`, just converted to log2 units:
/// `d(log2)/du = 1/(u ln2)`, so the correction term is `(c/u) * LOG2_E`
/// instead of plain `c/u` -- the extra multiply only scales the already-
/// small correction, not the (k-dominated) whole result, so it doesn't
/// reintroduce the double-rounding problem `ln`'s own doc comment
/// describes for the naive `log_2(x) * LN_2`. Same degenerate-correction
/// guard (`u == 0` or `u == inf` collapse `c/u` to NaN even though
/// `log_2(u)` alone is already the right answer there) and the same
/// trailing `x == 0.0` select for the opposite-signed-zero-addition trap
/// -- both copied from `log1p` verbatim. Avg ulp 0.102, max 3
/// (exhaustive).
#[allow(clippy::neg_cmp_op_on_partial_ord)]
#[inline(always)]
pub fn log2p1(x: f32) -> f32 {
    let u = 1.0 + x;
    let c = x - (u - 1.0);
    let corr = (c / u) * LOG2_E;
    let corr = if corr.is_finite() { corr } else { 0.0 };
    let normal = log_family_wrapper_no_denormal!(u, log_2_normal) + corr;
    if x == 0.0 { x } else { normal }
}

/// log10(1+x) (C23 `log10p1`), completing the C23 set next to `log2p1`
/// above: identical Sterbenz-exact-correction structure, just converted
/// to log10 units (`(c/u) * LOG10_E` instead of `* LOG2_E`) and calling
/// `log10` instead of `log_2` for the dominant term.
#[allow(clippy::neg_cmp_op_on_partial_ord)]
#[inline(always)]
pub fn log10p1(x: f32) -> f32 {
    let u = 1.0 + x;
    let c = x - (u - 1.0);
    let corr = (c / u) * std::f32::consts::LOG10_E;
    let corr = if corr.is_finite() { corr } else { 0.0 };
    let normal = log_family_wrapper_no_denormal!(u, log10_normal) + corr;
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
/// `r = x - k*ln2` done as an exact Cody-Waite reduction (`LN2_HI`/
/// `LN2_LO`, the same split `ln`/`log10` already use -- `k*LN2_HI` is
/// exact for this domain's k, and `x - k*LN2_HI` is exact by Sterbenz
/// since `k*ln2` tracks `x` closely), then a dedicated degree-5 minimax
/// poly for `e^r` on `[-ln2/2, ln2/2]`, scaled by `2^k`. Inherits exp2's
/// unchecked domain: only accurate while `x*log2(e)` stays in
/// `[-126, 128)`, i.e. roughly `x` in `[-87.3, 88.7)` -- outside that,
/// the exponent construction produces garbage rather than a clamped/
/// overflowed value. `expm1`, `sinh`, `cosh`, `sinh_throughput`, and
/// `cosh_throughput` inherit this poly and the same domain limit;
/// `powf`/`erf` route through `exp2_checked` and don't call this, and
/// `erfc`/`erfcx` route through [`exp_checked`], which is this poly and
/// this reduction with the input clamped first.
///
/// Scaling by `2^k` needs exp2_checked's k1/k2 split (not exp2's simpler
/// single-field trick), even though this function is otherwise
/// unchecked: `round` (unlike `floor`) can push `k` one integer past
/// where a single exponent-field construction stays valid -- e.g.
/// `x=88.37628` gives `x*log2e=127.50002`, which round pushes to `k=128`,
/// an exponent field reserved for inf/NaN (a real `exp(88.37628)=inf`
/// failure). Splitting into two representable halves sidesteps this by
/// construction.
///
/// `k`'s rounding uses the sin/cos-style magic-constant add
/// (`fma(x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC`) instead of `.round()`:
/// native `.round()` needs a multi-instruction ties-away emulation and
/// measured meaningfully slower here and in every caller. The swap to
/// the hardware's round-half-to-even differs only at exact half-integer
/// ties of `x*log2(e)`, verified zero accuracy difference on the
/// exhaustive sweep.
#[doc(alias = "expf")]
#[inline(always)]
pub fn exp(x: f32) -> f32 {
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let r = fma(-k, LN2_HI, x);
    let r = fma(-k, LN2_LO, r);
    // exp_r_poly!'s c0 and c1 are both pinned to exactly 1.0:
    // exp(r) = 1 + r + r^2*P(r) for tiny r, so a c1 off from 1.0 by even
    // ~6e-8 relative is a systematic bias right where exp(x) is most
    // commonly called (x near 0), and l0 = r + 1.0 needs no fma. c2..c5
    // are an ulp-weighted Chebyshev LP fit. Current: exp avg/max ulp
    // 0.074/3 (exhaustive).
    let p = exp_r_poly!(r);
    let (t1, t2) = exp2_field_split(k);
    p * t1 * t2
}

/// `e^x * 2^s` (backlog idea #117): a softmax/normalization building
/// block (rescaling a running exponential sum by a power of two costs
/// nothing extra here, instead of a separate multiply by `2^s` after an
/// ordinary `exp(x)`). Identical to [`exp`] except `s` folds directly
/// into the exponent-field split: `exp`'s own `k` already becomes two
/// exponent-field words via [`exp2_field_split`], and that split doesn't
/// care whether `k` came from `x`'s own reduction alone or has an extra
/// integer folded in first -- `exp2_field_split(k + s)` gives `2^(k+s)`
/// exactly as cheaply as `2^k`, so this is `exp`'s own body with one
/// `+ s as f32` added, not a separate multiply after the fact. Same
/// domain/accuracy as `exp` for `s=0`; for nonzero `s`, valid while
/// `k+s` itself stays within `exp2_field_split`'s own wide-but-not-
/// unlimited range (softmax-style rescaling needs nowhere near that).
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

/// exp(x), single-exponent-field tier (backlog ideas #23/#112): same
/// reduction and poly as [`exp`] (residual range unchanged, no refit
/// needed), but skips `exp2_field_split` entirely -- valid only while
/// `k=round(x*log2(e))` stays in `[-126,127]`, a single field's own
/// range, i.e. `x` in `[-87.68311, 88.37627]` (found by stepping one ulp
/// at a time through the real reduction to the exact boundary, same
/// method as `exp10_checked`'s clamp consolidation). Slightly narrower
/// than plain `exp`'s own unchecked `[-87.3, 88.7)` domain -- `exp`
/// needs the split specifically because `round` (unlike `floor`) can
/// push `k` one integer past a single field's valid range right at the
/// domain edge (see `exp`'s own doc comment); this tier's domain is
/// exactly small enough that `k` never reaches that edge in the first
/// place. No clamp: unchecked, like `exp` itself -- garbage outside the
/// documented domain, not saturated.
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

/// Full-range sibling of [`exp`]: `exp` already uses the k1/k2 split, so
/// the only change is clamping `x` *before* the reduction starts so `k`
/// never leaves the split's safe `[-151,128)` range, matching
/// `exp2_checked`'s early-clamp pattern. Clamp the input, not the
/// derived `k`: clamping `k` after the fact would desync it from an `r`
/// computed against the unclamped value. The bounds (`128/log2(e)`,
/// `-151/log2(e)`) are `exp2_checked`'s `k` boundary converted into
/// `x`'s units.
#[inline(always)]
pub fn exp_checked(x: f32) -> f32 {
    exp_reduce!(x.clamp(EXP_CLAMP_LO, EXP_CLAMP_HI))
}

/// A Pade approximant near 0 (where exp(x)-1 loses precision to
/// cancellation), exp(x)-1 directly elsewhere. See exp's doc comment for
/// the inherited unchecked-exp2 domain limit. The 5 Pade coefficients
/// were tuned as free parameters against f64::exp_m1 over |x| < 0.5:
/// max ulp 3, avg 0.109 on that branch.
#[doc(alias = "expm1f")]
#[inline(always)]
pub fn expm1(x: f32) -> f32 {
    let a = pade_expm1_ratio!(x, mul);
    // Deliberately a standalone copy of exp's reduction (not routed
    // through the public `exp` fn -- a shared-fn attempt regressed an
    // unrelated caller by +32%, see exp_r_poly!'s comment) ending in
    // fma(p*t1, t2, -1.0) instead of exp(x)-1.0: fuses the trailing
    // subtract into the last multiply, one rounding fewer. This branch
    // carries expm1's actual worst-case ulp (the Pade branch has
    // headroom). The poly is shared via exp_r_poly! and the field split
    // via exp2_field_split (an already-existing fn, verified via full
    // assembly diff to cost nothing here).
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let r = fma(-k, LN2_HI, x);
    let r = fma(-k, LN2_LO, r);
    let p = exp_r_poly!(r);
    let (t1, t2) = exp2_field_split(k);
    let b = fma(p * t1, t2, -1.0);
    if x.abs() < 0.5 { a } else { b }
}

/// expm1(x), single-exponent-field tier (backlog idea #201, the same
/// mechanism as [`exp_narrow`] applied here): identical Pade branch,
/// identical reduction/poly for the direct branch, but the direct
/// branch's combine drops straight to `fma(p, exp2int, -1.0)` (no `t1`
/// multiply at all, not just a narrower field) since there's only one
/// field to multiply by. Valid over the same `[-87.68311, 88.37627]`
/// domain as `exp_narrow` (identical `k=round(x*log2(e))` reduction, so
/// the same bit-level boundary applies). No clamp: unchecked, like
/// `expm1` itself.
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn expm1_narrow(x: f32) -> f32 {
    let a = pade_expm1_ratio!(x, mul);
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let r = fma(-k, LN2_HI, x);
    let r = fma(-k, LN2_LO, r);
    let p = exp_r_poly!(r);
    let exp2int = exp2int_field!(k);
    let b = fma(p, exp2int, -1.0);
    if x.abs() < 0.5 { a } else { b }
}

/// Full-range sibling of [`expm1`] (backlog idea #111): total over every
/// `f32`, at *better* throughput than `expm1` rather than the usual
/// checked-tier surcharge, because clamping the input first is what makes
/// a single exponent field legal -- so this pays for a clamp but drops the
/// whole `exp2_field_split` k1/k2 chain plus its `p * t1` multiply.
/// `expm1_checked(-inf) = -1`, `expm1_checked(inf) = inf`,
/// `expm1_checked(NaN) = NaN`, versus `expm1`'s garbage/NaN outside
/// roughly `[-87.3, 88.7)`.
///
/// Accuracy is not a tradeoff here: bit-identical to `expm1` everywhere
/// `expm1` is itself valid (verified over all 2^32 patterns), so max ulp
/// is `expm1`'s own 6, at the same worst point (`x ~ 0.9652`).
///
/// The cost split is real and goes both ways (mca): throughput
/// 1.695 -> 1.595 cyc/elem (-5.9%), latency 71.00 -> 75.06 cyc (+5.7%).
/// Fewer total ops buys the throughput, but the clamp lands at the *head*
/// of the dependency chain while the `exp2_field_split` work it replaces
/// was partly parallel to it, so the critical path gets longer even as the
/// op count falls -- the op-count-vs-path-depth distinction, in mirror
/// image. Throughput is the axis this crate optimizes, and for scale:
/// `exp` -> `exp_checked` pays +30% throughput for the same totality,
/// where this is *negative*.
///
/// The field is emitted at `k-1` (offset `382`, not [`expm1_narrow`]'s
/// `383`) with the missing factor of two folded into an exact `p + p`.
/// That indirection is the whole trick, and idea #111's own premise --
/// "the clamp caps `k` at 127 by construction" -- is what it corrects:
/// capping `k` at 127 would cap the output near `2.4e38`, so the top
/// third of an octave (`x` in `[88.376, 88.723)`, whose true results are
/// finite and representable up to `f32::MAX`) would come back short, and
/// everything above `ln(f32::MAX)` would *saturate finite* instead of
/// overflowing to `inf` -- the exact failure `edgecheck.rs` caught in the
/// rejected round-based `exp10_checked` reduction. Shifting the field
/// down one and doubling `p` instead lets `k` reach 128 while `k-1` stays
/// inside a single field's `[-126, 127]`, so the top of the range
/// overflows naturally, on its own, with no extra select.
///
/// Clamp bounds follow from that same one-field budget. Top is
/// `128/log2(e)`, [`exp_checked`]'s own bound: it sits a hair above
/// `ln(f32::MAX)`, so `x` at or past it overflows to `inf` as it should.
/// Bottom only has to be somewhere `k-1 >= -126` still holds while the
/// answer is already exactly `-1`: `e^x` is below half an ulp of 1 for
/// any `x < -17.4`, so `-86.0` (`k = -124`) clears the field's floor with
/// room to spare and every `x` below it correctly rounds to `-1.0`. The
/// Pade branch, its `|x| < 0.5` threshold, the Cody-Waite reduction and
/// `exp_r_poly!` are all shared with `expm1` unchanged.
#[inline(always)]
pub fn expm1_checked(x: f32) -> f32 {
    let a = pade_expm1_ratio!(x, mul);
    let xc = x.clamp(-86.0, 88.72283911167308);
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(xc, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let r = fma(-k, LN2_HI, xc);
    let r = fma(-k, LN2_LO, r);
    let p = exp_r_poly!(r);
    let exp2int = f32::from_bits(((k + 382_f32).to_bits() << 8) & EXPONENT_MASK);
    let b = fma(p + p, exp2int, -1.0);
    if x.abs() < 0.5 { a } else { b }
}

/// (e^x - 1)/x: the well-conditioned primitive behind financial
/// (continuously-compounded-rate) and ODE (exponential-integrator)
/// kernels, where callers otherwise write `expm1(x)/x` and hope `x`
/// never lands exactly on the removable singularity at 0. `expm1`'s Pade
/// branch is already this shape internally: `a = x * N(x)/D(x)`, so
/// `a/x = N(x)/D(x)` with the `x` factor cancelling algebraically before
/// any rounding -- no cancellation risk, not even at `x=0` itself
/// (`N(0)/D(0) = -120/-120 = 1.0` exactly, matching the true limit, so
/// no `x==0.0` select is needed at all). The direct branch (`|x|>=0.5`)
/// is `expm1`'s combine divided by `x`, a single extra rounding.
/// Duplicates `expm1`'s reduction/poly rather than routing through it
/// (same standalone-copy precedent). Inherits the unchecked-exp2 domain
/// limit: garbage outside roughly `x in [-87.3, 88.7)`.
#[inline(always)]
pub fn exp_m1_over_x(x: f32) -> f32 {
    let a = pade_expm1_ratio!(x);
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let r = fma(-k, LN2_HI, x);
    let r = fma(-k, LN2_LO, r);
    let p = exp_r_poly!(r);
    let (t1, t2) = exp2_field_split(k);
    let b = fma(p * t1, t2, -1.0) / x;
    if x.abs() < 0.5 { a } else { b }
}

/// (e^x - 1)/x, single-exponent-field tier (backlog idea #201, same
/// mechanism as [`exp_narrow`]/[`expm1_narrow`]): identical Pade branch,
/// identical reduction/poly for the direct branch, combine drops to
/// `fma(p, exp2int, -1.0) / x` (no `t1` multiply). Same
/// `[-87.68311, 88.37627]` domain as the others (identical
/// `k=round(x*log2(e))` reduction). No clamp: unchecked, like
/// `exp_m1_over_x` itself.
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn exp_m1_over_x_narrow(x: f32) -> f32 {
    let a = pade_expm1_ratio!(x);
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let r = fma(-k, LN2_HI, x);
    let r = fma(-k, LN2_LO, r);
    let p = exp_r_poly!(r);
    let exp2int = exp2int_field!(k);
    let b = fma(p, exp2int, -1.0) / x;
    if x.abs() < 0.5 { a } else { b }
}

/// 2^x - 1 (C23 `exp2m1`). Same cancellation problem as `expm1` (2^x is
/// close to 1 whenever x is close to 0, so computing 2^x first and
/// subtracting 1 loses low bits) and the same fix: `expm1`'s own Pade
/// approximant for e^y-1 is reused directly via the substitution
/// `y = x*LN_2` (2^x - 1 = e^{x ln2} - 1). This substitution would *not*
/// be safe for the direct branch: reducing through `exp2_checked(x*LOG2_E)`
/// the way exp's own doc comment warns against would reintroduce that
/// exact bug for large x, so the direct branch below instead duplicates
/// `exp2_checked`'s own k/f reduction and Q(f) poly verbatim (not routed
/// through the public `exp2_checked`, same reasoning as `expm1`'s own
/// standalone copy of `exp`'s reduction) and fuses the trailing `-1` into
/// the last multiply (`fma(p, t2, -1.0)`, one rounding instead of two).
///
/// Branch threshold is `|x| < 0.65` (backlog idea #7's "seam retunes not
/// yet done" list), not `expm1`'s own already-audited `0.5` -- checked
/// directly against the real exhaustive sweep (not just a refit), same
/// methodology as the crossover audit that shifted `asin`'s: `0.5` gives
/// avg ulp 0.0769, and every threshold tried between `0.55` and `0.75`
/// improves on that (avg bottoms out around `0.0766`-`0.0767` in
/// `0.65`-`0.7`) before `0.8`+ makes it worse again as the Pade branch's
/// own domain gets stretched. `0.65` keeps `|y| < 0.65*ln2 ≈ 0.4505`,
/// still comfortably inside the `|y|<0.5` domain `expm1` itself already
/// trusts this same approximant over -- no new coefficients, no new
/// domain risk, just using more of the already-valid range. Max ulp is
/// unaffected either way (4, exhaustive) -- the real worst point
/// (`x≈-0.3991`) sits well inside the Pade branch regardless of where
/// this threshold falls, so this is a pure avg-ulp win with no seam
/// discontinuity and no throughput/latency cost (branchless select, same
/// op count regardless of the literal). Inherits `exp2_checked`'s full
/// `[-151, 128)` clamp, so is total (never NaN/inf-producing outside its
/// true asymptotes): `exp2m1(-inf) = -1`, `exp2m1(inf) = inf`. The
/// round-domain Q(f) refit was tried and rejected here too (see IDEAS.md
/// §exp/exp2).
#[inline(always)]
#[allow(clippy::approx_constant)] // g0's constant term is a fitted minimax
// coefficient near ln(2), not ln(2) itself (bit pattern deliberately differs)
pub fn exp2m1(x: f32) -> f32 {
    let y = x * LN_2;
    let a = pade_expm1_ratio!(y, mul);

    let xs = x.clamp(-151.0, 128.0);
    let k = xs.floor();
    let f = xs - k;
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k1b = fma(xs, 0.5, ROUND_MAGIC) - (ROUND_MAGIC - 383.0);
    let k2b = (k + 766.0) - k1b;
    let t1 = f32::from_bits((k1b.to_bits() << 8) & EXPONENT_MASK);
    let t2 = f32::from_bits((k2b.to_bits() << 8) & EXPONENT_MASK);
    let q = exp2_q_poly!(f);
    let p = fma(q, t1 * f, t1);
    let b = fma(p, t2, -1.0);
    if x.abs() < 0.65 { a } else { b }
}

/// 10^x - 1 (C23 `exp10m1`), completing the C23 set next to `exp2m1`
/// above: same small-|x| Pade branch (`y = x*LN_10`, `10^x - 1 = e^{x
/// ln10} - 1`) plus [`exp10_checked`]'s own extreme-range reduction
/// (clamp, `exp10_reduction!`, `exp2_field_split`) for the direct branch,
/// with the trailing `-1` fused into the last multiply exactly like
/// `exp2m1`'s own `fma(p, t2, -1.0)`. Total over exp10_checked's full
/// domain: `exp10m1(-inf) = -1`, `exp10m1(inf) = inf`. Same clamp
/// consolidation as `exp10_checked` (idea #113, see its own doc comment):
/// the `[-45.154503, 38.53184]` bound makes a separate trailing `k`
/// clamp redundant here too (same `exp10_reduction!`, same boundary
/// math) -- verified bit-identical over the full exhaustive sweep.
///
/// Branch threshold is `|x| < 0.2`, not `exp2m1`'s `0.5`: the Pade
/// approximant (shared with `expm1`/`exp2m1` via `pade_expm1_ratio!`) is
/// only fitted/accurate for its argument `v` in `[-0.5, 0.5]` -- `expm1`
/// feeds it `v = x` directly (so its own `|x| < 0.5` threshold matches
/// exactly), and `exp2m1` feeds it `v = x*LN_2` (`|x| < 0.5` keeps `|v| <
/// 0.35`, still inside). `LN_10 ≈ 2.303` is more than 6x `LN_2`, so
/// reusing the same `0.5` threshold here would let `|v| = |x*LN_10|`
/// reach past 1.1 -- measured hundreds of ulp off near that seam. `0.2`
/// keeps `|v| < 0.47`, safely inside the fitted range.
#[inline(always)]
#[allow(clippy::approx_constant)] // g0's constant term is a fitted minimax
// coefficient near ln(2), not ln(2) itself (bit pattern deliberately differs)
pub fn exp10m1(x: f32) -> f32 {
    let y = x * std::f32::consts::LN_10;
    let a = pade_expm1_ratio!(y, mul);

    let xc = x.clamp(-45.154503, 38.53184);
    let (k, f) = exp10_reduction!(xc);
    let (t1, t2) = exp2_field_split(k);
    let q = exp2_q_poly!(f);
    let p = fma(q, t1 * f, t1);
    let b = fma(p, t2, -1.0);
    if x.abs() < 0.2 { a } else { b }
}

// exp2_checked's k1/k2 exponent-field split, factored out for
// exp/expm1/exp_pos_neg and friends. Any k1 + k2 == k works as long as
// both halves stay inside the exponent field, so k1 = round(k/2).
//
// Built on `exp2int_field!`'s magic constant (`1.5*2^23 + 127`) rather
// than the older `(k + 383) << 8 & EXPONENT_MASK` pair: folding the +127
// bias into the magic lands `k1 + 127` in the low 9 bits of a value in
// `[2^23, 2^24)`, so a single `<< 23` moves it into the exponent field
// with a zero sign bit and zero mantissa and *no mask* -- see that
// macro's own comment for why the older +383 form needed one. Two ops
// and three constants cheaper across every caller, and the `k1`/`k2`
// bookkeeping rides the same magic: `a - MAGIC` recovers `k1` exactly
// (the sum is an integer in the ulp-1 binade) and `(k - k1) + MAGIC` is
// exact for the same reason.
//
// The two forms pick opposite halves of a tie (`1.5*2^23` is even,
// `1.5*2^23 + 127` odd, and round-half-to-even resolves `k/2` against the
// magic's own parity), so k1 differs by 1 from the old form for odd k.
// That is not observable: `t1` is an exact power of two, so
// `fma(q, t1*f, t1)` is exactly `t1 * fma(q, f, 1)` with the *same*
// mantissa either way, and the single rounding into the denormal range
// still happens once, in the final `* t2`. Confirmed by a full exhaustive
// re-run of every affected function.
#[inline(always)]
fn exp2_field_split(k: f32) -> (f32, f32) {
    let a = fma(k, 0.5, EXP2INT_MAGIC);
    let k1 = a - EXP2INT_MAGIC;
    let b = (k - k1) + EXP2INT_MAGIC;
    let t1 = f32::from_bits(a.to_bits() << 23);
    let t2 = f32::from_bits(b.to_bits() << 23);
    (t1, t2)
}

// Shared exp(x)/exp(-x) for sinh/cosh: the Cody-Waite reduction only
// needs to happen once, since -x's reduction is exactly (-k, -r).
// exp's e^r poly splits into even/odd parts in r^2 (p(r) = e + r*o), so
// p(-r) = e - r*o reuses e/o at the cost of one more fma instead of a
// whole second poly; only the final exponent-field scaling (2^k vs 2^-k)
// is genuinely duplicated -- cheap bit-trick work, not fma-port
// pressure. Sharing one serial prefix is a large throughput win over two
// independent exp calls, at a small latency cost (the two calls used to
// overlap on the out-of-order CPU) -- kept, throughput is the metric
// this crate prioritizes. Coefficients are an ulp-weighted joint LP fit
// scoring p_pos/e^r and p_neg/e^-r simultaneously. Same unchecked-exp2
// domain limit as exp applies to both outputs.
#[inline(always)]
fn exp_pos_neg_half(x: f32) -> (f32, f32) {
    let (p_pos, p_neg, t1, t2, t1n, t2n) = exp_pos_neg_core!(x);
    (p_pos * t1 * t2, p_neg * t1n * t2n)
}

// exp_pos_neg_half, single-exponent-field tier (backlog idea #201, same
// mechanism as exp_narrow et al, one level further): both `+k` and `-k`
// must fit a single field's own valid range simultaneously here (unlike
// exp_narrow's one-sided k), which needs a domain a hair tighter than
// exp_narrow's own -- see sinh_narrow/cosh_narrow's own doc comment for
// the exact (symmetric) boundary. Standalone copy of exp_pos_neg_core!'s
// poly (not routed through the macro, which hardcodes the split) --
// same standalone-copy precedent as expm1/exp_checked's own reductions.
#[inline(always)]
fn exp_pos_neg_narrow_half(x: f32) -> (f32, f32) {
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let r = fma(-k, LN2_HI, x);
    let r = fma(-k, LN2_LO, r);
    // Pre-halved, exactly as in `exp_pos_neg_core!` -- see its comment for
    // why this is free and exact, and `sinh`'s for the premature-overflow
    // gap it closes.
    let c: [f32; 4] = [
        4.99993e-1 * 0.5,
        1.6667245e-1 * 0.5,
        4.188372e-2 * 0.5,
        8.300987e-3 * 0.5,
    ];
    let r2 = r * r;
    let e = fma(fma(c[2], r2, c[0]), r2, 0.5);
    let o = fma(fma(c[3], r2, c[1]), r2, 0.5);
    let p_pos = fma(r, o, e);
    let p_neg = fma(-r, o, e);
    let t = exp2int_field!(k);
    // reciprocal via bit-subtraction, same trick exp_pos_neg_core! uses
    // for t1n/t2n: exact for any power-of-two field.
    let tn = f32::from_bits(0x7F00_0000u32.wrapping_sub(t.to_bits()));
    (p_pos * t, p_neg * tn)
}

// sinh(x) ~ x * P(x^2), a two-fma odd approximation on |x| < 0.5. The
// leading coefficient is pinned to exactly 1.0 (tiny x returns x, its
// correctly-rounded sinh). c1/c2 are NOT the odd Taylor coefficients
// (1/6, 1/120): this is a degree-2 fit of sinh(x)/x over [0,0.5], one term
// shorter than the natural Taylor truncation. Dropping that term saves one
// fma on every sinh/sinh_throughput/sinh_checked call (all three evaluate
// this branch unconditionally, branchless-select) for a real ~6% throughput
// win; the cost is contained because sinh's worst case lives in the
// |x| >= 0.5 exp_pos_neg_half branch, not here -- this branch's max stays 3 ulp,
// under sinh's overall max, so the headline accuracy is unchanged (avg ulp
// rises slightly). c1/c2 are a least-squares fit (minimizes the branch's
// average ulp given that max headroom), not minimax. See IDEAS.md
// §hyperbolics for the same-degree refit that was rejected as a no-op first.
// Same role as expm1's Pade "a" branch: a cheap, cancellation-free small-x
// numerator.
#[inline(always)]
fn sinh_small(x: f32) -> f32 {
    let x2 = x * x;
    let c0 = 1.0f32;
    let c1 = 0.1666623055934906f32;
    let c2 = 0.00839646439999342f32;
    let p = fma(fma(c2, x2, c1), x2, c0);
    x * p
}

/// sinh(x) = 0.5*(exp(x) - exp(-x)) directly (via `exp_pos_neg_half`'s shared
/// reduction, see its own doc comment), except for |x| < 0.5 where exp(x)
/// and exp(-x) are both ~1 and the subtraction cancels almost all
/// precision (the same class of bug log1p/tanh had, see IDEAS.md) --
/// there, use the small-x poly form above instead, same branchless-select
/// pattern as expm1's Pade/exp split. See exp's doc comment for the
/// inherited unchecked-exp2 domain limit (only relevant on the `b` side,
/// unconditionally evaluated but only selected for |x| >= 0.5). cosh below
/// doesn't need this: it adds instead of subtracting, so it never cancels.
///
/// `exp_pos_neg_half` returns `0.5*e^(+-x)` rather than `e^(+-x)`, so
/// there is no trailing `0.5*` here -- and folding the halving into the
/// poly's own constants (free, exact) rather than paying a multiply also
/// hands plain `sinh`/`cosh` the premature-overflow fix `sinh_checked`
/// had to carry explicitly: the product never has to represent the
/// unhalved `e^x`, so `x` in roughly `[87.3, 89.4]` now returns its true
/// finite value instead of `inf` (e.g. `sinh(88.7228)` = 1.7014122e38,
/// which was `inf` before).
#[doc(alias = "sinhf")]
#[inline(always)]
pub fn sinh(x: f32) -> f32 {
    let a = sinh_small(x);
    let (ep, en) = exp_pos_neg_half(x);
    let b = ep - en;
    if x.abs() < 0.5 { a } else { b }
}

/// cosh(x) = 0.5*(exp(x) + exp(-x)) via `exp_pos_neg_half`'s shared
/// reduction (see its own doc comment) -- never cancels (adds instead of
/// subtracts), so unlike sinh needs no small-x branch. See exp's doc
/// comment for the inherited unchecked-exp2 domain limit, and `sinh`'s
/// for why the top of that range no longer overflows early.
#[doc(alias = "coshf")]
#[inline(always)]
pub fn cosh(x: f32) -> f32 {
    let (ep, en) = exp_pos_neg_half(x);
    ep + en
}

/// sinh(x), single-exponent-field tier (backlog idea #201, same
/// mechanism as `exp_narrow` et al, via [`exp_pos_neg_narrow`]): valid
/// over `[-87.68311, 87.68311]` -- symmetric and a hair tighter than
/// `exp_narrow`'s own `[-87.68311, 88.37627]`, because `exp_pos_neg_half`
/// needs *both* `k` and `-k` to fit a single field's `[-126,127]` range
/// simultaneously (found the same bit-level way: `k=126` is the last
/// safe value, since `k=127` would need `-k=-127`, one past the single
/// field's own lower edge).
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn sinh_narrow(x: f32) -> f32 {
    let a = sinh_small(x);
    let (ep, en) = exp_pos_neg_narrow_half(x);
    let b = ep - en;
    if x.abs() < 0.5 { a } else { b }
}

/// cosh(x), single-exponent-field tier -- see `sinh_narrow`'s own doc
/// comment for the domain and mechanism.
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn cosh_narrow(x: f32) -> f32 {
    let (ep, en) = exp_pos_neg_narrow_half(x);
    ep + en
}

/// `exp_pos_neg_half`, but with `x` clamped first so `k = round(x*log2e)`
/// never leaves the safe range for *both* `exp2_field_split(k)` and its
/// reciprocal-based negation, AND returning `0.5*exp(x)`/`0.5*exp(-x)`
/// directly instead of the raw pair:
///
/// 1) `exp_pos_neg_half` has no clamp, so for `|x|` large enough the
///    bit-trick exponent construction wraps around instead of saturating
///    (wrong-sign infinities, finite garbage, NaN for `+-inf` input).
///    The split (and its reciprocal negation) exactly matches
///    `2^k`/`2^-k` -- including correct saturation to `0`/`inf` -- for
///    `k` in `[-254, 254]`, with wraparound starting right outside
///    (verified by enumerating every integer `k`). Note 254, not 128:
///    128 is where `exp_checked` must stop because `exp(x)` alone
///    overflows, not where the split construction breaks. The clamp uses
///    `170.0` (`k` up to ~245.3), comfortably inside the proven ceiling
///    and far past the true `sinh`/`cosh` overflow threshold, so it only
///    ever discards inputs whose correct answer is already `+-inf`.
/// 2) The `0.5` MUST be applied inside the field split, not by the
///    caller: for `x` in roughly `[87.3, 89.4]`, `sinh(x)`/`cosh(x)` are
///    still finite (they're *half* of `exp(x)`, which overflows a bit
///    earlier), but `p_pos * t1 * t2` computes the full unscaled
///    `exp(x)` first, overflowing to `inf` before a caller could halve
///    it. Halving `t1`/`t1n` (exact for any power-of-two float down to
///    the denormal floor) *before* multiplying by `t2`/`t2n` means the
///    product only ever needs to represent `0.5*exp(x)` -- covering the
///    whole legitimately-finite window exactly, with `+-inf` correct
///    everywhere beyond it.
#[inline(always)]
fn exp_pos_neg_checked_half(x: f32) -> (f32, f32) {
    let x = x.clamp(-170.0, 170.0);
    let (p_pos, p_neg, t1, t2, t1n, t2n) = exp_pos_neg_core!(x);
    (p_pos * t1 * t2, p_neg * t1n * t2n)
}

/// Full-range sibling of [`sinh`] -- same construction, just built on
/// `exp_pos_neg_checked_half` instead of the unclamped
/// `exp_pos_neg_half` (both return the `0.5*exp(+-x)` halves, so neither
/// needs a separate `0.5*` multiply). See that function's own doc
/// comment for the correctness gaps this closes: `sinh`'s wrong-sign/NaN
/// behavior for large `|x|` (the premature-overflow gap that entry also
/// mentions is now closed for plain `sinh` too, see its own doc).
#[inline(always)]
pub fn sinh_checked(x: f32) -> f32 {
    let a = sinh_small(x);
    let (ep, en) = exp_pos_neg_checked_half(x);
    let b = ep - en;
    if x.abs() < 0.5 { a } else { b }
}

/// Full-range sibling of [`cosh`] -- same construction, just built on
/// `exp_pos_neg_checked_half` instead of the unchecked `exp_pos_neg`.
/// See `exp_pos_neg_checked_half`'s own doc comment for the
/// correctness gaps this closes.
#[inline(always)]
pub fn cosh_checked(x: f32) -> f32 {
    let (ep, en) = exp_pos_neg_checked_half(x);
    ep + en
}

// coshm1's own poly: `Q(u) = (cosh(sqrt(u)) - 1 - u/2) / u^2` on
// `u = x^2` in `[0, 4]`, i.e. the series with *both* leading terms peeled
// off, not just the constant. Minimax by LP (HiGHS) weighted by
// `u^2/(cosh(sqrt(u))-1)`, the factor that turns an absolute error here
// into the result's relative error. Degree 3, residual 0.177 ulp: degree 4
// halves that (0.088) and measures *identically* end to end, because what
// actually binds is the final `fma`'s own rounding plus `u = x*x`'s, not
// the fit -- so the extra term was dropped again.
const COSHM1_Q_COEFFS: [f32; 4] = [0.04166663, 0.0013889787, 2.4732362e-5, 2.963389e-7];

/// cosh(x) - 1 (backlog idea #144), the catenary/relativity primitive
/// where naive `cosh(x)-1` cancels badly for small `x` (`cosh(x)` is
/// `~1` there, the same class of cancellation `expm1`/`log1p` exist to
/// avoid). Two branches, each avoiding the *other*'s error mechanism:
///
/// - `|x| < 2`: `x^2/2 + (x^2)^2 * Q(x^2)`, with `x*(0.5*x)` -- a single
///   correctly-rounded multiply, exact all the way into the denormal
///   floor -- carrying the leading term. Peeling `x^2/2` out of the
///   polynomial is what makes this branch cheap *and* accurate: the
///   whole poly reaches the answer scaled by `u^2`, worth 28% of it at
///   `x = 2` and vanishing as `x -> 0`, so its own roundings are
///   demoted everywhere the cancellation would have mattered.
/// - `|x| >= 2`: `cosh_checked(x) - 1.0` directly. Subtracting 1 from a
///   value `>= 1` is *exact* in binary floating point (the result's
///   exponent drops by at most one, so `1` is always on the coarser
///   grid), so this branch is exactly as accurate as `cosh_checked`
///   itself, up to the conditioning factor `cosh/(cosh-1)` -- 1.36 at
///   `x = 2` and 1.0007 by `x = 4`, which is why the handover sits at
///   2 rather than lower. `cosh_checked` rather than plain `cosh` for
///   the same full-range reason `norm_pdf` uses `exp_checked`: `coshm1`
///   grows as fast as `cosh`, so a real caller reaches the unchecked
///   tier's domain edge easily.
///
/// The half-angle identity `2*sinh(x/2)^2` this replaces had no
/// cancellation either, but squaring doubles the relative error of
/// `sinh_checked` by construction -- verified rather than assumed, by
/// scoring `2 * relerr(sinh_checked(x/2))` per octave against the real
/// end-to-end error and finding it matched to three digits at every one.
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

/// tanh(x) = (e^2x - 1) / (e^2x + 1) = expm1(2x) / (expm1(2x) + 2), same
/// formula as before, reusing expm1's already-correct small-x handling
/// (its own Pade branch below |x|<0.5) instead of computing exp2(2x) and
/// cancelling `1.0 - (~1.0)` directly, which lost essentially all
/// precision for small x (a fuzz sweep found this the same
/// 300-million-ulp-average class of bug as log1p's, before that fix --
/// see IDEAS.md).
///
/// `2*x` is clamped to `[-87.0, 88.0]` before the reduction: the raw
/// `2*x` inherited exp's unchecked-domain garbage for |x| > ~44
/// (`tanh(50)`/`tanh(f32::MAX)` returned NaN where std saturates to 1).
/// An abs/mulsign restructuring was tried and rejected -- it still
/// needed a clamp to be safe, so it did no work the clamp alone doesn't
/// (see graveyard.md §hyperbolics). `f32::clamp` returns NaN unchanged
/// (unlike `.max`/`.min`), so NaN propagation isn't broken, and the
/// clamp is lossless even for in-range-but-large x: any x past the
/// boundary already has a true tanh value of exactly `+-1.0f32` many
/// orders of magnitude before reaching it, so clamping produces the
/// bit-identical correctly-rounded answer. The clamp's real perf cost
/// was accepted per this crate's usual "pay to fix wrong/NaN for
/// legitimate finite input" precedent.
#[doc(alias = "tanhf")]
#[inline(always)]
pub fn tanh(x: f32) -> f32 {
    // Standalone copy of expm1 (not a call through the public `expm1`
    // fn, same shared-helper scheduling risk as everywhere else) but
    // with a single exponent-field construction instead of exp's k1/k2
    // split: the clamp bound guarantees `k = round(x*2*log2e)` stays in
    // [-126, 127], comfortably short of the k=128 edge case the split
    // exists for. Same fma(p, exp2int, -1.0) tail fusion as expm1.
    //
    // Unlike the earlier form, `2*x` is never materialized as its own
    // value: every downstream constant is pre-scaled by the matching
    // power of two instead (2*LOG2_E, halved LN2_HI/LN2_LO, the poly's
    // c[n] coefficients each *2^(n+1), the Pade numerator/denominator
    // each /8 to match its two extra Horner multiplies). This is exact,
    // not an approximation: correctly-rounded arithmetic (every `fma`/
    // `*` step here) commutes exactly with power-of-2 scaling of all its
    // inputs, so each intermediate is bit-for-bit the old value at half
    // (or a smaller power-of-two fraction of) its former scale, all the
    // way through to the final `a`/`b` -- verified bit-identical
    // exhaustively, see IDEAS.md idea #119. Deletes the standalone
    // `2.0 * x` multiply from the critical path.
    let xc = x.clamp(-43.5, 44.0);

    const PADE_N_A: f32 = -1.9999927; // unscaled: shared with expm1/sinh_small's own copy
    const PADE_N_B: f32 = -120.0 / 4.0;
    const PADE_D_C1: f32 = 12.000030 / 2.0;
    const PADE_D_C2: f32 = 59.999996 / 4.0;
    const PADE_D_C3: f32 = -120.0 / 8.0;
    let a = xc * fma(PADE_N_A, xc * xc, PADE_N_B)
        / fma(xc, fma(xc, xc - PADE_D_C1, PADE_D_C2), PADE_D_C3);

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

    let e = if xc.abs() < 0.25 { a } else { b };
    e / (e + 2.0)
}

/// tanh'(x) = 1 - tanh(x)^2 (backlog idea #150), the gradient ML
/// backprop through a tanh activation needs. The obvious `1.0 -
/// tanh(x)*tanh(x)` composite was tried first and rejected as a real bug
/// (found by exhaustive fuzzing, not inspection): `tanh(x)` correctly
/// rounds to exactly `1.0f32` once `|x|` exceeds roughly `8.66` (`1 -
/// tanh(x)^2`'s true value there, `~2.5e-8`, is smaller than a whole ulp
/// of `tanh`'s own result near `1.0`, `~1.19e-7`), so the subtraction
/// cancels to exactly `0.0` for a wide band of `x` where the true
/// gradient, while small, is still many orders of magnitude away from
/// underflow -- max ulp in the billions, not a benign "near a true zero"
/// artifact.
///
/// Fixed with the standard stable form instead: `1 - tanh(x)^2 =
/// 4q/(1+q)^2` where `q = exp(-2|x|)` (derived from `tanh(x) =
/// (1-q)/(1+q)` for `x>=0`; squaring removes the sign, so using `|x|`
/// is exact, not an approximation, and this needs no clamp at all since
/// `-2|x| <= 0` always keeps `exp_checked`'s argument in its own
/// well-behaved, never-overflowing range). No cancellation anywhere: `q`
/// itself is the actual small quantity being tracked, never subtracted
/// from a same-magnitude value the way `1 - tanh(x)^2` is.
#[inline(always)]
pub fn tanh_grad(x: f32) -> f32 {
    let q = exp_checked(-2.0 * x.abs());
    4.0 * q / ((1.0 + q) * (1.0 + q))
}

/// logistic sigmoid, `1/(1+exp(-x))`, computed directly. The
/// algebraically-exact identity `sigmoid(x) = 0.5 + 0.5*tanh(x/2)` was
/// tried and rejected as a real bug: around `x = -17.3`, `tanh(x/2)`
/// correctly rounds to exactly `-1.0f32`, so `0.5 + 0.5*(-1.0)` gives
/// exactly `0.0` even though the true value (`~2.98e-8`) is nowhere near
/// f32's underflow threshold -- tanh's correct saturation discards
/// exactly the residual precision the identity needs (same class of bug
/// as `atanh`'s rejected single-log1p fusion). The direct form has no
/// such cancellation anywhere and gracefully saturates to exactly
/// `0.0`/`1.0` over the whole domain, never inf/nan for finite input.
#[inline(always)]
pub fn sigmoid(x: f32) -> f32 {
    // Standalone copy of exp's reduction (not routed through the public
    // `exp` fn, same pattern as expm1); the poly itself is shared via
    // `exp_r_poly!`. Single exponent-field construction, NOT
    // exp2_checked's k1/k2 split (tried: it works but doubles throughput
    // cost, an unjustified price here).
    //
    // The negative-side bound is `88.722839111673` (`128/log2(e)`, the
    // exact point where `k=round(y*log2e)` for `y=-x` reaches `128`),
    // deliberately wider than tanh's analogous `-87`: sigmoid's
    // asymptote is at `0`, which needs `e` to range all the way to where
    // it correctly *overflows* to `+inf` (so `1/(1+inf)=0` exactly) --
    // a narrower clamp froze the negative tail at a fixed wrong constant
    // for every `x` below it (sigmoid(-89..-inf) all ~6.05e-39). At
    // `k=128` the single-field trick naturally lands `exp2int` on the
    // `+inf` bit pattern, and `k` stays exactly `128` (never wrapping to
    // `129`) for every `y` up to ~89.05, comfortable margin around the
    // bound. The `x -> +inf` side only needs the already-representable
    // asymptote `1`, so `87.0` is fine there (same reasoning as tanh's
    // clamp). Known accepted gap: for `x` in roughly `(-104.7,-88.7)`
    // the true answer is a nonzero denormal but this returns exactly `0`
    // -- slightly early saturation on a sliver of denormal-scale
    // outputs, per the crate's "near a true zero, ulp isn't meaningful"
    // precedent.
    //
    // The negations are folded away rather than computed: clamp `x`
    // directly instead of `-x` (swap and negate the literal bounds),
    // fold the sign into `LOG2_E` for `k` (sign commutes exactly through
    // a multiply), carry `k` positive through both LN2_HI/LN2_LO fmas,
    // and apply one final negation to get `r`. Pure reassociation, bit-
    // exact (FMA rounds symmetrically under negation of all its inputs:
    // `fma(-a,b,-c) == -fma(a,b,c)` always); nets one fewer runtime
    // negation.
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

/// Hard-clamped piecewise-linear + one cubic correction (backlog idea
/// #191): an approx-tier `sigmoid` for ML-inference latency, no `exp`
/// call at all. `sigmoid(x) ~ 0.5 + x*(a + b*x^2)` for `|x| <=
/// 3.288051` (`a`/`b` a real minimax fit, scipy, not a hand-derived
/// Taylor truncation), clamped past that -- both the *input* (so the
/// cubic term, unbounded outside its fit domain, never gets a chance to
/// swing back the wrong way for large `|x|`; tried clamping only the
/// *output* first and it does exactly that, silently returning ~1.0 for
/// very negative `x` instead of ~0.0, found by checking the full range
/// rather than assuming the fit domain was wide enough) and the
/// *output* (belt-and-suspenders bound to `[0,1]`, though the input
/// clamp alone already keeps the poly's own range inside that here).
/// Max absolute error ~0.023 (minimax over the whole domain, jointly
/// fit with the clamp threshold itself, not just the polynomial) --
/// deliberately outside this crate's normal 0.5/2 ulp budget, the same
/// tier `exp2_approx`/`rsqrt_approx` occupy.
#[inline(always)]
pub fn sigmoid_fast(x: f32) -> f32 {
    let xc = x.clamp(-3.288051, 3.288051);
    let poly = fma(-0.006715598, xc * xc, 0.2178126);
    fma(poly, xc, 0.5).clamp(0.0, 1.0)
}

/// sigmoid'(x) = sigmoid(x) * (1 - sigmoid(x)) (backlog idea #150), the
/// gradient ML backprop through a sigmoid activation needs. Same class of
/// real cancellation bug as `tanh_grad`'s own doc comment, found the same
/// way: `sigmoid(x)` correctly rounds to exactly `1.0f32` once `x`
/// exceeds roughly `16.6`, so `1.0 - sigmoid(x)` (the true gradient's own
/// dominant factor there) cancels to exactly `0.0` -- max ulp in the
/// hundreds of millions, not underflow.
///
/// Fixed with `sigmoid(x)*(1-sigmoid(x)) = e/(1+e)^2`, `e = exp(-x)` --
/// but naively using `e = exp_checked(-x)` directly traded that bug for a
/// different one: for very negative `x`, `e` grows huge (not just up to
/// `sigmoid`'s own saturation point but arbitrarily large), and `(1+e)^2`
/// overflows to `inf` while `e` is still finite, silently giving `0.0`
/// for inputs whose true answer (e.g. `x=-44.36`, true value `~5.4e-20`)
/// is nowhere near underflow -- a real loss, not `atanh`'s "near a true
/// zero" class of benign artifact. Fixed like `tanh_grad`: this is an
/// even function (`sigmoid(-x)*(1-sigmoid(-x)) == sigmoid(x)*(1-
/// sigmoid(x))`, swap `s -> 1-s`), so folding onto `|x|` keeps `e` always
/// in `(0,1]`, never large enough for the square to overflow, no branch
/// needed at all.
#[inline(always)]
pub fn sigmoid_grad(x: f32) -> f32 {
    let e = exp_checked(-x.abs());
    e / ((1.0 + e) * (1.0 + e))
}

/// `log1p(e)` specialized for callers whose own domain already
/// guarantees `e` in `(0, 1]` (backlog idea #120: softplus/logaddexp's
/// own `e = exp(-something.min(87.0))`, `something >= 0`, so `e` never
/// leaves that range) -- `u = 1+e` then always lands in `(1, 2]`, safely
/// away from every special case `log1p`'s own wrapper exists for (never
/// zero, negative, denormal, inf, or nan), so this skips straight to
/// `ln_normal` + the Sterbenz correction, no wrapper, no
/// `corr.is_finite()` guard, no `x==0.0` select. Same domain-bypass
/// mechanism as `asinh`/`acosh`'s own `ln_normal` calls, distinct from
/// the already-rejected "log1p small-|x| branch" (which added a poly to
/// *every* general `log1p` call regardless of caller) -- this is a
/// separate callee only reachable from callers whose domain already
/// proves the skipped checks unreachable.
///
/// That bounded domain also makes the whole `ln` machinery unnecessary
/// (backlog idea #38): with `e` confined to `(0, 1]` there is no
/// reduction to do, so a single fitted polynomial covers the range
/// directly and the reduction, the exponent-field bit extraction, the
/// `k*LN2_HI + k*LN2_LO` recombine *and* the Sterbenz correction's
/// division all disappear. Written as `e + e^2*Q(e)` rather than
/// `e*P(e)`: the leading `e` is then exact and only the (much smaller)
/// correction carries the polynomial's rounding, which is also what
/// makes tiny `e` come back bit-exact -- `e*e` flushes to zero below
/// ~1e-19 and the result is literally `e`, matching the old form's own
/// `1.0 + e == 1.0` path.
///
/// Degree 9 in `Q`, evaluated Estrin with the last two coefficients
/// folded in at the `e^4` level (`ln_normal`'s own trick, so `e^8` never
/// has to be formed): idealized max relative error 0.078 ulp-equivalent
/// with the coefficients already rounded to f32.
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

/// softplus(x) = ln(1+e^x), the smooth approximation to `max(x,0)` ML
/// frameworks call `log1pexp`/`softplus`. Naive `(1.0+exp(x)).ln()`
/// overflows for large `x` (`exp(x)` alone does) and loses precision for
/// very negative `x` (`1.0+tiny` rounds to exactly `1.0`, the same
/// cancellation `log1p` exists to avoid) -- the standard numerically
/// stable form instead: `softplus(x) = max(x,0) + log1p(exp(-|x|))`
/// (verified algebraically: for `x>=0`, `x + ln(1+e^-x) =
/// ln(e^x) + ln(1+e^-x) = ln(e^x(1+e^-x)) = ln(e^x+1)`; for `x<0`,
/// `max(x,0)=0` and this reduces to `ln(1+e^x)` directly).
///
/// Two landmines in the obvious formulation, both avoided here:
/// 1. For `x` below the cutoff, `max(x,0)=0`, so the *entire* result is
///    the correction term -- clamping `exp`'s argument there would
///    replace the true (much smaller) correction with a fixed,
///    comparatively huge stand-in, and for negative `x` nothing hides
///    that error. So the correction is *selected* to exactly `0.0` once
///    `|x|` is far enough out that the true value is negligible at f32
///    precision, rather than feeding `exp` a clamped-but-wrong argument.
/// 2. `f32::max`/`min` follow IEEE `maxNum`/`minNum` semantics and
///    *discard* NaN rather than propagate it -- `x.max(0.0)` and the
///    exponent-clamping `min` would silently turn `softplus(NaN)` into
///    finite garbage. Guarded with an explicit trailing `is_nan` check.
#[inline(always)]
pub fn softplus(x: f32) -> f32 {
    let ax = x.abs();
    // `exp_narrow`, not `exp`: the `min(87.0)` above is already the guard
    // its single-exponent-field domain (`x` in `[-87.68311, 88.37627]`)
    // asks for, so the k1/k2 split is dead weight here -- `-ax.min(87.0)`
    // lands in `[-87, 0]` for every input, NaN included (`min` follows
    // IEEE `minNum` and returns `87.0`, and the trailing `is_nan` below
    // restores the NaN). Bit-identical, since `t1 * t2` and the single
    // field are the same exact power of two over this k range.
    let e = exp_narrow(-ax.min(87.0));
    let corr = if ax > 87.0 { 0.0 } else { log1p_unit(e) };
    let normal = x.max(0.0) + corr;
    if x.is_nan() { f32::NAN } else { normal }
}

/// logsigmoid(x) = ln(sigmoid(x)) = -softplus(-x) (backlog idea #118), the
/// numerically stable log-likelihood ML frameworks pair with `sigmoid`
/// (`ln(1/(1+e^-x))`, same cancellation trap as `softplus` itself for
/// large negative `x` -- reuses that fix directly instead of composing
/// `sigmoid(x).ln()`, which would underflow to `ln(0.0) = -inf` far
/// earlier than the true answer justifies).
#[inline(always)]
pub fn logsigmoid(x: f32) -> f32 {
    -softplus(-x)
}

/// logaddexp(a,b) = ln(e^a+e^b), the numerically stable "log of a sum of
/// exponentials" ML/statistics primitive (softmax/log-sum-exp's binary
/// building block -- in fact `softplus(x) == logaddexp(x, 0.0)`).
/// `max(a,b) + log1p(exp(-|a-b|))`, same derivation shape as `softplus`
/// (factor out `e^max(a,b)`, same algebra). Reuses `softplus`'s own two
/// fixes directly: the correction-term cutoff at `87.0` (here on
/// `|a-b|`, not `|x|`) instead of a clamped `exp` argument, and an
/// explicit trailing NaN guard (`a.max(b)` alone would silently discard
/// a NaN `a`/`b` the same way `softplus`'s own `x.max(0.0)` did).
///
/// Unlike `softplus`, `m` and `corr` here can be comparable-magnitude,
/// opposite-signed values that partially cancel (e.g. near
/// `logaddexp(-7e-5, -9.57)`): max ulp in fuzzing reaches the thousands
/// there while avg stays ~0.15 -- a real but narrow, accepted
/// cancellation cost, comparable in kind to `erfc`'s own outlier.
#[inline(always)]
pub fn logaddexp(a: f32, b: f32) -> f32 {
    let m = a.max(b);
    let d = (a - b).abs();
    // `exp_narrow` for the same reason as `softplus`'s -- see its comment.
    let e = exp_narrow(-d.min(87.0));
    let corr = if d > 87.0 { 0.0 } else { log1p_unit(e) };
    let normal = m + corr;
    if a.is_nan() || b.is_nan() { f32::NAN } else { normal }
}

/// `1/sqrt(2)` as a double-`f32` pair, same shape as
/// `FRAC_1_PI`/`RPI_LO` above: `RSQRT2_HI` is the nearest `f32` and
/// `RSQRT2_LO` the next
/// 24 bits of the remainder, together naming the constant to a relative
/// `2^-49`. Used by [`gelu`] to hand `erfc` an argument that is exact to
/// far more than `f32`, which its tail needs -- see [`gelu`]'s comment.
const RSQRT2_HI: f32 = std::f32::consts::FRAC_1_SQRT_2;
const RSQRT2_LO: f32 = 1.2101617485882343e-8;

/// GELU (Gaussian Error Linear Unit), the exact/erf-based form (as
/// opposed to the tanh approximation): `x * Phi(x)` where `Phi` is the
/// standard normal CDF. The de facto default activation in transformer
/// architectures (backlog idea #70).
///
/// `Phi(x) = 0.5*erfc(-x/sqrt2)`, *not* the algebraically-equivalent
/// `0.5*(1+erf(x/sqrt2))`: for negative `x`, `erf(x/sqrt2)` approaches
/// `-1`, so `1+erf(x/sqrt2)` cancels toward `0` and inherits `erf`'s own
/// small *absolute* error as a huge *relative* one (measured: >1e8 ulp
/// once `x` is a few units negative). `erfc(-x/sqrt2)` computes that
/// same near-zero tail value directly (its whole reason to exist, see
/// `erfc`'s own doc comment) instead of via subtractive cancellation, so
/// this form has no such blowup. `x = -inf` is the one input the
/// formula alone still mishandles: `erfc` saturates cleanly to `0.0`
/// there, but `0.0` times the *literal* `x = -inf` is an indeterminate
/// `0*inf`, even though the true limit (`x*Phi(x)` as `x -> -inf`) is
/// `0` -- same shape as `sqrt1pm1`'s own `x == inf` override below.
///
/// The argument `-x/sqrt2` is fed to `erfc` as a double-`f32` pair
/// rather than a single rounded `f32`, because `erfc` is violently
/// ill-conditioned in its argument out in the tail: the relative
/// sensitivity `|z * dln(erfc)/dz|` grows as `2z^2`, so at `x = -13`
/// (`z = 9.2`) the half-ulp already present in a rounded `z` reappears
/// as ~170 half-ulps of the result. `RSQRT2_HI + RSQRT2_LO` names
/// `1/sqrt2` to 2^-49, `dz` recovers the rest of the product exactly via
/// `fma`, and the first-order term `erfc(z+dz) = erfc(z)*(1 - 2*z*dz)`
/// (the `2z` asymptotic form of that same log-derivative) puts it back.
/// `np = max(-x, 0)` clamps the correction off for `x >= 0`, where the
/// true log-derivative decays like `exp(-z^2)` instead and `2z` would be
/// badly wrong; the clamp also keeps `dz` finite at `x = +inf`, and the
/// correction is applied to `erfc`'s own (bounded) result rather than to
/// `x*Phi(x)` so that `+inf` stays `+inf` instead of hitting `inf*0`.
#[inline(always)]
pub fn gelu(x: f32) -> f32 {
    let nx = -x;
    let np = nx.max(0.0);
    let zp = np * RSQRT2_HI;
    let dz = fma(np, RSQRT2_LO, fma(np, RSQRT2_HI, -zp));
    let e = erfc(nx * RSQRT2_HI);
    let normal = x * 0.5 * fma(-e, (zp + zp) * dz, e);
    if x == f32::NEG_INFINITY { 0.0 } else { normal }
}

/// SiLU / Swish: `x * sigmoid(x)` (backlog idea #70). Same `x = -inf`
/// `0*inf` indeterminate-form fix as `gelu` above (`sigmoid(-inf) = 0`
/// exactly, but `-inf * 0.0` alone is `NaN`, not the true limit `0`).
#[inline(always)]
pub fn silu(x: f32) -> f32 {
    let normal = x * sigmoid(x);
    if x == f32::NEG_INFINITY { 0.0 } else { normal }
}

/// Softsign: `x / (1 + |x|)` (backlog idea #70). `x = +-inf` is the one
/// input the formula alone mishandles (`inf/(inf+1)` is an
/// indeterminate `inf/inf`, even though the true limit is `+-1`) --
/// same shape as `sqrt1pm1`'s own `x == inf` override below.
#[inline(always)]
pub fn softsign(x: f32) -> f32 {
    let normal = x / (1.0 + x.abs());
    if x.is_infinite() { x.signum() } else { normal }
}

/// `sqrt(1+x) - 1`, the rationalized form that avoids the catastrophic
/// cancellation a caller writing `(1.0+x).sqrt() - 1.0` directly would
/// hit for small `|x|` (`1.0+x` rounds to exactly `1.0` well before `x`
/// itself underflows, so the naive form silently returns exactly `0`):
/// `sqrt(1+x) - 1 == x / (sqrt(1+x) + 1)` algebraically, and the
/// denominator's `+1` never cancels (`sqrt(1+x) >= 0` for any in-domain
/// `x`), so this is accurate across the *entire* domain with no branch
/// needed -- the same rationalization [`asinh`]/[`acosh`]/[`asin`] each
/// already re-derive inline for their own `sqrt(...) - 1`-shaped terms,
/// exposed here directly as its own function (backlog idea #132).
/// Verified against a high-precision (f64, itself rationalized the same
/// way to avoid the identical cancellation trap one level up) reference
/// over a 74M-sample fuzz: avg ulp 0.21, max ulp 2. Domain `x >= -1.0`
/// (else `NaN`, matching `1+x < 0`'s real domain error); `x == +inf` is
/// the one input the formula alone mishandles (`inf/(inf+1)` is an
/// indeterminate `inf/inf`), corrected with a trailing override.
#[inline(always)]
pub fn sqrt1pm1(x: f32) -> f32 {
    let normal = x / ((1.0 + x).sqrt() + 1.0);
    if x.is_infinite() { x } else { normal }
}

/// x^(3/2) (backlog idea #133): `x * sqrt(x)`, two correctly-rounded
/// hardware ops, always more accurate and faster than routing through
/// `powf`. Domain `x >= 0` (matching the real-valued convention);
/// `sqrt` alone already supplies every special case for free: `x=0`
/// gives `0*0=0`, `x=inf` gives `inf*inf=inf`, `x<0`/`NaN` give `NaN`
/// (IEEE `sqrt` of a negative number).
#[inline(always)]
pub fn pow_3_2(x: f32) -> f32 {
    x * x.sqrt()
}

/// Core of `x^(2/3)` for `a` positive and normal: `cbrt_normal`'s own
/// bit-trick seed `s` and residual `r = (s^3-a)/a`, but fitting
/// `(1+r)^(-2/3)` instead of `(1+r)^(-1/3)` -- `s^2 * (1+r)^(-2/3) ==
/// a^(2/3)` identically, exactly as `s * (1+r)^(-1/3) == a^(1/3)`, so the
/// two-thirds power is a *direct* fit rather than a cube root squared.
/// `scale`/`scale3` carry the caller's denormal rescale (`scale3` is
/// `scale/3`, see below); folding them in here rather than multiplying
/// the return value keeps them off the tail of the dependency chain,
/// which measures a real latency win. That is safe only because `scale`
/// is an exact power of two AND the poly stays in Horner form -- with the
/// Estrin schedule `cbrt_normal` uses, LLVM's vectorizer takes the scale
/// parameter as licence to duplicate this whole function for the tiny and
/// normal branches instead of computing once and blending, doubling every
/// fma/mul/div (see `cbrt_normal`'s own comment on the same trap).
///
/// The leading term is `s^2` rather than `s`, and `s*s` -- unlike the bit
/// pattern `s` -- is not exact, so its rounding `e2` reaches the result at
/// full weight and has to be added back. `e2` lands twice over:
///
/// - directly, since `a^(2/3) = (s2+e2)*(1+r)^(-2/3)` and the poly only
///   covers the `s2` part;
/// - through `r`, because `d = fma(s2, s, -a)` computes `s2*s - a`, short
///   of the intended `s^3 - a` by `e2*s ~ 2^-24 * a`. That shortfall
///   reaches the result multiplied by `-2/3` (against `cbrt`'s `-1/3`),
///   and `s2*s/a ~ 1`, so it is `-(2/3)*e2` to well within its own last
///   bit -- no second `fma` on the critical path needed to correct `d`.
///
/// The two together are `e2 - (2/3)*e2 = e2/3`, which is why the tail
/// addend is scaled by a third.
///
/// The poly is degree 4, not `cbrt_normal`'s degree 3: `(1+r)^(-2/3)`'s
/// series coefficients are ~3x `(1+r)^(-1/3)`'s, so degree 3 over the
/// seed's `r` range (`[-0.0999, 0.0894]`, exhaustive over all three of
/// the seed's exponent-mod-3 alignment classes) fits to only ~2 ulp where
/// cbrt's own degree 3 reaches ~0.5. Degree 4 fits to 0.13 ulp, leaving
/// the final rounding dominant.
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

/// x^(2/3) (backlog idea #133), fitted directly rather than as
/// `cbrt(x)^2`: squaring a cube root doubles its relative error before
/// the square's own rounding is even applied, and `cbrt`'s own budget is
/// ~0.28 avg ulp. See [`pow_2_3_normal`] for the kernel. `x^(2/3)` is
/// even, so this reads `|x|` and is defined for every real `x` (the
/// standard real-valued extension of a rational power via the odd root),
/// unlike `pow_3_2`. Denormals get `cbrt`'s rescale, with `(2^24)^(2/3) =
/// 2^16` undone on the way out; `+-0`, `+-inf` and `NaN` take the same
/// trailing select (`|x| + |x|` gives `+0`, `+inf` and `NaN` respectively,
/// all three the correct even-power answers).
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
    if ax == 0 || ax >= EXPONENT_MASK { a + a } else { r }
}

/// Hermite smoothstep (backlog idea #147), GLSL's `smoothstep(edge0,
/// edge1, x)`: `t^2*(3-2t)` on the normalized, clamped `t =
/// (x-edge0)/(edge1-edge0)` -- exact endpoints (`0` at `t=0`, `1` at
/// `t=1`) and zero slope at both, by construction of the Hermite basis,
/// not by any special-casing here. `edge0==edge1` is the one input this
/// leaves undefined (a `0/0` or `x/0`), matching every other
/// implementation of this GLSL primitive.
#[inline(always)]
pub fn smoothstep(edge0: f32, edge1: f32, x: f32) -> f32 {
    let t = ((x - edge0) / (edge1 - edge0)).clamp(0.0, 1.0);
    t * t * fma(-2.0, t, 3.0)
}

/// Perlin's "smootherstep" (backlog idea #147): `t^3*(6t^2-15t+10)`, the
/// degree-5 Hermite variant with zero *second* derivative at both
/// endpoints too (not just zero first derivative like [`smoothstep`]),
/// removing the curvature discontinuity animators call "banding" at the
/// seams. Same `edge0`/`edge1` normalization and domain caveat as
/// `smoothstep`.
#[inline(always)]
pub fn smootherstep(edge0: f32, edge1: f32, x: f32) -> f32 {
    let t = ((x - edge0) / (edge1 - edge0)).clamp(0.0, 1.0);
    let p = fma(t, fma(t, 6.0, -15.0), 10.0);
    t * t * t * p
}

/// asinh(x) = ln(x + sqrt(x^2+1)), with two fixes over the naive form:
///
/// 1. Small-x cliff: for |x| below ~6e-8, `x*x` is already too small to
///    survive being summed into `x^2+1` (it rounds to exactly 1.0 before
///    sqrt even runs), so `x + sqrt(x^2+1)` collapses to exactly 1.0 and
///    `ln(...)` returns exactly 0 instead of the correct tiny nonzero
///    value. Fixed the standard way: `sqrt(x^2+1) - 1 = x^2 / (sqrt(x^2+1)
///    + 1)` (rationalized, no cancellation -- `x^2` is computed as its own
///    multiply here, independent of the lossy `x^2+1` sum, so it keeps
///    full precision), then `asinh(x) = log1p(x + (sqrt(x^2+1) - 1))`.
/// 2. Large-negative-x cancellation (the actual worst case): for x very
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
    // runs for ax >= 2048, where the first dropped term is 1/(8*ax^4) <=
    // 7e-15 relative -- seven orders of magnitude under f32's own 6e-8, so
    // the two-term form is exact here. Costs one division and one fma
    // instead of a division, an add, a sqrt and a multiply.
    let inv_ax = 1.0 / ax;
    let rescaled_sq = fma(0.5, inv_ax, ax);
    let sq = if small { direct_sq } else { rescaled_sq };
    let sm1 = if small { ax2 / (sq + 1.0) } else { sq - 1.0 };
    let d = ax + sm1;
    // `log1p_finite(d)` and the overflow fallback `ln(ax) + LN_2` each
    // pay their own full `ln`-family poly + wrapper, unconditionally
    // (branchless), for every call -- but `ln(2*ax) = ln(ax) + ln(2)` is
    // exactly `ln_normal`'s own `koff` hook (it adds directly into the
    // pre-combine exponent field `k`, not as a post-hoc add onto an
    // already-rounded `ln(ax)`), so both branches reduce to one shared
    // `ln_normal` call on a selected (argument, koff) pair. `u = 1+d` is
    // always `>= 1` (finite-d branch, `d = ax+sm1 >= 0`) and `ax` is
    // always a genuine positive value here too, so neither ever needs
    // `log_family_wrapper!`'s zero/negative/denormal handling -- only
    // its inf/nan handling, which the trailing overrides below restore
    // (the raw `_normal` core doesn't propagate either, see its own doc
    // comment: `ax` is `NaN`/`+inf` exactly when `x` is, since `d`'s own
    // non-finiteness routes here).
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

/// ln(x + sqrt(x^2-1)), domain x >= 1 (NaN elsewhere). Four fixes over
/// the naive `x*x - 1.0` form:
///
/// 1. Sign loss: squaring erases x's sign, so sqrt(x^2-1) is the same
///    magnitude for +x and -x. Once |x| is large enough that ulp(x^2)
///    exceeds 1 (roughly |x| > 4096, far short of actual overflow) the
///    "-1" term vanishes entirely and sqrt(x^2-1) rounds to exactly |x|,
///    so `x + sqrt(x^2-1)` collapses to ~0 for negative x instead of
///    staying reliably negative -- ln of that silently returns finite
///    garbage (or +inf, once x^2 overflows) instead of the correct NaN,
///    for roughly the whole range x < -4096. Fixed with an explicit
///    domain select (cheap next to the sqrt+ln chain).
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
///    correction on top. Cut max ulp at the boundary from 1522 to 4.
#[doc(alias = "acoshf")]
#[inline(always)]
pub fn acosh(x: f32) -> f32 {
    // NOT the same "shared x2" opportunity as asinh: `direct` needs
    // x*x - 1.0 computed as a *single* rounding (the fma) because x is
    // near 1 at acosh's domain boundary, where it's a catastrophic-
    // cancellation subtraction -- a rounded-then-reused x2 loses exactly
    // the precision that cancellation needs (measured: max ulp 3 -> 700).
    // `direct` and `inv_x2` each need their own x*x in a different
    // rounding context.
    let direct = fma(x, x, -1.0).sqrt();
    // sqrt(x^2-1) = x - 1/(2*x) - 1/(8*x^3) - ..., same two-term expansion
    // (and same 7e-15 bound at this branch's own x >= 2048) as asinh's --
    // see its body comment.
    let inv_x = 1.0 / x;
    let rescaled = fma(-0.5, inv_x, x);
    let s = if x < 2048.0 { direct } else { rescaled };
    let d = (x - 1.0) + s;
    // Same shared-ln_normal merge as asinh (see its own doc comment for
    // the full mechanism): `log1p_finite(d)` and the `ln(x) + LN_2`
    // overflow fallback each paid a full ln-family poly + wrapper,
    // unconditionally, every call. `x` itself (not `ax`: acosh's domain
    // is `x >= 1`, no sign to strip) and `u = 1+d` are both always
    // positive here, so only inf/nan need restoring after the raw
    // `ln_normal` core -- `x < 1.0` (false for NaN) already runs last
    // and independently supplies the out-of-domain NaN, so the trailing
    // overrides only need to cover the in-domain `x >= 1` inf/nan cases.
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
    if x < 1.0 { f32::NAN } else { combined }
}

// atanh(x) ~ x*(1 + x^2/3 + x^4/5 + x^6/7 + ...), a degree-7 minimax
// refit of the odd Taylor series over |x| < 0.25 (idea #69): the
// leading coefficient is pinned to exactly 1.0 (same convention as
// asin_small/sinh_small), and only 3 non-leading terms are needed to
// stay near f32 precision over this narrow domain -- a fitted LP found
// the next two odd terms (x^8, x^10) converge to exactly 0, so this is
// cheaper than the idea's own "~5 odd terms" guess. Same role as
// asin_small/sinh_small: a cheap, cancellation-free small-x numerator
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

/// atanh(x) = 0.5*ln((1+x)/(1-x)), two branches (idea #69): a dedicated
/// small-x poly above for `|x| < 0.25`, and `0.5*log1p(2a/(1-a))` (a
/// single `log1p` call) for the rest, sign restored via `mulsign` --
/// half the `log1p` calls of a naive `0.5*(log1p(x)-log1p(-x))` form,
/// a real throughput win. See IDEAS.md idea #69 for why the single-
/// log1p form needs the dedicated small-x branch to be accurate (its
/// own worst case, un-split, lands inside `|x| < 0.25`). Domain `x` in
/// `[-1, 1]` falls out for free: `2a/(1-a)` is `+inf` at `a=1`
/// (`log1p(inf)=inf`, matching `atanh(+-1)=+-inf`), and is `< -1` for
/// any `a>1` (`log1p` already `NaN` there, matching `atanh`'s domain
/// edge). Current: avg/max ulp 0.0037/2 (exhaustive).
///
/// The `log1p` call is inlined rather than made, because `a = |x| >= 0`
/// makes `v = 2a/(1-a)` non-negative for the whole in-domain half and
/// therefore `u = 1+v >= 1`: `log1p`'s `x == 0.0` signed-zero select and
/// `ln`'s denormal rescale and zero-input `-inf` select are all
/// unreachable here (same "the guard is the licence" lever as
/// `asinh`/`acosh`'s shared-`ln_normal` merge). Only the out-of-domain
/// arms survive, and both are still exactly the cases `log_family_wrapper!`
/// covers: `u < 0` for `a > 1` and the `!(u < inf)` catch that turns
/// `a == 1`'s `u = +inf` into `+inf` and any `NaN` back into `NaN`.
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
    if a < 0.25 { small } else { big }
}

// degree-6 minimax poly (Horner via fma), fitted for acos's sqrt(1-|x|)
// factor; only acos uses it (asin has its own decoupled asin_poly
// below). The leading constant is acos_poly(0) and must equal exactly
// pi/2: `1.5707964` parses to the same bits as
// `std::f32::consts::FRAC_PI_2` (0x3fc90fdb) -- `1.5707963` is one digit
// short of round-trip precision and parses one ulp LOW, which alone cost
// acos most of its average error before being caught. The other six are a
// minimax fit of `acos(x)/sqrt(1-x)` weighted by `sqrt(1-x)/ulp(acos(x))`
// -- the weight that makes "minimax" mean the ulp error of the *shipped*
// combine, which peaks as x -> 1 -- then polished by coordinate descent
// scored on the real exhaustive sweep over all of [-1,1] (both signs; the
// `+PI` branch has its own error structure). Current: acos max ulp 4, avg
// 0.0555 (exhaustive, scored as the whole formula).
#[inline(always)]
fn acos_poly(x: f32) -> f32 {
    let u = 2.34696e-3f32;
    let u = fma(u, x, -1.1324962e-2);
    let u = fma(u, x, 2.7140902e-2);
    let u = fma(u, x, -4.895824e-2);
    let u = fma(u, x, 8.880409e-2);
    let u = fma(u, x, -2.145914e-1);
    fma(u, x, 1.5707964)
}

/// Dedicated asin-only copy of `acos_poly`'s shape (same 7-coefficient
/// Horner form, same `sqrt(1-x)*poly` combine, independently tuned),
/// decoupling asin from acos's protected coefficients: every joint refit
/// attempt died protecting acos's accuracy (see graveyard.md §asin/acos), so
/// not sharing coefficients means acos can't regress no matter what this
/// poly converges to. Fit against `acos(x)/sqrt(1-x)` restricted to
/// asin's actual domain for this branch, `x` in `[0.25, 1)` -- the
/// narrower domain is exactly the freedom the joint fits couldn't use.
#[inline(always)]
fn asin_poly(x: f32) -> f32 {
    let u = 1.3137129e-3f32;
    let u = fma(u, x, -7.664533e-3);
    let u = fma(u, x, 2.2003133e-2);
    let u = fma(u, x, -4.5330178e-2);
    let u = fma(u, x, 8.745548e-2);
    let u = fma(u, x, -2.1434246e-1);
    fma(u, x, 1.5707785)
}

/// acos(x), domain x in [-1,1] (result always in `[0,pi]`, never negative --
/// unlike sin/asinh/etc., acos isn't an odd function, so x=-0.0 has no
/// legitimate negative result the way it does for those). `mulsign`
/// (bit-based sign) and `x < 0.0` (value-based comparison) disagree on
/// exactly one input: `-0.0`, whose sign *bit* is set but whose *value*
/// equals `+0.0`. `mulsign` flipped `y`'s sign there (per the bit), while
/// `x < 0.0` correctly saw "not negative" and skipped the `+pi`
/// correction -- so `acos(-0.0)` came out `-pi/2` instead of the correct
/// `+pi/2`. Fixed by normalizing `-0.0` to `+0.0` before `mulsign` sees
/// it: `x + 0.0` is `-0.0 + 0.0 = +0.0` exactly (IEEE754's defined
/// round-to-nearest behavior for that one case) and a no-op for every
/// other `x`, genuinely negative or not.
#[doc(alias = "acosf")]
#[inline(always)]
pub fn acos(x: f32) -> f32 {
    const PI: f32 = std::f32::consts::PI;
    let a = x.abs();
    let y = (1.0 - a).sqrt() * acos_poly(a);
    mulsign(y, x + 0.0) + if x < 0.0 { PI } else { 0.0 }
}

/// acos(x) in degrees (backlog idea #123): plain composite, same
/// finding as `asind`'s own doc comment -- a rescaled-coefficient fold
/// was tried for `asin` and measured *worse* (real max ulp 16 vs the
/// naive composite's 11), so not attempted again here without new
/// evidence it would fare differently.
///
/// The single-word multiply here is deliberate, and is the one member of
/// this family where [`RAD_TO_DEG_HI`]'s two-word `fma` measures no
/// benefit: `acos`'s own error dominates, and `acosd`'s bulk mass sits
/// where the true answer is exactly `90` (representable, so no bias
/// crosses a rounding boundary). See graveyard.md.
#[inline(always)]
pub fn acosd(x: f32) -> f32 {
    acos(x) * RAD_TO_DEG_HI
}

// acos(x)/pi as its own degree-6 minimax poly (backlog idea #85 gives
// the half-turn form; the coefficients are not acos_poly's rescaled by
// 1/pi -- that rescale was the starting point and a dedicated fit beats
// it). Fitted as an ulp-weighted minimax (LP) of `acos(a)/(pi*sqrt(1-a))`
// against `a`, weighted by the `sqrt(1-a)*P(a)` combine's own sensitivity
// `sqrt(1-a)/ulp(acospi)` for whichever of the two halves (`y` and
// `1-y`) binds harder, then coordinate-descended over the f32
// quantisation. The trailing term stays pinned to exactly 0.5, which is
// what makes acospi(0) exact.
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

/// acos(x)/pi (backlog idea #85), the C23 half-turn convenience family.
/// Verified (exhaustive): a dedicated poly gives avg/max ulp 0.0438/3 vs
/// the naive `acos(x) * (1.0 / PI)` composite's 0.0595/5 -- a real win on
/// both axes with no mca cost (identical instruction stream to plain
/// acos), same verdict as asinpi's own 1/pi fold and the opposite of
/// asind's RAD_TO_DEG fold.
#[inline(always)]
pub fn acospi(x: f32) -> f32 {
    let a = x.abs();
    let y = (1.0 - a).sqrt() * acospi_poly(a);
    mulsign(y, x + 0.0) + if x < 0.0 { 1.0 } else { 0.0 }
}

// Odd approximation asin(x) ~ x * P(x^2), degree 3 in x^2, on |x| < 0.25.
// The leading coefficient is pinned to exactly 1.0 so tiny x returns x
// (its correctly-rounded asin). The other three are a minimax refit over
// [0, 0.25], *not* the odd Taylor series (1/6, 3/40, 15/336): Taylor is
// optimal only at x=0 and leaves the worst case at the 0.25 edge (where
// the next Taylor term is ~4.6e-7 relative), while an equal-degree minimax
// fit spreads that error for a lower max ulp at identical op count. Nearer
// |x|=1 the series would converge slowly (asin's sqrt singularity), so the
// other branch takes over there. See graveyard.md §asin/acos/atan.
#[inline(always)]
fn asin_small(x: f32) -> f32 {
    let x2 = x * x;
    let c0 = 1.0f32;
    let c1 = 0.166666746f32;
    let c2 = 0.074942857f32;
    let c3 = 0.0474379882f32;
    let p = fma(fma(fma(c3, x2, c2), x2, c1), x2, c0);
    x * p
}

/// Two branches, both computed unconditionally and selected (branchless,
/// auto-vectorizes): a dedicated odd minimax poly below
/// `|x| < 0.25` (see asin_small), and `asin(x) = pi/2 - acos(x)` above
/// it, via acos's own well-conditioned `sqrt(1-a) * poly(a)` formula (a
/// shrinking sqrt factor times a smooth bounded poly -- and
/// `pi/2 - acos(a)` doesn't cancel either, since acos(a) is small
/// exactly where pi/2 is O(1); near x=0 that difference IS catastrophic
/// cancellation, which is what the small branch exists to avoid). Sign
/// restored via `mulsign` (asin is odd). The poly is a dedicated
/// `asin_poly`, decoupled from acos's coefficients -- see its doc
/// comment. The crossover sits at 0.27, not 0.25 (asin_small's own fit
/// domain) -- after asin_small's minimax refit tightened its branch, the
/// worst case moved to just above the old 0.25 boundary, still inside
/// asin_small's error curve at that point, so shifting the boundary
/// (coefficients untouched) covers it with the small branch instead.
///
/// The big branch is one `fma`, not a multiply and a subtract: the
/// product `sqrt(1-a)*P(a)` is ~4x larger than the difference it feeds
/// just above the crossover (1.297 against 0.274), so rounding it to
/// f32 first costs ~2 ulp of the *result*. `fma` rounds once, at the
/// result's own magnitude, and is a whole instruction cheaper than the
/// `vmulps`/`vsubps` pair it replaces. Current: max ulp 5, avg 0.019
/// (exhaustive). An earlier three-branch design with a rational
/// mid-branch was strictly worse -- see IDEAS.md §asin/acos for that
/// history and the rejected refit variants.
#[doc(alias = "asinf")]
#[inline(always)]
pub fn asin(x: f32) -> f32 {
    let a = x.abs();
    let small = asin_small(x);
    let big = mulsign(fma(-(1.0 - a).sqrt(), asin_poly(a), FRAC_PI_2), x);
    if a < 0.27 { small } else { big }
}

// 180/pi as a double-f32, for the radians-to-degrees composites: HI+LO
// represents it to ~2^-49 relative, so `fma(y, HI, y*LO)` rounds once at
// the result's own magnitude where a single-word `y * K` rounds twice
// and carries whatever bias the f32 `K` has.
//
// That bias is not small here. `180.0 / std::f32::consts::PI` -- the
// obvious spelling -- computes `180/fl(pi)`, not `fl(180/pi)`, and lands
// one ulp low (`0x42652ee0` against `0x42652ee1`); as a lone multiplier
// that is a systematic **-0.46 ulp** relative error on every result.
// Even the correctly-rounded `fl(180/pi)` still biases every result by
// +0.098 ulp. (`1/pi` and `pi/180` have no such problem -- `1/fl(pi)`
// and `fl(pi)/180` each land on the correctly-rounded constant -- so
// this is specific to 180/pi, not a general rule about deriving
// constants from `PI`.)
//
// HI is deliberately the *low* neighbour rather than the correctly
// rounded one, which is what makes LO positive, which is what keeps
// `-0.0` a signed zero through the `fma`: with a negative LO, `y*LO`
// comes out `+0.0` for `y = -0.0`, and `-0.0 + 0.0` is `+0.0`, silently
// dropping the sign every odd function in this family has to preserve.
// The split point is otherwise free -- either neighbour as HI represents
// 180/pi to the same ~2^-49 once LO is added.
const RAD_TO_DEG_HI: f32 = 57.2957763671875;
const RAD_TO_DEG_LO: f32 = 3.1458948e-6;

// The same double-f32 treatment for 1/pi, for the half-turn composites.
// Here [`FRAC_1_PI`] is *already* the low neighbour -- the correctly
// rounded `fl(1/pi)` sits 0.43 ulp *below* the real 1/pi -- so it also
// serves as the HI word and only this positive tail has to be named.
// A lone `y * FRAC_1_PI` is therefore biased low by 0.34 to 0.68 ulp of
// the result depending on where in its binade the result lands, which is
// the whole of the gap between `atan`/`atan2` and their `pi` composites.
// LO being positive is what keeps `-0.0` signed through the `fma`, for
// exactly the reason spelled out for RAD_TO_DEG_LO above.
const FRAC_1_PI_LO: f32 = 1.28412765e-8;

/// asin(x) in degrees (backlog idea #123). The idea as originally
/// proposed -- "fold 180/pi into the poly/combine constants" -- was tried
/// and measured worse on max ulp than a post-multiply (see IDEAS.md
/// §asin/acos), so this stays a composite; what it does not stay is a
/// single-word multiply. `asin`'s result is already a rounded f32 and
/// `180/pi` is irrational, so a plain `asin(x) * K` rounds twice, and the
/// second rounding is biased by however far the f32 `K` sits from the
/// real `180/pi` (see [`RAD_TO_DEG_HI`]). The `fma` form rounds once, at
/// the result's own magnitude.
///
/// What is left is `asin`'s own error, amplified ~1.8x by the binade
/// shift from `asin`'s `[0.25,0.5)` to degrees' `[8,16)` at the worst
/// case -- `|x|` just above `asin`'s own 0.27 branch crossover, which is
/// exactly where `asin` itself is at its max. That is the whole of the
/// remaining max ulp, and nothing on this side of the composite can
/// reach it.
#[inline(always)]
pub fn asind(x: f32) -> f32 {
    let y = asin(x);
    fma(y, RAD_TO_DEG_HI, y * RAD_TO_DEG_LO)
}

// asin(x)/pi, seeded by rescaling each coefficient by 1/pi (backlog idea
// #85, the half-turn sibling of idea #123's degree fold) and then refit
// in half-turn space directly: NOT assumed safe just
// because idea #123's own RAD_TO_DEG fold measured worse for asind --
// verified separately since it's a different constant, not just a
// relabeling. Real exhaustive fuzz (pilot-tested here before extending
// to acospi/atanpi/atan2pi): folded gives avg/max ulp 0.0159/5 vs the
// plain `asin(x)/PI` composite's 0.2430/7 -- a real win on *both* axes
// here, unlike asind's fold (worse on max ulp there). Plausible reason
// for the opposite verdict: `1/pi` (~0.318) keeps rescaled coefficients
// in a similar-or-smaller magnitude range, while `180/pi` (~57.3)
// inflates them roughly 57x, apparently accumulating more absolute
// rounding per fma step than one final multiply costs. Not a general
// rule either way -- each site needs its own measurement.
#[inline(always)]
fn asinpi_small(x: f32) -> f32 {
    let x2 = x * x;
    let c1 = 0.053056117f32;
    let c2 = 0.023659006f32;
    let c3 = 0.017134449f32;
    // The leading coefficient is 1/pi, and unlike `asin_small`'s exact
    // 1.0 it is irrational: carrying it as the single word `FRAC_1_PI`
    // hands every result in this branch that constant's own 0.43-ulp low
    // bias (see [`FRAC_1_PI_LO`]) -- and this branch is where almost
    // every f32 in `asinpi`'s domain lands. So it is peeled out of the
    // Horner chain and carried as the double-f32 pair: `FRAC_1_PI_LO`
    // takes the trailing slot the leading coefficient used to hold (same
    // three fmas), the tail is scaled by `x` on its own, and the one
    // remaining fma adds `x*FRAC_1_PI` exactly, rounding once at the
    // result's own magnitude. `x*t`'s own rounding is ~0.004 ulp of the
    // result at the branch edge and vanishes from there down.
    //
    // `FRAC_1_PI_LO > 0` is also what keeps `asinpi(-0.0)` at `-0.0`:
    // `x*t` stays `-0.0` rather than becoming `+0.0`, so the fma adds two
    // negative zeros. See RAD_TO_DEG_LO's comment for the same argument.
    let t = fma(fma(fma(c3, x2, c2), x2, c1), x2, FRAC_1_PI_LO);
    fma(x, FRAC_1_PI, x * t)
}

// asin's `0.5 - sqrt(1-a)*P(a)` branch in half-turns, degree 6 over
// a in [0.27, 1). Seeded from `asin_poly`'s coefficients rescaled by
// 1/pi, then refit as an ulp-weighted minimax (LP) of
// `acos(a)/(pi*sqrt(1-a))` -- weight `sqrt(1-a)/ulp(asinpi)`, the
// combine's own sensitivity -- and coordinate-descended over the f32
// quantisation.
#[inline(always)]
fn asinpi_poly(x: f32) -> f32 {
    let u = 0.00047981358f32;
    let u = fma(u, x, -0.0026602724);
    let u = fma(u, x, 0.0073212385);
    let u = fma(u, x, -0.014663323);
    let u = fma(u, x, 0.027931212);
    let u = fma(u, x, -0.068246327);
    fma(u, x, 0.49999586)
}

/// asin(x)/pi (backlog idea #85), the C23 half-turn convenience family.
/// Rescaled dedicated coefficients (see `asinpi_small`'s own doc
/// comment) -- unlike `asind`'s analogous fold (rejected, see IDEAS.md),
/// this one measured as a real win, so it's used here instead of the
/// plain composite.
///
/// The big branch is one `fma`, not a multiply and a subtract, for the
/// same reason as [`asin`]'s -- and the cancellation here is sharper
/// (`0.5 - 0.412 = 0.088` just above the crossover, ~5.7x, against
/// asin's 4.7x), so the product's own rounding was worth even more.
/// Current: max ulp 5, avg 0.016 (exhaustive). The remaining max is the
/// big branch's, at `|x|` just above the crossover; the average is the
/// small branch's, and see `asinpi_small` for why it needs the two-word
/// `1/pi` that `acospi` (trailing constant exactly `0.5`) does not.
#[inline(always)]
pub fn asinpi(x: f32) -> f32 {
    let a = x.abs();
    let small = asinpi_small(x);
    let big = mulsign(fma(-(1.0 - a).sqrt(), asinpi_poly(a), 0.5), x);
    if a < 0.27 { small } else { big }
}

// 3/3 Pade-style rational approximation of atan on [0,1], seeded from a
// least-squares fit and coordinate-descent tuned. Current: atan avg/max
// ulp 0.068/4, atan2 0.069/3 (exhaustive). Numerator and denominator
// evaluate in parallel, so the depth cost over a lower-degree form is
// one fma, not two. See graveyard.md §asin/acos/atan for the fit history
// (including the "zero-move trap" of coordinate-descending a new
// coefficient from 0.0).
#[inline(always)]
fn atan_poly(x: f32) -> f32 {
    let a2 = 0.008830042167832291;
    let a1 = 0.2849778513254418;
    let a0 = 1.1271711055988247;
    let b2 = 5.0166193e-2;
    let b1 = 5.718157e-1;
    let b0 = 1.4605043e0;
    let x2 = x * x;
    let numer = fma(fma(fma(a2, x2, a1), x2, a0), x2, 1.0) * x;
    let denom = fma(fma(fma(b2, x2, b1), x2, b0), x2, 1.0);
    numer / denom
}

/// Straight port of jodiemath's atanf: reciprocates |x| > 1 into range
/// (atan(x) = pi/2 - atan(1/x)) before the poly, matching atan_poly's fit.
#[doc(alias = "atanf")]
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

/// atan(x) in degrees (backlog idea #123): composite -- see `asind`'s own
/// doc comment for why a rescaled-coefficient fold isn't attempted here
/// either, and [`RAD_TO_DEG_HI`] for why the `180/pi` multiply is a
/// two-word `fma` rather than the single constant it reads as.
#[inline(always)]
pub fn atand(x: f32) -> f32 {
    let y = atan(x);
    fma(y, RAD_TO_DEG_HI, y * RAD_TO_DEG_LO)
}

/// atan(x)/pi (backlog idea #85), the C23 half-turn convenience family.
/// A composite over `atan`, but not a single-word multiply: `atan`'s
/// result is already a rounded f32 and `1/pi` is irrational, so a plain
/// `atan(x) * K` rounds twice and the second rounding inherits the f32
/// `K`'s own 0.43-ulp low bias (see [`FRAC_1_PI_LO`]). The `fma` form
/// rounds once, at the result's own magnitude, which leaves essentially
/// only `atan`'s error behind.
///
/// A rescaled-coefficient fold (per asinpi/acospi's own precedent) is
/// deliberately not used here: `atan_poly` is a Pade rational whose
/// numerator and denominator share one unscaled trailing `+1.0`, so
/// folding 1/pi into the numerator's copy alone costs a second broadcast
/// and measures worse on both axes -- see graveyard.md.
#[inline(always)]
pub fn atanpi(x: f32) -> f32 {
    let y = atan(x);
    fma(y, FRAC_1_PI, y * FRAC_1_PI_LO)
}

/// atan(x), `|x| <= 1` contract (backlog idea #61): `atan_poly` alone is
/// already the whole answer over that domain (it's fitted directly
/// against atan on `[0,1]`), so this skips `atan`'s own `1/a`
/// reciprocal, `min`, and `FRAC_PI_2 - y` fold entirely -- for callers
/// who already know their input is bounded (e.g. already reduced via
/// some other identity). Garbage outside `[-1,1]`, same "unchecked
/// tier" convention as this crate's other `_unchecked`/narrow variants.
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn atan_bounded(x: f32) -> f32 {
    mulsign(atan_poly(x.abs()), x)
}

/// Latency-tier atan: a division-free odd degree-17 poly (fit directly
/// against atan(r) over r in `[0,1]`, not derived from atan_poly's own
/// rational) instead of atan_poly's 3/3 Pade form. Removes the division
/// atan_poly's own critical path pays (a division can't start until its
/// numerator/denominator resolve, unlike cbrt's early-starting `rcp`),
/// trading it for more fma depth -- opposite tradeoff to nearly every
/// other function in this crate, so a separate opt-in tier rather than a
/// replacement (same shape as `sinh_throughput`/`cosh_throughput`, just
/// favoring the other axis): slightly better latency, ~8% worse
/// throughput than `atan`. Accuracy is not a tradeoff (avg/max ulp
/// 0.052/3, slightly better than `atan`) -- use this over `atan` only
/// for a value on its own or a serial dependency chain.
///
/// The final combine applies `mulsign` to `p` and `FRAC_PI_2`
/// individually and selects between `sp` and `hpisignx - sp`, instead of
/// selecting on the unsigned value and applying one final `mulsign`:
/// `mulsign(a,x) - mulsign(b,x) == mulsign(a-b,x)` always (`mulsign` is
/// an exact sign-bit XOR, and IEEE754 subtraction commutes exactly with
/// negating both operands), so this is a pure reassociation -- verified
/// bit-identical, with a small real throughput win.
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
    if a < 1.0 { sp } else { hpisignx - sp }
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
///
/// `atan2(NaN, 0.0)`/`atan2(NaN, -0.0)` used to come out `+-FRAC_PI_2`
/// instead of the correct `NaN`: when
/// `nonzerox` is false, `r` collapses to bare `correction =
/// mulsign(FRAC_PI_2 - hpisignx, y)` with no `atan(y/x)` call at all --
/// but `mulsign` only ever reads `y`'s *sign bit*, it doesn't propagate
/// `y` being NaN, so a NaN `y` here silently degrades to a finite
/// `+-FRAC_PI_2` depending on which way that one bit happened to be set.
/// Every other input combination avoids this because `x != 0.0` routes
/// through `atan(y/x)`, and `y/x` is itself NaN whenever `y` is NaN
/// (`atan` propagates it correctly from there) -- only the `x == 0`
/// branch bypasses that path entirely. Fixed with an explicit trailing
/// override.
#[doc(alias = "atan2f")]
#[inline(always)]
pub fn atan2(y: f32, x: f32) -> f32 {
    let nonzerox = x != 0.0;
    let nonzeroy = y != 0.0;
    let bothzero = !nonzerox && !nonzeroy;
    let hpisignx = if nonzerox || bothzero { mulsign(FRAC_PI_2, x) } else { 0.0 };
    let correction = mulsign(FRAC_PI_2 - hpisignx, y);
    let r = if nonzerox { atan(y / x) + correction } else { correction };
    let r = if y.is_nan() { f32::NAN } else { r };
    // atan2(+-inf, +-inf): y/x is inf/inf, which is NaN, so the general
    // formula above can't produce an answer here at all. IEEE754/C99
    // define a canonical result by quadrant regardless (+-pi/4 or
    // +-3pi/4) -- not derived from any real ratio, since there isn't one
    // at true infinity, just a fixed convention.
    let bothinf = x.is_infinite() && y.is_infinite();
    let inf_result = mulsign(if x.is_sign_negative() { 3.0 * FRAC_PI_4 } else { FRAC_PI_4 }, y);
    if bothinf { inf_result } else { r }
}

/// atan2, latency tier (backlog idea #60): identical wrapper to [`atan2`]
/// (same zero/NaN/inf edge-case fixes, none of which depend on which
/// "atan"-shaped core fills in the ordinary case) but calls
/// [`atan_latency`] instead of [`atan`] for the `y/x` term. `atan2`
/// today stacks two divisions serially -- its own `y/x`, then
/// `atan_poly`'s internal rational division on top -- since
/// `atan_latency` is division-free past its own initial reciprocal,
/// this removes the second one from the critical path entirely. Same
/// tradeoff [`atan_latency`] documents on its own: better latency,
/// worse throughput than [`atan2`].
#[inline(always)]
pub fn atan2_latency(y: f32, x: f32) -> f32 {
    let nonzerox = x != 0.0;
    let nonzeroy = y != 0.0;
    let bothzero = !nonzerox && !nonzeroy;
    let hpisignx = if nonzerox || bothzero { mulsign(FRAC_PI_2, x) } else { 0.0 };
    let correction = mulsign(FRAC_PI_2 - hpisignx, y);
    let r = if nonzerox { atan_latency(y / x) + correction } else { correction };
    let r = if y.is_nan() { f32::NAN } else { r };
    let bothinf = x.is_infinite() && y.is_infinite();
    let inf_result = mulsign(if x.is_sign_negative() { 3.0 * FRAC_PI_4 } else { FRAC_PI_4 }, y);
    if bothinf { inf_result } else { r }
}

/// atan2(y,x) in degrees (backlog idea #123): the composite -- see
/// `asind`'s own doc comment for why a rescaled-coefficient fold isn't
/// attempted here either -- plus one branch the radian version doesn't
/// need.
///
/// `180/pi` is *greater than 1*, so unlike `atan2pi`'s `1/pi` it can
/// carry a denormal `atan2` result up into a binade with more mantissa
/// bits than the denormal ever had. Those bits are gone: `atan2` is
/// exactly right there (its own error at every one of these inputs
/// measures 0 ulp), it just cannot represent the answer, because `y/x`
/// underflowed inside it. Scaling afterwards cannot invent them back,
/// and the loss survives into results that are themselves perfectly
/// normal -- ~1.3% of the f32 plane scored over 8 ulp, up to 30, a
/// sixth of that with a normal result.
///
/// A denormal `atan2` means the angle is far inside `atan`'s linear
/// regime, where `atan2(y,x)` is just `y/x`, so `(180/pi) * atan2(y,x)`
/// is exactly `y * ((180/pi) / x)`. The association matters and is the
/// whole point: `(180/pi) / x` is always normal on this branch, so `y`
/// -- denormal or not -- is multiplied straight into the result's own
/// binade with a single rounding, and no intermediate is ever denormal.
/// Scaling `y` first instead (`(y * 180/pi) / x`) fixes the original
/// case but reintroduces the identical bug one step earlier for a
/// denormal `y`, where the scaled numerator is still denormal: measured
/// 77 max ulp that way, against 5 for this one.
///
/// `(180/pi) / x` is normal because this branch implies `|y| <
/// 1.2e-38 * |x|`, so a nonzero `y` forces `|x| >= 1.2e-7`, and `|x|`
/// can never exceed `f32::MAX`. It depends only on the arguments, so it
/// is off `atan2`'s dependency chain, but it is not free: throughput
/// +13%, latency +10%. `x == 0.0` is excluded because the quotient
/// degenerates there exactly where `atan2` already answers `+-0.0`
/// correctly; every other `x == 0` case gives `+-90` and never reaches
/// this branch at all.
#[inline(always)]
pub fn atan2d(y: f32, x: f32) -> f32 {
    let r = atan2(y, x);
    let normal = fma(r, RAD_TO_DEG_HI, r * RAD_TO_DEG_LO);
    // Single-word, and specifically the HI word rather than the
    // correctly-rounded `fl(180/pi)`: the two-word form would need a
    // second division, and swapping in the correctly-rounded single
    // constant -- which does measure better on this branch's average --
    // costs a 1-ulp regression at a pinned denormal edge case that is
    // currently correctly rounded. See graveyard.md; this branch exists
    // to recover *binade* bits from an underflowed `y/x`, not the last
    // ulp of a denormal.
    let tiny = y * (RAD_TO_DEG_HI / x);
    if r.abs() < f32::MIN_POSITIVE && x != 0.0 { tiny } else { normal }
}

/// atan2(y,x)/pi (backlog idea #85, C23 half-turn family): a composite
/// over `atan2`, with the same two-word `1/pi` multiply `atanpi` uses and
/// for the same reason (see [`FRAC_1_PI_LO`]). A rescaled-coefficient
/// fold is not used here either -- `atan2` composes the very `atan_poly`
/// Pade rational whose shared trailing constant makes that fold a loss in
/// `atanpi`, and `atan2` is the larger, more branch-heavy caller of it.
#[inline(always)]
pub fn atan2pi(y: f32, x: f32) -> f32 {
    let r = atan2(y, x);
    fma(r, FRAC_1_PI, r * FRAC_1_PI_LO)
}

/// atan2 without the x==0/both-zero/both-infinite special cases: contract
/// is x != 0.0 (and not both x and y infinite). Drops the nonzerox/
/// nonzeroy/bothzero selects and the whole bothinf branch atan2 pays on
/// every call -- see its own doc comment for exactly what those handle.
/// Bit-identical to atan2 whenever the contract holds.
#[inline(always)]
pub fn atan2_unchecked(y: f32, x: f32) -> f32 {
    let hpisignx = mulsign(FRAC_PI_2, x);
    let correction = mulsign(FRAC_PI_2 - hpisignx, y);
    atan(y / x) + correction
}

/// atan2 folded into a single positive turn (backlog idea #143), the
/// geo/graphics "bearing" convention: `atan2`'s own `(-pi, pi]` range
/// only needs a `+2*pi` fold on the negative half to land there,
/// branchless-select same as every other seam in this crate. NaN
/// propagates unchanged (the fold is a masked add and `NaN + 2*pi` is
/// NaN); `+pi` itself (the one boundary `atan2` can return) is already
/// in-range, no wraparound needed. The range is `[0, TAU]` *closed*,
/// not `[0, 2*pi)`: `f32::consts::TAU` is the nearest f32 to `2*pi` and
/// it sits just above it, so it is the correctly-rounded answer for
/// every angle in the last half-ulp of the turn and there is no way to
/// both round correctly and stay strictly under `2*pi`. `-0.0` is never
/// returned.
///
/// The fold is keyed on `y`'s sign bit, not on `atan2`'s own sign.
/// Cheaper: the mask comes straight off an argument, so it is ready long
/// before `atan2` is, instead of extending that result's dependency
/// chain the way `r < 0.0` does. And correct where `r < 0.0` is not:
/// `atan2` returns `-0.0` for a whole slab of genuinely negative `y`.
/// Once `|y/x|` falls below the smallest subnormal the quotient flushes
/// to `-0.0`, `atan(-0.0)` is `-0.0`, and the `x > 0` correction term is
/// `-0.0` as well, so nothing downstream still carries the fact that the
/// angle was ever nonzero; `r < 0.0` reads that as non-negative, skips
/// the fold, and answers `~0` where the true angle is `~2*pi` -- a full
/// turn out, and not a rare corner, since both operands uniform over all
/// f32 bit patterns puts ~2% of the plane in that slab. A sign bit
/// cannot underflow away.
///
/// The two disagree for a second reason at exactly one input, `y ==
/// -0.0` with `x >= +0.0`, where `atan2` is `-0.0` and this folds to
/// `TAU`. That is deliberate: it reads the sign bit as "the turn
/// approached from below", the same meaning it carries in the underflow
/// slab, and it is what keeps `-0.0` out of the output range of a
/// function whose whole job is to return a non-negative angle.
#[inline(always)]
pub fn atan2_pos(y: f32, x: f32) -> f32 {
    let r = atan2(y, x);
    if y.is_sign_negative() { r + std::f32::consts::TAU } else { r }
}

/// Straight port of jodiemath's tanf: sin(x)/cos(x), over `cos`'s domain
/// `|x| < 2^22 * pi` (~1.32e7), the narrower of the two (see their doc
/// comments).
#[doc(alias = "tanf")]
#[inline(always)]
pub fn tan(x: f32) -> f32 {
    sin_over_cos_domain(x) / cos(x)
}

// [`sin`] restricted to [`cos`]'s domain, and built from [`cos`]'s own
// `frac_x_over_pi!` grid rather than sin's wider multiple-of-4 one. Same
// `q`, same residual, same output as `sin` for every x inside that domain
// -- a ratio is only right as far out as its *denominator* is, so sin's
// extra two binades would buy `tan` nothing while costing it a second,
// non-shared reduction (mca: tan 105 -> 120 instrs, BlockRT 31 -> 36,
// latency 78 -> 89, for byte-identical output).
#[inline(always)]
fn sin_over_cos_domain(x: f32) -> f32 {
    let (nb, _, fc) = frac_x_over_pi!(x, ROUND_MAGIC);
    // `fc + nb` is a second magic round, of `x/pi` this time instead of
    // `x*FRAC_1_PI`: `nb` is already `ROUND_MAGIC + n` on the integer grid,
    // so adding the (corrected, |fc| <= 0.67) fraction back rounds to
    // `ROUND_MAGIC + round(x/pi)` and leaves the parity bit in `qb`'s low
    // mantissa bit. `sin` cannot reuse this shape -- its coarse grid is
    // too wide for `nb` to sit on the integer grid at all.
    let qb = fc + nb;
    let s = pi_reduce_and_poly!(x, qb - ROUND_MAGIC);
    f32::from_bits(s.to_bits() ^ (qb.to_bits() << 31))
}

// degree-6 minimax poly feeding erf's exp2-based tail (|x| >= 0.28),
// Estrin (3 fma's deep instead of Horner's 6, accuracy-neutral here --
// fma reassociation has to be checked per poly, not assumed either way).
// Coefficients are an ulp-weighted Chebyshev LP fit, weighted by the
// linearized sensitivity of erf's final `1 - 2^poly` combine. This is
// the `|x| >= 0.28` arm only; `erf`'s own overall numbers are on `erf`.
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

/// A Pade approximant near 0 (where the tail form loses precision to
/// cancellation), the exp2-based tail elsewhere. `|x|` MUST be clamped
/// before `erf_poly` sees it: it's a plain degree-6 polynomial, and its
/// positive leading coefficient means it eventually turns around and
/// grows to +inf for large |x| instead of staying deeply negative
/// (erf_poly(9) ~ -92, erf_poly(20) ~ +8698), so an unbounded
/// `exp2(erf_poly(|x|))` gave `erf(50) = -1.02e17` and `erf(+-inf) =
/// NaN` instead of +-1. The bound 10 matches erfc's clamp, comfortably
/// past where erf has saturated (erf_poly(10) = -83.8, safely negative);
/// `erf_poly`'s output over that whole clamped domain stays inside
/// `[-92, 0]` (never approaching `exp2`'s unchecked `[-126, 128)` bound,
/// let alone leaving it), so the extra `exp2_checked` insurance was never
/// reachable: `exp2` suffices. The one input that bypasses the `xa_bounded`
/// clamp is NaN itself (`NaN > 10.0` is false), but that propagates to NaN
/// through `erf_poly` before `exp2` ever sees it, and `exp2`'s bit-twiddled
/// exponent field only feeds a NaN-tainted multiply/fma from there, so the
/// result stays NaN regardless of that field's garbage value.
#[doc(alias = "erff")]
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
    // The Pade numerator's constant term is `erf'(0) = 2/sqrt(pi)`, pinned
    // by the approximant rather than fitted -- and it lands almost exactly
    // on an f32 tie, 0.49 ulp above the nearer representable value. A
    // single-word constant therefore carries that 0.49 ulp as *bias*, not
    // noise, into every result small enough that `A*x2` has vanished
    // beneath it: `x * fl(2/sqrt(pi))` is systematically one ulp high over
    // roughly half the f32 domain by bit pattern. Splitting it two-word and
    // folding LO into the same fma the rest of the numerator already needs
    // rounds once, at the result's own magnitude (same construction as
    // `atanpi`'s two-word `1/pi`).
    let numer = fma(x, f32::from_bits(0x3f906ebb), x * fma(f32::from_bits(0x3f174f6e), x2, f32::from_bits(0xb37bd649)));
    let denom = fma(fma(f32::from_bits(0x3e3e2be3), x2, f32::from_bits(0x3f5b6db7)), x2, 1.0);
    let a = numer / denom;
    let b = mulsign(1.0 - exp2(erf_poly(xa_bounded, x2)), x);
    if xa < 0.28 { a } else { b }
}

// The `|x|` clamps `erfc` and `erfcx` feed their `x*x` through. Both are
// set by one rule: make the exponent they hand `exp_reduce!` provably
// inside its valid range, so neither pays for a clamp `exp_checked`
// would have applied. The bounds that make them correct are not
// obvious, and each is one edit away from silently breaking, so they
// are proved at compile time rather than commented.
//
// `erfc` needs `-xs^2` at or above `EXP_CLAMP_LO` (so the reduction is
// valid) while still being far enough below zero that `e^-p` rounds to
// exactly 0 -- which is what makes clamping legal in the first place,
// since the true `erfc` has already reached exactly 0.0f32 by |x|~10.05.
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
// complementary error function, i.e. erfc(xa)*exp(xa^2). Reciprocal
// variable `v = 1/(2 + xa)` (so v in (0, 1/2], the whole half-line
// compressed into a finite interval), then a degree-10 minimax
// polynomial: `erfcx(xa) = v*P(v)`.
//
// The explicit `v` factor is what makes this exact on the *whole*
// domain instead of needing a clamp: `erfcx(xa) ~ 1/(xa*sqrt(pi))` as
// xa -> inf, and `v ~ 1/xa` there, so `P(0) = 1/sqrt(pi)` reproduces
// the asymptote by construction rather than having to be fitted (it is
// c0 below). Nothing here can overflow for any finite input, and `v` is
// 0 at xa = inf, so `erfcx` and `erfc` saturate correctly with no domain
// check at all. The old n/d rational this replaces was a degree-4/4 in
// `xa` itself, which needed a hard clamp to |xa| <= 10 to keep its
// denominator from overflowing -- and that clamp *froze* `erfcx` at
// `erfc_rational(10)` forever past x=10 (unbounded relative error).
//
// c0 carries that asymptote by itself, and it is held in *two* f32
// words: `c[0]` below is the low word, the high word (0x3f106ebb) is
// peeled out of the polynomial and applied by the final fma. One word
// cannot do it, because 1/sqrt(pi) lands almost exactly on an f32 tie --
// 0.49 ulp above the nearer representable value, 0.51 ulp below the
// other -- so *either* single-word choice leaves ~0.5 ulp of the
// constant behind. Out in the tail that leftover is the whole error, and
// it is bias rather than noise: the polynomial has collapsed to its
// constant term, the result is `c0*v` and nothing else, so a fixed
// relative offset in c0 arrives as a fixed, one-signed ~0.63 ulp offset
// in every result past |x| ~ 10. Two words and one fma remove it.
//
// Fitted by a relative-error Chebyshev LP over v in (0, 1/2] (an even
// sampling of v is a `1/x` sampling of the tail, so the fit is
// naturally weighted where the reciprocal variable resolves), then
// coordinate-descent polished on the f32 grid against the exact
// evaluation order below (`tune.rs`'s `erfcx` target), under the side
// constraint that `xa = 0` (where `2+xa` and `1/2` are both exact)
// reproduces `erfcx(0) = 1.0` bit-exactly. That constraint is load
// bearing and is what the polished low-order coefficients are holding:
// with an exact c0 the evaluation lands a hair the wrong side of the
// boundary on its own, so c1..c4 and c10 each sit an ulp off the LP's
// own values to pull it back. It is also the one property descent
// cannot be handed as an objective -- re-check it against edgecheck's
// pin after any refit. Estrin-grouped, 4 fma deep instead of Horner's
// 10: the division is already on the critical path ahead of it, and
// Horner measured only ~0.5 ulp better on max.
//
// What is left is dominated not by the fit (~0.5 ulp) but by forming
// `v` itself: rounding `2+xa` and then the reciprocal costs ~2.3 ulp
// that no polynomial can recover, since `erfcx` has a nonzero slope at
// 0 while `v` is stationary in relative terms there. Compensating that
// needs the exact residual of `2+xa` (4-6 more ops, correct ordering
// included) for ~1.5 ulp -- measured, rejected on cost.
#[inline(always)]
fn erfcx_pos(xa: f32) -> f32 {
    let v = 1.0 / (2.0 + xa);
    // c[0] is the *low* word of a two-word `1/sqrt(pi)`; the high word is
    // peeled out of the polynomial and applied in the final fma below.
    let c: [f32; 11] = [
        f32::from_bits(0xb2fbd649), 1.1283773, 1.974964, 2.8070478, 2.9756768, -3.7488432, 17.02367,
        -117.490135, 255.59447, -243.95302, 90.238,
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
    // `erfcx(x) -> 1/(x*sqrt(pi))` as x grows, so out in the tail this
    // whole polynomial has collapsed to its constant term and the result
    // is `v * 1/sqrt(pi)` and nothing else. 1/sqrt(pi) sits 0.49 ulp above
    // the nearer f32 (a near-tie: the other side is 0.51 ulp below), so
    // *either* single-word choice leaves ~0.5 ulp of the constant as pure
    // bias -- ~0.63 ulp of it in the result, in one direction, for every x
    // past ~10. Two words and one fma remove it outright, and the split is
    // free of the usual fit-perturbation worry because this coefficient is
    // the function's own asymptote, not a fitted parameter.
    fma(v, f32::from_bits(0x3f106ebb), v * fma(fma(t8, v4, hi), v4, lo))
}

/// `erfc(x) = exp(-x^2)*erfcx(|x|)` for `x >= 0`, reflected through
/// `erfc(x) = 2 - erfc(-x)` for `x < 0`. Both factors need care:
///
/// The Gaussian factor is the one that used to dominate this function's
/// error (max ulp 109). `exp(-x^2)` amplifies an *absolute* error in its
/// exponent into a *relative* error in the result, and `x*x` alone --
/// a single f32 rounding on a value up to 100 -- is already ~6e-6 of
/// absolute exponent error, i.e. dozens of ulp, before any exponential
/// runs. Fixed by keeping the square exactly: `x*x = p + pe` with
/// `pe = fma(x, x, -p)` the exact residual, then `exp(-(p+pe)) =
/// exp(-p)*exp(-pe) ~ exp(-p)*(1 - pe)` -- `|pe| <= 4e-6` here, so the
/// dropped second-order term is under 1e-11, and the whole correction
/// costs one extra fma folded into the `erfcx` factor. Routing through
/// [`exp_checked`] rather than `exp2_checked` is the other half: the
/// naive `exp2(-x*x*LOG2_E)` rounds a *second* time when it forms that
/// product, and f32's own `LOG2_E` is not even precise enough to carry
/// an exponent of magnitude ~144 (see `exp`'s doc comment) -- its
/// Cody-Waite reduction handles both.
///
/// The `erfcx` factor is [`erfcx_pos`] (see its comment): a degree-10
/// polynomial in `1/(2+|x|)`, ~0.5 ulp of fit error against the old
/// rational's ~15.
///
/// No clamp on `x` anywhere except [`ERFC_XS_CLAMP`] on the value being
/// squared -- far past where `erfc` has decayed under the smallest
/// denormal, so it can only ever fire on inputs whose result is exactly
/// `0.0`, and NaN passes through it untouched (`NaN > c` is false) and
/// propagates. That clamp does double duty: it also bounds `-x*x`
/// tightly enough that the exponential needs *no* clamp of its own (see
/// the const assertions on it), which is why this calls `exp_reduce!`
/// rather than [`exp_checked`].
///
/// Current: max ulp 6, avg 0.1949 (accuracy.rs's fuzz sweep restricted to
/// `|x| <= 10`, which is where the f64 reference stops being usable --
/// past ~10.05 the true `erfc` rounds to exactly `0.0` (or `2.0` for
/// `x < 0`) and this returns exactly that, pinned in edgecheck).
#[doc(alias = "erfcf")]
#[inline(always)]
pub fn erfc(x: f32) -> f32 {
    // The x<0 reflection with no compare and no select. `w` shifts x's
    // sign bit straight into the exponent field (0x4000_0000 is 2.0f),
    // and `mulsign` applies the same bit to `y`, so the two arms are
    // `y + 0` and `-y + 2` -- `2 - y` with the single rounding the old
    // `fma(y, z, w)` gave it, bit-identical, in integer ops on ports the
    // polynomial's fmas are not contending for.
    let w = f32::from_bits((x.to_bits() >> 1) & 0x4000_0000);
    let xa = x.abs();
    let xs = if xa > ERFC_XS_CLAMP { ERFC_XS_CLAMP } else { xa };
    let p = xs * xs;
    let pe = fma(xs, xs, -p);
    let r = erfcx_pos(xa);
    // No clamp on the exponent at all: `ERFC_XS_CLAMP` is chosen so
    // `-p` cannot leave `exp_reduce!`'s valid range, which the two const
    // assertions beside it check. NaN reaches here as NaN (`NaN > c` is
    // false), and comes back out through `r`.
    let e = exp_reduce!(-p);
    let y = e * fma(-r, pe, r);
    mulsign(y, x) + w
}

// erfinv's central branch (backlog idea #66): erfinv(x) = x*P(x^2) for
// |x| <= 0.7, degree 8 in u=x^2 (Estrin-grouped, same idiom as erf_poly).
// Coefficients are a real least-squares fit (scipy) of erfinv(x)/x
// against u over [0,0.7], not a transcription of any published
// algorithm's constants.
#[inline(always)]
fn erfinv_central_poly(u: f32) -> f32 {
    let c: [f32; 9] = [
        8.8622695e-1,
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

// erfinv's tail branch: `erfinv(x) = sign(x)*sqrt(w)*Q(w)`,
// `w = -ln(1-x^2)`, for `|x| > 0.7`. An ulp-weighted minimax (LP) fit of
// `erfinv(x)/sqrt(w)` against `w` over `w` in
// `[-ln(1-0.7^2), -ln(1-x_max^2)]` (`x_max` = the largest f32 below 1.0,
// so the fit's own domain exactly matches what an f32 caller can ever
// actually reach), weighted by the `sqrt(w)*Q` combine's own sensitivity
// `sqrt(w)/ulp(erfinv)` and then coordinate-descended over the f32
// quantisation -- a plain least-squares fit left ~3.4x more idealized
// error than the same degree can reach.
#[inline(always)]
fn erfinv_tail_poly(w: f32) -> f32 {
    let c: [f32; 9] = [
        8.862468e-1,
        1.0395482e-2,
        -2.2551943e-4,
        -1.02430866e-4,
        1.6865362e-5,
        -1.3347174e-6,
        6.0056976e-8,
        -1.4636747e-9,
        1.4956567e-11,
    ];
    let w2 = w * w;
    let w4 = w2 * w2;
    let l0 = fma(c[1], w, c[0]);
    let l1 = fma(c[3], w, c[2]);
    let l2 = fma(c[5], w, c[4]);
    let l3 = fma(c[7], w, c[6]);
    let r0 = fma(l1, w2, l0);
    let r1 = fma(l3, w2, l2);
    fma(c[8], w4 * w4, fma(r1, w4, r0))
}

/// Inverse error function (backlog idea #66): the sampling/ML staple
/// (inverse-CDF / Box-Muller-style transforms build on this). Two
/// branches, same shape as this crate's other `erf`/`erfc` splits:
/// `x*P(x^2)` for `|x| <= 0.7`, `sign(x)*sqrt(w)*Q(w)` (`w = -ln(1-x^2)`)
/// past it, where `erfinv` itself grows without bound as `|x| -> 1`.
/// `w` is computed as `-log1p(-x*x)`, reusing this crate's own
/// cancellation-safe `log1p` rather than `-(1.0-x*x).ln()`, which would
/// reintroduce exactly the precision loss `log1p` exists to avoid right
/// where it matters most (`x` close to `+-1`, `1-x*x` close to `0`).
/// `|x| > 1` needs no explicit domain-error handling: `-x*x < -1` there,
/// so `log1p`'s own existing domain guard already gives `NaN`, which
/// propagates through `sqrt`/the poly/`mulsign` unchanged.
///
/// `|x| == 1.0` exactly *does* need an explicit override, found by
/// fuzzing, not assumed: `w` correctly reaches `+inf` there
/// (`log1p(-1.0) == -inf`), but `erfinv_tail_poly`'s Estrin grouping
/// evaluates several partial sums independently before combining them,
/// and at `w=inf` different groups overflow to *opposite-signed*
/// infinities depending on their own local coefficient signs (unlike a
/// plain Horner chain, which stays consistently signed once the leading
/// term dominates) -- so the combine hits a genuine `-inf + inf = NaN`
/// instead of the correctly-signed `+-inf` erfinv actually has there.
#[inline(always)]
pub fn erfinv(x: f32) -> f32 {
    let u = x * x;
    // `-log1p(-u)`, with everything this call site can't reach removed
    // (the "guard is the licence" lever, see `atanh`/`log1p`). `t = 1-u`
    // is never denormal (`1 - fl(x*x)` is either `0` or `>= 2^-24`), and
    // `-u == 0` only at `x == 0`, where the central arm is selected -- so
    // both `ln`'s rescale and `log1p`'s signed-zero select are dead. The
    // `-inf` arm goes too: `t == 0` happens *exactly* when `|x| == 1`
    // (verified exhaustively, not argued), which the trailing override
    // below already owns. What is left collapses to a single select:
    // `t > 0.0` is false for `t <= 0` (`|x| > 1`, a real domain error)
    // and for NaN alike, and both want `NaN`. That does canonicalise the
    // NaN it returns rather than forwarding the input's payload, which is
    // fine by this crate's convention -- a NaN is a NaN, the same rule
    // `ulp_diff` and `worst_corpus` follow.
    let nu = -u;
    let t = 1.0 + nu;
    let c = nu - (t - 1.0);
    // No `corr.is_finite()` guard: `c` is the exact rounding error of
    // `1 + nu` whenever `|nu| <= 1`, so `corr` is finite for every `t`
    // this arm is selected for -- the only inputs that make it `NaN`
    // (`t <= 0`, i.e. `|x| >= 1`) already take the `NaN` arm below.
    let corr = c / t;
    let w = -if t > 0.0 { ln_normal(t, 0.0) + corr } else { f32::NAN };
    let central = x * erfinv_central_poly(u);
    let tail = mulsign(w.sqrt() * erfinv_tail_poly(w), x);
    let normal = if x.abs() <= 0.7 { central } else { tail };
    if x.abs() == 1.0 { f32::INFINITY.copysign(x) } else { normal }
}

// `norm_cdf`'s counterpart to `ERFC_XS_CLAMP`, in `x`'s units rather than
// `x/sqrt(2)`'s. Same two conditions, and they are asserted the same way:
// `x^2/2` must stay inside `exp_reduce!`'s valid range, and must already
// be far enough below zero that `e^-x^2/2` rounds to exactly 0 -- the
// true `norm_cdf` has reached exactly 0.0f32 by `x ~ -14.2`.
const NORM_CDF_XS_CLAMP: f32 = 14.44;
const _: () =
    assert!((NORM_CDF_XS_CLAMP as f64) * (NORM_CDF_XS_CLAMP as f64) * 0.5 <= -(EXP_CLAMP_LO as f64));
const _: () =
    assert!((NORM_CDF_XS_CLAMP as f64) * (NORM_CDF_XS_CLAMP as f64) * 0.5 > 103.97207708399179);

/// Standard normal CDF, `Φ(x) = 0.5*erfc(-x/sqrt(2))` (backlog idea
/// #67). Written out as `erfcx(|x|/sqrt(2)) * e^(-x^2/2)` rather than
/// composed on top of [`erfc`], because the composition's *argument* is
/// what limits it: `fl(x/sqrt(2))` carries a relative error of `2^-24`,
/// and `erfc`'s dominant `e^(-z^2)` factor amplifies that by `2*z^2`, so
/// a thin `0.5*erfc(-x*FRAC_1_SQRT_2)` reaches ~190 ulp near `x = -12.8`
/// no matter how accurate `erfc` itself is. Squaring `x` *before*
/// dividing by two removes the amplified rounding entirely -- `x^2/2` is
/// one rounding of `x^2` plus an exact halving -- and `pe` compensates
/// even that one. `erfcx` keeps its own argument in `x/sqrt(2)` units,
/// where the `2^-24` is harmless: `d(ln erfcx)/dz ~ -1/z`, so it stays a
/// `2^-24` relative error instead of being amplified.
///
/// The tail reflection is `erfc`'s own trick, one power of two down:
/// `w` shifts `x`'s sign bit into the exponent field, so the two arms
/// are `y + 0` and `-y + 2`, and the outer `0.5` (exact) turns those
/// into `Φ = y/2` and `Φ = 1 - y/2`. Saturation to exactly `0`/`1` at
/// both tails still comes for free.
#[inline(always)]
pub fn norm_cdf(x: f32) -> f32 {
    let xa = x.abs();
    // `erfcx_pos`, not `erfcx`: the argument is an absolute value, so
    // erfcx's own x<0 arm is dead, and LLVM does not prove that -- taking
    // the public wrapper left a whole second `exp_reduce!` (two
    // `exp2_field_split`s in the asm) in the region. x/sqrt(2) is formed
    // here and nowhere else, and erfcx is well-conditioned in it.
    let r = erfcx_pos(xa * std::f32::consts::FRAC_1_SQRT_2);
    let xs = if xa > NORM_CDF_XS_CLAMP { NORM_CDF_XS_CLAMP } else { xa };
    // p + pe == xs^2/2 exactly: the halving is exact, so this is just
    // `two_prod`'s error term on a single multiply. The clamp above is
    // what keeps `pe` from becoming `inf*inf - inf == NaN` at x = +-inf.
    let h = 0.5 * xs;
    let p = h * xs;
    let pe = fma(h, xs, -p);
    let e = exp_reduce!(-p);
    let y = e * fma(-r, pe, r);
    // `-x` as a sign carrier only (Φ(x) uses erfc(-x/sqrt(2))), so the
    // negation is a bit flip and stays exact at +-0.0: x = +0.0 takes the
    // `1 - y/2` arm and x = -0.0 the `y/2` arm, and both are exactly 0.5.
    let nx = -x;
    let w = f32::from_bits((nx.to_bits() >> 1) & 0x4000_0000);
    0.5 * (mulsign(y, nx) + w)
}

/// Inverse of `erfc` (backlog idea #139): `erfc(z) = 1 - erf(z)`, so
/// `erfc_inv(y) = erfinv(1 - y)` directly -- `1 - y` needs no
/// cancellation-safe handling the way `erfinv`'s own `w` does, since
/// `erfc_inv`'s intended use (small tail-probability `y` near `0`) lands
/// `1-y` close to `1`, already `erfinv`'s own well-conditioned regime
/// (Sterbenz-exact subtraction of two close, same-magnitude values, not
/// the "large cancellation" case that needs a dedicated correction).
#[inline(always)]
pub fn erfc_inv(y: f32) -> f32 {
    erfinv(1.0 - y)
}

/// Probit, the standard normal quantile function (backlog idea #139):
/// inverse of [`norm_cdf`], `probit(p) = sqrt(2)*erfinv(2p-1)` --
/// derived directly from `norm_cdf`'s own definition (`norm_cdf(x) =
/// 0.5*erfc(-x/sqrt(2))`, solved for `x` via `erfc_inv`/`erfinv`'s
/// oddness), not a separately-fit approximation. Completes the
/// sampling-stack trio with [`erfinv`]/[`erfc_inv`] (inverse-CDF
/// transforms, Box-Muller-style generators).
#[inline(always)]
pub fn probit(p: f32) -> f32 {
    std::f32::consts::SQRT_2 * erfinv(fma(2.0, p, -1.0))
}

/// Standard normal PDF, `φ(x) = exp(-x^2/2)/sqrt(2*pi)` (backlog idea
/// #67). The exponent is formed as a compensated square: `-0.5*x*x` is a
/// single rounding of magnitude `ulp(x^2/2)/2`, but it lands in an
/// *exponent*, where an absolute error `d` is a relative error `d` in
/// the result -- so at `|x| ~ 13` (`x^2/2 ~ 85`) that one rounding is
/// already ~70 ulp on its own, which is exactly where the naive form's
/// max sat. `p + pe == xs^2/2` exactly, and `fma(-pe, e, e)` applies
/// `e^-pe ~ 1 - pe` to put the missing part back.
///
/// The clamp is what keeps `pe` finite (`inf*inf - inf` is `NaN`), and
/// is placed where `e^-p` has already rounded to exactly 0, so it is a
/// no-op on the value. `exp_reduce!` rather than `exp_checked` for the
/// same reason [`erfc`] uses it: the clamp above already proves the
/// argument is in range, asserted next to the constant.
#[inline(always)]
pub fn norm_pdf(x: f32) -> f32 {
    const INV_SQRT_2PI: f32 = 0.3989422804014327;
    let xa = x.abs();
    let xs = if xa > NORM_CDF_XS_CLAMP { NORM_CDF_XS_CLAMP } else { xa };
    let h = 0.5 * xs;
    let p = h * xs;
    let pe = fma(h, xs, -p);
    let e = exp_reduce!(-p);
    INV_SQRT_2PI * fma(-pe, e, e)
}

// dawson's central branch (backlog idea #138): `x*P(u)/Q(u)`, `u=x^2`,
// degree 5/5, Estrin-grouped with each side's top pair folded in at the
// `u^2` level so `u^4` is never formed (exp_r_poly!'s own fold, applied
// to the numerator and the denominator alike). Real fit (scipy/HiGHS) of
// `dawsn(x)/x` against `u` over `|x| <= 4`, not a transcription of any
// published algorithm's constants.
//
// Minimax in *relative* error, not least squares: `dawsn(x)/x` falls by
// 30x across `u in [0,16]`, so an absolute-error objective spends its
// whole budget near `u=0` and leaves the top of the range 30x worse in
// the ulp that actually gets measured.
//
// `pc[0]` and `qc[0]` are both pinned to exactly `1.0`, which is what
// makes `dawson(x) == x` exact for small `|x|` -- the ratio is then
// exactly 1 there, and `x * 1.0` is exact. This is the one place where
// pinning a coefficient to its mathematically exact value pays for
// itself many times over rather than costing (cf. `erfc`): the pin is
// worth ~1.5 avg ulp across every octave below `2^-13`, which is most
// of the domain by sample count, and it costs only the fit freedom of a
// single constant term the minimax would otherwise have placed within
// its own error band of 1 anyway.
#[inline(always)]
fn dawson_central_ratio(u: f32) -> f32 {
    let pc: [f32; 6] = [
        1.0,
        -0.04486033,
        0.034414444,
        0.00092847738,
        0.00022258201,
        -4.9580024e-7,
    ];
    let qc: [f32; 6] = [1.0, 0.62182093, 0.1821982, 0.032993324, 0.0038007316, 0.00037382546];
    let u2 = u * u;
    let pl0 = fma(pc[1], u, pc[0]);
    let pl1 = fma(pc[3], u, pc[2]);
    let pl2 = fma(pc[5], u, pc[4]);
    let pr0 = fma(pl2, u2, pl1);
    let num = fma(pr0, u2, pl0);
    let ql0 = fma(qc[1], u, qc[0]);
    let ql1 = fma(qc[3], u, qc[2]);
    let ql2 = fma(qc[5], u, qc[4]);
    let qr0 = fma(ql2, u2, ql1);
    let den = fma(qr0, u2, ql0);
    num / den
}

// dawson's tail branch: `R(v)/(2x)`, `v=1/x^2`, degree 4, Estrin-grouped
// with the odd top coefficient folded in at the `v^2` level so `v^4` is
// never formed (exp_r_poly!'s own fold) -- a real least-squares fit of
// `2*x*dawsn(x)` against `v` over `|x| > 4` (`v` in `[0, 1/16]`),
// matching the `1 + v/2 + O(v^2)` asymptotic shape but fit directly
// rather than truncated from that series.
#[inline(always)]
fn dawson_tail_poly(v: f32) -> f32 {
    let c: [f32; 5] = [1.0, 0.4999613, 0.7583368, 1.4159758, 14.885198];
    let v2 = v * v;
    let l0 = fma(c[1], v, c[0]);
    let l1 = fma(c[3], v, c[2]);
    let r0 = fma(c[4], v2, l1);
    fma(r0, v2, l0)
}

/// Dawson's function `F(x) = exp(-x^2) * integral_0^x exp(t^2) dt`
/// (backlog idea #138): odd, `F(0)=0`, a single interior maximum
/// `F(x)~0.5410442` near `x~0.9241389`, decaying like `1/(2x)` for large
/// `|x|`. Two branches, same rational/asymptotic-tail shape as
/// [`erfcx`]: `x*P(u)/Q(u)` (`u=x^2`) for `|x| <= 4`, `R(v)/x` (halved
/// first, `v=1/x^2`) past it.
///
/// The tail branch halves before dividing by `x`, not after: `x` itself
/// never overflows (it's a finite input), but `2.0*x` does, for any
/// `|x|` above `f32::MAX/2` -- found by fuzzing, not assumed, since it
/// silently produces a wrong `0.0` (millions of ulp off) instead of the
/// correct tiny subnormal result there, rather than an obviously-wrong
/// `NaN`/`inf`. Halving first instead of last avoids this: halving never
/// overflows, and the following divide-by-`x` only shrinks the result
/// further.
///
/// Neither branch needs an explicit infinity override the way
/// `erfinv`'s tail does: the branch actually *returned* never sees its
/// own Estrin-overflow-to-NaN case here. The tail branch (selected for
/// huge `|x|`) computes `v=1/x^2`, which cleanly saturates to `0.0` (not
/// `NaN`) once `x*x` itself overflows to `+inf`, and `R(0)` is just its
/// own finite constant term -- no cancellation, no mixed-sign overflow
/// (`R`'s coefficients are all positive, so even `v=+inf` itself would
/// combine to a consistently-signed `+inf`, not `NaN`). The *central*
/// branch's `u=x^2` does overflow to `+inf` for that same huge `|x|`,
/// and its Estrin-grouped numerator (mixed coefficient signs) does hit
/// the same opposite-signed-infinities `NaN` erfinv's tail hit -- but
/// only in the discarded arm of the final `if`, for inputs where `tail`
/// is the one actually selected, so it never reaches the caller
/// (confirmed by edgecheck, not assumed).
#[doc(alias = "dawsn")]
#[inline(always)]
pub fn dawson(x: f32) -> f32 {
    let u = x * x;
    let central = x * dawson_central_ratio(u);
    let v = 1.0 / u;
    let tail = dawson_tail_poly(v) * 0.5 / x;
    if x.abs() <= 4.0 { central } else { tail }
}

/// logit(p) = ln(p/(1-p)), sigmoid's inverse (backlog idea #71).
///
/// The natural form, `ln(p) - log1p(-p)`, is accurate at both ends but
/// subtracts two nearly-equal logarithms near `p = 0.5`, where the true
/// result goes to zero: both terms are `~-ln(2)` there, and each carries
/// its own rounding of a quantity ~0.693 while their difference is only
/// ~`4*(p-0.5)`. That is a real relative error -- 1024 max ulp -- not
/// the "near a true zero, ulp isn't meaningful" artifact it was
/// documented as until 2026-07-31. The error grows smoothly as the
/// result shrinks (~`6e-8/ulp(result)` ulp), so there is no threshold
/// below which it is safely small; it just gets worse all the way in.
///
/// Fixed by giving the central band its own arm: `logit(p) =
/// 2*atanh(2p-1)`, evaluated with `atanh`'s own small-argument
/// polynomial over that poly's own `|x| < 0.25` domain. Nothing cancels
/// there -- the result is proportional to `2p-1`, which is *exact* for
/// any `p >= 0.25` (`2p` is an exact scaling, Sterbenz covers the
/// subtraction), and the doubling is exact too, so full relative
/// accuracy survives to the last f32 either side of 0.5.
///
/// Outside that band the difference form stays, unchanged: at the seam
/// the result is already ~0.51 against operands ~0.69, so it has
/// essentially nothing left to cancel, and it is what keeps denormal
/// `p`, the endpoints, and out-of-domain `p` correct. `p` outside
/// `[0,1]` comes out `NaN` (`ln(p)` for `p<0`, or `log1p(-p)`'s own
/// domain error for `p>1`); `logit(0)=-inf` and `logit(1)=inf` fall out
/// of `ln`/`log1p`'s existing `0`/`-1` special cases with no extra code.
///
/// A single `log1p((2p-1)/min(p,1-p))` covering the whole domain in one
/// arm was measured and rejected: better avg (0.2408 vs 0.2625) and the
/// same max, but its two divisions land in series where these two logs
/// run in parallel, for +34% latency. See graveyard.md.
#[inline(always)]
pub fn logit(p: f32) -> f32 {
    // `2p-1`, exact for every `p >= 0.25` (`2p` is an exact scaling,
    // Sterbenz covers the subtraction), which is what lets the central
    // arm keep full relative accuracy at the zero.
    let a = fma(2.0, p, -1.0);
    // logit(p) = 2*atanh(2p-1), reusing `atanh`'s own small-argument
    // polynomial on its own `|x| < 0.25` domain. No subtraction of
    // logarithms, so nothing cancels: the result is proportional to `a`,
    // and `a` is exact.
    let central = 2.0 * atanh_small(a);
    // Outside that band the difference form has nothing left to cancel
    // (at the seam the result is already ~0.51 against operands ~0.69)
    // and it is what keeps `p` denormal, `0`, `1` and out-of-domain
    // correct, so it stays exactly as it was.
    //
    // `log1p(-p)` inlined without the branches this call site can't
    // reach (same lever as `erfinv` above): `t = 1-p` is never denormal
    // (`>= 2^-24` for any `p < 1`), and `-p == 0` only at `p == 0`, where
    // `ln(p)`'s own `-inf` decides the result either way. `t == 0`
    // (`p == 1`) *is* live and must stay `-inf` — that is what makes
    // `logit(1) == +inf`. `t == +inf` (`p == -inf`) needs no arm of its
    // own: `ln(p)` is already `NaN` there.
    let np = -p;
    let t = 1.0 + np;
    let c = np - (t - 1.0);
    // No `corr.is_finite()` guard: the only `p` making it non-finite are
    // `p == 1` (`t == 0`, which takes the `spec` arm) and `p == -inf`
    // (where `ln(p)` is already `NaN`, so the subtraction below is `NaN`
    // either way).
    let corr = c / t;
    let spec = if t == 0.0 { f32::NEG_INFINITY } else { f32::NAN };
    let l = if t > 0.0 { ln_normal(t, 0.0) + corr } else { spec };
    let outer = ln(p) - l;
    if a.abs() < 0.25 { central } else { outer }
}

/// x*ln(y) (backlog idea #84), the entropy-sum kernel (`sum(p*ln(p))`
/// etc.): matches `scipy.special.xlogy`'s convention exactly -- `0` when
/// `x == 0`, *regardless of `y`* (even `y <= 0` or `y == NaN`), since
/// `0*ln(0)` is the indeterminate form entropy sums define away to `0`
/// by convention, and a caller summing many `x*ln(y)` terms wants every
/// zero-weight term to vanish from the sum without a NaN poisoning it,
/// not just the one exactly-`0*0` case. Every other `x` (including
/// `x < 0`) passes straight through to plain `x*ln(y)`, so `y <= 0`
/// still gives `+-inf`/`NaN` there exactly as bare multiplication would.
#[inline(always)]
pub fn xlogy(x: f32, y: f32) -> f32 {
    let normal = x * ln(y);
    if x == 0.0 { 0.0 } else { normal }
}

/// x*ln(1+y) (backlog idea #84), `xlogy`'s cancellation-safe sibling for
/// callers whose natural parameter is `1+y` (e.g. KL-divergence terms
/// written against a base rate) -- same `x == 0` override, and reuses
/// `log1p` instead of forming `1.0+y` and calling `xlogy`/`ln` directly,
/// avoiding the precision loss `log1p` itself exists to prevent for
/// small `y`.
#[inline(always)]
pub fn xlog1py(x: f32, y: f32) -> f32 {
    let normal = x * log1p(y);
    if x == 0.0 { 0.0 } else { normal }
}

/// (1+x)^n (backlog idea #72), the compound-interest/growth-rate
/// kernel: raising `1.0+x` to the `n`th power directly (via `powf`)
/// forms `1.0+x` as its own first step, losing exactly the low-order
/// bits of a small `x` the
/// same way a naive `exp(x)-1` loses them for `expm1` -- routing
/// through `n*log1p(x)` instead keeps `x`'s own precision intact all
/// the way through the exponent. `exp_checked` (not the unchecked
/// `exp`) since `n*log1p(x)` easily leaves the unchecked domain for
/// ordinary compounding inputs (many periods, or a large rate). Unlike
/// `powf`'s own dedicated `y==0`/`x==1` special cases, this composition
/// gives `NaN` for `compound(0.0, NaN)` (`log1p(0)=0`, `NaN*0=NaN`)
/// rather than `powf`'s C99-mandated `1.0` -- a real, deliberate
/// deviation for a thin composite, not something worth its own
/// override here.
///
/// Max ulp is in the hundreds, and that is real, not an artifact of
/// results near zero: `exp`'s error is `|n*log1p(x)| * relerr(log1p)`,
/// so a single ulp of the exponent is amplified by the exponent's own
/// magnitude, which reaches ~88 before the result overflows either way.
/// Both the `n * log1p(x)` product's own rounding and `log1p`'s ~1-2
/// ulp go through that multiplier, and both are already spent by the
/// time `exp_checked` sees them. Ordinary, fully-normal `(x, n)` pairs
/// land there -- it is the same amplification `powf` left the collapsed
/// single-f32 route to escape (see `powf_df_mag!`). [`compound_accurate`]
/// is that route here, at roughly twice the cost; this tier is the one
/// to reach for when `|n*log1p(x)|` stays modest, which for the
/// compound-interest reading of the arguments it usually does.
#[inline(always)]
pub fn compound(x: f32, n: f32) -> f32 {
    // `log1p` minus its trailing signed-zero select: that select only
    // changes `log1p(-0.0)` from `+0.0` to `-0.0`, and `exp_checked`
    // maps both zeros to exactly `1.0`, so the sign never reaches the
    // result.
    exp_checked(n * log1p_nonzero!(x))
}

// log2(1+x) as a double-float, `compound_accurate`'s exponent. Same
// shape as `log2p1`: `u = 1+x` with the bits that sum threw away
// recovered exactly as `c = x - (u-1)`, then `log2(u+c) = log2(u) +
// log2(1 + c/u)` with `|c/u| <= 2^-24`. Two differences, both forced by
// the fact that the result is about to be *multiplied* by an
// unrestricted `n` rather than returned:
//
//   - the leading term is `log2_df`, not `log_2`. Collapsing it to one
//     f32 first would throw away exactly the low bits `n` then
//     amplifies -- the same reason `powf` stopped doing that, see
//     `powf_df_mag!`.
//   - the correction has to carry its own low word too. `log2p1` can
//     round `(c/u) * LOG2_E` once because there it is a ulp-scale
//     addition onto a much bigger number. Here it is *not* small
//     compared to the leading term: whenever `1+x` rounds back to
//     exactly `1`, `log2_df(u)` is exactly `Df32(0, 0)` and the
//     correction is the entire answer, which is precisely the
//     small-`x` case `compound` exists for. So `LOG2_E` is split
//     hi/lo and the product kept as a pair.
//
// The quadratic Taylor term (`-e^2/2`) is not decoration either: at
// `u == 1` the answer is `log2(1+x)` for the full width of `x`, and
// dropping it leaves a relative `x/2` that `n` amplifies straight back
// into the tens of ulp.
//
// Degenerate `u` (`0`, negative, `+inf`, NaN) route around `log2_df`'s
// bit-level decomposition on the *high word only*, exactly as
// `powf_df_mag!` does and for the same reason: every path that reads
// the low word ends at `exp2_checked_df`'s own non-finite guard.
#[inline(always)]
fn log2p1_df(x: f32) -> Df32 {
    let u = 1.0 + x;
    let c = x - (u - 1.0);
    let l = log2_df(u);
    // e = c/u, as a pair. The residual is not optional: whenever `1+x`
    // rounds back to `1` the correction *is* the answer, so the single
    // rounding a bare `c/u` leaves gets amplified by `n` exactly like
    // the leading term's would be. One reciprocal serves both the
    // quotient and its refinement; `fma(-eh, u, c)` is exact.
    let rcp = 1.0 / u;
    let eh = c * rcp;
    let el = fma(-eh, u, c) * rcp;
    // log2(1+e) = (e - e^2/2) * log2(e) to well past double-float
    // precision at |e| <= 2^-24, as an exact two-product against the
    // split constant. `e^2` is already ~2^-48, so only its high word
    // is worth forming.
    let tl = fma(-0.5 * eh, eh, el);
    let ch = eh * LOG2_E;
    let cl = fma(eh, LOG2_E, -ch) + fma(eh, LOG2E_LO, tl * LOG2_E);
    // Full two-sum, not the quick form: `l.0` is exactly `0` for every
    // `x` small enough that `1+x` rounds to `1`, so `|l.0| >= |ch|` is
    // the one ordering that cannot be assumed here.
    let p = l.sloppy_add(ch);
    let q = Df32::from_quick_add(p.0, p.1 + cl);
    // Degenerate `u` replaces the *whole pair*, not just its high word
    // the way `powf_df_mag!` can get away with. There the low word is
    // only ever read by `exp2_checked_df`'s non-finite guard; here the
    // two-sum above would feed it back into the high word (`inf - inf`
    // in the residual, poisoning an otherwise correct `+-inf` into NaN).
    // `u == 0` is `x == -1`, exactly `0^n`; `u < 0` is `x < -1`, where
    // the real power does not exist; `u` non-finite is `x` non-finite.
    // `exp2_checked_df`'s clamp turns each into the right saturation.
    let deg = if u == 0.0 { f32::NEG_INFINITY } else { u };
    let deg = if u < 0.0 { f32::NAN } else { deg };
    if u > 0.0 && u < f32::INFINITY { q } else { Df32(deg, 0.0) }
}

/// (1+x)^n, accurate tier (see [`compound`] for the construction and
/// for why `1+x` is never formed in f32). [`compound`]'s exponent
/// `n * log1p(x)` is a single f32, and `exp`'s error is
/// `|n*log1p(x)| * relerr(log1p)`: one ulp in that product is amplified
/// by the exponent's own magnitude, which reaches ~88 before the result
/// overflows, so ordinary inputs land in the hundreds of ulp -- not the
/// near-zero-underflow artifact its doc comment used to claim (over 20M
/// samples, 5435 of the 5440 worst had a perfectly normal result).
/// This tier keeps the exponent a double-float from `log2p1_df` all the
/// way through the multiply by `n` and the `exp2` reconstruction,
/// collapsing to one f32 only at the very end -- the same fix, and the
/// same machinery, `powf` uses. Costs roughly twice [`compound`]; both
/// are kept because that is a real trade, not a free win.
#[inline(always)]
pub fn compound_accurate(x: f32, n: f32) -> f32 {
    exp2_checked_df(log2p1_df(x) * n)
}

/// erfcx(x) = e^(x^2)*erfc(x), the "scaled complementary error
/// function". For `x >= 0` this is [`erfcx_pos`] alone, with *no
/// exponential at all*: `erfc`'s own construction is
/// `exp(-x^2)*erfcx_pos(x)`, so multiplying by `exp(x^2)` cancels the
/// exponential exactly (algebraically, not numerically -- the exponent
/// is never computed, let alone cancelled). This is exactly what erfcx
/// is *for*: the naive `exp(x*x)*erfc(x)` a caller might otherwise
/// write breaks down before this function's domain gets interesting --
/// `exp(x*x)` alone overflows f32 for `|x| >~ 9.3`, while `erfcx`'s
/// true value there is still a small, well-behaved number
/// (`erfcx(x) ~ 1/(x*sqrt(pi))`).
///
/// Accurate over the *whole* positive half-line, not just `|x| <= 10`:
/// the reciprocal variable `1/(2+x)` reproduces that `1/(x*sqrt(pi))`
/// asymptote by construction and cannot overflow, so there is no domain
/// clamp to freeze against (see [`erfcx_pos`]). The rational this
/// replaced did freeze past `x = 10`, with relative error growing
/// without bound (~10% at 11, ~99% at 20, ~895% at 100); that is fixed,
/// and the separate `erfcx_checked` tier that used to paper over it with
/// an asymptotic branch is gone with it.
///
/// For x < 0, uses erfc's own reflection identity (`erfc(x) = 2 -
/// erfc(-x)` for x<0) to derive `erfcx(x) = 2*exp(x^2) - erfcx(-x)` --
/// unlike the x>=0 branch this does need one real exponential, because
/// `erfcx` genuinely diverges to `+inf` for sufficiently negative x
/// (`erfcx(-10) ~ 2*e^100`, far past `f32::MAX`): that is this
/// function's true mathematical behavior, not an implementation gap,
/// and `exp_checked`'s saturation makes it come out `+inf` correctly
/// rather than wrapping to garbage. `x*x` is split exactly the same way
/// [`erfc`] splits it (`exp(p+pe) ~ exp(p)*(1+pe)`, folded into the
/// already-needed `2*` scale as `2+2*pe`) -- that split, not the
/// polynomial, is what used to hold this branch at max ulp 126. With it
/// in place the negative branch is no longer where the worst case lives:
/// that moved to the positive branch near zero, where it is
/// [`erfcx_pos`]'s own `v`-formation floor. Written as a multiply rather
/// than `fma(g, 2.0*pe, ...)` on purpose: `g` saturates to `+inf`, and
/// `2+2*pe` is always positive, so `inf*(2+2*pe) = +inf` where a signed
/// `pe` in the multiplicand could have given `inf - inf`. The branch
/// diverges to `+inf` for `x < -9.382` (`x^2 > ~88.03`, where `2*e^(x^2)`
/// leaves f32); `exp_checked`'s saturation makes that come out `+inf`
/// rather than wrapping, and both sides of that boundary are pinned in
/// edgecheck.
///
/// Current: max ulp 6, avg 0.2142 (exhaustive sweep over `|x| <= 10`).
/// Was 126 / 0.3768. The rest of the domain is measured too, and this
/// is the part that used to not exist: `|x| <= 20` is 6 / 0.2154, and
/// `x >= 20` out to `f32::MAX` is **4 / 0.6478** over 1.04e9 samples,
/// scored against the asymptotic series (`exp(x^2)` overflows f64 past
/// x ~ 26.6, so the composed reference cannot reach there).
#[inline(always)]
pub fn erfcx(x: f32) -> f32 {
    let xa = x.abs();
    let r = erfcx_pos(xa);
    let xs = if xa > ERFCX_XS_CLAMP { ERFCX_XS_CLAMP } else { xa };
    let p = xs * xs;
    let pe = fma(xs, xs, -p);
    let g = exp_reduce!(p) * fma(pe, 2.0, 2.0);
    if x >= 0.0 { r } else { g - r }
}

/// 1/sqrt(x). Unlike most functions in this crate, no bit-trick seed or
/// fitted correction poly needed: `sqrt` and division are each already
/// correctly-rounded IEEE754 hardware operations (`x.sqrt()` isn't a
/// software approximation), so composing them directly costs at most
/// ~1 ulp (one rounding from each op) with zero special-casing --
/// `x <= 0.0` (including `-0.0`), `x.is_nan()`, and `x == inf` all
/// already give the right answer (`+inf`/`inf`, `NaN`, `0.0`
/// respectively) purely from IEEE754 semantics. (The crate's older
/// [`rsqrt_approx`] -- a Quake-style bit-trick seed with no correction,
/// ~4.8% max relative error -- is a deliberately-rough exploratory
/// function kept for the `_approx_plot` test suite, not a candidate
/// replacement.)
#[inline(always)]
pub fn rsqrt(x: f32) -> f32 {
    1.0 / x.sqrt()
}

/// hypot without the +-inf special case: contract is x, y both finite (or
/// both NaN-safe, since NaN propagates through fma/sqrt on its own) --
/// see hypot's own doc comment for the one case this drops (+-inf paired
/// with a NaN, where IEEE754/C99 defines +inf as the answer regardless).
/// Bit-identical to hypot whenever neither argument is infinite.
#[inline(always)]
pub fn hypot_unchecked(x: f32, y: f32) -> f32 {
    fma(x, x, y * y).sqrt()
}

/// Straight port of jodiemath's hypotf: naive sqrt(x^2+y^2), no anti-overflow
/// rescaling (unlike std's hypot) -- trades the overflow/underflow edge cases
/// for vectorizability, same tradeoff this crate makes for cbrt/sin/cos vs.
/// their std counterparts.
#[doc(alias = "hypotf")]
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

/// hypot with anti-overflow/underflow rescaling: extracts the larger
/// argument's exponent via bit tricks, rescales both arguments by an
/// exact power of two before squaring (so `x*x+y*y` can never overflow,
/// and never underflows the *dominant* term -- if the smaller term
/// underflows to exactly 0 after scaling, its true contribution was
/// already negligible at f32 precision, so this loses nothing real),
/// then scales the sqrt'd result back. Unlike the plain `hypot` above,
/// this recovers std-grade overflow/underflow behavior while staying
/// fully branchless (selects only, no early returns, so array loops
/// still auto-vectorize).
///
/// The scale exponent is rounded down to the nearest *even* value
/// (`es = 2*floor(e/2)`, `e >> 1` is Rust's arithmetic shift, i.e. floor
/// division for negative `e` too) rather than using `e` directly: `e`
/// alone can reach up to +-127, and 2^-127 has no *normal* single-word
/// representation (it's already past the denormal boundary), so building
/// it via this exponent-field bit trick would corrupt the value instead
/// of producing the reciprocal scale. Rounding to an even `es` keeps the
/// scale factor's own exponent safely within the representable normal
/// range for every valid input `m`, at the cost of the rescaled
/// magnitude landing in `[1,4)` instead of a tighter `[1,2)` -- still
/// small enough that squaring and summing never overflows.
///
/// `is_zero` is deliberately checked on `ax`/`ay` directly (`x.abs()`/
/// `y.abs()` are exactly zero), not on `m = ax.max(ay)`: `f32::max`
/// silently returns the *non-NaN* operand when only one argument is
/// NaN, so `hypot_checked(f32::NAN, 0.0)` would otherwise compute
/// `m = 0.0` and wrongly take the "both zero" branch, discarding NaN.
/// With the narrower check, `hypot_checked(NaN, 0.0)` still resolves to
/// exactly `NaN`: the exponent extraction degrades to a
/// garbage-but-finite scale either way, but `ax` (or `ay`) being NaN
/// itself propagates through the unconditional multiply regardless.
/// Max ulp 1 (sampled across the full exponent range, plus every
/// zero/NaN/inf combination directly).
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
    if x.is_infinite() || y.is_infinite() { f32::INFINITY } else { normal }
}

/// 1/hypot(x,y): normalizing a 2D vector (`(x,y) / hypot(x,y)`) is
/// hypot's single most common real use case, and computing the
/// reciprocal directly saves the caller their own separate division.
/// Same "compose already-correctly-rounded hardware ops" reasoning as
/// `rsqrt` -- `hypot_unchecked`'s `fma(x,x,y*y)` core plus one sqrt and
/// one division. The naive composition gets every zero/inf/nan special
/// case right *except one*, purely from IEEE754 semantics:
/// `rhypot(0,0)=inf`, `rhypot(x,inf)=0`, `rhypot(NaN,y)=NaN` all fall
/// out for free. The exception mirrors `hypot`'s own special case:
/// `+-inf` paired with `NaN` degrades to `1/sqrt(NaN)=NaN`, but
/// IEEE754/C99 defines `hypot(+-inf, NaN) = +inf` (infinity "wins"), so
/// the reciprocal should be `0` -- same override `hypot` uses. No
/// anti-overflow rescaling tier: `x*x+y*y` overflowing is already an
/// accepted tradeoff of the crate's default `hypot`.
#[inline(always)]
pub fn rhypot(x: f32, y: f32) -> f32 {
    let normal = 1.0 / fma(x, x, y * y).sqrt();
    if x.is_infinite() || y.is_infinite() { 0.0 } else { normal }
}

/// 3-arg Euclidean norm, `sqrt(x^2+y^2+z^2)` (backlog idea #55): a
/// graphics/physics staple (vector magnitude), same naive-fma-chain
/// tradeoff as [`hypot`] -- no anti-overflow rescaling, and the same
/// IEEE754/C99 "infinity wins over NaN" override for any argument being
/// `+-inf`.
#[inline(always)]
pub fn hypot3(x: f32, y: f32, z: f32) -> f32 {
    let normal = fma(x, x, fma(y, y, z * z)).sqrt();
    if x.is_infinite() || y.is_infinite() || z.is_infinite() { f32::INFINITY } else { normal }
}

/// Reciprocal of [`hypot3`], `1/sqrt(x^2+y^2+z^2)` -- the normalize-a-
/// vector building block ([`hypot3`]'s own point, per its doc comment).
/// Same `+-inf` override as [`rhypot`] and for the same reason: if one
/// argument is `+-inf` and another is `NaN`, the naive chain degrades to
/// `1/sqrt(inf + NaN) = 1/NaN = NaN` instead of the IEEE754/C99-defined
/// `0` (infinity "wins" over NaN in `hypot`'s own combine, so its
/// reciprocal should too).
#[inline(always)]
pub fn rnorm3(x: f32, y: f32, z: f32) -> f32 {
    let normal = 1.0 / fma(x, x, fma(y, y, z * z)).sqrt();
    if x.is_infinite() || y.is_infinite() || z.is_infinite() { 0.0 } else { normal }
}

/// 2D vector normalize (backlog idea #136, the operation users actually
/// want `rhypot` for): scales `(x,y)` by [`rhypot`] to unit magnitude.
/// Same zero-vector (`(NaN,NaN)`, no meaningful direction) and `+-inf`
/// override behavior as [`normalize4`], one level down.
#[inline(always)]
pub fn normalize2(x: f32, y: f32) -> (f32, f32) {
    let r = rhypot(x, y);
    (x * r, y * r)
}

/// 3D vector normalize (backlog idea #136): scales `(x,y,z)` by
/// [`rnorm3`] to unit magnitude. See [`normalize2`]/[`normalize4`] for
/// the shared zero-vector/`+-inf` behavior.
#[inline(always)]
pub fn normalize3(x: f32, y: f32, z: f32) -> (f32, f32, f32) {
    let r = rnorm3(x, y, z);
    (x * r, y * r, z * r)
}

/// 4-arg Euclidean norm (backlog idea #134, companion to [`hypot3`]):
/// same naive-fma-chain construction, tradeoff, and inf/NaN override,
/// one argument wider -- the quaternion-magnitude case.
#[inline(always)]
pub fn hypot4(w: f32, x: f32, y: f32, z: f32) -> f32 {
    let normal = fma(w, w, fma(x, x, fma(y, y, z * z))).sqrt();
    let any_inf = w.is_infinite() || x.is_infinite() || y.is_infinite() || z.is_infinite();
    if any_inf { f32::INFINITY } else { normal }
}

/// Reciprocal of [`hypot4`] -- see [`rnorm3`]'s own doc comment for why
/// the `+-inf` override is needed (same reason, one argument wider).
#[inline(always)]
pub fn rnorm4(w: f32, x: f32, y: f32, z: f32) -> f32 {
    let normal = 1.0 / fma(w, w, fma(x, x, fma(y, y, z * z))).sqrt();
    let any_inf = w.is_infinite() || x.is_infinite() || y.is_infinite() || z.is_infinite();
    if any_inf { 0.0 } else { normal }
}

/// Quaternion normalize (backlog idea #134): scales `(w,x,y,z)` by
/// [`rnorm4`] so the result has unit magnitude -- the actual operation
/// callers reach for `rnorm4` to build themselves, provided directly.
/// Inherits `rnorm4`'s own `+-inf`-in-any-component override (giving
/// every component `0` rather than `NaN`) and its zero-vector behavior
/// (`rnorm4(0,0,0,0)` is `+inf`, so `normalize4(0,0,0,0)` is
/// `(NaN,NaN,NaN,NaN)` via `0*inf` -- there's no meaningful unit
/// quaternion for the zero vector, so propagating `NaN` rather than
/// picking an arbitrary direction is the honest answer).
#[inline(always)]
pub fn normalize4(w: f32, x: f32, y: f32, z: f32) -> (f32, f32, f32, f32) {
    let r = rnorm4(w, x, y, z);
    (w * r, x * r, y * r, z * r)
}

/// `a*b - c*d`, computed via Kahan's compensated algorithm instead of the
/// naive two-multiply-one-subtract form (backlog idea #135): the naive
/// form's error is unbounded relative to the true result whenever `a*b`
/// and `c*d` are close in magnitude, since each product's own independent
/// rounding error survives the subtraction untouched. `w = c*d` (rounded
/// once), then `e = fma(-c,d,w)` recovers that rounding's exact error term
/// (a two-product, same construction `Df32::from_mul` uses elsewhere in
/// this crate, just inlined rather than returning a pair), and `f =
/// fma(a,b,-w)` folds `a*b`'s own rounding against the same `w` in one
/// step -- `f + e` then combines both error corrections, accurate to
/// within a couple ulp of the true value even at total cancellation, where
/// the naive form has no error bound there at all.
///
/// Same overflow tradeoff as `hypot`/`rhypot` above, one level removed:
/// if `a*b` or `c*d` individually overflows to `+-inf` in f32 even though
/// the true difference is finite, `w` becomes infinite and the correction
/// terms can degrade to `inf - inf = NaN` instead of the finite answer.
/// Only reachable once `|a*b|` or `|c*d|` approaches `f32::MAX`; well
/// inside that range (both products individually representable) this is
/// unaffected.
#[inline(always)]
pub fn diff_of_products(a: f32, b: f32, c: f32, d: f32) -> f32 {
    let w = c * d;
    let e = fma(-c, d, w);
    let f = fma(a, b, -w);
    f + e
}

/// 2D cross product (scalar "determinant" form, `ax*by - ay*bx`) via
/// [`diff_of_products`]: the standard signed-area/orientation predicate
/// in graphics and computational geometry, exactly the shape
/// [`diff_of_products`] exists to make accurate near cancellation (two
/// nearly-parallel or nearly-antiparallel vectors, where the true cross
/// product is small but each product term individually isn't).
#[inline(always)]
pub fn cross2(ax: f32, ay: f32, bx: f32, by: f32) -> f32 {
    diff_of_products(ax, by, ay, bx)
}

/// Complex modulus `|re+im*i|` (backlog idea #186): a named alias for
/// [`hypot_checked`], not the faster but overflow-prone [`hypot`] --
/// found by fuzzing, not assumed: `hypot`'s own doc comment already
/// documents "naive sqrt(x^2+y^2), no anti-overflow rescaling" as a
/// deliberate tradeoff for that function specifically, but a general-
/// purpose complex modulus silently going to `inf` for any two inputs
/// individually representable in f32 (found via `clog`'s own
/// composition: `hypot(-1413.68, 1.8447e19)` gives `inf`, not `≈im`) is
/// a worse default here, where nothing calls for `hypot`'s speed at the
/// cost of that gap.
#[doc(alias = "cabsf")]
#[inline(always)]
pub fn cabs(re: f32, im: f32) -> f32 {
    hypot_checked(re, im)
}

/// Complex argument (principal value, backlog idea #186): a named
/// alias for [`atan2`], same relationship to [`cabs`]/[`hypot`] as
/// `carg` has to `cabs` mathematically (`atan2`'s `(y, x)` argument
/// order already matches `carg`'s own `(im, re)`).
#[doc(alias = "cargf")]
#[inline(always)]
pub fn carg(re: f32, im: f32) -> f32 {
    atan2(im, re)
}

/// Complex exponential `e^(re+im*i) = e^re*(cos(im)+i*sin(im))`
/// (backlog idea #186), returned as `(re, im)`. Calls `cos`/`sin`
/// separately rather than a hand-fused "sincos": both are
/// `#[inline(always)]` over the same reduction, and this crate's own
/// `sinh_cosh` investigation (see IDEAS.md) already confirmed LLVM's
/// GVN shares that reduction across separate call sites for free, so a
/// combined function would cost real complexity for no real speedup.
#[inline(always)]
pub fn cexp(re: f32, im: f32) -> (f32, f32) {
    let m = exp(re);
    (m * cos(im), m * sin(im))
}

/// Complex natural log `ln(re+im*i) = ln(|re+im*i|) + i*arg(re+im*i)`
/// (backlog idea #186), returned as `(re, im)`. The imaginary part is
/// just [`carg`] (`atan2`'s principal-value convention is exactly what
/// the complex log's imaginary part is defined to be); the real part
/// needs more care than a direct `ln(cabs(re, im))`, for two distinct
/// reasons found by fuzzing, not assumed:
///
/// - `|z|` near `1` (e.g. `re=1.0000999, im=0.00242`) makes `ln(|z|)`
///   itself the near-zero-argument cancellation [`log1p`] exists to
///   avoid -- `ln` of a value that close to its own zero amplifies even
///   `cabs`'s own already-good ~1-ulp error into thousands of ulp in
///   the tiny output. `log1p(cabs-1.0)` recovers most of this (the
///   *only* fix available at this precision tier: the residual
///   thousands-of-ulp error left even after it is `cabs`'s own ~1-ulp
///   imprecision at magnitude ~1 having nowhere near enough absolute
///   precision left to support a correctly-rounded *much smaller*
///   output -- an information-theoretic floor of composing on top of
///   plain `f32` `cabs`, not something reachable without a compensated/
///   double-float magnitude the way `_accurate` tiers elsewhere in this
///   crate use, out of scope for this composite). `log1p(cabs-1.0)`
///   fixes what it can *only* in the narrow band `cabs` is actually
///   close to `1.0` (`|cabs-1.0| < 0.5` here) -- reaching for it
///   unconditionally was tried first and made things far worse (max
///   ulp over a billion)
///   for `cabs` far from `1`, e.g. very small: `log1p`'s own argument is
///   then close to `-1`, not `0`, none of the cancellation-avoidance
///   `log1p` provides actually applies there, and round-tripping an
///   already-computed `cabs` through `1.0+(cabs-1.0)` internally just
///   reintroduces the exact cancellation this was meant to avoid.
/// - `cabs` can overflow to `inf` even when both `re`/`im` are finite
///   (e.g. both individually near `f32::MAX`): the true mathematical
///   magnitude exceeds `f32::MAX` before `ln` of it would, so forming
///   `cabs` first throws away a real, still-representable answer.
///   Rescued the same way `hypot_checked` itself avoids overflow
///   internally: factor out the larger magnitude `mx` before squaring
///   (`ratio = mn/mx` stays in `[0,1]`, never overflows), so
///   `ln(cabs(re,im)) = ln(mx) + 0.5*log1p(ratio*ratio)` without ever
///   forming the too-large intermediate (`mx` is never close to `1`
///   here -- `cabs` wouldn't have overflowed if it were -- so this
///   branch uses plain `ln(mx)`, not `log1p`). Only taken when both
///   inputs are finite but `cabs` isn't -- a genuinely infinite/NaN
///   `re`/`im` instead falls through to plain `ln(cabs(re,im))`,
///   inheriting whatever convention `cabs`/`hypot_checked` already
///   establish there (e.g. infinity-wins-over-NaN) rather than
///   re-deriving it.
///
/// Neither `log1p` call is a real `log1p` call: both sit under a guard
/// that leaves `log1p` nothing to guard against ("the guard is the
/// licence", as in `atanh`/`erfinv`). `|mag-1| < 0.5` puts the first
/// one's `u = 1+v` in `(0.5, 1.5)` and `ratio in [0,1]` puts the
/// second's in `[1, 2]`, so on top of `ln`'s zero/negative/denormal/
/// non-finite arms, *both* of `log1p`'s own selects are dead too: `c/u`
/// cannot be non-finite over those ranges, and the `v == 0.0`
/// signed-zero guard has nothing to fix (`v` is `+0.0` at worst, and
/// `ln_normal(1.0, 0.0) + 0.0` is already `+0.0`). What is left is the
/// Sterbenz correction and the bare poly.
#[inline(always)]
pub fn clog(re: f32, im: f32) -> (f32, f32) {
    // `log1p(v)` for a `v` that is known to keep `1+v` positive, normal
    // and finite -- see this function's doc comment. Bit-identical to
    // `log1p` over every f32 in both call sites' licensed ranges
    // (`|v| < 0.5` and `[0, 1]`, checked exhaustively) with exactly one
    // exception, `v == -0.0`, where `log1p`'s signed-zero select returns
    // `-0.0` and this returns `+0.0`. Neither site can produce it, also
    // checked exhaustively rather than argued: `mag - 1.0` is `+0.0` for
    // every non-negative finite `mag` (IEEE `x - x` is `+0.0` under
    // round-to-nearest), and a square is never `-0.0`.
    #[inline(always)]
    fn log1p_guarded(v: f32) -> f32 {
        let u = 1.0 + v;
        let c = v - (u - 1.0);
        ln_normal(u, 0.0) + c / u
    }
    let mag = cabs(re, im);
    let log_mag = if mag.is_finite() {
        if (mag - 1.0).abs() < 0.5 { log1p_guarded(mag - 1.0) } else { ln(mag) }
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

/// `2*log2(e)` split into a double-float: the atanh form's leading
/// coefficient (see `log2_df`), which has to carry more than f32's own
/// 24 bits since it multiplies the *whole* result.
const LOG2E_2_HI: f32 = 2.885390043258667;
const LOG2E_2_LO: f32 = 3.851926067000022e-8;

/// Minimax fit of `1/(m+1)` over `m` in `[2^-0.5, 2^0.5]` -- only a seed:
/// one Newton step squares its error and the quotient refinement squares
/// it again, so this is nowhere near a 24-bit fit and does not need to
/// be. See `log2_df` for why it is degree 3 anyway.
///
/// In `m` rather than in `d = m + 1` (the same fit either way, since the
/// two differ by an affine change of variable) so that the seed does not
/// have to wait on the `m + 1` add: this poly is the head of the whole
/// function's dependency chain, and `d` is still wanted later, where it
/// is off the critical path.
const LOG2_ATANH_RCP: [f32; 4] = [
    0.9356114864349365,
    -0.6780579090118408,
    0.29953840374946594,
    -0.057135023176670074,
];

/// `2*log2(e) * (atanh(t) - t) / t^3` in `u = t^2`, over `u` in
/// `[0, (3-2*sqrt(2))^2]` -- the atanh series past its own leading term
/// (see `log2_df`).
const LOG2_ATANH_G: [f32; 4] = [
    0.9617967009544373,
    0.5770797729492188,
    0.4119008779525757,
    0.3365967571735382,
];

/// log2(x) as a double-float (Df32) instead of a collapsed f32, for
/// positive finite x only (same domain log_2_normal assumes -- callers
/// must guard zero/negative/inf/nan themselves). Denormal input is
/// handled the same way log_2's wrapper does (scale up, offset k).
///
/// The point of this function is *relative* accuracy well past f32's own
/// 24 bits, not just a spare low word: its only callers multiply the
/// result by `y` before exponentiating, so `powf`'s final
/// relative error is about `ln2 * |y*log2(x)| *` this function's own
/// relative error -- amplified by up to 128 (the largest `|y*log2(x)|`
/// with a finite result), which is ~7 bits. Landing within an ulp
/// therefore needs ~31 bits here.
///
/// That rules out `log_2_normal`'s shape (`k + s*P(s)` with `s = m - 1`),
/// however carefully its two roundings are captured as a low word: `P(s)`
/// is itself evaluated in f32, and for `x` near 1 -- exactly where `|y|`
/// is large enough to matter -- `k` is 0 and the result *is* `s*P(s)`, so
/// `P`'s own evaluation rounding lands directly on the result as a ~2^-24
/// relative error that no amount of low-word bookkeeping downstream can
/// recover. `|s|` reaches 0.414, so peeling terms off `P` doesn't fix it
/// either: each peeled term only buys ~2.3 bits, and 3+ of them would
/// each need their own exact double-float accumulation.
///
/// The atanh form reduces the argument instead. With `t = (m-1)/(m+1)`,
/// `log2(m) = 2*log2(e) * atanh(t)` and `|t| <= 3-2*sqrt(2) ~ 0.1716`, so
/// the series past the leading `t` term is only ~1% of the result: one
/// exactly-accumulated term (`2*log2(e) * t`, a two-product against a
/// double-float coefficient) plus a *plain f32* correction for the rest
/// clears 31 bits, where the `s*P(s)` form needs four. Everything the
/// correction's own rounding costs is scaled by that ~1%.
///
/// `t` has to be built to double-float accuracy without a division, and
/// is, in two squaring steps: `LOG2_ATANH_RCP` seeds `1/d` to ~2^-13 over
/// the narrow range `d = m + 1` covers, one Newton step squares that to
/// f32's own floor, and the usual quotient-refinement step (`tl`, from
/// the residual `s - th*d`) squares it again, past 2^-38.
///
/// That refinement needs the residual against the *exact* `d`, which
/// `m + 1` is not -- it drops exactly one bit of `m` over this range. The
/// way to get it is not to reconstruct the dropped bit but to never form
/// `d` at all: `d == s + 2` exactly as a real number (`s = m - 1` is
/// Sterbenz-exact), so
///
/// ```text
///     s - th*d  ==  s - th*(s + 2)  ==  (s - 2*th) - th*s
/// ```
///
/// and `s - 2*th` is *itself* Sterbenz-exact, because `s/(2*th)` is
/// `d/2` to within the seed's own error and `d/2` lies in
/// `[0.854, 1.207]`. So one fma produces it exactly and a second one
/// rounds the whole residual once, at the residual's own tiny scale. Two
/// ops, against four for carrying a `(dh, dl)` split of `d` around, and
/// `d` itself survives only as the Newton step's operand, where an
/// approximate value is all that was ever needed.
///
/// The seed is degree 3 rather than the degree 2 that squaring alone
/// would justify, because `tl`'s *size* matters as well as `t`'s
/// accuracy: `tl` sets how far the returned low word sits below the high
/// one, and `exp2_checked_df` folds that low word in as a first-order
/// term, so a low word one part in 2^-19 rather than ~1 ulp leaves a
/// visible second-order remainder. Cheaper than renormalizing after the
/// fact, which is the other way to fix it.
///
/// The correction term cubes its argument, which *triples* its relative
/// error, so it uses the refined `th + tl` and not `th` -- `th` alone
/// carries the seed-plus-Newton error, and 3x that against the term's ~1%
/// weight would dominate everything else in this function.
///
/// `x == 1` gives exactly `Df32(0, 0)`, which `powf` relies on.
#[inline(always)]
fn log2_df(x: f32) -> Df32 {
    let (xs, koff) = denormal_rescale!(x);
    // Same decomposition log_family_normal! does, spelled out (like
    // ln_normal/log10_normal's own copies) since nothing else here is
    // shared with the s = m - 1 poly form: m in [2^-0.5, 2^0.5), k exact.
    let e = (xs.to_bits() as i32).wrapping_sub(0x3f3504f3) >> 23;
    let m = f32::from_bits((xs.to_bits() as i32).wrapping_sub(e << 23) as u32);
    let k = e as f32 + koff;
    let s = m - 1.0;
    // 1/(m+1): minimax seed, one Newton step (r*(2 - d*r)). Both this and
    // the `tl` refinement below square the error they are given. Estrin,
    // not Horner: this sits at the head of the chain everything else
    // waits on, so a level of depth here is worth an extra multiply.
    // `d` is only ever the Newton step's operand, so the one bit `m + 1`
    // drops here costs nothing -- see this fn's doc comment.
    let rc = LOG2_ATANH_RCP;
    let d = 2.0 + s;
    let m2 = m * m;
    let r = fma(fma(rc[3], m, rc[2]), m2, fma(rc[1], m, rc[0]));
    let rcp = fma(-d, r, 2.0) * r;
    // t = s/d as (th, tl): th, then the residual s - th*d scaled back.
    // `fma(-2.0, th, s)` is exact and `fma(-th, s, ...)` rounds the whole
    // residual once, against the exact `s + 2` rather than the rounded
    // `m + 1`.
    let th = s * rcp;
    let rh = fma(-th, s, fma(-2.0, th, s));
    let tl = rh * rcp;
    // 2*log2(e) * t, exactly: two-product on the hi coefficient, both
    // remaining cross terms folded into the low word (tl*LO is ~2^-49
    // relative and is dropped).
    let hi = th * LOG2E_2_HI;
    let lo = fma(th, LOG2E_2_HI, -hi);
    let lo = fma(th, LOG2E_2_LO, lo);
    let lo = fma(tl, LOG2E_2_HI, lo);
    // 2*log2(e) * (atanh(t) - t) = t^3 * G(t^2), the ~1% correction --
    // built to run *beside* the `tl` chain above rather than after it,
    // which is where this function's latency used to go. Both halves
    // start from `th` alone:
    //   - G's argument may be sloppy. G varies by only ~2% across the
    //     whole range of t^2, so feeding it `th^2` instead of the refined
    //     t^2 moves it by ~2% of th's own error -- nothing. Only the t^3
    //     prefactor needs the refinement.
    //   - and it gets it without waiting for a squaring: t^3 is
    //     th^2 * (th + 3*tl) to within 3*th*tl^2, which is ~2^-44
    //     relative. So `uh` is ready early and only the one fma below is
    //     downstream of `tl`.
    let uh = th * th;
    let g = LOG2_ATANH_G;
    // Horner, unlike the seed poly above: off the critical path now (the
    // `tl` chain beside it is longer), so the extra multiply an Estrin
    // split would cost buys nothing back. Re-tested: an Estrin split plus
    // regrouping to `(uh*gp) * t3` does take a level off the chain and
    // buys 2.5% latency, but costs 2.2% throughput and 3 uOps -- the wrong
    // side of this crate's throughput-first trade.
    let gp = fma(fma(fma(g[3], uh, g[2]), uh, g[1]), uh, g[0]);
    let corr = (uh * fma(3.0, tl, th)) * gp;
    // Both sums are exact, so they can go in either order -- k first,
    // because `hi` is ready long before `corr` is and this puts that
    // two-sum beside the correction instead of after it. Both keep the
    // quick form's precondition: |k| >= 1 > |hi| whenever k != 0 (and
    // k == 0 makes its two-sum exact on its own), then |k + hi| is either
    // >= 0.5 or exactly `hi`, and `corr` is ~1% of `hi` in both cases.
    let p = Df32::from_quick_add(k, hi);
    let q = Df32::from_quick_add(p.0, corr);
    // corr is far larger than `lo` (~1% of hi, vs ~1 ulp), which is why
    // it joins as its own exact two-sum rather than being added into
    // `lo`: that would round the pair at corr's scale and cost ~7 bits.
    // Both residuals are ulp-scale by construction, so this last add is
    // lossless to well past what the result needs.
    Df32(q.0, q.1 + (p.1 + lo))
}

/// exp2 of a double-float argument, reusing exp2_checked's own clamp/
/// k1-k2-split/poly machinery on the hi component, then folding the lo
/// component in as a multiplicative correction:
/// `exp2(hi + lo) = exp2(hi) * exp2(lo) ~= exp2(hi) * (1 + lo*ln2)`
/// (first-order Taylor, valid since `lo` is always tiny relative to 1 by
/// Df32's own invariant). This is the step that actually captures the
/// precision `log2_df` preserved: without it, the lo component would just
/// be silently dropped and this would be no more accurate than the plain
/// `exp2_checked(v.to_f32())`.
///
/// The correction goes onto `p` -- the pre-`t2` partial product -- rather
/// than onto the finished result, which is free and strictly better on
/// both axes. `p` is `t1 * 2^f`, so it is always finite and *normal* by
/// construction (`t1` is bounded well inside the exponent field and
/// `2^f` is in `[1,2)`), where the finished result can be `0`, `+-inf`,
/// or denormal. That means: a denormal output now rounds exactly once, in
/// the final `* t2`, instead of having the correction applied *after* the
/// mantissa bits were already lost; and the non-finite guard below only
/// has to cover a genuinely non-finite `v.1`, not the ordinary overflow
/// case, which used to reach it as an `inf*negative + inf` NaN.
#[inline(always)]
#[allow(clippy::approx_constant)] // g0's constant term is a fitted minimax
// coefficient near ln(2), not ln(2) itself (bit pattern deliberately differs)
fn exp2_checked_df(v: Df32) -> f32 {
    // clamp *before* floor/subtract (matching exp2_checked exactly): if
    // v.0 itself is already +-inf (y large enough that y*log2(x)
    // overflows in the Df32 multiply), flooring first and clamping after
    // would compute `inf - floor(inf)` = NaN instead of correctly
    // saturating. (The round-domain refit was tried and rejected here
    // too -- accuracy win but real throughput cost, see IDEAS.md
    // §exp/exp2.)
    let xs = v.0.clamp(-151.0, 128.0);
    let k = xs.floor();
    let f = xs - k;
    let (t1, t2) = exp2_field_split(k);
    let q = exp2_q_poly!(f);
    let p = fma(q, t1 * f, t1);
    // Drop the correction whenever the clamp fired. `xs != v.0` is exactly
    // "the answer is saturated": above `128` every reachable `v.1` still
    // leaves `2^v` past `f32::MAX`, and below `-151` past the point where
    // anything rounds to a nonzero denormal, so there is no precision left
    // to refine either way -- and applying it anyway is *wrong*, not just
    // wasted. It was: `powf(f32::MAX, 1.0000001)` clamps `128.0000153` down
    // to `128`, discarding a factor of `1.0000106` that the correction (a
    // *negative* ~6e-8 here) then can't see, and the result comes back
    // `f32::MAX` where it must overflow.
    //
    // The same compare covers the non-finite cases for free, which is why
    // it replaces the `is_finite` guard the sibling `log1p`/`atanh`/`acosh`
    // functions carry: `v.1` is non-finite exactly when `v.0` was (an
    // `inf - inf` out of the Df32 multiply's residual step, reached when
    // `y*log2(x)` overflows on its own or `log2(x)` is `+-inf` for a zero
    // or infinite base), and every one of those has `v.0` outside the clamp
    // too. NaN takes the same path (`NaN == NaN` is false) and still
    // propagates, through `f`.
    //
    // Guarding the *correction* rather than the corrected value also keeps
    // it off the critical path -- `xs` and `v.0` are both ready long before
    // the poly is, where guarding `fma(p, c, p)` afterwards put the compare
    // and its select on the end of the chain. Zeroing `c` is exactly
    // equivalent, since `fma(p, 0.0, p)` is `p`.
    let c = if xs == v.0 { v.1 * LN_2 } else { 0.0 };
    fma(p, c, p) * t2
}

/// `exp2(log2(ax) * y)`, the magnitude half of the whole `powf` family.
///
/// `log2(ax)` stays a double-float (`Df32`) through the multiply by `y`
/// and the `exp2` reconstruction, collapsing to a single f32 only at the
/// very end (see `log2_df`/`exp2_checked_df`). That is the formula, not
/// an opt-in accuracy tier, because the single-f32 route has no way to be
/// merely *approximate* here: `exp2(log2(x)*y)`'s error is
/// `ln2 * |y*log2(x)| * relerr(log2)`, and collapsing `log2(x)` to one
/// f32 *before* the multiply throws away exactly the low bits `y` then
/// amplifies -- by up to 128, i.e. hundreds of ulp, at ordinary inputs
/// like `(1.21, 464.7)`. There is no cheap way to buy that back on the
/// collapsed route: a compensated two-product on the multiply alone
/// recovers ~13% of it, because the multiply is not where the error is.
///
/// `ax`'s degenerate values ride the same formula rather than a separate
/// fallback: `log2` of `+0` is `-inf` and of `+inf`/NaN is itself, and
/// `exp2_checked_df`'s own clamp then saturates each to the right
/// `0`/`+inf`/NaN. Two selects on the *high word only* is all that costs
/// -- the low word may be arbitrary (or NaN) there, because every path
/// that reads it ends at `exp2_checked_df`'s own non-finite guard, which
/// already exists for the overflow case. Deciding those four magnitudes
/// separately instead needs the sign of `y` as well and lands at roughly
/// three times the ops.
macro_rules! powf_df_mag {
    ($ax:expr, $y:expr) => {{
        let ax = $ax;
        let l = log2_df(ax);
        let lh = if ax == 0.0 { f32::NEG_INFINITY } else { l.0 };
        let lh = if !(ax < f32::INFINITY) { ax * ax } else { lh };
        exp2_checked_df(Df32(lh, l.1) * $y)
    }};
}

// Shared by powf/rootn: the negative-base/y-parity/y==0/x==+-1
// special-case combine, given each caller's own already-computed `mag`.
// Macro, not a fn -- see exp_r_poly!.
//
// For negative x, `exp2(log2(|x|)*y)` alone can't ever be negative (exp2
// of any real argument is positive), so routing straight through `mag`
// always gave NaN for x < 0.0 -- even for a well-defined case like
// `(-2.0)^3.0 = -8.0`. A real result only exists there when y is an
// integer: even y -> +mag, odd y -> -mag (reusing `parity`, the same
// integer-parity helper sin_checked/cos_checked already use),
// non-integer y -> NaN (correctly matches std, e.g. `(-8.0)^(1/3)` is
// NaN in f32 too -- real cube roots of negative numbers aren't picked by
// this branch).
//
// `x == -0.0` and `x == -inf` are C99-exempt from the "non-integer y ->
// NaN" rule above: unlike a
// genuinely negative *finite* real number (where a non-integer power
// really is undefined), `-0` and `-inf` are signed *boundary* values
// whose magnitude-only result (`mag`) is always well-defined -- only the
// *sign convention* depends on `y` being an odd integer specifically,
// for *any* `y`, integer or not (e.g. `(-0.0).powf(0.5) == 0.0`, not
// `NaN`, since only an odd-integer exponent would have kept `-0`'s
// sign). `x == 0.0` catches `x == -0.0` here since this whole branch
// only runs when `x.is_sign_negative()` is already true.
//
// `y` infinite is a third, independent special case: C99 defines
// `pow(x, +-inf)` purely by `|x|` relative to `1` (`mag` already is
// exactly that), never sign-flipped by `x`'s own sign regardless of
// integer-ness. Overrides the selection above (not folded into its own
// condition) since it must win even when the `y_int`/`x==0`/
// `x.is_infinite()` check above would have produced a sign-flipped
// answer.
//
// `x.is_sign_negative()` (bit-based), not `x < 0.0` (value-based): the
// latter disagrees with the former exactly at x = -0.0 (same class of
// bug as acos's own `-0.0` fix), which would silently
// route `(-0.0)^3.0` through the wrong (positive) branch instead of the
// correctly-signed `-0.0`.
//
// pow(1, y) = 1 for *any* y -- even inf, -inf, or NaN -- another
// dedicated IEEE754/C99 special case the log/exp2 formula can't derive
// on its own (log2_df(1) is exactly Df32(0,0), so mag is exp2 of 0*y; for
// y=inf/-inf/NaN that's a 0*inf or 0*NaN indeterminate form, degrading to
// NaN instead of the correct 1). pow(-1, +-inf) = 1 is a second, narrower
// C99 special case (unlike pow(1,y), it does *not* extend to pow(-1,NaN),
// which stays NaN).
//
// pow(x, 0) = 1 for *any* x -- even 0, negative, or NaN -- a dedicated
// IEEE754/C99 special case, not derivable from the log/exp2 formula
// (0*inf and NaN*0 both degrade to NaN above). Override last.
//
// Everything above is decided by a *sign multiplier* in `{+1, -1, NaN}`
// applied to `mag` with one multiply, rather than by a tree of selects
// over the magnitude itself. Two things fall out of that shape:
//
//   - `parity(y)` alone separates all three cases. It is `0` for an even
//     integer and `1` for an odd one, and for a non-integer `y` it is
//     neither -- which is exactly the domain-error condition -- so the
//     separate `y == y.trunc()` integer test the select-tree needed is
//     redundant here. Gating `par` on `x`'s sign bit (rather than gating
//     the final result) is what lets the same two compares serve the
//     positive-base case, where the answer is always `mag`.
//   - the `x == +-1` cases stop needing their own overrides once `mag`
//     itself is pinned: `|+-1|^y` is exactly `1` for *every* `y`, so
//     pinning it there fixes both `pow(1, +-inf)`, `pow(1, NaN)` and
//     `pow(-1, +-inf)` at once, and the sign multiplier handles
//     `(-1)^y`'s parity from there.
//
// The select-tree form this replaces cost ~28 vector ops, most of them
// from `y_int || x == 0.0 || x.is_infinite()`: a compound `||` compiles
// to a hand-assembled AVX-512 mask (`korb`/`kmovd`/`cmovnel`) instead of
// a compare-and-blend, the same codegen trap `powf`'s own `is_safe` test
// documents. Every compare here feeds exactly one blend.
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
        // what a negative base with a *non*-integer y gives: a domain
        // error, except where |x| alone decides the answer -- `ax + ax ==
        // ax` picks out exactly `+0` and `+inf` (the two C99-exempt
        // boundary magnitudes) in one compare, and infinite y is the
        // third, independent exemption.
        // Infinite y specifically, not `y - y != 0.0`: that spelling is two
        // constants cheaper and looked free, because a NaN y makes `mag`
        // NaN anyway -- except at `ax == 1`, where `mag` is *pinned* to 1
        // just above. `powf(-1, NaN)` came back `1.0` instead of NaN,
        // C99's one deliberate asymmetry with `powf(1, NaN)`.
        let spec = if ax + ax == ax { 1.0 } else { f32::NAN };
        let spec = if y.abs() == f32::INFINITY { 1.0 } else { spec };
        let sm = if par == 0.0 { 1.0 } else { spec };
        let sm = if par == 1.0 { -1.0 } else { sm };
        let r = mag * sm;
        if y == 0.0 { 1.0 } else { r }
    }};
}

/// `x^y`, C99 `pow` semantics (see `powf_sign_combine!` above for the
/// full special-case table and `powf_df_mag!` for the magnitude).
#[doc(alias = "pow")]
#[inline(always)]
#[allow(clippy::neg_cmp_op_on_partial_ord)] // `!(ax < inf)` catches NaN too
pub fn powf(x: f32, y: f32) -> f32 {
    let ax = x.abs();
    let mag = powf_df_mag!(ax, y);
    powf_sign_combine!(x, ax, y, mag)
}

/// `powf` restricted to `x > 0.0`, or `x` exactly `+0.0` (contract, not
/// asserted -- see below for the one value this excludes and why): drops
/// `powf_sign_combine!`'s entire negative-base tree (`par`/`spec`/the
/// sign multiplier and its multiply) -- none of it is reachable once `x`
/// can't be negative, so this keeps only the two overrides that still
/// apply for any `x >= 0.0`: `y == 0.0 -> 1.0` (the `log2_df(x)*0`
/// route degrades to a `0*inf`/`0*NaN` indeterminate form for
/// `y = +-inf`/`NaN`, so this can't be derived from the formula) and
/// `x == 1.0 -> 1.0` (same reasoning, `log2_df(1)*y` is `0*y`, degenerate
/// for `y = +-inf`).
///
/// `x == -0.0` is the one `x >= 0.0`-*valued* input this doesn't handle
/// (despite `-0.0 >= 0.0` being true): `powf(-0.0, y)` preserves `-0.0`'s
/// sign for odd-integer `y` (`(-0.0).powf(3.0) == -0.0`,
/// `(-0.0).powf(-1.0) == -inf`), and that sign restoration is exactly
/// the parity machinery this function exists to skip -- adding
/// it back just for `-0.0` would cost the same parity work on *every*
/// call this function is meant to avoid, defeating the point. Documented
/// out rather than silently wrong: `powf_pos(-0.0, y)` gives the
/// same *magnitude* as `powf` but always with `+0.0`'s sign, not `-0.0`'s
/// (verified: this is the only real behavioral gap versus `powf` over
/// `x >= 0.0`, confirmed by a 300M-sample fuzz with `-0.0` itself
/// excluded and two edgecheck pins documenting the gap directly).
///
/// This is deliberately *not* C23 `powr`/IEEE754-2008 `powr`, despite
/// computing the same thing for finite, non-edge-case inputs: `powr`'s
/// own special-case table is stricter than `pow`'s (e.g. IEEE754 defines
/// `powr(1, NaN) = NaN` and `powr(0, 0) = NaN`, both matching the
/// "totally connected exp(y*log(x))" composition literally, where this
/// function -- matching `powf`'s own C99 `pow` conventions instead --
/// gives `1.0` for both, same as `powf` does for `x >= 0`). Named for
/// what it verifiably does, not for a standard it doesn't fully
/// implement. Behavior for `x < 0.0` (true negatives, not `-0.0`) is
/// unspecified (not `NaN`-guaranteed like `powf`'s own domain error --
/// whatever `log2_df`'s own decomposition and the two overrides above
/// happen to produce).
#[inline(always)]
#[allow(clippy::neg_cmp_op_on_partial_ord)] // `!(x < inf)` catches NaN too
pub fn powf_pos(x: f32, y: f32) -> f32 {
    let mag = powf_df_mag!(x, y);
    let r = if x == 1.0 { 1.0 } else { mag };
    if y == 0.0 { 1.0 } else { r }
}

/// "Signed power", the graphics/shading convention for raising a
/// possibly-negative value to a power without `powf`'s NaN-for-
/// non-integer-exponent domain error: computes on `|x|` then reapplies
/// `x`'s own sign unconditionally (`mulsign`, not folded into the
/// integer-exponent parity `powf` uses), regardless of whether `y` is an
/// integer. Total for every finite `x`/`y` -- no domain error, unlike
/// `powf(-2.0, 0.5)` which is correctly `NaN` (real exponentiation of a
/// negative base to a non-integer power has no real result) but is
/// exactly the semantics some callers explicitly don't want (signed
/// gamma curves, symmetric shaping functions).
///
/// Built on [`powf_pos`] rather than duplicating its formula: `x.abs()`
/// is always `+0.0` for either zero input (never `-0.0`, unlike a raw
/// sign-bit-preserving negation would give for `x == -0.0`), so this
/// never hits `powf_pos`'s own documented `-0.0` gap -- verified bit-
/// identical to `mulsign(powf_pos(x.abs(), y), x)` by construction, not
/// separately fitted.
#[inline(always)]
pub fn signed_pow(x: f32, y: f32) -> f32 {
    mulsign(powf_pos(x.abs(), y), x)
}

/// [`powf`] without domain/sign checks: valid for `x` positive, normal,
/// and finite (the same domain [`log_2_unchecked`] requires) and
/// `y != 0.0`. No zero/inf/nan handling on `x` (`powf_df_mag!`'s two
/// high-word selects), no `y == 0.0` special case, no negative-base
/// parity handling -- those are exactly the branches [`powf`] pays on
/// every call regardless of whether they're ever hit. Mirrors `log_2`/
/// `log_2_unchecked` and `atan2`/`atan2_unchecked`'s own fast/full-safety
/// split. Bit-identical to [`powf`] on this domain, and the same
/// double-float accuracy: `exp2_checked_df` stays even though every
/// domain check is gone, because the precision it preserves is what makes
/// the answer right, not a safety net -- see `powf_df_mag!`.
#[inline(always)]
pub fn powf_unchecked(x: f32, y: f32) -> f32 {
    exp2_checked_df(log2_df(x) * y)
}

/// x^(1/n) for integer `n` (backlog idea #75, C23 `rootn`), the `cbrt`
/// generalization -- unlike `powf(x, 1.0/n as f32)`, which would reuse
/// `powf`'s own negative-base parity check on `1/n` itself (almost never
/// an integer, so that check would wrongly call `rootn(-8, 3)` domain
/// error instead of the real `-2`), this checks `n`'s own parity
/// directly, matching the actual C23 semantics: real for negative `x`
/// only when `n` is odd, domain error (`NaN`) for negative `x` with
/// even `n` or for `n == 0` regardless of `x` (verified against the
/// real glibc rootn implementation notes, not guessed).
///
/// The magnitude never forms `log2(|x|)` as a single `f32`. That value
/// carries `|x|`'s whole binary exponent (up to 149), so its own last
/// place is worth `ulp(149) ~ 1.5e-5` -- an absolute error `exp2` turns
/// straight back into a relative one, and one that dividing by `n`
/// only shrinks proportionally, so it is worst exactly at the small
/// `|n|` a caller is most likely to write. Instead the exponent leaves
/// the float path entirely: `|x| = m * 2^e` split on `log_2`'s own
/// `[sqrt(2)/2, sqrt(2))` window (so `log_2_unchecked` never pays its
/// `lm + k` rounding either), and the *exact* Euclidean split
/// `e = q*n + rr` with `0 <= rr < |n|` gives
///
/// ```text
/// log2(|x|)/n = q + (rr + log2(m))/n
/// ```
///
/// `q` is an integer handed to [`exp2_kf`] as its exponent field, and
/// the fractional argument `(rr + log2(m))/n` is self-normalising:
/// `rr < |n|` bounds the sum's own rounding by `|n|*2^-25`, which the
/// division by `n` scales right back down to `2^-25` no matter how
/// large `e` or `n` are. Accuracy is then flat in `n` rather than
/// degrading as `|n|` falls.
///
/// `|n| == 1` is excluded from that path rather than accommodated by
/// it, which is what keeps the `exp2_kf` contract (`q` in
/// `[-126, 128)`, normal result) satisfiable without a saturating
/// scale: at `|n| >= 2` the result exponent is at most `~149/2` in
/// magnitude, so it is always normal, while `|n| == 1` is the one case
/// that can overflow or go denormal. Both are exact closed forms
/// anyway -- `|x|` for `n == 1` and `1/|x|` for `n == -1`, correctly
/// rounded by hardware over the whole domain including the overflowing
/// and denormal results -- and that same pair is *also* the right
/// magnitude for zero, infinite and NaN `x` at any `n`, so one select
/// covers both and the shared sign/parity logic finishes each.
///
/// Not wired into `examples/mca.rs`/`mca_target.rs`: this function's own
/// multi-exit-path branching (`n==0`/`n==1`/negative-even-domain-error)
/// corrupts llvm-mca's inline-asm region markers for the whole assembly
/// file, the same documented harness limitation `ldexp`/`frexp`
/// already have. Use quickbench for this one too.
#[inline(always)]
pub fn rootn(x: f32, n: i32) -> f32 {
    let ax = x.abs();
    // |x| = m * 2^e with m in [sqrt(2)/2, sqrt(2)) -- log_2_normal's own
    // window and its own bit trick, so log_2_unchecked(m) sees k == 0 and
    // returns log2(m) in [-0.5, 0.5) with no exponent to add back and no
    // `lm + k` rounding. Reaching for frexp instead costs its [0.5, 1) ->
    // window shift and its zero/infinite selects, which the `degenerate`
    // arm below re-does anyway.
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
    // `|x|` for `n > 0` and `1/|x|` for `n < 0` is the whole answer for
    // three separate reasons at once: it is the exact `|n| == 1`
    // identity, and it is also the right magnitude for zero, infinite
    // and NaN `x` (where the mantissa split reads nothing meaningful and
    // only `n`'s sign matters). `n == 0` lands here too and is discarded
    // by the domain-error select at the end.
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
    let r = if neg_even_domain_error { f32::NAN } else { signed };
    if n == 0 { f32::NAN } else { r }
}

/// sRGB -> linear (backlog idea #146), IEC 61966-2-1's piecewise
/// transfer function: a linear "toe" near black (avoiding the power
/// curve's infinite slope at 0) below `0.04045`, `((c+0.055)/1.055)^2.4`
/// above it. `powf_pos` composition (not a dedicated poly) since this
/// crate's own `powf` family already reuses the same log2/exp2_checked
/// machinery every other power-law composite here does; a dedicated fit
/// is only worth it if this measures as a real hot path. Domain `c>=0`
/// (matching `powf_pos`'s own contract) -- the standard sRGB channel
/// range.
///
/// Spelled out rather than calling `powf_pos` so the toe's own guard can
/// pay for itself: the power arm's value is discarded for every `c` below
/// `0.04045`, which is exactly the range where its base could be
/// zero/negative/denormal, so `log_2`'s wrapper collapses to the inf/NaN
/// arm alone (see `log_family_wrapper_discarded_unless_normal!`).
/// `powf_pos`'s two overrides go too: `y == 0.0` is a constant here, and
/// `x == 1.0` is redundant because `log_2(1) * y` is exactly `0`.
#[inline(always)]
#[allow(clippy::neg_cmp_op_on_partial_ord)]
pub fn srgb_to_linear(c: f32) -> f32 {
    let low = c * (1.0 / 12.92);
    let b = (c + 0.055) * (1.0 / 1.055);
    let high = exp2_checked(log_family_wrapper_discarded_unless_normal!(b, log_2_normal) * 2.4);
    if c <= 0.04045 { low } else { high }
}

/// linear -> sRGB (backlog idea #146), the inverse transfer function:
/// linear below `0.0031308`, `1.055*l^(1/2.4) - 0.055` above it. See
/// `srgb_to_linear`'s own doc comment for the composition rationale, the
/// domain contract (`l>=0`), and why the power arm is spelled out instead
/// of calling `powf_pos`.
#[inline(always)]
#[allow(clippy::neg_cmp_op_on_partial_ord)]
pub fn linear_to_srgb(l: f32) -> f32 {
    let low = l * 12.92;
    let p = exp2_checked(log_family_wrapper_discarded_unless_normal!(l, log_2_normal) * (1.0 / 2.4));
    let high = fma(1.055, p, -0.055);
    if l <= 0.0031308 { low } else { high }
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
///
/// The ties-away-from-zero rounding here is an inherited port convention
/// (see this doc comment's first line), not a deliberate design goal --
/// prefer [`remainder_ieee`] for new code unless matching that specific
/// convention is what's actually wanted. It's the true IEEE754/C99
/// standard besides, and genuinely cheaper (see its own doc comment).
// Shared by remainder/remainder_ieee/fmod: given each caller's own `q`
// (`.round()`, `.round_ties_even()`, or `.trunc()` -- the one place they
// differ), the `x - q*y` combine plus the two special-case selects (the
// `x==0.0` sign-preservation guard and the `y` infinite/`x` finite
// no-reduction case, see `remainder`'s doc comment). Macro, not a fn --
// see exp_r_poly!.
macro_rules! remainder_style_combine {
    ($x:expr, $y:expr, $q:expr) => {{
        let normal = fma(-$q, $y, $x);
        // IEEE754 defines the sign of an exact-zero remainder/fmod
        // result to match x's sign (also fmod's own C99 "same sign as
        // x" contract) -- but when x is a nonzero exact multiple of y,
        // `-q*y` exactly cancels x, and IEEE754 subtraction of two
        // equal-magnitude values always gives *positive* zero regardless
        // of the "intended" sign, silently dropping it (same class of
        // bug as sinf_poly's own -0.0 note, one level further out: here
        // the cancellation is real, not a q that happens to be zero).
        // Confirmed against libm (Python's math.remainder/math.fmod, a
        // real C library, not a hand-derived assumption) for every
        // exact-multiple case checked. Guarded here instead of only at
        // the `x == 0.0` case below, which only protects x itself being
        // zero, not the computed result rounding to zero.
        let normal = if normal == 0.0 { normal.copysign($x) } else { normal };
        let r = if $x == 0.0 && !normal.is_nan() { $x } else { normal };
        if $y.is_infinite() && $x.is_finite() { $x } else { r }
    }};
}

#[doc(alias = "remainderf")]
#[inline(always)]
pub fn remainder(x: f32, y: f32) -> f32 {
    let q = (x / y).round();
    remainder_style_combine!(x, y, q)
}

/// [`remainder`], but with the quotient rounded ties-to-even instead of
/// ties-away-from-zero, matching true IEEE754 remainder semantics (the
/// two conventions differ only when `x/y` lands on an exact half-integer
/// tie -- `remainder`'s own doc comment documents this crate's ties-away
/// choice as a deliberate divergence, not an oversight; this variant is
/// for callers who need the standard instead).
///
/// This is genuinely *cheaper* than `remainder`, not the same cost:
/// x86's `vroundss` only has hardware support for round-to-nearest-even,
/// round-down, round-up, and truncate -- there is no native "round half
/// away from zero" mode. `.round_ties_even()` (this function) lowers to
/// one `vroundss`; `.round()` (`remainder`'s choice) needs LLVM to
/// emulate the away-from-zero tie-break in software first, 4 extra
/// serial instructions -- a real ~15% latency gap (see readme.md).
#[inline(always)]
pub fn remainder_ieee(x: f32, y: f32) -> f32 {
    let q = (x / y).round_ties_even();
    remainder_style_combine!(x, y, q)
}

/// remainder without domain checks: valid for `x != 0.0` and `y` finite
/// (not `+-inf`) -- skips the two special-case selects [`remainder`]'s own
/// doc comment describes (the `x == 0.0` sign-preservation guard and the
/// `y` infinite/`x` finite no-reduction case). Mirrors this crate's other
/// `_unchecked` cores; see [`remainder`] for the full-domain-safe version.
///
/// One more gap than that list implies: also skips [`remainder`]'s
/// exact-cancellation sign fix (see `remainder_style_combine!`'s own
/// comment) -- for nonzero `x` an exact multiple of `y`, this returns
/// `+0.0` unconditionally rather than `x`-signed zero, since the
/// `-q*y` cancellation loses the sign IEEE754 always drops on an exact
/// zero result unless something restores it. A nonzero-`x` counterpart
/// to the `x == 0.0` gap already documented above, not a new mechanism.
#[inline(always)]
pub fn remainder_unchecked(x: f32, y: f32) -> f32 {
    let q = (x / y).round();
    fma(-q, y, x)
}

/// Self-correcting variant of [`remainder`]: detects when `x/y`'s own
/// division rounding pushed `q` to the wrong side of a half-integer tie
/// (see [`remainder`]'s doc comment for the failure mode -- a rare but
/// real sign-flip bug, not just a large-ratio accuracy gap) and nudges `q`
/// by one integer to correct it. A correctly-rounded `q` always leaves
/// `|r0| <= |y|/2`, so that inequality failing is a direct signal to move
/// toward whichever side shrinks the residual and recompute -- both
/// branches are computed unconditionally and selected, matching this
/// crate's branchless style. 0 max ulp against an f64 reference for
/// `|x/y|` up to `2^24`, where `q` itself stops being an
/// exactly-representable f32 integer -- a separate, harder limit this
/// correction can't reach past (see [`remainder_wide`]). Costs a second
/// `fma` plus the correction's compare/select on every call, so kept as
/// an opt-in tier, matching sin/sin_checked and exp2/exp2_checked.
///
/// `r1 = fma(-adj, y, r0)` instead of `fma(-(q0+adj), y, x)`:
/// algebraically `x - (q0+adj)*y == r0 - adj*y`, so `r1` reuses the
/// already-computed `r0` -- one fewer add, verified bit-identical to the
/// un-simplified form on dense near-tie samples.
#[inline(always)]
pub fn remainder_checked(x: f32, y: f32) -> f32 {
    let q0 = (x / y).round();
    let r0 = fma(-q0, y, x);
    // The correction always moves `r0` toward zero, so it is
    // `r0 - copysign(|y|, r0)` -- no `+-1` multiplier to select and no
    // fma. (The two forms differ only at `r0 == +-0.0`, where the guard
    // below keeps `r0` anyway.) `ay` is shared with that guard.
    let ay = y.abs();
    let r1 = r0 - ay.copysign(r0);
    let normal = if r0.abs() > ay * 0.5 { r1 } else { r0 };
    // Same exact-cancellation sign bug as remainder_style_combine! (see
    // its own comment): a nonzero x that's an exact multiple of y
    // exactly cancels to +0.0 regardless of x's sign, silently dropping
    // it. IEEE754 defines an exact-zero remainder's sign to match x's.
    let normal = if normal == 0.0 { normal.copysign(x) } else { normal };
    let r = if x == 0.0 && !normal.is_nan() { x } else { normal };
    if y.is_infinite() && x.is_finite() { x } else { r }
}

/// [`remainder_checked`], but correct for `|x/y|` past `2^24` too (up to
/// roughly `2^48`). `remainder_checked`'s self-correction assumes
/// `q0 = (x/y).round()` differs from the true integer quotient by at
/// most one; past `2^24`, `q0` can only land on a coarse grid (gaps of
/// `2^(e-23)` at exponent `e`), so the true quotient can be tens of
/// integers away -- a severe gap, not a tail case (~85% of samples in
/// `[2^24, 1e9]` landed on a completely different multiple of `y`).
///
/// Fixed by computing the residual with `Df32` instead of a single
/// `fma`: `q0*y` via `Df32::from_mul` (an exact two-product) subtracted
/// from `x` gives the *exact* real-valued residual, unlimited by `q0`'s
/// coarse quantization -- collapsing that to f32 and dividing by `y`
/// recovers the correction `adj` exactly. A second exact Df32
/// subtraction applies `adj`, then `remainder_checked`'s near-tie logic
/// (`adj2`/`r2`) runs unchanged on the now-correct residual. This is
/// `remainder_checked`'s correction loop run twice -- once for the
/// coarse `q0`, once for a near-tie -- not a new algorithm. 0 avg/max
/// ulp and 0 gross errors up to `|x/y| ~ 2^48` against a double-f64
/// (106-bit) reference (plain f64 is not trustworthy that far out),
/// degrading past that where a single correction pass is no longer
/// enough. Substantial mca cost over `remainder_checked` -- a separate
/// opt-in tier so callers who don't need the range don't pay.
///
/// Two exceptions vs `remainder_checked` in its own `|x/y| < 2^24`
/// domain used to exist here (bit-identical everywhere else, and now
/// bit-identical in both these corners too, verified below):
///
/// 1. When `x/y` landed on an *exact* half-integer, `adj`'s blind
///    `.round()` (ties away from zero) treated the already-correctly-
///    resolved tie (`q0` itself resolves ties away from zero, landing
///    `r0` on exactly `+-ys/2`) as "one more whole `y` to remove",
///    flipping the sign -- fixed by switching `adj` to
///    `.round_ties_even()`: an exact `+-0.5` now rounds to `0` (no
///    spurious correction) while every non-tie multi-integer-
///    quantization gap `adj` exists to recover (`|adj| >= 1`, nowhere
///    near a `.5` boundary) rounds identically either way. Also faster,
///    same finding `remainder_ieee`'s own `q` already made:
///    `round_ties_even` lowers to a single native `vroundps` where
///    `.round()`'s ties-away needs extra emulation.
/// 2. `Df32::from_mul(q0, y)` rounds the intermediate product to a
///    single f32 before pairing it with its error term -- unlike a
///    hardware `fma`, it can overflow on an intermediate value: `q0*y`
///    can exceed `x` by up to `|y|/2`, so when `x` sits within a small
///    factor of `f32::MAX` the exact product can exceed `f32::MAX` even
///    though the true remainder is finite (gave NaN). Rescaling `x`
///    (and `y`, to keep the ratio `x/y` unchanged) by an exact power of
///    two (`0.125`) whenever `|x|` is within a 4x margin of `f32::MAX`
///    is exact, not approximate: remainder is homogeneous of degree 1
///    (`remainder(k*x,k*y) == k*remainder(x,y)` for `k>0`). The
///    original gate was `max(|x|,|y|) > f32::MAX/4`, rescaling even
///    when only `y` was huge and `x` was small/denormal-ish -- pushing
///    that already-tiny `x` into the denormal range where low mantissa
///    bits are unrepresentable for no reason: `q0*y` always tracks `x`
///    itself (within `|y|/2`), so the overflow risk is governed
///    entirely by `|x|`, regardless of `|y|`'s own magnitude (`q0`
///    rounds to `~0` whenever `|x|` is tiny relative to `|y|`, so
///    `q0*y` stays safely small too). Gating on `|x|` alone fixes this.
///
/// Verified zero mismatches against `remainder_checked` over 351M
/// in-domain samples for the tie fix (previously ~4 in that many, all
/// exact ties) and 179M targeted denormal-x/huge-y samples for the
/// rescale-gate fix (previously ~215K in that many, up to a few ulp) --
/// both exclusions removed from `examples/unchecked_parity.rs`'s
/// standing test now that they're fixed. The rescale-gate fix is also
/// faster on its own (one fewer operand in the gating comparison).
#[inline(always)]
pub fn remainder_wide(x: f32, y: f32) -> f32 {
    let big = x.abs() > f32::MAX * 0.25;
    let scale = if big { 0.125 } else { 1.0 };
    let unscale = if big { 8.0 } else { 1.0 };
    let xs = x * scale;
    let ys = y * scale;
    let q0 = (xs / ys).round();
    let r0_df = Df32::from_f32(xs) - Df32::from_mul(q0, ys);
    let r0 = r0_df.to_f32();
    let adj = (r0 / ys).round_ties_even();
    let r1_df = r0_df - Df32::from_mul(adj, ys);
    let r1 = r1_df.to_f32();
    // The last whole-`y` correction always moves `r1` toward zero, so it
    // is `r1 - copysign(|ys|, r1)` -- no `+-1` multiplier to select and
    // no fma. (The two forms differ only at `r1 == +-0.0`, where the
    // guard below keeps `r1` anyway.) `ays` is shared with that guard.
    let ays = ys.abs();
    let r2 = r1 - ays.copysign(r1);
    let normal = if r1.abs() > ays * 0.5 { r2 } else { r1 };
    // `unscale` is selected, not computed as `1.0 / scale`: both scales
    // are exact powers of two, so their reciprocals are exact literals
    // and a second blend replaces a full-width divide.
    let normal = normal * unscale;
    // Same exact-cancellation sign bug as remainder_style_combine! (see
    // its own comment): a nonzero x that's an exact multiple of y
    // exactly cancels to +0.0 regardless of x's sign, silently dropping
    // it. IEEE754 defines an exact-zero remainder's sign to match x's.
    let normal = if normal == 0.0 { normal.copysign(x) } else { normal };
    let r = if x == 0.0 && !normal.is_nan() { x } else { normal };
    if y.is_infinite() && x.is_finite() { x } else { r }
}

/// C's `fmod(x,y)`: truncated (round-toward-zero) division instead of
/// [`remainder`]'s round-to-nearest, so the result always has the same
/// sign as `x` (or is a correctly-signed zero) -- a real, defining
/// difference from `remainder`, not just a different tie-break (C99
/// `fmod` and IEEE754 `remainder` are two genuinely different
/// operations, not two conventions for the same one, unlike
/// `remainder`/`remainder_ieee`'s own ties-away-vs-ties-even split).
/// Same structure as `remainder` otherwise (same `x==0.0` sign guard,
/// same finite-x/infinite-y no-reduction special case), just `.trunc()`
/// instead of `.round()`. Matches Rust's own `%` operator on `f32` (`%`
/// implements C `fmod` semantics), aside from a low-probability (~3e-7)
/// failure mode -- the truncation analog of `remainder`'s documented
/// one: `x/y`'s division rounding can land the computed quotient on the
/// wrong side of an exact *integer* (`.trunc()`'s decision boundary),
/// putting `q` off by a whole 1 and the result off by exactly `y`. Not
/// corrected here for the same reason as `remainder`.
#[doc(alias = "fmodf")]
#[inline(always)]
pub fn fmod(x: f32, y: f32) -> f32 {
    let q = (x / y).trunc();
    remainder_style_combine!(x, y, q)
}

/// Self-correcting variant of [`fmod`]: detects when `x/y`'s own division
/// rounding pushed `q` to the wrong side of an *integer* boundary (see
/// [`fmod`]'s doc comment -- the truncation analog of
/// [`remainder_checked`]'s half-integer-tie fix) and nudges `q` by one to
/// correct it. Unlike `remainder_checked`'s single symmetric `|r0| >
/// |y|/2` test, a correctly-truncated `r0` must satisfy an *asymmetric*
/// condition (same sign as `x`, magnitude strictly less than `|y|`), so
/// an off-by-one `q` shows up as one of two different symptoms needing
/// opposite corrections: `r0`'s sign disagreeing with `x`'s (q0
/// over-subtracted, needs `+y` back) or `|r0| >= |y|` (q0
/// under-subtracted, needs one more `-y`) -- never both at once, since
/// the error is bounded to exactly one integer step. `adj`'s sign
/// picks the right one directly rather than trying both candidates.
/// Verified zero mismatches against `x as f64 % y as f64` over 160M
/// in-domain (`|x/y|<1000`) samples (previously ~4.9e-7 rate, matching
/// `fmod`'s own documented ~3e-7 characterization).
#[inline(always)]
pub fn fmod_checked(x: f32, y: f32) -> f32 {
    let q0 = (x / y).trunc();
    let r0 = fma(-q0, y, x);
    let wrong_sign = r0 != 0.0 && (r0 > 0.0) != (x > 0.0);
    // The correction's sign depends on *both* which failure mode fired
    // (overshoot/`wrong_sign` needs q0 nudged toward zero, undershoot
    // needs it nudged away) *and* whether x/y have the same sign --
    // verified against all 4 sign combinations crossed with both
    // failure modes by hand (8 cases) before trusting this, not derived
    // by inspection alone (an earlier y-sign-only version passed 4 of
    // those 8 and failed the rest).
    let same_sign = (x > 0.0) == (y > 0.0);
    let adj = if wrong_sign != same_sign { 1.0 } else { -1.0 };
    let r1 = fma(-adj, y, r0);
    let needs_fix = wrong_sign || r0.abs() >= y.abs();
    let normal = if needs_fix { r1 } else { r0 };
    // Same exact-cancellation sign bug as remainder_style_combine! (see
    // its own comment): a nonzero x that's an exact multiple of y
    // exactly cancels to +0.0 regardless of x's sign, silently dropping
    // it. IEEE754/C99 define an exact-zero fmod result's sign to match
    // x's. `wrong_sign` above deliberately excludes `r0 == 0.0` (it's
    // testing for a *nonzero* sign disagreement to detect the off-by-
    // one case), so this exact-zero case reaches here unflagged.
    let normal = if normal == 0.0 { normal.copysign(x) } else { normal };
    let r = if x == 0.0 && !normal.is_nan() { x } else { normal };
    if y.is_infinite() && x.is_finite() { x } else { r }
}

/// [`fmod`] without domain checks: valid for `x != 0.0` and `y` finite
/// (not `+-inf`) -- mirrors [`remainder_unchecked`]'s own contract and
/// reasoning exactly, just for `fmod`'s truncated convention.
#[inline(always)]
pub fn fmod_unchecked(x: f32, y: f32) -> f32 {
    let q = (x / y).trunc();
    fma(-q, y, x)
}

/// Rust-native `f32::rem_euclid` semantics (always-nonnegative
/// remainder, `0 <= result < |y|` for any finite nonzero `y`), built on
/// [`fmod`] instead of duplicating a general-purpose implementation:
/// `fmod`'s truncated remainder already has the right *magnitude*, so
/// this only needs to add `|y|` back whenever that remainder came out
/// negative (matching std's own formula exactly, just routed through
/// this crate's own already-optimized `fmod` instead of a libm call).
/// Real, substantial speed win over `f32::rem_euclid` in practice
/// (~2.3x in a direct wall-clock comparison) -- `fmod` is branchless and
/// vectorizes, where the general std path doesn't clearly do either.
///
/// Verified bit-identical to `f32::rem_euclid` over a 106M-sample fuzz
/// (`|x/y| < 1000`, matching `fmod`'s own established accuracy domain)
/// plus every special-value combination `f32::rem_euclid` itself
/// defines (zero, `+-0.0`, infinities, `y == 0.0`, `NaN`) -- correctly
/// benefits from the exact-cancellation sign fix `fmod` just received,
/// unlike a naive reimplementation would have without it. Non-obvious
/// bonus: this is *also* immune to `fmod`'s own documented rare (~3e-7)
/// off-by-a-whole-`y` quotient-boundary miss, confirmed over the same
/// 106M samples -- both that bug and this function's own `+|y|`
/// normalization step are "shift by exactly one `y`" adjustments on the
/// same quantity, so whenever `fmod`'s rare miss lands short by one `y`
/// in the negative direction, the normalization this function already
/// needs happens to correct it for free (see [`div_euclid`]'s own doc
/// comment for the case that *doesn't* get this same free correction).
#[inline(always)]
pub fn rem_euclid(x: f32, y: f32) -> f32 {
    let r = fmod(x, y);
    if r < 0.0 { r + y.abs() } else { r }
}

/// Rust-native `f32::div_euclid` semantics: the integer (well-defined
/// non-integer-typed here, since this is `f32`) quotient paired with
/// [`rem_euclid`] (`x == div_euclid(x,y) * y + rem_euclid(x,y)` for any
/// finite nonzero `y`, `rem_euclid` always in `[0, |y|)`). Same
/// `fmod`-truncation-plus-adjustment shape as `rem_euclid`, computing
/// its own `fmod` rather than sharing `rem_euclid`'s (no public two-
/// value return in this crate's style, and duplicating the one extra
/// `fma` this needs is cheaper than a tuple-returning detour would be
/// for callers who only want one half). Verified bit-identical to
/// `f32::div_euclid` over the same 106M-sample fuzz plus special-value
/// matrix as `rem_euclid`, *except* one gap `rem_euclid` doesn't share:
/// `q = (x/y).trunc()` here is computed independently of `fmod`'s own
/// (self-consistent) quotient, so unlike `rem_euclid` (which absorbs
/// `fmod`'s rare off-by-a-whole-`y` miss for free, see its own doc
/// comment), this one inherits it directly -- confirmed at essentially
/// the same rate `fmod` itself documents (47/106M samples, ~4.4e-7 vs
/// `fmod`'s own ~3e-7). Not chased further: matches an already-accepted,
/// already-documented tradeoff of the primitive this is built on, not a
/// new defect.
#[inline(always)]
pub fn div_euclid(x: f32, y: f32) -> f32 {
    let q = (x / y).trunc();
    let r = fmod(x, y);
    if r < 0.0 {
        if y > 0.0 { q - 1.0 } else { q + 1.0 }
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
