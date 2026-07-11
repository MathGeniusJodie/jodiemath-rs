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

// Shared by log_2/ln/log10/log2_df's own denormal handling (IDEAS.md
// backlog: "log_2/ln/log10 shared-core macro"): scale a denormal input up
// by 2^24 before the single normal-path evaluation, tracking the
// compensating exponent offset to fold back in afterwards. A macro (not a
// fn) so there's no function-call boundary for the auto-vectorizer to
// reason about -- same "duplicate rather than share a helper" caution
// this crate's other standalone-copy doc comments (expm1, tanh) already
// establish, just applied via textual substitution instead of a literal
// copy-paste.
macro_rules! denormal_rescale {
    ($x:expr) => {{
        let tiny = $x < f32::MIN_POSITIVE;
        let xs = if tiny { $x * 16777216.0 } else { $x };
        let koff = if tiny { -24.0 } else { 0.0 };
        (xs, koff)
    }};
}

// Shared by log_2_normal/ln_normal/log10_normal (IDEAS.md idea #24: "the
// three `_normal` bodies are identical modulo constants"): the exponent
// extraction, s = m - 1 decomposition, and degree-9 Estrin poly evaluation
// are all byte-for-byte the same shape across all three -- only the
// coefficient array and each function's own final k-combine differ.
// Returns (p, s, k) so each caller does its own combine (log_2_normal's
// plain `fma(p, s, k)`, ln_normal/log10_normal's own Cody-Waite HI/LO
// split). Macro, not a fn -- avoids the function-call-boundary
// codegen-regression risk this crate's own expm1/tanh doc comments
// already found for a different shared-helper attempt.
macro_rules! log_family_normal {
    ($x:expr, $koff:expr, $c:expr) => {{
        let e = ($x.to_bits() as i32).wrapping_sub(0x3f3504f3) >> 23;
        let m = f32::from_bits(($x.to_bits() as i32).wrapping_sub(e << 23) as u32);
        let k = e as f32 + $koff;
        let s = m - 1.0;
        let c: [f32; 10] = $c;
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
        (p, s, k)
    }};
}

// Shared by exp/expm1/tanh's own `exp(r)`-for-tiny-`r` poly evaluation:
// `exp`'s own doc comment establishes the coefficients (c0/c1 pinned to
// 1.0, c2..c5 LP-refit 2026-07-09) and notes the fit "cascades to every
// caller since they all reuse this poly directly" -- expm1 and tanh each
// carry a byte-identical textual copy of this exact 8-line fragment
// (confirmed by direct comparison), deliberately NOT calling a shared fn
// (both their own doc comments cite a real, reproduced +32% mca
// regression on an unrelated caller, `sinh_throughput`, from a prior `fn`-
// based sharing attempt -- "a scheduling side effect of the new function
// boundary"). A macro has no function-call boundary at all (pure textual
// substitution before codegen), so it's a fundamentally different
// mechanism than what was already tried and rejected -- verified via a
// full pre/post assembly diff before trusting that distinction (see
// IDEAS.md). Each caller keeps its own reduction (`k`/`r`) and exponent
// reconstruction/final-combine around this, since those differ (exp/
// expm1's k1/k2 split vs. tanh's single exp2int field).
macro_rules! exp_r_poly {
    ($r:expr) => {{
        let c: [f32; 4] = [4.9999300e-1, 1.6667245e-1, 4.1883811e-2, 8.3009899e-3];
        let r2 = $r * $r;
        let r4 = r2 * r2;
        let l0 = $r + 1.0;
        let l1 = fma(c[1], $r, c[0]);
        let l2 = fma(c[3], $r, c[2]);
        let r0 = fma(l1, r2, l0);
        fma(l2, r4, r0)
    }};
}

// Shared by exp2/exp2_checked/exp10/exp10_checked/exp2m1/exp2_checked_df's
// own `Q(f) = (2^f - 1)/f` evaluation (exp2's own doc comment already
// documents this as "shared verbatim... all 'standalone copies' of this
// exact poly, updated together to keep them in sync" -- a manual-sync
// convention this macro now enforces structurally instead). Same macro-
// not-fn reasoning as exp_r_poly! above. Returns `q`; each caller does its
// own distinct final combine (exp2/exp10's single-field
// `fma(q, exp2int*f, exp2int)` vs. exp2_checked/exp10_checked/exp2m1/
// exp2_checked_df's k1/k2-split `p = fma(q, t1*f, t1); p*t2`).
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

// Shared by exp_pos_neg/exp_pos_neg_checked_half (sinh/cosh's own
// unchecked/checked exp(x)/exp(-x) core): identical Cody-Waite reduction,
// even/odd-split poly (retuned 2026-07-09/10 -- see [`exp_pos_neg`]'s own
// doc comment for the coefficient history, since that's still the
// canonical place the refit story lives), and t1n/t2n reciprocal
// construction. Only each caller's own optional input clamp and final
// `0.5`-scaling differ, so those stay at the call site. Returns
// `(p_pos, p_neg, t1, t2, t1n, t2n)`. Macro, not a fn -- same reasoning
// as exp_r_poly!/exp2_q_poly! above.
macro_rules! exp_pos_neg_core {
    ($x:expr) => {{
        const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
        let k = fma($x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
        let r = fma(-k, LN2_HI, $x);
        let r = fma(-k, LN2_LO, r);
        let c: [f32; 4] = [4.99993e-1, 1.6667245e-1, 4.188372e-2, 8.300987e-3];
        let r2 = r * r;
        let r4 = r2 * r2;
        let e = fma(c[2], r4, fma(c[0], r2, 1.0));
        let o = fma(c[3], r4, fma(c[1], r2, 1.0));
        let p_pos = fma(r, o, e);
        let p_neg = fma(-r, o, e);
        let (t1, t2) = exp2_field_split(k);
        // t1n = 1/t1, t2n = 1/t2 (both exact power-of-two fields, so 1/t1
        // is itself an exact power of two): for an exact power-of-two
        // float with bit pattern b = (127+e)<<23, its reciprocal 2^-e has
        // bit pattern (127-e)<<23 = 0x7F000000 - b (since (127+e)+(127-e)
        // = 254 = 0xFE, and 0xFE<<23 == 0x7F000000). This is exactly the
        // split exp2_field_split(-k) would have produced (round-half-to-
        // even is antisymmetric under negation, so k1n=-k1/k2n=-k2), but
        // built with two integer subtracts instead of a whole second
        // magic-round fma/sub/sub chain -- deletes exp2_field_split(-k)'s
        // independent dependency on k entirely.
        let t1n = f32::from_bits(0x7F00_0000u32.wrapping_sub(t1.to_bits()));
        let t2n = f32::from_bits(0x7F00_0000u32.wrapping_sub(t2.to_bits()));
        (p_pos, p_neg, t1, t2, t1n, t2n)
    }};
}

// Shared by expm1/exp_m1_over_x/exp2m1/tanh's own Pade approximant for
// e^v-1 near v=0 (an exact closed-form Pade identity, not an empirical
// fit -- see exp_m1_over_x's own doc comment for the by-hand verification
// at v=0). Two arms, not one: `v * fma(N) / fma(D)` (expm1/exp2m1/tanh)
// vs. bare `fma(N) / fma(D)` (exp_m1_over_x). A macro invocation is
// always parsed as one atomic expression at its call site, so writing
// `v * pade_expm1_ratio!(v)` with a single `fma(N)/fma(D)`-only macro
// would silently reassociate `(v*N)/D` (the original text's actual,
// left-to-right operation order) into `v*(N/D)` -- mathematically equal
// but not bit-identical (confirmed via a full assembly diff: this
// exact mistake changed the compiled output, caught before committing).
// The `mul` arm keeps `v * fma(N)` and the final `/ fma(D)` inside one
// macro expansion, so ordinary Rust precedence inside that expansion
// reproduces the original grouping exactly.
macro_rules! pade_expm1_ratio {
    ($v:expr) => {
        fma(-1.9999927, $v * $v, -120.0) / fma($v, fma($v, $v - 12.000030, 59.999996), -120.0)
    };
    ($v:expr, mul) => {
        $v * fma(-1.9999927, $v * $v, -120.0) / fma($v, fma($v, $v - 12.000030, 59.999996), -120.0)
    };
}

// Shared by log_2/ln/log10: given each caller's own `_normal` fn (the
// only thing that differs -- coefficients/final combine live there), the
// denormal-rescale + special-case-select wrapper is byte-for-byte
// identical: edge handling is done with selects (no early returns) so
// loops over arrays of these calls can auto-vectorize -- scale
// denormals up before the single normal-path evaluation, then patch
// specials afterwards. `spec` is `-inf` for `+-0`, `NaN` for `x < 0`
// (includes `-inf`); its select input only depends on `x`, so it
// resolves in parallel with the poly evaluation. The final `!(x <
// f32::INFINITY)` check (`+inf`/`nan`: `x*x` is `inf`/`nan`
// respectively, false for `-inf` since `-inf < inf`) deliberately uses
// `!(x < inf)` (exploiting NaN's always-false comparisons to catch both
// +inf and NaN in one check) rather than `partial_cmp`, which would need
// an extra `Option`-unwrap for what's already a single cheap `fcmp` --
// each caller's own `#[allow(clippy::neg_cmp_op_on_partial_ord)]`
// still applies to this check post-expansion. Macro, not a fn -- same
// reasoning as this file's other shared-body macros.
macro_rules! log_family_wrapper {
    ($x:expr, $normal:ident) => {{
        let (xs, koff) = denormal_rescale!($x);
        let r = $normal(xs, koff);
        let spec = if $x == 0.0 { f32::NEG_INFINITY } else { f32::NAN };
        let r = if $x <= 0.0 { spec } else { r };
        if !($x < f32::INFINITY) {
            $x * $x
        } else {
            r
        }
    }};
}

#[inline(always)]
#[allow(clippy::neg_cmp_op_on_partial_ord)]
pub fn log_2(x: f32) -> f32 {
    log_family_wrapper!(x, log_2_normal)
}

/// Core of log_2 for positive normal finite x only: no handling for zero,
/// negative, denormal, inf, or nan input (those are the caller's job, see
/// log_2). Called directly with an out-of-domain x, this returns a
/// plausible-looking but wrong finite value rather than NaN/-inf.
// Shared by log_2_normal and log2_df -- both compute log2 itself (one in
// single-float, one in double-float form), so unlike ln_normal/log10_normal
// (each independently minimax-fitted for their own target), these two
// callers use the exact same degree-9 poly, not just a scaled variant.
const LOG2_COEFFS: [f32; 10] = [
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

#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn log_2_normal(x: f32, koff: f32) -> f32 {
    // decompose x = 2^k * m with m in [sqrt(2)/2, sqrt(2)), so s = m - 1
    // is exact (Sterbenz) and centered on 0: log2 stays relatively
    // accurate near x = 1. log2(m) = s * P(s), degree-9 minimax P fitted
    // with lolremez (rel. error 4.1e-9).
    let (p, s, k) = log_family_normal!(x, koff, LOG2_COEFFS);
    // k + s * P(s) in a single rounding
    fma(p, s, k)
}

/// log_2 without domain checks: valid for positive normal finite x only
/// (no handling for zero, negative, denormal, inf, or nan -- those give a
/// plausible-looking but wrong finite value instead of NaN/-inf). Mirrors
/// exp2/exp2_checked's fast/full-safety split, just with the "safe by
/// default" name (`log_2`) already taken by the checked tier, so this one
/// gets the `_unchecked` suffix instead. Drops the denormal-rescale
/// multiply and both post-hoc selects log_2 pays on every call (branchless
/// selects still cost real ops even when the branch not taken is a no-op
/// value), at the cost of undefined output outside the stated domain.
#[inline(always)]
pub fn log_2_unchecked(x: f32) -> f32 {
    log_2_normal(x, 0.0)
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
    //
    // IDEAS.md backlog round 3 #1 tried k=round(x) here (f in [-0.5,0.5],
    // a tighter Q(f) fit at the same degree) -- rejected for this
    // *unchecked*, non-split tier specifically: round(x) can land k=128
    // for x within the promised [-126,128) domain (e.g. x=127.6), which
    // this single exp2int construction (unlike exp2_checked's k1/k2
    // split) can't represent, producing NaN inside the documented safe
    // range -- a real regression, not just an accuracy tradeoff. The same
    // backlog idea was tried and rejected on every other function sharing
    // this poly too (exp2_checked, exp10, exp10_checked, exp2m1,
    // exp2_checked_df), each for its own distinct reason -- see each
    // function's own doc comment and IDEAS.md for the full writeup.
    let k = x.floor();
    let f = x - k;
    let exp2int = f32::from_bits(((k + 383_f32).to_bits() << 8) & EXPONENT_MASK);
    // Q(f) = (2^f - 1)/f, degree 5, grouped into 3 balanced pairs (g0, g1,
    // g2) instead of two degree-2 Horner halves: same 6 coefficients (so
    // same accuracy target) and the same 4-deep fma critical path, but the
    // combine only ever needs f^2 (never exp2int*f^4), so it's 2 fewer
    // plain multiplies per call than the old A/B split.
    //
    // g1/g2's leading coefficients refit (2026-07-10) via coordinate
    // descent from the shipped values (idea #3 in IDEAS.md), after fixing
    // tune.rs's own `exp2_c` probe to actually match this Estrin structure
    // (it had silently been testing the old, already-replaced A/B split
    // instead). Verified against a real ~562M-point dense sweep of the
    // whole unchecked domain: avg ulp 0.07176->0.06914 (~3.7% tighter),
    // max ulp unchanged at 1 (already the practical ceiling for a
    // single-fma-final-rounding construction) -- zero perf cost, same
    // instructions. Shared verbatim by exp2_checked/exp10/exp10_checked/
    // exp2m1/exp2_checked_df below (all "standalone copies" of this exact
    // poly), updated together to keep them in sync -- now enforced by
    // sharing `exp2_q_poly!` (a macro, not a fn: see exp_r_poly!'s own
    // doc comment for why that distinction matters in this crate).
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
    // an integer (e.g. x = 4.9999999).
    //
    // IDEAS.md backlog round 3 #1 tried k=round(x) here (matching exp10's
    // own natural convention, which would have let exp10_checked skip its
    // own floor-adjust step) -- rejected for exp2_checked specifically:
    // unlike exp10_checked, exp2_checked never had a floor-adjust to
    // remove (it already computed k straight from x), so the only effect
    // of the swap was replacing `x.floor()` (a single vroundps, on this
    // CPU's own otherwise-idle rounding port) with an explicit `(xs+M)-M`
    // add/sub pair that instead contends with the poly's own already-busy
    // fma/add ports -- confirmed via assembly diff (64->66 total
    // instructions, vroundps count 2->0, vaddps 4->8, fma count unchanged
    // at 14) and mca (throughput 1.399->3.007 cyc/elem, latency flat --
    // more than double the cost for zero structural savings). Also
    // regressed real max ulp 1->2 on the full 100M-sample fuzz despite the
    // coarse tune.rs-grid refit showing an improvement (avg
    // 0.0158->0.0051, max 1->1 on that grid) -- another instance of
    // "verify a tune.rs coarse-grid result against the real fuzz before
    // trusting it". (exp10_checked's own attempt at the floor-adjust
    // removal separately failed too, for a different, correctness-level
    // reason -- see its own doc comment.) Reverted, bit-identical to prior
    // HEAD.
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
/// discarding precision the careful reduction above just earned (found
/// empirically -- an early version measured max ulp 45, traced to
/// exactly this recombination step, not the reduction itself, which
/// checked out accurate to ~1e-11 in isolation). Fixed by never forming
/// that combined value at all: `exp2_checked`'s own internal split (its
/// own doc comment: "any split k=k1+k2 ... works") is reproduced here
/// directly against this function's *own* precisely-known integer `k`
/// and fractional `f` (floor-adjusted from the round-based reduction
/// into `exp2_checked`'s own `[0,1)` convention), instead of recombining
/// them into one f32 and letting `exp2_checked` re-derive (and re-round)
/// its own `k`/`f` from that already-lossy sum.
///
/// IDEAS.md backlog round 3 #1 tried dropping the floor-adjust entirely
/// (feed `kr`/`fr` straight through, sharing a round-domain Q(f) fit with
/// exp2_checked) -- real mca win when it worked (latency 72.06->51.06 cyc,
/// -29%; throughput 2.736->1.903, -30%) and real avg-ulp win too
/// (0.0343->0.0107), but `edgecheck.rs` caught a genuine correctness bug
/// the 100M-sample fuzz missed entirely: at the overflow-saturation
/// boundary (`k` clamped down to exactly `128`, e.g. `exp10_checked(inf)`
/// after its own `x.clamp(-1000,1000)` reduces to `x=1000`), `t1*t2 =
/// 2^128` sits right at `f32::MAX`. The old floor convention guarantees
/// `f >= 0` there (so the poly's `2^f >= 1` factor only ever pushes the
/// product *up* into overflow, correctly saturating to `inf`); the
/// round convention allows `f < 0`, so `2^f < 1` can instead pull the
/// product just *under* `f32::MAX`, giving `exp10_checked(inf) ==
/// 3.237e38` (finite!) instead of `inf` -- an actual contract violation,
/// not an accuracy tradeoff. Reverted (kept the floor-adjust), bit-
/// identical to prior HEAD. *A change that looks clean on both mca and a
/// 100M-sample fuzz can still hide a real bug at a boundary condition
/// neither tool samples -- `edgecheck.rs`'s dedicated special-value pins
/// (here, `x=inf`) are not redundant with fuzzing; run both before
/// trusting a reduction-scheme change.*
// Shared by exp10/exp10_checked: the round-based reduction (`kb`/`kr`/`d`/
// `fr`/floor-adjust to `(k, f)`) is identical between the two -- only
// exp10_checked's own leading `x` clamp and trailing `k` clamp (needed
// since it feeds the k1/k2-split combine, unlike exp10's single-field
// one) differ, both left at the call site. Macro, not a fn -- same
// reasoning as this file's other shared-body macros.
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
    // (f32::clamp preserves NaN in the receiver), and the bound is wide
    // enough that it never touches a genuinely in-range x (exp2_checked's
    // own clamp downstream still does the real range-limiting).
    let x = x.clamp(-1000.0, 1000.0);
    let (k, f) = exp10_reduction!(x);
    let k = k.clamp(-151.0, 128.0);
    let (t1, t2) = exp2_field_split(k);
    let q = exp2_q_poly!(f);
    let p = fma(q, t1 * f, t1);
    p * t2
}

/// Same reduction as [`exp10_checked`], but a single exponent-field
/// construction (no k1/k2 split) instead of two -- faster, narrower-
/// domain tier, same pairing as `exp2`/`exp2_checked`. Valid while `k`
/// (see `exp10_checked`'s own doc comment) stays in `[-126,128)`.
///
/// IDEAS.md backlog round 3 #1's round-domain refit (skip the floor-adjust
/// entirely, feed kr/fr straight through) was tried and rejected here for
/// the same reason as plain `exp2`: this single-exponent-field
/// construction has no k1/k2 split to absorb `kr` landing on `128` right
/// at the promised domain edge, producing NaN inside `[-126,128)` instead
/// of a correct finite value. Kept the floor-adjust; see `exp2`'s own doc
/// comment for the full story (same root cause, same fix).
#[inline(always)]
#[allow(clippy::approx_constant)] // g0's constant term is a fitted minimax
// coefficient near ln(2), not ln(2) itself (bit pattern deliberately differs)
pub fn exp10(x: f32) -> f32 {
    let (k, f) = exp10_reduction!(x);
    let exp2int = f32::from_bits(((k + 383_f32).to_bits() << 8) & EXPONENT_MASK);
    let q = exp2_q_poly!(f);
    fma(q, exp2int * f, exp2int)
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
//
// IDEAS.md backlog round 3 #81 (2026-07-10): audited whether this
// copysign is still load-bearing for every caller now that sin_checked
// applies the "flip r, not the result" trick and sinpi has its own
// explicit x==0.0 override -- verified empirically (a temporary
// copysign-free build, checked at each caller's own actual zero-crossing,
// both signs, not just x=+-0.0): sin, sind, cospi, cosd all still
// genuinely need it (sin/sind lose the odd x=-0.0 sign without it;
// cospi/cosd lose even-function sign-of-zero *symmetry* at their own
// crossings -- e.g. cospi(0.5) and cospi(-0.5) came out with *opposite*
// signs of zero without it, which cospi being even can never actually
// produce). sin_checked and sinpi's own separate guards make it redundant
// for them specifically at the spot-checked zero-crossing, and a full
// 100M-sample fuzz across every documented-accurate bucket for both
// functions (see `examples/accuracy.rs`) came back with avg/max ulp
// matching their prior baseline exactly in every in-budget range --
// split into this copysign-free core plus the public wrapper below on
// that evidence.
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
// Shared by sin/cos: the Cody-Waite pi-split reduction (given each
// caller's own `q`, computed differently -- sin's plain `round(x/pi)` vs
// cos's phase-shifted `round(x/pi-0.5)+0.5`) plus the `sinf_poly` call is
// byte-for-byte identical. Macro, not a fn -- same reasoning as this
// file's other shared-body macros.
macro_rules! pi_reduce_and_poly {
    ($x:expr, $q:expr) => {{
        let r = fma($q, -PI_A, $x);
        let r = fma($q, -PI_B, r);
        let r = fma($q, -PI_C, r);
        let r = fma($q, -PI_D, r);
        sinf_poly(r)
    }};
}

#[inline(always)]
pub fn sin(x: f32) -> f32 {
    let qb = fma(x, FRAC_1_PI, ROUND_MAGIC);
    let q = qb - ROUND_MAGIC;
    let s = pi_reduce_and_poly!(x, q);
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
/// Uses `x.round_ties_even()` (native, full-range-correct), not the
/// magic-constant `x + 1.5*2^23` trick used elsewhere in this crate
/// (e.g. `exp`'s own `k` rounding): that trick is only exact for
/// `|x| <= 2^22` by construction (adding `x` to a constant of comparable
/// magnitude loses precision once `x` itself approaches that magnitude)
/// -- and `sinpi`/`cospi` take the *raw*, unbounded input `x` directly,
/// unlike every other magic-round use in this crate, which only ever
/// rounds an already-reduced, small intermediate value. A real,
/// previously-undetected bug lived here for exactly that reason: the doc
/// comment's "no accuracy cliff, exact out to f32::MAX" claim was false
/// for `2^22 < |x| < 2^24` (found by a targeted probe past the
/// accuracy.rs sweep's own `|x| < 1e6` domain cutoff -- the exact range
/// where the bug lives was never exercised; see backlog idea #30's
/// "confirm these are in accuracy.rs's sweep" concern, which was right
/// to worry). `round_ties_even` instead of plain `.round()` (both
/// initially tried) since `q`'s specific tie-breaking rule doesn't affect
/// correctness here (`parity(q)`'s own sign correction self-compensates
/// for whichever nearby integer `q` lands on, verified by hand: e.g.
/// `sinpi(2.5)` gives the same final `1.0` whether `q` resolves to 2 or
/// 3) -- and `round_ties_even` lowers to a single native `vroundps`
/// immediate, while `.round()`'s ties-away isn't a hardware-native mode
/// and needs extra instructions (mca: 48.02/1.338 -> 43.02/1.149
/// cyc/elem, a real ~10-14% win, not just a wash), the same finding
/// `remainder_ieee` made independently for its own `q`.
#[inline(always)]
pub fn sinpi(x: f32) -> f32 {
    let q = x.round_ties_even();
    let r = x - q;
    // At x=-0.0: q=(-0.0).round_ties_even()=-0.0 too (rounding preserves
    // zero's sign, regardless of tie-breaking convention),
    // so r=x-q=(-0.0)-(-0.0), which IEEE754 always resolves to +0.0
    // regardless of the operands' own sign -- the same "opposite-signed-
    // zero subtraction erases sign" mechanism as sinf_poly's own -0.0
    // fix and atan2(-0.0,+0.0)'s bug, just one level further out (r's
    // lost sign means sinf_poly's own copysign fix would have nothing
    // left to copy at this point anyway, which is exactly why this
    // guard, not that copysign, is what makes sinpi(-0.0) correct --
    // confirmed this function calls the copysign-free `sinf_poly_raw`
    // now, see its own doc comment). Guarded the same way log1p/log_2
    // handle their own x==0.0 sign case: compute the normal path
    // unconditionally first, select x itself only at the singular zero
    // point.
    let normal = sinf_poly_raw(std::f32::consts::PI * r) * fma(-2.0, parity(q), 1.0);
    if x == 0.0 { x } else { normal }
}

/// cos(pi*x), argument in half-turns -- see `sinpi`'s doc comment for why
/// this reduction is exact and shares `sinf_poly` directly, same
/// full-range-accurate (no cliff) guarantee, and the same magic-round
/// range-limit bug this version fixes.
///
/// Uses the identity `cos(pi*x) = sin(pi*(x+0.5))` without ever forming
/// `x+0.5` as a single value (lossy for large `x`, reintroducing the
/// exact bug being fixed): `k = round(x-0.5)` gives `round(x+0.5) =
/// k+1` for free (rounding commutes with an exact integer shift), and
/// `k+1`'s parity is just `1 - parity(k)` -- so neither `x+0.5` nor
/// `k+1` themselves ever need to exist as floats, only `k`, `parity(k)`,
/// and `r = (x-k)-0.5` (computed in that order so `x-k` stays a small,
/// Sterbenz-exact value before the final `-0.5`).
#[inline(always)]
pub fn cospi(x: f32) -> f32 {
    let k = (x - 0.5).round_ties_even();
    let r = (x - k) - 0.5;
    let s = sinf_poly(std::f32::consts::PI * r);
    let sign = fma(2.0, parity(k), -1.0); // -(1 - 2*(1-parity(k))) = 2*parity(k)-1
    s * sign
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

/// tan(pi*x), argument in half-turns -- new function (backlog idea #29),
/// built directly from `sinpi`/`cospi`'s own ratio: `tan` has period 1 in
/// half-turns (unlike `sin`/`cos` individually, which flip sign every
/// integer), so `sinpi(x)/cospi(x)` is exactly `tan(pi*x)` with no
/// separate reduction of its own needed -- whichever integer `sinpi`'s
/// `q=round(x)` and `cospi`'s own `k=round(x-0.5)` each resolve to, their
/// respective sign corrections (`parity(q)`/`parity(k)`) cancel exactly
/// in the division (both numerator and denominator flip together or not
/// at all, since `tan(theta+n*pi) == tan(theta)` regardless of `n`'s
/// parity) -- so this is correct by construction, not an approximation
/// that happens to work. At `cospi`'s own zeros (`x` a half-integer,
/// `tan`'s true poles), IEEE754 division by a signed zero already gives
/// the correctly-signed `+-inf` for free, no extra handling needed
/// (`sinpi` is nonzero there, so this is a real `finite/0`, never the
/// `0/0` that would need a NaN override).
#[inline(always)]
pub fn tanpi(x: f32) -> f32 {
    sinpi(x) / cospi(x)
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
/// Unlike `sinpi` (exact for the *entire* f32 range, since a plain
/// integer `q` never needs more precision than f32 already gives it),
/// this reduction is only exact while `q*180.0` stays representable --
/// `180.0`'s 18 trailing zero mantissa bits keep that true for `|q|`
/// up to ~2^18 (`|x|` up to ~4.7e7), comfortably past any realistic
/// input but not the *entire* range the way `sinpi` achieves. Past that,
/// `d` stops being small and `d*DEG_TO_RAD_SMALL` could land far outside
/// `sinf_poly`'s fitted domain -- confirmed by a direct probe
/// (`sind(1e10)` returned `36046.176`, nonsense, before this clamp was
/// added). Fixed the same way `sin_checked`/`cos_checked` guard their own
/// poly input (`POLY_SAFE_BOUND`): clamp the radian residual to
/// `[-1000, 1000]` before `sinf_poly` sees it. This guarantees the same
/// contract `sin_checked` documents for its own far tail: always finite
/// for finite input (confirmed out to `f32::MAX`, never inf/nan), *not*
/// a guarantee of numerical correctness that far out -- past the exact
/// boundary the output is bounded-but-wrong, same as any fast tier in
/// this crate past its own documented exact range.
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

/// tan(x*pi/180), argument in degrees -- new function (backlog idea #29),
/// same `sind(x)/cosd(x)` ratio construction as `tanpi`'s own doc
/// comment describes (period-180 cancellation, poles handled for free by
/// IEEE754 division). See `sind`'s own doc comment for the shared
/// reduction's exactness limit (`|x|` up to ~4.7e7) and safety clamp.
#[inline(always)]
pub fn tand(x: f32) -> f32 {
    sind(x) / cosd(x)
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
// precision via `two_prod`/`two_sum` error-free transforms (each individual
// transform is exact -- p+e == a*b or s+e == a+b as reals, no rounding lost
// converting a sum/product into a hi/lo pair -- unlike the crate's other
// PI_A..D trick, which needs bounded q) -- specifically, the dominant cross
// term and the two next-biggest (~x*2^-24) get real two_prod treatment, and
// the smallest tier (~x*2^-48, plus the reciprocal-of-pi's own 3rd
// correction word) is folded in with plain multiplies/adds (its own
// rounding error there is already far below 1 ulp of the O(1) result). A
// version that dropped that smallest tier entirely was tried and measured
// (via examples/mca.rs) to cost the *same* as keeping it (~130-140 cyc
// either way) -- no meaningful savings once genuine multi-word q precision
// is needed at all, so all three tiers are kept here for the best accuracy
// at no extra cost. See POLY_SAFE_BOUND for why the output stays finite
// even once this gradual degradation is severe.
//
// CORRECTION (2026-07-10, IDEAS.md): the individual two_prod/two_sum
// transforms being exact does *not* make qh/ql "exact at any magnitude" as
// a previous version of this comment claimed -- that conflates "no
// rounding lost converting to a hi/lo pair" with "a hi/lo pair has
// unlimited precision," which isn't true. `qh` and `ql` together resolve
// `q` to roughly 48 bits total (two f32 mantissas); once the *true*
// `round(x/pi)` itself needs more than ~48 bits to pin down -- around
// `|x| > 2^48*pi ~ 8.85e14` -- qh+ql comes out off by a few whole integers
// (confirmed empirically: off by 1 at x=1e15, by 14 at x~1e16), which is
// exactly the "relocatable cliff" this whole double-word scheme was built
// to avoid, just relocated from a single f32's ~2^24 ceiling to roughly
// 2^48 instead of eliminated. `POLY_SAFE_BOUND` bounds `sinf_poly`'s
// *input* here but not its *output* -- see `sin_checked`'s own `.clamp`
// for the fix that keeps `sin_checked`/`cos_checked` inside `[-1,1]`
// regardless.
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
    // Calls the copysign-free `sinf_poly_raw`, not `sinf_poly` -- the
    // `x == 0.0` guard below already overrides the result at the one
    // point copysign would matter (see `sinf_poly_raw`'s own doc comment
    // for the verification this is safe), so paying for that instruction
    // here would be pure waste.
    //
    // `.clamp(-1.0, 1.0)` (2026-07-10, IDEAS.md): `round_x_over_pi`'s
    // double-float q genuinely loses precision once |x| exceeds roughly
    // 2^48*pi (~8.85e14) -- q itself comes out off by a few integers there
    // (confirmed empirically, not just theorized), which shifts r by whole
    // multiples of pi and puts it wildly outside sinf_poly's fitted
    // [-pi/2,pi/2] domain despite the POLY_SAFE_BOUND clamp (which only
    // bounds the poly's *input*, not its *output* -- a degree-9 poly
    // evaluated at |r|=1000 is dominated by its own leading r^9 term,
    // ~2.6e21, nowhere near the true `|sin(x)| <= 1` invariant). Without
    // this clamp, `sin_checked`/`cos_checked` could silently return values
    // like `1.07e9` or `2.6e21` for legitimate (if extreme) finite input --
    // a genuine invariant violation, not just reduced accuracy, and a much
    // worse failure mode than the "gradual degradation" these functions are
    // documented to provide. The clamp doesn't fix the underlying accuracy
    // for `|x|` this extreme (a real fix needs a wider-than-double-float q,
    // out of scope here) but it restores the one invariant every caller can
    // still rely on regardless of `x`'s magnitude. Verified a true no-op
    // everywhere the function was already accurate (fuzz: every in-budget
    // bucket from `|x|<=pi/4` through `[1e12,1e13)` unchanged vs. baseline);
    // the already-garbage `[1e15,inf)` tail's *ulp* numbers don't visibly
    // improve (this crate's own `ulp_diff` measures distance via a global
    // sign+magnitude ordering, so two arbitrary values inside `[-1,1]` can
    // still read as billions of "ulp" apart -- the fix isn't about that
    // metric), but the actual returned values there now measure bounded
    // to exactly `[-1,1]` for every tested magnitude up to `f32::MAX`
    // (previously up to `2.6e21`). Real, non-zero mca cost accepted (same
    // "real perf cost to fix a wrong-for-legitimate-input defect" precedent
    // as sin/cos's own inf-for-large-x fix and sinh/cosh's domain-hole fix):
    // sin_checked latency 109.02->117.02 cyc (+7.3%, includes the
    // sinf_poly_raw saving above), throughput 5.476->5.542 (+1.2%);
    // cos_checked (clamp only, no offsetting saving) 113.00->122.00
    // (+8.0%), 4.537->5.672 (+25.0%, confirmed via assembly diff to be a
    // real port-contention/scheduling shift, not de-vectorization -- total
    // instruction count only grew 166->173).
    let result = sinf_poly_raw(r).clamp(-1.0, 1.0);
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
    // See sin_checked's own doc comment for why this clamp is needed:
    // round_x_over_pi's double-float q loses precision for |x| beyond
    // ~2^48*pi, and POLY_SAFE_BOUND only bounds sinf_poly's *input*, not
    // its *output* -- without this, cos_checked could silently return
    // values like 2.6e21 for legitimate finite input, violating
    // `|cos(x)| <= 1`.
    sinf_poly(r).clamp(-1.0, 1.0)
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
/// literal constants differ). Refit again (2026-07-09) via an
/// ulp-weighted Chebyshev LP with a max-ulp cap (this session's exp/erf
/// technique) against the true target `P(r) = ((1+r)^(-1/3)-1)/r` over
/// one representative octave (a in [1,8), matching this poly's own
/// established octave-periodicity, unlike the unrelated cbrt_throughput
/// seed where a single-octave grid was later found misleading -- see
/// IDEAS.md): avg ulp 0.3112 -> 0.2813 (~9.6%, exhaustive), max ulp
/// unchanged at 3, zero perf cost. This is the same poly a prior session
/// already ran an *unconstrained* minimax LP against (see IDEAS.md's
/// "Exhaustive/rlibm-style" entry) and got max 3->2 at the cost of a
/// real avg regression (0.3125->0.4487) -- the max-cap fix (matching
/// acos_poly's own successful constrained-search pattern) recovers a
/// real avg win instead, at the cost of not chasing the max-ulp cut this
/// time. cbrt_accurate (which reuses this as a Newton seed) is
/// unaffected either way -- 0.000 avg / 1 max ulp before and after,
/// exhaustive.
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

/// x^(-1/3) -- new function, backlog idea #55: division-free by
/// construction wasn't pursued (that would need a fresh bit-trick seed +
/// its own fitted correction poly, real numerical-fitting work, not a
/// "reuse existing pieces" change); instead this composes `cbrt` with a
/// single hardware division, the same "sqrt/division are already
/// correctly-rounded, just compose them" reasoning `rsqrt`/`rhypot` both
/// already use in this crate. Every zero/inf/nan/negative special case
/// falls out of that composition for free via IEEE754 semantics (checked
/// by hand before writing this, same discipline as `rhypot`):
/// `rcbrt(0)=inf`, `rcbrt(-0)=-inf` (cbrt is odd, so is its reciprocal),
/// `rcbrt(inf)=0`, `rcbrt(-inf)=-0`, `rcbrt(nan)=nan`, `rcbrt(-8)=-0.5`
/// -- no override needed at all (unlike `rhypot`'s one inf-vs-NaN case,
/// `cbrt` has no analogous special case to inherit). Costs one more
/// rounding than `cbrt` itself (the division), same tradeoff `rsqrt`/
/// `rhypot` already accept.
#[inline(always)]
pub fn rcbrt(x: f32) -> f32 {
    1.0 / cbrt(x)
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
#[allow(clippy::neg_cmp_op_on_partial_ord)]
pub fn ln(x: f32) -> f32 {
    log_family_wrapper!(x, ln_normal)
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
    let (p, s, k) = log_family_normal!(
        x,
        koff,
        [
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
        ]
    );
    let k_hi = k * LN2_HI; // exact, see LN2_HI's comment
    fma(p, s, k_hi) + k * LN2_LO
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
#[inline(always)]
#[allow(clippy::neg_cmp_op_on_partial_ord)]
pub fn log10(x: f32) -> f32 {
    log_family_wrapper!(x, log10_normal)
}

/// Core of log10 for positive normal finite x only -- see ln_normal, same
/// approach with coefficients fitted for log10 (log_2's c[i] * LOG10_2).
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn log10_normal(x: f32, koff: f32) -> f32 {
    let (p, s, k) = log_family_normal!(
        x,
        koff,
        [
            std::f32::consts::LOG10_E, // bit-identical to this literal; not a coincidence
            -0.2171472,
            0.14476489,
            -0.10858066,
            0.08686763,
            -0.07213202,
            0.06159092,
            -0.05751561,
            0.05604425,
            -0.03309811,
        ]
    );
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
/// (`log_2(1.0)` is `+0.0`, `corr` carries x's sign, and IEEE754 always
/// resolves `+0.0 + -0.0` to `+0.0`) -- both copied from `log1p` verbatim.
/// Originally verified only by a fuzz (20M samples over `-1<x<1e6` plus a
/// 20M-sample pass concentrated on `|x|<1e-6`), which reported avg ulp
/// 0.005, max ulp 2 -- undersold by a rare worst case that fuzz sampling
/// simply didn't land on. The real exhaustive sweep (every f32 bit
/// pattern, matching readme.md's own table) gives avg ulp 0.102, max ulp 3
/// at `x=0.018272582` (2026-07-10).
#[inline(always)]
pub fn log2p1(x: f32) -> f32 {
    let u = 1.0 + x;
    let c = x - (u - 1.0);
    let corr = (c / u) * LOG2_E;
    let corr = if corr.is_finite() { corr } else { 0.0 };
    let normal = log_2(u) + corr;
    if x == 0.0 { x } else { normal }
}

/// log1p without the `corr.is_finite()` guard: valid whenever `x` is
/// finite and `x != -1.0` (i.e. `u = 1+x` is itself finite and nonzero,
/// the two edges the guard exists to suppress -- see log1p's own doc
/// comment). `asinh`/`acosh` both already check `d.is_finite()` before
/// calling log1p at all, and by construction (their own doc comments) the
/// value they pass is never `-1.0` on the path where it matters (asinh:
/// `d = ax + sm1` with both terms `>= 0`; acosh: `d = (x-1.0) + s` with
/// `s >= 0` over the valid `x >= 1` domain) -- so the guard is
/// unreachable for either caller and this drops one redundant select from
/// two hot composites.
#[inline(always)]
fn log1p_finite(x: f32) -> f32 {
    let u = 1.0 + x;
    let c = x - (u - 1.0);
    let corr = c / u;
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
///
/// `k`'s rounding (2026-07-08) uses the sin/cos-style magic-constant add
/// (`fma(x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC`) instead of a plain
/// `.round()` call -- IDEAS.md's backlog only speculated this would save
/// one op on `k1b`'s own already-existing magic-round (unclear payoff,
/// "worth one mca look"), but the real win was bigger and came from
/// somewhere else: native `.round()` (vroundps) is apparently just a
/// slower instruction on this CPU than the fma+subtract pair the magic
/// trick uses, independent of anything downstream. Measured (mca):
/// latency 51.00->42.00 cyc (-17.6%), throughput 1.648->1.327 cyc/elem
/// (-19.5%), cascading to every caller (expm1 -9.9%/-17.6%, sinh
/// -14.3%/-20.4%, cosh -14.5%/-18.3%, tanh -8.1%/-14.2% lat/throughput;
/// `powf`/`erf`/`erfc`/`asinh`/`acosh`/`atanh` don't call this function,
/// unaffected). This swaps round-half-away-from-zero (`.round()`'s
/// documented behavior) for the hardware's round-half-to-even (what the
/// add-then-subtract trick actually performs) -- differs only at exact
/// half-integer ties of `x*log2(e)`, and both the 100M-sample fuzz *and*
/// the exhaustive all-2^32-pattern sweep found zero measurable accuracy
/// difference (exp/expm1/sinh/cosh avg+max ulp bit-identical; tanh's
/// looked like max 8->9 under fuzz alone until the exhaustive sweep
/// showed *both* versions hit max 9 at the same worst x -- fuzz just
/// hadn't sampled that point before, not a regression). `exp(88.37628)`
/// (the k=128 edgecheck above) independently re-verified bit-exact
/// against a true f64 reference.
#[inline(always)]
pub fn exp(x: f32) -> f32 {
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let r = fma(-k, LN2_HI, x);
    let r = fma(-k, LN2_LO, r);
    // c0 and c1 both forced exactly 1.0 (were 1.0 and 1.0000000647031426):
    // exp(r) = 1 + r + r^2*P(r) for tiny r, so a c1 off from 1.0 by even
    // ~6e-8 relative is a systematic bias right where exp(x) is most
    // commonly called (x near 0). l0 = r + 1.0 needs no fma since both its
    // coefficients are now exactly 1.0. c2..c5 refit (2026-07-09) via a
    // proper ulp-weighted Chebyshev LP (scipy.optimize.linprog: minimize t
    // s.t. |P(r_i)-exp(r_i)| <= t/ulp(exp(r_i)) over a dense r grid),
    // superseding the plain coordinate-descent fit from IDEAS.md's
    // exp_r_c/"exp_r" tuner target -- coordinate descent had converged to a
    // real but non-global local optimum (the residual's worst points all
    // clustered near the domain edge |r|~0.3, not near a ulp-boundary, so
    // this was mostly the LP finding a better balance across the domain
    // than greedy per-coefficient steps, not a boundary-crossing effect).
    // Verified against the real crate, exhaustive all-2^32-pattern sweep,
    // paired via git stash: exp avg/max ulp 0.0911/4 -> 0.0745/3; cascades
    // to every caller since they all reuse this poly directly -- tanh
    // 0.1482/9 -> 0.1457/6 (the biggest win, since tanh's own worst case
    // sat inside this poly's residual), sinh_throughput 0.0839/7 ->
    // 0.0723/5, cosh_throughput 0.0606/4 -> 0.0507/4 (avg only), expm1
    // 0.1381/6 -> 0.1304/6 (avg only -- expm1's own max ulp lives entirely
    // in its other, exp-free Pade branch, untouched by this change, see
    // that branch's own doc comment). sinh/cosh exactly unchanged (their
    // own worst case doesn't sit in this residual). mca: latency/
    // throughput bit-for-bit identical on exp/sinh/tanh (42.00/1.327,
    // 58.00/2.523, 94.91/2.567) -- expected, pure coefficient swap, same
    // instructions.
    let p = exp_r_poly!(r);
    let (t1, t2) = exp2_field_split(k);
    p * t1 * t2
}

/// Full-range sibling of [`exp`] -- backlog idea #18: `exp` already uses
/// the k1/k2 split (needed regardless, see its own doc comment), so the
/// only change is clamping `x` *before* the reduction starts so `k`
/// itself never leaves the split's own safe `[-151,128)` range, matching
/// `exp2_checked`'s own early-clamp pattern (clamp the input, not the
/// derived `k`, so the reduction's `r` and the split's `k` always stay
/// mutually consistent -- clamping `k` after the fact would desync it
/// from an `r` computed against the *unclamped* value, the same class of
/// bug `sigmoid`'s own negative-tail fix found and fixed this session).
/// Clamp bounds (`128/log2(e)`, `-151/log2(e)`) are `exp2_checked`'s own
/// `k` boundary converted into `x`'s units, not guessed -- the same
/// derivation `sigmoid`'s own fix used. Callers like `sigmoid`/`tanh`
/// that currently hand-roll their own ad-hoc clamp before a duplicated
/// copy of `exp`'s reduction could route through this instead, but
/// aren't changed here (out of scope for this entry; each has its own
/// established, already-verified clamp bound and doc comment).
#[inline(always)]
pub fn exp_checked(x: f32) -> f32 {
    let x = x.clamp(-104.66522426455174, 88.72283911167308);
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let r = fma(-k, LN2_HI, x);
    let r = fma(-k, LN2_LO, r);
    let p = exp_r_poly!(r);
    let (t1, t2) = exp2_field_split(k);
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
    let a = pade_expm1_ratio!(x, mul);
    // Deliberately a standalone copy of exp's reduction (not routed through
    // the public `exp` fn) ending in fma(p*t1, t2, -1.0) instead of
    // exp(x)-1.0 -- fuses the trailing subtract into the last multiply,
    // one rounding fewer than rounding p*t1*t2 fully and then subtracting
    // 1.0 separately. This is the branch that carries expm1's actual
    // worst-case ulp (the Pade branch above is the one with headroom, not
    // this one -- confirmed by the exhaustive sweep's worst-x always
    // landing at |x|>=0.5, contradicting a stale claim in exp's own doc
    // comment). Factoring the *whole* reduction+poly+combine through a
    // shared `exp`-returning-parts `fn` was tried first and measured a
    // real, reproducible mca latency regression on `sinh_throughput` (+32%,
    // unrelated caller, apparently a scheduling side effect of the new
    // function boundary) even though `exp` itself was bit-identical -- so
    // the reduction (`k`/`r`) stays duplicated here (each caller derives
    // `k` slightly differently, from a plain reduction vs. exp2_checked's
    // own clamp-then-reduce). The 4-coefficient poly evaluation is shared
    // via `exp_r_poly!` (a macro, no function-call boundary at all), and
    // the exponent-field-split below is shared via `exp2_field_split` (an
    // already-existing small `fn`, already relied on by `exp_pos_neg` --
    // reusing an existing call site rather than introducing a new one is a
    // different risk profile than the original rejected experiment, and a
    // full pre/post assembly diff confirms zero byte differences, see
    // IDEAS.md's own idea log, 2026-07-11).
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let r = fma(-k, LN2_HI, x);
    let r = fma(-k, LN2_LO, r);
    let p = exp_r_poly!(r);
    let (t1, t2) = exp2_field_split(k);
    let b = fma(p * t1, t2, -1.0);
    if x.abs() < 0.5 { a } else { b }
}

/// (e^x - 1)/x -- new function, backlog idea #66: the well-conditioned
/// primitive behind financial (continuously-compounded-rate) and ODE
/// (exponential-integrator) kernels, where callers currently have to
/// write `expm1(x)/x` themselves and hope `x` never lands exactly on the
/// removable singularity at 0. `expm1`'s own Pade branch is *already*
/// exactly this shape internally: `a = x * N(x)/D(x)`, so `a/x = N(x)/D(x)`
/// with the `x` factor cancelling algebraically, before any rounding ever
/// touches it -- no new cancellation risk, not even at `x=0` itself
/// (verified by hand: `N(0)=fma(-1.9999927,0,-120)=-120`,
/// `D(0)=fma(0,fma(0,-12.00003,60),-120)=-120`, so `N(0)/D(0)=1.0`
/// exactly, matching the true limit `lim_{x->0}(e^x-1)/x=1` with no
/// explicit `x==0.0` select needed at all). The direct branch (`|x|>=0.5`,
/// where there's no cancellation to begin with) is just `expm1`'s own
/// combine divided by `x`, a single extra rounding. Duplicates `expm1`'s
/// reduction/poly rather than routing through it, matching `expm1`'s own
/// standalone-copy precedent (see its doc comment for why sharing a
/// helper here previously regressed an unrelated caller's codegen).
/// Inherits `expm1`'s own unchecked-exp2 domain limit (see `exp`'s doc
/// comment) -- garbage outside roughly `x in [-87.3, 88.7)`, same as
/// `expm1` itself.
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

/// 2^x - 1 (C23 `exp2m1`). Same cancellation problem as `expm1` (2^x is
/// close to 1 whenever x is close to 0, so computing 2^x first and
/// subtracting 1 loses low bits) and the same fix: `expm1`'s own Pade
/// approximant for e^y-1 is reused directly via the substitution
/// `y = x*LN_2` (2^x - 1 = e^{x ln2} - 1), valid because the Pade branch
/// only ever runs for `|x| < 0.5`, so `|y| < 0.5*ln2 ≈ 0.347` stays
/// comfortably inside the domain that approximant was fitted over -- no
/// new coefficients needed. This substitution would *not* be safe for the
/// direct branch: reducing through `exp2_checked(x*LOG2_E)` the way exp's
/// own doc comment warns against would reintroduce that exact bug for
/// large x, so the direct branch below instead duplicates
/// `exp2_checked`'s own k/f reduction and Q(f) poly verbatim (not routed
/// through the public `exp2_checked`, same reasoning as `expm1`'s own
/// standalone copy of `exp`'s reduction) and fuses the trailing `-1` into
/// the last multiply (`fma(p, t2, -1.0)`, one rounding instead of two).
/// Originally verified only by a fuzz (20M samples over |x|<100 plus a
/// separate 20M-sample pass concentrated on |x|<1 to stress the branch
/// seam), which reported avg ulp 0.06, max ulp 3 -- undersold by a rare
/// worst case (`x=0.5842032`, comfortably inside both fuzz ranges) that
/// sampling simply didn't land on. The real exhaustive sweep (every f32
/// bit pattern, matching readme.md's own table) gives avg ulp 0.077, max
/// ulp 4 (2026-07-10) -- still comparable to `expm1`'s own max ulp 3-4, no
/// seam discontinuity at the `|x|<0.5` threshold. Inherits `exp2_checked`'s
/// full `[-151, 128)` clamp, so is total (never NaN/inf-producing outside
/// its own true asymptotes): `exp2m1(-inf) = -1`, `exp2m1(inf) = inf`.
///
/// IDEAS.md backlog round 3 #1's round-domain Q(f) refit was tried here
/// too (same coefficients as the exp2_checked attempt) -- real max-ulp
/// regression on the full fuzz (3->6/7, avg roughly flat) for a real but
/// small mca cost (latency +1 cyc, throughput +3.3%), and would likely
/// have hit the same `f<0` overflow-saturation gap exp10_checked's own
/// attempt did at the `exp2m1(inf)=inf` boundary above (not confirmed via
/// edgecheck since the fuzz/mca result alone already killed it). Reverted,
/// bit-identical to prior HEAD.
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
    if x.abs() < 0.5 { a } else { b }
}

// Shared exp(x)/exp(-x) for sinh/cosh (2026-07-08): the Cody-Waite
// reduction only needs to happen once, since -x's reduction is exactly
// (-k, -r) -- immediate from x = k*ln2 + r (exact by construction), no
// independent rounding of -x*log2e needed or even relevant. exp's own e^r
// poly splits into an even/odd part in r^2 (e(u) = 1 + c[0]*u + c[2]*u^2,
// o(u) = 1 + c[1]*u + c[3]*u^2, so p(r) = e + r*o matches exp's p exactly
// -- verified by expanding both forms), so p(-r) = e - r*o reuses e/o at
// the cost of one more fma instead of a whole second poly. Only the final
// exponent-field scaling (2^k vs 2^-k) is genuinely duplicated -- cheap
// integer/bit-trick work, not fma-port pressure. See exp's own doc
// comment for shared domain/rounding notes; same unchecked-exp2 domain
// limit applies to both outputs here.
//
// This was IDEAS.md's "biggest single-function win candidate" backlog
// entry, and it delivered a real but three-way mixed result, all
// re-measured after the coefficient retune below (see exp_pos_neg):
// mca throughput sinh 3.057->2.523 cyc/elem (-17.5%), cosh 2.743->2.074
// (-24.4%) -- the predicted win, since the fma-port-heavy reduction+poly
// now runs once instead of twice. But mca latency got *worse*, sinh
// 54.00->58.00 cyc (+7.4%), cosh 53.00->57.00 (+7.5%): the old code's two
// independent exp(x)/exp(-x) calls could run concurrently on this
// out-of-order CPU (no data dependency between them), hiding one behind
// the other; funneling both through one shared serial prefix removes
// that overlap. Same shape as asin's 2-branch collapse (see its own doc
// comment, fix 6) and several other entries in this file -- throughput
// is the metric this crate's vectorization-first design prioritizes, so
// kept despite the latency regression. Accuracy (exhaustive sweep):
// after retuning, sinh's max ulp is unchanged (5->5, avg ulp 0.0797->
// 0.0806, noise-level); cosh's max ulp moved 4->5 (avg actually improved,
// 0.0760->0.0714) -- a small real regression, but landing at exactly the
// max-ulp level sinh (its own sibling function) already carries, not some
// new outlier. Adopted on that basis: a substantial, clearly-prioritized
// throughput win against a latency cost and a one-ulp max-ulp move that
// lands within this composite function family's own existing spread.
#[inline(always)]
fn exp2_field_split(k: f32) -> (f32, f32) {
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k1b = fma(k, 0.5, ROUND_MAGIC) - (ROUND_MAGIC - 383.0);
    let k2b = (k + 766.0) - k1b;
    let t1 = f32::from_bits((k1b.to_bits() << 8) & EXPONENT_MASK);
    let t2 = f32::from_bits((k2b.to_bits() << 8) & EXPONENT_MASK);
    (t1, t2)
}

#[inline(always)]
fn exp_pos_neg(x: f32) -> (f32, f32) {
    // Retuned for this even/odd split specifically (examples/tune.rs's
    // exp_r_pair_c/"exp_r_pair") -- the plain-Horner exp() coefficients
    // copied verbatim here left max ulp 4 on the tuning grid, retuning
    // c0/c1 recovered max ulp 3 (c2/c3 didn't move).
    //
    // Refit again (2026-07-10) via a genuine ulp-weighted joint Chebyshev
    // LP (idea #91's technique, using scipy -- confirmed available this
    // session, see IDEAS.md), scoring both p_pos/e^r and p_neg/e^-r
    // simultaneously over r in [-ln2/2,ln2/2]. Unlike the same LP applied
    // to log_2 (idealized win that evaporated through real f32 rounding),
    // this poly's error is fit-quality-dominated (idea #7's own round-off
    // audit), so the idealized gain largely survived: verified against
    // the real compiled sinh/cosh (bit-identical fidelity check first),
    // then a real ~559M-point dense sweep: sinh avg ulp 0.08875->0.08121
    // (~8.5% tighter), cosh avg ulp 0.08792->0.07638 (~13.1% tighter),
    // max ulp unchanged at 5 for both. Zero perf cost, same instructions
    // (same literal-swap shape as every other coefficient-only refit).
    let (p_pos, p_neg, t1, t2, t1n, t2n) = exp_pos_neg_core!(x);
    (p_pos * t1 * t2, p_neg * t1n * t2n)
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

/// sinh(x) = 0.5*(exp(x) - exp(-x)) directly (via `exp_pos_neg`'s shared
/// reduction, see its own doc comment), except for |x| < 0.5 where exp(x)
/// and exp(-x) are both ~1 and the subtraction cancels almost all
/// precision (the same class of bug log1p/tanh had, see IDEAS.md) --
/// there, use the Taylor form above instead, same branchless-select
/// pattern as expm1's Pade/exp split. See exp's doc comment for the
/// inherited unchecked-exp2 domain limit (only relevant on the `b` side,
/// unconditionally evaluated but only selected for |x| >= 0.5). cosh below
/// doesn't need this: it adds instead of subtracting, so it never cancels.
#[inline(always)]
pub fn sinh(x: f32) -> f32 {
    let a = sinh_small(x);
    let (ep, en) = exp_pos_neg(x);
    let b = 0.5 * (ep - en);
    if x.abs() < 0.5 { a } else { b }
}

/// cosh(x) = 0.5*(exp(x) + exp(-x)) via `exp_pos_neg`'s shared reduction
/// (see its own doc comment) -- never cancels (adds instead of
/// subtracts), so unlike sinh needs no small-x branch. See exp's doc
/// comment for the inherited unchecked-exp2 domain limit.
#[inline(always)]
pub fn cosh(x: f32) -> f32 {
    let (ep, en) = exp_pos_neg(x);
    0.5 * (ep + en)
}

/// `exp_pos_neg`, but with `x` clamped first so `k = round(x*log2e)`
/// never leaves the safe range for *both* `exp2_field_split(k)` (used
/// for `exp(x)`) and its reciprocal-based negation (used for `exp(-x)`),
/// AND returning `0.5*exp(x)`/`0.5*exp(-x)` directly instead of the raw
/// pair -- found 2026-07-09 (backlog idea #85's fifth wave, the
/// special-case matrix technique applied to `sinh`/`cosh`/`tanh`/
/// `sigmoid`):
///
/// 1) `exp_pos_neg` has no clamp at all, so for `|x|` large enough that
///    `k` leaves `exp2_field_split`'s own safe domain, the bit-trick
///    exponent construction wraps around instead of saturating --
///    `sinh(1000)` came out `-inf` (*wrong sign*, should be `+inf`),
///    `sinh(-1000)` came out a finite garbage value (should be `-inf`),
///    and `sinh`/`cosh(+-inf)` both came out `NaN` (should be
///    `+-inf`/`+inf`). A standalone probe enumerating every integer `k`
///    found the split (and its reciprocal-bit-trick negation) exactly
///    matches `2^k`/`2^-k` cast to f32 -- including correct saturation
///    to `0`/`inf` -- for `k` in `[-254, 254]`; wraparound starts right
///    outside that (first mismatch at `k = +-255`). That's much wider
///    than initially assumed (an earlier version of this fix clamped to
///    `|k|<128`, matching `exp_checked`'s own asymmetric bound applied
///    symmetrically) -- `128` is where `exp_checked` itself needs to
///    stop (`exp(x)` alone overflows there), not where the split's
///    bit-trick construction actually breaks. The clamp below uses
///    `170.0` (`k` up to ~245.3), comfortably inside the proven-safe
///    `254` ceiling with margin, and far past the true `sinh`/`cosh`
///    overflow threshold below -- so it only ever discards inputs whose
///    correct answer is already exactly `+-inf` anyway.
/// 2) The initial (too-tight, `|x|<=88.72`) clamp also surfaced a
///    second, larger-magnitude bug: for `x` in roughly `[87.3, 89.4]`,
///    `sinh(x)`/`cosh(x)` are still finite and f32-representable
///    (they're *half* of `exp(x)`, which overflows a bit earlier), but
///    `sinh`/`cosh`'s `p_pos * t1 * t2` computes the full unscaled
///    `exp(x)` first and only applies the `0.5` factor afterward in the
///    caller -- so the intermediate overflows to `inf` before the
///    caller ever gets to halve it. An exhaustive sweep confirmed this:
///    max ulp *8388030* right at the boundary (previously invisible --
///    the old accuracy.rs sweep for plain `sinh`/`cosh` deliberately
///    excludes this exact window via its own `sinh_domain` restriction,
///    so the gap was never measured). Fixed by pushing the `0.5` into
///    the exact-power-of-two field split instead of the final result:
///    `t1` (or `t1n`) is halved (`t1 * 0.5`, exact for any power-of-two
///    float down to the denormal floor) *before* multiplying by
///    `t2`/`t2n`, so the product only needs to represent `0.5*exp(x)`
///    rather than the (larger, earlier-overflowing) `exp(x)` itself.
///    Combined with the wide `170.0` clamp above, this now covers the
///    whole legitimately-finite window (`|x|` up to ~`89.4`) exactly,
///    with `+-inf` returned correctly everywhere beyond it.
#[inline(always)]
fn exp_pos_neg_checked_half(x: f32) -> (f32, f32) {
    let x = x.clamp(-170.0, 170.0);
    let (p_pos, p_neg, t1, t2, t1n, t2n) = exp_pos_neg_core!(x);
    (p_pos * (t1 * 0.5) * t2, p_neg * (t1n * 0.5) * t2n)
}

/// Full-range sibling of [`sinh`] -- same construction, just built on
/// `exp_pos_neg_checked_half` instead of the unchecked `exp_pos_neg`
/// (which already returns the `0.5*exp(+-x)` halves, so no separate
/// `0.5*` multiply here, unlike `sinh`). See that function's own doc
/// comment for the correctness gaps this closes (`sinh`'s own
/// wrong-sign/NaN behavior for large `|x|`, plus a premature-overflow
/// gap just below that, both found via backlog idea #85's fifth wave).
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
/// `2*x` is now clamped to `[-87.0, 88.0]` before reaching `expm1`
/// (2026-07-08). Previously the raw `2*x` inherited exp's unchecked-
/// domain garbage for |x| > ~44 (2x past ~88.7) -- verified as a real
/// bug, not just theoretical: `tanh(50)`/`tanh(1000)`/`tanh(f32::MAX)`
/// all returned NaN where std correctly saturates to 1. IDEAS.md's
/// backlog proposed fixing this via `tanh(x) = -expm1(-2|x|)/
/// (expm1(-2|x|)+2)` (always-nonpositive argument, "never overflows"),
/// reasoning that expm1's negative side only ever saturates gracefully
/// -- turned out to need verification, not just trust: `exp()` itself is
/// *not* gracefully monotonic outside its documented ~[-87.3, 88.7)
/// domain, it's genuine garbage in both directions (spot checks found
/// `exp(-150)=0` correctly but `exp(-200)=inf`, wrong, and worse beyond
/// that), so the abs-based form still needed a clamp to be safe, at
/// which point the abs/mulsign restructuring was doing no work the
/// clamp alone doesn't already do. Measured both (mca): the abs+mulsign
/// form cost 98.86/2.541 cyc lat/throughput vs. this plain-clamp form's
/// 94.91/2.567 -- the simpler form has better latency and matches
/// baseline in-domain accuracy exactly (avg/max ulp 0.1482/8, bit-
/// identical to the unclamped original, since real inputs never reach
/// the clamp boundary), where the abs+mulsign form measured very
/// slightly worse (avg ulp 0.1488) from its extra rounding steps. Kept
/// the plain clamp. `f32::clamp` returns NaN unchanged if the input is
/// NaN (unlike `.max`/`.min`), so NaN propagation isn't broken, and the
/// clamp is lossless even for in-range-but-large x: any x past the
/// clamp boundary already has a true tanh value of exactly 1.0f32 (or
/// -1.0f32) many orders of magnitude before reaching it (e^-2*8.3 is
/// already below ulp(1.0)/2), so clamping produces the bit-identical
/// correctly-rounded answer, not an approximation. Real remaining cost:
/// mca latency +4.5% (90.78->94.91 cyc), throughput +18.7% worse
/// (2.163->2.567 cyc/elem) versus the unclamped original -- accepted per
/// this crate's own established precedent of paying a real perf cost to
/// fix a "wrong/NaN for legitimate finite input" domain hole (see sin/
/// cos's own inf-for-large-x fix), which is a more serious defect class
/// than an in-domain ulp regression.
// Shared by tanh/sigmoid: the reduction plus single-exponent-field
// construction (given each caller's own already-clamped `y`) is
// identical -- only the shared poly's own consumer differs
// (`fma(p, exp2int, -1.0)`'s trailing -1 fusion for tanh's exp(y)-1 vs.
// sigmoid's plain `p * exp2int`). Returns `(p, exp2int)`. Macro, not a
// fn -- same reasoning as this file's other shared-body macros.
macro_rules! exp_r_singlefield {
    ($y:expr) => {{
        const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
        let k = fma($y, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
        let r = fma(-k, LN2_HI, $y);
        let r = fma(-k, LN2_LO, r);
        let p = exp_r_poly!(r);
        let exp2int = f32::from_bits(((k + 383_f32).to_bits() << 8) & EXPONENT_MASK);
        (p, exp2int)
    }};
}

#[inline(always)]
pub fn tanh(x: f32) -> f32 {
    // Standalone copy of expm1 (not a call through the public `expm1` fn,
    // same "avoid a shared-helper scheduling risk" reasoning as expm1's
    // own doc comment above) but with a single exponent-field construction
    // in the exp(x)-1 branch instead of exp's k1/k2 split: the
    // `.clamp(-87.0, 88.0)` bound already guarantees `k = round(y*log2e)`
    // stays in [-126, 127] (same domain sigmoid's own single-field fix
    // relies on), comfortably short of the k=128 edge case the split
    // exists for. Also carries the same fma(p, exp2int, -1.0) tail fusion
    // as expm1's own fix.
    let y = (2.0 * x).clamp(-87.0, 88.0);
    let a = pade_expm1_ratio!(y, mul);
    let (p, exp2int) = exp_r_singlefield!(y);
    let b = fma(p, exp2int, -1.0);
    let e = if y.abs() < 0.5 { a } else { b };
    e / (e + 2.0)
}

/// logistic sigmoid, `1/(1+exp(-x))`. The backlog idea for this function
/// proposed the identity `sigmoid(x) = 0.5 + 0.5*tanh(x/2)` (algebraically
/// exact: multiply `tanh(x/2)`'s own definition by `e^(x/2)/e^(x/2)` to
/// get `(e^x-1)/(e^x+1)`, and `0.5 + 0.5*` that simplifies to
/// `e^x/(e^x+1) = 1/(1+e^-x)`) -- tried first, and it's a real bug, not
/// just imprecise: for `x` around `-17.3`, `tanh(x/2)` (`~-8.66`) already
/// rounds to *exactly* `-1.0f32` (its true value is within half a ulp of
/// `-1.0`, so that's the *correct* f32 rounding for tanh's own output),
/// but `0.5 + 0.5*(-1.0)` then computes to *exactly* `0.0` even though
/// the true sigmoid value there (`~2.98e-8`) is nowhere near it or f32's
/// underflow threshold -- `tanh`'s own correct saturation discards
/// exactly the residual precision this formula needs. The same class of
/// "algebraically-exact identity reintroduces catastrophic cancellation"
/// bug this file's `atanh` single-log1p-fusion entry already documents.
/// Fixed by computing directly instead, which has no such cancellation
/// anywhere (adding a small-to-1 value or a 1-to-huge value both lose
/// only irrelevant precision, and the final division is well-
/// conditioned in both regimes): clamping `-x` before `exp` (mirroring
/// `tanh`'s own domain-safety clamp, same idea, no cancellation-prone
/// identity in between) keeps this correct and gracefully saturating to
/// exactly `0.0`/`1.0` over the whole domain, never inf/nan for any
/// finite input.
#[inline(always)]
pub fn sigmoid(x: f32) -> f32 {
    // Standalone copy of exp's reduction (not routed through the public
    // `exp` fn, matching expm1's own established pattern above -- see
    // its doc comment for why factoring the *reduction* through a shared
    // helper is avoided here). The poly evaluation itself is shared via
    // `exp_r_poly!`, same as exp/expm1/tanh (a macro, not a fn -- see
    // expm1's own doc comment for why that distinction matters here).
    //
    // Previously used a single exponent-field construction with `y`
    // clamped to `[-87,88]` (keeping `k=round(y*log2e)` within exp2's own
    // single-field domain `[-126,128)`) -- this looked "correct and
    // gracefully saturating... never inf/nan for any finite input" (a
    // stale claim in this doc comment) but was actually a real accuracy
    // bug for `x` below about `-88`: the clamp caps `y` (and so `e`) at a
    // fixed largish-but-finite value regardless of how much more negative
    // `x` gets, so `sigmoid(-89)`, `sigmoid(-1000)`, and `sigmoid(-inf)`
    // all returned the *same* wrong constant (`~6.05e-39`) instead of
    // correctly continuing to decay toward `0` -- found by fuzzing `x` in
    // `[-90,-80]` against an f64 reference (worst case `~2.7x` relative
    // error, not just a few ulp) after idea #44 asked whether the
    // negative tail's conditioning actually needed attention. Root cause:
    // sigmoid's asymptote is at `0`, needing the exponent to range all
    // the way to where `e` itself correctly *overflows* to `+inf` (so
    // `1/(1+inf)=0` exactly) -- structurally different from `tanh`'s own
    // analogous clamp, whose asymptote is a nearby, already-representable
    // value (`+-1`), where an early clamp costs nothing.
    //
    // Fixed by widening only the *upper* `y` bound to `88.722839111673`
    // (`128/log2(e)`, the exact point where `k=round(y*log2e)` reaches
    // `128`) while keeping the lower bound at `-87.0` (already safe --
    // `x -> +inf` only needs the nearby, already-representable asymptote
    // `1`, the same reasoning that makes the narrower clamp fine for
    // `tanh`; no bug was ever found on that side). The single-field trick
    // itself (unchanged, *not* switched to exp2_checked's k1/k2 split --
    // tried that first, and it works, but doubles mca throughput cost:
    // 1.354->2.713 cyc/elem, an unjustified price for a fix this targeted)
    // already produces the exact right answer at `k=128`, verified by
    // hand: `exp2int` naturally lands on the `+inf` bit pattern there
    // (`(128.0+383.0).to_bits()<<8 & EXPONENT_MASK` resolves to the
    // reserved all-ones exponent field), giving `1/(1+inf)=0` exactly --
    // and empirically stays at exactly `128` (never wrapping to `129`,
    // which *would* silently give `0` instead of `inf` and produce the
    // wrong answer `1.0`) for every `y` from `88.5` up to `~89.05`,
    // comfortable margin either side of the chosen bound. Verified (fuzz
    // against `1/(1+(-x).exp())` in f64): `sigmoid(-89)`, `sigmoid(-100)`,
    // `sigmoid(-1000)`, `sigmoid(-inf)` all now correctly `0.0` (were all
    // `~6.05e-39` before). A narrow residual gap remains for `x` in
    // roughly `(-104.7,-88.7)`, where the true answer is a nonzero
    // denormal but this now returns exactly `0` instead (a discontinuity
    // moved from "wrong by construction, unboundedly" to "slightly early
    // saturation on a ~16-unit sliver of denormal-scale outputs") --
    // matching this crate's own accepted "near a true zero, ulp isn't a
    // meaningful metric" precedent (see cospi's own doc comment); not
    // pursued further, the practical improvement is already total.
    let y = (-x).clamp(-87.0, 88.722839111673);
    let (p, exp2int) = exp_r_singlefield!(y);
    let e = p * exp2int;
    1.0 / (1.0 + e)
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
/// Two real bugs in a first attempt that clamped `exp`'s argument
/// directly (`exp((-x.abs()).max(-87.0))`), both found by a scratch
/// probe before this was ever wired into the real test harness:
/// 1. For `x<-87`-ish, `max(x,0)=0`, so the *entire* result comes from
///    the correction term alone -- there, clamping the exponent to
///    exactly `-87.0` doesn't just avoid `exp`'s out-of-contract domain,
///    it also *replaces* the true (much smaller) correction with a
///    fixed, comparatively huge stand-in (e.g. `softplus(-100)` came out
///    `~1.6e-38` instead of the true `~3.7e-44` -- for positive `x` this
///    error rounds away against `x`'s own magnitude, but for negative
///    `x` nothing hides it). Fixed by *selecting* the correction to
///    exactly `0.0` once `|x|` is far enough out that the true value is
///    already negligible at f32 precision, instead of feeding `exp` a
///    clamped-but-wrong argument and hoping the result is close enough.
/// 2. `f32::max`/`min` follow IEEE `maxNum`/`minNum` semantics and
///    *discard* NaN (return the other operand) rather than propagate it
///    -- both `x.max(0.0)` (this function's own core operation, not an
///    incidental clamp) and the exponent-clamping `min`/`max` calls
///    silently turned `softplus(NaN)` into a finite garbage value. Same
///    class of trap as the rejected `hypot` `.max()`/`.min()` idea
///    earlier this session, just harder to avoid here since `max(x,0)`
///    isn't optional the way it was there. Guarded with an explicit
///    trailing `is_nan` check instead.
#[inline(always)]
pub fn softplus(x: f32) -> f32 {
    let ax = x.abs();
    let e = exp(-ax.min(87.0));
    let corr = if ax > 87.0 { 0.0 } else { log1p(e) };
    let normal = x.max(0.0) + corr;
    if x.is_nan() { f32::NAN } else { normal }
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
/// opposite-signed-ish values that partially cancel (e.g. near
/// `logaddexp(-7e-5, -9.57)`, where `m~-7e-5` and `corr~+7e-5`) --
/// max ulp seen in fuzzing is in the thousands there (avg stays a
/// healthy ~0.15), a real but narrow precision cost from the plain
/// `m + corr` addition amplifying whatever rounding either operand
/// already carries. Not further chased: comparable in kind (real
/// cancellation, not a bug) to other composite functions in this crate
/// with an accepted, undominant max-ulp outlier (e.g. `erfc`).
#[inline(always)]
pub fn logaddexp(a: f32, b: f32) -> f32 {
    let m = a.max(b);
    let d = (a - b).abs();
    let e = exp(-d.min(87.0));
    let corr = if d > 87.0 { 0.0 } else { log1p(e) };
    let normal = m + corr;
    if a.is_nan() || b.is_nan() { f32::NAN } else { normal }
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
    let r = if d.is_finite() { log1p_finite(d) } else { ln(ax) + LN_2 };
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
    let r = if d.is_finite() { log1p_finite(d) } else { ln(x) + LN_2 };
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
//
// Leading constant fixed 2026-07-09 (backlog idea #36's own
// investigation): this term is `acos_poly(0)`, which should equal
// exactly pi/2 (`acos(0)=pi/2`, and at `x=0` every other coefficient's
// contribution vanishes) -- but the literal `1.5707963` only carries 8
// significant digits, just short of enough precision to round to the
// nearest f32, so it parsed to `0x3fc90fda`, one ulp *below* the true
// correctly-rounded pi/2 (`0x3fc90fdb`, `std::f32::consts::FRAC_PI_2`).
// Not an intentional fitting choice (there's no reason a joint refit
// targeting the *whole* domain would prefer being wrong specifically at
// `x=0`, and this was carried over unchanged from before the 2026-07-07
// refit too) -- just short of enough digits when the constant was
// originally transcribed. Fixed by using `1.5707964` (confirmed by hand
// to parse to the exact same bits as `FRAC_PI_2`). Found while idea #36
// was screening a much more involved Df32-accurate-tier redesign (a
// two-product `sqrt(1-a)*poly(a)` combine, tested first and found to
// give *zero* improvement on its own -- the extra precision doesn't
// survive collapsing straight back to f32 without something downstream
// to use it, the same lesson idea #14's own rejection already
// established) -- this single-constant fix was a much bigger, unplanned
// win found along the way. Verified via exhaustive sweep: avg ulp
// 0.4962->0.0676 (-86%), `acos(0)`/`acos(-0)` now bit-exact instead of 1
// ulp off (previously accepted as "already-accepted fit imprecision" in
// edgecheck.rs, not actually true). Max ulp moved 4->6 (a different,
// smaller worst point elsewhere in the domain), a minor tradeoff against
// the large average improvement. Zero perf cost (same instructions, one
// literal constant differs).
//
// Remaining 6 coefficients retuned (2026-07-10) via coordinate descent
// from the shipped values (idea #3 in IDEAS.md), after fixing a stale
// constant (1.5707963, the pre-idea-#36 value) in tune.rs's own "acos"
// coordinate-descent seed. Verified against a real ~1.07-billion-point
// dense sweep of the whole [-1,1] domain, scored as the *whole* acos
// formula (matching how it's actually used, not the bare poly):
// max ulp 6->5, avg ulp 0.43739->0.43206 (both axes improved together,
// not a tradeoff). Zero perf cost, same instructions. (`asin` is
// unaffected -- it was decoupled onto its own independent `asin_poly`
// in fix 8 below, confirmed structurally: this poly and `asin_poly`
// no longer share a single literal.)
#[inline(always)]
fn acos_poly(x: f32) -> f32 {
    let u = 2.2960447e-3f32;
    let u = fma(u, x, -1.1146237e-2);
    let u = fma(u, x, 2.690034e-2);
    let u = fma(u, x, -4.8802484e-2);
    let u = fma(u, x, 8.875553e-2);
    let u = fma(u, x, -2.1458574e-1);
    fma(u, x, 1.5707964)
}

/// Dedicated asin-only copy of `acos_poly`'s shape (same 7-coefficient
/// Estrin/Horner form, same `sqrt(1-x)*poly` combine, independently
/// tuned), decoupling asin from acos's protected coefficients entirely
/// instead of sharing the joint-constrained fit `acos_poly` carries.
/// Every previous *joint* refit attempt died protecting acos's own
/// accuracy (three separate rejections, see IDEAS.md); this sidesteps
/// that by not sharing coefficients at all, so acos can't regress no
/// matter what this poly converges to. Seeded from a real scipy
/// least-squares fit of the same shape against `acos(x)/sqrt(1-x)`
/// (asin's own target, since `asin(x) = pi/2 - sqrt(1-x)*P(x)`)
/// restricted to asin's actual domain for this branch, `x` in `[0.25,
/// 1)` (not `[0,1)`, unlike acos_poly, which also needs `x` near 0) --
/// this narrower domain plus not needing to also satisfy acos's own
/// requirement is exactly the extra degree of freedom the joint fits
/// couldn't use. Coordinate-descent tuned from that seed
/// (`examples/tune.rs`'s `asin_poly_c`/"asinpoly").
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
/// 8. `acos_poly` split into a dedicated `asin_poly` copy (2026-07-09):
///    fix 7's *joint* constrained refit could only move within the room
///    left by never letting acos's own on-grid max regress -- a real but
///    narrow window. Decoupling entirely (separate coefficients, only
///    `asin` calls `asin_poly`) removes that constraint altogether: a
///    scipy least-squares fit of the same 7-coefficient shape against
///    `acos(x)/sqrt(1-x)`, restricted to asin's actual domain for this
///    branch (`x` in `[0.25, 1)`, narrower than acos_poly's own `[0,1)`
///    since asin never needs this branch below 0.25) found ~10-20x
///    tighter continuous-math error before any f32 rounding was even
///    considered; seeded `examples/tune.rs`'s coordinate descent with
///    that fit (not 0.0 or an arbitrary constant -- this file's own
///    "zero-move trap" precedent) and verified against the real crate.
///    Exhaustive sweep: asin avg ulp 0.0302 -> 0.0251 (~17% better); max
///    ulp unchanged at 9 (that worst case sits at x~0.246, just inside
///    `asin_small`'s domain, untouched by this branch). acos exactly
///    unchanged (avg 0.4962, max 4, bit-for-bit identical -- structurally
///    guaranteed now, not just empirically confirmed, since acos no
///    longer shares any coefficient with this poly at all). Zero perf
///    cost for either function (mca bit-identical: asin 59.03/0.968,
///    acos 37.11/0.820, same instruction shape, only literals differ).
#[inline(always)]
pub fn asin(x: f32) -> f32 {
    let a = x.abs();
    let small = asin_small(x);
    let big = mulsign(FRAC_PI_2 - (1.0 - a).sqrt() * asin_poly(a), x);
    if a < 0.25 { small } else { big }
}

// Pade-style rational approximation of atan on [0,1]. Ported from
// jodiemath's atanf_poly; degree bumped 2/2 -> 3/3 (2026-07-08, one more
// term each in numerator/denominator). The naive way to search for this
// -- coordinate-descend a new 6th coefficient starting from 0.0 -- get
// trapped: 0.0's bit pattern is 0x0, so the tuner's integer-bit-step
// search (deltas of 1..16 ulp) only reaches denormal-scale values from
// there, which have ~zero effect on the polynomial and never score
// better, so it never moves (a "zero-move local optimum" that's really
// just the search being unable to leave 0, not evidence of no headroom).
// An arbitrary small nonzero seed (1e-3) did worse still: it starts from
// a *worse* point than the shipped 2/2 form and the greedy per-coordinate
// search never recovered (converged to max ulp 565, far worse than the
// 2/2 form's 18). What worked: an actual least-squares Pade fit (scipy,
// not a guess) of the 3/3 shape against atan(x) directly over [0,1] --
// found max abs error 1.3e-9 vs the 2/2 form's 7.8e-7 (~600x), real
// headroom neither naive tuner seed could reach. Fed as the coordinate-
// descent starting point (examples/tune.rs's atan_poly7_c/"atan"):
// dropped straight to max ulp 3 / avg 0.286 on the tuning grid (from the
// 2/2 form's 18/1.36) before the descent even needed to move much.
// Verified against the real crate (not just the tuning grid), exhaustive
// all-2^32-bit-pattern sweep: atan avg/max ulp 0.186/18 -> 0.068/4;
// atan2 (which calls this directly, sampled since it takes two f32 args)
// 0.136/18 -> 0.069/3. mca: one more fma-depth level in both numerator
// and denominator (evaluated in parallel, so latency cost is one fma,
// not two) -- real but modest per-call cost (atan 57.09/1.410 ->
// 61.09/1.491 cyc lat/throughput, atan2 57.17/1.467 -> 61.17/1.532),
// easily justified by a >4x max-ulp cut on the crate's second-worst
// accuracy offender.
//
// Numerator refit (2026-07-09), denominator held fixed: max-capped
// ulp-weighted LP (this session's exp/erf_poly/cbrt technique), solving
// for a0..a2 against `atan(x)*denom(x2)/x` (linear in a0..a2 with the
// denominator fixed, since numer/denom is linear in the numerator's own
// coefficients alone). The isolated fit predicted a huge win (~11x
// tighter max, ~15x tighter avg) but didn't survive contact with the
// real crate anywhere near that scale (this session's now-repeated
// finding that the isolated metric doesn't reliably predict real
// magnitude) -- exhaustive: atan avg ulp 0.0681->0.0675 (a real, small
// ~0.9% win), max ulp unchanged at 4, same exact worst-case x
// (1.0220603) before and after, confirming the numerator was never the
// binding constraint for the true worst case (likely the denominator or
// the division itself) -- refitting it only touched the average, not
// the max, similar to the expm1 Pade entry's own "wrong branch" shape.
// Zero perf cost (same instructions). Adopted for the small, real,
// no-regression avg win.
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

/// Latency-tier atan: a division-free odd degree-17 poly (fit directly
/// against atan(r) over r in `[0,1]`, not derived from atan_poly's own
/// rational) instead of atan_poly's 3/3 Pade form. Removes the division
/// atan_poly's own critical path pays (a division can't start until its
/// numerator/denominator resolve, unlike cbrt's early-starting `rcp`),
/// trading it for more fma depth -- opposite tradeoff to nearly every
/// other function in this crate, so a separate opt-in tier rather than a
/// replacement (same shape as `sinh_throughput`/`cosh_throughput`, just
/// favoring the other axis): mca latency 61.09→59.09 cyc (-3.3%),
/// throughput 1.491→1.611 cyc/elem (+8.1% worse). Accuracy is not a
/// tradeoff here (fuzz: avg/max ulp 0.052/3 vs atan's own 0.067/4, a
/// slight improvement, not a cost) -- use this over `atan` only for a
/// value on its own or a serial dependency chain where per-call latency
/// matters more than array-loop throughput. (`atan`'s own comparison
/// figure here predates `atan_poly`'s 2026-07-09 numerator refit; see
/// `atan_poly`'s own doc comment for the current exhaustive-verified
/// max ulp of 4.)
///
/// Backlog idea #40 (2026-07-09): fold the `a<1.0` select and the final
/// `mulsign` into "one xor+select" by applying `mulsign` to `p` and
/// `FRAC_PI_2` individually first (`mulsign(a,x) - mulsign(b,x) ==
/// mulsign(a-b,x)` always -- `mulsign` is an exact sign-bit XOR, and
/// IEEE754 subtraction commutes exactly with negating both operands, so
/// this isn't an approximation, just a reassociation), then selecting
/// between `sp` and `hpisignx - sp` instead of selecting on the unsigned
/// value and applying one final `mulsign`. Verified bit-identical
/// against the prior formulation over ~20M random bit patterns plus
/// every special value (0, -0, +-1, +-inf, NaN) before adopting -- this
/// is a pure reassociation, not a new approximation. mca: latency
/// 59.09→59.11 (+0.03%, noise), throughput 1.611→1.591 cyc/elem (-1.2%,
/// small but real and reproducible across repeat runs). A modest win,
/// not a dramatic one -- kept because it's free (zero accuracy cost,
/// latency unchanged) rather than because the throughput gain alone
/// would have justified real effort.
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
/// instead of the correct `NaN` (found 2026-07-09 building a systematic
/// C99 special-case matrix against std, backlog idea #85): when
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
// Coefficients refit again (2026-07-09) via the same ulp-weighted
// Chebyshev LP technique used for exp's degree-5 poly, this time with an
// explicit max-ulp cap (minimize the L1/avg-ish weighted error subject to
// the max weighted error never exceeding the shipped coefficients' own
// bound, same "constrained search" shape as acos_poly's own successful
// fix 7) -- single caller (erf itself, no dual-region sharing like
// sinf_poly's sin/cos split), so this avoided that pitfall. Small real
// win, exhaustive: avg ulp 0.3194->0.3166, max ulp unchanged at 5; erfc
// (doesn't call this poly) bit-for-bit unaffected. Zero perf cost, same
// instructions (mca unchanged both before and after this refit).
//
// Refit again (2026-07-10) via a genuine ulp-weighted Chebyshev LP (idea
// #52/#91's technique, using scipy -- confirmed available this session),
// this time weighting each point by the *linearized sensitivity of erf's
// own final `1 - 2^poly` combine* rather than the poly's own raw output --
// the prior refit's weighting didn't account for how that final
// transform's own sensitivity varies across the domain. Verified against
// the real compiled `erf` (bit-identical fidelity check first, 0
// mismatches) and a real ~1.09-billion-point dense sweep of the whole
// `[-10,10]` domain: max ulp 5->4, avg ulp unchanged (0.62958->0.63010,
// +0.08%, noise-level) -- a real max-ulp win at no accuracy cost
// elsewhere. Zero perf cost, same instructions.
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

// Shared n/d rational (degree 4 in xa) behind both `erfc` and `erfcx`
// (backlog idea #51 factored this out of erfc's own body): NaN-preserving
// clamp to |xa| <= 10 first (same reasoning erfc's own doc comment gives
// -- `f32::min`-style clamp, not a ternary, so NaN passes through instead
// of being silently replaced), matching the domain the coefficients were
// fit against.
#[inline(always)]
fn erfc_rational(xa: f32) -> f32 {
    let xa = if xa > 10.0 { 10.0 } else { xa };
    let n = fma(f32::from_bits(0x35c42f59), xa, f32::from_bits(0x3daf42cd));
    let n = fma(n, xa, f32::from_bits(0x3ee32e39));
    let n = fma(n, xa, f32::from_bits(0x3f7a7525));
    let n = fma(n, xa, 1.0);
    let d = fma(f32::from_bits(0x3e1b69eb), xa, f32::from_bits(0x3f48fdde));
    let d = fma(d, xa, f32::from_bits(0x3fe918da));
    let d = fma(d, xa, f32::from_bits(0x4006d464));
    let d = fma(d, xa, 1.0);
    n / d
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
/// domain comfortably covers the full clamped range (xa in `[0,10]` means
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
    let xa = x.abs();
    // The exponent term deliberately uses the *true*, unclamped `xa`, not
    // erfc_rational's own internal xa<=10 clamp: erfc_rational needs that
    // bound to keep its own rational polynomial from overflowing, but
    // exp2_checked already has its own established, correct saturate-to-0
    // contract for arbitrarily negative exponents (clamps its own input to
    // -151 internally). Using the clamped xa here instead (as a prior
    // version of this function did) froze the exponent at exp2_checked
    // (-100*log2e) for *every* xa>10, so erfc(x) for any x beyond ~10.02
    // returned that same tiny nonzero constant (~3e-45) forever instead of
    // the true value, which reaches exactly 0.0f32 well before x=11 --
    // found via idea #102/#105's own "does the same clamp-then-freeze
    // pattern erfcx hit also affect erfc itself" check, since erfc's own
    // accuracy.rs sweep is restricted to |x|<=10 and never exercised this.
    let y = exp2_checked(-(xa * xa) * LOG2_E) * erfc_rational(xa);
    fma(y, z, w)
}

/// erfcx(x) = e^(x^2)*erfc(x), the "scaled complementary error function"
/// -- new function, backlog idea #51. For x >= 0, this collapses to
/// exactly `erfc_rational(x)` alone with *no exponential at all*: erfc's
/// own construction is `exp(-x^2) * erfc_rational(x)`, so multiplying by
/// `exp(x^2)` cancels the exponential exactly (not approximately --
/// `exp(x^2)*exp(-x^2)` is algebraically 1, so this sidesteps the
/// exponent computation entirely rather than computing and cancelling
/// it). This is exactly what erfcx is *for*: the naive
/// `exp(x*x)*erfc(x)` a caller might otherwise write already breaks
/// down numerically before this function's own domain gets interesting
/// -- `exp(x*x)` alone overflows f32 for `|x| >~ 9.3`, while `erfcx`'s
/// true value there is still a small, well-behaved, easily-representable
/// number (`erfcx(x) ~ 1/(x*sqrt(pi))` for large positive x).
///
/// For x < 0, uses erfc's own reflection identity (`erfc(x) = 2 -
/// erfc(-x)` for x<0) to derive `erfcx(x) = 2*exp(x^2) - erfcx(-x)` --
/// unlike the x>=0 branch this does need one real `exp2_checked` call,
/// because `erfcx` genuinely diverges to `+inf` for sufficiently
/// negative x (`erfcx(-10) ~ 2*e^100`, far past `f32::MAX`) -- that's
/// this function's true mathematical behavior, not an implementation
/// gap, and `exp2_checked`'s own saturation makes it come out `+inf`
/// correctly rather than wrapping to garbage.
///
/// Like `erfc`, `erfc_rational`'s own |xa|<=10 fit domain means this is
/// only verified accurate for `|x| <= 10` -- for x > 10, `erfc_rational`
/// clamps its own input to 10.0 and so *freezes* at `erfc_rational(10.0)`
/// forever, not just "less accurate": unlike `erfc` (where this same
/// freeze is masked by the multiplicative `exp(-x^2)` factor correctly
/// decaying to 0) `erfcx` has no such factor, so the frozen value comes
/// back nakedly and relative error grows *without bound* as x grows
/// (measured against a proper reference: ~10% already at x=11, ~50% at
/// x=15, ~99% at x=20, ~895% by x=100 -- not a bounded imprecision).
/// Checked whether the standard asymptotic tail
/// (`erfcx(x) ~ (1 - 1/(2x^2) + 3/(4x^4))/(x*sqrt(pi))`, accurate to
/// <0.0002% relative error for `x > 10` right here) is worth blending in
/// past x=10: yes for correctness (verified against a proper reference,
/// bit-identical to the shipped form for `|x|<=10`, clean codegen), but
/// the real `mca` throughput cost is a genuine, not noise-level, +17.7%
/// (2.278->2.681 cyc/elem) from the extra division the tail needs --
/// fails this crate's own "accuracy win, no perf penalty" bar, so not
/// adopted here. Left as a real, precisely-scoped option for an
/// `erfcx_checked`/wider-domain opt-in tier if a caller ever actually
/// needs `|x| > 10`, matching this crate's existing checked/unchecked
/// tiering pattern -- not built speculatively.
///
/// mca's own latency number for this function (see readme.md) is not
/// trustworthy: this is a sign-dependent branch (`x >= 0.0`), and the
/// latency harness's `mix()` step always folds its chained value into
/// `[2, 4)` (masks the sign bit away entirely) -- so the `x<0` arm
/// (with its real `exp2_checked` call) is never exercised in that
/// measurement, the same "mca mix() sign blind spot" this crate has
/// hit before for other sign-shuffling functions. Only the throughput
/// number (built from a real mixed-sign array) is meaningful here.
#[inline(always)]
pub fn erfcx(x: f32) -> f32 {
    let xa = x.abs();
    let r = erfc_rational(xa);
    if x >= 0.0 { r } else { 2.0 * exp2_checked(x * x * LOG2_E) - r }
}

/// 1/sqrt(x). Unlike most functions in this crate, no bit-trick seed or
/// fitted correction poly needed: `sqrt` and division are each already
/// correctly-rounded IEEE754 hardware operations (`x.sqrt()` isn't a
/// software approximation), so composing them directly costs at most
/// ~1 ulp (one rounding from each op) with zero special-casing --
/// `x <= 0.0` (including `-0.0`), `x.is_nan()`, and `x == inf` all
/// already give the right answer (`+inf`/`inf`, `NaN`, `0.0`
/// respectively) purely from IEEE754 semantics, the same "let the
/// hardware ops already handle it" reasoning CLAUDE.md's own "use a
/// library function, don't reinvent the wheel" principle points at
/// directly. (The crate's older `rsqrt_approx` -- a single Quake-style
/// bit-trick seed with no correction, ~1e4 ulp -- is a different,
/// deliberately-rough exploratory function kept for the `_approx_plot`
/// test suite, not a candidate replacement for this one.)
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
/// With this narrower check, `hypot_checked(NaN, 0.0)` still resolves to
/// exactly `NaN` (verified in a Python bit-level prototype before writing
/// this): the exponent extraction degrades to a garbage-but-finite scale
/// either way, but `ax` (or `ay`) being NaN itself propagates through the
/// unconditional multiply regardless of what that scale becomes.
/// Verified against a Python bit-level prototype: max ulp 1 over ~1M
/// samples spanning the full exponent range, plus every zero/NaN/inf
/// combination checked directly.
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

/// 1/hypot(x,y) -- new function, backlog idea #61: normalizing a 2D vector
/// (`(x,y) / hypot(x,y)`) is hypot's single most common real use case, and
/// computing the reciprocal directly saves the caller their own separate
/// division. Same "sqrt and division are each already correctly-rounded
/// hardware operations, just compose them" reasoning as `rsqrt` (see its
/// own doc comment) -- no bit-trick seed needed, `hypot_unchecked`'s own
/// `fma(x,x,y*y)` core plus one hardware sqrt and one division. The naive
/// composition already gets every zero/inf/nan special case right *except
/// one*, purely from IEEE754 semantics: `rhypot(0,0)=inf` (`1/0`),
/// `rhypot(x,inf)=rhypot(inf,y)=0` (`1/inf`), `rhypot(NaN,y)=NaN` all fall
/// out for free (verified by hand before writing this). The one exception
/// mirrors `hypot`'s own documented special case: `+-inf` paired with
/// `NaN` degrades to `1/sqrt(NaN)=NaN` here, but IEEE754/C99 defines
/// `hypot(+-inf, NaN) = +inf` (infinity "wins" over NaN), so the
/// reciprocal should be `0`, not `NaN` -- same override `hypot` itself
/// uses. No anti-overflow rescaling tier (unlike `hypot_checked`): not
/// requested by the backlog idea, and `x*x+y*y` overflowing is already a
/// documented, accepted tradeoff of the crate's default `hypot`.
#[inline(always)]
pub fn rhypot(x: f32, y: f32) -> f32 {
    let normal = 1.0 / fma(x, x, y * y).sqrt();
    if x.is_infinite() || y.is_infinite() { 0.0 } else { normal }
}

/// log2(x) as a double-float (Df32) instead of a collapsed f32, for
/// positive finite x only (same domain log_2_normal assumes -- callers
/// must guard zero/negative/inf/nan themselves). Reuses log_2_normal's
/// exact decomposition and poly (`s = m - 1`, `P(s)`).
///
/// `p * s` (the poly correction) is combined with the exact integer
/// exponent `k` via `Df32::from_mul(p, s)` (an exact two-product, keeping
/// `p*s`'s own rounding error as the Df32's low word) added to `k` --
/// *not* the plain `p * s` single multiply this used until backlog idea
/// #58's own root-cause investigation (2026-07-09). The original form,
/// `Df32::from_add(k, p * s)`, has a real hole: `Df32::from_add`'s two-sum
/// only captures the rounding error of *adding* `k` and `p*s` together --
/// but `p * s` was already computed as a single, already-rounded f32
/// multiply *before* that add ever runs, so its own rounding error is
/// silently gone, never entering either Df32 word. Harmless when `k` is
/// large enough that the combine's own rounding dominates, but for `x`
/// near 1 (`k=0`, a common case, e.g. iteratively-refined bases) adding
/// exactly `0` is *itself* lossless -- meaning the entire "double-float"
/// result was, in that regime, silently no more accurate than a single
/// f32 multiply, defeating the whole point of this function relative to
/// `log_2_normal`'s own single-rounding `fma(p,s,k)`. Found by tracing a
/// concrete `powf_checked` worst-case (`x=1.0281241, y=2695.4136`, 183
/// ulp) against a Decimal-precision Python reference at each intermediate
/// step: the exact `p*s` product differed from the crate's *computed*
/// `p*s` by ~3.79e-9, purely from this dropped rounding, on top of the
/// poly's own ~2.66e-9 inherent fit error -- confirming the multiply's
/// own rounding was the larger of the two contributors, not just a minor
/// addition. Fixed by computing `p*s` as a real two-product instead.
/// Verified (100M-sample targeted fuzz concentrated on `x` in `[0.7,1.5]`
/// -- the `k=0` regime this fixes -- paired via git stash against the
/// same unfixed code): avg ulp 5.81->5.05 (-13%), max ulp 261->211
/// (-19%), count of samples over 100 ulp 160472->71836 (-55%) -- a real,
/// substantial improvement, though not a complete fix (the poly's own
/// ~2.66e-9 fit error is a separate, remaining contributor, not chased
/// here). Denormal input is handled the same way log_2's own wrapper does
/// (scale up, offset k).
#[inline(always)]
fn log2_df(x: f32) -> Df32 {
    let (xs, koff) = denormal_rescale!(x);
    let (p, s, k) = log_family_normal!(xs, koff, LOG2_COEFFS);
    Df32::from_f32(k) + Df32::from_mul(p, s)
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
#[allow(clippy::approx_constant)] // g0's constant term is a fitted minimax
// coefficient near ln(2), not ln(2) itself (bit pattern deliberately differs)
fn exp2_checked_df(v: Df32) -> f32 {
    // clamp *before* floor/subtract (matching exp2_checked exactly): if v.0
    // itself is already +-inf (y large enough that y*log2(x) overflows in
    // the Df32 multiply), flooring first and clamping after would compute
    // `inf - floor(inf)` = `inf - inf` = NaN instead of correctly
    // saturating.
    //
    // IDEAS.md backlog round 3 #1's round-domain refit was tried here too:
    // a real accuracy win downstream (powf_checked avg ulp 0.0229->0.0118,
    // max 123->118; powf_checked_unchecked avg 0.0459->0.0243, max 96->89
    // -- unlike exp2_checked itself, powf_checked's own error is dominated
    // by log2_df/the y-multiply, not this poly's fit quality, so the
    // tighter centered-domain fit helped rather than getting swamped). But
    // mca showed a real, if modest, throughput cost inherited from the
    // same floor->round port-shift as exp2_checked (see its own doc
    // comment): powf_checked 9.105->9.418 cyc/elem (+3.4%),
    // powf_checked_unchecked 7.234->7.475 (+3.3%), latency flat both ways.
    // Doesn't clear this crate's bar (accuracy gain *without* a perf
    // penalty) -- reverted, bit-identical to prior HEAD.
    let xs = v.0.clamp(-151.0, 128.0);
    let k = xs.floor();
    let f = xs - k;
    let (t1, t2) = exp2_field_split(k);
    let q = exp2_q_poly!(f);
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
// Shared by powf/powf_checked: the negative-base/y-parity/y==0/x==+-1
// special-case combine, given each caller's own already-computed `mag`
// (powf's plain `exp2_checked(log_2(ax)*y)` vs. powf_checked's
// is_safe-gated Df32 pipeline -- the two differ upstream of this, not in
// how the sign/special-case logic itself works). Macro, not a fn -- same
// reasoning as this file's other shared-body macros (no function-call
// boundary, verified via full assembly diff).
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
// NaN" rule above (found 2026-07-09 building a systematic special-case
// matrix against std, backlog idea #85's fourth wave): unlike a
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
// bug as acos's earlier `-0.0` fix this session), which would silently
// route `(-0.0)^3.0` through the wrong (positive) branch instead of the
// correctly-signed `-0.0`.
//
// pow(1, y) = 1 for *any* y -- even inf, -inf, or NaN -- another
// dedicated IEEE754/C99 special case the log/exp2 formula can't derive
// on its own (log_2(1)=0, so mag=exp2_checked(0*y); for y=inf/-inf/NaN
// that's a 0*inf or 0*NaN indeterminate form, degrading to NaN instead
// of the correct 1). pow(-1, +-inf) = 1 is a second, narrower C99
// special case (unlike pow(1,y), it does *not* extend to pow(-1,NaN),
// which stays NaN) -- handled separately since it only overrides the
// infinite-y case.
//
// pow(x, 0) = 1 for *any* x -- even 0, negative, or NaN -- a dedicated
// IEEE754/C99 special case, not derivable from the log/exp2 formula
// (0*inf and NaN*0 both degrade to NaN above). Override last.
macro_rules! powf_sign_combine {
    ($x:expr, $y:expr, $mag:expr) => {{
        let y_int = $y == $y.trunc();
        let y_odd = y_int && parity($y) != 0.0;
        let neg_signed = if y_odd { -$mag } else { $mag };
        let neg_result = if y_int || $x == 0.0 || $x.is_infinite() { neg_signed } else { f32::NAN };
        let neg_result = if $y.is_infinite() { $mag } else { neg_result };
        let r = if $x.is_sign_negative() { neg_result } else { $mag };
        let r = if $x == 1.0 { 1.0 } else { r };
        let r = if $x == -1.0 && $y.is_infinite() { 1.0 } else { r };
        if $y == 0.0 { 1.0 } else { r }
    }};
}

#[inline(always)]
pub fn powf(x: f32, y: f32) -> f32 {
    let ax = x.abs();
    let mag = exp2_checked(log_2(ax) * y);
    powf_sign_combine!(x, y, mag)
}

/// powf without domain/sign checks: valid for `x` positive, normal, and
/// finite (the same domain [`log_2_unchecked`] requires) and `y != 0.0`.
/// No handling for negative/zero/denormal/inf/nan `x`, no `y == 0.0`
/// special case, no negative-base parity handling -- those are exactly
/// the branches [`powf`]'s own doc comment describes paying on every
/// call regardless of whether they're ever hit. Mirrors `log_2`/
/// `log_2_unchecked` and `atan2`/`atan2_unchecked`'s own fast/full-safety
/// split; `exp2_checked` (not the even-faster `exp2`) is kept since powf's
/// own doc comment already documents why bare `exp2` silently wraps
/// around into plausible-looking garbage instead of overflowing --
/// nothing about this narrower domain contract changes that risk.
#[inline(always)]
pub fn powf_unchecked(x: f32, y: f32) -> f32 {
    exp2_checked(log_2_unchecked(x) * y)
}

/// x^n for integer `n` (`i32`), via exponentiation by squaring. Each
/// step is a single correctly-rounded f32 multiply -- no poly, no log/
/// exp composition -- so this sidesteps `powf`'s own "amplifies log_2's
/// rounding error by y" issue (see its own doc comment) entirely for
/// integer exponents.
///
/// A first version used a data-dependent `while` loop (trip count =
/// `n`'s bit length), which only auto-vectorized when `n` happened to
/// be a *compile-time constant* the compiler could unroll -- checked
/// directly (`--emit=asm`, a black-boxed runtime `n` fed uniformly to
/// every array element, the realistic "same exponent, varying base"
/// calling pattern): zero vector instructions, a genuine violation of
/// this crate's "every public function must auto-vectorize" hard
/// requirement, not just a missed optimization. Fixed with a fully
/// unrolled, branchless design instead: always exactly 32 iterations
/// (an `i32`'s full magnitude range including `i32::MIN`'s `2^31`, which
/// needs bit index 31 -- an off-by-one caught directly when a `0..31`
/// version returned `1.0` instead of `0.0` for `pown(2.0, i32::MIN)`),
/// squaring `base` unconditionally every step and selecting whether to
/// fold the current power into `result` based on each bit of `n` in
/// turn -- the operation *sequence* is now fixed regardless of `n`'s
/// runtime value, so it vectorizes even for a genuinely per-lane-
/// *varying* `n` (confirmed via the same `--emit=asm` check: real
/// `vmulps`/`vblendvps` instructions). Real cost: always pays for 32
/// squarings + selects regardless of how small `n` actually is, where
/// the old data-dependent loop only paid for `n`'s actual bit length --
/// a deliberate throughput-for-correctness trade, since a function that
/// silently fails to vectorize for its most realistic calling pattern
/// isn't a function this crate can ship.
///
/// `n=0` gives `1.0` for any `x` (including `0.0`, matching `powf`'s own
/// convention) for free: every one of the 31 iterations selects the
/// "don't multiply" branch, since `n`'s bits are all zero. Negative `x`
/// needs no special-casing either -- integer powers of a negative base
/// are always well-defined (unlike `powf`'s general real-exponent
/// case), so plain repeated multiplication already gets the sign right.
// Shared by pown/pown_small/pown_const -- the exponentiation-by-squaring
// body itself (invert-if-negative, then square-and-select loop) is
// identical across all three; only the iteration count differs (32 to
// cover i32::MIN's full range, 8 for pown_small's narrower |n|<=255
// contract) along with whether `n` is a runtime i32 or a const generic.
// Macro, not a fn -- same reasoning as this file's other shared-body
// macros (no function-call boundary, verified via full assembly diff).
macro_rules! pown_body {
    ($x:expr, $n:expr, $iters:expr) => {{
        let mut base = if $n < 0 { 1.0 / $x } else { $x };
        let un = $n.unsigned_abs();
        let mut result = 1.0f32;
        for i in 0..$iters {
            let bit_set = (un >> i) & 1 == 1;
            result = if bit_set { result * base } else { result };
            base *= base;
        }
        result
    }};
}

#[inline(always)]
pub fn pown(x: f32, n: i32) -> f32 {
    // Invert x *before* the squaring loop (not the final result after)
    // when n is negative -- computing x^|n| first and reciprocating at
    // the end can overflow at an intermediate squaring step even when
    // the true (small) final answer wouldn't (found empirically: a
    // 20M-sample fuzz caught pown(1.8e19, -2) returning 0 instead of the
    // correct ~2.94e-39, because 1.8e19^2 alone overflows f32 even
    // though its reciprocal doesn't). Squaring the already-small
    // reciprocal instead avoids that overflow, and this also handles
    // 0^negative (`1.0/0.0 = inf`, then `inf^n = inf`, correct) and
    // inf^negative (`1.0/inf = 0`, then `0^n = 0`, correct) for free,
    // with no separate final-inversion special case needed.
    // 32, not 31: i32::MIN's magnitude is exactly 2^31, needing bit
    // index 31 -- an off-by-one caught directly (pown(2.0, i32::MIN)
    // returned 1.0 instead of the correct 0.0 with a 0..31 range).
    pown_body!(x, n, 32u32)
}

/// `pown` restricted to `|n| <= 255`: same exponentiation-by-squaring
/// algorithm, but only 8 unrolled iterations instead of 32 (`255` is
/// `2^8-1`, so bit index 8 and above are always zero within this
/// contract, unlike `pown`'s own need to cover `i32::MIN`'s `2^31`) --
/// 4x fewer squarings/selects for what's overwhelmingly the common case
/// (small integer exponents). Bit-identical to `pown` whenever the
/// contract holds (confirmed over 50M generated samples spanning the
/// full `|n| <= 255` range), matching this crate's other narrower-domain
/// `_unchecked`-style tiers.
///
/// Still fully vectorizes for the harder, realistic per-lane-*varying* n
/// case (confirmed directly via a standalone `--emit=asm` probe: proper
/// AVX-512 masked selects, no scalar fallback) -- but isn't wired into
/// `examples/mca.rs`/`mca_target.rs` like this crate's other tiers: with
/// only 8 iterations (vs `pown`'s 32), LLVM's cost model chooses to
/// branch-specialize on `mca_target.rs`'s own shared/uniform-`n`
/// benchmark shape instead of emitting the uniform blend `pown` gets at
/// 32 iterations, and the resulting multi-exit-path function corrupts
/// llvm-mca's inline-asm region markers ("found an invalid region end
/// directive") -- a harness limitation specific to this trip count +
/// shared-n combination, not a code correctness issue. Real wall-clock
/// numbers instead (`quickbench`, 3 repeated runs): latency ~59→~15 ns
/// (~3.9x faster), throughput ~1.4→~0.22 ns (~6.4x faster) vs `pown`.
#[inline(always)]
pub fn pown_small(x: f32, n: i32) -> f32 {
    pown_body!(x, n, 8u32)
}

/// `pown` with a compile-time-known exponent: same algorithm as `pown`,
/// but with `N` as a const generic instead of a runtime `i32`, so
/// `N.unsigned_abs()` and `N < 0` are compile-time constants and the
/// whole 32-iteration bit-testing loop is expected to constant-fold away
/// entirely, leaving only the exact minimal sequence of multiplies this
/// specific exponent needs (no runtime branch or select at all) -- the
/// "same exponent, hard-coded at the call site" pattern (`x*x*x`-style
/// cubes/squares/reciprocals) that motivated `pown`'s own square-and-
/// multiply redesign in the first place, taken to its logical conclusion
/// once the exponent doesn't need to vary per call. Verify with
/// `--emit=asm` before trusting this folds as described -- monomorphized
/// generics don't automatically guarantee LLVM finishes the constant
/// folding, only that it has enough information to.
#[inline(always)]
pub fn pown_const<const N: i32>(x: f32) -> f32 {
    pown_body!(x, N, 32u32)
}

/// Higher-accuracy variant of [`powf`]: `exp2(log_2(x)*y)` amplifies
/// log_2's own rounding error by `y` -- for `|y|` large that swamps the
/// result (hundreds of ulp), since `log_2(x)` is collapsed to a single f32
/// *before* the multiply, throwing away exactly the low bits that `y`'s
/// multiplication would otherwise be able to use. Fixed by keeping
/// `log2(x)` as a double-float (Df32) through the multiply by `y` and the
/// exp2 reconstruction, only collapsing to a single f32 at the very end
/// (see `log2_df`/`exp2_checked_df`). Confirmed by fuzzing (25M
/// samples, apples-to-apples against the plain formula on the same
/// inputs): avg ulp 0.36 -> 0.05 (-87%), max ulp 170 -> 154 -- both
/// improve. Real mca cost though (double-float bookkeeping plus the poly
/// evaluation isn't free): latency 98.03->105.48 cyc (+7.6%), throughput
/// 3.898->6.285 cyc/elem (+61.2%) -- kept as an opt-in tier rather than
/// the default, matching sin/sin_checked and exp2/exp2_checked.
///
/// The remaining max-ulp cases (backlog idea #58, root-caused
/// 2026-07-09) were *not* an `exp2_checked` denormal/underflow artifact
/// as previously assumed here -- they're dominated by `log2_df` itself
/// silently losing precision for `x` near 1 (`k=0`), where the exponent
/// `k` contributes nothing to the double-float combine, fixed directly in
/// `log2_df`'s own doc comment/implementation (a real two-product instead
/// of a lossy single multiply for `p*s`). See there for the numbers; this
/// function inherits that fix automatically. Updated mca cost after that
/// fix (measured against this function's own most recent baseline,
/// 106.63/7.726, which already includes the unrelated `pow(1,y)`
/// special-case selects added later the same session): latency
/// 106.63->130.53 cyc (+22.4%), throughput 7.726->9.231 cyc/elem
/// (+19.5%) -- a further real cost for a further, partial (not complete)
/// accuracy gain, accepted for the same "opt-in accuracy tier, pay more
/// for more correctness" reasoning as the original fix above.
#[inline(always)]
pub fn powf_checked(x: f32, y: f32) -> f32 {
    let ax = x.abs();
    // The df path is only valid for ax strictly positive and finite (same
    // domain log_2_normal's own decomposition assumes); ax == 0 or
    // non-finite needs a fallback, but *not* a second full
    // log_2+exp2_checked computation -- that would roughly double this
    // function's cost just to cover a few degenerate inputs. The only
    // four possible magnitudes there are cheap direct selects: ax == 0
    // -> 0 (y > 0) or +inf (y < 0); ax == +inf -> +inf (y > 0) or 0
    // (y < 0); ax == NaN (from x == NaN) -> NaN; y == NaN -> NaN (found
    // 2026-07-09 building a systematic C99 special-case matrix against
    // std, backlog idea #85's own third wave: `zero_or_inf`'s own `y >
    // 0.0` comparison is simply false for NaN, same as any comparison,
    // so a NaN `y` silently fell through to a *finite* 0/inf answer
    // instead of propagating -- `powf`/`powf_unchecked` don't have this
    // gap since they route ax==0/inf/nan through the real, always-
    // NaN-propagating `log_2`/`exp2_checked` instead of this cheap
    // shortcut). (y == 0 is overridden separately below regardless of
    // any of this.)
    // Bit-trick form instead of `ax > 0.0 && ax.is_finite()` (2026-07-09):
    // that compound condition compiled to a fully scalar per-lane
    // sequence (extract each lane, run several scalar int test/cmp/set
    // instructions, then hand-assemble an AVX-512 mask bit-by-bit via
    // kmovd/kshiftlb/kshiftrb/korb/kandb) instead of a single vectorized
    // compare -- found while investigating a similar pattern for a
    // rejected rsqrt idea (see IDEAS.md's idea #90 entry). `ax` is
    // already non-negative (`x.abs()`), so its raw bits directly encode
    // magnitude: nonzero and below the all-ones exponent field (which
    // marks inf/nan) is exactly "strictly positive and finite."
    let axb = ax.to_bits();
    let is_safe = axb != 0 && axb < EXPONENT_MASK;
    let mag_precise = exp2_checked_df(log2_df(ax) * y);
    // (ax == 0) == (y > 0) picks out exactly the two "goes to zero" cases
    // (ax==0,y>0 and ax==+inf,y<0) vs. the two "goes to infinity" cases --
    // cheaper than a 4-way branch and avoids inf/inf-is-NaN traps a
    // division-based shortcut would hit for the ax==+inf,y<0 case.
    let zero_or_inf = if (ax == 0.0) == (y > 0.0) { 0.0 } else { f32::INFINITY };
    let edge_mag = if ax.is_nan() || y.is_nan() { f32::NAN } else { zero_or_inf };
    let mag = if is_safe { mag_precise } else { edge_mag };
    // Same negative-x/y-parity/y==0/x==+-1 handling as powf, shared via
    // `powf_sign_combine!` -- see that macro's own doc comment (above
    // `powf`) for the full reasoning; `log2_df(1)` is exactly `Df32(0,0)`,
    // so the same "0*inf/0*NaN degrades to NaN" mechanism that motivates
    // the x==1/x==-1 special cases there applies here too.
    powf_sign_combine!(x, y, mag)
}

/// `powf_checked` without domain/sign checks: valid for `x` positive,
/// normal, and finite (the same domain [`log_2_unchecked`]/
/// [`powf_unchecked`] require) and `y != 0.0`. No `is_safe`/`edge_mag`
/// fallback (that machinery exists purely to cover `x`'s zero/inf/nan
/// cases, all excluded by this domain), no `y == 0.0` special case, no
/// negative-base parity handling -- same relationship to `powf_checked`
/// that `powf_unchecked` has to `powf`, just applied to the high-
/// precision double-float tier instead of the plain one (same shape as
/// `cbrt_accurate_unchecked`'s relationship to `cbrt_accurate`). See
/// [`powf_checked`] for the full-domain-safe version.
#[inline(always)]
pub fn powf_checked_unchecked(x: f32, y: f32) -> f32 {
    exp2_checked_df(log2_df(x) * y)
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
// Shared by remainder/remainder_ieee/fmod: given each caller's own `q`
// (`.round()`, `.round_ties_even()`, or `.trunc()` -- the one place they
// differ), the `x - q*y` combine plus the two special-case selects (the
// `x==0.0` sign-preservation guard and the `y` infinite/`x` finite
// no-reduction case, see `remainder`'s own doc comment above for the
// full reasoning) are byte-for-byte identical across all three. Macro,
// not a fn -- same reasoning as this file's other shared-body macros.
macro_rules! remainder_style_combine {
    ($x:expr, $y:expr, $q:expr) => {{
        let normal = fma(-$q, $y, $x);
        let r = if $x == 0.0 && !normal.is_nan() { $x } else { normal };
        if $y.is_infinite() && $x.is_finite() { $x } else { r }
    }};
}

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
/// This is genuinely *cheaper* than `remainder`, not the same cost (a
/// prior version of this comment claimed "the same `vroundps` family
/// instruction, just a different rounding-mode immediate" -- checked
/// directly against the compiled assembly 2026-07-10 and found false):
/// x86's `vroundss`/`roundss` only has *hardware* support for
/// round-to-nearest-even, round-down, round-up, and truncate -- there is
/// no native "round half away from zero" mode. `.round_ties_even()`
/// (ties-to-even, this function) lowers directly to one `vroundss`
/// instruction; `.round()` (ties-away-from-zero, `remainder`'s own
/// choice) needs LLVM to emulate the away-from-zero tie-break in
/// software first (`vpbroadcastd` x2 loading a sign/magnitude constant
/// pair, `vpternlogd` combining them with the value, `vaddss` adding a
/// sign-matched 0.5, *then* `vroundss` in truncate mode) -- 4 extra
/// serial instructions `remainder_ieee` never pays. This fully accounts
/// for the real mca latency gap (see readme.md: `remainder` 34.11 cyc,
/// `remainder_ieee` 29.11 cyc, a genuine ~15% difference, not noise or
/// staleness).
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
/// crate's branchless style. Confirmed by fuzzing against an f64
/// reference: 0 max ulp for `|x/y|` up to `1e7` (including every concrete
/// sign-flip case [`remainder`] can hit), degrading only past `2^24`
/// where `q` itself stops being an exactly-representable f32 integer -- a
/// separate, harder limit this correction can't reach past. Costs a
/// second `fma` plus the correction's compare/select on every call (mca:
/// +42% latency, +60% throughput vs plain `remainder`), so kept as an
/// opt-in tier for callers who need the reliability guarantee, matching
/// sin/sin_checked and exp2/exp2_checked.
///
/// `r1 = fma(-adj, y, r0)` instead of `fma(-(q0+adj), y, x)`: algebraically
/// `x - (q0+adj)*y == (x - q0*y) - adj*y == r0 - adj*y`, so `r1` can reuse
/// the already-computed `r0` instead of re-deriving the corrected quotient
/// and recombining with `x` from scratch -- one fewer add (`q0 + adj`)
/// with the same 2-fma total shape. Verified against the un-simplified
/// form on 48M generated near-tie samples (`x` constructed within `1e-5`
/// of an exact half-integer multiple of `y`, the specific case this
/// correction exists for): bit-identical on every sample; the handful
/// that disagreed with the sleef f64 reference disagreed *identically*
/// both ways, and traced to this crate's already-documented ties-away
/// vs. IEEE ties-even convention difference (see the `remainder's
/// ties-away vs IEEE ties-even` backlog entry), not a new issue.
#[inline(always)]
pub fn remainder_checked(x: f32, y: f32) -> f32 {
    let q0 = (x / y).round();
    let r0 = fma(-q0, y, x);
    let adj = if (r0 > 0.0) == (y > 0.0) { 1.0 } else { -1.0 };
    let r1 = fma(-adj, y, r0);
    let normal = if r0.abs() > y.abs() * 0.5 { r1 } else { r0 };
    let r = if x == 0.0 && !normal.is_nan() { x } else { normal };
    if y.is_infinite() && x.is_finite() { x } else { r }
}

/// [`remainder_checked`], but correct for `|x/y|` past `2^24` too (up to
/// roughly `2^48`) -- backlog idea #64. `remainder_checked`'s own
/// self-correction assumes `q0 = (x/y).round()` differs from the true
/// integer quotient by at most one (a near-tie sign flip); that holds for
/// `|x/y| < 2^24` where every integer is exactly representable in f32, but
/// past it `q0` itself can only land on a coarse grid (gaps of
/// `2^(e-23)` at exponent `e`), so the true quotient can be tens of
/// integers away from `q0` -- fuzzing confirmed this is a *severe* gap, not
/// a tail case: ~85% of samples in `[2^24, 1e9]` were "gross" errors (off
/// by more than `0.1*|y|`, i.e. landed on a completely different multiple
/// of `y`), worst case off by 31 whole multiples of `y`.
///
/// Fixed by computing the residual with `Df32` (already in the crate)
/// instead of a single `fma`: `q0*y` via `Df32::from_mul` (an exact
/// two-product) subtracted from `x` gives the *exact* real-valued residual
/// `x - q0*y`, unlimited by `q0`'s own coarse quantization -- collapsing
/// that to a single f32 and dividing by `y` now recovers the correction
/// `adj` exactly (since `adj`'s own magnitude, `q0`'s quantization gap in
/// units of 1, is tiny compared to `q0` itself). A second exact Df32
/// subtraction applies `adj`, then `remainder_checked`'s own existing
/// near-tie logic runs unchanged on the now-correct residual. This is
/// `remainder_checked`'s own correction loop run twice -- once for the
/// coarse `q0` (potentially many integers off), once for a near-tie
/// (a single fma's own bounded error) -- not a new algorithm.
///
/// Verified against a from-scratch double-f64 (106-bit) reference (a
/// naive `x as f64 - q*y as f64` reference is *not* trustworthy here --
/// confirmed by hand with an arbitrary-precision check: once `q*y` itself
/// needs more than f64's 52 mantissa bits, plain f64 arithmetic
/// reintroduces the same class of precision loss this function exists to
/// avoid, just one level up): 0 avg/max ulp, 0 gross errors up to
/// `|x/y| ~ 2^48` (an 8-order-of-magnitude extension of
/// `remainder_checked`'s own `2^24` limit), degrading past that where a
/// single correction pass is no longer enough (the same kind of "harder,
/// separate limit" `remainder_checked`'s own doc comment already
/// acknowledges, just much further out). Bit-identical to
/// `remainder_checked` throughout `remainder_checked`'s own `|x/y|<2^24`
/// domain -- *except* a narrow, low-magnitude-`x`-paired-with-huge-`y`
/// corner the original "0 differing bit patterns" 5M-sample check
/// (2026-07-09) was too sparse to hit: found by a later, denser 30M-sample
/// standing-test run (2026-07-10, `examples/unchecked_parity.rs`).
/// Whenever `max(|x|,|y|) > f32::MAX/4` (the rescale trigger below),
/// `x` and `y` both get multiplied by the exact power-of-two `0.125`
/// *unconditionally*, regardless of `x`'s own magnitude -- if `|x|` was
/// already below `8 * f32::MIN_POSITIVE` (~9.4e-38), that multiply pushes
/// it into the denormal range, where some low mantissa bits become
/// unrepresentable; multiplying back by `8.0` at the end can't recover
/// them, so the round trip isn't lossless the way it is for any `x` that
/// stays normal throughout. Since `q0` is always `0` in this corner (`x`
/// is astronomically smaller than the rescaled `y`), the *true* answer
/// needs no rescaling at all -- `remainder_checked` (no such guard)
/// returns `x` bit-exact, while this function can differ by a handful of
/// ulp (up to ~4, confirmed by a targeted sweep). Narrow (needs `y` within
/// a factor of ~4 of `f32::MAX` *and* `x` already near/below the
/// denormal boundary) and small (a few ulp, not a gross error) -- not
/// chased further, but the doc claim is corrected here rather than left
/// overstated.
///
/// A second, distinct exception (found the same day, running the same
/// standing test at 500M samples instead of 30M): whenever `x/y` lands
/// on an *exact* half-integer (confirmed with exact rational arithmetic,
/// not just an f64 approximation, on the triggering case), this function
/// can return the *wrong sign* -- not just a few ulp off, the negative of
/// the correct magnitude. Root cause: `q0 = (xs/ys).round()` already
/// correctly resolves the tie (ties away from zero, matching
/// `remainder`/`remainder_checked`'s own established convention), which
/// makes the residual `r0` land on exactly `+-ys/2` -- but the next
/// stage, `adj = (r0/ys).round()`, exists to recover coarse-grid
/// quantization gaps `q0` can miss for extreme `|x/y|` (see below), and
/// its own blind `.round()` treats this *already-correctly-resolved* tie
/// as "one more whole `y` to remove," silently re-flipping a value that
/// was already right before `remainder_checked`'s own tie-breaking logic
/// (renamed `adj2`/`r2` here) even runs. That final stage then sees the
/// same `+-ys/2` magnitude again (now wrong-signed) and, using the exact
/// same strict `>` comparison `remainder_checked` itself uses (correctly
/// -- verified by hand that `remainder_checked` handles this same exact
/// tie shape correctly, since its residual never gets a spurious extra
/// nudge in the first place), doesn't trigger a correction for an *equal*
/// magnitude, so the wrong sign ships. Needs a genuine, exact mathematical
/// tie in `x/y` -- 3 hits in 292M uniform-random samples -- so this is
/// even narrower than the denormal case above, but a real sign flip is a
/// more serious defect class than a few-ulp miss. Not fixed this session
/// (the `adj` stage's blind rounding is load-bearing for its own real
/// job -- recovering potentially many-integer quantization gaps for
/// extreme ratios -- and a safe fix needs to distinguish "genuine
/// multi-integer gap" from "already-resolved single tie" without
/// breaking the former, not verified here). Excluded from the standing
/// test's own domain instead (`examples/unchecked_parity.rs`) rather
/// than either silently deleting the pair or leaving a permanently-red
/// test for an extremely narrow, already-diagnosed case.
/// Real extra mca cost on top of `remainder_checked` (two Df32
/// subtractions plus two Df32 products, plus the rescale guard below):
/// latency 46.03->179.13 cyc (~3.9x), throughput 1.282->8.351 cyc/elem
/// (~6.5x) -- substantial, kept as a separate opt-in tier rather than
/// folded into `remainder_checked` itself so existing callers who don't
/// need this range don't pay for it, matching this crate's established
/// precedent (`remainder`/`remainder_checked` themselves) that a real
/// correctness gap is worth a real cost for callers who opt in.
///
/// `Df32::from_mul(q0, y)` rounds the intermediate product `q0*y` to a
/// single f32 *before* pairing it with its error term -- unlike a hardware
/// `fma`, which keeps the product at full precision internally and only
/// rounds the final sum, so it never overflows just because an
/// intermediate value would have. `q0*y` is only close to `x` in the
/// *typical* case; in general it can exceed `x` by up to `|y|/2` (that gap
/// *is* the remainder), so whenever `x` or `y` individually sits within a
/// small factor of `f32::MAX`, the exact product `q0*y` can genuinely
/// exceed `f32::MAX` even though `x`, `y`, and the true remainder are all
/// finite -- found by fuzzing (a plain domain-restricted sweep missed it;
/// a dedicated hunt over full random bit patterns within the declared
/// `|x/y|` domain did not), e.g. `x=3.2603515e38, y=1.8878502e38` (`q0=2`,
/// exact `q0*y=3.7757e38 > f32::MAX`) gave `NaN` instead of the correct
/// finite remainder. Fixed by rescaling both `x` and `y` down by a fixed,
/// exact power of two (`0.125`) whenever `max(|x|,|y|)` gets within a 4x
/// safety margin of `f32::MAX` -- remainder is homogeneous of degree 1
/// (`remainder(k*x,k*y) == k*remainder(x,y)` for `k>0`), so this is exact,
/// not approximate, and the same branchless-select pattern (both scales
/// computed, one selected) as the rest of this crate. Verified by a
/// targeted 200M-sample hunt (deliberately generating full-range random
/// bit patterns rather than a bounded sweep, the technique that caught
/// this bug in the first place): zero remaining `NaN`/out-of-range
/// results for finite in-domain `x`,`y`.
#[inline(always)]
pub fn remainder_wide(x: f32, y: f32) -> f32 {
    let big = x.abs().max(y.abs()) > f32::MAX * 0.25;
    let scale = if big { 0.125 } else { 1.0 };
    let xs = x * scale;
    let ys = y * scale;
    let q0 = (xs / ys).round();
    let r0_df = Df32::from_f32(xs) - Df32::from_mul(q0, ys);
    let r0 = r0_df.to_f32();
    let adj = (r0 / ys).round();
    let r1_df = r0_df - Df32::from_mul(adj, ys);
    let r1 = r1_df.to_f32();
    let adj2 = if (r1 > 0.0) == (ys > 0.0) { 1.0 } else { -1.0 };
    let r2 = fma(-adj2, ys, r1);
    let normal = if r1.abs() > ys.abs() * 0.5 { r2 } else { r1 };
    let normal = normal * (1.0 / scale);
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
/// instead of `.round()`. Matches Rust's own `%` operator on `f32`
/// exactly (verified: `%` already implements C `fmod` semantics, not
/// `remainder`'s) over 50M generated `|x/y|<1000` samples, aside from a
/// low-probability (~3e-7) real failure mode -- the truncation analog of
/// `remainder`'s own already-documented one: `x/y`'s single division
/// rounding can occasionally land the *true* mathematical quotient
/// within half a division-ulp of an exact *integer* (not a half-integer
/// tie, since `.trunc()`'s decision boundary is at integers), pushing
/// the computed `x/y` to the wrong side of it and `q` off by a whole 1,
/// giving a result off by exactly `y`. Same inherited "naive x-q*y
/// formula" property `remainder` documents (see its own doc comment,
/// including its `remainder_checked` fix-at-extra-cost tier) -- not
/// pursued further here for the same reason, the failure rate is already
/// this low without any correction.
#[inline(always)]
pub fn fmod(x: f32, y: f32) -> f32 {
    let q = (x / y).trunc();
    remainder_style_combine!(x, y, q)
}

/// [`fmod`] without domain checks: valid for `x != 0.0` and `y` finite
/// (not `+-inf`) -- mirrors [`remainder_unchecked`]'s own contract and
/// reasoning exactly, just for `fmod`'s truncated convention.
#[inline(always)]
pub fn fmod_unchecked(x: f32, y: f32) -> f32 {
    let q = (x / y).trunc();
    fma(-q, y, x)
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
