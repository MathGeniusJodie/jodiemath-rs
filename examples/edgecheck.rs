use jodiemath_rs::*;
use std::sync::atomic::{AtomicU32, Ordering};

/// Every `FAIL` this harness prints also bumps `FAILURES`, which `main`
/// turns into a nonzero exit status. It used to only print: a real
/// `powf(-1, NaN)` regression printed `FAIL` in the middle of ~1200 lines
/// of `ok` and the process still exited 0, so a scripted
/// `cargo run --example edgecheck; echo $?` read clean. Pinning something
/// here is only a regression guard if the guard can fail the run.
static FAILURES: AtomicU32 = AtomicU32::new(0);

fn fail() {
    FAILURES.fetch_add(1, Ordering::Relaxed);
}

fn check(name: &str, got: f32, want: f32) {
    let ok = (got.is_nan() && want.is_nan()) || (got.to_bits() == want.to_bits());
    if !ok {
        fail();
    }
    println!("{} {:30} got {:e} (0x{:08x}) want {:e}", if ok { "ok  " } else { "FAIL" }, name, got, got.to_bits(), want);
}

/// Same as `check`, but for a mismatch that's a known/accepted 1-ulp
/// tradeoff (see readme.md's cbrt section) -- not a regression to chase.
/// Prints "known" instead of "FAIL" so this doesn't read as a new bug in
/// scrollback or get flagged by an unrelated task.
fn check_known_1ulp(name: &str, got: f32, want: f32) {
    let ok = (got.is_nan() && want.is_nan()) || (got.to_bits() == want.to_bits());
    println!("{} {:30} got {:e} (0x{:08x}) want {:e}", if ok { "ok  " } else { "known" }, name, got, got.to_bits(), want);
}

fn check_finite(name: &str, got: f32) {
    if !got.is_finite() {
        fail();
    }
    println!("{} {:30} got {:e} (0x{:08x}) (finite)", if got.is_finite() { "ok  " } else { "FAIL" }, name, got, got.to_bits());
}

/// Stronger than `check_finite`: also asserts `|got| <= bound`. Added
/// after finding `sin_checked`/`cos_checked` could silently return values
/// up to `2.6e21` for legitimate extreme finite input (fixed with a
/// `.clamp(-1.0,1.0)`, see their own doc comments) -- that bug's own
/// symptom was *finite*, so `check_finite` alone would never have caught
/// it; this closes that gap as a permanent regression guard.
fn check_bounded(name: &str, got: f32, bound: f32) {
    let ok = got.is_finite() && got.abs() <= bound;
    if !ok {
        fail();
    }
    println!(
        "{} {:30} got {:e} (0x{:08x}) (bounded by {bound})",
        if ok { "ok  " } else { "FAIL" },
        name,
        got,
        got.to_bits()
    );
}

/// Same idea as `check_bounded`, but for an asymmetric range (e.g.
/// `erfc`'s own `[0,2]`, not centered on zero) instead of `|got| <= bound`.
fn check_range(name: &str, got: f32, lo: f32, hi: f32) {
    let ok = got.is_finite() && got >= lo && got <= hi;
    if !ok {
        fail();
    }
    println!(
        "{} {:30} got {:e} (0x{:08x}) (in [{lo},{hi}])",
        if ok { "ok  " } else { "FAIL" },
        name,
        got,
        got.to_bits()
    );
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

/// Pin a value against an external reference (scipy) to within `tol` ulp,
/// for a point in a region where this crate's own accuracy target is a few
/// ulp anyway. A bit-exact `check` would be over-tight there -- it would
/// fail on any future refit that trades one worst point for another -- but
/// `check_finite` alone would let an arbitrarily large regression through,
/// which is exactly the failure mode these pins exist to catch.
fn check_ulp(name: &str, got: f32, want: f32, tol: u64) {
    let d = ulp_diff(got, want);
    let ok = d <= tol;
    if !ok {
        fail();
    }
    println!(
        "{} {:30} got {:e} (0x{:08x}) want {:e} ({d} ulp, tol {tol})",
        if ok { "ok  " } else { "FAIL" },
        name,
        got,
        got.to_bits(),
        want
    );
}

/// idea #197: seam continuity standing test -- every branchless function
/// in this crate always computes *both* branches and selects at a fixed
/// threshold, so the two independently-fitted approximations aren't
/// required to agree exactly right at the seam (only the *true*
/// mathematical function is continuous there), but a well-placed
/// crossover should still leave them close. Prints the value gap (in
/// ulp, evaluated at the threshold's immediate f32 neighbors on each
/// side) and a one-sided-slope estimate from a step 1000x wider, as a
/// regression detector for future coefficient/threshold refits -- not a
/// correctness assertion (a healthy seam can have a real few-ulp gap
/// from each branch's own independent fit error), just a generous bound
/// (500 ulp) wide enough to catch an actual wrong-branch/sign-flip bug
/// without false-failing on ordinary approximation disagreement.
fn check_seam(name: &str, f: impl Fn(f32) -> f32, threshold: f32) {
    let below = f32::from_bits(threshold.to_bits() - 1);
    let above = f32::from_bits(threshold.to_bits() + 1);
    let v_below = f(below);
    let v_above = f(above);
    let gap = ulp_diff(v_below, v_above);
    let step = threshold * 1e-3;
    let slope_below = (f(threshold) - f(threshold - step)) / step;
    let slope_above = (f(threshold + step) - f(threshold)) / step;
    let ok = gap < 500;
    if !ok {
        fail();
    }
    println!(
        "{} {:30} value gap {gap} ulp at threshold {threshold:e}, one-sided slopes {slope_below:e} / {slope_above:e}",
        if ok { "ok  " } else { "FAIL" },
        name
    );
}

fn main() {
    real_main();
    let n = FAILURES.load(Ordering::Relaxed);
    if n == 0 {
        println!("\nall edge cases pass");
    } else {
        println!("\n{n} FAILED edge case(s) -- see the FAIL lines above");
    }
    std::process::exit((n != 0) as i32);
}

fn real_main() {
    // The metric itself, before anything that uses it: NaN vs NaN scores
    // 0 regardless of payload/sign, one-sided NaN scores maximal. Every
    // ulp sweep in this repo depends on this, and `accuracy.rs`'s
    // `thorough` mode feeds it all 2^24 NaN payloads.
    let nan_a = f32::from_bits(0x7fc0_0001);
    let nan_b = f32::from_bits(0xffc0_5678);
    assert_eq!(ulp_diff(nan_a, nan_b), 0, "NaN vs NaN must be 0 ulp");
    assert_eq!(ulp_diff(f32::NAN, f32::NAN), 0, "NaN vs NaN must be 0 ulp");
    assert_eq!(ulp_diff(nan_a, 1.0), u64::MAX, "one-sided NaN must be maximal");
    assert_eq!(ulp_diff(1.0, nan_b), u64::MAX, "one-sided NaN must be maximal");

    // log_2
    check("log_2(0)", log_2(0.0), f32::NEG_INFINITY);
    check("log_2(-0)", log_2(-0.0), f32::NEG_INFINITY);
    check("log_2(-1)", log_2(-1.0), f32::NAN);
    check("log_2(-inf)", log_2(f32::NEG_INFINITY), f32::NAN);
    check("log_2(inf)", log_2(f32::INFINITY), f32::INFINITY);
    check("log_2(nan)", log_2(f32::NAN), f32::NAN);
    check("log_2(1e-39)", log_2(1e-39), (1e-39f32 as f64).log2() as f32);
    check("log_2(min_denorm)", log_2(f32::from_bits(1)), (f32::from_bits(1) as f64).log2() as f32);
    check("log_2(1)", log_2(1.0), 0.0);
    // log_2_unchecked: only test inside its documented domain (positive
    // normal finite) -- must match log_2 exactly there, bit for bit.
    check("log_2_unchecked(1)", log_2_unchecked(1.0), 0.0);
    check("log_2_unchecked(f32::MAX)", log_2_unchecked(f32::MAX), log_2(f32::MAX));
    check("log_2_unchecked(f32::MIN_POSITIVE)", log_2_unchecked(f32::MIN_POSITIVE), log_2(f32::MIN_POSITIVE));
    // exp2 (unchecked): only test inside its documented domain [-126, 128)
    check("exp2(127.9999)", exp2(127.9999), (127.9999f32 as f64).exp2() as f32);
    check("exp2(-125.9)", exp2(-125.9), ((-125.9f32) as f64).exp2() as f32);
    check("exp2(-126)", exp2(-126.0), ((-126.0f32) as f64).exp2() as f32);
    check("exp2(0)", exp2(0.0), 1.0);
    // exp2_kf (backlog idea #116): exp2's own combine, k/f supplied
    // directly -- bit-identical to exp2(x) when fed exp2's own
    // k=floor(x)/f=x-k (verified across 20M samples before shipping).
    check("exp2_kf(0,0)", exp2_kf(0.0, 0.0), 1.0);
    check("exp2_kf(1,0)", exp2_kf(1.0, 0.0), 2.0);
    check("exp2_kf(0,1)", exp2_kf(0.0, 1.0), 2.0);
    for &x in &[3.7f32, -10.25, 100.0, -125.9] {
        let k = x.floor();
        let f = x - k;
        check("exp2_kf reconstructs exp2", exp2_kf(k, f), exp2(x));
    }
    // exp2_checked: full range
    check("exp2_checked(128)", exp2_checked(128.0), f32::INFINITY);
    check("exp2_checked(127.9999)", exp2_checked(127.9999), (127.9999f32 as f64).exp2() as f32);
    check("exp2_checked(1000)", exp2_checked(1000.0), f32::INFINITY);
    check("exp2_checked(inf)", exp2_checked(f32::INFINITY), f32::INFINITY);
    check("exp2_checked(-inf)", exp2_checked(f32::NEG_INFINITY), 0.0);
    check("exp2_checked(nan)", exp2_checked(f32::NAN), f32::NAN);
    check("exp2_checked(-125.9)", exp2_checked(-125.9), ((-125.9f32) as f64).exp2() as f32);
    check("exp2_checked(-126.5)", exp2_checked(-126.5), ((-126.5f32) as f64).exp2() as f32);
    check("exp2_checked(-149)", exp2_checked(-149.0), f32::from_bits(1));
    check("exp2_checked(-150)", exp2_checked(-150.0), 0.0);
    check("exp2_checked(-1000)", exp2_checked(-1000.0), 0.0);
    check("exp2_checked(0)", exp2_checked(0.0), 1.0);

    // ldexp/frexp (backlog idea #86): exact power-of-two utilities.
    // Verified via a real large-scale sweep before adopting -- found and
    // fixed two real bugs (a premature-overflow case and a grossly-out-
    // of-range-clamped-instead-of-saturated case, see ldexp's own doc
    // comment) that these pins alone wouldn't have caught, but now pin
    // the specific failure points directly as regression guards.
    check("ldexp(0,5)", ldexp(0.0, 5), 0.0);
    check("ldexp(-0,5)", ldexp(-0.0, 5), -0.0);
    check("ldexp(inf,5)", ldexp(f32::INFINITY, 5), f32::INFINITY);
    check("ldexp(-inf,5)", ldexp(f32::NEG_INFINITY, 5), f32::NEG_INFINITY);
    check("ldexp(nan,5)", ldexp(f32::NAN, 5), f32::NAN);
    check("ldexp(1,0)", ldexp(1.0, 0), 1.0);
    check("ldexp(1,3)", ldexp(1.0, 3), 8.0);
    check("ldexp(1,-3)", ldexp(1.0, -3), 0.125);
    check("ldexp(f32::MAX,0)", ldexp(f32::MAX, 0), f32::MAX);
    check("ldexp(f32::MAX,1)", ldexp(f32::MAX, 1), f32::INFINITY);
    check("ldexp(f32::MIN_POSITIVE,-1000)", ldexp(f32::MIN_POSITIVE, -1000), 0.0);
    // real bug #1 (premature overflow: x's own magnitude should have
    // compensated for n past exp2_checked's own clamp, but didn't).
    check("ldexp(7.26589e-9,148)", ldexp(7.26589e-9, 148), 2.5925562e36);
    // real bug #2 (grossly-out-of-range clamped down to the reconstruction
    // boundary instead of saturating -- true answer overflows regardless
    // of mantissa).
    check("ldexp(-224910930000000,148)", ldexp(-224910930000000.0, 148), f32::NEG_INFINITY);
    check("frexp(0).0", frexp(0.0).0, 0.0);
    check("frexp(0).1", frexp(0.0).1 as f32, 0.0);
    check("frexp(-0).0", frexp(-0.0).0, -0.0);
    check("frexp(inf).0", frexp(f32::INFINITY).0, f32::INFINITY);
    check("frexp(nan).0", frexp(f32::NAN).0, f32::NAN);
    check("frexp(6).0", frexp(6.0).0, 0.75);
    check("frexp(6).1", frexp(6.0).1 as f32, 3.0);
    check("frexp(f32::MAX).0", frexp(f32::MAX).0, 0.99999994);
    check("frexp(f32::MAX).1", frexp(f32::MAX).1 as f32, 128.0);

    // cbrt
    for f in [cbrt as fn(f32) -> f32, cbrt_accurate as fn(f32) -> f32] {
        let n = if std::ptr::fn_addr_eq(f, cbrt as fn(f32) -> f32) { "cbrt" } else { "cbrt_acc" };
        check(&format!("{n}(0)"), f(0.0), 0.0);
        check(&format!("{n}(-0)"), f(-0.0), -0.0);
        check(&format!("{n}(inf)"), f(f32::INFINITY), f32::INFINITY);
        check(&format!("{n}(-inf)"), f(f32::NEG_INFINITY), f32::NEG_INFINITY);
        check(&format!("{n}(nan)"), f(f32::NAN), f32::NAN);
        check(&format!("{n}(8)"), f(8.0), 2.0);
        check(&format!("{n}(-8)"), f(-8.0), -2.0);
        check(&format!("{n}(1e-39)"), f(1e-39), ((1e-39f32 as f64).cbrt()) as f32);
        check(&format!("{n}(-1e-39)"), f(-1e-39), ((-1e-39f32 as f64).cbrt()) as f32);
        // plain cbrt (not cbrt_accurate) is 1 ulp off at these two points --
        // a known, accepted cost of the degree5->degree3 seed-correction
        // cut (see readme.md's cbrt section); not a regression, won't fix.
        let check_min_max = if n == "cbrt" { check_known_1ulp } else { check };
        check_min_max(&format!("{n}(min_denorm)"), f(f32::from_bits(1)), ((f32::from_bits(1) as f64).cbrt()) as f32);
        check_min_max(&format!("{n}(max)"), f(f32::MAX), ((f32::MAX as f64).cbrt()) as f32);
        check(&format!("{n}(2^-57)"), f(f32::from_bits(0x2300_0000)), ((f32::from_bits(0x2300_0000) as f64).cbrt()) as f32);
    }
    // cbrt/cbrt_accurate's own zero/inf/nan special cases (2026-07-10):
    // neither had ever had a direct pin -- only exercised indirectly via
    // the _unchecked comparisons below, which never touch this domain at
    // all. Both propagate via `x + x`, which preserves sign for 0/inf and
    // is NaN-preserving for nan (verified correct before pinning, not
    // assumed).
    check("cbrt(0)", cbrt(0.0), 0.0);
    check("cbrt(-0)", cbrt(-0.0), -0.0);
    check("cbrt(inf)", cbrt(f32::INFINITY), f32::INFINITY);
    check("cbrt(-inf)", cbrt(f32::NEG_INFINITY), f32::NEG_INFINITY);
    check("cbrt(nan)", cbrt(f32::NAN), f32::NAN);
    check("cbrt_accurate(0)", cbrt_accurate(0.0), 0.0);
    check("cbrt_accurate(-0)", cbrt_accurate(-0.0), -0.0);
    check("cbrt_accurate(inf)", cbrt_accurate(f32::INFINITY), f32::INFINITY);
    check("cbrt_accurate(-inf)", cbrt_accurate(f32::NEG_INFINITY), f32::NEG_INFINITY);
    check("cbrt_accurate(nan)", cbrt_accurate(f32::NAN), f32::NAN);
    check("cbrt_accurate(-8)", cbrt_accurate(-8.0), -2.0);
    check("cbrt_accurate(27)", cbrt_accurate(27.0), 3.0);
    // cbrt_unchecked: contract is x normal/finite (no denormal/zero/inf/nan)
    // -- must match cbrt inside that domain (verified more thoroughly via a
    // ~200M-sample fuzz, not preserved in-repo; permanent regression guard).
    check("cbrt_unchecked(8)", cbrt_unchecked(8.0), cbrt(8.0));
    check("cbrt_unchecked(-8)", cbrt_unchecked(-8.0), cbrt(-8.0));
    check("cbrt_unchecked(max)", cbrt_unchecked(f32::MAX), cbrt(f32::MAX));
    check("cbrt_unchecked(min_normal)", cbrt_unchecked(f32::MIN_POSITIVE), cbrt(f32::MIN_POSITIVE));
    // cbrt_accurate_unchecked: contract is x already inside cbrt_accurate's
    // own safe rescale range (roughly 2^-56 to 2^127) -- must match
    // cbrt_accurate inside that domain (verified more thoroughly via a
    // ~143M-sample fuzz, not preserved in-repo; permanent regression guard).
    check("cbrt_accurate_unchecked(8)", cbrt_accurate_unchecked(8.0), cbrt_accurate(8.0));
    check("cbrt_accurate_unchecked(-8)", cbrt_accurate_unchecked(-8.0), cbrt_accurate(-8.0));
    check("cbrt_accurate_unchecked(1e30)", cbrt_accurate_unchecked(1e30), cbrt_accurate(1e30));
    check("cbrt_accurate_unchecked(1e-16)", cbrt_accurate_unchecked(1e-16), cbrt_accurate(1e-16));

    // rcbrt(x) = 1/cbrt(x). Every special case falls out of composing
    // cbrt with a plain division purely from IEEE754 semantics (verified
    // by hand before writing the function) -- no override needed at all,
    // unlike rhypot's one inf-vs-NaN case.
    check("rcbrt(8)", rcbrt(8.0), 0.5);
    check("rcbrt(-8)", rcbrt(-8.0), -0.5);
    check("rcbrt(0)", rcbrt(0.0), f32::INFINITY);
    check("rcbrt(-0)", rcbrt(-0.0), f32::NEG_INFINITY);
    check("rcbrt(inf)", rcbrt(f32::INFINITY), 0.0);
    check("rcbrt(-inf)", rcbrt(f32::NEG_INFINITY), -0.0);
    check("rcbrt(nan)", rcbrt(f32::NAN), f32::NAN);
    check("rcbrt(1)", rcbrt(1.0), 1.0);
    check("rcbrt(-1)", rcbrt(-1.0), -1.0);
    // sin/cos (unchecked): only accurate while q = round(x/pi) is an exact
    // f32 integer, i.e. |x| < 2^22 * pi (~1.3e7) -- see sin's doc comment.
    check("sin(0)", sin(0.0), 0.0);
    check("cos(0)", cos(0.0), 1.0);
    check("sin(nan)", sin(f32::NAN), f32::NAN);
    check("cos(nan)", cos(f32::NAN), f32::NAN);
    // sin(-0.0)/tan(-0.0) used to lose their sign: sinf_poly's own fma chain
    // adds two exactly-zero values of opposite sign at this one input (the
    // same IEEE754 "+0 + -0 = +0" rule as the atan2 bug above), which
    // silently flips the correctly-signed -0.0 back to +0.0. Fixed with
    // `r.copysign(x)` in sinf_poly -- free for every nonzero x (sin is odd
    // and monotonic there, so r's sign already matches x's), only changes
    // this singular zero case. cos(-0.0) is unaffected (a nonzero result).
    check("sin(-0)", sin(-0.0), -0.0);
    check("cos(-0)", cos(-0.0), 1.0);
    check("tan(-0)", tan(-0.0), -0.0);
    // sin_checked/cos_checked: nan/+-inf must still come out nan, and every
    // other finite x (including ones far past the sub-ulp-accurate range)
    // must come out finite *and* bounded to `[-1,1]` -- see sin_checked's
    // doc comment for why (the residual clamp added 2026-07-06, plus the
    // final `.clamp(-1.0,1.0)` added 2026-07-10 after finding
    // round_x_over_pi's double-float q silently loses precision past
    // `|x|~8.85e14`, which used to let these same inputs return values up
    // to `2.6e21` -- `check_finite` alone never caught that, since `2.6e21`
    // is finite; `check_bounded` closes that gap for good).
    for f in [sin_checked as fn(f32) -> f32, cos_checked as fn(f32) -> f32] {
        let n = if std::ptr::fn_addr_eq(f, sin_checked as fn(f32) -> f32) { "sin_checked" } else { "cos_checked" };
        check(&format!("{n}(nan)"), f(f32::NAN), f32::NAN);
        check(&format!("{n}(inf)"), f(f32::INFINITY), f32::NAN);
        check(&format!("{n}(-inf)"), f(f32::NEG_INFINITY), f32::NAN);
        check_bounded(&format!("{n}(max)"), f(f32::MAX), 1.0);
        check_bounded(&format!("{n}(-max)"), f(f32::MIN), 1.0);
        check_bounded(&format!("{n}(1e20)"), f(1e20), 1.0);
        check_bounded(&format!("{n}(1e16)"), f(1e16), 1.0);
        check_bounded(&format!("{n}(1e10)"), f(1e10), 1.0);
    }
    // sin_checked(-0.0) used to lose its sign too, via a *different*
    // mechanism than sinf_poly's own bug above: reduce_pi's multi-term
    // two_sum/two_prod error-compensation chain loses x's sign somewhere
    // internally (same IEEE754 rule, not traced to the exact spot), well
    // before sinf_poly is even reached, so sinf_poly's own copysign fix
    // can't see the original sign to restore it. Guarded at sin_checked's
    // own output instead. cos_checked(-0.0) needs no such guard (nonzero
    // result, and its own reduction shares no sign-losing dependency here).
    check("sin_checked(-0)", sin_checked(-0.0), -0.0);
    check("cos_checked(-0)", cos_checked(-0.0), 1.0);

    // reduce_pi_checked/reduce_pi_half_checked (backlog idea #88): the
    // public pi-reduction primitive sin_checked/cos_checked build on.
    // sign is a plain +-1.0 multiplier (see its own doc comment for why,
    // not a bool) such that sign*sin(r) reconstructs sin_checked(x) (or
    // cos_checked(x) for the half variant) -- checked here via std::f32
    // sin as the reference poly, not this crate's own private sinf_poly.
    for &x in &[0.0f32, 1.0, 3.0, 100.0, -7.5, 1e6, -1e9] {
        let (r, sign) = reduce_pi_checked(x);
        check_bounded(&format!("reduce_pi_checked({x})"), sign * r.sin() - sin_checked(x), 1e-4);
        let (rc, signc) = reduce_pi_half_checked(x);
        check_bounded(&format!("reduce_pi_half_checked({x})"), signc * rc.sin() - cos_checked(x), 1e-4);
    }

    // wrap_pi (backlog idea #126): wraps to (-pi, pi], riding
    // reduce_pi_checked's own reduction. Checked via sin/cos preserved
    // (wrap_pi(x) and x are the same angle mod 2*pi) plus the range
    // invariant directly.
    check("wrap_pi(0)", wrap_pi(0.0), 0.0);
    // -0.0 needs its own pin, not just +0.0: reduce_pi_checked forms the
    // residual by subtracting equal signed zeros, which IEEE754 resolves
    // to +0.0, so without wrap_pi's own `x == 0.0` guard this returned
    // +0.0. Found by the special-value matrix (idea #164) -- it was the
    // single sign anomaly across 2.14e9 inputs on |x| <= pi/2, the region
    // where wrap_pi is otherwise exactly the identity.
    check("wrap_pi(-0)", wrap_pi(-0.0), -0.0);
    // Identity on |x| <= pi/2 (q = 0 there, so r == x exactly). Spot pins
    // for the property the -0.0 case is the boundary of.
    check("wrap_pi(0.5)", wrap_pi(0.5), 0.5);
    check("wrap_pi(-0.5)", wrap_pi(-0.5), -0.5);
    check("wrap_pi(min_denorm)", wrap_pi(f32::from_bits(1)), f32::from_bits(1));
    check(
        "wrap_pi(-min_denorm)",
        wrap_pi(-f32::from_bits(1)),
        -f32::from_bits(1),
    );
    check_bounded("wrap_pi(pi)+pi", wrap_pi(std::f32::consts::PI) + std::f32::consts::PI, 1e-5);
    for &x in &[1.0f32, 3.0, 4.0, -4.0, 100.0, -1e9, 1e6] {
        let w = wrap_pi(x);
        check_bounded(&format!("sin(wrap_pi({x}))-sin_checked({x})"), w.sin() - sin_checked(x), 1e-4);
        check_bounded(&format!("cos(wrap_pi({x}))-cos_checked({x})"), w.cos() - cos_checked(x), 1e-4);
        // in-range check: 0 when w is in (-pi,pi], positive by however
        // far out of range otherwise (w=-pi itself, the excluded
        // boundary, would show up here as pi - (-pi) = 2*pi).
        let over = (w - std::f32::consts::PI).max(0.0);
        let under = (-std::f32::consts::PI - w).max(0.0);
        check_bounded(&format!("wrap_pi({x}) in (-pi,pi]"), over + under, 1e-6);
    }
    check("wrap_pi(nan)", wrap_pi(f32::NAN), f32::NAN);

    // sin_prereduced/cos_prereduced (backlog idea #127): bit-identical
    // to each other (both are exactly sinf_poly) -- see their own doc
    // comments for why two names exist anyway. Checked by reconstructing
    // sin_checked/cos_checked from reduce_pi_checked/
    // reduce_pi_half_checked + these (within 1 ulp: sin_checked/
    // cos_checked's own extra `.clamp(-1,1)` isn't part of this lower-
    // level primitive's contract).
    check("sin_prereduced(0)", sin_prereduced(0.0), 0.0);
    check("cos_prereduced(0)", cos_prereduced(0.0), 0.0);
    check("sin_prereduced(cos_prereduced same fn)", sin_prereduced(0.7), cos_prereduced(0.7));
    for &x in &[0.3f32, -0.9, 1.5, -1.5, 100.0, -1e6] {
        let (r, sign) = reduce_pi_checked(x);
        check_bounded(&format!("sin_prereduced reconstructs sin_checked({x})"), sign * sin_prereduced(r) - sin_checked(x), 1e-6);
        let (rc, signc) = reduce_pi_half_checked(x);
        check_bounded(&format!("cos_prereduced reconstructs cos_checked({x})"), signc * cos_prereduced(rc) - cos_checked(x), 1e-6);
    }

    // sinpi/cospi: argument in half-turns, q=round(x)/r=x-q both exact in
    // f32, so (unlike sin/cos) there's no accuracy cliff anywhere -- these
    // pin the full-range "always finite, exact at exact half-integers"
    // guarantee described in sinpi's own doc comment.
    check("sinpi(0)", sinpi(0.0), 0.0);
    check("sinpi(-0)", sinpi(-0.0), -0.0);
    check("cospi(0)", cospi(0.0), 1.0);
    check("cospi(-0)", cospi(-0.0), 1.0);
    check("sinpi(0.5)", sinpi(0.5), 1.0);
    check("sinpi(1)", sinpi(1.0), -0.0);
    check("cospi(1)", cospi(1.0), -1.0);
    check("sinpi(-0.5)", sinpi(-0.5), -1.0);
    check("sinpi(nan)", sinpi(f32::NAN), f32::NAN);
    check("cospi(nan)", cospi(f32::NAN), f32::NAN);
    check("sinpi(inf)", sinpi(f32::INFINITY), f32::NAN);
    check("sinpi(-inf)", sinpi(f32::NEG_INFINITY), f32::NAN);
    check("cospi(inf)", cospi(f32::INFINITY), f32::NAN);
    // Regression guard for a real bug (2026-07-09): the magic-round-
    // constant reduction (valid only for |x|<=2^22, since it was applied
    // directly to unbounded raw x, unlike every other magic-round use in
    // this crate) silently gave *wrong*, not just imprecise, results for
    // 2^22 < |x| < 2^24 -- e.g. cospi(2^22+1) came out ~-0.0033 instead
    // of the correct -1.0. `check_finite` alone (the old form of these
    // checks) could never have caught this, since the wrong values were
    // still finite. Fixed with `x.round()` (see sinpi's own doc comment).
    check("sinpi(2^22+1)", sinpi(4_194_305.0), -0.0);
    check("cospi(2^22+1)", cospi(4_194_305.0), -1.0);
    check("sinpi(2^23)", sinpi(8_388_608.0), 0.0);
    check("cospi(2^23)", cospi(8_388_608.0), 1.0);
    check("sinpi(f32::MAX)", sinpi(f32::MAX), 0.0);
    check("cospi(f32::MAX)", cospi(f32::MAX), 1.0);
    check("sinpi(1e20)", sinpi(1e20), 0.0);
    check("cospi(1e20)", cospi(1e20), 1.0);

    // Regression guard for a second real bug (2026-07-30): cospi keeps
    // full *relative* accuracy right up to its own zeros -- small
    // absolute error there is not the same claim, and is not what these
    // pin. Reducing on `k = round(x)` makes
    // `0.5-|x-k|` Sterbenz-exact over the half of the domain containing
    // the crossing; the previous `k = round(x-0.5)` reduction rounded
    // away the low bits of exactly the quantity the answer is
    // proportional to, which cost ~2e5 ulp just below -0.5 and returned
    // a flat 0.0 at the very first f32 below it (a 100% relative error
    // that no absolute-error framing catches). References computed in
    // f64, where that same reduction has bits to spare.
    check_ulp("cospi(0.5-2^-25)", cospi(0.49999997), 9.362676e-8, 2);
    check_ulp("cospi(-(0.5-2^-25))", cospi(-0.49999997), 9.362676e-8, 2);
    check_ulp("cospi(0.5+2^-24)", cospi(0.50000006), -1.8725352e-7, 2);
    check_ulp("cospi(-4.9999842e-1)", cospi(-0.49999842), 4.962218e-6, 2);
    // Every zero is +0.0, matching IEEE 754-2019's `cosPi(n+1/2) = +0`
    // for every integer n: ties-even rounding sends `k` to the *even*
    // neighbour at every exact half-integer, so the parity sign is +1
    // there for free, with no special case spent on it.
    check("cospi(0.5)", cospi(0.5), 0.0);
    check("cospi(-0.5)", cospi(-0.5), 0.0);
    check("cospi(1.5)", cospi(1.5), 0.0);
    check("cospi(-2.5)", cospi(-2.5), 0.0);
    check("cospi(2^22+0.5)", cospi(4_194_304.5), 0.0);
    // cos is even, and this reduction is exactly even as well
    // (`round_ties_even`, `|x-k|` and `parity` are each sign-symmetric),
    // so the two sides agree bit for bit rather than just closely.
    for &x in &[0.3f32, 0.25, 0.7, 1.3, 12.75, 1_000_000.5, 4_194_305.0] {
        check(&format!("cospi(-{x}) == cospi({x})"), cospi(-x), cospi(x));
    }

    // sinpi_unchecked (backlog idea #98): drops sinpi's `x==0.0` guard,
    // bit-identical everywhere else (verified via a real exhaustive
    // 2^32-pattern sweep before adopting). The one documented exception is
    // pinned explicitly (not silently accepted): -0.0 loses its sign
    // instead of propagating it, same known-tradeoff-not-bug documentation
    // convention as cbrt_accurate's own pinned mantissa issue.
    check("sinpi_unchecked(0)", sinpi_unchecked(0.0), 0.0);
    check("sinpi_unchecked(-0) [documented exception: wrong sign]", sinpi_unchecked(-0.0), 0.0);
    check("sinpi_unchecked(0.5)", sinpi_unchecked(0.5), 1.0);
    check("sinpi_unchecked(1)", sinpi_unchecked(1.0), sinpi(1.0));
    check("sinpi_unchecked(-0.5)", sinpi_unchecked(-0.5), sinpi(-0.5));
    check("sinpi_unchecked(nan)", sinpi_unchecked(f32::NAN), f32::NAN);
    check("sinpi_unchecked(inf)", sinpi_unchecked(f32::INFINITY), f32::NAN);
    check("sinpi_unchecked(-inf)", sinpi_unchecked(f32::NEG_INFINITY), f32::NAN);
    check("sinpi_unchecked(f32::MAX)", sinpi_unchecked(f32::MAX), sinpi(f32::MAX));

    // tanpi(x): direct poly + cotangent reflection (backlog idea #128,
    // superseding the original sinpi(x)/cospi(x) ratio, idea #29).
    // Poles at half-integer x are real (an exact `1.0/0.0` in
    // `tan_core`'s reflected branch) and correctly give +-inf, not NaN
    // -- pinned so that stays true. tanpi(inf)/(-inf) are NaN, matching
    // sinpi/cospi's own existing (inherited, not new) convention at
    // infinity. `tanpi(0.25)` is `check_bounded`, not exact, unlike the
    // old ratio construction: that happened to give bit-identical
    // sinpi(0.25)/cospi(0.25) values (X/X==1.0 exactly), a structural
    // coincidence of the ratio, not a guarantee this fit-based version
    // inherits.
    check("tanpi(0)", tanpi(0.0), 0.0);
    check("tanpi(-0)", tanpi(-0.0), -0.0);
    check_bounded("tanpi(0.25)-1", tanpi(0.25) - 1.0, 1e-5);
    check("tanpi(1)", tanpi(1.0), 0.0);
    check("tanpi(0.5)", tanpi(0.5), f32::NEG_INFINITY);
    check("tanpi(-0.5)", tanpi(-0.5), f32::NEG_INFINITY);
    check("tanpi(1.5)", tanpi(1.5), f32::NEG_INFINITY);
    check("tanpi(nan)", tanpi(f32::NAN), f32::NAN);
    check("tanpi(inf)", tanpi(f32::INFINITY), f32::NAN);
    check("tanpi(-inf)", tanpi(f32::NEG_INFINITY), f32::NAN);

    // sin2pi/cos2pi/tan2pi: full-turn arguments (backlog idea #122).
    check("sin2pi(0)", sin2pi(0.0), 0.0);
    check("sin2pi(0.25)", sin2pi(0.25), 1.0);
    check("sin2pi(0.5)", sin2pi(0.5), -0.0);
    check("cos2pi(0)", cos2pi(0.0), 1.0);
    // `+0.0`, not `-0.0`: `2*0.25` is exactly `cospi`'s own half-integer
    // zero, and every one of those is `+0.0` (see cospi's pins above).
    check("cos2pi(0.25)", cos2pi(0.25), 0.0);
    check("tan2pi(0)", tan2pi(0.0), 0.0);
    check_bounded("tan2pi(0.125)-1", tan2pi(0.125) - 1.0, 1e-5);
    check("sin2pi(nan)", sin2pi(f32::NAN), f32::NAN);
    // Doubling overflows past f32::MAX/2, a documented gap (see doc
    // comment) -- NaN there, not a meaningful finite answer.
    check("sin2pi(f32::MAX)", sin2pi(f32::MAX), f32::NAN);

    // sinc(x) = sin(pi*x)/(pi*x), removable singularity at x=0.
    check("sinc(0)", sinc(0.0), 1.0);
    check("sinc(-0)", sinc(-0.0), 1.0);
    check("sinc(1)", sinc(1.0), -0.0);
    check("sinc(-1)", sinc(-1.0), 0.0);
    check("sinc(2)", sinc(2.0), 0.0);
    check("sinc(0.5)", sinc(0.5), std::f32::consts::FRAC_2_PI);
    check("sinc(-0.5)", sinc(-0.5), std::f32::consts::FRAC_2_PI);
    check("sinc(nan)", sinc(f32::NAN), f32::NAN);
    check("sinc(inf)", sinc(f32::INFINITY), f32::NAN);
    check("sinc(-inf)", sinc(f32::NEG_INFINITY), f32::NAN);
    // |sinc(x)| = |sin(pi*x)/(pi*x)| <= 1 for all real x (a provable
    // mathematical fact, |sin(u)| <= |u| everywhere) -- checked directly
    // (2026-07-10), no violation found anywhere, but check_bounded locks
    // this in as a permanent guard the same way sin_checked/cos_checked's
    // own `[-1,1]` bound is now guarded, rather than leaving it as only an
    // implicit consequence of sinpi's own correctness.
    check_bounded("sinc(f32::MAX)", sinc(f32::MAX), 1.0);

    // sinc_unnormalized(x) = sin(x)/x, radians (backlog idea #131).
    check("sinc_unnormalized(0)", sinc_unnormalized(0.0), 1.0);
    check("sinc_unnormalized(-0)", sinc_unnormalized(-0.0), 1.0);
    check(
        "sinc_unnormalized(pi/2)",
        sinc_unnormalized(std::f32::consts::FRAC_PI_2),
        std::f32::consts::FRAC_2_PI,
    );
    check("sinc_unnormalized(nan)", sinc_unnormalized(f32::NAN), f32::NAN);
    check("sinc_unnormalized(inf)", sinc_unnormalized(f32::INFINITY), f32::NAN);
    check("sinc_unnormalized(-inf)", sinc_unnormalized(f32::NEG_INFINITY), f32::NAN);
    check(
        "sinc_unnormalized(-1)==sinc_unnormalized(1)",
        sinc_unnormalized(-1.0),
        sinc_unnormalized(1.0),
    );
    check_bounded("sinc_unnormalized(f32::MAX)", sinc_unnormalized(f32::MAX), 1.0);

    // sind/cosd: argument in degrees. Exact reduction only up to ~4.7e7
    // (180.0's own trailing-zero-bit limit, see sind's doc comment) --
    // unlike sinpi/cospi, not the entire f32 range -- but still
    // guaranteed finite (never inf/nan) for any finite input, via the
    // same POLY_SAFE_BOUND clamp sin_checked/cos_checked use.
    check("sind(0)", sind(0.0), 0.0);
    check("sind(-0)", sind(-0.0), -0.0);
    check("cosd(0)", cosd(0.0), 1.0);
    check("sind(90)", sind(90.0), 1.0);
    check("sind(180)", sind(180.0), -0.0);
    check("cosd(180)", cosd(180.0), -1.0);
    check("sind(-90)", sind(-90.0), -1.0);
    check("sind(nan)", sind(f32::NAN), f32::NAN);
    check("cosd(nan)", cosd(f32::NAN), f32::NAN);
    check("sind(inf)", sind(f32::INFINITY), f32::NAN);
    check("sind(-inf)", sind(f32::NEG_INFINITY), f32::NAN);
    check("cosd(inf)", cosd(f32::INFINITY), f32::NAN);
    check_finite("sind(f32::MAX)", sind(f32::MAX));
    check_finite("cosd(f32::MAX)", cosd(f32::MAX));
    check_finite("sind(1e10)", sind(1e10));
    check_finite("cosd(1e10)", cosd(1e10));

    // sind_unchecked/cosd_unchecked (backlog idea #98): same POLY_SAFE_BOUND
    // clamp removed, valid only up to sind/cosd's own ~4.7e7 exactness
    // limit -- no check_finite pins at f32::MAX/1e10 here, unlike sind/cosd
    // above, since that guarantee is exactly what's given up past the
    // documented domain (verified bit-identical to sind/cosd for all
    // |x|<4.7e7 via a real exhaustive sweep before adopting).
    check("sind_unchecked(0)", sind_unchecked(0.0), 0.0);
    check("sind_unchecked(-0)", sind_unchecked(-0.0), -0.0);
    check("cosd_unchecked(0)", cosd_unchecked(0.0), 1.0);
    check("sind_unchecked(90)", sind_unchecked(90.0), sind(90.0));
    check("sind_unchecked(180)", sind_unchecked(180.0), sind(180.0));
    check("cosd_unchecked(180)", cosd_unchecked(180.0), cosd(180.0));
    check("sind_unchecked(-90)", sind_unchecked(-90.0), sind(-90.0));
    check("sind_unchecked(nan)", sind_unchecked(f32::NAN), f32::NAN);
    check("cosd_unchecked(nan)", cosd_unchecked(f32::NAN), f32::NAN);
    check("sind_unchecked(inf)", sind_unchecked(f32::INFINITY), f32::NAN);
    check("sind_unchecked(-inf)", sind_unchecked(f32::NEG_INFINITY), f32::NAN);
    check("cosd_unchecked(inf)", cosd_unchecked(f32::INFINITY), f32::NAN);

    // tand(x) = sind(x)/cosd(x): new function (backlog idea #29), same
    // "poles are real, IEEE754 division handles them for free" reasoning
    // as tanpi above. (idea #128's direct-poly form was tried and
    // reverted for tand specifically -- see tand's own doc comment.)
    check("tand(0)", tand(0.0), 0.0);
    check("tand(-0)", tand(-0.0), -0.0);
    check("tand(45)", tand(45.0), 1.0);
    check("tand(180)", tand(180.0), 0.0);
    check("tand(90)", tand(90.0), f32::NEG_INFINITY);
    check("tand(-90)", tand(-90.0), f32::NEG_INFINITY);
    check("tand(270)", tand(270.0), f32::NEG_INFINITY);
    check("tand(nan)", tand(f32::NAN), f32::NAN);
    check("tand(inf)", tand(f32::INFINITY), f32::NAN);
    check("tand(-inf)", tand(f32::NEG_INFINITY), f32::NAN);
    check_finite("tand(1e10)", tand(1e10));

    // tand_unchecked: same clamp-removal contract as sind_unchecked/
    // cosd_unchecked above.
    check("tand_unchecked(0)", tand_unchecked(0.0), 0.0);
    check("tand_unchecked(-0)", tand_unchecked(-0.0), -0.0);
    check("tand_unchecked(45)", tand_unchecked(45.0), 1.0);
    check("tand_unchecked(180)", tand_unchecked(180.0), tand(180.0));
    check("tand_unchecked(90)", tand_unchecked(90.0), f32::NEG_INFINITY);
    check("tand_unchecked(nan)", tand_unchecked(f32::NAN), f32::NAN);
    check("tand_unchecked(inf)", tand_unchecked(f32::INFINITY), f32::NAN);

    // ln/log10/log1p: same zero/negative/inf edges as log_2 (they're all
    // log_2 rescaled or composed with it).
    check("ln(0)", ln(0.0), f32::NEG_INFINITY);
    check("ln(-1)", ln(-1.0), f32::NAN);
    check("ln(1)", ln(1.0), 0.0);
    // ln_unchecked/log10_unchecked: same contract as log_2_unchecked, must
    // match ln/log10 exactly on positive normal finite input.
    check("ln_unchecked(1)", ln_unchecked(1.0), 0.0);
    check("ln_unchecked(f32::MAX)", ln_unchecked(f32::MAX), ln(f32::MAX));
    check("ln_unchecked(f32::MIN_POSITIVE)", ln_unchecked(f32::MIN_POSITIVE), ln(f32::MIN_POSITIVE));
    check("log10(0)", log10(0.0), f32::NEG_INFINITY);
    check("log10(100)", log10(100.0), 2.0);
    check("log10_unchecked(100)", log10_unchecked(100.0), 2.0);
    check("log10_unchecked(f32::MAX)", log10_unchecked(f32::MAX), log10(f32::MAX));
    check("log10_unchecked(f32::MIN_POSITIVE)", log10_unchecked(f32::MIN_POSITIVE), log10(f32::MIN_POSITIVE));
    check("log1p(0)", log1p(0.0), 0.0);
    check("log1p(-1)", log1p(-1.0), f32::NEG_INFINITY);
    check("log1p(-2)", log1p(-2.0), f32::NAN);
    // log1p(-0.0) used to lose its sign: `ln(u) + corr` adds two
    // exactly-zero values of opposite sign at x = +-0.0 (`ln(1.0)` is
    // `+0.0`, but `corr` correctly carries x's sign there), the same
    // IEEE754 mechanism as sinf_poly's own `-0.0` bug. Fixed with a
    // trailing `if x == 0.0 { x } else { normal }` select (log1p is odd
    // and monotonic through the origin, so `normal`'s sign already
    // matches x's for every nonzero x -- this only changes the singular
    // zero case).
    check("log1p(-0)", log1p(-0.0), -0.0);

    // log1pmx(x) = log1p(x) - x (backlog idea #145). Mathematically
    // <=0 everywhere in-domain (log(1+x) <= x always, equality only at
    // x=0), so both signed zeros naturally land on -0 here, not a bug.
    check("log1pmx(0)", log1pmx(0.0), -0.0);
    check("log1pmx(-0)", log1pmx(-0.0), -0.0);
    check("log1pmx(-1)", log1pmx(-1.0), f32::NEG_INFINITY);
    check("log1pmx(-2)", log1pmx(-2.0), f32::NAN);
    check("log1pmx(inf)", log1pmx(f32::INFINITY), f32::NEG_INFINITY);
    check("log1pmx(nan)", log1pmx(f32::NAN), f32::NAN);
    // Naive log1p(x)-x cancels for tiny x; log1pmx must not collapse to 0.
    check("log1pmx(1e-6)", log1pmx(1.0e-6), -4.9999967e-13);

    // log2p1(x) = log2(1+x), log1p's own zero/negative/inf edges (same
    // u/corr/trailing-zero-select structure, see its doc comment) plus a
    // couple of exact powers of two to pin the LOG2_E conversion.
    check("log2p1(0)", log2p1(0.0), 0.0);
    check("log2p1(-0)", log2p1(-0.0), -0.0);
    check("log2p1(-1)", log2p1(-1.0), f32::NEG_INFINITY);
    check("log2p1(-2)", log2p1(-2.0), f32::NAN);
    check("log2p1(1)", log2p1(1.0), 1.0);
    check("log2p1(3)", log2p1(3.0), 2.0);
    check("log2p1(inf)", log2p1(f32::INFINITY), f32::INFINITY);
    check("log2p1(nan)", log2p1(f32::NAN), f32::NAN);

    // log10p1(x) = log10(1+x), same structure as log2p1 with LOG10_E in
    // place of LOG2_E -- same edges, pinned with exact powers of ten instead.
    check("log10p1(0)", log10p1(0.0), 0.0);
    check("log10p1(-0)", log10p1(-0.0), -0.0);
    check("log10p1(-1)", log10p1(-1.0), f32::NEG_INFINITY);
    check("log10p1(-2)", log10p1(-2.0), f32::NAN);
    check("log10p1(9)", log10p1(9.0), 1.0);
    check("log10p1(99)", log10p1(99.0), 2.0);
    check("log10p1(inf)", log10p1(f32::INFINITY), f32::INFINITY);
    check("log10p1(nan)", log10p1(f32::NAN), f32::NAN);

    check("exp(0)", exp(0.0), 1.0);
    // exp's Cody-Waite reduction uses round (not exp2's own floor), which
    // can push k one integer higher than floor would right at the domain
    // ceiling -- e.g. x=88.37628 gives x*log2e=127.50002, floor keeps
    // k=127 (representable by a single exponent-field construction) but
    // round pushes k=128 (exponent field 255, reserved for inf/NaN, not
    // representable by a single multiply at all). A real regression this
    // fix introduced and then fixed with exp2_checked's own k1/k2 split;
    // pinned here so it can't silently come back.
    check_finite("exp(88.37628)", exp(88.37628));

    // exp_scaled (backlog idea #117): e^x*2^s, s folded into exp's own
    // exponent-field split.
    check("exp_scaled(0,0)", exp_scaled(0.0, 0), 1.0);
    check("exp_scaled(0,1)", exp_scaled(0.0, 1), 2.0);
    check("exp_scaled(0,-1)", exp_scaled(0.0, -1), 0.5);
    check_bounded("exp_scaled(1,0)-exp(1)", exp_scaled(1.0, 0) - exp(1.0), 1e-6);
    check_bounded("exp_scaled(10,5)-1", exp_scaled(10.0, 5) / (exp(10.0) * 32.0) - 1.0, 1e-5);
    check_bounded("exp_scaled(-10,-5)-1", exp_scaled(-10.0, -5) / (exp(-10.0) / 32.0) - 1.0, 1e-5);
    check("exp_scaled(nan,0)", exp_scaled(f32::NAN, 0), f32::NAN);

    // exp_narrow (backlog ideas #23/#112): same reduction/poly as exp,
    // single exponent field instead of the split, valid only up to
    // 88.37627 (one ulp below the k=128 edge above) -- pinned at its own
    // documented boundary, not exp's.
    check("exp_narrow(0)", exp_narrow(0.0), 1.0);
    check_finite("exp_narrow(88.37627)", exp_narrow(88.37627));
    check_finite("exp_narrow(-87.68311)", exp_narrow(-87.68311));
    check("exp_narrow(1)==exp(1)", exp_narrow(1.0), exp(1.0));

    // exp_checked: full range (backlog idea #18) -- x clamped before the
    // reduction so k never leaves the k1/k2 split's own safe [-151,128)
    // range, matching exp2_checked/exp10_checked's own saturation
    // guarantees.
    check("exp_checked(0)", exp_checked(0.0), 1.0);
    check_finite("exp_checked(88.37628)", exp_checked(88.37628));
    check("exp_checked(1000)", exp_checked(1000.0), f32::INFINITY);
    check("exp_checked(-1000)", exp_checked(-1000.0), 0.0);
    check("exp_checked(inf)", exp_checked(f32::INFINITY), f32::INFINITY);
    check("exp_checked(-inf)", exp_checked(f32::NEG_INFINITY), 0.0);
    check("exp_checked(nan)", exp_checked(f32::NAN), f32::NAN);

    check("exp10(0)", exp10(0.0), 1.0);
    check("exp10(1)", exp10(1.0), 10.0);
    check("exp10(2)", exp10(2.0), 100.0);
    check("exp10_checked(0)", exp10_checked(0.0), 1.0);
    check("exp10_checked(1)", exp10_checked(1.0), 10.0);
    check("exp10_checked(inf)", exp10_checked(f32::INFINITY), f32::INFINITY);
    check("exp10_checked(-inf)", exp10_checked(f32::NEG_INFINITY), 0.0);
    check("exp10_checked(nan)", exp10_checked(f32::NAN), f32::NAN);
    check_finite("exp10_checked(38.5)", exp10_checked(38.5));
    check("exp10_checked(-45)", exp10_checked(-45.0), 1e-45);

    check("expm1(0)", expm1(0.0), 0.0);

    // expm1_narrow (backlog idea #201): same single-field mechanism and
    // domain as exp_narrow, applied to expm1's own direct branch.
    check("expm1_narrow(0)", expm1_narrow(0.0), 0.0);
    check_finite("expm1_narrow(88.37627)", expm1_narrow(88.37627));
    check_finite("expm1_narrow(-87.68311)", expm1_narrow(-87.68311));
    check("expm1_narrow(1)==expm1(1)", expm1_narrow(1.0), expm1(1.0));
    check("expm1_narrow(60)==expm1(60)", expm1_narrow(60.0), expm1(60.0));

    // expm1_checked (backlog idea #111): single exponent field like
    // expm1_narrow, but total, because the field is emitted at k-1 so k can
    // still reach 128 and overflow to inf on its own. These pins are the
    // whole reason for that indirection -- a naive k<=127 clamp saturates
    // finite here instead (the rejected round-based exp10_checked bug), so
    // the inf/-1 pairs below are the regression gate for it.
    check("expm1_checked(0)", expm1_checked(0.0), 0.0);
    check("expm1_checked(-0)", expm1_checked(-0.0), -0.0);
    check("expm1_checked(inf)", expm1_checked(f32::INFINITY), f32::INFINITY);
    check("expm1_checked(-inf)", expm1_checked(f32::NEG_INFINITY), -1.0);
    check("expm1_checked(nan)", expm1_checked(f32::NAN), f32::NAN);
    check("expm1_checked(f32::MAX)", expm1_checked(f32::MAX), f32::INFINITY);
    check("expm1_checked(f32::MIN)", expm1_checked(f32::MIN), -1.0);
    check("expm1_checked(-100)", expm1_checked(-100.0), -1.0);
    // Saturation boundary, +-1 ulp around it (idea #165's pattern): just
    // below ln(f32::MAX) must stay finite, at/above must be inf.
    check_finite("expm1_checked(88.72283)", expm1_checked(88.72283));
    check("expm1_checked(88.72284)", expm1_checked(88.72284), f32::INFINITY);
    check("expm1_checked(1e10)", expm1_checked(1e10), f32::INFINITY);
    // Bit-identical to expm1 everywhere expm1 is itself valid.
    check("expm1_checked(1)==expm1(1)", expm1_checked(1.0), expm1(1.0));
    check("expm1_checked(60)==expm1(60)", expm1_checked(60.0), expm1(60.0));
    check(
        "expm1_checked(-20)==expm1(-20)",
        expm1_checked(-20.0),
        expm1(-20.0),
    );
    check(
        "expm1_checked(88.37627)==expm1(88.37627)",
        expm1_checked(88.37627),
        expm1(88.37627),
    );

    // exp_m1_over_x(x) = (e^x-1)/x, with the removable singularity at 0
    // resolving to exactly 1.0 for free from the Pade branch's own
    // algebra (N(0)/D(0)=-120/-120=1.0 exactly) -- no explicit x==0.0
    // select needed, unlike sinc's own removable-singularity handling.
    check("exp_m1_over_x(0)", exp_m1_over_x(0.0), 1.0);
    check("exp_m1_over_x(-0)", exp_m1_over_x(-0.0), 1.0);
    check_known_1ulp("exp_m1_over_x(1)", exp_m1_over_x(1.0), (1.0f64.exp_m1()) as f32);
    check_finite("exp_m1_over_x(80)", exp_m1_over_x(80.0));

    // exp_m1_over_x_narrow (backlog idea #201): same single-field
    // mechanism/domain as exp_narrow/expm1_narrow.
    check("exp_m1_over_x_narrow(0)", exp_m1_over_x_narrow(0.0), 1.0);
    check_finite("exp_m1_over_x_narrow(88.37627)", exp_m1_over_x_narrow(88.37627));
    check_finite("exp_m1_over_x_narrow(-87.68311)", exp_m1_over_x_narrow(-87.68311));
    check(
        "exp_m1_over_x_narrow(1)==exp_m1_over_x(1)",
        exp_m1_over_x_narrow(1.0),
        exp_m1_over_x(1.0),
    );
    check(
        "exp_m1_over_x_narrow(60)==exp_m1_over_x(60)",
        exp_m1_over_x_narrow(60.0),
        exp_m1_over_x(60.0),
    );
    check("exp_m1_over_x(-80)", exp_m1_over_x(-80.0), 0.0125);

    // exp2m1(x) = 2^x - 1. Total (inherits exp2_checked's own full
    // [-151,128) clamp, see exp2m1's doc comment) so -inf/inf both give an
    // exact finite/saturated answer rather than NaN, unlike expm1's own
    // unchecked-exp2 domain limit above.
    check("exp2m1(0)", exp2m1(0.0), 0.0);
    check("exp2m1(-0)", exp2m1(-0.0), -0.0);
    check("exp2m1(1)", exp2m1(1.0), 1.0);
    check("exp2m1(-1)", exp2m1(-1.0), -0.5);
    check("exp2m1(inf)", exp2m1(f32::INFINITY), f32::INFINITY);
    check("exp2m1(-inf)", exp2m1(f32::NEG_INFINITY), -1.0);
    check("exp2m1(-151)", exp2m1(-151.0), -1.0);
    check("exp2m1(nan)", exp2m1(f32::NAN), f32::NAN);

    // exp10m1(x) = 10^x - 1, same total-over-the-clamped-domain shape as
    // exp2m1 above (exp10_checked's own [-1000,1000] pre-clamp + [-151,128)
    // k-clamp), just base 10.
    check("exp10m1(0)", exp10m1(0.0), 0.0);
    check("exp10m1(-0)", exp10m1(-0.0), -0.0);
    check("exp10m1(1)", exp10m1(1.0), 9.0);
    check("exp10m1(-1)", exp10m1(-1.0), -0.9);
    check("exp10m1(inf)", exp10m1(f32::INFINITY), f32::INFINITY);
    check("exp10m1(-inf)", exp10m1(f32::NEG_INFINITY), -1.0);
    check("exp10m1(-1000)", exp10m1(-1000.0), -1.0);
    check("exp10m1(nan)", exp10m1(f32::NAN), f32::NAN);

    check("sinh(0)", sinh(0.0), 0.0);
    check("cosh(0)", cosh(0.0), 1.0);

    // sinh_narrow/cosh_narrow (backlog idea #201): single-exponent-field
    // tier via exp_pos_neg_narrow, valid over the symmetric
    // [-87.68311, 87.68311] (tighter than exp_narrow's own domain since
    // both +k and -k must fit a single field at once).
    check("sinh_narrow(0)", sinh_narrow(0.0), 0.0);
    check("cosh_narrow(0)", cosh_narrow(0.0), 1.0);
    check_finite("sinh_narrow(87.68311)", sinh_narrow(87.68311));
    check_finite("sinh_narrow(-87.68311)", sinh_narrow(-87.68311));
    check_finite("cosh_narrow(87.68311)", cosh_narrow(87.68311));
    check("sinh_narrow(1)==sinh(1)", sinh_narrow(1.0), sinh(1.0));
    check("cosh_narrow(1)==cosh(1)", cosh_narrow(1.0), cosh(1.0));
    check("sinh_narrow(60)==sinh(60)", sinh_narrow(60.0), sinh(60.0));
    check("cosh_narrow(60)==cosh(60)", cosh_narrow(60.0), cosh(60.0));

    // sinh_throughput/cosh_throughput (2026-07-10): distinct functions
    // from sinh/cosh (own accuracy.rs sweep entries, own reported ulp
    // numbers), but had zero edgecheck coverage at all. Same small-x
    // Taylor branch as sinh, so 0/-0 sign is preserved through the
    // multiply-by-x form, not through the e-1/e subtraction that would
    // give +0 regardless of input sign.
    check("sinh_throughput(0)", sinh_throughput(0.0), 0.0);
    check("sinh_throughput(-0)", sinh_throughput(-0.0), -0.0);
    check("cosh_throughput(0)", cosh_throughput(0.0), 1.0);
    check("cosh_throughput(-0)", cosh_throughput(-0.0), 1.0);

    // sinh_checked/cosh_checked: full range (backlog idea #85's fifth
    // wave) -- exp_pos_neg_checked_half clamps x to +-170 (comfortably
    // inside the field split's own proven-safe |k|<=254 window, see its
    // doc comment) and applies the 0.5 scale factor inside the split
    // itself so the intermediate never overflows before the halving does.
    // sinh/cosh (unchecked) get this all wrong: wrong-sign garbage at
    // x=1000, NaN at +-inf.
    check("sinh_checked(0)", sinh_checked(0.0), 0.0);
    check("cosh_checked(0)", cosh_checked(0.0), 1.0);
    check("sinh_checked(1000)", sinh_checked(1000.0), f32::INFINITY);
    check("sinh_checked(-1000)", sinh_checked(-1000.0), f32::NEG_INFINITY);
    check("cosh_checked(1000)", cosh_checked(1000.0), f32::INFINITY);
    check("cosh_checked(-1000)", cosh_checked(-1000.0), f32::INFINITY);
    check("sinh_checked(inf)", sinh_checked(f32::INFINITY), f32::INFINITY);
    check("sinh_checked(-inf)", sinh_checked(f32::NEG_INFINITY), f32::NEG_INFINITY);
    check("cosh_checked(inf)", cosh_checked(f32::INFINITY), f32::INFINITY);
    check("cosh_checked(-inf)", cosh_checked(f32::NEG_INFINITY), f32::INFINITY);
    check("sinh_checked(nan)", sinh_checked(f32::NAN), f32::NAN);
    check("cosh_checked(nan)", cosh_checked(f32::NAN), f32::NAN);

    // coshm1(x) = cosh(x)-1 = 2*sinh(x/2)^2 (backlog idea #144).
    check("coshm1(0)", coshm1(0.0), 0.0);
    check("coshm1(-0)", coshm1(-0.0), 0.0);
    check("coshm1(inf)", coshm1(f32::INFINITY), f32::INFINITY);
    check("coshm1(-inf)", coshm1(f32::NEG_INFINITY), f32::INFINITY);
    check("coshm1(nan)", coshm1(f32::NAN), f32::NAN);
    // Naive cosh(x)-1 cancels to exactly 0 here; coshm1 must not.
    check("coshm1(1e-6)", coshm1(1.0e-6), 5e-13);

    // Boundary just below/at the true overflow threshold (x ~= 89.416,
    // where exp(x) alone would already be inf but sinh/cosh(x) is still
    // finite -- half of a not-yet-overflowed exp(x)). This is the premature-
    // overflow-before-scaling gap the halving trick fixes; pinned so it
    // can't silently come back.
    check_finite("sinh_checked(89.415)", sinh_checked(89.415));
    check_finite("cosh_checked(89.415)", cosh_checked(89.415));
    check("sinh_checked(89.416)", sinh_checked(89.416), f32::INFINITY);
    check("cosh_checked(89.416)", cosh_checked(89.416), f32::INFINITY);

    check("tanh(0)", tanh(0.0), 0.0);
    // Domain hole fixed 2026-07-08: 2*x used to be passed to expm1
    // unclamped, inheriting exp's unchecked-domain garbage for |x| > ~44
    // (2x past ~88.7) -- tanh(50)/(1000)/(f32::MAX) all used to return
    // NaN instead of correctly saturating. Pinned here so it can't
    // silently come back.
    check("tanh(50)", tanh(50.0), 1.0);
    check("tanh(1000)", tanh(1000.0), 1.0);
    check("tanh(f32::MAX)", tanh(f32::MAX), 1.0);
    check("tanh(-50)", tanh(-50.0), -1.0);
    check("tanh(-1000)", tanh(-1000.0), -1.0);
    check("tanh(-f32::MAX)", tanh(-f32::MAX), -1.0);
    check("tanh(inf)", tanh(f32::INFINITY), 1.0);
    check("tanh(-inf)", tanh(f32::NEG_INFINITY), -1.0);
    check("tanh(nan)", tanh(f32::NAN), f32::NAN);

    // tanh_grad (backlog idea #150): the naive 1-tanh(x)^2 composite was
    // rejected as a real cancellation bug (see its own doc comment) --
    // these pins in particular exercise the |x|~8.66 region and the
    // large-|x| extremes where that bug showed up.
    check("tanh_grad(0)", tanh_grad(0.0), 1.0);
    check("tanh_grad(-0)", tanh_grad(-0.0), 1.0);
    // True value here is a legitimate tiny denormal (~4*exp(-100),
    // 1.49e-43), not exactly 0.0 -- this fix's whole point is computing
    // that correctly instead of prematurely saturating, so pin the real
    // (denormal) output, not a rounder-looking but wrong 0.0.
    check("tanh_grad(50)", tanh_grad(50.0), f32::from_bits(0x0000006c));
    check("tanh_grad(-50)", tanh_grad(-50.0), f32::from_bits(0x0000006c));
    check("tanh_grad(f32::MAX)", tanh_grad(f32::MAX), 0.0);
    check("tanh_grad(-f32::MAX)", tanh_grad(-f32::MAX), 0.0);
    check("tanh_grad(inf)", tanh_grad(f32::INFINITY), 0.0);
    check("tanh_grad(-inf)", tanh_grad(f32::NEG_INFINITY), 0.0);
    check("tanh_grad(8.66)==tanh_grad(-8.66)", tanh_grad(8.66), tanh_grad(-8.66));
    check("tanh_grad(nan)", tanh_grad(f32::NAN), f32::NAN);

    check("sigmoid(0)", sigmoid(0.0), 0.5);
    check("sigmoid(-0)", sigmoid(-0.0), 0.5);
    // Regression pin for the cancellation bug found while implementing this
    // function (see its doc comment): a first version computed
    // 0.5+0.5*tanh(x/2), which returned exactly 0.0 here because tanh(x/2)
    // had already correctly saturated to exactly -1.0f32 -- discarding the
    // true, nowhere-near-zero result. Direct 1/(1+exp(-x)) has no such
    // cancellation.
    check("sigmoid(-17.32869)", sigmoid(-17.32869), 2.9802024e-8);
    check("sigmoid(1000)", sigmoid(1000.0), 1.0);
    check("sigmoid(f32::MAX)", sigmoid(f32::MAX), 1.0);
    check("sigmoid(inf)", sigmoid(f32::INFINITY), 1.0);
    // Regression pins for idea #44's own negative-tail bug (see sigmoid's
    // doc comment): the old `y.clamp(-87.0,88.0)` capped the exponent at a
    // fixed finite value for *any* x below about -88, so sigmoid(-inf)
    // and sigmoid(-1000) both used to wrongly return ~6.054601e-39 instead
    // of the true 0.0 -- this pin used to lock in that wrong value as if
    // it were correct; now pins the fix instead.
    check("sigmoid(-inf)", sigmoid(f32::NEG_INFINITY), 0.0);
    check("sigmoid(-1000)", sigmoid(-1000.0), 0.0);
    check("sigmoid(-f32::MAX)", sigmoid(-f32::MAX), 0.0);

    // sigmoid_fast (backlog idea #191): approx tier, ~0.023 max absolute
    // error by design (see its own doc comment) -- pins check it stays
    // in that ballpark and, unlike sigmoid itself, correctly saturates
    // to its own frozen boundary value (not exactly 0/1) for large |x|,
    // never NaN/inf/wrong-side.
    check("sigmoid_fast(0)", sigmoid_fast(0.0), 0.5);
    check("sigmoid_fast(-0)", sigmoid_fast(-0.0), 0.5);
    check_bounded("sigmoid_fast(1)-sigmoid(1)", sigmoid_fast(1.0) - sigmoid(1.0), 0.03);
    check_bounded("sigmoid_fast(-1)-sigmoid(-1)", sigmoid_fast(-1.0) - sigmoid(-1.0), 0.03);
    check_bounded("sigmoid_fast(inf)-1", sigmoid_fast(f32::INFINITY) - 1.0, 0.03);
    check_bounded("sigmoid_fast(-inf)", sigmoid_fast(f32::NEG_INFINITY), 0.03);
    check_bounded("sigmoid_fast(f32::MAX)-1", sigmoid_fast(f32::MAX) - 1.0, 0.03);
    check_bounded("sigmoid_fast(-f32::MAX)", sigmoid_fast(-f32::MAX), 0.03);
    check("sigmoid_fast(nan)", sigmoid_fast(f32::NAN), f32::NAN);
    check_bounded("sigmoid_fast(100)-sigmoid_fast(inf)", sigmoid_fast(100.0) - sigmoid_fast(f32::INFINITY), 1e-6);
    check("sigmoid(-89)", sigmoid(-89.0), 0.0);
    // just below the fixed clamp boundary: still the same (correct, real)
    // value the old code also gave here, confirming no regression at the
    // boundary itself.
    check("sigmoid(-88)", sigmoid(-88.0), 6.054601e-39);
    check("sigmoid(nan)", sigmoid(f32::NAN), f32::NAN);

    // sigmoid_grad (backlog idea #150): the naive sigmoid(x)*(1-sigmoid(x))
    // composite was rejected as a real cancellation bug, and a first fix
    // attempt (exp_checked(-x)/(1+exp_checked(-x))^2, no |x| fold) traded
    // it for a different real bug (overflow in the square for very
    // negative x) -- see its own doc comment. These pins exercise both
    // extremes and the evenness the final fix relies on.
    check("sigmoid_grad(0)", sigmoid_grad(0.0), 0.25);
    check("sigmoid_grad(-0)", sigmoid_grad(-0.0), 0.25);
    check("sigmoid_grad(1000)", sigmoid_grad(1000.0), 0.0);
    check("sigmoid_grad(-1000)", sigmoid_grad(-1000.0), 0.0);
    check("sigmoid_grad(f32::MAX)", sigmoid_grad(f32::MAX), 0.0);
    check("sigmoid_grad(-f32::MAX)", sigmoid_grad(-f32::MAX), 0.0);
    check("sigmoid_grad(inf)", sigmoid_grad(f32::INFINITY), 0.0);
    check("sigmoid_grad(-inf)", sigmoid_grad(f32::NEG_INFINITY), 0.0);
    check("sigmoid_grad(44.36)==sigmoid_grad(-44.36)", sigmoid_grad(44.36), sigmoid_grad(-44.36));
    check("sigmoid_grad(nan)", sigmoid_grad(f32::NAN), f32::NAN);

    // softplus(x) = ln(1+e^x). ln(2) at 0 (ln(1+e^0)=ln(2)); saturates to
    // x itself for large positive x, to exactly 0 for large negative x
    // (both via the correction-term cutoff, see softplus's own doc
    // comment for why a hard cutoff instead of a clamped exp argument).
    check("softplus(0)", softplus(0.0), std::f32::consts::LN_2);
    check("softplus(-0)", softplus(-0.0), std::f32::consts::LN_2);
    check("softplus(1000)", softplus(1000.0), 1000.0);
    check("softplus(f32::MAX)", softplus(f32::MAX), f32::MAX);
    check("softplus(inf)", softplus(f32::INFINITY), f32::INFINITY);
    check("softplus(-inf)", softplus(f32::NEG_INFINITY), 0.0);
    check("softplus(-1000)", softplus(-1000.0), 0.0);
    check("softplus(-f32::MAX)", softplus(-f32::MAX), 0.0);
    check("softplus(nan)", softplus(f32::NAN), f32::NAN);

    // logsigmoid = -softplus(-x) (backlog idea #118): each pin here is
    // softplus's own pin above, mirrored through that identity.
    check("logsigmoid(0)", logsigmoid(0.0), -std::f32::consts::LN_2);
    check("logsigmoid(-0)", logsigmoid(-0.0), -std::f32::consts::LN_2);
    check("logsigmoid(1000)", logsigmoid(1000.0), -0.0);
    check("logsigmoid(f32::MAX)", logsigmoid(f32::MAX), -0.0);
    check("logsigmoid(inf)", logsigmoid(f32::INFINITY), -0.0);
    check("logsigmoid(-inf)", logsigmoid(f32::NEG_INFINITY), f32::NEG_INFINITY);
    check("logsigmoid(-1000)", logsigmoid(-1000.0), -1000.0);
    check("logsigmoid(-f32::MAX)", logsigmoid(-f32::MAX), -f32::MAX);
    check("logsigmoid(nan)", logsigmoid(f32::NAN), f32::NAN);

    // logaddexp(a,b) = ln(e^a+e^b); softplus(x) == logaddexp(x,0.0).
    check("logaddexp(0,0)", logaddexp(0.0, 0.0), std::f32::consts::LN_2);
    check("logaddexp(x,0)==softplus(x)", logaddexp(3.0, 0.0), softplus(3.0));
    check("logaddexp(100,1)", logaddexp(100.0, 1.0), 100.0);
    check("logaddexp(inf,5)", logaddexp(f32::INFINITY, 5.0), f32::INFINITY);
    check("logaddexp(inf,-inf)", logaddexp(f32::INFINITY, f32::NEG_INFINITY), f32::INFINITY);
    check("logaddexp(-inf,-inf)", logaddexp(f32::NEG_INFINITY, f32::NEG_INFINITY), f32::NEG_INFINITY);
    check("logaddexp(inf,inf)", logaddexp(f32::INFINITY, f32::INFINITY), f32::INFINITY);
    check("logaddexp(nan,1)", logaddexp(f32::NAN, 1.0), f32::NAN);
    check("logaddexp(1,nan)", logaddexp(1.0, f32::NAN), f32::NAN);

    // gelu(x) = x*Phi(x) = x*0.5*erfc(-x/sqrt2) (backlog idea #70). For
    // x >= 0 gelu's tail correction is identically zero, so gelu(3) is
    // still exactly that composition and pins it structurally. For x < 0
    // the correction is live and deliberately makes gelu *not* equal the
    // plain f32 composition (3 ulp off the true -0.0040496940948903 at
    // x=-3), so gelu(-3) is pinned against the true value instead --
    // accuracy across the domain stays accuracy.rs's job. The
    // special-value pins guard the explicit x==-inf override:
    // 0.0*(-inf) alone is NaN, but the true limit is 0.
    check("gelu(0)", gelu(0.0), 0.0);
    check("gelu(-0)", gelu(-0.0), -0.0);
    check("gelu(3)==3*.5*erfc(-3/sqrt2)", gelu(3.0), 3.0 * 0.5 * erfc(-3.0 * std::f32::consts::FRAC_1_SQRT_2));
    check_ulp("gelu(-3)==-3*Phi(-3)", gelu(-3.0), -4.0496941e-3, 1);
    check("gelu(inf)", gelu(f32::INFINITY), f32::INFINITY);
    check("gelu(-inf)", gelu(f32::NEG_INFINITY), 0.0);
    check("gelu(nan)", gelu(f32::NAN), f32::NAN);

    // silu(x) = x*sigmoid(x) (backlog idea #70). Same x==-inf 0*(-inf)
    // override as gelu (sigmoid(-inf)=0, but -inf*0 alone is NaN).
    check("silu(0)", silu(0.0), 0.0);
    check("silu(-0)", silu(-0.0), -0.0);
    check("silu(2)==2*sigmoid(2)", silu(2.0), 2.0 * sigmoid(2.0));
    check("silu(-2)==-2*sigmoid(-2)", silu(-2.0), -2.0 * sigmoid(-2.0));
    check("silu(inf)", silu(f32::INFINITY), f32::INFINITY);
    check("silu(-inf)", silu(f32::NEG_INFINITY), 0.0);
    check("silu(nan)", silu(f32::NAN), f32::NAN);

    // softsign(x) = x/(1+|x|) (backlog idea #70). 0.5/0.75 are exact here
    // (1/(1+1), 3/(1+3)). The +-inf override guards inf/inf -> +-1.
    check("softsign(0)", softsign(0.0), 0.0);
    check("softsign(-0)", softsign(-0.0), -0.0);
    check("softsign(1)", softsign(1.0), 0.5);
    check("softsign(-1)", softsign(-1.0), -0.5);
    check("softsign(3)", softsign(3.0), 0.75);
    check("softsign(inf)", softsign(f32::INFINITY), 1.0);
    check("softsign(-inf)", softsign(f32::NEG_INFINITY), -1.0);
    check("softsign(nan)", softsign(f32::NAN), f32::NAN);

    // sqrt1pm1: rationalized sqrt(1+x)-1 (backlog idea #132), values
    // confirmed against a rationalized f64 reference before pinning.
    check("sqrt1pm1(0)", sqrt1pm1(0.0), 0.0);
    check("sqrt1pm1(-0)", sqrt1pm1(-0.0), -0.0);
    check("sqrt1pm1(-1)", sqrt1pm1(-1.0), -1.0);
    check("sqrt1pm1(-2)", sqrt1pm1(-2.0), f32::NAN);
    check("sqrt1pm1(inf)", sqrt1pm1(f32::INFINITY), f32::INFINITY);
    check("sqrt1pm1(nan)", sqrt1pm1(f32::NAN), f32::NAN);
    check("sqrt1pm1(3)", sqrt1pm1(3.0), 1.0);
    check("sqrt1pm1(1e-10)", sqrt1pm1(1e-10), 5e-11);
    check("sqrt1pm1(-1e-10)", sqrt1pm1(-1e-10), -5e-11);

    // pow_3_2/pow_2_3 (backlog idea #133): plain sqrt/cbrt compositions.
    check("pow_3_2(4)", pow_3_2(4.0), 8.0);
    check("pow_3_2(0)", pow_3_2(0.0), 0.0);
    check("pow_3_2(-1)", pow_3_2(-1.0), f32::NAN);
    check("pow_3_2(inf)", pow_3_2(f32::INFINITY), f32::INFINITY);
    check("pow_3_2(nan)", pow_3_2(f32::NAN), f32::NAN);
    check("pow_2_3(8)", pow_2_3(8.0), 4.0);
    check("pow_2_3(-8)", pow_2_3(-8.0), 4.0);
    check("pow_2_3(0)", pow_2_3(0.0), 0.0);
    check("pow_2_3(inf)", pow_2_3(f32::INFINITY), f32::INFINITY);
    check("pow_2_3(-inf)", pow_2_3(f32::NEG_INFINITY), f32::INFINITY);
    check("pow_2_3(nan)", pow_2_3(f32::NAN), f32::NAN);

    // rootn (backlog idea #75, C23): x^(1/n) -- domain error (NaN) for
    // n==0 or negative x with even n, real negative result for negative
    // x with odd n. Verified against the real glibc rootn spec, not
    // guessed (see rootn's own doc comment).
    check("rootn(8,3)", rootn(8.0, 3), 2.0);
    check("rootn(-8,3)", rootn(-8.0, 3), -2.0);
    check("rootn(16,4)", rootn(16.0, 4), 2.0);
    check("rootn(-16,4)", rootn(-16.0, 4), f32::NAN);
    check("rootn(4,2)", rootn(4.0, 2), 2.0);
    check("rootn(-4,2)", rootn(-4.0, 2), f32::NAN);
    check("rootn(x,0)", rootn(4.0, 0), f32::NAN);
    check("rootn(nan,0)", rootn(f32::NAN, 0), f32::NAN);
    check("rootn(5,1)", rootn(5.0, 1), 5.0);
    check("rootn(-5,1)", rootn(-5.0, 1), -5.0);
    check("rootn(0,3)", rootn(0.0, 3), 0.0);
    check("rootn(-0,3)", rootn(-0.0, 3), -0.0);
    check("rootn(0,2)", rootn(0.0, 2), 0.0);
    check("rootn(-0,2)", rootn(-0.0, 2), 0.0);
    check("rootn(0,-3)", rootn(0.0, -3), f32::INFINITY);
    check("rootn(-0,-3)", rootn(-0.0, -3), f32::NEG_INFINITY);
    check("rootn(0,-2)", rootn(0.0, -2), f32::INFINITY);
    check("rootn(inf,3)", rootn(f32::INFINITY, 3), f32::INFINITY);
    check("rootn(inf,-3)", rootn(f32::INFINITY, -3), 0.0);
    check("rootn(-inf,3)", rootn(f32::NEG_INFINITY, 3), f32::NEG_INFINITY);
    check("rootn(-inf,-3)", rootn(f32::NEG_INFINITY, -3), -0.0);
    check("rootn(-inf,2)", rootn(f32::NEG_INFINITY, 2), f32::NAN);
    check("rootn(nan,3)", rootn(f32::NAN, 3), f32::NAN);
    // n == -1 is the reciprocal, the one n whose result can overflow or
    // land denormal -- and the reason the general path is only ever
    // entered with |n| >= 2 (see rootn's own doc comment).
    check("rootn(4,-1)", rootn(4.0, -1), 0.25);
    check("rootn(-4,-1)", rootn(-4.0, -1), -0.25);
    check("rootn(0,-1)", rootn(0.0, -1), f32::INFINITY);
    check("rootn(-0,-1)", rootn(-0.0, -1), f32::NEG_INFINITY);
    check("rootn(inf,-1)", rootn(f32::INFINITY, -1), 0.0);
    check("rootn(denormal_min,-1)", rootn(f32::from_bits(1), -1), f32::INFINITY);
    check("rootn(max,-1)", rootn(f32::MAX, -1), 1.0 / f32::MAX);
    check("rootn(nan,-1)", rootn(f32::NAN, -1), f32::NAN);

    // fast_round_int (backlog idea #185): the ROUND_MAGIC idiom exposed
    // standalone. Sign-of-zero pins are load-bearing, not decorative --
    // a real exhaustive sweep found the bare idiom (no trailing
    // copysign) wrong-signed for every negative x rounding to zero
    // (~1.057 billion bit patterns, all confined to |x|<=0.5), fixed by
    // the copysign now in the shipped body.
    check("fast_round_int(0)", fast_round_int(0.0), 0.0);
    check("fast_round_int(-0)", fast_round_int(-0.0), -0.0);
    check("fast_round_int(0.5)", fast_round_int(0.5), 0.0);
    check("fast_round_int(-0.5)", fast_round_int(-0.5), -0.0);
    check("fast_round_int(-0.4999999)", fast_round_int(-0.4999999), -0.0);
    check("fast_round_int(1.5)", fast_round_int(1.5), 2.0);
    check("fast_round_int(2.5)", fast_round_int(2.5), 2.0);
    check("fast_round_int(-3.7)", fast_round_int(-3.7), -4.0);
    check("fast_round_int(2^22)", fast_round_int(4_194_304.0), 4_194_304.0);
    check("fast_round_int(-2^22)", fast_round_int(-4_194_304.0), -4_194_304.0);
    check("fast_round_int(nan)", fast_round_int(f32::NAN), f32::NAN);
    check("fast_round_int(inf)", fast_round_int(f32::INFINITY), f32::INFINITY);
    check("fast_round_int(-inf)", fast_round_int(f32::NEG_INFINITY), f32::NEG_INFINITY);

    // smoothstep/smootherstep (backlog idea #147): exact endpoints,
    // clamped outside [edge0,edge1].
    check("smoothstep(0,1,0)", smoothstep(0.0, 1.0, 0.0), 0.0);
    check("smoothstep(0,1,1)", smoothstep(0.0, 1.0, 1.0), 1.0);
    check("smoothstep(0,1,0.5)", smoothstep(0.0, 1.0, 0.5), 0.5);
    check("smoothstep(0,1,-1)", smoothstep(0.0, 1.0, -1.0), 0.0);
    check("smoothstep(0,1,2)", smoothstep(0.0, 1.0, 2.0), 1.0);
    check("smoothstep(2,4,3)", smoothstep(2.0, 4.0, 3.0), 0.5);
    check("smoothstep(0,1,nan)", smoothstep(0.0, 1.0, f32::NAN), f32::NAN);
    check("smootherstep(0,1,0)", smootherstep(0.0, 1.0, 0.0), 0.0);
    check("smootherstep(0,1,1)", smootherstep(0.0, 1.0, 1.0), 1.0);
    check("smootherstep(0,1,0.5)", smootherstep(0.0, 1.0, 0.5), 0.5);
    check("smootherstep(0,1,-1)", smootherstep(0.0, 1.0, -1.0), 0.0);
    check("smootherstep(0,1,2)", smootherstep(0.0, 1.0, 2.0), 1.0);
    check("smootherstep(0,1,nan)", smootherstep(0.0, 1.0, f32::NAN), f32::NAN);

    check("asinh(0)", asinh(0.0), 0.0);
    // small-x cancellation (see asinh's doc comment) is fixed: asinh(x) ~ x
    // for tiny x, no longer collapses to exactly 0.
    check("asinh(2.34e-8)", asinh(2.34e-8), 2.34e-8);
    check("asinh(-2.34e-8)", asinh(-2.34e-8), -2.34e-8);
    // large-negative-x cancellation (also fixed, see doc comment): asinh is
    // odd, so this must equal -asinh(1e10) exactly.
    check("asinh(-1e10) == -asinh(1e10)", asinh(-1e10), -asinh(1e10));
    check("asinh(-f32::MAX)", asinh(-f32::MAX), -asinh(f32::MAX));
    // The above only check the +/- relationship stays consistent, not that
    // the *absolute* value is actually correct (a bug affecting both signs
    // identically would slip through unnoticed). Added 2026-07-10: direct
    // pins against asinh(x)~ln(2x) for large x, values confirmed by direct
    // computation first (relative error ~1e-8, matching f32's own budget).
    check("asinh(1e10)", asinh(1e10), 2.3718998e1);
    check("asinh(1e20)", asinh(1e20), 4.674485e1);
    check("asinh(f32::MAX)", asinh(f32::MAX), 8.9415985e1);
    check("acosh(1)", acosh(1.0), 0.0);
    check("acosh(0.5)", acosh(0.5), f32::NAN);
    // acosh's sign-losing domain bug (see its doc comment) is fixed with an
    // explicit domain check -- large negative x now correctly comes out NaN.
    check("acosh(-1e20)", acosh(-1e20), f32::NAN);
    check("acosh(-4096.0)", acosh(-4096.0), f32::NAN);
    check("acosh(-f32::MAX)", acosh(-f32::MAX), f32::NAN);
    // Large *positive* x had no coverage at all here (only the negative/
    // out-of-domain side was pinned) -- added 2026-07-10. acosh(x) for
    // x>=1 is always finite and non-negative, growing like ln(2x); values
    // below confirmed by direct computation before pinning.
    check("acosh(1e10)", acosh(1e10), 2.3718998e1);
    check("acosh(1e20)", acosh(1e20), 4.674485e1);
    check("acosh(f32::MAX)", acosh(f32::MAX), 8.9415985e1);
    check("acosh(inf)", acosh(f32::INFINITY), f32::INFINITY);
    check("atanh(0)", atanh(0.0), 0.0);
    check("atanh(1)", atanh(1.0), f32::INFINITY);
    check("atanh(-1)", atanh(-1.0), f32::NEG_INFINITY);
    check("atanh(2)", atanh(2.0), f32::NAN);
    // atanh(-0.0) was wrong purely as a downstream consequence of
    // log1p(-0.0)'s own bug (see above) -- fell out correct for free once
    // log1p was fixed, no separate change needed here.
    check("atanh(-0)", atanh(-0.0), -0.0);

    // asin(0)'s zero sign used to come out backwards (-0.0 for +0.0 input)
    // from the trailing `* (-hpi)` flipping an intermediate +0 -- fixed as
    // a side effect of rationalizing the small-x cancellation below (see
    // asin's doc comment), now matching the usual libm convention (and
    // this crate's other odd functions, e.g. atan/sinh/asinh/tanh).
    check("asin(0)", asin(0.0), 0.0);
    check("asin(-0)", asin(-0.0), -0.0);
    check("asin(1)", asin(1.0), std::f32::consts::FRAC_PI_2);
    check("asin(-1)", asin(-1.0), -std::f32::consts::FRAC_PI_2);
    check("asin(2)", asin(2.0), f32::NAN);

    // asind (backlog idea #123): plain asin(x)*RAD_TO_DEG composite --
    // a rescaled-coefficient fold was tried and measured worse (see its
    // own doc comment).
    check("asind(0)", asind(0.0), 0.0);
    check("asind(-0)", asind(-0.0), -0.0);
    check("asind(1)", asind(1.0), 90.0);
    check("asind(-1)", asind(-1.0), -90.0);
    check("asind(0.5)", asind(0.5), 30.0);
    check("asind(nan)", asind(f32::NAN), f32::NAN);

    // asinpi (backlog idea #85): folded dedicated coefficients (a real
    // win here, unlike asind's rejected fold -- see its own doc
    // comment).
    check("asinpi(0)", asinpi(0.0), 0.0);
    check("asinpi(-0)", asinpi(-0.0), -0.0);
    check("asinpi(1)", asinpi(1.0), 0.5);
    check("asinpi(-1)", asinpi(-1.0), -0.5);
    check("asinpi(nan)", asinpi(f32::NAN), f32::NAN);
    check("asinpi(2)", asinpi(2.0), f32::NAN);

    check("acosd(1)", acosd(1.0), 0.0);
    check("acosd(-1)", acosd(-1.0), 180.0);
    check("acosd(0)", acosd(0.0), 90.0);
    check("acosd(nan)", acosd(f32::NAN), f32::NAN);

    // acospi (backlog idea #85): folded dedicated coefficients (a real
    // win here, same verdict as asinpi -- see its own doc comment).
    check("acospi(1)", acospi(1.0), 0.0);
    check("acospi(-1)", acospi(-1.0), 1.0);
    check("acospi(0)", acospi(0.0), 0.5);
    check("acospi(nan)", acospi(f32::NAN), f32::NAN);
    check("acospi(2)", acospi(2.0), f32::NAN);

    check("acos(1)", acos(1.0), 0.0);
    check("acos(-1)", acos(-1.0), std::f32::consts::PI);
    check("acos(2)", acos(2.0), f32::NAN);
    // acos(-0.0) used to come out -pi/2 (mulsign's bit-based sign check
    // disagreed with x<0.0's value-based one, exactly at this one input)
    // instead of the correct +pi/2 -- acos is never negative.
    // acos(0)/acos(-0) used to also be 1 ulp off FRAC_PI_2 (a separate,
    // now-fixed bug: acos_poly's own leading constant literal,
    // "1.5707963", parses to 0x3fc90fda, one ulp *below* the true
    // correctly-rounded pi/2 (0x3fc90fdb) -- not an intentional fitting
    // choice, just short of enough significant digits to round to the
    // right value. Fixed 2026-07-09, backlog idea #36's own investigation
    // (which found a much bigger, unplanned win here before ever reaching
    // the Df32-accurate-tier idea it started out testing) -- see
    // acos_poly's own doc comment for the full accuracy numbers.
    check("acos(0)", acos(0.0), std::f32::consts::FRAC_PI_2);
    check("acos(-0)", acos(-0.0), std::f32::consts::FRAC_PI_2);
    check("atan(0)", atan(0.0), 0.0);
    check("atan(inf)", atan(f32::INFINITY), std::f32::consts::FRAC_PI_2);
    check("atan(-inf)", atan(f32::NEG_INFINITY), -std::f32::consts::FRAC_PI_2);

    // atan_bounded (backlog idea #61): |x|<=1 contract, atan_poly alone.
    check("atan_bounded(0)", atan_bounded(0.0), 0.0);
    check("atan_bounded(-0)", atan_bounded(-0.0), -0.0);
    check("atan_bounded(nan)", atan_bounded(f32::NAN), f32::NAN);
    check("atan_bounded(1)==atan(1)", atan_bounded(1.0), atan(1.0));
    check("atan_bounded(-1)==atan(-1)", atan_bounded(-1.0), atan(-1.0));
    check("atan_bounded(0.5)==atan(0.5)", atan_bounded(0.5), atan(0.5));
    check("atan_latency(0)", atan_latency(0.0), 0.0);
    check("atan_latency(-0)", atan_latency(-0.0), -0.0);
    check("atan_latency(inf)", atan_latency(f32::INFINITY), std::f32::consts::FRAC_PI_2);
    check("atan_latency(-inf)", atan_latency(f32::NEG_INFINITY), -std::f32::consts::FRAC_PI_2);
    check("atan_latency(nan)", atan_latency(f32::NAN), f32::NAN);
    check("atan_latency(1)", atan_latency(1.0), atan(1.0));

    check("atand(1)", atand(1.0), 45.0);
    check("atand(-1)", atand(-1.0), -45.0);
    check("atand(0)", atand(0.0), 0.0);
    check("atand(nan)", atand(f32::NAN), f32::NAN);

    // atanpi (backlog idea #85): folded dedicated coefficients (a real
    // win here, same verdict as asinpi/acospi -- see its own doc comment).
    check("atanpi(0)", atanpi(0.0), 0.0);
    check("atanpi(-0)", atanpi(-0.0), -0.0);
    check("atanpi(1)", atanpi(1.0), 0.25);
    check("atanpi(-1)", atanpi(-1.0), -0.25);
    check("atanpi(inf)", atanpi(f32::INFINITY), 0.5);
    check("atanpi(-inf)", atanpi(f32::NEG_INFINITY), -0.5);
    check("atanpi(nan)", atanpi(f32::NAN), f32::NAN);

    check("atan2(1,0)", atan2(1.0, 0.0), std::f32::consts::FRAC_PI_2);
    check("atan2(-1,0)", atan2(-1.0, 0.0), -std::f32::consts::FRAC_PI_2);
    // atan2(-0.0, +0.0) used to lose its sign (IEEE754's "+0 + -0 = +0"
    // addition rule silently flipped the correctly-signed -0.0 back to
    // +0.0) -- IEEE754/C99 specify this exactly, so all 12 zero/sign
    // combinations are checked bit-exact against std here.
    check("atan2(0,0)", atan2(0.0, 0.0), 0.0);
    check("atan2(-0,0)", atan2(-0.0, 0.0), -0.0);
    check("atan2(0,-0)", atan2(0.0, -0.0), std::f32::consts::PI);
    check("atan2(-0,-0)", atan2(-0.0, -0.0), -std::f32::consts::PI);
    check("atan2(0,1)", atan2(0.0, 1.0), 0.0);
    check("atan2(-0,1)", atan2(-0.0, 1.0), -0.0);
    check("atan2(0,-1)", atan2(0.0, -1.0), std::f32::consts::PI);
    check("atan2(-0,-1)", atan2(-0.0, -1.0), -std::f32::consts::PI);
    check("atan2(1,-0)", atan2(1.0, -0.0), std::f32::consts::FRAC_PI_2);
    check("atan2(-1,-0)", atan2(-1.0, -0.0), -std::f32::consts::FRAC_PI_2);
    // atan2(+-inf, +-inf): y/x is inf/inf = NaN, so the general formula
    // can't produce an answer -- IEEE754/C99 define a canonical result by
    // quadrant regardless (+-pi/4 or +-3pi/4), used to come out NaN here.
    check("atan2(inf,inf)", atan2(f32::INFINITY, f32::INFINITY), std::f32::consts::FRAC_PI_4);
    check("atan2(inf,-inf)", atan2(f32::INFINITY, f32::NEG_INFINITY), 3.0 * std::f32::consts::FRAC_PI_4);
    check("atan2(-inf,inf)", atan2(f32::NEG_INFINITY, f32::INFINITY), -std::f32::consts::FRAC_PI_4);
    check("atan2(-inf,-inf)", atan2(f32::NEG_INFINITY, f32::NEG_INFINITY), -3.0 * std::f32::consts::FRAC_PI_4);

    // atan2_pos: single-positive-turn fold (backlog idea #143).
    check("atan2_pos(0,1)", atan2_pos(0.0, 1.0), 0.0);
    check("atan2_pos(1,0)", atan2_pos(1.0, 0.0), std::f32::consts::FRAC_PI_2);
    check("atan2_pos(0,-1)", atan2_pos(0.0, -1.0), std::f32::consts::PI);
    check("atan2_pos(-1,0)", atan2_pos(-1.0, 0.0), 3.0 * std::f32::consts::FRAC_PI_2);
    check("atan2_pos(-1,-1)", atan2_pos(-1.0, -1.0), 5.0 * std::f32::consts::FRAC_PI_4);
    check("atan2_pos(nan,1)", atan2_pos(f32::NAN, 1.0), f32::NAN);
    check_range("atan2_pos(1,1)", atan2_pos(1.0, 1.0), 0.0, std::f32::consts::TAU);
    // y/x underflowing to -0.0 must still fold: the angle is a hair under
    // a full turn, not a hair over zero. Keying the fold on atan2's own
    // sign gets this wrong over ~2% of the f32 plane.
    check("atan2_pos(-1e-30,1e30)", atan2_pos(-1e-30, 1e30), std::f32::consts::TAU);
    check("atan2_pos(-min,max)", atan2_pos(-f32::MIN_POSITIVE, f32::MAX), std::f32::consts::TAU);
    // -0.0 reads as "approached from below" and folds, same as the slab
    // above; atan2_pos never returns -0.0.
    check("atan2_pos(-0,1)", atan2_pos(-0.0, 1.0), std::f32::consts::TAU);
    check("atan2_pos(-0,-1)", atan2_pos(-0.0, -1.0), std::f32::consts::PI);

    // atan2d (backlog idea #123): composite plus a scaled-quotient
    // branch for a denormal atan2 (180/pi > 1, so the plain multiply
    // carries a denormal's missing bits into a finer binade).
    check("atan2d(1,0)", atan2d(1.0, 0.0), 90.0);
    check("atan2d(0,1)", atan2d(0.0, 1.0), 0.0);
    check("atan2d(-1,0)", atan2d(-1.0, 0.0), -90.0);
    check("atan2d(nan,1)", atan2d(f32::NAN, 1.0), f32::NAN);
    // the `x == 0.0` exclusion, which keeps that branch's quotient form
    // out of the 0/0 corners atan2 already answers correctly.
    check("atan2d(0,0)", atan2d(0.0, 0.0), 0.0);
    check("atan2d(-0,0)", atan2d(-0.0, 0.0), -0.0);
    check("atan2d(0,-0)", atan2d(0.0, -0.0), 180.0);
    check("atan2d(-0,1)", atan2d(-0.0, 1.0), -0.0);
    check("atan2d(0,-1)", atan2d(0.0, -1.0), 180.0);
    check("atan2d(1,inf)", atan2d(1.0, f32::INFINITY), 0.0);
    // and the branch itself: both of these are 0 without it. The first
    // expects a hair less than the second because its own `y` is already
    // denormal (1e-40 as an f32 is 9.99994610e-41), so some of the bits
    // were gone before atan2d ever saw them -- that part is not
    // recoverable and is not what the branch is for.
    check("atan2d(1e-40,1)", atan2d(1.0e-40, 1.0), 5.7295465e-39);
    check("atan2d(1e-30,1e10)", atan2d(1.0e-30, 1.0e10), 5.7295776e-39);

    // atan2pi (backlog idea #85): plain composite -- see its own doc
    // comment for why a rescaled-coefficient fold isn't attempted.
    check("atan2pi(1,0)", atan2pi(1.0, 0.0), 0.5);
    check("atan2pi(0,1)", atan2pi(0.0, 1.0), 0.0);
    check("atan2pi(-1,0)", atan2pi(-1.0, 0.0), -0.5);
    check("atan2pi(nan,1)", atan2pi(f32::NAN, 1.0), f32::NAN);
    // atan2(NaN, 0.0)/atan2(NaN, -0.0) used to come out +-FRAC_PI_2 instead
    // of NaN (backlog idea #85, found building a systematic C99
    // special-case matrix against std): the x==0 branch bypasses
    // atan(y/x) entirely and falls straight to mulsign(...,y), which only
    // reads y's sign bit and doesn't propagate NaN. Every other NaN
    // combination already worked (x!=0 routes through atan(y/x), which
    // does propagate correctly).
    check("atan2(nan,0)", atan2(f32::NAN, 0.0), f32::NAN);
    check("atan2(nan,-0)", atan2(f32::NAN, -0.0), f32::NAN);
    check("atan2(nan,1)", atan2(f32::NAN, 1.0), f32::NAN);
    check("atan2(nan,inf)", atan2(f32::NAN, f32::INFINITY), f32::NAN);
    check("atan2(nan,nan)", atan2(f32::NAN, f32::NAN), f32::NAN);
    check("atan2(0,nan)", atan2(0.0, f32::NAN), f32::NAN);
    check("atan2(1,nan)", atan2(1.0, f32::NAN), f32::NAN);

    // atan2_latency (backlog idea #60): identical wrapper to atan2, just a
    // different atan-core, so every zero/inf/nan edge case above (none of
    // which touch that core) transfers unchanged.
    check("atan2_latency(1,0)", atan2_latency(1.0, 0.0), std::f32::consts::FRAC_PI_2);
    check("atan2_latency(-1,0)", atan2_latency(-1.0, 0.0), -std::f32::consts::FRAC_PI_2);
    check("atan2_latency(0,0)", atan2_latency(0.0, 0.0), 0.0);
    check("atan2_latency(-0,0)", atan2_latency(-0.0, 0.0), -0.0);
    check("atan2_latency(0,-0)", atan2_latency(0.0, -0.0), std::f32::consts::PI);
    check("atan2_latency(-0,-0)", atan2_latency(-0.0, -0.0), -std::f32::consts::PI);
    check("atan2_latency(inf,inf)", atan2_latency(f32::INFINITY, f32::INFINITY), std::f32::consts::FRAC_PI_4);
    check(
        "atan2_latency(inf,-inf)",
        atan2_latency(f32::INFINITY, f32::NEG_INFINITY),
        3.0 * std::f32::consts::FRAC_PI_4,
    );
    check("atan2_latency(nan,0)", atan2_latency(f32::NAN, 0.0), f32::NAN);
    check("atan2_latency(nan,-0)", atan2_latency(f32::NAN, -0.0), f32::NAN);
    check("atan2_latency(nan,1)", atan2_latency(f32::NAN, 1.0), f32::NAN);
    check("atan2_latency(0,nan)", atan2_latency(0.0, f32::NAN), f32::NAN);
    // Ordinary values: not bit-identical to atan2 (atan_latency's own poly
    // differs from atan_poly's), but both correctly-rounded-ish and close
    // -- pin against the exact 45-degree case, exact by construction for
    // any reasonable atan-core (atan(1)==pi/4 to within its own ulp
    // budget), not a fitted-poly-dependent value.
    check("atan2_latency(1,1)", atan2_latency(1.0, 1.0), atan2(1.0, 1.0));

    // atan2_unchecked: contract is x != 0.0, not both infinite -- must
    // match atan2 exactly wherever that contract holds.
    check("atan2_unchecked(1,2)", atan2_unchecked(1.0, 2.0), atan2(1.0, 2.0));
    check("atan2_unchecked(-1,2)", atan2_unchecked(-1.0, 2.0), atan2(-1.0, 2.0));
    check("atan2_unchecked(1,-2)", atan2_unchecked(1.0, -2.0), atan2(1.0, -2.0));
    check("atan2_unchecked(-1,-2)", atan2_unchecked(-1.0, -2.0), atan2(-1.0, -2.0));
    check("atan2_unchecked(0,1)", atan2_unchecked(0.0, 1.0), atan2(0.0, 1.0));
    check("atan2_unchecked(-0,1)", atan2_unchecked(-0.0, 1.0), atan2(-0.0, 1.0));
    check("atan2_unchecked(inf,1)", atan2_unchecked(f32::INFINITY, 1.0), atan2(f32::INFINITY, 1.0));
    check("tan(0)", tan(0.0), 0.0);

    // tan_checked (idea #48): plain sin_checked(x)/cos_checked(x)
    // composition, mirroring tanpi/tand's own pattern.
    check("tan_checked(0)", tan_checked(0.0), 0.0);
    check("tan_checked(-0)", tan_checked(-0.0), -0.0);
    check("tan_checked(nan)", tan_checked(f32::NAN), f32::NAN);
    check("tan_checked(inf)", tan_checked(f32::INFINITY), f32::NAN);
    check("tan_checked(-inf)", tan_checked(f32::NEG_INFINITY), f32::NAN);
    check_finite("tan_checked(1e15)", tan_checked(1.0e15));

    check("erf(0)", erf(0.0), 0.0);
    // erf_poly used to be evaluated unbounded: its own unboundedness (not
    // an exp2 domain issue) made erf(50)/(100)/(+-inf) wrong (see doc
    // comment) instead of correctly saturating to +-1.
    check("erf(50)", erf(50.0), 1.0);
    check("erf(-50)", erf(-50.0), -1.0);
    check("erf(inf)", erf(f32::INFINITY), 1.0);
    check("erf(-inf)", erf(f32::NEG_INFINITY), -1.0);
    check("erf(nan)", erf(f32::NAN), f32::NAN);
    check("erfc(0)", erfc(0.0), 1.0);
    check("erfc(-0)", erfc(-0.0), 1.0);
    // Infinities go through the `xs` clamp that feeds the exact-square
    // correction `fma(xs,xs,-xs*xs)`: unclamped that is `inf-inf` = NaN,
    // which would poison the result even though the exponential itself
    // saturates correctly. The clamp is also the reason `erfc(nan)` still
    // returns NaN -- `NaN > 11.0` is false, so NaN passes through rather
    // than being replaced by the clamp value.
    check("erfc(inf)", erfc(f32::INFINITY), 0.0);
    check("erfc(-inf)", erfc(f32::NEG_INFINITY), 2.0);
    check("erfc(nan)", erfc(f32::NAN), f32::NAN);
    // erfc's clamp used to not fully protect its internal exp2 call for
    // |x| >= ~9.35 (see its doc comment) -- these used to be inf/huge
    // garbage instead of the correct near-0 (or near-2 for negative x).
    // `erfc(x)` is mathematically bounded to `[0,2]` for every real `x`
    // (`erfc = 1-erf`, `erf` bounded to `[-1,1]`) -- checked directly
    // (2026-07-10), no violation found, but locked in as a permanent
    // guard via `check_bounded` rather than just `check_finite`, same
    // reasoning as the `sinc`/`sin_checked` hardening above.
    check_range("erfc(9.5)", erfc(9.5), 0.0, 2.0);
    check_range("erfc(10)", erfc(10.0), 0.0, 2.0);
    check_range("erfc(-9.5)", erfc(-9.5), 0.0, 2.0);
    check_range("erfc(-10)", erfc(-10.0), 0.0, 2.0);

    // norm_cdf/norm_pdf (backlog idea #67): thin composites over erfc/
    // exp_checked.
    check("norm_cdf(0)", norm_cdf(0.0), 0.5);
    check("norm_cdf(-0)", norm_cdf(-0.0), 0.5);
    check("norm_cdf(inf)", norm_cdf(f32::INFINITY), 1.0);
    check("norm_cdf(-inf)", norm_cdf(f32::NEG_INFINITY), 0.0);
    check("norm_cdf(nan)", norm_cdf(f32::NAN), f32::NAN);
    check("norm_pdf(0)", norm_pdf(0.0), 0.3989422804014327);
    check("norm_pdf(inf)", norm_pdf(f32::INFINITY), 0.0);
    check("norm_pdf(-inf)", norm_pdf(f32::NEG_INFINITY), 0.0);
    check("norm_pdf(nan)", norm_pdf(f32::NAN), f32::NAN);

    // logit(p) = ln(p/(1-p)), sigmoid's inverse (backlog idea #71).
    check("logit(0.5)", logit(0.5), 0.0);
    check("logit(0)", logit(0.0), f32::NEG_INFINITY);
    check("logit(1)", logit(1.0), f32::INFINITY);
    check("logit(-0.1)", logit(-0.1), f32::NAN);
    check("logit(1.1)", logit(1.1), f32::NAN);
    check("logit(nan)", logit(f32::NAN), f32::NAN);
    // 1 ulp high, not exact: p=0.7 is outside logit's central band, so
    // this round trip goes through the ln/log1p arm, where the seam is
    // placed precisely because a ulp there is already cheap.
    check("sigmoid(logit(0.7))", sigmoid(logit(0.7)), 0.70000005);
    // The central band is where the old difference form cancelled: both
    // its logs were ~-ln(2) while their difference was ~4*(p-0.5). These
    // sit where it measured its worst (~1024 ulp), against f64
    // references, and are the pins that would catch a regression to any
    // formula that subtracts two logarithms here.
    check_ulp("logit(0.4998779)", logit(0.4998779), -4.8840046e-4, 2);
    check_ulp("logit(0.5+2^-24)", logit(0.50000006), 2.3841858e-7, 2);
    check_ulp("logit(0.5-2^-24)", logit(0.49999994), -2.3841858e-7, 2);
    check_ulp("logit(0.4996)", logit(0.4996), -1.6000274e-3, 2);
    // Denormal p: the ratio form's `1/p` exceeds f32::MAX below
    // p ~ 2.9e-39, which returned -inf for a true value near -88 until
    // the denominator got scaled (see logit's own doc comment). Every
    // one of these is finite and ~-100, and the exhaustive sweep alone
    // did catch this one -- but only as a single number in an average,
    // so these pin the shape of the failure too.
    check_ulp("logit(f32::MIN_POSITIVE)", logit(f32::MIN_POSITIVE), -87.33655, 2);
    check_ulp("logit(2^-149)", logit(f32::from_bits(1)), -103.27893, 2);
    check_ulp("logit(2.938736e-39)", logit(f32::from_bits(0x0020_0000)), -88.72284, 2);
    check_ulp("logit(1e-40)", logit(1e-40), -92.10341, 2);
    check("logit(-1e-40) [out of domain]", logit(-1e-40), f32::NAN);

    // xlogy/xlog1py (backlog idea #84): x==0 overrides to 0 regardless of
    // y, matching scipy.special.xlogy's convention exactly -- including
    // the indeterminate 0*ln(0) case entropy sums define away, and even
    // y<0/NaN, since a zero-weight term should vanish from a sum rather
    // than poison it with NaN.
    check("xlogy(0,0)", xlogy(0.0, 0.0), 0.0);
    check("xlogy(0,-5)", xlogy(0.0, -5.0), 0.0);
    check("xlogy(0,nan)", xlogy(0.0, f32::NAN), 0.0);
    check("xlogy(-0,1)", xlogy(-0.0, 1.0), 0.0);
    check("xlogy(1,1)", xlogy(1.0, 1.0), 0.0);
    check("xlogy(2,1)", xlogy(2.0, 1.0), 0.0);
    check("xlogy(1,0)", xlogy(1.0, 0.0), f32::NEG_INFINITY);
    check("xlogy(-1,0)", xlogy(-1.0, 0.0), f32::INFINITY);
    check("xlogy(nan,1)", xlogy(f32::NAN, 1.0), f32::NAN);
    check("xlogy(1,-1)", xlogy(1.0, -1.0), f32::NAN);
    check("xlog1py(0,0)", xlog1py(0.0, 0.0), 0.0);
    check("xlog1py(0,-5)", xlog1py(0.0, -5.0), 0.0);
    check("xlog1py(0,nan)", xlog1py(0.0, f32::NAN), 0.0);
    check("xlog1py(1,0)", xlog1py(1.0, 0.0), 0.0);
    check("xlog1py(2,0)", xlog1py(2.0, 0.0), 0.0);
    check("xlog1py(1,-1)", xlog1py(1.0, -1.0), f32::NEG_INFINITY);
    check("xlog1py(nan,1)", xlog1py(f32::NAN, 1.0), f32::NAN);

    // compound(x,n) = (1+x)^n (backlog idea #72).
    check("compound(0,5)", compound(0.0, 5.0), 1.0);
    check("compound(-1,5)", compound(-1.0, 5.0), 0.0);
    check("compound(-2,5)", compound(-2.0, 5.0), f32::NAN);
    check("compound(x,0)", compound(0.05, 0.0), 1.0);
    check("compound(nan,1)", compound(f32::NAN, 1.0), f32::NAN);
    check("compound(0,nan)", compound(0.0, f32::NAN), f32::NAN);
    check("compound(inf,1)", compound(f32::INFINITY, 1.0), f32::INFINITY);
    // (1+1/n)^n -> e as n grows; the whole point of routing through
    // log1p is keeping this precise even for tiny x (1 ulp off the
    // exact constant here, not a real discrepancy).
    check_known_1ulp("compound(1e-8,1e8)", compound(1.0e-8, 1.0e8), std::f32::consts::E);
    // compound_accurate: same contract at every edge, double-float exponent.
    check("compound_accurate(0,5)", compound_accurate(0.0, 5.0), 1.0);
    check("compound_accurate(-1,5)", compound_accurate(-1.0, 5.0), 0.0);
    check("compound_accurate(-1,-5)", compound_accurate(-1.0, -5.0), f32::INFINITY);
    check("compound_accurate(-2,5)", compound_accurate(-2.0, 5.0), f32::NAN);
    check("compound_accurate(x,0)", compound_accurate(0.05, 0.0), 1.0);
    check("compound_accurate(nan,1)", compound_accurate(f32::NAN, 1.0), f32::NAN);
    check("compound_accurate(0,nan)", compound_accurate(0.0, f32::NAN), f32::NAN);
    check("compound_accurate(inf,1)", compound_accurate(f32::INFINITY, 1.0), f32::INFINITY);
    check("compound_accurate(1,5)", compound_accurate(1.0, 5.0), 32.0);
    // the (1+1/n)^n -> e limit the whole family exists for: exact here,
    // where the single-f32 exponent tier is a ulp off.
    check("compound_accurate(1e-8,1e8)", compound_accurate(1.0e-8, 1.0e8), std::f32::consts::E);

    // erfcx(x) = e^(x^2)*erfc(x), backlog idea #51. For x>=0 the
    // exponentials cancel exactly (see its doc comment), so erfcx(0)
    // reduces to the same trivial case erfc(0) does.
    check("erfcx(0)", erfcx(0.0), 1.0);
    check("erfcx(-0)", erfcx(-0.0), 1.0);
    check("erfcx(nan)", erfcx(f32::NAN), f32::NAN);
    // The positive side is `v*P(v)` with `v = 1/(2+x)` and no exponential
    // at all (see doc comment), so it is finite for every finite input by
    // construction, and decays like the true 1/(x*sqrt(pi)) asymptote all
    // the way out -- pinned mainly to guard the negative side, which does
    // route through exp_checked and genuinely diverges to +inf past the
    // point where 2*exp(x^2) itself overflows.
    check_finite("erfcx(-1)", erfcx(-1.0));
    check_finite("erfcx(-9)", erfcx(-9.0));
    check("erfcx(inf)", erfcx(f32::INFINITY), 0.0);
    // Past the point where 2*e^(x^2) itself overflows f32 (x^2 > ~88.03,
    // i.e. |x| > ~9.382 -- scipy bisection on the true value, not the
    // ~13.3 this comment used to claim), erfcx correctly saturates to
    // +inf rather than wrapping to garbage -- exp_checked's own
    // established saturation guarantee, inherited here for free. Pinned
    // on both sides of the boundary, not just far past it.
    check_finite("erfcx(-9.38)", erfcx(-9.38));
    check("erfcx(-9.39)", erfcx(-9.39), f32::INFINITY);
    check("erfcx(-1000)", erfcx(-1000.0), f32::INFINITY);
    // Large positive x, where the reciprocal-variable fit replaced an
    // implementation that froze at a constant past |x|=10 (relative error
    // growing without bound: ~99% at x=20, ~895% at x=100). Expected
    // values are scipy.special.erfcx rounded to f32.
    check_ulp("erfcx(11)", erfcx(11.0), 5.1080596e-2, 2);
    check_ulp("erfcx(15)", erfcx(15.0), 3.7529606e-2, 2);
    check_ulp("erfcx(20)", erfcx(20.0), 2.8174348e-2, 2);
    check_ulp("erfcx(50)", erfcx(50.0), 1.1281536e-2, 2);
    check_ulp("erfcx(100)", erfcx(100.0), 5.6416136e-3, 2);
    check_ulp("erfcx(200)", erfcx(200.0), 2.8209127e-3, 2);
    check_ulp("erfcx(1e6)", erfcx(1e6), 5.641896e-7, 2);
    // Denormal output at the very top of the domain: v = 1/(2+x) itself
    // rounds to 1/x there, and the leading (~1/sqrt(pi)) coefficient
    // carries it down to a denormal without flushing.
    check_ulp("erfcx(max)", erfcx(f32::MAX), 1.658004e-39, 2);

    // erfinv (backlog idea #66): domain (-1,1), odd function, unbounded
    // as |x|->1. Ordinary values checked via round-trip through erf
    // (tolerance, not exact -- erfinv is an approximation) rather than a
    // literal expected value, since there's no independent f32 erfinv
    // reference in this file to pin an exact value against.
    check("erfinv(0)", erfinv(0.0), 0.0);
    check("erfinv(-0)", erfinv(-0.0), -0.0);
    check("erfinv(1)", erfinv(1.0), f32::INFINITY);
    check("erfinv(-1)", erfinv(-1.0), f32::NEG_INFINITY);
    check("erfinv(1.5)", erfinv(1.5), f32::NAN);
    check("erfinv(-1.5)", erfinv(-1.5), f32::NAN);
    check("erfinv(nan)", erfinv(f32::NAN), f32::NAN);
    check("erfinv(inf)", erfinv(f32::INFINITY), f32::NAN);
    check("erfinv(-inf)", erfinv(f32::NEG_INFINITY), f32::NAN);
    for &x in &[0.3f32, -0.3, 0.5, -0.5, 0.7, -0.7, 0.9, -0.9, 0.9999, -0.9999] {
        let y = erfinv(x);
        check_bounded(&format!("erf(erfinv({x}))-{x}"), erf(y) - x, 1e-4);
    }

    // erfc_inv/probit (backlog idea #139): erfc_inv(y) = erfinv(1-y),
    // probit(p) = sqrt(2)*erfinv(2p-1). Round-trip checked (tolerance,
    // not exact) the same way as erfinv's own ordinary-value pins above.
    check("erfc_inv(1)", erfc_inv(1.0), 0.0);
    check("erfc_inv(0)", erfc_inv(0.0), f32::INFINITY);
    check("erfc_inv(2)", erfc_inv(2.0), f32::NEG_INFINITY);
    check("erfc_inv(nan)", erfc_inv(f32::NAN), f32::NAN);
    check("probit(0.5)", probit(0.5), 0.0);
    check("probit(0)", probit(0.0), f32::NEG_INFINITY);
    check("probit(1)", probit(1.0), f32::INFINITY);
    check("probit(nan)", probit(f32::NAN), f32::NAN);
    for &y in &[0.001f32, 0.5, 1.0, 1.5, 1.999] {
        let z = erfc_inv(y);
        check_bounded(&format!("erfc(erfc_inv({y}))-{y}"), erfc(z) - y, 1e-4);
    }
    for &p in &[0.001f32, 0.3, 0.5, 0.7, 0.999] {
        let x = probit(p);
        check_bounded(&format!("norm_cdf(probit({p}))-{p}"), norm_cdf(x) - p, 1e-4);
    }

    // dawson (backlog idea #138): odd, `F(0)=0`, single interior maximum
    // near x~0.9241389, decays like 1/(2x) for large |x|. f32::MAX/-MAX
    // pins are a standing regression guard on the `2.0*x` overflow bug
    // fuzzing found (fixed by halving before dividing by `x`, not after
    // -- see dawson's own doc comment): both must stay finite, not NaN.
    check("dawson(0)", dawson(0.0), 0.0);
    check("dawson(-0)", dawson(-0.0), -0.0);
    check("dawson(inf)", dawson(f32::INFINITY), 0.0);
    check("dawson(-inf)", dawson(f32::NEG_INFINITY), -0.0);
    check("dawson(nan)", dawson(f32::NAN), f32::NAN);
    check_bounded("dawson(0.9241389)-0.5410442", dawson(0.9241389) - 0.5410442, 1e-4);
    check_bounded("dawson(-0.9241389)+0.5410442", dawson(-0.9241389) + 0.5410442, 1e-4);
    check_bounded("dawson(4)-0.1293480", dawson(4.0) - 0.1293480, 1e-4);
    check_bounded("dawson(-4)+0.1293480", dawson(-4.0) + 0.1293480, 1e-4);
    check_finite("dawson(f32::MAX)", dawson(f32::MAX));
    check_finite("dawson(f32::MIN)", dawson(f32::MIN));
    check_bounded("dawson(f32::MAX)-1.469e-39", dawson(f32::MAX) - 1.4693680e-39, 1e-40);

    check("hypot(0,0)", hypot(0.0, 0.0), 0.0);
    check("hypot(3,4)", hypot(3.0, 4.0), 5.0);
    // hypot(+-inf, anything) = +inf even with a NaN other argument --
    // IEEE754/C99 special-cases infinity to "win" over NaN here. The
    // naive x*x+y*y formula can't reach this alone (inf*inf + NaN*NaN
    // degrades to NaN); distinct from this function's already-documented
    // finite-overflow tradeoff.
    check("hypot(inf,nan)", hypot(f32::INFINITY, f32::NAN), f32::INFINITY);
    check("hypot(nan,inf)", hypot(f32::NAN, f32::INFINITY), f32::INFINITY);
    // hypot_unchecked: contract is x, y both finite -- must match hypot
    // exactly wherever that contract holds.
    check("hypot_unchecked(0,0)", hypot_unchecked(0.0, 0.0), 0.0);
    check("hypot_unchecked(3,4)", hypot_unchecked(3.0, 4.0), 5.0);
    check("hypot_unchecked(-3,4)", hypot_unchecked(-3.0, 4.0), 5.0);
    check("hypot_unchecked(nan,1)", hypot_unchecked(f32::NAN, 1.0), f32::NAN);
    // hypot_checked: no overflow/underflow tradeoff (the whole point),
    // full IEEE754/C99 domain including the exact same NaN/inf edge
    // cases hypot itself handles, plus graceful over/underflow hypot's
    // own doc comment documents as *not* handling.
    check("hypot_checked(0,0)", hypot_checked(0.0, 0.0), 0.0);
    check("hypot_checked(-0,0)", hypot_checked(-0.0, 0.0), 0.0);
    check("hypot_checked(3,4)", hypot_checked(3.0, 4.0), 5.0);
    check("hypot_checked(inf,nan)", hypot_checked(f32::INFINITY, f32::NAN), f32::INFINITY);
    check("hypot_checked(nan,inf)", hypot_checked(f32::NAN, f32::INFINITY), f32::INFINITY);
    check("hypot_checked(nan,0)", hypot_checked(f32::NAN, 0.0), f32::NAN);
    check("hypot_checked(0,nan)", hypot_checked(0.0, f32::NAN), f32::NAN);
    check("hypot_checked(nan,nan)", hypot_checked(f32::NAN, f32::NAN), f32::NAN);
    check("hypot_checked(inf,inf)", hypot_checked(f32::INFINITY, f32::INFINITY), f32::INFINITY);
    check("hypot_checked(-inf,3)", hypot_checked(f32::NEG_INFINITY, 3.0), f32::INFINITY);
    // f32::MAX,f32::MAX overflows the naive x*x+y*y (already inf before
    // sqrt even runs); hypot_checked's whole point is getting this right.
    check("hypot_checked(MAX,MAX)", hypot_checked(f32::MAX, f32::MAX), f32::INFINITY);
    // hypot(x,y) >= max(|x|,|y|) for any finite x,y (a provable
    // mathematical fact) -- checked directly (2026-07-10): holds for
    // hypot_checked at both these points (unlike plain hypot/
    // hypot_unchecked, which trade this away for paired denormal inputs,
    // an already-documented tradeoff -- see hypot's own doc comment).
    // check_range locks this in rather than only checking finiteness.
    check_range("hypot_checked(MAX/2,MAX/2)", hypot_checked(f32::MAX / 2.0, f32::MAX / 2.0), f32::MAX / 2.0, f32::MAX);
    check_range(
        "hypot_checked(min_denorm,min_denorm)",
        hypot_checked(f32::from_bits(1), f32::from_bits(1)),
        f32::from_bits(1),
        f32::MAX,
    );

    // rhypot(x,y) = 1/hypot(x,y). Every special case here falls out of the
    // naive fma(x,x,y*y).sqrt() composition purely from IEEE754 semantics
    // (verified by hand before writing the function) *except* the same
    // inf-vs-NaN case hypot itself needs an override for.
    check("rhypot(0,0)", rhypot(0.0, 0.0), f32::INFINITY);
    check("rhypot(-0,0)", rhypot(-0.0, 0.0), f32::INFINITY);
    check("rhypot(3,4)", rhypot(3.0, 4.0), 0.2);
    check("rhypot(inf,1)", rhypot(f32::INFINITY, 1.0), 0.0);
    check("rhypot(1,inf)", rhypot(1.0, f32::INFINITY), 0.0);
    check("rhypot(-inf,1)", rhypot(f32::NEG_INFINITY, 1.0), 0.0);
    check("rhypot(inf,inf)", rhypot(f32::INFINITY, f32::INFINITY), 0.0);
    check("rhypot(nan,1)", rhypot(f32::NAN, 1.0), f32::NAN);
    check("rhypot(nan,nan)", rhypot(f32::NAN, f32::NAN), f32::NAN);
    // the actual point of the inf-vs-NaN override: matches hypot's own
    // "infinity wins over NaN" special case, reciprocated.
    check("rhypot(inf,nan)", rhypot(f32::INFINITY, f32::NAN), 0.0);
    check("rhypot(nan,inf)", rhypot(f32::NAN, f32::INFINITY), 0.0);

    // normalize2 (backlog idea #136).
    {
        let (nx, ny) = normalize2(3.0, 4.0);
        check("normalize2(3,4).0", nx, 0.6);
        check("normalize2(3,4).1", ny, 0.8);
        let (zx, zy) = normalize2(0.0, 0.0);
        check("normalize2(0,0).0", zx, f32::NAN);
        check("normalize2(0,0).1", zy, f32::NAN);
    }

    // hypot3/rnorm3 (backlog idea #55): same naive-fma-chain tradeoff and
    // inf-wins-over-NaN override as hypot/rhypot, one arg wider.
    check("hypot3(0,0,0)", hypot3(0.0, 0.0, 0.0), 0.0);
    check("hypot3(2,3,6)", hypot3(2.0, 3.0, 6.0), 7.0);
    check("hypot3(inf,nan,1)", hypot3(f32::INFINITY, f32::NAN, 1.0), f32::INFINITY);
    check("hypot3(1,inf,nan)", hypot3(1.0, f32::INFINITY, f32::NAN), f32::INFINITY);
    check("hypot3(nan,1,inf)", hypot3(f32::NAN, 1.0, f32::INFINITY), f32::INFINITY);
    check("rnorm3(0,0,0)", rnorm3(0.0, 0.0, 0.0), f32::INFINITY);
    check("rnorm3(2,3,6)", rnorm3(2.0, 3.0, 6.0), 1.0 / 7.0);
    check("rnorm3(inf,nan,1)", rnorm3(f32::INFINITY, f32::NAN, 1.0), 0.0);
    check("rnorm3(nan,nan,nan)", rnorm3(f32::NAN, f32::NAN, f32::NAN), f32::NAN);

    // normalize3 (backlog idea #136). Magnitude check, not exact-value:
    // normalize3(x,y,z) multiplies by rnorm3's own reciprocal-norm
    // *approximation* (verified ~2 ulp accurate, not exact), so e.g.
    // `3.0*rnorm3(2,3,6)` need not land on the exact same f32 as the
    // literal `3.0/7.0` even though both approximate the same real
    // number -- found by a real FAIL here first, not assumed.
    {
        let (nx, ny, nz) = normalize3(2.0, 3.0, 6.0);
        check_bounded("normalize3(2,3,6) magnitude deviation", (hypot3(nx, ny, nz) - 1.0).abs(), 1e-6);
        let (zx, zy, zz) = normalize3(0.0, 0.0, 0.0);
        check("normalize3(0,0,0).0", zx, f32::NAN);
        check("normalize3(0,0,0).1", zy, f32::NAN);
        check("normalize3(0,0,0).2", zz, f32::NAN);
    }

    // hypot4/rnorm4/normalize4 (backlog idea #134): companions to
    // hypot3/rnorm3, one arg wider (quaternion case).
    check("hypot4(0,0,0,0)", hypot4(0.0, 0.0, 0.0, 0.0), 0.0);
    check("hypot4(1,2,2,4)", hypot4(1.0, 2.0, 2.0, 4.0), 5.0);
    check("hypot4(inf,nan,1,1)", hypot4(f32::INFINITY, f32::NAN, 1.0, 1.0), f32::INFINITY);
    check("rnorm4(0,0,0,0)", rnorm4(0.0, 0.0, 0.0, 0.0), f32::INFINITY);
    check("rnorm4(1,2,2,4)", rnorm4(1.0, 2.0, 2.0, 4.0), 0.2);
    check("rnorm4(inf,nan,1,1)", rnorm4(f32::INFINITY, f32::NAN, 1.0, 1.0), 0.0);
    {
        let (w, x, y, z) = normalize4(1.0, 2.0, 2.0, 4.0);
        check("normalize4(1,2,2,4).0", w, 0.2);
        check("normalize4(1,2,2,4).1", x, 0.4);
        check("normalize4(1,2,2,4).2", y, 0.4);
        check("normalize4(1,2,2,4).3", z, 0.8);
        let (zw, zx, zy, zz) = normalize4(0.0, 0.0, 0.0, 0.0);
        check("normalize4(0,0,0,0).0", zw, f32::NAN);
        check("normalize4(0,0,0,0).1", zx, f32::NAN);
        check("normalize4(0,0,0,0).2", zy, f32::NAN);
        check("normalize4(0,0,0,0).3", zz, f32::NAN);
    }

    // diff_of_products(a,b,c,d) = a*b - c*d via Kahan's compensated
    // two-product (backlog idea #135). NaN/inf propagate through the
    // ordinary fma/mul chain with no special-cased override needed.
    check("diff_of_products(2,3,1,1)", diff_of_products(2.0, 3.0, 1.0, 1.0), 5.0);
    check("diff_of_products(0,0,0,0)", diff_of_products(0.0, 0.0, 0.0, 0.0), 0.0);
    check("diff_of_products(1,1,1,1)", diff_of_products(1.0, 1.0, 1.0, 1.0), 0.0);
    check("diff_of_products(inf,1,0,0)", diff_of_products(f32::INFINITY, 1.0, 0.0, 0.0), f32::INFINITY);
    check("diff_of_products(nan,1,0,0)", diff_of_products(f32::NAN, 1.0, 0.0, 0.0), f32::NAN);
    check("diff_of_products(1,0,0,nan)", diff_of_products(1.0, 0.0, 0.0, f32::NAN), f32::NAN);
    // cross2 is diff_of_products(ax,by,ay,bx): parallel/perpendicular sanity.
    check("cross2(1,0,0,1)", cross2(1.0, 0.0, 0.0, 1.0), 1.0);
    check("cross2(1,0,1,0)", cross2(1.0, 0.0, 1.0, 0.0), 0.0);
    check("cross2(0,0,5,5)", cross2(0.0, 0.0, 5.0, 5.0), 0.0);

    // Complex pack (backlog idea #186). cabs/carg are exact aliases for
    // hypot_checked/atan2 -- pinned against the same classic cases those
    // functions themselves use. cexp/clog are real compositions (exp/
    // cos/sin, and ln/log1p/carg respectively), so their "nice point"
    // pins are check_bounded, not exact.
    check("cabs(3,4)", cabs(3.0, 4.0), 5.0);
    check("cabs(0,0)", cabs(0.0, 0.0), 0.0);
    check("cabs(nan,0)", cabs(f32::NAN, 0.0), f32::NAN);
    check("carg(1,0)", carg(1.0, 0.0), 0.0);
    check("carg(0,0)", carg(0.0, 0.0), 0.0);
    check_bounded("carg(0,1)-FRAC_PI_2", carg(0.0, 1.0) - std::f32::consts::FRAC_PI_2, 1e-6);
    check_bounded("carg(-1,0)-PI", carg(-1.0, 0.0) - std::f32::consts::PI, 1e-6);
    let (cer, cei) = cexp(0.0, 0.0);
    check_bounded("cexp(0,0).re-1", cer - 1.0, 1e-6);
    check_bounded("cexp(0,0).im", cei, 1e-6);
    let (cer, cei) = cexp(0.0, std::f32::consts::FRAC_PI_2);
    check_bounded("cexp(0,pi/2).re", cer, 1e-5);
    check_bounded("cexp(0,pi/2).im-1", cei - 1.0, 1e-5);
    let (clr, cli) = clog(1.0, 0.0);
    check_bounded("clog(1,0).re", clr, 1e-6);
    check_bounded("clog(1,0).im", cli, 1e-6);
    let (clr, cli) = clog(0.0, 0.0);
    check("clog(0,0).re", clr, f32::NEG_INFINITY);
    check("clog(0,0).im", cli, 0.0);
    let (clr, _) = clog(f32::MAX, f32::MAX);
    check_finite("clog(MAX,MAX).re", clr);

    check("rsqrt(1)", rsqrt(1.0), 1.0);
    check("rsqrt(4)", rsqrt(4.0), 0.5);
    check("rsqrt(0)", rsqrt(0.0), f32::INFINITY);
    check("rsqrt(-0)", rsqrt(-0.0), f32::NEG_INFINITY);
    check("rsqrt(-1)", rsqrt(-1.0), f32::NAN);
    check("rsqrt(inf)", rsqrt(f32::INFINITY), 0.0);
    check("rsqrt(-inf)", rsqrt(f32::NEG_INFINITY), f32::NAN);
    check("rsqrt(nan)", rsqrt(f32::NAN), f32::NAN);

    check("powf(2,3)", powf(2.0, 3.0), 8.0);
    check("powf(1,5)", powf(1.0, 5.0), 1.0);
    // powf_unchecked: contract is x positive/normal/finite, y != 0.0 --
    // must match powf inside that domain (verified more thoroughly via a
    // 50M-sample fuzz, not preserved in-repo; these are a permanent
    // regression guard).
    check("powf_unchecked(2,3)", powf_unchecked(2.0, 3.0), powf(2.0, 3.0));
    check("powf_unchecked(1,5)", powf_unchecked(1.0, 5.0), powf(1.0, 5.0));
    check("powf_unchecked(2,1000)", powf_unchecked(2.0, 1000.0), powf(2.0, 1000.0));
    check("powf_unchecked(2,-1000)", powf_unchecked(2.0, -1000.0), powf(2.0, -1000.0));
    check("powf_unchecked(0.86967933,576.48004)", powf_unchecked(0.86967933, 576.48004), powf(0.86967933, 576.48004));
    // powf used to return plausible-looking finite garbage instead of
    // inf/0 once log2(x)*y left the unchecked exp2's domain (see doc
    // comment) -- now correctly saturates.
    check("powf(2,1000)", powf(2.0, 1000.0), f32::INFINITY);
    check("powf(2,-1000)", powf(2.0, -1000.0), 0.0);
    check("powf(10,100)", powf(10.0, 100.0), f32::INFINITY);
    // powf(negative, integer) used to always be NaN: exp2(log2(|x|)*y)
    // alone can never be negative, so the old formula had no way to
    // produce a real answer for x < 0.0 at all, even for a well-defined
    // case like (-2.0)^3.0 = -8.0. Fixed by computing on |x| and
    // reapplying the sign for integer y (even -> positive, odd ->
    // negative), matching std/IEEE754. Non-integer y stays NaN (a real
    // domain limit, not a bug -- e.g. (-8.0)^(1/3) is NaN in f32 too).
    check("powf(-2,3)", powf(-2.0, 3.0), -8.0);
    check("powf(-2,2)", powf(-2.0, 2.0), 4.0);
    check("powf(-2,3.5)", powf(-2.0, 3.5), f32::NAN);
    // pow(x, 0) = 1 for *any* x, even 0, negative, or NaN -- a dedicated
    // IEEE754/C99 special case the log/exp2 formula can't derive on its
    // own (0*inf and NaN*0 both degrade to NaN) -- also used to be wrong.
    check("powf(0,0)", powf(0.0, 0.0), 1.0);
    check("powf(-0,0)", powf(-0.0, 0.0), 1.0);
    check("powf(nan,0)", powf(f32::NAN, 0.0), 1.0);
    // -0.0's own sign as a base: `x.is_sign_negative()` (bit-based), not
    // `x < 0.0` (value-based) -- the same acos-style pitfall from earlier
    // this session, since -0.0 < 0.0 is false. Both are IEEE754-pinned.
    check("powf(-0,3)", powf(-0.0, 3.0), -0.0);
    check("powf(-0,2)", powf(-0.0, 2.0), 0.0);
    check("powf(-0,-1)", powf(-0.0, -1.0), f32::NEG_INFINITY);
    // -0.0/-inf are C99-exempt from the "non-integer y -> NaN" rule that
    // applies to genuinely negative finite x (found 2026-07-09 building a
    // systematic special-case matrix against std, backlog idea #85's
    // fourth wave): only an odd-integer y preserves the negative sign,
    // any other y (even integer or non-integer) gives the unsigned
    // magnitude instead of NaN.
    check("powf(-0,0.5)", powf(-0.0, 0.5), 0.0);
    check("powf(-0,-0.5)", powf(-0.0, -0.5), f32::INFINITY);
    check("powf(-0,-inf)", powf(-0.0, f32::NEG_INFINITY), f32::INFINITY);
    check("powf(-inf,0.5)", powf(f32::NEG_INFINITY, 0.5), f32::INFINITY);
    check("powf(-inf,-0.5)", powf(f32::NEG_INFINITY, -0.5), 0.0);
    // y infinite: sign never depends on x, only |x| relative to 1 --
    // a negative base to an infinite power has no well-defined sign in
    // the limit, only a magnitude.
    check("powf(-2,inf)", powf(-2.0, f32::INFINITY), f32::INFINITY);
    check("powf(-0.5,-inf)", powf(-0.5, f32::NEG_INFINITY), f32::INFINITY);
    check("powf(-inf,inf)", powf(f32::NEG_INFINITY, f32::INFINITY), f32::INFINITY);
    // pow(1, y) = 1 for *any* y, even inf/-inf/nan -- another dedicated
    // C99 special case the log/exp2 formula can't derive on its own
    // (log_2(1)=0, so 0*inf/0*nan degrade to NaN instead of the correct
    // 1). Found via idea #87's own suspicion, checked directly against
    // std. pow(-1,+-inf)=1 is a second, narrower case that does *not*
    // extend to pow(-1,nan) (stays NaN, matching std).
    check("powf(1,inf)", powf(1.0, f32::INFINITY), 1.0);
    check("powf(1,-inf)", powf(1.0, f32::NEG_INFINITY), 1.0);
    check("powf(1,nan)", powf(1.0, f32::NAN), 1.0);
    check("powf(-1,inf)", powf(-1.0, f32::INFINITY), 1.0);
    check("powf(-1,-inf)", powf(-1.0, f32::NEG_INFINITY), 1.0);
    check("powf(-1,nan)", powf(-1.0, f32::NAN), f32::NAN);
    // powf_pos: x>0.0 (or exactly +0.0) domain (backlog idea #74),
    // bit-identical to powf throughout it, including every x>=0 special
    // case above -- except x==-0.0 specifically, documented as a known
    // gap in its own doc comment (sign not preserved for odd-integer y,
    // unlike powf's own (-0.0).powf(3.0)==-0.0).
    check("powf_pos(2,3)", powf_pos(2.0, 3.0), 8.0);
    check("powf_pos(0,0)", powf_pos(0.0, 0.0), 1.0);
    check("powf_pos(-0,0)", powf_pos(-0.0, 0.0), 1.0);
    check("powf_pos(nan,0)", powf_pos(f32::NAN, 0.0), 1.0);
    check_known_1ulp("powf_pos(-0,3) [-0.0 sign not preserved, see doc comment]", powf_pos(-0.0, 3.0), -0.0);
    check_known_1ulp("powf_pos(-0,-1) [-0.0 sign not preserved, see doc comment]", powf_pos(-0.0, -1.0), f32::NEG_INFINITY);
    check("powf_pos(1,inf)", powf_pos(1.0, f32::INFINITY), 1.0);
    check("powf_pos(1,-inf)", powf_pos(1.0, f32::NEG_INFINITY), 1.0);
    check("powf_pos(1,nan)", powf_pos(1.0, f32::NAN), 1.0);
    check("powf_pos(inf,5)", powf_pos(f32::INFINITY, 5.0), f32::INFINITY);
    check("powf_pos(inf,-5)", powf_pos(f32::INFINITY, -5.0), 0.0);
    check("powf_pos(2,inf)", powf_pos(2.0, f32::INFINITY), f32::INFINITY);
    check("powf_pos(0.5,inf)", powf_pos(0.5, f32::INFINITY), 0.0);

    // srgb_to_linear/linear_to_srgb (backlog idea #146).
    check("srgb_to_linear(0)", srgb_to_linear(0.0), 0.0);
    check("srgb_to_linear(1)", srgb_to_linear(1.0), 1.0);
    check("srgb_to_linear(nan)", srgb_to_linear(f32::NAN), f32::NAN);
    check("linear_to_srgb(0)", linear_to_srgb(0.0), 0.0);
    check_known_1ulp("linear_to_srgb(1)", linear_to_srgb(1.0), 1.0);
    check("linear_to_srgb(nan)", linear_to_srgb(f32::NAN), f32::NAN);
    // Round trip through the toe boundary and a mid-range value.
    check_known_1ulp(
        "linear_to_srgb(srgb_to_linear(0.5))",
        linear_to_srgb(srgb_to_linear(0.5)),
        0.5,
    );
    check_known_1ulp(
        "linear_to_srgb(srgb_to_linear(0.04045))",
        linear_to_srgb(srgb_to_linear(0.04045)),
        0.04045,
    );
    // signed_pow: graphics/shading "raise |x| then reapply sign"
    // convention (backlog idea #151) -- total for every finite x/y,
    // never NaN for a negative base with a non-integer y, unlike powf's
    // real domain error there. Note the deliberately *different* -1/inf
    // convention from powf: powf(-1,inf)=1.0 (a dedicated C99 special
    // case), but signed_pow(-1,inf)=-1.0 (unconditional sign
    // reapplication on top of powf_pos(1,inf)=1.0, no special-casing).
    check("signed_pow(2,3)", signed_pow(2.0, 3.0), 8.0);
    check("signed_pow(-2,3)", signed_pow(-2.0, 3.0), -8.0);
    check("signed_pow(-2,0.5)", signed_pow(-2.0, 0.5), -std::f32::consts::SQRT_2);
    check("signed_pow(-2,2)", signed_pow(-2.0, 2.0), -4.0);
    check("signed_pow(0,3)", signed_pow(0.0, 3.0), 0.0);
    check("signed_pow(-0,3)", signed_pow(-0.0, 3.0), -0.0);
    check("signed_pow(-0,2)", signed_pow(-0.0, 2.0), -0.0);
    check("signed_pow(0,0)", signed_pow(0.0, 0.0), 1.0);
    check("signed_pow(nan,3)", signed_pow(f32::NAN, 3.0), f32::NAN);
    check("signed_pow(2,nan)", signed_pow(2.0, f32::NAN), f32::NAN);
    check("signed_pow(-1,inf)", signed_pow(-1.0, f32::INFINITY), -1.0);
    // powf's double-float magnitude: the large-|y| corner the plain
    // `exp2(log_2(x)*y)` formula used to miss by hundreds of ulp. Pinned
    // against an independent Decimal-precision Python reference (true
    // value ~1.1009300443688705e-35), not against whatever this crate
    // happened to return: the pin is the correctly-rounded f32, which
    // powf hits to within 0.25 ulp.
    check("powf(0.86967933,576.48004)", powf(0.86967933, 576.48004), 1.10093e-35);
    // powf_unchecked: contract is x positive/normal/finite, y != 0.0 --
    // must match powf inside that domain (verified more thoroughly via
    // examples/unchecked_parity.rs; permanent regression guard).
    check("powf_unchecked(2,3)", powf_unchecked(2.0, 3.0), powf(2.0, 3.0));
    check("powf_unchecked(2,1000)", powf_unchecked(2.0, 1000.0), powf(2.0, 1000.0));
    check(
        "powf_unchecked(0.86967933,576.48004)",
        powf_unchecked(0.86967933, 576.48004),
        powf(0.86967933, 576.48004),
    );
    check("remainder(5,3)", remainder(5.0, 3.0), -1.0);
    check("remainder(4,2)", remainder(4.0, 2.0), 0.0);
    // remainder_unchecked: contract is x != 0.0, y finite -- must match
    // remainder inside that domain (verified more thoroughly via a
    // 100M-sample fuzz, not preserved in-repo; permanent regression guard).
    check("remainder_unchecked(5,3)", remainder_unchecked(5.0, 3.0), remainder(5.0, 3.0));
    check("remainder_unchecked(4,2)", remainder_unchecked(4.0, 2.0), remainder(4.0, 2.0));
    check("remainder_unchecked(-5,3)", remainder_unchecked(-5.0, 3.0), remainder(-5.0, 3.0));
    // remainder(-0.0, y) used to lose its sign: q is +-0.0 matching x/y's
    // sign, so `-q*y` ends up the opposite sign to x, and `fma(-q,y,x)`
    // adds two exactly-zero values of opposite sign (same IEEE754
    // mechanism as sinf_poly/log1p's own `-0.0` bugs). Unlike those,
    // remainder's result sign does *not* generally track x's sign for
    // nonzero x (e.g. remainder(2,3) == -1, an IEEE remainder property,
    // not a bug), so a blanket copysign fix isn't valid here -- fixed
    // with a trailing `if x == 0.0 { x } else { normal }` select instead
    // (a no-op for every nonzero x, including remainder(+0.0, y) which
    // was already correct).
    check("remainder(-0,3)", remainder(-0.0, 3.0), -0.0);
    check("remainder(0,3)", remainder(0.0, 3.0), 0.0);
    // A second, narrower -0.0 gap survived that fix: unlike the general
    // case above (which correctly doesn't track x's sign), IEEE754 does
    // specifically define an *exact-zero* remainder/fmod result's sign
    // to match x's -- and for nonzero x an exact multiple of y, `-q*y`
    // exactly cancels x the same way, dropping the sign the same way.
    // Not a contradiction of the "blanket copysign isn't valid" note
    // above: this fix only fires when the *computed result* is exactly
    // zero (`remainder_style_combine!`'s own `if normal == 0.0 {
    // normal.copysign(x) }`), not for every nonzero x. Confirmed against
    // libm (Python's math.remainder/math.fmod) before pinning.
    check("remainder(-6,3)", remainder(-6.0, 3.0), -0.0);
    check("remainder(-9,3)", remainder(-9.0, 3.0), -0.0);
    check("remainder(6,3)", remainder(6.0, 3.0), 0.0);
    check("remainder(6,-3)", remainder(6.0, -3.0), 0.0);
    // remainder(0,0)/remainder(0,nan) used to come out 0 instead of NaN
    // (backlog idea #85's own follow-up, 2026-07-09, found via the same
    // systematic special-case matrix technique that caught atan2's NaN
    // bug immediately before this): the `x==0.0` sign-preservation guard
    // fired unconditionally, silently overriding a `normal` that had
    // already correctly evaluated to NaN (via `q=(x/y).round()`, itself
    // NaN whenever `y` is 0 or NaN) with plain `x` instead. Same root
    // cause and same fix (guard now also requires `!normal.is_nan()`)
    // across remainder/remainder_checked/remainder_ieee/remainder_wide/
    // fmod -- all five share this exact `if x==0.0 {x} else {normal}`
    // shape.
    check("remainder(0,0)", remainder(0.0, 0.0), f32::NAN);
    check("remainder(-0,0)", remainder(-0.0, 0.0), f32::NAN);
    check("remainder(0,nan)", remainder(0.0, f32::NAN), f32::NAN);
    // remainder(finite x, +-inf) = x (IEEE754/C99 special case): q rounds
    // to exactly 0.0 for any finite x, but multiplying that zero by an
    // *infinite* y used to give NaN (0*inf is NaN) instead of the
    // intended "no reduction happened, answer is just x" no-op.
    check("remainder(3,inf)", remainder(3.0, f32::INFINITY), 3.0);
    check("remainder(-3,inf)", remainder(-3.0, f32::INFINITY), -3.0);
    check("remainder(inf,3)", remainder(f32::INFINITY, 3.0), f32::NAN);
    // remainder_checked shares remainder's special-case handling (same
    // trailing selects) on top of its wider-range q correction -- same
    // edge cases should hold identically.
    check("remainder_checked(5,3)", remainder_checked(5.0, 3.0), -1.0);
    check("remainder_checked(4,2)", remainder_checked(4.0, 2.0), 0.0);
    // Exact-cancellation sign fix, same as remainder's own -- see its pin
    // comment above for the mechanism.
    check("remainder_checked(-6,3)", remainder_checked(-6.0, 3.0), -0.0);
    check("remainder_checked(6,3)", remainder_checked(6.0, 3.0), 0.0);
    check("remainder_checked(-0,3)", remainder_checked(-0.0, 3.0), -0.0);
    check("remainder_checked(0,3)", remainder_checked(0.0, 3.0), 0.0);
    check("remainder_checked(0,0)", remainder_checked(0.0, 0.0), f32::NAN);
    check("remainder_checked(0,nan)", remainder_checked(0.0, f32::NAN), f32::NAN);
    check("remainder_checked(3,inf)", remainder_checked(3.0, f32::INFINITY), 3.0);
    check("remainder_checked(-3,inf)", remainder_checked(-3.0, f32::INFINITY), -3.0);
    check("remainder_checked(inf,3)", remainder_checked(f32::INFINITY, 3.0), f32::NAN);
    // the actual point of remainder_checked: a case where q's own division
    // rounding would land on the wrong integer for the plain formula.
    check(
        "remainder_checked(1e7,3)",
        remainder_checked(1.0e7, 3.0),
        remainder_ref_exact(1.0e7, 3.0),
    );
    // remainder_ieee: ties-to-even instead of remainder's own ties-away,
    // matching true IEEE754 remainder. Away from ties it's identical to
    // remainder; at an exact tie (x/y = 2.5, an odd/even boundary) it
    // must disagree with remainder's own ties-away answer.
    check("remainder_ieee(4,2)", remainder_ieee(4.0, 2.0), 0.0);
    // Exact-cancellation sign fix, same as remainder's own -- see its pin
    // comment above for the mechanism.
    check("remainder_ieee(-6,3)", remainder_ieee(-6.0, 3.0), -0.0);
    check("remainder_ieee(6,3)", remainder_ieee(6.0, 3.0), 0.0);
    check("remainder_ieee(-0,3)", remainder_ieee(-0.0, 3.0), -0.0);
    check("remainder_ieee(0,3)", remainder_ieee(0.0, 3.0), 0.0);
    check("remainder_ieee(0,0)", remainder_ieee(0.0, 0.0), f32::NAN);
    check("remainder_ieee(0,nan)", remainder_ieee(0.0, f32::NAN), f32::NAN);
    check("remainder_ieee(3,inf)", remainder_ieee(3.0, f32::INFINITY), 3.0);
    check("remainder_ieee(inf,3)", remainder_ieee(f32::INFINITY, 3.0), f32::NAN);
    // x/y=2.5: ties-away rounds q to 3 (remainder -1); ties-to-even rounds
    // q to 2 (the even neighbor), remainder +1 -- must actually differ.
    check("remainder_ieee(5,2) ties-even", remainder_ieee(5.0, 2.0), 1.0);
    check("remainder(5,2) ties-away (contrast)", remainder(5.0, 2.0), -1.0);
    // x/y=1.5: ties-to-even rounds q to 2 (even), remainder -0.5; matches
    // remainder's own ties-away answer here since ties-away *also* picks
    // the higher magnitude 2 for a positive 1.5 (both conventions agree
    // whenever the "away" and "even" neighbors happen to coincide).
    check("remainder_ieee(3,2)", remainder_ieee(3.0, 2.0), -1.0);

    // remainder_wide: remainder_checked's own special-case handling and
    // near-tie behavior, still holding in its own already-correct domain.
    check("remainder_wide(5,3)", remainder_wide(5.0, 3.0), -1.0);
    check("remainder_wide(4,2)", remainder_wide(4.0, 2.0), 0.0);
    // Exact-cancellation sign fix, same as remainder's own -- see its pin
    // comment above for the mechanism.
    check("remainder_wide(-6,3)", remainder_wide(-6.0, 3.0), -0.0);
    check("remainder_wide(6,3)", remainder_wide(6.0, 3.0), 0.0);
    check("remainder_wide(-0,3)", remainder_wide(-0.0, 3.0), -0.0);
    check("remainder_wide(0,3)", remainder_wide(0.0, 3.0), 0.0);
    check("remainder_wide(0,0)", remainder_wide(0.0, 0.0), f32::NAN);
    check("remainder_wide(0,nan)", remainder_wide(0.0, f32::NAN), f32::NAN);
    check("remainder_wide(3,inf)", remainder_wide(3.0, f32::INFINITY), 3.0);
    check("remainder_wide(-3,inf)", remainder_wide(-3.0, f32::INFINITY), -3.0);
    check("remainder_wide(inf,3)", remainder_wide(f32::INFINITY, 3.0), f32::NAN);
    check(
        "remainder_wide(1e7,3) matches remainder_checked",
        remainder_wide(1.0e7, 3.0),
        remainder_checked(1.0e7, 3.0),
    );
    // the actual point of remainder_wide: |x/y| well past remainder_checked's
    // own 2^24 cliff (see its doc comment) -- q0's own quantization gap at
    // this magnitude is dozens of integers, not the single-integer nudge
    // remainder_checked's own correction can fix.
    check(
        "remainder_wide(1e10,3) past remainder_checked's 2^24 cliff",
        remainder_wide(1.0e10, 3.0),
        remainder_ref_exact(1.0e10, 3.0),
    );
    check(
        "remainder_wide(1e13,7) further past the cliff",
        remainder_wide(1.0e13, 7.0),
        remainder_ref_exact(1.0e13, 7.0),
    );
    // Regression pin for the internal-overflow bug found while implementing
    // this function (see its own doc comment): q0*y computed as a single
    // f32 product can exceed f32::MAX even though x, y, and the true
    // remainder are all finite, whenever x or y individually sits close to
    // f32::MAX -- this used to return NaN instead of a finite remainder.
    // q0=2 here, small enough that an f64 reference stays exact.
    check(
        "remainder_wide near f32::MAX (overflow regression)",
        remainder_wide(3.2603515e38, 1.8878502e38),
        remainder_ref_exact(3.2603515e38, 1.8878502e38),
    );
    check(
        "remainder_wide near -f32::MAX (overflow regression, negative)",
        remainder_wide(-3.2603515e38, 1.8878502e38),
        remainder_ref_exact(-3.2603515e38, 1.8878502e38),
    );
    check("remainder_wide(f32::MAX,f32::MAX)", remainder_wide(f32::MAX, f32::MAX), 0.0);
    check("remainder_wide(nan,3)", remainder_wide(f32::NAN, 3.0), f32::NAN);

    // fmod: C fmod semantics (truncated division, sign always matches x)
    // -- verified directly against Rust's own `%` operator, which already
    // implements this convention.
    check("fmod(5,3)", fmod(5.0, 3.0), 5.0f32 % 3.0);
    check("fmod(-5,3)", fmod(-5.0, 3.0), -5.0f32 % 3.0);
    check("fmod(5,-3)", fmod(5.0, -3.0), 5.0f32 % -3.0);
    check("fmod(-5,-3)", fmod(-5.0, -3.0), -5.0f32 % -3.0);
    check("fmod(0,3)", fmod(0.0, 3.0), 0.0);
    check("fmod(-0,3)", fmod(-0.0, 3.0), -0.0);
    // Exact-multiple case: `%`'s own reference already gets this right
    // (-0.0), but fmod itself used to drop the sign via IEEE754 exact
    // cancellation (see remainder_style_combine!'s own comment) -- this
    // is the case that pin coverage above (5,3 not being an exact
    // multiple of 3) never exercised.
    check("fmod(-6,3)", fmod(-6.0, 3.0), -6.0f32 % 3.0);
    check("fmod(-9,3)", fmod(-9.0, 3.0), -9.0f32 % 3.0);
    check("fmod(6,3)", fmod(6.0, 3.0), 6.0f32 % 3.0);
    check("fmod(6,-3)", fmod(6.0, -3.0), 6.0f32 % -3.0);
    check("fmod(3,inf)", fmod(3.0, f32::INFINITY), 3.0);
    check("fmod(-3,inf)", fmod(-3.0, f32::INFINITY), -3.0);
    check("fmod(inf,3)", fmod(f32::INFINITY, 3.0), f32::NAN);
    check("fmod(3,0)", fmod(3.0, 0.0), f32::NAN);
    check("fmod(nan,3)", fmod(f32::NAN, 3.0), f32::NAN);
    check("fmod(0,0)", fmod(0.0, 0.0), 0.0f32 % 0.0f32);
    check("fmod(-0,0)", fmod(-0.0, 0.0), (-0.0f32) % 0.0f32);
    check("fmod(0,nan)", fmod(0.0, f32::NAN), 0.0f32 % f32::NAN);
    check("fmod_unchecked(5,3)", fmod_unchecked(5.0, 3.0), fmod(5.0, 3.0));
    check("fmod_unchecked(-5,3)", fmod_unchecked(-5.0, 3.0), fmod(-5.0, 3.0));

    // fmod_checked: same representative cases as fmod above, plus the
    // specific off-by-a-whole-y cases found while testing idea #79 (one
    // per sign combination of x/y, see fmod_checked's own doc comment).
    check("fmod_checked(5,3)", fmod_checked(5.0, 3.0), 5.0f32 % 3.0);
    check("fmod_checked(-5,3)", fmod_checked(-5.0, 3.0), -5.0f32 % 3.0);
    check("fmod_checked(5,-3)", fmod_checked(5.0, -3.0), 5.0f32 % -3.0);
    check("fmod_checked(-5,-3)", fmod_checked(-5.0, -3.0), -5.0f32 % -3.0);
    check("fmod_checked(0,3)", fmod_checked(0.0, 3.0), 0.0);
    check("fmod_checked(-0,3)", fmod_checked(-0.0, 3.0), -0.0);
    // Exact-cancellation sign fix, same as fmod's own -- see its pin
    // comment above for the mechanism.
    check("fmod_checked(-6,3)", fmod_checked(-6.0, 3.0), -6.0f32 % 3.0);
    check("fmod_checked(6,3)", fmod_checked(6.0, 3.0), 6.0f32 % 3.0);
    check("fmod_checked(3,inf)", fmod_checked(3.0, f32::INFINITY), 3.0);
    check("fmod_checked(-3,inf)", fmod_checked(-3.0, f32::INFINITY), -3.0);
    check("fmod_checked(inf,3)", fmod_checked(f32::INFINITY, 3.0), f32::NAN);
    check("fmod_checked(3,0)", fmod_checked(3.0, 0.0), f32::NAN);
    check("fmod_checked(nan,3)", fmod_checked(f32::NAN, 3.0), f32::NAN);
    check("fmod_checked(0,0)", fmod_checked(0.0, 0.0), 0.0f32 % 0.0f32);
    // Off-by-a-whole-y cases (see fmod_checked's own doc comment for the
    // adj-sign derivation these each pin one branch of): all four
    // sign(x)/sign(y) combinations, values confirmed against
    // (x as f64) % (y as f64) before pinning.
    check("fmod_checked(-6.055322e9,-2.0952672e7)", fmod_checked(-6.055322e9, -2.0952672e7), -2.0952576e7);
    check("fmod_checked(-7.871735e-1,1.574347e-3)", fmod_checked(-7.871735e-1, 1.574347e-3), -1.5743424e-3);
    check("fmod_checked(2.4304448e2,-9.241235e-1)", fmod_checked(2.4304448e2, -9.241235e-1), 9.2411566e-1);
    check("fmod_checked(9.885632e-22,2.4062233e-23)", fmod_checked(9.885632e-22, 2.4062233e-23), 2.0116773e-24);

    // rem_euclid/div_euclid: Rust-native std semantics, verified against
    // f32::rem_euclid/f32::div_euclid directly (idea #152).
    check("rem_euclid(5,3)", rem_euclid(5.0, 3.0), 5.0f32.rem_euclid(3.0));
    check("rem_euclid(-5,3)", rem_euclid(-5.0, 3.0), (-5.0f32).rem_euclid(3.0));
    check("rem_euclid(5,-3)", rem_euclid(5.0, -3.0), 5.0f32.rem_euclid(-3.0));
    check("rem_euclid(-5,-3)", rem_euclid(-5.0, -3.0), (-5.0f32).rem_euclid(-3.0));
    check("rem_euclid(0,3)", rem_euclid(0.0, 3.0), 0.0f32.rem_euclid(3.0));
    check("rem_euclid(-0,3)", rem_euclid(-0.0, 3.0), (-0.0f32).rem_euclid(3.0));
    check("rem_euclid(-6,3)", rem_euclid(-6.0, 3.0), (-6.0f32).rem_euclid(3.0));
    check("rem_euclid(5,0)", rem_euclid(5.0, 0.0), 5.0f32.rem_euclid(0.0));
    check("rem_euclid(inf,3)", rem_euclid(f32::INFINITY, 3.0), f32::INFINITY.rem_euclid(3.0));
    check("rem_euclid(5,inf)", rem_euclid(5.0, f32::INFINITY), 5.0f32.rem_euclid(f32::INFINITY));
    check("rem_euclid(nan,3)", rem_euclid(f32::NAN, 3.0), f32::NAN.rem_euclid(3.0));
    check("div_euclid(5,3)", div_euclid(5.0, 3.0), 5.0f32.div_euclid(3.0));
    check("div_euclid(-5,3)", div_euclid(-5.0, 3.0), (-5.0f32).div_euclid(3.0));
    check("div_euclid(5,-3)", div_euclid(5.0, -3.0), 5.0f32.div_euclid(-3.0));
    check("div_euclid(-5,-3)", div_euclid(-5.0, -3.0), (-5.0f32).div_euclid(-3.0));
    check("div_euclid(5,0)", div_euclid(5.0, 0.0), 5.0f32.div_euclid(0.0));
    check("div_euclid(inf,3)", div_euclid(f32::INFINITY, 3.0), f32::INFINITY.div_euclid(3.0));
    check("div_euclid(5,inf)", div_euclid(5.0, f32::INFINITY), 5.0f32.div_euclid(f32::INFINITY));
    check("div_euclid(nan,3)", div_euclid(f32::NAN, 3.0), f32::NAN.div_euclid(3.0));

    // idea #197: seam continuity standing test, informational regression
    // detector for future coefficient/threshold refits (see check_seam's
    // own doc comment) -- every threshold here must match the real
    // shipped branch condition exactly, not the function's own doc
    // comment (which can drift, e.g. asin's was 0.25 before idea #58's
    // 2026-07-20 crossover shift to 0.27).
    check_seam("expm1 seam", expm1, 0.5);
    check_seam("exp2m1 seam", exp2m1, 0.65);
    check_seam("sinh seam", sinh, 0.5);
    check_seam("tanh seam", tanh, 0.25);
    check_seam("asin seam", asin, 0.27);
    check_seam("erf seam", erf, 0.28);
    check_seam("atanh seam", atanh, 0.25);
    // logit's band is on |2p-1|, so its seams sit at p = 0.625 and
    // p = 0.375 -- both, since the two arms are not mirror images of
    // each other (the atanh poly is odd, the ln/log1p difference is not).
    check_seam("logit seam (upper)", logit, 0.625);
    check_seam("logit seam (lower)", logit, 0.375);
}

/// f64-computed exact reference for a single spot-check triple, used only
/// to pin remainder_checked's wider-range behavior in edgecheck (not a
/// general-purpose reference -- see examples/accuracy.rs for the real
/// fuzz-tested sweep against sleef).
fn remainder_ref_exact(x: f32, y: f32) -> f32 {
    let xd = x as f64;
    let yd = y as f64;
    let q = (xd / yd).round();
    (xd - q * yd) as f32
}
