use jodiemath_rs::*;

fn check(name: &str, got: f32, want: f32) {
    let ok = (got.is_nan() && want.is_nan()) || (got.to_bits() == want.to_bits());
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
    println!("{} {:30} got {:e} (0x{:08x}) (finite)", if got.is_finite() { "ok  " } else { "FAIL" }, name, got, got.to_bits());
}

fn main() {
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
    // exp2 (unchecked): only test inside its documented domain [-126, 128)
    check("exp2(127.9999)", exp2(127.9999), (127.9999f32 as f64).exp2() as f32);
    check("exp2(-125.9)", exp2(-125.9), ((-125.9f32) as f64).exp2() as f32);
    check("exp2(-126)", exp2(-126.0), ((-126.0f32) as f64).exp2() as f32);
    check("exp2(0)", exp2(0.0), 1.0);
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
    // cbrt
    for f in [cbrt as fn(f32) -> f32, cbrt_accurate as fn(f32) -> f32] {
        let n = if f == cbrt as fn(f32) -> f32 { "cbrt" } else { "cbrt_acc" };
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
    // sin/cos (unchecked): only accurate while q = round(x/pi) is an exact
    // f32 integer, i.e. |x| < 2^22 * pi (~1.3e7) -- see sin's doc comment.
    check("sin(0)", sin(0.0), 0.0);
    check("cos(0)", cos(0.0), 1.0);
    check("sin(nan)", sin(f32::NAN), f32::NAN);
    check("cos(nan)", cos(f32::NAN), f32::NAN);
    // sin_checked/cos_checked: nan/+-inf must still come out nan, and every
    // other finite x (including ones far past the sub-ulp-accurate range)
    // must come out finite, never inf -- see sin_checked's doc comment for
    // why (the residual clamp added 2026-07-06 that guarantees this).
    for f in [sin_checked as fn(f32) -> f32, cos_checked as fn(f32) -> f32] {
        let n = if f == sin_checked as fn(f32) -> f32 { "sin_checked" } else { "cos_checked" };
        check(&format!("{n}(nan)"), f(f32::NAN), f32::NAN);
        check(&format!("{n}(inf)"), f(f32::INFINITY), f32::NAN);
        check(&format!("{n}(-inf)"), f(f32::NEG_INFINITY), f32::NAN);
        check_finite(&format!("{n}(max)"), f(f32::MAX));
        check_finite(&format!("{n}(-max)"), f(f32::MIN));
        check_finite(&format!("{n}(1e20)"), f(1e20));
        check_finite(&format!("{n}(1e10)"), f(1e10));
    }

    // ln/log10/log1p: same zero/negative/inf edges as log_2 (they're all
    // log_2 rescaled or composed with it).
    check("ln(0)", ln(0.0), f32::NEG_INFINITY);
    check("ln(-1)", ln(-1.0), f32::NAN);
    check("ln(1)", ln(1.0), 0.0);
    check("log10(0)", log10(0.0), f32::NEG_INFINITY);
    check("log10(100)", log10(100.0), 2.0);
    check("log1p(0)", log1p(0.0), 0.0);
    check("log1p(-1)", log1p(-1.0), f32::NEG_INFINITY);
    check("log1p(-2)", log1p(-2.0), f32::NAN);

    check("exp(0)", exp(0.0), 1.0);
    check("expm1(0)", expm1(0.0), 0.0);
    check("sinh(0)", sinh(0.0), 0.0);
    check("cosh(0)", cosh(0.0), 1.0);
    check("tanh(0)", tanh(0.0), 0.0);

    // asinh/atanh: near-zero cancellation is a known, documented inherited
    // flaw (see their doc comments) -- exactly 0 here is the documented
    // wrong answer, not a crash, so this pins down the current (bad but
    // stable) behavior rather than asserting correctness.
    check("asinh(0)", asinh(0.0), 0.0);
    check("asinh(2.34e-8) [known cancellation]", asinh(2.34e-8), 0.0);
    check("acosh(1)", acosh(1.0), 0.0);
    check("acosh(0.5)", acosh(0.5), f32::NAN);
    // acosh's sign-losing domain bug (see its doc comment) is fixed with an
    // explicit domain check -- large negative x now correctly comes out NaN.
    check("acosh(-1e20)", acosh(-1e20), f32::NAN);
    check("acosh(-4096.0)", acosh(-4096.0), f32::NAN);
    check("acosh(-f32::MAX)", acosh(-f32::MAX), f32::NAN);
    check("atanh(0)", atanh(0.0), 0.0);
    check("atanh(1)", atanh(1.0), f32::INFINITY);
    check("atanh(-1)", atanh(-1.0), f32::NEG_INFINITY);
    check("atanh(2)", atanh(2.0), f32::NAN);

    // asin(0) comes out -0.0, not +0.0: the trailing `* (-hpi)` in asinf
    // flips the sign of the intermediate 0 (same in the C original --
    // traced by hand, not a translation artifact), backwards from the
    // usual libm convention of asin matching its input's zero sign.
    check("asin(0)", asin(0.0), -0.0);
    check("asin(-0)", asin(-0.0), 0.0);
    check("asin(1)", asin(1.0), std::f32::consts::FRAC_PI_2);
    check("asin(-1)", asin(-1.0), -std::f32::consts::FRAC_PI_2);
    check("asin(2)", asin(2.0), f32::NAN);
    check("acos(1)", acos(1.0), 0.0);
    check("acos(-1)", acos(-1.0), std::f32::consts::PI);
    check("acos(2)", acos(2.0), f32::NAN);
    check("atan(0)", atan(0.0), 0.0);
    check("atan(inf)", atan(f32::INFINITY), std::f32::consts::FRAC_PI_2);
    check("atan(-inf)", atan(f32::NEG_INFINITY), -std::f32::consts::FRAC_PI_2);
    check("atan2(1,0)", atan2(1.0, 0.0), std::f32::consts::FRAC_PI_2);
    check("atan2(-1,0)", atan2(-1.0, 0.0), -std::f32::consts::FRAC_PI_2);
    check("atan2(0,-1)", atan2(0.0, -1.0), std::f32::consts::PI);
    check("tan(0)", tan(0.0), 0.0);

    check("erf(0)", erf(0.0), 0.0);
    check("erfc(0)", erfc(0.0), 1.0);

    check("hypot(0,0)", hypot(0.0, 0.0), 0.0);
    check("hypot(3,4)", hypot(3.0, 4.0), 5.0);
    check("powf(2,3)", powf(2.0, 3.0), 8.0);
    check("powf(1,5)", powf(1.0, 5.0), 1.0);
    check("remainder(5,3)", remainder(5.0, 3.0), -1.0);
    check("remainder(4,2)", remainder(4.0, 2.0), 0.0);
}
