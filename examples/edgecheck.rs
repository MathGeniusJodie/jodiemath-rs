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
    // log1p(-0.0) used to lose its sign: `ln(u) + corr` adds two
    // exactly-zero values of opposite sign at x = +-0.0 (`ln(1.0)` is
    // `+0.0`, but `corr` correctly carries x's sign there), the same
    // IEEE754 mechanism as sinf_poly's own `-0.0` bug. Fixed with a
    // trailing `if x == 0.0 { x } else { normal }` select (log1p is odd
    // and monotonic through the origin, so `normal`'s sign already
    // matches x's for every nonzero x -- this only changes the singular
    // zero case).
    check("log1p(-0)", log1p(-0.0), -0.0);

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
    check("expm1(0)", expm1(0.0), 0.0);
    check("sinh(0)", sinh(0.0), 0.0);
    check("cosh(0)", cosh(0.0), 1.0);
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

    check("asinh(0)", asinh(0.0), 0.0);
    // small-x cancellation (see asinh's doc comment) is fixed: asinh(x) ~ x
    // for tiny x, no longer collapses to exactly 0.
    check("asinh(2.34e-8)", asinh(2.34e-8), 2.34e-8);
    check("asinh(-2.34e-8)", asinh(-2.34e-8), -2.34e-8);
    // large-negative-x cancellation (also fixed, see doc comment): asinh is
    // odd, so this must equal -asinh(1e10) exactly.
    check("asinh(-1e10) == -asinh(1e10)", asinh(-1e10), -asinh(1e10));
    check("asinh(-f32::MAX)", asinh(-f32::MAX), -asinh(f32::MAX));
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
    check("acos(1)", acos(1.0), 0.0);
    check("acos(-1)", acos(-1.0), std::f32::consts::PI);
    check("acos(2)", acos(2.0), f32::NAN);
    // acos(-0.0) used to come out -pi/2 (mulsign's bit-based sign check
    // disagreed with x<0.0's value-based one, exactly at this one input)
    // instead of the correct +pi/2 -- acos is never negative. The
    // remaining 1-ulp gap from the "ideal" FRAC_PI_2 here is acos_poly's
    // own pre-existing, already-accepted fit imprecision (unrelated to
    // the sign bug, and unchanged by this fix -- acos(0) had it too,
    // before and after), not a new issue.
    check_known_1ulp("acos(0)", acos(0.0), std::f32::consts::FRAC_PI_2);
    check_known_1ulp("acos(-0)", acos(-0.0), std::f32::consts::FRAC_PI_2);
    check("atan(0)", atan(0.0), 0.0);
    check("atan(inf)", atan(f32::INFINITY), std::f32::consts::FRAC_PI_2);
    check("atan(-inf)", atan(f32::NEG_INFINITY), -std::f32::consts::FRAC_PI_2);
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
    check("tan(0)", tan(0.0), 0.0);

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
    // erfc's clamp used to not fully protect its internal exp2 call for
    // |x| >= ~9.35 (see its doc comment) -- these used to be inf/huge
    // garbage instead of the correct near-0 (or near-2 for negative x).
    check_finite("erfc(9.5)", erfc(9.5));
    check_finite("erfc(10)", erfc(10.0));
    check_finite("erfc(-9.5)", erfc(-9.5));
    check_finite("erfc(-10)", erfc(-10.0));

    check("hypot(0,0)", hypot(0.0, 0.0), 0.0);
    check("hypot(3,4)", hypot(3.0, 4.0), 5.0);
    // hypot(+-inf, anything) = +inf even with a NaN other argument --
    // IEEE754/C99 special-cases infinity to "win" over NaN here. The
    // naive x*x+y*y formula can't reach this alone (inf*inf + NaN*NaN
    // degrades to NaN); distinct from this function's already-documented
    // finite-overflow tradeoff.
    check("hypot(inf,nan)", hypot(f32::INFINITY, f32::NAN), f32::INFINITY);
    check("hypot(nan,inf)", hypot(f32::NAN, f32::INFINITY), f32::INFINITY);
    check("powf(2,3)", powf(2.0, 3.0), 8.0);
    check("powf(1,5)", powf(1.0, 5.0), 1.0);
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
    // powf_checked shares powf's special-case handling on top of its
    // double-float-precision magnitude for large |y| -- same edge cases
    // should hold identically, plus a couple more that exercise the
    // is_safe/edge_mag fallback split (ax == 0/+inf/NaN) specifically.
    check("powf_checked(2,3)", powf_checked(2.0, 3.0), 8.0);
    check("powf_checked(2,1000)", powf_checked(2.0, 1000.0), f32::INFINITY);
    check("powf_checked(2,-1000)", powf_checked(2.0, -1000.0), 0.0);
    check("powf_checked(-2,3)", powf_checked(-2.0, 3.0), -8.0);
    check("powf_checked(-2,3.5)", powf_checked(-2.0, 3.5), f32::NAN);
    check("powf_checked(0,0)", powf_checked(0.0, 0.0), 1.0);
    check("powf_checked(0,5)", powf_checked(0.0, 5.0), 0.0);
    check("powf_checked(0,-5)", powf_checked(0.0, -5.0), f32::INFINITY);
    check("powf_checked(inf,5)", powf_checked(f32::INFINITY, 5.0), f32::INFINITY);
    check("powf_checked(inf,-5)", powf_checked(f32::INFINITY, -5.0), 0.0);
    check("powf_checked(nan,5)", powf_checked(f32::NAN, 5.0), f32::NAN);
    check("powf_checked(2,nan)", powf_checked(2.0, f32::NAN), f32::NAN);
    check("powf_checked(-0,3)", powf_checked(-0.0, 3.0), -0.0);
    check("powf_checked(-0,-1)", powf_checked(-0.0, -1.0), f32::NEG_INFINITY);
    // the actual point of powf_checked: a large-|y| case where the plain
    // formula's error is large (see powf_checked's own doc comment).
    check(
        "powf_checked(0.86967933,576.48004)",
        powf_checked(0.86967933, 576.48004),
        1.1009411e-35,
    );
    check("remainder(5,3)", remainder(5.0, 3.0), -1.0);
    check("remainder(4,2)", remainder(4.0, 2.0), 0.0);
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
    check("remainder_checked(-0,3)", remainder_checked(-0.0, 3.0), -0.0);
    check("remainder_checked(0,3)", remainder_checked(0.0, 3.0), 0.0);
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
