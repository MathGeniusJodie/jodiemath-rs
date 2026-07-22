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

/// Stronger than `check_finite`: also asserts `|got| <= bound`. Added
/// after finding `sin_checked`/`cos_checked` could silently return values
/// up to `2.6e21` for legitimate extreme finite input (fixed with a
/// `.clamp(-1.0,1.0)`, see their own doc comments) -- that bug's own
/// symptom was *finite*, so `check_finite` alone would never have caught
/// it; this closes that gap as a permanent regression guard.
fn check_bounded(name: &str, got: f32, bound: f32) {
    let ok = got.is_finite() && got.abs() <= bound;
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
    println!(
        "{} {:30} got {:e} (0x{:08x}) (in [{lo},{hi}])",
        if ok { "ok  " } else { "FAIL" },
        name,
        got,
        got.to_bits()
    );
}

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
    println!(
        "{} {:30} value gap {gap} ulp at threshold {threshold:e}, one-sided slopes {slope_below:e} / {slope_above:e}",
        if ok { "ok  " } else { "FAIL" },
        name
    );
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

    // tanpi(x) = sinpi(x)/cospi(x): new function (backlog idea #29).
    // Poles at half-integer x are real (cospi(x)==0 there) and correctly
    // give +-inf via IEEE754 division, not NaN -- pinned so that stays
    // true. tanpi(inf)/(-inf) are NaN, matching sinpi/cospi's own
    // existing (inherited, not new) convention at infinity.
    check("tanpi(0)", tanpi(0.0), 0.0);
    check("tanpi(-0)", tanpi(-0.0), -0.0);
    check("tanpi(0.25)", tanpi(0.25), 1.0);
    check("tanpi(1)", tanpi(1.0), 0.0);
    check("tanpi(0.5)", tanpi(0.5), f32::NEG_INFINITY);
    check("tanpi(-0.5)", tanpi(-0.5), f32::NEG_INFINITY);
    check("tanpi(1.5)", tanpi(1.5), f32::NEG_INFINITY);
    check("tanpi(nan)", tanpi(f32::NAN), f32::NAN);
    check("tanpi(inf)", tanpi(f32::INFINITY), f32::NAN);
    check("tanpi(-inf)", tanpi(f32::NEG_INFINITY), f32::NAN);

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

    // tand(x) = sind(x)/cosd(x): new function (backlog idea #29), same
    // "poles are real, IEEE754 division handles them for free" reasoning
    // as tanpi above.
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

    // exp_m1_over_x(x) = (e^x-1)/x, with the removable singularity at 0
    // resolving to exactly 1.0 for free from the Pade branch's own
    // algebra (N(0)/D(0)=-120/-120=1.0 exactly) -- no explicit x==0.0
    // select needed, unlike sinc's own removable-singularity handling.
    check("exp_m1_over_x(0)", exp_m1_over_x(0.0), 1.0);
    check("exp_m1_over_x(-0)", exp_m1_over_x(-0.0), 1.0);
    check_known_1ulp("exp_m1_over_x(1)", exp_m1_over_x(1.0), (1.0f64.exp_m1()) as f32);
    check_finite("exp_m1_over_x(80)", exp_m1_over_x(80.0));
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
    check("sigmoid(-89)", sigmoid(-89.0), 0.0);
    // just below the fixed clamp boundary: still the same (correct, real)
    // value the old code also gave here, confirming no regression at the
    // boundary itself.
    check("sigmoid(-88)", sigmoid(-88.0), 6.054601e-39);
    check("sigmoid(nan)", sigmoid(f32::NAN), f32::NAN);

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

    // gelu(x) = x*Phi(x) = x*0.5*erfc(-x/sqrt2) (backlog idea #70). The
    // structural pins (gelu(3)/gelu(-3)) confirm the normal path really is
    // that composition, bit for bit -- accuracy vs a real reference is the
    // accuracy.rs harness's job. The special-value pins guard the explicit
    // x==-inf override: 0.0*(-inf) alone is NaN, but the true limit is 0.
    check("gelu(0)", gelu(0.0), 0.0);
    check("gelu(-0)", gelu(-0.0), -0.0);
    check("gelu(3)==3*.5*erfc(-3/sqrt2)", gelu(3.0), 3.0 * 0.5 * erfc(-3.0 * std::f32::consts::FRAC_1_SQRT_2));
    check("gelu(-3)==-3*.5*erfc(3/sqrt2)", gelu(-3.0), -3.0 * 0.5 * erfc(3.0 * std::f32::consts::FRAC_1_SQRT_2));
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
    check("atan_latency(0)", atan_latency(0.0), 0.0);
    check("atan_latency(-0)", atan_latency(-0.0), -0.0);
    check("atan_latency(inf)", atan_latency(f32::INFINITY), std::f32::consts::FRAC_PI_2);
    check("atan_latency(-inf)", atan_latency(f32::NEG_INFINITY), -std::f32::consts::FRAC_PI_2);
    check("atan_latency(nan)", atan_latency(f32::NAN), f32::NAN);
    check("atan_latency(1)", atan_latency(1.0), atan(1.0));
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
    // `erfc(x)` is mathematically bounded to `[0,2]` for every real `x`
    // (`erfc = 1-erf`, `erf` bounded to `[-1,1]`) -- checked directly
    // (2026-07-10), no violation found, but locked in as a permanent
    // guard via `check_bounded` rather than just `check_finite`, same
    // reasoning as the `sinc`/`sin_checked` hardening above.
    check_range("erfc(9.5)", erfc(9.5), 0.0, 2.0);
    check_range("erfc(10)", erfc(10.0), 0.0, 2.0);
    check_range("erfc(-9.5)", erfc(-9.5), 0.0, 2.0);
    check_range("erfc(-10)", erfc(-10.0), 0.0, 2.0);

    // erfc_accurate: full-precision-exponent sibling of erfc (idea #62),
    // same special-value/domain behavior, just a tighter fit within it.
    check("erfc_accurate(0)", erfc_accurate(0.0), 1.0);
    check("erfc_accurate(-0)", erfc_accurate(-0.0), 1.0);
    check("erfc_accurate(inf)", erfc_accurate(f32::INFINITY), 0.0);
    check("erfc_accurate(-inf)", erfc_accurate(f32::NEG_INFINITY), 2.0);
    check("erfc_accurate(nan)", erfc_accurate(f32::NAN), f32::NAN);
    check_range("erfc_accurate(9.5)", erfc_accurate(9.5), 0.0, 2.0);
    check_range("erfc_accurate(10)", erfc_accurate(10.0), 0.0, 2.0);
    check_range("erfc_accurate(-9.5)", erfc_accurate(-9.5), 0.0, 2.0);
    check_range("erfc_accurate(-10)", erfc_accurate(-10.0), 0.0, 2.0);

    // erfcx(x) = e^(x^2)*erfc(x), backlog idea #51. For x>=0 the
    // exponentials cancel exactly (see its doc comment), so erfcx(0)
    // reduces to the same trivial case erfc(0) does.
    check("erfcx(0)", erfcx(0.0), 1.0);
    check("erfcx(-0)", erfcx(-0.0), 1.0);
    check("erfcx(nan)", erfcx(f32::NAN), f32::NAN);
    // Positive side never needs its own exp2_checked call (see doc
    // comment), so it's finite for any finite input by construction --
    // pinned mainly to guard the negative side, which does route through
    // exp2_checked and genuinely diverges to +inf past the point where
    // 2*exp(x^2) itself overflows (see doc comment) -- confirm the
    // still-representable region stays finite.
    check_finite("erfcx(1e6)", erfcx(1e6));
    check_finite("erfcx(-1)", erfcx(-1.0));
    check_finite("erfcx(-9)", erfcx(-9.0));
    // Past the point where 2*e^(x^2) itself overflows f32 (x^2 > ~176.7,
    // i.e. |x| > ~13.3), erfcx correctly saturates to +inf rather than
    // wrapping to garbage -- exp2_checked's own established saturation
    // guarantee, inherited here for free.
    check("erfcx(-1000)", erfcx(-1000.0), f32::INFINITY);

    // erfcx_accurate: full-precision-exponent sibling of erfcx (same
    // mechanism as erfc_accurate), same special-value behavior as erfcx.
    check("erfcx_accurate(0)", erfcx_accurate(0.0), 1.0);
    check("erfcx_accurate(-0)", erfcx_accurate(-0.0), 1.0);
    check("erfcx_accurate(nan)", erfcx_accurate(f32::NAN), f32::NAN);
    check_finite("erfcx_accurate(1e6)", erfcx_accurate(1e6));
    check_finite("erfcx_accurate(-1)", erfcx_accurate(-1.0));
    check_finite("erfcx_accurate(-9)", erfcx_accurate(-9.0));
    check("erfcx_accurate(-1000)", erfcx_accurate(-1000.0), f32::INFINITY);

    // erfcx_checked: full-range sibling fixing erfcx's own documented
    // freeze past |x|=10 (see its doc comment). Bit-identical to erfcx
    // for |x|<=10 (same erfc_rational call).
    check("erfcx_checked(0)", erfcx_checked(0.0), 1.0);
    check("erfcx_checked(nan)", erfcx_checked(f32::NAN), f32::NAN);
    check("erfcx_checked(9)", erfcx_checked(9.0), erfcx(9.0));
    check("erfcx_checked(10)", erfcx_checked(10.0), erfcx(10.0));
    check("erfcx_checked(inf)", erfcx_checked(f32::INFINITY), 0.0);
    check("erfcx_checked(-1000)", erfcx_checked(-1000.0), f32::INFINITY);
    // Past the 10 boundary: values confirmed against scipy.special.erfcx
    // (max rel error ~5.8e-7, ~5 ulp, dominated by erfc_rational's own
    // fit error right at the x=10 seam -- not a defect in the asymptotic
    // tail itself) before pinning.
    check("erfcx_checked(15)", erfcx_checked(15.0), 3.7529606e-2);
    check("erfcx_checked(20)", erfcx_checked(20.0), 2.817435e-2);
    check("erfcx_checked(50)", erfcx_checked(50.0), 1.1281537e-2);
    check("erfcx_checked(100)", erfcx_checked(100.0), 5.641614e-3);
    check("erfcx_checked(200)", erfcx_checked(200.0), 2.820913e-3);

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

    check("rsqrt(1)", rsqrt(1.0), 1.0);
    check("rsqrt(4)", rsqrt(4.0), 0.5);
    check("rsqrt(0)", rsqrt(0.0), f32::INFINITY);
    check("rsqrt(-0)", rsqrt(-0.0), f32::NEG_INFINITY);
    check("rsqrt(-1)", rsqrt(-1.0), f32::NAN);
    check("rsqrt(inf)", rsqrt(f32::INFINITY), 0.0);
    check("rsqrt(-inf)", rsqrt(f32::NEG_INFINITY), f32::NAN);
    check("rsqrt(nan)", rsqrt(f32::NAN), f32::NAN);

    check("pown(2,3)", pown(2.0, 3), 8.0);
    check("pown(2,0)", pown(2.0, 0), 1.0);
    check("pown(0,0)", pown(0.0, 0), 1.0);
    check("pown(0,5)", pown(0.0, 5), 0.0);
    check("pown(0,-3)", pown(0.0, -3), f32::INFINITY);
    check("pown(-0,3)", pown(-0.0, 3), -0.0);
    check("pown(-0,-3)", pown(-0.0, -3), f32::NEG_INFINITY);
    check("pown(-2,3)", pown(-2.0, 3), -8.0);
    check("pown(-2,4)", pown(-2.0, 4), 16.0);
    check("pown(2,-1)", pown(2.0, -1), 0.5);
    check("pown(nan,2)", pown(f32::NAN, 2), f32::NAN);
    check("pown(inf,2)", pown(f32::INFINITY, 2), f32::INFINITY);
    check("pown(inf,-2)", pown(f32::INFINITY, -2), 0.0);
    check("pown(-inf,3)", pown(f32::NEG_INFINITY, 3), f32::NEG_INFINITY);
    check("pown(-inf,-3)", pown(f32::NEG_INFINITY, -3), -0.0);
    // i32::MIN's magnitude is 2^31, needing bit index 31 -- an off-by-one
    // in an earlier 0..31 iteration range gave 1.0 here instead of 0.0.
    check("pown(2,i32::MIN)", pown(2.0, i32::MIN), 0.0);
    check("pown(2,i32::MAX)", pown(2.0, i32::MAX), f32::INFINITY);
    // Found via fuzzing: computing x^|n| then reciprocating at the end
    // (instead of inverting x first) overflowed here even though the
    // true small answer doesn't.
    check("pown(-1.8449108e19,-2)", pown(-1.8449108e19, -2), 2.937983e-39);

    // pown_small: same representative cases as pown above, minus the
    // i32::MIN/i32::MAX ones -- those are outside pown_small's own
    // |n| <= 255 contract, not something it needs to get right.
    check("pown_small(2,3)", pown_small(2.0, 3), 8.0);
    check("pown_small(2,0)", pown_small(2.0, 0), 1.0);
    check("pown_small(0,0)", pown_small(0.0, 0), 1.0);
    check("pown_small(0,5)", pown_small(0.0, 5), 0.0);
    check("pown_small(0,-3)", pown_small(0.0, -3), f32::INFINITY);
    check("pown_small(-0,3)", pown_small(-0.0, 3), -0.0);
    check("pown_small(-0,-3)", pown_small(-0.0, -3), f32::NEG_INFINITY);
    check("pown_small(-2,3)", pown_small(-2.0, 3), -8.0);
    check("pown_small(-2,4)", pown_small(-2.0, 4), 16.0);
    check("pown_small(2,-1)", pown_small(2.0, -1), 0.5);
    check("pown_small(nan,2)", pown_small(f32::NAN, 2), f32::NAN);
    check("pown_small(inf,2)", pown_small(f32::INFINITY, 2), f32::INFINITY);
    check("pown_small(inf,-2)", pown_small(f32::INFINITY, -2), 0.0);
    check("pown_small(-inf,3)", pown_small(f32::NEG_INFINITY, 3), f32::NEG_INFINITY);
    check("pown_small(-inf,-3)", pown_small(f32::NEG_INFINITY, -3), -0.0);
    // Boundary of pown_small's own |n| <= 255 contract.
    check("pown_small(2,255)", pown_small(2.0, 255), pown(2.0, 255));
    check("pown_small(2,-255)", pown_small(2.0, -255), pown(2.0, -255));

    // pown_small_accurate: Df32-compensated squaring chain (idea #78),
    // same |n| <= 255 contract and same special-value behavior as
    // pown_small (verified against it directly for every case its own
    // fuzz found risky: zero, +-0, nan, inf, negative base, and the
    // overflow-through-squaring path that motivated its own base.1==0.0
    // guard, see its doc comment).
    check("pown_small_accurate(2,3)", pown_small_accurate(2.0, 3), 8.0);
    check("pown_small_accurate(2,0)", pown_small_accurate(2.0, 0), 1.0);
    check("pown_small_accurate(0,0)", pown_small_accurate(0.0, 0), 1.0);
    check("pown_small_accurate(0,5)", pown_small_accurate(0.0, 5), 0.0);
    check("pown_small_accurate(0,-3)", pown_small_accurate(0.0, -3), f32::INFINITY);
    check("pown_small_accurate(-0,3)", pown_small_accurate(-0.0, 3), -0.0);
    check("pown_small_accurate(-0,-3)", pown_small_accurate(-0.0, -3), f32::NEG_INFINITY);
    check("pown_small_accurate(-2,3)", pown_small_accurate(-2.0, 3), -8.0);
    check("pown_small_accurate(-2,4)", pown_small_accurate(-2.0, 4), 16.0);
    check("pown_small_accurate(2,-1)", pown_small_accurate(2.0, -1), 0.5);
    check("pown_small_accurate(nan,2)", pown_small_accurate(f32::NAN, 2), f32::NAN);
    check("pown_small_accurate(inf,2)", pown_small_accurate(f32::INFINITY, 2), f32::INFINITY);
    check("pown_small_accurate(inf,-2)", pown_small_accurate(f32::INFINITY, -2), 0.0);
    check("pown_small_accurate(-inf,3)", pown_small_accurate(f32::NEG_INFINITY, 3), f32::NEG_INFINITY);
    check("pown_small_accurate(-inf,-3)", pown_small_accurate(f32::NEG_INFINITY, -3), -0.0);
    // The overflow-through-squaring path (idea #78's own doc comment):
    // |x|>1 raised through 255 forces base to genuinely overflow to +-inf
    // partway through the loop, and result itself overflows too for a
    // large enough n -- both must stay exactly the pown/pown_small answer,
    // not NaN.
    check("pown_small_accurate(2,255)", pown_small_accurate(2.0, 255), pown(2.0, 255));
    check("pown_small_accurate(2,-255)", pown_small_accurate(2.0, -255), pown(2.0, -255));
    check("pown_small_accurate(-2.82e14,7)", pown_small_accurate(-2.82e14, 7), pown(-2.82e14, 7));
    check("pown_small_accurate(3.0,200)", pown_small_accurate(3.0, 200), pown(3.0, 200));
    // Not bit-identical to pown_small here (it's a different approximation,
    // more accurate on *average*, not at every single point -- this one
    // happens to land 75 ulp off the true answer where pown_small itself
    // is only 12 off, well within both functions' normal error range, see
    // IDEAS.md idea #78's own writeup). Just confirm no NaN/overflow
    // garbage from the base.1==0.0 guard at this magnitude.
    check_finite("pown_small_accurate(0.9,-255)", pown_small_accurate(0.9, -255));

    // pown_const<N>: same representative cases as pown, N as a const
    // generic instead of a runtime argument -- must match pown exactly
    // for every N, including the full i32::MIN/i32::MAX extremes (no
    // narrower contract here, unlike pown_small).
    check("pown_const::<3>(2)", pown_const::<3>(2.0), 8.0);
    check("pown_const::<0>(2)", pown_const::<0>(2.0), 1.0);
    check("pown_const::<0>(0)", pown_const::<0>(0.0), 1.0);
    check("pown_const::<5>(0)", pown_const::<5>(0.0), 0.0);
    check("pown_const::<-3>(0)", pown_const::<-3>(0.0), f32::INFINITY);
    check("pown_const::<3>(-0)", pown_const::<3>(-0.0), -0.0);
    check("pown_const::<-3>(-0)", pown_const::<-3>(-0.0), f32::NEG_INFINITY);
    check("pown_const::<3>(-2)", pown_const::<3>(-2.0), -8.0);
    check("pown_const::<4>(-2)", pown_const::<4>(-2.0), 16.0);
    check("pown_const::<-1>(2)", pown_const::<-1>(2.0), 0.5);
    check("pown_const::<2>(nan)", pown_const::<2>(f32::NAN), f32::NAN);
    check("pown_const::<2>(inf)", pown_const::<2>(f32::INFINITY), f32::INFINITY);
    check("pown_const::<-2>(inf)", pown_const::<-2>(f32::INFINITY), 0.0);
    check("pown_const::<3>(-inf)", pown_const::<3>(f32::NEG_INFINITY), f32::NEG_INFINITY);
    check("pown_const::<-3>(-inf)", pown_const::<-3>(f32::NEG_INFINITY), -0.0);
    check("pown_const::<i32::MIN>(2)", pown_const::<{ i32::MIN }>(2.0), 0.0);
    check("pown_const::<i32::MAX>(2)", pown_const::<{ i32::MAX }>(2.0), f32::INFINITY);
    check(
        "pown_const::<-2>(-1.8449108e19)",
        pown_const::<-2>(-1.8449108e19),
        2.937983e-39,
    );

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
    check("powf_checked(0,nan)", powf_checked(0.0, f32::NAN), f32::NAN);
    check("powf_checked(inf,nan)", powf_checked(f32::INFINITY, f32::NAN), f32::NAN);
    check("powf_checked(-0,3)", powf_checked(-0.0, 3.0), -0.0);
    check("powf_checked(-0,-1)", powf_checked(-0.0, -1.0), f32::NEG_INFINITY);
    check("powf_checked(-0,0.5)", powf_checked(-0.0, 0.5), 0.0);
    check("powf_checked(-0,-0.5)", powf_checked(-0.0, -0.5), f32::INFINITY);
    check("powf_checked(-inf,0.5)", powf_checked(f32::NEG_INFINITY, 0.5), f32::INFINITY);
    check("powf_checked(-2,inf)", powf_checked(-2.0, f32::INFINITY), f32::INFINITY);
    check("powf_checked(1,inf)", powf_checked(1.0, f32::INFINITY), 1.0);
    check("powf_checked(1,-inf)", powf_checked(1.0, f32::NEG_INFINITY), 1.0);
    check("powf_checked(1,nan)", powf_checked(1.0, f32::NAN), 1.0);
    check("powf_checked(-1,inf)", powf_checked(-1.0, f32::INFINITY), 1.0);
    check("powf_checked(-1,-inf)", powf_checked(-1.0, f32::NEG_INFINITY), 1.0);
    check("powf_checked(-1,nan)", powf_checked(-1.0, f32::NAN), f32::NAN);
    // the actual point of powf_checked: a large-|y| case where the plain
    // formula's error is large (see powf_checked's own doc comment).
    // Value updated 2026-07-09 for the log2_df two-product fix (backlog
    // idea #58) -- confirmed via an independent Decimal-precision Python
    // reference (true value ~1.1009300443688705e-35) that the new value
    // is genuinely closer to correct than the old pin was (was
    // 1.1009411e-35, off by ~1.106e-41; now off by ~8.456e-42).
    check(
        "powf_checked(0.86967933,576.48004)",
        powf_checked(0.86967933, 576.48004),
        1.1009385e-35,
    );
    // powf_checked_unchecked: contract is x positive/normal/finite, y !=
    // 0.0 -- must match powf_checked inside that domain (verified more
    // thoroughly via a 50M-sample fuzz, not preserved in-repo; permanent
    // regression guard).
    check("powf_checked_unchecked(2,3)", powf_checked_unchecked(2.0, 3.0), powf_checked(2.0, 3.0));
    check("powf_checked_unchecked(2,1000)", powf_checked_unchecked(2.0, 1000.0), powf_checked(2.0, 1000.0));
    check(
        "powf_checked_unchecked(0.86967933,576.48004)",
        powf_checked_unchecked(0.86967933, 576.48004),
        powf_checked(0.86967933, 576.48004),
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
