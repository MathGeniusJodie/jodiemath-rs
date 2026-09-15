// Idea #166: which functions produce correctly-rounded denormal outputs and
// which flush or saturate early? `exp2_checked` documents its behaviour;
// most others were unaudited. The saturation-pin sweep (idea #165) already
// caught `sigmoid` returning exactly 0 across a sliver of denormal-scale
// outputs, which is what prompted finishing this properly.
//
// Two distinct cases, audited separately because they fail for different
// reasons:
//
//   (A) normal input -> denormal output. The function's own machinery has to
//       carry a result below `MIN_POSITIVE` without flushing: an exponent
//       field built by a bit trick simply cannot represent it, so this is
//       where early saturation lives.
//   (B) denormal input -> denormal output. Mostly the near-zero identity
//       region (sin(x) ~ x etc.), where the risk is the reverse: a
//       reduction or rescale mangling a subnormal argument.
//
// Metric is relative error, not ulp: at denormal magnitudes one ulp is
// 1.4e-45, so raw ulp reads in the millions for an answer that is merely
// off in its last bits (see idea #165's entry). A result that is exactly 0
// where the truth is a nonzero denormal is reported separately as a flush,
// since relative error is 1.0 for every such case and that hides how early
// the flush starts.
use jodiemath_rs::*;

const MIN_NORM: f64 = f32::MIN_POSITIVE as f64;

struct Row {
    name: &'static str,
    worst_rel: f64,
    worst_x: f32,
    n_denorm: u64,
    n_flushed: u64,
    first_flush_x: f32,
    last_ok_x: f32,
}

fn audit(name: &'static str, f: impl Fn(f32) -> f32, r: impl Fn(f64) -> f64, xs: &[f32]) -> Row {
    let mut row = Row {
        name,
        worst_rel: 0.0,
        worst_x: 0.0,
        n_denorm: 0,
        n_flushed: 0,
        first_flush_x: 0.0,
        last_ok_x: 0.0,
    };
    for &x in xs {
        let want = r(x as f64);
        let a = want.abs();
        // Only interested where the *correctly rounded f32* answer is a
        // nonzero denormal. Testing `a < MIN_NORM` in f64 alone is wrong: it
        // also admits values below ~7e-46, which round to f32 zero, so
        // returning 0 there is correct and counting it as a flush inflates
        // the figure (an earlier version did, reporting erfc at "39%
        // flushed" when it is not premature at all).
        if !want.is_finite() || !(a > 0.0 && a < MIN_NORM) || (a as f32) == 0.0 {
            continue;
        }
        row.n_denorm += 1;
        let got = f(x);
        if got == 0.0 {
            if row.n_flushed == 0 {
                row.first_flush_x = x;
            }
            row.n_flushed += 1;
            continue;
        }
        row.last_ok_x = x;
        let rel = ((got as f64) - want).abs() / a;
        if rel > row.worst_rel {
            row.worst_rel = rel;
            row.worst_x = x;
        }
    }
    row
}

fn report(title: &str, rows: &[Row]) {
    println!("\n=== {title} ===");
    println!(
        "{:<16} {:>9} {:>9} {:>12} {:>14}  {}",
        "function", "denorm", "flushed", "worst rel", "worst x", "note"
    );
    for r in rows {
        if r.n_denorm == 0 {
            println!(
                "{:<16} {:>9} {:>9} {:>12} {:>14}  no denormal outputs in range",
                r.name, 0, "-", "-", "-"
            );
            continue;
        }
        let note = if r.n_flushed == 0 {
            "carries denormals".to_string()
        } else if r.n_flushed == r.n_denorm {
            format!("FLUSHES ALL (from x={:e})", r.first_flush_x)
        } else {
            format!(
                "flushes {:.0}% (first at x={:e})",
                100.0 * r.n_flushed as f64 / r.n_denorm as f64,
                r.first_flush_x
            )
        };
        println!(
            "{:<16} {:>9} {:>9} {:>12.3e} {:>14e}  {}",
            r.name, r.n_denorm, r.n_flushed, r.worst_rel, r.worst_x, note
        );
    }
}

fn main() {
    // (A) normal inputs whose true result is denormal: dense sweeps over the
    // far tails where each function's output collapses toward zero.
    let lin = |lo: f32, hi: f32, n: usize| -> Vec<f32> {
        (0..=n)
            .map(|i| lo + (hi - lo) * (i as f32) / (n as f32))
            .collect()
    };

    let a_rows = vec![
        audit(
            "exp2_checked",
            exp2_checked,
            f64::exp2,
            &lin(-150.0, -125.0, 300_000),
        ),
        audit(
            "exp_checked",
            exp_checked,
            f64::exp,
            &lin(-104.0, -87.0, 300_000),
        ),
        audit(
            "exp10_checked",
            exp10_checked,
            |v| 10f64.powf(v),
            &lin(-45.2, -37.8, 300_000),
        ),
        // NOTE: unchecked `exp2` is deliberately absent. Its doc scopes it to
        // "(non-denormal, finite, nonzero) results only", so it has no
        // denormal outputs *within its own domain* -- sweeping it past -126
        // measures documented garbage (an early version of this audit did,
        // and reported a meaningless 5.8e76 relative error).
        audit(
            "sigmoid",
            sigmoid,
            |v| 1.0 / (1.0 + (-v).exp()),
            &lin(-105.0, -87.0, 300_000),
        ),
        audit("erfc", erfc, |v| libm_erfc(v), &lin(9.0, 10.6, 300_000)),
        audit(
            "norm_pdf",
            norm_pdf,
            |v| (-0.5 * v * v).exp() / (2.0 * std::f64::consts::PI).sqrt(),
            &lin(12.0, 14.5, 300_000),
        ),
        audit(
            "tanh_grad",
            tanh_grad,
            |v| {
                let t = v.tanh();
                1.0 - t * t
            },
            &lin(20.0, 45.0, 300_000),
        ),
        audit(
            "sigmoid_grad",
            sigmoid_grad,
            |v| {
                let s = 1.0 / (1.0 + (-v).exp());
                s * (1.0 - s)
            },
            &lin(60.0, 92.0, 300_000),
        ),
        // The list above was hand-picked and had never been checked for
        // completeness; these five are the rest of the crate's
        // exponential-tailed 1-arg functions, added after `softplus` turned
        // out to flush 17 units of `x` early with nothing reporting it.
        audit(
            "softplus",
            softplus,
            |v| v.exp().ln_1p(),
            &lin(-110.0, -85.0, 300_000),
        ),
        audit(
            "softplus_checked",
            softplus_checked,
            |v| v.exp().ln_1p(),
            &lin(-110.0, -85.0, 300_000),
        ),
        audit(
            "logsigmoid",
            logsigmoid,
            |v| -((-v).exp().ln_1p()),
            &lin(85.0, 110.0, 300_000),
        ),
        audit(
            "logsigmoid_checked",
            logsigmoid_checked,
            |v| -((-v).exp().ln_1p()),
            &lin(85.0, 110.0, 300_000),
        ),
        audit(
            "silu",
            silu,
            |v| v / (1.0 + (-v).exp()),
            &lin(-110.0, -85.0, 300_000),
        ),
        audit(
            "silu_checked",
            silu_checked,
            |v| v / (1.0 + (-v).exp()),
            &lin(-110.0, -85.0, 300_000),
        ),
        // First 2-arg coverage in this file. Both lists here take
        // `fn(f32) -> f32`, so `logaddexp` and every other 2-arg function
        // were simply absent -- and `logaddexp` carries `softplus`'s exact
        // cutoff, just on `|a-b|` instead of `|x|`. Currying `b = 0` is
        // not an approximation of that: `logaddexp(x, 0) == softplus(x)`
        // identically, so the curried sweep is the flush band itself.
        audit(
            "logaddexp(x,0)",
            |x| logaddexp(x, 0.0),
            |v| v.exp().ln_1p(),
            &lin(-110.0, -85.0, 300_000),
        ),
        audit(
            "logaddexp_checked(x,0)",
            |x| logaddexp_checked(x, 0.0),
            |v| v.exp().ln_1p(),
            &lin(-110.0, -85.0, 300_000),
        ),
        audit(
            "logaddexp_accurate(x,0)",
            |x| logaddexp_accurate(x, 0.0),
            |v| v.exp().ln_1p(),
            &lin(-110.0, -85.0, 300_000),
        ),
        audit(
            "gelu",
            gelu,
            |v| v * 0.5 * libm_erfc(-v / std::f64::consts::SQRT_2),
            &lin(-15.5, -12.0, 300_000),
        ),
        audit(
            "norm_cdf",
            norm_cdf,
            |v| 0.5 * libm_erfc(-v / std::f64::consts::SQRT_2),
            &lin(-15.0, -13.0, 300_000),
        ),
    ];
    report("(A) normal input -> denormal output", &a_rows);

    // (B) denormal inputs, near-zero identity region. Sample the denormal
    // bit patterns directly (there are only 2^23 per sign).
    let mut dn: Vec<f32> = Vec::new();
    let mut b: u32 = 1;
    while b < 0x0080_0000 {
        dn.push(f32::from_bits(b));
        dn.push(-f32::from_bits(b));
        b = b.wrapping_add(97); // stride, ~86k samples per sign
    }

    let b_rows = vec![
        audit("sin", sin, f64::sin, &dn),
        audit("sin_wide", sin_wide, f64::sin, &dn),
        audit("tan", tan, f64::tan, &dn),
        audit("asin", asin, f64::asin, &dn),
        audit("atan", atan, f64::atan, &dn),
        audit("sinh", sinh, f64::sinh, &dn),
        audit("asinh", asinh, f64::asinh, &dn),
        audit("atanh", atanh, f64::atanh, &dn),
        audit("expm1", expm1, f64::exp_m1, &dn),
        audit("expm1_checked", expm1_checked, f64::exp_m1, &dn),
        audit("log1p", log1p, f64::ln_1p, &dn),
        audit("erf", erf, |v| libm_erf(v), &dn),
        audit("dawson", dawson, |v| v - v * v * v * 2.0 / 3.0, &dn),
        audit("sinpi", sinpi, |v| (v * std::f64::consts::PI).sin(), &dn),
        audit("sqrt1pm1", sqrt1pm1, |v| v / ((1.0 + v).sqrt() + 1.0), &dn),
        audit(
            "gelu",
            gelu,
            |v| 0.5 * v * (1.0 + libm_erf(v / 2f64.sqrt())),
            &dn,
        ),
        audit("silu", silu, |v| v / (1.0 + (-v).exp()), &dn),
        audit(
            "silu_checked",
            silu_checked,
            |v| v / (1.0 + (-v).exp()),
            &dn,
        ),
        audit("softsign", softsign, |v| v / (1.0 + v.abs()), &dn),
        audit("wrap_pi", wrap_pi, |v| v, &dn),
    ];
    report("(B) denormal input -> denormal output", &b_rows);

    // The percentage of flushed samples is a weak signal on its own -- what
    // matters is how *early* a function flushes relative to where the true
    // result genuinely rounds to zero. A function flushing only in the last
    // representable denormal unit is doing essentially the right thing; one
    // flushing an octave early is losing real answers.
    println!("\n=== how early does each flush start? ===");
    println!(
        "{:<16} {:>16} {:>16} {:>12}",
        "function", "first flush at x", "true zero at x", "premature by"
    );
    let width = |name: &str, f: &dyn Fn(f32) -> f32, r: &dyn Fn(f64) -> f64, lo: f32, hi: f32| {
        // scan in the direction of decreasing |result|
        let n = 400_000;
        let mut first_flush: Option<f32> = None;
        let mut true_zero: Option<f32> = None;
        for i in 0..=n {
            let x = lo + (hi - lo) * (i as f32) / (n as f32);
            let want = r(x as f64);
            let rounds_to_zero = (want.abs() as f32) == 0.0;
            if first_flush.is_none() && f(x) == 0.0 && !rounds_to_zero {
                first_flush = Some(x);
            }
            if true_zero.is_none() && rounds_to_zero {
                true_zero = Some(x);
            }
        }
        match (first_flush, true_zero) {
            (Some(a), Some(b)) => {
                println!("{name:<16} {a:>16e} {b:>16e} {:>11.4} in x", (b - a).abs())
            }
            (Some(a), None) => println!("{name:<16} {a:>16e} {:>16} {:>12}", "(beyond range)", "-"),
            (None, _) => println!("{name:<16} {:>16} {:>16} {:>12}", "never", "-", "none"),
        }
    };
    width("exp2_checked", &exp2_checked, &f64::exp2, -125.0, -152.0);
    width("exp_checked", &exp_checked, &f64::exp, -87.0, -105.0);
    width(
        "exp10_checked",
        &exp10_checked,
        &|v: f64| 10f64.powf(v),
        -37.8,
        -45.5,
    );
    width(
        "sigmoid",
        &sigmoid,
        &|v: f64| 1.0 / (1.0 + (-v).exp()),
        -87.0,
        -106.0,
    );
    width("erfc", &erfc, &libm_erfc, 9.0, 10.8);
    width(
        "norm_pdf",
        &norm_pdf,
        &|v: f64| (-0.5 * v * v).exp() / (2.0 * std::f64::consts::PI).sqrt(),
        12.0,
        14.8,
    );

    // Same five as the (A) table above. Note this helper carries its *own*
    // list -- it is the one that produces the actionable "premature by"
    // figure, and it silently omitted every function that flushes its whole
    // denormal range, which is exactly the set worth looking at. Each scan
    // has to *start* where the function is still correct, so these ranges
    // begin outside the denormal band and walk into it.
    width(
        "softplus",
        &softplus,
        &|v: f64| v.exp().ln_1p(),
        -80.0,
        -110.0,
    );
    width(
        "softplus_checked",
        &softplus_checked,
        &|v: f64| v.exp().ln_1p(),
        -80.0,
        -110.0,
    );
    width(
        "logsigmoid",
        &logsigmoid,
        &|v: f64| -((-v).exp().ln_1p()),
        80.0,
        110.0,
    );
    width(
        "logsigmoid_checked",
        &logsigmoid_checked,
        &|v: f64| -((-v).exp().ln_1p()),
        80.0,
        110.0,
    );
    width(
        "silu",
        &silu,
        &|v: f64| v / (1.0 + (-v).exp()),
        -80.0,
        -112.0,
    );
    width(
        "silu_checked",
        &silu_checked,
        &|v: f64| v / (1.0 + (-v).exp()),
        -80.0,
        -112.0,
    );
    width(
        "gelu",
        &gelu,
        &|v: f64| v * 0.5 * libm_erfc(-v / std::f64::consts::SQRT_2),
        -12.0,
        -16.0,
    );
    width(
        "norm_cdf",
        &norm_cdf,
        &|v: f64| 0.5 * libm_erfc(-v / std::f64::consts::SQRT_2),
        -13.0,
        -15.5,
    );

    let flushers: Vec<&str> = a_rows
        .iter()
        .chain(b_rows.iter())
        .filter(|r| r.n_flushed > 0)
        .map(|r| r.name)
        .collect();
    println!("\n=== summary ===");
    println!("functions that flush some denormal output: {flushers:?}");
}

// erf/erfc in f64 via a series/continued fraction good to ~1e-15 in the
// ranges used above -- std has no erf, and the tails here are exactly where
// a naive 1-erf(x) would cancel to nothing.
fn libm_erfc(x: f64) -> f64 {
    if x < 0.0 {
        return 2.0 - libm_erfc(-x);
    }
    if x < 2.0 {
        return 1.0 - libm_erf(x);
    }
    // continued fraction for the scaled complementary error function
    let mut f = 0.0f64;
    for k in (1..=300).rev() {
        f = (k as f64) / 2.0 / (x + f);
    }
    (-x * x).exp() / (x + f) / std::f64::consts::PI.sqrt()
}

fn libm_erf(x: f64) -> f64 {
    if x.abs() > 2.0 {
        return x.signum() * (1.0 - libm_erfc(x.abs()));
    }
    // Taylor/Maclaurin: erf(x) = 2/sqrt(pi) * sum (-1)^n x^(2n+1)/(n!(2n+1))
    let mut term = x;
    let mut sum = x;
    for n in 1..200 {
        term *= -x * x / (n as f64);
        sum += term / (2 * n + 1) as f64;
        if term.abs() < 1e-300 {
            break;
        }
    }
    2.0 / std::f64::consts::PI.sqrt() * sum
}
