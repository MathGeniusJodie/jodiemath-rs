// Idea #165: systematic saturation-boundary pins. Every input clamp and
// overflow threshold gets swept at +-N ulp around it and compared against
// an f64 reference, generalizing the exp10_checked bug where a round-based
// reduction returned a *finite* 3.237e38 for exp10_checked(inf) instead of
// inf -- a boundary-only failure that neither fuzzing nor mca samples.
//
// Why a window rather than a single pin: a misplaced clamp is usually off
// by a small number of ulp, so the failure sits just inside or just
// outside the constant, not at it. Sweeping +-64 ulp catches "the clamp is
// one ulp too generous" as well as "the clamp saturates too early".
//
// Reference is f64: at these magnitudes the true result is either well
// within f32 range or overflows, and f64 has ~8 orders of headroom either
// way, so `(ref as f32)` is the correctly-rounded answer including the
// overflow-to-inf transition.
#![allow(clippy::approx_constant)]
use jodiemath_rs::*;

fn ulps_away(x: f32, n: i32) -> f32 {
    if x == 0.0 || !x.is_finite() {
        return x;
    }
    let b = x.to_bits() as i64;
    // move along the magnitude axis: away from zero for positive n
    let signed = if x > 0.0 { b + n as i64 } else { b - n as i64 };
    f32::from_bits(signed as u32)
}

fn same(a: f32, b: f32) -> bool {
    if a.is_nan() && b.is_nan() {
        return true;
    }
    a.to_bits() == b.to_bits()
}

fn ulp_err(got: f32, want_f64: f64) -> f64 {
    let want = want_f64 as f32;
    if same(got, want) {
        return 0.0;
    }
    if !got.is_finite() || !want.is_finite() {
        // one saturated and the other didn't: a boundary failure
        return f64::INFINITY;
    }
    let a = want.abs();
    let ulp = if a == 0.0 {
        f32::from_bits(1) as f64
    } else {
        (f32::from_bits(a.to_bits() + 1) - a) as f64
    };
    ((got as f64) - (want as f64)).abs() / ulp
}

struct Target {
    name: &'static str,
    f: fn(f32) -> f32,
    r: fn(f64) -> f64,
    // clamp/threshold constants to sweep around
    bounds: &'static [f32],
    // tolerated ulp inside the domain (each function's own documented max)
    tol: f64,
    // Skip points whose *true* result is a denormal. Only set where the
    // function documents early saturation on a denormal-output sliver --
    // `sigmoid` does, verbatim: "for x in roughly (-104.7,-88.7) the true
    // answer is a nonzero denormal but this returns exactly 0". Scoped to
    // denormal-magnitude outputs rather than to an x range, so a new
    // failure at normal output magnitudes still fails the gate.
    skip_denormal_out: bool,
}

fn main() {
    let targets: Vec<Target> = vec![
        Target {
            name: "exp2_checked",
            f: exp2_checked,
            r: f64::exp2,
            bounds: &[-151.0, 128.0, -126.0, -149.0],
            tol: 1.0,
            skip_denormal_out: false,
        },
        Target {
            name: "exp10_checked",
            f: exp10_checked,
            r: |x| 10f64.powf(x),
            bounds: &[-45.154503, 38.53184],
            tol: 1.0,
            skip_denormal_out: false,
        },
        Target {
            name: "exp_checked",
            f: exp_checked,
            r: f64::exp,
            bounds: &[-104.66522426455174, 88.72283911167308],
            tol: 3.0,
            skip_denormal_out: false,
        },
        Target {
            name: "expm1_checked",
            f: expm1_checked,
            r: f64::exp_m1,
            bounds: &[-86.0, 88.72283911167308],
            tol: 6.0,
            skip_denormal_out: false,
        },
        Target {
            name: "exp2m1",
            f: exp2m1,
            r: |x| (x * std::f64::consts::LN_2).exp_m1(),
            bounds: &[-151.0, 128.0],
            tol: 4.0,
            skip_denormal_out: false,
        },
        Target {
            name: "exp10m1",
            f: exp10m1,
            r: |x| (x * std::f64::consts::LN_10).exp_m1(),
            bounds: &[-37.0, 38.53184],
            tol: 4.0,
            skip_denormal_out: false,
        },
        Target {
            name: "tanh",
            f: tanh,
            r: f64::tanh,
            bounds: &[-43.5, 44.0],
            tol: 3.0,
            skip_denormal_out: false,
        },
        Target {
            name: "sigmoid",
            f: sigmoid,
            r: |x| 1.0 / (1.0 + (-x).exp()),
            bounds: &[-88.722839111673, 87.0],
            tol: 3.0,
            skip_denormal_out: true,
        },
        Target {
            name: "sinh_checked",
            f: sinh_checked,
            r: f64::sinh,
            bounds: &[-89.4, 89.4, -170.0, 170.0],
            tol: 5.0,
            skip_denormal_out: false,
        },
        Target {
            name: "cosh_checked",
            f: cosh_checked,
            r: f64::cosh,
            bounds: &[-89.4, 89.4, -170.0, 170.0],
            tol: 5.0,
            skip_denormal_out: false,
        },
        Target {
            name: "erfcx",
            f: erfcx,
            r: |x| {
                // erfcx(x) = exp(x^2)*erfc(x); use the asymptotic form for
                // large x where exp(x^2) overflows f64.
                if x > 30.0 {
                    let t = 1.0 / (x * x);
                    1.0 / (x * std::f64::consts::PI.sqrt())
                        * (1.0 - 0.5 * t + 0.75 * t * t)
                } else {
                    f64::NAN // not asserted below 30
                }
            },
            bounds: &[],
            tol: f64::INFINITY,
            skip_denormal_out: false,
        },
    ];

    const WINDOW: i32 = 64;
    let mut total_bad = 0usize;

    for t in &targets {
        if t.bounds.is_empty() {
            continue;
        }
        let mut worst = 0.0f64;
        let mut worst_x = 0.0f32;
        let mut sat_bad: Vec<String> = Vec::new();
        let mut n = 0u64;
        let mut denorm_skipped = 0u64;

        for &b in t.bounds {
            for k in -WINDOW..=WINDOW {
                let x = ulps_away(b, k);
                if !x.is_finite() {
                    continue;
                }
                n += 1;
                let got = (t.f)(x);
                let want64 = (t.r)(x as f64);
                if want64.is_nan() {
                    continue;
                }
                if t.skip_denormal_out
                    && want64.abs() < f32::MIN_POSITIVE as f64
                    && want64 != 0.0
                {
                    denorm_skipped += 1;
                    continue;
                }
                let e = ulp_err(got, want64);
                if e.is_infinite() {
                    // saturation disagreement: one side inf/finite, other not
                    sat_bad.push(format!(
                        "{}({:e}) = {:e}, f64 ref -> {:e}",
                        t.name, x, got, want64 as f32
                    ));
                } else if e > worst {
                    worst = e;
                    worst_x = x;
                }
            }
        }

        let over = worst > t.tol;
        let bad = !sat_bad.is_empty() || over;
        println!(
            "{:<16} {:>5} pts around {:?}  worst {:.2} ulp (tol {}) at {:e}{}",
            t.name,
            n,
            t.bounds,
            worst,
            t.tol,
            worst_x,
            if bad { "   <-- CHECK" } else { "" }
        );
        if denorm_skipped > 0 {
            println!(
                "      ({denorm_skipped} pts skipped: true result is a denormal, documented early-saturation gap)"
            );
        }
        for s in sat_bad.iter().take(6) {
            println!("      saturation mismatch: {s}");
        }
        if sat_bad.len() > 6 {
            println!("      ... and {} more", sat_bad.len() - 6);
        }
        if bad {
            total_bad += 1;
        }
    }

    // Direct saturation pins: past every upper bound the answer must be
    // the mathematical limit, not a finite leftover.
    println!("\n=== saturation limits far outside the clamp ===");
    let limits: &[(&str, fn(f32) -> f32, f32, f32)] = &[
        ("exp2_checked(+big)", exp2_checked, 1e30, f32::INFINITY),
        ("exp2_checked(-big)", exp2_checked, -1e30, 0.0),
        ("exp10_checked(+big)", exp10_checked, 1e30, f32::INFINITY),
        ("exp10_checked(-big)", exp10_checked, -1e30, 0.0),
        ("exp_checked(+big)", exp_checked, 1e30, f32::INFINITY),
        ("exp_checked(-big)", exp_checked, -1e30, 0.0),
        ("expm1_checked(+big)", expm1_checked, 1e30, f32::INFINITY),
        ("expm1_checked(-big)", expm1_checked, -1e30, -1.0),
        ("exp2m1(+big)", exp2m1, 1e30, f32::INFINITY),
        ("exp2m1(-big)", exp2m1, -1e30, -1.0),
        ("exp10m1(+big)", exp10m1, 1e30, f32::INFINITY),
        ("exp10m1(-big)", exp10m1, -1e30, -1.0),
        ("tanh(+big)", tanh, 1e30, 1.0),
        ("tanh(-big)", tanh, -1e30, -1.0),
        ("sigmoid(+big)", sigmoid, 1e30, 1.0),
        ("sigmoid(-big)", sigmoid, -1e30, 0.0),
        ("sinh_checked(+big)", sinh_checked, 1e30, f32::INFINITY),
        ("sinh_checked(-big)", sinh_checked, -1e30, f32::NEG_INFINITY),
        ("cosh_checked(+big)", cosh_checked, 1e30, f32::INFINITY),
        ("cosh_checked(-big)", cosh_checked, -1e30, f32::INFINITY),
    ];
    for &(name, f, x, want) in limits {
        let got = f(x);
        let ok = same(got, want);
        println!(
            "  {:<22} {:>12e}  want {:>12e}  {}",
            name,
            got,
            want,
            if ok { "ok" } else { "MISMATCH" }
        );
        if !ok {
            total_bad += 1;
        }
    }

    println!("\n{}", if total_bad == 0 { "PASS" } else { "FAIL" });
    if total_bad != 0 {
        std::process::exit(1);
    }
}
