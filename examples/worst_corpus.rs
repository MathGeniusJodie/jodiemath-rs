// Idea #168: worst-case corpus regression gate. The hours-long exhaustive
// sweeps are the real accuracy authority, but they are far too slow to run
// per commit. This is the fast counterpart: evaluate every public 1-arg
// function at a fixed corpus of historically-hard inputs and assert the
// results are *bit-identical* to a blessed golden file.
//
// Deliberately reference-free. A ulp gate would need an f64 reference per
// function (100+ of them, the reason #168 sat in the backlog); comparing
// against blessed bits instead detects any behavioural change at the inputs
// that have historically mattered, in milliseconds, with no reference
// machinery. It answers "did anything move?", not "is it correct?" -- the
// sweeps answer the second question.
//
// The corpus is the special values, a spread of magnitudes across the whole
// exponent range, the branch seams this crate actually has (0.25/0.27/0.28/
// 0.5/0.65/2048/...), and the specific worst-x values recorded in IDEAS.md
// and the readme.
//
// Usage:
//   cargo run --release --example worst_corpus            # check
//   cargo run --release --example worst_corpus -- --bless # regenerate
// Any intentional accuracy change is expected to fail this gate; re-bless
// and eyeball the diff, which is the point.
//
// KNOWN LIMITATION, measured rather than assumed -- do not read a pass as
// "nothing changed". The gate was validated by deliberately perturbing the
// shared Pade coefficient, and it catches a decisive change (0.5% on that
// coefficient moves 61 of the 9720 entries, exit 1). But two *smaller* real
// changes slipped through:
//   - a 1-ulp change to that same coefficient (-1.9999927 -> -1.9999925):
//     invisible, because it multiplies `v*v` against a -120.0 term, so at
//     these inputs it shifts the sum by ~1e-10 relative, far under f32's
//     ~6e-8 resolution.
//   - moving exp10m1's branch seam 0.2 -> 0.21: invisible, because the two
//     branches agree bit-for-bit at the corpus point x=0.2 (a well-placed
//     seam is *supposed* to, which is exactly why seam points make weak
//     canaries).
// So this detects behavioural changes at these inputs, not arbitrarily small
// ones, and a corpus gate is only ever as good as its input list. The
// exhaustive sweeps in accuracy.rs remain the authority; this is a cheap
// early-warning net between them.
//
// One deliberate exemption from bit-exactness: a NaN result compared against
// a NaN golden always passes, whatever the payload/sign/quiet bits are. A
// NaN is a NaN -- payload bits carry no numeric meaning, so a payload-only
// move is not a behavioural change. This matches the rule every ulp metric
// in the repo follows (see accuracy.rs's `ulp_diff`). Values are still
// compared bit-for-bit. Validated the same way the gate itself was, by
// deliberately perturbing the golden file rather than by inspection: a
// payload-only edit (7fc00000 -> ffc0dead) passes, the same entry edited to
// a finite value (3f800000) still fails, and an ordinary entry moved by one
// ulp still fails. Note the shipped corpus only feeds canonical NaNs in, so
// this exemption is about what a function *returns*, not what it is given.
#![allow(clippy::approx_constant)]
use jodiemath_rs::*;
use std::io::Write;

const GOLDEN: &str = "examples/support/worst_corpus.golden";

fn corpus() -> Vec<f32> {
    let mut v: Vec<f32> = vec![
        // specials
        0.0, -0.0, f32::INFINITY, f32::NEG_INFINITY, f32::NAN,
        f32::MIN_POSITIVE, -f32::MIN_POSITIVE,
        f32::from_bits(1), -f32::from_bits(1),
        f32::from_bits(0x007f_ffff), f32::MAX, f32::MIN,
        // branch seams present in this crate
        0.25, -0.25, 0.27, -0.27, 0.28, -0.28, 0.5, -0.5, 0.65, -0.65,
        0.2, -0.2, 1.0, -1.0, 2048.0, -2048.0, 3.288051, -3.288051,
        // clamp boundaries
        88.72283911167308, -104.66522426455174, -86.0, 128.0, -151.0,
        38.53184, -45.154503, 44.0, -43.5, 87.0, -88.722839111673,
        170.0, -170.0, 89.4, -89.4,
        // recorded worst-x values (IDEAS.md / readme / this session)
        0.9652361, -0.39914432, -0.057932023, 4.1295314e-7, 8953.539,
        0.25323957, 0.1505, 2.480704e-1, 9.184352e5, -1.3169037e7,
        1.031, 0.0155, 0.111, 0.965, -2.3,
        // magnitude spread across the exponent range
        1e-38, 1e-30, 1e-20, 1e-10, 1e-5, 1e-2, 3.0, 10.0, 100.0,
        1e5, 1e10, 1e20, 1e30, 1e38,
        -1e-38, -1e-30, -1e-10, -1e-2, -3.0, -100.0, -1e10, -1e30, -1e38,
        // pi-related, where the trig reductions are hardest
        std::f32::consts::PI, -std::f32::consts::PI,
        std::f32::consts::FRAC_PI_2, -std::f32::consts::FRAC_PI_2,
        std::f32::consts::TAU, 1e6, 1.3176794e7, 2.6e7,
    ];
    v.dedup();
    v
}

fn main() {
    let bless = std::env::args().any(|a| a == "--bless");
    let fns: Vec<(&str, fn(f32) -> f32)> = vec![
        ("acos", acos as fn(f32) -> f32),
        ("acosd", acosd as fn(f32) -> f32),
        ("acosh", acosh as fn(f32) -> f32),
        ("acospi", acospi as fn(f32) -> f32),
        ("asin", asin as fn(f32) -> f32),
        ("asind", asind as fn(f32) -> f32),
        ("asinh", asinh as fn(f32) -> f32),
        ("asinpi", asinpi as fn(f32) -> f32),
        ("atan", atan as fn(f32) -> f32),
        ("atan_bounded", atan_bounded as fn(f32) -> f32),
        ("atand", atand as fn(f32) -> f32),
        ("atanh", atanh as fn(f32) -> f32),
        ("atan_latency", atan_latency as fn(f32) -> f32),
        ("atanpi", atanpi as fn(f32) -> f32),
        ("cbrt", cbrt as fn(f32) -> f32),
        ("cbrt_accurate", cbrt_accurate as fn(f32) -> f32),
        ("cbrt_accurate_unchecked", cbrt_accurate_unchecked as fn(f32) -> f32),
        ("cbrt_approx", cbrt_approx as fn(f32) -> f32),
        ("cbrt_fast", cbrt_fast as fn(f32) -> f32),
        ("cbrt_normal", cbrt_normal as fn(f32) -> f32),
        ("cbrt_throughput", cbrt_throughput as fn(f32) -> f32),
        ("cbrt_unchecked", cbrt_unchecked as fn(f32) -> f32),
        ("cos", cos as fn(f32) -> f32),
        ("cos2pi", cos2pi as fn(f32) -> f32),
        ("cos_checked", cos_checked as fn(f32) -> f32),
        ("cosd", cosd as fn(f32) -> f32),
        ("cosd_unchecked", cosd_unchecked as fn(f32) -> f32),
        ("cosh", cosh as fn(f32) -> f32),
        ("cosh_checked", cosh_checked as fn(f32) -> f32),
        ("coshm1", coshm1 as fn(f32) -> f32),
        ("cosh_narrow", cosh_narrow as fn(f32) -> f32),
        ("cosh_throughput", cosh_throughput as fn(f32) -> f32),
        ("cospi", cospi as fn(f32) -> f32),
        ("dawson", dawson as fn(f32) -> f32),
        ("erf", erf as fn(f32) -> f32),
        ("erfc", erfc as fn(f32) -> f32),
        ("erfcx", erfcx as fn(f32) -> f32),
        ("erfinv", erfinv as fn(f32) -> f32),
        ("exp", exp as fn(f32) -> f32),
        ("exp10", exp10 as fn(f32) -> f32),
        ("exp10_checked", exp10_checked as fn(f32) -> f32),
        ("exp10m1", exp10m1 as fn(f32) -> f32),
        ("exp2", exp2 as fn(f32) -> f32),
        ("exp2_approx", exp2_approx as fn(f32) -> f32),
        ("exp2_checked", exp2_checked as fn(f32) -> f32),
        ("exp2m1", exp2m1 as fn(f32) -> f32),
        ("exp_checked", exp_checked as fn(f32) -> f32),
        ("expm1", expm1 as fn(f32) -> f32),
        ("expm1_checked", expm1_checked as fn(f32) -> f32),
        ("expm1_narrow", expm1_narrow as fn(f32) -> f32),
        ("exp_m1_over_x", exp_m1_over_x as fn(f32) -> f32),
        ("exp_m1_over_x_narrow", exp_m1_over_x_narrow as fn(f32) -> f32),
        ("exp_narrow", exp_narrow as fn(f32) -> f32),
        ("fast_round_int", fast_round_int as fn(f32) -> f32),
        ("gelu", gelu as fn(f32) -> f32),
        ("ln", ln as fn(f32) -> f32),
        ("ln_unchecked", ln_unchecked as fn(f32) -> f32),
        ("log10", log10 as fn(f32) -> f32),
        ("log10p1", log10p1 as fn(f32) -> f32),
        ("log10_unchecked", log10_unchecked as fn(f32) -> f32),
        ("log1p", log1p as fn(f32) -> f32),
        ("log1pmx", log1pmx as fn(f32) -> f32),
        ("log_2", log_2 as fn(f32) -> f32),
        ("log2_approx", log2_approx as fn(f32) -> f32),
        ("log2p1", log2p1 as fn(f32) -> f32),
        ("log_2_unchecked", log_2_unchecked as fn(f32) -> f32),
        ("logsigmoid", logsigmoid as fn(f32) -> f32),
        ("norm_cdf", norm_cdf as fn(f32) -> f32),
        ("norm_pdf", norm_pdf as fn(f32) -> f32),
        ("pow_2_3", pow_2_3 as fn(f32) -> f32),
        ("pow_3_2", pow_3_2 as fn(f32) -> f32),
        ("rcbrt", rcbrt as fn(f32) -> f32),
        ("rcp_approx", rcp_approx as fn(f32) -> f32),
        ("rsqrt", rsqrt as fn(f32) -> f32),
        ("rsqrt_approx", rsqrt_approx as fn(f32) -> f32),
        ("sigmoid", sigmoid as fn(f32) -> f32),
        ("sigmoid_fast", sigmoid_fast as fn(f32) -> f32),
        ("sigmoid_grad", sigmoid_grad as fn(f32) -> f32),
        ("silu", silu as fn(f32) -> f32),
        ("sin", sin as fn(f32) -> f32),
        ("sin2pi", sin2pi as fn(f32) -> f32),
        ("sinc", sinc as fn(f32) -> f32),
        ("sin_checked", sin_checked as fn(f32) -> f32),
        ("sinc_unnormalized", sinc_unnormalized as fn(f32) -> f32),
        ("sind", sind as fn(f32) -> f32),
        ("sind_unchecked", sind_unchecked as fn(f32) -> f32),
        ("sinh", sinh as fn(f32) -> f32),
        ("sinh_checked", sinh_checked as fn(f32) -> f32),
        ("sinh_narrow", sinh_narrow as fn(f32) -> f32),
        ("sinh_throughput", sinh_throughput as fn(f32) -> f32),
        ("sinpi", sinpi as fn(f32) -> f32),
        ("sinpi_unchecked", sinpi_unchecked as fn(f32) -> f32),
        ("softplus", softplus as fn(f32) -> f32),
        ("softsign", softsign as fn(f32) -> f32),
        ("sqrt1pm1", sqrt1pm1 as fn(f32) -> f32),
        ("sqrt_approx", sqrt_approx as fn(f32) -> f32),
        ("tan", tan as fn(f32) -> f32),
        ("tan2pi", tan2pi as fn(f32) -> f32),
        ("tan_checked", tan_checked as fn(f32) -> f32),
        ("tand", tand as fn(f32) -> f32),
        ("tand_unchecked", tand_unchecked as fn(f32) -> f32),
        ("tanh", tanh as fn(f32) -> f32),
        ("tanh_grad", tanh_grad as fn(f32) -> f32),
        ("tanpi", tanpi as fn(f32) -> f32),
        ("wrap_pi", wrap_pi as fn(f32) -> f32),
    ];

    let xs = corpus();
    let mut lines: Vec<String> = Vec::new();
    for (name, f) in &fns {
        for &x in &xs {
            lines.push(format!("{} {:08x} {:08x}", name, x.to_bits(), f(x).to_bits()));
        }
    }

    if bless {
        std::fs::create_dir_all("examples/support").ok();
        let mut out = std::fs::File::create(GOLDEN).expect("couldn't write golden file");
        writeln!(out, "# worst_corpus golden values: <fn> <x_bits> <result_bits>").unwrap();
        writeln!(out, "# {} functions x {} inputs", fns.len(), xs.len()).unwrap();
        for l in &lines {
            writeln!(out, "{l}").unwrap();
        }
        println!("blessed {} entries ({} functions x {} inputs) -> {GOLDEN}", lines.len(), fns.len(), xs.len());
        return;
    }

    let golden = match std::fs::read_to_string(GOLDEN) {
        Ok(s) => s,
        Err(_) => {
            eprintln!("no golden file at {GOLDEN}; run with --bless first");
            std::process::exit(2);
        }
    };
    let expected: Vec<&str> = golden.lines().filter(|l| !l.starts_with('#')).collect();

    if expected.len() != lines.len() {
        println!(
            "corpus size changed: golden has {} entries, this build produces {} \
             (a function or input was added/removed -- re-bless and review)",
            expected.len(),
            lines.len()
        );
        std::process::exit(1);
    }

    // A NaN is a NaN: two NaN results are the same result whatever their
    // payload, sign or quiet bits say, so a payload-only move is not a
    // behavioural change and must not fail the gate. This is the same rule
    // every ulp metric in the repo follows (see accuracy.rs's `ulp_diff`);
    // the gate is bit-exact about *values*, not about NaN bookkeeping.
    // Everything else stays a raw bit comparison.
    fn same_result(a: &str, b: &str) -> bool {
        if a == b {
            return true;
        }
        let (ap, bp): (Vec<&str>, Vec<&str>) = (a.split(' ').collect(), b.split(' ').collect());
        // same function and same input, differing only in a NaN result
        ap.len() == 3
            && bp.len() == 3
            && ap[0] == bp[0]
            && ap[1] == bp[1]
            && [ap[2], bp[2]].iter().all(|h| {
                u32::from_str_radix(h, 16).map(|b| f32::from_bits(b).is_nan()).unwrap_or(false)
            })
    }

    let mut diffs: Vec<String> = Vec::new();
    for (got, want) in lines.iter().zip(expected.iter()) {
        if !same_result(got, want) {
            let gp: Vec<&str> = got.split(' ').collect();
            let wp: Vec<&str> = want.split(' ').collect();
            diffs.push(format!(
                "  {} at x={} : {} -> {}   ({:e} -> {:e})",
                gp[0], gp[1], wp[2], gp[2],
                f32::from_bits(u32::from_str_radix(wp[2], 16).unwrap()),
                f32::from_bits(u32::from_str_radix(gp[2], 16).unwrap()),
            ));
        }
    }

    println!(
        "worst_corpus: {} entries ({} functions x {} inputs)",
        lines.len(),
        fns.len(),
        xs.len()
    );
    if diffs.is_empty() {
        println!("all bit-identical to the golden corpus");
        return;
    }
    println!("{} value(s) changed:", diffs.len());
    for d in diffs.iter().take(40) {
        println!("{d}");
    }
    if diffs.len() > 40 {
        println!("  ... and {} more", diffs.len() - 40);
    }
    println!("\nIf these changes are intentional, re-run with --bless.");
    std::process::exit(1);
}
