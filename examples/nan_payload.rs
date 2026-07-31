// Idea #175: NaN payload / quietness propagation. `special_matrix.rs`
// already asserts that every function returns a *quiet* NaN for a NaN input;
// this asks the finer question it doesn't: does the input NaN's *payload*
// (the low 22 mantissa bits) and sign survive, or does the function
// canonicalize to the default NaN?
//
// IEEE754 permits either -- a function may return one of its input NaNs
// unchanged, or a fresh canonical quiet NaN -- so neither answer is a bug.
// This is documentation-grade: the point is to record which behaviour each
// function actually has, since callers who smuggle diagnostic tags in NaN
// payloads need to know where they survive.
//
// Mechanism note for reading the results: payloads survive an operation that
// merely *propagates* an operand (a multiply, an fma, a select), and are lost
// wherever hardware produces a NaN from scratch -- e.g. `inf * 0`, `0/0`, or
// a comparison-driven select that picks a computed constant instead of x.
#![allow(clippy::approx_constant)]
use jodiemath_rs::*;

const CANON: u32 = 0x7fc0_0000; // default quiet NaN, empty payload

fn main() {
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

    // distinctive payloads, plus a negative-signed NaN
    let inputs: [(&str, u32); 4] = [
        ("qNaN payload 0x2a", 0x7fc0_002a),
        ("qNaN payload 0x155555", 0x7fd5_5555),
        ("-qNaN payload 0x2a", 0xffc0_002a),
        ("qNaN canonical", CANON),
    ];

    let mut preserve = 0usize;
    let mut canonical = 0usize;
    let mut other = 0usize;
    let mut off_domain = 0usize;
    let mut rows: Vec<String> = Vec::new();

    for (name, f) in &fns {
        let mut verdicts: Vec<&str> = Vec::new();
        for &(_, bits) in &inputs {
            let x = f32::from_bits(bits);
            let y = f(x);
            if !y.is_nan() {
                verdicts.push("NOT-NAN");
            } else if y.to_bits() == bits {
                verdicts.push("exact");
            } else if y.to_bits() & 0x003f_ffff == bits & 0x003f_ffff {
                verdicts.push("payload-kept");
            } else if y.to_bits() & 0x7fff_ffff == CANON {
                verdicts.push("canonical");
            } else {
                verdicts.push("other");
            }
        }
        // classify the function by its behaviour on the tagged payloads
        let all_exact = verdicts.iter().take(3).all(|v| *v == "exact");
        let all_canon = verdicts.iter().take(3).all(|v| *v == "canonical");
        let no_nan = verdicts.iter().all(|v| *v == "NOT-NAN");
        let class = if no_nan {
            // The `_unchecked`/`_approx` tiers: their docs promise nothing
            // for non-finite input, so they return a number rather than any
            // NaN at all. Not a payload question -- broken out so it doesn't
            // pollute the "mixed" count. (special_matrix.rs exempts exactly
            // this set by name for the same reason.)
            off_domain += 1;
            "off-domain tier (no NaN)"
        } else if all_exact {
            preserve += 1;
            "preserves payload+sign"
        } else if all_canon {
            canonical += 1;
            "canonicalizes"
        } else {
            other += 1;
            "payload kept, sign varies"
        };
        rows.push(format!("{name:<26} {class:<24} {verdicts:?}"));
    }

    println!("{} public 1-arg functions\n", fns.len());
    println!("{:<26} {:<24} {}", "function", "class", "per-input verdicts");
    println!("{}", "-".repeat(100));
    for r in &rows {
        println!("{r}");
    }
    println!(
        "\nsummary: {preserve} preserve payload+sign exactly, {canonical} canonicalize, \
{other} keep the payload but vary the sign, {off_domain} are off-domain tiers returning no NaN at all"
    );
    println!("\nNeither behaviour is an IEEE754 violation -- this is a record, not a gate.");
    println!("`special_matrix.rs` is what asserts the part that IS required: NaN in -> quiet NaN out.");
}
