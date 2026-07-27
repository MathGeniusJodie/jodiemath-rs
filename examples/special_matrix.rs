// Idea #164: systematic special-value matrix over every public 1-arg
// f32 -> f32 function. Five +-0 bugs (acos, atan2, sinf_poly, sinpi,
// remainder) were found ad hoc before this existed; this makes the sweep
// standing and exhaustive over the function list instead of incidental.
//
// Two classes of check run automatically, needing no per-function
// reference:
//   (1) NaN propagation: f(NaN) must be NaN. A function that turns NaN
//       into a number is silently wrong for every caller downstream.
//   (2) NaN payload quietness: the result must be a *quiet* NaN.
// Everything else is printed as a matrix for review, since the correct
// value at +-0/+-inf is function-specific (and for the _unchecked/_approx
// tiers, deliberately unspecified outside their domains).
#![allow(clippy::approx_constant)]
use jodiemath_rs::*;

// Functions whose own doc comments already promise nothing for non-finite
// input: the `_unchecked` cores (no special-value selects at all, that
// being the point) and the `_approx` bit-trick tiers (pure exponent-field
// arithmetic, so a NaN's exponent field is just read as a large number).
// Listed explicitly rather than pattern-matched on the name so that a new
// function silently inheriting garbage NaN behaviour still fails the gate.
const NAN_EXEMPT: &[&str] = &[
    "exp2_approx",
    "ln_unchecked",
    "log10_unchecked",
    "log2_approx",
    "log_2_unchecked",
    "rcp_approx",
    "rsqrt_approx",
    "sqrt_approx",
];

fn cls(v: f32) -> String {
    if v.is_nan() {
        // quiet NaN has the top mantissa bit set
        let quiet = v.to_bits() & 0x0040_0000 != 0;
        format!("NaN{}", if quiet { "" } else { "(SIGNALING)" })
    } else if v == 0.0 {
        if v.is_sign_negative() { "-0".to_string() } else { "+0".to_string() }
    } else if v.is_infinite() {
        if v < 0.0 { "-inf".to_string() } else { "+inf".to_string() }
    } else {
        format!("{v:.6e}")
    }
}

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
        ("erfc_accurate", erfc_accurate as fn(f32) -> f32),
        ("erfcx", erfcx as fn(f32) -> f32),
        ("erfcx_accurate", erfcx_accurate as fn(f32) -> f32),
        ("erfcx_checked", erfcx_checked as fn(f32) -> f32),
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

    println!("{} public 1-arg f32->f32 functions\n", fns.len());
    println!(
        "{:<26} {:>14} {:>14} {:>14} {:>14} {:>14}",
        "function", "f(+0)", "f(-0)", "f(+inf)", "f(-inf)", "f(NaN)"
    );
    println!("{}", "-".repeat(102));

    let mut nan_broken: Vec<&str> = Vec::new();
    let mut signaling: Vec<&str> = Vec::new();
    let mut zero_sign_differs: Vec<&str> = Vec::new();

    for (name, f) in &fns {
        let p0 = f(0.0);
        let n0 = f(-0.0);
        let pi = f(f32::INFINITY);
        let ni = f(f32::NEG_INFINITY);
        let nn = f(f32::NAN);

        if !nn.is_nan() {
            if !NAN_EXEMPT.contains(name) {
                nan_broken.push(name);
            }
        } else if nn.to_bits() & 0x0040_0000 == 0 {
            signaling.push(name);
        }
        // Record whether +-0 map to different signs -- expected for odd
        // functions, a red flag for even ones. Reviewed, not auto-failed.
        if p0.to_bits() != n0.to_bits() {
            zero_sign_differs.push(name);
        }

        println!(
            "{:<26} {:>14} {:>14} {:>14} {:>14} {:>14}",
            name, cls(p0), cls(n0), cls(pi), cls(ni), cls(nn)
        );
    }

    println!("\n=== automatic checks ===");
    println!("NaN not propagated ({}): {:?}", nan_broken.len(), nan_broken);
    println!("signaling NaN returned ({}): {:?}", signaling.len(), signaling);
    println!(
        "\nf(+0) and f(-0) differ in bits ({}) -- expected for odd fns, review for even ones:\n  {:?}",
        zero_sign_differs.len(),
        zero_sign_differs
    );

    if !nan_broken.is_empty() || !signaling.is_empty() {
        println!("\nFAIL");
        std::process::exit(1);
    }
    println!("\nautomatic checks PASS");
}
