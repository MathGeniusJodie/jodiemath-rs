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
#![allow(clippy::approx_constant, clippy::type_complexity)]
use jodiemath_rs::*;

include!("support/unary_fns.rs");

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
        if v.is_sign_negative() {
            "-0".to_string()
        } else {
            "+0".to_string()
        }
    } else if v.is_infinite() {
        if v < 0.0 {
            "-inf".to_string()
        } else {
            "+inf".to_string()
        }
    } else {
        format!("{v:.6e}")
    }
}

fn main() {
    assert_unary_fns_complete();
    let fns = UNARY_FNS;

    println!("{} public 1-arg f32->f32 functions\n", fns.len());
    println!(
        "{:<26} {:>14} {:>14} {:>14} {:>14} {:>14}",
        "function", "f(+0)", "f(-0)", "f(+inf)", "f(-inf)", "f(NaN)"
    );
    println!("{}", "-".repeat(102));

    let mut nan_broken: Vec<&str> = Vec::new();
    let mut signaling: Vec<&str> = Vec::new();
    let mut zero_sign_differs: Vec<&str> = Vec::new();

    for (name, f) in fns {
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
            name,
            cls(p0),
            cls(n0),
            cls(pi),
            cls(ni),
            cls(nn)
        );
    }

    println!("\n=== automatic checks ===");
    println!(
        "NaN not propagated ({}): {:?}",
        nan_broken.len(),
        nan_broken
    );
    println!(
        "signaling NaN returned ({}): {:?}",
        signaling.len(),
        signaling
    );
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
