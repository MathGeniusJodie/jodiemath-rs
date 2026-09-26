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
#![allow(clippy::approx_constant, clippy::type_complexity)]
use jodiemath_rs::*;

include!("support/unary_fns.rs");

const CANON: u32 = 0x7fc0_0000; // default quiet NaN, empty payload

fn main() {
    assert_unary_fns_complete();
    let fns = UNARY_FNS;

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

    for (name, f) in fns {
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
    println!("{:<26} {:<24} per-input verdicts", "function", "class");
    println!("{}", "-".repeat(100));
    for r in &rows {
        println!("{r}");
    }
    println!(
        "\nsummary: {preserve} preserve payload+sign exactly, {canonical} canonicalize, \
{other} keep the payload but vary the sign, {off_domain} are off-domain tiers returning no NaN at all"
    );
    println!("\nNeither behaviour is an IEEE754 violation -- this is a record, not a gate.");
    println!(
        "`special_matrix.rs` is what asserts the part that IS required: NaN in -> quiet NaN out."
    );
}
