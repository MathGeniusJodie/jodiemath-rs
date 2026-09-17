// Idea #164, second half: the +-0/+-inf/NaN *cross product* for 2-arg
// functions. This is where the original atan2 and remainder +-0 bugs
// lived, and the 1-arg sweep (examples/special_matrix.rs) deliberately
// didn't cover it.
//
// For the functions with an f64 std counterpart, std implements the IEEE
// 754 / C99 Annex F special-value rules, and at these inputs every correct
// result is exactly representable in f32 -- so a bitwise comparison of
// jodie(x,y) against (std_f64(x,y) as f32) is a real check, not a
// precision-tolerance one. Functions without a counterpart are printed
// for review instead.
//
// Known and deliberate convention divergences are allowlisted with the
// reason, so this exits nonzero only on something new.
#![allow(clippy::approx_constant)]
use jodiemath_rs::*;

const VALS: &[(&str, f32)] = &[
    ("-inf", f32::NEG_INFINITY),
    ("-2", -2.0),
    ("-1", -1.0),
    ("-0", -0.0),
    ("+0", 0.0),
    ("1", 1.0),
    ("2", 2.0),
    ("+inf", f32::INFINITY),
    ("NaN", f32::NAN),
];

fn same(a: f32, b: f32) -> bool {
    if a.is_nan() && b.is_nan() {
        return true;
    }
    a.to_bits() == b.to_bits()
}

fn show(v: f32) -> String {
    if v.is_nan() {
        "NaN".into()
    } else if v == 0.0 {
        if v.is_sign_negative() {
            "-0".into()
        } else {
            "+0".into()
        }
    } else if v.is_infinite() {
        if v < 0.0 {
            "-inf".into()
        } else {
            "+inf".into()
        }
    } else {
        format!("{v}")
    }
}

struct Case {
    name: &'static str,
    jodie: fn(f32, f32) -> f32,
    reference: fn(f64, f64) -> f64,
    // (x, y) pairs where a divergence from f64 std is known and intended.
    // Empty for every function as of the first run (2026-07-27) -- all
    // nine match std exactly. Document the reason inline when adding one.
    exempt: &'static [(&'static str, &'static str)],
}

fn main() {
    let cases: Vec<Case> = vec![
        Case {
            name: "atan2",
            jodie: atan2,
            reference: f64::atan2,
            exempt: &[],
        },
        Case {
            name: "powf",
            jodie: powf,
            reference: f64::powf,
            exempt: &[],
        },
        Case {
            name: "fmod",
            jodie: fmod,
            reference: |x, y| x % y,
            exempt: &[],
        },
        Case {
            name: "fmod_checked",
            jodie: fmod_checked,
            reference: |x, y| x % y,
            exempt: &[],
        },
        Case {
            name: "div_euclid",
            jodie: div_euclid,
            reference: f64::div_euclid,
            exempt: &[],
        },
        Case {
            name: "rem_euclid",
            jodie: rem_euclid,
            reference: f64::rem_euclid,
            exempt: &[],
        },
    ];

    let mut total_mismatch = 0usize;
    let mut findings: Vec<String> = Vec::new();

    for c in &cases {
        let mut rows: Vec<String> = Vec::new();
        for &(xn, x) in VALS {
            for &(yn, y) in VALS {
                let got = (c.jodie)(x, y);
                let want64 = (c.reference)(x as f64, y as f64);
                let want = want64 as f32;
                if same(got, want) {
                    continue;
                }
                if c.exempt.contains(&(xn, yn)) {
                    continue;
                }
                rows.push(format!(
                    "    {}({:>4}, {:>4}) = {:>6}   std f64 -> {:>6}",
                    c.name,
                    xn,
                    yn,
                    show(got),
                    show(want)
                ));
            }
        }
        if rows.is_empty() {
            println!("{:<16} all 81 combos match f64 std", c.name);
        } else {
            println!("{:<16} {} mismatch(es):", c.name, rows.len());
            for r in &rows {
                println!("{r}");
            }
            total_mismatch += rows.len();
            findings.push(format!("{} ({})", c.name, rows.len()));
        }
    }

    // Print-only: no std counterpart to assert against.
    println!("\n=== print-only (no f64 std counterpart) ===");
    let printonly: &[(&str, fn(f32, f32) -> f32)] = &[
        ("xlogy", xlogy),
        ("xlog1py", xlog1py),
        ("signed_pow", signed_pow),
        ("mulsign", mulsign),
    ];
    for (name, f) in printonly {
        // Only show the rows involving a zero, an infinity or a NaN in
        // either slot -- that's the whole point here, and printing all 81
        // per function would bury it.
        let mut shown = 0;
        let mut buf: Vec<String> = Vec::new();
        for &(xn, x) in VALS {
            for &(yn, y) in VALS {
                let interesting = x == 0.0 || y == 0.0 || !x.is_finite() || !y.is_finite();
                if !interesting {
                    continue;
                }
                buf.push(format!("{}({},{})={}", name, xn, yn, show(f(x, y))));
                shown += 1;
            }
        }
        println!("  {name} ({shown} special combos):");
        for chunk in buf.chunks(5) {
            println!("      {}", chunk.join("  "));
        }
    }

    println!("\n=== summary ===");
    if total_mismatch == 0 {
        println!("all asserted functions match f64 std on the full 9x9 cross product");
    } else {
        println!("{total_mismatch} mismatch(es) across: {findings:?}");
        std::process::exit(1);
    }
}
