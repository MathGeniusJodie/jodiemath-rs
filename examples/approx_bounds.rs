// Verifies the _approx tier's documented error bounds (idea #188 gave
// exp2_approx/log2_approx/rsqrt_approx real doc-comment bounds, but nothing
// checked them), and measures the three that carry no doc comment at all
// so their bounds can be written down too.
//
// These functions are deliberately outside the crate's 0.5/2 ulp budget --
// they are bit-trick seeds, not accurate kernels -- so the metric is the
// one each doc actually claims (max *relative* error, or max *absolute*
// error where the true value crosses zero), not ulp.
//
// Positive normals are swept on a stride rather than exhaustively: these
// are smooth per-octave bit-trick approximations, so a stride of 16 over
// 2^31 patterns (~134M samples) locates the max reliably while keeping the
// run to seconds. cbrt_approx's own comment notes its error does not repeat
// across octaves, which the stride still covers since it visits every
// octave densely.
use jodiemath_rs::*;

const STRIDE: u32 = 16;

fn sweep_positive_normals(
    name: &str,
    f: impl Fn(f32) -> f32,
    r: impl Fn(f64) -> f64,
    relative: bool,
) -> (f64, f32) {
    let lo = f32::MIN_POSITIVE.to_bits(); // first positive normal
    let hi = f32::MAX.to_bits();
    let mut worst = 0.0f64;
    let mut worst_x = 0.0f32;
    let mut bits = lo;
    while bits <= hi {
        let x = f32::from_bits(bits);
        let got = f(x) as f64;
        let want = r(x as f64);
        if want.is_finite() && got.is_finite() {
            let e = if relative {
                if want != 0.0 { ((got - want) / want).abs() } else { 0.0 }
            } else {
                (got - want).abs()
            };
            if e > worst {
                worst = e;
                worst_x = x;
            }
        }
        bits = bits.wrapping_add(STRIDE);
    }
    let unit = if relative { "relative" } else { "absolute" };
    println!("  {name:<14} max {unit} error {worst:.6}  at x = {worst_x:e}");
    (worst, worst_x)
}

fn main() {
    println!("=== documented bounds (idea #188) ===");

    // exp2_approx: doc claims max relative error ~6.1% over |x| < 120.
    // Domain is a value range, not a bit range, so sample x directly.
    let mut e_worst = 0.0f64;
    let mut e_worst_x = 0.0f32;
    const N: i64 = 40_000_000;
    for i in 0..=N {
        let x = (-120.0 + 240.0 * (i as f64) / (N as f64)) as f32;
        let got = exp2_approx(x) as f64;
        let want = (x as f64).exp2();
        if want != 0.0 && want.is_finite() && got.is_finite() {
            let e = ((got - want) / want).abs();
            if e > e_worst {
                e_worst = e;
                e_worst_x = x;
            }
        }
    }
    println!(
        "  exp2_approx    max relative error {:.6} ({:.3}%)  at x = {:e}   [doc: ~6.1% over |x|<120]",
        e_worst,
        e_worst * 100.0,
        e_worst_x
    );
    let exp2_ok = e_worst <= 0.062;

    // log2_approx: doc claims max absolute error ~0.086 over positive normals.
    let (l2, _) = sweep_positive_normals("log2_approx", log2_approx, f64::log2, false);
    println!("                 [doc: max absolute error ~0.086 over positive normal x]");
    let log2_ok = l2 <= 0.09;

    // rsqrt_approx: doc claims max relative error ~4.8% over positive normals.
    let (rs, _) = sweep_positive_normals(
        "rsqrt_approx",
        rsqrt_approx,
        |v| 1.0 / v.sqrt(),
        true,
    );
    println!("                 [doc: max relative error ~4.8% over positive normal x]");
    let rsqrt_ok = rs <= 0.049;

    println!("\n=== undocumented _approx members (no doc comment at all) ===");
    println!("  full positive-normal range:");
    sweep_positive_normals("sqrt_approx", sqrt_approx, f64::sqrt, true);
    sweep_positive_normals("rcp_approx", rcp_approx, |v| 1.0 / v, true);
    sweep_positive_normals("cbrt_approx", cbrt_approx, f64::cbrt, true);

    // The full-range numbers above are dominated by the extreme octaves,
    // where these bit tricks have no range handling at all (rcp_approx's
    // reciprocal underflows toward denormals, cbrt_approx's seed addition
    // leaves the normal exponent field). Re-measure over a mid-range domain
    // so a usable bound can actually be written down for each.
    println!("\n  restricted to |x| in [1e-30, 1e30] (mid-range, where a bound is meaningful):");
    for (name, f, r) in [
        ("sqrt_approx", sqrt_approx as fn(f32) -> f32, f64::sqrt as fn(f64) -> f64),
        ("rcp_approx", rcp_approx, |v: f64| 1.0 / v),
        ("cbrt_approx", cbrt_approx, f64::cbrt),
    ] {
        let mut worst = 0.0f64;
        let mut worst_x = 0.0f32;
        let mut bits = 1e-30f32.to_bits();
        let hi = 1e30f32.to_bits();
        while bits <= hi {
            let x = f32::from_bits(bits);
            let got = f(x) as f64;
            let want = r(x as f64);
            if want != 0.0 && want.is_finite() && got.is_finite() {
                let e = ((got - want) / want).abs();
                if e > worst {
                    worst = e;
                    worst_x = x;
                }
            }
            bits = bits.wrapping_add(STRIDE);
        }
        println!(
            "  {name:<14} max relative error {worst:.6} ({:.3}%)  at x = {worst_x:e}",
            worst * 100.0
        );
    }

    println!("\n=== verdict ===");
    for (name, ok, claim) in [
        ("exp2_approx", exp2_ok, "~6.1% relative"),
        ("log2_approx", log2_ok, "~0.086 absolute"),
        ("rsqrt_approx", rsqrt_ok, "~4.8% relative"),
    ] {
        println!("  {name:<14} documented {claim:<18} {}", if ok { "HOLDS" } else { "VIOLATED" });
    }
    if !(exp2_ok && log2_ok && rsqrt_ok) {
        println!("\nFAIL: a documented bound does not hold");
        std::process::exit(1);
    }
    println!("\nall documented _approx bounds hold");
}
