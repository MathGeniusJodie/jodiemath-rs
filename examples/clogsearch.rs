//! Adversarial worst case for `clog`'s real part, where the blind 2-arg
//! fuzz in `accuracy.rs` is structurally weak in *both* of its halves.
//!
//! `Re clog(z) = ln|z|`, and its hard region is the unit circle: there
//! `ln|z| -> 0` while `re^2 + im^2 -> 1`, so the answer is a near-total
//! cancellation and its *relative* error is whatever relative error
//! `v = re^2 + im^2 - 1` carries. Two separate things hide that from the
//! standing sweep:
//!
//! 1. **Sampling.** `accuracy.rs` draws `re` and `im` as independent
//!    uniform f32 bit patterns, which lands within a few ulp of `|z| = 1`
//!    essentially never. This walks the manifold on purpose: pick `re`,
//!    solve for `im`, then step `im` by a few ulp.
//! 2. **The reference.** `accuracy.rs` uses `(re as f64).hypot(im as
//!    f64).ln()`. `hypot` is correctly rounded, so `|z|` carries a
//!    relative `2^-53` -- but `ln` of it is `~|z|-1`, so that becomes a
//!    relative `2^-53 / |v|`, i.e. ~2 f32 ulp once `|v| ~ 1e-9` and
//!    unbounded below that. The reference cannot score the region it is
//!    being asked about. Here `re^2 + im^2` is instead built *exactly* as
//!    an f64 pair (both squares are exact in f64, and `two_sum` catches
//!    the one rounding in their sum), so `v` is exact and `log1p` is the
//!    only rounding in the whole reference.
//!
//! Exits nonzero if the measured max exceeds `MAX_ULP_BUDGET`, so this can
//! be run as a standing gate alongside `edgecheck`/`worst_corpus`.

use jodiemath_rs::clog;

/// Budget for `Re clog` on the unit-circle manifold.
const MAX_ULP_BUDGET: u64 = 4;

/// Error-free sum of two f64s (Knuth), used only by the reference.
fn two_sum_f64(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let v = s - a;
    (s, (a - (s - v)) + (b - v))
}

/// `ln|z|` to ~`2^-52` relative even when `|z| -> 1`. `re*re` and `im*im`
/// are each exact in f64 (24 bits squared is 48), so `s + es` is exactly
/// `re^2 + im^2`, and `s - 1.0` is Sterbenz-exact for every `s` in
/// `[0.5, 2]` -- which is every `s` this search generates.
fn ln_abs_ref(re: f32, im: f32) -> f64 {
    let (re, im) = (re as f64, im as f64);
    let (s, es) = two_sum_f64(re * re, im * im);
    0.5 * ((s - 1.0) + es).ln_1p()
}

fn ulp_diff(got: f32, want: f32) -> u64 {
    if got.is_nan() && want.is_nan() {
        return 0;
    }
    if got == want {
        return 0;
    }
    if !got.is_finite() || !want.is_finite() {
        return u64::MAX;
    }
    let key = |v: f32| {
        let b = v.to_bits() as i64;
        if b < 0 { i64::MIN.wrapping_sub(b).wrapping_neg() ^ (1i64 << 63) } else { b }
    };
    let (a, b) = (key(got), key(want));
    a.abs_diff(b)
}

fn main() {
    // Walk the manifold: for each `re`, the `im` that puts `|z|` closest to
    // 1, then a window of neighbouring f32s around it. The window is what
    // sweeps `|v|` down through zero -- at offset 0 the residual is
    // whatever f32 rounding leaves, and each step moves it by ~1 ulp of
    // `im^2`.
    let mut worst = 0u64;
    let mut worst_at = (0.0f32, 0.0f32);
    let mut worst_v = 0.0f64;
    let mut buckets = [(0u64, 0usize); 5]; // by |v| decade
    let mut n = 0usize;

    let steps = 200_000u32;
    for i in 0..steps {
        // re in (0, 1/sqrt(2)]: past that the roles swap and the same
        // points are revisited with re/im exchanged.
        let re = (i as f32 + 0.5) / steps as f32 * std::f32::consts::FRAC_1_SQRT_2;
        if re <= 0.0 {
            continue;
        }
        let im0 = (1.0f64 - re as f64 * re as f64).sqrt() as f32;
        for k in -4i32..=4 {
            let im = if k == 0 {
                im0
            } else {
                let mut v = im0;
                for _ in 0..k.abs() {
                    v = if k > 0 { f32::from_bits(v.to_bits() + 1) } else { f32::from_bits(v.to_bits() - 1) };
                }
                v
            };
            if im <= 0.0 || !im.is_finite() {
                continue;
            }
            let want = ln_abs_ref(re, im);
            if want == 0.0 || !want.is_finite() {
                continue;
            }
            let (got_re, _) = clog(re, im);
            let d = ulp_diff(got_re, want as f32);
            n += 1;
            let (s, es) = two_sum_f64(re as f64 * re as f64, im as f64 * im as f64);
            let v = ((s - 1.0) + es).abs();
            let b = match v {
                v if v >= 1e-6 => 0,
                v if v >= 1e-7 => 1,
                v if v >= 1e-8 => 2,
                v if v >= 1e-9 => 3,
                _ => 4,
            };
            if d > buckets[b].0 {
                buckets[b].0 = d;
            }
            buckets[b].1 += 1;
            if d > worst {
                worst = d;
                worst_at = (re, im);
                worst_v = v;
            }
        }
    }

    println!("clog Re, unit-circle manifold ({n} samples)");
    let names = ["|v| >= 1e-6", "1e-7..1e-6", "1e-8..1e-7", "1e-9..1e-8", "|v| < 1e-9"];
    for (b, name) in names.iter().enumerate() {
        println!("  {name:<14} max ulp {:>12}   ({} samples)", buckets[b].0, buckets[b].1);
    }
    println!("  worst {worst} ulp at re={:e} im={:e}  (|v| = {worst_v:e})", worst_at.0, worst_at.1);

    if worst > MAX_ULP_BUDGET {
        eprintln!("FAIL: Re clog max {worst} ulp exceeds budget {MAX_ULP_BUDGET}");
        std::process::exit(1);
    }
    println!("PASS (budget {MAX_ULP_BUDGET})");
}
