// Adversarial worst-case search for the powf family (backlog idea #93's
// 2-arg importance sampling, specialised to powf), where a blind 2-arg
// fuzz is structurally weak: powf's error is
//
//     ln2 * |y*log2(x)| * (relative error of the log2 it feeds exp2)
//
// so the max is reached only where *both* factors are extreme at once --
// |y*log2(x)| pinned just inside the largest exponent with a finite
// result, and x inside the one octave where log2's own relative error is
// undiluted by the integer exponent (m in [2^-0.5, 2^0.5), i.e. k == 0;
// outside it |log2(x)| >= 0.5 divides the same absolute error down).
// Random (x, y) pairs hit that corner with vanishing probability -- the
// standard fuzz reports powf at 2 ulp where the real worst case
// is twice that -- so this picks y from x instead of drawing it, which
// makes the search exhaustive over the band rather than statistical.
//
// The reference is f64 powf: ~2^-52 relative, i.e. ~2^-28 of an f32 ulp,
// so it is a sound judge at the scale being measured here.
use jodiemath_rs::*;

fn ulp_diff(a: f32, b: f32) -> u64 {
    fn ord(x: f32) -> i64 {
        let b = x.to_bits();
        if b & 0x8000_0000 != 0 { -((b & 0x7fff_ffff) as i64) } else { b as i64 }
    }
    if a.is_nan() || b.is_nan() {
        return if a.is_nan() == b.is_nan() { 0 } else { u64::MAX };
    }
    (ord(a) - ord(b)).unsigned_abs()
}

/// Sweep `[lo, hi)` by `stride`, choosing `y` per `x` so that
/// `y*log2(x) ~ target`, and score against the f64 reference.
fn sweep(name: &str, f: impl Fn(f32, f32) -> f32, target: f64, lo: u32, hi: u32, stride: u32) -> u64 {
    let (mut max, mut sum, mut n, mut wx, mut wy) = (0u64, 0u64, 0u64, 0f32, 0f32);
    let mut b = lo;
    while b < hi {
        let x = f32::from_bits(b);
        b += stride;
        let l = (x as f64).log2();
        if l == 0.0 || !x.is_normal() {
            continue;
        }
        let y = (target / l) as f32;
        let r = (x as f64).powf(y as f64);
        // skip anything that lands on a denormal/zero/inf, where an ulp
        // count stops being a meaningful relative measure
        if !r.is_normal() || !(r as f32).is_normal() {
            continue;
        }
        let d = ulp_diff(f(x, y), r as f32);
        sum += d;
        n += 1;
        if d > max {
            (max, wx, wy) = (d, x, y);
        }
    }
    println!(
        "{name:>22}  y*log2(x)~{target:>7.1}  avg ulp {:9.4}  max ulp {max:>6}  at x={wx:e} y={wy:e}",
        sum as f64 / n as f64
    );
    max
}

fn main() {
    // the k == 0 octave, exhaustively: every f32 in [2^-0.5, 2^0.5)
    let (k0lo, k0hi) = (0x3f3504f3u32, 0x3fb504f3u32);
    // every normal positive f32, strided -- confirms the k == 0 octave
    // really is the worst case rather than assuming it
    let (alllo, allhi) = (0x0080_0000u32, 0x7f80_0000u32);
    let mut worst = [0u64; 2];
    for target in [127.9f64, 64.0, -64.0, -125.9] {
        for (i, (name, f)) in [
            ("powf", &powf as &dyn Fn(f32, f32) -> f32),
            ("powf_unchecked", &powf_unchecked),
        ]
        .into_iter()
        .enumerate()
        {
            let a = sweep(name, f, target, k0lo, k0hi, 1);
            let b = sweep(&format!("{name} (all x)"), f, target, alllo, allhi, 61);
            worst[i] = worst[i].max(a).max(b);
        }
        println!();
    }
    println!("worst over all sweeps:  powf {}  powf_unchecked {}", worst[0], worst[1]);
}
