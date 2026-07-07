// Coordinate-descent ULP tuner for polynomial coefficients.
use jodiemath_rs::*;

#[inline(always)]
fn fma(a: f32, b: f32, c: f32) -> f32 {
    a.mul_add(b, c)
}

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

const EXPONENT_MASK: u32 = 0x7f800000;

#[inline(always)]
fn exp2_c(x: f32, c: &[f32]) -> f32 {
    let k = x.floor();
    let f = x - k;
    let exp2int = f32::from_bits(((k + 383_f32).to_bits() << 8) & EXPONENT_MASK);
    fma(
        fma(fma(c[0], f, c[1]), f, c[2]),
        exp2int * (f * f) * (f * f),
        fma(fma(fma(c[3], f, c[4]), f, c[5]), exp2int * f, exp2int),
    )
}

#[inline(always)]
fn log2_c(x: f32, c: &[f32]) -> f32 {
    let e = (x.to_bits() as i32).wrapping_sub(0x3f3504f3) >> 23;
    let m = f32::from_bits((x.to_bits() as i32).wrapping_sub(e << 23) as u32);
    let k = e as f32;
    let s = m - 1.0;
    let s2 = s * s;
    let s4 = s2 * s2;
    let l0 = fma(c[1], s, c[0]);
    let l1 = fma(c[3], s, c[2]);
    let l2 = fma(c[5], s, c[4]);
    let l3 = fma(c[7], s, c[6]);
    let l4 = fma(c[9], s, c[8]);
    let r0 = fma(l1, s2, l0);
    let r1 = fma(l3, s2, l2);
    let r2 = fma(l4, s4, r1);
    let p = fma(r2, s4, r0);
    fma(p, s, k)
}

fn score(
    f: &dyn Fn(f32, &[f32]) -> f32,
    reference: &dyn Fn(f64) -> f64,
    grid: &[f32],
    c: &[f32],
) -> (u64, u64) {
    let mut sum = 0u64;
    let mut max = 0u64;
    for &x in grid {
        let r = reference(x as f64) as f32;
        let d = ulp_diff(f(x, c), r);
        sum += d;
        max = max.max(d);
    }
    (max, sum)
}

#[inline(always)]
fn mulsign_c(x: f32, y: f32) -> f32 {
    f32::from_bits(x.to_bits() ^ (y.to_bits() & 0x8000_0000))
}

// asin's mid-branch rational correction (see src/lib.rs's asin doc
// comment): a2 = (a^2-a)/d + a where d is this 3-coefficient poly (plus
// the leading constant folded in as c[3]), then the already-fixed
// rationalized sqrt step. Only the correction's 4 coefficients are being
// tuned here -- the rationalization and mulsign/branch structure are
// fixed, matching what's actually shipped in src/lib.rs.
#[inline(always)]
fn asin_mid_c(x: f32, c: &[f32]) -> f32 {
    let a = x.abs();
    let d = fma(-c[0], a, c[1]);
    let d = fma(-a, d, c[2]);
    let d = fma(-a, d, c[3]);
    let a2 = (a * a - a) / d + a;
    let sq = (1.0 - a2).sqrt();
    let sm1 = -a2 / (sq + 1.0);
    mulsign_c(sm1, x) * (-std::f32::consts::FRAC_PI_2)
}

fn tune(
    name: &str,
    f: &dyn Fn(f32, &[f32]) -> f32,
    reference: &dyn Fn(f64) -> f64,
    grid: &[f32],
    init: &[f32],
) {
    let mut c: Vec<f32> = init.to_vec();
    let mut best = score(f, reference, grid, &c);
    println!("{name}: start max {} avg {:.5}", best.0, best.1 as f64 / grid.len() as f64);
    let mut improved = true;
    while improved {
        improved = false;
        for i in 0..c.len() {
            for delta in [1i32, -1, 2, -2, 4, -4, 8, -8, 16, -16] {
                let mut trial = c.clone();
                trial[i] = f32::from_bits((trial[i].to_bits() as i32 + delta) as u32);
                let s = score(f, reference, grid, &trial);
                if s < best {
                    best = s;
                    c = trial;
                    improved = true;
                }
            }
        }
    }
    println!(
        "{name}: tuned max {} avg {:.5}  coeffs: {:?}",
        best.0,
        best.1 as f64 / grid.len() as f64,
        c.iter().map(|v| format!("{v:e}")).collect::<Vec<_>>()
    );
}

fn main() {
    let which = std::env::args().nth(1).unwrap_or_default();
    if which.contains("exp2") || which.is_empty() {
        // grid over (-126, 128)
        let mut grid = vec![];
        let mut b = 1e-6f32.to_bits();
        while b <= 126.0f32.to_bits() {
            grid.push(f32::from_bits(b));
            grid.push(-f32::from_bits(b));
            b += 997;
        }
        let init = [2.1702255e-4, 1.2439688e-3, 9.678841e-3, 5.5483342e-2, 2.4022984e-1, 6.9314698e-1];
        tune("exp2", &exp2_c, &|x| x.exp2(), &grid, &init);
    }
    if which.contains("log2") || which.is_empty() {
        let mut grid = vec![];
        let mut b = 0x0080_0000u32;
        while b < 0x7f80_0000 {
            grid.push(f32::from_bits(b));
            b += 1499;
        }
        let init = [
            1.4426950, -0.72134735, 0.48089824, -0.36069664, 0.28856741,
            -0.23961740, 0.20460062, -0.19106275, 0.18617496, -0.10994955,
        ];
        tune("log2", &log2_c, &|x| x.log2(), &grid, &init);
    }
    if which.contains("asin") {
        // asin's mid branch is only ever evaluated for a = |x| in
        // [0.1, 0.9) in the shipped code (see src/lib.rs's asin) -- tune
        // against exactly that range, not the whole [0,1] domain, so the
        // objective matches what's actually on the hot path.
        let mut grid = vec![];
        let mut b = 0.1f32.to_bits();
        while b < 0.9f32.to_bits() {
            grid.push(f32::from_bits(b));
            grid.push(-f32::from_bits(b));
            b += 37;
        }
        let init = [0.0392588, 0.179323, 1.75866, -3.66063];
        tune("asin_mid", &asin_mid_c, &|x| x.asin(), &grid, &init);
    }
}
