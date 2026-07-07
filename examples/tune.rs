// Coordinate-descent ULP tuner for polynomial coefficients.
use jodiemath_rs::*;
// Scalar (non-SIMD, no portable_simd/nightly needed) 1.0-ulp reference,
// same crate accuracy.rs uses for its vectorized ground truth.
use sleef::f64::erf_u10 as erf_ref;
use sleef::f64::erfc_u15 as erfc_ref;

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

// atan_poly (see src/lib.rs): a Pade form approximating atan(x) directly
// for x in [0,1] -- atan() calls this on a.min(1/a), always in that
// range, so that's exactly the grid to tune against.
#[inline(always)]
fn atan_poly_c(x: f32, c: &[f32]) -> f32 {
    let x2 = x * x;
    (fma(fma(c[0], x2, c[1]), x2, 1.0) * x) / fma(fma(x2, c[2], c[3]), x2, 1.0)
}

// acos_poly (see src/lib.rs), scored as the *whole* acos(x) formula for
// x >= 0 (sqrt(1-x)*acos_poly(x)) rather than the bare poly value -- the
// sqrt factor's own rounding interacts with the poly, so tuning the poly
// in isolation could mistune it relative to what actually ships.
// acos_poly is also reused by asin's near-1 branch, so an improvement
// here benefits both callers.
#[inline(always)]
fn acos_poly_c(x: f32, c: &[f32]) -> f32 {
    let u = fma(c[0], x, c[1]);
    let u = fma(u, x, c[2]);
    let u = fma(u, x, c[3]);
    let u = fma(u, x, c[4]);
    let u = fma(u, x, c[5]);
    let poly = fma(u, x, c[6]);
    (1.0 - x).sqrt() * poly
}

// erf's tail branch (see src/lib.rs's erf): scored as the whole
// mulsign(1.0 - exp2(erf_poly(xa)), x) formula, xa in [0.28, 10] (exactly
// where this branch is used in the shipped code; the Pade near-zero
// branch below 0.28 is untouched). Uses f32::exp2 (std) as a stand-in for
// the shipped exp2_checked -- within this bounded range erf_poly never
// leaves exp2's safe domain, so the checked/unchecked distinction doesn't
// matter for tuning purposes.
#[inline(always)]
fn erf_tail_c(x: f32, c: &[f32]) -> f32 {
    let xa = x.abs().min(10.0);
    let u = fma(c[0], xa, c[1]);
    let u = fma(u, xa, c[2]);
    let u = fma(u, xa, c[3]);
    let u = fma(u, xa, c[4]);
    let u = fma(u, xa, c[5]);
    let poly = fma(u, xa, c[6]);
    mulsign_c(1.0 - poly.exp2(), x)
}

// erf's near-zero Pade branch (see src/lib.rs's erf), used for |x| < 0.28
// -- that's exactly the grid to tune against; the tail branch above 0.28
// is untouched.
#[inline(always)]
fn erf_near0_c(x: f32, c: &[f32]) -> f32 {
    let x2 = x * x;
    let numer = x * fma(c[0], x2, c[1]);
    let denom = fma(fma(c[2], x2, c[3]), x2, 1.0);
    numer / denom
}

// cbrt_normal (see src/lib.rs): the degree-3 correction poly's 4
// coefficients, against the *existing* ax/3-based seed (a separate,
// bigger question -- swapping to the cheaper (bits>>16)*0x5556 seed --
// is deferred, see IDEAS.md's cbrt seed entry).
#[inline(always)]
fn cbrt_normal_c(x: f32, c: &[f32]) -> f32 {
    const SIGN_MASK: u32 = 0x8000_0000;
    let ax = x.to_bits() & !SIGN_MASK;
    let a = f32::from_bits(ax);
    let rcp = 1.0 / a;
    let s = f32::from_bits(ax / 3 + 0x2a509a07u32);
    let s2 = s * s;
    let d = fma(s2, s, -a);
    let r = d * rcp;
    let r2 = r * r;
    let a1 = fma(c[1], r, c[0]);
    let b1 = fma(c[3], r, c[2]);
    let p = fma(b1, r2, a1);
    let ss = f32::from_bits(s.to_bits() | (x.to_bits() & SIGN_MASK));
    let sr = ss * r;
    fma(sr, p, ss)
}

// erfc's rational*gaussian tail (see src/lib.rs's erfc): the 8 named
// coefficients (4 for n, 4 for d) are tuned; the two Horner chains'
// trailing "+1.0" leading terms are left fixed, matching the shipped
// structure exactly (not changing the algebraic shape, only retuning
// what's already parameterized).
#[inline(always)]
fn erfc_c(x: f32, c: &[f32]) -> f32 {
    let z = if x < 0.0 { -1.0 } else { 1.0 };
    let w = if x < 0.0 { 2.0 } else { 0.0 };
    let xa = x.abs().min(10.0);
    let n = fma(c[0], xa, c[1]);
    let n = fma(n, xa, c[2]);
    let n = fma(n, xa, c[3]);
    let n = fma(n, xa, 1.0);
    let d = fma(c[4], xa, c[5]);
    let d = fma(d, xa, c[6]);
    let d = fma(d, xa, c[7]);
    let d = fma(d, xa, 1.0);
    let y = (-(xa * xa) * std::f32::consts::LOG2_E).exp2() * n / d;
    fma(y, z, w)
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
    if which.contains("atan") {
        // atan_poly is only ever called on a.min(1/a), i.e. x in [0,1].
        let mut grid = vec![];
        let mut b = 1e-7f32.to_bits();
        while b <= 1.0f32.to_bits() {
            grid.push(f32::from_bits(b));
            b += 37;
        }
        let init = [0.040634338, 0.65748954, 0.17133473, 0.9907859];
        tune("atan_poly", &atan_poly_c, &|x| x.atan(), &grid, &init);
    }
    if which.contains("acos") {
        // acos/asin's near-1 branch both evaluate this for x = |input| in
        // [0,1] (acos directly; asin only above a > 0.9, but tuning the
        // whole [0,1] range keeps the poly consistent for both callers).
        let mut grid = vec![];
        let mut b = 0.0f32.to_bits();
        while b < 1.0f32.to_bits() {
            grid.push(f32::from_bits(b));
            b += 10000;
        }
        let init = [2.2960134e-3, -1.1146357e-2, 2.6900099e-2, -4.8802612e-2, 8.875567e-2, -2.1458527e-1, 1.5707962];
        tune("acos_poly", &acos_poly_c, &|x| x.acos(), &grid, &init);
    }
    if which.contains("erf") {
        // erf's tail branch is only ever used for xa in [0.28, 10] (see
        // erf's doc comment).
        let mut grid = vec![];
        let mut b = 0.28f32.to_bits();
        while b < 10.0f32.to_bits() {
            grid.push(f32::from_bits(b));
            grid.push(-f32::from_bits(b));
            b += 3000;
        }
        let init = [3.118769e-4, -4.67225e-3, 3.3162573e-2, -1.5214339e-1, -9.1684705e-1, -1.6282598, 3.1332566e-5];
        tune("erf_tail", &erf_tail_c, &erf_ref, &grid, &init);

        // erf's near-zero Pade branch, |x| < 0.28.
        let mut grid = vec![];
        let mut b = 0f32.to_bits();
        while b < 0.28f32.to_bits() {
            grid.push(f32::from_bits(b));
            grid.push(-f32::from_bits(b));
            b += 2000;
        }
        let init = [0.5910557508468628, 1.128379225730896, 0.18571428954601288, 0.8571428656578064];
        tune("erf_near0", &erf_near0_c, &erf_ref, &grid, &init);
    }
    if which == "erfc" {
        // erfc's whole domain is xa in [0, 10] (clamped inside the
        // function itself).
        let mut grid = vec![];
        let mut b = 0.0f32.to_bits();
        while b < 10.0f32.to_bits() {
            grid.push(f32::from_bits(b));
            grid.push(-f32::from_bits(b));
            b += 3000;
        }
        let init = [
            1.461691795157094e-6, 0.08557674288749695, 0.44371211528778076, 0.9783496856689453,
            0.15177123248577118, 0.7851238250732422, 1.8210692405700684, 2.1067135334014893,
        ];
        tune("erfc", &erfc_c, &erfc_ref, &grid, &init);
    }
    if which.contains("cbrt") {
        // one octave [1,2) is representative: the bit-trick seed's
        // relative error pattern repeats across octaves (see
        // cbrt_normal's doc comment -- fitted against the seed's error
        // range, not a specific magnitude range).
        let mut grid = vec![];
        let mut b = 1.0f32.to_bits();
        while b < 2.0f32.to_bits() {
            grid.push(f32::from_bits(b));
            grid.push(-f32::from_bits(b));
            b += 5;
        }
        let init = [-0.33333164, 0.22220786, -0.17394418, 0.1482371];
        tune("cbrt_normal", &cbrt_normal_c, &|x| x.cbrt(), &grid, &init);
    }
}
