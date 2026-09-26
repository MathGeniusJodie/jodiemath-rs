// Exhaustive (or strided) ULP sweep of the radian trig family against glibc's
// f64 sin/cos/tan, which reduce huge arguments exactly, rounded to f32.
//
//   cargo run --release --example trig_sweep [stride] [filter]
//
// stride 1 (default) scores all 2^32 bit patterns; stride 97 is a ~1s screen.
// Scoring follows accuracy.rs's `ulp_diff` (NaN vs NaN is 0).
use jodiemath_rs::*;
use std::sync::atomic::{AtomicU64, Ordering};

const BLOCK: usize = 4096;

struct Case {
    name: &'static str,
    f: fn(&[f32; BLOCK], &mut [f32; BLOCK]),
    reference: fn(f64) -> f64,
    limit: f32,
}

macro_rules! case {
    ($name:literal, $f:expr, $r:expr, $limit:expr) => {
        Case {
            name: $name,
            f: |i, o| {
                for (o, &x) in o.iter_mut().zip(i.iter()) {
                    *o = $f(x);
                }
            },
            reference: $r,
            limit: $limit,
        }
    };
}

/// sin(pi x) or cos(pi x) from the exact split x = n + r, |r| <= 1/2.
fn half_turns(x: f64, cos: bool) -> f64 {
    let n = x.round_ties_even();
    let r = std::f64::consts::PI * (x - n);
    // cos is exactly zero at the half-integers (tanpi's convention there is -inf).
    let v = if cos && (x - n).abs() == 0.5 { 0.0 } else if cos { r.cos() } else { r.sin() };
    if n.rem_euclid(2.0) == 1.0 {
        -v
    } else {
        v
    }
}

fn ulp_diff(a: f32, b: f32) -> u64 {
    fn ord(x: f32) -> i64 {
        let b = x.to_bits();
        if b & 0x8000_0000 != 0 {
            -((b & 0x7fff_ffff) as i64)
        } else {
            b as i64
        }
    }
    if a.is_nan() || b.is_nan() {
        return if a.is_nan() == b.is_nan() { 0 } else { u64::MAX };
    }
    (ord(a) - ord(b)).unsigned_abs()
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let stride: u64 = args.get(1).map_or(1, |s| s.parse().expect("stride"));
    let filter = args.get(2).cloned().unwrap_or_default();
    let pi = std::f32::consts::PI;
    let cases = [
        case!("sin_wide", sin_wide, f64::sin, f32::INFINITY),
        case!("cos_wide", cos_wide, f64::cos, f32::INFINITY),
        case!("tan_wide", tan_wide, f64::tan, f32::INFINITY),
        case!("sin", sin, f64::sin, 16777216.0 * pi),
        case!("cos", cos, f64::cos, 4194304.0 * pi),
        case!("tan", tan, f64::tan, 8388608.0 * pi),
        case!("tan22", tan, f64::tan, 4194304.0 * pi),
        case!("sinpi", sinpi, |x: f64| half_turns(x, false), f32::INFINITY),
        case!("cospi", cospi, |x: f64| half_turns(x, true), f32::INFINITY),
        case!("tanpi", tanpi, |x: f64| {
            let c = half_turns(x, true);
            if c == 0.0 { f64::NEG_INFINITY } else { half_turns(x, false) / c }
        }, f32::INFINITY),
    ];
    for (name, f) in [("sin", sin as fn(f32) -> f32), ("sin_wide", sin_wide), ("tan", tan), ("tan_wide", tan_wide), ("sinpi", sinpi), ("tanpi", tanpi)] {
        for z in [0.0f32, -0.0] {
            assert_eq!(f(z).to_bits(), z.to_bits(), "{name}({z}) lost the sign of zero");
        }
    }
    let threads = std::thread::available_parallelism().map_or(8, |n| n.get());
    let total = (1u64 << 32).div_ceil(stride);
    for case in cases.iter().filter(|c| c.name.contains(&filter)) {
        let start = std::time::Instant::now();
        let next = AtomicU64::new(0);
        let results: Vec<(u64, u64, u64, u32)> = std::thread::scope(|s| {
            let handles: Vec<_> = (0..threads)
                .map(|_| {
                    s.spawn(|| {
                        let (mut sum, mut n, mut max, mut worst) = (0u64, 0u64, 0u64, 0u32);
                        let mut input = [0f32; BLOCK];
                        let mut output = [0f32; BLOCK];
                        loop {
                            let base = next.fetch_add(BLOCK as u64, Ordering::Relaxed);
                            if base >= total {
                                break;
                            }
                            let len = (total - base).min(BLOCK as u64) as usize;
                            for (k, v) in input.iter_mut().enumerate() {
                                let idx = (base + (k as u64).min(len as u64 - 1)) * stride;
                                *v = f32::from_bits(idx as u32);
                            }
                            (case.f)(&input, &mut output);
                            for k in 0..len {
                                let x = input[k];
                                if !(x.abs() < case.limit) && x.is_finite() {
                                    continue;
                                }
                                let r = (case.reference)(x as f64) as f32;
                                let d = ulp_diff(output[k], r);
                                sum += d.min(1 << 32);
                                n += 1;
                                if d > max {
                                    max = d;
                                    worst = x.to_bits();
                                }
                            }
                        }
                        (sum, n, max, worst)
                    })
                })
                .collect();
            handles.into_iter().map(|h| h.join().unwrap()).collect()
        });
        let sum: u64 = results.iter().map(|r| r.0).sum();
        let n: u64 = results.iter().map(|r| r.1).sum();
        let (max, worst) = results.iter().map(|r| (r.2, r.3)).max().unwrap();
        println!(
            "{:<10} avg {:.4}  max {:>3}  worst x={:e} ({:#010x})  n={}  {:.1}s",
            case.name,
            sum as f64 / n as f64,
            max,
            f32::from_bits(worst),
            worst,
            n,
            start.elapsed().as_secs_f64()
        );
    }
}
