// Latency and throughput of the radian trig family with the input exponent
// unknown to the compiler. quickbench's `Band::mix` pins the exponent, which
// lets LLVM constant-fold everything exponent-driven (for the wide tier, the
// whole 1/pi window), so its latency column flatters them.
//
//   cargo run --release --example trig_bench [filter]
//
// Latency: a dependency chain whose every step ORs a pseudo-random exponent
// (loaded from a black-boxed table, off the chain) into the previous result.
// Throughput: independent evaluations over inputs with random exponents in
// the same range.
use jodiemath_rs::*;
use std::hint::black_box;
use std::time::Instant;

const REPS: usize = 7;
const LAT_ITERS: usize = 1 << 22;
const TP_LEN: usize = 4096;
const TP_PASSES: usize = 512;

/// Exponent fields for each tier's domain: all finite magnitudes for the wide
/// tier, |x| < 2^22 for the narrow one, both from 2^-4 up.
fn exponents(max_biased: u32, seed: u64) -> Vec<u32> {
    let mut s = seed;
    (0..256)
        .map(|_| {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let lo = 123u32;
            (lo + ((s >> 33) as u32) % (max_biased - lo)) << 23
        })
        .collect()
}

fn bench(name: &str, max_biased: u32, f: impl Fn(f32) -> f32) {
    let exps = black_box(exponents(max_biased, 7));
    let mut lat = f64::INFINITY;
    for _ in 0..REPS {
        let mut x = f32::from_bits(0x3fc0_0000);
        let t = Instant::now();
        for i in 0..LAT_ITERS {
            x = f32::from_bits((f(x).to_bits() & 0x807f_ffff) | exps[i & 255]);
        }
        black_box(x);
        lat = lat.min(t.elapsed().as_nanos() as f64 / LAT_ITERS as f64);
    }
    let mut s = 99u64;
    let mut input = [0f32; TP_LEN];
    for v in input.iter_mut() {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        *v = f32::from_bits(((s >> 32) as u32 & 0x807f_ffff) | exps[(s >> 24) as usize & 255]);
    }
    let input = black_box(input);
    let mut out = [0f32; TP_LEN];
    let mut thr = f64::INFINITY;
    for _ in 0..REPS {
        let t = Instant::now();
        for _ in 0..TP_PASSES {
            for (o, &x) in out.iter_mut().zip(input.iter()) {
                *o = f(x);
            }
            black_box(&mut out);
        }
        thr = thr.min(t.elapsed().as_nanos() as f64 / (TP_LEN * TP_PASSES) as f64);
    }
    println!("{name:10} latency {lat:6.2} ns   throughput {thr:6.3} ns/op");
}

fn main() {
    let filter = std::env::args().nth(1).unwrap_or_default();
    macro_rules! run {
        ($name:literal, $max:expr, $f:expr) => {
            if $name.contains(filter.as_str()) {
                bench($name, $max, $f);
            }
        };
    }
    run!("sin_wide", 255, sin_wide);
    run!("cos_wide", 255, cos_wide);
    run!("tan_wide", 255, tan_wide);
    run!("sin", 149, sin);
    run!("cos", 149, cos);
    run!("tan", 149, tan);
}
