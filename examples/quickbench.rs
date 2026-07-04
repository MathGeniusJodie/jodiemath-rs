// Quick latency + throughput bench. Much faster to iterate than criterion.
// Latency: serial dependency chain (output feeds next input).
// Throughput: independent evaluations over an array.
use jodiemath_rs::*;
use std::hint::black_box;
use std::time::Instant;

include!("support/mca_common.rs");

const LAT_ITERS: u64 = 4_000_000;
const TP_ARR: usize = 4096;
const TP_PASSES: usize = 1024;
const REPS: usize = 7;

fn bench_latency(name: &str, f: impl Fn(f32) -> f32) {
    let mut best = f64::INFINITY;
    for _ in 0..REPS {
        // dependency chain; keep value in a sane domain by mixing back toward ~2.0
        let mut x = 1.234_f32;
        let start = Instant::now();
        for _ in 0..LAT_ITERS {
            x = mix(f(x));
        }
        black_box(x);
        let ns = start.elapsed().as_nanos() as f64 / LAT_ITERS as f64;
        best = best.min(ns);
    }
    println!("{:22} latency    {:6.2} ns/op (min of {REPS}, incl chain overhead)", name, best);
}

fn bench_throughput(name: &str, f: impl Fn(f32) -> f32) {
    // fixed-size arrays: no bounds checks, so the loop can auto-vectorize
    let mut input = [0f32; TP_ARR];
    for i in 0..TP_ARR {
        input[i] = 2.0 + (i as f32) * (2.0 / TP_ARR as f32);
    }
    let input = black_box(input);
    let mut out = [0f32; TP_ARR];
    let mut best = f64::INFINITY;
    for _ in 0..REPS {
        let start = Instant::now();
        for _ in 0..TP_PASSES {
            for (o, &x) in out.iter_mut().zip(input.iter()) {
                *o = f(x);
            }
            black_box(&mut out);
        }
        let total = (TP_ARR * TP_PASSES) as f64;
        let ns = start.elapsed().as_nanos() as f64 / total;
        best = best.min(ns);
    }
    println!("{:22} throughput {:6.3} ns/op (min of {REPS})", name, best);
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let filter = args.get(1).map(|s| s.as_str()).unwrap_or("");
    let run = |n: &str| filter.is_empty() || n.contains(filter);

    macro_rules! bench {
        ($name:expr, $f:expr) => {
            if run($name) {
                bench_latency($name, $f);
                bench_throughput($name, $f);
            }
        };
    }

    bench!("nop", |x: f32| x);
    bench!("cbrt", cbrt);
    bench!("std cbrt", |x: f32| x.cbrt());
    bench!("cbrt_accurate", cbrt_accurate);
    bench!("cbrt_throughput", cbrt_throughput);
    bench!("cbrt_fast", cbrt_fast);
    bench!("exp2", exp2);
    bench!("exp2_checked", exp2_checked);
    bench!("std exp2", |x: f32| x.exp2());
    bench!("log_2", log_2);
    bench!("std log2", |x: f32| x.log2());
    bench!("sin", sin);
    bench!("std sin", |x: f32| x.sin());
    bench!("cos", cos);
    bench!("std cos", |x: f32| x.cos());
}
