// Accuracy harness: ULP sweep of each function against an f64-computed
// reference, rounded to f32.
//
// Two modes:
//   - quick (default): fuzz-style uniform random bit patterns, fast enough
//     to run on every iteration.
//   - thorough: exhaustively evaluates every one of the 2^32 f32 bit
//     patterns (both signs, every denormal, every NaN payload) -- run this
//     before trusting an accuracy number in the readme, not on every change.
// Both modes share one engine (`sweep`): only how a trial index maps to a
// bit pattern differs. Work is split across all cores with std::thread,
// and within a thread the jodie function is evaluated over a fixed-size
// array (same idiom as examples/quickbench.rs's throughput loop) so it
// auto-vectorizes instead of paying scalar call overhead per element --
// per the readme, these should run close to memcpy speed once vectorized.
// The reference call is f64 libm, unavoidably scalar, and dominates wall
// time regardless.
//
// Usage: cargo run --release --example accuracy [thorough|quick] [filter]
use jodiemath_rs::*;
use rand::RngExt;
use std::time::Instant;

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

struct Stats {
    sum: u64,
    max: u64,
    worst_x: f32,
    n: u64,
}

impl Stats {
    fn zero() -> Stats {
        Stats { sum: 0, max: 0, worst_x: 0.0, n: 0 }
    }
    fn combine(self, other: Stats) -> Stats {
        Stats {
            sum: self.sum + other.sum,
            max: self.max.max(other.max),
            worst_x: if self.max >= other.max { self.worst_x } else { other.worst_x },
            n: self.n + other.n,
        }
    }
}

fn report(name: &str, s: &Stats, since: Instant) {
    println!(
        "{:24} avg ulp {:>10.4}  max ulp {:>10}  worst x {:e} ({:>12} samples, {:>7.2}s elapsed)",
        name,
        s.sum as f64 / s.n as f64,
        s.max,
        s.worst_x,
        s.n,
        since.elapsed().as_secs_f64(),
    );
}

const BATCH: usize = 4096;

/// Runs `total` trials split across all cores. `bit_at(i)` maps a trial
/// index to the f32 bit pattern to test. `in_domain` filters which decoded
/// values count towards the stats; out-of-domain values are still fed
/// through `f` (see below) but not scored. NaN inputs are not filtered out
/// specially: ulp_diff treats any NaN-vs-NaN pair as a 0-ulp match, so a
/// full sweep also checks that every one of the ~2^25 NaN payloads still
/// propagates to *some* NaN instead of silently producing a finite value.
///
/// `f` is evaluated over a whole fixed-size array per batch, not one
/// element at a time, so the loop auto-vectorizes (bounds-check-free,
/// no branch inside) -- the same reason out-of-domain inputs aren't
/// filtered before this call: a per-element skip would reintroduce a
/// branch and block vectorization.
fn sweep(
    total: u64,
    bit_at: impl Fn(u64) -> u32 + Sync,
    in_domain: impl Fn(f32) -> bool + Sync,
    f: impl Fn(f32) -> f32 + Sync,
    reference: impl Fn(f64) -> f64 + Sync,
) -> Stats {
    let threads = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1) as u64;
    let chunk = total.div_ceil(threads);
    std::thread::scope(|scope| {
        (0..threads)
            .map(|t| {
                let bit_at = &bit_at;
                let in_domain = &in_domain;
                let f = &f;
                let reference = &reference;
                let start = t * chunk;
                let end = ((t + 1) * chunk).min(total);
                scope.spawn(move || {
                    let mut s = Stats::zero();
                    let mut xs = [0f32; BATCH];
                    let mut ys = [0f32; BATCH];
                    let mut i = start;
                    while i < end {
                        let n = (end - i).min(BATCH as u64) as usize;
                        for (j, x) in xs.iter_mut().enumerate().take(n) {
                            *x = f32::from_bits(bit_at(i + j as u64));
                        }
                        for (y, &x) in ys.iter_mut().zip(xs.iter()).take(n) {
                            *y = f(x);
                        }
                        for j in 0..n {
                            let x = xs[j];
                            if in_domain(x) {
                                let r = reference(x as f64) as f32;
                                let d = ulp_diff(ys[j], r);
                                s.sum += d;
                                if d > s.max {
                                    s.max = d;
                                    s.worst_x = x;
                                }
                                s.n += 1;
                            }
                        }
                        i += n as u64;
                    }
                    s
                })
            })
            .collect::<Vec<_>>()
            .into_iter()
            .map(|h| h.join().unwrap())
            .reduce(Stats::combine)
            .unwrap()
    })
}

fn exhaustive(
    in_domain: impl Fn(f32) -> bool + Sync,
    f: impl Fn(f32) -> f32 + Sync,
    reference: impl Fn(f64) -> f64 + Sync,
) -> Stats {
    sweep(1u64 << 32, |i| i as u32, in_domain, f, reference)
}

fn fuzz(
    samples: u64,
    in_domain: impl Fn(f32) -> bool + Sync,
    f: impl Fn(f32) -> f32 + Sync,
    reference: impl Fn(f64) -> f64 + Sync,
) -> Stats {
    sweep(samples, |_| rand::rng().random::<u32>(), in_domain, f, reference)
}

fn main() {
    if cfg!(debug_assertions) {
        eprintln!(
            "accuracy: refusing to run in a debug build -- thorough mode alone is ~4.3 billion \
             evaluations per function. Run with cargo run --release --example accuracy"
        );
        std::process::exit(1);
    }

    let args: Vec<String> = std::env::args().skip(1).collect();
    let thorough = args.iter().any(|a| a == "thorough" || a == "--thorough");
    let filter = args
        .iter()
        .find(|a| !matches!(a.as_str(), "thorough" | "--thorough" | "quick" | "--quick"))
        .cloned()
        .unwrap_or_default();
    let run = |n: &str| filter.is_empty() || n.contains(filter.as_str());

    // fuzz-mode sample count: chosen so the whole suite finishes in a few
    // seconds (reference calls dominate at ~20ns/sample single-threaded,
    // parallelized across all cores here)
    const QUICK_SAMPLES: u64 = 100_000_000;

    println!(
        "mode: {} ({} cores)",
        if thorough { "thorough (exhaustive, every f32 bit pattern)" } else { "quick (fuzz, random bit patterns)" },
        std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1)
    );
    let t0 = Instant::now();

    macro_rules! measure {
        ($domain:expr, $f:expr, $reference:expr) => {
            if thorough {
                exhaustive($domain, $f, $reference)
            } else {
                fuzz(QUICK_SAMPLES, $domain, $f, $reference)
            }
        };
    }

    let everywhere = |_: f32| true;

    if run("cbrt") {
        let s = measure!(everywhere, cbrt, |x: f64| x.cbrt());
        report("cbrt", &s, t0);
        let s = measure!(everywhere, cbrt_accurate, |x: f64| x.cbrt());
        report("cbrt_accurate", &s, t0);
        // bit-trick-only experiments: only ever designed/tuned for
        // positive normal x (their bit tricks assume a normal exponent
        // field), so including denormals just measures "how wrong is a
        // function fed inputs it was never meant for" (blows up to ~1e8
        // ulp there, swamping the number that's actually informative)
        let positive_normal = |x: f32| x >= f32::MIN_POSITIVE && x.is_finite();
        let s = measure!(positive_normal, cbrt_throughput, |x: f64| x.cbrt());
        report("cbrt_throughput (+)", &s, t0);
        let s = measure!(positive_normal, cbrt_fast, |x: f64| x.cbrt());
        report("cbrt_fast (+)", &s, t0);
        let s = measure!(everywhere, |x: f32| x.cbrt(), |x: f64| x.cbrt());
        report("std cbrt", &s, t0);
    }
    if run("log") {
        let s = measure!(everywhere, log_2, |x: f64| x.log2());
        report("log_2", &s, t0);
        let s = measure!(everywhere, |x: f32| x.log2(), |x: f64| x.log2());
        report("std log2", &s, t0);
    }
    if run("exp") {
        // unchecked exp2's documented domain: [-126, 128) (normal results
        // only). Filtering both jodie's and std's inputs to the same set
        // keeps the comparison apples-to-apples.
        let exp2_domain = |x: f32| (-126.0..128.0).contains(&x);
        let s = measure!(exp2_domain, exp2, |x: f64| x.exp2());
        report("exp2", &s, t0);
        let s = measure!(everywhere, exp2_checked, |x: f64| x.exp2());
        report("exp2_checked", &s, t0);
        let s = measure!(exp2_domain, |x: f32| x.exp2(), |x: f64| x.exp2());
        report("std exp2", &s, t0);
    }
    if run("sin") {
        for (name, hi) in [("sin |x|<=pi/4", 0.785398_f32), ("sin |x|<=10", 10.0), ("sin |x|<=1000", 1000.0)] {
            let domain = move |x: f32| x.abs() <= hi;
            let s = measure!(domain, sin, |x: f64| x.sin());
            report(name, &s, t0);
        }
        let s = measure!(everywhere, sin, |x: f64| x.sin());
        report("sin (all f32)", &s, t0);
        let s = measure!(everywhere, |x: f32| x.sin(), |x: f64| x.sin());
        report("std sin (all f32)", &s, t0);
    }
    if run("cos") {
        for (name, hi) in [("cos |x|<=pi/4", 0.785398_f32), ("cos |x|<=10", 10.0), ("cos |x|<=1000", 1000.0)] {
            let domain = move |x: f32| x.abs() <= hi;
            let s = measure!(domain, cos, |x: f64| x.cos());
            report(name, &s, t0);
        }
        let s = measure!(everywhere, cos, |x: f64| x.cos());
        report("cos (all f32)", &s, t0);
        let s = measure!(everywhere, |x: f32| x.cos(), |x: f64| x.cos());
        report("std cos (all f32)", &s, t0);
    }

    println!("total: {:.2}s", t0.elapsed().as_secs_f64());
}
