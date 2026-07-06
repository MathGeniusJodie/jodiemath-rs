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

/// Fuzz-only sweep for two-argument functions (no exhaustive mode: 2^64
/// pairs is infeasible). Otherwise mirrors `sweep` -- random bit patterns
/// for both args, split across all cores.
fn fuzz2(
    samples: u64,
    in_domain: impl Fn(f32, f32) -> bool + Sync,
    f: impl Fn(f32, f32) -> f32 + Sync,
    reference: impl Fn(f64, f64) -> f64 + Sync,
) -> Stats {
    let threads = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1) as u64;
    let chunk = samples.div_ceil(threads);
    std::thread::scope(|scope| {
        (0..threads)
            .map(|_| {
                let in_domain = &in_domain;
                let f = &f;
                let reference = &reference;
                scope.spawn(move || {
                    let mut s = Stats::zero();
                    let mut xs = [0f32; BATCH];
                    let mut ys = [0f32; BATCH];
                    let mut zs = [0f32; BATCH];
                    let mut i = 0u64;
                    while i < chunk {
                        let n = (chunk - i).min(BATCH as u64) as usize;
                        for j in 0..n {
                            xs[j] = f32::from_bits(rand::rng().random::<u32>());
                            ys[j] = f32::from_bits(rand::rng().random::<u32>());
                        }
                        for j in 0..n {
                            zs[j] = f(xs[j], ys[j]);
                        }
                        for j in 0..n {
                            if in_domain(xs[j], ys[j]) {
                                let r = reference(xs[j] as f64, ys[j] as f64) as f32;
                                let d = ulp_diff(zs[j], r);
                                s.sum += d;
                                if d > s.max {
                                    s.max = d;
                                    s.worst_x = xs[j];
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
        // unchecked sin's documented exact-integer range: |x| < 2^22 * pi
        // (the round-via-fma magic-constant trick's exact range). Filtering
        // both jodie's and std's inputs to the same set keeps the
        // comparison apples-to-apples, same as exp2/exp2_checked above.
        let sin_domain = |x: f32| x.abs() < (1u32 << 22) as f32 * std::f32::consts::PI;
        for (name, hi) in [
            ("sin |x|<=pi/4", 0.785398_f32),
            ("sin |x|<=10", 10.0),
            ("sin |x|<=1000", 1000.0),
            ("sin |x|<=1e6", 1e6),
        ] {
            let domain = move |x: f32| x.abs() <= hi;
            let s = measure!(domain, sin, |x: f64| x.sin());
            report(name, &s, t0);
        }
        let s = measure!(sin_domain, sin, |x: f64| x.sin());
        report("sin (in-domain)", &s, t0);
        let s = measure!(sin_domain, |x: f32| x.sin(), |x: f64| x.sin());
        report("std sin (in-domain)", &s, t0);
        for (name, hi) in [
            ("sin_checked |x|<=pi/4", 0.785398_f32),
            ("sin_checked |x|<=10", 10.0),
            ("sin_checked |x|<=1000", 1000.0),
            ("sin_checked |x|<=1e6", 1e6),
        ] {
            let domain = move |x: f32| x.abs() <= hi;
            let s = measure!(domain, sin_checked, |x: f64| x.sin());
            report(name, &s, t0);
        }
        // magnitude buckets past the 1-ulp-average guarantee: verifies the
        // reduction degrades gradually (not a cliff) well beyond 1e6, per
        // its doc comment.
        for (name, lo, hi) in [
            ("sin_checked [1e7,1e8)", 1e7, 1e8),
            ("sin_checked [1e9,1e10)", 1e9, 1e10),
            ("sin_checked [1e12,1e13)", 1e12, 1e13),
            ("sin_checked [1e15,1e16)", 1e15, 1e16),
        ] {
            let domain = move |x: f32| x.abs() >= lo && x.abs() < hi;
            let s = measure!(domain, sin_checked, |x: f64| x.sin());
            report(name, &s, t0);
        }
        let s = measure!(everywhere, sin_checked, |x: f64| x.sin());
        report("sin_checked (all f32)", &s, t0);
        let s = measure!(everywhere, |x: f32| x.sin(), |x: f64| x.sin());
        report("std sin (all f32)", &s, t0);
    }
    if run("cos") {
        let cos_domain = |x: f32| x.abs() < (1u32 << 22) as f32 * std::f32::consts::PI;
        for (name, hi) in [
            ("cos |x|<=pi/4", 0.785398_f32),
            ("cos |x|<=10", 10.0),
            ("cos |x|<=1000", 1000.0),
            ("cos |x|<=1e6", 1e6),
        ] {
            let domain = move |x: f32| x.abs() <= hi;
            let s = measure!(domain, cos, |x: f64| x.cos());
            report(name, &s, t0);
        }
        let s = measure!(cos_domain, cos, |x: f64| x.cos());
        report("cos (in-domain)", &s, t0);
        let s = measure!(cos_domain, |x: f32| x.cos(), |x: f64| x.cos());
        report("std cos (in-domain)", &s, t0);
        for (name, hi) in [
            ("cos_checked |x|<=pi/4", 0.785398_f32),
            ("cos_checked |x|<=10", 10.0),
            ("cos_checked |x|<=1000", 1000.0),
            ("cos_checked |x|<=1e6", 1e6),
        ] {
            let domain = move |x: f32| x.abs() <= hi;
            let s = measure!(domain, cos_checked, |x: f64| x.cos());
            report(name, &s, t0);
        }
        for (name, lo, hi) in [
            ("cos_checked [1e7,1e8)", 1e7, 1e8),
            ("cos_checked [1e9,1e10)", 1e9, 1e10),
            ("cos_checked [1e12,1e13)", 1e12, 1e13),
            ("cos_checked [1e15,1e16)", 1e15, 1e16),
        ] {
            let domain = move |x: f32| x.abs() >= lo && x.abs() < hi;
            let s = measure!(domain, cos_checked, |x: f64| x.cos());
            report(name, &s, t0);
        }
        let s = measure!(everywhere, cos_checked, |x: f64| x.cos());
        report("cos_checked (all f32)", &s, t0);
        let s = measure!(everywhere, |x: f32| x.cos(), |x: f64| x.cos());
        report("std cos (all f32)", &s, t0);
    }

    if run("ln") {
        let s = measure!(everywhere, ln, |x: f64| x.ln());
        report("ln", &s, t0);
        let s = measure!(everywhere, |x: f32| x.ln(), |x: f64| x.ln());
        report("std ln", &s, t0);
    }
    if run("log10") {
        let s = measure!(everywhere, log10, |x: f64| x.log10());
        report("log10", &s, t0);
        let s = measure!(everywhere, |x: f32| x.log10(), |x: f64| x.log10());
        report("std log10", &s, t0);
    }
    if run("log1p") {
        let s = measure!(everywhere, log1p, |x: f64| x.ln_1p());
        report("log1p", &s, t0);
        let s = measure!(everywhere, |x: f32| x.ln_1p(), |x: f64| x.ln_1p());
        report("std log1p", &s, t0);
    }
    if run("expm1") {
        // expm1's large-|x| branch (and exp itself) calls exp2, only
        // accurate while x*log2(e) stays inside exp2's unchecked domain
        // (see the "exp" block above for exp2 itself).
        let exp_domain = |x: f32| (-126.0..128.0).contains(&(x * std::f32::consts::LOG2_E));
        let s = measure!(exp_domain, exp, |x: f64| x.exp());
        report("exp", &s, t0);
        let s = measure!(exp_domain, |x: f32| x.exp(), |x: f64| x.exp());
        report("std exp", &s, t0);
        let s = measure!(exp_domain, expm1, |x: f64| x.exp_m1());
        report("expm1", &s, t0);
        let s = measure!(exp_domain, |x: f32| x.exp_m1(), |x: f64| x.exp_m1());
        report("std expm1", &s, t0);
    }
    if run("sinh") {
        // sinh/cosh use both exp(x) and exp(-x): restrict to where both
        // stay inside exp2's unchecked domain.
        let sinh_domain = |x: f32| {
            let e = x * std::f32::consts::LOG2_E;
            e > -126.0 && e < 126.0
        };
        let s = measure!(sinh_domain, sinh, |x: f64| x.sinh());
        report("sinh", &s, t0);
        let s = measure!(sinh_domain, |x: f32| x.sinh(), |x: f64| x.sinh());
        report("std sinh", &s, t0);
        let s = measure!(sinh_domain, cosh, |x: f64| x.cosh());
        report("cosh", &s, t0);
        let s = measure!(sinh_domain, |x: f32| x.cosh(), |x: f64| x.cosh());
        report("std cosh", &s, t0);
    }
    if run("tanh") {
        // tanh uses exp(2x): same reasoning as sinh/cosh, halved.
        let tanh_domain = |x: f32| {
            let e = 2.0 * x * std::f32::consts::LOG2_E;
            e > -126.0 && e < 126.0
        };
        let s = measure!(tanh_domain, tanh, |x: f64| x.tanh());
        report("tanh", &s, t0);
        let s = measure!(tanh_domain, |x: f32| x.tanh(), |x: f64| x.tanh());
        report("std tanh", &s, t0);
    }
    if run("asinh") {
        let s = measure!(everywhere, asinh, |x: f64| x.asinh());
        report("asinh", &s, t0);
        let s = measure!(everywhere, |x: f32| x.asinh(), |x: f64| x.asinh());
        report("std asinh", &s, t0);
    }
    if run("acosh") {
        let s = measure!(everywhere, acosh, |x: f64| x.acosh());
        report("acosh", &s, t0);
        let s = measure!(everywhere, |x: f32| x.acosh(), |x: f64| x.acosh());
        report("std acosh", &s, t0);
    }
    if run("atanh") {
        let s = measure!(everywhere, atanh, |x: f64| x.atanh());
        report("atanh", &s, t0);
        let s = measure!(everywhere, |x: f32| x.atanh(), |x: f64| x.atanh());
        report("std atanh", &s, t0);
    }
    if run("asin") {
        let s = measure!(everywhere, asin, |x: f64| x.asin());
        report("asin", &s, t0);
        let s = measure!(everywhere, |x: f32| x.asin(), |x: f64| x.asin());
        report("std asin", &s, t0);
    }
    if run("acos") {
        let s = measure!(everywhere, acos, |x: f64| x.acos());
        report("acos", &s, t0);
        let s = measure!(everywhere, |x: f32| x.acos(), |x: f64| x.acos());
        report("std acos", &s, t0);
    }
    if run("atan") {
        let s = measure!(everywhere, atan, |x: f64| x.atan());
        report("atan", &s, t0);
        let s = measure!(everywhere, |x: f32| x.atan(), |x: f64| x.atan());
        report("std atan", &s, t0);
    }
    if run("tan") {
        let tan_domain = |x: f32| x.abs() < (1u32 << 22) as f32 * std::f32::consts::PI;
        let s = measure!(tan_domain, tan, |x: f64| x.tan());
        report("tan (in-domain)", &s, t0);
        let s = measure!(tan_domain, |x: f32| x.tan(), |x: f64| x.tan());
        report("std tan (in-domain)", &s, t0);
    }
    if run("erf") {
        // erf's tail branch calls exp2(erf_poly(|x|)), which overflows
        // exp2's unchecked domain once |x| is beyond ~8.5 -- well past
        // where erf has already saturated to +-1 at f32 precision (see
        // erf's doc comment), so this bound only excludes the range where
        // the "true" answer is indistinguishable from the saturated value.
        let erf_domain = |x: f32| x.abs() < 6.0;
        let s = measure!(erf_domain, erf, libm::erf);
        report("erf", &s, t0);
        // erfc clamps to |x|<=10 before its own exp2 call, but that clamp
        // doesn't fully protect exp2's domain -- see erfc's doc comment for
        // the ~9.35 threshold this mirrors.
        let erfc_domain = |x: f32| x.abs() < 9.3;
        let s = measure!(erfc_domain, erfc, libm::erfc);
        report("erfc", &s, t0);
    }

    // two-argument functions: fuzz-only (exhaustive over 2^64 pairs isn't
    // feasible), smaller sample count since each trial needs two RNG draws.
    const TWOARG_SAMPLES: u64 = 10_000_000;
    if run("atan2") {
        let s = fuzz2(TWOARG_SAMPLES, |_, _| true, atan2, |y: f64, x: f64| y.atan2(x));
        report("atan2", &s, t0);
        let s = fuzz2(TWOARG_SAMPLES, |_, _| true, |y: f32, x: f32| y.atan2(x), |y: f64, x: f64| y.atan2(x));
        report("std atan2", &s, t0);
    }
    if run("hypot") {
        // naive x*x+y*y overflows f32 once |x| or |y| exceeds ~sqrt(f32::MAX)
        // (~1.8e19), and underflows (or flushes clean to 0, losing the
        // input's magnitude entirely) once |x|,|y| drop below ~sqrt of the
        // smallest denormal (~3.7e-23) -- hypot's doc comment calls both out
        // as the accepted tradeoff for avoiding std::hypot's anti-overflow
        // rescaling, so (like exp2/erf/erfc above) restrict to a range
        // clear of both, for a meaningful ulp number.
        let hypot_domain = |x: f32, y: f32| {
            let ok = |v: f32| v == 0.0 || (v.abs() > 1e-15 && v.abs() < 1e18);
            ok(x) && ok(y)
        };
        let s = fuzz2(TWOARG_SAMPLES, hypot_domain, hypot, |x: f64, y: f64| x.hypot(y));
        report("hypot", &s, t0);
        let s = fuzz2(TWOARG_SAMPLES, hypot_domain, |x: f32, y: f32| x.hypot(y), |x: f64, y: f64| x.hypot(y));
        report("std hypot", &s, t0);
    }
    if run("powf") {
        // restrict to jodie powf's actual domain: x > 0 (log_2's domain),
        // and the exponent log2(x)*y kept inside exp2's unchecked range.
        let pow_domain = |x: f32, y: f32| x > 0.0 && (-126.0..128.0).contains(&(x.log2() * y));
        let s = fuzz2(TWOARG_SAMPLES, pow_domain, powf, |x: f64, y: f64| x.powf(y));
        report("powf", &s, t0);
        let s = fuzz2(TWOARG_SAMPLES, pow_domain, |x: f32, y: f32| x.powf(y), |x: f64, y: f64| x.powf(y));
        report("std powf", &s, t0);
    }
    if run("remainder") {
        // x - round(x/y)*y loses precision to cancellation once |x/y| is
        // large: round(x/y)*y's absolute error scales with ulp(x), which
        // swamps the true remainder (at most |y|/2) once x/y is big enough
        // -- an inherited property of the naive formula (same in the C
        // original), not specific to this port. Bound |x/y| to stay in the
        // formula's reliable range.
        let remainder_domain = |x: f32, y: f32| y != 0.0 && (x / y).abs() < 1000.0;
        let s = fuzz2(TWOARG_SAMPLES, remainder_domain, remainder, libm::remainder);
        report("remainder", &s, t0);
    }

    println!("total: {:.2}s", t0.elapsed().as_secs_f64());
}
