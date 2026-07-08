// Accuracy harness: ULP sweep of each function against an f64 reference,
// rounded to f32.
//
// Two modes:
//   - quick (default): fuzz-style uniform random bit patterns, fast enough
//     to run on every iteration.
//   - thorough: exhaustively evaluates every one of the 2^32 f32 bit
//     patterns (both signs, every denormal, every NaN payload) -- run this
//     before trusting an accuracy number in the readme, not on every change.
// Both modes share one engine (`sweep`): only how a trial index maps to a
// bit pattern differs. Work is split across half the machine's cores
// (`worker_threads`), and the process nices itself (`nice_self`) so the rest
// stays responsive for other work while a sweep runs -- both cost nothing
// when the machine is otherwise idle, since niceness only matters under
// contention. Within a thread the jodie function is evaluated over a
// fixed-size array (same idiom as examples/quickbench.rs's throughput loop)
// so it auto-vectorizes instead of paying scalar call overhead per element --
// per the readme, these should run close to memcpy speed once vectorized.
// The f64 reference is computed the same way now: batched through the
// `sleef` crate's SIMD functions instead of one scalar libm call per
// element, so it's vectorized too and no longer dominates wall time the way
// a scalar reference used to. Uses the cheapest ULP bucket sleef offers per
// function (u35 where available, else u10/u15) -- even 3.5 ULP of *f64*
// error is ~1e8x tighter than f32 ever needs, so it's plenty of headroom for
// a reference we're rounding down to f32 anyway, and the coarser buckets
// need less internal precision (no double-float correction), so they're
// faster. Out-of-domain inputs are still
// fed through both `f` and `reference` (see `sweep`'s doc comment below) but
// not scored -- a per-element skip would reintroduce a branch and block
// vectorization in either one.
//
// Requires nightly (sleef depends on the unstable `portable_simd` feature):
//   cargo +nightly run --release --example accuracy [thorough|quick] [filter]
#![feature(portable_simd)]
use jodiemath_rs::*;
use rand::RngExt;
// Precision buckets picked for speed, not accuracy: even u35 (3.5 ULP of
// *f64*) is ~1e8x tighter than f32 ever needs, so it's plenty for a ground
// truth we're going to round down to f32 anyway. Only functions with no u35
// variant fall back to u10 (or u15 for erfc); pow/remainder have no ULP
// bucket at all (exact-ish by construction).
use sleef::f64x::{
    acos_u35, acosh_u10, asin_u35, asinh_u10, atan2_u35, atan_u35, atanh_u10, cbrt_u35, cos_u35,
    cosh_u35, erf_u10, erfc_u15, exp10_u35, exp2_u35, exp_u10, expm1_u10, hypot_u35, log10_u10,
    log1p_u10, log2_u35, log_u35, pow_u10, remainder as remainder_ref, sin_u35, sinh_u35, tan_u35,
    tanh_u35,
};
use std::simd::num::SimdFloat;
use std::simd::{Select, Simd, StdFloat};
use std::time::Instant;

const LANES: usize = 8;
type F64xN = Simd<f64, LANES>;

/// Workaround for a sleef-rs 0.3.3 bug: sin_u35/cos_u35/tan_u35's (and the
/// u10 variants', same underlying code) large-argument path (`rempi`)
/// derives a lookup-table index from the input's exponent, and
/// NaN/+-Infinity's sentinel exponent overflows that table (panics) instead
/// of producing a result. sin/cos/tan of any non-finite input is NaN
/// regardless (matches std), so route those lanes around sleef entirely
/// instead of feeding it a value it can't handle.
fn trig_safe(v: F64xN, f: impl Fn(F64xN) -> F64xN) -> F64xN {
    let finite = v.is_finite();
    let safe = finite.select(v, F64xN::splat(0.0));
    finite.select(f(safe), F64xN::splat(f64::NAN))
}
fn sin_ref(v: F64xN) -> F64xN {
    trig_safe(v, sin_u35)
}
fn cos_ref(v: F64xN) -> F64xN {
    trig_safe(v, cos_u35)
}
fn tan_ref(v: F64xN) -> F64xN {
    trig_safe(v, tan_u35)
}
// sinpi/cospi's own point: q=round(x), r=x-q is *exact*, so pi*r is a
// small, precisely-representable angle -- computing pi*x directly (the
// naive reference) reintroduces the argument-reduction imprecision
// sinpi/cospi exist to avoid, and would give a *wrong* reference right
// at sin's own zeros (any tiny rounding error in a large pi*x is a huge
// *relative* error exactly where the true value is ~0). Mirror the same
// reduction here, just in f64, so the reference is trustworthy at any
// magnitude the exact-in-f32 reduction is (all of f32, in principle;
// this harness still only trusts it up to where f64 keeps q exact,
// i.e. far beyond f32's own 2^24 limit).
fn parity_f64(q: F64xN) -> F64xN {
    q - F64xN::splat(2.0) * (q * F64xN::splat(0.5)).floor()
}
fn sinpi_ref(v: F64xN) -> F64xN {
    trig_safe(v, |x: F64xN| {
        let q = x.round();
        let r = x - q;
        let s = sin_u35(r * F64xN::splat(std::f64::consts::PI));
        let sign = F64xN::splat(1.0) - F64xN::splat(2.0) * parity_f64(q);
        s * sign
    })
}
fn cospi_ref(v: F64xN) -> F64xN {
    trig_safe(v, |x: F64xN| {
        let q = (x - F64xN::splat(0.5)).round() + F64xN::splat(0.5);
        let r = x - q;
        let s = sin_u35(r * F64xN::splat(std::f64::consts::PI));
        let sign = F64xN::splat(2.0) * parity_f64(q - F64xN::splat(0.5)) - F64xN::splat(1.0);
        s * sign
    })
}
// sind/cosd's own point (see their doc comments): q=round(x/180),
// d=x-q*180 keeps the residual small and precise, so d*pi/180 is a
// small, accurate angle -- computing x*pi/180 directly (the naive
// reference) reintroduces the same imprecision-at-large-x problem
// sinpi_ref's own comment describes. Mirror the reduction here too.
fn sind_ref(v: F64xN) -> F64xN {
    trig_safe(v, |x: F64xN| {
        let q = (x / F64xN::splat(180.0)).round();
        let d = x - q * F64xN::splat(180.0);
        let s = sin_u35(d * F64xN::splat(std::f64::consts::PI / 180.0));
        let sign = F64xN::splat(1.0) - F64xN::splat(2.0) * parity_f64(q);
        s * sign
    })
}
fn cosd_ref(v: F64xN) -> F64xN {
    trig_safe(v, |x: F64xN| {
        let q = (x / F64xN::splat(180.0) - F64xN::splat(0.5)).round() + F64xN::splat(0.5);
        let d = x - q * F64xN::splat(180.0);
        let s = sin_u35(d * F64xN::splat(std::f64::consts::PI / 180.0));
        let sign = F64xN::splat(2.0) * parity_f64(q - F64xN::splat(0.5)) - F64xN::splat(1.0);
        s * sign
    })
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

/// Lowers this process's scheduling priority so a sweep doesn't compete with
/// foreground work -- niceing only matters under CPU contention, so it costs
/// nothing when the machine is otherwise idle. Threads spawned later inherit
/// the nice value set here, so this only needs to run once, before `sweep`
/// or `fuzz2` spawn their workers.
#[cfg(unix)]
fn nice_self() {
    // SAFETY: setpriority(PRIO_PROCESS, 0, _) only ever affects the calling
    // process's own niceness; failure just leaves the default priority.
    if unsafe { libc::setpriority(libc::PRIO_PROCESS, 0, 19) } != 0 {
        eprintln!("accuracy: couldn't lower process priority (continuing anyway)");
    }
}
#[cfg(not(unix))]
fn nice_self() {}

/// Half the machine's logical cores (rounded down, minimum 1) -- leaves the
/// other half free for whatever else is running, alongside `nice_self`, to
/// keep the machine responsive while a sweep runs.
fn worker_threads() -> u64 {
    let total = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(2);
    (total / 2).max(1) as u64
}

const BATCH: usize = 4096;

/// Runs `total` trials split across half the machine's cores. `bit_at(i)`
/// maps a trial index to the f32 bit pattern to test. `in_domain` filters
/// which decoded values count towards the stats; out-of-domain values are
/// still fed through `f` and `reference` (see below) but not scored. NaN
/// inputs are not filtered out specially: ulp_diff treats any NaN-vs-NaN
/// pair as a 0-ulp match, so a full sweep also checks that every one of the
/// ~2^25 NaN payloads still propagates to *some* NaN instead of silently
/// producing a finite value.
///
/// `f` is evaluated over a whole fixed-size array per batch, not one
/// element at a time, so the loop auto-vectorizes (bounds-check-free,
/// no branch inside) -- the same reason out-of-domain inputs aren't
/// filtered before this call: a per-element skip would reintroduce a
/// branch and block vectorization. `reference` is likewise evaluated over
/// `LANES` f64s at a time via the `sleef` crate, instead of one scalar libm
/// call per element.
fn sweep(
    total: u64,
    bit_at: impl Fn(u64) -> u32 + Sync,
    in_domain: impl Fn(f32) -> bool + Sync,
    f: impl Fn(f32) -> f32 + Sync,
    reference: impl Fn(F64xN) -> F64xN + Sync,
) -> Stats {
    let threads = worker_threads();
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
                        let mut k = 0;
                        while k < n {
                            let m = (n - k).min(LANES);
                            let mut xd = [0f64; LANES];
                            for idx in 0..m {
                                xd[idx] = xs[k + idx] as f64;
                            }
                            let r = reference(F64xN::from_array(xd)).to_array();
                            for idx in 0..m {
                                let x = xs[k + idx];
                                if in_domain(x) {
                                    let rf = r[idx] as f32;
                                    let d = ulp_diff(ys[k + idx], rf);
                                    s.sum += d;
                                    if d > s.max {
                                        s.max = d;
                                        s.worst_x = x;
                                    }
                                    s.n += 1;
                                }
                            }
                            k += m;
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
/// for both args, split across half the machine's cores, reference batched
/// `LANES` at a time via `sleef`.
fn fuzz2(
    samples: u64,
    in_domain: impl Fn(f32, f32) -> bool + Sync,
    f: impl Fn(f32, f32) -> f32 + Sync,
    reference: impl Fn(F64xN, F64xN) -> F64xN + Sync,
) -> Stats {
    let threads = worker_threads();
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
                        let mut k = 0;
                        while k < n {
                            let m = (n - k).min(LANES);
                            let mut xd = [0f64; LANES];
                            let mut yd = [0f64; LANES];
                            for idx in 0..m {
                                xd[idx] = xs[k + idx] as f64;
                                yd[idx] = ys[k + idx] as f64;
                            }
                            let r =
                                reference(F64xN::from_array(xd), F64xN::from_array(yd)).to_array();
                            for idx in 0..m {
                                if in_domain(xs[k + idx], ys[k + idx]) {
                                    let rf = r[idx] as f32;
                                    let d = ulp_diff(zs[k + idx], rf);
                                    s.sum += d;
                                    if d > s.max {
                                        s.max = d;
                                        s.worst_x = xs[k + idx];
                                    }
                                    s.n += 1;
                                }
                            }
                            k += m;
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
    reference: impl Fn(F64xN) -> F64xN + Sync,
) -> Stats {
    sweep(1u64 << 32, |i| i as u32, in_domain, f, reference)
}

fn fuzz(
    samples: u64,
    in_domain: impl Fn(f32) -> bool + Sync,
    f: impl Fn(f32) -> f32 + Sync,
    reference: impl Fn(F64xN) -> F64xN + Sync,
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

    nice_self();

    let args: Vec<String> = std::env::args().skip(1).collect();
    let thorough = args.iter().any(|a| a == "thorough" || a == "--thorough");
    let filter = args
        .iter()
        .find(|a| !matches!(a.as_str(), "thorough" | "--thorough" | "quick" | "--quick"))
        .cloned()
        .unwrap_or_default();
    let run = |n: &str| filter.is_empty() || n.contains(filter.as_str());

    // fuzz-mode sample count: chosen so the whole suite finishes in a few
    // seconds (reference is vectorized via sleef now, so it no longer
    // dominates; still split across half the machine's cores, niced, so
    // iterating doesn't compete with foreground work)
    const QUICK_SAMPLES: u64 = 100_000_000;

    let total_cores = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
    println!(
        "mode: {} ({} of {} cores, niced)",
        if thorough { "thorough (exhaustive, every f32 bit pattern)" } else { "quick (fuzz, random bit patterns)" },
        worker_threads(),
        total_cores,
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
        let s = measure!(everywhere, cbrt, cbrt_u35);
        report("cbrt", &s, t0);
        let s = measure!(everywhere, cbrt_accurate, cbrt_u35);
        report("cbrt_accurate", &s, t0);
        // bit-trick-only experiments: only ever designed/tuned for
        // positive normal x (their bit tricks assume a normal exponent
        // field), so including denormals just measures "how wrong is a
        // function fed inputs it was never meant for" (blows up to ~1e8
        // ulp there, swamping the number that's actually informative)
        let positive_normal = |x: f32| x >= f32::MIN_POSITIVE && x.is_finite();
        let s = measure!(positive_normal, cbrt_throughput, cbrt_u35);
        report("cbrt_throughput (+)", &s, t0);
        let s = measure!(positive_normal, cbrt_fast, cbrt_u35);
        report("cbrt_fast (+)", &s, t0);
        let s = measure!(everywhere, |x: f32| x.cbrt(), cbrt_u35);
        report("std cbrt", &s, t0);
    }
    if run("log") {
        let s = measure!(everywhere, log_2, log2_u35);
        report("log_2", &s, t0);
        // positive_normal: log_2_unchecked's documented contract (see its
        // doc comment) -- undefined for zero/negative/denormal/inf/nan.
        let positive_normal = |x: f32| x >= f32::MIN_POSITIVE && x.is_finite();
        let s = measure!(positive_normal, log_2_unchecked, log2_u35);
        report("log_2_unchecked (+)", &s, t0);
        let s = measure!(everywhere, |x: f32| x.log2(), log2_u35);
        report("std log2", &s, t0);
    }
    if run("exp") {
        // unchecked exp2's documented domain: [-126, 128) (normal results
        // only). Filtering both jodie's and std's inputs to the same set
        // keeps the comparison apples-to-apples.
        let exp2_domain = |x: f32| (-126.0..128.0).contains(&x);
        let s = measure!(exp2_domain, exp2, exp2_u35);
        report("exp2", &s, t0);
        let s = measure!(everywhere, exp2_checked, exp2_u35);
        report("exp2_checked", &s, t0);
        let s = measure!(exp2_domain, |x: f32| x.exp2(), exp2_u35);
        report("std exp2", &s, t0);
        // exp10/exp10_checked's own domains, mirroring exp2/exp2_checked's
        // unchecked-vs-checked split: x*log2(10) must stay in exp2's own
        // [-126,128) (unchecked) or exp2_checked's wider [-151,128).
        let exp10_domain = |x: f32| (-126.0..128.0).contains(&(x * std::f32::consts::LOG2_10));
        let exp10_checked_domain = |x: f32| (-151.0..128.0).contains(&(x * std::f32::consts::LOG2_10));
        let s = measure!(exp10_domain, exp10, exp10_u35);
        report("exp10", &s, t0);
        let s = measure!(exp10_checked_domain, exp10_checked, exp10_u35);
        report("exp10_checked", &s, t0);
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
            let s = measure!(domain, sin, sin_ref);
            report(name, &s, t0);
        }
        let s = measure!(sin_domain, sin, sin_ref);
        report("sin (in-domain)", &s, t0);
        let s = measure!(sin_domain, |x: f32| x.sin(), sin_ref);
        report("std sin (in-domain)", &s, t0);
        for (name, hi) in [
            ("sin_checked |x|<=pi/4", 0.785398_f32),
            ("sin_checked |x|<=10", 10.0),
            ("sin_checked |x|<=1000", 1000.0),
            ("sin_checked |x|<=1e6", 1e6),
        ] {
            let domain = move |x: f32| x.abs() <= hi;
            let s = measure!(domain, sin_checked, sin_ref);
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
            let s = measure!(domain, sin_checked, sin_ref);
            report(name, &s, t0);
        }
        let s = measure!(everywhere, sin_checked, sin_ref);
        report("sin_checked (all f32)", &s, t0);
        let s = measure!(everywhere, |x: f32| x.sin(), sin_ref);
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
            let s = measure!(domain, cos, cos_ref);
            report(name, &s, t0);
        }
        let s = measure!(cos_domain, cos, cos_ref);
        report("cos (in-domain)", &s, t0);
        let s = measure!(cos_domain, |x: f32| x.cos(), cos_ref);
        report("std cos (in-domain)", &s, t0);
        for (name, hi) in [
            ("cos_checked |x|<=pi/4", 0.785398_f32),
            ("cos_checked |x|<=10", 10.0),
            ("cos_checked |x|<=1000", 1000.0),
            ("cos_checked |x|<=1e6", 1e6),
        ] {
            let domain = move |x: f32| x.abs() <= hi;
            let s = measure!(domain, cos_checked, cos_ref);
            report(name, &s, t0);
        }
        for (name, lo, hi) in [
            ("cos_checked [1e7,1e8)", 1e7, 1e8),
            ("cos_checked [1e9,1e10)", 1e9, 1e10),
            ("cos_checked [1e12,1e13)", 1e12, 1e13),
            ("cos_checked [1e15,1e16)", 1e15, 1e16),
        ] {
            let domain = move |x: f32| x.abs() >= lo && x.abs() < hi;
            let s = measure!(domain, cos_checked, cos_ref);
            report(name, &s, t0);
        }
        let s = measure!(everywhere, cos_checked, cos_ref);
        report("cos_checked (all f32)", &s, t0);
        let s = measure!(everywhere, |x: f32| x.cos(), cos_ref);
        report("std cos (all f32)", &s, t0);
    }
    if run("sinpi") {
        // sinpi/cospi's own reduction (q=round(x), r=x-q, both exact in
        // f32) is valid over the *entire* f32 range -- but this test's
        // reference (x*pi computed directly in f64, then a standard
        // sin/cos) reintroduces exactly the argument-reduction precision
        // problem sinpi/cospi exist to avoid, once x is large enough that
        // even f64 can't represent x*pi to within a fraction of a half-
        // turn. Restricting to |x|<1e6 keeps the reference itself
        // trustworthy (f64 has 29 more mantissa bits than f32, plenty of
        // headroom at this scale); edgecheck.rs separately verifies the
        // large-x "exactly 0, never inf/nan" tail behavior without
        // relying on a precise reference there.
        let sinpi_domain = |x: f32| x.abs() < 1e6;
        let s = measure!(sinpi_domain, sinpi, sinpi_ref);
        report("sinpi (|x|<1e6)", &s, t0);
        let s = measure!(sinpi_domain, cospi, cospi_ref);
        report("cospi (|x|<1e6)", &s, t0);
    }
    if run("sind") {
        // sind/cosd's own exact-reduction range is ~4.7e7 (see their doc
        // comments, limited by 180.0's trailing zero bits, unlike
        // sinpi/cospi's full-f32-range exactness) -- |x|<1e6 is
        // comfortably inside that and keeps this test's own f64-based
        // reference trustworthy too.
        let sind_domain = |x: f32| x.abs() < 1e6;
        let s = measure!(sind_domain, sind, sind_ref);
        report("sind (|x|<1e6)", &s, t0);
        let s = measure!(sind_domain, cosd, cosd_ref);
        report("cosd (|x|<1e6)", &s, t0);
    }

    if run("ln") {
        let s = measure!(everywhere, ln, log_u35);
        report("ln", &s, t0);
        let positive_normal = |x: f32| x >= f32::MIN_POSITIVE && x.is_finite();
        let s = measure!(positive_normal, ln_unchecked, log_u35);
        report("ln_unchecked (+)", &s, t0);
        let s = measure!(everywhere, |x: f32| x.ln(), log_u35);
        report("std ln", &s, t0);
    }
    if run("log10") {
        let s = measure!(everywhere, log10, log10_u10);
        report("log10", &s, t0);
        let positive_normal = |x: f32| x >= f32::MIN_POSITIVE && x.is_finite();
        let s = measure!(positive_normal, log10_unchecked, log10_u10);
        report("log10_unchecked (+)", &s, t0);
        let s = measure!(everywhere, |x: f32| x.log10(), log10_u10);
        report("std log10", &s, t0);
    }
    if run("log1p") {
        let s = measure!(everywhere, log1p, log1p_u10);
        report("log1p", &s, t0);
        let s = measure!(everywhere, |x: f32| x.ln_1p(), log1p_u10);
        report("std log1p", &s, t0);
    }
    if run("expm1") {
        // expm1's large-|x| branch (and exp itself) calls exp2, only
        // accurate while x*log2(e) stays inside exp2's unchecked domain
        // (see the "exp" block above for exp2 itself).
        let exp_domain = |x: f32| (-126.0..128.0).contains(&(x * std::f32::consts::LOG2_E));
        let s = measure!(exp_domain, exp, exp_u10);
        report("exp", &s, t0);
        let s = measure!(exp_domain, |x: f32| x.exp(), exp_u10);
        report("std exp", &s, t0);
        let s = measure!(exp_domain, expm1, expm1_u10);
        report("expm1", &s, t0);
        let s = measure!(exp_domain, |x: f32| x.exp_m1(), expm1_u10);
        report("std expm1", &s, t0);
    }
    if run("sinh") {
        // sinh/cosh use both exp(x) and exp(-x): restrict to where both
        // stay inside exp2's unchecked domain.
        let sinh_domain = |x: f32| {
            let e = x * std::f32::consts::LOG2_E;
            e > -126.0 && e < 126.0
        };
        let s = measure!(sinh_domain, sinh, sinh_u35);
        report("sinh", &s, t0);
        let s = measure!(sinh_domain, |x: f32| x.sinh(), sinh_u35);
        report("std sinh", &s, t0);
        let s = measure!(sinh_domain, cosh, cosh_u35);
        report("cosh", &s, t0);
        let s = measure!(sinh_domain, |x: f32| x.cosh(), cosh_u35);
        report("std cosh", &s, t0);
        let s = measure!(sinh_domain, sinh_throughput, sinh_u35);
        report("sinh_throughput", &s, t0);
        let s = measure!(sinh_domain, cosh_throughput, cosh_u35);
        report("cosh_throughput", &s, t0);
    }
    if run("tanh") {
        // tanh uses exp(2x): same reasoning as sinh/cosh, halved.
        let tanh_domain = |x: f32| {
            let e = 2.0 * x * std::f32::consts::LOG2_E;
            e > -126.0 && e < 126.0
        };
        let s = measure!(tanh_domain, tanh, tanh_u35);
        report("tanh", &s, t0);
        let s = measure!(tanh_domain, |x: f32| x.tanh(), tanh_u35);
        report("std tanh", &s, t0);
    }
    if run("sigmoid") {
        // Domain matches sigmoid's own exp(-x) clamp (see its doc
        // comment): accurate while -x stays in exp's own good range.
        let sigmoid_domain = |x: f32| {
            let e = x * std::f32::consts::LOG2_E;
            e > -126.0 && e < 126.0
        };
        // Direct 1/(1+exp(-x)) in f64, *not* the 0.5+0.5*tanh(x/2)
        // identity -- that identity has exactly the cancellation bug
        // sigmoid's own doc comment describes, just pushed out to a
        // larger |x| in f64 (tanh saturates to exactly -1.0 once
        // |x/2| exceeds ~18.7 in f64, well inside this domain), so it
        // would silently give a *wrong* reference for part of the swept
        // range instead of a merely imprecise one.
        let sigmoid_ref = |v: F64xN| {
            F64xN::splat(1.0) / (F64xN::splat(1.0) + exp_u10(-v))
        };
        let s = measure!(sigmoid_domain, sigmoid, sigmoid_ref);
        report("sigmoid", &s, t0);
    }
    if run("asinh") {
        let s = measure!(everywhere, asinh, asinh_u10);
        report("asinh", &s, t0);
        let s = measure!(everywhere, |x: f32| x.asinh(), asinh_u10);
        report("std asinh", &s, t0);
    }
    if run("acosh") {
        let s = measure!(everywhere, acosh, acosh_u10);
        report("acosh", &s, t0);
        let s = measure!(everywhere, |x: f32| x.acosh(), acosh_u10);
        report("std acosh", &s, t0);
    }
    if run("atanh") {
        let s = measure!(everywhere, atanh, atanh_u10);
        report("atanh", &s, t0);
        let s = measure!(everywhere, |x: f32| x.atanh(), atanh_u10);
        report("std atanh", &s, t0);
    }
    if run("asin") {
        let s = measure!(everywhere, asin, asin_u35);
        report("asin", &s, t0);
        let s = measure!(everywhere, |x: f32| x.asin(), asin_u35);
        report("std asin", &s, t0);
    }
    if run("acos") {
        let s = measure!(everywhere, acos, acos_u35);
        report("acos", &s, t0);
        let s = measure!(everywhere, |x: f32| x.acos(), acos_u35);
        report("std acos", &s, t0);
    }
    if run("atan") {
        let s = measure!(everywhere, atan, atan_u35);
        report("atan", &s, t0);
        let s = measure!(everywhere, |x: f32| x.atan(), atan_u35);
        report("std atan", &s, t0);
    }
    if run("tan") {
        let tan_domain = |x: f32| x.abs() < (1u32 << 22) as f32 * std::f32::consts::PI;
        let s = measure!(tan_domain, tan, tan_ref);
        report("tan (in-domain)", &s, t0);
        let s = measure!(tan_domain, |x: f32| x.tan(), tan_ref);
        report("std tan (in-domain)", &s, t0);
    }
    if run("erf") {
        // erf_poly used to be evaluated unbounded on |x|, which was wrong
        // (not just imprecise) well before this bound -- see erf's doc
        // comment. Fixed now, so this measures the whole domain.
        let s = measure!(everywhere, erf, erf_u10);
        report("erf", &s, t0);
        // erfc clamps to |x|<=10 before its own exp2 call; that clamp used
        // to not fully protect exp2's unchecked domain (see erfc's doc
        // comment), fixed by routing through exp2_checked instead -- the
        // domain restriction here now only matches erfc's own clamp
        // (previously it stopped short at 9.3 specifically to dodge the
        // since-fixed bug).
        let erfc_domain = |x: f32| x.abs() <= 10.0;
        let s = measure!(erfc_domain, erfc, erfc_u15);
        report("erfc", &s, t0);
    }

    // two-argument functions: fuzz-only (exhaustive over 2^64 pairs isn't
    // feasible), smaller sample count since each trial needs two RNG draws.
    const TWOARG_SAMPLES: u64 = 10_000_000;
    if run("atan2") {
        let s = fuzz2(TWOARG_SAMPLES, |_, _| true, atan2, atan2_u35);
        report("atan2", &s, t0);
        let s = fuzz2(TWOARG_SAMPLES, |_, _| true, |y: f32, x: f32| y.atan2(x), atan2_u35);
        report("std atan2", &s, t0);
    }
    if run("rsqrt") {
        // x > 0.0 only (0/negative/nan/inf are all correct "for free" via
        // plain IEEE754 semantics, see rsqrt's own doc comment -- not a
        // fuzz-density target). f64's own sqrt (a real op, not a sleef
        // reference) is precise enough ground truth for verifying f32-level
        // rsqrt accuracy.
        let rsqrt_domain = |x: f32| x > 0.0 && x.is_finite();
        let s = measure!(rsqrt_domain, rsqrt, |v: F64xN| F64xN::splat(1.0) / v.sqrt());
        report("rsqrt", &s, t0);
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
        let s = fuzz2(TWOARG_SAMPLES, hypot_domain, hypot, hypot_u35);
        report("hypot", &s, t0);
        let s = fuzz2(TWOARG_SAMPLES, hypot_domain, |x: f32, y: f32| x.hypot(y), hypot_u35);
        report("std hypot", &s, t0);
    }
    if run("pown") {
        // pown(x, n) takes an i32 exponent, not the f32/f64 pair shape
        // the rest of this harness is built around (SIMD reference via
        // sleef) -- it's also a plain repeated-multiply algorithm with
        // no fitted approximation to verify at SIMD width, so a simple
        // scalar random sweep against f64::powi is sufficient ground
        // truth here, no vectorized reference needed.
        //
        // Two ranges, since exponentiation-by-squaring's accumulated
        // rounding grows with |n| (each squaring/multiply step is its
        // own correctly-rounded op, but there are ~log2(|n|) of them) --
        // small n (the overwhelmingly common real use: squares, cubes,
        // reciprocals) is near-perfect; large n degrades gracefully,
        // same characteristic any repeated-squaring algorithm has, not a
        // bug to chase further for a utility function like this one.
        let pown_sweep = |lo: i32, hi: i32, label: &str| {
            let mut sum = 0u64;
            let mut max = 0u64;
            let mut worst = (0.0f32, 0i32);
            let n_samples = 20_000_000u64;
            for _ in 0..n_samples {
                let x = f32::from_bits(rand::rng().random::<u32>());
                let n: i32 = rand::rng().random_range(lo..=hi);
                if !x.is_finite() || x == 0.0 {
                    continue;
                }
                let got = pown(x, n);
                let want = (x as f64).powi(n) as f32;
                let d = ulp_diff(got, want);
                sum += d;
                if d > max {
                    max = d;
                    worst = (x, n);
                }
            }
            println!(
                "{:24} avg ulp {:>10.4}  max ulp {:>10}  worst x={:e},n={} ({:>12} samples, {:>7.2}s elapsed)",
                label,
                sum as f64 / n_samples as f64,
                max,
                worst.0,
                worst.1,
                n_samples,
                t0.elapsed().as_secs_f64(),
            );
        };
        pown_sweep(-8, 8, "pown (|n|<=8)");
        pown_sweep(-64, 64, "pown (|n|<=64)");
    }
    if run("powf") {
        // x != 0 (x == 0 is its own exact case, not a fuzz-density target)
        // and the exponent log2(|x|)*y kept inside exp2's unchecked range.
        // Negative x is now in-domain too (see powf's own doc comment: a
        // real result exists whenever y is an integer) -- this used to be
        // filtered out entirely (`x > 0.0`), which dodged the bug instead
        // of exercising it, the same pattern erf/erfc's filters had.
        let pow_domain =
            |x: f32, y: f32| x != 0.0 && (-126.0..128.0).contains(&(x.abs().log2() * y));
        let s = fuzz2(TWOARG_SAMPLES, pow_domain, powf, pow_u10);
        report("powf", &s, t0);
        let s = fuzz2(TWOARG_SAMPLES, pow_domain, |x: f32, y: f32| x.powf(y), pow_u10);
        report("std powf", &s, t0);
    }
    if run("powf_checked") {
        // Same domain as powf's own sweep -- see powf_checked's doc
        // comment for why this variant is substantially more accurate
        // for large |y| (avg ulp -87% in a controlled comparison).
        let pow_domain =
            |x: f32, y: f32| x != 0.0 && (-126.0..128.0).contains(&(x.abs().log2() * y));
        let s = fuzz2(TWOARG_SAMPLES, pow_domain, powf_checked, pow_u10);
        report("powf_checked", &s, t0);
    }
    // Both remainder variants use *ties-away-from-zero* rounding for q (see
    // remainder's own doc comment) but sleef's `remainder_ref` implements
    // true IEEE754 remainder, which is ties-to-*even* -- a different, also
    // "correct," convention that disagrees with ties-away by a full `y`
    // whenever x/y lands acceptably close to an exact half-integer tie.
    // That's not a bug in either implementation, just two valid conventions
    // disagreeing at their boundary, but it shows up in this fuzz sweep as
    // an occasional spurious billions-of-ulp reading unrelated to either
    // function's real accuracy (confirmed by hand: sleef's remainder(2.5,
    // 1.0) = 0.5, ties-to-even, vs jodiemath's ties-away q=3 giving -0.5).
    // Excluded here so the sweep measures real accuracy, not convention
    // disagreement.
    let near_tie = |x: f32, y: f32| ((x / y).abs().fract() - 0.5).abs() < 1e-4;
    if run("remainder") {
        // x - round(x/y)*y loses precision to cancellation once |x/y| is
        // large: round(x/y)*y's absolute error scales with ulp(x), which
        // swamps the true remainder (at most |y|/2) once x/y is big enough
        // -- an inherited property of the naive formula (same in the C
        // original), not specific to this port. Bound |x/y| to stay in the
        // formula's reliable range. Also: even well inside that range, a
        // low-probability but real bug exists (see remainder's own doc
        // comment) where x/y's f32 division rounding crosses a tie boundary
        // that a more precise division wouldn't have -- remainder_checked
        // below fixes it; excluded from `remainder`'s own sweep via
        // near_tie since it's a documented, known limitation, not something
        // this sweep is meant to catch.
        let remainder_domain =
            |x: f32, y: f32| y != 0.0 && (x / y).abs() < 1000.0 && !near_tie(x, y);
        let s = fuzz2(TWOARG_SAMPLES, remainder_domain, remainder, remainder_ref);
        report("remainder", &s, t0);
    }
    if run("remainder_checked") {
        // remainder_checked() self-corrects q by one when x/y's own
        // division rounding pushed it to the wrong integer, which holds up
        // cleanly (0 max ulp, fuzz-tested) all the way up to where q itself
        // stops being an exactly-representable f32 integer (2^24) -- past
        // that point q's own rounding is the limit, not the division, and
        // no amount of one-integer nudging can fix it (confirmed:
        // bounding at 2e7 instead of 1e7 immediately produces billions of
        // ulp of error). Bound |x/y| comfortably under that 2^24 cliff.
        let remainder_domain =
            |x: f32, y: f32| y != 0.0 && (x / y).abs() < 10000000.0 && !near_tie(x, y);
        let s = fuzz2(TWOARG_SAMPLES, remainder_domain, remainder_checked, remainder_ref);
        report("remainder_checked", &s, t0);
    }

    println!("total: {:.2}s", t0.elapsed().as_secs_f64());
}
