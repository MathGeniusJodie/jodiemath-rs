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
use std::simd::cmp::SimdPartialEq;
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
        // Never forms `k+0.5` as a single f64 (loses precision once `x`
        // exceeds f64's own ~2^53 exact-integer range, the same bug
        // class the real cospi's own fix avoids for f32's ~2^24 -- see
        // its doc comment) -- k, r=(x-k)-0.5, and parity(k) stay
        // separate the whole way through instead.
        let k = (x - F64xN::splat(0.5)).round();
        let r = (x - k) - F64xN::splat(0.5);
        let s = sin_u35(r * F64xN::splat(std::f64::consts::PI));
        let sign = F64xN::splat(2.0) * parity_f64(k) - F64xN::splat(1.0);
        s * sign
    })
}
fn sinc_ref(v: F64xN) -> F64xN {
    trig_safe(v, |x: F64xN| {
        let is_zero = x.simd_eq(F64xN::splat(0.0));
        let safe_x = is_zero.select(F64xN::splat(1.0), x);
        let normal = sinpi_ref(x) / (F64xN::splat(std::f64::consts::PI) * safe_x);
        is_zero.select(F64xN::splat(1.0), normal)
    })
}
// tanpi/tand's own references: same ratio construction as the real
// functions (see their doc comments) -- reuses sinpi_ref/cospi_ref
// (resp. sind_ref/cosd_ref below) directly rather than a naive
// tan(pi*x)/tan(x*pi/180), for the same "don't reintroduce the
// large-x reduction imprecision the real function exists to avoid"
// reasoning as sinc_ref above. sinpi_ref/cospi_ref already each handle
// non-finite input via their own trig_safe wrapper, so no extra
// wrapping needed here.
fn tanpi_ref(v: F64xN) -> F64xN {
    sinpi_ref(v) / cospi_ref(v)
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
fn tand_ref(v: F64xN) -> F64xN {
    sind_ref(v) / cosd_ref(v)
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

// 4-arg analog of fuzz2, for diff_of_products/cross2 (backlog idea #135).
// Unlike fuzz2's transcendental references (sleef, batched via SIMD purely
// to amortize their own call cost), a*b-c*d's f64 reference is a few plain
// arithmetic ops -- no SIMD batching needed, plain scalar per-sample is
// already fast enough. No exhaustive/"thorough" mode either: a real 4-arg
// sweep is 2^128 combinations, not a coherent concept the way 1-arg
// exhaustive (2^32) is.
fn fuzz4(
    samples: u64,
    f: impl Fn(f32, f32, f32, f32) -> f32 + Sync,
    reference: impl Fn(f64, f64, f64, f64) -> f64 + Sync,
) -> Stats {
    let threads = worker_threads();
    let chunk = samples.div_ceil(threads);
    std::thread::scope(|scope| {
        (0..threads)
            .map(|_| {
                let f = &f;
                let reference = &reference;
                scope.spawn(move || {
                    let mut s = Stats::zero();
                    let mut i = 0u64;
                    while i < chunk {
                        i += 1;
                        let a = f32::from_bits(rand::rng().random::<u32>());
                        let b = f32::from_bits(rand::rng().random::<u32>());
                        let c = f32::from_bits(rand::rng().random::<u32>());
                        let d = f32::from_bits(rand::rng().random::<u32>());
                        if !(a.is_finite() && b.is_finite() && c.is_finite() && d.is_finite()) {
                            continue;
                        }
                        let refv32 = reference(a as f64, b as f64, c as f64, d as f64) as f32;
                        if !refv32.is_finite() {
                            continue;
                        }
                        let got = f(a, b, c, d);
                        let dd = ulp_diff(got, refv32);
                        s.sum += dd;
                        if dd > s.max {
                            s.max = dd;
                            s.worst_x = a;
                        }
                        s.n += 1;
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
        // cbrt_unchecked's own contract: x normal/finite, either sign
        // (unlike log_2_unchecked's positive-only domain -- cbrt_normal
        // already reapplies x's own sign bit internally).
        let normal_finite = |x: f32| x.abs() >= f32::MIN_POSITIVE && x.is_finite();
        let s = measure!(normal_finite, cbrt_unchecked, cbrt_u35);
        report("cbrt_unchecked (+)", &s, t0);
        // cbrt_accurate_unchecked's own contract: x already inside
        // cbrt_accurate's own safe rescale range (matches its `!small &&
        // !big` domain exactly, see cbrt_accurate's own doc comment).
        let accurate_safe_range = |x: f32| {
            let ax = x.to_bits() & 0x7fff_ffff;
            ax >= 0x2380_0000 && ax < 0x7f00_0000
        };
        let s = measure!(accurate_safe_range, cbrt_accurate_unchecked, cbrt_u35);
        report("cbrt_accurate_unchecked (+)", &s, t0);
    }
    if run("rcbrt") {
        let rcbrt_ref = |v: F64xN| F64xN::splat(1.0) / cbrt_u35(v);
        let s = measure!(everywhere, rcbrt, rcbrt_ref);
        report("rcbrt", &s, t0);
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
            ("sin |x|<=pi/4", std::f32::consts::FRAC_PI_4),
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
            ("sin_checked |x|<=pi/4", std::f32::consts::FRAC_PI_4),
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
            ("cos |x|<=pi/4", std::f32::consts::FRAC_PI_4),
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
            ("cos_checked |x|<=pi/4", std::f32::consts::FRAC_PI_4),
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
        // sinpi_ref/cospi_ref (above) already mirror sinpi/cospi's own
        // exact reduction in f64 rather than computing x*pi directly, so
        // the reference stays trustworthy over the *entire* f32 range
        // (f64's 2^53 exact-integer limit is far beyond f32's own 2^24,
        // and past that every f32 is already even by representability,
        // so no extra reference precision is even needed there) -- no
        // domain restriction needed. This full-range sweep is exactly
        // what caught a real bug (2026-07-09): the previous |x|<1e6
        // restriction happened to sit entirely below 2^22, silently
        // missing that the old magic-constant-based reduction (valid
        // only for |x|<=2^22, since it applied the trick directly to
        // unbounded raw x, unlike every other magic-round use in this
        // crate) was wrong for 2^22 < |x| < 2^24 despite the doc
        // comment's "exact out to f32::MAX" claim. Fixed with `x.round()`
        // (see sinpi's own doc comment); this sweep is what would have
        // caught it originally.
        let s = measure!(everywhere, sinpi, sinpi_ref);
        report("sinpi (all f32)", &s, t0);
        let s = measure!(everywhere, cospi, cospi_ref);
        report("cospi (all f32)", &s, t0);
    }
    if run("tanpi") {
        // Same full-range exactness as sinpi/cospi (tanpi is built
        // directly on their own reduction, see its doc comment) -- no
        // domain restriction needed for the reduction itself, but real
        // poles (cospi(x)==0, at half-integer x) mean huge ulp right at
        // those points is expected and harmless (matching cospi's own
        // "near a zero" caveat, just amplified by the division) rather
        // than a bug.
        let s = measure!(everywhere, tanpi, tanpi_ref);
        report("tanpi (all f32)", &s, t0);
    }
    if run("sinc") {
        // sinc(x) = sin(pi*x)/(pi*x) via sinpi, so sinc_ref reuses
        // sinpi_ref directly rather than a naive sin(pi*x)/(pi*x) (which
        // would reintroduce sinpi's own large-x imprecision into the
        // reference -- the same trap already hit once with sind_ref).
        // Restricted to |x|<1e6: well past that, sinc(x)'s true value is
        // already indistinguishable from 0 at f32 precision (|sinc(x)|
        // <= 1/(pi*|x|)), so ulp comparisons there measure noise near a
        // genuine zero, not real accuracy (same class of artifact as
        // cospi's own near-zero-crossing ulp blowup).
        let sinc_domain = |x: f32| x.abs() < 1e6;
        let s = measure!(sinc_domain, sinc, sinc_ref);
        report("sinc (|x|<1e6)", &s, t0);
    }
    if run("sind") {
        // sind/cosd's own exact-reduction range is ~4.7e7 (see their doc
        // comments, limited by 180.0's trailing zero bits, unlike
        // sinpi/cospi's full-f32-range exactness) -- widened to match
        // that documented boundary exactly (was |x|<1e6, leaving the
        // 1e6-4.7e7 "documented accurate but never actually tested" gap
        // unexercised -- the same shape of coverage hole that hid a real
        // bug in sinpi/cospi; checked directly with a throwaway scratch
        // harness using a proper reduction-based reference before
        // trusting this domain change, and confirmed clean: max ulp 2
        // throughout, no bug here, just closing the gap). This test's own
        // f64-based reference (`sind_ref`/`cosd_ref` above) stays
        // trustworthy well past this boundary too.
        let sind_domain = |x: f32| x.abs() < 4.7e7;
        let s = measure!(sind_domain, sind, sind_ref);
        report("sind (|x|<4.7e7)", &s, t0);
        let s = measure!(sind_domain, cosd, cosd_ref);
        report("cosd (|x|<4.7e7)", &s, t0);
        // tand: same domain as sind/cosd (built directly on their own
        // reduction); real poles at cosd(x)==0 give expected huge ulp
        // there, same caveat as tanpi above.
        let s = measure!(sind_domain, tand, tand_ref);
        report("tand (|x|<4.7e7)", &s, t0);
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
    if run("log2p1") {
        // log2(1+x) via log1p_u10(x)/ln(2), not the naive log2_u35(1.0+v):
        // forming 1.0+x directly in f64 hits the same cancellation trap
        // (for the many f32 x with |x| below f64's own ~2.22e-16 relative
        // precision near 1.0) that log1p's own real-function design exists
        // to avoid, just at a different, still-real-for-small-f32-x
        // threshold -- reusing sleef's own cancellation-safe log1p_u10
        // keeps the reference honest instead of reintroducing that bug on
        // the test side (the third time this session a naive f64
        // reference has needed the same fix, see sinc_ref/sind's own gap
        // check).
        let log2p1_ref = |v: F64xN| log1p_u10(v) / F64xN::splat(std::f64::consts::LN_2);
        let s = measure!(everywhere, log2p1, log2p1_ref);
        report("log2p1", &s, t0);
    }
    if run("log10p1") {
        // Same cancellation-safe log1p_u10 reference as log2p1 above, just
        // divided by ln(10) instead of ln(2).
        let log10p1_ref = |v: F64xN| log1p_u10(v) / F64xN::splat(std::f64::consts::LN_10);
        let s = measure!(everywhere, log10p1, log10p1_ref);
        report("log10p1", &s, t0);
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
    if run("exp_checked") {
        // Full range (backlog idea #18): properly saturates to 0/inf, so
        // unlike exp's own exp_domain restriction above, this gets the
        // whole f32 domain.
        let s = measure!(everywhere, exp_checked, exp_u10);
        report("exp_checked", &s, t0);
    }
    if run("exp_m1_over_x") {
        // Same inherited unchecked-exp2 domain as expm1 itself (see its
        // own doc comment).
        let exp_domain = |x: f32| (-126.0..128.0).contains(&(x * std::f32::consts::LOG2_E));
        // removable singularity at 0: (e^x-1)/x -> 1 as x -> 0.
        let exp_m1_over_x_ref = |v: F64xN| {
            let is_zero = v.simd_eq(F64xN::splat(0.0));
            is_zero.select(F64xN::splat(1.0), expm1_u10(v) / v)
        };
        let s = measure!(exp_domain, exp_m1_over_x, exp_m1_over_x_ref);
        report("exp_m1_over_x", &s, t0);
    }
    if run("exp2m1") {
        // 2^x - 1 via expm1_u10(x*ln2), not naive exp2_u35(v)-1.0: same
        // cancellation-trap reasoning as log2p1_ref just above, mirrored
        // (exp2_u35(v) rounds to exactly 1.0 in f64 for any |x| below
        // ~3.2e-16, which many f32 denormals/small-normals are, silently
        // erasing the true small-but-f32-representable result).
        // exp2m1 itself is total (inherits exp2_checked's own full clamp,
        // see its doc comment), so no domain restriction needed here,
        // unlike expm1 above.
        let exp2m1_ref = |v: F64xN| expm1_u10(v * F64xN::splat(std::f64::consts::LN_2));
        let s = measure!(everywhere, exp2m1, exp2m1_ref);
        report("exp2m1", &s, t0);
    }
    if run("exp10m1") {
        // Same expm1_u10-via-change-of-base reference as exp2m1 above, and
        // likewise total over exp10_checked's own clamped domain, so no
        // domain restriction needed.
        let exp10m1_ref = |v: F64xN| expm1_u10(v * F64xN::splat(std::f64::consts::LN_10));
        let s = measure!(everywhere, exp10m1, exp10m1_ref);
        report("exp10m1", &s, t0);
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
    if run("sinh_checked") {
        // Full range (backlog idea #85's fifth wave): properly saturates
        // to +-inf instead of sinh/cosh's own wraparound-to-garbage or
        // NaN-at-infinity, so unlike sinh/cosh's own sinh_domain
        // restriction above, this gets the whole f32 domain.
        let s = measure!(everywhere, sinh_checked, sinh_u35);
        report("sinh_checked", &s, t0);
        let s = measure!(everywhere, cosh_checked, cosh_u35);
        report("cosh_checked", &s, t0);
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
    if run("softplus") {
        // Restricted to |x|<80: softplus's own correction-term cutoff
        // (see its doc comment) creates a real, deliberate discontinuity
        // right around |x|=87 -- true value ~1.6e-38 on one side, exactly
        // 0.0 on the other, both equally "correct" in the sense that
        // neither is distinguishable from the other at any scale that
        // matters, but a raw ulp comparison right at that seam reports
        // millions of ulp for what's actually a sub-denormal-scale
        // difference (the same "near a value too small to matter"
        // artifact already documented for cospi elsewhere in this file).
        // |x|<80 stays comfortably clear of the seam on both sides.
        let softplus_domain = |x: f32| x.abs() < 80.0;
        let softplus_ref = |v: F64xN| {
            v.simd_max(F64xN::splat(0.0)) + log1p_u10(exp_u10(-v.abs()))
        };
        let s = measure!(softplus_domain, softplus, softplus_ref);
        report("softplus (|x|<80)", &s, t0);
    }
    if run("logaddexp") {
        // Same |x|<80 reasoning as softplus (its own doc comment) --
        // logaddexp shares the identical correction-term cutoff shape,
        // just on |a-b| instead of |x|.
        let logaddexp_domain = |a: f32, b: f32| a.abs() < 80.0 && b.abs() < 80.0;
        let logaddexp_ref = |a: F64xN, b: F64xN| {
            let m = a.simd_max(b);
            let d = (a - b).abs();
            m + log1p_u10(exp_u10(-d))
        };
        let s = fuzz2(TWOARG_SAMPLES, logaddexp_domain, logaddexp, logaddexp_ref);
        report("logaddexp (|a|,|b|<80)", &s, t0);
    }
    if run("gelu") {
        // x * Phi(x) via erfc_u15 (see gelu's own doc comment for why not
        // erf_u10 -- 1+erf(z) cancels toward 0 for negative x). Same
        // |erfc's argument|<=10 domain restriction as erfc's own block
        // above, translated through gelu's `-x/sqrt2` argument (so
        // |x|<=10*sqrt2); outside it, gelu still returns a sane saturated
        // value via erfc's own saturation (see gelu's doc comment), just
        // not one this fuzz screen claims ulp accuracy for.
        let gelu_domain = |x: f32| x.abs() <= 10.0 * std::f32::consts::SQRT_2;
        let gelu_ref = |v: F64xN| {
            v * F64xN::splat(0.5) * erfc_u15(-v * F64xN::splat(std::f64::consts::FRAC_1_SQRT_2))
        };
        let s = measure!(gelu_domain, gelu, gelu_ref);
        report("gelu", &s, t0);
    }
    if run("silu") {
        // x * sigmoid(x), same sigmoid_ref/domain as sigmoid's own block
        // above (this crate's sigmoid is only calibrated accurate inside
        // that exp(-x)-doesn't-overflow range; outside it, silu still
        // returns a sane saturated value via sigmoid's own clamp, just not
        // one this fuzz screen claims ulp accuracy for).
        let silu_domain = |x: f32| {
            let e = x * std::f32::consts::LOG2_E;
            e > -126.0 && e < 126.0
        };
        let silu_ref = |v: F64xN| v / (F64xN::splat(1.0) + exp_u10(-v));
        let s = measure!(silu_domain, silu, silu_ref);
        report("silu", &s, t0);
    }
    if run("softsign") {
        // `x = +-inf` excluded: same reference-side `inf/inf` indeterminate
        // form as gelu's `-inf` case above (softsign itself special-cases
        // it; edgecheck.rs pins the actual behavior).
        let finite = |x: f32| x.is_finite();
        let softsign_ref = |v: F64xN| v / (F64xN::splat(1.0) + v.abs());
        let s = measure!(finite, softsign, softsign_ref);
        report("softsign", &s, t0);
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
        let s = measure!(everywhere, atan_latency, atan_u35);
        report("atan_latency", &s, t0);
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
        let s = measure!(erfc_domain, erfc_accurate, erfc_u15);
        report("erfc_accurate", &s, t0);
        // erfcx_ref: no sleef erfcx bucket, so compose exp(x^2)*erfc_u15(x)
        // directly in f64 -- safe over this domain (x^2 <= 100 is nowhere
        // near f64's own ~709 exp overflow point) and multiplication
        // doesn't lose relative precision the way addition/subtraction
        // would, so a tiny erfc(x) times a huge exp(x^2) is still an
        // accurate f64 product.
        let erfcx_ref = |v: F64xN| exp_u10(v * v) * erfc_u15(v);
        let s = measure!(erfc_domain, erfcx, erfcx_ref);
        report("erfcx", &s, t0);
        let s = measure!(erfc_domain, erfcx_accurate, erfcx_ref);
        report("erfcx_accurate", &s, t0);
        // erfcx_checked's own wider domain, extended past erfcx's |x|<=10
        // fit boundary (see its doc comment) -- still safely inside where
        // erfcx_ref's f64 `exp_u10(v*v)` doesn't itself overflow (v*v <=
        // 400, nowhere near f64's ~709 exp overflow point), so the same
        // reference stays trustworthy this far out; the asymptotic tail's
        // own verified range (up to x=200, see erfcx_checked's doc
        // comment) was checked separately against scipy.special.erfcx
        // since sleef has no f64 exp headroom left to compose a reference
        // that far.
        let erfcx_checked_domain = |x: f32| x.abs() <= 20.0;
        let s = measure!(erfcx_checked_domain, erfcx_checked, erfcx_ref);
        report("erfcx_checked", &s, t0);
    }

    // two-argument functions: fuzz-only (exhaustive over 2^64 pairs isn't
    // feasible), smaller sample count since each trial needs two RNG draws.
    const TWOARG_SAMPLES: u64 = 10_000_000;
    if run("atan2") {
        let s = fuzz2(TWOARG_SAMPLES, |_, _| true, atan2, atan2_u35);
        report("atan2", &s, t0);
        let s = fuzz2(TWOARG_SAMPLES, |_, _| true, |y: f32, x: f32| y.atan2(x), atan2_u35);
        report("std atan2", &s, t0);
        // atan2_unchecked's documented contract: x != 0.0, not both infinite.
        let atan2_domain = |_: f32, x: f32| x != 0.0;
        let s = fuzz2(TWOARG_SAMPLES, atan2_domain, atan2_unchecked, atan2_u35);
        report("atan2_unchecked (+)", &s, t0);
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
        // hypot_unchecked's documented contract: x, y both finite. Same
        // overflow/underflow-avoidance domain restriction as hypot itself.
        let s = fuzz2(TWOARG_SAMPLES, hypot_domain, hypot_unchecked, hypot_u35);
        report("hypot_unchecked (+)", &s, t0);
        // hypot_checked has no overflow/underflow tradeoff to work around
        // (that's the whole point), so it gets the full domain -- every
        // finite magnitude, zero, inf, and nan.
        let s = fuzz2(TWOARG_SAMPLES, |_, _| true, hypot_checked, hypot_u35);
        report("hypot_checked", &s, t0);
    }
    if run("rhypot") {
        // Same overflow/underflow tradeoff as hypot/hypot_unchecked (see
        // hypot's own domain comment above) -- rhypot shares the identical
        // fma(x,x,y*y) core.
        let hypot_domain = |x: f32, y: f32| {
            let ok = |v: f32| v == 0.0 || (v.abs() > 1e-15 && v.abs() < 1e18);
            ok(x) && ok(y)
        };
        let rhypot_ref = |v: F64xN, w: F64xN| F64xN::splat(1.0) / hypot_u35(v, w);
        let s = fuzz2(TWOARG_SAMPLES, hypot_domain, rhypot, rhypot_ref);
        report("rhypot", &s, t0);
    }
    if run("diff_of_products") {
        // Same overflow tradeoff as hypot_domain above, applied to the two
        // *products* rather than the raw operands directly: fuzz4 has no
        // domain-filter parameter (unlike fuzz2/measure!), so the bound is
        // baked into the reference closure itself -- reject samples where
        // either product would leave f32's representable range, matching
        // diff_of_products' own doc comment.
        let bounded = |a: f64, b: f64, c: f64, d: f64| -> f64 {
            let ok = |p: f64| p.abs() < 1e30;
            if ok(a * b) && ok(c * d) { a * b - c * d } else { f64::NAN }
        };
        let s = fuzz4(QUICK_SAMPLES, diff_of_products, bounded);
        report("diff_of_products", &s, t0);
    }
    if run("cross2") {
        let bounded = |ax: f64, ay: f64, bx: f64, by: f64| -> f64 {
            let ok = |p: f64| p.abs() < 1e30;
            if ok(ax * by) && ok(ay * bx) { ax * by - ay * bx } else { f64::NAN }
        };
        let s = fuzz4(QUICK_SAMPLES, cross2, bounded);
        report("cross2", &s, t0);
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
        let pown_sweep = |f: &dyn Fn(f32, i32) -> f32, lo: i32, hi: i32, label: &str| {
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
                let got = f(x, n);
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
        pown_sweep(&pown, -8, 8, "pown (|n|<=8)");
        pown_sweep(&pown, -64, 64, "pown (|n|<=64)");
        // pown_small's own contract (|n| <= 255) -- bit-identical to pown
        // whenever both are in-domain, so this also doubles as a
        // regression check on the two shared |n|<=8/|n|<=64 buckets.
        pown_sweep(&pown_small, -8, 8, "pown_small (|n|<=8)");
        pown_sweep(&pown_small, -64, 64, "pown_small (|n|<=64)");
        pown_sweep(&pown_small, -255, 255, "pown_small (|n|<=255)");
        // pown_small_accurate: Df32-compensated squaring chain (see its
        // own doc comment), same |n|<=255 contract as pown_small -- real
        // improvement expected (avg/max ulp roughly halved on the
        // dense-|n| bucket), not a full fix.
        pown_sweep(&pown_small_accurate, -8, 8, "pown_small_accurate (|n|<=8)");
        pown_sweep(&pown_small_accurate, -64, 64, "pown_small_accurate (|n|<=64)");
        pown_sweep(&pown_small_accurate, -255, 255, "pown_small_accurate (|n|<=255)");
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
    if run("powf_pos") {
        // powf_pos's own contract: x >= 0.0 (see its doc comment for the
        // one excluded value, x == -0.0, negligible density in a random
        // fuzz). Same exponent-range restriction as powf's own sweep.
        let pow_domain =
            |x: f32, y: f32| x >= 0.0 && (-126.0..128.0).contains(&(x.log2() * y));
        let s = fuzz2(TWOARG_SAMPLES, pow_domain, powf_pos, pow_u10);
        report("powf_pos", &s, t0);
    }
    if run("powf_unchecked") {
        // powf_unchecked's own contract: x positive/normal/finite (same as
        // log_2_unchecked), y != 0.0. Narrower than powf's own domain above
        // (no negative x), same exponent-range restriction.
        let pow_domain = |x: f32, y: f32| {
            x >= f32::MIN_POSITIVE && x.is_finite() && y != 0.0 && (-126.0..128.0).contains(&(x.log2() * y))
        };
        let s = fuzz2(TWOARG_SAMPLES, pow_domain, powf_unchecked, pow_u10);
        report("powf_unchecked (+)", &s, t0);
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
    if run("powf_checked_unchecked") {
        // powf_checked_unchecked's own contract: x positive/normal/finite
        // (same as powf_unchecked), y != 0.0.
        let pow_domain = |x: f32, y: f32| {
            x >= f32::MIN_POSITIVE && x.is_finite() && y != 0.0 && (-126.0..128.0).contains(&(x.log2() * y))
        };
        let s = fuzz2(TWOARG_SAMPLES, pow_domain, powf_checked_unchecked, pow_u10);
        report("powf_checked_unchecked (+)", &s, t0);
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
    if run("remainder_unchecked") {
        // remainder_unchecked's own contract: x != 0.0, y finite. Same
        // reliable-range/near-tie exclusions as remainder's own sweep above.
        let remainder_domain =
            |x: f32, y: f32| x != 0.0 && y != 0.0 && y.is_finite() && (x / y).abs() < 1000.0 && !near_tie(x, y);
        let s = fuzz2(TWOARG_SAMPLES, remainder_domain, remainder_unchecked, remainder_ref);
        report("remainder_unchecked (+)", &s, t0);
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
    if run("remainder_wide") {
        // remainder_wide extends remainder_checked's own correction past
        // the 2^24 cliff (see its own doc comment) using an exact Df32
        // residual instead of a single fma, verified by hand to hold
        // cleanly up to |x/y| ~ 2^48 (2^48 == 2.81e14) -- bound comfortably
        // under that, matching remainder_checked's own "stay well clear of
        // the cliff" convention. sleef's remainder_ref is a real IEEE754
        // remainder implementation (its own internal argument reduction,
        // not a naive single-f64-pass formula), so it stays trustworthy as
        // a reference at this magnitude, unlike a naive `xf - q*yf` f64
        // reference would be (confirmed by hand with an arbitrary-
        // precision Python check during development -- plain f64
        // arithmetic loses precision once q*y needs more than f64's own
        // 52 mantissa bits, which happens well before this domain's edge).
        let remainder_domain =
            |x: f32, y: f32| y != 0.0 && (x / y).abs() < 2.0e14 && !near_tie(x, y);
        let s = fuzz2(TWOARG_SAMPLES, remainder_domain, remainder_wide, remainder_ref);
        report("remainder_wide", &s, t0);
    }
    if run("remainder_ieee") {
        // remainder_ieee rounds q ties-to-even instead of remainder's own
        // ties-away, matching sleef's true-IEEE754 reference exactly at
        // ties (verified separately in edgecheck.rs with pinned exact-tie
        // cases, e.g. remainder_ieee(5,2)==1.0 vs remainder(5,2)==-1.0).
        // This sweep still excludes near_tie, same as remainder/
        // remainder_checked above: that exclusion is about a different,
        // already-known issue (x/y's own division rounding flipping which
        // *integer* q lands on near, but not exactly at, a tie -- the
        // problem remainder_checked exists to fix), which affects
        // remainder_ieee identically to remainder since both share the
        // same plain `x/y` division, only differing in the final
        // rounding-mode convention.
        let remainder_domain =
            |x: f32, y: f32| y != 0.0 && (x / y).abs() < 1000.0 && !near_tie(x, y);
        let s = fuzz2(TWOARG_SAMPLES, remainder_domain, remainder_ieee, remainder_ref);
        report("remainder_ieee", &s, t0);
    }
    if run("fmod") {
        // fmod's own failure mode (see its doc comment) is the
        // truncation analog of remainder's near_tie exclusion, but at
        // integer boundaries instead of half-integer ones (`.trunc()`'s
        // decision changes at integers, not half-integers) -- exclude
        // the same way, `near_int` instead of `near_tie`.
        let near_int = |x: f32, y: f32| {
            let f = (x / y).fract().abs();
            f < 1e-4 || f > 1.0 - 1e-4
        };
        let fmod_domain =
            |x: f32, y: f32| y != 0.0 && (x / y).abs() < 1000.0 && !near_int(x, y);
        let fmod_ref = |a: F64xN, b: F64xN| a % b;
        let s = fuzz2(TWOARG_SAMPLES, fmod_domain, fmod, fmod_ref);
        report("fmod", &s, t0);
        let fmod_domain_unchecked = |x: f32, y: f32| {
            x != 0.0 && y != 0.0 && y.is_finite() && (x / y).abs() < 1000.0 && !near_int(x, y)
        };
        let s = fuzz2(TWOARG_SAMPLES, fmod_domain_unchecked, fmod_unchecked, fmod_ref);
        report("fmod_unchecked (+)", &s, t0);
        // fmod_checked: no near_int exclusion -- that's exactly the
        // domain it fixes (see its own doc comment), so this doubles as
        // a standing regression check on the off-by-a-whole-y bug.
        let fmod_checked_domain = |x: f32, y: f32| y != 0.0 && (x / y).abs() < 1000.0;
        let s = fuzz2(TWOARG_SAMPLES, fmod_checked_domain, fmod_checked, fmod_ref);
        report("fmod_checked", &s, t0);
    }

    println!("total: {:.2}s", t0.elapsed().as_secs_f64());
}
