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
use std::simd::cmp::{SimdPartialEq, SimdPartialOrd};
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
        //
        // Deliberately keeps this `round(x-0.5)` shape rather than
        // mirroring the real cospi's newer `round(x)` one: in f64 both
        // are exact for an f32 input (24 bits of input, ~53 available,
        // so the rounding that made the shape matter in f32 cannot
        // happen here), and this one signs its zeros the way `tanpi_ref`
        // below needs -- `sinpi_ref/cospi_ref` there lands on a constant
        // `-inf` at every pole only because these zeros alternate. The
        // real `cospi`'s `+0.0`-everywhere zeros are invisible to
        // `ulp_diff`, which ranks `+0.0` and `-0.0` equal.
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
// sinc_unnormalized(x) = sin(x)/x, radians (backlog idea #131). Unlike
// sinc_ref (which routes through sinpi_ref's own exact round(x)
// reduction, never handing sin_u35 an extreme argument directly), there
// is no analogous exact-reduction trick for plain radians -- sleef's
// sin_u35 itself panics (internal table index out of bounds) on a
// sufficiently huge finite argument, so this also has to guard *large*
// |v|, not just non-finite v the way trig_safe does. Fine here: the
// accuracy.rs `sinc_domain` restriction (|x|<1e6) already excludes
// these from scoring, this is purely about not crashing the *call*
// on the unrestricted lanes "thorough" mode still feeds through.
fn sinc_unnormalized_ref(v: F64xN) -> F64xN {
    let safe = v.abs().simd_lt(F64xN::splat(1e15)) & v.is_finite();
    let safe_v = safe.select(v, F64xN::splat(1.0));
    let normal = sin_u35(safe_v) / safe_v;
    let is_zero = v.simd_eq(F64xN::splat(0.0));
    is_zero.select(F64xN::splat(1.0), safe.select(normal, F64xN::splat(f64::NAN)))
}
// xlogy/xlog1py (backlog idea #84): x==0 overrides to 0 regardless of y
// (matching scipy.special.xlogy's own convention, see the real
// functions' doc comments), otherwise plain x*ln(y)/x*ln(1+y). log_u35
// itself already returns NaN for y<=0 and propagates NaN/inf correctly,
// so no extra domain guarding needed beyond the x==0 override.
fn xlogy_ref(x: F64xN, y: F64xN) -> F64xN {
    let is_zero = x.simd_eq(F64xN::splat(0.0));
    is_zero.select(F64xN::splat(0.0), x * log_u35(y))
}
fn xlog1py_ref(x: F64xN, y: F64xN) -> F64xN {
    let is_zero = x.simd_eq(F64xN::splat(0.0));
    is_zero.select(F64xN::splat(0.0), x * log1p_u10(y))
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

/// ULP distance along the monotonic ordering of f32 bit patterns.
///
/// **NaN vs NaN is always exactly 0, whatever the two bit patterns are.**
/// Payload, sign and quiet bits carry no numeric meaning, so a difference
/// there is not an accuracy difference -- and `thorough` mode sweeps all
/// 2^24 NaN payloads, so scoring them would swamp the max-ulp column with
/// something that isn't error. Exactly one side being NaN is still a
/// maximal error. Any bespoke scoring harness must copy this rule; the
/// one deliberate exception in this repo is `worst_corpus.rs`, which is a
/// bit-exactness regression gate rather than an accuracy test.
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
    if run("pow_3_2") {
        let pos = |x: f32| x >= 0.0;
        let pow_3_2_ref = |v: F64xN| v * v.sqrt();
        let s = measure!(pos, pow_3_2, pow_3_2_ref);
        report("pow_3_2", &s, t0);
    }
    if run("pow_2_3") {
        let pow_2_3_ref = |v: F64xN| {
            let c = cbrt_u35(v);
            c * c
        };
        let s = measure!(everywhere, pow_2_3, pow_2_3_ref);
        report("pow_2_3", &s, t0);
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
        // std on the same |x|<=1e6 restriction the readme's row quotes, so
        // both columns of that row come from the same input set.
        let s = measure!(|x: f32| x.abs() <= 1e6, |x: f32| x.sin(), sin_ref);
        report("std sin |x|<=1e6", &s, t0);
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
            ("sin_checked [1e13,1e14)", 1e13, 1e14),
            ("sin_checked [1e14,1e15)", 1e14, 1e15),
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
        // std on the same |x|<=1e6 restriction the readme's row quotes.
        let s = measure!(|x: f32| x.abs() <= 1e6, |x: f32| x.cos(), cos_ref);
        report("std cos |x|<=1e6", &s, t0);
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
            ("cos_checked [1e13,1e14)", 1e13, 1e14),
            ("cos_checked [1e14,1e15)", 1e14, 1e15),
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
        // sinpi_unchecked's one documented difference from sinpi (wrong
        // sign of zero at x=-0.0 only) is invisible to ulp_diff here --
        // it treats +0.0/-0.0 as equal (see its own `ord` helper) -- so
        // `everywhere` is safe; edgecheck.rs's exact-bits check is what
        // actually pins the -0.0 exception.
        let s = measure!(everywhere, sinpi_unchecked, sinpi_ref);
        report("sinpi_unchecked (+)", &s, t0);
        let s = measure!(everywhere, cospi, cospi_ref);
        report("cospi (all f32)", &s, t0);
    }
    if run("tanpi") {
        // Same full-range exactness as sinpi/cospi (tanpi is built
        // directly on their own reduction, see its doc comment) -- no
        // domain restriction needed for the reduction itself, but real
        // poles (cospi(x)==0, at half-integer x) mean a huge ulp right at
        // those points would be a pole artifact rather than a bug. (Not
        // a licence to assume it: `tan_core` pins the poles to an exact
        // `-inf` the reference matches, so this measures clean anyway.
        // "Near a true zero/pole, ulp isn't meaningful" only holds when
        // the reduction is exact there -- assuming it for `cospi` is
        // what hid a real 8.7e8-ulp bug until 2026-07-30.)
        let s = measure!(everywhere, tanpi, tanpi_ref);
        report("tanpi (all f32)", &s, t0);
    }
    if run("2pi") {
        // sin2pi/cos2pi/tan2pi (backlog idea #122): 2*v is exact in f64
        // for any f32 v (nowhere near f64's own overflow), so the
        // reference just doubles before sinpi_ref/etc, mirroring the
        // real function exactly. Restricted to |x|<f32::MAX/2 -- past
        // that the real function's own `2.0*x` overflows to +-inf, a
        // documented gap (see sin2pi's own doc comment), not something
        // to score here.
        let half_max = |x: f32| x.abs() < f32::MAX / 2.0;
        let sin2pi_ref = |v: F64xN| sinpi_ref(v * F64xN::splat(2.0));
        let s = measure!(half_max, sin2pi, sin2pi_ref);
        report("sin2pi", &s, t0);
        let cos2pi_ref = |v: F64xN| cospi_ref(v * F64xN::splat(2.0));
        let s = measure!(half_max, cos2pi, cos2pi_ref);
        report("cos2pi", &s, t0);
        let tan2pi_ref = |v: F64xN| tanpi_ref(v * F64xN::splat(2.0));
        let s = measure!(half_max, tan2pi, tan2pi_ref);
        report("tan2pi", &s, t0);
    }
    if run("sinc") {
        // sinc(x) = sin(pi*x)/(pi*x) via sinpi, so sinc_ref reuses
        // sinpi_ref directly rather than a naive sin(pi*x)/(pi*x) (which
        // would reintroduce sinpi's own large-x imprecision into the
        // reference -- the same trap already hit once with sind_ref).
        // Restricted to |x|<1e6: well past that, sinc(x)'s true value is
        // already indistinguishable from 0 at f32 precision (|sinc(x)|
        // <= 1/(pi*|x|)), so ulp comparisons there measure noise near a
        // value too small for f32 to resolve, not real accuracy. (The
        // true value being unresolvable is the claim; a *crossing* alone
        // would not justify this -- see tanpi's note above.)
        let sinc_domain = |x: f32| x.abs() < 1e6;
        let s = measure!(sinc_domain, sinc, sinc_ref);
        report("sinc (|x|<1e6)", &s, t0);
        // sinc_unnormalized (backlog idea #131): same near-zero-crossing
        // caveat as sinc itself, same |x|<1e6 restriction.
        let s = measure!(sinc_domain, sinc_unnormalized, sinc_unnormalized_ref);
        report("sinc_unnormalized (|x|<1e6)", &s, t0);
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
        let s = measure!(sind_domain, sind_unchecked, sind_ref);
        report("sind_unchecked (+)", &s, t0);
        let s = measure!(sind_domain, cosd, cosd_ref);
        report("cosd (|x|<4.7e7)", &s, t0);
        let s = measure!(sind_domain, cosd_unchecked, cosd_ref);
        report("cosd_unchecked (+)", &s, t0);
        // tand: same domain as sind/cosd (built directly on their own
        // reduction); real poles at cosd(x)==0 give expected huge ulp
        // there, same caveat as tanpi above.
        let s = measure!(sind_domain, tand, tand_ref);
        report("tand (|x|<4.7e7)", &s, t0);
        let s = measure!(sind_domain, tand_unchecked, tand_ref);
        report("tand_unchecked (+)", &s, t0);
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
        // log1pmx (backlog idea #145): direct f64 log1p_u10(v)-v cancels
        // for tiny v the same way the naive f32 form does, just at f64's
        // own (much smaller) precision floor -- rationalized for |v| below
        // that floor via the leading-order term (higher Taylor terms are
        // utterly negligible there), same "reference must itself avoid
        // the cancellation trap" precedent as sqrt1pm1's own reference.
        let log1pmx_ref = |v: F64xN| {
            let tiny = v.abs().simd_lt(F64xN::splat(1e-6));
            let small_ref = v * v * F64xN::splat(-0.5);
            let big_ref = log1p_u10(v) - v;
            let is_pos_inf = v.simd_eq(F64xN::splat(f64::INFINITY));
            let normal = tiny.select(small_ref, big_ref);
            is_pos_inf.select(F64xN::splat(f64::NEG_INFINITY), normal)
        };
        let s = measure!(everywhere, log1pmx, log1pmx_ref);
        report("log1pmx", &s, t0);
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
        // exp_narrow's own documented domain (backlog ideas #23/#112):
        // k=round(x*log2(e)) in [-126,127], a single exponent field's own
        // range -- slightly narrower than exp_domain above.
        let exp_narrow_domain = |x: f32| (-87.68311..=88.37627).contains(&x);
        let s = measure!(exp_narrow_domain, exp_narrow, exp_u10);
        report("exp_narrow", &s, t0);
        let s = measure!(exp_domain, expm1, expm1_u10);
        report("expm1", &s, t0);
        let s = measure!(exp_domain, |x: f32| x.exp_m1(), expm1_u10);
        report("std expm1", &s, t0);
        // expm1_narrow's own documented domain (backlog idea #201, same
        // mechanism as exp_narrow).
        let exp_narrow_domain = |x: f32| (-87.68311..=88.37627).contains(&x);
        let s = measure!(exp_narrow_domain, expm1_narrow, expm1_u10);
        report("expm1_narrow", &s, t0);
        // expm1_checked (backlog idea #111) is total, so it gets the whole
        // domain rather than exp_domain -- the same `everywhere` treatment
        // exp_checked gets below, and the point of the function. std is
        // measured over that same full domain too, so the readme's row can
        // compare like with like (the "std expm1" run above is scored over
        // exp_domain only, so its numbers are not interchangeable).
        let s = measure!(everywhere, expm1_checked, expm1_u10);
        report("expm1_checked", &s, t0);
        let s = measure!(everywhere, |x: f32| x.exp_m1(), expm1_u10);
        report("std expm1 (everywhere)", &s, t0);
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
        // exp_m1_over_x_narrow's own documented domain (backlog idea
        // #201, same mechanism as exp_narrow/expm1_narrow).
        let exp_narrow_domain = |x: f32| (-87.68311..=88.37627).contains(&x);
        let s = measure!(exp_narrow_domain, exp_m1_over_x_narrow, exp_m1_over_x_ref);
        report("exp_m1_over_x_narrow", &s, t0);
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
        // sinh_narrow/cosh_narrow (backlog idea #201): same domain
        // concern as sinh_domain above (both +k/-k must stay in exp2's
        // unchecked range), just tighter in x -- sinh_domain already
        // covers it.
        let s = measure!(sinh_domain, sinh_narrow, sinh_u35);
        report("sinh_narrow", &s, t0);
        let s = measure!(sinh_domain, cosh_narrow, cosh_u35);
        report("cosh_narrow", &s, t0);
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
        // coshm1 (backlog idea #144): same half-angle identity as the
        // reference, computed in f64 to avoid the reference itself
        // cancelling for small x.
        let coshm1_ref = |v: F64xN| {
            let s = sinh_u35(v * F64xN::splat(0.5));
            F64xN::splat(2.0) * s * s
        };
        let s = measure!(everywhere, coshm1, coshm1_ref);
        report("coshm1", &s, t0);
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
        // tanh_grad (backlog idea #150): `1 - tanh_u35(v)^2` in f64 was
        // tried first and rejected as a *reference* bug, not a
        // tanh_grad bug -- tanh_u35 itself correctly rounds to exactly
        // 1.0 in f64 well before this domain's edge, so the subtraction
        // silently gives a wrong 0.0 reference over a wide band (same
        // "pushed out to a larger |x| in f64" trap as sigmoid_ref's own
        // doc comment below). The stable `4q/(1+q)^2` form (`q =
        // exp(-2|x|)`, tanh_grad's own doc comment derives it) has no
        // such cancellation at any scale.
        let tanh_grad_ref = |v: F64xN| {
            let q = exp_u10(v.abs() * F64xN::splat(-2.0));
            F64xN::splat(4.0) * q / ((F64xN::splat(1.0) + q) * (F64xN::splat(1.0) + q))
        };
        let s = measure!(tanh_domain, tanh_grad, tanh_grad_ref);
        report("tanh_grad", &s, t0);
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
        // sigmoid_grad (backlog idea #150): `sigmoid_ref(v)*(1-sigmoid_ref(v))`
        // was tried first and rejected as a *reference* bug (same "pushed
        // out to a larger |x| in f64" trap as sigmoid_ref's own doc
        // comment above) -- sigmoid_ref itself rounds to exactly 1.0 in
        // f64 well before this domain's edge, silently giving a wrong 0.0
        // reference. The stable `e/(1+e)^2` form (sigmoid_grad's own doc
        // comment derives it) has no such cancellation.
        let sigmoid_grad_ref = |v: F64xN| {
            let e = exp_u10(-v);
            e / ((F64xN::splat(1.0) + e) * (F64xN::splat(1.0) + e))
        };
        let s = measure!(sigmoid_domain, sigmoid_grad, sigmoid_grad_ref);
        report("sigmoid_grad", &s, t0);
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
        // artifact documented for sinc above).
        // |x|<80 stays comfortably clear of the seam on both sides.
        let softplus_domain = |x: f32| x.abs() < 80.0;
        let softplus_ref = |v: F64xN| {
            v.simd_max(F64xN::splat(0.0)) + log1p_u10(exp_u10(-v.abs()))
        };
        let s = measure!(softplus_domain, softplus, softplus_ref);
        report("softplus (|x|<80)", &s, t0);
    }
    if run("logsigmoid") {
        // Same |x|<80 reasoning as softplus's own doc comment (this is
        // -softplus(-x), so the identical seam sits at the same |x|=87).
        let logsigmoid_domain = |x: f32| x.abs() < 80.0;
        let logsigmoid_ref = |v: F64xN| {
            -((-v).simd_max(F64xN::splat(0.0)) + log1p_u10(exp_u10(-v.abs())))
        };
        let s = measure!(logsigmoid_domain, logsigmoid, logsigmoid_ref);
        report("logsigmoid (|x|<80)", &s, t0);
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
        // asind (backlog idea #123): plain composite -- a rescaled-
        // coefficient fold was tried and measured worse, see its own
        // doc comment.
        let rad_to_deg = 180.0 / std::f64::consts::PI;
        let asind_ref = |v: F64xN| asin_u35(v) * F64xN::splat(rad_to_deg);
        let s = measure!(everywhere, asind, asind_ref);
        report("asind", &s, t0);
        // asinpi (backlog idea #85): folded dedicated coefficients,
        // compared directly against the naive asin(x)/PI post-multiply
        // to see whether the fold is actually worth it here (idea #123's
        // own RAD_TO_DEG fold measured worse for asind -- checking
        // separately, not assuming the same verdict transfers).
        let inv_pi = 1.0 / std::f64::consts::PI;
        let asinpi_ref = |v: F64xN| asin_u35(v) * F64xN::splat(inv_pi);
        let s = measure!(everywhere, asinpi, asinpi_ref);
        report("asinpi", &s, t0);
        let inv_pi_f32 = 1.0f32 / std::f32::consts::PI;
        let s = measure!(everywhere, move |x: f32| asin(x) * inv_pi_f32, asinpi_ref);
        report("asin(x)/PI (naive)", &s, t0);
    }
    if run("acos") {
        let s = measure!(everywhere, acos, acos_u35);
        report("acos", &s, t0);
        let s = measure!(everywhere, |x: f32| x.acos(), acos_u35);
        report("std acos", &s, t0);
        let rad_to_deg = 180.0 / std::f64::consts::PI;
        let acosd_ref = |v: F64xN| acos_u35(v) * F64xN::splat(rad_to_deg);
        let s = measure!(everywhere, acosd, acosd_ref);
        report("acosd", &s, t0);
        // acospi (backlog idea #85): folded vs naive, tested separately
        // per asinpi's own "each site needs its own measurement" finding.
        let inv_pi = 1.0 / std::f64::consts::PI;
        let acospi_ref = |v: F64xN| acos_u35(v) * F64xN::splat(inv_pi);
        let s = measure!(everywhere, acospi, acospi_ref);
        report("acospi", &s, t0);
        let inv_pi_f32 = 1.0f32 / std::f32::consts::PI;
        let s = measure!(everywhere, move |x: f32| acos(x) * inv_pi_f32, acospi_ref);
        report("acos(x)/PI (naive)", &s, t0);
    }
    if run("atan") {
        let s = measure!(everywhere, atan, atan_u35);
        report("atan", &s, t0);
        let s = measure!(everywhere, |x: f32| x.atan(), atan_u35);
        report("std atan", &s, t0);
        let s = measure!(everywhere, atan_latency, atan_u35);
        report("atan_latency", &s, t0);
        // atan_bounded's own documented contract (backlog idea #61).
        let bounded_domain = |x: f32| x.abs() <= 1.0;
        let s = measure!(bounded_domain, atan_bounded, atan_u35);
        report("atan_bounded", &s, t0);
        let rad_to_deg = 180.0 / std::f64::consts::PI;
        let atand_ref = |v: F64xN| atan_u35(v) * F64xN::splat(rad_to_deg);
        let s = measure!(everywhere, atand, atand_ref);
        report("atand", &s, t0);
        // atanpi (backlog idea #85): plain composite -- see its own doc
        // comment for why a rescaled-coefficient fold isn't kept here.
        let inv_pi = 1.0 / std::f64::consts::PI;
        let atanpi_ref = |v: F64xN| atan_u35(v) * F64xN::splat(inv_pi);
        let s = measure!(everywhere, atanpi, atanpi_ref);
        report("atanpi", &s, t0);
    }
    if run("tan") {
        let tan_domain = |x: f32| x.abs() < (1u32 << 22) as f32 * std::f32::consts::PI;
        let s = measure!(tan_domain, tan, tan_ref);
        report("tan (in-domain)", &s, t0);
        let s = measure!(tan_domain, |x: f32| x.tan(), tan_ref);
        report("std tan (in-domain)", &s, t0);
    }
    if run("tan_checked") {
        // Full range gradual degradation, like sin_checked/cos_checked
        // themselves (idea #48) -- no domain restriction *for validity*,
        // but expect the reported avg/max ulp to look alarming: tan has
        // a genuine pole every pi, and at large |x| those poles sit
        // closer together than the local float spacing, so *any*
        // correctly-behaving tan implementation shows unbounded relative
        // error near them (e.g. x=-4.4230258e15 sits at x/pi ==
        // -1407892830220377.5, essentially exactly a pole -- verified
        // cos_checked(x)=-7.88e-6, correctly near zero, not a bug). Same
        // "ulp isn't meaningful near a true zero/pole" class of artifact
        // as cosh's own near-zero case elsewhere in this crate, just at
        // infinity instead of zero -- and here, unlike `cospi`'s former
        // version of this excuse, the reduction genuinely has run out of
        // bits by that magnitude, which is what makes it an excuse.
        let s = measure!(everywhere, tan_checked, tan_ref);
        report("tan_checked", &s, t0);
    }
    if run("erf") {
        // erf_poly used to be evaluated unbounded on |x|, which was wrong
        // (not just imprecise) well before this bound -- see erf's doc
        // comment. Fixed now, so this measures the whole domain.
        let s = measure!(everywhere, erf, erf_u10);
        report("erf", &s, t0);
        // erfc/erfcx are now accurate over the whole real line (no fit
        // clamp anywhere -- see erfcx_pos), but the f64 reference is the
        // limit: sleef's erfc_u15 underflows to 0 well before f32's own
        // domain ends, and erfcx_ref below composes exp(x^2), so both are
        // only trustworthy while x^2 stays far from f64's ~709 exp
        // overflow. |x| <= 10 keeps erfc_u15 above its own underflow, and
        // is already past x ~ 10.05 where the true erfc rounds to exactly
        // 0.0 (2.0 for x < 0) and this crate returns exactly that.
        let erfc_domain = |x: f32| x.abs() <= 10.0;
        let s = measure!(erfc_domain, erfc, erfc_u15);
        report("erfc", &s, t0);
        // erfcx_ref: no sleef erfcx bucket, so compose exp(x^2)*erfc_u15(x)
        // directly in f64 -- safe over this domain (x^2 <= 400 is nowhere
        // near f64's own ~709 exp overflow point) and multiplication
        // doesn't lose relative precision the way addition/subtraction
        // would, so a tiny erfc(x) times a huge exp(x^2) is still an
        // accurate f64 product.
        let erfcx_ref = |v: F64xN| exp_u10(v * v) * erfc_u15(v);
        let s = measure!(erfc_domain, erfcx, erfcx_ref);
        report("erfcx", &s, t0);
        // erfcx's own wider domain: past |x|=10 the old rational froze,
        // so this range used to belong to a separate `erfcx_checked` tier;
        // the reciprocal-variable fit covers it directly now.
        let erfcx_wide = |x: f32| x.abs() <= 20.0;
        let s = measure!(erfcx_wide, erfcx, erfcx_ref);
        report("erfcx (|x|<=20)", &s, t0);
        // The rest of the domain, all the way to f32::MAX. `exp(x^2)` is
        // unusable as a reference here (it overflows f64 past x~26.6), but
        // that is exactly where the standard asymptotic series becomes an
        // excellent reference in its own right:
        //   erfcx(x) ~ 1/(x*sqrt(pi)) * sum (-1)^n (2n-1)!! / (2x^2)^n
        // Truncating after the t^4 term leaves a relative error bounded by
        // the first dropped term, 945*t^5 with t = 1/(2x^2): at the x=20
        // left edge that is ~2.9e-12, ~5 orders of magnitude under f32's
        // own ~6e-8 resolution, and it only shrinks as x grows. So this row
        // measures this crate's error, not the reference's. Positive side
        // only -- the negative side is `2*e^(x^2) - erfcx(|x|)`, which has
        // genuinely overflowed to +inf for every x < -9.382 (pinned on
        // both sides of that boundary in edgecheck) and carries no
        // accuracy information out here.
        let erfcx_tail_ref = |v: F64xN| {
            let t = F64xN::splat(0.5) / (v * v);
            let p = F64xN::splat(1.0)
                - t * (F64xN::splat(1.0)
                    - t * (F64xN::splat(3.0)
                        - t * (F64xN::splat(15.0) - t * F64xN::splat(105.0))));
            p / (v * F64xN::splat(std::f64::consts::PI.sqrt()))
        };
        let erfcx_tail = |x: f32| x >= 20.0 && x.is_finite();
        let s = measure!(erfcx_tail, erfcx, erfcx_tail_ref);
        report("erfcx (x>=20)", &s, t0);
    }
    if run("erfinv") {
        // No sleef erfinv bucket, so verify via round-trip through erf_u10
        // instead of a direct reference: erf and erfinv are computed via
        // completely different mechanisms (erf's own poly/exp2 combine vs
        // erfinv's central/tail fit), so erf_u10(erfinv(x)) landing back
        // on x is real, independent evidence, not circular.
        let n_samples = 5_000_000u64;
        let mut max_dev = 0.0f64;
        let mut worst_x = 0.0f32;
        for _ in 0..n_samples {
            let x = f32::from_bits(rand::rng().random::<u32>());
            if !(x.abs() < 1.0) {
                continue;
            }
            let y = erfinv(x) as f64;
            let back = erf_u10(F64xN::splat(y)).to_array()[0];
            let dev = (back - x as f64).abs();
            if dev > max_dev {
                max_dev = dev;
                worst_x = x;
            }
        }
        println!(
            "{:24} max |erf(erfinv(x))-x| {:>10.6e}  worst x={:e} ({:>12} samples, {:>7.2}s elapsed)",
            "erfinv",
            max_dev,
            worst_x,
            n_samples,
            t0.elapsed().as_secs_f64(),
        );
        // probit/erfc_inv (backlog idea #139): same round-trip approach,
        // one level further out (through norm_cdf/erfc's own already-
        // verified accuracy) -- no sleef bucket for either exists.
        let mut max_dev_p = 0.0f64;
        let mut worst_p = 0.0f32;
        for _ in 0..n_samples {
            let p = f32::from_bits(rand::rng().random::<u32>());
            if !(p > 0.0 && p < 1.0) {
                continue;
            }
            let x = probit(p) as f64;
            let norm_cdf_ref = |v: F64xN| {
                let z = F64xN::splat(-std::f64::consts::FRAC_1_SQRT_2) * v;
                F64xN::splat(0.5) * erfc_u15(z)
            };
            let back = norm_cdf_ref(F64xN::splat(x)).to_array()[0];
            let dev = (back - p as f64).abs();
            if dev > max_dev_p {
                max_dev_p = dev;
                worst_p = p;
            }
        }
        println!(
            "{:24} max |norm_cdf(probit(p))-p| {:>10.6e}  worst p={:e} ({:>12} samples, {:>7.2}s elapsed)",
            "probit",
            max_dev_p,
            worst_p,
            n_samples,
            t0.elapsed().as_secs_f64(),
        );
        let mut max_dev_y = 0.0f64;
        let mut worst_y = 0.0f32;
        for _ in 0..n_samples {
            let y = f32::from_bits(rand::rng().random::<u32>());
            if !(y > 0.0 && y < 2.0) {
                continue;
            }
            let z = erfc_inv(y) as f64;
            let back = erfc_u15(F64xN::splat(z)).to_array()[0];
            let dev = (back - y as f64).abs();
            if dev > max_dev_y {
                max_dev_y = dev;
                worst_y = y;
            }
        }
        println!(
            "{:24} max |erfc(erfc_inv(y))-y| {:>10.6e}  worst y={:e} ({:>12} samples, {:>7.2}s elapsed)",
            "erfc_inv",
            max_dev_y,
            worst_y,
            n_samples,
            t0.elapsed().as_secs_f64(),
        );
    }
    if run("norm_cdf") {
        // Both compose already-full-range primitives (erfc/exp_checked),
        // so no domain restriction needed (backlog idea #67).
        let norm_cdf_ref = |v: F64xN| {
            F64xN::splat(0.5) * erfc_u15(-v * F64xN::splat(std::f64::consts::FRAC_1_SQRT_2))
        };
        let s = measure!(everywhere, norm_cdf, norm_cdf_ref);
        report("norm_cdf", &s, t0);
        let norm_pdf_ref =
            |v: F64xN| exp_u10(v * v * F64xN::splat(-0.5)) * F64xN::splat(0.3989422804014326779399460599);
        let s = measure!(everywhere, norm_pdf, norm_pdf_ref);
        report("norm_pdf", &s, t0);
    }
    if run("logit") {
        // Domain (0,1) (backlog idea #71) -- log1p_u10(-v) is
        // cancellation-safe for v near 1 the same way the real
        // implementation is.
        //
        // Deliberately keeps the difference form the real logit now only
        // uses *outside* its central band: in f64 the cancellation near
        // p=0.5 that made it a bug in f32 costs ~1e-16 absolute against
        // a result ~4e-4, i.e. ~1e-13 relative -- orders below an f32
        // ulp, so this stays a trustworthy reference for the very region
        // whose f32 version it is no longer good enough to compute.
        let unit_open = |x: f32| x > 0.0 && x < 1.0;
        let logit_ref = |v: F64xN| log_u35(v) - log1p_u10(-v);
        let s = measure!(unit_open, logit, logit_ref);
        report("logit", &s, t0);
    }
    if run("srgb") {
        // Standard sRGB channel range (backlog idea #146); powf_pos's
        // own x>=0 contract.
        let unit_range = |x: f32| (0.0..=1.0).contains(&x);
        let srgb_to_linear_ref = |v: F64xN| {
            let low = v * F64xN::splat(1.0 / 12.92);
            let high = pow_u10(
                (v + F64xN::splat(0.055)) * F64xN::splat(1.0 / 1.055),
                F64xN::splat(2.4),
            );
            v.simd_le(F64xN::splat(0.04045)).select(low, high)
        };
        let s = measure!(unit_range, srgb_to_linear, srgb_to_linear_ref);
        report("srgb_to_linear", &s, t0);
        let linear_to_srgb_ref = |v: F64xN| {
            let low = v * F64xN::splat(12.92);
            let high = F64xN::splat(1.055) * pow_u10(v, F64xN::splat(1.0 / 2.4))
                - F64xN::splat(0.055);
            v.simd_le(F64xN::splat(0.0031308)).select(low, high)
        };
        let s = measure!(unit_range, linear_to_srgb, linear_to_srgb_ref);
        report("linear_to_srgb", &s, t0);
    }

    // two-argument functions: fuzz-only (exhaustive over 2^64 pairs isn't
    // feasible), smaller sample count since each trial needs two RNG draws.
    const TWOARG_SAMPLES: u64 = 10_000_000;
    if run("atan2") {
        let s = fuzz2(TWOARG_SAMPLES, |_, _| true, atan2, atan2_u35);
        report("atan2", &s, t0);
        let s = fuzz2(TWOARG_SAMPLES, |_, _| true, atan2_latency, atan2_u35);
        report("atan2_latency", &s, t0);
        let s = fuzz2(TWOARG_SAMPLES, |_, _| true, |y: f32, x: f32| y.atan2(x), atan2_u35);
        report("std atan2", &s, t0);
        // atan2_unchecked's documented contract: x != 0.0, not both infinite.
        let atan2_domain = |_: f32, x: f32| x != 0.0;
        let s = fuzz2(TWOARG_SAMPLES, atan2_domain, atan2_unchecked, atan2_u35);
        report("atan2_unchecked (+)", &s, t0);
        // atan2_pos: single-positive-turn fold (backlog idea #143). Keyed
        // on y's sign bit, same as the shipped function -- an f64 atan2
        // never underflows to -0.0 on f32 inputs, so `r < 0.0` here would
        // score atan2_pos's deliberate `y == -0.0` fold as a full turn of
        // error at the one input where the two conventions differ.
        let atan2_pos_ref = |y: F64xN, x: F64xN| {
            let r = atan2_u35(y, x);
            y.is_sign_negative().select(r + F64xN::splat(std::f64::consts::TAU), r)
        };
        let s = fuzz2(TWOARG_SAMPLES, |_, _| true, atan2_pos, atan2_pos_ref);
        report("atan2_pos", &s, t0);
        // atan2d (backlog idea #123): plain composite.
        let atan2d_ref =
            |y: F64xN, x: F64xN| atan2_u35(y, x) * F64xN::splat(180.0 / std::f64::consts::PI);
        let s = fuzz2(TWOARG_SAMPLES, |_, _| true, atan2d, atan2d_ref);
        report("atan2d", &s, t0);
        // atan2pi (backlog idea #85): plain composite -- see its own doc
        // comment for why a rescaled-coefficient fold isn't attempted.
        let atan2pi_ref = |y: F64xN, x: F64xN| atan2_u35(y, x) * F64xN::splat(1.0 / std::f64::consts::PI);
        let s = fuzz2(TWOARG_SAMPLES, |_, _| true, atan2pi, atan2pi_ref);
        report("atan2pi", &s, t0);
    }
    if run("xlogy") {
        let s = fuzz2(TWOARG_SAMPLES, |_, _| true, xlogy, xlogy_ref);
        report("xlogy", &s, t0);
        let s = fuzz2(TWOARG_SAMPLES, |_, _| true, xlog1py, xlog1py_ref);
        report("xlog1py", &s, t0);
    }
    if run("compound") {
        // Domain x > -1 (backlog idea #72), log1p's own real-domain
        // limit; exp_checked handles any resulting exponent magnitude.
        let compound_domain = |x: f32, _: f32| x > -1.0;
        let compound_ref = |x: F64xN, n: F64xN| exp_u10(n * log1p_u10(x));
        let s = fuzz2(TWOARG_SAMPLES, compound_domain, compound, compound_ref);
        report("compound", &s, t0);
        let s = fuzz2(TWOARG_SAMPLES, compound_domain, compound_accurate, compound_ref);
        report("compound_accurate", &s, t0);
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
    if run("sqrt1pm1") {
        // First ulp sweep for this function -- it had edgecheck pins but no
        // accuracy coverage. Reference is `x / (sqrt(1+x) + 1)`, the
        // algebraically-equal form that does *not* repeat the cancellation
        // sqrt1pm1 itself exists to avoid: computing `sqrt(1+x) - 1`
        // directly in f64 would lose exactly the low bits being measured
        // for small |x| (same "the reference must not repeat the
        // cancellation trap" rule as log1pmx's own reference).
        let sqrt1pm1_domain = |x: f32| x >= -1.0 && x.is_finite();
        let s = measure!(sqrt1pm1_domain, sqrt1pm1, |v: F64xN| {
            v / ((F64xN::splat(1.0) + v).sqrt() + F64xN::splat(1.0))
        });
        report("sqrt1pm1", &s, t0);
    }
    if run("wrap_pi") {
        // First ulp sweep for wrap_pi too (edgecheck had the range and
        // sin/cos-preservation invariants, but no ulp measurement).
        //
        // Restricted to |x| <= 1e4 on purpose: wrap_pi rides
        // reduce_pi_checked's double-float reduction, which stays accurate
        // far past what a single f64 word can reference. Beyond this range
        // edgecheck's sin(wrap_pi(x)) == sin_checked(x) invariant is what
        // covers the reduction (out to 1e9).
        //
        // TAU is carried as two words (TAU_HI + TAU_LO) rather than one.
        // That is not pedantry: with a single-word TAU the reference itself
        // is off by ~25 ulp at the worst points, i.e. it would be the less
        // accurate of the two things being compared. Measured directly --
        // at x = 8953.539 a 1-word and 2-word f64 reference disagree by
        // 24.56 ulp.
        //
        // Even with the 2-word reference the reported max stays large, and
        // it is the crate's usual near-a-true-zero artifact rather than a
        // real defect: the worst inputs are the ones sitting almost exactly
        // on a multiple of 2*pi, where the answer is ~1e-7 formed by near
        // total cancellation of operands ~1e4, so one f32 ulp of the
        // *result* is a vanishingly small absolute quantity. Away from
        // those points wrap_pi measures <= 0.41 ulp. Same class as
        // compound's documented near-zero case.
        let wrap_domain = |x: f32| x.abs() <= 1e4;
        let s = measure!(wrap_domain, wrap_pi, |v: F64xN| {
            // 2*pi split into two f64 words: TAU_LO holds the part that the
            // nearest-f64 TAU_HI drops.
            let tau_hi = F64xN::splat(6.283185307179586);
            let tau_lo = F64xN::splat(2.4492935982947064e-16);
            let pi = F64xN::splat(std::f64::consts::PI);
            let tau = tau_hi + tau_lo;
            let q = (v / tau_hi).round();
            let r = (v - q * tau_hi) - q * tau_lo;
            // fold into (-pi, pi]: r lands in [-tau/2, tau/2] but the
            // half-open convention needs r > pi pulled down and
            // r <= -pi pushed up.
            let r = r.simd_gt(pi).select(r - tau, r);
            r.simd_le(-pi).select(r + tau, r)
        });
        report("wrap_pi", &s, t0);
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
        // normalize2 (backlog idea #136): check the result actually has
        // unit magnitude, not a ulp comparison against a single reference.
        let mut max_dev = 0.0f64;
        for _ in 0..1_000_000u64 {
            let x = f32::from_bits(rand::rng().random::<u32>());
            let y = f32::from_bits(rand::rng().random::<u32>());
            if !(hypot_domain(x, y)) || (x == 0.0 && y == 0.0) {
                continue;
            }
            let (nx, ny) = normalize2(x, y);
            let mag = hypot(nx, ny) as f64;
            max_dev = max_dev.max((mag - 1.0).abs());
        }
        println!("{:24} max |magnitude-1| {:>10.6e}", "normalize2", max_dev);
    }
    if run("hypot3") {
        // Same naive-fma-chain overflow/underflow tradeoff as hypot (see
        // its own domain comment) -- no sleef hypot3, so a direct f64
        // sqrt(x^2+y^2+z^2) is ground truth (plenty of headroom rounding
        // down to f32).
        let ok = |v: f32| v == 0.0 || (v.abs() > 1e-15 && v.abs() < 1e18);
        let n_samples = 10_000_000u64;
        for (label, f) in [
            ("hypot3", hypot3 as fn(f32, f32, f32) -> f32),
            ("rnorm3", rnorm3 as fn(f32, f32, f32) -> f32),
        ] {
            let mut sum = 0u64;
            let mut max = 0u64;
            let mut worst = (0.0f32, 0.0f32, 0.0f32);
            for _ in 0..n_samples {
                let x = f32::from_bits(rand::rng().random::<u32>());
                let y = f32::from_bits(rand::rng().random::<u32>());
                let z = f32::from_bits(rand::rng().random::<u32>());
                if !(ok(x) && ok(y) && ok(z)) {
                    continue;
                }
                let got = f(x, y, z);
                let norm = ((x as f64).powi(2) + (y as f64).powi(2) + (z as f64).powi(2)).sqrt();
                let want = (if label == "hypot3" { norm } else { 1.0 / norm }) as f32;
                let d = ulp_diff(got, want);
                sum += d;
                if d > max {
                    max = d;
                    worst = (x, y, z);
                }
            }
            println!(
                "{:24} avg ulp {:>10.4}  max ulp {:>10}  worst x={:e},y={:e},z={:e} ({:>12} samples, {:>7.2}s elapsed)",
                label,
                sum as f64 / n_samples as f64,
                max,
                worst.0,
                worst.1,
                worst.2,
                n_samples,
                t0.elapsed().as_secs_f64(),
            );
        }
        // normalize3 (backlog idea #136): same magnitude check as
        // normalize2/normalize4.
        let ok = |v: f32| v == 0.0 || (v.abs() > 1e-15 && v.abs() < 1e18);
        let mut max_dev = 0.0f64;
        for _ in 0..1_000_000u64 {
            let x = f32::from_bits(rand::rng().random::<u32>());
            let y = f32::from_bits(rand::rng().random::<u32>());
            let z = f32::from_bits(rand::rng().random::<u32>());
            if !(ok(x) && ok(y) && ok(z)) || (x == 0.0 && y == 0.0 && z == 0.0) {
                continue;
            }
            let (nx, ny, nz) = normalize3(x, y, z);
            let mag = hypot3(nx, ny, nz) as f64;
            max_dev = max_dev.max((mag - 1.0).abs());
        }
        println!("{:24} max |magnitude-1| {:>10.6e}", "normalize3", max_dev);
    }
    if run("hypot4") {
        // Same tradeoff and ground-truth approach as hypot3 above, one
        // argument wider (backlog idea #134).
        let ok = |v: f32| v == 0.0 || (v.abs() > 1e-15 && v.abs() < 1e18);
        let n_samples = 10_000_000u64;
        for (label, f) in [
            ("hypot4", hypot4 as fn(f32, f32, f32, f32) -> f32),
            ("rnorm4", rnorm4 as fn(f32, f32, f32, f32) -> f32),
        ] {
            let mut sum = 0u64;
            let mut max = 0u64;
            let mut worst = (0.0f32, 0.0f32, 0.0f32, 0.0f32);
            for _ in 0..n_samples {
                let w = f32::from_bits(rand::rng().random::<u32>());
                let x = f32::from_bits(rand::rng().random::<u32>());
                let y = f32::from_bits(rand::rng().random::<u32>());
                let z = f32::from_bits(rand::rng().random::<u32>());
                if !(ok(w) && ok(x) && ok(y) && ok(z)) {
                    continue;
                }
                let got = f(w, x, y, z);
                let norm = ((w as f64).powi(2)
                    + (x as f64).powi(2)
                    + (y as f64).powi(2)
                    + (z as f64).powi(2))
                .sqrt();
                let want = (if label == "hypot4" { norm } else { 1.0 / norm }) as f32;
                let d = ulp_diff(got, want);
                sum += d;
                if d > max {
                    max = d;
                    worst = (w, x, y, z);
                }
            }
            println!(
                "{:24} avg ulp {:>10.4}  max ulp {:>10}  worst w={:e},x={:e},y={:e},z={:e} ({:>12} samples, {:>7.2}s elapsed)",
                label,
                sum as f64 / n_samples as f64,
                max,
                worst.0,
                worst.1,
                worst.2,
                worst.3,
                n_samples,
                t0.elapsed().as_secs_f64(),
            );
        }
        // normalize4: check the result actually has unit magnitude
        // (hypot4 of the normalized components ~= 1), not a ulp
        // comparison against a single reference value.
        let mut max_dev = 0.0f64;
        for _ in 0..1_000_000u64 {
            let w = f32::from_bits(rand::rng().random::<u32>());
            let x = f32::from_bits(rand::rng().random::<u32>());
            let y = f32::from_bits(rand::rng().random::<u32>());
            let z = f32::from_bits(rand::rng().random::<u32>());
            if !(ok(w) && ok(x) && ok(y) && ok(z)) || (w == 0.0 && x == 0.0 && y == 0.0 && z == 0.0)
            {
                continue;
            }
            let (nw, nx, ny, nz) = normalize4(w, x, y, z);
            let mag = hypot4(nw, nx, ny, nz) as f64;
            let dev = (mag - 1.0).abs();
            if dev > max_dev {
                max_dev = dev;
            }
        }
        println!("{:24} max |magnitude-1| {:>10.6e}", "normalize4", max_dev);
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
    if run("complex") {
        // cexp/clog (backlog idea #186) return a (f32,f32) pair, not the
        // single f32 fuzz2/fuzz4 expect, so this is a custom loop, same
        // shape as erfinv's/rootn's own. cabs/carg are exact aliases for
        // hypot_checked/atan2 (bit-identical by construction, verified
        // by their own extensive fuzzing already), so no separate check
        // for those two -- only the genuinely new compositions get
        // fuzzed.
        let n_samples = 5_000_000u64;
        let mag_domain = |x: f32| x.abs() < 80.0; // exp(re) overflows past here
        // cexp's own im has no such bound mathematically (cos/sin are
        // total), but sin/cos's *own* accuracy already degrades past
        // their documented large-argument range -- not a new cexp gap,
        // so this just avoids re-measuring an already-known limit here.
        let im_domain = |x: f32| x.abs() < 1e4;
        let mut max_ulp_re = 0u64;
        let mut max_ulp_im = 0u64;
        let mut worst = (0.0f32, 0.0f32);
        for _ in 0..n_samples {
            let re = f32::from_bits(rand::rng().random::<u32>());
            let im = f32::from_bits(rand::rng().random::<u32>());
            if !mag_domain(re) || !im_domain(im) {
                continue;
            }
            let (got_re, got_im) = cexp(re, im);
            let want_re = (re as f64).exp() * (im as f64).cos();
            let want_im = (re as f64).exp() * (im as f64).sin();
            let d_re = ulp_diff(got_re, want_re as f32);
            let d_im = ulp_diff(got_im, want_im as f32);
            if d_re > max_ulp_re {
                max_ulp_re = d_re;
                worst = (re, im);
            }
            if d_im > max_ulp_im {
                max_ulp_im = d_im;
            }
        }
        println!(
            "{:24} max ulp re {:>10} im {:>10}  worst (re,im)={:e},{:e} ({:>12} samples, {:>7.2}s elapsed)",
            "cexp",
            max_ulp_re,
            max_ulp_im,
            worst.0,
            worst.1,
            n_samples,
            t0.elapsed().as_secs_f64(),
        );
        let mut max_ulp_re = 0u64;
        let mut max_ulp_im = 0u64;
        let mut worst = (0.0f32, 0.0f32);
        for _ in 0..n_samples {
            let re = f32::from_bits(rand::rng().random::<u32>());
            let im = f32::from_bits(rand::rng().random::<u32>());
            if !re.is_finite() || !im.is_finite() || (re == 0.0 && im == 0.0) {
                continue;
            }
            let (got_re, got_im) = clog(re, im);
            let want_re = (re as f64).hypot(im as f64).ln();
            let want_im = (im as f64).atan2(re as f64);
            let d_re = ulp_diff(got_re, want_re as f32);
            let d_im = ulp_diff(got_im, want_im as f32);
            if d_re > max_ulp_re {
                max_ulp_re = d_re;
                worst = (re, im);
            }
            if d_im > max_ulp_im {
                max_ulp_im = d_im;
            }
        }
        println!(
            "{:24} max ulp re {:>10} im {:>10}  worst (re,im)={:e},{:e} ({:>12} samples, {:>7.2}s elapsed)",
            "clog",
            max_ulp_re,
            max_ulp_im,
            worst.0,
            worst.1,
            n_samples,
            t0.elapsed().as_secs_f64(),
        );
        // Round-trip: clog(cexp(re,im)) should recover re exactly-ish and
        // im modulo 2*pi, wrapped to carg's own principal range -- check
        // against a bounded re (cexp's own domain above) and any finite im.
        let mut max_dev_re = 0.0f64;
        let mut max_dev_im = 0.0f64;
        for _ in 0..n_samples {
            let re = f32::from_bits(rand::rng().random::<u32>());
            let im = f32::from_bits(rand::rng().random::<u32>());
            if !mag_domain(re) || !im_domain(im) {
                continue;
            }
            let (a, b) = cexp(re, im);
            let (back_re, back_im) = clog(a, b);
            let dev_re = (back_re as f64 - re as f64).abs();
            // wrap im to (-pi,pi] the same way carg's atan2 does before comparing
            let two_pi = std::f64::consts::TAU;
            let im_wrapped = im as f64 - two_pi * ((im as f64 + std::f64::consts::PI) / two_pi).floor();
            let dev_im = (back_im as f64 - im_wrapped).abs();
            max_dev_re = max_dev_re.max(dev_re);
            max_dev_im = max_dev_im.max(dev_im);
        }
        println!(
            "{:24} max |re-back| {:>10.4e} max |im-back| {:>10.4e} ({:>12} samples, {:>7.2}s elapsed)",
            "clog(cexp(.))",
            max_dev_re,
            max_dev_im,
            n_samples,
            t0.elapsed().as_secs_f64(),
        );
    }
    if run("rootn") {
        // x^(1/n) (backlog idea #75): reference matches the documented
        // C23 domain-error/sign rules (see rootn's own doc comment),
        // not a naive x.powf(1.0/n) -- that would mishandle negative x
        // with odd n the exact same way a naive implementation would.
        let rootn_ref = |x: f64, n: i32| -> f64 {
            if n == 0 {
                return f64::NAN;
            }
            if x == 0.0 {
                let n_odd = n % 2 != 0;
                let mag = if n > 0 { 0.0 } else { f64::INFINITY };
                return if n_odd { mag.copysign(x) } else { mag };
            }
            if x < 0.0 {
                return if n % 2 == 0 { f64::NAN } else { -((-x).powf(1.0 / n as f64)) };
            }
            x.powf(1.0 / n as f64)
        };
        let ns: [i32; 24] = [
            1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 7, -7, 8, -8, 10, -10, 31, -31, 1000, -1000,
            1 << 24, -(1 << 24), i32::MAX, i32::MIN,
        ];
        let n_samples = 2_000_000u64;
        for &n in &ns {
            let mut sum = 0u64;
            let mut max = 0u64;
            let mut worst = 0.0f32;
            for _ in 0..n_samples {
                let x = f32::from_bits(rand::rng().random::<u32>());
                if !x.is_finite() {
                    continue;
                }
                let want_f64 = rootn_ref(x as f64, n);
                if !want_f64.is_finite() {
                    continue;
                }
                let got = rootn(x, n);
                let want = want_f64 as f32;
                let d = ulp_diff(got, want);
                sum += d;
                if d > max {
                    max = d;
                    worst = x;
                }
            }
            println!(
                "{:24} avg ulp {:>10.4}  max ulp {:>10}  worst x={:e} ({:>12} samples, {:>7.2}s elapsed)",
                format!("rootn(x,{n})"),
                sum as f64 / n_samples as f64,
                max,
                worst,
                n_samples,
                t0.elapsed().as_secs_f64(),
            );
        }
    }
    if run("ldexp") {
        // ldexp/frexp (backlog idea #86): exact bit manipulations, not
        // approximations -- expected max ulp is always exactly 0. Two
        // real bugs were found and fixed during development via a much
        // larger dedicated scratch sweep (see ldexp's own doc comment);
        // this is a permanent, smaller-scale regression guard, not a
        // substitute for that sweep.
        let n_samples = 20_000_000u64;
        let mut sum = 0u64;
        let mut max = 0u64;
        let mut worst = (0.0f32, 0i32);
        for _ in 0..n_samples {
            let x = f32::from_bits(rand::rng().random::<u32>());
            let n: i32 = rand::rng().random_range(-2000..=2000);
            if !x.is_finite() {
                continue;
            }
            let got = ldexp(x, n);
            let want = (x as f64 * 2f64.powi(n.clamp(-1100, 1100))) as f32;
            let d = ulp_diff(got, want);
            sum += d;
            if d > max {
                max = d;
                worst = (x, n);
            }
        }
        println!(
            "{:24} avg ulp {:>10.4}  max ulp {:>10}  worst x={:e},n={} ({:>12} samples, {:>7.2}s elapsed)",
            "ldexp",
            sum as f64 / n_samples as f64,
            max,
            worst.0,
            worst.1,
            n_samples,
            t0.elapsed().as_secs_f64(),
        );
        // frexp: exact reconstruction + mantissa-range check, not a ulp
        // measurement (it's a decomposition, not an approximation) --
        // reported as a bad-count, not avg/max ulp.
        let mut fbad = 0u64;
        for _ in 0..n_samples {
            let x = f32::from_bits(rand::rng().random::<u32>());
            let (m, e) = frexp(x);
            let ok = if x == 0.0 {
                m.to_bits() == x.to_bits() && e == 0
            } else if !x.is_finite() {
                (m.is_nan() && x.is_nan()) || m.to_bits() == x.to_bits()
            } else {
                let recon = (m as f64 * 2f64.powi(e)) as f32;
                recon.to_bits() == x.to_bits() && m.abs() >= 0.5 && m.abs() < 1.0
            };
            if !ok {
                fbad += 1;
            }
        }
        println!(
            "{:24} bad reconstructions {:>10} / {:<12} ({:>7.2}s elapsed)",
            "frexp",
            fbad,
            n_samples,
            t0.elapsed().as_secs_f64(),
        );
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
    if run("dawson") {
        // No sleef bucket for Dawson's function, so the reference here is
        // a real independent computation, not a round-trip: Simpson's-rule
        // quadrature of the defining integral D(x) = x * integral_0^1
        // exp(-x^2*(1-s^2)) ds for |x| <= 5 (N=800; measured against
        // scipy.special.dawsn, its own relative error is <1e-11 out to
        // x=2, 1.5e-8 (0.13 f32 ulp) at x=4 and 8.8e-8 (0.74 ulp) at the
        // x=5 handover -- so it is a real reference below x~4, and adds
        // up to ~1 ulp of its own noise in [4,5]),
        // the literal double-factorial asymptotic series for |x| > 5 (15
        // terms, converges to near f64 precision well before the series'
        // own eventual divergence past its optimal truncation point).
        let dawson_ref = |x: f64| -> f64 {
            let ax = x.abs();
            let mag = if ax <= 5.0 {
                const N: usize = 800;
                let h = 1.0 / N as f64;
                let mut sum = (-ax * ax).exp() + 1.0;
                for i in 1..N {
                    let s = i as f64 * h;
                    let f = (-ax * ax * (1.0 - s * s)).exp();
                    sum += if i % 2 == 1 { 4.0 * f } else { 2.0 * f };
                }
                ax * (h / 3.0) * sum
            } else {
                let v = 1.0 / (ax * ax);
                let mut term = 1.0;
                let mut acc = 1.0;
                for k in 1..=15 {
                    term *= (2.0 * k as f64 - 1.0) * v * 0.5;
                    acc += term;
                }
                acc / (2.0 * ax)
            };
            mag.copysign(x)
        };
        let n_samples = 2_000_000u64;
        let mut sum = 0u64;
        let mut max = 0u64;
        let mut worst = 0.0f32;
        for _ in 0..n_samples {
            let x = f32::from_bits(rand::rng().random::<u32>());
            if !x.is_finite() {
                continue;
            }
            let want = dawson_ref(x as f64) as f32;
            let got = dawson(x);
            let d = ulp_diff(got, want);
            sum += d;
            if d > max {
                max = d;
                worst = x;
            }
        }
        println!(
            "{:24} avg ulp {:>10.4}  max ulp {:>10}  worst x={:e} ({:>12} samples, {:>7.2}s elapsed)",
            "dawson",
            sum as f64 / n_samples as f64,
            max,
            worst,
            n_samples,
            t0.elapsed().as_secs_f64(),
        );
    }

    if run("identities") {
        // Cross-function algebraic identity fuzz (backlog idea #167): a
        // cheap bug detector independent of any single function's own
        // f64 reference. Domain-restricted per identity to where the
        // check itself stays well-conditioned -- e.g. cosh(x)^2-sinh(x)^2
        // and exp(x)-1 vs expm1(x) both look like real failures if
        // checked by *absolute* difference at large x (both sides are
        // individually huge, so the identity's true near-zero residual
        // is swamped by rounding in forming the huge intermediates
        // themselves -- a real "does this test even make sense" trap,
        // confirmed by checking relative error instead: exactly 0 at
        // every magnitude tried). Every identity here already passed a
        // 30M-sample sweep with no real finding -- kept as a standing
        // regression gate, not because a bug was found.
        let n_samples = 20_000_000u64;
        let identity = |name: &str, domain: &dyn Fn(f32) -> bool, resid: &dyn Fn(f32) -> f64, tol: f64| {
            let mut max_dev = 0.0f64;
            let mut worst = 0.0f32;
            let mut count = 0u64;
            for _ in 0..n_samples {
                let x = f32::from_bits(rand::rng().random::<u32>());
                if !domain(x) {
                    continue;
                }
                let d = resid(x).abs();
                count += 1;
                if d > max_dev {
                    max_dev = d;
                    worst = x;
                }
            }
            let flag = if max_dev > tol { "FLAG" } else { "ok  " };
            println!(
                "{flag} {name:34} max |residual| {:>12.4e} worst x={:e} (n={count}, {:.2}s elapsed)",
                max_dev,
                worst,
                t0.elapsed().as_secs_f64(),
            );
        };
        identity(
            "sin_checked^2+cos_checked^2=1",
            &|x| x.is_finite() && x.abs() < 8.85e14, // sin_checked's own documented exact-reduction limit
            &|x| {
                let s = sin_checked(x) as f64;
                let c = cos_checked(x) as f64;
                s * s + c * c - 1.0
            },
            1e-4,
        );
        identity(
            "tanh(x)=sinh(x)/cosh(x)",
            &|x| x.is_finite() && x.abs() < 80.0,
            &|x| (tanh(x) as f64) - (sinh(x) as f64) / (cosh(x) as f64),
            1e-4,
        );
        identity("exp(ln(x))=x", &|x| x > 0.0 && x.is_finite(), &|x| (exp(ln(x)) as f64) / (x as f64) - 1.0, 1e-4);
        identity(
            "ln(exp(x))=x",
            &|x| x.is_finite() && x.abs() < 80.0,
            &|x| (ln(exp(x)) as f64) - (x as f64),
            1e-2,
        );
        identity(
            "sigmoid(x)+sigmoid(-x)=1",
            &|x| x.is_finite(),
            &|x| (sigmoid(x) as f64) + (sigmoid(-x) as f64) - 1.0,
            1e-4,
        );
        identity("erf(x)+erfc(x)=1", &|x| x.is_finite(), &|x| (erf(x) as f64) + (erfc(x) as f64) - 1.0, 1e-4);
        identity("erf(-x)=-erf(x)", &|x| x.is_finite(), &|x| (erf(-x) as f64) + (erf(x) as f64), 1e-6);
        identity(
            "atan2(sin(x),cos(x))=x [|x|<pi]",
            &|x| x.abs() < 3.0,
            &|x| (atan2(sin(x), cos(x)) as f64) - (x as f64),
            1e-2,
        );
        identity(
            "log1p(x)=ln(1+x)",
            &|x| x > -1.0 && x.is_finite() && x.abs() < 1e6,
            &|x| (log1p(x) as f64) - (ln(1.0 + x) as f64),
            1e-2,
        );
    }

    println!("total: {:.2}s", t0.elapsed().as_secs_f64());
}
