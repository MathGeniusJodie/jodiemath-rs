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
/// Chain seed. Only its mantissa survives the first `Band::mix`, so the same
/// literal works for every band.
const SEED: f32 = 1.234;
/// How far `check_in_domain` walks the chain before declaring the band sound.
const CHECK_ITERS: usize = 4096;

/// Untimed pre-pass: walk the same chain `bench_latency` will time and shout
/// if the function ever goes non-finite. A NaN means the band is outside the
/// function's domain, which silently turns the row into a measurement of the
/// *domain check* rather than the function -- for std that is an early `ret`
/// a few ns long, so the comparison inverts. See `Band`'s own doc comment for
/// the rows this actually happened to; `log1pmx` is the one this check caught
/// on its own, and the one worth keeping the check around for: its magnitude
/// was never the problem, its own negative output walked the chain past the
/// `x > -1` edge, which is not something you spot by reading a domain off a
/// signature.
fn check_in_domain(name: &str, band: Band, f: impl Fn(f32) -> f32) {
    let mut x = band.mix(SEED);
    for _ in 0..CHECK_ITERS {
        let y = f(x);
        if !y.is_finite() {
            eprintln!(
                "!! {name}: f({x:e}) = {y} -- band is outside this function's domain, \
                 the numbers below do not measure it. Pick a different Band."
            );
            return;
        }
        x = band.mix(y);
    }
}

fn bench_latency(name: &str, band: Band, f: impl Fn(f32) -> f32) {
    let mut best = f64::INFINITY;
    for _ in 0..REPS {
        // dependency chain; `mix` keeps the value inside the band (and so
        // inside the function's domain) without branching
        let mut x = band.mix(SEED);
        let start = Instant::now();
        for _ in 0..LAT_ITERS {
            x = band.mix(f(x));
        }
        black_box(x);
        let ns = start.elapsed().as_nanos() as f64 / LAT_ITERS as f64;
        best = best.min(ns);
    }
    println!("{:22} latency    {:6.2} ns/op (min of {REPS}, incl chain overhead)", name, best);
}

// idea #78: N independent serial dependency chains interleaved in the
// same loop, instead of bench_latency's single chain. A single chain's
// own measured latency can be masked or exaggerated by how much
// cross-call ILP the CPU's out-of-order engine can extract when there's
// only one outstanding dependency chain to fill scheduling gaps with --
// 4 independent chains approximates a caller that has some real
// concurrent work available, closer to how throughput callers (many
// independent elements) sit between this and the single-chain extreme.
const N_STREAMS: usize = 4;
fn bench_latency_n(name: &str, band: Band, f: impl Fn(f32) -> f32) {
    let mut best = f64::INFINITY;
    for _ in 0..REPS {
        let mut xs = [1.234_f32, 1.876, 1.456, 1.987].map(|s| band.mix(s));
        let start = Instant::now();
        for _ in 0..(LAT_ITERS / N_STREAMS as u64) {
            for x in xs.iter_mut() {
                *x = band.mix(f(*x));
            }
        }
        black_box(xs);
        let ns = start.elapsed().as_nanos() as f64 / (LAT_ITERS / N_STREAMS as u64 * N_STREAMS as u64) as f64;
        best = best.min(ns);
    }
    println!(
        "{:22} latency(x{N_STREAMS}) {:6.2} ns/op (min of {REPS}, incl chain overhead)",
        name, best
    );
}

fn bench_throughput(name: &str, band: Band, f: impl Fn(f32) -> f32) {
    // fixed-size arrays: no bounds checks, so the loop can auto-vectorize.
    // Same band as the latency chain, for the same domain reason -- this
    // array was hardcoded to [2,4) too.
    let (lo, hi) = band.range();
    let mut input = [0f32; TP_ARR];
    for i in 0..TP_ARR {
        input[i] = lo + (i as f32) * ((hi - lo) / TP_ARR as f32);
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
    // Real wall-clock benchmarking, same contention concern as accuracy.rs's
    // own nice_self -- lower priority before any real work runs.
    nice_self();
    let args: Vec<String> = std::env::args().collect();
    if args.get(1).map(|s| s.as_str()) == Some("latencyn") {
        // idea #78 spot-check: single-chain vs 4-independent-chain latency
        // for a spread of functions with different codegen shapes (exp2's
        // short balanced-Estrin chain, sinh's branch-selected two-branch
        // combine, acos's plain Horner).
        bench_latency("exp2", Band::Two, exp2);
        bench_latency_n("exp2", Band::Two, exp2);
        bench_latency("sinh", Band::Two, sinh);
        bench_latency_n("sinh", Band::Two, sinh);
        bench_latency("acos", Band::Half, acos);
        bench_latency_n("acos", Band::Half, acos);
        return;
    }
    let filter = args.get(1).map(|s| s.as_str()).unwrap_or("");
    let run = |n: &str| filter.is_empty() || n.contains(filter);

    // Third argument is the input `Band`. It defaults to `Band::Two` (|x| in
    // [2,4)), which is only valid for functions defined there -- anything
    // restricted to [-1,1], (0,1) or [0,1] must name its own band, or the row
    // measures a domain check instead of the function. `check_in_domain`
    // enforces that at runtime.
    macro_rules! bench {
        ($name:expr, $f:expr) => {
            bench!($name, $f, Band::Two)
        };
        ($name:expr, $f:expr, $band:expr) => {
            if run($name) {
                check_in_domain($name, $band, $f);
                bench_latency($name, $band, $f);
                bench_throughput($name, $band, $f);
            }
        };
    }

    bench!("nop", |x: f32| x);
    bench!("cbrt", cbrt);
    bench!("std cbrt", |x: f32| x.cbrt());
    bench!("cbrt_unchecked", cbrt_unchecked);
    bench!("cbrt_accurate", cbrt_accurate);
    bench!("cbrt_accurate_unchecked", cbrt_accurate_unchecked);
    bench!("cbrt_fast", cbrt_fast);
    bench!("rcbrt", rcbrt);
    bench!("pow_3_2", pow_3_2);
    bench!("pow_2_3", pow_2_3);
    bench!("smoothstep", move |x: f32| smoothstep(0.0, 1.0, x), Band::Half);
    bench!("smootherstep", move |x: f32| smootherstep(0.0, 1.0, x), Band::Half);
    bench!("exp2", exp2);
    bench!("exp2_kf", |f: f32| exp2_kf(3.0, f));
    bench!("exp2_checked", exp2_checked);
    bench!("std exp2", |x: f32| x.exp2());
    bench!("exp10", exp10);
    bench!("exp10_checked", exp10_checked);
    bench!("log_2", log_2);
    bench!("log_2_unchecked", log_2_unchecked);
    bench!("std log2", |x: f32| x.log2());
    bench!("sin", sin);
    bench!("sin_checked", sin_checked);
    bench!("std sin", |x: f32| x.sin());
    bench!("cos", cos);
    bench!("cos_checked", cos_checked);
    bench!("sin_wide", sin_wide);
    bench!("cos_wide", cos_wide);
    bench!("std cos", |x: f32| x.cos());
    // reduce_pi_checked/reduce_pi_half_checked (idea #88): (f32,f32),
    // same tuple-adapter reasoning as cexp/clog above.
    bench!("reduce_pi_checked", |x: f32| {
        let (r, s) = reduce_pi_checked(x);
        r + s
    });
    bench!("reduce_pi_half_checked", |x: f32| {
        let (r, s) = reduce_pi_half_checked(x);
        r + s
    });
    bench!("wrap_pi", wrap_pi);
    bench!("sin_prereduced", sin_prereduced);
    bench!("cos_prereduced", cos_prereduced);
    bench!("sinpi", sinpi);
    bench!("sinpi_unchecked", sinpi_unchecked);
    bench!("cospi", cospi);
    bench!("tanpi", tanpi);
    bench!("sin2pi", sin2pi);
    bench!("cos2pi", cos2pi);
    bench!("tan2pi", tan2pi);
    bench!("sinc", sinc);
    bench!("sinc_unnormalized", sinc_unnormalized);
    bench!("sind", sind);
    bench!("sind_unchecked", sind_unchecked);
    bench!("cosd", cosd);
    bench!("cosd_unchecked", cosd_unchecked);
    bench!("tand", tand);
    bench!("tand_unchecked", tand_unchecked);

    bench!("ln", ln);
    bench!("ln_unchecked", ln_unchecked);
    bench!("std ln", |x: f32| x.ln());
    bench!("log10", log10);
    bench!("log10_unchecked", log10_unchecked);
    bench!("std log10", |x: f32| x.log10());
    bench!("log1p", log1p);
    // domain is x > -1, and log1pmx(x) is negative for every x, so the chain
    // walks itself out of domain under Band::Two (caught by check_in_domain).
    bench!("log1pmx", log1pmx, Band::Half);
    bench!("std log1p", |x: f32| x.ln_1p());
    bench!("log2p1", log2p1);
    bench!("exp", exp);
    bench!("exp_scaled", |x: f32| exp_scaled(x, 3));
    bench!("std exp", |x: f32| x.exp());
    bench!("exp_narrow", exp_narrow);
    bench!("exp_checked", exp_checked);
    bench!("expm1", expm1);
    bench!("expm1_narrow", expm1_narrow);
    bench!("std expm1", |x: f32| x.exp_m1());
    bench!("exp_m1_over_x", exp_m1_over_x);
    bench!("exp_m1_over_x_narrow", exp_m1_over_x_narrow);
    bench!("exp2m1", exp2m1);
    bench!("sinh", sinh);
    bench!("std sinh", |x: f32| x.sinh());
    bench!("sinh_narrow", sinh_narrow);
    bench!("cosh", cosh);
    bench!("std cosh", |x: f32| x.cosh());
    bench!("cosh_narrow", cosh_narrow);
    bench!("sinh_throughput", sinh_throughput);
    bench!("cosh_throughput", cosh_throughput);
    bench!("sinh_checked", sinh_checked);
    bench!("cosh_checked", cosh_checked);
    bench!("coshm1", coshm1);
    bench!("tanh", tanh);
    bench!("std tanh", |x: f32| x.tanh());
    bench!("tanh_grad", tanh_grad);
    bench!("sigmoid", sigmoid);
    bench!("sigmoid_fast", sigmoid_fast);
    bench!("sigmoid_grad", sigmoid_grad);
    bench!("softplus", softplus);
    bench!("logsigmoid", logsigmoid);
    bench!("logaddexp", |x: f32| logaddexp(x, 0.0));
    bench!("logaddexp_accurate", |x: f32| logaddexp_accurate(x, 0.0));
    bench!("asinh", asinh);
    bench!("std asinh", |x: f32| x.asinh());
    bench!("acosh", acosh);
    bench!("std acosh", |x: f32| x.acosh());
    bench!("atanh", atanh, Band::Half);
    bench!("std atanh", |x: f32| x.atanh(), Band::Half);
    bench!("asin", asin, Band::Half);
    bench!("asind", asind, Band::Half);
    bench!("asinpi", asinpi, Band::Half);
    bench!("std asin", |x: f32| x.asin(), Band::Half);
    bench!("acos", acos, Band::Half);
    bench!("std acos", |x: f32| x.acos(), Band::Half);
    bench!("acosd", acosd, Band::Half);
    bench!("acospi", acospi, Band::Half);
    bench!("atan", atan);
    bench!("std atan", |x: f32| x.atan());
    bench!("atan_latency", atan_latency);
    bench!("atan_bounded", atan_bounded);
    bench!("atand", atand);
    bench!("atanpi", atanpi);
    // black_box'd 2nd arg (not a literal 1.0): a compile-time-constant 2nd
    // arg lets LLVM fold away atan2's own special-case branches entirely,
    // silently hiding their real cost -- this matters here specifically
    // because it would otherwise make atan2_unchecked look like a wash
    // instead of the real win it is. See readme.md's own todo note on
    // this exact fixed-argument limitation.
    let atan2_x2 = std::hint::black_box(1.0);
    bench!("atan2", move |x: f32| atan2(x, atan2_x2));
    bench!("atan2_latency", move |x: f32| atan2_latency(x, atan2_x2));
    bench!("std atan2", move |x: f32| x.atan2(atan2_x2));
    bench!("atan2_unchecked", move |x: f32| atan2_unchecked(x, atan2_x2));
    bench!("atan2_pos", move |x: f32| atan2_pos(x, atan2_x2));
    bench!("atan2d", move |x: f32| atan2d(x, atan2_x2));
    bench!("atan2pi", move |x: f32| atan2pi(x, atan2_x2));
    bench!("tan", tan);
    bench!("tan_checked", tan_checked);
    bench!("std tan", |x: f32| x.tan());
    bench!("erf", erf);
    bench!("erfc", erfc);
    bench!("norm_cdf", norm_cdf);
    bench!("norm_pdf", norm_pdf);
    bench!("logit", logit, Band::Half);
    bench!("compound", move |x: f32| compound(x, 5.0));
    bench!("xlogy", move |x: f32| xlogy(x, 2.0));
    bench!("xlog1py", move |x: f32| xlog1py(x, 1.0));
    bench!("ldexp", move |x: f32| ldexp(x, 5));
    bench!("frexp", |x: f32| {
        let (m, e) = frexp(x);
        m + e as f32
    });
    bench!("erfcx", erfcx);
    bench!("erfinv", erfinv, Band::Half);
    bench!("erfc_inv", erfc_inv, Band::Half);
    bench!("probit", probit, Band::Half);
    bench!("dawson", dawson);
    // black_box'd 2nd arg, same reasoning as atan2 above.
    let hypot_y = std::hint::black_box(1.0);
    bench!("hypot", move |x: f32| hypot(x, hypot_y));
    bench!("std hypot", move |x: f32| x.hypot(hypot_y));
    bench!("hypot_unchecked", move |x: f32| hypot_unchecked(x, hypot_y));
    bench!("hypot_checked", move |x: f32| hypot_checked(x, hypot_y));
    bench!("rhypot", move |x: f32| rhypot(x, hypot_y));
    // Complex pack (idea #186): cabs/carg fix the 2nd arg like hypot/
    // atan2 above. cexp/clog return (f32,f32) -- bench!/mca need a
    // single f32, so this sums the pair as a cheap adapter (dominated
    // by the same ops either way; not used for mca given clog's own
    // branching, see mca_target.rs's own note).
    bench!("cabs", move |x: f32| cabs(x, hypot_y));
    bench!("carg", move |x: f32| carg(x, hypot_y));
    bench!("cexp", move |x: f32| {
        let (a, b) = cexp(x, hypot_y);
        a + b
    });
    bench!("clog", move |x: f32| {
        let (a, b) = clog(x, hypot_y);
        a + b
    });
    bench!("normalize2", move |x: f32| {
        let (a, b) = normalize2(x, hypot_y);
        a + b
    });
    bench!("hypot3", move |x: f32| hypot3(x, hypot_y, 2.0));
    bench!("rnorm3", move |x: f32| rnorm3(x, hypot_y, 2.0));
    bench!("normalize3", move |x: f32| {
        let (a, b, c) = normalize3(x, hypot_y, 2.0);
        a + b + c
    });
    bench!("hypot4", move |x: f32| hypot4(x, hypot_y, 2.0, 3.0));
    bench!("rnorm4", move |x: f32| rnorm4(x, hypot_y, 2.0, 3.0));
    bench!("normalize4", move |x: f32| {
        let (a, b, c, d) = normalize4(x, hypot_y, 2.0, 3.0);
        a + b + c + d
    });
    bench!("rsqrt", rsqrt);
    // black_box'd 2nd arg, same reasoning as atan2/hypot above -- a literal
    // exponent lets LLVM constant-fold powf's y==0.0/y_int/y_odd branches
    // away entirely (2.0 is a compile-time-known even integer), understating
    // the real branchy cost. See readme.md's own todo note on this.
    let powf_y = std::hint::black_box(2.0);
    bench!("powf", move |x: f32| powf(x, powf_y));
    bench!("srgb_to_linear", srgb_to_linear, Band::Half);
    bench!("linear_to_srgb", linear_to_srgb, Band::Half);
    bench!("std powf", move |x: f32| x.powf(powf_y));
    bench!("powf_unchecked", move |x: f32| powf_unchecked(x, powf_y));
    bench!("rootn", |x: f32| rootn(x, 3));
    // black_box'd 2nd arg, same reasoning as powf just above.
    let remainder_y = std::hint::black_box(3.0);
    bench!("remainder", move |x: f32| remainder(x, remainder_y));
    bench!("remainder_unchecked", move |x: f32| remainder_unchecked(x, remainder_y));
    bench!("remainder_checked", move |x: f32| remainder_checked(x, remainder_y));
    bench!("remainder_ieee", move |x: f32| remainder_ieee(x, remainder_y));
    bench!("remainder_wide", move |x: f32| remainder_wide(x, remainder_y));
    bench!("fmod", move |x: f32| fmod(x, remainder_y));
    bench!("fmod_unchecked", move |x: f32| fmod_unchecked(x, remainder_y));
}
