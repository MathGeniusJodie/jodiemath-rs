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

// idea #78: N independent serial dependency chains interleaved in the
// same loop, instead of bench_latency's single chain. A single chain's
// own measured latency can be masked or exaggerated by how much
// cross-call ILP the CPU's out-of-order engine can extract when there's
// only one outstanding dependency chain to fill scheduling gaps with --
// 4 independent chains approximates a caller that has some real
// concurrent work available, closer to how throughput callers (many
// independent elements) sit between this and the single-chain extreme.
const N_STREAMS: usize = 4;
fn bench_latency_n(name: &str, f: impl Fn(f32) -> f32) {
    let mut best = f64::INFINITY;
    for _ in 0..REPS {
        let mut xs = [1.234_f32, 1.876, 1.456, 1.987];
        let start = Instant::now();
        for _ in 0..(LAT_ITERS / N_STREAMS as u64) {
            for x in xs.iter_mut() {
                *x = mix(f(*x));
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
    if args.get(1).map(|s| s.as_str()) == Some("latencyn") {
        // idea #78 spot-check: single-chain vs 4-independent-chain latency
        // for a spread of functions with different codegen shapes (pown's
        // fully-unrolled branchy loop, exp2's short balanced-Estrin chain,
        // sinh's branch-selected two-branch combine, acos's plain Horner).
        bench_latency("pown", |x: f32| pown(x, black_box(5)));
        bench_latency_n("pown", |x: f32| pown(x, black_box(5)));
        bench_latency("exp2", exp2);
        bench_latency_n("exp2", exp2);
        bench_latency("sinh", sinh);
        bench_latency_n("sinh", sinh);
        bench_latency("acos", acos);
        bench_latency_n("acos", acos);
        return;
    }
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
    bench!("cbrt_unchecked", cbrt_unchecked);
    bench!("cbrt_accurate", cbrt_accurate);
    bench!("cbrt_accurate_unchecked", cbrt_accurate_unchecked);
    bench!("cbrt_throughput", cbrt_throughput);
    bench!("cbrt_fast", cbrt_fast);
    bench!("rcbrt", rcbrt);
    bench!("pow_3_2", pow_3_2);
    bench!("pow_2_3", pow_2_3);
    bench!("smoothstep", move |x: f32| smoothstep(0.0, 1.0, x));
    bench!("smootherstep", move |x: f32| smootherstep(0.0, 1.0, x));
    bench!("exp2", exp2);
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
    bench!("std cos", |x: f32| x.cos());
    bench!("sinpi", sinpi);
    bench!("cospi", cospi);
    bench!("tanpi", tanpi);
    bench!("sin2pi", sin2pi);
    bench!("cos2pi", cos2pi);
    bench!("tan2pi", tan2pi);
    bench!("sinc", sinc);
    bench!("sinc_unnormalized", sinc_unnormalized);
    bench!("sind", sind);
    bench!("cosd", cosd);
    bench!("tand", tand);

    bench!("ln", ln);
    bench!("ln_unchecked", ln_unchecked);
    bench!("std ln", |x: f32| x.ln());
    bench!("log10", log10);
    bench!("log10_unchecked", log10_unchecked);
    bench!("std log10", |x: f32| x.log10());
    bench!("log1p", log1p);
    bench!("log1pmx", log1pmx);
    bench!("std log1p", |x: f32| x.ln_1p());
    bench!("log2p1", log2p1);
    bench!("exp", exp);
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
    bench!("sigmoid", sigmoid);
    bench!("softplus", softplus);
    bench!("logaddexp", |x: f32| logaddexp(x, 0.0));
    bench!("asinh", asinh);
    bench!("std asinh", |x: f32| x.asinh());
    bench!("acosh", acosh);
    bench!("std acosh", |x: f32| x.acosh());
    bench!("atanh", atanh);
    bench!("std atanh", |x: f32| x.atanh());
    bench!("asin", asin);
    bench!("asind", asind);
    bench!("std asin", |x: f32| x.asin());
    bench!("acos", acos);
    bench!("std acos", |x: f32| x.acos());
    bench!("acosd", acosd);
    bench!("atan", atan);
    bench!("std atan", |x: f32| x.atan());
    bench!("atan_latency", atan_latency);
    bench!("atan_bounded", atan_bounded);
    bench!("atand", atand);
    // black_box'd 2nd arg (not a literal 1.0): a compile-time-constant 2nd
    // arg lets LLVM fold away atan2's own special-case branches entirely,
    // silently hiding their real cost -- this matters here specifically
    // because it would otherwise make atan2_unchecked look like a wash
    // instead of the real win it is. See readme.md's own todo note on
    // this exact fixed-argument limitation.
    let atan2_x2 = std::hint::black_box(1.0);
    bench!("atan2", move |x: f32| atan2(x, atan2_x2));
    bench!("std atan2", move |x: f32| x.atan2(atan2_x2));
    bench!("atan2_unchecked", move |x: f32| atan2_unchecked(x, atan2_x2));
    bench!("atan2_pos", move |x: f32| atan2_pos(x, atan2_x2));
    bench!("atan2d", move |x: f32| atan2d(x, atan2_x2));
    bench!("tan", tan);
    bench!("tan_checked", tan_checked);
    bench!("std tan", |x: f32| x.tan());
    bench!("erf", erf);
    bench!("erfc", erfc);
    bench!("norm_cdf", norm_cdf);
    bench!("norm_pdf", norm_pdf);
    bench!("logit", logit);
    bench!("compound", move |x: f32| compound(x, 5.0));
    bench!("erfcx", erfcx);
    // black_box'd 2nd arg, same reasoning as atan2 above.
    let hypot_y = std::hint::black_box(1.0);
    bench!("hypot", move |x: f32| hypot(x, hypot_y));
    bench!("std hypot", move |x: f32| x.hypot(hypot_y));
    bench!("hypot_unchecked", move |x: f32| hypot_unchecked(x, hypot_y));
    bench!("hypot_checked", move |x: f32| hypot_checked(x, hypot_y));
    bench!("rhypot", move |x: f32| rhypot(x, hypot_y));
    bench!("rsqrt", rsqrt);
    // black_box'd 2nd arg, same reasoning as atan2/hypot above -- a literal
    // exponent lets LLVM constant-fold powf's y==0.0/y_int/y_odd branches
    // away entirely (2.0 is a compile-time-known even integer), understating
    // the real branchy cost. See readme.md's own todo note on this.
    let powf_y = std::hint::black_box(2.0);
    bench!("powf", move |x: f32| powf(x, powf_y));
    bench!("srgb_to_linear", srgb_to_linear);
    bench!("linear_to_srgb", linear_to_srgb);
    bench!("std powf", move |x: f32| x.powf(powf_y));
    bench!("powf_unchecked", move |x: f32| powf_unchecked(x, powf_y));
    bench!("powf_checked", move |x: f32| powf_checked(x, powf_y));
    bench!("powf_checked_unchecked", move |x: f32| powf_checked_unchecked(x, powf_y));
    {
        // black_box'd once, not per-call -- see pown's own mca_target.rs
        // comment for why that placement matters.
        let n = black_box(5);
        bench!("pown", |x: f32| pown(x, n));
        bench!("pown_small", |x: f32| pown_small(x, n));
        bench!("pown_small_accurate", |x: f32| pown_small_accurate(x, n));
    }
    // N baked in at compile time (not black_box'd -- that's the whole
    // point of pown_const, unlike pown/pown_small above).
    bench!("pown_const<5>", |x: f32| pown_const::<5>(x));
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
