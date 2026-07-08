// Marked functions for llvm-mca: each wraps a function under test in a
// `# LLVM-MCA-BEGIN <name>`/`# LLVM-MCA-END` comment pair, which llvm-mca
// recognizes as a named code region when fed the compiler's `--emit=asm`
// output (examples/mca.rs drives that pipeline). This file is never
// meaningfully "run" -- it only needs to compile and keep every function
// alive so the region markers survive into the .s file.
//
// llvm-mca has no branch predictor: it treats a marked region as one
// literal instruction stream, so a real conditional branch inside it gets
// both sides counted as if they always execute. Two edge cases in this
// crate (cbrt's denormal rescale, log_2's denormal/zero/inf checks) compile
// to actual jumps in *scalar* code even though the crate's whole design
// goal is branchless selects -- that branchless guarantee only holds once
// LLVM if-converts a real vectorized loop (see readme.md). So:
//   - latency (serial dependency chain, necessarily scalar): call the
//     branchless `*_normal` core directly, skipping the rare-input edge
//     branch entirely. This matches quickbench's own methodology, which
//     also only ever feeds in-domain inputs so the branch predictor learns
//     one direction and the untaken path costs ~0 cycles on real hardware.
//   - throughput (array of independent values): use the real public
//     function inside an actual `for (o, &x) in ...zip...` loop, same idiom
//     as quickbench, so LLVM's loop vectorizer if-converts the edge checks
//     into masked selects -- genuinely branchless, no special-casing needed.
// std::f32 functions aren't included here: they call out to external libm
// symbols llvm-mca can't see inside; use quickbench for std comparisons.
use jodiemath_rs::*;
use seq_macro::seq;
use std::hint::black_box;

include!("support/mca_common.rs");
// seq! needs a literal range, not a const, so the `0..64` below can't read
// CHAIN_LEN directly -- this ties them together so drift is a compile error
// instead of silently wrong numbers.
const _: () = assert!(CHAIN_LEN == 64);

macro_rules! latency_fn {
    ($fn_name:ident, $mca_name:literal, $f:expr) => {
        #[inline(never)]
        pub fn $fn_name(mut x: f32) -> f32 {
            let f = $f;
            // A bare marker comment has no data dependency, so LLVM is free
            // to hoist/sink pure register arithmetic across it (verified:
            // without the xmm0 operand below, most of the chain leaked
            // outside the region, undercounting latency by ~30x). Routing
            // x through an explicit xmm0 in/out operand forces a real SSA
            // dependency edge at exactly the BEGIN/END boundary.
            unsafe {
                core::arch::asm!(concat!("# LLVM-MCA-BEGIN ", $mca_name), inout("xmm0") x, options(nostack, preserves_flags));
            }
            seq!(_i in 0..64 {
                x = mix(f(x));
            });
            unsafe {
                core::arch::asm!("# LLVM-MCA-END", inout("xmm0") x, options(nostack, preserves_flags));
            }
            x
        }
    };
}

macro_rules! throughput_fn {
    ($fn_name:ident, $mca_name:literal, $f:expr) => {
        #[inline(never)]
        pub fn $fn_name(input: &[f32; ARR_LEN], output: &mut [f32; ARR_LEN]) {
            let f = $f;
            unsafe { core::arch::asm!(concat!("# LLVM-MCA-BEGIN ", $mca_name)) };
            for (o, &x) in output.iter_mut().zip(input.iter()) {
                *o = f(x);
            }
            unsafe { core::arch::asm!("# LLVM-MCA-END") };
        }
    };
}

// black_box, not a bare identity: mix() is idempotent once a value is
// already in [2,4), so `mix(identity(x))` repeated 64x collapses to a
// single mix at compile time (LLVM proves the chain is a no-op after the
// first step) and undercounts the very bookkeeping overhead this baseline
// is supposed to measure. black_box makes each step opaque so all 64 mixes
// actually execute.
latency_fn!(lat_nop, "nop_latency", |x: f32| black_box(x));
throughput_fn!(thr_nop, "nop_throughput", |x: f32| x);

latency_fn!(lat_cbrt, "cbrt_latency", cbrt_normal);
throughput_fn!(thr_cbrt, "cbrt_throughput", cbrt);

// full public cbrt (with its tiny-input select), not just the _normal core --
// this is the fair "real function, all edge-case handling included" latency
// number, since lat_cbrt above intentionally skips the rare-input branch.
latency_fn!(lat_cbrt_wrapped, "cbrt_wrapped_latency", cbrt);

latency_fn!(lat_cbrt_accurate, "cbrt_accurate_latency", |x: f32| cbrt_accurate_normal(x, 1.0));
throughput_fn!(thr_cbrt_accurate, "cbrt_accurate_throughput", cbrt_accurate);

latency_fn!(lat_cbrt_throughput_fn, "cbrt_throughput_fn_latency", cbrt_throughput);
throughput_fn!(thr_cbrt_throughput_fn, "cbrt_throughput_fn_throughput", cbrt_throughput);

latency_fn!(lat_cbrt_fast, "cbrt_fast_latency", cbrt_fast);
throughput_fn!(thr_cbrt_fast, "cbrt_fast_throughput", cbrt_fast);

latency_fn!(lat_exp2, "exp2_latency", exp2);
throughput_fn!(thr_exp2, "exp2_throughput", exp2);

latency_fn!(lat_exp2_checked, "exp2_checked_latency", exp2_checked);
throughput_fn!(thr_exp2_checked, "exp2_checked_throughput", exp2_checked);

latency_fn!(lat_exp10, "exp10_latency", exp10);
throughput_fn!(thr_exp10, "exp10_throughput", exp10);

latency_fn!(lat_exp10_checked, "exp10_checked_latency", exp10_checked);
throughput_fn!(thr_exp10_checked, "exp10_checked_throughput", exp10_checked);

latency_fn!(lat_log2, "log2_latency", |x: f32| log_2_normal(x, 0.0));
throughput_fn!(thr_log2, "log2_throughput", log_2);

latency_fn!(lat_sin, "sin_latency", sin);
throughput_fn!(thr_sin, "sin_throughput", sin);

latency_fn!(lat_sin_checked, "sin_checked_latency", sin_checked);
throughput_fn!(thr_sin_checked, "sin_checked_throughput", sin_checked);

latency_fn!(lat_cos, "cos_latency", cos);
throughput_fn!(thr_cos, "cos_throughput", cos);

latency_fn!(lat_cos_checked, "cos_checked_latency", cos_checked);
throughput_fn!(thr_cos_checked, "cos_checked_throughput", cos_checked);

latency_fn!(lat_sinpi, "sinpi_latency", sinpi);
throughput_fn!(thr_sinpi, "sinpi_throughput", sinpi);

latency_fn!(lat_cospi, "cospi_latency", cospi);
throughput_fn!(thr_cospi, "cospi_throughput", cospi);

latency_fn!(lat_sind, "sind_latency", sind);
throughput_fn!(thr_sind, "sind_throughput", sind);

latency_fn!(lat_cosd, "cosd_latency", cosd);
throughput_fn!(thr_cosd, "cosd_throughput", cosd);

latency_fn!(lat_ln, "ln_latency", ln);
throughput_fn!(thr_ln, "ln_throughput", ln);

latency_fn!(lat_log10, "log10_latency", log10);
throughput_fn!(thr_log10, "log10_throughput", log10);

latency_fn!(lat_log1p, "log1p_latency", log1p);
throughput_fn!(thr_log1p, "log1p_throughput", log1p);

latency_fn!(lat_exp, "exp_latency", exp);
throughput_fn!(thr_exp, "exp_throughput", exp);

latency_fn!(lat_expm1, "expm1_latency", expm1);
throughput_fn!(thr_expm1, "expm1_throughput", expm1);

latency_fn!(lat_sinh, "sinh_latency", sinh);
throughput_fn!(thr_sinh, "sinh_throughput", sinh);

latency_fn!(lat_cosh, "cosh_latency", cosh);
throughput_fn!(thr_cosh, "cosh_throughput", cosh);

// "_fn" disambiguates these region names from sinh/cosh's own "_throughput"
// mode above, same convention already used for cbrt_throughput below.
latency_fn!(lat_sinh_throughput_fn, "sinh_throughput_fn_latency", sinh_throughput);
throughput_fn!(thr_sinh_throughput_fn, "sinh_throughput_fn_throughput", sinh_throughput);

latency_fn!(lat_cosh_throughput_fn, "cosh_throughput_fn_latency", cosh_throughput);
throughput_fn!(thr_cosh_throughput_fn, "cosh_throughput_fn_throughput", cosh_throughput);

latency_fn!(lat_tanh, "tanh_latency", tanh);
throughput_fn!(thr_tanh, "tanh_throughput", tanh);

latency_fn!(lat_asinh, "asinh_latency", asinh);
throughput_fn!(thr_asinh, "asinh_throughput", asinh);

latency_fn!(lat_acosh, "acosh_latency", acosh);
throughput_fn!(thr_acosh, "acosh_throughput", acosh);

latency_fn!(lat_atanh, "atanh_latency", atanh);
throughput_fn!(thr_atanh, "atanh_throughput", atanh);

latency_fn!(lat_asin, "asin_latency", asin);
throughput_fn!(thr_asin, "asin_throughput", asin);

latency_fn!(lat_acos, "acos_latency", acos);
throughput_fn!(thr_acos, "acos_throughput", acos);

latency_fn!(lat_atan, "atan_latency", atan);
throughput_fn!(thr_atan, "atan_throughput", atan);

latency_fn!(lat_atan2, "atan2_latency", |x: f32| atan2(x, 1.0));
throughput_fn!(thr_atan2, "atan2_throughput", |x: f32| atan2(x, 1.0));

latency_fn!(lat_tan, "tan_latency", tan);
throughput_fn!(thr_tan, "tan_throughput", tan);

latency_fn!(lat_erf, "erf_latency", erf);
throughput_fn!(thr_erf, "erf_throughput", erf);

latency_fn!(lat_erfc, "erfc_latency", erfc);
throughput_fn!(thr_erfc, "erfc_throughput", erfc);

latency_fn!(lat_hypot, "hypot_latency", |x: f32| hypot(x, 1.0));
throughput_fn!(thr_hypot, "hypot_throughput", |x: f32| hypot(x, 1.0));

latency_fn!(lat_powf, "powf_latency", |x: f32| powf(x, 2.0));
throughput_fn!(thr_powf, "powf_throughput", |x: f32| powf(x, 2.0));

latency_fn!(lat_powf_checked, "powf_checked_latency", |x: f32| powf_checked(x, 2.0));
throughput_fn!(thr_powf_checked, "powf_checked_throughput", |x: f32| powf_checked(
    x, 2.0
));

latency_fn!(lat_remainder, "remainder_latency", |x: f32| remainder(x, 3.0));
throughput_fn!(thr_remainder, "remainder_throughput", |x: f32| remainder(x, 3.0));

latency_fn!(lat_remainder_checked, "remainder_checked_latency", |x: f32| remainder_checked(
    x, 3.0
));
throughput_fn!(thr_remainder_checked, "remainder_checked_throughput", |x: f32| {
    remainder_checked(x, 3.0)
});

fn main() {
    // smoke test only: exercises every marked function once so `cargo run
    // --release --example mca_target` succeeds on its own. The interesting
    // artifact is the assembly, produced separately by examples/mca.rs.
    let arr_in: [f32; ARR_LEN] = std::array::from_fn(|i| 2.0 + i as f32 * (2.0 / ARR_LEN as f32));
    let mut arr_out = [0f32; ARR_LEN];

    macro_rules! run_all {
        ($($lat:ident, $thr:ident);* $(;)?) => {
            $(
                black_box($lat(black_box(1.234)));
                $thr(black_box(&arr_in), black_box(&mut arr_out));
            )*
        };
    }

    black_box(lat_cbrt_wrapped(black_box(1.234)));

    run_all!(
        lat_nop, thr_nop;
        lat_cbrt, thr_cbrt;
        lat_cbrt_accurate, thr_cbrt_accurate;
        lat_cbrt_throughput_fn, thr_cbrt_throughput_fn;
        lat_cbrt_fast, thr_cbrt_fast;
        lat_exp2, thr_exp2;
        lat_exp2_checked, thr_exp2_checked;
        lat_exp10, thr_exp10;
        lat_exp10_checked, thr_exp10_checked;
        lat_log2, thr_log2;
        lat_sin, thr_sin;
        lat_sin_checked, thr_sin_checked;
        lat_cos, thr_cos;
        lat_cos_checked, thr_cos_checked;
        lat_sinpi, thr_sinpi;
        lat_cospi, thr_cospi;
        lat_sind, thr_sind;
        lat_cosd, thr_cosd;
        lat_ln, thr_ln;
        lat_log10, thr_log10;
        lat_log1p, thr_log1p;
        lat_exp, thr_exp;
        lat_expm1, thr_expm1;
        lat_sinh, thr_sinh;
        lat_cosh, thr_cosh;
        lat_sinh_throughput_fn, thr_sinh_throughput_fn;
        lat_cosh_throughput_fn, thr_cosh_throughput_fn;
        lat_tanh, thr_tanh;
        lat_asinh, thr_asinh;
        lat_acosh, thr_acosh;
        lat_atanh, thr_atanh;
        lat_asin, thr_asin;
        lat_acos, thr_acos;
        lat_atan, thr_atan;
        lat_atan2, thr_atan2;
        lat_tan, thr_tan;
        lat_erf, thr_erf;
        lat_erfc, thr_erfc;
        lat_hypot, thr_hypot;
        lat_powf, thr_powf;
        lat_powf_checked, thr_powf_checked;
        lat_remainder, thr_remainder;
        lat_remainder_checked, thr_remainder_checked;
    );
    black_box(&arr_out);
}
