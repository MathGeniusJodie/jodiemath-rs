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

// cbrt_unchecked == cbrt_normal, so latency here is expected to match
// lat_cbrt above exactly (this crate's own convention: latency calls the
// branchless *_normal core directly already) -- measured explicitly
// anyway for a complete row, same as log2_unchecked's own pair below.
latency_fn!(lat_cbrt_unchecked, "cbrt_unchecked_latency", cbrt_unchecked);
throughput_fn!(thr_cbrt_unchecked, "cbrt_unchecked_throughput", cbrt_unchecked);

latency_fn!(lat_cbrt_accurate, "cbrt_accurate_latency", |x: f32| cbrt_accurate_normal(x, 1.0));
throughput_fn!(thr_cbrt_accurate, "cbrt_accurate_throughput", cbrt_accurate);

// cbrt_accurate_unchecked == cbrt_accurate_normal(x, 1.0), so latency here
// is expected to match lat_cbrt_accurate above exactly (same reasoning as
// cbrt_unchecked) -- measured explicitly anyway for a complete row.
latency_fn!(lat_cbrt_accurate_unchecked, "cbrt_accurate_unchecked_latency", cbrt_accurate_unchecked);
throughput_fn!(thr_cbrt_accurate_unchecked, "cbrt_accurate_unchecked_throughput", cbrt_accurate_unchecked);

latency_fn!(lat_cbrt_throughput_fn, "cbrt_throughput_fn_latency", cbrt_throughput);
throughput_fn!(thr_cbrt_throughput_fn, "cbrt_throughput_fn_throughput", cbrt_throughput);

latency_fn!(lat_cbrt_fast, "cbrt_fast_latency", cbrt_fast);
throughput_fn!(thr_cbrt_fast, "cbrt_fast_throughput", cbrt_fast);

latency_fn!(lat_rcbrt, "rcbrt_latency", rcbrt);
throughput_fn!(thr_rcbrt, "rcbrt_throughput", rcbrt);

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

latency_fn!(lat_log2_unchecked, "log2_unchecked_latency", log_2_unchecked);
throughput_fn!(thr_log2_unchecked, "log2_unchecked_throughput", log_2_unchecked);

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

latency_fn!(lat_sinc, "sinc_latency", sinc);
throughput_fn!(thr_sinc, "sinc_throughput", sinc);

latency_fn!(lat_sind, "sind_latency", sind);
throughput_fn!(thr_sind, "sind_throughput", sind);

latency_fn!(lat_cosd, "cosd_latency", cosd);
throughput_fn!(thr_cosd, "cosd_throughput", cosd);

latency_fn!(lat_ln, "ln_latency", ln);
throughput_fn!(thr_ln, "ln_throughput", ln);

latency_fn!(lat_ln_unchecked, "ln_unchecked_latency", ln_unchecked);
throughput_fn!(thr_ln_unchecked, "ln_unchecked_throughput", ln_unchecked);

latency_fn!(lat_log10, "log10_latency", log10);
throughput_fn!(thr_log10, "log10_throughput", log10);

latency_fn!(lat_log10_unchecked, "log10_unchecked_latency", log10_unchecked);
throughput_fn!(thr_log10_unchecked, "log10_unchecked_throughput", log10_unchecked);

latency_fn!(lat_log1p, "log1p_latency", log1p);
throughput_fn!(thr_log1p, "log1p_throughput", log1p);

latency_fn!(lat_log2p1, "log2p1_latency", log2p1);
throughput_fn!(thr_log2p1, "log2p1_throughput", log2p1);

latency_fn!(lat_exp, "exp_latency", exp);
throughput_fn!(thr_exp, "exp_throughput", exp);

latency_fn!(lat_expm1, "expm1_latency", expm1);
throughput_fn!(thr_expm1, "expm1_throughput", expm1);

latency_fn!(lat_exp_m1_over_x, "exp_m1_over_x_latency", exp_m1_over_x);
throughput_fn!(thr_exp_m1_over_x, "exp_m1_over_x_throughput", exp_m1_over_x);

latency_fn!(lat_exp2m1, "exp2m1_latency", exp2m1);
throughput_fn!(thr_exp2m1, "exp2m1_throughput", exp2m1);

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

latency_fn!(lat_sigmoid, "sigmoid_latency", sigmoid);
throughput_fn!(thr_sigmoid, "sigmoid_throughput", sigmoid);

latency_fn!(lat_softplus, "softplus_latency", softplus);
throughput_fn!(thr_softplus, "softplus_throughput", softplus);

latency_fn!(lat_logaddexp, "logaddexp_latency", |x: f32| logaddexp(x, 0.0));
throughput_fn!(thr_logaddexp, "logaddexp_throughput", |x: f32| logaddexp(x, 0.0));

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

latency_fn!(lat_atan_latency, "atan_latency_latency", atan_latency);
throughput_fn!(thr_atan_latency, "atan_latency_throughput", atan_latency);

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

latency_fn!(lat_hypot_checked, "hypot_checked_latency", |x: f32| hypot_checked(x, 1.0));
throughput_fn!(thr_hypot_checked, "hypot_checked_throughput", |x: f32| hypot_checked(x, 1.0));

latency_fn!(lat_rhypot, "rhypot_latency", |x: f32| rhypot(x, 1.0));
throughput_fn!(thr_rhypot, "rhypot_throughput", |x: f32| rhypot(x, 1.0));

latency_fn!(lat_rsqrt, "rsqrt_latency", rsqrt);
throughput_fn!(thr_rsqrt, "rsqrt_throughput", rsqrt);

// black_box'd 2nd arg, same reasoning as pown's `n` below and atan2/hypot
// above: a literal `2.0` exponent is a compile-time-known even integer,
// letting LLVM fold away powf's y==0.0/y_int/y_odd branches entirely and
// understating its real branchy cost.
latency_fn!(lat_powf, "powf_latency", {
    let y = black_box(2.0);
    move |x: f32| powf(x, y)
});
throughput_fn!(thr_powf, "powf_throughput", {
    let y = black_box(2.0);
    move |x: f32| powf(x, y)
});

latency_fn!(lat_powf_unchecked, "powf_unchecked_latency", {
    let y = black_box(2.0);
    move |x: f32| powf_unchecked(x, y)
});
throughput_fn!(thr_powf_unchecked, "powf_unchecked_throughput", {
    let y = black_box(2.0);
    move |x: f32| powf_unchecked(x, y)
});

// black_box'd n, computed *once* before the loop/chain (not per-call
// inside the closure -- that placement let LLVM see through it and
// simplify pown's fixed 32-iteration loop away, undercounting the
// realistic "runtime exponent, unknown at compile time" cost this
// function is actually built for, see its own doc comment).
latency_fn!(lat_pown, "pown_latency", {
    let n = black_box(5);
    move |x: f32| pown(x, n)
});
throughput_fn!(thr_pown, "pown_throughput", {
    let n = black_box(5);
    move |x: f32| pown(x, n)
});
// pown_small deliberately NOT wired up here: with only 8 unrolled
// iterations (vs pown's 32), LLVM's cost model chooses to branch-
// specialize on the shared black_box'd `n` this harness uses (cheap
// enough to be worth it at this trip count, unlike pown's 32) instead of
// emitting the uniform blend/select pown gets -- confirmed correct and
// still fully vectorized for the harder, realistic per-lane-varying-`n`
// case (checked directly via a standalone --emit=asm probe: proper
// AVX-512 masked selects, no scalar fallback), but the branch-specialized
// shared-n form has multiple return paths, each carrying its own copy of
// this macro's inline-asm END marker, which corrupts llvm-mca's region
// parser ("found an invalid region end directive"). A harness limitation
// specific to this trip count + shared-n combination, not a code
// correctness issue -- see quickbench.rs for pown_small's real wall-clock
// numbers instead.

latency_fn!(lat_powf_checked, "powf_checked_latency", {
    let y = black_box(2.0);
    move |x: f32| powf_checked(x, y)
});
throughput_fn!(thr_powf_checked, "powf_checked_throughput", {
    let y = black_box(2.0);
    move |x: f32| powf_checked(x, y)
});

latency_fn!(lat_powf_checked_unchecked, "powf_checked_unchecked_latency", {
    let y = black_box(2.0);
    move |x: f32| powf_checked_unchecked(x, y)
});
throughput_fn!(thr_powf_checked_unchecked, "powf_checked_unchecked_throughput", {
    let y = black_box(2.0);
    move |x: f32| powf_checked_unchecked(x, y)
});

// black_box'd 2nd arg, same reasoning as powf above.
latency_fn!(lat_remainder, "remainder_latency", {
    let y = black_box(3.0);
    move |x: f32| remainder(x, y)
});
throughput_fn!(thr_remainder, "remainder_throughput", {
    let y = black_box(3.0);
    move |x: f32| remainder(x, y)
});

latency_fn!(lat_remainder_unchecked, "remainder_unchecked_latency", {
    let y = black_box(3.0);
    move |x: f32| remainder_unchecked(x, y)
});
throughput_fn!(thr_remainder_unchecked, "remainder_unchecked_throughput", {
    let y = black_box(3.0);
    move |x: f32| remainder_unchecked(x, y)
});

latency_fn!(lat_remainder_checked, "remainder_checked_latency", {
    let y = black_box(3.0);
    move |x: f32| remainder_checked(x, y)
});
throughput_fn!(thr_remainder_checked, "remainder_checked_throughput", {
    let y = black_box(3.0);
    move |x: f32| remainder_checked(x, y)
});
latency_fn!(lat_remainder_ieee, "remainder_ieee_latency", {
    let y = black_box(3.0);
    move |x: f32| remainder_ieee(x, y)
});
throughput_fn!(thr_remainder_ieee, "remainder_ieee_throughput", {
    let y = black_box(3.0);
    move |x: f32| remainder_ieee(x, y)
});
latency_fn!(lat_remainder_wide, "remainder_wide_latency", {
    let y = black_box(3.0);
    move |x: f32| remainder_wide(x, y)
});
throughput_fn!(thr_remainder_wide, "remainder_wide_throughput", {
    let y = black_box(3.0);
    move |x: f32| remainder_wide(x, y)
});

latency_fn!(lat_fmod, "fmod_latency", {
    let y = black_box(3.0);
    move |x: f32| fmod(x, y)
});
throughput_fn!(thr_fmod, "fmod_throughput", {
    let y = black_box(3.0);
    move |x: f32| fmod(x, y)
});
latency_fn!(lat_fmod_unchecked, "fmod_unchecked_latency", {
    let y = black_box(3.0);
    move |x: f32| fmod_unchecked(x, y)
});
throughput_fn!(thr_fmod_unchecked, "fmod_unchecked_throughput", {
    let y = black_box(3.0);
    move |x: f32| fmod_unchecked(x, y)
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
        lat_cbrt_unchecked, thr_cbrt_unchecked;
        lat_cbrt_accurate, thr_cbrt_accurate;
        lat_cbrt_accurate_unchecked, thr_cbrt_accurate_unchecked;
        lat_cbrt_throughput_fn, thr_cbrt_throughput_fn;
        lat_cbrt_fast, thr_cbrt_fast;
        lat_rcbrt, thr_rcbrt;
        lat_exp2, thr_exp2;
        lat_exp2_checked, thr_exp2_checked;
        lat_exp10, thr_exp10;
        lat_exp10_checked, thr_exp10_checked;
        lat_log2, thr_log2;
        lat_log2_unchecked, thr_log2_unchecked;
        lat_sin, thr_sin;
        lat_sin_checked, thr_sin_checked;
        lat_cos, thr_cos;
        lat_cos_checked, thr_cos_checked;
        lat_sinpi, thr_sinpi;
        lat_cospi, thr_cospi;
        lat_sinc, thr_sinc;
        lat_sind, thr_sind;
        lat_cosd, thr_cosd;
        lat_ln, thr_ln;
        lat_ln_unchecked, thr_ln_unchecked;
        lat_log10, thr_log10;
        lat_log10_unchecked, thr_log10_unchecked;
        lat_log1p, thr_log1p;
        lat_log2p1, thr_log2p1;
        lat_exp, thr_exp;
        lat_expm1, thr_expm1;
        lat_exp_m1_over_x, thr_exp_m1_over_x;
        lat_exp2m1, thr_exp2m1;
        lat_sinh, thr_sinh;
        lat_cosh, thr_cosh;
        lat_sinh_throughput_fn, thr_sinh_throughput_fn;
        lat_cosh_throughput_fn, thr_cosh_throughput_fn;
        lat_tanh, thr_tanh;
        lat_sigmoid, thr_sigmoid;
        lat_softplus, thr_softplus;
        lat_logaddexp, thr_logaddexp;
        lat_asinh, thr_asinh;
        lat_acosh, thr_acosh;
        lat_atanh, thr_atanh;
        lat_asin, thr_asin;
        lat_acos, thr_acos;
        lat_atan, thr_atan;
        lat_atan_latency, thr_atan_latency;
        lat_atan2, thr_atan2;
        lat_tan, thr_tan;
        lat_erf, thr_erf;
        lat_erfc, thr_erfc;
        lat_hypot, thr_hypot;
        lat_hypot_checked, thr_hypot_checked;
        lat_rhypot, thr_rhypot;
        lat_rsqrt, thr_rsqrt;
        lat_powf, thr_powf;
        lat_powf_unchecked, thr_powf_unchecked;
        lat_pown, thr_pown;
        lat_powf_checked, thr_powf_checked;
        lat_powf_checked_unchecked, thr_powf_checked_unchecked;
        lat_remainder, thr_remainder;
        lat_remainder_unchecked, thr_remainder_unchecked;
        lat_fmod, thr_fmod;
        lat_fmod_unchecked, thr_fmod_unchecked;
        lat_remainder_checked, thr_remainder_checked;
        lat_remainder_ieee, thr_remainder_ieee;
        lat_remainder_wide, thr_remainder_wide;
    );
    black_box(&arr_out);
}
