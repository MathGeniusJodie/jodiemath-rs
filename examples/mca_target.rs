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

// Optional 4th argument is the input `Band`, defaulting to `Band::Two`
// (|x| in [2,4)) exactly as this macro used to hardcode. Unlike quickbench's,
// this file's band does not currently change any number: llvm-mca never sees
// a runtime value, and switching the domain-restricted rows to `Band::Half`
// left all ~150 regions **bit-identical**. The plausible mechanism for it to
// matter -- the band is a compile-time constant, so LLVM's known-bits
// analysis could in principle prove the chain's magnitude and fold an
// `|x|`-threshold select (asin's `a < 0.27`, atanh's `a < 0.25`) out of the
// measured region -- was measured and does not happen. The band is threaded
// through anyway to keep this file's chain identical to quickbench's, per
// mca_common.rs, and so that this stays true by construction rather than by
// luck if LLVM's reasoning ever sharpens.
macro_rules! latency_fn {
    ($fn_name:ident, $mca_name:literal, $f:expr) => {
        latency_fn!($fn_name, $mca_name, $f, Band::Two);
    };
    ($fn_name:ident, $mca_name:literal, $f:expr, $band:expr) => {
        #[inline(never)]
        pub fn $fn_name(mut x: f32) -> f32 {
            let f = $f;
            const BAND: Band = $band;
            // Seed into the band *before* the region marker, so the first of
            // the 64 chain steps sees an in-domain input like the other 63
            // without this costing an instruction inside the measured region.
            x = BAND.mix(x);
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
                x = BAND.mix(f(x));
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

// black_box, not a bare identity: Band::mix is idempotent once a value is
// already in the band, so `mix(identity(x))` repeated 64x collapses to a
// single mix at compile time (LLVM proves the chain is a no-op after the
// first step) and undercounts the very bookkeeping overhead this baseline
// is supposed to measure. black_box makes each step opaque so all 64 mixes
// actually execute.
latency_fn!(lat_nop, "nop_latency", |x: f32| black_box(x));
throughput_fn!(thr_nop, "nop_throughput", |x: f32| x);

latency_fn!(lat_fast_round_int, "fast_round_int_latency", fast_round_int);
throughput_fn!(
    thr_fast_round_int,
    "fast_round_int_throughput",
    fast_round_int
);
latency_fn!(lat_std_round, "std_round_latency", |x: f32| x.round());
throughput_fn!(thr_std_round, "std_round_throughput", |x: f32| x.round());

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
throughput_fn!(
    thr_cbrt_unchecked,
    "cbrt_unchecked_throughput",
    cbrt_unchecked
);

latency_fn!(lat_cbrt_accurate, "cbrt_accurate_latency", |x: f32| {
    cbrt_accurate_normal(x, 1.0)
});
throughput_fn!(thr_cbrt_accurate, "cbrt_accurate_throughput", cbrt_accurate);

// cbrt_accurate_unchecked == cbrt_accurate_normal(x, 1.0), so latency here
// is expected to match lat_cbrt_accurate above exactly (same reasoning as
// cbrt_unchecked) -- measured explicitly anyway for a complete row.
latency_fn!(
    lat_cbrt_accurate_unchecked,
    "cbrt_accurate_unchecked_latency",
    cbrt_accurate_unchecked
);
throughput_fn!(
    thr_cbrt_accurate_unchecked,
    "cbrt_accurate_unchecked_throughput",
    cbrt_accurate_unchecked
);

latency_fn!(lat_cbrt_fast, "cbrt_fast_latency", cbrt_fast);
throughput_fn!(thr_cbrt_fast, "cbrt_fast_throughput", cbrt_fast);

latency_fn!(lat_rcbrt, "rcbrt_latency", rcbrt);
throughput_fn!(thr_rcbrt, "rcbrt_throughput", rcbrt);

latency_fn!(lat_pow_3_2, "pow_3_2_latency", pow_3_2);
throughput_fn!(thr_pow_3_2, "pow_3_2_throughput", pow_3_2);

latency_fn!(lat_pow_2_3, "pow_2_3_latency", pow_2_3);
throughput_fn!(thr_pow_2_3, "pow_2_3_throughput", pow_2_3);

latency_fn!(
    lat_smoothstep,
    "smoothstep_latency",
    {
        let e0 = black_box(0.0);
        let e1 = black_box(1.0);
        move |x: f32| smoothstep(e0, e1, x)
    },
    Band::Half
);
throughput_fn!(thr_smoothstep, "smoothstep_throughput", {
    let e0 = black_box(0.0);
    let e1 = black_box(1.0);
    move |x: f32| smoothstep(e0, e1, x)
});

latency_fn!(
    lat_smootherstep,
    "smootherstep_latency",
    {
        let e0 = black_box(0.0);
        let e1 = black_box(1.0);
        move |x: f32| smootherstep(e0, e1, x)
    },
    Band::Half
);
throughput_fn!(thr_smootherstep, "smootherstep_throughput", {
    let e0 = black_box(0.0);
    let e1 = black_box(1.0);
    move |x: f32| smootherstep(e0, e1, x)
});

latency_fn!(lat_exp2, "exp2_latency", exp2);
throughput_fn!(thr_exp2, "exp2_throughput", exp2);

latency_fn!(lat_exp2_kf, "exp2_kf_latency", |f: f32| exp2_kf(3.0, f));
throughput_fn!(thr_exp2_kf, "exp2_kf_throughput", |f: f32| exp2_kf(3.0, f));

latency_fn!(lat_exp2_checked, "exp2_checked_latency", exp2_checked);
throughput_fn!(thr_exp2_checked, "exp2_checked_throughput", exp2_checked);

latency_fn!(lat_exp10, "exp10_latency", exp10);
throughput_fn!(thr_exp10, "exp10_throughput", exp10);

latency_fn!(lat_exp10_checked, "exp10_checked_latency", exp10_checked);
throughput_fn!(thr_exp10_checked, "exp10_checked_throughput", exp10_checked);

latency_fn!(lat_log2, "log2_latency", |x: f32| log_2_normal(x, 0.0));
throughput_fn!(thr_log2, "log2_throughput", log_2);

latency_fn!(
    lat_log2_unchecked,
    "log2_unchecked_latency",
    log_2_unchecked
);
throughput_fn!(
    thr_log2_unchecked,
    "log2_unchecked_throughput",
    log_2_unchecked
);

latency_fn!(lat_sin, "sin_latency", sin);
throughput_fn!(thr_sin, "sin_throughput", sin);

latency_fn!(lat_sin_wide, "sin_wide_latency", sin_wide);
throughput_fn!(thr_sin_wide, "sin_wide_throughput", sin_wide);

latency_fn!(lat_cos_wide, "cos_wide_latency", cos_wide);
throughput_fn!(thr_cos_wide, "cos_wide_throughput", cos_wide);

latency_fn!(lat_cos, "cos_latency", cos);
throughput_fn!(thr_cos, "cos_throughput", cos);

latency_fn!(lat_wrap_pi, "wrap_pi_latency", wrap_pi);
throughput_fn!(thr_wrap_pi, "wrap_pi_throughput", wrap_pi);

latency_fn!(lat_sin_prereduced, "sin_prereduced_latency", sin_prereduced);
throughput_fn!(
    thr_sin_prereduced,
    "sin_prereduced_throughput",
    sin_prereduced
);
latency_fn!(lat_cos_prereduced, "cos_prereduced_latency", cos_prereduced);
throughput_fn!(
    thr_cos_prereduced,
    "cos_prereduced_throughput",
    cos_prereduced
);

latency_fn!(lat_sinpi, "sinpi_latency", sinpi);
throughput_fn!(thr_sinpi, "sinpi_throughput", sinpi);

latency_fn!(lat_cospi, "cospi_latency", cospi);
throughput_fn!(thr_cospi, "cospi_throughput", cospi);

latency_fn!(lat_tanpi, "tanpi_latency", tanpi);
throughput_fn!(thr_tanpi, "tanpi_throughput", tanpi);

latency_fn!(lat_sin2pi, "sin2pi_latency", sin2pi);
throughput_fn!(thr_sin2pi, "sin2pi_throughput", sin2pi);

latency_fn!(lat_cos2pi, "cos2pi_latency", cos2pi);
throughput_fn!(thr_cos2pi, "cos2pi_throughput", cos2pi);

latency_fn!(lat_tan2pi, "tan2pi_latency", tan2pi);
throughput_fn!(thr_tan2pi, "tan2pi_throughput", tan2pi);

latency_fn!(
    lat_sinc_unnormalized,
    "sinc_unnormalized_latency",
    sinc_unnormalized
);
throughput_fn!(
    thr_sinc_unnormalized,
    "sinc_unnormalized_throughput",
    sinc_unnormalized
);

latency_fn!(lat_sind_unchecked, "sind_unchecked_latency", sind_unchecked);
throughput_fn!(
    thr_sind_unchecked,
    "sind_unchecked_throughput",
    sind_unchecked
);

latency_fn!(lat_cosd_unchecked, "cosd_unchecked_latency", cosd_unchecked);
throughput_fn!(
    thr_cosd_unchecked,
    "cosd_unchecked_throughput",
    cosd_unchecked
);

latency_fn!(lat_tand_unchecked, "tand_unchecked_latency", tand_unchecked);
throughput_fn!(
    thr_tand_unchecked,
    "tand_unchecked_throughput",
    tand_unchecked
);

latency_fn!(lat_ln, "ln_latency", ln);
throughput_fn!(thr_ln, "ln_throughput", ln);

latency_fn!(lat_ln_unchecked, "ln_unchecked_latency", ln_unchecked);
throughput_fn!(thr_ln_unchecked, "ln_unchecked_throughput", ln_unchecked);

latency_fn!(lat_log10, "log10_latency", log10);
throughput_fn!(thr_log10, "log10_throughput", log10);

latency_fn!(
    lat_log10_unchecked,
    "log10_unchecked_latency",
    log10_unchecked
);
throughput_fn!(
    thr_log10_unchecked,
    "log10_unchecked_throughput",
    log10_unchecked
);

latency_fn!(lat_log1p, "log1p_latency", log1p);
throughput_fn!(thr_log1p, "log1p_throughput", log1p);

latency_fn!(lat_log1pmx, "log1pmx_latency", log1pmx, Band::Half);
throughput_fn!(thr_log1pmx, "log1pmx_throughput", log1pmx);

latency_fn!(lat_log2p1, "log2p1_latency", log2p1);
throughput_fn!(thr_log2p1, "log2p1_throughput", log2p1);

latency_fn!(lat_log10p1, "log10p1_latency", log10p1);
throughput_fn!(thr_log10p1, "log10p1_throughput", log10p1);

latency_fn!(lat_exp, "exp_latency", exp);
throughput_fn!(thr_exp, "exp_throughput", exp);

latency_fn!(lat_exp_scaled, "exp_scaled_latency", |x: f32| exp_scaled(
    x, 3
));
throughput_fn!(
    thr_exp_scaled,
    "exp_scaled_throughput",
    |x: f32| exp_scaled(x, 3)
);

latency_fn!(lat_exp_narrow, "exp_narrow_latency", exp_narrow);
throughput_fn!(thr_exp_narrow, "exp_narrow_throughput", exp_narrow);

latency_fn!(lat_exp_checked, "exp_checked_latency", exp_checked);
throughput_fn!(thr_exp_checked, "exp_checked_throughput", exp_checked);

latency_fn!(lat_expm1, "expm1_latency", expm1);
throughput_fn!(thr_expm1, "expm1_throughput", expm1);

latency_fn!(lat_expm1_narrow, "expm1_narrow_latency", expm1_narrow);
throughput_fn!(thr_expm1_narrow, "expm1_narrow_throughput", expm1_narrow);

latency_fn!(lat_expm1_checked, "expm1_checked_latency", expm1_checked);
throughput_fn!(thr_expm1_checked, "expm1_checked_throughput", expm1_checked);

latency_fn!(
    lat_exp_m1_over_x_narrow,
    "exp_m1_over_x_narrow_latency",
    exp_m1_over_x_narrow
);
throughput_fn!(
    thr_exp_m1_over_x_narrow,
    "exp_m1_over_x_narrow_throughput",
    exp_m1_over_x_narrow
);

latency_fn!(lat_exp2m1, "exp2m1_latency", exp2m1);
throughput_fn!(thr_exp2m1, "exp2m1_throughput", exp2m1);

latency_fn!(lat_exp10m1, "exp10m1_latency", exp10m1);
throughput_fn!(thr_exp10m1, "exp10m1_throughput", exp10m1);

latency_fn!(lat_sinh, "sinh_latency", sinh);
throughput_fn!(thr_sinh, "sinh_throughput", sinh);

latency_fn!(lat_sinh_narrow, "sinh_narrow_latency", sinh_narrow);
throughput_fn!(thr_sinh_narrow, "sinh_narrow_throughput", sinh_narrow);

latency_fn!(lat_cosh, "cosh_latency", cosh);
throughput_fn!(thr_cosh, "cosh_throughput", cosh);

latency_fn!(lat_cosh_narrow, "cosh_narrow_latency", cosh_narrow);
throughput_fn!(thr_cosh_narrow, "cosh_narrow_throughput", cosh_narrow);

// "_fn" disambiguates these region names from sinh/cosh's own "_throughput"
// mode above, same convention already used for cbrt_throughput below.
latency_fn!(
    lat_sinh_throughput_fn,
    "sinh_throughput_fn_latency",
    sinh_throughput
);
throughput_fn!(
    thr_sinh_throughput_fn,
    "sinh_throughput_fn_throughput",
    sinh_throughput
);

latency_fn!(
    lat_cosh_throughput_fn,
    "cosh_throughput_fn_latency",
    cosh_throughput
);
throughput_fn!(
    thr_cosh_throughput_fn,
    "cosh_throughput_fn_throughput",
    cosh_throughput
);

latency_fn!(lat_sinh_checked, "sinh_checked_latency", sinh_checked);
throughput_fn!(thr_sinh_checked, "sinh_checked_throughput", sinh_checked);

latency_fn!(lat_cosh_checked, "cosh_checked_latency", cosh_checked);
throughput_fn!(thr_cosh_checked, "cosh_checked_throughput", cosh_checked);

latency_fn!(lat_coshm1, "coshm1_latency", coshm1);
throughput_fn!(thr_coshm1, "coshm1_throughput", coshm1);

latency_fn!(lat_tanh, "tanh_latency", tanh);
throughput_fn!(thr_tanh, "tanh_throughput", tanh);

latency_fn!(lat_tanh_grad, "tanh_grad_latency", tanh_grad);
throughput_fn!(thr_tanh_grad, "tanh_grad_throughput", tanh_grad);

latency_fn!(lat_sigmoid, "sigmoid_latency", sigmoid);
throughput_fn!(thr_sigmoid, "sigmoid_throughput", sigmoid);

latency_fn!(lat_sigmoid_fast, "sigmoid_fast_latency", sigmoid_fast);
throughput_fn!(thr_sigmoid_fast, "sigmoid_fast_throughput", sigmoid_fast);

latency_fn!(lat_sigmoid_grad, "sigmoid_grad_latency", sigmoid_grad);
throughput_fn!(thr_sigmoid_grad, "sigmoid_grad_throughput", sigmoid_grad);

latency_fn!(lat_logsigmoid, "logsigmoid_latency", logsigmoid);
throughput_fn!(thr_logsigmoid, "logsigmoid_throughput", logsigmoid);

latency_fn!(
    lat_logsigmoid_checked,
    "logsigmoid_checked_latency",
    logsigmoid_checked
);
throughput_fn!(
    thr_logsigmoid_checked,
    "logsigmoid_checked_throughput",
    logsigmoid_checked
);

latency_fn!(lat_gelu, "gelu_latency", gelu);
throughput_fn!(thr_gelu, "gelu_throughput", gelu);

latency_fn!(lat_silu, "silu_latency", silu);
throughput_fn!(thr_silu, "silu_throughput", silu);

latency_fn!(lat_silu_checked, "silu_checked_latency", silu_checked);
throughput_fn!(thr_silu_checked, "silu_checked_throughput", silu_checked);

latency_fn!(lat_softsign, "softsign_latency", softsign);
throughput_fn!(thr_softsign, "softsign_throughput", softsign);

latency_fn!(lat_sqrt1pm1, "sqrt1pm1_latency", sqrt1pm1);
throughput_fn!(thr_sqrt1pm1, "sqrt1pm1_throughput", sqrt1pm1);

latency_fn!(lat_asinh, "asinh_latency", asinh);
throughput_fn!(thr_asinh, "asinh_throughput", asinh);

latency_fn!(lat_acosh, "acosh_latency", acosh);
throughput_fn!(thr_acosh, "acosh_throughput", acosh);

latency_fn!(lat_atanh, "atanh_latency", atanh, Band::Half);
throughput_fn!(thr_atanh, "atanh_throughput", atanh);

latency_fn!(lat_asin, "asin_latency", asin, Band::Half);
throughput_fn!(thr_asin, "asin_throughput", asin);

latency_fn!(lat_asind, "asind_latency", asind, Band::Half);
throughput_fn!(thr_asind, "asind_throughput", asind);

latency_fn!(lat_asinpi, "asinpi_latency", asinpi, Band::Half);
throughput_fn!(thr_asinpi, "asinpi_throughput", asinpi);

latency_fn!(lat_acos, "acos_latency", acos, Band::Half);
throughput_fn!(thr_acos, "acos_throughput", acos);

latency_fn!(lat_acosd, "acosd_latency", acosd, Band::Half);
throughput_fn!(thr_acosd, "acosd_throughput", acosd);

latency_fn!(lat_acospi, "acospi_latency", acospi, Band::Half);
throughput_fn!(thr_acospi, "acospi_throughput", acospi);

latency_fn!(lat_atan, "atan_latency", atan);
throughput_fn!(thr_atan, "atan_throughput", atan);

latency_fn!(lat_atan_latency, "atan_latency_latency", atan_latency);
throughput_fn!(thr_atan_latency, "atan_latency_throughput", atan_latency);

latency_fn!(lat_atan_bounded, "atan_bounded_latency", atan_bounded);
throughput_fn!(thr_atan_bounded, "atan_bounded_throughput", atan_bounded);

latency_fn!(lat_atand, "atand_latency", atand);
throughput_fn!(thr_atand, "atand_throughput", atand);

latency_fn!(lat_atan2, "atan2_latency", |x: f32| atan2(x, 1.0));
throughput_fn!(thr_atan2, "atan2_throughput", |x: f32| atan2(x, 1.0));

latency_fn!(lat_atan2_unchecked, "atan2_unchecked_latency", |x: f32| {
    atan2_unchecked(x, 1.0)
});
throughput_fn!(
    thr_atan2_unchecked,
    "atan2_unchecked_throughput",
    |x: f32| atan2_unchecked(x, 1.0)
);

latency_fn!(lat_atan2_latency, "atan2_latency_latency", |x: f32| {
    atan2_latency(x, 1.0)
});
throughput_fn!(thr_atan2_latency, "atan2_latency_throughput", |x: f32| {
    atan2_latency(x, 1.0)
});

latency_fn!(lat_atan2_pos, "atan2_pos_latency", |x: f32| atan2_pos(
    x, 1.0
));
throughput_fn!(thr_atan2_pos, "atan2_pos_throughput", |x: f32| atan2_pos(
    x, 1.0
));

latency_fn!(lat_tan, "tan_latency", tan);
throughput_fn!(thr_tan, "tan_throughput", tan);

latency_fn!(lat_tan_wide, "tan_wide_latency", tan_wide);
throughput_fn!(thr_tan_wide, "tan_wide_throughput", tan_wide);

latency_fn!(lat_erf, "erf_latency", erf);
throughput_fn!(thr_erf, "erf_throughput", erf);

latency_fn!(lat_erfc, "erfc_latency", erfc);
throughput_fn!(thr_erfc, "erfc_throughput", erfc);

latency_fn!(lat_logit, "logit_latency", logit, Band::Half);
throughput_fn!(thr_logit, "logit_throughput", logit);

// xlogy/xlog1py must take a *varying* y, unlike atan2/hypot/powf above,
// which are fine with a loop-invariant second argument because their own
// expensive work still depends on x. Here it doesn't: xlogy(x,y) = x*ln(y)
// puts the entire transcendental on y alone, so a `black_box`ed invariant
// y let LLVM hoist the whole ln/log1p out of the loop and the regions
// measured little more than `x * precomputed_constant`. That is optimal
// codegen, not a de-vectorization bug -- but it made the numbers
// meaningless: xlog1py reported 9.03 cyc latency while log1p by itself is
// 52.19, and it also tripped codegen_check, whose scalar-vdivss check was
// seeing the hoisted loop-invariant divide rather than a per-element one.
// Passing x for both operands keeps every operand varying. No CSE hazard:
// x is a plain multiplier and y is the log's argument, so the two uses
// share no subexpression.
latency_fn!(lat_xlogy, "xlogy_latency", |x: f32| xlogy(x, x));
throughput_fn!(thr_xlogy, "xlogy_throughput", |x: f32| xlogy(x, x));
latency_fn!(lat_xlog1py, "xlog1py_latency", |x: f32| xlog1py(x, x));
throughput_fn!(thr_xlog1py, "xlog1py_throughput", |x: f32| xlog1py(x, x));

latency_fn!(lat_erfinv, "erfinv_latency", erfinv, Band::Half);
throughput_fn!(thr_erfinv, "erfinv_throughput", erfinv);

latency_fn!(lat_erfc_inv, "erfc_inv_latency", erfc_inv, Band::Half);
throughput_fn!(thr_erfc_inv, "erfc_inv_throughput", erfc_inv);

// cabs/carg (idea #186): exact aliases for hypot_checked/atan2, wired
// the same way those already are. cexp/clog return (f32,f32), which
// doesn't fit latency_fn!/throughput_fn!'s Fn(f32)->f32 shape -- see
// quickbench.rs for those instead (an adapter closure works there;
// clog's own branching also risks the multi-exit-path region-marker
// corruption already documented for ldexp/frexp/rootn, not
// worth the risk for a function whose cost is just its already-measured
// constituents).
// clog's *real* part only: `carg` is a plain `atan2` alias already measured
// on its own, and the near-1 branch (`0.5*log1p(re^2+im^2-1)`) is the part
// with any algebra in it. Both arms are computed and blended in the
// vectorized loop, so this region covers the whole value path. The second
// argument is derived from `x` rather than a constant so it cannot hoist.
latency_fn!(lat_clog_re, "clog_re_latency", |x: f32| clog(x, 1.0 - x).0);
throughput_fn!(thr_clog_re, "clog_re_throughput", |x: f32| clog(x, 1.0 - x)
    .0);

latency_fn!(lat_cabs, "cabs_latency", |x: f32| cabs(x, 1.0));
throughput_fn!(thr_cabs, "cabs_throughput", |x: f32| cabs(x, 1.0));
latency_fn!(lat_carg, "carg_latency", |x: f32| carg(x, 1.0));
throughput_fn!(thr_carg, "carg_throughput", |x: f32| carg(x, 1.0));

latency_fn!(lat_normalize2, "normalize2_latency", |x: f32| {
    let (a, b) = normalize2(x, 1.0);
    a + b
});
throughput_fn!(thr_normalize2, "normalize2_throughput", |x: f32| {
    let (a, b) = normalize2(x, 1.0);
    a + b
});

latency_fn!(lat_hypot3, "hypot3_latency", |x: f32| hypot3(x, 1.0, 2.0));
throughput_fn!(thr_hypot3, "hypot3_throughput", |x: f32| hypot3(
    x, 1.0, 2.0
));
latency_fn!(lat_rnorm3, "rnorm3_latency", |x: f32| rnorm3(x, 1.0, 2.0));
throughput_fn!(thr_rnorm3, "rnorm3_throughput", |x: f32| rnorm3(
    x, 1.0, 2.0
));

latency_fn!(lat_normalize3, "normalize3_latency", |x: f32| {
    let (a, b, c) = normalize3(x, 1.0, 2.0);
    a + b + c
});
throughput_fn!(thr_normalize3, "normalize3_throughput", |x: f32| {
    let (a, b, c) = normalize3(x, 1.0, 2.0);
    a + b + c
});

latency_fn!(lat_hypot4, "hypot4_latency", |x: f32| hypot4(
    x, 1.0, 2.0, 3.0
));
throughput_fn!(thr_hypot4, "hypot4_throughput", |x: f32| hypot4(
    x, 1.0, 2.0, 3.0
));
latency_fn!(lat_rnorm4, "rnorm4_latency", |x: f32| rnorm4(
    x, 1.0, 2.0, 3.0
));
throughput_fn!(thr_rnorm4, "rnorm4_throughput", |x: f32| rnorm4(
    x, 1.0, 2.0, 3.0
));
latency_fn!(lat_normalize4, "normalize4_latency", |x: f32| {
    let (a, b, c, d) = normalize4(x, 1.0, 2.0, 3.0);
    a + b + c + d
});
throughput_fn!(thr_normalize4, "normalize4_throughput", |x: f32| {
    let (a, b, c, d) = normalize4(x, 1.0, 2.0, 3.0);
    a + b + c + d
});

// black_box'd b/c/d: unlike hypot/rhypot's plain `1.0` above (fine there --
// no branch or constant-foldable sub-expression depends on it), diff_of_products
// computes w=c*d and e=fma(-c,d,w) from those operands alone -- with literal
// constants LLVM folds both to compile-time values and the chain collapses to
// one fma + one add, understating the real 1 mul + 2 fma + 1 add cost. Same
// reasoning as powf's `y` above.
latency_fn!(lat_diff_of_products, "diff_of_products_latency", {
    let b = black_box(1.7);
    let c = black_box(2.3);
    let d = black_box(0.9);
    move |x: f32| diff_of_products(x, b, c, d)
});
throughput_fn!(thr_diff_of_products, "diff_of_products_throughput", {
    let b = black_box(1.7);
    let c = black_box(2.3);
    let d = black_box(0.9);
    move |x: f32| diff_of_products(x, b, c, d)
});

latency_fn!(lat_cross2, "cross2_latency", {
    let ay = black_box(1.7);
    let bx = black_box(2.3);
    let by = black_box(0.9);
    move |ax: f32| cross2(ax, ay, bx, by)
});
throughput_fn!(thr_cross2, "cross2_throughput", {
    let ay = black_box(1.7);
    let bx = black_box(2.3);
    let by = black_box(0.9);
    move |ax: f32| cross2(ax, ay, bx, by)
});

latency_fn!(lat_rsqrt, "rsqrt_latency", rsqrt);
throughput_fn!(thr_rsqrt, "rsqrt_throughput", rsqrt);

// black_box'd 2nd arg, same reasoning as powf's `y` above and atan2/hypot
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

latency_fn!(lat_powf_pos, "powf_pos_latency", {
    let y = black_box(2.0);
    move |x: f32| powf_pos(x, y)
});
throughput_fn!(thr_powf_pos, "powf_pos_throughput", {
    let y = black_box(2.0);
    move |x: f32| powf_pos(x, y)
});

latency_fn!(
    lat_srgb_to_linear,
    "srgb_to_linear_latency",
    srgb_to_linear,
    Band::Half
);
throughput_fn!(
    thr_srgb_to_linear,
    "srgb_to_linear_throughput",
    srgb_to_linear
);

latency_fn!(
    lat_linear_to_srgb,
    "linear_to_srgb_latency",
    linear_to_srgb,
    Band::Half
);
throughput_fn!(
    thr_linear_to_srgb,
    "linear_to_srgb_throughput",
    linear_to_srgb
);

// signed_pow's mulsign is a branchless bit operation (not a runtime
// select on x's sign), so unlike erfcx's own mix()-sign-blindness
// caveat, this number is representative regardless of mix() erasing
// the sign each chain hop.
latency_fn!(lat_signed_pow, "signed_pow_latency", {
    let y = black_box(2.0);
    move |x: f32| signed_pow(x, y)
});
throughput_fn!(thr_signed_pow, "signed_pow_throughput", {
    let y = black_box(2.0);
    move |x: f32| signed_pow(x, y)
});

latency_fn!(lat_powf_unchecked, "powf_unchecked_latency", {
    let y = black_box(2.0);
    move |x: f32| powf_unchecked(x, y)
});
throughput_fn!(thr_powf_unchecked, "powf_unchecked_throughput", {
    let y = black_box(2.0);
    move |x: f32| powf_unchecked(x, y)
});

latency_fn!(lat_fmod, "fmod_latency", {
    let y = black_box(3.0);
    move |x: f32| fmod(x, y)
});
// fmod_throughput deliberately NOT wired up here -- same
// branch-specialization issue as remainder_throughput above (identical
// short-body shape, same fix), see that comment for the full mechanism.
latency_fn!(lat_fmod_checked, "fmod_checked_latency", {
    let y = black_box(3.0);
    move |x: f32| fmod_checked(x, y)
});
// fmod_checked_throughput deliberately NOT wired up either, same reason.
latency_fn!(lat_fmod_unchecked, "fmod_unchecked_latency", {
    let y = black_box(3.0);
    move |x: f32| fmod_unchecked(x, y)
});
throughput_fn!(thr_fmod_unchecked, "fmod_unchecked_throughput", {
    let y = black_box(3.0);
    move |x: f32| fmod_unchecked(x, y)
});

// rem_euclid/div_euclid: same short-body-over-fmod shape as
// fmod/fmod_checked above, so throughput deliberately NOT wired up
// either, same branch-specialization reasoning.
latency_fn!(lat_rem_euclid, "rem_euclid_latency", {
    let y = black_box(3.0);
    move |x: f32| rem_euclid(x, y)
});
latency_fn!(lat_div_euclid, "div_euclid_latency", {
    let y = black_box(3.0);
    move |x: f32| div_euclid(x, y)
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
        lat_fast_round_int, thr_fast_round_int;
        lat_std_round, thr_std_round;
        lat_cbrt, thr_cbrt;
        lat_cbrt_unchecked, thr_cbrt_unchecked;
        lat_cbrt_accurate, thr_cbrt_accurate;
        lat_cbrt_accurate_unchecked, thr_cbrt_accurate_unchecked;
        lat_cbrt_fast, thr_cbrt_fast;
        lat_rcbrt, thr_rcbrt;
        lat_pow_3_2, thr_pow_3_2;
        lat_pow_2_3, thr_pow_2_3;
        lat_smoothstep, thr_smoothstep;
        lat_smootherstep, thr_smootherstep;
        lat_exp2, thr_exp2;
        lat_exp2_kf, thr_exp2_kf;
        lat_exp2_checked, thr_exp2_checked;
        lat_exp10, thr_exp10;
        lat_exp10_checked, thr_exp10_checked;
        lat_log2, thr_log2;
        lat_log2_unchecked, thr_log2_unchecked;
        lat_sin, thr_sin;
        lat_sin_wide, thr_sin_wide;
        lat_cos_wide, thr_cos_wide;
        lat_cos, thr_cos;
        lat_wrap_pi, thr_wrap_pi;
        lat_sin_prereduced, thr_sin_prereduced;
        lat_cos_prereduced, thr_cos_prereduced;
        lat_sinpi, thr_sinpi;
        lat_cospi, thr_cospi;
        lat_tanpi, thr_tanpi;
        lat_sin2pi, thr_sin2pi;
        lat_cos2pi, thr_cos2pi;
        lat_tan2pi, thr_tan2pi;
        lat_sinc_unnormalized, thr_sinc_unnormalized;
        lat_sind_unchecked, thr_sind_unchecked;
        lat_cosd_unchecked, thr_cosd_unchecked;
        lat_tand_unchecked, thr_tand_unchecked;
        lat_ln, thr_ln;
        lat_ln_unchecked, thr_ln_unchecked;
        lat_log10, thr_log10;
        lat_log10_unchecked, thr_log10_unchecked;
        lat_log1p, thr_log1p;
        lat_log1pmx, thr_log1pmx;
        lat_log2p1, thr_log2p1;
        lat_log10p1, thr_log10p1;
        lat_exp, thr_exp;
        lat_exp_scaled, thr_exp_scaled;
        lat_exp_narrow, thr_exp_narrow;
        lat_exp_checked, thr_exp_checked;
        lat_expm1, thr_expm1;
        lat_expm1_narrow, thr_expm1_narrow;
        lat_expm1_checked, thr_expm1_checked;
        lat_exp_m1_over_x_narrow, thr_exp_m1_over_x_narrow;
        lat_exp2m1, thr_exp2m1;
        lat_exp10m1, thr_exp10m1;
        lat_sinh, thr_sinh;
        lat_sinh_narrow, thr_sinh_narrow;
        lat_cosh, thr_cosh;
        lat_cosh_narrow, thr_cosh_narrow;
        lat_sinh_throughput_fn, thr_sinh_throughput_fn;
        lat_cosh_throughput_fn, thr_cosh_throughput_fn;
        lat_sinh_checked, thr_sinh_checked;
        lat_cosh_checked, thr_cosh_checked;
        lat_coshm1, thr_coshm1;
        lat_tanh, thr_tanh;
        lat_tanh_grad, thr_tanh_grad;
        lat_sigmoid, thr_sigmoid;
        lat_sigmoid_fast, thr_sigmoid_fast;
        lat_sigmoid_grad, thr_sigmoid_grad;
        lat_logsigmoid, thr_logsigmoid;
        lat_logsigmoid_checked, thr_logsigmoid_checked;
        lat_gelu, thr_gelu;
        lat_silu, thr_silu;
        lat_silu_checked, thr_silu_checked;
        lat_softsign, thr_softsign;
        lat_sqrt1pm1, thr_sqrt1pm1;
        lat_asinh, thr_asinh;
        lat_acosh, thr_acosh;
        lat_atanh, thr_atanh;
        lat_asin, thr_asin;
        lat_asind, thr_asind;
        lat_asinpi, thr_asinpi;
        lat_acos, thr_acos;
        lat_acosd, thr_acosd;
        lat_acospi, thr_acospi;
        lat_atan, thr_atan;
        lat_atan_latency, thr_atan_latency;
        lat_atan_bounded, thr_atan_bounded;
        lat_atand, thr_atand;
        lat_atan2, thr_atan2;
        lat_atan2_unchecked, thr_atan2_unchecked;
        lat_atan2_latency, thr_atan2_latency;
        lat_atan2_pos, thr_atan2_pos;
        lat_tan, thr_tan;
        lat_tan_wide, thr_tan_wide;
        lat_erf, thr_erf;
        lat_erfc, thr_erfc;
        lat_logit, thr_logit;
        lat_xlogy, thr_xlogy;
        lat_xlog1py, thr_xlog1py;
        lat_erfinv, thr_erfinv;
        lat_erfc_inv, thr_erfc_inv;
        lat_clog_re, thr_clog_re;
        lat_cabs, thr_cabs;
        lat_carg, thr_carg;
        lat_normalize2, thr_normalize2;
        lat_hypot3, thr_hypot3;
        lat_rnorm3, thr_rnorm3;
        lat_normalize3, thr_normalize3;
        lat_hypot4, thr_hypot4;
        lat_rnorm4, thr_rnorm4;
        lat_normalize4, thr_normalize4;
        lat_diff_of_products, thr_diff_of_products;
        lat_cross2, thr_cross2;
        lat_rsqrt, thr_rsqrt;
        lat_powf, thr_powf;
        lat_powf_pos, thr_powf_pos;
        lat_srgb_to_linear, thr_srgb_to_linear;
        lat_linear_to_srgb, thr_linear_to_srgb;
        lat_signed_pow, thr_signed_pow;
        lat_powf_unchecked, thr_powf_unchecked;
        lat_fmod_unchecked, thr_fmod_unchecked;
    );
    // fmod: latency-only (see its own
    // throughput_fn! omission comments above for why), same standalone
    black_box(lat_fmod(black_box(1.234)));
    black_box(lat_fmod_checked(black_box(1.234)));
    black_box(lat_rem_euclid(black_box(1.234)));
    black_box(lat_div_euclid(black_box(1.234)));
    black_box(&arr_out);
}
