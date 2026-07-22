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

latency_fn!(lat_pow_3_2, "pow_3_2_latency", pow_3_2);
throughput_fn!(thr_pow_3_2, "pow_3_2_throughput", pow_3_2);

latency_fn!(lat_pow_2_3, "pow_2_3_latency", pow_2_3);
throughput_fn!(thr_pow_2_3, "pow_2_3_throughput", pow_2_3);

latency_fn!(lat_smoothstep, "smoothstep_latency", {
    let e0 = black_box(0.0);
    let e1 = black_box(1.0);
    move |x: f32| smoothstep(e0, e1, x)
});
throughput_fn!(thr_smoothstep, "smoothstep_throughput", {
    let e0 = black_box(0.0);
    let e1 = black_box(1.0);
    move |x: f32| smoothstep(e0, e1, x)
});

latency_fn!(lat_smootherstep, "smootherstep_latency", {
    let e0 = black_box(0.0);
    let e1 = black_box(1.0);
    move |x: f32| smootherstep(e0, e1, x)
});
throughput_fn!(thr_smootherstep, "smootherstep_throughput", {
    let e0 = black_box(0.0);
    let e1 = black_box(1.0);
    move |x: f32| smootherstep(e0, e1, x)
});

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

latency_fn!(lat_tanpi, "tanpi_latency", tanpi);
throughput_fn!(thr_tanpi, "tanpi_throughput", tanpi);

latency_fn!(lat_sinc, "sinc_latency", sinc);
throughput_fn!(thr_sinc, "sinc_throughput", sinc);

latency_fn!(lat_sinc_unnormalized, "sinc_unnormalized_latency", sinc_unnormalized);
throughput_fn!(thr_sinc_unnormalized, "sinc_unnormalized_throughput", sinc_unnormalized);

latency_fn!(lat_sind, "sind_latency", sind);
throughput_fn!(thr_sind, "sind_throughput", sind);

latency_fn!(lat_cosd, "cosd_latency", cosd);
throughput_fn!(thr_cosd, "cosd_throughput", cosd);

latency_fn!(lat_tand, "tand_latency", tand);
throughput_fn!(thr_tand, "tand_throughput", tand);

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

latency_fn!(lat_log1pmx, "log1pmx_latency", log1pmx);
throughput_fn!(thr_log1pmx, "log1pmx_throughput", log1pmx);

latency_fn!(lat_log2p1, "log2p1_latency", log2p1);
throughput_fn!(thr_log2p1, "log2p1_throughput", log2p1);

latency_fn!(lat_log10p1, "log10p1_latency", log10p1);
throughput_fn!(thr_log10p1, "log10p1_throughput", log10p1);

latency_fn!(lat_exp, "exp_latency", exp);
throughput_fn!(thr_exp, "exp_throughput", exp);

latency_fn!(lat_exp_narrow, "exp_narrow_latency", exp_narrow);
throughput_fn!(thr_exp_narrow, "exp_narrow_throughput", exp_narrow);

latency_fn!(lat_exp_checked, "exp_checked_latency", exp_checked);
throughput_fn!(thr_exp_checked, "exp_checked_throughput", exp_checked);

latency_fn!(lat_expm1, "expm1_latency", expm1);
throughput_fn!(thr_expm1, "expm1_throughput", expm1);

latency_fn!(lat_expm1_narrow, "expm1_narrow_latency", expm1_narrow);
throughput_fn!(thr_expm1_narrow, "expm1_narrow_throughput", expm1_narrow);

latency_fn!(lat_exp_m1_over_x, "exp_m1_over_x_latency", exp_m1_over_x);
throughput_fn!(thr_exp_m1_over_x, "exp_m1_over_x_throughput", exp_m1_over_x);

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
latency_fn!(lat_sinh_throughput_fn, "sinh_throughput_fn_latency", sinh_throughput);
throughput_fn!(thr_sinh_throughput_fn, "sinh_throughput_fn_throughput", sinh_throughput);

latency_fn!(lat_cosh_throughput_fn, "cosh_throughput_fn_latency", cosh_throughput);
throughput_fn!(thr_cosh_throughput_fn, "cosh_throughput_fn_throughput", cosh_throughput);

latency_fn!(lat_sinh_checked, "sinh_checked_latency", sinh_checked);
throughput_fn!(thr_sinh_checked, "sinh_checked_throughput", sinh_checked);

latency_fn!(lat_cosh_checked, "cosh_checked_latency", cosh_checked);
throughput_fn!(thr_cosh_checked, "cosh_checked_throughput", cosh_checked);

latency_fn!(lat_coshm1, "coshm1_latency", coshm1);
throughput_fn!(thr_coshm1, "coshm1_throughput", coshm1);

latency_fn!(lat_tanh, "tanh_latency", tanh);
throughput_fn!(thr_tanh, "tanh_throughput", tanh);

latency_fn!(lat_sigmoid, "sigmoid_latency", sigmoid);
throughput_fn!(thr_sigmoid, "sigmoid_throughput", sigmoid);

latency_fn!(lat_softplus, "softplus_latency", softplus);
throughput_fn!(thr_softplus, "softplus_throughput", softplus);

latency_fn!(lat_logaddexp, "logaddexp_latency", |x: f32| logaddexp(x, 0.0));
throughput_fn!(thr_logaddexp, "logaddexp_throughput", |x: f32| logaddexp(x, 0.0));

latency_fn!(lat_gelu, "gelu_latency", gelu);
throughput_fn!(thr_gelu, "gelu_throughput", gelu);

latency_fn!(lat_silu, "silu_latency", silu);
throughput_fn!(thr_silu, "silu_throughput", silu);

latency_fn!(lat_softsign, "softsign_latency", softsign);
throughput_fn!(thr_softsign, "softsign_throughput", softsign);

latency_fn!(lat_sqrt1pm1, "sqrt1pm1_latency", sqrt1pm1);
throughput_fn!(thr_sqrt1pm1, "sqrt1pm1_throughput", sqrt1pm1);

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

latency_fn!(lat_atan_bounded, "atan_bounded_latency", atan_bounded);
throughput_fn!(thr_atan_bounded, "atan_bounded_throughput", atan_bounded);

latency_fn!(lat_atan2, "atan2_latency", |x: f32| atan2(x, 1.0));
throughput_fn!(thr_atan2, "atan2_throughput", |x: f32| atan2(x, 1.0));

latency_fn!(lat_atan2_pos, "atan2_pos_latency", |x: f32| atan2_pos(x, 1.0));
throughput_fn!(thr_atan2_pos, "atan2_pos_throughput", |x: f32| atan2_pos(x, 1.0));

latency_fn!(lat_tan, "tan_latency", tan);
throughput_fn!(thr_tan, "tan_throughput", tan);

latency_fn!(lat_tan_checked, "tan_checked_latency", tan_checked);
throughput_fn!(thr_tan_checked, "tan_checked_throughput", tan_checked);

latency_fn!(lat_erf, "erf_latency", erf);
throughput_fn!(thr_erf, "erf_throughput", erf);

latency_fn!(lat_erfc, "erfc_latency", erfc);
throughput_fn!(thr_erfc, "erfc_throughput", erfc);

latency_fn!(lat_norm_cdf, "norm_cdf_latency", norm_cdf);
throughput_fn!(thr_norm_cdf, "norm_cdf_throughput", norm_cdf);

latency_fn!(lat_norm_pdf, "norm_pdf_latency", norm_pdf);
throughput_fn!(thr_norm_pdf, "norm_pdf_throughput", norm_pdf);

latency_fn!(lat_logit, "logit_latency", logit);
throughput_fn!(thr_logit, "logit_throughput", logit);

latency_fn!(lat_compound, "compound_latency", {
    let n = black_box(5.0);
    move |x: f32| compound(x, n)
});
throughput_fn!(thr_compound, "compound_throughput", {
    let n = black_box(5.0);
    move |x: f32| compound(x, n)
});

latency_fn!(lat_erfc_accurate, "erfc_accurate_latency", erfc_accurate);
throughput_fn!(thr_erfc_accurate, "erfc_accurate_throughput", erfc_accurate);

latency_fn!(lat_erfcx, "erfcx_latency", erfcx);
throughput_fn!(thr_erfcx, "erfcx_throughput", erfcx);

latency_fn!(lat_erfcx_accurate, "erfcx_accurate_latency", erfcx_accurate);
throughput_fn!(thr_erfcx_accurate, "erfcx_accurate_throughput", erfcx_accurate);

latency_fn!(lat_erfcx_checked, "erfcx_checked_latency", erfcx_checked);
throughput_fn!(thr_erfcx_checked, "erfcx_checked_throughput", erfcx_checked);

latency_fn!(lat_hypot, "hypot_latency", |x: f32| hypot(x, 1.0));
throughput_fn!(thr_hypot, "hypot_throughput", |x: f32| hypot(x, 1.0));

latency_fn!(lat_hypot_checked, "hypot_checked_latency", |x: f32| hypot_checked(x, 1.0));
throughput_fn!(thr_hypot_checked, "hypot_checked_throughput", |x: f32| hypot_checked(x, 1.0));

latency_fn!(lat_rhypot, "rhypot_latency", |x: f32| rhypot(x, 1.0));
throughput_fn!(thr_rhypot, "rhypot_throughput", |x: f32| rhypot(x, 1.0));

// black_box'd b/c/d: unlike hypot/rhypot's plain `1.0` above (fine there --
// no branch or constant-foldable sub-expression depends on it), diff_of_products
// computes w=c*d and e=fma(-c,d,w) from those operands alone -- with literal
// constants LLVM folds both to compile-time values and the chain collapses to
// one fma + one add, understating the real 1 mul + 2 fma + 1 add cost. Same
// reasoning as pown's `n`/powf's `y` above.
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

latency_fn!(lat_powf_pos, "powf_pos_latency", {
    let y = black_box(2.0);
    move |x: f32| powf_pos(x, y)
});
throughput_fn!(thr_powf_pos, "powf_pos_throughput", {
    let y = black_box(2.0);
    move |x: f32| powf_pos(x, y)
});

latency_fn!(lat_srgb_to_linear, "srgb_to_linear_latency", srgb_to_linear);
throughput_fn!(thr_srgb_to_linear, "srgb_to_linear_throughput", srgb_to_linear);

latency_fn!(lat_linear_to_srgb, "linear_to_srgb_latency", linear_to_srgb);
throughput_fn!(thr_linear_to_srgb, "linear_to_srgb_throughput", linear_to_srgb);

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

// pown_const<N>: idea #153, standing codegen verification (not just
// doc-comment "advice") that the const-generic exponent really does
// constant-fold the whole 32-iteration bit-testing loop away, for a
// representative small positive N (few bits set), a multi-bit positive N,
// and a negative N (reciprocal path). codegen_check.rs asserts these
// regions contain no branch/loop-back instruction at all.
latency_fn!(lat_pown_const3, "pown_const3_latency", pown_const::<3>);
latency_fn!(lat_pown_const7, "pown_const7_latency", pown_const::<7>);
latency_fn!(lat_pown_const_neg5, "pown_const_neg5_latency", pown_const::<-5>);

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
// remainder_throughput deliberately NOT wired up here (same reasoning as
// pown_small's own precedent above): after the x==0.0 zero/nan fix
// (backlog idea #85's own follow-up, 2026-07-09), LLVM branch-specializes
// this short function's vectorized loop on the shared black_box'd `y`,
// producing multiple physical return paths that each carry their own
// copy of this macro's inline-asm END marker, corrupting llvm-mca's
// region parser ("found an invalid region end directive"). codegen_check
// still confirms the real function has no scalar-fallback signatures --
// a harness limitation, not a code correctness issue. remainder_checked/
// remainder_wide (longer, already more complex bodies) don't hit this
// threshold and stay wired up normally. See quickbench.rs for
// remainder's own real wall-clock numbers instead.

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
// remainder_ieee_throughput deliberately NOT wired up here -- same
// branch-specialization issue as remainder_throughput above (identical
// short-body shape, same fix), see that comment for the full mechanism.
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
        lat_cbrt, thr_cbrt;
        lat_cbrt_unchecked, thr_cbrt_unchecked;
        lat_cbrt_accurate, thr_cbrt_accurate;
        lat_cbrt_accurate_unchecked, thr_cbrt_accurate_unchecked;
        lat_cbrt_throughput_fn, thr_cbrt_throughput_fn;
        lat_cbrt_fast, thr_cbrt_fast;
        lat_rcbrt, thr_rcbrt;
        lat_pow_3_2, thr_pow_3_2;
        lat_pow_2_3, thr_pow_2_3;
        lat_smoothstep, thr_smoothstep;
        lat_smootherstep, thr_smootherstep;
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
        lat_tanpi, thr_tanpi;
        lat_sinc, thr_sinc;
        lat_sinc_unnormalized, thr_sinc_unnormalized;
        lat_sind, thr_sind;
        lat_cosd, thr_cosd;
        lat_tand, thr_tand;
        lat_ln, thr_ln;
        lat_ln_unchecked, thr_ln_unchecked;
        lat_log10, thr_log10;
        lat_log10_unchecked, thr_log10_unchecked;
        lat_log1p, thr_log1p;
        lat_log1pmx, thr_log1pmx;
        lat_log2p1, thr_log2p1;
        lat_log10p1, thr_log10p1;
        lat_exp, thr_exp;
        lat_exp_narrow, thr_exp_narrow;
        lat_exp_checked, thr_exp_checked;
        lat_expm1, thr_expm1;
        lat_expm1_narrow, thr_expm1_narrow;
        lat_exp_m1_over_x, thr_exp_m1_over_x;
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
        lat_sigmoid, thr_sigmoid;
        lat_softplus, thr_softplus;
        lat_logaddexp, thr_logaddexp;
        lat_gelu, thr_gelu;
        lat_silu, thr_silu;
        lat_softsign, thr_softsign;
        lat_sqrt1pm1, thr_sqrt1pm1;
        lat_asinh, thr_asinh;
        lat_acosh, thr_acosh;
        lat_atanh, thr_atanh;
        lat_asin, thr_asin;
        lat_acos, thr_acos;
        lat_atan, thr_atan;
        lat_atan_latency, thr_atan_latency;
        lat_atan_bounded, thr_atan_bounded;
        lat_atan2, thr_atan2;
        lat_atan2_pos, thr_atan2_pos;
        lat_tan, thr_tan;
        lat_tan_checked, thr_tan_checked;
        lat_erf, thr_erf;
        lat_erfc, thr_erfc;
        lat_norm_cdf, thr_norm_cdf;
        lat_norm_pdf, thr_norm_pdf;
        lat_logit, thr_logit;
        lat_compound, thr_compound;
        lat_erfc_accurate, thr_erfc_accurate;
        lat_erfcx, thr_erfcx;
        lat_erfcx_accurate, thr_erfcx_accurate;
        lat_erfcx_checked, thr_erfcx_checked;
        lat_hypot, thr_hypot;
        lat_hypot_checked, thr_hypot_checked;
        lat_rhypot, thr_rhypot;
        lat_diff_of_products, thr_diff_of_products;
        lat_cross2, thr_cross2;
        lat_rsqrt, thr_rsqrt;
        lat_powf, thr_powf;
        lat_powf_pos, thr_powf_pos;
        lat_srgb_to_linear, thr_srgb_to_linear;
        lat_linear_to_srgb, thr_linear_to_srgb;
        lat_signed_pow, thr_signed_pow;
        lat_powf_unchecked, thr_powf_unchecked;
        lat_pown, thr_pown;
        lat_powf_checked, thr_powf_checked;
        lat_powf_checked_unchecked, thr_powf_checked_unchecked;
        lat_remainder_unchecked, thr_remainder_unchecked;
        lat_fmod_unchecked, thr_fmod_unchecked;
        lat_remainder_checked, thr_remainder_checked;
        lat_remainder_wide, thr_remainder_wide;
    );
    // remainder/remainder_ieee/fmod: latency-only (see their own
    // throughput_fn! omission comments above for why), same standalone
    // call shape as lat_cbrt_wrapped above.
    black_box(lat_remainder(black_box(1.234)));
    black_box(lat_remainder_ieee(black_box(1.234)));
    black_box(lat_fmod(black_box(1.234)));
    black_box(lat_fmod_checked(black_box(1.234)));
    black_box(lat_rem_euclid(black_box(1.234)));
    black_box(lat_div_euclid(black_box(1.234)));
    black_box(lat_pown_const3(black_box(1.234)));
    black_box(lat_pown_const7(black_box(1.234)));
    black_box(lat_pown_const_neg5(black_box(1.234)));
    black_box(&arr_out);
}
