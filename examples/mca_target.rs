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

latency_fn!(lat_cbrt, "cbrt_latency", |x: f32| cbrt_normal(x, 1.0));
throughput_fn!(thr_cbrt, "cbrt_throughput", cbrt);

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

latency_fn!(lat_log2, "log2_latency", |x: f32| log_2_normal(x, 0.0));
throughput_fn!(thr_log2, "log2_throughput", log_2);

latency_fn!(lat_sin, "sin_latency", sin);
throughput_fn!(thr_sin, "sin_throughput", sin);

latency_fn!(lat_cos, "cos_latency", cos);
throughput_fn!(thr_cos, "cos_throughput", cos);

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

    run_all!(
        lat_nop, thr_nop;
        lat_cbrt, thr_cbrt;
        lat_cbrt_accurate, thr_cbrt_accurate;
        lat_cbrt_throughput_fn, thr_cbrt_throughput_fn;
        lat_cbrt_fast, thr_cbrt_fast;
        lat_exp2, thr_exp2;
        lat_exp2_checked, thr_exp2_checked;
        lat_log2, thr_log2;
        lat_sin, thr_sin;
        lat_cos, thr_cos;
    );
    black_box(&arr_out);
}
