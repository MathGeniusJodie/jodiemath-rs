// Standalone scalar wrappers for `--emit=asm` + llvm-mca latency timelines.
use jodiemath_rs::*;
#[no_mangle]
#[inline(never)]
pub fn scalar_sin_wide(x: f32) -> f32 {
    sin_wide(x)
}
#[no_mangle]
#[inline(never)]
pub fn scalar_cos_wide(x: f32) -> f32 {
    cos_wide(x)
}
#[no_mangle]
#[inline(never)]
pub fn scalar_tan_wide(x: f32) -> f32 {
    tan_wide(x)
}
fn main() {
    let x = std::hint::black_box(1.5f32);
    println!(
        "{}",
        scalar_sin_wide(x) + scalar_cos_wide(x) + scalar_tan_wide(x)
    );
}
