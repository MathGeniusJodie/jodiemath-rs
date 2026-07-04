// Shared between mca.rs and mca_target.rs via include!() -- they compile as
// separate binaries, so this is the one source of truth for the constants
// that tie the two together (the divisors in mca.rs must match the actual
// chain depth / array width baked into mca_target.rs's asm).
#[allow(dead_code)] // quickbench.rs includes this file for mix() but doesn't use these
pub const CHAIN_LEN: usize = 64;
#[allow(dead_code)]
pub const ARR_LEN: usize = 16;

// Keeps a serial dependency chain in a safe domain ([2, 4)) without a real
// branch. Shared by quickbench.rs and mca_target.rs so both harnesses clamp
// the same way.
#[inline(always)]
#[allow(dead_code)] // mca.rs includes this file for CHAIN_LEN/ARR_LEN but doesn't call mix()
pub fn mix(y: f32) -> f32 {
    f32::from_bits((y.to_bits() & 0x007f_ffff) | 0x4000_0000)
}
