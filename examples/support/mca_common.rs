// Shared between mca.rs and mca_target.rs via include!() -- they compile as
// separate binaries, so this is the one source of truth for the constants
// that tie the two together (the divisors in mca.rs must match the actual
// chain depth / array width baked into mca_target.rs's asm).
#[allow(dead_code)] // quickbench.rs includes this file for mix() but doesn't use these
pub const CHAIN_LEN: usize = 64;
#[allow(dead_code)]
pub const ARR_LEN: usize = 16;

// Keeps a serial dependency chain in a safe magnitude range (|x| in [2, 4))
// without a real branch. Shared by quickbench.rs and mca_target.rs so both
// harnesses clamp the same way. Total: the mask-and-or produces a normal
// value for *any* input bits, NaN and inf included (inf and zero both have a
// zero mantissa, so they come out as +-2.0).
//
// The sign bit is deliberately **carried through** rather than cleared
// (backlog idea #198). Clearing it let LLVM prove the chain's value was
// non-negative, and it then folded away the sign-handling work in every
// function whose tail does any -- `copysign`/`mulsign`, `abs`, sign selects
// -- silently deleting those instructions from the measured region. Carrying
// the sign makes it genuinely data-dependent and unknowable at compile time,
// so that work stays in.
//
// Note this is stronger than idea #198's own suggestion of injecting an
// *alternating* sign per chain step: the chain is `seq!`-unrolled, so a
// per-iteration constant would still be compile-time known and would fold
// exactly the same way. It has to depend on the data.
#[inline(always)]
#[allow(dead_code)] // mca.rs includes this file for CHAIN_LEN/ARR_LEN but doesn't call mix()
pub fn mix(y: f32) -> f32 {
    f32::from_bits((y.to_bits() & 0x807f_ffff) | 0x4000_0000)
}

/// Lowers this process's scheduling priority so a benchmark run doesn't
/// compete with foreground work on whatever machine it's run on -- mirrors
/// accuracy.rs's own `nice_self`. Nice values are inherited across
/// fork/exec, so calling this before mca.rs spawns `cargo rustc`/`llvm-mca`
/// nices those children too, not just mca.rs's own (otherwise idle) process
/// -- llvm-mca itself is the real CPU/memory hog in that pipeline.
#[cfg(unix)]
#[allow(dead_code)] // mca_target.rs includes this file but never runs standalone
pub fn nice_self() {
    // SAFETY: setpriority(PRIO_PROCESS, 0, _) only ever affects the calling
    // process's own niceness; failure just leaves the default priority.
    if unsafe { libc::setpriority(libc::PRIO_PROCESS, 0, 19) } != 0 {
        eprintln!("couldn't lower process priority (continuing anyway)");
    }
}
#[cfg(not(unix))]
#[allow(dead_code)]
pub fn nice_self() {}
