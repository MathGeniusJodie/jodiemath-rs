// Shared between mca.rs and mca_target.rs via include!() -- they compile as
// separate binaries, so this is the one source of truth for the constants
// that tie the two together (the divisors in mca.rs must match the actual
// chain depth / array width baked into mca_target.rs's asm).
#[allow(dead_code)] // quickbench.rs includes this file for mix() but doesn't use these
pub const CHAIN_LEN: usize = 64;
#[allow(dead_code)]
pub const ARR_LEN: usize = 16;

/// The magnitude band a benchmark pins its inputs to -- the dependency
/// chain's `mix` step for latency, the input array for throughput.
///
/// This exists because a latency chain feeds a function its own output, so
/// the band has to sit **inside that function's domain**. It didn't, for a
/// while: the band was hardcoded to `Two` and `asin`/`acos`/`atanh` were
/// benched entirely on `|x| in [2,4)`, where they are undefined. libm's own
/// asinf/acosf/atanhf take their domain check and return NaN in a couple of
/// nanoseconds, while this crate's branchless versions compute the full
/// result regardless -- so readme.md published those three as "0.2x vs std"
/// when in-domain they are 1.1-1.4x *ahead* of std. Same for the throughput
/// array, in the other direction (it undersold them by ~2x). Pick the band
/// from the function's domain, not by default.
///
/// Domain here means the domain of the *chain*, not just of one call:
/// `log1pmx` fits in `[2,4)` fine but returns a negative value for every
/// input, and `mix` carries the sign, so it walked itself out past its own
/// `x > -1` edge in one hop. quickbench's `check_in_domain` exists because
/// that is not visible in a signature.
///
/// Tell-tale that a row has this problem: std's latency comes out roughly
/// equal to std's own throughput. A real dependency chain costs several
/// times its throughput; equality means the chain is short-circuiting.
#[derive(Clone, Copy, PartialEq, Eq)]
#[allow(dead_code)] // each harness only names the variants it actually needs
pub enum Band {
    /// `|x|` in `[2, 4)`. The default, and correct for any function that is
    /// defined on all reals or on `x >= 1` (`acosh`, the logs, the exps...).
    Two,
    /// `|x|` in `[0.5, 1)`. For functions defined only on `[-1, 1]` /
    /// `(-1, 1)` / `(0, 1)`, and for `erfc_inv`'s `(0, 2)`. Lands on the
    /// main branch of the ones that split on `|x|`.
    Half,
    /// `|x|` in `[0.125, 0.25)`. Same domains as `Half`, but lands on the
    /// small-`|x|` branch instead (`asin_small`, `atanh_small`, ...).
    Eighth,
    /// `|x|` in `[2^30, 2^31)`. For wide-range functions (`sin_wide`, `cos_wide`,
    /// `tan_wide`) whose purpose is large-magnitude reduction where standard
    /// Cody-Waite reduction breaks down and full Payne-Hanek reduction is required.
    Large,
}

impl Band {
    /// The f32 exponent field `mix` ORs in. The mantissa is taken from the
    /// data and the sign is carried, so only the exponent is fixed.
    #[inline(always)]
    const fn exponent_bits(self) -> u32 {
        match self {
            Band::Two => 0x4000_0000,
            Band::Half => 0x3f00_0000,
            Band::Eighth => 0x3e00_0000,
            Band::Large => 0x4e80_0000,
        }
    }

    /// Half-open magnitude range `[lo, hi)`, for building a throughput input
    /// array that covers the same band the latency chain runs in.
    #[inline(always)]
    #[allow(dead_code)] // mca_target.rs's throughput regions take an opaque array
    pub const fn range(self) -> (f32, f32) {
        match self {
            Band::Two => (2.0, 4.0),
            Band::Half => (0.5, 1.0),
            Band::Eighth => (0.125, 0.25),
            Band::Large => (1073741824.0, 2147483648.0),
        }
    }

    /// Keeps a serial dependency chain in this band without a real branch.
    /// Shared by quickbench.rs and mca_target.rs so both harnesses clamp the
    /// same way. Total: the mask-and-or produces a normal value for *any*
    /// input bits, NaN and inf included (inf and zero both have a zero
    /// mantissa, so they come out as +-the band's low end).
    ///
    /// The sign bit is deliberately **carried through** rather than cleared
    /// (backlog idea #198). Clearing it let LLVM prove the chain's value was
    /// non-negative, and it then folded away the sign-handling work in every
    /// function whose tail does any -- `copysign`/`mulsign`, `abs`, sign
    /// selects -- silently deleting those instructions from the measured
    /// region. Carrying the sign makes it genuinely data-dependent and
    /// unknowable at compile time, so that work stays in.
    ///
    /// Note this is stronger than idea #198's own suggestion of injecting an
    /// *alternating* sign per chain step: the chain is `seq!`-unrolled, so a
    /// per-iteration constant would still be compile-time known and would
    /// fold exactly the same way. It has to depend on the data.
    ///
    /// The carried sign is also why `Half`/`Eighth` are safe for the
    /// one-sided domains (`logit`, `probit`, `erfc_inv`): each of those maps
    /// its positive half to a positive result, so a chain seeded positive
    /// stays positive. `quickbench`'s `check_in_domain` verifies that rather
    /// than trusting it.
    #[inline(always)]
    #[allow(dead_code)] // mca.rs includes this file for CHAIN_LEN/ARR_LEN only
    pub fn mix(self, y: f32) -> f32 {
        f32::from_bits((y.to_bits() & 0x807f_ffff) | self.exponent_bits())
    }
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
