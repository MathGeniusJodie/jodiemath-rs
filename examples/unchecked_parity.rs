// IDEAS.md backlog round 3 #12: "Standing test: every `_unchecked` bit-
// matches its checked sibling on the documented domain (several doc
// comments promise this; nothing enforces it)." `edgecheck.rs` only pins
// 2-3 spot values per pair, not a systematic sweep -- this fuzzes each
// pair across the exact domain closures `examples/accuracy.rs` already
// uses for its own ulp reporting (reused verbatim to stay consistent with
// the crate's established domain definitions), asserting bit-for-bit
// equality rather than measuring accuracy. Run with:
//   cargo run --release --example unchecked_parity
// (debug builds are slow for 50M samples/pair but not refused outright,
// unlike accuracy.rs's exhaustive mode -- this is fuzz-only, not 2^32).
use jodiemath_rs::*;

struct Rng(u64);
impl Rng {
    fn next_u32(&mut self) -> u32 {
        // xorshift64
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 32) as u32
    }
    fn next_f32(&mut self) -> f32 {
        f32::from_bits(self.next_u32())
    }
}

#[must_use]
fn check1(name: &str, n: u64, domain: impl Fn(f32) -> bool, checked: impl Fn(f32) -> f32, unchecked: impl Fn(f32) -> f32) -> bool {
    let mut rng = Rng(0x9E3779B97F4A7C15 ^ (name.len() as u64 + 1));
    let mut checked_count = 0u64;
    let mut mismatches = 0u64;
    let mut first_mismatch: Option<(f32, f32, f32)> = None;
    for _ in 0..n {
        let x = rng.next_f32();
        if !domain(x) {
            continue;
        }
        checked_count += 1;
        let a = checked(x);
        let b = unchecked(x);
        if a.to_bits() != b.to_bits() && !(a.is_nan() && b.is_nan()) {
            mismatches += 1;
            if first_mismatch.is_none() {
                first_mismatch = Some((x, a, b));
            }
        }
    }
    if mismatches > 0 {
        let (x, a, b) = first_mismatch.unwrap();
        println!(
            "MISMATCH {name}: {mismatches}/{checked_count} in-domain samples differ; first at x={x:e} checked={a:e} (0x{:08x}) unchecked={b:e} (0x{:08x})",
            a.to_bits(), b.to_bits()
        );
        false
    } else {
        println!("ok       {name}: {checked_count} in-domain samples, bit-identical");
        true
    }
}

#[must_use]
fn check2(
    name: &str,
    n: u64,
    domain: impl Fn(f32, f32) -> bool,
    checked: impl Fn(f32, f32) -> f32,
    unchecked: impl Fn(f32, f32) -> f32,
) -> bool {
    let mut rng = Rng(0x2545F4914F6CDD1D ^ (name.len() as u64 + 1));
    let mut checked_count = 0u64;
    let mut mismatches = 0u64;
    let mut first_mismatch: Option<(f32, f32, f32, f32)> = None;
    for _ in 0..n {
        let x = rng.next_f32();
        let y = rng.next_f32();
        if !domain(x, y) {
            continue;
        }
        checked_count += 1;
        let a = checked(x, y);
        let b = unchecked(x, y);
        if a.to_bits() != b.to_bits() && !(a.is_nan() && b.is_nan()) {
            mismatches += 1;
            if first_mismatch.is_none() {
                first_mismatch = Some((x, y, a, b));
            }
        }
    }
    if mismatches > 0 {
        let (x, y, a, b) = first_mismatch.unwrap();
        println!(
            "MISMATCH {name}: {mismatches}/{checked_count} in-domain samples differ; first at x={x:e} y={y:e} checked={a:e} (0x{:08x}) unchecked={b:e} (0x{:08x})",
            a.to_bits(), b.to_bits()
        );
        false
    } else {
        println!("ok       {name}: {checked_count} in-domain samples, bit-identical");
        true
    }
}

fn main() {
    const N: u64 = 50_000_000;
    let mut ok = true;

    let positive_normal = |x: f32| x >= f32::MIN_POSITIVE && x.is_finite();
    ok &= check1("log_2 / log_2_unchecked", N, positive_normal, log_2, log_2_unchecked);
    ok &= check1("ln / ln_unchecked", N, positive_normal, ln, ln_unchecked);
    ok &= check1("log10 / log10_unchecked", N, positive_normal, log10, log10_unchecked);

    let normal_finite = |x: f32| x.abs() >= f32::MIN_POSITIVE && x.is_finite();
    ok &= check1("cbrt / cbrt_unchecked", N, normal_finite, cbrt, cbrt_unchecked);

    let accurate_safe_range = |x: f32| {
        let ax = x.to_bits() & 0x7fff_ffff;
        ax >= 0x2380_0000 && ax < 0x7f00_0000
    };
    ok &= check1(
        "cbrt_accurate / cbrt_accurate_unchecked",
        N,
        accurate_safe_range,
        cbrt_accurate,
        cbrt_accurate_unchecked,
    );

    let atan2_domain = |x: f32, y: f32| x != 0.0 && !(x.is_infinite() && y.is_infinite());
    ok &= check2("atan2 / atan2_unchecked", N, atan2_domain, atan2, atan2_unchecked);

    // fmod/remainder correct an exact-cancellation sign bug their own
    // _unchecked twins don't (see remainder_style_combine!'s own
    // comment): when nonzero x is an exact multiple of y, `-q*y` exactly
    // cancels x, and IEEE754 exact-cancellation always gives +0.0
    // regardless of the "should be x's sign" convention -- fmod/
    // remainder patch this with an explicit copysign, but the
    // _unchecked twins skip it (same class of omission their own doc
    // comments already document for x==0.0 itself, just one case wider:
    // a nonzero x whose *result* happens to land on exactly zero).
    let exact_multiple = |x: f32, y: f32| (-(x / y).round()).mul_add(y, x) == 0.0;
    let rem_domain = |x: f32, y: f32| x != 0.0 && y.is_finite() && !exact_multiple(x, y);
    ok &= check2("fmod / fmod_unchecked", N, rem_domain, fmod, fmod_unchecked);
    ok &= check2("remainder / remainder_unchecked", N, rem_domain, remainder, remainder_unchecked);

    let pow_domain = |x: f32, y: f32| {
        x >= f32::MIN_POSITIVE && x.is_finite() && y != 0.0 && (-126.0..128.0).contains(&(x.log2() * y))
    };
    ok &= check2("powf / powf_unchecked", N, pow_domain, powf, powf_unchecked);

    // remainder_wide's own contract: bit-identical to remainder_checked
    // (not to plain remainder) throughout remainder_checked's own
    // |x/y|<2^24 domain -- a different pair than remainder/
    // remainder_unchecked already above, never added to this standing
    // test either. No near_tie exclusion needed here (unlike accuracy.rs's
    // own remainder_checked sweep): that exclusion is about comparing
    // against a *reference*, not about whether these two implementations
    // agree with *each other*. Two exceptions used to be excluded here:
    // an exact half-integer x/y tie flipping sign (found at 500M
    // samples), and a rescale-near-f32::MAX guard pushing an
    // already-tiny x into denormal-underflow territory (found at 30M
    // samples, up to ~4 ulp). Both are now fixed (round_ties_even
    // instead of round for adj; gating the rescale on |x| alone instead
    // of max(|x|,|y|), since q0*y tracks x regardless of y's own
    // magnitude -- see remainder_wide's own doc comment for both) and
    // verified with no exclusion needed -- removed from here accordingly.
    let remainder_wide_domain = |x: f32, y: f32| y != 0.0 && (x / y).abs() < 16777216.0;
    ok &= check2(
        "remainder_checked / remainder_wide",
        N,
        remainder_wide_domain,
        remainder_checked,
        remainder_wide,
    );

    if !ok {
        std::process::exit(1);
    }
}
