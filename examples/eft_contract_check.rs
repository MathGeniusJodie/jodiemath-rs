// Standing verification of the exactness contracts documented on the
// public EFT toolkit (idea #184): two_sum, quick_two_sum, two_prod,
// mulsign. Ground truth is f64, which represents every f32 sum and
// product exactly (24-bit operands -> at most a 48-bit product, well
// inside f64's 53), so "exact" here is a real check, not a tolerance.
//
// These four were internal helpers used only on bounded, well-behaved
// values; going public means the contracts have to hold for arbitrary
// user input, which is strictly more than the internal callers ever
// exercised. This is what caught two_prod's real range precondition.
//
// Exits nonzero if any documented contract is violated.

use jodiemath_rs::{mulsign, quick_two_sum, two_prod, two_sum};

// two_prod's documented safe band: exact iff 2^-102 <= |a*b| <= f32::MAX.
fn two_prod_lo() -> f64 {
    (-102f64).exp2()
}

struct Rng(u64);
impl Rng {
    fn next_u32(&mut self) -> u32 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        (x >> 32) as u32
    }
    // Uniform over every f32 bit pattern: all magnitudes, denormals and
    // specials included. Deliberately harsher than any real workload.
    fn wild(&mut self) -> f32 {
        f32::from_bits(self.next_u32())
    }
    // Uniform over finite normals with exponents in [-60, 60] -- the
    // regime the crate's own internal callers actually live in.
    fn tame(&mut self) -> f32 {
        let m = self.next_u32() & 0x7f_ffff;
        let e = ((self.next_u32() % 121) as i32 - 60 + 127) as u32;
        let s = self.next_u32() & 1;
        f32::from_bits((s << 31) | (e << 23) | m)
    }
}

fn main() {
    let lo = two_prod_lo();
    let hi = f32::MAX as f64;
    let mut rng = Rng(0x243f_6a88_85a3_08d3);
    let mut failures = 0u64;

    // ---- two_sum / quick_two_sum / two_prod over random pairs ----
    let mut ts_fail = 0u64;
    let mut ts_first = (0f32, 0f32);
    let mut qts_fail = 0u64;
    let mut qts_first = (0f32, 0f32);
    let mut tp_fail = 0u64;
    let mut tp_first = (0f32, 0f32);
    let mut tp_inrange = 0u64;

    const N: u64 = 40_000_000;
    for i in 0..N {
        let (a, b) = if i % 2 == 0 {
            (rng.tame(), rng.tame())
        } else {
            (rng.wild(), rng.wild())
        };
        if !a.is_finite() || !b.is_finite() {
            continue;
        }

        // two_sum: exact for any finite a, b whose sum does not overflow.
        let (s, e) = two_sum(a, b);
        if s.is_finite() && (s as f64 + e as f64) != (a as f64 + b as f64) {
            if ts_fail == 0 {
                ts_first = (a, b);
            }
            ts_fail += 1;
        }

        // quick_two_sum: exact only where |a| >= |b|, so order first.
        let (qa, qb) = if a.abs() >= b.abs() { (a, b) } else { (b, a) };
        let (qs, qe) = quick_two_sum(qa, qb);
        if qs.is_finite() && (qs as f64 + qe as f64) != (qa as f64 + qb as f64) {
            if qts_fail == 0 {
                qts_first = (qa, qb);
            }
            qts_fail += 1;
        }

        // two_prod: exact within the documented band only.
        let exact = a as f64 * b as f64;
        let m = exact.abs();
        if m >= lo && m <= hi {
            tp_inrange += 1;
            let (p, pe) = two_prod(a, b);
            if !p.is_finite() || (p as f64 + pe as f64) != exact {
                if tp_fail == 0 {
                    tp_first = (a, b);
                }
                tp_fail += 1;
            }
        }
    }

    println!("=== random pairs ({N} drawn, mixed tame/wild) ===");
    println!("two_sum        inexact (finite sums): {ts_fail}  first {ts_first:?}");
    println!("quick_two_sum  inexact (|a|>=|b|):    {qts_fail}  first {qts_first:?}");
    println!("two_prod       inexact in-band:       {tp_fail} of {tp_inrange}  first {tp_first:?}");
    failures += ts_fail + qts_fail + tp_fail;

    // ---- two_prod: confirm the band edges really are where they are ----
    // Below 2^-103 the error term is unrepresentable, so failures are the
    // rule, not the exception. This asserts the doc's claim that the
    // precondition is real rather than defensive boilerplate.
    let mut below = 0u64;
    let mut below_fail = 0u64;
    for _ in 0..20_000_000 {
        let a = rng.wild();
        let b = rng.wild();
        if !a.is_finite() || !b.is_finite() || a == 0.0 || b == 0.0 {
            continue;
        }
        let exact = a as f64 * b as f64;
        let m = exact.abs();
        if m == 0.0 || m >= (-103f64).exp2() {
            continue;
        }
        below += 1;
        let (p, pe) = two_prod(a, b);
        if !p.is_finite() || (p as f64 + pe as f64) != exact {
            below_fail += 1;
        }
    }
    println!("\n=== two_prod below the band (2^-103) ===");
    println!("  {below_fail} of {below} inexact -- the precondition is load-bearing");
    if below > 1000 && below_fail == 0 {
        println!("  UNEXPECTED: no failures below the band; the documented bound may be too conservative");
        failures += 1;
    }

    // ---- overflow behaviour, as documented ----
    println!("\n=== overflow (documented: two_sum e=NaN, two_prod p,e opposite infinities) ===");
    let (os, oe) = two_sum(f32::MAX, f32::MAX);
    let (op, ope) = two_prod(f32::MAX, f32::MAX);
    println!("  two_sum(MAX,MAX)  = ({os:e}, {oe:e})");
    println!("  two_prod(MAX,MAX) = ({op:e}, {ope:e})");
    if !(os.is_infinite() && oe.is_nan()) {
        println!("  MISMATCH: two_sum overflow no longer matches its doc");
        failures += 1;
    }
    if !(op.is_infinite() && ope.is_infinite() && op.signum() != ope.signum()) {
        println!("  MISMATCH: two_prod overflow no longer matches its doc");
        failures += 1;
    }

    // ---- two_sum exactness through the subnormal range ----
    let mut sub_fail = 0u64;
    for i in 1..3_000_000u32 {
        let a = f32::from_bits(i);
        let b = f32::from_bits((i.wrapping_mul(2_654_435_761) & 0x007f_ffff) | 1);
        let (s, e) = two_sum(a, b);
        if (s as f64 + e as f64) != (a as f64 + b as f64) {
            sub_fail += 1;
        }
    }
    println!("\n=== two_sum in the subnormal band ===");
    println!("  inexact: {sub_fail}");
    failures += sub_fail;

    // ---- mulsign: x * sign(y), exhaustively ----
    // The documented contract is an xor of sign bits, so -0.0 counts as a
    // negative y and the flip is unconditional on y's sign bit alone.
    let mut ms_fail = 0u64;
    for &y in &[1.0f32, -1.0, 0.0, -0.0, f32::NAN, -f32::NAN, f32::MAX, f32::MIN] {
        let y_neg = y.to_bits() & 0x8000_0000 != 0;
        for bits in 0..=u32::MAX {
            let want = if y_neg { bits ^ 0x8000_0000 } else { bits };
            if mulsign(f32::from_bits(bits), y).to_bits() != want {
                ms_fail += 1;
            }
        }
    }
    println!("\n=== mulsign (exhaustive: 8 y values x all 2^32 x) ===");
    println!("  mismatches vs xor-of-sign-bits: {ms_fail}");
    failures += ms_fail;

    // The documented divergence from copysign must actually hold, or the
    // whole reason this function exists alongside copysign is wrong.
    let mut agree_when_should_differ = 0u64;
    for &x in &[-2.0f32, -1.0, -f32::MAX, -f32::MIN_POSITIVE] {
        for &y in &[-3.0f32, -1.0, -f32::MAX] {
            if mulsign(x, y).to_bits() == x.copysign(y).to_bits() {
                agree_when_should_differ += 1;
            }
        }
    }
    println!("  x<0,y<0 cases where mulsign==copysign (must be 0): {agree_when_should_differ}");
    failures += agree_when_should_differ;

    println!("\n{}", if failures == 0 { "ALL CONTRACTS HOLD" } else { "CONTRACT VIOLATIONS FOUND" });
    if failures != 0 {
        std::process::exit(1);
    }
}
