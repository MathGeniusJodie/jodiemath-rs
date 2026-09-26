// Hardware FMA is required; software fallback changes rounding and performance.
#[cfg(all(not(target_feature = "fma"), not(doctest)))]
compile_error!(
    "jodiemath-rs requires hardware FMA (target-feature=+fma or target-cpu=native) -- \
     without it, f32::mul_add falls back to a slower, differently-rounded software path \
     and every accuracy/perf figure in this crate's docs is invalid. Build with \
     `RUSTFLAGS=\"-C target-cpu=native\"` or ensure .cargo/config.toml's rustflags \
     aren't being overridden by an environment RUSTFLAGS variable."
);

mod doublefloat;
mod pitable;
use doublefloat::Df32;

const SIGN_MASK: u32 = 0x80000000;
const EXPONENT_MASK: u32 = 0x7f800000;

#[inline(always)]
fn fma(a: f32, b: f32, c: f32) -> f32 {
    a.mul_add(b, c)
}

// 2^k as bare exponent field for integer k in [-127, 128]. k=128 yields +inf, k=-127 yields +0.0.
const EXP2INT_MAGIC: f32 = 12583039.0; // 1.5 * 2^23 + 127
macro_rules! exp2int_field {
    ($k:expr) => {
        f32::from_bits(($k + EXP2INT_MAGIC).to_bits() << 23)
    };
}

/// Round to the nearest integer (ties-to-even) for `|x| <= 2^22`.
/// Uses the magic constant `1.5 * 2^23`. Preserves the sign of zero for `x` in `[-0.5, 0)`.
#[inline(always)]
pub fn fast_round_int(x: f32) -> f32 {
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    (fma(x, 1.0, ROUND_MAGIC) - ROUND_MAGIC).copysign(x)
}

macro_rules! denormal_rescale {
    ($x:expr) => {{
        let tiny = $x < f32::MIN_POSITIVE;
        let xs = if tiny { $x * 16777216.0 } else { $x };
        let koff = if tiny { -24.0 } else { 0.0 };
        (xs, koff)
    }};
}

// exp(r) for tiny r.
macro_rules! exp_r_poly {
    ($r:expr) => {{
        let c: [f32; 4] = [4.999_93e-1, 1.6667245e-1, 4.188_381e-2, 8.300_99e-3];
        let r2 = $r * $r;
        let l0 = $r + 1.0;
        let l1 = fma(c[1], $r, c[0]);
        let l2 = fma(c[3], $r, c[2]);
        let m = fma(l2, r2, l1);
        fma(m, r2, l0)
    }};
}

// Cody-Waite reduction and poly; caller must ensure x is in EXP_CLAMP_LO..=EXP_CLAMP_HI.
macro_rules! exp_reduce {
    ($x:expr) => {{
        let x = $x;
        const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
        let k = fma(x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
        let r = fma(-k, LN2_HI, x);
        let r = fma(-k, LN2_LO, r);
        let c: [f32; 5] = [
            0.50000006,
            0.16666451,
            0.041665636,
            0.008_374_871,
            0.0013946877,
        ];
        let r2 = r * r;
        let l1 = fma(c[1], r, c[0]);
        let l2 = fma(c[3], r, c[2]);
        let l3 = fma(c[4], r2, l2);
        let m = fma(l3, r2, l1);
        let s = fma(m, r2, r);
        let (t1, t2) = exp2_field_split(k);
        fma(s, t1, t1) * t2
    }};
}

/// Clamp bounds for `exp_checked` in `x` units.
const EXP_CLAMP_LO: f32 = -104.665_22;
const EXP_CLAMP_HI: f32 = 88.722_84;

// Q(f) = (2^f - 1)/f for f in [0, 1).
macro_rules! exp2_q_poly {
    ($f:expr) => {{
        let f2 = $f * $f;
        let g0 = fma(2.4022985e-1, $f, 6.93147e-1);
        let g1 = fma(9.678817e-3, $f, 5.548333e-2);
        let g2 = fma(2.1702255e-4, $f, 1.2439643e-3);
        let h = fma(g2, f2, g1);
        fma(h, f2, g0)
    }};
}

// Q(f) = (2^f - 1)/f refit for f in [-0.5, 0.5].
macro_rules! exp2_q_poly_centered {
    ($f:expr) => {{
        let f2 = $f * $f;
        let g0 = fma(2.402265e-1, $f, 6.931_472e-1);
        let g1 = fma(9.618_238e-3, $f, 5.5503574e-2);
        let g2 = fma(1.5403504e-4, $f, 1.3390731e-3);
        let h = fma(g2, f2, g1);
        fma(h, f2, g0)
    }};
}

macro_rules! exp_pos_neg_core {
    ($x:expr) => {{
        const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
        let k = fma($x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
        let r = fma(-k, LN2_HI, $x);
        let r = fma(-k, LN2_LO, r);
        let c: [f32; 5] = [
            0.49999994 * 0.5,
            0.16666521 * 0.5,
            0.041_668_33 * 0.5,
            8.3687045e-3 * 0.5,
            1.3814511e-3 * 0.5,
        ];
        let r2 = r * r;
        let e = fma(fma(fma(c[4], r2, c[2]), r2, c[0]), r2, 0.5);
        let o = fma(fma(c[3], r2, c[1]), r2, 0.5);
        let p_pos = fma(r, o, e);
        let p_neg = fma(-r, o, e);
        let (t1, t2) = exp2_field_split(k);
        // Exact reciprocal for power-of-two float: 2^-e bits are 0x7F000000 - b.
        let t1n = f32::from_bits(0x7F00_0000u32.wrapping_sub(t1.to_bits()));
        let t2n = f32::from_bits(0x7F00_0000u32.wrapping_sub(t2.to_bits()));
        (p_pos, p_neg, t1, t2, t1n, t2n)
    }};
}

macro_rules! log_family_edges {
    ($x:expr, $r:expr) => {{
        let r = $r;
        let spec = if $x == 0.0 {
            f32::NEG_INFINITY
        } else {
            f32::NAN
        };
        let r = if $x <= 0.0 { spec } else { r };
        if !($x < f32::INFINITY) {
            $x * $x
        } else {
            r
        }
    }};
}

macro_rules! log_family_wrapper {
    ($x:expr, $normal:ident) => {
        log_family_edges!($x, {
            let (xs, koff) = denormal_rescale!($x);
            $normal(xs, koff)
        })
    };
}

// log wrapper for arguments known not to be positive denormals.
macro_rules! log_family_wrapper_no_denormal {
    ($x:expr, $normal:ident) => {
        log_family_edges!($x, $normal($x, 0.0))
    };
}

// log wrapper when caller discards the result unless x is positive normal.
macro_rules! log_family_wrapper_discarded_unless_normal {
    ($x:expr, $normal:ident) => {
        if !($x < f32::INFINITY) {
            $x * $x
        } else {
            $normal($x, 0.0)
        }
    };
}

#[doc(alias = "log2f")]
#[doc(alias = "log2")]
#[inline(always)]
#[allow(clippy::neg_cmp_op_on_partial_ord)]
pub fn log_2(x: f32) -> f32 {
    log_family_wrapper!(x, log_2_normal)
}

/// Core log_2 for positive normal finite x only (caller handles edge cases).
const LOG2_Q_COEFFS: [f32; 9] = [
    -0.7213475,
    0.48089963,
    -0.36067435,
    0.28850868,
    -0.24009936,
    0.2058956,
    -0.18871288,
    0.17711402,
    -0.10358754,
];

#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn log_2_normal(x: f32, koff: f32) -> f32 {
    // Decompose x = 2^k * m with m in [sqrt(2)/2, sqrt(2)), centered at 1.
    let e = (x.to_bits() as i32).wrapping_sub(0x3f3504f3) >> 23;
    let m = f32::from_bits((x.to_bits() as i32).wrapping_sub(e << 23) as u32);
    let k = e as f32 + koff;
    let s = m - 1.0;
    let c = LOG2_Q_COEFFS;
    let s2 = s * s;
    let s4 = s2 * s2;
    let l0 = fma(c[1], s, c[0]);
    let l1 = fma(c[3], s, c[2]);
    let l2 = fma(c[5], s, c[4]);
    let l3 = fma(c[7], s, c[6]);
    let a = s2 * l0;
    let w0 = fma(l2, s2, l1);
    let w1 = fma(c[8], s2, l3);
    let v = fma(w1, s4, w0);
    let sq = fma(v, s4, a);
    // k joins last to avoid double-rounding at k's scale.
    let lm = fma(s, std::f32::consts::LOG2_E, sq);
    lm + k
}

/// log_2 without domain checks: valid for positive normal finite x only.
#[inline(always)]
pub fn log_2_unchecked(x: f32) -> f32 {
    log_2_normal(x, 0.0)
}

/// exp2 without domain checks: valid for x in [-126, 128), i.e. normal
/// (non-denormal, finite, nonzero) results only.
#[doc(alias = "exp2f")]
#[inline(always)]
#[allow(clippy::approx_constant)] // g0's constant term is a fitted minimax
pub fn exp2(x: f32) -> f32 {
    // exp2int must come from the same floor(x) as f to avoid double-counting.
    let k = x.floor();
    let f = x - k;
    let exp2int = exp2int_field!(k);
    let q = exp2_q_poly!(f);
    fma(q, exp2int * f, exp2int)
}

/// Computes `2^(k+f)` for integer-valued `k` in `[-126, 128)` and `f` in `[0, 1)`.
#[inline(always)]
#[allow(clippy::approx_constant)]
pub fn exp2_kf(k: f32, f: f32) -> f32 {
    let exp2int = exp2int_field!(k);
    let q = exp2_q_poly!(f);
    fma(q, exp2int * f, exp2int)
}

#[inline(always)]
#[allow(clippy::approx_constant)] // g0's constant term is a fitted minimax
pub fn exp2_checked(x: f32) -> f32 {
    // Splitting k = k1 + k2 keeps both powers of 2 representable across clamp range.
    let xs = x.clamp(-151.0, 128.0);
    let k = xs.floor();
    let f = xs - k;
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k1b = fma(xs, 0.5, ROUND_MAGIC) - (ROUND_MAGIC - 383.0); // k1 + 383
    let k2b = (k + 766.0) - k1b; // k2 + 383, exact: all integers
    let t1 = f32::from_bits((k1b.to_bits() << 8) & EXPONENT_MASK);
    let t2 = f32::from_bits((k2b.to_bits() << 8) & EXPONENT_MASK);
    let q = exp2_q_poly!(f);
    let p = fma(q, t1 * f, t1);
    p * t2
}

/// Computes `x * 2^n` (`ldexp`/`scalbn`).
/// Decomposes `x` via [`frexp`], adds `n` to the exponent, and reconstructs the float.
/// Overflows to `+/-inf`, underflows to `+/-0.0`.
#[inline(always)]
pub fn ldexp(x: f32, n: i32) -> f32 {
    let (mantissa, e) = frexp(x);
    let target_exp_wide = e as i64 + n as i64;
    let overflow = target_exp_wide > 128;
    let underflow = target_exp_wide < -151;
    let target_exp = target_exp_wide.clamp(-151, 128) as f32;
    let (t1, t2) = exp2_field_split(target_exp);
    let reconstructed = mantissa * t1 * t2;
    let saturated = if overflow {
        f32::INFINITY.copysign(x)
    } else if underflow {
        0.0f32.copysign(x)
    } else {
        reconstructed
    };
    if x == 0.0 || !x.is_finite() {
        x
    } else {
        saturated
    }
}

/// Decomposes `x` into `(mantissa, exponent)` such that `x == mantissa * 2^exponent`,
/// with `mantissa` in `[0.5, 1)`. If `x` is zero or non-finite, returns `(x, 0)`.
#[inline(always)]
pub fn frexp(x: f32) -> (f32, i32) {
    let ax = x.abs();
    let (xs, koff) = denormal_rescale!(ax);
    let bits = xs.to_bits();
    let raw_exp = (bits >> 23) & 0xFF;
    let mantissa_bits = (bits & 0x807FFFFF) | (126u32 << 23);
    let mantissa = f32::from_bits(mantissa_bits).copysign(x);
    let exponent = raw_exp as i32 - 126 + koff as i32;
    let is_special = x == 0.0 || !x.is_finite();
    (
        if is_special { x } else { mantissa },
        if is_special { 0 } else { exponent },
    )
}

/// Cody-Waite reduction for `10^x`.
macro_rules! exp10_reduction {
    ($x:expr) => {{
        const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
        let kb = fma($x, std::f32::consts::LOG2_10, ROUND_MAGIC);
        let kr = kb - ROUND_MAGIC; // round(x*log2(10)), coarse multiply is fine
        let d = fma(-kr, LOG10_2_HI, $x);
        let d = fma(-kr, LOG10_2_LO, d);
        let fr = d * std::f32::consts::LOG2_10; // small, precise correction in log2 units, in [-0.5, 0.5]
        // Adjust from round's [-0.5, 0.5] to exp2's [0, 1) convention.
        let adjust = if fr < 0.0 { 1.0 } else { 0.0 };
        let k = kr - adjust;
        let f = fr + adjust;
        (k, f)
    }};
}

#[inline(always)]
#[allow(clippy::approx_constant)] // g0's constant term is a fitted minimax
pub fn exp10_checked(x: f32) -> f32 {
    // Clamp prevents +-inf from producing an inf-inf NaN in reduction.
    let x = x.clamp(-45.154503, 38.53184);
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let kb = fma(x, std::f32::consts::LOG2_10, ROUND_MAGIC);
    let k = kb - ROUND_MAGIC; // round(x*log2(10))
    let d = fma(-k, LOG10_2_HI, x);
    let d = fma(-k, LOG10_2_LO, d);
    let f = d * std::f32::consts::LOG2_10;
    let (t1, t2) = exp2_field_split(k);
    let q = exp2_q_poly_centered!(f);
    let p = fma(q, t1 * f, t1);
    p * t2
}

/// `exp10` via a single exponent field. Valid for `k` in `[-126, 128)`.
#[inline(always)]
#[allow(clippy::approx_constant)] // g0's constant term is a fitted minimax
pub fn exp10(x: f32) -> f32 {
    let (k, f) = exp10_reduction!(x);
    let exp2int = exp2int_field!(k);
    let q = exp2_q_poly!(f);
    fma(q, exp2int * f, exp2int)
}

// sin(x) ~= x + x^3*p(x^2) on [-pi/2, pi/2]; odd, and keeps the sign of zero.
#[inline(always)]
fn sinf_poly(x: f32) -> f32 {
    let c0 = -0.166_666_6_f32;
    let c1 = 8.3330662e-3f32;
    let c2 = -1.9809603e-4f32;
    let c3 = 2.6057806e-6f32;
    let y = x * x;
    let y2 = y * y;
    // `+ 0.0` turns x3 = -0 into +0, so p*x3 = -0 and x + p*x3 keeps x's zero.
    let x3 = fma(y, x, 0.0);
    let a = fma(c1, y, c0);
    let b = fma(c3, y, c2);
    let p = fma(b, y2, a);
    fma(p, x3, x)
}

// Cody-Waite pi split; trailing zeros keep leading reduction steps exact.
const PI_A: f32 = 3.140625;
const PI_B: f32 = 0.000_967_025_76;
const PI_C: f32 = 6.278_329_5e-7;
const PI_D: f32 = 1.078_060_6e-14;
const FRAC_1_PI: f32 = std::f32::consts::FRAC_1_PI;

// Magic constant for round-to-nearest-even on |v| < 2^22.
const ROUND_MAGIC: f32 = 12582912.0;

// Rounds |v| < 2^24 to multiples of 4.
const ROUND_MAGIC_4: f32 = 50331648.0;

/// Computes `sin(x)` via single-f32 Cody-Waite range reduction.
/// Accurate for `|x| < 2^24 * pi` (~5.27e7).
macro_rules! pi_reduce_and_poly {
    ($x:expr, $q:expr) => {{
        let r = fma($q, -PI_A, $x);
        let r = fma($q, -PI_B, r);
        let r = fma($q, -PI_C, r);
        let r = fma($q, -PI_D, r);
        sinf_poly(r)
    }};
}

macro_rules! frac_x_over_pi {
    ($x:expr, $magic:expr) => {{
        let nb = fma($x, FRAC_1_PI, $magic);
        let n = nb - $magic;
        let f = fma($x, FRAC_1_PI, -n);
        (nb, n, fma($x, RPI_LO, f))
    }};
}

/// Computes `sin(x)` (radians) via single-f32 Cody-Waite range reduction.
/// Accurate for `|x| < 2^24 * pi` (~5.27e7).
#[doc(alias = "sinf")]
#[inline(always)]
pub fn sin(x: f32) -> f32 {
    let (_, n, fc) = frac_x_over_pi!(x, ROUND_MAGIC_4);
    // Fine round of fc to complete round(x/pi) without losing integer precision.
    let nm = n - ROUND_MAGIC;
    let qb = fc + ROUND_MAGIC;
    let q = nm + qb;
    let s = pi_reduce_and_poly!(x, q);
    // q's parity is in qb's lowest mantissa bit since n is a multiple of 4.
    let parity = qb.to_bits() << 31;
    f32::from_bits(s.to_bits() ^ parity)
}
/// Computes `cos(x)` (radians) via single-f32 Cody-Waite range reduction.
/// Accurate for `|x| < 2^22 * pi` (~1.32e7).
#[doc(alias = "cosf")]
#[inline(always)]
pub fn cos(x: f32) -> f32 {
    // Nearest half-odd-integer is n + copysign(0.5, fc).
    let (nb, n, fc) = frac_x_over_pi!(x, ROUND_MAGIC);
    let q = n + 0.5f32.copysign(fc);
    let s = pi_reduce_and_poly!(x, q);
    let parity = (nb.to_bits() << 31) ^ (!fc.to_bits() & SIGN_MASK);
    f32::from_bits(s.to_bits() ^ parity)
}

/// Computes `sin(pi * x)`, argument in half-turns. Exact at integers and total over all finite f32.
#[inline(always)]
pub fn sinpi(x: f32) -> f32 {
    let q = x.round_ties_even();
    let r = x - q;
    // Preserve -0.0: subtraction of equal zeros yields +0.0 in IEEE 754.
    let normal = sinf_poly(std::f32::consts::PI * r) * fma(-2.0, parity(q), 1.0);
    if x == 0.0 {
        x
    } else {
        normal
    }
}

/// Computes `cos(pi * x)`, argument in half-turns. Exact at half-integers and total over all finite f32.
#[inline(always)]
pub fn cospi(x: f32) -> f32 {
    let k = x.round_ties_even();
    let r = x - k;
    let s = sinf_poly(std::f32::consts::PI * (0.5 - r.abs()));
    s * fma(-2.0, parity(k), 1.0)
}

/// Unnormalized sinc function: `sin(x) / x` in radians, with `sinc(0) = 1.0`.
#[inline(always)]
pub fn sinc_unnormalized(x: f32) -> f32 {
    let normal = sin_wide(x) / x;
    if x == 0.0 {
        1.0
    } else {
        normal
    }
}

// B(u) = tan(pi*w)/w - pi for w in [0, 0.25], u = w*w.
#[inline(always)]
fn tan_poly(u: f32) -> f32 {
    let c: [f32; 7] = [
        -8.742278e-8,
        10.335385,
        40.82169,
        160.9828,
        741.586_5,
        701.914_2,
        28496.229,
    ];
    let u2 = u * u;
    let u4 = u2 * u2;
    let l0 = fma(c[1], u, c[0]);
    let l1 = fma(c[3], u, c[2]);
    let l2 = fma(c[6], u2, fma(c[5], u, c[4]));
    fma(u4, l2, fma(u2, l1, l0))
}

#[inline(always)]
fn tan_core(r: f32, s: f32) -> f32 {
    let ar = r.abs();
    let direct = fma(r, std::f32::consts::PI, r * tan_poly(r * r));
    let reflected = mulsign(1.0 / fma(s, std::f32::consts::PI, s * tan_poly(s * s)), r);
    let normal = if ar <= 0.25 { direct } else { reflected };
    if s == 0.0 {
        f32::NEG_INFINITY
    } else {
        normal
    }
}

/// Computes `tan(pi * x)`, argument in half-turns. Total over all finite f32.
#[inline(always)]
pub fn tanpi(x: f32) -> f32 {
    let q = x.round_ties_even();
    let r = x - q;
    let ar = r.abs();
    let normal = tan_core(r, 0.5 - ar);
    if x == 0.0 {
        x
    } else {
        normal
    }
}

// q = round(x/pi) must be an exact integer for x - q*pi to land in [-pi/2, pi/2].
const RPI_LO: f32 = 1.284_127_65e-8;
// 1.5 * 2^52: rounds f64 to integer; parity lands in bit 0.
const ROUND_MAGIC64: f64 = 6755399441055744.0;

/// Error-free transformation (Knuth two-sum): returns `(s, e)` such that `s + e == a + b` exactly.
#[inline(always)]
pub fn two_sum(a: f32, b: f32) -> (f32, f32) {
    let s = a + b;
    let v = s - a;
    let e = (a - (s - v)) + (b - v);
    (s, e)
}

/// Error-free transformation (Fast2Sum): returns `(s, e)` such that `s + e == a + b` exactly.
/// Requires `|a| >= |b|`.
#[inline(always)]
pub fn quick_two_sum(a: f32, b: f32) -> (f32, f32) {
    let s = a + b;
    let e = b - (s - a);
    (s, e)
}

/// Error-free product: returns `(p, e)` such that `p + e == a * b` exactly.
#[inline(always)]
pub fn two_prod(a: f32, b: f32) -> (f32, f32) {
    let p = a * b;
    let e = fma(a, b, -p);
    (p, e)
}

// Parity of exact-integer float q without saturating float-to-int cast.
#[inline(always)]
fn parity(q: f32) -> f32 {
    fma(-2.0, (q * 0.5).floor(), q)
}

const PI_UNIT: f32 = std::f32::consts::PI * (1.0 / 2147483648.0);
const PI_UNIT_LO: f32 = -8.742278e-8 * (1.0 / 2147483648.0);

/// Payne-Hanek reduction of `|x|` modulo `2*pi`, in 32-bit lanes only.
///
/// `|x|/pi + GRID/2^32 = (h - 2^30 - 64 + rem) * 2^-31 + tail/pi (mod 2)`:
/// `h` is exact fixed point whose wraparound is the mod 2, `|rem| <= 1`,
/// `tail` is tiny and already in radians. `GRID = 0` is the sin grid, `2^30`
/// (a quarter period) the cos grid. Meaningless where [`off_window`].
#[inline(always)]
fn reduce_pi_wide<const GRID: u32>(b: u32) -> (u32, f32, f32) {
    // ulp 2 on [2^24, 2^25): m*g1 < 2^24 always rounds inside one binade.
    const M: f32 = 16777216.0;
    const K: [u32; 7] = pitable::WORDS;
    // The 96 bits of beta(e) = (2^(e-150)/pi) mod 2 that matter start at bit
    // t = e - CUT of K: words q..q+3 (q = t >> 5 <= 3) funnel-shifted by
    // t & 31. The words are picked by a select tree over immediates, never an
    // indexed load: that vectorizes to gathers (or scalar loads on targets
    // without them), which cost 3/4 of the old table version's time. Plain
    // `if`s get folded back into a lookup table, hence select_unpredictable.
    // Only bits 0..6 of t are read, so the sign bit above e is harmless.
    let t = (b >> 23).wrapping_sub(pitable::CUT);
    let (b0, b1) = (t & 32 != 0, t & 64 != 0);
    let sel = core::hint::select_unpredictable::<u32>;
    let u = |i: usize| sel(b0, K[i + 1], K[i]);
    let (u0, u1, u2, u3, u4, u5) = (u(0), u(1), u(2), u(3), u(4), u(5));
    let s = t & 31;
    let funnel = |hi: u32, lo: u32| (hi << s) | ((lo >> 1) >> (31 - s));
    let (w0, w1) = (sel(b1, u2, u0), sel(b1, u3, u1));
    let (w2, w3) = (sel(b1, u4, u2), sel(b1, u5, u3));
    let (x0, x1, x2) = (funnel(w0, w1), funnel(w1, w2), funnel(w2, w3));
    // x0 = floor(beta * 2^31) mod 2^32; next 24 bits exact in n1, 24 more in n2.
    let n1 = (x1 >> 8) as f32;
    let n2 = (((x1 << 16) & 0x00ff_0000) | (x2 >> 16)) as f32;
    let mf = f32::from_bits((b & 0x007f_ffff) | 0x3f00_0000); // m * 2^-24
    let m = (b & 0x007f_ffff) | 0x0080_0000;
    let kb = fma(mf, n1, M);
    let rem = fma(mf, n1, M - kb);
    let tail = n2 * (mf * (PI_UNIT * (1.0 / 16777216.0)));
    let h = m
        .wrapping_mul(x0)
        .wrapping_add(
            GRID.wrapping_add((1 << 30) + 64)
                .wrapping_sub(M.to_bits() << 1),
        )
        .wrapping_add(kb.to_bits() << 1);
    (h, rem, tail)
}

/// Lanes `reduce_pi_wide` does not cover: `|x| < 1`, inf and NaN.
#[inline(always)]
fn off_window(b: u32) -> bool {
    // e - CUT lands in [0, 128) exactly for the windowed exponents, for
    // either sign (the sign bit adds 256).
    (b >> 23).wrapping_sub(pitable::CUT) & 128 != 0
}

/// `pi` times the signed distance from a reduction word to the nearest
/// integer, and that integer's parity as a sign mask.
#[inline(always)]
fn reduced_angle(h: u32, rem: f32, tail: f32) -> (f32, u32) {
    (angle(hi_part(h), h, rem, tail), h & SIGN_MASK)
}

/// The centred part of a reduction word that is a multiple of 128, so its
/// conversion to f32 is exact.
#[inline(always)]
fn hi_part(h: u32) -> i32 {
    (h & 0x7fff_ff80) as i32 - (1 << 30)
}

/// `pi * (hi + lo + rem) * 2^-31 + tail`, where `lo` is the rest of `h`,
/// rounded once. `lo + rem` is exact exactly when it cancels.
#[inline(always)]
fn angle(hi: i32, h: u32, rem: f32, tail: f32) -> f32 {
    let lo = (h & 0x7f) as i32 - 64;
    let f = hi as f32;
    let d = fma(lo as f32 + rem, PI_UNIT, fma(f, PI_UNIT_LO, tail));
    fma(f, PI_UNIT, d)
}

/// `clamp(-1, 1)` that keeps NaN.
#[inline(always)]
fn clamp_unit(s: f32) -> f32 {
    let s = if s < -1.0 { -1.0 } else { s };
    if s > 1.0 {
        1.0
    } else {
        s
    }
}

/// Computes `sin(x)` with no magnitude limit across all finite f32.
#[inline(always)]
pub fn sin_wide(x: f32) -> f32 {
    let b = x.to_bits();
    let (h, rem, tail) = reduce_pi_wide::<0>(b);
    let (r, parity) = reduced_angle(h, rem, tail);
    let r = f32::from_bits(r.to_bits() ^ parity ^ (b & SIGN_MASK));
    // x - (x - x) is x, or NaN for inf.
    let r = if off_window(b) { x - (x - x) } else { r };
    clamp_unit(sinf_poly(r))
}

/// Computes `cos(x)` with no magnitude limit across all finite f32.
#[inline(always)]
pub fn cos_wide(x: f32) -> f32 {
    let b = x.to_bits();
    let (h, rem, tail) = reduce_pi_wide::<{ 1 << 30 }>(b);
    let (r, parity) = reduced_angle(h, rem, tail);
    let c = clamp_unit(sinf_poly(f32::from_bits(r.to_bits() ^ parity)));
    // Near 1 the odd poly at pi/2 - |x| is noisy by half an ulp; 1 + y*C(y)
    // is not. x - (x - x) is x, or NaN for inf.
    if off_window(b) {
        cos_unit(x - (x - x))
    } else {
        c
    }
}

/// Largest magnitude [`wrap_pi`] can return: the largest `f32` whose *exact*
/// value is below `pi`, one ulp under `f32::consts::PI`.
pub const WRAP_PI_MAX: f32 = f32::from_bits(0x40490fda);

/// `(sin r, cos r)` for `|r| <= pi/4`: relative error 2^-28 and 2^-33
/// (weighted Remez, `tools/remez_trig.py`). Keeps the sign of zero like
/// `sinf_poly`.
#[inline(always)]
fn sincos_quarter(r: f32) -> (f32, f32) {
    let s: [f32; 3] = [-0.16666655, 0.00833216, -0.00019515282];
    let c: [f32; 4] = [-0.5, 0.04166662, -0.0013886682, 2.4383566e-5];
    let y = r * r;
    let y2 = y * y;
    let sp = fma(s[2], y2, fma(s[1], y, s[0]));
    let cp = fma(fma(c[3], y, c[2]), y2, fma(c[1], y, c[0]));
    (fma(sp, fma(y, r, 0.0), r), fma(cp, y, 1.0))
}

/// `cos r` for `|r| <= 1`: relative error 2^-30 (weighted Remez,
/// `tools/remez_trig.py`).
#[inline(always)]
fn cos_unit(r: f32) -> f32 {
    let c: [f32; 4] = [-0.49999997, 0.04166646, -0.0013882959, 2.4118643e-5];
    let y = r * r;
    let cp = fma(fma(c[3], y, c[2]), y * y, fma(c[1], y, c[0]));
    fma(cp, y, 1.0)
}

/// `(sin r, cos r)` for `|r| <= 1`: relative error 2^-34 and 2^-30.
/// Keeps the sign of zero.
#[inline(always)]
fn sincos_unit(r: f32) -> (f32, f32) {
    let s: [f32; 4] = [-0.16666667, 0.008333316, -0.00019836116, 2.69481e-6];
    let y = r * r;
    let y2 = y * y;
    let sp = fma(fma(s[3], y, s[2]), y2, fma(s[1], y, s[0]));
    (fma(sp, fma(y, r, 0.0), r), cos_unit(r))
}

/// Computes `tan(x)` with no magnitude limit across all finite f32.
#[inline(always)]
pub fn tan_wide(x: f32) -> f32 {
    let b = x.to_bits();
    let (h, rem, tail) = reduce_pi_wide::<0>(b);
    let off = off_window(b);
    // Beyond a quarter turn, reduce onto the cos grid instead and use
    // tan(t) = -cos(t - pi/2) / sin(t - pi/2): |r| stays within pi/4.
    // Off the window |x| < 1 already, and r = x.
    let near = (h ^ (h << 1)) & (1 << 30) != 0 || off;
    // Selects, not branches: `near` is a coin flip for arbitrary inputs.
    let sel = core::hint::select_unpredictable::<u32>;
    let hi = sel(
        near,
        hi_part(h) as u32,
        hi_part(h.wrapping_add(1 << 30)) as u32,
    );
    let r = angle(hi as i32, h, rem, tail);
    let r = f32::from_bits(r.to_bits() ^ (b & SIGN_MASK) ^ sel(near, 0, SIGN_MASK));
    let r = if off { x - (x - x) } else { r };
    let (s, c) = sincos_unit(r);
    let (n, d) = core::hint::select_unpredictable(near, (s, c), (c, s));
    n / d
}

/// Core of cbrt for normal finite x: bit-trick seed followed by degree-3 polynomial correction.
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
fn cbrt_normal(x: f32) -> f32 {
    let ax = x.to_bits() & !SIGN_MASK;
    let a = f32::from_bits(ax);
    let rcp = 1.0 / a; // independent of the seed chain, starts immediately
    let s = f32::from_bits(ax / 3 + 0x2a509a07u32);
    let s2 = s * s;
    let d = fma(s2, s, -a);
    let r = d * rcp;
    let c1 = -0.333_331_4_f32;
    let c2 = 0.222_213_36_f32;
    let c3 = -0.173_940_27_f32;
    let c4 = 0.147_204_53_f32;
    let r2 = r * r;
    let a1 = fma(c2, r, c1);
    let b1 = fma(c4, r, c3);
    let p = fma(b1, r2, a1);
    let ss = f32::from_bits(s.to_bits() | (x.to_bits() & SIGN_MASK));
    let sr = ss * r;
    fma(sr, p, ss)
}

#[doc(alias = "cbrtf")]
#[inline(always)]
pub fn cbrt(x: f32) -> f32 {
    let ax = x.to_bits() & !SIGN_MASK;
    let tiny = ax < 0x0080_0000; // denormal or zero: rescale by 2^24 = (2^8)^3
    let xs = if tiny { x * 16777216.0 } else { x };
    let scale = if tiny { 0.00390625 } else { 1.0 };
    let r = cbrt_normal(xs) * scale;
    if ax == 0 || ax >= EXPONENT_MASK {
        x + x
    } else {
        r
    }
}

/// `cbrt` without domain checks: valid for normal finite `x`.
#[inline(always)]
pub fn cbrt_unchecked(x: f32) -> f32 {
    cbrt_normal(x)
}

/// cbrt to ~0.5 ulp. Caller must ensure x is in roughly [2^-56, 2^127] to prevent underflow/overflow.
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn cbrt_accurate_normal(x: f32, scale: f32) -> f32 {
    let ax = x.to_bits() & !SIGN_MASK;
    let a = f32::from_bits(ax);
    let neg_rcp3_scale = -((1.0 / a) * (1.0 / 3.0)) * scale;
    let y = cbrt_normal(x);
    let y2 = Df32::from_mul(y, y);
    let y3 = y2 * y;
    let e = (y3.0 - x) + y3.1;
    // Approximates 1/(3y^2) as |y|/(3a) using y^3 ~ x, avoiding a division.
    let neg_den_recip_scale = y.abs() * neg_rcp3_scale;
    fma(e, neg_den_recip_scale, y * scale)
}

#[inline(always)]
pub fn cbrt_accurate(x: f32) -> f32 {
    let ax = x.to_bits() & !SIGN_MASK;
    // Rescale to prevent intermediate overflow (above 2^127) or underflow (below 2^-56).
    const SCALE_UP: f32 = f32::from_bits(0x7e80_0000); // 2^126
    const SCALE_UP_OUT: f32 = f32::from_bits(0x2a80_0000); // 2^-42
    const SCALE_DN: f32 = f32::from_bits(0x0080_0000); // 2^-126
    const SCALE_DN_OUT: f32 = f32::from_bits(0x5480_0000); // 2^42
    let small = ax < 0x2380_0000;
    let big = ax >= 0x7f00_0000; // 2^127; inf/nan land here too, fixed up below
    let xs = if small {
        x * SCALE_UP
    } else if big {
        x * SCALE_DN
    } else {
        x
    };
    let scale = if small {
        SCALE_UP_OUT
    } else if big {
        SCALE_DN_OUT
    } else {
        1.0
    };
    let r = cbrt_accurate_normal(xs, scale);
    if ax == 0 || ax >= EXPONENT_MASK {
        x + x
    } else {
        r
    }
}

/// `cbrt_accurate` without domain checks: valid for normal finite `x`.
#[inline(always)]
pub fn cbrt_accurate_unchecked(x: f32) -> f32 {
    cbrt_accurate_normal(x, 1.0)
}

/// Core of rcbrt for normal finite x: bit-trick seed and degree-3 correction.
#[inline(always)]
fn rcbrt_normal(x: f32) -> f32 {
    let ax = x.to_bits() & !SIGN_MASK;
    let a = f32::from_bits(ax);
    let t = f32::from_bits(0x54a20d0eu32.wrapping_sub(ax / 3));
    let t2 = t * t;
    let e = fma(a * t2, t, -1.0);
    let d1 = -0.333_333_3_f32;
    let d2 = 0.222_217_9_f32;
    let d3 = -0.172_844_5_f32;
    let d4 = 0.145_519_35_f32;
    let d5 = -0.125_217_36_f32;
    let e2 = e * e;
    let a1 = fma(d2, e, d1);
    let b1 = fma(d4, e, d3);
    let p = fma(fma(d5, e2, b1), e2, a1);
    let ts = f32::from_bits(t.to_bits() | (x.to_bits() & SIGN_MASK));
    fma(ts * e, p, ts)
}

/// Computes `x^(-1/3)` (reciprocal cube root).
#[inline(always)]
pub fn rcbrt(x: f32) -> f32 {
    let ax = x.to_bits() & !SIGN_MASK;
    let tiny = ax < 0x0080_0000; // denormal or zero
    let xs = if tiny { x * 16777216.0 } else { x };
    let scale = if tiny { 256.0 } else { 1.0 };
    let r = rcbrt_normal(xs) * scale;
    let spec = f32::from_bits(EXPONENT_MASK.wrapping_sub(ax) | (x.to_bits() & SIGN_MASK));
    if ax == 0 || ax >= EXPONENT_MASK {
        spec
    } else {
        r
    }
}

/// Bit-trick approximation of `cbrt(x)` with two rational refinement steps.
pub fn cbrt_approx(x: f32) -> f32 {
    let y = f32::from_bits(0x2a509849u32 + (x.to_bits() / 3));
    let y = (x + 2. * (y * y) * y) / (3. * (y * y));
    (2. * x * y + (y * y) * (y * y)) / (x + 2. * (y * y) * y)
}
/// Single-bit-trick seed for `sqrt(x)`.
pub fn sqrt_approx(x: f32) -> f32 {
    f32::from_bits(0x1FBD22DF + (x.to_bits() >> 1))
}
/// Single-bit-trick seed for `1/x`.
pub fn rcp_approx(x: f32) -> f32 {
    f32::from_bits(0x7EEF370B - x.to_bits())
}
/// Single-bit-trick seed for `2^x`.
pub fn exp2_approx(x: f32) -> f32 {
    -f32::from_bits((x + 383.).to_bits() << 8)
}

/// Single-bit-trick seed for `log2(x)`.
pub fn log2_approx(x: f32) -> f32 {
    f32::from_bits((x).to_bits() >> 8 | 256_f32.to_bits()) - 383.
}

/// Quake-style bit-trick seed for `1/sqrt(x)`.
pub fn rsqrt_approx(x: f32) -> f32 {
    f32::from_bits(0x5F33E79F - (x.to_bits() >> 1))
}

/// Latency-optimal `cbrt` approximation using two bit-trick seeds.
#[inline(always)]
pub fn cbrt_fast(x: f32) -> f32 {
    let s = f32::from_bits(0x2a4d_def1u32.wrapping_add((x.to_bits() >> 16) * 0x5556u32));
    let r = f32::from_bits(0x68ff_2381u32.wrapping_sub((x.to_bits() >> 16) * 0xaaac));
    let s = fma(s * s, s * -r, fma(r, x, s));
    fma(s * s, s * -r, fma(r, x, s))
}

/// Computes `x * sign(y)` via an XOR of sign bits (preserves magnitude of `x`).
#[inline(always)]
pub fn mulsign(x: f32, y: f32) -> f32 {
    f32::from_bits(x.to_bits() ^ (y.to_bits() & SIGN_MASK))
}

const LN_2: f32 = std::f32::consts::LN_2;
const LOG2_E: f32 = std::f32::consts::LOG2_E;

// Low words of two-word log2(e) and log10(e) for tiny |x| < 2^-24.
const LOG2_E_LO: f32 = f32::from_bits(0x32a5_7060);
const LOG10_E_LO: f32 = f32::from_bits(0xb22d_91af);
const FRAC_PI_2: f32 = std::f32::consts::FRAC_PI_2;
const FRAC_PI_4: f32 = std::f32::consts::FRAC_PI_4;

// Cody-Waite split of ln(2). LN2_HI has 9 trailing zero mantissa bits so k*LN2_HI is exact.
const LN2_HI: f32 = 0.693_145_75;
const LN2_LO: f32 = 1.428_606_8e-6;
const LOG10_2_HI: f32 = 0.301_025_4;
const LOG10_2_LO: f32 = 4.605_039e-6;

/// Computes the natural logarithm of `x`.
#[doc(alias = "logf")]
#[doc(alias = "log")]
#[inline(always)]
#[allow(clippy::neg_cmp_op_on_partial_ord)]
pub fn ln(x: f32) -> f32 {
    log_family_wrapper!(x, ln_normal)
}

/// Core of ln for positive normal finite x only.
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn ln_normal(x: f32, koff: f32) -> f32 {
    let e = (x.to_bits() as i32).wrapping_sub(0x3f3504f3) >> 23;
    let m = f32::from_bits((x.to_bits() as i32).wrapping_sub(e << 23) as u32);
    let k = e as f32 + koff;
    let s = m - 1.0;
    let c: [f32; 8] = [
        -0.4999999,
        0.33333948,
        -0.2500179,
        0.1996218,
        -0.16569935,
        0.14916396,
        -0.1431012,
        0.08741673,
    ];
    let s2 = s * s;
    let s4 = s2 * s2;
    let l0 = fma(c[1], s, c[0]);
    let l1 = fma(c[3], s, c[2]);
    let l2 = fma(c[5], s, c[4]);
    let l3 = fma(c[7], s, c[6]);
    let a = s2 * l0;
    let u = fma(l3, s2, l2);
    let w = fma(u, s2, l1);
    let sq = fma(w, s4, a);
    let base = fma(k, LN2_LO, s);
    fma(k, LN2_HI, base + sq)
}

/// `ln` without domain checks: valid for positive normal finite `x`.
#[inline(always)]
pub fn ln_unchecked(x: f32) -> f32 {
    ln_normal(x, 0.0)
}

/// Computes the base-10 logarithm of `x`.
#[doc(alias = "log10f")]
#[inline(always)]
#[allow(clippy::neg_cmp_op_on_partial_ord)]
pub fn log10(x: f32) -> f32 {
    log_family_wrapper!(x, log10_normal)
}

/// Core of log10 for positive normal finite x only.
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn log10_normal(x: f32, koff: f32) -> f32 {
    let e = (x.to_bits() as i32).wrapping_sub(0x3f3504f3) >> 23;
    let m = f32::from_bits((x.to_bits() as i32).wrapping_sub(e << 23) as u32);
    let k = e as f32 + koff;
    let s = m - 1.0;
    let c: [f32; 8] = [
        -0.21714722,
        0.14476636,
        -0.10857988,
        0.086721875,
        -0.07200416,
        0.06459543,
        -0.06182998,
        0.038040668,
    ];
    let s2 = s * s;
    let s4 = s2 * s2;
    let l0 = fma(c[1], s, c[0]);
    let l1 = fma(c[3], s, c[2]);
    let l2 = fma(c[5], s, c[4]);
    let l3 = fma(c[7], s, c[6]);
    let a = fma(s2, l0, k * LOG10_2_LO);
    let u = fma(l3, s2, l2);
    let w = fma(u, s2, l1);
    let sq = fma(w, s4, a);
    // s*LOG10_E must stay inside the closing fma to avoid double rounding.
    fma(k, LOG10_2_HI, fma(s, std::f32::consts::LOG10_E, sq))
}

/// `log10` without domain checks: valid for positive normal finite `x`.
#[inline(always)]
pub fn log10_unchecked(x: f32) -> f32 {
    log10_normal(x, 0.0)
}

/// Computes `ln(1 + x)`, accurate for small `|x|`.
macro_rules! log1p_nonzero {
    ($x:expr) => {{
        let u = 1.0 + $x;
        let c = $x - (u - 1.0);
        let corr = c / u;
        let corr = if corr.is_finite() { corr } else { 0.0 };
        log_family_wrapper_no_denormal!(u, ln_normal) + corr
    }};
}

#[doc(alias = "log1pf")]
#[allow(clippy::neg_cmp_op_on_partial_ord)]
#[inline(always)]
pub fn log1p(x: f32) -> f32 {
    let normal = log1p_nonzero!(x);
    if x == 0.0 {
        x
    } else {
        normal
    }
}

/// Computes `ln(1 + x) - x`, accurate for small `|x|`.
#[inline(always)]
#[allow(clippy::neg_cmp_op_on_partial_ord)]
pub fn log1pmx(x: f32) -> f32 {
    let u = 1.0 + x;
    let c = x - (u - 1.0);
    let corr = c / u;
    let corr = if corr.is_finite() { corr } else { 0.0 };
    let e = (u.to_bits() as i32).wrapping_sub(0x3f3504f3) >> 23;
    let m = f32::from_bits((u.to_bits() as i32).wrapping_sub(e << 23) as u32);
    let k = e as f32;
    let w = m - 1.0;
    let z = if x.abs() < 0.5 { x } else { w };
    let z2 = z * z;
    const C1: f32 = -0.666_666_5;
    const C2: f32 = 0.499_999_52;
    const C3: f32 = -0.400_019_9;
    const C4: f32 = 0.333_370_42;
    const C5: f32 = -0.285_057_43;
    const C6: f32 = 0.249_024_99;
    const C7: f32 = -0.231_316_34;
    const C8: f32 = 0.211_575_52;
    const C9: f32 = -0.122_651_83;
    const C10: f32 = 0.099_456_76;
    const C11: f32 = -0.324_865_7;
    const C12: f32 = 0.320_051_94;
    let z4 = z2 * z2;
    let z8 = z4 * z4;
    let a0 = fma(C2, z, C1);
    let a1 = fma(C4, z, C3);
    let a2 = fma(C6, z, C5);
    let a3 = fma(C8, z, C7);
    let a4 = fma(C10, z, C9);
    let a5 = fma(C12, z, C11);
    let b0 = fma(a1, z2, a0);
    let b1 = fma(a3, z2, a2);
    let b2 = fma(a5, z2, a4);
    let c0 = fma(b1, z4, b0);
    let g = fma(b2, z8, c0);
    let q = fma(g, z, 1.0);
    let p = -0.5 * z2 * q;
    let big = log_family_edges!(u, {
        let t = fma(k, LN2_HI, -x) + w;
        t + (p + fma(k, LN2_LO, corr))
    });
    let normal = if x.abs() < 0.5 { p } else { big };
    if x == f32::INFINITY {
        f32::NEG_INFINITY
    } else {
        normal
    }
}

/// Computes `log2(1 + x)` (C23 `log2p1`).
#[allow(clippy::neg_cmp_op_on_partial_ord)]
#[inline(always)]
pub fn log2p1(x: f32) -> f32 {
    let u = 1.0 + x;
    let c = x - (u - 1.0);
    // Two-word log2(e) avoids bias for tiny |x| < 2^-24.
    let cu = c / u;
    let corr = fma(cu, LOG2_E, cu * LOG2_E_LO);
    let corr = if corr.is_finite() { corr } else { 0.0 };
    let normal = log_family_wrapper_no_denormal!(u, log_2_normal) + corr;
    if x == 0.0 {
        x
    } else {
        normal
    }
}

/// Computes `log10(1 + x)` (C23 `log10p1`).
#[allow(clippy::neg_cmp_op_on_partial_ord)]
#[inline(always)]
pub fn log10p1(x: f32) -> f32 {
    let u = 1.0 + x;
    let c = x - (u - 1.0);
    // Two-word log10(e) avoids bias for tiny |x| < 2^-24.
    let cu = c / u;
    let corr = fma(cu, std::f32::consts::LOG10_E, cu * LOG10_E_LO);
    let corr = if corr.is_finite() { corr } else { 0.0 };
    let normal = log_family_wrapper_no_denormal!(u, log10_normal) + corr;
    if x == 0.0 {
        x
    } else {
        normal
    }
}

/// Computes `e^x` via Cody-Waite range reduction.
#[doc(alias = "expf")]
#[inline(always)]
pub fn exp(x: f32) -> f32 {
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let r = fma(-k, LN2_HI, x);
    let r = fma(-k, LN2_LO, r);
    // c0 and c1 are pinned to 1.0 to preserve correct behavior near zero.
    let p = exp_r_poly!(r);
    let (t1, t2) = exp2_field_split(k);
    p * t1 * t2
}

/// Computes `e^x * 2^s`.
#[inline(always)]
pub fn exp_scaled(x: f32, s: i32) -> f32 {
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let r = fma(-k, LN2_HI, x);
    let r = fma(-k, LN2_LO, r);
    let p = exp_r_poly!(r);
    let (t1, t2) = exp2_field_split(k + s as f32);
    p * t1 * t2
}

/// Computes `e^x` via a single exponent field. Valid for `x` in `[-87.3, 88.7]`.
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn exp_narrow(x: f32) -> f32 {
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let r = fma(-k, LN2_HI, x);
    let r = fma(-k, LN2_LO, r);
    let p = exp_r_poly!(r);
    let exp2int = exp2int_field!(k);
    p * exp2int
}

/// Computes `e^x` across the full f32 domain with saturation.
#[inline(always)]
pub fn exp_checked(x: f32) -> f32 {
    exp_reduce!(x.clamp(EXP_CLAMP_LO, EXP_CLAMP_HI))
}

// Evaluates e^r - 1 as r + r^2*P(r) to avoid cancellation near zero.
macro_rules! expm1_p_poly {
    ($r:expr, $r2:expr) => {{
        let c: [f32; 5] = [0.5, 1.6666504e-1, 4.166_678e-2, 8.370_725e-3, 1.3916677e-3];
        let l1 = fma(c[1], $r, c[0]);
        let l2 = fma(c[3], $r, c[2]);
        let m = fma(c[4], $r2, l2);
        fma(m, $r2, l1)
    }};
}

macro_rules! expm1_r_poly {
    ($r:expr) => {{
        let r = $r;
        let r2 = r * r;
        fma(expm1_p_poly!(r, r2), r2, r)
    }};
}

// Exponent field for 2^(k-1).
const EXPM1_HALF_MAGIC: f32 = 12583038.0;

// Below 2^-125, expm1(x) is x and intermediate 2^(k-1) would denormalize.
const EXPM1_LINEAR: f32 = 2.0 * f32::MIN_POSITIVE;

/// Computes `e^x - 1`, avoiding cancellation near zero.
#[doc(alias = "expm1f")]
#[inline(always)]
pub fn expm1(x: f32) -> f32 {
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let r = fma(-k, LN2_HI, x);
    let r = fma(-k, LN2_LO, r);
    let e = expm1_r_poly!(r);
    let t = f32::from_bits((k + EXPM1_HALF_MAGIC).to_bits() << 23);
    let b = fma(e, t, t - 0.5);
    let b = b + b;
    if x.abs() < EXPM1_LINEAR {
        x
    } else {
        b
    }
}

/// `expm1` via a single exponent field. Valid for `x` in `[-87.3, 88.7]`.
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn expm1_narrow(x: f32) -> f32 {
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let r = fma(-k, LN2_HI, x);
    let r = fma(-k, LN2_LO, r);
    let e = expm1_r_poly!(r);
    let t = exp2int_field!(k);
    let b = fma(e, t, t - 1.0);
    f32::from_bits(b.to_bits() | (x.to_bits() & SIGN_MASK))
}

/// `expm1` across the full f32 domain with saturation.
#[inline(always)]
pub fn expm1_checked(x: f32) -> f32 {
    let xc = x.clamp(-86.0, 88.722_84);
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(xc, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let r = fma(-k, LN2_HI, xc);
    let r = fma(-k, LN2_LO, r);
    let e = expm1_r_poly!(r);
    let t = f32::from_bits((k + EXPM1_HALF_MAGIC).to_bits() << 23);
    let b = fma(e, t, t - 0.5);
    let b = b + b;
    if x.abs() < EXPM1_LINEAR {
        x
    } else {
        b
    }
}

/// `exp_m1_over_x` via a single exponent field. Valid for `x` in `[-87.3, 88.7]`.
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn exp_m1_over_x_narrow(x: f32) -> f32 {
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let r = fma(-k, LN2_HI, x);
    let r = fma(-k, LN2_LO, r);
    let r2 = r * r;
    let p = expm1_p_poly!(r, r2);
    let q = fma(r, p, 1.0);
    let e = fma(r2, p, r);
    let t = exp2int_field!(k);
    let b = fma(e, t, t - 1.0);
    if k == 0.0 {
        q
    } else {
        b / x
    }
}

macro_rules! exp2m1_f_poly {
    ($f:expr, $f2:expr) => {{
        let c: [f32; 5] = [
            2.402265e-1,
            5.55035e-2,
            9.618533e-3,
            1.3395752e-3,
            1.526698e-4,
        ];
        let l1 = fma(c[1], $f, c[0]);
        let l2 = fma(c[3], $f, c[2]);
        let m = fma(c[4], $f2, l2);
        fma(m, $f2, l1)
    }};
}

// Below 2^-124, intermediate 2^(k-1) would denormalize.
const EXP2M1_LINEAR: f32 = 4.0 * f32::MIN_POSITIVE;

/// Computes `2^x - 1`, avoiding cancellation near zero.
#[inline(always)]
pub fn exp2m1(x: f32) -> f32 {
    let xs = x.clamp(-126.0, 128.0);
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = (xs + ROUND_MAGIC) - ROUND_MAGIC;
    let f = xs - k;
    let f2 = f * f;
    let fl = f * LN_2;
    let big = fma(f2, exp2m1_f_poly!(f, f2), fl);
    let t = f32::from_bits((k + EXPM1_HALF_MAGIC).to_bits() << 23);
    let b = fma(big, t, t - 0.5);
    let b = b + b;
    if x.abs() < EXP2M1_LINEAR {
        fl
    } else {
        b
    }
}

macro_rules! exp10m1_d_poly {
    ($d:expr, $d2:expr) => {{
        let c: [f32; 5] = [2.650949, 2.0346525, 1.1712452, 0.5420898, 0.20779254];
        let l1 = fma(c[1], $d, c[0]);
        let l2 = fma(c[3], $d, c[2]);
        let m = fma(c[4], $d2, l2);
        fma(m, $d2, l1)
    }};
}

// Below 2^-125, intermediate 2^(k-1) would denormalize.
const EXP10M1_LINEAR: f32 = 2.0 * f32::MIN_POSITIVE;

/// Computes `10^x - 1`, avoiding cancellation near zero.
#[inline(always)]
pub fn exp10m1(x: f32) -> f32 {
    let xs = x.clamp(-37.0, 38.53184);
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(xs, std::f32::consts::LOG2_10, ROUND_MAGIC) - ROUND_MAGIC;
    let d = fma(-k, LOG10_2_HI, xs);
    let d = fma(-k, LOG10_2_LO, d);
    let d2 = d * d;
    let dl = d * std::f32::consts::LN_10;
    let big = fma(d2, exp10m1_d_poly!(d, d2), dl);
    let t = f32::from_bits((k + EXPM1_HALF_MAGIC).to_bits() << 23);
    let b = fma(big, t, t - 0.5);
    let b = b + b;
    if x.abs() < EXP10M1_LINEAR {
        dl
    } else {
        b
    }
}

// Split k = k1 + k2 into two representable exponent fields.
#[inline(always)]
fn exp2_field_split(k: f32) -> (f32, f32) {
    let a = fma(k, 0.5, EXP2INT_MAGIC);
    let k1 = a - EXP2INT_MAGIC;
    let b = (k - k1) + EXP2INT_MAGIC;
    let t1 = f32::from_bits(a.to_bits() << 23);
    let t2 = f32::from_bits(b.to_bits() << 23);
    (t1, t2)
}

// Evaluates exp(x) and exp(-x) simultaneously via even/odd poly decomposition.
#[inline(always)]
fn exp_pos_neg_half(x: f32) -> (f32, f32) {
    let (p_pos, p_neg, t1, t2, t1n, t2n) = exp_pos_neg_core!(x);
    (p_pos * t1 * t2, p_neg * t1n * t2n)
}

#[inline(always)]
fn exp_pos_neg_narrow_half(x: f32) -> (f32, f32) {
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(x, LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let r = fma(-k, LN2_HI, x);
    let r = fma(-k, LN2_LO, r);
    let c: [f32; 5] = [
        0.49999994 * 0.5,
        0.16666521 * 0.5,
        0.041_668_33 * 0.5,
        8.3687045e-3 * 0.5,
        1.3814511e-3 * 0.5,
    ];
    let r2 = r * r;
    let e = fma(fma(fma(c[4], r2, c[2]), r2, c[0]), r2, 0.5);
    let o = fma(fma(c[3], r2, c[1]), r2, 0.5);
    let p_pos = fma(r, o, e);
    let p_neg = fma(-r, o, e);
    let t = exp2int_field!(k);
    // Exact power-of-two reciprocal via bit subtraction.
    let tn = f32::from_bits(0x7F00_0000u32.wrapping_sub(t.to_bits()));
    (p_pos * t, p_neg * tn)
}

// sinh(x) ~ x * P(x^2) for |x| < 0.5.
#[inline(always)]
fn sinh_small(x: f32) -> f32 {
    let x2 = x * x;
    let c0 = 1.0f32;
    let c1 = 0.166_662_3_f32;
    let c2 = 0.008_396_464_f32;
    let p = fma(fma(c2, x2, c1), x2, c0);
    x * p
}

/// Computes the hyperbolic sine `sinh(x)`.
#[doc(alias = "sinhf")]
#[inline(always)]
pub fn sinh(x: f32) -> f32 {
    let a = sinh_small(x);
    let (ep, en) = exp_pos_neg_half(x);
    let b = ep - en;
    if x.abs() < 0.5 {
        a
    } else {
        b
    }
}

/// Computes the hyperbolic cosine `cosh(x)`.
#[doc(alias = "coshf")]
#[inline(always)]
pub fn cosh(x: f32) -> f32 {
    let (ep, en) = exp_pos_neg_half(x);
    ep + en
}

/// `sinh` via a single exponent field. Valid for `|x| <= 88.7`.
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn sinh_narrow(x: f32) -> f32 {
    let a = sinh_small(x);
    let (ep, en) = exp_pos_neg_narrow_half(x);
    let b = ep - en;
    if x.abs() < 0.5 {
        a
    } else {
        b
    }
}

/// `cosh` via a single exponent field. Valid for `|x| <= 88.7`.
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn cosh_narrow(x: f32) -> f32 {
    let (ep, en) = exp_pos_neg_narrow_half(x);
    ep + en
}

/// `exp_pos_neg_half` with clamped input, returning `(0.5*exp(x), 0.5*exp(-x))`.
#[inline(always)]
fn exp_pos_neg_checked_half(x: f32) -> (f32, f32) {
    let x = x.clamp(-170.0, 170.0);
    let (p_pos, p_neg, t1, t2, t1n, t2n) = exp_pos_neg_core!(x);
    (p_pos * t1 * t2, p_neg * t1n * t2n)
}

/// `sinh` across the full f32 domain with saturation.
#[inline(always)]
pub fn sinh_checked(x: f32) -> f32 {
    let a = sinh_small(x);
    let (ep, en) = exp_pos_neg_checked_half(x);
    let b = ep - en;
    if x.abs() < 0.5 {
        a
    } else {
        b
    }
}

/// `cosh` across the full f32 domain with saturation.
#[inline(always)]
pub fn cosh_checked(x: f32) -> f32 {
    let (ep, en) = exp_pos_neg_checked_half(x);
    ep + en
}

// Q(u) = (cosh(sqrt(u)) - 1 - u/2) / u^2 for u in [0, 4].
const COSHM1_Q_COEFFS: [f32; 4] = [0.04166663, 0.0013889787, 2.4732362e-5, 2.963389e-7];

/// Computes `cosh(x) - 1`, avoiding cancellation near zero.
#[inline(always)]
pub fn coshm1(x: f32) -> f32 {
    let u = x * x;
    let u2 = u * u;
    let c = COSHM1_Q_COEFFS;
    let a = fma(c[1], u, c[0]);
    let b = fma(c[3], u, c[2]);
    let q = fma(b, u2, a);
    // x*(0.5*x) prevents underflow when x*x would be denormal.
    let small = fma(u2, q, x * (0.5 * x));
    let big = cosh_checked(x) - 1.0;
    if x.abs() < 2.0 {
        small
    } else {
        big
    }
}

/// Throughput-optimized `sinh(x)`.
#[inline(always)]
pub fn sinh_throughput(x: f32) -> f32 {
    let a = sinh_small(x);
    let e = exp(x);
    let b = 0.5 * (e - 1.0 / e);
    if x.abs() < 0.5 {
        a
    } else {
        b
    }
}

/// Throughput-optimized `cosh(x)`.
#[inline(always)]
pub fn cosh_throughput(x: f32) -> f32 {
    let e = exp(x);
    0.5 * (e + 1.0 / e)
}

/// Computes the hyperbolic tangent `tanh(x)`.
#[doc(alias = "tanhf")]
#[inline(always)]
pub fn tanh(x: f32) -> f32 {
    let xc = x.clamp(-43.5, 44.0);

    // Small arm denominator: x*coth(x) over [0, 0.8].
    const COTH1: f32 = 0.33333313;
    const COTH2: f32 = -0.02221908;
    const COTH3: f32 = 0.0021010686;
    const COTH4: f32 = -0.00018176674;
    let x2 = xc * xc;
    let dp = fma(fma(fma(COTH4, x2, COTH3), x2, COTH2), x2, COTH1);
    let ds = fma(dp, x2, 1.0);

    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    const LOG2_E_X2: f32 = 2.0 * LOG2_E;
    let k = fma(xc, LOG2_E_X2, ROUND_MAGIC) - ROUND_MAGIC;
    const LN2_HI_HALF: f32 = LN2_HI / 2.0;
    const LN2_LO_HALF: f32 = LN2_LO / 2.0;
    let rh = fma(-k, LN2_HI_HALF, xc);
    let rh = fma(-k, LN2_LO_HALF, rh);
    const D0: f32 = 4.0 * 4.999_93e-1;
    const D1: f32 = 8.0 * 1.6667245e-1;
    const D2: f32 = 16.0 * 4.188_381e-2;
    const D3: f32 = 32.0 * 8.300_99e-3;
    let rh2 = rh * rh;
    let l0 = fma(2.0, rh, 1.0);
    let l1 = fma(D1, rh, D0);
    let l2 = fma(D3, rh, D2);
    let m = fma(l2, rh2, l1);
    let p = fma(m, rh2, l0);
    let exp2int = exp2int_field!(k);
    let b = fma(p, exp2int, -1.0);

    // Select before division to execute divide only once.
    let small = xc.abs() < 0.8;
    let n = if small { xc } else { b };
    let d = if small { ds } else { b + 2.0 };
    n / d
}

/// Computes the derivative of `tanh`: `1 - tanh(x)^2`.
#[inline(always)]
pub fn tanh_grad(x: f32) -> f32 {
    let q = exp_checked(-2.0 * x.abs());
    4.0 * q / ((1.0 + q) * (1.0 + q))
}

/// Logistic sigmoid: `1 / (1 + exp(-x))`.
#[inline(always)]
pub fn sigmoid(x: f32) -> f32 {
    let xc = x.clamp(-88.722_84, 87.0);
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    let k = fma(xc, -LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let t1 = fma(k, LN2_HI, xc);
    let t2 = fma(k, LN2_LO, t1);
    let r = -t2;
    let p = exp_r_poly!(r);
    let exp2int = exp2int_field!(k);
    let e = p * exp2int;
    1.0 / (1.0 + e)
}

/// Fast piecewise approximation of sigmoid.
#[inline(always)]
pub fn sigmoid_fast(x: f32) -> f32 {
    let xc = x.clamp(-3.288051, 3.288051);
    let poly = fma(-0.006715598, xc * xc, 0.2178126);
    fma(poly, xc, 0.5).clamp(0.0, 1.0)
}

/// Derivative of sigmoid: `sigmoid(x) * (1 - sigmoid(x))`.
#[inline(always)]
pub fn sigmoid_grad(x: f32) -> f32 {
    let e = exp_checked(-x.abs());
    e / ((1.0 + e) * (1.0 + e))
}

/// `log1p` specialized for callers guaranteeing `e` in `(0, 1]`, bypassing domain checks.
#[inline(always)]
fn log1p_unit(e: f32) -> f32 {
    let c: [f32; 10] = [
        -0.499_999_88,
        0.333_326_9,
        -0.249_885_8,
        0.198_979_3,
        -0.161_293_28,
        0.124_671_06,
        -0.083_073_795,
        0.041_981_1,
        -0.013_631_347,
        0.002_072_919_4,
    ];
    let e2 = e * e;
    let e4 = e2 * e2;
    let l0 = fma(c[1], e, c[0]);
    let l1 = fma(c[3], e, c[2]);
    let l2 = fma(c[5], e, c[4]);
    let l3 = fma(c[7], e, c[6]);
    let l4 = fma(c[9], e, c[8]);
    let r0 = fma(l1, e2, l0);
    let r1 = fma(l3, e2, l2);
    let r1b = fma(l4, e4, r1);
    let q = fma(r1b, e4, r0);
    fma(e2, q, e)
}

#[inline(always)]
fn softplus_impl(x: f32) -> f32 {
    let ax = x.abs();
    let e = exp_narrow(-ax.min(87.0));
    let corr = if ax > 87.0 { 0.0 } else { log1p_unit(e) };
    let normal = x.max(0.0) + corr;
    if x.is_nan() {
        f32::NAN
    } else {
        normal
    }
}

// Computes exp(-a) * 2^64 to keep intermediate normal across denormal range.
macro_rules! exp_neg_scaled64 {
    ($a:expr) => {{
        const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
        let a = $a;
        let k = fma(a, -LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
        let t1 = fma(k, LN2_HI, a);
        let t2 = fma(k, LN2_LO, t1);
        let p = exp_r_poly!(-t2);
        p * exp2int_field!(k + 64.0)
    }};
}

#[inline(always)]
fn softplus_checked_impl(x: f32) -> f32 {
    const P64: f32 = 5.421_011e-20; // 2^-64, exact
    let e = exp_neg_scaled64!(x.abs().min(105.0)) * P64;
    let normal = x.max(0.0) + log1p_unit(e);
    if x.is_nan() {
        f32::NAN
    } else {
        normal
    }
}

/// Log-sigmoid: `ln(sigmoid(x)) = -softplus(-x)`.
#[inline(always)]
pub fn logsigmoid(x: f32) -> f32 {
    -softplus_impl(-x)
}

/// Log-sigmoid across the full f32 domain.
#[inline(always)]
pub fn logsigmoid_checked(x: f32) -> f32 {
    -softplus_checked_impl(-x)
}

/// Gaussian Error Linear Unit (GELU): `x * Phi(x)`.
#[inline(always)]
pub fn gelu(x: f32) -> f32 {
    let xa = x.abs();
    let r = erfcx_pos(xa * std::f32::consts::FRAC_1_SQRT_2);
    let xs = if xa > NORM_CDF_XS_CLAMP {
        NORM_CDF_XS_CLAMP
    } else {
        xa
    };
    let h = 0.5 * xs;
    let p = h * xs;
    let pe = fma(h, xs, -p);
    let e = exp_reduce!(-p);
    // Preserves sign of -0.0 in closing fma.
    let s = (x.to_bits() as i32 >> 31) as u32;
    let addend = f32::from_bits(x.to_bits() & (!s | 0x8000_0000));
    fma(-(h * e), fma(-r, pe, r), addend)
}

/// SiLU / Swish activation: `x * sigmoid(x)`.
#[inline(always)]
pub fn silu(x: f32) -> f32 {
    let normal = x * sigmoid(x);
    if x == f32::NEG_INFINITY {
        0.0
    } else {
        normal
    }
}

/// SiLU across the full f32 domain.
#[inline(always)]
pub fn silu_checked(x: f32) -> f32 {
    const ROUND_MAGIC: f32 = 12582912.0; // 1.5 * 2^23
    const TWO64: f32 = 18446744073709551616.0; // 2^64, exact
    let ax = x.abs();
    let k = fma(ax, -LOG2_E, ROUND_MAGIC) - ROUND_MAGIC;
    let t1 = fma(k, LN2_HI, ax);
    let t2 = fma(k, LN2_LO, t1);
    let r = -t2;
    let p = exp_r_poly!(r);
    let e2 = p * exp2int_field!(k + 64.0); // exp(-|x|) * 2^64
    let num = x * if x < 0.0 { e2 } else { TWO64 };
    let normal = num / (TWO64 + e2);
    if ax > 110.0 {
        x.max(-0.0)
    } else {
        normal
    }
}

/// Softsign: `x / (1 + |x|)`.
#[inline(always)]
pub fn softsign(x: f32) -> f32 {
    let normal = x / (1.0 + x.abs());
    if x.is_infinite() {
        x.signum()
    } else {
        normal
    }
}

/// Computes `sqrt(1 + x) - 1`, avoiding catastrophic cancellation near zero.
/// Domain: `x >= -1.0`.
#[inline(always)]
pub fn sqrt1pm1(x: f32) -> f32 {
    let normal = x / ((1.0 + x).sqrt() + 1.0);
    if x.is_infinite() {
        x
    } else {
        normal
    }
}

/// Computes `x^(3/2) = x * sqrt(x)` for `x >= 0`.
#[inline(always)]
pub fn pow_3_2(x: f32) -> f32 {
    x * x.sqrt()
}

/// Core of `x^(2/3)` for positive normal `x` via direct minimax fit of `(1+r)^(-2/3)`.
#[inline(always)]
fn pow_2_3_normal(a: f32, scale: f32, scale3: f32) -> f32 {
    let ax = a.to_bits();
    let rcp = 1.0 / a; // independent of the seed chain, starts immediately
    let s = f32::from_bits(ax / 3 + 0x2a509a07u32);
    let s2u = s * s;
    let e2 = fma(s, s, -s2u) * scale3;
    let d = fma(s2u, s, -a);
    let r = d * rcp;
    let s2 = s2u * scale; // exact: scale is a power of two
    let c1 = -0.666_666_8_f32;
    let c2 = 0.555_541_16_f32;
    let c3 = -0.493_700_95_f32;
    let c4 = 0.457_742_93_f32;
    let c5 = -0.440_396_8_f32;
    let p = fma(fma(fma(fma(c5, r, c4), r, c3), r, c2), r, c1);
    s2 + fma(s2 * r, p, e2)
}

/// Computes `x^(2/3)` for `x >= 0`.
#[inline(always)]
pub fn pow_2_3(x: f32) -> f32 {
    // Rescale denormals: scales x by 2^24, result comes back scaled by 2^16.
    const OUT: f32 = 1.525_878_9e-5; // 2^-16
    const OUT3: f32 = OUT * (1.0 / 3.0);
    let ax = x.to_bits() & !SIGN_MASK;
    let a = f32::from_bits(ax);
    let tiny = ax < 0x0080_0000;
    let ascaled = if tiny { a * 16777216.0 } else { a };
    let scale = if tiny { OUT } else { 1.0 };
    let scale3 = if tiny { OUT3 } else { 1.0 / 3.0 };
    let r = pow_2_3_normal(ascaled, scale, scale3);
    if ax == 0 || ax >= EXPONENT_MASK {
        a + a
    } else {
        r
    }
}

/// Hermite smoothstep on `[edge0, edge1]`.
#[inline(always)]
pub fn smoothstep(edge0: f32, edge1: f32, x: f32) -> f32 {
    let t = ((x - edge0) / (edge1 - edge0)).clamp(0.0, 1.0);
    t * t * fma(-2.0, t, 3.0)
}

/// Perlin smootherstep on `[edge0, edge1]`.
#[inline(always)]
pub fn smootherstep(edge0: f32, edge1: f32, x: f32) -> f32 {
    let t = ((x - edge0) / (edge1 - edge0)).clamp(0.0, 1.0);
    let p = fma(t, fma(t, 6.0, -15.0), 10.0);
    t * t * t * p
}

/// Computes the inverse hyperbolic sine `asinh(x)`.
#[doc(alias = "asinhf")]
#[inline(always)]
pub fn asinh(x: f32) -> f32 {
    let ax = x.abs();
    let small = ax < 2048.0;
    let ax2 = ax * ax;
    let direct_sq = (ax2 + 1.0).sqrt();
    // For ax >= 2048, sqrt(ax^2+1) ~= ax + 1/(2*ax) is exact to f32 precision.
    let inv_ax = 1.0 / ax;
    let rescaled_sq = fma(0.5, inv_ax, ax);
    let sq = if small { direct_sq } else { rescaled_sq };
    let sm1 = if small { ax2 / (sq + 1.0) } else { sq - 1.0 };
    let d = ax + sm1;
    let finite_d = d.is_finite();
    let u = 1.0 + d;
    let c = d - (u - 1.0);
    let corr = c / u;
    let arg = if finite_d { u } else { ax };
    let koff = if finite_d { 0.0 } else { 1.0 };
    let shared = ln_normal(arg, koff);
    let combined = if finite_d { shared + corr } else { shared };
    let combined = if x.is_finite() { combined } else { ax };
    mulsign(combined, x)
}

/// Computes the inverse hyperbolic cosine `acosh(x)` for `x >= 1`.
#[doc(alias = "acoshf")]
#[inline(always)]
pub fn acosh(x: f32) -> f32 {
    // x*x - 1.0 must be a single FMA to avoid catastrophic cancellation near x=1.
    let direct = fma(x, x, -1.0).sqrt();
    // For x >= 2048, sqrt(x^2-1) ~= x - 1/(2*x) is exact to f32 precision.
    let inv_x = 1.0 / x;
    let rescaled = fma(-0.5, inv_x, x);
    let s = if x < 2048.0 { direct } else { rescaled };
    let d = (x - 1.0) + s;
    let finite_d = d.is_finite();
    let u = 1.0 + d;
    let c = d - (u - 1.0);
    let corr = c / u;
    let arg = if finite_d { u } else { x };
    let koff = if finite_d { 0.0 } else { 1.0 };
    let shared = ln_normal(arg, koff);
    let combined = if finite_d { shared + corr } else { shared };
    let combined = if x.is_finite() { combined } else { x };
    if x < 1.0 {
        f32::NAN
    } else {
        combined
    }
}

// atanh(x) ~ x * P(x^2) on |x| < 0.25.
#[inline(always)]
fn atanh_small(x: f32) -> f32 {
    let x2 = x * x;
    let c0 = 0.3333338f32;
    let c1 = 0.19981473f32;
    let c2 = 0.15260197f32;
    let p = fma(fma(c2, x2, c1), x2, c0);
    fma(x * x2, p, x)
}

/// Computes the inverse hyperbolic tangent `atanh(x)` for `|x| < 1`.
#[doc(alias = "atanhf")]
#[inline(always)]
#[allow(clippy::neg_cmp_op_on_partial_ord)]
pub fn atanh(x: f32) -> f32 {
    let a = x.abs();
    let small = atanh_small(x);
    let v = 2.0 * a / (1.0 - a);
    let u = 1.0 + v;
    let c = v - (u - 1.0);
    let corr = c / u;
    let corr = if corr.is_finite() { corr } else { 0.0 };
    let l = ln_normal(u, 0.0) + corr;
    let l = if u <= 0.0 { f32::NAN } else { l };
    let l = if !(u < f32::INFINITY) { u * u } else { l };
    let big = mulsign(0.5 * l, x);
    if a < 0.25 {
        small
    } else {
        big
    }
}

// asin(sqrt(t))/sqrt(t) on t in [0, 0.25].
#[inline(always)]
fn acos_poly(t: f32) -> f32 {
    let u = 4.2285666e-2f32;
    let u = fma(u, t, 2.407_541e-2);
    let u = fma(u, t, 4.5502156e-2);
    let u = fma(u, t, 7.494872e-2);
    let u = fma(u, t, 1.6666777e-1);
    fma(u, t, 1.0)
}

/// Dedicated asin.
#[inline(always)]
fn asin_poly(x: f32) -> f32 {
    let u = -2.0342327e-3f32;
    let u = fma(u, x, 1.20692495e-2);
    let u = fma(u, x, -3.609505e-2);
    let u = fma(u, x, 8.2684554e-2);
    let u = fma(u, x, -2.1304381e-1);
    fma(u, x, 1.5706329)
}

/// Computes `acos(x)` in radians for `x` in `[-1, 1]`. Result is in `[0, pi]`.
#[doc(alias = "acosf")]
#[inline(always)]
pub fn acos(x: f32) -> f32 {
    const PI: f32 = std::f32::consts::PI;
    const LO_RATIO: f32 = -2.7827534e-8;
    let na = f32::from_bits(x.to_bits() | SIGN_MASK); // -|x|
    let small = na > -0.5;
    // t = x^2 below crossover, (1-|x|)/2 above it.
    let t = if small { x * x } else { fma(na, 0.5, 0.5) };
    let y = t.sqrt();
    let m = mulsign(if small { na } else { y + y }, x);
    let c = if small {
        FRAC_PI_2
    } else if x < 0.0 {
        PI
    } else {
        0.0
    };
    fma(m, acos_poly(t), c * LO_RATIO) + c
}

/// Computes `acos(x)` in degrees for `x` in `[-1, 1]`.
#[inline(always)]
pub fn acosd(x: f32) -> f32 {
    acos(x) * RAD_TO_DEG_HI
}

#[inline(always)]
fn acospi_poly(x: f32) -> f32 {
    let u = 7.541_4e-4_f32;
    let u = fma(u, x, -3.6262998e-3);
    let u = fma(u, x, 8.662372e-3);
    let u = fma(u, x, -1.55939115e-2);
    let u = fma(u, x, 2.8268332e-2);
    let u = fma(u, x, -6.830643e-2);
    fma(u, x, 5e-1)
}

/// Computes `acos(x) / pi` in half-turns for `x` in `[-1, 1]`.
#[inline(always)]
pub fn acospi(x: f32) -> f32 {
    let a = x.abs();
    let y = (1.0 - a).sqrt() * acospi_poly(a);
    mulsign(y, x + 0.0) + if x < 0.0 { 1.0 } else { 0.0 }
}

// asin(x) ~ x * P(x^2) on |x| < 0.5.
#[inline(always)]
fn asin_small(x: f32) -> f32 {
    let x2 = x * x;
    let c0 = 1.0f32;
    let c1 = 0.16666752f32;
    let c2 = 0.074952975f32;
    let c3 = 0.04547038f32;
    let c4 = 0.02417949f32;
    let c5 = 0.042166352f32;
    let p = fma(
        fma(fma(fma(fma(c5, x2, c4), x2, c3), x2, c2), x2, c1),
        x2,
        c0,
    );
    x * p
}

/// Computes `asin(x)` in radians for `x` in `[-1, 1]`.
#[doc(alias = "asinf")]
#[inline(always)]
pub fn asin(x: f32) -> f32 {
    let a = x.abs();
    let small = asin_small(x);
    let big = mulsign(fma(-(1.0 - a).sqrt(), asin_poly(a), FRAC_PI_2), x);
    if a < 0.5 {
        small
    } else {
        big
    }
}

// 180/pi split into double-f32 to prevent rounding bias in degree conversions.
const RAD_TO_DEG_HI: f32 = 57.295_776;
const RAD_TO_DEG_LO: f32 = 3.1458948e-6;

// 1/pi split into double-f32 for half-turn conversions.
const FRAC_1_PI_LO: f32 = 1.28412765e-8;

/// Computes `asin(x)` in degrees for `x` in `[-1, 1]`.
#[inline(always)]
pub fn asind(x: f32) -> f32 {
    let y = asin(x);
    fma(y, RAD_TO_DEG_HI, y * RAD_TO_DEG_LO)
}

#[inline(always)]
fn asinpi_small(x: f32) -> f32 {
    let x2 = x * x;
    let c1 = 0.053051922f32;
    let c2 = 0.023858273f32;
    let c3 = 0.014473671f32;
    let c4 = 0.0076965746f32;
    let c5 = 0.013421959f32;
    let t = fma(
        fma(fma(fma(fma(c5, x2, c4), x2, c3), x2, c2), x2, c1),
        x2,
        FRAC_1_PI_LO,
    );
    fma(x, FRAC_1_PI, x * t)
}

#[inline(always)]
fn asinpi_poly(x: f32) -> f32 {
    let u = -0.0006497958f32;
    let u = fma(u, x, 0.0038506954);
    let u = fma(u, x, -0.011503078);
    let u = fma(u, x, 0.026329458);
    let u = fma(u, x, -0.06781764);
    fma(u, x, 0.4999485)
}

/// Computes `asin(x) / pi` in half-turns for `x` in `[-1, 1]`.
#[inline(always)]
pub fn asinpi(x: f32) -> f32 {
    let a = x.abs();
    let small = asinpi_small(x);
    let big = mulsign(fma(-(1.0 - a).sqrt(), asinpi_poly(a), 0.5), x);
    if a < 0.5 {
        small
    } else {
        big
    }
}

// 3/3 Pade rational approximation of atan on [0, 1].
#[inline(always)]
fn atan_poly(x: f32) -> f32 {
    let a2 = 0.008_830_043;
    let a1 = 0.284_977_85;
    let a0 = 1.127_171_2;
    let b2 = 5.0166193e-2;
    let b1 = 5.718157e-1;
    let b0 = 1.4605043e0;
    let x2 = x * x;
    let numer = fma(x * x2, fma(fma(a2, x2, a1), x2, a0), x);
    let denom = fma(fma(fma(b2, x2, b1), x2, b0), x2, 1.0);
    numer / denom
}

/// Computes `atan(x)` in radians.
#[doc(alias = "atanf")]
#[inline(always)]
pub fn atan(x: f32) -> f32 {
    let a = x.abs();
    let y = a.min(1.0 / a);
    let y = atan_poly(y);
    let y = if a < 1.0 { y } else { FRAC_PI_2 - y };
    mulsign(y, x)
}

/// Computes `atan(x)` in degrees.
#[inline(always)]
pub fn atand(x: f32) -> f32 {
    let y = atan(x);
    fma(y, RAD_TO_DEG_HI, y * RAD_TO_DEG_LO)
}

/// `atan(x)` for `|x| <= 1`.
#[doc(hidden)] // pub only so examples/mca_target.rs can benchmark it directly
#[inline(always)]
pub fn atan_bounded(x: f32) -> f32 {
    mulsign(atan_poly(x.abs()), x)
}

/// Latency-optimized division-free `atan(x)`.
#[inline(always)]
pub fn atan_latency(x: f32) -> f32 {
    let a = x.abs();
    let r = a.min(1.0 / a);
    let r2 = r * r;
    let r4 = r2 * r2;
    let c0 = -0.33333167;
    let c1 = 0.19994265;
    let c2 = -0.14216055;
    let c3 = 0.10689225;
    let c4 = -0.07608681;
    let c5 = 0.04395558;
    let c6 = -0.01687014;
    let c7 = 0.0030569038;
    let lo = fma(fma(fma(c2, r2, c1), r2, c0), r2, 1.0);
    let hi = fma(fma(fma(fma(c7, r2, c6), r2, c5), r2, c4), r2, c3);
    let p = fma(hi, r4 * r4, lo) * r;
    let sp = mulsign(p, x);
    let hpisignx = mulsign(FRAC_PI_2, x);
    if a < 1.0 {
        sp
    } else {
        hpisignx - sp
    }
}

/// Computes the four-quadrant arctangent `atan2(y, x)` in radians.
#[doc(alias = "atan2f")]
#[inline(always)]
pub fn atan2(y: f32, x: f32) -> f32 {
    let nonzerox = x != 0.0;
    let nonzeroy = y != 0.0;
    let bothzero = !nonzerox && !nonzeroy;
    let hpisignx = if nonzerox || bothzero {
        mulsign(FRAC_PI_2, x)
    } else {
        0.0
    };
    let correction = mulsign(FRAC_PI_2 - hpisignx, y);
    let r = if nonzerox {
        atan(y / x) + correction
    } else {
        correction
    };
    let r = if y.is_nan() { f32::NAN } else { r };
    // atan2(+-inf, +-inf): quadrant-dependent result required by IEEE 754 / C99.
    let bothinf = x.is_infinite() && y.is_infinite();
    let inf_result = mulsign(
        if x.is_sign_negative() {
            3.0 * FRAC_PI_4
        } else {
            FRAC_PI_4
        },
        y,
    );
    if bothinf {
        inf_result
    } else {
        r
    }
}

/// Latency-optimized `atan2(y, x)`.
#[inline(always)]
pub fn atan2_latency(y: f32, x: f32) -> f32 {
    let nonzerox = x != 0.0;
    let nonzeroy = y != 0.0;
    let bothzero = !nonzerox && !nonzeroy;
    let hpisignx = if nonzerox || bothzero {
        mulsign(FRAC_PI_2, x)
    } else {
        0.0
    };
    let correction = mulsign(FRAC_PI_2 - hpisignx, y);
    let r = if nonzerox {
        atan_latency(y / x) + correction
    } else {
        correction
    };
    let r = if y.is_nan() { f32::NAN } else { r };
    let bothinf = x.is_infinite() && y.is_infinite();
    let inf_result = mulsign(
        if x.is_sign_negative() {
            3.0 * FRAC_PI_4
        } else {
            FRAC_PI_4
        },
        y,
    );
    if bothinf {
        inf_result
    } else {
        r
    }
}

/// `atan2` without special zero/infinite case handling.
#[inline(always)]
pub fn atan2_unchecked(y: f32, x: f32) -> f32 {
    let hpisignx = mulsign(FRAC_PI_2, x);
    let correction = mulsign(FRAC_PI_2 - hpisignx, y);
    atan(y / x) + correction
}

/// `atan2` folded into `[0, 2*pi)` (single positive turn).
#[inline(always)]
pub fn atan2_pos(y: f32, x: f32) -> f32 {
    let r = atan2(y, x);
    if y.is_sign_negative() {
        r + std::f32::consts::TAU
    } else {
        r
    }
}

/// Computes `tan(x)` in radians for `|x| < 2^23 * pi`.
#[doc(alias = "tanf")]
#[inline(always)]
pub fn tan(x: f32) -> f32 {
    // q = round(2x/pi) exactly: a coarse round to multiples of 4, then a fine
    // one of the remainder, as in `sin`.
    let nb = fma(x, std::f32::consts::FRAC_2_PI, ROUND_MAGIC_4);
    let n = nb - ROUND_MAGIC_4;
    let f = fma(x, std::f32::consts::FRAC_2_PI, -n);
    let fc = fma(x, 2.0 * RPI_LO, f);
    let qb = fc + ROUND_MAGIC;
    let q = (n - ROUND_MAGIC) + qb;
    let r = fma(q, -0.5 * PI_A, x);
    let r = fma(q, -0.5 * PI_B, r);
    let r = fma(q, -0.5 * PI_C, r);
    let r = fma(q, -0.5 * PI_D, r);
    // Odd q: tan(x) = -cos(r) / sin(r).
    let odd = qb.to_bits() & 1 != 0;
    let flip = core::hint::select_unpredictable(odd, SIGN_MASK, 0);
    let (s, c) = sincos_quarter(f32::from_bits(r.to_bits() ^ flip));
    let (n, d) = core::hint::select_unpredictable(odd, (c, s), (s, c));
    n / d
}

#[inline(always)]
fn erf_poly(x: f32, x2: f32) -> f32 {
    let a6 = 2.838_853e-4_f32;
    let a5 = -4.4954885e-3f32;
    let a4 = 3.273_625e-2_f32;
    let a3 = -1.5164591e-1f32;
    let a2 = -9.171_398e-1_f32;
    let a1 = -1.6281782f32;
    let a0 = 2.2989703e-5f32;
    let x4 = x2 * x2;
    let b0 = fma(a1, x, a0);
    let b1 = fma(a3, x, a2);
    let b2 = fma(a5, x, a4);
    let c0 = fma(b1, x2, b0);
    let c1 = fma(a6, x2, b2);
    fma(c1, x4, c0)
}

/// Error function `erf(x)`.
#[doc(alias = "erff")]
#[inline(always)]
pub fn erf(x: f32) -> f32 {
    let xa = x.abs();
    let xa_bounded = if xa > 10.0 { 10.0 } else { xa };
    let x2 = xa_bounded * xa_bounded;
    // Pade constant term 2/sqrt(pi) split into two words to avoid rounding bias.
    let numer = fma(
        x,
        f32::from_bits(0x3f906ebb),
        x * fma(f32::from_bits(0x3f174f6e), x2, f32::from_bits(0xb37bd649)),
    );
    let denom = fma(
        fma(f32::from_bits(0x3e3e2be3), x2, f32::from_bits(0x3f5b6db7)),
        x2,
        1.0,
    );
    let a = numer / denom;
    let b = mulsign(1.0 - exp2(erf_poly(xa_bounded, x2)), x);
    if xa < 0.28 {
        a
    } else {
        b
    }
}

// Clamps on |x| ensure exponents stay within exp_reduce!'s valid range.
const ERFC_XS_CLAMP: f32 = 10.21;
const _: () = assert!((ERFC_XS_CLAMP as f64) * (ERFC_XS_CLAMP as f64) <= -(EXP_CLAMP_LO as f64));
// Below 150*ln2, e^-p underflows to 0.
const _: () = assert!((ERFC_XS_CLAMP as f64) * (ERFC_XS_CLAMP as f64) > 103.97207708399179);

// Upper bound ensuring 2*e^(xs^2) overflows to +inf within f32 range.
const ERFCX_XS_CLAMP: f32 = 9.41;
const _: () = assert!((ERFCX_XS_CLAMP as f64) * (ERFCX_XS_CLAMP as f64) <= EXP_CLAMP_HI as f64);
const _: () = assert!((ERFCX_XS_CLAMP as f64) * (ERFCX_XS_CLAMP as f64) > 88.02969187150839);

// erfcx(xa) = erfc(xa) * exp(xa^2) for xa >= 0.
#[inline(always)]
fn erfcx_pos(xa: f32) -> f32 {
    let v0 = 1.0 / (2.0 + xa);
    let v = if xa <= 2.0 {
        fma(-0.5 * xa, v0, 0.5)
    } else {
        v0
    };
    // c[0] is the low word of 1/sqrt(pi); high word is applied in the final fma.
    let c: [f32; 11] = [
        f32::from_bits(0xb2fbd649),
        1.1283773,
        1.974964,
        2.8070478,
        2.9756768,
        -3.7488432,
        17.02367,
        -117.490135,
        255.59447,
        -243.95302,
        90.238,
    ];
    let v2 = v * v;
    let v4 = v2 * v2;
    let p01 = fma(c[1], v, c[0]);
    let p23 = fma(c[3], v, c[2]);
    let p45 = fma(c[5], v, c[4]);
    let p67 = fma(c[7], v, c[6]);
    let t9 = fma(c[10], v, c[9]);
    let t8 = fma(t9, v, c[8]);
    let lo = fma(p23, v2, p01);
    let hi = fma(p67, v2, p45);
    fma(
        v,
        f32::from_bits(0x3f106ebb),
        v * fma(fma(t8, v4, hi), v4, lo),
    )
}

/// Complementary error function `erfc(x) = 1 - erf(x)`.
#[doc(alias = "erfcf")]
#[inline(always)]
pub fn erfc(x: f32) -> f32 {
    // Reflection for x < 0: erfc(-x) = 2 - erfc(x).
    let w = f32::from_bits((x.to_bits() >> 1) & 0x4000_0000);
    let xa = x.abs();
    let xs = if xa > ERFC_XS_CLAMP {
        ERFC_XS_CLAMP
    } else {
        xa
    };
    let p = xs * xs;
    let pe = fma(xs, xs, -p);
    let r = erfcx_pos(xa);
    let e = exp_reduce!(-p);
    fma(mulsign(e, x), fma(-r, pe, r), w)
}

// erfinv(x) = x*P(x^2) for |x| <= 0.7.
#[inline(always)]
fn erfinv_central_poly_m1(u: f32) -> f32 {
    let c: [f32; 9] = [
        -0.11377308,
        2.3201263e-1,
        1.2761366e-1,
        8.534858e-2,
        7.737176e-2,
        -1.8498806e-2,
        2.6768243e-1,
        -3.5444248e-1,
        3.377774e-1,
    ];
    let u2 = u * u;
    let u4 = u2 * u2;
    let l0 = fma(c[1], u, c[0]);
    let l1 = fma(c[3], u, c[2]);
    let l2 = fma(c[5], u, c[4]);
    let l3 = fma(c[7], u, c[6]);
    let r0 = fma(l1, u2, l0);
    let r1 = fma(l3, u2, l2);
    fma(c[8], u4 * u4, fma(r1, u4, r0))
}

// erfinv(x) = sign(x)*sqrt(w)*Q with w = -ln(1-x^2) for |x| > 0.7.
#[inline(always)]
fn erfinv_tail_poly_m1(t: f32) -> f32 {
    let c: [f32; 12] = [
        -0.103669584,
        0.019397417,
        0.007896575,
        -0.0021255405,
        -0.0007704421,
        -1.969348e-5,
        -2.4305053e-5,
        0.00033242995,
        -0.00024495937,
        7.880177e-5,
        -1.2436317e-5,
        7.908929e-7,
    ];
    let t2 = t * t;
    let t4 = t2 * t2;
    let l0 = fma(c[1], t, c[0]);
    let l1 = fma(c[3], t, c[2]);
    let l2 = fma(c[5], t, c[4]);
    let l3 = fma(c[7], t, c[6]);
    let l4 = fma(c[9], t, c[8]);
    let l5 = fma(c[11], t, c[10]);
    let r0 = fma(l1, t2, l0);
    let r1 = fma(l3, t2, l2);
    let r2 = fma(l5, t2, l4);
    fma(fma(r2, t4, r1), t4, r0)
}

// Far-tail polynomial in t = sqrt(w) - 7 for erfc_inv/probit.
#[inline(always)]
fn erfinv_far_poly_m1(t: f32) -> f32 {
    let c: [f32; 8] = [
        -0.018711485,
        0.0039010558,
        -0.00063050824,
        9.089283e-5,
        -1.2126031e-5,
        1.5352017e-6,
        -1.7609956e-7,
        1.2450847e-8,
    ];
    let t2 = t * t;
    let t4 = t2 * t2;
    let l0 = fma(c[1], t, c[0]);
    let l1 = fma(c[3], t, c[2]);
    let l2 = fma(c[5], t, c[4]);
    let l3 = fma(c[7], t, c[6]);
    let r0 = fma(l1, t2, l0);
    let r1 = fma(l3, t2, l2);
    fma(r1, t4, r0)
}

const ERFC_INV_W_FAR: f32 = 16.0;

// |erfc_inv(n)| for n in (0, 1].
#[inline(always)]
fn erfc_inv_half(n: f32) -> f32 {
    let x = 1.0 - n;
    let central = fma(x, erfinv_central_poly_m1(x * x), x);
    let s = fma(-n, n, n + n);
    let (ss, koff) = denormal_rescale!(s);
    let w = -ln_normal(ss, koff);
    let v = w.sqrt();
    let q = if w > ERFC_INV_W_FAR {
        erfinv_far_poly_m1(v - 7.0)
    } else {
        erfinv_tail_poly_m1(v - 1.0)
    };
    let mag = if x <= 0.7 { central } else { fma(v, q, v) };
    // n == 0 is the pole (+inf); n < 0 is NaN.
    let edge = if n == 0.0 { f32::INFINITY } else { f32::NAN };
    if n > 0.0 {
        mag
    } else {
        edge
    }
}

/// Inverse error function for `x` in `(-1, 1)`.
#[inline(always)]
pub fn erfinv(x: f32) -> f32 {
    let ax = x.abs();
    // Evaluates 1 - x^2 as (1-|x|)*(1+|x|) to avoid cancellation.
    let n = 1.0 - ax;
    let s = fma(-n, n, n + n);
    let w = -if s > 0.0 { ln_normal(s, 0.0) } else { f32::NAN };
    // Both arms apply sign at the end to preserve -0.0.
    let central = fma(ax, erfinv_central_poly_m1(x * x), ax);
    let v = w.sqrt();
    let tail = fma(v, erfinv_tail_poly_m1(v - 1.0), v);
    let normal = mulsign(if ax <= 0.7 { central } else { tail }, x);
    if ax == 1.0 {
        f32::INFINITY.copysign(x)
    } else {
        normal
    }
}

const NORM_CDF_XS_CLAMP: f32 = 14.44;
const _: () = assert!(
    (NORM_CDF_XS_CLAMP as f64) * (NORM_CDF_XS_CLAMP as f64) * 0.5 <= -(EXP_CLAMP_LO as f64)
);
const _: () =
    assert!((NORM_CDF_XS_CLAMP as f64) * (NORM_CDF_XS_CLAMP as f64) * 0.5 > 103.97207708399179);

/// Inverse complementary error function for `x` in `(0, 2)`.
#[inline(always)]
pub fn erfc_inv(y: f32) -> f32 {
    let n = if y < 1.0 { y } else { 2.0 - y };
    mulsign(erfc_inv_half(n), 1.0 - y)
}

/// Logit function: `ln(p / (1 - p))` for `p` in `(0, 1)`.
#[inline(always)]
pub fn logit(p: f32) -> f32 {
    let a = fma(2.0, p, -1.0);
    // logit(p) = 2*atanh(2p-1) near p=0.5 to avoid cancellation.
    let central = 2.0 * atanh_small(a);
    let np = -p;
    let t = 1.0 + np;
    let c = np - (t - 1.0);
    let corr = c / t;
    let spec = if t == 0.0 {
        f32::NEG_INFINITY
    } else {
        f32::NAN
    };
    let l = if t > 0.0 {
        ln_normal(t, 0.0) + corr
    } else {
        spec
    };
    let outer = ln(p) - l;
    if a.abs() < 0.25 {
        central
    } else {
        outer
    }
}

/// Computes `x * ln(y)`, with `0 * ln(y) = 0`.
#[inline(always)]
pub fn xlogy(x: f32, y: f32) -> f32 {
    let normal = x * ln(y);
    if x == 0.0 {
        0.0
    } else {
        normal
    }
}

/// Computes `x * ln(1 + y)`, with `0 * ln(1 + y) = 0`.
#[inline(always)]
pub fn xlog1py(x: f32, y: f32) -> f32 {
    let normal = x * log1p(y);
    if x == 0.0 {
        0.0
    } else {
        normal
    }
}

/// Computes `1 / sqrt(x)`.
#[inline(always)]
pub fn rsqrt(x: f32) -> f32 {
    1.0 / x.sqrt()
}

#[inline(always)]
fn hypot_checked(x: f32, y: f32) -> f32 {
    let ax = x.abs();
    let ay = y.abs();
    let m = ax.max(ay);
    let is_zero = ax == 0.0 && ay == 0.0;
    let m_safe = if is_zero { 1.0 } else { m };
    let tiny = m_safe < f32::MIN_POSITIVE;
    let pre = if tiny { 16777216.0 } else { 1.0 }; // 2^24
    let post = if tiny { 1.0 / 16777216.0 } else { 1.0 };
    let ax = ax * pre;
    let ay = ay * pre;
    let m_safe = m_safe * pre;
    let e = ((m_safe.to_bits() >> 23) as i32) - 127;
    let es = 2 * (e >> 1);
    let scale = f32::from_bits(((127 - es) as u32) << 23);
    let descale = f32::from_bits(((127 + es) as u32) << 23);
    let xs = ax * scale;
    let ys = ay * scale;
    let normal = fma(xs, xs, ys * ys).sqrt() * descale * post;
    let normal = if is_zero { 0.0 } else { normal };
    if x.is_infinite() || y.is_infinite() {
        f32::INFINITY
    } else {
        normal
    }
}

#[inline(always)]
fn rhypot(x: f32, y: f32) -> f32 {
    let normal = 1.0 / fma(x, x, y * y).sqrt();
    if x.is_infinite() || y.is_infinite() {
        0.0
    } else {
        normal
    }
}

/// Computes `sqrt(x^2 + y^2 + z^2)`.
#[inline(always)]
pub fn hypot3(x: f32, y: f32, z: f32) -> f32 {
    let normal = fma(x, x, fma(y, y, z * z)).sqrt();
    if x.is_infinite() || y.is_infinite() || z.is_infinite() {
        f32::INFINITY
    } else {
        normal
    }
}

/// Computes `1 / sqrt(x^2 + y^2 + z^2)`.
#[inline(always)]
pub fn rnorm3(x: f32, y: f32, z: f32) -> f32 {
    let normal = 1.0 / fma(x, x, fma(y, y, z * z)).sqrt();
    if x.is_infinite() || y.is_infinite() || z.is_infinite() {
        0.0
    } else {
        normal
    }
}

/// Normalizes a 2D vector `(x, y)`.
#[inline(always)]
pub fn normalize2(x: f32, y: f32) -> (f32, f32) {
    let r = rhypot(x, y);
    (x * r, y * r)
}

/// Normalizes a 3D vector `(x, y, z)`.
#[inline(always)]
pub fn normalize3(x: f32, y: f32, z: f32) -> (f32, f32, f32) {
    let r = rnorm3(x, y, z);
    (x * r, y * r, z * r)
}

/// Computes `sqrt(w^2 + x^2 + y^2 + z^2)`.
#[inline(always)]
pub fn hypot4(w: f32, x: f32, y: f32, z: f32) -> f32 {
    let normal = fma(w, w, fma(x, x, fma(y, y, z * z))).sqrt();
    let any_inf = w.is_infinite() || x.is_infinite() || y.is_infinite() || z.is_infinite();
    if any_inf {
        f32::INFINITY
    } else {
        normal
    }
}

/// Computes `1 / sqrt(w^2 + x^2 + y^2 + z^2)`.
#[inline(always)]
pub fn rnorm4(w: f32, x: f32, y: f32, z: f32) -> f32 {
    let normal = 1.0 / fma(w, w, fma(x, x, fma(y, y, z * z))).sqrt();
    let any_inf = w.is_infinite() || x.is_infinite() || y.is_infinite() || z.is_infinite();
    if any_inf {
        0.0
    } else {
        normal
    }
}

/// Normalizes a 4D quaternion/vector `(w, x, y, z)`.
#[inline(always)]
pub fn normalize4(w: f32, x: f32, y: f32, z: f32) -> (f32, f32, f32, f32) {
    let r = rnorm4(w, x, y, z);
    (w * r, x * r, y * r, z * r)
}

/// Computes `a*b - c*d` using Kahan's compensated algorithm.
#[inline(always)]
pub fn diff_of_products(a: f32, b: f32, c: f32, d: f32) -> f32 {
    let w = c * d;
    let e = fma(-c, d, w);
    let f = fma(a, b, -w);
    f + e
}

/// 2D cross product (`ax*by - ay*bx`) via Kahan's algorithm.
#[inline(always)]
pub fn cross2(ax: f32, ay: f32, bx: f32, by: f32) -> f32 {
    diff_of_products(ax, by, ay, bx)
}

/// Complex modulus `|re + im * i|`.
#[doc(alias = "cabsf")]
#[inline(always)]
pub fn cabs(re: f32, im: f32) -> f32 {
    hypot_checked(re, im)
}

/// Complex argument (principal value in `(-pi, pi]`).
#[doc(alias = "cargf")]
#[inline(always)]
pub fn carg(re: f32, im: f32) -> f32 {
    atan2(im, re)
}

/// Complex exponential `e^(re + im * i) = e^re * (cos(im) + i * sin(im))`.
#[inline(always)]
pub fn cexp(re: f32, im: f32) -> (f32, f32) {
    let m = exp(re);
    (m * cos(im), m * sin(im))
}

/// Complex natural log `ln(re+im*i) = ln(|re+im*i|) + i*arg(re+im*i)`, returned
/// as `(re, im)`.
#[inline(always)]
pub fn clog(re: f32, im: f32) -> (f32, f32) {
    // log1p for inputs where 1+v is known to be positive and normal.
    #[inline(always)]
    fn log1p_guarded(v: f32) -> f32 {
        let u = 1.0 + v;
        let c = v - (u - 1.0);
        ln_normal(u, 0.0) + c / u
    }
    let mag = cabs(re, im);
    let log_mag = if mag.is_finite() {
        if (mag - 1.0).abs() < 0.5 {
            let are = re.abs();
            let aim = im.abs();
            let a = are.max(aim) as f64;
            let b = are.min(aim) as f64;
            let v = f64::mul_add(a, a, -1.0) + b * b;
            0.5 * log1p_guarded(v as f32)
        } else {
            ln(mag)
        }
    } else if re.is_finite() && im.is_finite() {
        let are = re.abs();
        let aim = im.abs();
        let (mx, mn) = if are > aim { (are, aim) } else { (aim, are) };
        let ratio = mn / mx;
        ln(mx) + 0.5 * log1p_guarded(ratio * ratio)
    } else {
        ln(mag)
    };
    (log_mag, carg(re, im))
}

/// `(log2(1+s) - s*log2(e)) / s^2` over `s` in `[2^-0.5 - 1, 2^0.5 - 1]`.
const LOG2_Q_F64: [f64; 12] = [
    -7.213_475_205_084_016e-1,
    4.808_983_510_959_339e-1,
    -3.606_737_270_105_828e-1,
    2.885_381_586_685_297e-1,
    -2.404_514_367_623_340_5e-1,
    2.061_490_767_697_333_7e-1,
    -1.803_115_697_859_057_6e-1,
    1.591_062_250_458_101_8e-1,
    -1.433_471_501_632_151_1e-1,
    1.431_704_659_215_636_6e-1,
    -1.409_765_414_534_232_3e-1,
    7.781_246_573_228_973e-2,
];

/// `log2(x)` in f64, for positive finite `x` (denormals included; callers must
/// guard zero/negative/inf/nan themselves). The `powf` family's log half.
#[inline(always)]
fn log2_f64(x: f32) -> f64 {
    let (xs, koff) = denormal_rescale!(x);
    // Decompose x = 2^k * m with m in [2^-0.5, 2^0.5).
    let e = (xs.to_bits() as i32).wrapping_sub(0x3f3504f3) >> 23;
    let m = f32::from_bits((xs.to_bits() as i32).wrapping_sub(e << 23) as u32);
    let k = (e as f32 + koff) as f64;
    let s = (m as f64) - 1.0;
    let c = LOG2_Q_F64;
    let s2 = s * s;
    let s4 = s2 * s2;
    let s8 = s4 * s4;
    let p0 = f64::mul_add(c[1], s, c[0]);
    let p1 = f64::mul_add(c[3], s, c[2]);
    let p2 = f64::mul_add(c[5], s, c[4]);
    let p3 = f64::mul_add(c[7], s, c[6]);
    let p4 = f64::mul_add(c[9], s, c[8]);
    let p5 = f64::mul_add(c[11], s, c[10]);
    let q0 = f64::mul_add(p1, s2, p0);
    let q1 = f64::mul_add(p3, s2, p2);
    let q2 = f64::mul_add(p5, s2, p4);
    let r0 = f64::mul_add(q1, s4, q0);
    let qs = f64::mul_add(q2, s8, r0);
    let sq = s2 * qs;
    let lm = f64::mul_add(s, std::f64::consts::LOG2_E, sq);
    lm + k
}

/// `(2^f - 1)/f` over `f` in `[-0.5, 0.5]`.
const EXP2_F64_E: [f64; 6] = [
    0.6931472028549269,
    0.24022647913384074,
    0.05550332471225973,
    0.009618437395496837,
    0.0013398874430087457,
    0.0001535334944368378,
];

/// `2^v` for an f64 `v`, narrowed to f32. The `powf` family's exp half.
#[inline(always)]
fn exp2_f64_to_f32(v: f64) -> f32 {
    let vc = v.clamp(-200.0, 200.0);
    // Uses magic round to extract exponent field without saturating float-to-int cast.
    let nm = vc + ROUND_MAGIC64;
    let n = nm - ROUND_MAGIC64;
    let f = vc - n;
    let c = EXP2_F64_E;
    let f2 = f * f;
    let l0 = f64::mul_add(c[1], f, c[0]);
    let l1 = f64::mul_add(c[3], f, c[2]);
    let l2 = f64::mul_add(c[5], f, c[4]);
    let r0 = f64::mul_add(l1, f2, l0);
    let r1 = f64::mul_add(l2, f2 * f2, r0);
    let p = f64::mul_add(r1, f, 1.0);
    let scale = f64::from_bits(nm.to_bits().wrapping_add(1023) << 52);
    (p * scale) as f32
}

/// `exp2(log2(ax) * y)`, the magnitude half of the whole `powf` family.
macro_rules! powf_f64_mag {
    ($ax:expr, $y:expr) => {{
        let ax = $ax;
        let l = log2_f64(ax);
        let l = if ax == 0.0 { f64::NEG_INFINITY } else { l };
        let l = if !(ax < f32::INFINITY) {
            (ax * ax) as f64
        } else {
            l
        };
        exp2_f64_to_f32(l * ($y as f64))
    }};
}

macro_rules! powf_sign_combine {
    ($x:expr, $ax:expr, $y:expr, $mag:expr) => {{
        let x = $x;
        let y = $y;
        let ax = $ax;
        let mag = if ax == 1.0 { 1.0 } else { $mag };
        let par = fma(-2.0, (y * 0.5).floor(), y);
        let par = if x.is_sign_negative() { par } else { 0.0 };
        let spec = if ax + ax == ax { 1.0 } else { f32::NAN };
        let spec = if y.abs() == f32::INFINITY { 1.0 } else { spec };
        let sm = if par == 0.0 { 1.0 } else { spec };
        let sm = if par == 1.0 { -1.0 } else { sm };
        let r = mag * sm;
        if y == 0.0 {
            1.0
        } else {
            r
        }
    }};
}

/// Computes `x^y` (C99 `pow` semantics).
#[doc(alias = "pow")]
#[inline(always)]
#[allow(clippy::neg_cmp_op_on_partial_ord)] // `!(ax < inf)` catches NaN too
pub fn powf(x: f32, y: f32) -> f32 {
    let ax = x.abs();
    let mag = powf_f64_mag!(ax, y);
    powf_sign_combine!(x, ax, y, mag)
}

/// `powf` for `x > 0` (or `x == +0.0`).
#[inline(always)]
#[allow(clippy::neg_cmp_op_on_partial_ord)] // `!(x < inf)` catches NaN too
pub fn powf_pos(x: f32, y: f32) -> f32 {
    let mag = powf_f64_mag!(x, y);
    let r = if x == 1.0 { 1.0 } else { mag };
    if y == 0.0 {
        1.0
    } else {
        r
    }
}

/// Signed power: `copysign(|x|^y, x)`.
#[inline(always)]
pub fn signed_pow(x: f32, y: f32) -> f32 {
    mulsign(powf_pos(x.abs(), y), x)
}

/// `powf` without domain or sign checks.
#[inline(always)]
pub fn powf_unchecked(x: f32, y: f32) -> f32 {
    exp2_f64_to_f32(log2_f64(x) * y as f64)
}

// Exact constants for (c + 0.055) / 1.055 in sRGB conversion.
const SRGB_INV_1055: f32 = 0.947_867_3;
const SRGB_OFF_1055: f32 = 0.052_132_7;

/// Converts an sRGB color component in `[0, 1]` to linear (IEC 61966-2-1).
#[inline(always)]
#[allow(clippy::neg_cmp_op_on_partial_ord)]
pub fn srgb_to_linear(c: f32) -> f32 {
    let low = c * (1.0 / 12.92);
    let b = fma(c, SRGB_INV_1055, SRGB_OFF_1055);
    let p = exp2_checked(log_family_wrapper_discarded_unless_normal!(b, log_2_normal) * 0.4);
    let high = b * b * p;
    if c <= 0.04045 {
        low
    } else {
        high
    }
}

/// Converts a linear color component in `[0, 1]` to sRGB.
#[inline(always)]
#[allow(clippy::neg_cmp_op_on_partial_ord)]
pub fn linear_to_srgb(l: f32) -> f32 {
    let low = l * 12.92;
    let p =
        exp2_checked(log_family_wrapper_discarded_unless_normal!(l, log_2_normal) * (1.0 / 2.4));
    let high = fma(1.055, p, -0.055);
    if l <= 0.0031308 {
        low
    } else {
        high
    }
}

/// `x - round(x/y)*y` with ties away from zero (via `f32::round`, not IEEE 754 ties-to-even).
/// Only reliable for moderate `|x/y|`.
macro_rules! remainder_style_combine {
    ($x:expr, $y:expr, $q:expr) => {{
        let normal = fma(-$q, $y, $x);
        // Subtraction of equal values cancels to +0.0 in IEEE 754; restore sign of x.
        let normal = if normal == 0.0 {
            normal.copysign($x)
        } else {
            normal
        };
        let r = if $x == 0.0 && !normal.is_nan() {
            $x
        } else {
            normal
        };
        if $y.is_infinite() && $x.is_finite() {
            $x
        } else {
            r
        }
    }};
}

/// Truncated floating-point remainder `x - trunc(x/y) * y`.
#[doc(alias = "fmodf")]
#[inline(always)]
pub fn fmod(x: f32, y: f32) -> f32 {
    let q = (x / y).trunc();
    remainder_style_combine!(x, y, q)
}

/// Self-correcting `fmod` for `|x/y| <= 2^24`.
#[inline(always)]
pub fn fmod_checked(x: f32, y: f32) -> f32 {
    let q0 = (x / y).trunc();
    let r0 = fma(-q0, y, x);
    let wrong_sign = r0 != 0.0 && (r0 > 0.0) != (x > 0.0);
    let same_sign = (x > 0.0) == (y > 0.0);
    let adj = if wrong_sign != same_sign { 1.0 } else { -1.0 };
    let r1 = fma(-adj, y, r0);
    let needs_fix = wrong_sign || r0.abs() >= y.abs();
    let normal = if needs_fix { r1 } else { r0 };
    // Subtraction of equal values cancels to +0.0 in IEEE 754; restore sign of x.
    let normal = if normal == 0.0 {
        normal.copysign(x)
    } else {
        normal
    };
    let r = if x == 0.0 && !normal.is_nan() {
        x
    } else {
        normal
    };
    if y.is_infinite() && x.is_finite() {
        x
    } else {
        r
    }
}

/// `fmod` without domain checks.
#[inline(always)]
pub fn fmod_unchecked(x: f32, y: f32) -> f32 {
    let q = (x / y).trunc();
    fma(-q, y, x)
}

/// Computes the Euclidean remainder `x - y * floor(x / y)` (`0 <= result < |y|`).
#[inline(always)]
pub fn rem_euclid(x: f32, y: f32) -> f32 {
    let r = fmod(x, y);
    if r < 0.0 {
        r + y.abs()
    } else {
        r
    }
}

/// Computes the Euclidean quotient `floor(x / y)`.
#[inline(always)]
pub fn div_euclid(x: f32, y: f32) -> f32 {
    let q = (x / y).trunc();
    let r = fmod(x, y);
    if r < 0.0 {
        if y > 0.0 {
            q - 1.0
        } else {
            q + 1.0
        }
    } else {
        q
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        assert_eq!(log_2(1.0), 0.0);
        assert_eq!(log_2(2.0), 1.0);
        assert_eq!(log_2(4.0), 2.0);
        assert_eq!(log_2(8.0), 3.0);
    }

    fn ulp_error(
        range: std::ops::Range<i32>,
        scale: f64,
        f: impl Fn(f32) -> f32,
        reference: impl Fn(f64) -> f64,
    ) -> f32 {
        let count = (range.end - range.start) as f64;
        let err: u64 = range
            .map(|x| {
                let xf = x as f64 * scale;
                let ref_val = reference(xf) as f32;
                ref_val.to_bits().abs_diff(f(xf as f32).to_bits()) as u64
            })
            .sum();
        err as f32 / count as f32
    }

    #[test]
    fn cbrt_precision() {
        println!(
            "jodie cbrt error: {}",
            ulp_error(1..10000, 1.0, cbrt, |x| x.cbrt())
        );
        println!(
            "std   cbrt error: {}",
            ulp_error(1..10000, 1.0, |x| x.cbrt(), |x| x.cbrt())
        );
    }
    #[test]
    fn cbrt_accurate_precision() {
        println!(
            "jodie cbrt accurate error: {}",
            ulp_error(1..10000, 1.0, cbrt_accurate, |x| x.cbrt())
        );
        println!(
            "std   cbrt error: {}",
            ulp_error(1..10000, 1.0, |x| x.cbrt(), |x| x.cbrt())
        );
    }
    #[test]
    fn exp2_precision() {
        println!(
            "jodie exp2 error: {}",
            ulp_error(-100..100, 0.01, exp2, |x| x.exp2())
        );
        println!(
            "std   exp2 error: {}",
            ulp_error(-100..100, 0.01, |x| x.exp2(), |x| x.exp2())
        );
    }
    #[test]
    fn log2_precision() {
        println!(
            "jodie log2 error: {}",
            ulp_error(2..1000, 0.1, log_2, |x| x.log2())
        );
        println!(
            "std   log2 error: {}",
            ulp_error(2..1000, 0.1, |x| x.log2(), |x| x.log2())
        );
    }
    #[test]
    fn sin_precision() {
        println!(
            "jodie sin error: {}",
            ulp_error(-100..100, 0.01, sin, |x| x.sin())
        );
        println!(
            "std   sin error: {}",
            ulp_error(-100..100, 0.01, |x| x.sin(), |x| x.sin())
        );
    }
    #[test]
    fn cos_precision() {
        println!(
            "jodie cos error: {}",
            ulp_error(-100..100, 0.01, cos, |x| x.cos())
        );
        println!(
            "std   cos error: {}",
            ulp_error(-100..100, 0.01, |x| x.cos(), |x| x.cos())
        );
    }

    fn plot_approx(
        path: &str,
        x_start: f32,
        x_end: f32,
        approx: impl Fn(f32) -> f32,
        truth: impl Fn(f32) -> f32,
    ) {
        use plotters::prelude::*;
        let xs: Vec<f32> = (0..1000)
            .map(|i| x_start + (x_end - x_start) * i as f32 / 999.0)
            .collect();
        let all_y: Vec<f32> = xs.iter().flat_map(|&x| [approx(x), truth(x)]).collect();
        let y_min = all_y
            .iter()
            .cloned()
            .filter(|y| y.is_finite())
            .fold(f32::INFINITY, f32::min);
        let y_max = all_y
            .iter()
            .cloned()
            .filter(|y| y.is_finite())
            .fold(f32::NEG_INFINITY, f32::max);
        let root = BitMapBackend::new(path, (480, 480)).into_drawing_area();
        root.fill(&WHITE).unwrap();
        let mut chart = ChartBuilder::on(&root)
            .margin(5)
            .x_label_area_size(30)
            .y_label_area_size(30)
            .build_cartesian_2d(x_start..x_end, y_min..y_max)
            .unwrap();
        chart.configure_mesh().draw().unwrap();
        chart
            .draw_series(LineSeries::new(xs.iter().map(|&x| (x, approx(x))), &BLACK))
            .unwrap()
            .label("integer approx")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLACK));
        chart
            .draw_series(LineSeries::new(xs.iter().map(|&x| (x, truth(x))), &RED))
            .unwrap()
            .label("ground truth")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], RED));
        chart.configure_series_labels().draw().unwrap();
        root.present().expect("Unable to write result to file");
        std::process::Command::new("pngquant")
            .args(["--force", "--ext", ".png", "16", "--", path])
            .status()
            .unwrap();
    }

    fn plot_error(path: &str, x_start: f32, x_end: f32, f: impl Fn(f32) -> f32) {
        use plotters::prelude::*;
        let samples: Vec<(f32, f32)> = (0..1000)
            .map(|i| x_start + (x_end - x_start) * i as f32 / 999.0)
            .map(|x| (x, f(x)))
            .filter(|&(_, y)| y.is_finite())
            .collect();
        let y_min = samples
            .iter()
            .map(|&(_, y)| y)
            .fold(f32::INFINITY, f32::min);
        let y_max = samples
            .iter()
            .map(|&(_, y)| y)
            .fold(f32::NEG_INFINITY, f32::max);
        let root = BitMapBackend::new(path, (480, 480)).into_drawing_area();
        root.fill(&WHITE).unwrap();
        let mut chart = ChartBuilder::on(&root)
            .margin(5)
            .x_label_area_size(30)
            .y_label_area_size(50)
            .build_cartesian_2d(x_start..x_end, y_min..y_max)
            .unwrap();
        chart
            .configure_mesh()
            .y_label_formatter(&|y| format!("{:.2e}", y))
            .draw()
            .unwrap();
        chart.draw_series(LineSeries::new(samples, &BLACK)).unwrap();
        root.present().expect("Unable to write result to file");
        std::process::Command::new("pngquant")
            .args(["--force", "--ext", ".png", "16", "--", path])
            .status()
            .unwrap();
    }

    #[test]
    fn cbrt_approx_plot() {
        plot_approx("cbrt_approx.png", 1., 128., cbrt_approx, |x| x.cbrt());
    }
    #[test]
    fn sqrt_approx_plot() {
        plot_approx("sqrt_approx.png", 1., 128., sqrt_approx, |x| x.sqrt());
    }
    #[test]
    fn rcp_approx_plot() {
        plot_approx("rcp_approx.png", 1., 10., rcp_approx, |x| 1.0 / x);
    }
    #[test]
    fn exp2_approx_plot() {
        plot_approx("exp2_approx.png", 0., 10., exp2_approx, |x| x.exp2());
    }
    #[test]
    fn log2_approx_plot() {
        plot_approx("log2_approx.png", 1., 128., log2_approx, |x| x.log2());
    }
    #[test]
    fn sin_plot() {
        plot_approx("sin.png", -20., 20., sin, |x| x.sin());
    }
    #[test]
    fn cos_plot() {
        plot_approx("cos.png", -20., 20., cos, |x| x.cos());
    }
    #[test]
    fn rsqrt_approx_plot() {
        plot_approx("rsqrt_approx.png", 1., 128., rsqrt_approx, |x| {
            1.0 / x.sqrt()
        });
    }

    #[test]
    fn log_2_error() {
        plot_error("log_2_error.png", 1., 128., |x| {
            log_2(x) / (x as f64).log2() as f32 - 1.0
        });
    }
    #[test]
    fn exp2_error() {
        plot_error("exp2_error.png", 0., 10., |x| {
            exp2(x) / (x as f64).exp2() as f32 - 1.0
        });
    }
    #[test]
    fn sin_error() {
        plot_error("sin_error.png", -20., 20., |x| {
            sin(x) / (x as f64).sin() as f32 - 1.0
        });
    }
    #[test]
    fn cos_error() {
        plot_error("cos_error.png", -20., 20., |x| {
            cos(x) / (x as f64).cos() as f32 - 1.0
        });
    }
    #[test]
    fn cbrt_error() {
        plot_error("cbrt_error.png", 1., 128., |x| {
            cbrt(x) / (x as f64).cbrt() as f32 - 1.0
        });
    }
    #[test]
    fn cbrt_accurate_error() {
        plot_error("cbrt_accurate_error.png", 1., 128., |x| {
            cbrt_accurate(x) / (x as f64).cbrt() as f32 - 1.0
        });
    }
    #[test]
    fn cbrt_approx_error() {
        plot_error("cbrt_approx_error.png", 1., 128., |x| {
            cbrt_approx(x) / (x as f64).cbrt() as f32 - 1.0
        });
    }

    #[test]
    fn wide_trig_special_values() {
        // Zero and signed zero
        assert_eq!(sin_wide(0.0).to_bits(), 0.0f32.to_bits());
        assert_eq!(sin_wide(-0.0).to_bits(), (-0.0f32).to_bits());
        assert_eq!(cos_wide(0.0), 1.0);
        assert_eq!(cos_wide(-0.0), 1.0);
        assert_eq!(tan_wide(0.0).to_bits(), 0.0f32.to_bits());
        assert_eq!(tan_wide(-0.0).to_bits(), (-0.0f32).to_bits());

        // Infinities and NaNs
        assert!(sin_wide(f32::INFINITY).is_nan());
        assert!(sin_wide(f32::NEG_INFINITY).is_nan());
        assert!(sin_wide(f32::NAN).is_nan());
        assert!(cos_wide(f32::INFINITY).is_nan());
        assert!(cos_wide(f32::NEG_INFINITY).is_nan());
        assert!(cos_wide(f32::NAN).is_nan());
        assert!(tan_wide(f32::INFINITY).is_nan());
        assert!(tan_wide(f32::NEG_INFINITY).is_nan());
        assert!(tan_wide(f32::NAN).is_nan());

        // Subnormals (e < pitable::CUT)
        let subnorm = f32::from_bits(0x0000_0001); // min positive subnormal
        assert_eq!(sin_wide(subnorm), subnorm);
        assert_eq!(cos_wide(subnorm), 1.0);
        assert_eq!(tan_wide(subnorm), subnorm);
    }

    #[test]
    fn wide_trig_pythagorean_identity() {
        let mut rng = 987654321u64;
        for _ in 0..10_000 {
            rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1);
            let b = (rng >> 32) as u32;
            let x = f32::from_bits(b);
            if !x.is_finite() {
                continue;
            }
            let s = sin_wide(x) as f64;
            let c = cos_wide(x) as f64;
            let pyth = s * s + c * c;
            assert!(
                (pyth - 1.0).abs() < 1e-4,
                "Pythagorean identity failed for x={x:e}: sin={s}, cos={c}, sum={pyth}"
            );
        }
    }
}
