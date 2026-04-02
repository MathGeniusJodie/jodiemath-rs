use std::ops::{Add, Div, Mul, Neg, Sub};

/// Double-float: an unevaluated sum (.0, .1) where |.1| <= ulp(.0)/2.
/// Provides ~2× f32 precision (~14 decimal digits) using pairs of f32s.
#[derive(Clone, Copy, Debug)]
pub struct Df32(pub f32, pub f32);

// ── Helpers ──────────────────────────────────────────────────────────

#[inline(always)]
fn fma(a: f32, b: f32, c: f32) -> f32 {
    a.mul_add(b, c)
}

/// Error-free addition (Knuth two-sum).
#[inline(always)]
fn two_sum(a: f32, b: f32) -> (f32, f32) {
    let s = a + b;
    let v = s - a;
    let e = (a - (s - v)) + (b - v);
    (s, e)
}

/// Fast path when |a| >= |b| is guaranteed.
#[inline(always)]
fn quick_two_sum(a: f32, b: f32) -> (f32, f32) {
    let s = a + b;
    let e = b - (s - a);
    (s, e)
}

// ── Constructors ─────────────────────────────────────────────────────

impl Df32 {
    /// Exact conversion from a single f32.
    #[inline(always)]
    pub const fn from_f32(v: f32) -> Self {
        Self(v, 0.0)
    }

    /// Build from the product of two f32s (error-free).
    #[inline(always)]
    pub fn from_mul(a: f32, b: f32) -> Self {
        let p = a * b;
        let e = fma(a, b, -p);
        Self(p, e)
    }

    /// Build from the sum of two f32s (error-free).
    #[inline(always)]
    pub fn from_add(a: f32, b: f32) -> Self {
        let (s, e) = two_sum(a, b);
        Self(s, e)
    }

    #[inline(always)]
    pub fn from_quick_add(a: f32, b: f32) -> Self {
        let (s, e) = quick_two_sum(a, b);
        Self(s, e)
    }

    /// Collapse back to a single f32.
    #[inline(always)]
    pub fn to_f32(self) -> f32 {
        self.0 + self.1
    }

    /// Fast division returning a single f32 (the original div_ff_ff_f32).
    #[inline(always)]
    pub fn div_to_f32(self, rhs: Self) -> f32 {
        let rcp = 1.0 / rhs.0;
        let q1 = self.0 * rcp;
        let rh = fma(-q1, rhs.0, self.0) + fma(-q1, rhs.1, self.1);
        fma(rh, rcp, q1)
    }

    #[inline(always)]
    pub fn quick_add(self, b: f32) -> Self {
        let s = self.0 + b;
        let lo = self.1 + (b - (s - self.0));
        Self(s, lo)
    }
    #[inline(always)]
    pub fn square(self) -> Self {
        let p = self.0 * self.0;
        let e = fma(self.0, self.0, -p);
        //let lo = fma(self.1, self.0, e) + fma(self.0, self.1, self.1 * self.1);
        let lo = fma(self.0, self.1 * 2., e);
        Self(p, lo)
    }
    #[inline(always)]
    pub fn quick_add_df(self, rhs: Self) -> Self {
        let (s, e) = quick_two_sum(self.0, rhs.0);
        //let (s, e) = quick_two_sum(s, e + self.1 + rhs.1);
        //Self(s,e)
        Self(s, e + self.1 + rhs.1)
    }
    #[inline(always)]
    pub fn abs(self) -> Self {
        if self.0 > 0.{
            self
        } else {
            Self(-self.0,-self.1)
        }
    }
}

impl From<Df32> for f32 {
    #[inline(always)]
    fn from(d: Df32) -> f32 {
        d.0 + d.1
    }
}

impl From<f32> for Df32 {
    #[inline(always)]
    fn from(v: f32) -> Self {
        Self::from_f32(v)
    }
}

// ── Addition ─────────────────────────────────────────────────────────

impl Add for Df32 {
    type Output = Self;
    #[inline(always)]
    fn add(self, rhs: Self) -> Self {
        let (s, e) = two_sum(self.0, rhs.0);
        let e = e + self.1 + rhs.1;
        let (s, e) = quick_two_sum(s, e);
        Self(s, e)
    }
}

impl Add<f32> for Df32 {
    type Output = Self;
    #[inline(always)]
    fn add(self, rhs: f32) -> Self {
        let (s, e) = two_sum(self.0, rhs);
        let e = e + self.1;
        let (s, e) = quick_two_sum(s, e);
        Self(s, e)
    }
}

impl Add<Df32> for f32 {
    type Output = Df32;
    #[inline(always)]
    fn add(self, rhs: Df32) -> Df32 {
        rhs + self
    }
}

// ── Subtraction ──────────────────────────────────────────────────────

impl Neg for Df32 {
    type Output = Self;
    #[inline(always)]
    fn neg(self) -> Self {
        Self(-self.0, -self.1)
    }
}

impl Sub for Df32 {
    type Output = Self;
    #[inline(always)]
    fn sub(self, rhs: Self) -> Self {
        self + (-rhs)
    }
}

impl Sub<f32> for Df32 {
    type Output = Self;
    #[inline(always)]
    fn sub(self, rhs: f32) -> Self {
        self + (-rhs)
    }
}

impl Sub<Df32> for f32 {
    type Output = Df32;
    #[inline(always)]
    fn sub(self, rhs: Df32) -> Df32 {
        (-rhs) + self
    }
}

// ── Multiplication ───────────────────────────────────────────────────

impl Mul for Df32 {
    type Output = Self;
    #[inline(always)]
    fn mul(self, rhs: Self) -> Self {
        let p = self.0 * rhs.0;
        let e = fma(self.0, rhs.0, -p);
        let lo = fma(self.1, rhs.0, e) + fma(self.0, rhs.1, rhs.1 * self.1);
        //let lo = fma(self.1, rhs.0, fma(self.0, rhs.1, e));
        Self(p, lo)
    }
}

impl Mul<f32> for Df32 {
    type Output = Self;
    #[inline(always)]
    fn mul(self, rhs: f32) -> Self {
        let p = self.0 * rhs;
        let e = fma(self.0, rhs, -p);
        let lo = fma(self.1, rhs, e);
        Self(p, lo)
    }
}

impl Mul<Df32> for f32 {
    type Output = Df32;
    #[inline(always)]
    fn mul(self, rhs: Df32) -> Df32 {
        rhs * self
    }
}

// ── Division ─────────────────────────────────────────────────────────

impl Div for Df32 {
    type Output = Self;
    #[inline(always)]
    fn div(self, rhs: Self) -> Self {
        let rcp = 1.0 / rhs.0;
        let q1 = self.0 * rcp;
        let rh = fma(-q1, rhs.0, self.0) + fma(-q1, rhs.1, self.1);
        let q2 = rh * rcp;
        let (hi, lo) = quick_two_sum(q1, q2);
        Self(hi, lo)
    }
}

impl Div<f32> for Df32 {
    type Output = Self;
    #[inline(always)]
    fn div(self, rhs: f32) -> Self {
        self / Df32::from_f32(rhs)
    }
}

impl Div<Df32> for f32 {
    type Output = Df32;
    #[inline(always)]
    fn div(self, rhs: Df32) -> Df32 {
        Df32::from_f32(self) / rhs
    }
}

// ── Display ──────────────────────────────────────────────────────────

impl std::fmt::Display for Df32 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let val = self.0 as f64 + self.1 as f64;
        write!(f, "{val}")
    }
}

// ── Tests ────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tuple_access() {
        let x = Df32(3.14, 1e-10);
        assert_eq!(x.0, 3.14);
        assert_eq!(x.1, 1e-10);
    }

    #[test]
    fn test_destructure() {
        let Df32(hi, lo) = Df32::from_mul(3.0, 7.0);
        assert_eq!(hi, 3.0 * 7.0);
        assert!(lo.abs() < 1e-6);
    }

    #[test]
    fn test_sum_precision() {
        let a = Df32::from_f32(1.0);
        let tiny = Df32::from_f32(1e-10);
        let sum = a + tiny;
        let result = sum.0 as f64 + sum.1 as f64;
        let expected = 1.0_f64 + 1e-10_f64;
        assert!(
            (result - expected).abs() < 1e-14,
            "got {result}, expected {expected}"
        );
    }

    #[test]
    fn test_cross_type_ops() {
        let d = Df32::from_f32(3.0);
        let _ = d * 2.0_f32;
        let _ = 2.0_f32 * d;
        let _ = d + 1.0_f32;
        let _ = 1.0_f32 + d;
        let _ = d / 2.0_f32;
        let _ = Df32::from_mul(2.0, 3.0);
    }
}
