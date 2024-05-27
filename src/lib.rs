const SIGN_MASK: u32 = 0x80000000;
const EXPONENT_MASK: u32 = 0x7f800000;
const MANTISSA_MASK: u32 = 0x007fffff;

fn fmaf(a: f32, b: f32, c: f32) -> f32 {
    return a.mul_add(b, c);
}

fn log_2_mantissa(x: f32) -> f32 {
    return fmaf(
        fmaf(f32::from_bits(0x40153ebb), x, f32::from_bits(0x413b8af9)),
        x,
        f32::from_bits(0x409c1a68),
    ) * (x - 1.)
        / fmaf(
            fmaf(
                fmaf(f32::from_bits(0x3ecfca47), x, f32::from_bits(0x409f8156)),
                x,
                f32::from_bits(0x40d76ca4),
            ),
            x,
            1.,
        );
}

pub fn log_2(x: f32) -> f32 {
    // log2(x*y) == log2(x)+log2(y)
    let mantissa = f32::from_bits(1_f32.to_bits() | (x.to_bits() & MANTISSA_MASK));
    let log2exponent = f32::from_bits(256_f32.to_bits() | ((x.to_bits() & EXPONENT_MASK) >> 8)) - 383.;
    return log2exponent + log_2_mantissa(mantissa);
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
}
