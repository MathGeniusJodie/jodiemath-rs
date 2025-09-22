use std::f32::consts::E;

use plotters::data::Quartiles;

const SIGN_MASK: u32 = 0x80000000;
const EXPONENT_MASK: u32 = 0x7f800000;
const MANTISSA_MASK: u32 = 0x007fffff;

fn fma(a: f32, b: f32, c: f32) -> f32 {
    return a.mul_add(b, c);
}

fn mulsign(x: f32, y: f32) -> f32 {
    return f32::from_bits(x.to_bits() ^ (y.to_bits() & SIGN_MASK));
}

fn log_2_mantissa(x: f32) -> f32 {
    let a = f32::from_bits(0x40153ebb);
    let b = f32::from_bits(0x413b8af9);
    let c = f32::from_bits(0x409c1a68);
    let d = f32::from_bits(0x3ecfca47);
    let e = f32::from_bits(0x409f8156);
    let f = f32::from_bits(0x40d76ca4);
    return fma(fma(a, x, b), x, c) * (x - 1.) / fma(fma(fma(d, x, e), x, f), x, 1.);
}

pub fn log_2(x: f32) -> f32 {
    // log2(x*y) == log2(x)+log2(y)
    let mantissa = f32::from_bits(1_f32.to_bits() | (x.to_bits() & MANTISSA_MASK));
    let log2exponent =
        f32::from_bits(256_f32.to_bits() | ((x.to_bits() & EXPONENT_MASK) >> 8)) - 383.;
    return log2exponent + log_2_mantissa(mantissa);
}

fn exp2_fract(x: f32) -> f32 {
    let a = f32::from_bits(0x3af71c15);
    let b = f32::from_bits(0x3c130514);
    let c = f32::from_bits(0x3d64b437);
    let d = f32::from_bits(0x3e75ea9e);
    let e = f32::from_bits(0x3f317271);
    fma(fma(fma(fma(fma(a, x, b), x, c), x, d), x, e), x, 1_f32)
}

pub fn exp2(x: f32) -> f32 {
    // exp2(floor(x))*exp2(fract(x)) == exp2(x)
    let exp2int = f32::from_bits(((x + 383_f32).to_bits() << 8) & EXPONENT_MASK);
    let fract = x - x.floor();
    exp2int * exp2_fract(fract)
}

fn sinf_poly(x: f32) -> f32 {
    let a = f32::from_bits(0xb2cc0ff1);
    let b = f32::from_bits(0x3638a80e);
    let c = f32::from_bits(0xb9500b44);
    let d = f32::from_bits(0x3c088883);
    let e = f32::from_bits(0xbe2aaaaa);
    let x2 = x * x;
    fma(fma(fma(fma(fma(a, x2, b), x2, c), x2, d), x2, e), x2, 1.) * x
}

pub fn sin(x: f32) -> f32 {
    let tau: f64 = 6.2831853071795864769252867665590057683943387987502116419498891846;
    let tauh: f32 = tau as f32;
    let taul: f32 = (tau - (tauh as f64)) as f32;
    let rtau: f64 = 0.1591549430918953357688837633725143620344596457404564487476673440;
    let rtauh: f32 = rtau as f32;
    let rtaul: f32 = (rtau - (rtauh as f64)) as f32;
    let pi = 3.14159265358979323846;
    let pilo = f32::from_bits(0x33bbbd2e);
    let rpi = 0.3183098861837907;
    let z = x - (x * rtauh).round() * tauh + (x * rtauh).round() * taul;
    let y = (x - (x * rpi).round() * pi + (x * rpi).round() * pilo).abs();
    return sinf_poly(mulsign(y, z));
}

// todo: make less convoluted
pub fn cos(x: f32) -> f32 {
    let tau: f64 = 6.2831853071795864769252867665590057683943387987502116419498891846;
    let tauh: f32 = tau as f32;
    let taul: f32 = (tau - (tauh as f64)) as f32;
    let rtau: f64 = 0.1591549430918953357688837633725143620344596457404564487476673440;
    let rtauh: f32 = rtau as f32;
    let rtaul: f32 = (rtau - (rtauh as f64)) as f32;

    let m = (fma(x, (rtauh * 2f32), x * (rtaul * 2f32))).floor();
    let s = (((m as i32) & 1) * 2 - 1) as f32;

    let high = fma(m, tauh / 2f32, tauh / 4f32);
    let errorhigh = tauh / 4f32 + fma(m, tauh / 2f32, -high);
    let low = fma(m, taul / 2f32, taul / 4f32);

    let mut x = x;
    x -= high;
    x -= errorhigh;
    x -= low;

    sinf_poly(x * s)
}

//f0 = (s^3-x)
//f1 = (3 s^2)
//f2 = (6s)
//f3 = 6
//
// (6.0*f0*f1*f1-3.0*f0*f0*f2) / (6.0*f1*f1*f1-6.0*f0*f1*f2+f0*f0*f3)
// (3 s (2 s^3 + x) (s^3 - x))/(10 s^6 + 16 s^3 x + x^2)

pub fn cbrt(x: f32) -> f32 {
    let s = f32::from_bits(x.to_bits() / 3 + 709982100);
    let s3 = s * s * s;
    fma(
        fma(6., s3, 3. * x),
        s * (x - s3) / fma(fma(10., s3, 16. * x), s3, x * x),
        s,
    )
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
    #[test]
    fn cbrt_precision() {
        let mut err: u64 = 0;
        for x in 1..1000 {
            let reference = (x as f64).cbrt() as f32;
            let result = cbrt(x as f32);
            err += reference.to_bits().abs_diff(result.to_bits()) as u64;
        }
        println!("total error: {}", (err as f32) / 1000.);
    }
    #[test]
    fn exp2_precision() {
        let mut err = 0;
        for x in -100..100 {
            let reference = ((x as f32) * 0.1).exp2();
            let result = exp2((x as f32) * 0.1);
            err += reference.to_bits().abs_diff(result.to_bits());
        }
        println!("total error: {}", (err as f32) / 200.);
    }
    #[test]
    fn log2_precision() {
        let mut err = 0;
        for x in 2..1000 {
            let reference = (x as f64).log2() as f32;
            let result = log_2(x as f32);
            err += reference.to_bits().abs_diff(result.to_bits());
        }
        println!("total error: {}", (err as f32) / 1000.);
    }
    #[test]
    fn sin_precision() {
        let mut err = 0;
        for x in -1000..1000 {
            let reference = ((x as f32) * 0.1).sin();
            let result = sin((x as f32) * 0.1);
            err += reference.to_bits().abs_diff(result.to_bits());
        }
        println!("total error: {}", (err as f32) / 2000.);
    }
    #[test]
    fn cos_precision() {
        let mut err = 0;
        for x in -1000..1000 {
            let reference = ((x as f32) * 0.1f32).cos() as f32;
            let result = cos((x as f32) * 0.1);
            err += reference.to_bits().abs_diff(result.to_bits()) as u64;
        }
        println!("total error: {}", (err as f32) / 2000.);
    }
    #[test]
    fn cbrt_plot() {
        use plotters::prelude::*;
        let root = BitMapBackend::new("cbrt.png", (640, 480)).into_drawing_area();
        root.fill(&WHITE).unwrap();
        let mut chart = ChartBuilder::on(&root)
            .caption("Cube Root", ("sans-serif", 50).into_font())
            .margin(5)
            .x_label_area_size(30)
            .y_label_area_size(30)
            .build_cartesian_2d(0f32..128f32, -0.1f32..0.1f32)
            .unwrap();

        chart.configure_mesh().draw().unwrap();

        chart
            .draw_series(LineSeries::new(
                (0..1000)
                    .map(|x| x as f32)
                    .map(|x| (x, cbrt(x as f32) / (x as f32).cbrt() - 1.0)),
                &RED,
            ))
            .unwrap()
            .label("std")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));
    }

    // plot cos and builtin cos in 2 different colors
    #[test]
    fn cos_plot() {
        use plotters::prelude::*;
        let root = BitMapBackend::new("cos.png", (640, 480)).into_drawing_area();
        root.fill(&WHITE).unwrap();
        let mut chart = ChartBuilder::on(&root)
            .caption("Cosine", ("sans-serif", 50).into_font())
            .margin(5)
            .x_label_area_size(30)
            .y_label_area_size(30)
            .build_cartesian_2d(-20f32..20f32, -1.1f32..1.1f32)
            .unwrap();
        chart.configure_mesh().draw().unwrap();
        chart
            .draw_series(LineSeries::new(
                (-200..200).map(|x| x as f32 * 0.1).map(|x| (x, cos(x))),
                &RED,
            ))
            .unwrap()
            .label("approx")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));
        chart
            .draw_series(LineSeries::new(
                (-200..200).map(|x| x as f32 * 0.1).map(|x| (x, x.cos())),
                &BLUE,
            ))
            .unwrap()
            .label("std")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));
        chart.configure_series_labels().draw().unwrap();
        chart
            .draw_series(LineSeries::new(
                (-200..200).map(|x| x as f32 * 0.1).map(|x| (x, x.cos())),
                &BLUE,
            ))
            .unwrap()
            .label("std")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));
        chart.configure_series_labels().draw().unwrap();

        root.present().expect("Unable to write result to file, please make sure 'plotters' crate is in your Cargo.toml");
    }
}
