use std::f32::consts::E;

use plotters::data::Quartiles;

const SIGN_MASK: u32 = 0x80000000;
const EXPONENT_MASK: u32 = 0x7f800000;
const MANTISSA_MASK: u32 = 0x007fffff;

fn fmaf(a: f32, b: f32, c: f32) -> f32 {
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
    return fmaf(fmaf(a, x, b), x, c) * (x - 1.) / fmaf(fmaf(fmaf(d, x, e), x, f), x, 1.);
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
    fmaf(fmaf(fmaf(fmaf(fmaf(a, x, b), x, c), x, d), x, e), x, 1_f32)
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
    fmaf(
        fmaf(fmaf(fmaf(fmaf(a, x2, b), x2, c), x2, d), x2, e),
        x2,
        1.,
    ) * x
}

pub fn sin(x: f32) -> f32 {
    let tau = 6.28318530717958647692;
    let taulo = f32::from_bits(0x343bbd2e);
    let rtau = 0.15915494309189533576;
    let pi = 3.14159265358979323846;
    let pilo = f32::from_bits(0x33bbbd2e);
    let rpi = 0.3183098861837907;
    let z = x - (x * rtau).round() * tau + (x * rtau).round() * taulo;
    let y = (x - (x * rpi).round() * pi + (x * rpi).round() * pilo).abs();
    return sinf_poly(mulsign(y, z));
}

// todo: make less convoluted
pub fn cos(x: f32) -> f32 {
    let tau: f32 = f32::from_bits(0x40C90FDB);
    let taulo: f32 = f32::from_bits(0x50824525);
    let rtau: f32 = f32::from_bits(0x3E22F983);
    let rtaulo: f32 = f32::from_bits(0x31DC9C88);
    let pi = tau / 2.;
    let rpi = rtau * 2.;
    let rpilo = rtaulo * 2.;
    let hpi = tau / 4.;
    let hpilo = taulo / 4.;

    let m = fmaf(x, rpi, x * rpilo).floor();
    let s = (((m as u32) & 1) as f32) * 2. - 1.;

    let high = fmaf(m, pi, hpi);
    let errorhigh = hpi + fmaf(m, pi, -high);
    let low = fmaf(m, pi, hpilo);

    let mut x = x;
    x -= high;
    x -= errorhigh;
    x -= low;

    sinf_poly(x * s)
}

pub fn cbrt(x: f32) -> f32 {
    let mut s = f32::from_bits(
        unsafe{fmaf(
            x.to_bits() as f32,
            0.33333333,
            709983100.
        ).to_int_unchecked()}
    );
    let s3 = s * s * s;
    s -= (3.0*s*fmaf(2.0,s3,x)*(s3-x))
        / fmaf(fmaf(10.0,s3,16.0*x),s3,x*x);
    return s;
}




//f0 = (s^3-x)
//f1 = (3 s^2)
//f2 = (6s)
//f3 = 6
//
// (6.0*f0*f1*f1-3.0*f0*f0*f2) / (6.0*f1*f1*f1-6.0*f0*f1*f2+f0*f0*f3)
// (3 s (2 s^3 + x) (s^3 - x))/(10 s^6 + 16 s^3 x + x^2)





/*
pub fn cbrt(x: f32) -> f32 {

    let third_x = x * -0.33333333;

    let mut y = 
        f32::from_bits(
        unsafe{
        fmaf(x.to_bits() as f32 , -0.33333333, 1418472267.) // initial guess for inverse cube root
        .to_int_unchecked()});

    let k1 = 1.752319676;
    let k2 = -1.2509524245;
    let k3 = 0.5093818292;
    let a = x*(y*y);
    y = fmaf((y*y)*a,fmaf(a,k3*y,k2),k1*y);

    y = fmaf((y*y)*(y*y), third_x, 1.3333333*y);
    x * y * y
}*/

/*
1:  float InvCbrt21 (float x){

5:     int i = *(int*) &x;
6:     i = 0x548c2b4b − i/3;
7:     float y = *(float*) &i;
8:     float c = x*y*y*y;
9:     y = y*(k1 − c*(k2 − k3*c));
10:      c = 1.0f − x*y*y*y;//fmaf
11:      y = y*(1.0f + 0.333333333333f*c);//fmaf
12:      return y;
13:  }
 */


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
        let mut err = 0.;
        for x in 1..1000 {
            let reference = (x as f32).cbrt();
            let result = cbrt(x as f32);
            err += (result/reference - 1.).abs();
        }
        println!("total error: {}", err/1000.);
    }
    #[test]
    fn exp2_precision() {
        let mut err = 0.;
        for x in -100..100 {
            let reference = ((x as f32)*0.1).exp2();
            let result = exp2((x as f32)*0.1);
            err += (result/reference - 1.).abs();
        }
        println!("total error: {}", err/200.);
    }
    #[test]
    fn log2_precision() {
        let mut err = 0.;
        for x in 2..1000 {
            let reference = (x as f32).log2();
            let result = log_2(x as f32);
            err += (result/reference - 1.).abs();
        }
        println!("total error: {}", err/1000.);
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
            .build_cartesian_2d(0f32..128f32, -0.001f32..0.001f32)
            .unwrap();

        chart.configure_mesh().draw().unwrap();

        chart
            .draw_series(LineSeries::new(
                (0..1000).map(|x| x as f32).map(|x| (x, cbrt(x as f32)-1.0)),
                &RED,
            ))
            .unwrap()
            .label("std")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));
    }
}
