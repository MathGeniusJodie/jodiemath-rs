// godbolt flags -C opt-level=3 -C target_feature=+fma

mod doublefloat;
use doublefloat::Df32;

const SIGN_MASK: u32 = 0x80000000;
const EXPONENT_MASK: u32 = 0x7f800000;
const MANTISSA_MASK: u32 = 0x007fffff;
use std::f32::consts::TAU as TAU32;
use std::f64::consts::TAU as TAU64;
const TAUDF: Df32 = Df32(TAU32,(TAU64 - (TAU32 as f64)) as f32);
const RTAUDF:Df32 = Df32(
    (1./TAU64) as f32,
    ((1./TAU64) - (((1./TAU64) as f32) as f64)) as f32
);
const RPIDF:Df32 = Df32(RTAUDF.0*2.,RTAUDF.1*2.);
const HPIDF: Df32 = Df32(TAUDF.0/4.,TAUDF.1/4.);
const PIDF: Df32 = Df32(TAUDF.0/2.,TAUDF.1/2.);

#[inline(always)]
fn fma(a: f32, b: f32, c: f32) -> f32 {
    a.mul_add(b, c)
}
#[inline(always)]
fn mulsign(x: f32, y: f32) -> f32 {
    f32::from_bits(x.to_bits() ^ (y.to_bits() & SIGN_MASK))
}

#[inline(always)]
pub fn log_2(x: f32) -> f32 {
    let a = f32::from_bits(0x40153ebb);
    let b = f32::from_bits(0x41163b4a);
    let c = f32::from_bits(0xc09c1a68);
    let d = f32::from_bits(0x3ecfca47);
    let e = f32::from_bits(0x409f8156);
    let f = f32::from_bits(0x40d76ca4);
    let g = f32::from_bits(0xc0dafb8a);
    // log2(x*y) == log2(x)+log2(y)
    let m = f32::from_bits(1_f32.to_bits() | (x.to_bits() & MANTISSA_MASK));
    let log2exponent =
        f32::from_bits(256_f32.to_bits() | ((x.to_bits() & EXPONENT_MASK) >> 8)) - 383.;
    log2exponent + fma(m * m, fma(a, m, b), fma(g, m, c)) / fma(m * m, fma(d, m, e), fma(f, m, 1.))
}

#[inline(always)]
pub fn exp2(x: f32) -> f32 {
    let a = f32::from_bits(0x3af71c15);
    let b = f32::from_bits(0x3c130514);
    let c = f32::from_bits(0x3d64b437);
    let d = f32::from_bits(0x3e75ea9e);
    let e = f32::from_bits(0x3f317271);
    // exp2(floor(x))*exp2(fract(x)) == exp2(x)
    let exp2int = f32::from_bits(((x + 383_f32).to_bits() << 8) & EXPONENT_MASK);
    let f = x - x.floor();
    exp2int * fma(fma(fma(fma(fma(a, f, b), f, c), f, d), f, e), f, 1_f32)
}
#[inline(always)]
fn sinf_poly(x: f32) -> f32 {
    let a = f32::from_bits(0xb2cc0ff1);
    let b = f32::from_bits(0x3638a80e);
    let c = f32::from_bits(0xb9500b44);
    let d = f32::from_bits(0x3c088883);
    let e = f32::from_bits(0xbe2aaaaa);
    let x2 = x * x;
    fma(fma(fma(fma(fma(a, x2, b), x2, c), x2, d), x2, e), x2, 1.) * x
}
#[inline(always)]
pub fn sin(x: f32) -> f32 {
    let q = fma(x,RTAUDF.0,0.25).round();
    let high = fma(-q, TAUDF.0, HPIDF.0);
    let errorhigh = HPIDF.0 + fma(-q, TAUDF.0, -high);
    let low = fma(-q, TAUDF.1, HPIDF.1);
    let y = x + high;
    let z = (-HPIDF).quick_add_to_f32(Df32(y,errorhigh+low).abs());
    sinf_poly(z)
}
#[inline(always)]
pub fn cos(x: f32) -> f32 {
    let q = (x*RTAUDF.0).round();
    let e = q*TAUDF.1;
    let y = fma(q, TAUDF.0, -x);
    let z = HPIDF.quick_add_to_f32(-Df32(y,e).abs());
    sinf_poly(z)
}

#[inline(always)]
pub fn cbrt(x: f32) -> f32 {
    let s = f32::from_bits(x.to_bits() / 3 + 0x2a509a07u32);
    let s2 = s * s;
    fma(
        fma(0.6 * s, s2, 0.3 * x),
        fma(s2, -s2, x * s) / fma(fma(s, s2, 1.6 * x), s * s2, x * x * 0.1),
        s,
    )
}

#[inline(always)]
pub fn cbrt_accurate(x: f32) -> f32 {
    let s = f32::from_bits(0x2a4ddef1u32.wrapping_add(x.to_bits() / 3));
    let r = f32::from_bits(0x68ff2381u32.wrapping_sub((x.to_bits() / 3) << 1));
    let s = fma(s * s, s * -r, fma(r, x, s));
    let s = fma(s * s, s * -r, fma(r, x, s));
    let s2 = Df32::from_mul(s, s);
    let s32x = {
        let b = fma(s2.0, s * 2., x);
        let p = x - b;
        let e = fma(s2.0, s * 2., p) - (p + b - x);
        let lo = fma(s2.1, s * 2., e);
        Df32(b, lo)
    };
    let s2xps4 = {
        let s40 = s2.0 * s2.0;
        let e = fma(s2.0, s2.0, -s40);
        let s41 = fma(s2.0 * 2., s2.1, fma(s2.1, s2.1, e));
        let p = s * 2. * x;
        let e = fma(s * 2., x, -p);
        let s = p + s40;
        Df32(s, s40 - (s - p) + e + s41)
    };
    return s2xps4.div_to_f32(s32x);
}

// higher throughput cbrt experiment, 7 ulp average error
fn cbrt_throughput(x: f32) -> f32 {
    let r = f32::from_bits(0xd461ff81u32.wrapping_sub(x.to_bits() / 3));
    let r = fma(r * r, (r * r) * x, r * f32::from_bits(0x3fb6e3d7));
    let r = fma(r * r, (r * r) * x, r * f32::from_bits(0x3fe09c2a));
    r * r * x
}

fn cbrt_blazing(x: f32) -> f32 {
    f32::from_bits(0x2a4ddef1u32.wrapping_add(x.to_bits() / 3))
}
fn sqrt_blazing(x: f32) -> f32 {
    f32::from_bits(0x1FBD22DFu32.wrapping_add(x.to_bits() >> 1))
}
fn rcp_blazing(x: f32) -> f32 {
    f32::from_bits(0x7EEF370Bu32.wrapping_sub(x.to_bits()))
}
fn exp2_blazing(x: f32) -> f32 {
    f32::from_bits((x + 383.).to_bits() << 8 & 0x7FFFFFFF)
}
fn log2_blazing(x: f32) -> f32 {
    f32::from_bits((x).to_bits() >> 8 | 0x43800000) - 383.
}
fn rsqrt_blazing(x: f32) -> f32 {
    f32::from_bits(0x5F33E79Fu32.wrapping_sub(x.to_bits() >> 1))
}

// 50 average ulp error 32 cycle latency 5.5 cycle rthroughput
fn cbrt_fast(x: f32) -> f32 {
    let s = f32::from_bits(0x2a4ddef1u32.wrapping_add(x.to_bits() / 3));
    let r = f32::from_bits(0x68ff2381u32.wrapping_sub((x.to_bits() / 3) << 1));
    let s = fma(s * s, s * -r, fma(r, x, s));
    fma(s * s, s * -r, fma(r, x, s))
}

pub fn cbrt_constant(x: f32, c: &[u32]) -> f32 {
    let s = f32::from_bits(c[0].wrapping_add(x.to_bits() / 3));
    let r = f32::from_bits(c[1].wrapping_sub((x.to_bits() / 3) << 1));
    return fma(s * s, s * -r, fma(r, x, s));

    let s2 = Df32::from_mul(s, s);
    let s3 = s2 * s;
    let s3x = s3 * x;
    let s6 = s3 * s3;
    let x2 = Df32::from_mul(x, x);
    // h3
    fma(
        s,
        (
            Df32(fma(-2., s6.0, s3x.0) + x2.0, -2. * s6.1 + s3x.1 + x2.1)
            //(Df32(s6.0*-2.,s6.1*-2.).quick_add_df(s3x) + x2) * 3.
        )
        .div_to_f32(
            Df32(
                fma(s3x.0, 16. / 3., fma(s6.0, 10. / 3., x2.0 * (1. / 3.))),
                fma(s3x.1, 16. / 3., fma(s6.1, 10. / 3., x2.1 * (1. / 3.))),
            ), //Df32(s3x.0*16.,s3x.1*16.) + s6 * 10. + x2
        ),
        s,
    )

    //h4
    //s - ((s3 - x)*(Df32(s3x.0*16.,s3x.1*16.).quick_add_df(10.*s6) + x2)).div_to_f32(
    //    (15.*s6 + 15.*x2 + 51.*s3x)*s2
    //)

    /*
    let s = f32::from_bits((x.to_bits() / 3).wrapping_add(
        c0
        //0x2a53c472
    )
    );
    let rs = f32::from_bits(
        c1
        //0x543846f5
         .wrapping_sub( x.to_bits()/ 3));
    let s = fma(fma(s*s,-s,x),rs*rs,s);

    let s32 = Df32::from_mul(s,s) * (s*2.);
    s * 2f32 + ((s32*-1.5)*s).div_to_f32(s32.quick_add(x))*/
    //s
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;
    use rand::RngExt;

    fn run_descent(
        f: impl Fn(f32, &[u32]) -> f32,
        reference: impl Fn(f32) -> f32,
        initial_consts: &[u32],
    ) {
        let mut consts: Vec<u32> = initial_consts.to_vec();
        let iters = 100;

        let mut best_err: u64 = 0;
        let mut steps: u64 = 0;
        for _ in 1..iters {
            let x: f32 = rand::rng().random::<f32>().abs();
            let result = f(x, &consts);
            best_err += reference(x).to_bits().abs_diff(result.to_bits()) as u64;
        }
        let mut tryagain = false;
        println!(
            "optimizing! starting error: {}",
            best_err as f64 / iters as f64
        );
        let mut deltas: Vec<u32> = consts
            .iter()
            .map(|&_| {
                let n: u32 = rand::rng().random();
                f32::from_bits(n) as i32 as u32
            })
            .collect();
        loop {
            let mut new_err: u64 = 0;
            let mut old_err: u64 = 0;
            if !tryagain {
                deltas = consts
                    .iter()
                    .map(|&_| {
                        let n: u32 = rand::rng().random();
                        return f32::from_bits(n) as i32 as u32;
                        if steps < 10 {
                            n
                        } else {
                            if steps & 1 == 0 {
                                n
                            } else {
                                f32::from_bits(n) as i32 as u32
                            }
                        }
                    })
                    .collect();
            }
            steps += 1;
            let new_consts: Vec<u32> = consts
                .iter()
                .zip(deltas.iter())
                .map(|(&c, &d)| c.wrapping_add(d))
                .collect();
            if new_consts.iter().zip(consts.iter()).all(|(&x, &y)| x == y) {
                continue;
            }
            for _ in 1..iters {
                let x: f32 = rand::rng().random::<f32>().abs();
                let ref_val = reference(x);
                let old_result = f(x, &consts);
                let new_result = f(x, &new_consts);
                new_err += ref_val.to_bits().abs_diff(new_result.to_bits()) as u64;
                old_err += ref_val.to_bits().abs_diff(old_result.to_bits()) as u64;
            }
            if new_err <= old_err {
                consts = new_consts.clone();
                if new_err < best_err {
                    let mut new_err_nomul: u64 = 0;
                    // no mul check
                    for _ in 1..iters {
                        let x: f32 = rand::rng().random::<f32>().abs();
                        let new_result = f(x, &new_consts);
                        new_err_nomul +=
                            reference(x).to_bits().abs_diff(new_result.to_bits()) as u64;
                    }
                    if new_err_nomul < best_err {
                        best_err = new_err_nomul;
                        let const_strs: Vec<String> =
                            consts.iter().map(|c| format!("0x{:x}u32", c)).collect();
                        println!(
                            "new best consts {} with error {} nomul: {} step:{}",
                            const_strs.join(", "),
                            new_err as f64 / iters as f64,
                            new_err_nomul as f64 / iters as f64,
                            steps
                        );
                    }
                    tryagain = true;
                } else {
                    tryagain = false;
                }
            } else {
                tryagain = false;
            }
            if new_err == best_err {
                best_err = new_err;
                consts = new_consts;
            }
        }
    }
/*
    #[test]
    fn descent2() {
        run_descent(
            |x, consts| cbrt_constant(x, consts),
            |x| (x as f64).cbrt() as f32,
            &[0x2a4dc461u32, 0x6900943cu32],
        );
    }*/

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
            ulp_error(1..1000, 1.0, cbrt, |x| x.cbrt())
        );
        println!(
            "std   cbrt error: {}",
            ulp_error(1..1000, 1.0, |x| x.cbrt(), |x| x.cbrt())
        );
    }
    #[test]
    fn cbrt_accurate_precision() {
        println!(
            "jodie cbrt accurate error: {}",
            ulp_error(1..1000, 1.0, cbrt_accurate, |x| x.cbrt())
        );
        println!(
            "std   cbrt error: {}",
            ulp_error(1..1000, 1.0, |x| x.cbrt(), |x| x.cbrt())
        );
    }
    #[test]
    fn exp2_precision() {
        println!(
            "jodie exp2 error: {}",
            ulp_error(-100..100, 0.1, exp2, |x| x.exp2())
        );
        println!(
            "std   exp2 error: {}",
            ulp_error(-100..100, 0.1, |x| x.exp2(), |x| x.exp2())
        );
    }
    #[test]
    fn log2_precision() {
        println!(
            "jodie log2 error: {}",
            ulp_error(2..1000, 1.0, log_2, |x| x.log2())
        );
        println!(
            "std   log2 error: {}",
            ulp_error(2..1000, 1.0, |x| x.log2(), |x| x.log2())
        );
    }
    #[test]
    fn sin_precision() {
        println!(
            "jodie sin error: {}",
            ulp_error(-100..100, 0.1, sin, |x| x.sin())
        );
        println!(
            "std   sin error: {}",
            ulp_error(-100..100, 0.1, |x| x.sin(), |x| x.sin())
        );
    }
    #[test]
    fn cos_precision() {
        println!(
            "jodie cos error: {}",
            ulp_error(-100..100, 0.1, cos, |x| x.cos())
        );
        println!(
            "std   cos error: {}",
            ulp_error(-100..100, 0.1, |x| x.cos(), |x| x.cos())
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
        let y_min = all_y.iter().cloned().filter(|y| y.is_finite()).fold(f32::INFINITY, f32::min);
        let y_max = all_y.iter().cloned().filter(|y| y.is_finite()).fold(f32::NEG_INFINITY, f32::max);
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
        plot_approx("cbrt_approx.png", 1., 128., cbrt_blazing, |x| x.cbrt());
    }
    #[test]
    fn sqrt_approx_plot() {
        plot_approx("sqrt_approx.png", 1., 128., sqrt_blazing, |x| x.sqrt());
    }
    #[test]
    fn rcp_approx_plot() {
        plot_approx("rcp_approx.png", 1., 10., rcp_blazing, |x| 1.0 / x);
    }
    #[test]
    fn exp2_approx_plot() {
        plot_approx("exp2_approx.png", 0., 10., exp2_blazing, |x| x.exp2());
    }
    #[test]
    fn log2_approx_plot() {
        plot_approx("log2_approx.png", 1., 128., log2_blazing, |x| x.log2());
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
        plot_approx("rsqrt_approx.png", 1., 128., rsqrt_blazing, |x| {
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
}
