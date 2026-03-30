// godbolt flags -C opt-level=3 -C target_feature=+fma

mod doublefloat;
use doublefloat::Df32;

const SIGN_MASK: u32 = 0x80000000;
const EXPONENT_MASK: u32 = 0x7f800000;
const MANTISSA_MASK: u32 = 0x007fffff;

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
fn exp2_fract(x: f32) -> f32 {
    let a = f32::from_bits(0x3af71c15);
    let b = f32::from_bits(0x3c130514);
    let c = f32::from_bits(0x3d64b437);
    let d = f32::from_bits(0x3e75ea9e);
    let e = f32::from_bits(0x3f317271);
    fma(fma(fma(fma(fma(a, x, b), x, c), x, d), x, e), x, 1_f32)
    //alternative for less latency
    //fma(x*x,fma(fma(a,x,b),x*x,fma(c,x,d)),fma(e,x,1.))
}
#[inline(always)]
pub fn exp2(x: f32) -> f32 {
    // exp2(floor(x))*exp2(fract(x)) == exp2(x)
    let exp2int = f32::from_bits(((x + 383_f32).to_bits() << 8) & EXPONENT_MASK);
    let fract = x - x.floor();
    exp2int * exp2_fract(fract)
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
    let tau: f64 = std::f64::consts::TAU;
    let tauh: f32 = tau as f32;
    let taul: f32 = (tau - (tauh as f64)) as f32;
    let rtau: f64 = 0.159_154_943_091_895_35;
    let rtauh: f32 = rtau as f32;
    //let rtaul: f32 = (rtau - (rtauh as f64)) as f32;
    let pi = std::f32::consts::PI;
    let pilo = f32::from_bits(0x33bbbd2e);
    let rpi = 0.318_309_87;
    let z = x - (x * rtauh).round() * tauh + (x * rtauh).round() * taul;
    let y = (x - (x * rpi).round() * pi + (x * rpi).round() * pilo).abs();
    sinf_poly(mulsign(y, z))
}
#[inline(always)]
// todo: make less convoluted
pub fn cos(x: f32) -> f32 {
    let tau: f64 = std::f64::consts::TAU;
    let tauh: f32 = tau as f32;
    let taul: f32 = (tau - (tauh as f64)) as f32;
    let rtau: f64 = 0.159_154_943_091_895_35;
    let rtauh: f32 = rtau as f32;
    let rtaul: f32 = (rtau - (rtauh as f64)) as f32;

    let m = (fma(x, rtauh * 2f32, x * (rtaul * 2f32))).floor();
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

#[inline(always)]
pub fn cbrt(x: f32) -> f32 {
    let s = f32::from_bits(x.to_bits() / 3 + 709982100);
    let s2 = s * s;
    fma(
        0.3 * fma(2.*s2, s2, s*x),
        fma(s2, -s, x) / fma(fma(s, s2, 1.6*x), s*s2, x*x*0.1),
        s,
    )
}

#[inline(always)]
pub fn cbrt_accurate(x: f32) -> f32 {
    let s = f32::from_bits(
        0x2a4f536eu32.wrapping_add(x.to_bits() / 3)
    );
    let r = f32::from_bits(
        0x69043c30u32.wrapping_sub((x.to_bits()/3)<<1)
    );
    let s = fma(s*s,s*-r,fma(r,x,s));
    let s = fma(s*s,s*-r,fma(r,x,s));
    let s3 = Df32::from_mul(s,s) * s;
    return s * 2f32 - ((s3*1.5)*s).div_to_f32(s3.quick_add(x*0.5));
}
/*
-1/3*((s^3 - x)*(10*s^6 + 16*s^3*x + x^2))/
(s^2*(5*s^6 + 17*s^3*x + 5*x^2))
*/

pub fn cbrt_constant(x: f32, c:u32) -> f32 {
    let s = f32::from_bits(x.to_bits() / 3 + 
    c
    //710137493
    );
    let rs = f32::from_bits(
        //c
        1412840964
         - x.to_bits() / 3);
    let s = fma(fma(s*s,-s,x),rs*rs,s);

    let s32 = Df32::from_mul(s,s) * (s*2.);
    s * 2f32 + ((s32*-1.5)*s).div_to_f32(s32.quick_add(x))
}


pub fn cbrt_constant_2(x: f32, c0:u32, c1:u32) -> f32 {
    let s = f32::from_bits(
        c0.wrapping_add(x.to_bits() / 3)
    );
    let r = f32::from_bits(
        c1.wrapping_sub((x.to_bits()/3)<<1)
    );
    let s = fma(s*s,s*-r,fma(r,x,s));
    let s = fma(s*s,s*-r,fma(r,x,s));
    //let s2 = s * s;
    //let s = fma(s2*s2, -1.5 / fma(s, s2, x*0.5), s * 2.);
    //return s;
    let s3 = Df32::from_mul(s,s) * s;
    return s * 2f32 - ((s3*1.5)*s).div_to_f32(s3.quick_add(x*0.5));

    let s2 = Df32::from_mul(s,s);
    let s3 = s2 * s;
    let s3x = s3 * x;
    let s6 = s3*s3;
    let x2 = Df32::from_mul(x,x);
    // h3
    fma(s,(
            Df32(fma(-2.,s6.0,s3x.0)+x2.0,-2.*s6.1+s3x.1+x2.1)
            //(Df32(s6.0*-2.,s6.1*-2.).quick_add_df(s3x) + x2) * 3.
        ).div_to_f32(
            Df32(
                fma(s3x.0,16./3.,fma(s6.0,10./3.,x2.0*(1./3.))),
                fma(s3x.1,16./3.,fma(s6.1,10./3.,x2.1*(1./3.))))
            //Df32(s3x.0*16.,s3x.1*16.) + s6 * 10. + x2
        )
    ,s)

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
    
    #[test]
    fn descent2() {
        let mut c0:u32 = 0x2a4f536e;
        let mut c1:u32 = 0x69043c30;

        let iters = 10_000_000;

        let mut best_err: u64 = 0;
        for x in 1..iters {
            let x:f32 = rand::rng().random::<f32>().abs();
            let reference = (x as f64).cbrt() as f32;
            let result = cbrt_constant_2(x as f32,c0,c1);
            best_err += reference.to_bits().abs_diff(result.to_bits()) as u64;
        }
        println!("optimizing! starting error: {}",best_err as f64 / iters as f64);
        loop{
            let mut new_err: u64 = 0;
            let mut old_err: u64 = 0;
            let n0: u32 =  rand::rng().random();
            let n1: u32 =  rand::rng().random();
            let new_const0 = c0.wrapping_add(f32::from_bits(n0) as i32 as u32);
            let new_const1 = c1.wrapping_add(f32::from_bits(n1) as i32 as u32);
            for x in 1..iters {
                let x:f32 = rand::rng().random::<f32>().abs();
                let reference = (((x as f32)) as f64).cbrt() as f32;
                let old_result = cbrt_constant_2((x as f32), c0,c1);
                let new_result = cbrt_constant_2((x as f32), new_const0, new_const1);
                new_err += reference.to_bits().abs_diff(new_result.to_bits()) as u64;
                old_err += reference.to_bits().abs_diff(old_result.to_bits()) as u64;
            }
            if new_err <= old_err{
                c0 = new_const0;
                c1 = new_const1;
                if new_err < best_err{

                    let mut new_err_nomul: u64 = 0;
                    // no mul check
                    for x in 1..iters {
                        let x:f32 = rand::rng().random::<f32>().abs();
                        let reference = (((x as f32)) as f64).cbrt() as f32;
                        let new_result = cbrt_constant_2((x as f32), new_const0, new_const1);
                        new_err_nomul += reference.to_bits().abs_diff(new_result.to_bits()) as u64;
                    }
                    if new_err_nomul < best_err{
                        best_err = new_err_nomul;
                        println!("new best consts {:x} {:x} with error {} nomul: {}",c0,c1,new_err as f64/iters as f64, new_err_nomul as f64/iters as f64);
                    }
                }
            }
            if new_err == best_err{
                best_err = new_err;
                c0 = new_const0;
                c1 = new_const1;
                //print!(".");
            }
        }
    }

    /*
    #[test]
    fn descent() {
        let mut c = 709982101;
        let mut step = 1;

        let mut best_err: u64 = 0;
        for x in 1..1000 {
            let reference = (x as f64).cbrt() as f32;
            let result = cbrt_constant(x as f32,c);
            best_err += reference.to_bits().abs_diff(result.to_bits()) as u64;
        }
        println!("optimizing! starting error: {}",best_err);
        let mut u = c;
        let mut d = c; 
        loop{
            u += step;
            d -= step;

            let mut u_err: u64 = 0;
            let mut d_err: u64 = 0;
            for x in 1..1000 {
                let reference = (x as f64).cbrt() as f32;
                let u_result = cbrt_constant(x as f32, u);
                let d_result = cbrt_constant(x as f32, d);
                u_err += reference.to_bits().abs_diff(u_result.to_bits()) as u64;
                d_err += reference.to_bits().abs_diff(d_result.to_bits()) as u64;
            }
            if u_err > best_err && d_err > best_err {
                println!("new best const {} with error {}",c,best_err);
                break;
            }
            if d_err < best_err{
                best_err = d_err;
                c = d;
                u = d;
                //println!("new best const {} with error {}",c,best_err);
            }
            if u_err < best_err{
                best_err = u_err;
                c = u;
                d = u;
                //println!("new best const {} with error {}",c,best_err);
            }
            if u_err >= best_err || d_err >= best_err {
                //println!("step {} too small, increasing",step);
                //step+=1;
            }
        }
    }*/

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
        println!("jodie cbrt error: {}", (err as f32) / 1000.);
        let mut err: u64 = 0;
        for x in 1..1000 {
            let reference = (x as f64).cbrt() as f32;
            let result = (x as f32).cbrt();
            err += reference.to_bits().abs_diff(result.to_bits()) as u64;
        }
        println!("std cbrt error: {}", (err as f32) / 1000.);
    }
    #[test]
    fn cbrt_accurate_precision() {
        let mut err: u64 = 0;
        for x in 1..1000 {
            let reference = (x as f64).cbrt() as f32;
            let result = cbrt_accurate(x as f32);
            err += reference.to_bits().abs_diff(result.to_bits()) as u64;
        }
        println!("jodie accurate cbrt error: {}", (err as f32) / 1000.);
        let mut err: u64 = 0;
        for x in 1..1000 {
            let reference = (x as f64).cbrt() as f32;
            let result = (x as f32).cbrt();
            err += reference.to_bits().abs_diff(result.to_bits()) as u64;
        }
        println!("std cbrt error: {}", (err as f32) / 1000.);
    }

    #[test]
    fn exp2_precision() {
        let mut err = 0;
        for x in -100..100 {
            let reference = ((x as f64) * 0.1).exp2() as f32;
            let result = exp2((x as f32) * 0.1);
            err += reference.to_bits().abs_diff(result.to_bits());
        }
        println!("jodie exp2 error: {}", (err as f32) / 200.);
        let mut err = 0;
        for x in -100..100 {
            let reference = ((x as f64) * 0.1).exp2() as f32;
            let result = (x as f32 * 0.1).exp2();
            err += reference.to_bits().abs_diff(result.to_bits());
        }
        println!("std exp2 error: {}", (err as f32) / 200.);
    }
    #[test]
    fn log2_precision() {
        let mut err = 0;
        for x in 2..1000 {
            let reference = (x as f64).log2() as f32;
            let result = log_2(x as f32);
            err += reference.to_bits().abs_diff(result.to_bits());
        }
        println!("jodie log2 error: {}", (err as f32) / 1000.);
        let mut err = 0;
        for x in 2..1000 {
            let reference = (x as f64).log2() as f32;
            let result = (x as f32).log2();
            err += reference.to_bits().abs_diff(result.to_bits());
        }
        println!("std log2 error: {}", (err as f32) / 1000.);
    }
    #[test]
    fn sin_precision() {
        let mut err = 0;
        for x in -100..100 {
            let reference = ((x as f64) * 0.1).sin() as f32;
            let result = sin((x as f32) * 0.1);
            err += reference.to_bits().abs_diff(result.to_bits());
        }
        println!("jodie sin error: {}", (err as f32) / 200.);
        let mut err = 0;
        for x in -100..100 {
            let reference = ((x as f64) * 0.1).sin() as f32;
            let result = (x as f32 * 0.1).sin();
            err += reference.to_bits().abs_diff(result.to_bits());
        }
        println!("std sin error: {}", (err as f32) / 200.);
    }
    #[test]
    fn cos_precision() {
        let mut err = 0;
        for x in -100..100 {
            let reference = ((x as f64) * 0.1).cos() as f32;
            let result = cos((x as f32) * 0.1);
            err += reference.to_bits().abs_diff(result.to_bits()) as u64;
        }
        println!("jodie cos error: {}", (err as f32) / 200.);
        let mut err = 0;
        for x in -100..100 {
            let reference = ((x as f64) * 0.1).cos() as f32;
            let result = (x as f32 * 0.1).cos();
            err += reference.to_bits().abs_diff(result.to_bits()) as u64;
        }
        println!("std cos error: {}", (err as f32) / 200.);
    }
    #[test]
    fn cbrt_plot() {
        println!("{}", cbrt_accurate(9.0));
        println!("{}", cbrt_accurate(17.0));
        println!("{}", cbrt_accurate(33.0));
        println!("{}", cbrt_accurate(65.0));
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
                    .map(|x| (x, cbrt_accurate(x) / x.cbrt() - 1.0)),
                &RED,
            ))
            .unwrap()
            .label("std")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], RED));
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
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], RED));
        chart
            .draw_series(LineSeries::new(
                (-200..200).map(|x| x as f32 * 0.1).map(|x| (x, x.cos())),
                &BLUE,
            ))
            .unwrap()
            .label("std")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLUE));
        chart.configure_series_labels().draw().unwrap();
        chart
            .draw_series(LineSeries::new(
                (-200..200).map(|x| x as f32 * 0.1).map(|x| (x, x.cos())),
                &BLUE,
            ))
            .unwrap()
            .label("std")
            .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLUE));
        chart.configure_series_labels().draw().unwrap();

        root.present().expect("Unable to write result to file, please make sure 'plotters' crate is in your Cargo.toml");
    }
}
