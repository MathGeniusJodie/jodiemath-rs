// SCRATCH (w2): exhaustive ulp probe of dawson against an f64 reference
// independent of the harness's Simpson quadrature. Delete before committing.
//
// |x| <= 7: D(x) = exp(-x^2) * sum_{n>=0} x^(2n+1)/(n!(2n+1)).  The sum is
// all-positive, so there is no cancellation anywhere -- the exp(-x^2) that
// cancels the sum's growth is applied once, at the end.  Validated against
// scipy.special.dawsn to <= 6.2e-16 relative at x = 1e-5..10.
// |x| > 7: the double-factorial asymptotic series, exact to f64 there
// (0.0 relative against scipy at x = 7, 8, 10, 20).
use jodiemath_rs::*;

fn dawson_ref(x: f64) -> f64 {
    let ax = x.abs();
    let mag = if ax <= 7.0 {
        let mut s = 0.0f64;
        let mut t = ax;
        for n in 0..200 {
            let d = t / (2 * n + 1) as f64;
            s += d;
            if d < s * 1e-19 {
                break;
            }
            t = t * ax * ax / (n + 1) as f64;
        }
        (-ax * ax).exp() * s
    } else {
        let v = 1.0 / (2.0 * ax * ax);
        let mut acc = 1.0f64;
        let mut term = 1.0f64;
        for k in 1..=20 {
            term *= (2 * k - 1) as f64 * v;
            acc += term;
        }
        acc / (2.0 * ax)
    };
    if x < 0.0 { -mag } else { mag }
}

const PC: [f64; 7] = [
    1.0, -0.085751414, 0.037434783, -0.0004054072, 0.00019858626, 5.6392253e-8, 2.8313497e-8,
];
const QC: [f64; 7] = [
    1.0, 0.5809171, 0.15803601, 0.026255792, 0.002856104, 0.00021274923, 5.002549e-6,
];

/// the shipped central rational, evaluated in f64 (Horner, no rounding of
/// consequence) -- isolates the f32 evaluation chain from the fit itself
fn central_f64(u: f64) -> f64 {
    let mut n = PC[6];
    let mut d = QC[6];
    for k in (0..6).rev() {
        n = n * u + PC[k];
        d = d * u + QC[k];
    }
    n / d
}

/// mode 0 shipped; 1 = f64 rational on fl(x*x); 2 = f64 rational on exact
/// x^2; 3 = oracle x * fl(dawsn(x)/x)
fn dawson_mode(x: f32, mode: u32) -> f32 {
    if mode == 0 || x.abs() > 4.0 {
        return dawson(x);
    }
    let xd = x as f64;
    let r = match mode {
        1 => central_f64((x * x) as f64),
        2 => central_f64(xd * xd),
        _ => {
            let t = dawson_ref(xd) / xd;
            (t as f32) as f64
        }
    };
    x * (r as f32)
}

fn ulp_diff(got: f32, want_f64: f64) -> f64 {
    let want = want_f64 as f32;
    if got.is_nan() && want.is_nan() {
        return 0.0;
    }
    if got == want {
        return 0.0;
    }
    if !got.is_finite() || !want.is_finite() {
        return f64::INFINITY;
    }
    let b = want.abs().to_bits().max(1);
    let u = (f32::from_bits(b + 1) - f32::from_bits(b)) as f64;
    ((got as f64) - want_f64).abs() / u
}

type Slot = (f64, f64, u64, f32);

fn sweep(lo_b: u32, hi_b: u32, stride: u32, mode: u32) -> (f64, f32, f64, u64, Vec<Slot>) {
    let mut worst = 0.0f64;
    let mut worst_x = 0.0f32;
    let mut sum = 0.0f64;
    let mut n = 0u64;
    let mut oct: Vec<Slot> = vec![(0.0, 0.0, 0, 0.0); 256];
    let mut b = lo_b;
    while b < hi_b {
        let x = f32::from_bits(b);
        let d = ulp_diff(dawson_mode(x, mode), dawson_ref(x as f64));
        sum += d;
        n += 1;
        if d > worst {
            worst = d;
            worst_x = x;
        }
        let slot = &mut oct[(b >> 23) as usize];
        slot.0 += d;
        slot.2 += 1;
        if d > slot.1 {
            slot.1 = d;
            slot.3 = x;
        }
        b += stride;
    }
    (worst, worst_x, sum, n, oct)
}

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let lo: f32 = a.first().map(|s| s.parse().unwrap()).unwrap_or(1e-3);
    let hi: f32 = a.get(1).map(|s| s.parse().unwrap()).unwrap_or(32.0);
    let stride: u32 = a.get(2).map(|s| s.parse().unwrap()).unwrap_or(1);
    let mode: u32 = a.get(3).map(|s| s.parse().unwrap()).unwrap_or(0);
    let lo_b = lo.to_bits();
    let hi_b = hi.to_bits();
    println!("dawson over [{lo:e}, {hi:e}) : {} patterns, stride {stride}", hi_b - lo_b);

    const NT: u32 = 4;
    let chunk = ((hi_b - lo_b) as u64).div_ceil(NT as u64) as u32;
    let mut hs = Vec::new();
    for i in 0..NT {
        let s0 = lo_b + ((i * chunk) / stride) * stride;
        let s1 = ((s0 as u64).max((lo_b + (i + 1) * chunk) as u64).min(hi_b as u64)) as u32;
        hs.push(std::thread::spawn(move || sweep(s0, s1, stride, mode)));
    }
    let mut worst = 0.0f64;
    let mut worst_x = 0.0f32;
    let mut sum = 0.0f64;
    let mut n = 0u64;
    let mut oct: Vec<Slot> = vec![(0.0, 0.0, 0, 0.0); 256];
    for h in hs {
        let (w, wx, s, c, o) = h.join().unwrap();
        if w > worst {
            worst = w;
            worst_x = wx;
        }
        sum += s;
        n += c;
        for (i, e) in o.iter().enumerate() {
            oct[i].0 += e.0;
            oct[i].2 += e.2;
            if e.1 > oct[i].1 {
                oct[i].1 = e.1;
                oct[i].3 = e.3;
            }
        }
    }
    for (i, (s, m, c, wx)) in oct.iter().enumerate() {
        if *c == 0 {
            continue;
        }
        println!(
            "  2^{:<5} [{:e},{:e})  avg {:8.4}  max {:9.3} at {:e}  ({} pts)",
            i as i32 - 127,
            f32::from_bits((i as u32) << 23),
            f32::from_bits(((i + 1) as u32) << 23),
            s / *c as f64,
            m,
            wx,
            c
        );
    }
    println!("mode {mode}  dawson avg {:.4}  max {:.3} at x = {:e} ({:#010x})", sum / n as f64, worst, worst_x, worst_x.to_bits());
}
