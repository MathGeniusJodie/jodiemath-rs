// Idea #169: per-function ULP error *structure*, not just avg/max. A single
// max tells you nothing about whether the error is one isolated pocket, a
// broad plateau, or a step at a branch seam -- and that distinction is what
// decides whether a refit, a seam move, or a new sub-branch is the right
// lever. Locating an error concentration is a proven approach in this crate:
// it is what rescued atanh (idea #69), where a rejected "real accuracy loss"
// turned out to be confined to a narrow region around x~0.111 and a
// dedicated small-x branch there converted it into a win on every axis.
//
// Reports, per function: a histogram over ulp buckets, and a per-region
// breakdown so a concentration is visible as a region whose max dominates
// the global one. Regions are the crate's own branch structure where it is
// known (the |x| < seam split), plus a log-magnitude sweep either side.
use jodiemath_rs::*;

fn ulp_err(got: f32, want64: f64) -> f64 {
    if got.is_nan() && want64.is_nan() {
        return 0.0;
    }
    let want = want64 as f32;
    if got == want {
        return 0.0;
    }
    if !got.is_finite() || !want.is_finite() {
        return f64::INFINITY;
    }
    let a = want.abs();
    let ulp = if a == 0.0 {
        f32::from_bits(1) as f64
    } else {
        (f32::from_bits(a.to_bits() + 1) - a) as f64
    };
    ((got as f64) - (want as f64)).abs() / ulp
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
}

fn profile(
    name: &str,
    f: impl Fn(f32) -> f32,
    r: impl Fn(f64) -> f64,
    lo: f32,
    hi: f32,
    seam: f32,
) {
    println!("\n=== {name}  over |x| in [{lo:e}, {hi:e}], known seam at {seam} ===");
    let mut rng = Rng(0x51ed_2701_a3bf_9d14);

    // ulp histogram
    let mut hist = [0u64; 12];
    let mut total = 0u64;
    // per-decade breakdown
    const NB: usize = 24;
    let mut b_max = [0.0f64; NB];
    let mut b_sum = [0.0f64; NB];
    let mut b_n = [0u64; NB];
    let mut b_worst = [0.0f32; NB];
    let l_lo = (lo as f64).log10();
    let l_hi = (hi as f64).log10();

    const N: u64 = 24_000_000;
    for _ in 0..N {
        // log-uniform magnitude, random sign
        let u = rng.next_u32() as f64 / u32::MAX as f64;
        let mag = 10f64.powf(l_lo + (l_hi - l_lo) * u) as f32;
        let x = if rng.next_u32() & 1 == 0 { mag } else { -mag };
        if !x.is_finite() || x == 0.0 {
            continue;
        }
        let e = ulp_err(f(x), r(x as f64));
        if !e.is_finite() {
            continue;
        }
        total += 1;
        let hb = (e.ceil() as usize).min(11);
        hist[hb] += 1;

        let li = (((mag as f64).log10() - l_lo) / (l_hi - l_lo) * NB as f64) as usize;
        let li = li.min(NB - 1);
        b_n[li] += 1;
        b_sum[li] += e;
        if e > b_max[li] {
            b_max[li] = e;
            b_worst[li] = x;
        }
    }

    print!("  ulp histogram: ");
    for (i, c) in hist.iter().enumerate() {
        if *c > 0 {
            let pct = 100.0 * *c as f64 / total as f64;
            let lbl = if i == 11 {
                ">10".to_string()
            } else {
                format!("<={i}")
            };
            print!("{lbl}:{pct:.2}%  ");
        }
    }
    println!();

    let gmax = b_max.iter().cloned().fold(0.0f64, f64::max);
    println!("  per-magnitude-band (global max {gmax:.2}):");
    for i in 0..NB {
        if b_n[i] == 0 {
            continue;
        }
        let bl = 10f64.powf(l_lo + (l_hi - l_lo) * i as f64 / NB as f64);
        let bh = 10f64.powf(l_lo + (l_hi - l_lo) * (i + 1) as f64 / NB as f64);
        let mark = if b_max[i] >= gmax - 0.001 {
            "  <-- carries the max"
        } else {
            ""
        };
        let side = if bh <= seam as f64 {
            "below seam"
        } else if bl >= seam as f64 {
            "above seam"
        } else {
            "STRADDLES seam"
        };
        println!(
            "    [{:>9.2e},{:>9.2e}) {:>14}  avg {:.4}  max {:>6.2}  worst x {:>13e}{}",
            bl,
            bh,
            side,
            b_sum[i] / b_n[i] as f64,
            b_max[i],
            b_worst[i],
            mark
        );
    }
}

fn main() {
    // The three joint-worst budgeted functions (max ulp 6 each).
    profile("expm1", expm1, f64::exp_m1, 1e-3, 80.0, 0.5);
    profile(
        "exp_m1_over_x",
        exp_m1_over_x,
        |v| if v == 0.0 { 1.0 } else { v.exp_m1() / v },
        1e-3,
        80.0,
        0.5,
    );
    profile("tanh", tanh, f64::tanh, 1e-3, 20.0, 0.25);
}
