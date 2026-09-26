// Scores sin(x) ~= x + x^3*(c0 + c1 y + c2 y^2 + c3 y^3), y = x^2, over every
// f32 in [0, 1.5707965] (the reach of every caller's reduced argument): max
// and mean ulp against f64 sin rounded to f32, and how many results exceed 1.
//   cargo run --release --example sinpoly_tune [c0 c1 c2 c3]   (score)
//   cargo run --release --example sinpoly_tune search c0 c1 c2 c3
use std::sync::atomic::{AtomicU64, Ordering};

const BLOCK: usize = 4096;
const TOP: u32 = 0x3fc9_0fdc; // 1.5707965

#[inline(always)]
fn poly(x: f32, c: [f32; 4]) -> f32 {
    let y = x * x;
    let y2 = y * y;
    let x3 = y * x;
    let a = c[1].mul_add(y, c[0]);
    let b = c[3].mul_add(y, c[2]);
    let p = b.mul_add(y2, a);
    p.mul_add(x3, x)
}

fn score(c: [f32; 4]) -> (u64, f64, u64) {
    score_from(c, 0)
}

/// Score over [from, 1.5707965] only.
fn score_from(c: [f32; 4], from: u32) -> (u64, f64, u64) {
    let threads = std::thread::available_parallelism().map_or(8, |n| n.get());
    let next = AtomicU64::new(from as u64);
    let parts: Vec<(u64, u64, u64, u64)> = std::thread::scope(|s| {
        (0..threads)
            .map(|_| {
                s.spawn(|| {
                    let (mut max, mut sum, mut n, mut over) = (0u64, 0u64, 0u64, 0u64);
                    let mut xs = [0f32; BLOCK];
                    let mut ys = [0f32; BLOCK];
                    loop {
                        let base = next.fetch_add(BLOCK as u64, Ordering::Relaxed);
                        if base > TOP as u64 {
                            break;
                        }
                        for (k, v) in xs.iter_mut().enumerate() {
                            *v = f32::from_bits((base as u32 + k as u32).min(TOP));
                        }
                        for (y, &x) in ys.iter_mut().zip(xs.iter()) {
                            *y = poly(x, c);
                        }
                        for k in 0..BLOCK {
                            if base + k as u64 > TOP as u64 {
                                break;
                            }
                            let r = (xs[k] as f64).sin() as f32;
                            let d = (ys[k].to_bits() as i64 - r.to_bits() as i64).unsigned_abs();
                            max = max.max(d);
                            sum += d;
                            n += 1;
                            over += (ys[k] > 1.0) as u64;
                        }
                    }
                    (max, sum, n, over)
                })
            })
            .collect::<Vec<_>>()
            .into_iter()
            .map(|h| h.join().unwrap())
            .collect()
    });
    let max = parts.iter().map(|p| p.0).max().unwrap();
    let sum: u64 = parts.iter().map(|p| p.1).sum();
    let n: u64 = parts.iter().map(|p| p.2).sum();
    let over: u64 = parts.iter().map(|p| p.3).sum();
    (max, sum as f64 / n as f64, over)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let default = [-0.166_666_6_f32, 8.3330662e-3, -1.9809603e-4, 2.6057806e-6];
    let search = args.get(1).is_some_and(|a| a == "search");
    let off = if search { 2 } else { 1 };
    let mut c = default;
    if args.len() >= off + 4 {
        for i in 0..4 {
            c[i] = args[off + i].parse().unwrap();
        }
    }
    let (m, a, o) = score(c);
    let top = score_from(c, 0x3f80_0000);
    println!(
        "{:?}: max {m} avg {a:.5} over1 {o}  [1, pi/2]: max {} avg {:.5}",
        c, top.0, top.1
    );
    if args.get(1).is_some_and(|a| a == "scan") {
        // Grid over c2/c3 ulp offsets: print every point with no result above 1.
        for d3 in (-6000i32..=0).step_by(250) {
            for d2 in (-400i32..=400).step_by(50) {
                let mut t = default;
                t[2] = f32::from_bits((t[2].to_bits() as i32 + d2) as u32);
                t[3] = f32::from_bits((t[3].to_bits() as i32 + d3) as u32);
                let s = score(t);
                if s.2 == 0 {
                    println!("d2 {d2} d3 {d3} {:?}: max {} avg {:.5}", t, s.0, s.1);
                }
            }
        }
        return;
    }
    if !search {
        return;
    }
    // Coordinate descent in ulp steps from a feasible start: minimise avg
    // subject to no result above 1 and max <= 2.
    // Weight the top of the range: cos of small arguments all lands near pi/2.
    let full = |t: [f32; 4]| {
        let a = score(t);
        let b = score_from(t, 0x3f80_0000);
        (a.0, a.1 + 20.0 * b.1, a.2)
    };
    let key = |s: (u64, f64, u64)| (s.2 > 0 || s.0 > 2, (s.1 * 1e7) as u64);
    let mut best = (c, full(c));
    for step in [64i32, 16, 4, 1] {
        let mut improved = true;
        while improved {
            improved = false;
            for i in 0..4 {
                for d in [-step, step] {
                    let mut t = best.0;
                    t[i] = f32::from_bits((t[i].to_bits() as i32 + d) as u32);
                    let s = full(t);
                    if key(s) < key(best.1) {
                        println!("  {:?}: max {} avg {:.5} over1 {}", t, s.0, s.1, s.2);
                        best = (t, s);
                        improved = true;
                    }
                }
            }
        }
    }
    println!(
        "best {:?}: max {} avg {:.5} over1 {}",
        best.0, best.1 .0, best.1 .1, best.1 .2
    );
}
