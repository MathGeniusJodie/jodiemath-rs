//! Differential test + first-principles perf check for the gather-free
//! register-permute Payne-Hanek (`reduce_pi_wide_x8`) against the shipped
//! autovectorized scalar path (`sin_wide`/`cos_wide`).
//!
//! Correctness claim under test: the x8 path is BIT-IDENTICAL to the scalar
//! path on every input pattern, because its three 29-bit chunks are
//! bit-identical to the gathered table values (verified exhaustively over
//! exponents in Python against pitable.rs) and every downstream op is the
//! same operation on the same lane values.
//!
//! Run: cargo run --release --example wide_x8

use jodiemath_rs::*;
use std::hint::black_box;

fn fill(xs: &mut [f32], seed: u64) {
    let mut s = seed | 1;
    for x in xs.iter_mut() {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        *x = f32::from_bits((s & 0xffff_ffff) as u32);
    }
}

fn main() {
    println!("== differential sweep: x8 path vs scalar sin_wide/cos_wide ==");
    let n0 = 1 << 21;
    let mut inputs = vec![0f32; n0];
    fill(&mut inputs, 42);
    // pin coverage of every exponent incl. denormals and inf/nan edges
    for e in 0..=255u32 {
        inputs[e as usize] = f32::from_bits(e << 23);
        inputs[256 + e as usize] = f32::from_bits((e << 23) | 1);
        inputs[512 + e as usize] = f32::from_bits((e << 23) | 0x007f_ffff);
        // adversarial: mantissas around pi-related near-wrap points
        inputs[768 + e as usize] = f32::from_bits((e << 23) | (e * 131));
    }
    let mut out_sin = vec![0f32; inputs.len()];
    let mut out_cos = vec![0f32; inputs.len()];
    unsafe {
        sin_wide_x8_slice(&inputs, &mut out_sin);
        cos_wide_x8_slice(&inputs, &mut out_cos);
    }
    let mut fails = 0usize;
    for (i, (&x, &s)) in inputs.iter().zip(out_sin.iter()).enumerate() {
        let want = sin_wide(x);
        if s.to_bits() != want.to_bits() {
            if fails < 8 {
                println!(
                    "SIN MISMATCH idx={i} x={x:e} ({:#010x}): got {s:e} ({:#010x}) want {want:e} ({:#010x})",
                    x.to_bits(),
                    s.to_bits(),
                    want.to_bits()
                );
            }
            fails += 1;
        }
    }
    for (i, (&x, &c)) in inputs.iter().zip(out_cos.iter()).enumerate() {
        let want = cos_wide(x);
        if c.to_bits() != want.to_bits() {
            if fails < 8 {
                println!("COS MISMATCH idx={i} x={x:e}: got {c:e} want {want:e}");
            }
            fails += 1;
        }
    }
    println!(
        "bit-exactness over {} patterns (+256 pinned exponents): {}",
        inputs.len(),
        if fails == 0 { "PASS" } else { "FAIL" }
    );
    if fails > 0 {
        return;
    }

    println!("\n== timing (this process only; run serialized via jm bench) ==");
    timed_compare();
}

fn timed_compare() {
    // Mirror quickbench's throughput methodology: L1-resident array, min over
    // reps. Three distributions: quickbench-style Band::Two magnitudes,
    // chain-range huge inputs (the wide tier's actual reason to exist), and
    // everything-including-tiny.
    for (label, lo, hi) in [
        ("band-two ", 4u64, 5u64),
        ("huge     ", 96u64, 255u64),
        ("all-exps ", 0u64, 260u64),
    ] {
        let n = 1 << 14;
        let mut xs = vec![0f32; n];
        fill(&mut xs, 99);
        for (i, x) in xs.iter_mut().enumerate() {
            let e = lo + ((i as u64).wrapping_mul(2654435761) % (hi - lo));
            *x = f32::from_bits((e as u32) << 23 | (x.to_bits() & 0x007f_ffff));
        }
        let mut o1 = vec![0f32; n];
        let mut o2 = vec![0f32; n];

        let mut ts_scalar = Vec::new();
        let mut ts_x8 = Vec::new();
        const REPS: usize = 200;
        for _ in 0..REPS {
            let t = std::time::Instant::now();
            for (o, &x) in o1.iter_mut().zip(xs.iter()) {
                *o = sin_wide(x);
            }
            black_box(&o1);
            ts_scalar.push(t.elapsed());
            let t = std::time::Instant::now();
            unsafe {
                sin_wide_x8_slice(&xs, &mut o2);
            }
            black_box(&o2);
            ts_x8.push(t.elapsed());
        }
        report("sin_wide", &xs, &o1, &o2);
        let bs = ts_scalar.iter().min().unwrap().as_nanos() as f64 / n as f64;
        let bx = ts_x8.iter().min().unwrap().as_nanos() as f64 / n as f64;
        println!(
            "{label} sin_wide autovec {bs:.2} ns/elem   x8 {bx:.2} ns/elem   speedup {:.2}x",
            bs / bx
        );

        ts_scalar.clear();
        ts_x8.clear();
        for _ in 0..REPS {
            let t = std::time::Instant::now();
            for (o, &x) in o1.iter_mut().zip(xs.iter()) {
                *o = cos_wide(x);
            }
            black_box(&o1);
            ts_scalar.push(t.elapsed());
            let t = std::time::Instant::now();
            unsafe {
                cos_wide_x8_slice(&xs, &mut o2);
            }
            black_box(&o2);
            ts_x8.push(t.elapsed());
        }
        report("cos_wide", &xs, &o1, &o2);
        let bs = ts_scalar.iter().min().unwrap().as_nanos() as f64 / n as f64;
        let bx = ts_x8.iter().min().unwrap().as_nanos() as f64 / n as f64;
        println!(
            "{label} cos_wide autovec {bs:.2} ns/elem   x8 {bx:.2} ns/elem   speedup {:.2}x",
            bs / bx
        );
    }
}

fn report(name: &str, xs: &[f32], a: &[f32], b: &[f32]) {
    let mut first = None;
    let mut cnt = 0usize;
    for (i, (x, y)) in a.iter().zip(b.iter()).enumerate() {
        if x.to_bits() != y.to_bits() {
            cnt += 1;
            if first.is_none() {
                first = Some(i);
            }
        }
    }
    match first {
        None => println!("{name}: outputs bit-identical"),
        Some(i) => {
            let j = i;
            println!(
                "{name}: {} DIFFS, first at {j}: x={:e} ({:#x}) scalar={:e} x8={:e}",
                cnt, xs[j], xs[j].to_bits(), a[j], b[j]
            );
        }
    }
}
