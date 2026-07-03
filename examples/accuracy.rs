// Accuracy harness: avg + max ULP over dense strided sampling of each
// function's domain. Reference values computed in f64 and rounded to f32.
use jodiemath_rs::*;

fn ulp_diff(a: f32, b: f32) -> u64 {
    fn ord(x: f32) -> i64 {
        let b = x.to_bits();
        if b & 0x8000_0000 != 0 {
            -((b & 0x7fff_ffff) as i64)
        } else {
            b as i64
        }
    }
    if a.is_nan() || b.is_nan() {
        return if a.is_nan() == b.is_nan() { 0 } else { u64::MAX };
    }
    (ord(a) - ord(b)).unsigned_abs()
}

struct Stats {
    sum: u64,
    max: u64,
    worst_x: f32,
    n: u64,
}

fn measure(
    lo: f32,
    hi: f32,
    stride: u32,
    f: impl Fn(f32) -> f32,
    reference: impl Fn(f64) -> f64,
) -> Stats {
    // bit-pattern walk only goes upward: negative ranges must be swept as
    // mirrored positives (negate inside the closures), or they'd silently
    // contribute zero samples
    assert!(lo > 0.0 && lo <= hi, "measure() needs 0 < lo <= hi");
    let (blo, bhi) = (lo.to_bits(), hi.to_bits());
    let mut s = Stats { sum: 0, max: 0, worst_x: 0.0, n: 0 };
    let mut process = |x: f32| {
        let r = reference(x as f64) as f32;
        let d = ulp_diff(f(x), r);
        s.sum += d;
        if d > s.max {
            s.max = d;
            s.worst_x = x;
        }
        s.n += 1;
    };
    // walk bit patterns from lo..hi (same sign assumed), plus mirrored negatives if lo<0 handled by caller
    let mut b = blo;
    while b <= bhi {
        process(f32::from_bits(b));
        b = b.wrapping_add(stride);
        if b < stride { break; }
    }
    s
}

fn report(name: &str, s: &Stats) {
    println!(
        "{:24} avg ulp {:>10.4}  max ulp {:>8}  worst x {:e} ({} samples)",
        name,
        s.sum as f64 / s.n as f64,
        s.max,
        s.worst_x,
        s.n
    );
}

fn combine(a: Stats, b: Stats) -> Stats {
    Stats {
        sum: a.sum + b.sum,
        max: a.max.max(b.max),
        worst_x: if a.max >= b.max { a.worst_x } else { b.worst_x },
        n: a.n + b.n,
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let filter = args.get(1).map(|s| s.as_str()).unwrap_or("");
    let run = |n: &str| filter.is_empty() || n.contains(filter);

    // positive-domain functions: sample all positive normal floats with stride
    if run("cbrt") {
        // full positive normal range, ~2M samples
        let s = measure(1e-38, 3e38, 1024, cbrt, |x| x.cbrt());
        report("cbrt (+)", &s);
        let s = measure(1e-38, 3e38, 1024, |x| cbrt(-x), |x| (-x).cbrt());
        report("cbrt (-)", &s);
        let s = measure(1e-38, 3e38, 1024, cbrt_accurate, |x| x.cbrt());
        report("cbrt_accurate (+)", &s);
        let s = measure(1e-38, 3e38, 1024, cbrt_throughput, |x| x.cbrt());
        report("cbrt_throughput (+)", &s);
        let s = measure(1e-38, 3e38, 1024, cbrt_fast, |x| x.cbrt());
        report("cbrt_fast (+)", &s);
        let s = measure(1e-38, 3e38, 1024, |x| x.cbrt(), |x| x.cbrt());
        report("std cbrt (+)", &s);
    }
    if run("log") {
        let s = measure(1e-38, 3e38, 1024, log_2, |x| x.log2());
        report("log_2", &s);
        let s = measure(1e-38, 3e38, 1024, |x| x.log2(), |x| x.log2());
        report("std log2", &s);
    }
    if run("exp") {
        // negative ranges must be swept as mirrored positives: measure()
        // walks bit patterns upward, so a (negative, negative) range has
        // blo > bhi and silently contributes zero samples
        // unchecked exp2 domain: [-126, 128) (normal results only)
        let s = combine(
            measure(1e-30, 127.9, 512, exp2, |x| x.exp2()),
            measure(1e-30, 126.0, 512, |x| exp2(-x), |x| (-x).exp2()),
        );
        report("exp2", &s);
        // checked variant additionally covers the denormal-result range
        let s = combine(
            combine(
                measure(1e-30, 127.9, 512, exp2_checked, |x| x.exp2()),
                measure(1e-30, 126.0, 512, |x| exp2_checked(-x), |x| (-x).exp2()),
            ),
            measure(126.0, 149.9, 64, |x| exp2_checked(-x), |x| (-x).exp2()),
        );
        report("exp2_checked", &s);
        let s = combine(
            measure(1e-30, 127.9, 512, |x| x.exp2(), |x| x.exp2()),
            measure(1e-30, 126.0, 512, |x| (-x).exp2(), |x| (-x).exp2()),
        );
        report("std exp2", &s);
    }
    if run("sin") {
        for (name, hi, stride) in [("sin [0,pi/4]", 0.785398_f32, 64u32), ("sin [0,10]", 10.0, 64), ("sin [0,1000]", 1000.0, 64)] {
            let s = combine(
                measure(1e-30, hi, stride, sin, |x| x.sin()),
                measure(1e-30, hi, stride, |x| sin(-x), |x| (-x).sin()),
            );
            report(name, &s);
        }
        let s = combine(
            measure(1e-30, 1000.0, 64, |x| x.sin(), |x| x.sin()),
            measure(1e-30, 1000.0, 64, |x| (-x).sin(), |x| (-x).sin()),
        );
        report("std sin [0,1000]", &s);
    }
    if run("cos") {
        for (name, hi, stride) in [("cos [0,pi/4]", 0.785398_f32, 64u32), ("cos [0,10]", 10.0, 64), ("cos [0,1000]", 1000.0, 64)] {
            let s = combine(
                measure(1e-30, hi, stride, cos, |x| x.cos()),
                measure(1e-30, hi, stride, |x| cos(-x), |x| (-x).cos()),
            );
            report(name, &s);
        }
        let s = combine(
            measure(1e-30, 1000.0, 64, |x| x.cos(), |x| x.cos()),
            measure(1e-30, 1000.0, 64, |x| (-x).cos(), |x| (-x).cos()),
        );
        report("std cos [0,1000]", &s);
    }
}
