use criterion::{criterion_group, criterion_main, Criterion};
use jodiemath_rs::{cbrt, cos, exp2, log_2, sin};
use std::hint::black_box;

fn bench_cbrt(c: &mut Criterion) {
    c.bench_function("jodie exp2", |b| {
        b.iter(|| {
            let input: [f32; 8] = black_box([12345.0; 8]);
            let mut out: [f32; 8] = [0.0; 8];
            for i in 0..8 {
                out[i] = exp2(input[i]);
            }
            for i in 0..8 {
                black_box(out[i]);
            }
        })
    });
    c.bench_function("std exp2", |b| {
        b.iter(|| {
            let input: [f32; 8] = black_box([12345.0; 8]);
            let mut out: [f32; 8] = [0.0; 8];
            for i in 0..8 {
                out[i] = input[i].exp2();
            }
            for i in 0..8 {
                black_box(out[i]);
            }
        })
    });
    c.bench_function("jodie cbrt", |b| {
        b.iter(|| {
            let input: [f32; 8] = black_box([12345.0; 8]);
            let mut out: [f32; 8] = [0.0; 8];
            for i in 0..8 {
                out[i] = cbrt(input[i]);
            }
            for i in 0..8 {
                black_box(out[i]);
            }
        })
    });
    c.bench_function("std cbrt", |b| {
        b.iter(|| {
            let input: [f32; 8] = black_box([12345.0; 8]);
            let mut out: [f32; 8] = [0.0; 8];
            for i in 0..8 {
                out[i] = input[i].cbrt();
            }
            for i in 0..8 {
                black_box(out[i]);
            }
        })
    });
    c.bench_function("jodie cos", |b| {
        b.iter(|| {
            let input: [f32; 8] = black_box([12345.0; 8]);
            let mut out: [f32; 8] = [0.0; 8];
            for i in 0..8 {
                out[i] = cos(input[i]);
            }
            for i in 0..8 {
                black_box(out[i]);
            }
        })
    });
    c.bench_function("std cos", |b| {
        b.iter(|| {
            let input: [f32; 8] = black_box([12345.0; 8]);
            let mut out: [f32; 8] = [0.0; 8];
            for i in 0..8 {
                out[i] = input[i].cos();
            }
            for i in 0..8 {
                black_box(out[i]);
            }
        })
    });
    c.bench_function("jodie log_2", |b| {
        b.iter(|| {
            let input: [f32; 8] = black_box([12345.0; 8]);
            let mut out: [f32; 8] = [0.0; 8];
            for i in 0..8 {
                out[i] = log_2(input[i]);
            }
            for i in 0..8 {
                black_box(out[i]);
            }
        })
    });
    c.bench_function("std log_2", |b| {
        b.iter(|| {
            let input: [f32; 8] = black_box([12345.0; 8]);
            let mut out: [f32; 8] = [0.0; 8];
            for i in 0..8 {
                out[i] = input[i].log2();
            }
            for i in 0..8 {
                black_box(out[i]);
            }
        })
    });
    c.bench_function("jodie sin", |b| {
        b.iter(|| {
            let input: [f32; 8] = black_box([12345.0; 8]);
            let mut out: [f32; 8] = [0.0; 8];
            for i in 0..8 {
                out[i] = sin(input[i]);
            }
            for i in 0..8 {
                black_box(out[i]);
            }
        })
    });
    c.bench_function("std sin", |b| {
        b.iter(|| {
            let input: [f32; 8] = black_box([12345.0; 8]);
            let mut out: [f32; 8] = [0.0; 8];
            for i in 0..8 {
                out[i] = input[i].sin();
            }
            for i in 0..8 {
                black_box(out[i]);
            }
        })
    });
}

criterion_group!(benches, bench_cbrt);
criterion_main!(benches);
