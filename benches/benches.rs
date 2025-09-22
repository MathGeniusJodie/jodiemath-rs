use criterion::{criterion_group, criterion_main, Criterion};
//use std::hint::black_box;
use criterion::black_box;
use jodiemath_rs::{cbrt, cos, exp2, log_2, sin};

fn bench_cbrt(c: &mut Criterion) {
    c.bench_function("jodie cbrt", |b| b.iter(|| cbrt(black_box(12345.0))));
    c.bench_function("std cbrt", |b| b.iter(|| black_box(12345.0f32).cbrt()));
    c.bench_function("jodie cos", |b| b.iter(|| cos(black_box(12345.0))));
    c.bench_function("std cos", |b| b.iter(|| black_box(12345.0f32).cos()));
    c.bench_function("jodie exp2", |b| b.iter(|| exp2(black_box(12345.0))));
    c.bench_function("std exp2", |b| b.iter(|| black_box(12345.0f32).exp2()));
    c.bench_function("jodie log_2", |b| b.iter(|| log_2(black_box(12345.0))));
    c.bench_function("std log_2", |b| b.iter(|| black_box(12345.0f32).log2()));
    c.bench_function("jodie sin", |b| b.iter(|| sin(black_box(12345.0))));
    c.bench_function("std sin", |b| b.iter(|| black_box(12345.0f32).sin()));
}

criterion_group!(benches, bench_cbrt);
criterion_main!(benches);
