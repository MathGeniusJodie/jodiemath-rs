
use criterion::{criterion_group, criterion_main, Criterion};
use std::hint::black_box;
use jodiemath_rs::cbrt;

fn cbrt_bench() {
    fn bench_cbrt(c: &mut Criterion) {
        c.bench_function("jodie cbrt", |b| b.iter(|| cbrt(black_box(12345.0))));
        c.bench_function("std cbrt", |b| b.iter(|| black_box(12345.0f32).cbrt()));
    }
    criterion_group!(benches, bench_cbrt);
    criterion_main!(benches);
}