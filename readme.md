# faster? f32 math functions
Attempting to provide faster implementations of common f32 math functions with a similar level of accuracy to the standard library.
In the future I want to also provide perfectly rounded versions and less accurate but faster versions.
Currently only the cbrt and exp2 functions are significantly faster, but I have a lot of ideas for improvements.

# precision (average ulp)
```
         jodie cbrt error: 0.83
           std cbrt error: 0
jodie accurate cbrt error: 0

jodie exp2 error: 0.315
  std exp2 error: 0.105

jodie  cos error: 0.45
  std  cos error: 0.08

jodie  sin error: 0.3
  std  sin error: 0.16

jodie log2 error: 0.16
  std log2 error: 0.09
```

# benchmarks
```
Run on i5-1145G7

Latency (lower is better)
              | jodie  | std    | improvement
--------------|--------|--------|------------
         cbrt | 2.4 ns | 5.4 ns | 2.3x
cbrt_accurate | 7.3 ns | 5.4 ns | 0.7x
          cos | 4.7 ns | 6.7 ns | 1.4x
         exp2 | 2.1 ns | 5.4 ns | 2.6x
         log2 | 1.7 ns | 3.8 ns | 2.2x
          sin | 5.3 ns | 6.5 ns | 1.2x

Throughput (higher is better) (single thread)
              | jodie    | std      | improvement
--------------|----------|----------|------------
         cbrt | 4.7 GB/s | 0.7 GB/s | 6.7x
cbrt_accurate | 2.7 GB/s | 0.7 GB/s | 3.8x
          cos | 4.1 GB/s | 0.6 GB/s | 6.8x
         exp2 | 6.6 GB/s | 0.7 GB/s | 3.9x
         log2 | 6.0 GB/s | 1.0 GB/s | 6.0x
          sin | 3.7 GB/s | 0.6 GB/s | 6.1x

```

# todo:
- do principled and thourough analysis of dependency chains and rounding errors to find optimizations
- investigate higher order for last iter on accurate cbrt
- add inverse trig functions
- add tan()
- perfectly rounded versions