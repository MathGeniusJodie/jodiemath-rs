# faster? f32 math functions
Attempting to provide faster implementations of common f32 math functions with a similar level of accuracy to the standard library.
In the future I want to also provide perfectly rounded versions and less accurate but faster versions.
Currently only the cbrt and exp2 functions are significantly faster, but I have a lot of ideas for improvements.

# precision (average ulp)
```
         jodie cbrt error: 0.574
           std cbrt error: 0.094
jodie accurate cbrt error: 0

jodie exp2 error: 1.27
  std exp2 error: 0.985

jodie  cos error: 6.37
  std  cos error: 6.13

jodie  sin error: 9.715
  std  sin error: 7.33

jodie log2 error: 0.055
  std log2 error: 0
```

# benchmarks
```
Run on i7-12700H

Latency (lower is better)
              | jodie  | std    | perf improvement
--------------|--------|--------|-----------------
         cbrt | 2.4 ns | 5.4 ns | 2.3x
accurate cbrt | 7.3 ns | 5.4 ns | 0.7x
          cos | 4.7 ns | 6.7 ns | 1.4x
         exp2 | 1.5 ns | 4.9 ns | 3.3x
         log2 | 1.7 ns | 3.8 ns | 2.2x
          sin | 6.6 ns | 6.9 ns | 1.0x

Throughput (higher is better) (single thread)
              | jodie    | std      | perf improvement
--------------|----------|----------|-----------------
         cbrt | 3.2 GB/s | 0.5 GB/s | 6.3x
accurate cbrt | 1.8 GB/s | 0.5 GB/s | 3.6x
          cos | 7.5 GB/s | 1.1 GB/s | 6.8x
         exp2 | 6.8 GB/s | 1.4 GB/s | 4.9x
         log2 | 6.6 GB/s | 1.9 GB/s | 3.5x
          sin | 5.8 GB/s | 1.1 GB/s | 5.2x

```

# todo:
- add inverse trig functions
- add tan()
- less accurate but faster versions
- perfectly rounded versions