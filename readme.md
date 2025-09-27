# faster? f32 math functions
Attempting to provide faster implementations of common f32 math functions with a similar level of accuracy to the standard library.
In the future I want to also provide perfectly rounded versions and less accurate but faster versions.
Currently only the cbrt and exp2 functions are significantly faster, but I have a lot of ideas for improvements.

# precision (average ulp)
```
jodie cbrt error: 0.646
  std cbrt error: 0.094

jodie exp2 error: 1.27
  std exp2 error: 0.985

jodie  cos error: 6.37
  std  cos error: 6.13

jodie  sin error: 9.715
  std  sin error: 7.33

jodie log2 error: 0.037
  std log2 error: 0
```

# benchmarks
```
    | jodie   | std     | perf improvement
----|---------|---------|-----------------
cbrt| 11.5 ns | 60.3 ns | 5.2x
 cos| 21.2 ns | 27.0 ns | 1.3x
exp2|  8.4 ns | 22.6 ns | 2.7x
log2|  8.8 ns | 14.3 ns | 1.6x
 sin| 38.8 ns | 25.3 ns | 0.7x (slower)
```

# todo:
- fix sin() performance
- less accurate but faster versions
- add inverse trig functions
- add tan()
- perfectly rounded versions