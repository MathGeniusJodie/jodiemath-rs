# faster? f32 math functions
Attempting to provide faster implementations of common f32 math functions with a similar level of accuracy to the standard library.
There are also perfectly rounded variants and faster-but-sloppier variants.

cbrt, cbrt_accurate, log2 and exp2_checked are full-range correct: negatives,
denormals, zero, inf and nan are all handled. exp2's default is the fast
unchecked version, only valid for x in [-126, 128); exp2_checked handles the
full range (overflow to inf, denormal underflow, nan) for ~2.5 ns extra
latency. cbrt_accurate is within 1 ulp on every input tested (exhaustively,
see below) -- **not** perfectly rounded on every input, despite an earlier
claim here: the full 2^32-pattern sweep found a systematic 1-ulp miss at one
specific mantissa, recurring at every octave from ~2^-126 to ~2^127.
**sin/cos: gradual degradation, not a cliff -- sub-ulp (avg) out to ~1e13,
decaying smoothly from there, and always finite for every finite input.**
q = round(x/pi) must be an exact integer for the reduction to land in
[-pi/2, pi/2]; a single f32 q is a binary either-or (exactly right, or off
by a whole integer once |x| crosses q's exact-integer ceiling), which
shifts the residual by a whole multiple of pi and puts the degree-9
polynomial (fit only for [-pi/2, pi/2]) hopelessly outside its domain -- a
relocatable *cliff*, not a slope, no matter how q is rounded (tried both
the classic magic-constant trick and a native hardware `.round()`; the
latter only moves the cliff from ~1.3e7 to ~2.6e7, same shape). Fixed by
giving q a second, small f32 word (qh, ql) -- genuine double-float
precision via `two_prod`/`two_sum` error-free transforms (exact at any
magnitude, unlike the crate's other PI_A..D trick, which needs bounded q).
A cheap clamp on the *reduced residual* (before it reaches the polynomial)
still guarantees the output can never overflow to `inf`, however bad the
reduction gets. Verified
exhaustively over all 2^32 bit patterns (examples/accuracy.rs,
examples/edgecheck.rs): avg ulp stays flat (~0.25-0.35) from small x out
through ~1e13, then climbs gradually (not a jump) through ~1e19 before
plateauing at essentially-uncorrelated-but-finite output for anything
larger. **Real bug found and fixed along the way**: folding cos's -0.5
phase-shift offset directly into the two_prod's dominant term breaks once
that term's own ulp exceeds 1 (|x| ~ 1.68e7) -- adding a fixed 0.5 to an
already-coarse-ulp float just rounds it away, silently corrupting cos's
residual by a whole pi. This is the same class of bug as everything else
in this reduction (a small quantity lost against a big-ulp value) and
likely also affected an earlier, since-removed version of this same
double-float design that was never checked with a cos-specific
magnitude-bucketed sweep (only sin's was; the two aren't symmetric here).
Fixed by folding the offset into the small correction word instead, where
it survives regardless of the dominant term's own precision. A middle
ground that keeps *only* the two dominant tiers (dropping the third,
smallest one) was tried and measured to cost the same as keeping all three
(~130-140 cyc either way) -- there isn't a cheaper partial version once you
need genuine multi-word q precision at all, so the version here just keeps
all three tiers for the best accuracy at no extra cost.

All functions auto-vectorize, it's a hard requirement

# precision (see examples/accuracy.rs)
Fuzz-mode (100M random f32 bit patterns/function); pass `thorough` for an
exhaustive sweep of all 2^32 patterns instead (few minutes, needs --release).
```
                  | jodie avg  | jodie max | std avg | std max
------------------|------------|-----------|---------|--------
             cbrt |    0.326   |     3     |    0    |    0
    cbrt_accurate |    0.000   |     1     |    0    |    0
             exp2 |    0.030   |     1     |  0.000  |    1
     exp2_checked |    0.016   |     1     |  0.000  |    1
             log2 |    0.003   |     3     |  0.000  |    1
    sin (|x|<=1e6)|    0.036   |     2     |  0.000  |    1
    cos (|x|<=1e6)|    0.081   |     3     |  0.000  |    1
    sin (all f32) | (degrades gradually past |x| ~ 1e13 -- see note above; 0.007/1 for std)
    cos (all f32) | (degrades gradually past |x| ~ 1e13 -- see note above; 0.007/1 for std)
```
sin/cos rows are exhaustive (all 2^32 bit patterns); the rest of this table
is the default 100M-sample fuzz mode. Magnitude-bucketed avg/max ulp shows
the actual shape of the degradation -- flat and low for a very long
stretch, then a real but gradual climb, never a sudden jump to garbage and
never `inf` (confirmed by an exhaustive all-2^32-pattern check that no
finite input produces a non-finite output). These buckets are now a
permanent part of `examples/accuracy.rs` (quick-fuzz mode, ~1.3M
samples/bucket) instead of an ad hoc uncommitted script, so the shape claim
stays checkable after future changes instead of just asserted:
```
        range | sin avg |  sin max | cos avg |  cos max
--------------|---------|----------|---------|----------
     |x|<=1e6 |   0.036 |        2 |   0.081 |        3
   [1e7,1e8)  |   0.253 |        2 |   0.262 |        6
  [1e9,1e10)  |   0.252 |       24 |   0.252 |        2
 [1e12,1e13)  |   0.285 |    10535 |   0.282 |     3037
 [1e15,1e16)  |   3.2e8 |    2.4e9 |   9.7e8 |     2.4e9
```

# benchmarks
Run on i5-1145G7, -C target-cpu=native (now set in .cargo/config.toml)
```
Serial latency (dependency chain, examples/quickbench.rs; lower is better)
              | jodie   | std     | improvement
--------------|---------|---------|------------
         cbrt | 12.9 ns | 22.0 ns | 1.7x
cbrt_accurate | 17.0 ns | 22.0 ns | 1.3x
          cos | 30.1 ns | 12.6 ns | 0.4x
         exp2 |  8.5 ns | 13.3 ns | 1.6x
 exp2_checked | 13.3 ns | 13.3 ns | 1.0x
         log2 | 13.1 ns | 14.9 ns | 1.1x
          sin | 29.1 ns | 12.8 ns | 0.4x
```
sin/cos's latency is worse than std again (2026-07-06, later same day:
reinstated a double-float reduction -- see the accuracy note above for why
the cheap single-word version was reverted a second time). Deliberate
trade, same shape as the very first version of this reduction: gradual,
predictable degradation instead of a cliff costs real latency. (Improved
slightly from 32.9/31.7 ns by the same-day critical-path shortening
described below -- still 0.4x at this rounding.)
```
Throughput (independent array evals over [f32; 4096], examples/quickbench.rs; lower is better)
              | jodie    | std     | improvement
--------------|----------|---------|------------
         cbrt | 0.37 ns  | 3.93 ns | 10.6x
cbrt_accurate | 0.67 ns  | 3.93 ns | 5.9x
          cos | 1.52 ns  | 2.92 ns | 1.9x
         exp2 | 0.23 ns  | 3.13 ns | 13.7x
 exp2_checked | 0.48 ns  | 3.13 ns | 6.6x
         log2 | 0.53 ns  | 4.10 ns | 7.8x
          sin | 1.51 ns  | 3.05 ns | 2.0x
```
sin/cos's throughput is back down to ~1.6x std (was ~9-11x with the cheap
single-word-plus-clamp version, ~13-14x with the original
single-word-but-inf-prone version). This is the direct cost of genuine
gradual degradation: a double-float q needs several `two_prod`/`two_sum`
error-free transforms (each 2-6 ops) to stay accurate well past a single
f32's exact-integer range, instead of one cheap magic-constant rounding
trick. Tried to find a cheaper partial version (only the two biggest
cross-term tiers, dropping the smallest) -- measured to cost the *same* as
keeping all three (~130-140 cyc either way, see the mca table below), so
there's no meaningfully cheaper middle ground once multi-word q is needed
at all; kept all three tiers since the third one is free once you're
already paying for the other two. Two vectorization/perf pitfalls to
remember if this reduction is ever revisited (see jodiemath-workflow memory
for the full trail, including the cos pre_offset bug described above):
- a `q as i64` cast for the parity bit doesn't vectorize at all -- Rust's
  float-to-int cast is saturating, so LLVM falls back to a scalar
  convert-with-NaN/range-check per lane (~9x slower than std). Parity here
  uses a `floor`-based "mod 2" instead (`q - 2*(q*0.5).floor()`), the same
  instruction family as the `.round()` already used elsewhere, so it stays
  vectorized (confirmed via `--emit=asm`: 0 scalar convert instructions).
- combining several small correction terms via plain adds first (instead
  of feeding all of them through a sequential compensated-sum loop) cuts
  the summation loop from 7 iterations to 4 with no accuracy cost --
  verified bit-for-bit identical output over the full accuracy.rs sweep.
Two further micro-optimizations landed later the same day, found by
re-running `llvm-mca --bottleneck-analysis` on the *current* code (it now
skews more dependency-chain-bound than the ~39%-resource-pressure region an
earlier attempt at this measured, so shortening the critical path pays off
here where it didn't before -- always re-check, don't assume an old
bottleneck-analysis finding still applies after the surrounding code
changes):
- `reduce_pi`'s first correction merge (`two_sum(s, -e1)`, combining the
  residual right after subtracting the dominant term with that
  subtraction's own two_prod rounding-error) uses the cheaper 3-op
  `quick_two_sum` (Fast2Sum) instead of full 6-op `two_sum`. **Correction to
  an initial claim here**: a first check (only near exact multiples of pi,
  a narrow and unrepresentatively well-conditioned slice) seemed to show the
  `|s|>=|e1|` ordering Fast2Sum needs for *exactness* holds everywhere; a
  broader recheck (uniform random bit patterns, not just near-exact
  multiples of pi) found this is false -- real violations starting around
  `|x| ~ 1e3` and exceeding 80% of samples by `|x| ~ 1e8`+. Kept anyway,
  because Fast2Sum's failure mode when misordered is *bounded* (`e` off by
  up to ~1 ulp of `s`, not unbounded), and re-verified against the metric
  that actually matters -- the full accuracy.rs exhaustive sweep plus a
  magnitude-bucketed sweep against std out to f32::MAX -- shows no
  measurable difference from the full-`two_sum` version at any magnitude.
  Lesson: checking a theoretical invariant (an ordering assumption) is not
  the same as checking the thing that matters (final ulp error); it can be
  technically false while still being practically harmless. The other three
  terms in that loop (`p2`, `p3`, `tier2`) violate this ordering far more
  severely (`p3` in particular is usually exactly 0 but occasionally ~pi,
  comparable to the whole residual, whenever `ql != 0`) and were left on
  full `two_sum`, untested whether quick_two_sum would be harmless there too.
  Tried going further and pre-combining more of the 4 correction terms into
  a shallower tree (a single merge instead of 4 sequential ones, or even
  just pairing off `p3`+`tier2` in parallel with the dominant subtraction)
  -- both made latency *and* throughput worse (e.g. the full-tree version:
  sin 140→149 cyc latency, 7.609→8.545 cyc/elem), despite being exact and
  despite the region being dependency-bound: more total register-live-range
  pressure in the 16-wide vectorized loop apparently outweighs the shorter
  chain. Reverted both; kept only the one verified-safe single swap.
- the final parity combine (`parity(parity(qh) + parity(ql))`) doesn't need
  its 3rd `floor`-based `parity()` call: `parity(qh)` and `parity(ql)` are
  each exactly 0.0 or 1.0, so their sum mod 2 is just whether they differ,
  i.e. `if pq == pl { 0.0 } else { 1.0 }` -- a compare+select instead of a
  4-op floor chain, still branchless/vectorized (confirmed via `--emit=asm`:
  0 scalar convert instructions).
Both changes verified against the full accuracy.rs exhaustive sweep (`sin
|x|<=1e6` 0.036→0.036 avg / 2→2 max ulp, `cos` 0.077→0.077 avg / 2→2 max
ulp -- unchanged) and a magnitude-bucketed sweep confirming the same
gradual (non-cliff) shape past 1e6. Net: sin 140.00→136.00 cyc latency,
7.609→7.109 cyc/elem throughput; cos 144.00→140.00 cyc latency,
6.657→5.850 cyc/elem throughput.

**Further pass, same reduction, later session: revisits the "pairing off
p3+tier2" idea rejected just above, this time getting a real win from it --
the earlier attempt most likely had the same sign bug this session found
and fixed, though its code wasn't kept to confirm directly.** `round_x_over_pi`'s
`two_sum(e0, x*RPI_LO)` (comparable-magnitude operands, same
ordering-not-guaranteed-but-bounded situation as the swap above) moved to
`quick_two_sum` -- free, no measurable accuracy change. Then, in
`reduce_pi`: `p3` and `tier2` are the only two of the four correction terms
that depend on `ql` (round_x_over_pi's last-ready output, needing its whole
chain), while `p1`/`p2`/the dominant `x - p1` subtraction only need `qh`,
ready much earlier -- so combining `p3+tier2` into one value via `two_sum`
runs parallel to the `qh`-only work instead of stacking as two more
sequential merges after it, and the merge of that combined value into the
running residual was further downgraded to `quick_two_sum`. **Real bug hit
while building this**: `two_sum(p3, tier2)` guarantees `p3t + e3t == p3 +
tier2` exactly, so subtracting `(p3+tier2)` from the residual means
subtracting *both* `p3t` and `e3t` -- the first attempt added `e3t` into the
running low-order correction instead of subtracting it, which is a sign
error, not a rounding one. It passed small-`x` spot checks but broke `cos`
badly on the full sweep (`|x|<=1e6` avg 2.51 ulp, max 37M ulp, worst case
right at `x ≈ 2.5π` and `3.5π` -- `cos`'s zero crossings, exactly where `p3`
and `tier2` partially cancel, making `e3t` unexpectedly large instead of the
negligible correction it usually is). Fixed by flipping the sign; re-ran the
full exhaustive sweep before trusting it further. This is very likely the
same class of mistake behind the "made things worse" result recorded just
above for the identical-sounding idea -- a sign error wouldn't necessarily
look like a correctness bug in an mca run (mca has no notion of numerical
correctness, only codegen), so a broken version could plausibly still
compile to *legitimately worse* code for unrelated reasons and get reverted
without the sign bug itself ever being noticed. Not provable without the old
code, which wasn't kept, but worth recording as a caution: an idea "already
tried and found not to help" is only as trustworthy as the correctness of
the attempt, and llvm-mca's cycle counts can't tell you whether the code
being measured was actually right.
`llvm-mca --bottleneck-analysis` on `cos_throughput` (this session, before
these changes) showed a roughly even split -- 53.5% resource pressure
(`ICXPort0`/`ICXPort1`, the fma/mul ports) vs. 61.7% register/data
dependencies -- slightly skewed towards the dependency chain, consistent
with a chain-shortening restructure (rather than an op-count cut) paying
off here. Verified against the full exhaustive accuracy.rs sweep: `sin
|x|<=1e6` unchanged (0.036 avg / 2 max ulp); `cos |x|<=1e6` 0.077→0.081 avg
/ 2→3 max ulp -- a small, bounded cost from the two additional
`quick_two_sum` downgrades, still >10x under the 1-ulp-average budget. A
magnitude-bucketed sweep (see the precision table above, now a permanent
part of examples/accuracy.rs) confirms the same gradual, non-cliff shape
survives out past 1e15. Net: sin 136.00→124.00 cyc latency (-8.8%),
7.109→6.421 cyc/elem throughput (-9.7%); cos 140.00→128.00 cyc latency
(-8.6%), 5.850→5.354 cyc/elem throughput (-8.5%).
The throughput gap vs std comes almost entirely from vectorization: std's
functions have branches, so LLVM can't vectorize loops that call them.
Absolute ns swing session-to-session with CPU thermal state (the laptop
throttles up to ~2.5x mid-session) -- only trust jodie-vs-std ratios measured
in the same run. examples/mca.rs gives a thermal-noise-free second opinion in
cycles instead of ns:
```
theoretical cost from llvm-mca (-mcpu=native, 100 iterations)
                    | latency (cyc) | throughput (cyc)
--------------------|----------------|------------------
cbrt                |          35.06 |             1.629
cbrt_accurate       |          59.06 |             3.129
exp2                |          35.00 |             0.841
exp2_checked        |          43.06 |             1.399
log2                |          34.23 |             1.556
sin                 |         124.00 |             6.421
cos                 |         128.00 |             5.354
```

# tools
- `cargo run --release --example accuracy [thorough] [filter]` - avg/max ulp against an f64 reference.
  Default mode fuzzes 100M random f32 bit patterns per function (a few seconds); `thorough` exhaustively
  sweeps all 2^32 bit patterns instead (every denormal, every NaN payload, both signs -- a few minutes,
  multi-threaded, refuses to run in a debug build)
- `cargo run --release --example quickbench [filter]` - latency (serial dependency chain) + throughput, min of 7 reps
- `cargo run --release --example edgecheck` - bit-exact checks of edge cases (0, -0, denormals, inf, nan, domain boundaries)
- `cargo run --release --example tune` - coordinate-descent ulp tuning of polynomial coefficients
- `cargo run --release --example mca` - theoretical latency/throughput straight from llvm-mca's scheduler
  model for the host CPU (requires `llvm-mca` on PATH). No wall-clock timing, so no thermal-throttling
  noise, and much faster to iterate on than quickbench; see examples/mca_target.rs for the marker
  functions it analyzes and why each region is built the way it is (llvm-mca has no branch predictor,
  so branchy edge-case handling has to be routed around, not just measured through).

# todo:
- do principled and thourough analysis of dependency chains and rounding errors to find optimizations
- add inverse trig functions
- add tan()
- perfectly rounded versions
