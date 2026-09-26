# graveyard.md

Everything that has been tried: shipped, rejected, or closed, with the numbers that decided it. Open ideas
are in `IDEAS.md`. Check here before starting on a function; re-running an old negative result is the most
common way to waste a session. "Rejected" means measured and rejected, not guessed.

Speed figures are from the Intel i5-1145G7 unless an entry says otherwise (the crate now runs on Zen 5).
Each entry starts with **Shipped.**, **Rejected.**, **Closed.** (nothing to do) or **Finding.** (a
measurement or diagnosis). Within a section, entries are grouped by function
or topic, oldest first.

Two rules that come up again and again:

- **A rejection may only reject one implementation.** Three ideas here were rejected, re-derived from
  scratch later, and shipped (`exp10_checked`'s round-based reduction, the `ln_normal` fma fold, the
  `pre_offset` dead-add removal). If a failure came from one edge case or a scheduling artifact and the
  prize was large, re-derive it rather than resurrecting the old code.
- **Disproving an idea's reason doesn't disprove the change.** One transform sat unused for a year because
  its stated reason ("moves work to integer ports") was measured false; the transform itself was a real
  3-to-2-op win at nine sites for a different reason.

## Method and tooling

### Accuracy gates: check a composite over its own domain
**Shipped.** `gelu` inherited `erfc`'s `|x| <= 10√2` accuracy screen even though its saturated tails can be checked against the f64 reference. A second screen over all finite f32 inputs exposed non-finite results in 10.03% of samples from the old implementation (worst sampled input `-8.74e22`); `edgecheck` had covered infinities and zeros but not large finite inputs. An oracle-kernel substitution measures how much the kernel contributes, not where a composite's worst error lies. `worst_corpus` checks movement, not correctness: its old golden data contained three `gelu` NaNs, so inspect surprising values before blessing them.

### accuracy harness: probe outside chosen bounds with a sound reference
**Finding.** Uniform independent inputs missed `clog`'s near-unit-circle failure: the standing sweep reported max 3 ulp, while a manifold probe found 4096. Its `hypot(...).ln()` reference rounded the magnitude before taking a logarithm, losing the precision needed near magnitude one; the replacement formed squared magnitude exactly with f64 squares and `two_sum`. A separate bounds audit found `sinc`, `wrap_pi`, and `cexp` still within their recorded maxima when widened, while degree trig and unscaled norm bounds were genuinely load-bearing. A domain gate, sampling distribution, or reference must be checked in the difficult region itself; wider sampling alone cannot repair an inaccurate reference.

### Accuracy harness: stale example binaries and sampled two-argument maxima
**Finding.** `cargo build --release` does not rebuild `examples/`; an unchanged `accuracy` binary once made a changed polynomial appear bit-identical. Build with `cargo build --release --example accuracy` or `cargo run`, and check the binary timestamp when a change that must affect results appears to do nothing. Two-argument functions have no exhaustive sweep: `atan2` and `atan2_unchecked` each reached max 4 rather than the published 3 in ordinary 10M-sample repeats. Record their maxima over at least three runs, not a single lucky sample.

### Accuracy screening: sample near the failure manifold
**Finding.** Uniform random sampling missed the largest errors in wide-range argument reduction: `sin_checked`'s `[1e8,1e10)` max varied between 2, 3 and 10 ulp in quick runs, while exhaustive enumeration found 47; three bands had true maxima of 6, 47 and 67 where sampling reported 2, 3 and 21. A targeted 11.5M-point probe within four f32 ulp of trigonometric zeros also exposed parity errors that four million random samples per band missed. Averages were much more stable; treat sampled maxima as lower bounds, and enumerate narrow bands or target the cancellation manifold.

### Accuracy sweeps: distinguish sample noise from changed code
**Finding.** Quick fuzz is unseeded: `cos2pi` appeared to regress from average/max 0.0737/51472 to 0.1172/3294199 despite byte-identical assembly; unchanged paths for `norm_cdf` and `gelu` also moved. The quick average can drift about 0.0001, and near-zero cancellation can move maxima dramatically. Diff the affected assembly and use exhaustive or reproducible samples before attributing an accuracy change to an edit.

### Accuracy sweeps: rare maxima need exhaustive evidence
**Finding.** A 400k-samples-per-octave screen missed all seven 6-ulp inputs among roughly 185 million central-branch inputs of the now-removed `dawson`; it reported max 5 where exhaustive measurement found 6. Similarly, sampled `atan` had max 3, but exhaustive measurement found max 4 at `x = 1.0220603` inside the reflected `|x| >= 1` branch; after the numerator peel the worst input moved to `0.9360248`. Sampled two-argument maxima need repeated runs and should not replace published maxima without stronger evidence.

### Assembly measurements: regenerate the `.s` file before scoring it
**Finding.** `cargo build --release --example mca_target` does not emit assembly, and `cargo rustc` can skip an unchanged example because Cargo does not fingerprint the raw `--emit=asm` argument. Use `touch examples/mca_target.rs` followed by `cargo rustc --release --example mca_target -- --emit=asm -C debuginfo=0`, and check the resulting `.s` timestamp before running `mca_region.py`. Stale assembly once made a behavior-changing edit appear byte-identical on all five measurement rungs.

### Binary mca targets: keep the expensive operand live in the loop
**Shipped.** `xlogy(x,y)` and `xlog1py(x,y)` benchmarks held `y` invariant, letting LLVM hoist the entire logarithm; passing the changing `x` as both operands fixed them. `xlogy`'s reported latency/throughput changed from 21.57/3.506 to 58.94/1.612 cycles, near `ln`'s 56.86/1.614; `xlog1py` changed from 9.03/1.858 to 57.23/2.301, near `log1p`'s 56.09/2.335. Packed/scalar instructions in the latter region changed from 12/61 to 73/18, clearing a false `codegen_check` failure caused by the hoisted scalar divide. A composite apparently cheaper than its expensive primitive warrants inspection for hoisting or elimination; the other invariant-argument targets were checked and unaffected.

### Branch cutoffs: compare errors on both sides after changing a branch
**Finding.** The original cutoffs for `sinh` (0.5), `tanh` (0.25), `expm1` (0.5), `asin` (0.25), and `erf` (0.28) matched branch-error crossovers. An `asin_small` refit later justified moving asin's cutoff to 0.27, so the earlier result was not permanent. After `sinh_small` changed, sinh's worst inputs still lay far from 0.5 (about 3.16, 0.85, and -7.32 for its variants), leaving no reason to move its boundary.

### Branchy benchmarks: check the input band before trusting a cost estimate
**Finding.** In the `tand` experiment, llvm-mca's instruction and uOp counts predicted a throughput regression, while its cycle estimate predicted a gain; paired wall-clock runs found throughput improvements of 2.6% on the stock small-input band and 9.5% over a full period. The stock band never reached the reflected arm, so a representative band was essential. Scalar llvm-mca latency was not comparable with branch-free rows because the new region contained 448 jumps; measure both arms and use paired runs when branch selection matters.

### Coefficient probes: keep evaluation order and hard regions faithful
**Finding.** A cbrt coefficient search reduced max 3→2 on a 12-octave sample, but full `cbrt` retained max 3 and worsened avg 0.281→0.303 because its hardest input was at the tiny/denormal rescale boundary (`1.3057394e-38`). An apparent `erf_poly` max 5→4 gain vanished when the probe used the shipped Estrin chain rather than Horner; dense max stayed 5 and avg moved by at most 0.02%. Max-capped weighted fitting is preferable when average error matters: a plain minimax cbrt correction reduced max 3→2 but increased avg error 43%.

### Combine-weighted LP fits: avoid weights that vanish at the hard inputs
**Rejected.** Weighting coefficients by sensitivity of the final combine helped `erf_poly`, but worsened `asin_poly` (grid max 5→6, avg 0.887→1.390; dense avg 0.05681→0.06438, max unchanged) and the `erfc` rational numerator (max 109→129, avg -4.6%). Their weights vanish near the difficult regions, respectively `a→1` and large `xa`, allowing error to move there. A sigmoid fit's theoretical 13.7% gain became only avg 0.11482→0.11465 with max 4 unchanged. Tanh's `exp2int` range spans ~76 orders of magnitude and defeated normalized HiGHS; its error is also split about evenly between poly and final combine.

### Coverage: add real sweeps where edge tests are insufficient
**Finding.** An audit found `sqrt1pm1` had only edge pins: an exhaustive 3.20-billion-input sweep measured average 0.2053, max 2 ulp, using `x/(sqrt(1+x)+1)` as the cancellation-free reference. The now-removed `wrap_pi` needed a two-word f64 tau reference; a single-word reference falsely reported max 66 rather than 41 ulp. Check reference precision when the implementation itself performs compensated reduction.

### Denormal audit: distinguish subnormal inputs from subnormal results
**Finding.** The denormal-input identity region was clean: `sin`, `tan`, `asin`, `atan`, `sinh`, `asinh`, `atanh`, `expm1`, `expm1_checked`, `log1p` and `softsign` had zero relative error; `erf`, `sinpi`, `gelu` and `silu` were within one representable step at the bottom. For normal inputs yielding subnormal results, `exp2_checked`, `exp_checked`, `exp10_checked` and `erfc` did not flush prematurely; `sigmoid` flushed from about −88.38 although true zero begins at −103.97, a 15.6-wide gap. Count a flush only when the correctly rounded f32 reference is a *nonzero* subnormal: an f64 result below `MIN_POSITIVE` may already round to zero. Do not sweep unchecked functions outside their documented domains.

### denormal audit: hand-picked function lists omitted whole ranges
**Finding.** Extending `denormal_audit`'s hand-maintained case and width lists exposed complete premature denormal flushing in `silu`, `softplus`, and `logsigmoid`: respectively 20.28, 16.97, and 16.97 units of input early. The existing audit had reported only `sigmoid` and `norm_pdf`; its width measurements also omitted functions with no `last_ok_x` inside the sweep. Coverage lists are part of the test's contract and need an explicit completeness check.

### EFT toolkit: audit bounds before publishing internal helpers
**Shipped.** `two_prod`, `two_sum`, `quick_two_sum` and `mulsign` became public; the pre/post `mca_target.s` diff found bit-identical inlined assembly. `two_prod` is exact only when `2^-102 <= |a*b| <= f32::MAX`: interpreting its old “no overflow” note as a guarantee for every input would have been wrong. The weaker `2^-103` lower bound still allowed about 21,000 failures among 42 million in-range pairs. `two_sum` handles finite non-overflowing pairs including subnormals; `mulsign` matched xor-of-sign-bits exhaustively. `eft_contract_check.rs` tests 40 million mixed pairs and all 2^32 bit patterns of one `mulsign` operand against eight other values.

### erfcx throughput: arbitrate contradictory simulator metrics with alternating binaries
**Finding.** Removing an `erfcx` clamp reduced instructions 124→123, uOps 143→140, Block RThroughput 39→38 and latency 4.3%, yet simulated throughput worsened 3.7%. Four alternating A/B rounds of prebuilt binaries favored the clamp-free version in all four throughput comparisons (0.858/0.911/0.889/0.942 versus 0.952/0.952/0.984/0.953 ns/op). The recorded simulator column was retained but is not evidence of a real slowdown; alternating runs matter because an unchanged control drifted substantially within the session.

### Error profiles: locate cancellation instead of reading only max ulp
**Shipped.** `error_profile.rs` records ulp histograms and per-magnitude-band avg/max/worst input. `expm1` averaged 1.116 ulp across its seam and 0.877 just above (max 5 and 6), versus about 0.35 elsewhere; `exp_m1_over_x` was 1.033/1.020 (max 5/6) versus about 0.39; `tanh` was 0.977/0.831 (max 6/4) versus about 0.40–0.50. The direct `exp(x)-1` branch amplifies error by `exp(x)/(exp(x)-1)`: 2.541 at 0.5, 1.582 at 1 and 1.157 at 2. `tanh` feeds `expm1(2x)`, so its seam at 0.25 has the same 2.541 factor. Profiling identifies a narrow cancellation band that a whole-domain avg/max obscures.

### Exponent fields: vector integer conversion versus magic rounding
**Rejected.** `to_int_unchecked` yielded packed conversions and integer operations, with all 148 codegen regions clean, but an `exp2_checked` integer exponent-field path increased throughput cost from 1.399 to 1.584 cycles/element (+13.2%): it needed nine operations versus eight, including an explicit exponent bias that magic rounding embeds. `to_int_unchecked` on NaN is undefined behavior, and `.clamp()` does not remove NaN. A proposed one-FMA merge of the magic-round and field-extraction steps was mathematically wrong in 9,999,993 of 10,000,001 probed cases; the intervening de-bias changes the bit representation and cannot just be deleted. Range hints had no remaining float-to-int conversion to improve.

### Fit screens: measure the final f32 construction, not just the polynomial residual
**Finding.** Several lower isolated residuals became regressions once reduction, fma order, branch crossover and output rounding were included. For example, a three-and-a-half-word pi reduction improved bare residual maximum 1.19e-7→1.04e-7 but made sin max ulp 412→1,824,546 and cos 2780→386,929 on 100M inputs near zeros. Use same-seed, same-precision baseline and candidate runs, including boundary and exhaustive sweeps when a rare maximum matters: acos's quick-fuzz baseline max 4 was actually 5 exhaustively, and atan2's 10M-sample max fluctuated 3–4 even with unchanged code.

### Global LLVM flags: mixed regressions prevent crate-wide adoption
**Rejected.** `-force-vector-interleave=4/8` was bit-identical to default; 2 helped roughly 60% of ~70 functions but caused a +165.6% outlier. Forcing zmm AVX-512 helped 60/70 functions (median -19.4%) but slowed `exp_checked` 88.3% and another function 84.5% on the i5-1145G7: its single 512-bit FMA port exposed serial polynomial dependencies that the default double-unrolled ymm code hides. Stable Cargo could not scope these flags per function.

### Guarded kernels: follow reachability, not just direct calls
**Finding.** A caller's branch condition or construction can make a general kernel's exceptional-case handling unreachable, even when a branchless discarded arm still evaluates it. The sweep found such sites in hyperbolic functions, logs, inverse error functions and power-based color conversions; the latter were missed by a search for direct `log_2` calls because they reach it through `powf_pos`. Retain handling for NaN and infinity if the outer select can expose those results, and verify the licensed range, including signed zero, exhaustively when practical. In the later re-screen, unconstrained arguments to `rootn`, `xlog1py`, `compound`, `powf` and `signed_pow` gave no such license; the first three of those functions have since been removed.

### Inlined code: check whether LLVM already shares work
**Finding.** A proposed `sinh_cosh` pair computed `exp_pos_neg` once and matched separate calls bit-for-bit over all f32 patterns. After inlining, LLVM emitted the same 102-instruction sequence for the pair and for `sinh(x) + cosh(x)`; both measured 60.00 cycles latency and 2.219 cycles/element throughput in mca. No public pair was added: composition already gets the common-subexpression elimination at one call site, though not necessarily across non-inlinable boundaries.

### Integer polynomial evaluation: vector integer multiplication did not free the FMA ports
**Rejected.** On Tiger Lake, `vpmulld` costs two uOps and 1.00 reciprocal throughput versus FMA's one uOp and 0.50. A widening multiply-high/shift/add chain would use about 2.6× the port pressure of one FMA per term on the same execution ports. Resource screening ruled it out before fitting.

### inverse-function tests: round trips can conceal catastrophic outputs
**Finding.** The old `accuracy.rs` round trips reported residuals around 1.93e-7 for `probit` and `erfc_inv`, even when about 80% of uniformly sampled bit patterns produced infinity: applying the inverse map to infinity returned zero, close to tiny inputs in absolute error. Direct ulp rows exposed the problem and non-finite counters made it unmistakable. `edgecheck` also used round trips and sampled `probit` only from 0.001 upward. Score the function directly when an ill-conditioned inverse can erase the error under test.

### Isolated worst-case patches: check whether errors cluster
**Rejected.** `erfc` had about 4,900 bad points spread continuously across the domain, making a compare-select patch list impractical. The single bad mantissa in `cbrt_accurate` is patchable but remains an accepted won't-fix.

### jm check: inserted lines can cause false ownership warnings
**Finding.** `cmd_check` compares post-image diff line numbers with function ranges from master's pre-image `lib.rs`. A 52-line net insertion confined to `norm_cdf`/`norm_pdf` falsely reported later functions as “NOT YOURS”; inspect removed lines, whose coordinates have not shifted, before treating such a warning as an ownership violation.

### Leading-term peels: check correction size and signed zero
**Finding.** Rewriting `x*P(u)` as `fma(x, P(u)-1, x)` attenuates polynomial roundings only when the correction is small relative to the result; the `dawson` tail's correction was at most 3.3%, while its central branch's correction reached 30 times the answer at `x=4`. Such a peel can also turn `-0.0` into `+0.0`: when `P-1` is negative, both the zero-product sign and the final zero-addition matter. An ulp sweep cannot distinguish signed zeros; the `erfinv(-0)` edge check caught this, and applying the sign after evaluation on `|x|` preserved it without added throughput cost.

### Literal and NaN audits: check live hazards rather than adding generic helpers
**Closed.** The erroneous hand-typed pi/2 coefficient in `acos_poly` already had an `acos(±0)` regression pin; other mathematical constants used named constants, while fitted literals have no independently known intended bits. All seven `.max()`/`.min()` sites were either guarded by final NaN propagation, received NaN in both operands, or were already fixed; no further live NaN-discard bug justified a shared helper. By contrast, min/max is not automatically a safe replacement for a compare/select because IEEE maxNum/minNum can discard NaN.

### LLVM and build: test the emitted instructions
**Closed.** `codegen-units=1` plus LTO left the bench-profile assembly unchanged. Five LLVM flags (`-enable-unroll-and-jam`, `-extra-vectorizer-passes`, `-slp-vectorize-hor`, `-vectorizer-maximize-bandwidth`, `-enable-masked-interleaved-mem-accesses`) each left the 383662-instruction assembly unchanged; a control, `-force-vector-width=2`, changed it to 393186. When supplying `RUSTFLAGS`, include `-C target-cpu=native` because the environment replaces this repository's configured rustflags; check build errors, refresh the assembly even when Cargo considers the target current, and select the newest emitted `.s` rather than an old hash-named file.

### LLVM canonicalization: source rewrites without machine-code changes
**Rejected.** Combining zero/nonfinite tests with `wrapping_sub(1)` in cbrt and powf, and moving the asinh/acosh branch before their square root, produced byte-identical assembly. The latter was also bit-exact across all 2^32 inputs of each function and unchanged in llvm-mca latency and throughput. LLVM already performs these transformations; retaining more obscure source buys nothing.

### llvm-mca: branch layout and scheduler pressure can reverse the apparent result
**Finding.** Removing an `asinh` branch reduced instructions 156→145, uOps 163→151 and divider instructions 6→4, yet llvm-mca throughput worsened 4.136→4.760 cycles/element and its latency estimate rose 69.41→123.88 cycles. The old scalar latency region followed the cheaper, last-laid-out branch rather than the actual in-domain dependency chain; its 69.41-cycle figure was misleading. Throughput was limited by llvm-mca's modeled 60-entry scheduler (87–99% `SCHEDQ` stalls across comparable kernels), versus the machine's larger real window. Twenty alternating hardware rounds, normalized by unchanged `acosh`, differed only 1.4%, below the control's 4.7% drift. Check branch layout and scheduler statistics when instruction and cycle signals conflict.

### llvm-mca: check throughput against hardware and the resource counts
**Finding.** `sigmoid_grad` ranked seventh most expensive among 144 throughput regions at 4.386 cycles/element, versus 2.076 for the structurally similar `tanh_grad`. Yet it had fewer instructions (67 vs 74), uOps (73 vs 81), and Block RThroughput (22 vs 24), with the same two divides and 16 FMAs. Four hardware runs found `sigmoid_grad` faster every time, by 11–15% in the three clean runs. Its mca IPC was only 0.95 and 100 simulated iterations took 7017 cycles, essentially one 70-cycle latency chain per iteration: mca failed to overlap independent elements. For a large throughput region, inspect IPC and the instruction/uOp/resource ladder before trusting mca's cycles column.

### llvm-mca: isolate conditional arms before citing latency
**Shipped.** `tools/mca_arms.py` removes the arm an in-domain input does not take and measures the remaining chain. Of 151 latency regions, 39 contained branches and 24 published latencies exceeded *both* independently measured arms: `exp2m1` published 80 cycles versus arms of 37/48, `erf` 87 versus 44/65.98, and `acosh` 89.08 versus 46.08/77.02. Concatenated arms can create false dependencies; conversely `asinh` published 69.41 while its real mixed-arm chain was 88.02, because the cheap last-written register won. The throughput counterparts were branchless masked selections and did not have this problem. The apparent `srgb_to_linear` latency improvement of 110.02 → 104.30 was also an artifact: its in-domain arms measured 69.02 → 69.08.

### llvm-mca: confirm table rows before using them
**Finding.** A full remeasurement found 10 stale rows among 79 published rows: `erf` latency 83.98 → 87.00, `sinh_throughput` 1.943 → 1.689 cycles/element, `cosh_throughput` 1.616 → 1.498, and `erfcx` 2.896 → 2.827. The `remainder`, `remainder_ieee`, and `fmod` latency rows were also stale (34.11 → 40.13, 29.11 → 35.13, and 29.11 → 35.13); `sin`, `asinh`, and `erfc` drifted by less than 1%. Remeasure a row before using it to justify a change: stale numbers can conceal either a gain or a regression.

### llvm-mca regions: temporary wiring can corrupt unrelated measurements
**Finding.** Temporarily wiring `clog` into `mca_target.rs` yielded a useful same-wiring A/B, but changed unrelated `remainder_wide` region extraction from 4056 to 1495 instructions for latency and 134 to 67 for throughput. The multi-exit marker problem remains; compare the entire region list before retaining new instrumentation, not just the new region and two controls.

### llvm-mca throughput: a retained loop invalidates division by 16
**Shipped.** The throughput harness divided a region's simulated cycles by 16 assuming LLVM fully unrolled its 16-element loop. `tan_wide` retained a loop stepping by four, making its apparent 3.479 cycles/element actually 13.915; `clog_re` similarly measured 15.489 rather than 7.744 with a step of eight, though that row was unpublished. Both mca tools now read the induction step and report the retained loop; forward branches such as `powf`'s do not imply a loop. Check the assembly before trusting a surprisingly cheap composite.

### Lookup tables: price the lookup before fitting a polynomial
**Rejected.** On the i5-1145G7, `vgatherdps` has 4.0-cycle block throughput per ymm operation versus 0.5 for an FMA. One gather alone costs more than the five-FMA, 2.5-cycle exp2 polynomial it would replace, before adding the LUT's residual polynomial; it also cannot profitably replace enough operations in log_2. An in-register 16-entry `vpermi2ps` lookup costs 1.0 cycle, so this result does not reject a separately implemented SIMD tier using register tables.

### mca assembly custody and focused measurements
**Finding.** A concurrent `--emit=asm` build of another source state can overwrite the assembly between `examples/mca` building it and `llvm-mca` reading it, silently measuring the wrong tree. Serialize builds or use separate worktrees and target directories. Extracting only the relevant `LLVM-MCA-BEGIN`/`END` regions reproduced full-harness results in seconds instead of roughly 45 minutes and enabled bottleneck analysis. The latency chain's `mix()` sign-bit masking was also fixed: sign-dependent work had been absent from latency measurements.

### mca throughput: inspect the instructions and resource bound when the schedule looks implausible
**Finding.** Removing four vector multiplies and one broadcast from `coshm1` produced a reported +61.6% mca slowdown, despite an otherwise identical packed instruction mix. Its absolute 4.700 cycles/element estimate was also twice `sinh_checked`'s 2.340 for near-identical regions (102 versus 99 instructions); that row is not reliable for cross-function comparisons. For an `acosh` select merge, mca likewise reported +7.3% despite six fewer instructions, 16700→16100 uOps, unchanged 40.0 block resource throughput, and fewer expensive operations. Check opcode mix, uOps and `Block RThroughput`, not simulated cycles alone; conversely an `erfinv` fold decreased instructions 177→164 but increased uOps 19400→20300, a credible regression. A shared reduction in the removed `sincos_checked` also looked 19% faster in mca but was 1.4% slower on hardware after bench-shape corrections.

### Monotonicity: small local reversals are inherent to these approximations
**Finding.** An exhaustive adjacent-float walk over [-4,4] plus 50M wider-domain samples found 1–2-ulp drops in `sigmoid`, `tanh`, `erf`, and `atan` (and the then-present `softplus`), across 0.015–0.19% of steps. The checker's initial report of catastrophic sigmoid drops was a print-label bug; adjacent outputs at 0.18189578/0.18189579 were 0.54534906/0.54534894. Repair would require a different monotonicity-preserving construction, not a cheap local patch.

### NaN payload survey: quietness is a contract, payload preservation is not
**Shipped.** `nan_payload.rs` feeds three tagged NaNs to 108 unary functions; the later recorded classification was 76 preserving payload and sign, 21 varying sign, three canonicalizing (`erfinv`, `softplus`, `logsigmoid`) and eight yielding no NaN on inputs outside their promised domains. Comparisons, sign operations and hardware sqrt/select paths explain why the payload can change. The survey is a record rather than a gate; `special_matrix.rs` separately asserts quiet NaN behavior for functions promising it.

### PGO and static throughput: use the right measurement
**Finding.** A PGO/BOLT bench comparison was inconclusive under thermal noise: the same binary varied from 6.183 to 11.935 ns/op. `llvm-mca` cannot assess PGO's whole-program inlining and layout effects. Static throughput can also disagree with controlled wall time: one later `probit` change was reported +15.1% by mca but measured 1.599→1.527 ns/op in repeated quickbench runs.

### Polynomial evaluation: fold singly used fourth powers
**Shipped.** Rewriting `l0 + l1*x2 + l2*x4` as `l0 + (l1 + l2*x2)*x2` removes one multiplication and one rounding without adding dependency depth. Folding exp, paired-exp, tan and other callers reduced whole-file instructions 360429→358643; 29 throughput regions shrank and none grew. A separately copied tanh polynomial also dropped 79→77 instructions and exhaustive average/max from 0.1457/6 to 0.1452/5. Search by computation shape, not only macro name. Do not apply where the power is shared or the fold lengthens the critical path: erf throughput worsened 3.2%; a resource-bound erfinv variant reduced instructions 177→164 but increased uops 19400→20300 and regressed callers 21.8–25.3%.

### Polynomial fitting: measure the wrapped function, not just the fit
**Finding.** Degree reductions failed for `log_2` (9→8 exceeded its 2-ulp cap at 3–5), `exp2` (Q degree 5→4 estimated 40× worse relative error), and `sinf_poly` (9→7 estimated 178× worse). Local coefficient descent found no useful movement in `log_2`, `exp2`, `sinf_poly`, `expm1_near0`, `exp`, `asin_poly`, `ln`, `log10`, `erf_near0`, and two atan polys (typically <0.5%). A zero coefficient is a bad seed for bit-step search: `atan_poly` appeared stuck at max 18 ulp, whereas a least-squares seed reached 3 on the same grid; an arbitrary 1e-3 seed reached 565. Continuous-fit improvements can disappear under f32 evaluation: `log_2`'s weighted LP improved the ideal residual 8.4× but its actual Estrin chain scored avg 0.25128 vs 0.25137 over ~533M inputs, and `ln_normal`'s predicted 75% improvement yielded ~0.09% exhaustively.

### polynomial fitting: scale tiny LP objectives above solver tolerance
**Finding.** The `ln_normal` minimax LP minimized raw relative error around 1e-8 while HiGHS's default primal feasibility tolerance was 1e-7. It returned success and an objective of `-0.0` with arbitrary feasible coefficients; scaling weights by 2^24 made the degree-eight optimum 0.0685 ulp-equivalent. Sequentially fixing each coefficient to f32 and re-solving improved degree-seven error from 0.5832 to 0.5059 ulp-equivalent (continuous optimum 0.5033). Scale objectives to order one and check plausibility independently of solver status.

### Polynomial fitting: scale LP constraints and score the actual evaluation
**Finding.** On the `tanh`, `log10`, and sinh/cosh fits, unscaled small residuals fell below HiGHS's feasibility tolerance, sometimes producing a reported zero error. Scaling constraints to ulp units and verifying independently on a dense grid gave meaningful fits. Idealized minimax error alone was not predictive when coefficient cancellation or f32 evaluation dominated: use a simulated or measured real chain, and inspect `sum |c_k x^k| / |P(x)|` before increasing degree.

### Polynomial screening: separate approximation error from rounding-chain error
**Finding.** Compare exact-arithmetic, ulp-weighted error against a same-degree minimax fit before changing coefficients. `acos_poly` had 3.18 ulp-equivalent error against an LP optimum of 1.49 (1.54 after f32 quantization), while `log_2` was already at its 0.3092 optimum, pinned near zero by the rounded `LOG2_E`. A fitting grid must cover the binding endpoint and both signed branches: the old positive-only `acos_poly` grid stopped at 0.9998 and reported max 2.717 while the error reached 3.185 nearer 1. A stride-8 scan of 266M inputs still missed a true max of 5 and reported 4; use exhaustive bit-step-1 sweeps for max-ulp claims.

### Polynomial searches: score the actual rounded chain
**Finding.** An ideal-arithmetic, ulp-weighted LP screen compared shipped f32 coefficients with re-quantised minimax fits. Large fit ratios did not necessarily become real improvements: `asin_poly` had 9.8× apparent headroom, but over all 16,106,127 f32 values in `[0.27,1)` its candidate changed max/average from 6/0.87462 to 8/1.17783. Its worst input, 0.27004012, lies at a cancelling branch crossover. `atan_poly` had 14.2× relative headroom but only 0.228 ulp of absolute idealised error, so its real max of 4 was mostly outside the fit. Test a better polynomial through the complete f32 computation before investing in coefficient search.

### Polynomial searches: exhaustive domains and search objectives
**Finding.** An `asin_poly` coordinate descent trained on stride-32 samples improved its training max/average from 6/0.930 to 5/0.815, but regressed on the full 16.1M-value domain from 6/0.875 to 7/0.911. The original max-first basin hopping on `acos_poly` likewise produced a coarse max of 2 but real average 0.961 versus shipped 0.496. An average-first search with a maximum-error cap improved its coarse average 0.14447→0.12758 while holding max at 4; the exhaustive result was only 0.0650→0.0647 with max 5 unchanged, so the coefficient change was reverted and the search tool retained. Enumerate small discrete domains during search, not just final verification.

### Precision tables: stale log and cosine rows
**Shipped.** Exhaustive measurement corrected `log2p1` from 0.102/3 to 0.092/2 avg/max ulp, and `cos_checked (|x|<=1e6)` from 0.081/3 to 0.076/2. A thorough filtered run needs `accuracy thorough cos`, not the row label `cos_checked`, because the harness filters on its gate name. `cos_fast (|x|<2^22*pi)` averaged 0.288 rather than 0.291; its 2780 max reproduced.

### Probes: isolate signs, binades and factors outside the shared harness
**Finding.** An out-of-tree probe crate depending on the worktree and `sleef = "0.3.3"`, with `-C target-cpu=native`, allowed two compiled versions to score identical inputs and references without changing `accuracy.rs`. Per-sign and per-binade scans located error windows; scoring a composite against an exact-in-f64 factor and its actual f32 other factor separated their contributions. Without the native CPU flag the crate rejects a missing FMA, and a previous executable can otherwise be mistaken for a fresh measurement.

### Rational minimax: verify solver feasibility at sub-1e-7 scales
**Finding.** A joint six-coefficient nonlinear fit of `atan_poly` reduced idealized error from 1.949e-9 (0.033 ulp-equivalent) to 1.567e-9 (0.026), too small to explain the real 4-ulp maximum. Rounded-to-f32 coefficients changed the idealized score from the shipped 0.291 to 0.145 ulp-equivalent, but the shipped coefficients were optimized against the real chain, so that is not an observed improvement. HiGHS's default ~1e-7 primal tolerance incorrectly called demonstrably feasible ~1e-9 constraints infeasible; rescale around the incumbent or use a seeded nonlinear optimizer, and always check that the incumbent is reported feasible.

### Saturation pins: test beyond the immediate clamp boundary
**Shipped.** `saturation_pins.rs` compares roughly 3,400 values in ±64-ulp windows around ten functions' bounds against f64, plus 20 far-outside saturation pins. `exp2_checked`, `exp10_checked`, `exp_checked`, `expm1_checked`, `exp2m1`, `exp10m1`, `tanh`, `sinh_checked` and `cosh_checked` stayed within 1 ulp in-window, and all 20 limits were exact. `sigmoid`'s known early flush produced 2,098,176 ulp near −88.72235; the gate exempts 129 points by *denormal true result*, rather than by input range, to keep normal-output regressions visible. Raw ulp is a poor severity measure for subnormal outputs.

### Special-value tests: sweep the public unary and binary surfaces
**Shipped.** `special_matrix.rs` checks ±0, ±inf and NaN for 108 unary functions; it found a sixth signed-zero bug after five had been found ad hoc. Its NaN quietness/propagation gate names eight off-domain `_unchecked`/`_approx` exemptions explicitly, so newly added functions do not inherit an exemption by naming convention. `special_matrix2.rs` covers the 9×9 cross product for binary functions: all nine with an f64 standard-library counterpart matched bitwise in 729 comparisons, with no new bugs. The remaining binary conventions, including `xlogy(±0,y)=+0` even for NaN `y`, were checked against their contracts; off-domain unchecked tiers were not asserted.

### Strided sweeps: useful for deltas, not absolute maxima
**Finding.** The same every-128th-f32 sweep measured `erfinv`/`erfc_inv`/`probit` maxima of 3.67/4.50/4.98 against scipy, whereas exhaustive maxima were 4/5/5. Compare variants on identical strided points, but publish exhaustive maxima where available.

### Toolchain changes: remeasure scheduler-dependent rejections
**Finding.** On rustc 1.98.0-nightly, the previously rejected ln/log10 trailing-FMA association and `pre_offset` dead-add removal became wins; the latter measured `tan_checked` -6.3% and `sin_checked` -1.5% at the time. The `reduce_pi` depth-two rebalance still cost 3 cycles. Remeasure scheduling-only rejections after toolchain changes, including both associations of an FMA fold; structural dependency-depth losses are less likely to reverse.

### Trig invariant: measure bounds separately from ulp accuracy
**Finding.** An exhaustive scan of `|result| <= 1` needs no reference function and caught millions of one-ulp overshoots that the accuracy rows could not distinguish from ordinary approximation error. Removing the clamp from `sin_wide`/`cos_wide` produced 2,726,588 violations, worst 1.0000001, despite exact argument reduction and passing edge pins. The polynomial itself overshoots near its maximum; a proof about reduction does not prove the output bound. The clamp was retained and can only improve accuracy against a true sine/cosine bounded by one.

### Vectorization check: recognize double-precision arithmetic
**Shipped.** `codegen_check` formerly required a packed `ps` arithmetic mnemonic, falsely rejecting an eight-lane f64-internal `remainder_wide` loop containing `vdivpd`, `vfnmadd213pd` and `vrndscalepd`. It now accepts `pd` arithmetic too. Such checks must recognize the actual internal arithmetic type, not infer it from the f32 input and output.

### Wide-trig benchmarks: measure real cycles and avoid constant exponents
**Finding.** On the i5-1145G7, `mca` priced gathers too cheaply: the gather-free x8 experiment moved Block RThroughput only 38 → 36 while measured `sin_wide` throughput moved 9.69 → 1.56 ns/element. On the Ryzen AI Max+ 395, `reduce_pi_wide` instead filled the FP scheduler (`fp_sch_rsrc_stall` around 65% of cycles); neither Block RThroughput nor port counts predicted the result. The old latency regions also constant-folded the exponent-driven window. `examples/trig_bench.rs` mixes a black-boxed pseudo-random exponent into latency chains; interleaved A/B runs with `tools/ab.sh` and exhaustive `trig_sweep` were used for the final pass.

### Worst-case corpus: a fast change detector, not an accuracy oracle
**Shipped.** `worst_corpus.rs` checks 9,720 blessed output bits (108 unary functions × 90 inputs) in 0.058 s without building f64 reference implementations; full sweeps still establish correctness. A 0.5% shared Padé coefficient perturbation moved 61 outputs, but a one-ulp coefficient perturbation and an `exp10m1` seam move from 0.2 to 0.21 both escaped the corpus. A passing corpus does not prove outputs unchanged: test a new gate with perturbations at several magnitudes, including inputs on both sides of seams rather than only the seam itself.

## f64 and double-float

### cbrt_accurate: no useful f64 port of a small residual block
**Closed.** Its `Df32` residual accounts for 7 of 77 instructions in the measured unchecked-throughput region, while accuracy already measured avg 0.000/max 1 ulp. f64 could form the residual more precisely, but there was no demonstrated accuracy gain to buy with a port.

### clog: compensate both squares near the unit circle
**Shipped.** Computing the real part as `0.5*log1p(re²+im²−1)` with error terms for *both* squares and their sum avoids rounding `cabs` at magnitude one. Standing fuzz max fell 4717→3 ulp, latency 170.63→158.93 cycles and throughput changed 8.166→8.227 cycles/element (+0.75%). Compensating only one square left max 1466. A targeted 400k-point unit-circle probe instead found old 1.68e7 versus new 2048 ulp: neither has a uniform ulp bound as the true log approaches zero, but compensation moves the cancellation threshold from roughly 2^-24 to 2^-48.

### clog: replace the f32 error-free sum near the unit circle
**Shipped.** The f32 decomposition of `re²+im²−1` was exact algebraically but summed its correction words in f32, introducing absolute error around 8e-15. On the `|z|=1` manifold, max error was 4096 ulp for `|v|<1e-9`; using `a=max(|re|,|im|)`, `b=min(...)`, and f64 `v=mul_add(a,a,-1)+b*b` brought that band to zero and the `1e-9..1e-8` band from 14 to 1. Subtracting one from the larger square is essential: it makes that FMA exact in the near-unit-circle branch; subtracting from the smaller square does not. The measured region also improved from 263 to 258 instructions, Block RThroughput 67.5 → 61, throughput 8.057 → 7.728, and latency 295.31 → 281.88.

### compound: preserving log1p's low word before exponentiation
**Finding.** A claimed underflow explanation for ~200 max ulp was false: 5,435 of 5,440 samples above 50 ulp had normal outputs, with `|n*log1p(x)|` between 30 and 86. The former double-float `compound_accurate` reduced max 236→4–5 and avg 0.195→0.019, but throughput rose 3.724→9.239 cycles/element and latency 96.75→153.56. Its correction required both a division residual and a quadratic Taylor term: omitting the former raised max to 44 because multiplication by n amplifies its lost low bits too. Replacing the original function with `powf(1+x,n)` was much worse (max 1.07e9/avg 5.0e6 versus 240/0.29), since forming `1+x` loses the increment.

### compound_accurate: replace double-float bookkeeping with f64 reduction
**Shipped.** The historical function's f64 atanh-form `log2(1+x)` kept `x` as the numerator, avoiding the cancellation that forced `log2p1_df` to reconstruct lost bits of `1+x`. Instructions fell 255 → 159, uOps 302 → 211, Block RThroughput 84.5 → 76, throughput 9.239 → 6.173 cycles/element (-33.2%), latency 153.56 → 120.77; quick accuracy went from avg 0.0194/max 5–6 to 0.0007/max 1 ulp. This completed removal of the unused double-float log/exp helpers (115 lines); the public `Df32` module remained.

### Double-float log scheduling: overlap correction with refinement
**Finding.** The first higher-precision powf log implementation cost 14% latency and 24% throughput. Refitting a reciprocal seed in the un-offset variable, evaluating a slowly varying correction polynomial from an unrefined argument while retaining the refined prefactor, and ordering exact two-sums by operand readiness eliminated the latency cost without changing coefficients; throughput cost fell to 7–9%. An Estrin split after that correction left the critical path unchanged, while Horner evaluation of the seed worsened both metrics.

### f64 ports: count bookkeeping, not just operations
**Finding.** The historical `remainder_wide`, `powf`, and `compound_accurate` ports changed throughput by 7.688 → 2.266, 8.459 → 6.996, and 9.239 → 6.173 cycles/element respectively, with max error improving from a 2^48 to a 2^53 contract, 3 → 1 ulp, and 5 → 1 ulp. On the measured machine f64 vector FMA costs twice as much per eight lanes, so a port pays off when it removes enough f32 error-free-transform work. `compound_accurate` cut 38% of its instructions; leaner `powf` cut only 10–18% and its Block RThroughput worsened.

### f64 trigonometric reduction: one f64 pi word was insufficient
**Rejected.** A newly built f64-reduction sin/cos path reached max 637 ulp even in its smallest bucket; two-product compensation reduced that to 70 but could not repair error already present in the f64 pi constant. Extra precision in arithmetic does not restore bits absent from the reduction constant.

### logaddexp_accurate: widen the cancellation chain
**Shipped.** The removed function ran `m + log1p(exp(-d))` wholly in f64, including `d = |a-b|`, and narrowed only once; its correction used `2*atanh(s)` with `s = 1/(1+2*exp(d))`. On 20,000 high-precision pairs near the zero curve, avg/max ulp was 0.2675/72.2 versus 16,917,850/34,372,019,720 for the f32 tiers. The remaining worst sample had absolute error 6.3e-17 at a true result of 1.065e-11: cancellation magnifies even f64 error. Its 7.616 throughput and 121.517 latency cycles versus 2.449/74.110 for `logaddexp` bought the accuracy; f32 algebraic rearrangements cannot remove an absolute error of about 2^-24 in the correction.

### logaddexp_accurate: avoid an unnecessary denormal polynomial tail
**Finding.** Its atanh variable reaches about 1.3e-56 over the useful domain; powers in the Estrin tail could become subnormal near `d = 44`. Flooring `u = s*s` at 1e-30 kept the tail normal while changing it by only about 2^-101 relative to the leading `2s`. For wide-domain polynomials, check intermediate powers as well as the final value.

### powf: the error was in log precision, not the final multiplication
**Finding.** At a point with −312 ulp error, a compensated two-product multiplication recovered only about 13% (to −270); the single-f32 log_2 error was magnified by the exponent. A double-float pair containing only the exact rounding residue of the final polynomial multiply and add is not a higher-precision logarithm: an f32-evaluated polynomial still limits relative precision near x=1. In the historical checked tier, changing the log's evaluation form cut max from at least 203 to 3 ulp and average 0.046→0.019, at unchanged latency and 7–9% more throughput cost; the plain fast powf remained around 292 ulp.

### powf: replace an inaccurate fast route with the double-float route
**Shipped.** `powf` absorbed `powf_checked` and uses `exp2_checked_df(log2_df(abs(x))*y)`: max ulp fell from at least 292 to 3 and avg from 0.181 to 0.019. Against the former checked tier, latency fell 152.89→125.17 cycles (−18.1%) and throughput 9.918→8.459 (−14.7%); against the inaccurate formula, cost rose 19.3%/49.7%. Exponent amplification made the fast formula hundreds of ulp wrong at ordinary inputs, not a useful accuracy tier. Route zero/inf/NaN through log-high-word selects and existing exp clamps rather than a separate fallback tree; a compact parity-based sign multiplier replaced roughly 28 vector ops with 20.

### powf: guard the exponent correction before the dependency-chain tail
**Shipped.** Moving `exp2_checked_df`'s correction onto its finite, normal pre-`t2` product and guarding the correction rather than the final result reduced simulated region throughput 10.3% on its own. It also preserved denormal rounding and fixed an overflow case: `powf(f32::MAX, 1.0000001)` had incorrectly returned `f32::MAX` after a negative correction to an exponent clamped from 128.0000153 to 128. In `log2_df`, rewriting `s-th*(s+2)` as `(s-2*th)-th*s` made its exact-denominator residual two fmas instead of a four-operation split, with bit-identical output.

### powf: cheaper double-float log corrections did not retain accuracy
**Rejected.** A multiplicative `log2_df` correction (`hi*(1+w)`) saved two ops but raised `powf` max ulp 3→10: it omitted the low word's scaling and exposed the full relative error of `u`, where the absolute correction exposed only the polynomial's roughly 0.0176 log-derivative. A degree-5 seed without the Newton step saved 3.6% latency and 2.2% throughput but raised max ulp 3→4. Estrin regrouping saved 2.5% latency at 2.2% more throughput and three more uOps. Compensated arithmetic is only useful when the kernel and its conditioning preserve the extra bits.

### powf: move the logarithm–exponent chain into f64
**Shipped.** Replacing `log2_df`, double-float multiplication and `exp2_checked_df` with an f64 chain brought the search worst from 3 to 1 ulp, worst per-band average 0.3027→0.0106, and quick average 0.0193→0.0007. `powf` throughput fell 8.459→6.996 cycles/element (−17.3%) with essentially unchanged latency, 125.17→125.81; `powf_unchecked` fell 6.509→6.171 (−5.2%). An f64 eight-lane FMA costs twice an f32 eight-lane FMA on the i5-1145G7, so f64 wins only where it removes enough double-float bookkeeping. A degree-5 reciprocal seed and one Newton step avoided the roughly 5% latency penalty of the degree-3/two-step candidate; shorter polynomials also reduced f64 constant broadcasts.

### remainder_wide: two exceptions to a claimed bitwise equivalence
**Finding.** Compared with historical `remainder_checked`, `remainder_wide` differed by up to about 4 ulp on 2798/29M extreme-magnitude denormal-rescale cases, and flipped sign at 3/292M exact half-integer ties because an intermediate `.round()` corrected the tie twice. Neither exception was changed; they were documented and excluded from the equivalence test's domain. The tiers were not interchangeable: 8.186 versus 1.544 cycles/element and 142 versus 47 instructions for wide versus checked, reflecting wide's large-ratio contract.

### remainder_wide: replace double-float reduction with f64
**Shipped.** The f64 rewrite cut 134→67 instructions, 141→85 uOps, throughput 7.688→2.266 cycles/element (−70.5%) and latency 171.17→58.13 (−66.0%). Its quotient is exact as an integer through 2^53, avoiding the coarse-quotient correction, second division and double-float residual; f64's exponent range also avoids the overflow rescale. In simulation the old chain first had 1532 inexact results among 6000 samples in `[2^50,2^52]`, versus zero for f64 through `[2^52,2^53]`; widening the real harness to `|x/y|<4e15` gave old 2.55e6 avg / 3.42e9 max ulp against new 0 / 0. The final single correction remains necessary because rounded f64 division can cross a half-integer near quotient 2^52.

### sRGB transfer functions: the log2/exp2 round trip sets the error
**Closed.** A 4.9M-point scan of `srgb_to_linear`'s power arm found 1.496 avg / 13 max ulp overall, 1.670 / 6 from rounding `(c+0.055)/1.055`, and 1.730 / 13 from `log_2_normal`, multiplication by 2.4, and `exp2_checked`. At the 13-ulp worst point the argument-rounding contribution is zero. Carrying the argument in double-float could improve average but not max; `log2_df` is limited by its f32 kernel's evaluation, and changing the shared log kernel for this bounded `[0,1]` application is not warranted. `linear_to_srgb` reaches 7 max ulp with exponent `1/2.4`, which amplifies log error less.

## exp family

### exp: weaving the scale factor into the polynomial
**Rejected.** It improved exhaustive avg/max ulp 0.0745/3→0.0522/2, but throughput worsened 1.327→1.393 cycles/element: it trades a parallel add for an additional FMA on busy ports (4 mul + 4 FMA + 1 add versus 4 mul + 5 FMA). Weaving either scale factor has the same cost. Pre-multiplying the two factors is incorrect near overflow: at x=88.376274 their product alone is `inf`, while the interleaved expression returns finite `2.4061984e38`. Scratch functions coexisting in one mca target misleadingly looked faster (1.274/1.303); replacing the real function in isolation reproduced 1.393 three times.

### exp and exp_checked: retain the round-based reduction
**Closed.** A floor-based reduction would need a same-shape degree-five exponential polynomial with about 4.9 ulp-equivalent idealized error, versus roughly 2.0 for the existing polynomial, or a higher degree. Narrowing the existing reduction's domain to avoid a split field achieved the intended single-field path in `exp_narrow` without refitting the polynomial.

### exp kernels: centering and rational refits
**Rejected.** Centering exp2's remainder produced larger, not smaller, coefficients; the incumbent monotonic fit already uses its edge minimum. Rational refits for sin and erf kernels either fitted far worse or timed out; do not assume a different basis or rational form helps without checking its critical-path division and realised fit.

### exp10: adjusting the reduction and its low word
**Rejected.** Replacing compare/select/add/sub with `fr.floor()` unexpectedly increased max 1→2 and `exp10_checked` throughput 2.736→4.389 cycles/element (+60.4%). Varying `LOG10_2_LO` by ±1–8 ulp left `exp10` and the then-current checked variant at about 0.0306–0.0307 avg and max 1, supporting the decision not to add a third Cody–Waite word (the earlier exp10 score was 0.0343 avg/max 2 on >2.2B samples). Source-level equivalence and fewer-looking operations do not guarantee equivalent f32 output or scheduling.

### exp10_checked: centered Q and round reduction
**Shipped.** A re-derived centered degree-5 Q removes the floor-adjust compare/select and two adds/subtracts without moving the shipped saturation clamps. Exhaustive avg ulp fell 0.0307→0.0082, max remained 1, throughput 2.361→1.897 cycles/element (-19.7%), and latency 63.56→51.06 (-19.7%); `exp10_checked(inf)` remains `inf`. An earlier centered implementation had returned finite `3.237e38` for infinity and max 2, so edge and saturation pins are essential; that particular failure did not recur in the re-derived version.

### exp10m1: the 0.2 crossover is already on a plateau
**Rejected.** Exhaustive 4.29-billion-sample sweeps at 0.18, 0.205, 0.21 and 0.215 matched the shipped 0.2 crossover at average/max 0.1279/4 and worst x≈−0.057932023; 0.15 worsened to 0.1283/5 and 0.2171 to 0.1280 average. Since |x·ln(10)| must remain below 0.5, widening stops near 0.217 (+8.5%), unlike exp2m1's 30% widening. The max lies well inside the Padé branch and cannot be moved by this seam.

### exp10m1: fit the reduced decimal residue directly
**Shipped.** The near-zero Padé and separate `exp10_checked`-style arm became one degree-4 approximation of `10^d - 1`, combined as `2^k D + (2^k - 1)`. Fitting the centred reduction residue `d` avoids the `d*LOG2_10` rounding and removes the seam, division and floor-adjust; a low-magnitude arm returns `d*LN_10` to preserve negative zero. Exhaustive avg/max improved from 0.1279/4 to 0.0946/3 ulp; throughput fell from 2.629 to 1.342 cycles/element, with 106→64 instructions and 116→69 uOps. The lower clamp tightened from -45.154503 to -37.0 because results below -7.53 already round to -1; the upper 38.53184 remains the overflow boundary.

### exp10m1: a higher degree and a two-word logarithm did not repay their cost
**Rejected.** With the shipped one-word `LN_10`, degrees 4 and 5 both had a 0.2316-ulp idealized margin; an exact-`ln10` oracle reduced degree 5 to 0.0032 but degree 4 only to 0.1830. The constant, not degree-5 fit error, would require a two-word peel and another term, while reduction and closing-fma roundings remain around 0.5 ulp each. The new branch-shaped latency is genuinely 55 cycles on either arm, versus old arms of 37/80 cycles (whose published 112 was an artifact); the former near-zero arm thus pays 18 cycles despite the throughput and accuracy wins.

### exp2: centered reduction needs whole-function and edge checks
Rejected for `exp2`, `exp10`, `exp2_checked`, and `exp2m1`. Replacing floor with round centers the polynomial interval at [-0.5,0.5], but a degree-4 Q reached max 13 ulp; degree 5 made unchecked `exp2`/`exp10` produce NaN inside their documented domains when round gave k=128 (for example x=127.6). `exp2_checked` throughput grew 1.399→3.007 cycles/element and max 1→2; `exp2m1` max grew 3→6/7 and throughput +3.3%. The double-float variant improved `powf_checked` avg 0.0229→0.0118 and max 123→118, but cost ~3.4% throughput; it too was reverted.

### exp2: refit its Q polynomial with a max-capped LP
**Rejected.** The weighted continuous error improved 0.043→0.013, but actual `exp2`, `exp2_checked`, `exp10`, and `exp10_checked` all regressed from max 1→2 ulp. The original fit had no worst-case margin for a continuous objective to spend.

### exp2: fit 2^f directly with a pinned constant
**Rejected.** Replacing `1+f·Q(f)` by `P(f)` with coefficient zero pinned to 1 would save a multiply but loses a degree of freedom across all of [0,1). The converged fit reached max 463/avg 232 ulp versus shipped max 2/avg 0.203.

### exp2: select-tree table and polynomial
**Rejected.** A proposed two-level lookup with a degree-3 polynomial was about 24 times less accurate than the existing fit in the isolated screen. Even a three-level degree-8 version worsened real avg ulp 0.20→0.43 and max 2→3.

### exp2 and expm1: folding scaling constants into dedicated fits
Rejected after screening. Folding `LN_2` into `exp2m1`'s rational fit offered only ~5.6× idealized margin, comparable to constant-folded pi/degree trig fits that regressed in real f32 evaluation (2.6× and 4.9× margins). The shared `exp2_q_poly` and `pade_expm1_ratio` would affect multiple callers; no implementation was made.

### exp2 exponent-field split: simplify shared scaling
**Shipped.** Replacing per-word `(k+383)<<8 & EXPONENT_MASK` with `exp2_field_split` removed two ops and three constants across roughly 15 functions. Throughput improved 13.2% for `sinh_checked`, 12.7% for `exp10_checked`, 8.8% for `exp_checked`, 4.7% for `exp_m1_over_x` and 1.0% for `expm1`, with no row slower. Although the split chooses a different intermediate exponent for odd `k`, every reachable `k` against 4,096 mantissas gave bit-identical results after final scaling.

### exp2_checked: integer exponent-field construction
**Rejected.** An ordinary `k as i32` cast de-vectorized the region despite an earlier runtime clamp, because LLVM could not prove its bound. `to_int_unchecked` restored packed conversion (two `vcvttps2dq`, zero scalar conversions; all 148 codegen regions passed), but needed a NaN-quashed copy through `max(-151).min(128)` to avoid undefined conversion while preserving NaN output. Throughput still worsened 1.399→1.584 (+13.2%) for `exp2_checked`, 2.437→2.651 for `erfc`, and 5.651→5.714 for `powf`. Magic rounding embeds the exponent bias in eight operations; the integer route takes nine and its conversion competes for the same ports.

### exp2m1: widen the Padé crossover
**Shipped.** Moving the branchless-select crossover from |x|=0.5 to 0.65 reduced exhaustive average ulp 0.0769→0.0766, held max at 4, and added no llvm-mca cost. This succeeded because the changed branch covers a substantial fraction of inputs.

### exp2m1: a narrower caller-specific Padé fit had no useful target
**Rejected.** After the crossover moved to 0.65, its Padé argument reaches |y|≈0.4505 rather than the formerly assumed 0.347. The incumbent's idealised error peaks at an interior y≈0.334 and remains 0.617 ulp-equivalent under either narrower domain; the actual max 4 at x≈−0.3991 has only 0.158 ulp idealised fit error. The error is primarily evaluation rounding, not the shared fit's domain width.

### exp2m1: fit in the original argument with a rounded reduction
**Shipped.** The old near-zero Padé first rounded `x*ln(2)`, and a floor-based reduction would send small negative inputs to a remainder near one, losing the subtraction (`x=-1e-8` could return zero). Rounding `k` instead keeps `f=x-k` in `[-0.5,0.5]`; a direct fit of `2^f-1` with the leading `f*ln(2)` peeled into an fma removes the Padé and division. Exhaustive average/max improved from 0.0769/4 to 0.0400/2 ulp, and throughput from 1.843 to 1.278 cycles/element (30.7% faster); the main old scalar arm was 48 cycles versus 52 after. The denormal select also returns the signed leading term at negative zero, which the peeled fma alone would turn positive.

### exp_m1_over_x: factor out the removable singularity
**Shipped.** With `e^r-1=r+r²P(r)`, its near-zero quotient becomes `fma(r,P(r),1)` and is exactly one at zero without division; only the `k!=0` arm divides. Exhaustive average/max fell from 0.0709/5 to 0.0170/2 ulp and throughput from 1.639 to 1.279 cycles/element. The narrow variant's throughput fell from 1.347 to 1.276. The old near-zero Padé's division had been solving a singularity that an explicit polynomial factor removes algebraically.

### exp_pos_neg: adding terms to the even/odd polynomial
Rejected after screening. Sinh's worst-point budget had 3.055 ulp from rounding and 1.670 from truncation; its ~1.3e-7 relative polynomial error is magnified about 8× by reconstruction. A degree-7 even/odd fit across `|r|≤ln(2)/2` gave only ~5.9× idealized margin and its new odd coefficient converged to zero. Splitting coefficients for sinh and cosh offered no narrower residual domain: both use [-ln(2)/2,ln(2)/2] and had the same 1.325e-7 maximum fit error.

### exp_r_poly: a degree bump helps accuracy but costs every caller
**Rejected.** Degree 5→6 had ~33× idealized fitting margin and improved average ulp 9–37% for `exp`, `exp_checked`, `expm1`, `exp_m1_over_x`, and `sigmoid`: respectively 0.0744→0.0473, 0.0389→0.0246, 0.1304→0.1189, 0.0729→0.0626, and 0.0925→0.0835. Max improved 3→2, 3→2, 6→4, and 6→5 for the first four; sigmoid stayed at 4. But all five paid 4–5 cycles latency (+6–9.5%) and +5–23% throughput, with sigmoid +23.2%. The shared macro spreads both gains and cost.

### exp_reduce!: halving the reduction interval costs more than it saves
**Rejected.** Reducing by `ln2/2` rather than `ln2` would shrink the degree-5 minimax error from 7.5e-8 (1.26 ulp) to 1.2e-9 (0.02 ulp). Reconstructing odd half-powers needs a parity/floor calculation, blend, and a two-part `sqrt(2)` scale: roughly six added instructions, versus five fmas and 2.5 Block RThroughput cycles for the entire existing polynomial. Dropping to degree 4 at half range gives 1.37 ulp, essentially the original 1.26, so it cannot fund that reconstruction. The shared reduction serves `exp_checked` and the Gaussian callers, making this cost widespread.

### exp_reduce!: peel the leading one and use degree six for its callers
**Shipped.** A peeled degree-6 approximation of `e^r - 1` reconstructs with `fma(s, t1, t1)`, where `t1` is a power of two, instead of rounding `1+r` inside the polynomial. On `|r| <= ln(2)/2`, isolated max/avg fell from 2.392/0.7165 to 0.830/0.2654 ulp; exhaustive `exp_checked` went 0.0370/3 → 0.0042/1, `erfcx` 0.1439/6 → 0.1333/4, and `erfc` 0.1289/7 → 0.1217/6. The peel saves an add and turns a multiply into an fma, paying for the extra coefficient: one fewer arithmetic instruction per call and unchanged Block RThroughput across the five measured regions, though latency gains an fma level (50→54 cycles for `exp_checked`). A distributed degree-6 grouping kept latency at 50.001 cycles but raised Block RThroughput in every region (for example 41→42); the folded form was chosen for throughput. A rejection measured on a broader shared kernel need not apply to a narrower caller set.

### expm1: raising the near-zero Padé degree
**Rejected.** Degree 3→5 reduced near-zero avg ulp 0.138→0.135 but left the function's max 6 in the other `exp(x)-1` branch and cost +7% throughput. At its worst point, rounding and truncation/coefficient errors were comparable (2.946 and 2.579 ulp), not a single isolated bottleneck.

### expm1: half-argument doubling buys a max ulp at too much cost
**Rejected.** In the cancellation band, `expm1(2u)=a*(a+2)` for `a=expm1(u)` has amplification 1.124 at `x=0.5`, versus 2.541 for direct `exp(x)-1`. Doubling throughout `|x|<1` changed max 6→5 but avg 0.1304→0.1506. A third branch limited to `[0.5,1)` improved both avg 0.1304→0.1268 and max 6→5, but throughput worsened 1.695→2.366 cycles/element (+39.6%) because it evaluates another Padé including a divide. The apparent latency improvement was a harness artifact, not grounds to accept the cost.

### expm1: widening the shared Padé fit cannot cheaply move the seam
**Rejected.** An idealized same-degree minimax refit yielded max f32-ulp-equivalent fit error 0.65 on `|v|<0.5`, 3.15 on `<0.65`, 10.98 on `<0.8` and 42.04 on `<1`, before f32 evaluation rounding. Thus moving the seam from 0.5 to 1 with today's [2/3] rational cannot beat the whole function's 6-ulp max; moving the seam alone also evaluates the current fit outside its domain. The shipped fit's idealized 1.37 versus attainable 0.65 on its existing domain is not a demonstrated real-chain win, and the worst errors are on the other side of the seam.

### expm1: approximate the difference rather than subtracting one
**Shipped.** The reduction now approximates `E=e^r-1` directly and reconstructs `2^k E+(2^k-1)` with an fma, avoiding the cancellation in forming `e^r-1`. This removed the Padé, its division and the near-zero seam: exhaustive average/max fell from 0.1291/5 to 0.0166/2 ulp; throughput fell from 1.604 to 1.087 cycles/element (32.2%). The scalar near-zero arm's latency rose from 32 to 47 cycles, but the vectorized old loop paid for the division even on lanes selecting the other arm. A tiny-input select is still needed below `2^-125`: halving the intermediate before doubling lost all 2^24 denormal inputs and also lost negative zero without that select.

### expm1_checked: offset the field and double the power
**Shipped.** Clamping `k` to 127 would incorrectly saturate finite outputs near the top of f32's range and suppress overflow. Forming a single exponent field at `k-1` and then adding the resulting power to itself permits `k=128`; the clamp is [-86.0, 128/log2(e)], and exhaustive checking found `k-1` in [-125,127]. The result is bit-identical to `expm1` throughout its valid domain (max 6 ulp); throughput improved 1.695→1.595 cycles/element (-5.9%), though latency grew 71.00→75.06 cycles because the clamp precedes the dependency chain.

### Integer exponent fields: use a maskless shifted magic constant
**Shipped.** `exp2int_field!` uses `f32::from_bits((k + 12583039.0).to_bits() << 23)` instead of shifting by 8 and masking the sign bit. Folding the exponent bias into an exact integer addition in the 2^23 binade saves one operation (two versus three); this is an operation-count win, not the hypothesized integer-port migration. The forms are bit-identical for every integer `k` in [-127,128], including underflow to zero and overflow to infinity. Nine sites were consolidated into the macro, shrinking 14 public functions: `exp2` 41→38 instructions, `sigmoid` 55→51, `exp10` 63→60. Unclamped out-of-domain outputs can now be negative nonsense (`exp2(-200)` changed from +7.6e-6 to -7.2e16); all 135 changed worst-corpus entries were outside the affected functions' documented domains, and clamped callers did not change.

### Shared exponential reduction: two small erfc-side changes did not pay
**Rejected.** Carrying the low Cody–Waite term multiplicatively would eliminate a rounding worth at most 0.25 ulp in `exp`, but adds an operation and a dependency level to five callers. Moving `erfc`'s `pe` adjustment into the exponent is an operation-count wash and introduces another reduction rounding up to 1.5e-8; the existing linearization error is below 1e-11. Changing the reduction does not require changing the shared polynomial, but neither change improved the callers enough to justify its cost.

### Split exponent fields: smaller code was reverted over unchecked outputs
**Rejected.** A maskless `exp2_field_split` took six operations instead of eight and shrank 23 throughput regions without growing any: `exp` 60→55 instructions and `powf_checked` 251→245. It was bit-identical for `k` in [-200,200), with identical in-domain fuzz across 16 groups, but the unchecked `cosh(2048)` changed from infinity to a plausible finite 1.0666397e35; all 55 changed corpus entries were unchecked and out of domain. It was reverted. Applying the variant only to callers that clamp `k` remains a possible zero-behavior-change optimization; a verified bit-trick range must include every caller's *reachable* values, not merely normal test values.

## log family

### Denormal normalization: do not replace the rescale with leading-zero shifts
**Rejected.** AVX-512CD does lower vector `leading_zeros` to packed `vplzcntd`, but `denormal_rescale!` already takes one compare, one multiplication and two selects. Exact shift normalization needs leading-zero count, shift-count arithmetic, shift, masking, exponent reconstruction and gating: seven or eight operations instead of four. The exponent-offset select cannot be moved after the rounded log result without changing accuracy; the hypothesized port-shifting benefit had also lost in a related exponent-field measurement.

### ln and callers: exhaustive accuracy after the kernel change
**Finding.** Exhaustive full-pattern checks measured `ln` avg/max 0.117/3 → 0.0036/1, `ln_unchecked` 0.235/3 → 0.0073/1, `log1p` 0.097/4 → 0.0249/2, `asinh` 0.149/3 → 0.0340/2, and `acosh` 0.060/4 → 0.0031/3; `atanh` remained 0.0037/2 and `log1pmx` 0.0632/3. Quick fuzz had underestimated `acosh`'s new maximum as 2 rather than 3. A better polynomial fit was not the main remedy: where the split constant is combined and where the leading term is rounded mattered more.

### ln and log10: a trailing FMA fusion cost a cycle
**Rejected.** A fresh baseline measurement showed the proposed fusion reproducibly added one cycle; an earlier apparent win had compared against a stale baseline. Re-measure before relying on a recorded cycle count.

### ln_normal: fold the low exponent word before the polynomial
**Shipped.** `fma(p,s,fma(k,LN2_LO,k_hi))` reduced `ln_unchecked` throughput 1.113→1.018 (-8.5%) and latency 38.22→34.06 cycles, with downstream gains including `acosh` -3.9%, `log1pmx` -3.9%, and `asinh` -3.3%, without accuracy loss. The opposite association, `fma(k,LN2_LO,fma(p,s,k_hi))`, still cost a cycle; it is scheduling, not merely operation count. Applying the successful association to `log10_normal` worsened exhaustive `log10` avg/max 0.1265/3→0.1280/4 and `log10p1` 0.2319/3→0.2347/4, so only ln's half shipped.

### ln_normal: better coefficients did not beat evaluation rounding
**Rejected.** Across all 8,388,608 mantissas, the shipped degree-8 polynomial scored idealised 1.101 and real max/average 3/0.4531; an old degree-9 polynomial improved idealised error about 22× to 0.05 but measured 3/0.4644. A correctly rounded polynomial oracle reached 1/0.2723, showing the remaining gap arises in the evaluation chain. A minimax degree-8 candidate measured 3/0.4670; real-chain descent improved aggregate ln average only 0.235184→0.235132 (0.022%) while raising asinh max 3→4, missed by 100M quick samples. Preserve the existing coefficient tradeoff.

### ln_normal: peel the exact leading term and delay the large rounding
**Shipped.** Write `ln(m)=s+s²Q(s)` and combine the split `k*ln(2)` as `base=fma(k,LN2_LO,s)`, then `fma(k,LN2_HI,base+sq)`. The old nested FMA rounded a nearly full-size `k*ln(2)` before rounding the output; the new form keeps `k*LN2_HI` exact and rounds the full-size result once. Exhaustive positive-normal enumeration moved the normal kernel from avg 0.234319/max 3 to 0.007262/max 1 ulp; restructuring alone gave 0.009211/max 3, while peeling alone gave 0.232666/max 1. In the `k=0` octave the peel, not the combine, moved avg 0.453148 → 0.216130 and max 3 → 1. Score both the reduced octave and the aggregate: the mechanisms are invisible to each other's preferred sample set.

### ln_normal: degree eight was too expensive for its gain
**Rejected.** A peeled degree-eight polynomial improved normal-kernel avg from 0.007262 to 0.006315 and `k=0` avg from 0.216130 to 0.059757, with max still 1. It added three instructions, four uOps, and raised `ln_unchecked` Block RThroughput 16 → 17 (also `log1p` 21 → 22 and `asinh` 42 → 44). On public quick-fuzz rows the improvements were small (`ln` 0.0036 → 0.0031, `log1p` 0.0249 → 0.0239), with unchanged maxima. The shipped degree-seven chain kept throughput essentially flat but added about one serial FMA: `ln_unchecked` latency 34.06 → 38.06 cycles.

### log wrappers: special-case comparisons answer different questions
**Closed.** One comparison decides which exceptional value to produce (−infinity for zero or NaN otherwise); another decides whether the normal result may be used. They cannot simply be consolidated as redundant tests, and prior compound-comparison rewrites compiled to identical code.

### log10: peel the leading term and move the low exponent word into the polynomial
**Shipped.** Fit `log10(m) = s*LOG10_E + s²*Q(s)` with degree-7 `Q`, rather than an unpeeled degree-8 polynomial. Fold `k*LOG10_2_LO` into the polynomial's low group before the final `fma(k, HI, fma(s, LOG10_E, sq))`; the low word is ready early, removing the peel's apparent extra dependency. Across all positive normal f32 inputs the unchecked core improved from average 0.254004 / max 3 ulp to 0.006904 / max 1. Exhaustive `log10` improved 0.127/3 → 0.0034/1, `log10_unchecked` 0.255/3 → 0.0069/1, and `log10p1` 0.2319/3 → 0.1458/2; unchecked throughput improved 1.113 → 1.022 cycles/element (62 → 60 instructions), with no material latency cost. The leading `LOG10_E` rounding imposes a 0.39029-ulp-equivalent fit floor; degree 8 reaches that floor but buys no useful improvement. This placement does not transfer to `log2` (no low word) or `ln` (its leading `s` cannot be hoisted into that group without worsening max error).

### log10_normal: test the new combine independently of an older failed fold
**Finding.** The previous `log10` nested-FMA fold was worse, but the `ln_normal` style delayed combine is different. A stride-251 sweep of all positive normals measured the historical `log10_normal` chain at avg 0.254374/max 2, the restructured combine at 0.009362/max 2, and a peeled degree-seven fit at 0.007301/max 1. It was not shipped in this investigation because its `exp10_checked` dependency domain still needed testing.

### log10_normal: a split leading constant required a refit and cost too much
**Rejected.** Adding `LOG10_E_LO` without refitting the polynomial made sampled average error worse (0.006941 to 0.007195 ulp), because the existing fit already compensates for the one-word constant where it can. Refitting with two words and increasing degree from 7 to 8 improved exhaustive `log10p1` from 0.0387/2 to 0.0304/1 ulp, `log10` from 0.0034/1 to 0.0029/1 and `log10_unchecked` from 0.0069/1 to 0.0058/1. It cost two fmas: throughput rose 8.7% for `log10p1`, 11.7% for `log10`, and 12.4% for `log10_unchecked`; the degree-7 split alone still left `log10p1` at max 2 and cost 4.6–6.2%. A two-word constant is not a drop-in replacement beside a polynomial fitted against its one-word value.

### log1p: a new small branch and coefficient refits
**Rejected.** A branchless small-|x| polynomial added evaluation on every call for avg 0.073→0.068 and max 4→3, but increased throughput 48.1%. A joint output-correction refit looked promising on a coarse grid (-1.6% avg, ln max 3→2), then regressed on 100M real samples from avg/max 0.0966/4→0.1085/7. Direct minimax refits of ln/log10 and an integer `koff` fold offered no gain; LLVM already performed the latter reordering.

### log1p: early reciprocal instead of final correction division
**Rejected.** Rewriting `c/u` as `c*(1/u)` issued early worsened latency 52.19→53.09 cycles (+1.7%) and throughput 2.337→2.374 (+1.6%); avg ulp 0.0966→0.0969 with max 4 unchanged. The longer `ln(u)` evaluation already gave the scheduler slack to hide the division. The same reordering in the now-removed `exp_m1_over_x` did reduce latency, but at accuracy and throughput cost.

### log1p and log2p1: two_sum without using its error term
**Closed.** `two_sum(a,b).0` is literally the same `a+b` rounding, so substituting it for `ln(u)+corr` cannot change output. Preserving more accuracy would require calculating and using a low word, not discarding `two_sum`'s second result.

### log1p, log2p1 and log10p1: omit impossible denormal handling
**Shipped.** A positive rounded `1+x` is never subnormal: the smallest positive value is 2^-24, confirmed over all 2^32 inputs. A wrapper without denormal rescaling keeps zero, negative, infinity and NaN handling and is bit-identical to the old wrapper for all three functions over every input. Throughput fell from 2.289 to 1.857 cycles/element for `log1p` (-18.9%) and 2.328 to 1.886 for `log2p1` (-19.0%); `xlog1py` fell 2.324→1.857. Apparent mca latency regressions for three downstream inverse-error functions were a branch-trace artifact: the scalar baseline's conditional rescale skips work that mca does not model with a predictor; interleaved wall-clock runs favored the new code.

### log1pmx: discard the unreachable signed-zero selection
**Shipped.** Its direct `log1p(x)-x` arm is chosen only when |x|≥0.5, so the `x==0` signed-zero selection cannot affect its output. Sharing the stripped body through `log1p_nonzero!` reduced the throughput region from 155 to 151 instructions (-2.6%) without copying a kernel. The zero, negative and infinity cases of `1+x` remain reachable and are retained.

### log1pmx: correct the reference and reuse its polynomial for the large arm
**Shipped.** The old reference's 11-ulp worst case at `9.53e-7` came from replacing `log1p(v)-v` by only `-v²/2` below `1e-6`; the omitted cubic contributes 10.6 ulp near that cutoff. A series through `v^7` with handover at `1e-3` fixes the reference. Against it the real old max was 8 ulp just past the `|x|=0.5` polynomial boundary. For the large arm, reduce `1+x=2^k m`, set the exact `w=m-1`, and reuse the existing degree-12 `log1pmx(w)` polynomial to reconstruct `ln(1+x)-x`; this removes the separate ln polynomial. Estrin grouping with the leading `1.0` peeled preserves 2-ulp fuzz max while avoiding a 12-deep Horner critical path. Exhaustive avg changed 0.0541 → 0.0632 and max 8 → 3 (fuzz reports 2 after); throughput-region instructions 151 → 133, Block RThroughput 37 → 29, latency cycles 387804 → 371712. Widening the original fit to `x=1` would require roughly four more fma operations; merely reordering the direct arm left max at 8.

### log2p1 and log10p1: split the correction's scaling constant
**Shipped.** For `|x|<2^-24`, the log kernel contributes zero and the scaled correction is the entire answer, so the one-word `LOG2_E` or `LOG10_E` bias is not attenuated. Computing the big product inside an fma with the low-word product as addend cut exhaustive average error from 0.0919 to 0.0324 ulp for `log2p1` and from 0.1458 to 0.0387 for `log10p1`; both maxima stayed 2. Throughput rose 3.1% and 2.1%, respectively. Reversing the fma operands rounds the large product first and is worse; a constant that usually scales a small correction still matters wherever the other term vanishes.

### log_2: Horner and atanh reduction are not free precision
**Rejected.** Full Horner improved avg/max ulp 0.0031/3→0.0019/2 with fewer operations but increased mca latency 55% and throughput 42% due to its serial nine-step chain; a shorter Horner tail retained the same operations and longer critical path than Estrin. The atanh reduction `t=(m-1)/(m+1)` tightened ideal math ~1000× but slightly worsened the single-f32 result and cost 43% latency with a direct division. The form later became useful in double-float `log2_df`, where downstream multiplication amplifies low bits and the division can be replaced by a reciprocal seed and refinement; rejection here applies only to single-f32 `log_2`.

### log_2: denormal bit tricks and unchecked fast path
**Closed.** A constant adjustment to the exponent bit trick cannot replace ×2^24 normalization of denormals: the required shift varies from 1 to 23 bits. Removing a supposedly unnecessary `koff=0` add from unchecked logs was a codegen no-op because LLVM already specializes it away. A non-negative `log_family_edges` variant did remove two instructions from a then-present `powf_throughput` region (217→215), but worsened throughput 1.1% with flat latency.

### log_2: peel the leading term to reduce full-weight rounding
**Shipped.** Separating `log2(e)*s` from the degree-9 polynomial and fitting a degree-8 residual with ulp weighting cut exhaustive max 3→1 and avg 0.0031→0.0028, with one fewer instruction and uOp. Throughput changed −0.1% for `log2` and +6.6% for `log2_unchecked`; latency rose 11.2%. The integer `k` must join last in a separate rounding: folding it into the peeled combine raised avg to 0.1249 (about 40×), while rounding the leading term before the combine defeated the peel. Incorporating `s²` into the low polynomial group retained three Estrin levels; multiplying afterwards added a level and 11% latency.

### logit: use atanh near its central zero
**Shipped.** `ln(p)-log1p(-p)` loses accuracy near p=0.5 because two values near -ln(2) cancel. Using `2*atanh(2p-1)` for `|2p-1|<0.25` reduced exhaustive avg/max ulp 0.2758/1024→0.2625/3, at latency 59.03→59.91 cycles and throughput 3.176→3.551 (+11.8%). A single-domain signed `log1p` quotient scored better avg 0.2408 with max 3, but two serial divisions cost latency 79.19 and throughput 4.028; it also required careful denormal and out-of-domain handling. A near-zero max-ulp blowup can expose real cancellation, not merely a meaningless metric.

## sin, cos, tan and the pi-scaled trig

### cospi: reduce around the nearest integer, then reflect
**Shipped.** Exhaustive avg/max ulp fell 0.2813/868814811→0.0578/2; latency 51→47 cycles, throughput 1.283→1.226 cycles/element, and instructions 51→48. The old `round(x-0.5)` route rounded away low bits in `x-k` near half-integer zeros: at -0.49999997 it returned zero rather than about 9.36e-8. With `r=x-round(x)` exact, the reflected `0.5-|r|` is exact in the half-domain containing those zeros. `cos2pi` inherited avg/max 0.2835/868814811→0.0583/2; the new zeros are +0 and the function is bitwise even. Near-zero ulp catastrophes warrant checking reduction exactness before dismissal.

### Fast sin and cos: allocate more bits to the final pi words
**Shipped.** Keeping `PI_A` and `PI_B` narrow for exact early products while widening `PI_C` to `0x3528885a` and `PI_D` to `0x284234c5` reduced exhaustive |x|≤1e6 maximum error from 58 to 3 ulp for `sin` and 88 to 3 for `cos`; average error fell 12.6% and 6.2%. Only constants changed, with byte-identical instruction streams. The old Sleef split was designed for a much smaller range: its unsummed pi residual was -2.435e-18 versus -1.906e-22 after the change, creating an absolute error proportional to the quadrant near zeros. In-domain maxima of 219 (`sin`) and 2769 (`cos`) still come from incorrect quadrant selection past about 1.3e7, not the pi split. A 100M random fuzz run reported the old |x|≤1e6 max as 4 instead of the exhaustive 58; use exhaustive screening for fast-trig maxima.

### Fast sin and cos: two-product reduction is too expensive for the tier
**Rejected.** Replacing the four-FMA pi reduction by `two_prod(q, PI_HI)` plus low words improved `sin` over |x|≤1e6 from 0.0409/7 to 0.0356/2 avg/max ulp, but cost 1.151→1.278 cycles/element (+11.0%); `cos` cost 1.406→1.651 (+17.4%). It did not move the ~1.3e7 cliff, which comes from wrong `q`, not the multiplication by pi. `PI_TINY` cannot be dropped: the five-operation variant reached max 8758 ulp over |x|≤1e6 and 126119 in-domain because a small absolute error dominates near zeros. An intermediate tier was plausible at 1.278 cycles/element, but not added; score the final trig output rather than isolated residual-relative error near cancellation.

### Quadrant parity: keep the three-operation floating-point formula
**Rejected.** `round_x_over_pi` produces rounded `qh` and `ql` with `vroundps`, not a magic-add result with conveniently placed parity bits. `fma(-2, floor(q*0.5), q)` is three operations and correct even for |q|≥2^24, where every representable integer is even; general exponent-and-mantissa parity extraction needs about seven integer operations. The cheap magic shortcut is valid only below 2^22, but `qh` reaches that magnitude within the former checked-sine range. XOR-based sign combination also costs the same two operations as the current compare-and-select.

### reduce_pi: other shorter forms lost correctness or scheduling
**Rejected.** Balancing its four-deep sum into a two-deep tree added three cycles to both callers, including on rustc 1.98.0-nightly, because the last-ready `e3t` enters sooner rather than last. Replacing e2's `two_prod` with a plain multiply increased sin's in-domain max 2→51,054 ulp at |x|≤1e6; doing that for e3 gave divergent throughput by caller and much worse off-contract error. `round_ties_even` for qh raised cos's max 2→6. Integer bit-op parity was bit-exact but slower; clamping only the polynomial left the raw residual unbounded and restored infinity for finite input.

### reduce_pi_wide: wider gathered elements lost vector width
**Rejected.** Packing the chunks into two `u64` table words was arithmetically exact but changed LLVM's vector factor from 8 to 4; `sin_wide` Block RThroughput rose to 90. The `vpgatherdq` indices were still dwords, so the 8-byte *element* width, not index width, caused the cliff. A fewer-gather table must not assume that a wider element preserves vectorization.

### reduce_pi_wide: instruction reductions that did not reduce cycles
**Shipped one scheduling change; rejected the others.** Combining the three table statics into one and interleaving their lookups removed two instructions and two uOps per region; i5-1145G7 `mca` throughput estimates moved 4.441 → 4.393 (sin), 5.134 → 5.111 (cos), 6.424 → 6.421 (tan), with bit-identical accuracy. Issuing all six long-latency gathers together instead caused queue stalls. Exponent-field scaling reduced Block RThroughput 42 → 36 but slowed sin 4.441 → 4.513; splitting parity added 433 cycles/100 iterations; sharing gather masks reduced Block RThroughput 42 → 37 but increased simulated cycles to 9841/100 iterations. Separate gather-mask rematerialization breaks dependencies, and this loop had scheduling slack: compare total cycles, not just instruction or port counts.

### reduce_pi_wide: replace gathered planes with a 32-bit select-tree window
**Shipped.** On the Ryzen AI Max+ 395 (Zen 5, AVX-512), the preceding two `u64` products per lane made AVX2 and NEON scalarize. A 32-bit reduction word now wraps the integer product modulo two, while f32 FMA recovers its rounded integer part and residual; splitting `pi` into two words and converting the centred word in an exactly representable high part improved the measured sin error (max 3 → 2, average 0.20 → 0.14). Raising the small-input cutoff from 96 to 115 brought the average to 0.127. The gathered intermediate ran sin at 0.62 ns/op versus 0.757 before; freezing the gather index ran at 0.16, identifying the three gathers as roughly three-quarters of the cost.

**Shipped.** Each exponent's table row is a shifted window of the same `1/pi` bit string. Selecting immediate words by exponent bits, then funnel-shifting them in 32-bit lanes, removed the gathers and brought sin/cos to about 0.355 ns/op without changing exhaustive accuracy. `select_unpredictable` was needed: conditional indexing through a const array became gathers, ordinary selects over constants became a lookup table (1.23 ns/op), and `array::from_fn` de-vectorized the loop (4.7 ns/op). This throughput gain increased scalar sin/cos latency from 11.1 to 14.2 ns and tan from 14.1 to 16.6 ns; the select/cmov chain and float-to-integer crossing cost more than a scalar table load.

### reduce_pi_wide: later small instruction experiments
**Rejected.** With the final 32-bit window, changing the low seven word bits to mantissa insertion yielded the same instruction count; alternative mantissa or AVX2 blend-mask bit tricks canonicalized to identical assembly. An AVX2 checked shift worsened the funnel-shift lowering, and dropping the tail word was not accurate near `k*pi` (remainders can be around `2^-29`). Reordering the FMA residual chain increased scalar latency 14.2 → 14.7 ns. A proposed three-word Cody–Waite narrow-sin reduction was not implemented: its second rounding is at the magnitude of `q*P3`, not of the small final residual.

### reduce_pi_wide: three older table-era micro-changes
**Rejected.** On the i5-1145G7 table-based reducer, a shared `2^-28` scale for all `u32` planes made deep planes contribute `2^29` and `2^58` too much (exhaustive average around `7e8` ulp). Selecting the small-input bypass residual in f32 with an explicit NaN select cut instructions 12% but added 2% sin cycles. Replacing `reduce_pi64`'s two `vrndscale` operations with magic-add rounding made no meaningful difference (±0.7%). These were not rerun against the gather-free reducer.

### sin: use a coarse magic grid and round the small correction again
**Shipped.** A multiple-of-four coarse quotient grid followed by a second magic round of its small correction extended the documented domain from `2^22*pi` to `2^24*pi` (5.2707178e7) for one extra add. Throughput barely moved, 1.776→1.778 cycles/element, with unchanged 64-cycle latency; exhaustive scoring of all 2,559,713,206 in-domain f32 patterns gave 0.0457 avg / 2 max ulp. The quotient must remain an exactly representable f32 integer, which makes 2^24 the natural limit for this architecture. A `cvtps2dq` alternative reached the same range but cost 1.906 cycles/element and 75-cycle latency, plus an out-of-range guard. The analogous coarse grid for `cos` was wrong even inside its old domain (0.2806 avg / 2 max to 0.4523 / 205), because its half-integer quotient needs a different correction; `tan` retained its shared reduction with `cos` rather than paying 105→120 instructions for byte-identical results.

### sin and cos: shorter fast reduction and polynomial refits
**Rejected.** Shortening the PI_A..D chain from four to three steps changed avg/max ulp 0.06/220→1.48/866M near sin's zeros. Folding cos's q into a single f32 constant failed representability or rounded near its zeros. A quantized `sinf_poly` refit found no headroom; LP variants worsened cos_checked avg 0.081→0.79 and 0.0495→0.0542 despite isolated improvement, because its residual distribution differs from sin's. Per-caller polynomial copies had the same actual worst residual (2.169e-08 absolute error) and only ~1.25–3× idealized margin.

### sin and cos: replace the PI_D reduction word with a fitted split
**Rejected.** Although the bare reduction residual improved from max 1.19e-7 to 1.04e-7, a 100M-sample real fuzz worsened sin avg/max ulp 0.0645/412→0.1399/1,824,546 and cos 0.2917/2780→0.3320/386,929. Near zeros, relative error exposes the missing precision of the removed pi word.

### sin and cos: leave a one-ulp bound overshoot in the fast tiers
**Rejected.** Within their documented domain, `sin`/`cos` exceed magnitude one on 660/2,720,382 f32 patterns, and `sin_fast`/`cos_fast` on 670/2,720,361; every excess is only 1.0000001. A clamp raised `sin` throughput 1.778 -> 2.033 (+14.3%) and latency 64 -> 72 (+12.5%), and `cos` 1.654 -> 1.883 (+13.8%) and 61 -> 69 (+13.1%). Their doc comments instead direct callers requiring a strict bound to the clamped checked or wide tiers.

### sin, cos and tan: two-word reciprocal pi fixes large-angle reduction
**Shipped.** Across the documented `|x| < 2^22*pi` domain, a single f32 `1/pi` perturbed `x/pi` by as much as 0.17 near the upper limit, sometimes rounding the quotient to the adjacent integer and evaluating the degree-9 sine polynomial outside its fitted `[-pi/2,pi/2]` interval. The magic-round mechanism itself was not the fault; `cos` also lost a rounding by materializing `x/pi-0.5`. A high/low reciprocal-pi product recovers the fraction before rounding the quotient and parity. Exhaustive max ulp changed `sin` 219 → 2, `cos` 2762 → 2, `tan` 3019 → 4; averages 0.0594/0.2902/0.3231 → 0.0422/0.0833/0.1177. On the i5-1145G7, throughput rose 1.151 → 1.776, 1.406 → 1.654, and 2.532 → 3.153 cycles/element respectively; this is still much cheaper than the then-available wide checked reductions. The four pi words in the remainder were already adequate, and a third reciprocal-pi word did not help.

### sin_checked and cos_checked: simplify the wide reduction carefully
Shipped in the then-existing functions. Removing sin's `+0.0` pre-offset through a const-generic path improved `tan_checked` throughput 8.787→8.235 (-6.3%) and `sin_checked` 5.311→5.232 (-1.5%) with unchanged accuracy; the add had been semantically observable at -0.0 and was not eliminated by LLVM. Replacing the final `two_sum(p3,tier2)` with `quick_two_sum` later became bit-identical over all 2^32 inputs and improved throughput 4.854→4.698 for sin and 4.349→4.079 for cos after their safety clamp was removed. Before that neighboring change it had helped sin 1.6% but hurt cos 13.7%, illustrating why scheduling results must be rechecked after surrounding code changes.

### sin_wide and cos_wide: exact wide-range argument reduction
**Shipped.** An f32 is a 24-bit integer times a power of two, so a 256-row, three-f64-word table of `(2^E/pi) mod 2` can reduce its exact value without forming huge `x/pi`. The first word has 29 significant bits so its product with the significand is exact in f64; significant-bit splitting also covers small inputs and subnormals without a bypass. Exhaustive avg/max ulp became 0.1269/2 for `sin_wide` and 0.1501/2 for `cos_wide`; the old `_checked` reduction could return essentially unrelated values at huge magnitudes (max 2,130,706,432 ulp). Against exact-rational probes spanning every exponent, parity had zero mismatches in 120,000 samples and the worst reduced-fraction relative error was 2^-53.

### sin_wide and cos_wide: the f64 gather halves vector width
**Finding.** `sin_wide` throughput/latency measured 8.045/91.05 cycles versus `sin_checked`'s 2.495/82.00; `cos_wide` measured 9.299/99.06 versus 3.157/87.00. The three f64 gathers matter, but LLVM also drops the vectorization factor from eight to four, doubling the cost of other arithmetic. A clamp-less probe with no gather measured 3.286 cycles/element, 32-bit gathers kept width eight at 5.022, and the three-f64-gather version measured 7.540; reducing to two f64 gathers still measured 7.539. A 32-bit-table redesign might avoid the width cliff, but needs extra reconstruction and was not implemented. The wide functions remain separate tiers because their roughly 3x throughput price buys a genuinely unrestricted domain.

### sin_wide, cos_wide, tan_wide: exact u32 chunks replaced f64 tables
**Shipped, then superseded by the gather-free reducer.** Four 27-bit `u32` planes replaced `[f64; 256]` tables: a 24-bit mantissa times a 27-bit chunk is exact in f64, eliminating FMA rounding recovery and `two_sum`. A low-exponent bypass had to select the residual before shared grid/sign processing; fixed-position truncation alone had 22% relative error at exponent 45, and selecting the final output would make small-input `cos_wide` wrong. On the i5-1145G7, vector-factor-8 `mca` throughput estimates fell from 8.045/9.299/13.915 to 4.920/5.434/6.982 for sin/cos/tan; exhaustive averages and maxima stayed 0.1269/2, 0.1501/2, 0.2495/4. Widening chunks to 29 bits cut four planes to three and eight gathers to six per 16 elements; moving the bypass cutoff from 76 to 96 kept seam error around 0.03 ulp (cutoff 90 gave around 0.5 ulp). Throughput then fell to 4.441/5.134/6.424 without changing exhaustive accuracy.

### sin_wide, cos_wide, tan_wide: cutoff 127 and direct small inputs
**Shipped.** Raising the window cutoff from 115 to 127 left only two select levels. Inputs with `|x| < 1` take direct polynomial paths; the same off-window test catches infinities/NaNs, for which `fma(x, 0, x)` supplies NaN while preserving finite `x`, including `-0`. `cos_wide` uses a separate even polynomial near one because the old `pi/2 - |x|` odd-polynomial route was often one ulp low. On Zen 5 AVX-512, sin/cos/tan throughput moved 0.357 → 0.301, 0.358 → 0.309, 0.440 → 0.413 ns/op; AVX2 moved 0.91 → 0.73, 0.90 → 0.79, 1.06 → 0.91. Exhaustive average/max ulps were 0.1269/2 (sin, unchanged), 0.1270/2 (cos, from 0.1432/2), and 0.2307/3 (tan, from 0.2303/3). Latencies were 15.26 → 15.23, 14.76 → 14.52, and 17.00 → 17.17 ns respectively.

### sinf_poly and cos: alternative fitting shapes
**Rejected.** A sinf_poly rational fit was five orders of magnitude worse than the existing polynomial; its idealised LP result was also 0.206 versus the shipped 0.166 ulp-equivalent. A dedicated even cosine polynomial has unbounded relative error near cosine's zeros. Reducing absolute polynomial error cannot meet a relative-accuracy contract at a zero.

### sinpi and cospi: recover argument rounding with two_prod
**Rejected.** A derivative correction left sinpi avg ulp bit-identical and improved cospi avg 0.0861→0.0768 (~11%), but cost throughput +21.8% and +28.8%, respectively. Sinpi's polynomial fit error already dominated its rounding budget; cospi had more recoverable argument error, but not for free.

### sinpi and cospi: fit π into a dedicated polynomial
**Rejected.** A direct `sin(πr)` fit had ~2.6× idealized residual margin yet worsened sinpi avg/max 0.1969/2→0.2065/3 and cospi avg 0.1079→0.1105. Its apparent reduction in cospi's huge near-half-integer outlier did not fix the actual reduction bug, since corrected separately. Rounding of the intermediate construction matters more than a continuous fit's isolated score.

### tan: share reduction with sincos
**Rejected.** Two fused-sincos/direct-tan attempts increased tan avg ulp 0.33→1.05 and max 3,000→32M near poles. Additive combines cancel badly near a residual of π/2; an accurate version would need multiword reduction.

### tan: reduce modulo pi/2 and evaluate a dedicated polynomial
**Rejected.** The original sin/cos ratio was already within 0–1 ulp at actual poles; its roughly 3000-ulp worst cases came instead from the large-x sin/cos domain cliff. The new reduction was worse across its smaller safe domain. Computing paired even/odd polynomials after the same reduction saved only one coefficient degree and increased total operations about 25%.

### tan: sharing sin and cos reductions across inlining
**Finding.** mca measured tan throughput at 2.532 cycles/element, nearly the sum of sin's 1.151 and cos's 1.406, so little common work is shared. Their reductions fuse multiplication with different additive operands, leaving no standalone multiply for LLVM to eliminate; explicit sharing needs different reduction logic and careful rounding-tie handling. The analogous hand-shared `sincos_checked` reduction lost 1.4% on hardware despite mca's predicted 19% win, so mca alone cannot establish a payoff here.

### tan: one reduction loses accuracy at the poles
**Closed.** Computing tan from a single sin/cos reduction could remove one `round_x_over_pi`/`reduce_pi` pair and parity work, but a single-f32 reduced argument loses relative information when cosine approaches zero. Re-expressing cosine as a shifted sine does not recover those lost bits; a double-float remainder would be needed, and the shared checked reduction had previously measured slower. This is only a possible different, explicitly less-accurate tier, not a replacement for tan's contract.

### tan: combine the numerator and denominator signs once
**Shipped.** The sine and cosine reductions share `n`; combine their parity masks on the quotient instead of applying two signs separately. Exhaustive comparison against the old quotient was bit-identical over all f32 patterns. The throughput region fell from 105 to 101 instructions, 107 to 104 uOps, and 3.153 to 2.985 cycles/element (−5.3%); latency stayed 78 cycles.

### tan: dropping the sign mask fails at reduction ties
**Rejected.** Replacing the denominator with `|sin(r_c)|` looked cheaper (105 → 92 instructions, 3.153 → 2.658 cycles/element) but flipped the sign at 12 inputs, all near poles. At `f32(pi/2)`, the f32 reduction rounds `fc` to exactly 0.5; tie-to-even chooses the other sine quotient, so the residuals have the same rather than opposite signs. The result becomes +2.2877334e7 instead of −2.2877334e7. Exhaustive bit comparison, not quick fuzz, is needed to test identities that presume a reduction chooses the mathematically nearest interval.

### tan and tan_wide: reduce once onto a quarter-period pair
**Shipped.** For angles past a quarter turn, select the cos grid and compute the ratio from one small sin/cos polynomial pair rather than evaluating two full-range sine polynomials. `tan_wide` average/max improved 0.243/4 → 0.230/3 ulp and throughput 0.468 → 0.442 ns/op in the first Zen 5 pass; `tan` improved 0.117 → 0.094 average ulp (max 4), 0.254 → 0.190 ns/op, and latency 15.9 → 15.0 ns, while its magnitude range grew from `2^22*pi` to `2^23*pi`. The quadrant choice needs `select_unpredictable`: ordinary branches duplicated the polynomial paths and mispredicted roughly half the time, adding 1.2 ns. A three-term small-path sine over `[-1,1]` raised `tan_wide`'s max from 3 to 4 ulp; four terms retained 3.

### tan_wide: shared table work still costs 3.2x throughput
**Shipped.** Using the wide reducer made `tan_wide` exhaustive avg/max ulp 0.2495/4 across all f32 patterns, versus `tan_checked`'s roughly 406 million avg and 2.3 billion max in quick fuzz. The two reduction calls share their three gathers through compiler CSE rather than doing six, but the checked calls share work too. Corrected mca throughput/latency is 13.915/116.06 versus 4.289/101.00 cycles for `tan_checked`; f64 gathers again reduce vector width. A near-pole argument had previously excused the huge errors, but the wide reducer's worst input was 1.3138148 and its max only 4: the old reduction, not an inherent pole limitation, caused them.

### tanpi: preserve both parity corrections in a shared reduction
**Rejected.** Dropping sinpi's and cospi's parity signs before taking their ratio flips the answer on about half the domain: at x=0.1 the correct ratio is +0.3249, while the unsigned-polynomial ratio is -0.3249. Sharing the rounding operations while preserving the varying parity relationship remains possible, but needs tie-case analysis for a modest expected saving.

### tanpi: dedicated polynomial and pole-distance reflection
**Shipped.** The direct polynomial construction required subtracting the distance to a pole before scaling to radians; subtracting after scaling lost the precision needed near the pole. This ordering was caught by real fuzzing.

### tanpi: peel the leading pi term out of the polynomial
**Shipped.** Fitting the remainder `B(w²)=tan(pi*w)/w-fl(pi)` and closing with `fma(w, PI, w*B(w²))` reduces the fraction of the answer carried by the polynomial; it is not the previously rejected operation of merely folding pi into coefficients. Exhaustive average/max improved from 0.2267/5 to 0.0327/2 ulp, while throughput went from 2.223 to 2.155 cycles/element and instructions from 94 to 93. In the direct arm the remainder carries at most 21.5% of the result, so its evaluation error is attenuated; an ideal-polynomial screen predicted the improvement before refitting. The degree-6 Horner alternative (0.2309/4 with worse latency) is dominated.

### tanpi and tan2pi: add a degree to the shared tangent polynomial
**Shipped.** Correctly rounding the polynomial's target in a screen gave 3.66 max / 0.61 avg ulp against degree 5's 7.85 / 2.05 in the tested chain; a degree-6 fit reached 4.28 / 0.71, while changing the rounded `PI*r` argument was not the main error. Regrouping the new terms in two halves (`lo + u^4*hi`) shipped: exhaustive `tanpi` avg/max 0.3856/8 → 0.2267/5 and `tan2pi` 0.2273/5; throughput 2.031 → 2.223 cycles/element (+9.4%), latency 78.88 → 76.00. A Horner-spine variant reached 4 max but 0.2309 avg and 88.72 latency, too costly for one max ulp. A polynomial degree increase pays off when its fit is actually binding and the callers are limited.

### Wide trig: register-permute x8 prototype and forced autovector width
**Shipped the default vector-width flag; the explicit x8 prototype was not retained.** An experimental AVX-512 x8 reducer reconstructed shifted windows of `1/pi` from one 48-byte constant with register permutes, eliminating gathers; on the i5-1145G7, L1-resident sin/cos throughput was 9.69 → 1.56 and 7.29 → 1.67 ns/element, with bit-identical differential results over 2 million random patterns and pinned exponents. Intrinsics did not scalarize into the portable API. For the then-table-based portable path, `.cargo/config.toml` instead forced vector factor 16: min-of-seven sin/cos/tan throughput went 4.520 → 2.529, 4.730 → 2.746, 5.138 → 2.824 ns/op, with 125/142 throughput rows improved and ten regressions at most 10%. The LLVM flag is toolchain-dependent and global to vectorizable loops in the build; these timings predate the later gather-free scalar-source rewrite.

### Wide trig: sign, zero and polynomial-fitting pitfalls
**Shipped sign fixes; rejected the polynomial refit.** The 2026-09 reducer had returned `-cos(x)` for negative inputs; a Pythagorean-identity test squared away the sign error. `sinf_poly` now preserves negative zero in its own cubic term, removing an external sign-copy operation. Refitting it to eliminate the wide-tier clamp saved two operations and 0.7 ns latency, but biased values near `pi/2`: cos average rose 0.083 → 0.102 ulp, cospi 0.058 → 0.069, and cos_wide 0.143 → 0.151. Optimizing one polynomial consumer cannot be judged without the others.

## asin, acos, atan, atan2

### acos: replace the trailing select with sign algebra
**Rejected.** A reassociation introducing `FRAC_PI_2 - y` differed at 60,142 of 49.6M samples and raised max ulp 6→121 through cancellation near x=0. A variant subtracting only signed constants was bit-identical but raised mca throughput 0.820→0.834 cycles/element (+1.7%). Fewer selects do not guarantee a win.

### acos: refit the polynomial against output ulps
**Shipped.** Weighting approximation error by `sqrt(1-x)/ulp(acos(x))` and accounting for the negative branch gave `acos` max 5→4 ulp and average 0.0650→0.0555 (-14.6%); `acosd` average fell 0.0634→0.0545. Only constants changed, so throughput and latency instruction streams are unchanged. The weighted LP seed provided nearly all the gain; optimizing the real chain moved its average only 0.112046→0.111928 further. Local ±64-ulp coefficient scans could not find this fit (the winning coefficient moved about 220000 ulp); `acospi` has a separate polynomial, with screened idealized error 3.00 against a 1.61 LP optimum (1.62 quantized), but was not changed here. The refit redistributed errors across input bands: the curated worst corpus got 3 ulp worse in aggregate while global average and max improved.

### acos: the asin crossover improvement does not apply
**Closed.** A strided scan found max 3 ulp throughout the positive domain, rather than a concentrated window; the negative side was generally better. Its single `sqrt(1-a)*P(a)` expression has no crossover to move. Earlier polynomial rejections caused by worsening `asin` are no longer binding because `asin` now has its own polynomial, but the best old degree-6-to-7 result improved `acos` avg/max 0.496/4 -> 0.490/3 at +7–11% mca, not enough to take here.

### acos: share the small-argument polynomial across both branches
**Shipped.** For small `|x|`, use `pi/2 - asin(x)`; for large positive `x`, use `2*asin(sqrt((1-x)/2))`. Both evaluate one degree-5 polynomial approximating `asin(sqrt(t))/sqrt(t)` on `t in [0,1/4]`, with different multipliers and addends. A single-fma combine alone had barely changed accuracy (0.1119→0.1118 average, max 4) but reduced throughput cost from 0.820 to 0.771 cycles/element; the shared-poly rewrite lowered in-domain avg/max from 0.1118/4 to 0.00435/2, at 0.961 cycles/element, 56 instructions and 59 uOps versus 43/45 before. Without low words the shared form cost 0.881 cycles/element and yielded 0.0768/2. Degree 6 barely changed the final average (0.00435→0.00433); degree 4 regressed to max 4.

### acos: put the constant's low word inside the final rounding
**Finding.** Adding a low word after `fma(-x, P, FRAC_PI_2)` did nothing: it was below half an ulp of the already-rounded result. Computing `fma(-x, P, PI_2_LO) + FRAC_PI_2` instead moved the shared-poly average from 0.0768 to 0.00435 ulp; the matching `pi` and `pi/2` low words can be generated from one ratio. Adding a two-word `PI` to the *old* polynomial had instead worsened negative-branch average from 0.1191 to 0.9352, because that fit relied on `fl(pi) = 2*fl(pi/2)` near zero. Refit before changing a constant against which the coefficients were tuned.

### acos: changing sqrt scaling did not improve the binding cost
**Rejected.** Replacing `2*sqrt(t)` with the bit-equivalent `sqrt(4t)` shortened the apparent dependency chain, but the polynomial still needed `t`, so computing `4t` duplicated work. Instructions/uOps rose from 56/59 to 58/62, while binding latency only moved 45.98→45.00 cycles and throughput worsened.

### acos_poly: Horner versus Estrin
**Rejected.** Estrin reduced latency but changed fma rounding, raising asin max ulp 9→12 and acos 4→5; retuning raised acos to 6. The shared polynomial must be judged on both callers.

### acos_poly: higher degree and constrained joint refits
**Rejected.** A degree-7 fit improved acos avg/max ulp 0.496/4→0.490/3 but left asin unchanged and cost 7–11% mca throughput on both functions. Joint objectives and LP attempts improved fitted scores but raised acos max 4→5 or asin max 9→12; allowing the pi/2 constant to drift raised acos avg 0.496→1.905. A later coordinate descent improved its 106k-point grid avg 0.14045→0.12982, but exhaustive production evaluation was effectively unchanged (0.0650→0.0649 avg, max 5 in both). The quick-fuzz baseline's apparent max 4 had missed a rare point.

### acos_poly: double-float leading term and pi/2 split
**Rejected.** Splitting the leading term reduced acos avg ulp 0.496→0.068, but increased max 4→5 and worsened asin max 9→12. Retuning recovered asin but left acos max at 6 with extra mca cost. A separate `acos_accurate` combining Df32 pi/2 with a two-product sqrt/poly multiplication brought no improvement after conversion to f32; a transcription error in the shared pi/2 constant was fixed instead.

### acos_poly: rational fitting did not converge
**Rejected.** Two rational degrees failed to converge within 60 seconds; the tested average-first basin hop also produced no material exhaustive improvement. The existing real-chain-tuned coefficients are a better baseline than idealised fit ratios.

### acospi: no asin-style crossover to move
**Closed.** `acospi` uses one expression over the full domain, `0.5 - sqrt(1-a)*P(a)`, with max 3 / avg 0.0438 ulp. Unlike `asinpi`, it has no small branch or crossover window to retune; `acos` has the same structural distinction.

### asin: improve the small branch with one Taylor term
**Rejected.** An extra term genuinely improved avg/max ulp 0.0251/9→0.0202/7, but cost +7.9% throughput with Horner or +12.3% with Estrin. The branch is evaluated unconditionally, so its cost cannot be limited to inputs that need it.

### asin: fit the [0.25,0.5) branch for ulp error
**Rejected.** This band's avg ulp was about 1.4 with max 6–7, but an ulp-weighted Chebyshev LP optimizing the isolated residual (weighted max 2.13→0.10) worsened the real band avg/max 1.426/7→2.046/8 and whole-domain 0.731/7→1.165/8. At the worst point the original isolated fit residual was only ~0.53 ulp-equivalent, not the source of the real error. A combine-sensitivity refit likewise increased avg ulp 13.3%, because its sensitivity vanished near the hard a→1 region.

### asin: route through atan2 and sqrt
**Rejected.** `atan2(x, sqrt((1-x)(1+x)))` improved max ulp 9→4 but nearly doubled avg 0.0506→0.0953. mca latency rose 69.8% and throughput 274%; a binary function's full reduction cost dwarfed the hoped-for gain.

### asin: moving the seam does not cure cancellation
**Rejected.** A prefix/suffix scan of both arms found the same max 6 ulp for every seam in [0.2608490,0.2994175], including the shipped 0.27; the current band average 0.25473 is within 0.0002 of the best. There were 2875 six-ulp inputs between 0.27 and 1, forming a broad big-arm plateau rather than a seam defect. Subtracting roughly 1.266 from 1.5708 around x=0.3 amplifies rounding in the big arm. A direct midrange polynomial would add an unconditionally evaluated branch to a 0.968-cycle/element function, so it was not built.

### asin: fit the constant actually used by the reconstruction
**Finding.** On `[0.27,1)`, replacing the big polynomial with rounded `acos(a)/sqrt(1−a)` still gave max 4 ulp, whereas fitting `(fl(pi/2)−asin(a))/sqrt(1−a)` gave max 2: the stored `pi/2` differs from the exact value by 4.371e−8. Adding an explicit `PI2_LO` to the already tuned chain made max/average worse, 5/0.722 → 6/1.135. Fit against the *implemented* reconstruction; a corrected target also has a 0.3658-ulp-equivalent floor at the endpoint, so extra polynomial degree does not remove this bias.

### asin and asind: move the branch boundary to 0.5
**Shipped.** Every `asin` input with at least 3 ulp error lay in `[0.27, 0.4997]`, where the large branch amplifies subtraction error and `1−|x|` is not guaranteed exact. Move the crossover 0.27 → 0.5, raise the small polynomial from degree 3 to 5 in `x²`, and narrow the large polynomial enough to drop degree 6 → 5. Exhaustive `asin` improved average/max 0.0188/5 → 0.0158/2 ulp and `asind` 0.0222/9 → 0.0195/4; throughput costs 0.900 → 0.961 cycles/element (+6.8%), while the binding large arm's latency fell 40.99 → 36.99. The former estimate that reaching 0.5 required about 12 terms used Taylor convergence rather than an ulp-weighted minimax fit.

### asin and asinpi: fuse the subtractive big-branch product
**Shipped.** Explicit `fma(-s, P, FRAC_PI_2)` avoids rounding `s*P` before subtracting it near `asin`'s crossover; Rust without fast-math does not contract the written subtraction automatically. Exhaustive `asin` max 6 → 5, avg 0.0199 → 0.0188; instructions 57 → 55, throughput 0.968 → 0.900 cycles/element, latency 60.99 → 56.74. The same change to `asinpi`'s `0.5-s*P` gave max 7 → 5, avg 0.2375 → 0.2364, throughput 0.974 → 0.901. A `two_prod` correction matched `asin` accuracy but cost 61 instructions and 1.150 cycles/element; use one fma when the only purpose is to undo a product rounding before addition.

### asin_small: FMA and Estrin rewrites
**Rejected.** Replacing `a*a-a` with `fma(a,a,-a)` was bit-identical but lowered throughput; a small-polynomial Estrin rewrite was not bit-exact and raised asin's max from 9 to 10. Fewer apparent source operations do not establish a machine-level win.

### asind: use a two-word degrees multiplier
**Shipped.** `180.0 / f32::consts::PI` yields `0x42652ee0`, one f32 step below the correctly rounded `180/pi` (`0x42652ee1`); its -0.4606-ulp relative bias accounted for nearly all of `asind`'s average error. Exhaustive avg 0.3378 → 0.0222 with a two-word fma multiplier, while max stayed 9 because the radian `asin` error is amplified by binade geometry near the crossover. Instructions rose 58 → 61, Block RThroughput 14 → 15, latency +6.4%. The high word must be the lower neighbor so the low word is positive: choosing the higher word makes `asind(-0.0)` become `+0.0`.

### asinpi: split the small polynomial's leading reciprocal-pi coefficient
**Shipped.** The `|x|<0.27` polynomial had a lone rounded `1/pi` coefficient, unlike `acospi`'s coefficient-folded polynomial. Peeling its low part into `fma(x, FRAC_1_PI, x*t)` changed exhaustive average 0.2375→0.0159 ulp; max remained 5 at the untouched large-argument branch. Throughput rose 0.901→0.981 cycles/element (+8.9%), 56→61 instructions and 62→68 uOps. Dropping a polynomial degree instead would incur roughly 9 ulp at the branch edge.

### asinpi: the same crossover and degree exchange
**Shipped.** Its ≥3-ulp population also lay in `[0.27,0.5)`. Moving the crossover 0.27 → 0.5, raising the small polynomial degree 3 → 5 and lowering the large degree 6 → 5 improved exhaustive average/max 0.0159/5 → 0.0057/2 ulp. Throughput rose 0.981 → 1.059 cycles/element (+8.0%); the binding latency moved from the large arm at 40.99 cycles to the small arm at 39.99. Unlike `asin`, the pi-scaled small arm now binds, so another term there would directly cost latency.

### atan: compensate reciprocal rounding or use three intervals
**Rejected.** Correcting the `1/a` fold left avg ulp 0.0675 unchanged and worsened max 3→4 at a=1, while throughput rose 1.491→4.196 cycles/element (+181%; atan2 +57%). A three-interval transform confirmed an accuracy opportunity but added a second serial division for half the domain, raising latency 72% and throughput 99%. Screen serial divider dependencies before investing in fit work.

### atan: port atan_latency's sign reassociation
**Rejected.** Applying `mulsign` to the polynomial and pi/2 terms separately was bit-exact but worsened atan throughput 1.491→1.529 cycles/element (+2.5%) and atan2 1.694→1.731 (+2.2%) on re-screening, with flat latency. The same algebra's benefit in division-free `atan_latency` does not transfer to the ordinary reciprocal path.

### atan: joint rational refit does not address the real error
**Closed.** The shipped six free rational coefficients have 0.033 ulp-equivalent exact-arithmetic error; a converged joint fit reached 0.026, against a measured 4-ulp maximum. Division and f32 rounding dominate, and optimizing an idealized rational rather than the real execution chain would not establish a gain.

### atan2: divide in both directions up front
**Rejected.** Computing `y/x` and `x/y` independently passed edge cases and marginally improved same-seed avg ulp 0.068735→0.068683 over 99M samples; max remained 4. Latency stayed ~61.2 cycles but throughput worsened 1.662→1.729 cycles/element (+4.0%) under increased divider pressure.

### atan2: compensate the division residual
**Rejected.** After correcting a NaN-sentinel measurement bug, the apparent tenfold improvement vanished into noise; polynomial error dominated. The extra correction cost +14% latency and +43% throughput.

### atan2: simplify both-zero and pi/2-sign handling
**Closed.** The alternative boolean expression compiled byte-for-byte identically; LLVM's InstCombine already performs it.

### atan2_pos: choose the full-turn fold from the input sign
**Shipped.** For negative y with `y/x` underflowing to -0.0, testing the rounded angle (`r < 0`) missed an entire turn. This affected 420,587/20M uniform-bit-pattern samples, not merely a reference seam: max/avg ulp fell 1,086,918,619/2.29e7→3/0.063 after using y's sign bit; throughput improved 1.993→1.856 cycles/element (-6.9%), latency unchanged at 72.19. The output range is `[0, TAU]`, since correctly rounded f32 TAU can exceed exact 2π. Testing `y<0 || r<0` also preserved one -0.0 convention but cost +2.2% throughput, so was rejected.

### atan_latency: refit or shorten its polynomial
**Rejected.** An LP refit predicted ~6% isolated improvement but measured avg ulp 0.0516 versus 0.0517, indistinguishable. Dropping degree 8→7 idealized to 1.41 ulp-equivalent against real max 3, too little margin for a safe degree reduction.

### atan_poly: refit the rational polynomial
**Rejected.** A denominator-only LP fit barely changed the coefficients (<1% improvement); a numerator-only isolated fit predicted 11–15× but delivered only ~0.9% in the real function, with the same worst point. A seeded 4/4 Padé fit reached ~5e-11–1e-10 continuous maximum absolute error, yet the best real f32 attempt tied max ulp 3 while increasing avg 0.0678→0.0688; other seeds reached max 4. Higher degree and a superior continuous fit do not substitute for tuning the actual fma chain.

### atan_poly: Horner and Chebyshev alternatives
**Rejected.** Switching atan_poly from Horner to Estrin changed average/max 0.068/3→0.072/5 and worsened throughput 156%. Clenshaw evaluation of the latency-oriented degree-17 atan polynomial was screened out: roughly 17 serial recurrence steps would replace a roughly five-deep parallel split. A change of basis must be judged on both dependency depth and final rounded accuracy.

### atan_poly: peel its rational numerator's leading one
**Shipped.** Replacing `x*(1+x²*A(x²))/D(x²)` with `fma(x*x², A(x²), x)/D(x²)` kept the operation count unchanged while saving a full-weight rounding and shortening the dependency path into division. Exhaustive `atan` improved from 0.0675/4 to 0.0627/3 ulp and `atan_bounded` from 0.0526/3 to 0.0431/3; sampled `atan2` improved from 0.0683/3 to 0.0660/3 and `atan2_pos` from 0.0632/3 to 0.0621/3. Throughput instruction/uOp counts stayed unchanged; modeled latency dropped three cycles for `atan`, `atan2` and `atan_bounded`. The same numerator also improved `atand` from 0.0628/4 to 0.0582/3.

### atand and acosd: price the degrees multiplier separately
Shipped for `atand`: the two-word `180/pi` multiplier changed avg/max 0.3682/5 → 0.0628/4 for three more instructions and unchanged Block RThroughput. Rejected for `acosd`: 0.0545/5 → 0.0546/5, because the many tiny inputs round to exactly 90 degrees and `acos`'s own error dominates elsewhere. Do not infer composite gain merely from improving a shared constant.

### inverse-trig combines: add low words of pi constants
**Rejected.** A low-word correction to atan's pi/2 fold left max ulp 3 and avg ~0.0675 unchanged but raised throughput 3.2% for atan and 11.9% for atan2. Asin's large branch regressed avg 0.0199→0.0227 and max 6→7; acos's `+PI` regressed avg 0.0650→0.2722 at an ordinary interior input. An apparent atan2 max 4→3 was only sampling fluctuation across repeated 10M-sample runs. All four sites were reverted.

## Hyperbolic functions and activations

### asinh: removing the large-argument rescale did not improve hardware speed
**Rejected.** Above 2048, `ln(x)+LN_2` approximates `asinh(x)` within 0.0625 ulp at the threshold, falling to 0.0039 at 8192, so the intermediate rescale arm could be removed without changing quick accuracy (0.1493 avg / 3 max). Although this eliminated 11 instructions, 12 uOps and two vector divides, normalized hardware speed changed by only −1.4%, below run-to-run control drift. The actual in-domain square-root/divide/log chain was untouched; removal offered no established speed or accuracy benefit.

### asinh and acosh: the 2048 rescale seam is insensitive
**Rejected.** Sweeping thresholds from 100 to 8000 left average/max ulp (to four decimal places) and worst inputs unchanged for both functions. Their worst inputs, about 0.0155 for asinh and 1.031 for acosh, lie far from the seam; both branches have comparable accuracy throughout the tested range.

### asinh and acosh: replace guarded large-input square roots
**Shipped.** For |x|≥2048, replacing `x*sqrt(1±1/x²)` with `fma(±0.5,1/x,x)` drops a relative next term of at most 7e-15, far below f32 precision. Throughput improved 5.406→4.169 cycles/element for `asinh` (-22.9%) and 4.718→3.569 for `acosh` (-24.4%); latency was effectively flat. Exhaustive accuracy was identical for `asinh` (0.1493 avg / 3 max ulp) and slightly better for `acosh` (avg 0.0597→0.0596, max 4 unchanged). The branch guard, not a new general-purpose sqrt, supplies the accuracy guarantee.

### atanh: use one log1p with a small-input branch
**Shipped.** A signed-x one-log1p formula lost up to 31,303 ulp near -1 because rounding before the singular log amplified error; using `|x|` avoided that failure but alone worsened avg/max 0.0313/3→0.0352/4. Adding a dedicated `|x|<0.25` polynomial corrected its worst point near x=0.111 and improved avg/max to 0.0037/2 with 25.5% lower throughput cost, at +5.9% latency.

### atanh: inline only the reachable log1p behavior
**Shipped.** In its large arm, `v=2|x|/(1-|x|)` is nonnegative for in-domain input and `1+v≥1`, making denormal and zero handling unnecessary; the zero-input arm is discarded by `atanh`'s own select. Out-of-domain and NaN/infinity handling remain. Throughput improved 3.325→2.903 cycles/element (-12.7%) and latency 75.47→69.42 (-8.0%); exhaustive average/max error stayed 0.0037/2 ulp.

### coshm1: use a local small-input polynomial and subtract one for large inputs
**Shipped.** The former `2*sinh_checked(x/2)^2` really did double sinh's relative error: at `[2,4)` its predicted 1.6581 avg / 9 max ulp matched actual 1.6587 / 9. Below `|x|=2`, a degree-3 fit of the residual after `x*(x/2)` avoids cancellation and preserves tiny-subnormal rounding; above it, `cosh_checked(x)-1` has no harmful cancellation. Exhaustive avg/max 0.0864/12 → 0.0287/4 (100M fuzz 0.0832/10 → 0.0287/4). Throughput-region instructions rose 97 → 108 (+11.3%), latency-region cycles 3468 → 3600 (+3.8%); a degree-4 fit halved ideal residual but did not change end-to-end 0.0287/4.

### exp_pos_neg: halve the polynomial constants rather than its outputs
**Shipped.** Scaling four fitted coefficients and the two leading ones by 0.5 removed output multiplies for sinh, cosh and their checked forms while keeping ordinary accuracy effectively unchanged. Applying the factor before exponent-field multiplication also fixed finite large sinh values previously returned as infinity: sinh(88.7228) became 1.7014122e38 and sinh(89.4) became 3.348863e38. Latency fell 2–4 cycles for the affected paths; cosh_checked throughput improved 7.3%, sinh 1.1%, sinh_checked 0.2%, while cosh worsened 0.9%. The narrow variants also improved throughput 6.5% (sinh_narrow, 1.630→1.524) and 5.6% (cosh_narrow, 1.350→1.274), with unchanged corpus and avg/max ulp. An earlier attempt that merely moved a runtime multiply improved latency but worsened sinh/cosh throughput 4.6%/1.4%; compile-time coefficient scaling is the operation-removing version.

### exp_pos_neg: derive the negative half from a reciprocal identity
**Rejected.** Replacing a parallel fma with two multiplies, an fma and a division depending on the positive result increased sinh latency 56→70 cycles and throughput 1.971→2.286 (+16%); sinh_checked rose 58.06→72.06 and 2.345→2.721 (+16%). The mca cost was decisive before accuracy testing.

### exp_pos_neg_checked_half: replace half-multiplies with exponent-bit decrements
**Rejected.** The bit trick was exhaustively bit-exact, but mca throughput improved sinh_checked 2.527→2.346 cycles/element (-7.2%) while worsening cosh_checked 2.212→3.203 (+44.8%). Inlining reschedules the shared kernel separately in each caller; operation counts alone missed the large asymmetry.

### gelu: compensate the rounded argument of erfc, not erfc's kernel
**Shipped.** `gelu(x)=x*0.5*erfc(-x/sqrt2)` previously reached 199 max ulp near `x=-13.09` despite `erfc`'s 7: in the `|x|=8–16` octave, up to 197 ulp came from rounding the argument and only two from the kernel. A split `1/sqrt2` constant plus an fma residual and first-order `erfc(z+dz)` correction reduced exhaustive avg 0.3477→0.2095 and max 199→10, at latency 70.61→74.85 (+6.0%) and throughput 3.217→3.470 (+7.9%). Gate correction to negative `x`, where the tail's sensitivity grows roughly as `2z²`; applying the asymptotic correction on positive `x` added roughly 84 ulp at 13. Correct the bounded `erfc` output before multiplying by `x`, avoiding `inf*0` at positive infinity.

### gelu: its remaining maximum follows erfc
**Closed.** `gelu` measures 0.209 avg / 9 max ulp; replacing only its `erfc` call with a correctly rounded result while retaining all other f32 operations leaves 2.86 max. `erfc`'s relative error passes through the multiplication, sometimes costing twice as many ulps across a binade; improving the surrounding `gelu` arithmetic cannot remove the binding max.

### gelu: keep the Gaussian intermediate normal and square before scaling
**Shipped.** The old `x*norm_cdf(x)` arrangement passed through a denormal probability near `x=-13.10` even though the final product was normal. Reassociating `0.5*|x|` with the Gaussian before multiplication by the `erfcx` factor reduced that band's avg/max from 1.8218/9 to 0.3525/4 ulp. Computing the Gaussian exponent from an accurately formed `x²/2` also removed a rounded `x/sqrt(2)` and its residual-correction machinery; on `[-8,-1]` this independently changed avg/max from 0.8822/7 to 0.8492/6. Exhaustive whole-domain avg/max became 0.1309/6 from 0.1374/9, and throughput improved 3.979→3.434 cycles/unit (-13.7%, 164→140 instructions). A correctly rounded `erfc` alone did not fix the old worst point: check intermediate underflow before blaming the callee.

### gelu: eliminate large-finite NaNs and correct the negative tail
**Shipped.** In the old correction, `(zp + zp)*dz` overflowed to infinity for large `|x|`; multiplying it by an erfc factor of zero produced NaN. The first observed NaN was at `-7.850934e22`, and 10.03% of finite f32 inputs returned NaN instead of zero. The revised form has no such product, returns `-0.0` throughout the saturated negative tail, and deliberately returns `-0.0` at `-inf` as the continuous limit of `x*Phi(x)`; the old infinity override was removed. For composites `x*f(x)` with exponentially decaying `f`, check where the factor becomes denormal separately from where the product underflows.

### logsigmoid: use the already-licensed narrow exponential
**Shipped.** The `softplus`-style correction clamps its exponent argument into [-87,0], inside `exp_narrow`'s single-field domain, so the split field in `exp` was unnecessary. For the then-public softplus and logaddexp, this took throughput 3.006→2.534 cycles/element (-15.7%), while `logsigmoid` fell 2.969→2.663 (-10.3%); their argument constructions were checked for bit-identical behavior. Public `softplus` and `logaddexp` have subsequently been removed, although internal softplus kernels remain.

### logsigmoid and bounded log1p correction: fit the caller's range
**Shipped.** The shared `log1p_unit` correction for `softplus`-style callers receives `e` in (0,1], so `e+e²Q(e)` replaces logarithm reduction, Sterbenz correction and division without changing `exp`. Degree-nine `Q` takes two multiplies and ten FMAs; the rounded-coefficient idealized maximum is 0.078 ulp-equivalent. At this stage `logsigmoid` throughput fell 4.224→2.969 cycles/element (-29.7%); the then-public `softplus` fell 4.100→3.006 (-26.7%), with exhaustive average 0.0833→0.0768 ulp and max 4 unchanged. Weight fits by `e²/ln(1+e)`, not the absolute error of `Q`: for degree eight that changed the estimated result error from 0.541 to 0.086 ulp-equivalent.

### logsigmoid_checked: retain the subnormal exponential tail
**Shipped.** The ordinary function truncates its correction beyond roughly `|x|=87`, while its true output persists through the denormal band to about `x=-103.97`. Replacing `exp_narrow` with `exp_checked` would have raised modeled throughput cost 18.9% (2.604 → 3.095 cycles/element, Block RThroughput 28 → 34). The checked tier uses an exponent offset of 64 and a final exact power-of-two scaling instead; the same approach to `softplus` cost only 2.3%, rather than 21.3% for its full checked exponential.

### sigmoid: one-sided exponent evaluation cost more
**Rejected.** Evaluating exp(+|x|) and choosing `s` or `1-s` preserved special-value pins, including sigmoid(−88), but latency increased 57.06→61.11 cycles (+7.1%) and throughput cost 1.283→1.468 cycles/element (+14.4%). Avoiding negative exponent-field values did not pay for abs and sign-dependent result selection; the candidate was reverted before a full accuracy sweep.

### sigmoid: reciprocal correction targets the wrong side
**Rejected.** The identity `w = 1-e*w` for `w = 1/(1+e)` attenuates reciprocal error when `e < 1`, or `x > 0`. That side's exhaustive max is 1.71 ulp against 3.35 for `x < 0`; even an exact reciprocal of the same f32 exponential attributes at most 1.5 ulp to division. Improving the positive side would not address the limiting exponential error.

### sigmoid_grad: the apparent throughput hotspot was a simulator artifact
**Closed.** Hardware found it faster than `tanh_grad` in all four paired runs (`sigmoid_grad` 0.840–0.903 ns/op against 0.970–1.602), despite mca reporting 4.386 versus 2.076 cycles/element. The function needed no code change; the IPC failure and resource-count cross-check are recorded under Method and tooling.

### silu: retain the average-accuracy checked tier alongside the max-accuracy default
**Finding.** Over the common `|x|≤87` domain, `silu` measured average/max 0.0904/3.85 ulp and `silu_checked` 0.0673/4.69. Their opposite exponential orientations differ on 232,916,327 of 2,237,399,042 scanned patterns, so neither dominates. A tail saturation must select the *result*, not clamp the exponent: otherwise unbounded `x` can multiply a fixed tiny exponential back into range. Use `x.max(-0.0)` where the negative tail must retain signed zero; `max(0.0)` is not signed-zero-safe.

### silu_checked: preserve the tail using a scaled quotient
**Shipped.** `silu`'s true tail persists to about `x=-108.6`, while its sigmoid-based implementation flushes early and can multiply a denormal sigmoid into a normal answer with only 17–22 bits left. With `e2=exp(-|x|)*2^64`, evaluate negative inputs as `(x*e2)/(2^64+e2)` and nonnegative inputs as `(x*2^64)/(2^64+e2)`. The checked region grew 57 → 67 instructions and 60 → 71 uOps, kept Block RThroughput at 17, and changed modeled throughput 1.407 → 1.466 cycles/element and in-domain latency 65.02 → 60.97. Explicitly multiplying by 2^-64 instead was bit-identical but cost four more instructions and raised Block RThroughput to 19.

### sinh and cosh: reassociate their shared kernel
**Rejected.** Dedicated reassociations cut latency around 11–12%, but cosh throughput worsened and sinh avg ulp worsened about twelvefold. Sharing a kernel does not imply both callers have the same accuracy or scheduling outcome.

### sinh and cosh: raise the shared even polynomial by one degree
**Shipped.** A correctly rounded reduced `exp(±r)` oracle cut sampled `cosh` average/max from 0.6649/3.629 to 0.1031/1.164 ulp, identifying the shared polynomial rather than the final sum as the main error. Treat the even/odd halves as one approximation of `exp(r)/2`; degree 5 → 6 needs one extra fma in the even half and reduces idealized error 1.76 → 0.053 ulp. `sinh` exhaustive average/max improved 0.078/4 → 0.061/3 while throughput rose 1.720 → 1.824 cycles/element and latency 51 → 54 cycles; `cosh` throughput rose 1.606 → 1.695. The separate `sinh_throughput` and `cosh_throughput` paths did not change, retaining their speed tier.

### sinh_checked: keep the 0.5 small-argument crossover
**Rejected.** Moving it to 0.3/0.4 changed average 0.0429→0.0439/0.0431 and max 5→7; 0.6/0.7 changed average to 0.0509/0.1019 and max to 28/124. Cosh has no analogous cancellation branch because its exponentials are added, not subtracted.

### sinh_small: same-degree minimax refit and shorter polynomial
**Shipped.** Replacing three Taylor coefficients with a same-degree minimax fit changed exhaustive small-branch avg ulp only 0.02952→0.02951; its max 2→1 did not move sinh's headline max 4–5 in the other branch, so the refit was rejected. Instead, dropping a term to degree 2 spent that branch's spare precision: branch max 2→3, sinh avg 0.073→0.082, headline max still 5, and throughput fell 5.6% for sinh, 6.1% for sinh_throughput and 7.2% for sinh_checked.

### Small hyperbolic polynomial: Estrin was not a free speedup
**Rejected.** Estrin evaluation of `sinh_small` worsened latency, throughput and accuracy. A shorter-looking polynomial dependency structure needs actual codegen and whole-function testing.

### tanh: replace the expm1 formula with one rational
Rejected after fitting. A full-domain `P(x²)/Q(x²)` needed 13 free coefficients to converge; a two-domain split needed 14. Neither offered a simpler route than the existing expm1-based implementation.

### tanh: evaluate a positive/negative exponential ratio
**Rejected.** `(ep-en)/(ep+en)` plus a small-x branch improved avg ulp 0.146→0.024 with unchanged max, but required two full polynomial evaluations rather than the existing single expm1 path. Latency rose 2.7% and throughput 25.8%.

### tanh: correct the final division residual
**Rejected.** A reciprocal-and-fma correction was stopped at the mca screen: latency 82.72→94.36 cycles (+14.1%) and throughput 1.859→2.329 cycles/element (+25.3%). Accuracy was not measured because the extra throughput cost already failed the criterion.

### tanh: changing the quarter-unit seam does not help
**Rejected.** Every seam in [0.2499507,0.3039100], including 0.25, retains max 6 ulp; moving to the average-minimizing 0.2716 changes exhaustive average only 0.1457→0.1456 and worsens two curated edge inputs. Both arms independently reach roughly 6 ulp. A true-value-relative scan initially suggested max 6.084→5.557, but both still round to a six-ulp bit distance in the harness; use the actual `ulp_diff` metric for seam decisions.

### tanh: fit `x/tanh(x)` and divide once
**Shipped.** The old `expm1(2x)/(expm1(2x)+2)` arm fixed the crossover at 0.25, exactly where the direct arm's sensitivity to relative exp error, `1/sinh(2x)`, is 1.919. Fit the small arm as `x/D(x²)` (degree 4) and move its seam to 0.8, where that sensitivity is 0.421; select numerator and denominator before a single division. Exhaustive average/max improved 0.1452/5 → 0.0440/2 ulp over 2,220,710,048 scored patterns. Instructions stayed 77, uOps fell 89 → 88, BlockRT 22 → 21, division instructions 4 → 2 and published latency 85.64 → 63 cycles; mca's contradictory 1.731 → 1.759 cycles/element alone was not accepted as a regression. Fitting the reciprocal preserves the leading `x` exactly, provided the reciprocal is polynomial-friendly over the fit domain; changing an arm's form also invalidates an earlier optimal-seam measurement.

## erf family

### erf: refit near-zero and tail polynomial branches
**Rejected.** Direct refits left max unchanged and shifted avg less than 0.3%. A numerator LP predicted 84% isolated improvement but worsened real avg ulp 0.317→0.325 around the x=0.28 branch crossover. Degree 6→5 in the tail idealized to 6.35 ulp-equivalent, already worse than its own ~4–5 ulp, while erf's worst input was in the near-zero Padé branch. Check crossover neighborhoods and the final construction rather than isolated fits.

### erf: replace the tail subtraction with exp2m1
**Rejected.** The proposed `1-exp2_checked` cancellation fix changed zero bits among 711k samples. Latency fell 18%, but throughput rose 13.6% because the unused exp2m1 Padé branch is paid for on every call.

### erf: changing the 0.28 boundary alongside coefficients
**Closed.** A cheap screen of eight boundary candidates with existing coefficients found a flat plateau around 0.28, consistent with fixed-boundary refits. No full joint optimizer was built, so this is evidence against easy seam headroom rather than proof that every joint fit loses.

### erf_poly: remove its small leading coefficient
**Rejected.** Zeroing a0≈3.4e-5 raised exhaustive avg/max ulp 0.3166/5→2.3367/539. The coefficient is essential, not a negligible fit artifact.

### erfc: rational and exponent refits
**Rejected.** Estrin evaluation of the rational numerator/denominator was speed-neutral and raised avg error 2.7%; a centered variable was ~73 times worse. Exact exponent multiplication via two-product improved avg/max ulp 0.311/109→0.306/93 but cost 5–7% mca throughput. Splitting the exponent by a mantissa mask was exhaustively bit-identical because the subsequent multiply by LOG2_E remained uncompensated. A single log-space degree-15 polynomial converged poorly; borrowing erf_poly yielded about 297,000 avg ulp.

### erfc: target the [0.03,0.25] error bump
**Rejected.** At the worst point fit truncation was 3.96e-7 versus combined rounding 2.22e-7, but a band-biased refit's local grid max 80→79 became a full-domain regression: max/avg 109/0.321→111/0.342. An unbiased coordinate descent against the actual formula produced no meaningful movement. The positive side, avg/max 0.4815/109 versus the negative side's 0.1397/6, is harder largely because the same absolute error corresponds to more ulps near zero output; a sign-specific refit does not repair that.

### erfc: split the domain and narrow saturation
**Rejected.** Splitting [0,10] at 2 made the underlying approximations 100–380 times tighter and improved avg ulp 0.311→0.189, but hardly moved max 106→105: downstream error dominates. True positive output reaches f32 zero only around x=10.05, leaving no useful room inside the current clamp; negative output rounds to 2 around x=3.83, but that side's max is already only 6.

### erfc: combine-sensitivity and compensated evaluation
**Rejected.** A sensitivity-weighted numerator refit reduced avg 4.6% but worsened max 109→129 as the sensitivity vanished at difficult inputs. Compensated Horner initially lowered avg about 14% without moving max because the exponent expression alone could err by 87 ulp. Even after an exponent fix, compensated numerator/denominator evaluation improved an accurate variant only 0.245→0.227 avg and 12→10 max while throughput more than doubled 3.033→6.189 cycles/element (latency 74.00→87.11). Evaluation precision is not a substitute for a better fit.

### erfc: degree-5/5 rational fit
**Rejected.** An LP fit with coordinate-descent polishing improved exhaustive avg ulp across five consumers by about 13–20% (erfc 0.3055→0.2486; erfc_accurate 0.2452→0.1948), but four maxima stayed flat and the targeted erfc_accurate max worsened 13→14. Two extra fmas increased throughput 3.9–6.5% across callers; erfcx_checked latency jumped 41.00→60.22 cycles. The new numerator term was nearly unused; changing denominator degree alone may be a better future experiment.

### erfc: reciprocal-variable fit and exact-square correction
**Shipped.** The old degree-4/4 rational had about 15 ulp *fit* error, so improving its evaluation could not resolve a max near 109. A shared degree-10 `erfcx_pos(x)=v*P(v)`, `v=1/(2+x)`, maps the nonnegative half-line into `(0,1/2]` and builds its inverse-`x` tail into the leading factor. For `erfc`, `p=x*x` and `pe=fma(x,x,-p)` retain the square's low word, and `exp_checked(-p)*(1-pe)` avoids exponent error that would otherwise become relative output error. Exhaustive `erfc` avg/max improved 0.3055/109→0.1993/7; the initial cost was throughput 2.437→3.284 (+34.8%) and latency 64.00→70.36 (+9.9%). A double-float exponent passed through the f32 `exp2_checked_df` kernel instead reached max 23 versus 8 for Cody–Waite on the same coefficients.

### erfc: remove a redundant exponential clamp without changing output bits
**Shipped.** Bounding the squared input so `exp_checked`'s own clamp is provably inactive, with compile-time window assertions, and folding the sign using bits and `mulsign` kept corpus and exhaustive accuracy bit-identical at avg/max 0.1993/7. Latency improved 70.36→62.28 cycles (−11.5%) and throughput 3.284→2.899 (−11.7%); versus before the accuracy rewrite, latency is now lower than 64.00 but throughput remains 19.0% higher. The clamp-free reduction was factored into a macro, and a full assembly diff verified other `exp_checked` callers unchanged.

### erfc: degree and evaluation order cannot be chosen on instruction count alone
**Rejected.** Dropping from degree 10 to 9 raised ideal fit error from 0.542 to at least 2.9 ulp over 12 reciprocal offsets; real-chain max went 4.049→6.058→10.400 for degrees 10→9→8. An even/odd polynomial split saved an op and improved mca latency/throughput 2.0%/3.9%, but exhaustive `erfc` max rose 7→8 and its paired tail kernel rose 6→8 even after retuning coefficients. Horner saved 6.5% throughput but added 13.5% latency through a ten-fma dependency chain.

### erfc: pinning the asymptotic coefficient constrains the entire fit
**Rejected.** For the reciprocal-variable polynomial, forcing its constant to correctly rounded `1/sqrt(pi)` increased `erfc` avg/max from 0.1993/7 to 0.2003/8. The free coefficient was one ulp low and let the other coefficients reduce whole-domain error; its far-tail avg cost was only 0.6245→0.6478 ulp for `x>=20`, with max 4 unchanged. The residual max of the positive-tail kernel was about 5.9 ulp, only about 0.5 attributable to fit: rounding `2+x` and then its reciprocal costs about 2.3 ulp that a refit cannot recover. Exact-residual compensation gained about 1.5 ulp at four to six more ops and was rejected.

### erfc and erfcx: fold the Gaussian reflection into fma
**Shipped.** Moving the sign to the Gaussian factor fused `erfc`'s negative-arm reflection into one fma, improving that arm's max/avg from 2.87/0.216 to 2.61/0.193 with one fewer instruction. `erfcx`'s negative arm likewise fused `exp*(2+2*pe)-r` into `fma(exp, 2+2*pe, -r)`, preserving its infinity behavior. These folds accompanied the reciprocal change; the positive `erfc` arm remained bit-identical.

### erfc and erfcx: the remaining maxima come from both factors
**Finding.** After the reciprocal change, `erfc` over `[1,4]` had max 6.79 ulp, with up to 4.19 attributable to its erfcx factor and 4.19 to its Gaussian path. `erfcx`'s negative-arm worst sample near -0.514 had total max 5.90, with 1.74 from the positive erfcx factor and 5.23 from the Gaussian; its positive arm was max 4.03. The public maxima cannot be explained or removed by tuning `erfcx_pos` alone. A two-branch erfcx fit needed about 18 polynomial degrees versus the shipped 10, both branches evaluated per lane, for limited benefit; reciprocal offsets below two degraded the fit (for example offset 1 gave max 23.6 on a 5M-point grid versus 6.05 at two).

### erfc rational: Chebyshev recurrence is a poor fit for degree four
**Closed.** The existing numerator/denominator are degree 4/4 Horner chains; Clenshaw would not rescue parallelism, and the conditioning advantage of Chebyshev recurrence is primarily relevant at much higher degrees. A separate minimax candidate changed erfc average 0.3055→0.2258 but increased max 109→121, so it was not adopted.

### erfc_inv: avoid forming a tail argument that rounds to one
**Shipped.** Historically `erfc_inv(y)=erfinv(1-y)` returned infinity for `y<2^-24`; a direct ulp sweep found 79.68% non-finite results, avg 837462091/max 1057303482 ulp. Reducing `n=min(y,2-y)` and forming `w=-ln(n*(2-n))` directly avoids the cancellation; an additional polynomial in `sqrt(w)-7` covers `w=16..102.8` beyond the old `erfinv` reach. The new direct sweep measured avg 1.5212/max 16 with no non-finite results. Block RThroughput rose 46 → 58, instructions 171 → 209 and uOps 186 → 270; the former infinity-producing implementation was not retained as a speed tier.

### erfc_inv reference: refine the logarithm of the tail probability
**Shipped.** The new f64 direct-reference path uses Newton refinement against `erf_u10`/`erfc_u15` and accepts separately supplied erfc and erf targets, preserving whichever one is exact in its branch. Newton on `erfc` crawls in the far tail; Newton on `ln(erfc)` converges in about three steps, provided the residual is computed with `log1p` of the relative difference rather than subtracting two numbers near -100. It was validated to 1e-15 relative against SciPy over target values from 1e-45 to 1.

### erfcx: split the asymptotic leading coefficient
**Shipped.** A single f32 `1/sqrt(pi)` coefficient is a near tie between adjacent representable words, so either choice causes a persistent tail bias; the original tail signed mean was -0.63 ulp. A two-word coefficient and final fma brought x≥20 avg/max ulp 0.6478/4→0.2684/2 and signed mean to approximately zero. The fit needed subsequent small coefficient adjustments to retain `erfcx_pos(0)==1.0`.

### erfcx: repair the rational's frozen tail above x=10
**Rejected.** Clamping `erfc_rational` at 10 freezes the value used by erfcx after its exponential factor cancels: relative error grows from ~10% at 11 to ~99% at 20 and ~895% at 100. An asymptotic branch was bit-identical below 10 and limited relative error to 0.00015% up to 200, but an unconditional extra division cost 17.7% mca throughput. The documentation was corrected to describe the freeze.

### erfcx: compute the exponent with double-float precision
**Rejected.** Feeding Df32-precision `x²·LOG2_E` into exp2_checked_df lowered worst-point error 121→20 ulp and domain avg/max 0.20/122→0.15/20, but cost +18.0% throughput and still missed the single-digit-max target.

### erfcx: changing the reciprocal offset does not cure evaluation error
**Rejected.** Refitting `erfcx_pos` with `v=1/(a+x)` improved endpoint conditioning at `a=3`, but its real f32 Estrin chain worsened from max/average 6.16/1.093 ulp at `a=2` to 7.04/1.382 on 3M points. Intermediate coefficients grew and cancelled even though the endpoint condition number fell. Horner improved the original coefficients to 5.35/1.069 but adds ten dependent fmas after a division; an `a=2` LP refit gave 5.73/1.076, too weak a margin against coefficients already tuned on the real chain. Check conditioning throughout the interval and the actual evaluation order, not just at its endpoint.

### erfcx_pos: evaluation, not the fit, limits accuracy
**Finding.** On a dense `[1e-4,12]` scan, its f32 Estrin chain measured 0.730 avg / 5.08 max ulp; an exact residual for the reduced argument gave 0.672 / 4.20, an exact f64 argument still 0.639 / 3.75, and the fit in f64 arithmetic only 0.290 / 0.72. More accurate argument formation or a refit at the same degree cannot remove the polynomial evaluation error. Horner with the plain argument measured 0.660 / 4.36, but loses on max after argument correction and has greater dependency depth.

### erfcx_pos: Horner reordering breaks an evaluation-specific fit
**Rejected.** With shipped coefficients, replacing the degree-10 Estrin tree by Horner made exhaustive `erfc` avg/max 0.1948/7 -> 1.3899/7 and `erfcx` 0.2080/6 -> 1.4028/6; it also broke exact pins `erfc(0) = erfcx(0) = 1`. The coefficients were polished for the Estrin instruction sequence, so reassociation requires a new fit rather than reusing them. Despite fewer instructions, Horner increased erfc/erfcx latency 62.08 -> 83.33 (+34.2%) and 66.99 -> 79.99 (+19.4%) by replacing a four-level tree with ten dependent fmas. The unchanged baseline's `erfc` exhaustive max was 7, not the old readme's 6.

### erfcx_pos: peeling leading coefficients scarcely helps
**Rejected.** Regrouping `v*HI + v*P(v)` to peel `c0` and `c1` improved a simulated shipped-coefficient grid from max/avg 5.336/0.8154 to 5.241/0.8051; a fresh LP fit gave 4.980/0.7932 -> 4.942/0.7829. Even with exact `v`, max only moved 3.371 -> 3.111. The final `v*P` product and addition to `v*HI` still round at full result scale; peeling only attenuates earlier polynomial roundings, unlike the log-family use case. The grid misses exhaustive worst points, so these are relative screening figures, not published maxima.

### erfcx_pos: re-substitute the unrounded argument into its reciprocal
**Shipped.** For `v0 = fl(1/fl(2+xa))`, use `fma(-0.5*xa, v0, 0.5)` at `xa <= 2`, otherwise `v0`; the identity `1/(2+x) = 0.5 - (x/2)/(2+x)` suppresses the denominator's lost low bits while preserving `erfcx_pos(0)` exactly. Exhaustive avg/max ulp changed `erfc` 0.1948/7 -> 0.1289/7, `erfcx` 0.2080/6 -> 0.1439/6, and `gelu` 0.2029/9 -> 0.1433/9. A real-chain simulation gave max/avg 6.223/0.9997 -> 4.187/0.6897, reaching the exact-reciprocal max floor of 4.187. Throughput rose 2.7–6.7% across the measured callers while latency fell 1.1–4.4%; unlike a denominator error-free transform, this needs only two arithmetic instructions plus selection. Above `xa = 2` the identity amplifies error, and clamping its multiplier instead of selecting is not algebraically valid.

### erfcx_pos: unrefitted Horner offers no useful accuracy gain
**Rejected.** With current degree-10 coefficients and the two-word `1/sqrt(pi)` constant, changing only Estrin evaluation to Horner moved isolated max/avg from 3.740/0.6247 to 3.739/0.6209 ulp; end-to-end `erfc` moved 5.252/0.4478 to 5.252/0.4461. The old roughly 0.5-ulp max benefit no longer reproduces, while Horner has 19–34% more latency and depth 10 rather than 4. A meaningful new Horner experiment would require refitting coefficients against its own evaluation order.

### erfinv: inline log1p under a signed argument
**Shipped.** Its `log1p(-x*x)` input guarantees that positive `1-x*x` is not denormal; the explicit |x|=1 override handles the zero case, and the other arm handles x=0. Removing dead guards was bit-identical on all 4,278,190,083 non-NaN patterns; the 16,777,213 NaN patterns now return canonical NaN instead of preserving payload. Combined with the no-denormal wrapper change, `erfinv` throughput fell 4.185→3.321 cycles/element (-20.6%); `erfc_inv` fell 4.372→3.514 (-19.6%).

### erfinv: eliminate its remaining cancellation before refitting
**Finding.** The revised inverse sweep left `erfinv` near avg 0.3853/max 69 ulp, essentially unchanged from 0.3855/max 71. Its worst tail error comes from `1-fl(x*x)`, which incurs absolute error around 2^-25 and relative error about 1.7e-4 at `x=0.99983`. Forming `(1-|x|)*(1+|x|)` instead is the untried route; the polynomial was not shown to be the cause.

### erfinv: preserve negative zero after a leading-term peel
**Shipped.** The central polynomial peel changed `erfinv(-0.0)` to positive zero despite passing ulp sweeps, because a negative correction multiplied by negative zero gives positive zero and `+0 + -0` rounds to positive zero. Evaluating both arms at `|x|` and applying `mulsign` after selecting the arm restores the edge and is bit-identical for nonzero inputs; throughput uOps and block throughput stayed flat. Keep the explicit signed-zero edge test for this and other odd-function peels.

### erfinv: retune the central leading coefficient on the peeled grid
**Shipped.** Writing the leading term as `fma(x, c0-1, x)` puts `c0-1` on an eight-times-finer f32 grid. Moving that stored coefficient four grid steps from its old value reduced exhaustive average from 0.3615 to 0.0472 ulp across 2,130,706,432 patterns; max remained 3. The predicted exponent-uniform averages were 0.3589 before and 0.0478 after. Only `.rodata` changed: the measured assembly regions were byte-identical. A coefficient's mathematically nearest representable value is not necessarily the best value for the fitted full chain.

### erfinv: do not replace the tail fit with an ill-conditioned LP
**Rejected.** Tail LP candidates improved sampled average only from 0.475 to 0.41–0.43 while raising max from 3 to 5–7. In the monomial basis, degree 11 with `t` up to 3 made HiGHS report a 5.11-ulp optimum even though shipped coefficients measured 0.875: that is a conditioning failure, not evidence against the shipped polynomial. A retry needs a better-conditioned basis such as Chebyshev.

### erfinv and erfc_inv: avoid a squared-input rounding in the tail reduction
**Shipped.** Form `1−x²` as `n*(2−n)` with `n=1−|x|`, using one fma, rather than `1−fl(x*x)` plus a correction for the *later* subtraction. On 5,033,165 positive-tail inputs `erfinv` average/max improved 5.0129/71.670 → 4.9566/11.312 ulp, while instructions fell 171 → 166, uOps 186 → 183 and throughput 3.251 → 3.212 cycles/element. The previous correction could not recover information lost in `x*x`.

### erfinv and erfc_inv: recenter the shared tail polynomial
**Shipped.** The degree-8 polynomial in `w=−ln(1−x²)` ranged over `[0.673,16]`; increasing degree in `w` amplified f32 coefficient error through cancellation (up to 15.6× at degree 10). A degree-9 fit in `t=sqrt(w)−1` keeps the amplifier at 1.6–3.4×, with only an exact subtract because the square root was already computed. `erfinv` positive-tail average/max improved 4.9566/11.312 → 1.5565/5.765 ulp; `erfc_inv` sampled tail improved 5.2835/15.555 → 1.4128/5.853. The old fit stopped at `w=15.9424`, the `erfinv` endpoint, but `erfc_inv` reaches 16; a shared fit must cover the union of caller domains. This costs two instructions per call and one BlockRT (about 1.7–2.3%).

### erfinv and erfc_inv: a twelfth tail coefficient crosses the fit cliff
**Shipped.** Idealized tail error for 10/11/12/13 coefficients was 3.32/2.38/0.78/0.84 ulp. Adding two fmas improved strided `erfinv` max 5.0727 -> 3.6714 and `erfc_inv` 5.8861 -> 4.5009; exhaustive `erfinv` avg/max changed 0.3692/6 -> 0.3642/4, with new exhaustive maxima 5 for `erfc_inv` and the removed `probit`. Block RThroughput rose by two (+3.3–4.4%). The central polynomial's shipped nine coefficients already measured 0.464 idealized ulp versus 0.504 for refits at nine or eleven, so central-arm chain rounding, not its fit, was the remaining target.

### erfinv and erfc_inv: refit the far polynomial over its actual range
**Shipped.** The peeled far branch had a fit contribution of about 0.35 average ulp, versus 0.262 from `sqrt(w)` rounding and 0.130 from rounding `w`; its polynomial was not below the chain's error floor. An ulp-weighted L1 fit with a capped maximum, quantizing coefficients sequentially and re-solving the LP, replaced the old fit. Scored across all 10,620,503 f32 values of its input in `[4,10.1285]` with perfect `w`, polynomial-plus-combine avg fell 0.5974→0.2586 ulp and max stayed 2, with the same eight coefficients and instructions. The fitting grid must include the real branch seam near 4.06; ending at 4.0 had produced an 8-ulp extrapolation error despite a 0.43-ulp idealized fit.

## cbrt, sqrt, powers, roots and sRGB

### cbrt: seeds and correction degree
**Rejected.** A degree-2 correction left max 112 ulp even after 41 seed candidates; an integer-division-free seed saved two instructions (5→3) but left max 33. Joint degree-3 seed/constant searches changed results by at most about 0.3%, while a 2/2 rational correction tightened the approximation in isolation about 100× but increased latency 28.5% and slightly worsened real average. The existing correction is limited by final f32 rounding as well as fit quality.

### cbrt: degree-four correction improved accuracy but cost every caller
**Rejected.** An exhaustive 2^32-value check measured cbrt average/max 0.28/3→0.0884/1 and cbrt_unchecked 0.28/3→0.0887/1. Yet latency rose 35.06→39.06 cycles for cbrt and 59.06→63.08 for accurate callers; rcbrt rose 46.30→50.30, with throughput penalties of 2.6–19.6%. The candidate was reverted despite clearing the accuracy budget; its fitting scaffold remained in `examples/tune.rs`.

### cbrt_accurate: cheaper seed plus Halley did not shorten the chain
**Rejected.** After correction of overflow ordering and the large-value rescale threshold, the candidate matched the Newton version bit-for-bit, but latency rose 59.06→75.16 cycles (+27.3%) and throughput cost about 1.2%. The added division and sign reconstruction outweighed the seed-polynomial work saved.

### cbrt_fast: constant tuning does not fix its approximation floor
**Rejected.** A second reciprocal-seed offset improved avg/max 57.4/554 → 53.7/464 with unchanged 28.047-cycle latency, but throughput worsened 0.839 → 0.910 cycles/element, roughly `cbrt_unchecked`'s 0.906. Exact `bits/3` as well gave 52.8/383 but latency 29.047 and throughput 0.854. Fixed-reciprocal Newton steps retain errors tied to a roughly 3% seed; neither version approaches an ulp, so the speed-focused tier kept its original constants.

### Color power curves: strip log guards whose result is discarded
**Shipped.** The power arm of `srgb_to_linear` and `linear_to_srgb` is selected only above the linear toe, where its base is positive and normal. The outer select discards bad-base results below the seam; the infinity/NaN path remains live. Inlining the licensed `powf_pos`/`log_2` composition and removing redundant power overrides was bit-identical over all 2^32 inputs for both functions. Throughput-region instructions fell 163→132 and 164→134; mca throughput fell 4.772→3.666 (-23.2%) and 4.285→3.411 (-20.4%) cycles/element respectively, with no other region changing.

### linear_to_srgb: correcting the outer constant made accuracy worse
**Rejected.** Folding `1.055` into the exponent moved quick-fuzz avg/max from 0.1209/8 to 0.1207/9; adding a two-word correction to `1.055` moved them to 0.1228/9. Although the low f32 constant seems to cost over an ulp near the toe, its error cancels some of the exponent chain's error. The input reduction and `exp2_checked` errors are amplified roughly 3.06 times there; correct composite constants jointly and measure the whole output rather than assuming a more exact isolated constant helps.

### pow_2_3: fit its correction directly instead of squaring cbrt
**Shipped.** Squaring `cbrt(x)` doubled its relative error and added a rounding: exhaustive avg/max 0.7626/7 versus `cbrt`'s 0.281/3. The existing bit seed's residual `r=(s³-a)/a` also permits `s²(1+r)^(-2/3)`; a degree-4 correction and a one-multiply adjustment for the rounded `s²` delivered exhaustive 0.1032/1. The degree-3 fit had 2.03 ideal ulp-equivalent error against degree 4's 0.125; dropping the `s²` correction gave 0.190/2. Throughput-region instructions rose 7500 → 7600, Block RThroughput 16 → 18 (+12.5%), latency cycles 339506 → 352804 (+3.9%). Horner avoided LLVM duplicating the rescaled kernel for denormal/normal branches, which Estrin made jump from 7600 to 12200 instructions.

### powf: selecting an exact-square branch taxes every call
**Rejected.** A branchless `y==2` selection gave exact 0/0 average/max there instead of 5.42/45, but increased general throughput cost 4.3% because both paths execute regardless of exponent. Callers needing a square can write `x*x`.

### rcbrt: inverse seed needs its own three-octave fit
**Rejected.** A direct inverse-cube-root seed eliminated two divisions and improved latency 46.30→31.50 cycles (−32.0%) and throughput cost 1.666→1.511 (−9.3%). Reusing cbrt's polynomial instead changed average/max 0.418/5→8.94/62: calibrating on `[1,2)` and exact powers tested only one of three exponent-mod-3 rounding classes. Across `[1,8)`, even the best exactness-preserving magic constant left max about 23.2; fit a new correction over all three classes before revisiting this architecture.

### rcbrt: seed the reciprocal cube root directly
**Shipped.** `1/cbrt(x)` paid two divisions (one already in `cbrt_normal`) and moved cbrt's relative error to an unfavorable binade, yielding exhaustive avg/max 0.418/5. A downward bit seed for `a^(-1/3)` makes `a*t³-1` an fma rather than a division; a degree-4 refit over its actual residual range `[-0.100965,0.103657]` gave 0.117/1. A degree-3 refit reached 2 max but 0.768 avg, and reusing cbrt's coefficients had about 11 ideal ulp-equivalent error on the new range. Throughput instructions 75 → 76, uOps +6%, simulated cycles +1.4%, latency cycles +1.5%; four vector divides disappeared and Block RThroughput improved 20 → 15. The seed residual repeats every three exponent classes, so its range can be enumerated without scanning all f32 values.

### srgb_to_linear: peel the integer part of the exponent
**Finding.** Computing `b^2.4` as `b*b * b^0.4`, together with an FMA using exact-value rounded constants for `b=(c+0.055)/1.055`, improved exhaustive `[0,1]` accuracy from 0.1046 avg / 14 max to 0.0823 / 7 ulp. Latency fell 110.02→104.30, while throughput rose 3.831→4.289 cycles/element (+12%). Correcting the base constants alone regressed to 0.1339 / 15 because old reciprocal and exponent constant errors canceled; reducing the logarithmic exponent by sixfold and correcting the constants together was necessary. Check for compensating systematic biases before changing one constant in a composite.

### srgb_to_linear_fast: modeled throughput did not survive hardware testing
**Rejected.** The pre-split form looked 10.7% faster in mca (3.831 versus 4.289 cycles/element), but eight alternating hardware rounds found minima of 0.915 ns versus 0.898 for the shipped split form; the in-binary `atan` control was 0.319 for both. Block RThroughput stayed 41 and arm-isolated latency was flat at 69.02 versus 69.08 cycles. With no measured speed benefit and worse accuracy than the split form (which changed max 14 → 7, avg 0.1046 → 0.0823), no `_fast` API was added.

## Complex, geometry and fmod

### clog: strip both guarded log1p wrappers
**Shipped.** In its near-one branch `1+v` is in (0.5,1.5), and in its rescaled branch it is in [1,2], so the general `log1p` exceptional cases and zero-sign fix are dead; Sterbenz correction and `ln_normal` remain. Exhaustive range checks found only a `v=-0.0` difference in the stripped expression, which neither caller can produce. Interleaved quickbench throughput improved from 5.92/5.62/5.73 to 5.66/5.35/5.25 ns over three runs (roughly 5–8%); latency was flat.

### rsqrt and rhypot: Newton-free residual correction
**Rejected.** `e=fma(r,r*x,-1); r_new=fma(-0.5*r,e,r)` improved rsqrt avg ulp 0.2599→0.1226 with max still 1, but raised latency 43% and throughput 8.8%. Applying the same idea to rhypot worsened avg 0.065→0.175: refining a reciprocal against its already-rounded `x²+y²` can erase beneficial error cancellation.

## Removed functions

### asind: absorb degrees conversion into polynomial coefficients
**Rejected.** Folding RAD_TO_DEG into asin coefficients yielded exhaustive avg ulp 0.3376 versus 0.3387 for the simple composite, but raised max 11→16. Scaling output coefficients, despite being algebraically linear, changes fma rounding; the simple composite was retained for the degree inverse-trig functions.

### atan2d: avoid denormal intermediate angles
**Shipped.** Its original max ulp 30 came from underflow in atan2's internal ratio, then magnification by 180/pi; every sample over 8 ulp had a denormal intermediate even though some final outputs were normal. For that branch, `y*((180/pi)/x)` kept intermediates normal and improved avg/max 0.458/30→0.175/5 at +13.1% throughput (1.881→2.128 cycles/element) and +9.7% latency. `(y*(180/pi))/x` instead reached max 77 because the scaled numerator could remain denormal. This association lesson transfers to any conversion multiplying a tiny result by a factor above one.

### atan2d: improve its degree conversion without changing its denormal divisor
Shipped at the time. A two-word `180/pi` multiplier changed sampled avg/max 0.1749/5 → 0.1072/4 for three more instructions and unchanged Block RThroughput. Its denormal division arm retained the old divisor: a rounded single constant improved that arm's sample average 0.0258 → 0.0110 but regressed a pinned edge case by one ulp.

### atanpi and atan2pi: split the final reciprocal-pi multiply
**Shipped.** A high/low `1/pi` conversion lowered `atanpi` average 0.2783→0.0775 ulp and `atan2pi` 0.1383→0.1137, each retaining max 4; the positive low word preserves signed negative zero. Both cost three instructions and three uOps, with one dependent FMA adding four latency cycles. Folding the scale into `atanpi`'s Padé numerator instead gave 0.2432 avg at 1.612 cycles/element, worse in accuracy than the split form's 0.0775 at 1.657: the rational normally shares an unscaled trailing constant between numerator and denominator. `acospi` already folds `1/pi` into its dedicated polynomial and had no terminal multiply to fix.

### atanpi and atan2pi: fold in half-turns before rescaling
**Finding.** These functions are gone, but their measurements show why rescaling a callee's answer cannot remove an inexact quadrant constant already added inside it. Moving `atanpi`'s reflection from `fl(pi/2)-t` to `0.5-h` lowered exhaustive avg from 0.0733 to 0.0308 ulp (max 4 unchanged), at 1.643→1.724 cycles/element. For `atan2pi`, using exact half-turn quadrant constants lowered sampled avg from 0.1119 to 0.0702 (max 3 unchanged), at 1.867→1.907 cycles/unit; the positive arm remained mostly limited by `y/x` rounding (0.1011 of its 0.1177 average). Branch-specific measurements overturned a supposed twofold error floor inferred from `atan2`.

### cbrt_throughput: single-octave tuning did not generalise
**Rejected.** A single-octave grid suggested max 16→12, but real average worsened 6.73→7.01. A 60-octave retry exposed max-first search bias, moving average 6.83→22.73. Seed errors may differ among exponent-alignment classes even if a related cbrt kernel repeats across octaves.

### cbrt_throughput: delete a tier dominated by the ordinary kernel
Shipped at the time. It measured 6.74 avg / 74 max ulp, 42.05 latency and 0.917 throughput cycles/element, versus `cbrt_unchecked`'s 0.282/3, 35.06 and 0.906 on a broader domain; it was deleted. Its inverse-cube-root iteration has four dependent floating levels per step versus two in `cbrt_fast`'s coupled update, so retuning constants could not make it competitive.

### compound: a downstream consumer can erase signed-zero differences
**Finding.** The historical `compound` fed log1p into exp, which maps both signed zeros to 1, so replacing log1p with its nonzero variant was bit-identical and reduced throughput-region instructions 159→155. When considering a wrapper removal, examine distinctions discarded by consumers as well as preconditions established by callers.

### cos_checked: converting the f64 sign to f32 was not free
**Rejected.** `fc as f32` would replace a shift/extract pair with one conversion, saving two instructions, but uOps and BlockRT stayed flat and the critical-path latency rose 87 → 90 cycles, despite a claimed throughput change 3.157 → 3.099. The sign is applied before the polynomial, so this conversion is on its dependency chain. The function has since been removed.

### cos_fast: fixing a second rounding cost too much
**Rejected.** Recomputing the cosine quadrant from a shared nearest-integer reduction rather than `round(x/pi−0.5)` would recover information lost in the latter's double rounding, but increased instructions 57 → 64 and throughput 1.406 → 1.526 cycles/element (+8.5%). The accuracy of that candidate was not measured; at this cost it no longer served its speed-tier purpose. `cos_fast` has since been removed.

### dawson: pin exact small-input behavior before refitting
Shipped at the time. Its central rational's constant term was `1+2^-23`, making every sufficiently small input about 1–2 ulp high; pinning it to 1 and fitting relative rather than absolute error changed harness avg/max 0.805/60 → 0.154/15 at no instruction cost (78 → 77). A tail relative-minimax refit worsened whole-function avg/max 0.1537/15 → 0.1692/16 because most tail samples lie near zero in its transformed coordinate. An exact coefficient is valuable when it makes an entire region exact, but a uniform minimax objective need not optimize a bit-pattern-weighted average.

### dawson: increase both rational and tail degrees
**Finding.** Raising the central rational from `[5/5]` to `[6/6]` and replacing the tail's degree-4 least-squares fit with a degree-5, log-uniform-weighted L1 fit under a max-error cap lowered bit-pattern average 0.1527→0.0561 and max 15→5 ulp. Modeled throughput rose 1.938→1.992 cycles/element (+2.8%) and latency 65.89→68.97. The old `[5/5]` was already near its degree's minimax optimum (13.12 versus 13.07 ulp-equivalent); the old tail's seam alone reached 21.9. Rescaling the rational fit variable to `u/16` avoided an ill-conditioned monomial linear program; matching numerator and denominator degrees preserved favorable SIMD codegen.

### dawson: derive inverse-square from the existing reciprocal
**Shipped.** Its tail used one divide for `1/x²` and another for `1/x`; taking `w=0.5/x`, then evaluating the polynomial in `w²` with rescaled coefficients, removed one divide. The tail region's Block RThroughput fell 30 → 23, modeled throughput 1.992 → 1.701 cycles/element and latency 68.97 → 62.11; whole-function avg error rose only 0.0561 → 0.0585 ulp with max still 5. In divide-bound regions, a coefficient scaling can be much cheaper than a second reciprocal.

### dawson: reciprocal-polynomial orientation requires distant complex zeros
**Rejected.** Although `dawson(x)=x*P(x²)/Q(x²)` seems a candidate for `x/E(x²)`, approximating `E=x/dawson(x)` on `|x|≤4` gave 68 ulp relative error even at polynomial degree 14, versus the existing 12-coefficient rational's 2.8 ulp-equivalent fit. Nearby complex zeros of `dawson` become poles of `E`; a smooth positive real-axis graph does not establish polynomial approximability. The function has since been removed.

### dawson: replace a noisy sampled quadrature reference
**Finding.** An 800-point Simpson reference was expensive and itself noisy by up to 0.74 f32 ulp near its handover. An all-positive series for `|x|≤7` and a 20-term asymptotic expansion beyond it agreed with SciPy within 6.2e−16 relative over the checked range; exhaustive measurement then gave average/max 0.0581/6 ulp, versus the old 2M-sample 0.059/5. When vectorizing such a reference, handle NaN explicitly (`simd_min/max` may discard it), clamp the unselected arm against overflow, and hoist constant reciprocals out of per-lane divisions (quick run 77 → 18 seconds). The function has since been removed.

### dawson: coefficient descent traded average accuracy for a small max gain
**Rejected.** Real-chain coordinate descent over 12 free rational coefficients improved max 5.481 → 5.103 ulp (−6.9%) but worsened average 0.5536 → 0.5753 (+3.9%). A decomposition of the remaining ~5.7-ulp max attributed roughly 2.2 to the fit, 0.8 to rounded `x²`, 1.4 to f32 evaluation and 1.3 to the final multiply/rounding. The function has since been removed.

### dawson: argument compensation needs a derivative
**Rejected.** Its worst errors spread across `[1,4]`, not a crossover-sized window: a strided scan reached 5.50 and 5.54 max ulp in the `2^0` and `2^1` binades. Recovering the residual of `fl(x*x)` takes one fma, but applying it through `x*P(u)/Q(u)` requires derivatives of both polynomials. An exact square would lower the measured max only 4.322 -> 3.517, too little for those extra evaluations.

### dawson: peel the tail's small correction, not its central rational
**Shipped.** The tail rewrite `fma(w*z,T(z),w)` reduced tail average/max from 0.4637/5 to 0.2842/3 ulp, matching its fit-only floor of 0.284190/3; whole-function exhaustive average fell from 0.0581 to 0.0502, with max still 6 in the central branch. It removed a multiply, and throughput improved from 1.701 to 1.645 cycles/element. Rejected for the central branch: its correction grows to 30 times the answer at `x=4`, and even an exact-correction oracle reached 11.90/32 versus the original 0.81/5 in that octave; attenuation must be checked across the entire proposed domain.

### dawson: peel constants inside a rational, not the ratio itself
**Finding.** Peeling each pinned `1.0` in numerator and denominator of its former rational changed exhaustive central-branch avg from 0.43740 to 0.41352 ulp, with max 6 on both sides and unchanged 22.00 Block RThroughput. Peeling the ratio as `x + x*(P/Q-1)` instead caused cancellation because `P/Q` shrank about 30-fold across the branch. Seven inputs still had 6-ulp error after the useful peel, versus 21 before.

### erfcx: reciprocal fit removed a frozen far tail
**Shipped.** The now-removed `erfcx` had frozen its positive tail at the rational's `x=10` clamp, giving about 10% relative error at 11, 99% at 20 and 895% at 100. The reciprocal-variable degree-10 fit removed that clamp and changed exhaustive avg/max 0.3768/126→0.2142/6; an asymptotic reference over `x>=20` gave max 4 and avg 0.6478 across 1.04 billion samples. Initial latency rose 62.02→69.97 and throughput 2.278→2.792 (+22.6%); later clamp removal reduced latency to 66.99, with its throughput simulator reading contradicted by alternating wall-clock trials. A frozen approximation can appear bounded in ordinary sweeps while having unbounded relative tail error.

### exp_m1_over_x: issue reciprocal before the polynomial
**Rejected.** Replacing final `/x` with an early `*(1/x)` improved latency 83→75 cycles (-9.6%) but worsened throughput 1.798→1.842 (+2.4%) and avg/max ulp 0.0729/6→0.0757/7. A divider can overlap other work without making the extra rounding or multiplication free.

### exp_m1_over_x: moving its crossover was not free
**Rejected.** Moving the historical 0.5 seam to 0.52–0.55 improved average 0.0729→0.0728 with max 6 unchanged, but throughput cost rose 1.798→1.822 cycles/element (+1.3%). Wider moves reached average 0.0850 and max 19; a tiny accuracy gain did not justify real execution cost.

### exp_m1_over_x_checked: numerator overflow precedes quotient overflow
**Closed.** A proposed checked tier was not built, and the original `exp_m1_over_x` has since been removed. `exp(x)` overflows at about 88.7228 while `(exp(x)-1)/x` remains representable to about 93.2582 (`x=93` gives roughly 2.64e38), so clamping the numerator cannot make the quotient total. A selected large-x arm could compute `exp(x-ln(x))`; over the narrow [88.7228,93.2582] selection band, `ln(x)` spans only [4.4856,4.5354], making a small fitted polynomial (~4 FMAs, estimated +0.3 cycles/element) more plausible than evaluating a full log on every call. This was an estimate, not a measured implementation.

### hypot: compensate the square-root result
**Rejected.** `e=fma(r,-r,s); r+e/(2r)` made no measurable accuracy difference for an already nearly rounded result, while latency increased 109% and throughput 87%.

### hypot: ordering operands before FMA serialized the path
**Rejected.** Max-first pairing improved average ulp about 15% with max still 1, but hypot latency/throughput costs rose 24%/4% and hypot_checked 1.8%/8.5%. Using `.max()`/`.min()` raised latency 38% and also changed NaN semantics. Sorting operands can cost more than the rounding improvement, especially before a critical-path FMA.

### logaddexp: cancellation makes maximum-ulp A/B comparisons unreliable
**Finding.** Its former fused `log1p_unit` kernel improved average error from 0.155 to 0.141 ulp, and the eight-run average ranges did not overlap (old minimum 0.1509, new maximum 0.1451). Maxima ranged 992–38183 before and 1013–15556 after, with a separate new-code run reaching 249305: a heavy-tailed cancellation maximum does not support a directional claim from a few repeats. For such metrics, report no measurable max change unless repeated ranges separate.

### logaddexp: cancellation needs a more accurate exponential
**Finding.** Near `e^a+e^b=1`, `m+log1p(exp(-|a-b|))` cancels and sampled max swung 1549–27783 ulp according to proximity to the zero curve. Correcting subtraction of `|a-b|` removes only one error term; the roughly `2^-24` relative error of the exponential remains, whereas observed cancellation depth `|m|/|result|` of `1e4–1e5` calls for roughly `2^-40` precision. Equivalent algebraic rewrites and a Newton correction do not remove that requirement. The then-related `softplus` and `logsigmoid` had no analogous cancellation because their summed terms were nonnegative.

### logaddexp_checked: retain representable tails when the larger input is near zero
**Shipped historically.** The ordinary correction cut off beyond `|a−b|=87`; `logaddexp_checked` carried denormals down to about 103.972, restoring 2,224,564 representable f32 corrections in `(87,104]`. It cost 100 → 103 instructions, throughput 2.449 → 2.506 cycles/element (+2.3%) and improved latency 74.110 → 73.345; the tiers were bit-identical throughout `[0,87]`. For a reference near `(0,−88)`, use `m+log1p(exp(−d))`, not `(exp(a)+exp(b)).ln()`: even f64 rounds the latter to zero. The function and its base counterpart have since been removed.

### norm_cdf: correcting the argument was not worth its cost
**Rejected.** At its max-7 worst point, argument rounding contributed only 0.415 ulp, versus 2.931 from `erfcx_pos` evaluation and at most 3 from the exponential factor. Fully carrying a two-word `x/sqrt(2)` residual into `erfcx_pos` changed exhaustive average from 0.0657 to 0.0656 ulp, left max at 7, and increased throughput from 3.436 to 3.811 cycles/unit (10.9%). The callee's sensitivity determines whether a split helps: `erfcx`'s logarithmic condition was -0.729 there, unlike the steep `erfc` used by `gelu`. A residual computed at infinity also produces `inf-inf=NaN`, so such splits require explicit edge handling.

### norm_cdf and norm_pdf: composite argument rounding dominated kernel error
**Shipped.** These removed composites improved exhaustive max ulp 295→8 and 67→4 respectively by changing where their arguments were rounded. For `norm_cdf`, the rounded `x/sqrt2` fed `erfc` an error amplified by roughly `2z²` (about 160 ulp at `x=-12.8`); for `norm_pdf`, rounding `-x*x/2` inside an exponent contributed roughly 70 ulp near `|x|=13`. Square first and halve exactly where possible, retain the multiplication residual with fma, and assess an outer function's log-derivative before attempting to improve its already-accurate kernel.

### norm_cdf and norm_pdf: avoid an expensive public erfcx wrapper for positive arguments
**Finding.** Calling public `erfcx` on an argument known positive increased throughput cost by 37.6% (132 → 198 instructions; `vpslld` count 4 → 8) because the unused negative arm kept a second exponential reduction alive. Calling `erfcx_pos` directly recovered the cost. Inspect generated opcodes for doubled expensive work when composing a branchy public wrapper.

### norm_pdf: split the output scale constant
**Shipped.** A two-word `1/sqrt(2π)` reduced avg/max error 0.0320/4 → 0.0269/3 ulp, at four extra instructions/uOps, Block RThroughput 24 → 25, and latency 62.09 → 66.09 cycles. The single-word constant had a 0.4767-ulp bias that survived at outputs away from zero and the correctly rounded peak. By contrast, the removed `sinc` used the same rounded π in numerator and denominator, so correcting only its denominator would break cancellation (its avg 0.094 was already below `sinpi`'s 0.197).

### normalize2/normalize3 slices: batch API infrastructure is missing
**Closed.** Single-call `normalize2`, `normalize3`, `normalize4`, `hypot4` and `rnorm4` shipped, but the proposed slice-batched normalization kernel needs a slice API tier. No slice variant or measured speedup resulted.

### probit: preserve tiny probabilities in inverse-tail reduction
**Shipped.** Forming `2p-1` rounded to -1 for `p<2^-25`, producing `-inf` for 79.53% of uniform-bit samples despite finite true answers such as `probit(1e-45)≈-14.12`. Routing through the shared half-domain inverse with exact `2p` changed direct avg/max from 832718338/1053729553 to 1.5929/15 ulp and removed non-finite results. Its throughput-region Block RThroughput rose 47 → 60, instructions 176 → 213 and uOps 193 → 275. Near an inverse endpoint, preserve the small input itself rather than subtracting it from one.

### probit: the shared inverse-error tail fit helped its removed caller
**Shipped historically.** The `sqrt(w)−1` tail refit also improved `probit` sampled average/max 5.1376/15.367 → 1.3483/7.194 ulp. Its range reaching `w=16`, rather than `erfinv`'s 15.9424, was one reason the shared fit had to be extended. `probit` has since been removed.

### probit: a tail-only sqrt(2) fold became obsolete
**Rejected.** Folding `sqrt(2)*sqrt(w)` into `sqrt(2*w)` reduced pre-degree-bump `probit` avg/max from 0.9444/7.0831 to 0.8825/6.7963 at +1.6% Block RThroughput. Once the twelfth inverse-erf tail coefficient landed, `probit`'s limiting max moved to the central arm (4.98 strided versus 4.59 in the tail), so the fold bought no max improvement for four added instructions. Recheck a proposed win after changes to a shared dependency; the old baseline may no longer be relevant.

### rootn: keep integer exponents out of the rounded logarithm
Shipped at the time. Forming `log2(|x|)` as f32 lost up to `7.6e-6` absolute near exponent 149, which `exp2` amplified to about 44 ulp. Splitting `|x|=m*2^e` and `e=q*n+rr` kept `q` integral and fed only `(rr+log2(m))/n` to the floating exponent path; exhaustive `n=3` avg/max became 0.15088/2, versus sampled 5.07/43 before, and `n=-1` sampled 10.10/45 → 0/1. Constant-`n` instructions rose 91 → 101 (`n=3`) or 111 (`n=-3`), while dynamic `n` rose 104 → 143 due to integer division. This applies when the outer exponent's integer part survives exactly, not to an arbitrary power such as 2.4.

### sin_checked: narrow-range accuracy is not a reason to slow the fast tier
**Finding.** The two-product fast reduction matched the former `sin_checked` tier's 0.0356 avg / 2 max ulp over |x|≤1e6 at 1.278 versus 5.311 cycles/element, but it did not fix wrong quadrant selection beyond ~1.3e7. The checked function no longer exists; a new intermediate tier would need to justify its extra API surface and 11–17% cost over the fast tier.

### sin_checked and cos_checked: wide-reduction lessons
**Finding.** Their former reductions showed that exact residuals near zeros matter more than a cheaper approximate product: shortening the PI chain raised max 220→866M and dropping e2's error term raised in-domain max 2→51,054. A quick-fuzz `wrap_pi` max varied 1/1/2 versus 12/2/12 across repeats, while exhaustive runs agreed at avg 0.0244/max 41; rare cancellation outliers require matching input sets, not merely repeated random samples.

### sin_checked and cos_checked: magic-round the low quotient parity
**Rejected.** The supposed small remainder actually includes terms proportional to x and reaches about 9.4e30 near f32::MAX, far outside magic rounding's ~2²² range. When the magic constant disappears in the sum, the extracted mantissa bit ceases to represent integer parity: full-range avg error rose to about 3.3e8 ulp and max 2,130,706,432. A large-magnitude path would be necessary before using this trick.

### sin_checked and cos_checked: simplify reduction without changing bits
**Finding.** Three operation-level changes passed an exhaustive 2^32-input bitwise diff, including both reduction helpers: deleting a residual pre-clamp already subsumed by the output `[-1,1]` clamp, using `quick_two_sum` at the final merge, and spelling subtraction as Fast2Diff to avoid a materialized sign flip. Throughput fell 5.232→4.698 (−10.2%) and 4.603→4.079 (−11.4%) cycles/element respectively, with latencies 117.02→108.02 and 122→113. A degree-function path lacking an output clamp could not safely lose its residual pre-clamp; merging two error terms before reduction also rounded away needed precision near large multiples of pi.

### sin_checked and cos_checked: two reciprocal-pi words suffice for quotient rounding
**Finding.** Removing a third `1/pi` correction word and fusing the remaining correction into one FMA lowered throughput 4.698→4.546 (−3.2%) and 4.079→4.037 (−1.0%) cycles/element. Exhaustive results below 1e13 kept the same maxima and worst inputs; an 11.5M-point near-zero probe was bit-identical. The quotient's tiny third word can only alter a half-integer tie, where either quotient/parity representation is consistent until magnitudes already beyond the two-word quotient's useful range. This does *not* justify dropping a third word of pi from the output residual, where its error is amplified near zeros.

### sin_checked and cos_checked: forcing even coarse-quotient parity was not free
**Rejected.** Rounding the coarse quotient to an even integer removed nine instructions and eleven uOps but raised modeled throughput 4.546→4.701 for sine and 4.037→5.020 cycles/element for cosine; sine near-zero max also rose 1.516→2.199 ulp. Moving pi's tiny correction to the full quotient was essential: omitting it gave 5.801 max for both. Independently substituting ties-to-even for the existing round raised cosine's targeted near-zero max 3.460→5.801 while ordinary random-band samples missed it.

### sinc: compensate division after sinpi
**Rejected.** Direct probes found roughly 2 ulp from division even with a perfect sinpi numerator, but replacing `s/denom` with an explicit reciprocal and residual correction cost +12.3% throughput and +1.9% latency. Near denormal x, the reciprocal overflowed to infinity and `tiny_s * inf` produced NaN, whereas direct division remained valid. A compensated reciprocal is not a drop-in replacement across exponent extremes.

### sinc: shared argument rounding beat improving either operand alone
**Rejected.** On `|x|` from 1e-3 to 1e6, the shipped `sinpi(x)/(PI*x)` averaged 0.394 ulp (max 3); replacing only its numerator with a correctly rounded one averaged 0.472, and replacing only its denominator averaged 0.510. Below 0.5, `sinpi` and the denominator share the same rounded `PI*x`, so the ratio behaves like `sin(t)/t` at a common `t`; exhaustive average there was 0.3585 versus 0.4636 even when both operands were independently correctly rounded. Check separate operand oracles before improving one side of a ratio whose errors may cancel.

### sincos_checked: share the angle reduction
**Rejected.** A bit-exact shared reduction passed exhaustive checks and mca predicted roughly 19% throughput improvement, but corrected wall-clock benchmarks showed it about 1.4% slower. This is a concrete counterexample to treating mca throughput as hardware throughput.

### sind and cosd: absorb degrees-to-radians scaling in a polynomial
**Rejected.** Directly fitting the degree-9 polynomial in degree units had ~4.9× idealized margin but worsened sind avg 0.1236→0.1453 and cosd avg/max 0.0725/2→0.1297/3. A stronger continuous fit did not survive f32 evaluation with the irrational scaling constant folded in.

### sind and cosd: improve the degree-to-radian reduction
**Rejected.** Two-product correction gave no measurable accuracy gain and cost +27%/+27.9% throughput. A distinct hi/lo literal split improved exhaustive sind avg ulp 0.1237→0.0675 but cost 11% throughput; cosd gained nothing and cost 13.2%. Test distinct error mechanisms independently rather than generalizing from one failed correction.

### sinpi_unchecked and hypot_unchecked: an unguarded tier needs a measurable win
**Finding.** The deleted variants saved vector-body instructions (48→43 and 18→9), but Block RThroughput stayed at 17 and 12; throughput improved only 1.133→1.117 and 0.766→0.763 cycles/element. In exchange they mishandled `sinpi(-0.0)` and `hypot(±inf, NaN)`. Guards issuing beside the binding polynomial or sqrt chain were effectively free. The former `expm1_checked` documentation also became stale after `expm1` adopted the same field split: its clamp cost 1.087→1.320 cycles/element (+21%), not a throughput improvement.

### softplus: scale the exponential instead of splitting its exponent
**Shipped.** The original tail flushed output at about `x=-87.3` while `softplus(-87.3)` is still a normal 1.2192433e-38, losing the entire denormal output range. Substituting `exp_checked` cost `softplus` 21.3% modeled throughput (2.449 → 2.970 cycles/element), so that route was declined. A checked tier instead offsets the exponent by 64, computes a normal scaled exponential, then multiplies by exact 2^-64: `softplus_checked` costs 103 versus 100 instructions, Block RThroughput 30 versus 28, throughput 2.506 versus 2.449, and latency 73.36 versus 74.11. It was bit-identical to `softplus` for all 2,237,399,042 tested patterns with `|x|<=87`, and scored zero ulp at integer inputs -88 through -103 where the old function flushed.

### softplus and logaddexp: final sum and NaN guard
**Closed.** `two_sum(m,corr).0` is bit-identical to `m+corr`; its unused low word cannot fix the upstream correction error that reaches ~1e4 ulp after near-cancellation. Dropping the NaN guard was rejected: Rust `f32::max/min` ignore a NaN operand (`NaN.max(0)=0`), so a missing operand would silently turn into a plausible finite result rather than propagate NaN.

### softplus and logaddexp: 87 is a negligible-correction cutoff, not a crossover
**Closed.** The accuracy domain ended at |x|<80, since measuring ulps near 87 can turn sub-denormal-scale correction differences into misleading enormous counts. There was no valid seam measurement to optimise; unlike a true crossover, the cutoff marks where the correction is no longer meaningful.

### softplus and logaddexp: guarded kernels can remove whole general-purpose paths
**Finding.** The former public functions' `e∈(0,1]` correction allowed a direct polynomial instead of general `ln` and division; using `exp_narrow` later reduced their then-measured throughput another 15.7%. Public versions have been removed, while the bounded correction and internal softplus implementations remain. The transferable technique is to fit and price a kernel on the caller's guaranteed domain, not on the general function's domain.

### tan_checked: shared reduction parity could be cancelled
**Shipped historically.** Its two checked halves used the same reduction parity; moving the remaining sign to the quotient cut 119 → 108 instructions, 138 → 126 uOps, throughput 4.598 → 4.289 cycles/element (−6.7%) and latency 103.095 → 101, with exhaustive bit identity. Unlike `tan`, this parity cancelled outright. The function has since been removed.

### tan_checked and cos_checked: large-input reduction had a half-integer limit
**Closed.** The removed `tan_checked` inherited the checked sine/cosine reduction's large-input failure: quick-sweep avg/max were 3.3e8/2.34e9 ulp, versus `sin_checked`'s 3.1e8/2.13e9. `cos_checked` lost its half-integer quotient at about `2^47*pi` when a residual's ulp reached one, a binade before sine's roughly `2^48*pi` limit; over `[1e14,8.9e14)`, cosine max was 2.7e11 against sine's 2.3e8, and the two residuals could collapse to the same value. A renormalized two-word residual improved tan's `[1e14,8.9e14)` max from 1.0e11 to 9.5e4 but worsened `[1e9,1e12)` from 535 to 1.97e4, while the dominating `|x|>1e19` errors remained. Correct all-f32 reduction at that scale needs exponent-indexed bits of `1/pi` and more than the roughly 72 bits in the existing pi split, not a cheap quotient or tan-only patch.

### tand: reuse the direct tanpi polynomial
**Rejected.** Reusing sind's reduced argument to measure pole distance caused about 15% relative error and more than a million ulp near odd multiples of 90. A separate cosd-style reduction fixed that (max ulp 3→12), but removed the speed advantage over the original ratio.

### tand: folding independent signs does not remove shared work
**Rejected.** Its sine and cosine parities came from separately rounded degree grids, so folding signs across the quotient removed only three instructions (95 → 92), not a common parity. BlockRT fell 30 → 29.5 but mca throughput rose 2.533 → 2.655 cycles/element; inlining both reductions would duplicate code for little benefit. The function has since been removed.

### tand: use signed pole distance after an inexact period reduction
**Shipped.** A direct degree-6 peeled tangent polynomial improved average/max from 0.1773/3 to 0.0383/2 ulp. Because `round(x/180)` can yield a remainder slightly beyond 90 degrees, `90-|d|` loses the reflected answer's sign; `mulsign(90,d)-d` preserves it without another reduction. Paired wall-clock runs over a full period found throughput 1.960 to 1.774 and latency 51.62 to 41.13 in the measured units, despite increases in instruction and uOp counts. The signed pole distance generalizes to reflected functions fed by an inexact quotient.

### wrap_pi: near-zero ulp maxima need a wider reference
**Finding.** Over |x|≤10000, historical wrap_pi measured average/max 0.0244/41 ulp; the worst x≈8953.539 reduced nearly to zero and values away from multiples of 2π stayed within 0.41 ulp. A one-word f64 tau reference instead reported max 66; its reduction differed from a two-word reference by 24.56 ulp at that input.

### wrap_pi: signed zero was lost in reduction
**Shipped.** The now-removed `wrap_pi` needed an `x==0` guard because subtracting signed zeros in its reduction turned `-0` into `+0`; throughput changed 4.103→4.280 cycles/element (+4.3%), latency remained about 92. Within the exact-identity range of 2.14 billion tested inputs there were three anomalies: `-0` and both rounded `±pi/2` branch boundaries. Preserve signed zero before reduction when the function's contract requires it, rather than trying to reconstruct the sign afterwards.

### wrap_pi: exact products matter in a reference reduction
**Finding.** Its old f64 reference rounded `q*tau_hi` before cancellation and falsely reported 41 max ulp; splitting tau so each product was exact reduced the reported max to 1. A two-word pi fold lowered the misrounding rate 0.0244 → 0.0043 but raised latency 92.02 → 115.02 with reported max still 1; reducing modulo the full period instead put a rounding tie on the branch cut and produced 2.16e9 max ulp. Reference splits must ensure exact products, and reduction ties must not coincide with discontinuities.

### wrap_pi: clamp to the representable interval bound
**Finding.** Exhaustive inspection found 1,273,675,032 finite f32 inputs outside the promised interval `(-pi,pi]`, including every input beyond about 2.83e22, and sixteen results exactly `-f32::consts::PI` at the excluded endpoint. A clamp to `WRAP_PI_MAX` (`0x40490fda`, the largest f32 below pi), with a cheaper copysign-based half-turn fold, restored the range; exhaustive accuracy stayed 0.0244 avg / 1 max ulp. Throughput changed 3.756→3.854 cycles/element (+2.6%) and latency 91.111→96.017 (+5.4%). Half-open mathematical intervals need an explicit representable endpoint; edge tests limited to 1e9 and a 1e-6 tolerance missed both failures.
