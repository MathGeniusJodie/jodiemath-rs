# jodiemath-rs

## You are one of several Claude instances working on this crate

Others are editing `src/lib.rs` right now, in their own git worktrees. Two
instances tuning the same polynomial silently overwrite each other, and two
instances tuning functions that *share* a kernel corrupt each other's
measurements without either noticing. `tools/jm` prevents both.

**Before editing `src/lib.rs`:**

```sh
./tools/jm next            # what is free
./tools/jm claim <fn>      # e.g. ./tools/jm claim erfc
```

`claim` takes any function name and locks the whole **domain** that function
belongs to -- every public function sharing its private helpers. Claiming
`erfc` also gets you `erfcx` and `exp_checked`, because all three compile
`erfcx_pos`. That is the point: they cannot be worked on independently.

If a domain is taken, **pick another one**. Do not wait, and do not edit it
anyway.

**Before every commit:**

```sh
./tools/jm check           # non-zero if your diff leaves your claim
```

**When you are done with a function** (not when the session ends):

```sh
./tools/jm release
```

### Shared kernels are frozen unless you hold the core lock

`fma`, `sinf_poly`, `sinf_poly_raw`, `parity`, `exp2_field_split`,
`exp2int_field`, `exp2_q_poly`, `exp_r_poly`, `log_family_edges`,
`denormal_rescale` are used by 6-94 public functions each. Nobody owns them.
To edit one: `./tools/jm core claim`, which locks out everyone. A core edit
moves every function downstream of it, so re-run the accuracy sweep for the
**whole crate** before landing, not just your domain.

This prose list goes stale; `./tools/jm core` computes the real one from the
source. Trust the tool over this paragraph.

### f64 in hot paths: hard rules (measured 2026-08)

This crate computes in f32. Wider intermediates come in two flavors:
f64, and double-float/`Df32`+fma (error-free transforms built from f32,
see `src/doublefloat.rs`).

**Default for extra precision is `Df32`/fma, not f64.** Reaching for f64
requires a recorded head-to-head measurement against a *good* attempt at
the `Df32`/fma version -- fitted, passing its accuracy gate, not a straw
man. Every f64 site currently in the crate carries such a test in
graveyard.md; do not re-run them, and do not add a new f64 site without
one:

| f64 site | Df32/fma attempt tested | outcome |
|---|---|---|
| `reduce_pi64` (sin/cos/tan_checked, reduce_pi_checked, wrap_pi) | yes: double-f32 EFT reduction | Df32 ~2x slower, 3 decades less range |
| `powf`/`powf_pos`/`powf_unchecked` (`log2_f64`, `exp2_f64_to_f32`) | yes: `log2_df` + `Df32*f32` + `exp2_checked_df` | f64 -5..17% cyc AND max ulp 3 -> 1 |
| `compound_accurate` (`log2p1_f64`) | yes: `log2p1_df` chain | f64 -33% cyc, -21% latency, max ulp 5 -> 1 |
| `remainder_wide` | yes: Df32 chain | f64 ~-70% instructions ("Df32 doing f64's job") |
| `logaddexp_accurate` | yes: Df32 sketch analyzed | f64 more accurate (2^-53 vs ~2^-47) and cheaper to write |
| `clog` near-unit `v` | yes: f32 predecessor failed; Df32 form priced | f64 fragment kept (same throughput, fewer uops) |
| `reduce_pi_wide` | n/a | gathers dominate; arithmetic is not the lever |

Why f64 keeps winning here despite costing 2x per lane (measured directly:
`vfmadd213pd` RT 1.0 vs `vfmadd213ps` RT 0.5, both 8 lanes -- see
graveyard "powf: the Df32 chain was doing f64's job"): an f64 rewrite only
wins when it *shortens the algorithm*, and dropping error-free transforms
does exactly that. `Df32` stays the right default when the wide step is
short (one or two ops where conversions would dominate), when inputs are
already split (f32 inputs enter `Df32` for free), or when f64's exponent
range would mask an overflow that the contract wants caught.

Additional hazards measured on this machine:
- Packed/vectorized f64 arithmetic performs **as mca predicts** (real ≈
  0.93-1.27x simulated across powf, logaddexp_accurate,
  compound_accurate, remainder_wide, sin_checked). Do not "optimize away"
  f64 math on the theory that it is secretly slow -- it is not.
- The hazards that ARE real:
  1. **f64/qword tables and f64 gathers** collapse surrounding loops from
     VF=8 to VF=4 and then run ~4.5x worse than simulated (sin_wide;
     graveyard "the gather: 3.2x throughput"). Dword `u32` chunk tables
     are the only gather shape allowed.
  2. **Scalar f64** costs ~15x the packed form (powf measured both ways).
     Any per-lane branch can cause this; keep f64 paths branchless
     (selects/blends), and split `_checked`/`_unchecked` tiers rather
     than branching per call.
  3. **f64 divides** are effectively unpipelined; use them only where the
     contract needs exact division (remainder_wide), never as a shortcut
     reciprocal.

Therefore: (a) new precision extensions default to `Df32`/fma; shipping
f64 instead requires the head-to-head table above to gain a row with f64
winning on cycles AND accuracy; (b) never widen a table's element type
past dword; (c) keep every f64 hot path branchless/select-based so it
stays packed, splitting `_checked`/`_unchecked` tiers rather than
branching per lane; (d) f64 divides only where exact division is the
contract (remainder_wide); (e) any change touching an f64 path needs a
perf-stat cycle count against its mca region before landing (recipe in
IDEAS.md "mca vs reality"; wall-clock ns cannot resolve <10% on this
machine, and llvm-mca alone cannot be trusted for the wide/gather tiers).

### Measurement

- **`llvm-mca` needs no lock.** It is static analysis, and each worktree has
  its own `target/`, so nobody clobbers your `mca_target-*.s`.
- **Iterate with `tools/mca_region.py`, not the full harness.** It extracts
  just the `LLVM-MCA-BEGIN`/`END` regions you name and reproduces
  `examples/mca`'s numbers to the digit in seconds instead of ~20 minutes, and
  it prints instrs / uOps / Block RThroughput together — the whole escalation
  ladder, which you need because mca's throughput column alone has been wrong
  in both directions here.

      python3 tools/mca_region.py \
        "$(ls -t target/release/examples/mca_target-*.s | head -1)" \
        asin_throughput asin_latency
- **Never set `CARGO_TARGET_DIR`.** A shared target dir puts every instance's
  `--emit=asm` output at the same path and silently corrupts asm comparisons.
  `./tools/jm doctor` checks this.
- **Wall-clock timing must be serialised:** `./tools/jm bench -- cargo run
  --release --example quickbench`. Another instance compiling on the same 8
  cores makes a timing run meaningless.
- Accuracy sweeps use half the cores and self-nice; two at once is fine, three
  thrashes. This machine also runs unrelated multi-core jobs for hours at a
  time, so wall-clock duration tells you nothing — one more reason perf
  conclusions come from `llvm-mca` only.
- **Do not poll for a background job with `until ! pgrep -f "<cmd>"`.** The
  pattern matches the waiting shell's own command line, so the loop waits on
  itself forever — it has already hung an instance for hours. Run the job with
  `run_in_background` and let the completion notification arrive.

### If this worktree holds work that is not yours, stop

One worktree is meant to have one instance in it, but a session that was
waiting on a long background job can wake up after a newer instance has taken
over its worktree. If you find uncommitted changes you did not make, or a
claim held by this worktree that you did not take:

**Report it and change nothing.** Do not commit, `git checkout`, stash,
revert, or `jm release` any of it — every one of those destroys in-progress
work belonging to the current occupant. `jm status` shows which claims this
worktree holds.

### Never use `git stash`

`refs/stash` lives in the shared `.git` dir, **not** per worktree. Every
instance pushes onto one stack, so a mistimed `git stash pop` silently applies
another instance's work into your tree. Use `./tools/jm stash` /
`./tools/jm unstash`, which keep the patch in your worktree-private git dir.
(Tracked files only; untracked files are left in place and reported.)

### Landing your work

```sh
./tools/jm land            # rebase onto master, then fast-forward master
```

`land` serialises on a lock, so simultaneous landings queue instead of racing.
Never `git merge` into master by hand and never push -- this repo is local-only
and far ahead of `origin`.

`IDEAS.md`, `readme.md` and `examples/*` are shared by everyone and are the
usual source of rebase conflicts. Keep edits to them minimal and land often.
`graveyard.md` is set to union-merge, so appending to it never conflicts.

## Working on this crate

Read `IDEAS.md` for what is open and **`graveyard.md` before proposing
anything** -- this crate has a long history of measured negative results, and
re-running one is the most common way to waste a session. Both files also
document the screening order that has saved the most time (oracle screen before
refitting a poly; instruction count -> uOps -> `Block RThroughput` before
believing an mca throughput delta).
