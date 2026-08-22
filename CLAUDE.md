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

This crate computes in f32. f64 appears only inside accuracy tiers
(`*_checked`, `*_wide`, `*_accurate`, `powf`) whose documented contracts
require >24-bit intermediates. Perf-stat cycle counts on this machine
(i5-1145G7) against llvm-mca established:

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

Therefore: (a) never introduce f64 into an f32 hot path unless the
function's documented contract needs >24-bit intermediates, and say in a
comment exactly which step needs the width; (b) never widen a table's
element type past dword; (c) a change touching an f64 path needs a
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
