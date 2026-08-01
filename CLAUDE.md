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
`denormal_rescale`, `pade_expm1_ratio` are used by 6-94 public functions each.
Nobody owns them. To edit one: `./tools/jm core claim`, which locks out
everyone. A core edit moves every function downstream of it, so re-run the
accuracy sweep for the **whole crate** before landing, not just your domain.

### Measurement

- **`llvm-mca` needs no lock.** It is static analysis, and each worktree has
  its own `target/`, so nobody clobbers your `mca_target-*.s`.
- **Never set `CARGO_TARGET_DIR`.** A shared target dir puts every instance's
  `--emit=asm` output at the same path and silently corrupts asm comparisons.
  `./tools/jm doctor` checks this.
- **Wall-clock timing must be serialised:** `./tools/jm bench -- cargo run
  --release --example quickbench`. Another instance compiling on the same 8
  cores makes a timing run meaningless.
- Accuracy sweeps use half the cores and self-nice; two at once is fine, three
  thrashes.

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
