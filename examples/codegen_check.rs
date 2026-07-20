// Permanent codegen regression guard, from IDEAS.md's "CI asm grep" idea:
// compiles examples/mca_target.rs to assembly (same pipeline as mca.rs) and
// greps each `*_throughput` region -- the vectorized, real-public-function
// loop -- for the two specific classes of silent de-vectorization this
// crate has actually been bitten by before (see jodiemath-workflow memory):
//   - a scalar `call` instruction (libm fallback, or a de-vectorized
//     per-element function call the auto-vectorizer gave up on)
//   - `cvttsd2si`/`cvttss2si` (the saturating-cast incident: a `as` numeric
//     cast lowered to a scalar per-lane extract+convert instead of a single
//     packed instruction, silently de-vectorizing an otherwise-packed loop)
// and confirms at least one packed (non-scalar-suffixed) SIMD arithmetic
// instruction is actually present, so an empty/optimized-away region
// doesn't pass by vacuous truth.
//
// This is a coarse, best-effort net, not a full codegen verifier: it
// doesn't (and can't cheaply) prove the *whole* region vectorized cleanly,
// only that these specific known failure signatures aren't present. Treat
// a clean run as "no known regression class detected," not "codegen is
// optimal."
//
// Also checks the separate `pown_const*_latency` regions (idea #153):
// `pown_const`'s own doc comment says to verify with `--emit=asm` that its
// const-generic exponent fully constant-folds the 32-iteration bit-testing
// loop away instead of trusting it as advice -- this makes that a standing
// assertion (no branch/call/loop instruction in the region) instead.
use std::path::PathBuf;
use std::process::Command;

fn main() {
    let manifest_dir = env!("CARGO_MANIFEST_DIR");

    let target_src = PathBuf::from(manifest_dir).join("examples/mca_target.rs");
    std::fs::File::open(&target_src)
        .and_then(|f| f.set_modified(std::time::SystemTime::now()))
        .expect("couldn't touch examples/mca_target.rs to force a rebuild");

    eprintln!("compiling examples/mca_target.rs to assembly...");
    let status = Command::new("cargo")
        .current_dir(manifest_dir)
        .args(["rustc", "--release", "--example", "mca_target", "--", "--emit=asm", "-C", "debuginfo=0"])
        .status()
        .expect("failed to run `cargo rustc` -- is cargo on PATH?");
    if !status.success() {
        eprintln!("cargo rustc failed, aborting");
        std::process::exit(1);
    }

    let examples_dir = PathBuf::from(manifest_dir).join("target/release/examples");
    let asm_path = std::fs::read_dir(&examples_dir)
        .expect("couldn't read target/release/examples")
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| {
            p.file_name().and_then(|n| n.to_str()).is_some_and(|n| n.starts_with("mca_target-") && n.ends_with(".s"))
        })
        .filter_map(|p| std::fs::metadata(&p).and_then(|m| m.modified()).ok().map(|t| (t, p)))
        .max_by_key(|(t, _)| *t)
        .map(|(_, p)| p)
        .expect("no mca_target-*.s found after `cargo rustc --emit=asm` -- did the example build?");

    let text = std::fs::read_to_string(&asm_path).expect("couldn't read the generated .s file");
    let lines: Vec<&str> = text.lines().collect();

    // Collect (region_name, [instruction lines]) for every *_throughput
    // region, the same BEGIN/END marker convention mca.rs's own parser
    // relies on.
    let mut regions: Vec<(String, Vec<&str>)> = Vec::new();
    // idea #153: pown_const<N>'s own doc comment says "verify with
    // --emit=asm before trusting" its 32-iteration bit-testing loop fully
    // constant-folds away for a given N -- standing test instead of manual
    // advice, so a future toolchain/LLVM change that silently breaks the
    // fold gets caught. Collected separately (different naming pattern and
    // a different check below) from the *_throughput regions.
    let mut pown_const_regions: Vec<(String, Vec<&str>)> = Vec::new();
    let mut current: Option<(String, Vec<&str>)> = None;
    let mut current_is_pown_const = false;
    for &line in &lines {
        let trimmed = line.trim();
        if let Some(name) = trimmed.strip_prefix("# LLVM-MCA-BEGIN ") {
            if name.ends_with("_throughput") {
                current = Some((name.to_string(), Vec::new()));
                current_is_pown_const = false;
            } else if name.starts_with("pown_const") {
                current = Some((name.to_string(), Vec::new()));
                current_is_pown_const = true;
            }
        } else if trimmed == "# LLVM-MCA-END" {
            if let Some(region) = current.take() {
                if current_is_pown_const {
                    pown_const_regions.push(region);
                } else {
                    regions.push(region);
                }
            }
        } else if let Some((_, body)) = current.as_mut() {
            body.push(line);
        }
    }
    assert!(!regions.is_empty(), "found no *_throughput regions -- did the BEGIN/END marker format change?");
    assert!(
        !pown_const_regions.is_empty(),
        "found no pown_const* regions -- did examples/mca_target.rs drop its lat_pown_const* wiring?"
    );

    let packed_simd = ["ymm", "zmm", "xmm"]; // xmm still packed (128-bit); scalar forms use an "ss"/"sd" mnemonic suffix, not just xmm registers
    let mut failures = Vec::new();
    for (name, body) in &pown_const_regions {
        let has_branch_or_call =
            body.iter().any(|l| { let t = l.trim_start(); t.starts_with('j') || t.starts_with("call") || t.starts_with("loop") });
        if has_branch_or_call {
            failures.push(format!(
                "{name}: contains a branch/call/loop instruction -- pown_const's const-generic \
                 exponent loop may not have fully constant-folded away (see its own doc comment)"
            ));
        }
    }
    for (name, body) in &regions {
        let has_call = body.iter().any(|l| l.trim_start().starts_with("call"));
        // "contains", not "starts_with": the AVX-encoded mnemonic has a
        // "v" prefix (vcvttss2si/vcvttsd2si), not a bare cvtt... form --
        // starts_with alone silently misses it (caught by deliberately
        // reintroducing the historical saturating-cast pattern into a
        // scratch copy of hypot and confirming this check failed to catch
        // it before this fix, see IDEAS.md).
        let has_saturating_cast = body.iter().any(|l| l.contains("cvttsd2si") || l.contains("cvttss2si"));
        // IDEAS.md idea #75's own explicit ask ("no vsqrtss/vdivss") isn't
        // fully covered by has_packed_arith below, which only confirms *at
        // least one* packed op is present -- a scalar sqrt/div could still
        // hide alongside otherwise-packed code (a partial de-vectorization
        // of just one sub-computation) without tripping that check at all.
        // Checked directly: currently zero regions have this (verified
        // 2026-07-10), but nothing was actually asserting it.
        let has_scalar_sqrt_or_div = body.iter().any(|l| {
            let t = l.trim_start();
            t.starts_with("vsqrtss") || t.starts_with("vdivss") || t.starts_with("vsqrtsd") || t.starts_with("vdivsd")
        });
        let has_packed_arith = body.iter().any(|l| {
            let t = l.trim_start();
            let is_arith = t.starts_with("vadd")
                || t.starts_with("vsub")
                || t.starts_with("vmul")
                || t.starts_with("vdiv")
                || t.starts_with("vfmadd")
                || t.starts_with("vfnmadd")
                || t.starts_with("vfmsub");
            // packed forms end in "ps" (single-precision packed); scalar
            // forms end in "ss" -- reject scalar-only, require at least one
            // genuinely packed op somewhere using a wide register.
            is_arith && t.contains("ps\t") && packed_simd.iter().any(|r| t.contains(r))
        });
        if has_call {
            failures.push(format!("{name}: contains a `call` instruction (libm fallback / de-vectorized loop)"));
        }
        if has_saturating_cast {
            failures.push(format!("{name}: contains cvttsd2si/cvttss2si (saturating-cast de-vectorization, see jodiemath-workflow memory)"));
        }
        if has_scalar_sqrt_or_div {
            failures.push(format!("{name}: contains a scalar vsqrtss/vdivss/vsqrtsd/vdivsd (partial de-vectorization -- should be the packed vXXXps/vXXXpd form)"));
        }
        // nop_throughput is the harness's own deliberately-trivial identity
        // baseline (`|x: f32| x`, see mca_target.rs) -- no arithmetic by
        // design, not a real function to check.
        if !has_packed_arith && name != "nop_throughput" {
            failures.push(format!("{name}: no packed (ymm/zmm) arithmetic instruction found -- loop may not have vectorized at all"));
        }
    }

    println!(
        "checked {} *_throughput region(s) and {} pown_const* region(s) in {}",
        regions.len(),
        pown_const_regions.len(),
        asm_path.display()
    );
    if failures.is_empty() {
        println!("ok: no known de-vectorization signatures found in any region");
    } else {
        println!("FAIL: {} issue(s) found:", failures.len());
        for f in &failures {
            println!("  {f}");
        }
        std::process::exit(1);
    }
}
