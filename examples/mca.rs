// Theoretical latency/throughput via llvm-mca -- no wall-clock timing, so
// no thermal-throttling noise (the i5-1145G7 in this readme throttles up to
// ~2.5x mid-session; see quickbench's nop-baseline caveat). This drives a
// two-step pipeline:
//   1. compile examples/mca_target.rs to assembly (`cargo rustc --emit=asm`)
//   2. run `llvm-mca -mcpu=native --json` over it and parse the per-region
//      summaries (see mca_target.rs for why each region is built the way it
//      is, and why the numbers are cycles, not nanoseconds: llvm-mca models
//      an LLVM scheduler description, not a specific clock frequency).
//
// Requires `llvm-mca` on PATH (Debian/Arch: part of the `llvm` package).
use serde_json::Value;
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::process::Command;

include!("support/mca_common.rs");
const MCA_ITERATIONS: u32 = 100;

fn main() {
    let manifest_dir = env!("CARGO_MANIFEST_DIR");

    // Cargo's fingerprint cache doesn't account for the `--emit=asm` passed
    // below as a raw rustc arg, so if mca_target.rs hasn't changed since the
    // last (non-asm) build, `cargo rustc` treats this as a no-op and never
    // regenerates the .s file. Bumping the mtime forces a real rebuild.
    let target_src = PathBuf::from(manifest_dir).join("examples/mca_target.rs");
    std::fs::File::open(&target_src)
        .and_then(|f| f.set_modified(std::time::SystemTime::now()))
        .expect("couldn't touch examples/mca_target.rs to force a rebuild");

    eprintln!("compiling examples/mca_target.rs to assembly...");
    let status = Command::new("cargo")
        .current_dir(manifest_dir)
        .args([
            "rustc",
            "--release",
            "--example",
            "mca_target",
            "--",
            "--emit=asm",
            "-C",
            "debuginfo=0",
        ])
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
            p.file_name()
                .and_then(|n| n.to_str())
                .is_some_and(|n| n.starts_with("mca_target-") && n.ends_with(".s"))
        })
        .filter_map(|p| std::fs::metadata(&p).and_then(|m| m.modified()).ok().map(|t| (t, p)))
        .max_by_key(|(t, _)| *t)
        .map(|(_, p)| p)
        .expect("no mca_target-*.s found after `cargo rustc --emit=asm` -- did the example build?");

    eprintln!("running llvm-mca on {}...", asm_path.display());
    let output = Command::new("llvm-mca")
        .arg("-mcpu=native")
        .arg(format!("--iterations={MCA_ITERATIONS}"))
        .arg("--json")
        .arg(&asm_path)
        .output()
        .expect("failed to run llvm-mca -- is it installed and on PATH? (Debian/Arch: `llvm` package)");
    if !output.status.success() {
        eprintln!("llvm-mca failed:\n{}", String::from_utf8_lossy(&output.stderr));
        std::process::exit(1);
    }

    let json: Value =
        serde_json::from_slice(&output.stdout).expect("llvm-mca did not produce valid JSON");
    let regions = json["CodeRegions"]
        .as_array()
        .expect("no CodeRegions in llvm-mca output");

    let mut latency: BTreeMap<String, f64> = BTreeMap::new();
    let mut throughput: BTreeMap<String, f64> = BTreeMap::new();
    for region in regions {
        let name = region["Name"].as_str().unwrap_or("");
        let summary = &region["SummaryView"];
        let total_cycles = summary["TotalCycles"].as_f64().expect("missing TotalCycles");
        let iterations = summary["Iterations"].as_f64().expect("missing Iterations");
        if let Some(key) = name.strip_suffix("_latency") {
            latency.insert(key.to_string(), total_cycles / (iterations * CHAIN_LEN as f64));
        } else if let Some(key) = name.strip_suffix("_throughput") {
            throughput.insert(key.to_string(), total_cycles / (iterations * ARR_LEN as f64));
        }
    }

    let order = [
        "nop",
        "cbrt",
        "cbrt_wrapped",
        "cbrt_accurate",
        "cbrt_throughput_fn",
        "cbrt_fast",
        "exp2",
        "exp2_checked",
        "exp10",
        "exp10_checked",
        "log2",
        "log2_unchecked",
        "sin",
        "sin_checked",
        "cos",
        "cos_checked",
        "sinpi",
        "cospi",
        "sind",
        "cosd",
        "ln",
        "ln_unchecked",
        "log10",
        "log10_unchecked",
        "log1p",
        "exp",
        "expm1",
        "sinh",
        "cosh",
        "sinh_throughput_fn",
        "cosh_throughput_fn",
        "tanh",
        "sigmoid",
        "asinh",
        "acosh",
        "atanh",
        "asin",
        "acos",
        "atan",
        "atan_latency",
        "atan2",
        "tan",
        "erf",
        "erfc",
        "hypot",
        "hypot_checked",
        "rsqrt",
        "powf",
        "pown",
        "powf_checked",
        "remainder",
        "remainder_checked",
    ];

    println!();
    println!("theoretical cost from llvm-mca (-mcpu=native, {MCA_ITERATIONS} iterations)");
    println!("latency: branchless *_normal core, 64-deep serial dependency chain, cycles/call");
    println!("throughput: real public function, a 16-wide (two AVX2 vectors) auto-vectorized block, cycles/element");
    println!();
    println!("{:19} | latency (cyc) | throughput (cyc)", "");
    println!("{:-<19}-|-{:->14}-|-{:->17}", "", "", "");
    for key in order {
        println!(
            "{:19} | {:>14} | {:>17}",
            key,
            latency
                .get(key)
                .map(|v| format!("{v:.2}"))
                .unwrap_or_else(|| "?".into()),
            throughput
                .get(key)
                .map(|v| format!("{v:.3}"))
                .unwrap_or_else(|| "?".into()),
        );
    }
    println!();
    println!(
        "nop is the harness's own overhead (chain mixing / array bookkeeping) -- \
         compare other rows against it, same as quickbench's readme caveat."
    );
}
