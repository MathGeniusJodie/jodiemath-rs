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
use std::process::Command;

include!("support/mca_common.rs");
include!("support/mca_asm.rs");
const MCA_ITERATIONS: u32 = 100;

/// Elements one simulated pass through each `*_throughput` region really
/// covers, keyed by region name; absent means the region is straight-line
/// and `ARR_LEN` is right.
///
/// A region counts as "loop kept" when it contains a backward branch to
/// one of its own labels. The element step is then the `addq $N, %r..`
/// that drives the `cmpq $ARR_LEN` beside it -- LLVM's own induction
/// variable, so it is the unroll factor by construction rather than a
/// guess. Text-scanning the asm rather than asking llvm-mca, which
/// reports cycles for the block it was given and has no idea the block
/// was supposed to be 16 elements wide.
fn loop_steps(asm: &str) -> BTreeMap<String, usize> {
    let mut out = BTreeMap::new();
    let mut name: Option<String> = None;
    let mut body: Vec<&str> = Vec::new();
    for line in asm.lines() {
        if let Some(rest) = line.split("# LLVM-MCA-BEGIN ").nth(1) {
            name = Some(rest.split_whitespace().next().unwrap_or("").to_string());
            body.clear();
            continue;
        }
        let Some(n) = name.clone() else { continue };
        if line.contains("# LLVM-MCA-END") {
            if n.ends_with("_throughput") {
                let labels: Vec<&str> = body
                    .iter()
                    .filter_map(|l| l.trim().strip_suffix(':'))
                    .filter(|l| l.starts_with(".LBB"))
                    .collect();
                let mut seen: Vec<&str> = Vec::new();
                let mut backward = false;
                for l in &body {
                    let t = l.trim();
                    if let Some(lab) = t.strip_suffix(':') {
                        seen.push(lab);
                    }
                    if t.starts_with('j') {
                        if let Some(tgt) = t.split_whitespace().nth(1) {
                            if labels.contains(&tgt) && seen.contains(&tgt) {
                                backward = true;
                            }
                        }
                    }
                }
                if backward {
                    for l in &body {
                        let t = l.trim();
                        if let Some(rest) = t.strip_prefix("addq\t$") {
                            if let Some((num, reg)) = rest.split_once(", %r") {
                                if !reg.is_empty() {
                                    if let Ok(v) = num.parse::<usize>() {
                                        if v > 0 && v <= ARR_LEN {
                                            out.insert(n.clone(), v);
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            name = None;
            continue;
        }
        body.push(line);
    }
    out
}

fn main() {
    // Both this process and the `cargo rustc`/`llvm-mca` children it spawns
    // below inherit this niceness (nice values survive fork/exec) -- llvm-mca
    // itself is the real CPU/memory hog here, easily pegging a full core for
    // over a minute on this crate's larger assembly files.
    nice_self();
    // Unlike quickbench, this used to silently ignore any argument (no
    // env::args() call at all) -- `./mca powf` looked like it filtered but
    // actually always printed the full ~52-row table. Real filter support:
    // still recompiles/reruns llvm-mca over everything (region names aren't
    // known until the asm exists), just narrows the printed rows.
    let args: Vec<String> = std::env::args().collect();
    let filter = args.get(1).map(|s| s.as_str()).unwrap_or("");

    let (asm_path, asm_text) = emit_mca_target_asm();

    // rustc's LLVM can be newer than the llvm-mca on PATH; `.prefalign`
    // (function-level alignment, LLVM 23+) is the directive that older
    // assemblers reject, and it never sits inside a measured region.
    let mca_input: String = asm_text
        .lines()
        .filter(|l| !l.trim_start().starts_with(".prefalign"))
        .flat_map(|l| [l, "\n"])
        .collect();

    eprintln!("running llvm-mca on {}...", asm_path.display());
    let mut child = Command::new("llvm-mca")
        .arg("-mcpu=native")
        .arg(format!("--iterations={MCA_ITERATIONS}"))
        .arg("--json")
        .arg("-")
        .stdin(std::process::Stdio::piped())
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .expect(
            "failed to run llvm-mca -- is it installed and on PATH? (Debian/Arch: `llvm` package)",
        );
    let mut stdin = child.stdin.take().expect("llvm-mca stdin");
    let writer = std::thread::spawn(move || {
        use std::io::Write;
        stdin.write_all(mca_input.as_bytes())
    });
    let output = child.wait_with_output().expect("llvm-mca did not run");
    writer
        .join()
        .expect("llvm-mca stdin writer panicked")
        .expect("couldn't write asm to llvm-mca");
    if !output.status.success() {
        eprintln!(
            "llvm-mca failed:\n{}",
            String::from_utf8_lossy(&output.stderr)
        );
        std::process::exit(1);
    }

    let json: Value =
        serde_json::from_slice(&output.stdout).expect("llvm-mca did not produce valid JSON");
    let regions = json["CodeRegions"]
        .as_array()
        .expect("no CodeRegions in llvm-mca output");

    // A throughput region is *meant* to be the fully-unrolled body of the
    // 16-element loop in `throughput_fn!`, which is what makes
    // `TotalCycles / (iterations * ARR_LEN)` cycles-per-element. When the
    // region gets big enough that LLVM keeps the loop instead of unrolling
    // it, llvm-mca simulates one *iteration*, and dividing by 16 is then
    // wrong by exactly the unroll factor -- silently, and in the
    // flattering direction. `loop_steps` reads the real step out of the
    // asm; see its own comment.
    let steps = loop_steps(&asm_text);
    let mut latency: BTreeMap<String, f64> = BTreeMap::new();
    let mut throughput: BTreeMap<String, f64> = BTreeMap::new();
    let mut looped: Vec<String> = Vec::new();
    let mut order: Vec<String> = Vec::new();
    for region in regions {
        let name = region["Name"].as_str().unwrap_or("");
        let summary = &region["SummaryView"];
        let total_cycles = summary["TotalCycles"]
            .as_f64()
            .expect("missing TotalCycles");
        let iterations = summary["Iterations"].as_f64().expect("missing Iterations");
        let key = name
            .strip_suffix("_latency")
            .or_else(|| name.strip_suffix("_throughput"));
        if let Some(key) = key.filter(|k| !order.iter().any(|o| o == k)) {
            order.push(key.to_string());
        }
        if let Some(key) = name.strip_suffix("_latency") {
            latency.insert(
                key.to_string(),
                total_cycles / (iterations * CHAIN_LEN as f64),
            );
        } else if let Some(key) = name.strip_suffix("_throughput") {
            let n = steps.get(name).copied().unwrap_or(ARR_LEN);
            if n != ARR_LEN {
                looped.push(format!("{key} (/{n})"));
            }
            throughput.insert(key.to_string(), total_cycles / (iterations * n as f64));
        }
    }
    if !looped.is_empty() {
        eprintln!(
            "note: loop kept (divided by the real step, not {ARR_LEN}): {}",
            looped.join(", ")
        );
    }


    // llvm-mca reports regions in asm layout order, which is meaningless;
    // list them in the order mca_target.rs declares them instead.
    let target_text = std::fs::read_to_string(concat!(env!("CARGO_MANIFEST_DIR"), "/examples/mca_target.rs")).expect("couldn't read mca_target.rs");
    let declared_at = |key: &str| {
        ["latency", "throughput"]
            .iter()
            .filter_map(|kind| target_text.find(&format!("\"{key}_{kind}\"")))
            .min()
            .unwrap_or(usize::MAX)
    };
    order.sort_by_key(|key| declared_at(key));

    println!();
    println!("theoretical cost from llvm-mca (-mcpu=native, {MCA_ITERATIONS} iterations)");
    println!("latency: branchless *_normal core, 64-deep serial dependency chain, cycles/call");
    println!("throughput: real public function, a 16-element auto-vectorized block, cycles/element");
    println!();
    println!("{:19} | latency (cyc) | throughput (cyc)", "");
    println!("{:-<19}-|-{:->14}-|-{:->17}", "", "", "");
    for key in &order {
        if !filter.is_empty() && !key.contains(filter) {
            continue;
        }
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
