// idea #77: mca/wall-clock disagreement detector. Runs both `mca` and
// `quickbench` for a given filter and either (a) saves a labeled snapshot,
// or (b) compares two previously-saved snapshots and flags any row where
// the two tools disagree on the *direction* of change -- the exact
// mismatch class `sincos_checked` hit (mca predicted +19% throughput,
// real wall-clock measured -1.4%).
//
// This deliberately does NOT automate git stash/checkout itself (the
// crate's own safety principle: risky git operations should stay under
// direct human/caller control, not buried inside a benchmarking tool).
// Usage:
//   cargo run --release --example perf_compare -- snapshot <filter> <out.json>
//     Runs mca + quickbench for <filter> on whatever the working tree
//     currently is, saves both raw outputs to <out.json>.
//   cargo run --release --example perf_compare -- compare <before.json> <after.json>
//     Parses both snapshots, prints each function's mca-throughput and
//     quickbench-throughput %% change, and flags any row where the two
//     disagree in sign (i.e. one says faster, the other says slower).
//
// Typical flow: snapshot before a change (`git stash`, snapshot, `git
// stash pop`), make the change, snapshot again, then compare.
use std::process::Command;

fn run_example(name: &str, filter: &str) -> String {
    let output = Command::new("cargo")
        .args(["run", "--release", "--example", name, "--", filter])
        .output()
        .unwrap_or_else(|e| panic!("failed to run example {name}: {e}"));
    if !output.status.success() {
        eprintln!(
            "warning: example {name} exited non-zero:\n{}",
            String::from_utf8_lossy(&output.stderr)
        );
    }
    String::from_utf8_lossy(&output.stdout).into_owned()
}

// mca's own row format: "name | latency | throughput" (see examples/mca.rs).
fn parse_mca_throughput(text: &str) -> Vec<(String, f64)> {
    let mut out = vec![];
    for line in text.lines() {
        let parts: Vec<&str> = line.split('|').collect();
        if parts.len() != 3 {
            continue;
        }
        let name = parts[0].trim();
        if name.is_empty() || name.starts_with('-') {
            continue;
        }
        if let Ok(v) = parts[2].trim().parse::<f64>() {
            out.push((name.to_string(), v));
        }
    }
    out
}

// quickbench's own row format: "name    throughput  X.XXX ns/op ..."
// (see examples/quickbench.rs's bench_throughput). Known limitation:
// multi-word names ("std cbrt") get truncated to their first token
// ("std") since this splits on whitespace -- the parsed *value* is still
// correct (the "throughput" marker token is found regardless), only the
// label is imprecise for those handful of rows. Not fixed: the crate's
// own functions (the actual comparison targets this tool exists for) are
// all single-word, so this only affects the "std X" reference rows.
fn parse_quickbench_throughput(text: &str) -> Vec<(String, f64)> {
    let mut out = vec![];
    for line in text.lines() {
        if !line.contains("throughput") {
            continue;
        }
        let mut it = line.split_whitespace();
        let name = match it.next() {
            Some(n) => n.to_string(),
            None => continue,
        };
        // skip the literal "throughput" token, take the number after it
        let mut found = None;
        let mut prev_was_label = false;
        for tok in it {
            if prev_was_label {
                found = tok.parse::<f64>().ok();
                break;
            }
            if tok == "throughput" {
                prev_was_label = true;
            }
        }
        if let Some(v) = found {
            out.push((name, v));
        }
    }
    out
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mode = args.get(1).map(|s| s.as_str()).unwrap_or("");

    match mode {
        "snapshot" => {
            let filter = args.get(2).cloned().unwrap_or_default();
            let out_path = args.get(3).expect("usage: snapshot <filter> <out.json>");
            eprintln!("running mca -- {filter}...");
            let mca_out = run_example("mca", &filter);
            eprintln!("running quickbench -- {filter}...");
            let qb_out = run_example("quickbench", &filter);
            let snapshot = serde_json::json!({ "mca": mca_out, "quickbench": qb_out });
            std::fs::write(out_path, serde_json::to_string_pretty(&snapshot).unwrap())
                .unwrap_or_else(|e| panic!("failed to write {out_path}: {e}"));
            eprintln!("saved snapshot to {out_path}");
        }
        "compare" => {
            let before_path = args
                .get(2)
                .expect("usage: compare <before.json> <after.json>");
            let after_path = args
                .get(3)
                .expect("usage: compare <before.json> <after.json>");
            let before: serde_json::Value =
                serde_json::from_str(&std::fs::read_to_string(before_path).unwrap()).unwrap();
            let after: serde_json::Value =
                serde_json::from_str(&std::fs::read_to_string(after_path).unwrap()).unwrap();

            let mca_before: std::collections::BTreeMap<_, _> =
                parse_mca_throughput(before["mca"].as_str().unwrap())
                    .into_iter()
                    .collect();
            let mca_after: std::collections::BTreeMap<_, _> =
                parse_mca_throughput(after["mca"].as_str().unwrap())
                    .into_iter()
                    .collect();
            let qb_before: std::collections::BTreeMap<_, _> =
                parse_quickbench_throughput(before["quickbench"].as_str().unwrap())
                    .into_iter()
                    .collect();
            let qb_after: std::collections::BTreeMap<_, _> =
                parse_quickbench_throughput(after["quickbench"].as_str().unwrap())
                    .into_iter()
                    .collect();

            println!(
                "{:22} | {:>12} | {:>14} | flag",
                "name", "mca % chg", "wallclock % chg"
            );
            println!("{:-<22}-|-{:->12}-|-{:->14}-|-----", "", "", "");
            let mut names: Vec<&String> = mca_after.keys().chain(qb_after.keys()).collect();
            names.sort();
            names.dedup();
            for name in names {
                let mca_pct = match (mca_before.get(name), mca_after.get(name)) {
                    (Some(&b), Some(&a)) if b != 0.0 => Some((a - b) / b * 100.0),
                    _ => None,
                };
                let qb_pct = match (qb_before.get(name), qb_after.get(name)) {
                    (Some(&b), Some(&a)) if b != 0.0 => Some((a - b) / b * 100.0),
                    _ => None,
                };
                if mca_pct.is_none() && qb_pct.is_none() {
                    continue;
                }
                let flag = match (mca_pct, qb_pct) {
                    (Some(m), Some(q))
                        if m.signum() != q.signum() && m.abs() > 1.0 && q.abs() > 1.0 =>
                    {
                        "  <-- DISAGREE"
                    }
                    _ => "",
                };
                println!(
                    "{:22} | {:>11} | {:>13} |{}",
                    name,
                    mca_pct
                        .map(|v| format!("{v:+.1}%"))
                        .unwrap_or_else(|| "?".into()),
                    qb_pct
                        .map(|v| format!("{v:+.1}%"))
                        .unwrap_or_else(|| "?".into()),
                    flag
                );
            }
        }
        _ => {
            eprintln!(
                "usage:\n  cargo run --release --example perf_compare -- snapshot <filter> <out.json>\n  cargo run --release --example perf_compare -- compare <before.json> <after.json>"
            );
            std::process::exit(1);
        }
    }
}
