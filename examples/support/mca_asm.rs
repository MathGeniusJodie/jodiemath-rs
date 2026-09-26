// Shared by mca.rs and codegen_check.rs: compile examples/mca_target.rs to
// assembly and return (path, text) of the fresh .s file.
fn emit_mca_target_asm() -> (std::path::PathBuf, String) {
    use std::path::{Path, PathBuf};
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
    let status = std::process::Command::new("cargo")
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

    // Older cargo writes `--emit=asm` output to target/release/examples/
    // mca_target-<hash>.s, newer cargo to target/release/build/<pkg>/<hash>/
    // out/mca_target.s; take the newest of either.
    fn asm_files(dir: &Path, out: &mut Vec<PathBuf>) {
        let Ok(entries) = std::fs::read_dir(dir) else {
            return;
        };
        for path in entries.filter_map(|e| e.ok()).map(|e| e.path()) {
            if path.is_dir() {
                asm_files(&path, out);
            } else if path
                .file_name()
                .and_then(|n| n.to_str())
                .is_some_and(|n| n.starts_with("mca_target") && n.ends_with(".s"))
            {
                out.push(path);
            }
        }
    }
    let mut candidates = Vec::new();
    asm_files(&PathBuf::from(manifest_dir).join("target/release"), &mut candidates);
    let asm_path = candidates
        .into_iter()
        .filter_map(|p| std::fs::metadata(&p).and_then(|m| m.modified()).ok().map(|t| (t, p)))
        .max_by_key(|(t, _)| *t)
        .map(|(_, p)| p)
        .expect("no mca_target*.s found after `cargo rustc --emit=asm` -- did the example build?");
    let text = std::fs::read_to_string(&asm_path).expect("couldn't read the emitted asm");
    (asm_path, text)
}
