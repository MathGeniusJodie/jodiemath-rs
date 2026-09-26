#!/usr/bin/env python3
"""Extract named LLVM-MCA regions from a .s file and run llvm-mca on each in
isolation.  Same numbers as running llvm-mca over the whole file (verified:
asin_throughput 1840 cycles both ways), in seconds instead of ~20 minutes.

usage: tools/mca_region.py [<file.s>] <region> [<region> ...]
       (without a file, the newest mca_target*.s under target/)
prints: region, instrs, uOps, Block RThroughput, TotalCycles, cyc/unit

The columns are deliberately the whole escalation ladder in one call, because
mca's throughput number alone is not trustworthy on this crate -- it has been
wrong in both directions (overstating sin/cos/tan's cost by 54% where instrs,
uOps and RThroughput agreed on ~+20%, and understating norm_cdf's). Read
instrs + uOps + BlockRT together before believing a delta.

Find the newest .s with:
  ls -t target/release/examples/mca_target-*.s | head -1
and note that `cargo rustc --emit=asm` must have run since your last edit --
examples/mca.rs bumps mca_target.rs's mtime for exactly that reason.
"""
import json
import re
import subprocess
import sys
import tempfile
import os

CHAIN_LEN = 64   # mca.rs: latency regions divide by iterations*CHAIN_LEN
ARR_LEN = 16     # mca.rs: throughput regions divide by iterations*ARR_LEN
ITERS = 100

def newest_asm():
    """Newest mca_target .s: cargo puts --emit=asm output in
    target/release/examples/ on older toolchains and under
    target/release/build/<pkg>/<hash>/out/ on newer ones."""
    found = []
    for root, _, files in os.walk("target"):
        found += [os.path.join(root, f) for f in files
                  if f.startswith("mca_target") and f.endswith(".s")]
    if not found:
        sys.exit("no mca_target*.s under target/ -- run "
                 "`cargo rustc --release --example mca_target -- --emit=asm` first")
    return max(found, key=os.path.getmtime)


if sys.argv[1].endswith(".s") and os.path.isfile(sys.argv[1]):
    path, wanted = sys.argv[1], sys.argv[2:]
else:
    path, wanted = newest_asm(), sys.argv[1:]

regions, cur, buf = {}, None, []
for ln in open(path).read().split("\n"):
    m = re.search(r"# LLVM-MCA-BEGIN (\S+)", ln)
    if m:
        cur, buf = m.group(1), [ln]
        continue
    if cur is not None:
        buf.append(ln)
        if "# LLVM-MCA-END" in ln:
            regions[cur] = buf
            cur = None

missing = [w for w in wanted if w not in regions]
if missing:
    print("regions not found:", missing, file=sys.stderr)
    sys.exit(1)


def loop_step(body, default):
    """Elements one pass through a *throughput* region actually covers.

    A throughput region is meant to be the fully-unrolled body of a
    16-element loop, so dividing TotalCycles by ARR_LEN gives cycles per
    element. When the region is big enough that LLVM keeps the loop
    instead of unrolling it, llvm-mca simulates *one iteration* -- and the
    /16 is then wrong by exactly the unroll factor. Detect it: a backward
    branch to a label defined in the region, plus the induction step of
    the `addq $N, %r..` / `cmpq $16, %r..` pair that drives it.

    Real cases at the time of writing: tan_wide (N=4, 4x low) and clog_re
    (N=8, 2x low). Silent, and in the flattering direction.
    """
    pos = {}
    for i, l in enumerate(body):
        m = re.match(r"^(\.LBB\d+_\d+):", l.strip())
        if m:
            pos[m.group(1)] = i
    back = False
    for i, l in enumerate(body):
        m = re.match(r"^\s+j[a-z]+\s+(\.LBB\d+_\d+)", l)
        if m and m.group(1) in pos and pos[m.group(1)] < i:
            back = True
    if not back:
        return default, False
    for l in body:
        m = re.match(r"^\s+addq\s+\$(\d+), %r", l)
        if m and 0 < int(m.group(1)) <= default:
            return int(m.group(1)), True
    return default, True

print(f"{'region':<28} {'instrs':>7} {'uOps':>7} {'BlockRT':>8} {'cycles':>8} {'cyc/unit':>9}")
for w in wanted:
    with tempfile.NamedTemporaryFile("w", suffix=".s", delete=False) as f:
        f.write(".text\n" + "\n".join(regions[w]) + "\n")
        tmp = f.name
    r = subprocess.run(
        ["llvm-mca", "-mcpu=native", f"--iterations={ITERS}", "--json", tmp],
        capture_output=True, text=True)
    os.unlink(tmp)
    if r.returncode != 0:
        print(w, "FAILED:", r.stderr.strip()[:300], file=sys.stderr)
        continue
    s = json.loads(r.stdout)["CodeRegions"][0]["SummaryView"]
    if w.endswith("_latency"):
        div, looped = CHAIN_LEN, False
    else:
        div, looped = loop_step(regions[w], ARR_LEN)
    per = s["TotalCycles"] / (s["Iterations"] * div)
    mark = f"  (loop kept, /{div} not /{ARR_LEN})" if looped else ""
    print(f"{w:<28} {s['Instructions']//ITERS:>7} {s['TotaluOps']//ITERS:>7} "
          f"{s['BlockRThroughput']:>8.2f} {s['TotalCycles']:>8} {per:>9.3f}{mark}")
