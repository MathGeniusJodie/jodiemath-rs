#!/usr/bin/env python3
"""Extract named LLVM-MCA regions from a .s file and run llvm-mca on each in
isolation.  Same numbers as running llvm-mca over the whole file (verified:
asin_throughput 1840 cycles both ways), in seconds instead of ~20 minutes.

usage: tools/mca_region.py <file.s> <region> [<region> ...]
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

path = sys.argv[1]
wanted = sys.argv[2:]

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
    div = CHAIN_LEN if w.endswith("_latency") else ARR_LEN
    per = s["TotalCycles"] / (s["Iterations"] * div)
    print(f"{w:<28} {s['Instructions']//ITERS:>7} {s['TotaluOps']//ITERS:>7} "
          f"{s['BlockRThroughput']:>8.2f} {s['TotalCycles']:>8} {per:>9.3f}")
