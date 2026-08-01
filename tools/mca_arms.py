#!/usr/bin/env python3
"""Isolate each arm of a branch-shaped llvm-mca *latency* region.

llvm-mca has no branch predictor: it simulates a region as one straight-line
instruction stream.  For the diamond `jcc L1 / <A> / jmp L2 / L1: <B> / L2:`
it therefore executes A *and* B back to back, so the number it reports is
neither arm's.  Two distinct failure modes follow, and this crate has one
measured example of each:

  - **concatenation (overstates).**  If B reads a register A clobbered, the
    two arms fuse into one artificially long dependency chain.  `acosh`
    publishes 89.08 against arms of 46.08 and 77.02 -- above *both*.  24 of
    the 39 branch-shaped latency regions in `mca_target.rs` are like this,
    the exp/expm1 family worst (`exp2m1` 80.00 against arms 37.00/48.00).
  - **last-writer wins (understates).**  If both arms write the same
    register and the *cheap* one is laid out last, mca times the cheap one.
    `asinh` publishes 69.41; its real in-domain chain is 88.02.

Throughput regions are not affected -- they are vectorized and LLVM
if-converts them to masked selects, so there is no branch to mis-simulate.

Usage:
    tools/mca_arms.py <file.s> <region_latency> [...]

prints the as-published figure next to each arm in isolation.  A published
value outside the interval its own two arms span is proof mca is simulating
neither.  Which arm is the honest one is a property of the *function*, not of
the asm -- read the guard.  For a region with several branches, pass a
comma-separated pattern instead of a single region name choice by editing the
call to `strip` (e.g. `strip(v, "fall,target")` for `asinh`, whose first
branch takes the fall-through in domain and whose second takes the target).
"""
import re, subprocess, tempfile, os, sys

def load(path):
    out={}; cur=None; buf=[]
    for ln in open(path):
        m=re.search(r"# LLVM-MCA-BEGIN (\S+)", ln)
        if m: cur=m.group(1); buf=[]; continue
        if cur is not None:
            if "# LLVM-MCA-END" in ln: out[cur]=buf; cur=None; continue
            buf.append(ln.rstrip())
    return out

JCC = re.compile(r'^\s+j(?!mp\b)\w+\s+(\.LBB\S+)')
JMP = re.compile(r'^\s+jmp\s+(\.LBB\S+)')
LBL = re.compile(r'^(\.LBB\S+):')

def strip(lines, keep):
    """keep is 'fall'/'target', or a comma list cycled over branches in order."""
    pat = keep.split(',') if ',' in keep else [keep]
    out=[]; i=0; n=len(lines); b=0
    while i < n:
        m = JCC.match(lines[i])
        if not m:
            out.append(lines[i]); i+=1; continue
        l1 = m.group(1)
        # locate jmp L2, label L1, label L2 after i
        j = next((k for k in range(i+1,n) if JMP.match(lines[k])), None)
        k1 = next((k for k in range(i+1,n) if LBL.match(lines[k]) and LBL.match(lines[k]).group(1)==l1), None)
        if j is None or k1 is None or not (i < j < k1):
            out.append(lines[i]); i+=1; continue
        l2 = JMP.match(lines[j]).group(1)
        k2 = next((k for k in range(k1+1,n) if LBL.match(lines[k]) and LBL.match(lines[k]).group(1)==l2), None)
        if k2 is None:
            out.append(lines[i]); i+=1; continue
        A = lines[i+1:j]      # fall-through arm
        B = lines[k1+1:k2]    # jump-target arm
        out.extend(A if pat[b % len(pat)]=='fall' else B); b+=1
        i = k2+1              # drop the L2 label too
    return out

def run(lines, iters=100):
    with tempfile.NamedTemporaryFile('w', suffix='.s', delete=False) as f:
        f.write(".text\n" + "\n".join(lines) + "\n"); t=f.name
    r = subprocess.run(["llvm-mca","-mcpu=native",f"--iterations={iters}","--json",t],
                       capture_output=True, text=True)
    os.unlink(t)
    if r.returncode != 0: return None
    import json
    d = json.loads(r.stdout)["CodeRegions"][0]["SummaryView"]
    return d["TotalCycles"]/(d["Iterations"]*64.0)

if __name__ == "__main__":
    regs = load(sys.argv[1])
    print(f"{'region':<32} {'as-published':>13} {'fall-arm':>10} {'target-arm':>11}")
    for name in sys.argv[2:]:
        v = regs[name]
        base = run(v)
        fa = run(strip(v,'fall')); ta = run(strip(v,'target'))
        f = lambda x: f"{x:.2f}" if x else "  --  "
        print(f"{name:<32} {f(base):>13} {f(fa):>10} {f(ta):>11}")
