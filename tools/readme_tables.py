#!/usr/bin/env python3
"""Rewrite readme.md's generated tables from harness output.

usage: tools/readme_tables.py <stats-dir> <readme.md>

<stats-dir> holds the raw output of tools/readme_stats.sh: accuracy.txt
(`accuracy thorough`), powfsearch.txt, clogsearch.txt, mca.txt and
quickbench.<round>.txt (min taken per row), plus machine.txt. Each table replaces the text between its
`<!-- BEGIN generated:NAME -->` / `<!-- END generated:NAME -->` markers.
"""
import re
import sys
from pathlib import Path

VARIANT_SUFFIXES = (
    "_unchecked", "_checked", "_narrow", "_throughput", "_latency",
    "_pos", "_accurate", "_bounded", "_fast", "_wide",
)
ALIASES = {"log_2": "log2"}
BAND = re.compile(r"^\S+ [\[|]|\(naive\)")


def root(name):
    changed = True
    while changed:
        changed = False
        for s in VARIANT_SUFFIXES:
            if name.endswith(s):
                name, changed = name[: -len(s)], True
    return ALIASES.get(name, name)


def std_for(name, std):
    """The std row a jodie row compares against, or None."""
    name = name.replace(" (+)", "")
    head, _, qual = name.partition(" ")
    r = root(head)
    candidates = [name, ALIASES.get(head, head)]
    if qual:
        candidates.append(f"{r} {qual}")
    if head.endswith("_checked"):
        candidates.append(f"{r} (everywhere)")
    candidates += [r, f"{r} (in-domain)"]
    return next((std[c] for c in candidates if c in std), None)


def table(header, rows):
    widths = [max(len(str(r[i])) for r in [header] + rows) for i in range(len(header))]
    rows = [[str(c).replace("|", "\\|") for c in r] for r in rows]
    header = [str(c).replace("|", "\\|") for c in header]
    fmt = lambda r: "| " + " | ".join(
        str(c).ljust(w) if i == 0 else str(c).rjust(w) for i, (c, w) in enumerate(zip(r, widths))
    ) + " |"
    sep = "|" + "|".join(":" + "-" * (w + 1) if i == 0 else "-" * (w + 1) + ":" for i, w in enumerate(widths)) + "|"
    return "\n".join([fmt(header), sep] + [fmt(r) for r in rows])


def samples(n):
    return "2^32" if n == 2**32 else f"{n:.2e}"


ULP_ROW = re.compile(
    r"^(?P<name>.+?)\s+avg ulp\s+(?P<avg>[\d.]+)\s+max ulp\s+(?P<max>\d+)\s+worst (?P<worst>.*?)"
    r"\s\(\s*(?P<n>\d+) samples"
)
ELAPSED = re.compile(r",?\s*\(?\s*[\d.]+s elapsed\)?")


def accuracy(stats):
    lines = (stats / "accuracy.txt").read_text().splitlines()
    main, bands, other = [], [], []
    for line in lines[1:]:
        m = ULP_ROW.match(line)
        if not m:
            if line.strip() and not line.startswith("total:"):
                line = re.sub(r"\(\s+", "(", ELAPSED.sub("", line).rstrip())
                other.append(line + ")" * (line.count("(") - line.count(")")))
            continue
        row = (m["name"], float(m["avg"]), int(m["max"]), int(m["n"]), m["worst"].strip().removeprefix("x "))
        (bands if BAND.search(row[0]) else main).append(row)
    std = {r[0][4:]: r for r in main if r[0].startswith("std ")}
    ours = [r for r in main if not r[0].startswith("std ")]

    out = []
    for name, avg, mx, n, worst in ours:
        s = std_for(name, std)
        out.append([
            f"`{name}`", f"{avg:.4f}", mx,
            f"{s[1]:.4f}" if s else "-", s[2] if s else "-",
            samples(n), worst,
        ])
    body = table(["function", "avg ulp", "max ulp", "std avg", "std max", "inputs", "worst x"], out)
    band_body = table(
        ["row", "avg ulp", "max ulp", "inputs"],
        [[f"`{r[0]}`", f"{r[1]:.4f}", r[2], samples(r[3])] for r in bands],
    )

    powf = [l for l in (stats / "powfsearch.txt").read_text().splitlines() if l.startswith("worst over")]
    clog = [l.strip() for l in (stats / "clogsearch.txt").read_text().splitlines()
            if l.strip().startswith(("worst", "PASS", "FAIL"))]
    other += ["", "powfsearch: " + " ".join(powf)] + ["clogsearch: " + l for l in clog]
    return {
        "accuracy": body,
        "accuracy-bands": band_body,
        "accuracy-other": "```\n" + "\n".join(other).strip() + "\n```",
    }


BENCH_ROW = re.compile(r"^(?P<name>.+?)\s+(?P<kind>latency|throughput)\s+(?P<ns>[\d.]+) ns/op")


def bench(stats):
    rows = {}
    for run in sorted(stats.glob("quickbench.*.txt")):
        for line in run.read_text().splitlines():
            m = BENCH_ROW.match(line)
            if m:
                kinds = rows.setdefault(m["name"], {})
                kinds[m["kind"]] = min(kinds.get(m["kind"], float("inf")), float(m["ns"]))
    std = {k[4:]: v for k, v in rows.items() if k.startswith("std ")}

    def paired(name):
        head = name.split(" ")[0]
        for c in (name, ALIASES.get(head, head), root(head)):
            if c in std:
                return std[c]
        return None

    out = {}
    for kind, digits in (("latency", 2), ("throughput", 3)):
        body = []
        for name, v in rows.items():
            if name.startswith("std "):
                continue
            s = paired(name)
            ours, theirs = v[kind], s[kind] if s else None
            body.append([
                f"`{name}`", f"{ours:.{digits}f}",
                f"{theirs:.{digits}f}" if theirs is not None else "-",
                f"{theirs / ours:.1f}x" if theirs is not None and ours > 0 else "-",
            ])
        out[f"bench-{kind}"] = table(["function", "jodie ns", "std ns", "speedup"], body)
    return out


MCA_ROW = re.compile(r"^(?P<name>\S+)\s*\|\s*(?P<lat>[\d.?]+)\s*\|\s*(?P<thr>[\d.?]+)\s*$")


def mca(stats):
    body = [
        [f"`{m['name']}`", m["lat"], m["thr"]]
        for m in map(MCA_ROW.match, (stats / "mca.txt").read_text().splitlines())
        if m
    ]
    return {"mca": table(["function", "latency (cyc)", "throughput (cyc/elem)"], body)}


def machine(stats):
    return {"machine": "```\n" + (stats / "machine.txt").read_text().strip() + "\n```"}


def main():
    stats, readme = Path(sys.argv[1]), Path(sys.argv[2])
    sections = {**machine(stats), **accuracy(stats), **bench(stats), **mca(stats)}
    text = readme.read_text()
    for name, body in sections.items():
        pat = re.compile(
            rf"(<!-- BEGIN generated:{re.escape(name)} -->\n).*?(<!-- END generated:{re.escape(name)} -->)",
            re.S,
        )
        if not pat.search(text):
            sys.exit(f"readme has no generated:{name} markers")
        text = pat.sub(lambda m: m[1] + body + "\n" + m[2], text)
    readme.write_text(text)


if __name__ == "__main__":
    main()
