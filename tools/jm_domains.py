#!/usr/bin/env python3
"""Derive ownership domains for src/lib.rs from its own call graph.

Every function in this crate lives in one flat file and most of them share
kernels: `sin`, `cospi`, `tan2pi` and `sind` all bottom out in the same
`pi_reduce_and_poly!` macro and `sinf_poly`. So "two Claudes must not work on
the same function" is the wrong unit -- editing `sinf_poly` for `sin` silently
moves `cospi`'s error too. The unit that actually has to be exclusive is the
*cluster of public functions that share private code*.

This computes those clusters:

  1. Parse top-level `fn` and `macro_rules!` definitions with line ranges.
     Comments are stripped before edges are read -- doc comments in this crate
     name other functions constantly ("see `sin`"), and taking those as calls
     collapses everything into one blob.
  2. Edge u -> v when u's body mentions v (call or macro invocation).
  3. Units reached by a large share of the API are CORE infrastructure
     (`fma`, `mulsign`, `two_sum`, ...). Those are not owned by anyone;
     touching one needs the global lock, because it moves every function.
  4. Domain = connected component of public functions joined by any shared
     non-core unit.

Output is JSON on stdout: domains, their members, line ranges, and core units.
"""

import json
import re
import sys
from pathlib import Path

# A unit reached by more than this fraction of the public API is treated as
# core infrastructure rather than as evidence that two functions are related.
CORE_REACH_FRACTION = 0.03
CORE_REACH_MIN = 6

DEF_RE = re.compile(
    r"^(?P<vis>pub(?:\([a-z]+\))?\s+)?(?:const\s+|unsafe\s+)*fn\s+(?P<name>\w+)"
)
MACRO_RE = re.compile(r"^macro_rules!\s+(?P<name>\w+)")
# Top-level consts, chiefly the per-function polynomial coefficient arrays.
# These are where most of this crate's real edits land, so they have to be
# attributable to an owner. Ownership is by *reference* rather than position:
# a coefficient array belongs to whichever domain names it. `const fn` is
# excluded by requiring a SCREAMING_CASE identifier.
# The name must start with a letter: this crate uses `const _: () = assert!(..)`
# for compile-time checks, and every one of those is called `_`. They collapse
# into a single unit, and a bare `_` is a token in half the bodies in the file,
# so admitting them welds unrelated domains together.
CONST_RE = re.compile(
    r"^(?P<vis>pub(?:\([a-z]+\))?\s+)?(?:const|static)\s+(?P<name>[A-Z][A-Z_0-9]*)\s*:"
)
MOD_RE = re.compile(r"^(?:pub\s+)?mod\s+(?P<name>\w+)\s*\{")
ATTR_RE = re.compile(r"^\s*(?://|#\[|#!\[)")


def strip_comments(text: str) -> str:
    """Remove // and /* */ comments. Keeps line count stable."""
    out = []
    i = 0
    n = len(text)
    while i < n:
        c = text[i]
        if c == '"':  # skip string literal
            out.append(c)
            i += 1
            while i < n:
                out.append(text[i])
                if text[i] == "\\":
                    i += 2
                    if i - 1 < n:
                        out.append(text[i - 1] if i - 1 < n else "")
                    continue
                if text[i] == '"':
                    i += 1
                    break
                i += 1
            continue
        if c == "/" and i + 1 < n and text[i + 1] == "/":
            while i < n and text[i] != "\n":
                i += 1
            continue
        if c == "/" and i + 1 < n and text[i + 1] == "*":
            i += 2
            while i + 1 < n and not (text[i] == "*" and text[i + 1] == "/"):
                if text[i] == "\n":
                    out.append("\n")
                i += 1
            i += 2
            continue
        out.append(c)
        i += 1
    return "".join(out)


def parse_units(lines):
    """Top-level fn / macro_rules! definitions with inclusive line ranges.

    A unit's range is extended upward over its doc comment and attributes so
    that editing the docs counts as touching the function.
    """
    units = []
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        m = DEF_RE.match(line) or MACRO_RE.match(line) or MOD_RE.match(line)
        cm = None if m else CONST_RE.match(line)
        if not m and not cm:
            i += 1
            continue
        if cm is not None:
            # A const ends where its brackets balance on a line carrying the
            # terminating semicolon -- covers both the one-liner and the
            # multi-line coefficient array.
            name = cm.group("name")
            depth = 0
            end = i
            for j in range(i, n):
                depth += lines[j].count("[") - lines[j].count("]")
                depth += lines[j].count("(") - lines[j].count(")")
                if depth <= 0 and lines[j].rstrip().endswith(";"):
                    end = j
                    break
            else:
                end = i
            units.append(
                {
                    "name": name,
                    "kind": "const",
                    "pub": False,  # shared unit, attributed to its referrers
                    "start": i + 1,
                    "end": end + 1,
                    "def_line": i + 1,
                }
            )
            i = end + 1
            continue
        kind = (
            "fn"
            if DEF_RE.match(line)
            else ("macro" if MACRO_RE.match(line) else "mod")
        )
        name = m.group("name")
        is_pub = bool(DEF_RE.match(line) and m.group("vis"))
        # find the closing brace at column 0
        end = i
        j = i
        while j < n:
            if lines[j] == "}" or lines[j].startswith("} "):
                end = j
                break
            j += 1
        else:
            end = n - 1
        # extend upward over attributes and doc comments
        start = i
        k = i - 1
        while k >= 0 and ATTR_RE.match(lines[k]) and lines[k].strip():
            start = k
            k -= 1
        units.append(
            {
                "name": name,
                "kind": kind,
                "pub": is_pub,
                "start": start + 1,  # 1-indexed, inclusive
                "end": end + 1,
                "def_line": i + 1,
            }
        )
        i = end + 1
    return units


def build_graph(units, clean_lines):
    names = {u["name"] for u in units}
    edges = {}
    for u in units:
        body = "\n".join(clean_lines[u["def_line"] : u["end"]])
        found = set()
        for tok in re.findall(r"\b(\w+)\b", body):
            if tok in names and tok != u["name"]:
                found.add(tok)
        edges[u["name"]] = found
    return edges


def reachable(start, edges, unit_by_name):
    """Private units reachable from `start`, stopping at public boundaries.

    Traversal does not expand another public function's body. If `sinh` calls
    `exp2` and `exp2` calls `exp2_field_split`, that helper belongs to `exp2`,
    not to `sinh` -- `exp2` is a published, separately tested contract, and
    expanding through it welds every exp-descended function into one blob.
    Cross-boundary calls are reported separately as `depends_on`.
    """
    seen = set()
    stack = list(edges.get(start, ()))
    while stack:
        v = stack.pop()
        if v in seen:
            continue
        seen.add(v)
        u = unit_by_name.get(v)
        if u is not None and u["pub"]:
            continue  # stop at the public contract
        stack.extend(edges.get(v, ()))
    return seen


class DSU:
    def __init__(self):
        self.p = {}

    def find(self, x):
        self.p.setdefault(x, x)
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[rb] = ra


def main():
    # Accepts a repo root or a lib.rs path directly. `jm` passes a snapshot of
    # *master's* lib.rs so every worktree derives the same map -- an instance's
    # own uncommitted edits must not reshape ownership under other instances.
    arg = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
    src = arg if arg.is_file() else arg / "src" / "lib.rs"
    text = src.read_text()
    raw_lines = text.split("\n")
    clean_lines = strip_comments(text).split("\n")

    units = parse_units(raw_lines)
    unit_by_name = {u["name"]: u for u in units}
    # `mod tests` and anything inside it is not part of the API surface
    units = [u for u in units if u["kind"] != "mod"]
    unit_by_name = {u["name"]: u for u in units}

    edges = build_graph(units, clean_lines)
    pubs = [u["name"] for u in units if u["pub"]]
    npub = len(pubs)

    # what each public function transitively depends on
    reach = {p: reachable(p, edges, unit_by_name) for p in pubs}

    # how many public functions reach each unit
    reachers = {}
    for p in pubs:
        for v in reach[p]:
            reachers.setdefault(v, set()).add(p)

    core_cut = max(CORE_REACH_MIN, int(npub * CORE_REACH_FRACTION))
    core = {
        v
        for v, rs in reachers.items()
        if len(rs) >= core_cut and v in unit_by_name and not unit_by_name[v]["pub"]
    }

    # Union public functions that share a non-core PRIVATE unit. Sharing a
    # public dependency deliberately does not union: `sinh` and `erf` both
    # reaching `exp2` is a contract dependency, not joint ownership, and
    # unioning on it collapses 123 of 171 functions into one domain that
    # nobody could ever claim.
    dsu = DSU()
    for p in pubs:
        dsu.find(p)
    for v, rs in reachers.items():
        if v in core:
            continue
        u = unit_by_name.get(v)
        if u is None or u["pub"]:
            continue
        rs = sorted(rs)
        for other in rs[1:]:
            dsu.union(rs[0], other)

    groups = {}
    for p in pubs:
        groups.setdefault(dsu.find(p), []).append(p)

    # name each domain after its lowest-line public function, and collect the
    # private units it exclusively owns
    domains = {}
    for members in groups.values():
        members = sorted(members, key=lambda m: unit_by_name[m]["def_line"])
        label = members[0]
        owned = set()
        for m in members:
            for v in reach[m]:
                if v in core:
                    continue
                if v in unit_by_name and not unit_by_name[v]["pub"]:
                    owned.add(v)
        ranges = []
        for name in list(members) + sorted(owned):
            u = unit_by_name[name]
            ranges.append([u["start"], u["end"], name])
        domains[label] = {
            "public": members,
            "private": sorted(owned),
            "ranges": sorted(ranges),
            "lines": sum(r[1] - r[0] + 1 for r in ranges),
        }

    # Advisory cross-domain edges: D calls a public function owned by E. Not a
    # claim conflict (E's contract is stable), but if you change E's accuracy
    # or shape, D's numbers move -- so `jm claim` reports it.
    domain_of = {}
    for label, d in domains.items():
        for m in d["public"]:
            domain_of[m] = label
    for label, d in domains.items():
        calls, called_by = set(), set()
        for m in d["public"]:
            for v in edges.get(m, ()):
                o = domain_of.get(v)
                if o and o != label:
                    calls.add(o)
        for other, od in domains.items():
            if other == label:
                continue
            for m in od["public"]:
                for v in edges.get(m, ()):
                    if domain_of.get(v) == label:
                        called_by.add(other)
        d["calls"] = sorted(calls)
        d["called_by"] = sorted(called_by)

    out = {
        "source": str(src),
        "n_public": npub,
        "core": sorted(core),
        "core_ranges": sorted(
            [unit_by_name[c]["start"], unit_by_name[c]["end"], c] for c in core
        ),
        "domains": domains,
    }
    json.dump(out, sys.stdout, indent=1)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
