#!/usr/bin/env python3
"""Weighted Remez fits for the [-pi/4, pi/4] sin/cos pair used by tan_wide.

  sin(r) ~= r + r^3 * S(r^2),   cos(r) ~= 1 + r^2 * C(r^2)

minimising the relative error of sin/cos over |r| <= R. Prints f32
coefficients (nearest) in Rust syntax.
"""
import sys
import mpmath as mp

mp.mp.dps = 50
R = mp.mpf(sys.argv[1]) if len(sys.argv) > 1 else mp.pi / 4 * (1 + mp.mpf(2) ** -20)
NS = int(sys.argv[2]) if len(sys.argv) > 2 else 3
NC = int(sys.argv[3]) if len(sys.argv) > 3 else 4

def remez(g, w, n, a, b, iters=30):
    """min over deg-(n-1) p of max |w(y) (g(y) - p(y))| on [a, b]."""
    m = n + 1
    xs = [(a + b) / 2 + (b - a) / 2 * mp.cos(mp.pi * (m - 1 - i) / (m - 1)) for i in range(m)]
    for _ in range(iters):
        A = mp.matrix(m, m)
        rhs = mp.matrix(m, 1)
        for i, x in enumerate(xs):
            for j in range(n):
                A[i, j] = x ** j
            A[i, n] = (-1) ** i / w(x)
            rhs[i] = g(x)
        sol = mp.lu_solve(A, rhs)
        c = [sol[j] for j in range(n)]
        err = lambda x: w(x) * (g(x) - sum(c[j] * x ** j for j in range(n)))
        # locate extrema on a fine grid, then refine
        grid = [a + (b - a) * k / 4000 for k in range(4001)]
        vals = [err(x) for x in grid]
        ext = [grid[0]]
        for k in range(1, len(grid) - 1):
            if (vals[k] - vals[k - 1]) * (vals[k + 1] - vals[k]) <= 0:
                ext.append(grid[k])
        ext.append(grid[-1])
        # keep m alternating extrema of largest magnitude
        while len(ext) > m:
            es = [abs(err(x)) for x in ext]
            k = es.index(min(es))
            ext.pop(k)
        if len(ext) == m:
            xs = ext
    return c, max(abs(err(x)) for x in [a + (b - a) * k / 20000 for k in range(20001)])

def f32(v):
    import struct
    return struct.unpack("<f", struct.pack("<f", float(v)))[0]

Y = R * R
eps = mp.mpf(10) ** -40
gs = lambda y: (mp.sin(mp.sqrt(y)) / mp.sqrt(y) - 1) / y if y > eps else mp.mpf(-1) / 6 + y / 120
ws = lambda y: y / (mp.sin(mp.sqrt(y)) / mp.sqrt(y)) if y > eps else y
gc = lambda y: (mp.cos(mp.sqrt(y)) - 1) / y if y > eps else mp.mpf(-1) / 2 + y / 24
wc = lambda y: y / mp.cos(mp.sqrt(y)) if y > eps else y
cs, es = remez(gs, ws, NS, Y * mp.mpf(10) ** -6, Y)
cc, ec = remez(gc, wc, NC, Y * mp.mpf(10) ** -6, Y)
print(f"// R = {mp.nstr(R, 10)}: sin rel err {mp.nstr(es, 3)} (2^{mp.nstr(mp.log(es, 2), 4)}), cos rel err {mp.nstr(ec, 3)} (2^{mp.nstr(mp.log(ec, 2), 4)})")
print("let s: [f32; %d] = [%s];" % (NS, ", ".join(repr(f32(v)) for v in cs)))
print("let c: [f32; %d] = [%s];" % (NC, ", ".join(repr(f32(v)) for v in cc)))
