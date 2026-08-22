#!/usr/bin/env python3
"""Validate + regenerate src/pitable.rs tables for reduce_pi_wide.

beta(e) = (2^(e-150)/pi) mod 2, sliced into u32 chunks at fixed bit
positions. Validates against the shipped layout (units 2^-28 / 2^-29 /
2^-58, bits down to 2^-86), then emits the proposed uniform-scale layout.

Exact integer arithmetic throughout: pi to S bits via Machin's formula,
so every comparison against the shipped tables is exact, not float-y.
"""
from fractions import Fraction

S = 600   # pi known to 2^-S; swamps every N we slice at
N = 100   # bits of beta carried in `full`

def arctan_inv(x, one):
    """atan(1/x) * one, integer."""
    total = term = one // x
    x2 = x * x
    k = 3
    sign = -1
    while term != 0:
        term = term // x2
        total += sign * (term // k)
        sign = -sign
        k += 2
    return total

def pi_scaled():
    one = 1 << S
    return 16 * arctan_inv(5, one) - 4 * arctan_inv(239, one)  # ~= pi * 2^S

PI_S = pi_scaled()

def beta_full(e):
    """floor(beta(e) * 2^N), beta = (2^(e-150)/pi) mod 2, exact.

    beta < 2, so this needs N+1 bits -- the top bit IS beta's integer
    part (the parity bit). Reducing mod 2^N would chop it off.
    """
    # num // PI_S ~= 2^(e-150+N)/pi; mod 2^(N+1) implements the mod 2
    num = 1 << (e - 150 + N + S)
    return (num // PI_S) % (1 << (N + 1))

# ---- shipped layout ----
# W0 = floor(V*2^28)                      bits  2^0 .. 2^-28 (unit 2^-28)
# W1 = floor(V*2^57) & (2^29-1)           bits 2^-29 .. 2^-57 (unit 2^-29 -> 2^-57)
# W2 = floor(V*2^86) & (2^29-1)           bits 2^-58 .. 2^-86 (unit 2^-58 -> 2^-86)
import re
shipped = open("src/pitable.rs").read()
tab = shipped[shipped.rindex("= ["):]  # only the literal
# strip comments so header digits don't pollute the scan
tab = re.sub(r"//[^\n]*", "", tab)
nums = re.findall(r"\d+", tab)
assert len(nums) == 768, len(nums)
vals = [int(x) for x in nums]
W0o, W1o, W2o = vals[0:256], vals[256:512], vals[512:768]
bad = 0
for e in range(256):
    f = beta_full(e)
    if (f >> (N - 28)) != W0o[e]: bad += 1; print("W0", e, f >> (N-28), W0o[e])
    if ((f >> (N - 57)) & (2**29 - 1)) != W1o[e]: bad += 1; print("W1", e)
    if ((f >> (N - 86)) & (2**29 - 1)) != W2o[e]: bad += 1; print("W2", e)
print("shipped-layout validation mismatches:", bad)
assert bad == 0

# ---- proposed uniform-scale layout ----
# All three planes consumed as p_i = (m * 2^-28) * W_i, so every chunk is in
# units of 2^-28:
#   W0' = floor(V*2^28)                    bits  2^0 .. 2^-28  (< 2^29)
#   W1' = floor(V*2^56) & (2^28-1)         bits 2^-29 .. 2^-56 (< 2^28)
#   W2' = floor(V*2^85) & (2^29-1)         bits 2^-57 .. 2^-85 (< 2^29)
# Exactness with 24-bit m: m*W0' < 2^53, m*W1' < 2^52, m*W2' < 2^53.
for e in range(256):
    f = beta_full(e)
    w1n = (f >> (N - 56)) & (2**28 - 1)
    w2n = (f >> (N - 85)) & (2**29 - 1)
    assert w1n < 2**28 and w2n < 2**29

def emit(name, rows, comment):
    out = [f"    [ // {comment}"]
    for r in range(0, 256, 8):
        out.append("    " + ", ".join(str(v) for v in rows[r:r+8]) + ",")
    out.append("    ],")
    return "\n".join(out)

planes = []
for idx, (shift, mask) in enumerate([(N-28, None), (N-56, 2**28-1), (N-85, 2**29-1)]):
    rows = []
    for e in range(256):
        f = beta_full(e)
        v = f >> shift
        if mask: v &= mask
        rows.append(v)
    planes.append(rows)

open("/tmp/pitable_new.txt", "w").write(
    emit("P0", planes[0], "plane 0: bits 2^0 .. 2^-28 of beta (top bit = parity)") + "\n" +
    emit("P1", planes[1], "plane 1: bits 2^-29 .. 2^-56 of beta") + "\n" +
    emit("P2", planes[2], "plane 2: bits 2^-57 .. 2^-85 of beta") + "\n")
print("new tables written to /tmp/pitable_new.txt")

# cross-check: new W0' must equal shipped W0
for e in range(256):
    assert planes[0][e] == W0o[e]
print("W0' == shipped W0 for all e")
