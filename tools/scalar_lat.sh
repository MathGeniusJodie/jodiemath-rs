#!/usr/bin/env bash
# Scalar latency of the wrappers in examples/scalar_asm.rs under llvm-mca.
# xmm0 feeds back into itself, so cycles/iteration is the latency with the
# input's exponent unknown -- unlike mca_target's latency regions, whose
# band-clamped chains let LLVM constant-fold everything exponent-driven
# (for the wide tier: the whole 1/pi window).
#   tools/scalar_lat.sh [scalar_sin_wide ...]
set -euo pipefail
cd "$(dirname "$0")/.."
touch examples/scalar_asm.rs
cargo rustc --release --example scalar_asm -q -- --emit=asm -C debuginfo=0 2>&1 | grep -E "^error" -A5 || true
S=$(ls -t $(find target/release -name 'scalar_asm*.s') | head -1)
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
for f in "${@:-scalar_sin_wide scalar_cos_wide scalar_tan_wide}"; do
  for g in $f; do
    awk -v f="$g" '$0 ~ "^"f":" {p=1; next} p && /^\s*retq/ {exit} p' "$S" | grep -v '^\s*\.' | grep -v '^\.L' > "$tmp/$g.s"
    printf "%-18s " "$g"
    llvm-mca -mcpu=native -iterations=100 "$tmp/$g.s" 2>/dev/null |
      awk '/^Instructions:/{i=$2} /^Total Cycles:/{c=$3} END{printf "instrs %d  cycles %.1f\n", i/100, c/100}'
  done
done
