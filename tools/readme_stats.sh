#!/usr/bin/env bash
# Re-run every readme number and rewrite readme.md's generated tables:
# accuracy (exhaustive 2^32 sweep + powfsearch + clogsearch), llvm-mca, and
# quickbench wall-clock (3 rounds under the jm bench lock). ~10 minutes on the 32-core
# Zen 5 box.
#   tools/readme_stats.sh            run everything, then rewrite readme.md
#   tools/readme_stats.sh --tables   only rewrite readme.md from the last run
set -euo pipefail
cd "$(dirname "$0")/.."
if [ -z "${IN_NIX_SHELL:-}" ]; then exec nix-shell --run "$(printf '%q ' "$0" "$@")"; fi
out=target/readme-stats

if [ "${1:-}" != "--tables" ]; then
  mkdir -p "$out"
  cargo build --release --example accuracy --example powfsearch --example clogsearch \
    --example mca --example quickbench
  ex=target/release/examples
  {
    echo "$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2 | xargs), $(nproc) threads"
    rustc --version
    echo "llvm-mca $(llvm-mca --version | sed -n 's/.*LLVM version \(.*\)/\1/p' | xargs), -mcpu=$(llvm-mca --version | sed -n 's/.*Host CPU: //p')"
    echo "RUSTFLAGS from .cargo/config.toml: $(sed -n 's/^rustflags = //p' .cargo/config.toml)"
    echo "commit $(git rev-parse --short HEAD)$(git diff --quiet HEAD -- src || echo " + uncommitted src changes"), $(date -u +%Y-%m-%d)"
  } > "$out/machine.txt"
  "$ex/accuracy" thorough > "$out/accuracy.txt"
  "$ex/powfsearch" > "$out/powfsearch.txt"
  "$ex/clogsearch" > "$out/clogsearch.txt"
  "$ex/mca" > "$out/mca.txt"
  # Min over rounds: a single min-of-7 still swings up to ~1.5x on a few rows.
  rm -f "$out"/quickbench*.txt
  for round in 1 2 3; do
    ./tools/jm bench -- "$ex/quickbench" > "$out/quickbench.$round.txt"
  done
fi

python3 tools/readme_tables.py "$out" readme.md
