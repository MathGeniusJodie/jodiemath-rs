#!/usr/bin/env bash
# Interleaved A/B of examples/trig_bench.rs: the committed HEAD (built in a
# scratch worktree) against the working tree. Prints each side's minimum
# over the rounds.  tools/ab.sh [rounds] [filter]
set -euo pipefail
cd "$(dirname "$0")/.."
rounds=${1:-3}; filter=${2:-}
base=$(git rev-parse --show-toplevel)/../.jm-ab-base
if [ ! -d "$base" ]; then git worktree add -f --detach "$base" HEAD >/dev/null 2>&1; fi
git -C "$base" checkout -q --detach "$(git rev-parse HEAD)"
cp examples/trig_bench.rs "$base/examples/"
(cd "$base" && cargo build --release -q --example trig_bench 2>&1 | grep -E "^error" -A5 || true)
cargo build --release -q --example trig_bench 2>&1 | grep -E "^error" -A5 || true
for i in $(seq "$rounds"); do
  "$base/target/release/examples/trig_bench" "$filter" | sed 's/^/A /'
  ./target/release/examples/trig_bench "$filter" | sed 's/^/B /'
done | awk '{k=$1" "$2; if(!(k in l)||$4<l[k])l[k]=$4; if(!(k in t)||$7<t[k])t[k]=$7; if(!(k in o)){o[k]=n++; ks[n-1]=k}}
  END{for(i=0;i<n;i++){k=ks[i]; printf "%-12s lat %6.2f  thr %6.3f\n", k, l[k], t[k]}}' | sort -k2,2 -k1,1
