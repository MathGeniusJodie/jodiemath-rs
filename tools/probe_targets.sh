#!/usr/bin/env bash
# Codegen of the trig throughput loops across SIMD targets: counts gathers,
# 64-bit-lane ops (u64 multiplies and their widen/narrow shuffles) and scalar
# FP ops in each loop. A scratch crate includes src/lib.rs by path, so cross
# targets need only their std (rustup target add aarch64-unknown-linux-gnu).
#   tools/probe_targets.sh
set -euo pipefail
repo=$(cd "$(dirname "$0")/.." && pwd)
P=$(mktemp -d)
trap 'rm -rf "$P"' EXIT
mkdir -p "$P/src"
cat > "$P/Cargo.toml" <<TOML
[package]
name = "probe"
version = "0.1.0"
edition = "2021"
TOML
{
  echo '#![allow(dead_code, unused_imports)]'
  echo "#[path = \"$repo/src/lib.rs\"] mod jm;"
  for f in sin_wide cos_wide tan_wide sin cos tan; do
    echo "#[no_mangle] pub fn thr_$f(i: &[f32; 1024], o: &mut [f32; 1024]) { for (o, &x) in o.iter_mut().zip(i.iter()) { *o = jm::$f(x); } }"
  done
} > "$P/src/lib.rs"
cp "$repo/rust-toolchain.toml" "$P/"
cd "$P"
emit() {
  local tag=$1 target=$2; shift 2
  CARGO_TARGET_DIR="$P/target-$tag" RUSTFLAGS="$*" cargo rustc --release --lib -q ${target:+--target $target} -- --emit=asm 2>&1 | grep -E "^error" -A5 || true
  local s
  s=$(find "$P/target-$tag" -name 'probe-*.s' | head -1)
  for f in sin_wide cos_wide tan_wide sin cos tan; do
    body=$(awk -v f="thr_$f" '$0 ~ "^"f":" {p=1; next} p && /^\s*(retq|ret)\s*$/ {exit} p' "$s" | grep -v '^\s*[.#]' | grep -v '^\.L\|^\s*$' || true)
    printf "%-10s %-9s instrs %4d  gathers %2d  64-bit-lane %2d  scalar-fp %3d\n" "$tag" "$f" \
      "$(printf "%s\n" "$body" | wc -l)" \
      "$(printf "%s\n" "$body" | grep -ciE 'gather' || true)" \
      "$(printf "%s\n" "$body" | grep -cE 'pmuludq|vpmovzxdq|vpmovqd|umull|uzp1.*\.4s' || true)" \
      "$(printf "%s\n" "$body" | grep -cE '^\s*(vmovss|vmulss|vfmadd...ss|fmadd\s+s|fmul\s+s)' || true)"
  done
}
emit avx512 "" -C target-cpu=native
emit avx2 "" -C target-cpu=x86-64-v3
# aarch64 has no `fma` target feature; `doctest` skips lib.rs's x86 FMA guard.
emit neon aarch64-unknown-linux-gnu -C target-feature=+neon --cfg doctest
