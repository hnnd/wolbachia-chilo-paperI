#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CFG="$ROOT/config/config.yaml"
log(){ echo "[$(date '+%F %T')] $*" >&2; }
die(){ echo "[FATAL] $*" >&2; exit 1; }
cfg(){ awk -F': *' -v k="$1" '$1==k{sub(/[ \t]+#.*$/,"",$2); gsub(/^[ \t]+|[ \t]+$/,"",$2); print $2; exit}' "$CFG"; }
need(){ for f in "$@"; do [[ -s $f ]] || die "missing/empty: $f"; done; }
mark_done(){ touch "$1.done"; }
n_samples(){ echo $(($(wc -l < "$ROOT/config/samples.tsv")-1)); }
# 取第 i 个样本（1-based，跳过表头）的第 col 列
sample_field(){ awk -F'\t' -v i="$1" -v c="$2" 'NR==i+1{print $c}' "$ROOT/config/samples.tsv"; }
