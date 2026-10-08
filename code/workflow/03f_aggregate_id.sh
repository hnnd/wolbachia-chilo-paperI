#!/usr/bin/env bash
# Stage 03f aggregator (login node) — collate de novo BLAST IDs.
set -euo pipefail
source "$(dirname "$0")/lib/common.sh"
source "$HOME/opt/miniforge3/etc/profile.d/conda.sh" && conda activate wolb-r
cd "$ROOT"
out="results/06_denovo"
ls "$out"/*.blast.tsv >/dev/null 2>&1 || die "no *.blast.tsv — run 'bash workflow/03f_assemble_id.sh' first"
# Prefix each BLAST row with its sample id.
: > "$out/denovo_blast_long.tsv"
for f in "$out"/*.blast.tsv; do
  sid=$(basename "$f" .blast.tsv)
  [[ -s $f ]] && awk -v s="$sid" -F'\t' '{print s"\t"$0}' "$f" >> "$out/denovo_blast_long.tsv"
done
need "$out/denovo_blast_long.tsv"
Rscript workflow/03f_id_analysis.R
need reports/stage03f_denovo_id.md results/06_denovo/per_sample_denovo_id.tsv
