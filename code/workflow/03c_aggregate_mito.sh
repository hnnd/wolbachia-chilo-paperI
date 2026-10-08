#!/usr/bin/env bash
# Stage 03c aggregator (login node, no SLURM) — concatenate per-sample mito
# counts and run the A5-full specificity analysis.
set -euo pipefail
source "$(dirname "$0")/lib/common.sh"
source "$HOME/opt/miniforge3/etc/profile.d/conda.sh" && conda activate wolb-r
cd "$ROOT"
out="results/03_profile"
ls "$out"/*.mito.tsv >/dev/null 2>&1 || die "no *.mito.tsv found — run 'bash workflow/03c_mito_dualsource.sh' (array) first"
cat "$out"/*.mito.tsv > "$out/mito_long.tsv"
need "$out/mito_long.tsv"
Rscript workflow/03c_mito_analysis.R
need reports/stage03c_mito_dualsource.md \
     results/03_profile/per_sample_mito.tsv \
     results/03_profile/mito_specificity_stats.tsv \
     reports/figs/mito_dualsource.png
