#!/usr/bin/env bash
# Stage 03e aggregator (login node) — concatenate per-sample census + analyse.
set -euo pipefail
source "$(dirname "$0")/lib/common.sh"
source "$HOME/opt/miniforge3/etc/profile.d/conda.sh" && conda activate wolb-r
cd "$ROOT"
out="results/03_profile"
ls "$out"/*.census.tsv >/dev/null 2>&1 || die "no *.census.tsv — run 'bash workflow/03e_census.sh' (array) first"
cat "$out"/*.census.tsv > "$out/census_long.tsv"
need "$out/census_long.tsv"
Rscript workflow/03e_census_analysis.R
need reports/stage03e_census.md \
     results/03_profile/per_sample_census.tsv \
     results/03_profile/genus_wolbachia.tsv \
     reports/figs/census_composition.png
