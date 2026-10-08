#!/usr/bin/env bash
# Stage 3b — A3 dose-response + A5-lite source specificity (login node, no SLURM).
# Pure reanalysis of results/03_profile/per_sample_profile.tsv produced by Stage 2b.
set -euo pipefail
source "$(dirname "$0")/lib/common.sh"
source "$HOME/opt/miniforge3/etc/profile.d/conda.sh" && conda activate wolb-r   # R lives in wolb-r env
cd "$ROOT"; need results/03_profile/per_sample_profile.tsv
Rscript workflow/03b_dose_source.R
need reports/stage03b_dose_source.md \
     results/03_profile/dose_response_stats.tsv \
     results/03_profile/source_specificity_stats.tsv \
     reports/figs/dose_response.png reports/figs/source_specificity.png
