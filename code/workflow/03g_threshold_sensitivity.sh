#!/usr/bin/env bash
# Stage 03g (A6) — threshold sensitivity of the necessity result (login node, no SLURM).
# Sweeps Wolbachia breadth x rpm and parasitoid para_frac thresholds and shows the
# Wolbachia+/parasitoid- counterexample count stays 0 across a wide region.
# Pure reanalysis of results/03_profile/per_sample_profile.tsv.
set -euo pipefail
source "$(dirname "$0")/lib/common.sh"
source "$HOME/opt/miniforge3/etc/profile.d/conda.sh" && conda activate wolb-r   # R lives in wolb-r env
cd "$ROOT"; need results/03_profile/per_sample_profile.tsv
Rscript workflow/03g_threshold_sensitivity.R
need reports/stage03g_threshold_sensitivity.md \
     results/03_profile/threshold_sweep.tsv \
     reports/figs/threshold_sensitivity.png
