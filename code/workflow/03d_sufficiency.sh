#!/usr/bin/env bash
# Stage 03d — necessity + detection-limited sufficiency (login node, no SLURM).
# Core evidence for evidence line A: parasitoid presence is necessary and
# (detection-limited) sufficient for Wolbachia detection. Pure reanalysis of
# results/03_profile/per_sample_profile.tsv.
set -euo pipefail
source "$(dirname "$0")/lib/common.sh"
source "$HOME/opt/miniforge3/etc/profile.d/conda.sh" && conda activate wolb-r   # R lives in wolb-r env
cd "$ROOT"; need results/03_profile/per_sample_profile.tsv
Rscript workflow/03d_sufficiency.R
need reports/stage03d_sufficiency.md \
     results/03_profile/sufficiency_stats.tsv \
     results/03_profile/positivity_by_load.tsv \
     reports/figs/positivity_curve.png reports/figs/detection_threshold_scatter.png
