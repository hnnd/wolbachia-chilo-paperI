#!/usr/bin/env bash
# Stage 03h (A4) — multivariable logistic controlling pan-abundance confounders
# (login node, no SLURM). Tests whether parasitoid load remains the dominant
# predictor of Wolbachia positivity after adjusting for sequencing depth, non-moth
# biomass fraction and host plant. Pure reanalysis of per_sample_profile.tsv.
set -euo pipefail
source "$(dirname "$0")/lib/common.sh"
source "$HOME/opt/miniforge3/etc/profile.d/conda.sh" && conda activate wolb-r   # R lives in wolb-r env
cd "$ROOT"; need results/03_profile/per_sample_profile.tsv
Rscript workflow/03h_multivariate.R
need reports/stage03h_multivariate.md \
     results/03_profile/multivariate_coef.tsv \
     reports/figs/multivariate_or.png
