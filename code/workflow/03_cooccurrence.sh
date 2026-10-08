#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/lib/common.sh"
source "$HOME/opt/miniforge3/etc/profile.d/conda.sh" && conda activate wolb-r   # R lives in wolb-r env
cd "$ROOT"; need results/03_profile/per_sample_profile.tsv
Rscript workflow/03_cooccurrence.R
need reports/stage03_cooccurrence.md results/03_profile/cooccurrence_stats.tsv \
     reports/figs/cooccurrence_box.png reports/figs/cooccurrence_freq_host.png
