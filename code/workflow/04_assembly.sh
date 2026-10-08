#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/lib/common.sh"
source "$HOME/opt/miniforge3/etc/profile.d/conda.sh" && conda activate wolb
cd "$ROOT"; mkdir -p results/05_typing
# Write pos_samples.txt (may be empty) before any die so state is inspectable
csvtk -t filter2 -f '$wolb_pos==1' results/03_profile/per_sample_profile.tsv | csvtk -t cut -f sample | csvtk -t del-header > results/05_typing/pos_samples.txt
M=$(wc -l < results/05_typing/pos_samples.txt); log "wolb-positive samples: $M"
[[ $M -ge 1 ]] || die "no positive samples"
sbatch --array="1-$M%6" -J wolb_asm workflow/jobs/assemble_one.sh
