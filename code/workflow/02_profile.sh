#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/lib/common.sh"
N=$(n_samples); RANGE=${1:-1-$N}
sbatch --array="$RANGE%12" -J profile "$ROOT/workflow/jobs/profile_one.sh"
