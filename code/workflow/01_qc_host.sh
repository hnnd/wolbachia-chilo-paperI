#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/lib/common.sh"
N=$(n_samples); RANGE=${1:-1-$N}   # pass "1-6" to run only the pilot
log "submitting QC+host array $RANGE (of $N)"
sbatch --array="$RANGE%12" -J qc_host "$ROOT/workflow/jobs/qc_host_one.sh"
