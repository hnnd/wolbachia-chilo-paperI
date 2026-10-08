#!/usr/bin/env bash
# Stage 03c submitter — A5-full dual-source mitochondrial mapping over the cohort.
# Readiness gate: build the combined {Chilo + parasitoid} mito index before
# launching the array, so we don't fan out 247 jobs that all fail on a missing
# reference (matches the fail-fast discipline of 02c_kraken.sh).
set -euo pipefail
source "$(dirname "$0")/lib/common.sh"

ref="$ROOT/refs/mito_dualsource/mito_combined"
if [[ ! -s $ref.l2b || ! -s $ref.mbw ]]; then
  log "combined mito index missing — building it now"
  bash "$ROOT/workflow/lib/build_mito_ref.sh"
fi
need "$ref.l2b" "$ref.mbw"

N=$(n_samples); RANGE=${1:-1-$N}
sbatch --array="$RANGE%12" -J mito "$ROOT/workflow/jobs/mito_one.sh"
