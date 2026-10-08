#!/usr/bin/env bash
# Stage 03e submitter — per-sample parasitoid census (broad mito reference).
# Readiness gate: build the broad reference before fanning out the array.
set -euo pipefail
source "$(dirname "$0")/lib/common.sh"

ref="$ROOT/refs/parasitoid_broad/parasitoid_mito_broad"
if [[ ! -s $ref.l2b || ! -s $ref.mbw ]]; then
  log "broad parasitoid mito index missing — building it now"
  bash "$ROOT/workflow/lib/fetch_broad_parasitoid_mito.sh"
fi
need "$ref.l2b" "$ref.mbw"

N=$(n_samples); RANGE=${1:-1-$N}
sbatch --array="$RANGE%12" -J census "$ROOT/workflow/jobs/census_one.sh"
