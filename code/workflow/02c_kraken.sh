#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/lib/common.sh"
# Stage 2c submitter: Kraken2 taxonomic profile over the cohort.
# Refuses to launch unless the DB is fully present (avoids 247 failing jobs).
kdb="$(cfg kraken_db)"
[[ -f "$kdb/hash.k2d" && -f "$kdb/opts.k2d" && -f "$kdb/taxo.k2d" ]] \
  || die "kraken2 DB not fully present at $kdb (need hash/opts/taxo.k2d); wait for copy to finish"
N=$(n_samples); RANGE=${1:-1-$N}
log "submitting Kraken2 array $RANGE (of $N), %3 (each task ~150G RAM)"
sbatch --array="$RANGE%3" -J kraken "$ROOT/workflow/jobs/kraken_one.sh"
