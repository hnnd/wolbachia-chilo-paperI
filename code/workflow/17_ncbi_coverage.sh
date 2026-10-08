#!/usr/bin/env bash
# Stage 17 — NCBI submission coverage.  Submitter only (login node).
#
# Before running, the three submission FASTAs must be in results/17_ncbi_cov/
# (docs/submission/ncbi/genomes/*.fsa mirrored to the cluster):
#   wCchiA.fsa  wCchiB.fsa  wCchiA_unplaced.fsa
#
# Result (results/17_ncbi_cov/coverage_summary.tsv) supplies the mandatory
# "Genome coverage" field of the NCBI Genome Info table.
set -euo pipefail
source "$(dirname "$0")/lib/common.sh"
out="$ROOT/results/17_ncbi_cov"
for f in wCchiA.fsa wCchiB.fsa wCchiA_unplaced.fsa; do need "$out/$f"; done
sbatch -J ncbi_cov "$ROOT/workflow/jobs/ncbi_cov_one.sh"
