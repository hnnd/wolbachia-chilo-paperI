#!/usr/bin/env bash
# Stage 03f submitter (T2) — de novo mito assembly + BLAST species ID.
# Targets: non-Cotesia samples with enough parasitoid-mito signal to assemble
# (where 03e best-hit was unreliable), plus a few Cotesia positive controls.
set -euo pipefail
source "$(dirname "$0")/lib/common.sh"
source "$HOME/opt/miniforge3/etc/profile.d/conda.sh" && conda activate wolb

out="$ROOT/results/06_denovo"; mkdir -p "$out"
ref="$ROOT/refs/parasitoid_broad/parasitoid_mito_broad"
db="$ROOT/refs/parasitoid_broad/mito_blastdb"
pset="$ROOT/refs/parasitoid_broad/parasitoid_only.fa"

# Gate 1: broad mito reference (bait).
if [[ ! -s $ref.l2b ]]; then bash "$ROOT/workflow/lib/fetch_broad_parasitoid_mito.sh"; fi
need "$ref.l2b" "$ref.mbw" "$pset"
# Gate 2: local BLAST DB from the parasitoid mito set.
if [[ ! -s $db.nin && ! -s $db.ndb ]]; then
  log "building local mito BLAST DB"
  makeblastdb -in "$pset" -dbtype nucl -out "$db" -title parasitoid_mito >/dev/null
fi

# Build target list (idempotent): non-Cotesia & total_para>=500, + 3 Cotesia controls.
targets="$out/assemble_targets.txt"
need "$ROOT/results/03_profile/per_sample_census.tsv"
conda activate wolb-r
Rscript -e '
  suppressMessages(library(readr))
  d <- read_tsv("results/03_profile/per_sample_census.tsv", show_col_types=FALSE)
  noncot <- d$sample[d$dom_genus!="Cotesia" & d$total_para>=500]
  ctrl   <- head(d$sample[d$dom_genus=="Cotesia"][order(-d$total_para[d$dom_genus=="Cotesia"])], 3)
  writeLines(unique(c(noncot, ctrl)), "results/06_denovo/assemble_targets.txt")
' 2>/dev/null
need "$targets"
N=$(wc -l < "$targets")
log "Stage 03f targets: $N samples"
RANGE=${1:-1-$N}
sbatch --array="$RANGE%12" -J asmid "$ROOT/workflow/jobs/assemble_mito_one.sh"
