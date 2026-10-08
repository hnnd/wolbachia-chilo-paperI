#!/usr/bin/env bash
#SBATCH -p cpu -c 16 --mem=150G -o logs/kraken_%A_%a.out -e logs/kraken_%A_%a.err
set -euo pipefail
# Stage 2c sidecar: Kraken2 taxonomic profile on non-host reads.
# Separate from Stage 2 because the standard DB (~104GB) needs ~150G RAM.
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}"   # SLURM copies script to spool; locate repo via submit dir
source "$HERE/workflow/lib/common.sh"
source "$HOME/opt/miniforge3/etc/profile.d/conda.sh" && conda activate wolb
i=${SLURM_ARRAY_TASK_ID:?}; T=$(cfg threads)
sid=$(sample_field "$i" 1)
nh="$ROOT/results/02_nonhost"; out="$ROOT/results/03_profile"; mkdir -p "$out"
r1="$nh/${sid}_R1.fq.gz"; r2="$nh/${sid}_R2.fq.gz"; need "$r1" "$r2"
[[ -f $out/$sid.kraken.done ]] && { log "$sid kraken done, skip"; exit 0; }
kdb="$(cfg kraken_db)"
# fail loud if DB incomplete — all three index files must be present
[[ -f "$kdb/hash.k2d" && -f "$kdb/opts.k2d" && -f "$kdb/taxo.k2d" ]] \
  || die "kraken2 DB incomplete at $kdb (need hash/opts/taxo.k2d)"
kraken2 --db "$kdb" --threads "$T" --paired "$r1" "$r2" \
  --report "$out/$sid.kraken.report" --output /dev/null
need "$out/$sid.kraken.report"
mark_done "$out/$sid.kraken"
