#!/usr/bin/env bash
#SBATCH -p cpu -c 8 --mem=16G -o logs/census_%A_%a.out -e logs/census_%A_%a.err
# Stage 03e worker — per-sample parasitoid census against a BROAD mitogenome
# reference (122 parasitoid mitos across families + Chilo). Maps non-host reads,
# records primary-best-hit read count per contig. Aggregation to genus (done in
# 03e_census_analysis.R) is robust to within-genus cross-mapping; the dominant
# taxon is reliable, the low tail may be conserved-region leakage.
set -euo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}"
source "$HERE/workflow/lib/common.sh"
source "$HOME/opt/miniforge3/etc/profile.d/conda.sh" && conda activate wolb
i=${SLURM_ARRAY_TASK_ID:?}; T=$(cfg threads)
sid=$(sample_field "$i" 1)
nh="$ROOT/results/02_nonhost"; out="$ROOT/results/03_profile"; mkdir -p "$out"
r1="$nh/${sid}_R1.fq.gz"; r2="$nh/${sid}_R2.fq.gz"; need "$r1" "$r2"
[[ -f $out/$sid.census.done ]] && exit 0

ref="$ROOT/refs/parasitoid_broad/parasitoid_mito_broad"; need "$ref.l2b" "$ref.mbw"
bam="$out/$sid.census.bam"
minibwa map -x sr -t "$T" "$ref" "$r1" "$r2" \
  | samtools sort -@ "$T" -o "$bam" - && samtools index "$bam"

# Primary best-hit counts per contig; getline-loaded so zero-hit samples still
# emit one row per contig (see slurm-array-empty-file-awk-bug).
samtools view -F 0x904 "$bam" | cut -f3 | sort | uniq -c \
  | awk '{print $2"\t"$1}' > "$out/$sid.census.cnt"
samtools idxstats "$bam" | awk -F'\t' -v s="$sid" -v cf="$out/$sid.census.cnt" '
  BEGIN{ while((getline l < cf) > 0){ split(l, b, "\t"); cnt[b[1]] = b[2] } }
  $1!="*" && $3+$4>=0 { split($1, a, "|"); r=($1 in cnt)?cnt[$1]:0;
    if (r > 0) print s"\t"a[1]"\t"$1"\t"r }   # keep only contigs with reads (broad ref is large)
' > "$out/$sid.census.tsv"
# Guarantee a non-empty file even when nothing maps (avoids need failure).
[[ -s $out/$sid.census.tsv ]] || printf "%s\tNONE\tNONE\t0\n" "$sid" > "$out/$sid.census.tsv"

rm -f "$bam"* "$out/$sid.census.cnt"
mark_done "$out/$sid.census"
