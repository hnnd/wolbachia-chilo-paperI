#!/usr/bin/env bash
#SBATCH -p cpu -c 8 --mem=16G -o logs/mito_%A_%a.out -e logs/mito_%A_%a.err
# Stage 03c worker — dual-source mitochondrial mapping (A5-full).
# Maps non-host reads (which RETAIN Chilo mito because the host assembly has no
# mito contig — verified: refs/host/chilo_suppressalis.ann is chromosomes only)
# to the combined {Chilo + parasitoid} mitogenome reference, and records the
# primary-best-hit read count per mitogenome. Symmetric, mutually-exclusive
# insect-mito counts let Stage 03c test Wolbachia ~ parasitoid-mito controlling
# for Chilo-mito and depth — the clean specificity test A5-lite could not do.
set -euo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}"   # SLURM spools the script; locate repo via submit dir
source "$HERE/workflow/lib/common.sh"
source "$HOME/opt/miniforge3/etc/profile.d/conda.sh" && conda activate wolb
i=${SLURM_ARRAY_TASK_ID:?}; T=$(cfg threads)
sid=$(sample_field "$i" 1)
nh="$ROOT/results/02_nonhost"; out="$ROOT/results/03_profile"; mkdir -p "$out"
r1="$nh/${sid}_R1.fq.gz"; r2="$nh/${sid}_R2.fq.gz"; need "$r1" "$r2"
[[ -f $out/$sid.mito.done ]] && exit 0

ref="$ROOT/refs/mito_dualsource/mito_combined"; need "$ref.l2b" "$ref.mbw"
bam="$out/$sid.mito.bam"
minibwa map -x sr -t "$T" "$ref" "$r1" "$r2" \
  | samtools sort -@ "$T" -o "$bam" - && samtools index "$bam"

# Per-contig primary read counts (exclude secondary/supplementary/unmapped:
# -F 0x904). Counts are loaded via getline in BEGIN so a zero-mito sample (empty
# count file) still yields one row per contig — the NR==FNR two-file idiom breaks
# when the first file is empty, silently emitting nothing. idxstats is the main
# input (always lists all 12 contigs), so the row-per-contig output is guaranteed.
samtools view -F 0x904 "$bam" | cut -f3 | sort | uniq -c \
  | awk '{print $2"\t"$1}' > "$out/$sid.mito.cnt"     # contig<TAB>count (may be empty)
samtools idxstats "$bam" | awk -F'\t' -v s="$sid" -v cf="$out/$sid.mito.cnt" '
  BEGIN{ while((getline l < cf) > 0){ split(l, b, "\t"); cnt[b[1]] = b[2] } }
  $1!="*"{ split($1, a, "|"); print s"\t"a[1]"\t"$1"\t"$2"\t"(($1 in cnt)?cnt[$1]:0) }
' > "$out/$sid.mito.tsv"
need "$out/$sid.mito.tsv"

rm -f "$bam"* "$out/$sid.mito.cnt"
mark_done "$out/$sid.mito"
