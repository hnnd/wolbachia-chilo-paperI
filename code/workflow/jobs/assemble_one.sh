#!/usr/bin/env bash
#SBATCH -p fat -c 32 --mem=200G -o logs/asm_%A_%a.out -e logs/asm_%A_%a.err
set -euo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}"   # SLURM copies script to spool; locate repo via submit dir
source "$HERE/workflow/lib/common.sh"
source "$HOME/opt/miniforge3/etc/profile.d/conda.sh" && conda activate wolb
T=32
sid=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ROOT/results/05_typing/pos_samples.txt")
[[ -n "$sid" ]] || die "empty sid for SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"
nh="$ROOT/results/02_nonhost"; ad="$ROOT/results/04_assembly/$sid"; td="$ROOT/results/05_typing"
mkdir -p "$ad" "$td"; [[ -f $ad/assembly.done ]] && exit 0
# enrich reads that map to Wolbachia, then assemble
minibwa map -x sr -t "$T" refs/wolbachia/wolbachia_refs "$nh/${sid}_R1.fq.gz" "$nh/${sid}_R2.fq.gz" \
  | samtools view -@ "$T" -bF 12 - | samtools sort -n -@ "$T" -o "$ad/wolb.bam" -
samtools fastq -@ "$T" -1 "$ad/w_R1.fq.gz" -2 "$ad/w_R2.fq.gz" -0 /dev/null -s /dev/null "$ad/wolb.bam"
spades.py --isolate -t "$T" -m 200 -1 "$ad/w_R1.fq.gz" -2 "$ad/w_R2.fq.gz" -o "$ad/spades"
if [[ ! -s "$ad/spades/scaffolds.fasta" ]]; then
  log "SPAdes produced no scaffolds for $sid (insufficient Wolbachia reads?) — not marking done"
  exit 0
fi
ln -sf "$ad/spades/scaffolds.fasta" "$ad/assembly.fasta"
# quality + typing (all best-effort; failures do not abort)
quast.py -o "$ad/quast" "$ad/assembly.fasta" >/dev/null 2>&1 || true
mlst "$ad/assembly.fasta" > "$td/$sid.mlst.tsv" || true
# wsp: blastn against wsp reference DB
blastn -query "$ad/assembly.fasta" -db refs/wolbachia/wsp_db -outfmt 6 -max_target_seqs 1 > "$td/$sid.wsp.tsv" 2>/dev/null || true
mark_done "$ad/assembly"
