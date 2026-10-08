#!/usr/bin/env bash
#SBATCH -p cpu -c 16 --mem=32G -o logs/qc_%A_%a.out -e logs/qc_%A_%a.err
set -euo pipefail
# Under SLURM the script is copied to a spool dir, so $0 can't locate the repo;
# use the submit dir (=ROOT), with a non-SLURM fallback.
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}"
source "$HERE/workflow/lib/common.sh"
source "$HOME/opt/miniforge3/etc/profile.d/conda.sh" && conda activate wolb
i=${SLURM_ARRAY_TASK_ID:?}; T=$(cfg threads)
sid=$(sample_field "$i" 1); fq1=$(sample_field "$i" 6); fq2=$(sample_field "$i" 7)
out="$ROOT/results/02_nonhost"; qc="$ROOT/results/01_qc"; mkdir -p "$out" "$qc"
[[ -f $out/$sid.done ]] && { log "$sid done, skip"; exit 0; }
# 1) fastp
fastp -i "$fq1" -I "$fq2" -o "$qc/${sid}_clean_R1.fq.gz" -O "$qc/${sid}_clean_R2.fq.gz" \
  -j "$qc/$sid.fastp.json" -h "$qc/$sid.fastp.html" -w "$T" 2> "$ROOT/logs/$sid.fastp.log"
# 2) map to host; keep pairs both-unmapped (-f 12 -F 256) => non-host
minibwa map -x sr -t "$T" "$(cfg host_idx)" "$qc/${sid}_clean_R1.fq.gz" "$qc/${sid}_clean_R2.fq.gz" \
  | samtools view -@ "$T" -bf 12 -F 256 - \
  | samtools sort -@ "$T" -n -o "$out/$sid.nonhost.bam" -
samtools fastq -@ "$T" -1 "$out/${sid}_R1.fq.gz" -2 "$out/${sid}_R2.fq.gz" -0 /dev/null -s /dev/null "$out/$sid.nonhost.bam"
# guard: minibwa can exit 0 yet emit nothing (e.g. bad index) -> fail loud on empty output
need "$out/${sid}_R1.fq.gz" "$out/${sid}_R2.fq.gz"
# 3) stats
total=$(($(zcat "$qc/${sid}_clean_R1.fq.gz" | wc -l)/4*2))
nonhost=$(($(zcat "$out/${sid}_R1.fq.gz" | wc -l)/4*2))
hostm=$((total-nonhost)); pct=$(awk -v a=$hostm -v b=$total 'BEGIN{printf "%.4f",b?a/b:0}')
printf "sample\ttotal\thost_mapped\tnonhost\thost_pct\n%s\t%d\t%d\t%d\t%s\n" "$sid" "$total" "$hostm" "$nonhost" "$pct" > "$out/$sid.hoststat.txt"
rm -f "$qc/${sid}_clean_R1.fq.gz" "$qc/${sid}_clean_R2.fq.gz" "$out/$sid.nonhost.bam"
mark_done "$out/$sid"
