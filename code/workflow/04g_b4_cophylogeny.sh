#!/bin/bash
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=60G
#SBATCH -J b4_cophylo
#SBATCH -o logs/b4_%j.out
#SBATCH -e logs/b4_%j.err
set -euo pipefail
# Stage 04g — B4 co-phylogeny (SLURM; mapping full non-host is too slow for login node).
#
# Tests co-transmission: does the parasitoid (Cotesia flavipes) mito haplotype track
# the Wolbachia strain composition across samples? Under parasitoid HT + vertical
# co-transmission of Wolbachia within the parasitoid, samples sharing a Cotesia mito
# haplotype should share the same Wolbachia (A:B) signature.
#
# Per high-load co-infection sample (breadth>=0.7):
#   1. Cotesia mito consensus  — map non-host → combined mito ref, consensus on the
#      Cotesia_flavipes contig (parasitoid haplotype).
#   2. Wolbachia A:B dosage     — map enriched wolb reads to wMel(A) & wPip(B), mean
#      depth on each; B-fraction = depthB/(depthA+depthB) (approx; conserved core
#      cross-maps, so this is relative not absolute).
# Then: align Cotesia consensuses, count variable sites (is there ANY haplotype
# variation to test?), NJ/ML haplotype tree, and tabulate haplotype vs B-fraction.
#
# Output: reports/stage04g_cophylogeny.md, results/07_cophylo/
HERE="${SLURM_SUBMIT_DIR:-/mnt/inaisfs/home/wangyunsheng/wangyunsheng/work/chilo/reseq/01.wolb}"
cd "$HERE"
source /mnt/inaisfs/home/wangyunsheng/wangyunsheng/miniforge3/etc/profile.d/conda.sh
conda activate wolb

T="${SLURM_CPUS_PER_TASK:-16}"
prof="results/03_profile/per_sample_profile.tsv"
pos="results/05_typing/pos_samples.txt"
od="results/07_cophylo"; mkdir -p "$od"
nh="results/02_nonhost"
mitoref="refs/mito_dualsource/mito_combined"
ctg="Cotesia_flavipes|NC_063945.1"
refA="refs/wolbachia/strains/NC_002978.6"   # wMel A
refB="refs/wolbachia/strains/NC_010981.1"   # wPip B
MINBREADTH=0.70

catfa="$od/cotesia_mito.all.fa"; : > "$catfa"
abtsv="$od/wolb_ab.tsv"; echo -e "sample\tdepthA\tdepthB\tBfrac\tcot_cov_bp" > "$abtsv"

meandepth() { # bam -> mean depth over covered+uncovered ref positions of a region ($2 optional)
  local bam="$1" reg="${2:-}"
  if [[ -n "$reg" ]]; then samtools depth -a -r "$reg" "$bam" 2>/dev/null; else samtools depth -a "$bam" 2>/dev/null; fi \
    | awk '{s+=$3; n++} END{printf "%.2f", (n?s/n:0)}'
}

while read -r s; do
  br=$(awk -F'\t' -v x="$s" '$1==x{print $3}' "$prof")
  awk -v b="${br:-0}" -v m="$MINBREADTH" 'BEGIN{exit !(b>=m)}' || continue
  r1="$nh/${s}_R1.fq.gz"; r2="$nh/${s}_R2.fq.gz"
  wr1="results/04_assembly/$s/w_R1.fq.gz"; wr2="results/04_assembly/$s/w_R2.fq.gz"
  [[ -s "$r1" && -s "$r2" && -s "$wr1" && -s "$wr2" ]] || continue
  echo "[$(date +%T)] $s : Cotesia mito consensus"
  cb="$od/.$s.mito.bam"
  minibwa map -x sr -t "$T" "$mitoref" "$r1" "$r2" 2>/dev/null | samtools sort -@ "$T" -o "$cb" - 2>/dev/null
  samtools index "$cb" 2>/dev/null || true
  cons="$od/$s.cot_mito.fa"
  samtools consensus -m simple -d 5 -c 0.6 --min-MQ 20 -r "$ctg" -o "$cons" "$cb" 2>/dev/null || : > "$cons"
  rm -f "$cb" "$cb".bai
  cov=$(seqkit seq -s -w0 "$cons" 2>/dev/null | tr -cd 'ACGTacgt' | wc -c)
  awk -v s="$s" '/^>/{print ">"s; next}{print}' "$cons" >> "$catfa"

  echo "[$(date +%T)] $s : Wolbachia A:B dosage"
  ba="$od/.$s.A.bam"; bb="$od/.$s.B.bam"
  minibwa map -x sr -t "$T" "$refA" "$wr1" "$wr2" 2>/dev/null | samtools sort -@ "$T" -o "$ba" - 2>/dev/null
  minibwa map -x sr -t "$T" "$refB" "$wr1" "$wr2" 2>/dev/null | samtools sort -@ "$T" -o "$bb" - 2>/dev/null
  dA=$(meandepth "$ba"); dB=$(meandepth "$bb"); rm -f "$ba" "$bb"
  bfrac=$(awk -v a="$dA" -v b="$dB" 'BEGIN{printf "%.3f", ((a+b)>0? b/(a+b):0)}')
  printf "%s\t%s\t%s\t%s\t%s\n" "$s" "$dA" "$dB" "$bfrac" "$cov" >> "$abtsv"
done < "$pos"

nseq=$(grep -c '^>' "$catfa" || true)
echo "[$(date +%T)] Cotesia mito consensuses: $nseq"

# align + variable sites + haplotype tree
aln="$od/cotesia_mito.aln"
mafft --auto --thread "$T" "$catfa" > "$aln" 2>/dev/null || cp "$catfa" "$aln"
( cd "$od" && rm -f cotesia_ml.* && iqtree -s "cotesia_mito.aln" -m TEST -bb 1000 -nt AUTO -ntmax "$T" -pre cotesia_ml -redo >/dev/null 2>&1 ) || echo "WARN iqtree"

conda run -n wolb-r Rscript workflow/04g_cophylo_report.R 2>&1 | tail -5 || \
  echo "WARN: report R failed (will still have raw tables)"
echo "[$(date +%T)] DONE. tables: $abtsv ; aln: $aln"
