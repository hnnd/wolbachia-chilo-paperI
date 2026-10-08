#!/usr/bin/env bash
# Stage 04h (B5) — objective same-strain criterion (login node, no SLURM):
#   (1) MLST allele-profile concordance from results/05_typing/*.mlst.tsv
#   (2) whole-genome fastANI among contiguous assemblies + reference strains
# fastANI runs in the `wolb` env; the report/figures in `wolb-r`.
set -euo pipefail
source "$(dirname "$0")/lib/common.sh"
source "$HOME/opt/miniforge3/etc/profile.d/conda.sh"
cd "$ROOT"

WORK=results/05_typing/strain_id; mkdir -p "$WORK"

# reference strains: wMel(A) wRi(A) wPip(B) wTpre(B) wBm(D outgroup)
: > "$WORK/ref_list.txt"
for r in refs/wolbachia/strains/NC_002978.6.fa refs/wolbachia/strains/NC_012416.1.fa \
         refs/wolbachia/strains/NC_010981.1.fa refs/wolbachia/strains/NZ_CM003641.1.fa \
         refs/wolb_phylo/genomes/wBm.fna; do
  [[ -s $r ]] && echo "$r" >> "$WORK/ref_list.txt"
done
need "$WORK/ref_list.txt"

# query: Chilo Wolbachia assemblies >= 500 kb (small ones are too fragmented for ANI)
: > "$WORK/query_list.txt"
for d in results/04_assembly/*/; do
  f="$d/assembly.fasta"; [[ -s $f ]] || continue
  bp=$(grep -v '^>' "$f" | tr -d '\n' | wc -c)
  (( bp >= 500000 )) && echo "$f" >> "$WORK/query_list.txt"
done
need "$WORK/query_list.txt"
cat "$WORK/query_list.txt" "$WORK/ref_list.txt" > "$WORK/all_list.txt"

# fastANI all-vs-all (idempotent: skip if already computed)
if [[ ! -s $WORK/ani_all.tsv ]]; then
  conda activate wolb
  command -v fastANI >/dev/null || die "fastANI not found in env wolb (Stage 0a: mamba install -n wolb -c bioconda fastani)"
  fastANI --ql "$WORK/all_list.txt" --rl "$WORK/all_list.txt" -o "$WORK/ani_all.tsv" -t 8
fi
need "$WORK/ani_all.tsv"

conda activate wolb-r   # R lives in wolb-r env
Rscript workflow/04h_strain_identity.R
need reports/stage04h_strain_identity.md \
     results/05_typing/strain_id/mlst_profiles.tsv \
     results/05_typing/strain_id/mlst_recurrence.tsv \
     reports/figs/strain_mlst_concordance.png reports/figs/strain_ani_heatmap.png
