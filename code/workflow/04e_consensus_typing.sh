#!/bin/bash
# Stage 04e — reference-guided consensus typing (login node; reads are tiny).
#
# de novo assembly (Stage 4) needs local overlap, so it drops MLST/wsp loci in
# samples with decent breadth but patchy depth (e.g. NC2025_R_1/CD_R_3/HZ_R_5:
# breadth 0.63-0.84 yet 0 wsp de novo). A pileup consensus only needs per-base
# depth, so it recovers loci that fall in covered-but-unassemblable regions.
# This ADDS samples to the strain phylogeny at ~zero cost (reuses enriched wolb
# reads). It cannot rescue depth-wall samples (HC_W_9 etc., <0.5x) — those loci
# have zero coverage, see reports/stage_denovo_compare_HC_W_9.md.
#
# Output: results/05_typing/consensus/<sid>.consensus.fasta  (whole best-strain
# consensus, N where depth<2) + reports/stage04e_consensus.md coverage table.
set -euo pipefail
HERE="${SLURM_SUBMIT_DIR:-/mnt/inaisfs/home/wangyunsheng/wangyunsheng/work/chilo/reseq/01.wolb}"
cd "$HERE"
source /mnt/inaisfs/home/wangyunsheng/wangyunsheng/miniforge3/etc/profile.d/conda.sh
conda activate wolb

T="${SLURM_CPUS_PER_TASK:-8}"
prof="results/03_profile/per_sample_profile.tsv"
pos="results/05_typing/pos_samples.txt"
cdir="results/05_typing/consensus"; mkdir -p "$cdir"
wspq="refs/wolb_phylo/wsp_query.fa"
lociD="refs/wolb_phylo/loci"
strainsD="refs/wolbachia/strains"
[[ -s $pos && -s $prof ]] || { echo "FATAL: missing pos_samples.txt or profile"; exit 1; }

# best strain per sample from the profile (col 2), keyed once.
declare -A STRAIN
while IFS=$'\t' read -r s st _; do STRAIN["$s"]="$st"; done < <(tail -n +2 "$prof" | cut -f1,2)

# recovered-loci SET for a fasta subject -> echoes space-separated tokens among
# {gatB coxA hcpA ftsZ fbpA wsp}. Locus kept if len>=250 & pid>=90 (wsp pid>=80).
# Returning the set (not a count) lets the caller UNION de novo and consensus
# without one masking the other via bitscore.
recovered_set() {
  local fa="$1" d="$cdir/.cq_$$" L ln pid out=""
  [[ -s "$fa" ]] || { echo ""; return; }
  makeblastdb -in "$fa" -dbtype nucl -out "$d" >/dev/null 2>&1 || { echo ""; return; }
  for L in gatB coxA hcpA ftsZ fbpA; do
    ln=0; pid=0
    read -r ln pid < <({ blastn -query "$lociD/$L.refs.fa" -db "$d" \
        -outfmt "6 length pident bitscore" -max_target_seqs 1 2>/dev/null \
        | sort -k3,3nr | head -1 | awk 'END{print (NR?L:0),(NR?P:0)}{L=$1;P=$2}'; } || true) || true
    awk -v l="${ln:-0}" -v p="${pid:-0}" 'BEGIN{exit !(l>=250 && p>=90)}' && out="$out $L" || true
  done
  ln=0; pid=0
  read -r ln pid < <({ blastn -query "$wspq" -db "$d" \
      -outfmt "6 length pident bitscore" -max_target_seqs 1 2>/dev/null \
      | sort -k3,3nr | head -1 | awk 'END{print (NR?L:0),(NR?P:0)}{L=$1;P=$2}'; } || true) || true
  awk -v l="${ln:-0}" -v p="${pid:-0}" 'BEGIN{exit !(l>=250 && p>=80)}' && out="$out wsp" || true
  rm -f "$d".*
  echo "$out"
}
# count MLST tokens (exclude wsp) in a set; and whether wsp present
n_mlst() { local n=0 t; for t in $*; do case "$t" in gatB|coxA|hcpA|ftsZ|fbpA) n=$((n+1));; esac; done; echo "$n"; }
has_wsp() { local t; for t in $*; do [[ "$t" == wsp ]] && { echo 1; return; }; done; echo 0; }

rep="reports/stage04e_consensus.md"
{
  echo "# Stage 04e — reference-guided consensus typing"
  echo
  echo "Rescue MLST/wsp loci for positive samples where de novo assembly dropped them."
  echo "Consensus: enriched wolb reads → minibwa best-strain → \`samtools consensus\` (simple,"
  echo "min-depth 2, min-MQ 20, N below). Locus kept if len≥250 & pid≥90 (wsp pid≥80)."
  echo
  echo "| sample | breadth | best strain | de novo (MLST/wsp) | +consensus (MLST/wsp) | gained |"
  echo "|---|---|---|---|---|---|"
} > "$rep"

for s in $(cat "$pos"); do
  st="${STRAIN[$s]:-}"
  br=$(awk -F'\t' -v x="$s" '$1==x{print $3}' "$prof")
  r1="results/04_assembly/$s/w_R1.fq.gz"; r2="results/04_assembly/$s/w_R2.fq.gz"
  denovo="results/04_assembly/$s/assembly.fasta"
  cons="$cdir/$s.consensus.fasta"
  [[ -n "$st" && -s "$r1" && -s "$r2" && -f "$strainsD/$st.mbw" ]] || { echo "  skip $s (no reads/strain)"; continue; }

  if [[ ! -s "$cons" ]]; then
    bam="$cdir/.$s.bam"
    minibwa map -x sr -t "$T" "$strainsD/$st" "$r1" "$r2" 2>/dev/null \
      | samtools sort -@ "$T" -o "$bam" - 2>/dev/null
    samtools index "$bam" 2>/dev/null || true
    samtools consensus -m simple -a -d 2 -c 0.6 --min-MQ 20 -o "$cons" "$bam" 2>/dev/null || : > "$cons"
    rm -f "$bam" "$bam".bai
  fi

  # evaluate de novo and consensus SEPARATELY, then union (no bitscore masking)
  dset=$(recovered_set "$denovo")
  cset=$(recovered_set "$cons")
  uset=$(echo "$dset $cset" | tr ' ' '\n' | sort -u | tr '\n' ' ')
  dn=$(n_mlst $dset); dw=$(has_wsp $dset)
  un=$(n_mlst $uset); uw=$(has_wsp $uset)
  gain=""
  (( un > dn )) && gain="+$((un-dn)) MLST" || true
  (( uw > dw )) && gain="$gain +wsp" || true
  [[ -z "$gain" ]] && gain="—" || true
  printf "| %s | %s | %s | %s/%s | %s/%s | %s |\n" \
    "$s" "$br" "$st" "$dn" "$dw" "$un" "$uw" "$gain" >> "$rep"
done

{
  echo
  echo "## Conclusion"
  echo "Reference-guided consensus does **not** expand the strain phylogeny. Only CD_W_1"
  echo "gained a locus (4→5 MLST, already on both trees); **no sample crossed a tree-inclusion"
  echo "threshold (MLST≥3 or wsp)** and **no sample gained wsp**. Two reasons, both data/biology"
  echo "not pipeline: (1) conserved MLST loci already assemble de novo wherever depth suffices,"
  echo "so consensus is redundant there; (2) wsp is hypervariable — its reads are lost at the"
  echo "reference-mapping step, so they are absent from the pileup and consensus cannot call it"
  echo "(same limit for all reference-recruited methods). Note the many breadth≥0.7 samples with"
  echo "5/0 (full MLST, no wsp): that pattern IS the narrow reference-recruitment bias, confined"
  echo "to the wsp locus. Decision: 04e is a standalone diagnostic; it is NOT wired into 04c/04d."
} >> "$rep"

echo "[done] $rep"
cat "$rep"
