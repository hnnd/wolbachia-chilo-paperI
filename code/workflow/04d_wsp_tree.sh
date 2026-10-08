#!/usr/bin/env bash
# Stage 04d (evidence line B3, decisive test) — Wolbachia wsp gene tree INCLUDING
# Cotesia flavipes Wolbachia (the dominant Chilo parasitoid, whose Wolbachia has no
# genome but does have wsp on GenBank). The MLST tree (04c) lacked Cotesia and was
# inconclusive; wsp lets us ask directly: does the Chilo Wolbachia wsp group with
# Cotesia flavipes wsp (H1, parasitoid-origin) or elsewhere?
# Caveat: wsp is single-locus and recombination-prone — read alongside MLST (04c).
set -euo pipefail
source "$(dirname "$0")/lib/common.sh"
source "$HOME/opt/miniforge3/etc/profile.d/conda.sh" && conda activate wolb
cd "$ROOT"
T=$(cfg threads)
ph="refs/wolb_phylo"; work="results/06_phylo"; mkdir -p "$work"
need "$ph/panel.tsv" "$ph/cotesia/cotesia_wsp.fa"

# wsp query for BLAST extraction (Wolbachia WalbB wsp; conserved enough to seed any strain).
qf="$ph/wsp_query.fa"
if [[ ! -s $qf ]]; then efetch -db nuccore -id AF020059.1 -format fasta 2>/dev/null > "$qf"; fi
need "$qf"

# extract_wsp <subject.fna> <label>  -> appends a wsp record to $out (empty if no hit)
out="$work/wsp.comb.fa"; : > "$out"
extract_wsp() {
  local subj="$1" label="$2" db hit sid ss se bs ln pid region seq rcf
  [[ -s $subj ]] || return 0
  db="$work/.wdb_$$"; makeblastdb -in "$subj" -dbtype nucl -out "$db" >/dev/null 2>&1 || return 0
  samtools faidx "$subj" 2>/dev/null || true
  hit=$(blastn -query "$qf" -db "$db" -outfmt "6 sseqid sstart send bitscore length pident" \
          -max_target_seqs 1 2>/dev/null | sort -k4,4nr | head -1 || true)
  rm -f "$db".*
  [[ -z "$hit" ]] && return 0
  read -r sid ss se bs ln pid <<<"$hit"
  if (( se >= ss )); then region="$sid:$ss-$se"; rcf=0; else region="$sid:$se-$ss"; rcf=1; fi
  seq=$(samtools faidx "$subj" "$region" 2>/dev/null | seqkit seq -s -w0 2>/dev/null || true)
  [[ $rcf -eq 1 ]] && seq=$(printf "%s\n" "$seq" | rev | tr 'ACGTacgtNn' 'TGCAtgcaNn')
  [[ -z "$seq" ]] && return 0
  printf ">%s\n%s\n" "$label" "$seq" >> "$out"
}

# 1. Cotesia wsp (GenBank) — sanitize '|' for IQ-TREE.
sed 's/|/_/g' "$ph/cotesia/cotesia_wsp.fa" >> "$out"
# 2. wsp from each panel reference genome.
while IFS=$'\t' read -r strain host htype sg acc; do
  extract_wsp "$ph/genomes/$strain.fna" "${strain}_${htype}"
done < <(tail -n +2 "$ph/panel.tsv")
# 3. wsp from each wolb+ Chilo assembly.
for s in $(cat results/05_typing/pos_samples.txt); do
  extract_wsp "results/04_assembly/$s/assembly.fasta" "Chilo_${s}_CHILO"
done

n=$(grep -c '^>' "$out"); log "wsp records: $n"; need "$out"
nchilo=$(grep -c '^>Chilo_' "$out" || true); log "Chilo wsp recovered: $nchilo"
[[ $nchilo -ge 1 ]] || die "no Chilo sample yielded wsp — cannot test"

# 4. Align + ML tree (midpoint-rooted in report; wsp too short/recombinant for fixed outgroup).
mafft --auto --thread "$T" "$out" > "$work/wsp.aln" 2>/dev/null
( cd "$work" && rm -f wsp.aln.* wsp_iqtree.* \
  && iqtree -s wsp.aln -m TEST -bb 1000 -alrt 1000 -nt AUTO -ntmax "$T" -pre wsp_iqtree -redo >/dev/null 2>&1 ) || \
  log "WARN: wsp iqtree non-zero"
need "$work/wsp_iqtree.treefile"
conda run -n wolb-r Rscript workflow/04d_wsp_report.R
