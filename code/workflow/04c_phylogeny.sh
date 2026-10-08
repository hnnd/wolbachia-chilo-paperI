#!/usr/bin/env bash
# Stage 04c (evidence line B3) — Wolbachia MLST strain phylogeny.
# Extracts the 5 pubMLST loci from our wolb+ sample assemblies, aligns them with the
# host-annotated reference panel (parasitoid wasps vs Lepidoptera vs flies; Stage
# fetch_wolb_phylo_refs.sh), and builds a partitioned ML tree. The discriminating
# read-out: do the Chilo Wolbachia sequences nest among PARASITOID Wolbachia (H1) or
# with the closely-related corn-borer / other Lepidoptera Wolbachia (H0)?
set -euo pipefail
source "$(dirname "$0")/lib/common.sh"
source "$HOME/opt/miniforge3/etc/profile.d/conda.sh" && conda activate wolb
cd "$ROOT"
T=$(cfg threads)
loci=(gatB coxA hcpA ftsZ fbpA)
mlstdb="$HOME/opt/miniforge3/envs/wolb/db/pubmlst/wolbachia"
ph="refs/wolb_phylo"; work="results/06_phylo"; mkdir -p "$work/loci"

# Gate: reference panel
if [[ ! -s $ph/loci/gatB.refs.fa ]]; then bash workflow/lib/fetch_wolb_phylo_refs.sh; fi
need "$ph/panel.tsv" "$ph/loci/gatB.refs.fa"

# 1. Extract each MLST locus from every wolb+ sample assembly (best BLAST hit region).
#    A sample is kept if it yields >=3 loci with length>=250 & pident>=90.
pos=results/05_typing/pos_samples.txt; need "$pos"
: > "$work/sample_loci_count.tsv"
declare -A SEQ   # SEQ[locus|sample]=seq
for s in $(cat "$pos"); do
  fa="results/04_assembly/$s/assembly.fasta"
  [[ -s $fa ]] || continue
  db="$work/.db_$s"; makeblastdb -in "$fa" -dbtype nucl -out "$db" >/dev/null 2>&1
  samtools faidx "$fa" 2>/dev/null || true
  ngood=0
  for L in "${loci[@]}"; do
    hit=$(blastn -query "$mlstdb/$L.tfa" -db "$db" \
            -outfmt "6 sseqid sstart send bitscore length pident" -max_target_seqs 1 2>/dev/null \
          | sort -k4,4nr | head -1 || true)
    [[ -z "$hit" ]] && continue
    read -r sid ss se bs ln pid <<<"$hit"
    awk_ok=$(awk -v l="$ln" -v p="$pid" 'BEGIN{print (l>=250 && p>=90)?1:0}')
    [[ "$awk_ok" == "1" ]] || continue
    if (( se >= ss )); then region="$sid:$ss-$se"; rcf=0; else region="$sid:$se-$ss"; rcf=1; fi
    seq=$(samtools faidx "$fa" "$region" 2>/dev/null | seqkit seq -s -w0 2>/dev/null || true)
    [[ $rcf -eq 1 ]] && seq=$(printf "%s\n" "$seq" | rev | tr 'ACGTacgtNn' 'TGCAtgcaNn')
    [[ -z "$seq" ]] && continue
    SEQ["$L|$s"]="$seq"; ngood=$((ngood+1))
  done
  printf "%s\t%s\n" "$s" "$ngood" >> "$work/sample_loci_count.tsv"
  rm -f "$db".*
done
keep=$(awk -F'\t' '$2>=3{print $1}' "$work/sample_loci_count.tsv")
nkeep=$(printf "%s\n" "$keep" | grep -c . || true)
log "samples with >=3 MLST loci recovered: $nkeep"
[[ $nkeep -ge 1 ]] || die "no sample recovered >=3 loci — phylogeny underpowered"

# 2. Per locus: refs + kept samples -> mafft align. Sanitize '|' to '_' for IQ-TREE.
alns=()
for L in "${loci[@]}"; do
  comb="$work/loci/$L.comb.fa"; sed 's/|/_/g' "$ph/loci/$L.refs.fa" > "$comb"
  for s in $keep; do
    [[ -n "${SEQ[$L|$s]:-}" ]] && printf ">Chilo_%s_CHILO\n%s\n" "$s" "${SEQ[$L|$s]}" >> "$comb"
  done
  mafft --auto --thread "$T" "$comb" > "$work/loci/$L.aln" 2>/dev/null
  alns+=("$work/loci/$L.aln")
done

# 3. Concatenate into partitioned supermatrix (missing loci -> gaps).
python3 workflow/lib/concat_aln.py "$work/supermatrix.fasta" "$work/partition.txt" "${alns[@]}"
need "$work/supermatrix.fasta" "$work/partition.txt"

# 4. ML tree (partitioned, ultrafast bootstrap + SH-aLRT), outgroup = Brugia (wBm).
og="wBm_NEMATODE"
( cd "$work" && rm -f supermatrix.fasta.* iqtree.* \
  && iqtree -s supermatrix.fasta -spp partition.txt -m TEST -bb 1000 -alrt 1000 \
       -nt AUTO -ntmax "$T" -o "$og" -pre iqtree -redo >/dev/null 2>&1 ) || \
  log "WARN: iqtree returned non-zero (check $work/iqtree.log)"
need "$work/iqtree.treefile"
log "tree: $work/iqtree.treefile"
conda run -n wolb-r Rscript workflow/04c_phylo_report.R
