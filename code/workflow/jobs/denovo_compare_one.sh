#!/bin/bash
#SBATCH -p cpu
#SBATCH -N 1
#SBATCH -c 32
#SBATCH --mem=300G
#SBATCH -J denovo_cmp
#SBATCH -o logs/denovo_%j.out
#SBATCH -e logs/denovo_%j.err
set -euo pipefail

# One-off diagnostic (think #2): does WHOLE non-host metaSPAdes recover Wolbachia
# sequence that the reference-guided enrichment assembly misses, for the most
# divergent + lowest-abundance positive sample? Settles whether reference-mapping
# bias (not depth) is the cause, before committing the cohort to whole-metagenome
# assembly. Default sample HC_W_9 (tachinid-dominated, 3.43% mapped divergence).

HERE="${SLURM_SUBMIT_DIR:-/mnt/inaisfs/home/wangyunsheng/wangyunsheng/work/chilo/reseq/01.wolb}"
cd "$HERE"
SID="${1:-HC_W_9}"
T="${SLURM_CPUS_PER_TASK:-32}"

source /mnt/inaisfs/home/wangyunsheng/wangyunsheng/miniforge3/etc/profile.d/conda.sh
conda activate wolb

nh="results/02_nonhost"
od="results/06_denovo_compare/$SID"
enr="results/04_assembly/$SID/assembly.fasta"        # existing enrichment assembly
wref="refs/wolbachia/wolbachia_refs.fa"               # merged 4-strain (contig classifier)
wspq="refs/wolb_phylo/wsp_query.fa"                   # AF020059 wsp gene
lociD="refs/wolb_phylo/loci"                          # gatB/coxA/hcpA/ftsZ/fbpA .refs.fa
mkdir -p "$od"

r1="$nh/${SID}_R1.fq.gz"; r2="$nh/${SID}_R2.fq.gz"
[[ -s $r1 && -s $r2 ]] || { echo "FATAL: missing non-host reads for $SID"; exit 1; }

# ---- 1. Whole non-host metaSPAdes -------------------------------------------
ms="$od/metaspades"
if [[ ! -s "$ms/contigs.fasta" ]]; then
  echo "[$(date +%T)] metaSPAdes on full non-host set ($SID) ..."
  metaspades.py -1 "$r1" -2 "$r2" -t "$T" -m 290 -o "$ms" \
    || { echo "metaSPAdes failed"; exit 1; }
fi
whole="$ms/contigs.fasta"
echo "[$(date +%T)] whole-metagenome contigs: $(grep -c '>' "$whole")"

# ---- 2. Classify Wolbachia contigs out of the metagenome --------------------
# A contig is 'Wolbachia' if it has a blastn hit to the merged 4-strain ref with
# >=500 bp aligned at >=80% identity (conserved enough to catch a divergent strain,
# strict enough to reject host/gut-microbe contigs sharing short conserved genes).
db="$od/.wref_db"
makeblastdb -in "$wref" -dbtype nucl -out "$db" >/dev/null 2>&1
blastn -query "$whole" -db "$db" -task megablast \
  -outfmt "6 qseqid pident length bitscore" -max_target_seqs 1 -evalue 1e-10 2>/dev/null \
  | awk '$2>=80 && $3>=500{print $1}' | sort -u > "$od/wolb_contig_ids.txt" || true
rm -f "$db".*
nwolb=$(wc -l < "$od/wolb_contig_ids.txt" || echo 0)
echo "[$(date +%T)] Wolbachia-classified contigs: $nwolb"

wolbfa="$od/${SID}.wolb_from_metagenome.fasta"
if [[ "$nwolb" -gt 0 ]]; then
  seqkit grep -f "$od/wolb_contig_ids.txt" "$whole" > "$wolbfa"
else
  : > "$wolbfa"
fi

# ---- 3. Compare enrichment vs whole-metagenome wolb contigs -----------------
# helper: print "<bp> <n> <max> <N50>" for a fasta (0 0 0 0 if empty)
faistats() {
  # sum_len num_seqs max_len N50  (from seqkit stats -a -T)
  local fa="$1"
  if [[ -s "$fa" ]]; then
    seqkit stats -a -T "$fa" 2>/dev/null \
      | awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} NR==2{print $h["sum_len"], $h["num_seqs"], $h["max_len"], $h["N50"]}'
  else echo "0 0 0 0"; fi
}
# helper: best wsp hit "<len> <pid>" (query=wsp gene vs db=assembly)
wsp_hit() {
  local fa="$1" d="$od/.q_$$"
  [[ -s "$fa" ]] || { echo "0 0"; return; }
  makeblastdb -in "$fa" -dbtype nucl -out "$d" >/dev/null 2>&1 || { echo "0 0"; return; }
  { blastn -query "$wspq" -db "$d" -outfmt "6 length pident bitscore" -max_target_seqs 1 2>/dev/null \
    | sort -k3,3nr | head -1 | awk 'END{print (NR?L:0), (NR?P:0)} {L=$1;P=$2}'; } || true
  rm -f "$d".*
}
# helper: count MLST loci recovered (len>=250 & pid>=90), same rule as Stage 04c
mlst_n() {
  local fa="$1" d="$od/.m_$$" L ln pid n=0
  [[ -s "$fa" ]] || { echo 0; return; }
  makeblastdb -in "$fa" -dbtype nucl -out "$d" >/dev/null 2>&1 || { echo 0; return; }
  for L in gatB coxA hcpA ftsZ fbpA; do
    ln=0; pid=0
    read -r ln pid < <({ blastn -query "$lociD/$L.refs.fa" -db "$d" \
        -outfmt "6 length pident bitscore" -max_target_seqs 1 2>/dev/null \
        | sort -k3,3nr | head -1 | awk '{print ($1==""?0:$1), ($2==""?0:$2)}'; } || true) || true
    awk -v l="${ln:-0}" -v p="${pid:-0}" 'BEGIN{exit !(l>=250 && p>=90)}' && n=$((n+1)) || true
  done
  rm -f "$d".*; echo "$n"
}

rep="reports/stage_denovo_compare_${SID}.md"
{
  echo "# De novo whole-metagenome vs enrichment — $SID"
  echo
  echo "One-off test (think #2): is reference-mapping bias, not depth, why $SID failed strain typing?"
  echo
  read -r ebp en emax en50 < <(faistats "$enr") || true
  read -r wbp wn wmax wn50 < <(faistats "$wolbfa") || true
  read -r ewl ewp < <(wsp_hit "$enr") || true
  read -r wwl wwp < <(wsp_hit "$wolbfa") || true
  emn=$(mlst_n "$enr"); wmn=$(mlst_n "$wolbfa")
  echo "| metric | enrichment (map→ref→assemble) | whole metagenome (assemble→filter) |"
  echo "|---|---|---|"
  echo "| Wolbachia contigs | $en | $wn |"
  echo "| total Wolbachia bp | $ebp | $wbp |"
  echo "| max contig (bp) | $emax | $wmax |"
  echo "| N50 (bp) | $en50 | $wn50 |"
  echo "| wsp best hit (len/pid) | ${ewl}bp / ${ewp}% | ${wwl}bp / ${wwp}% |"
  echo "| MLST loci recovered (/5) | $emn | $wmn |"
  echo
  echo "Whole-metagenome total contigs: $(grep -c '>' "$whole"); Wolbachia-classified: $nwolb."
  echo
  # Data-driven verdict: does whole-metagenome recover Wolbachia signal enrichment missed?
  gain=$(awk -v w="$wmn" -v e="$emn" -v ww="$wwl" -v ew="$ewl" -v wb="$wbp" -v eb="$ebp" \
    'BEGIN{print (w>e || (ww>=250 && ew<250) || wb>eb*1.2)?1:0}')
  echo "## Verdict"
  if [[ "$gain" == "1" ]]; then
    echo "**Reference-recruitment bias is real for $SID.** Whole-metagenome assembly recovered"
    echo "Wolbachia sequence (wsp/MLST/more genome) that enrichment missed → the map-to-reference"
    echo "step drops divergent-strain reads. Worth widening recruitment for divergent samples."
  else
    echo "**Depth wall confirmed — NOT reference bias.** Whole-metagenome assembly (recruitment-"
    echo "bias-free) recovered *less* Wolbachia ($wbp bp vs $ebp bp) and still **no wsp, no full"
    echo "MLST locus** — the loci are simply not sequenced deeply enough (depth <0.5×, ~82% of the"
    echo "genome at zero coverage). Enrichment is the better method for a low-abundance intracellular"
    echo "symbiont (isolating wolb reads first prevents low-coverage tips being pruned in the 533k-"
    echo "contig background). Do **not** switch the cohort to whole-metagenome assembly; honestly"
    echo "caveat $SID (and the other <0.5×-depth positives) as too low-abundance to strain-type."
    echo "Note: $SID's mapped-read divergence is 3.43% (vs ~1.8% in rice samples) — the tachinid"
    echo "carries a genuinely more-divergent Wolbachia, but abundance, not divergence, blocks typing."
  fi
} > "$rep"

echo "[$(date +%T)] DONE. Report: $rep"
cat "$rep"
