#!/usr/bin/env bash
#SBATCH -p cpu -c 16 --mem=32G -o logs/asm_%A_%a.out -e logs/asm_%A_%a.err
# Stage 03f worker (T2) — de novo mitogenome assembly + BLAST ID for samples whose
# parasitoid could not be reliably named by best-hit (true species absent from the
# reference -> read scatter). Baits parasitoid-mito reads, assembles them de novo
# (so identity does NOT depend on the reference being the right species), then
# BLASTs the assembled mitogenome against a local arthropod-mito DB for the nearest
# relative + percent identity (confidence). Exact species for low-pident hits needs
# a follow-up remote nt BLAST on the assembled contig.
set -euo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}"
source "$HERE/workflow/lib/common.sh"
source "$HOME/opt/miniforge3/etc/profile.d/conda.sh" && conda activate wolb
i=${SLURM_ARRAY_TASK_ID:?}; T=$(cfg threads)
targets="$ROOT/results/06_denovo/assemble_targets.txt"; need "$targets"
sid=$(sed -n "${i}p" "$targets"); [[ -n $sid ]] || die "no sample at target line $i"
nh="$ROOT/results/02_nonhost"; out="$ROOT/results/06_denovo"; mkdir -p "$out"
r1="$nh/${sid}_R1.fq.gz"; r2="$nh/${sid}_R2.fq.gz"; need "$r1" "$r2"
[[ -f $out/$sid.done ]] && exit 0

ref="$ROOT/refs/parasitoid_broad/parasitoid_mito_broad"; need "$ref.l2b" "$ref.mbw"
db="$ROOT/refs/parasitoid_broad/mito_blastdb"
wd="$out/$sid.work"; rm -rf "$wd"; mkdir -p "$wd"

# 1. Bait: collect read names with ANY mate mapping to the parasitoid mito ref,
#    then pull the FULL pairs from the original fastqs. (Mito reads are mostly
#    single-mate hits against the small ref, so samtools fastq -s would discard
#    them — id-based pair extraction keeps both mates for assembly.)
minibwa map -x sr -t "$T" "$ref" "$r1" "$r2" \
  | samtools view -F 4 - | cut -f1 | sort -u > "$wd/ids.txt"
nb=$(wc -l < "$wd/ids.txt")
echo "[asm] $sid baited read-ids=$nb" >&2

# 2. Assemble baited pairs de novo (megahit: fast, low memory).
asm=""
if [[ $nb -ge 50 ]]; then
  seqkit grep -f "$wd/ids.txt" "$r1" -o "$wd/R1.fq" 2>/dev/null
  seqkit grep -f "$wd/ids.txt" "$r2" -o "$wd/R2.fq" 2>/dev/null
  megahit -1 "$wd/R1.fq" -2 "$wd/R2.fq" -t "$T" --min-contig-len 1000 \
    -o "$wd/megahit" >/dev/null 2>&1 || true
  [[ -s $wd/megahit/final.contigs.fa ]] && asm="$wd/megahit/final.contigs.fa"
fi

# 3. BLAST contigs vs local arthropod-mito DB; keep best hit per contig.
: > "$out/$sid.blast.tsv"
if [[ -n $asm ]]; then
  seqkit replace -p '\s.*$' -r '' "$asm" > "$out/$sid.contigs.fa"   # clean contig names
  blastn -query "$out/$sid.contigs.fa" -db "$db" -max_target_seqs 3 -evalue 1e-20 \
    -num_threads "$T" -outfmt "6 qseqid qlen sseqid pident length bitscore stitle" 2>/dev/null \
    | sort -k1,1 -k6,6nr | awk -F'\t' '!seen[$1]++' > "$out/$sid.blast.tsv" || true
fi

rm -rf "$wd"
mark_done "$out/$sid"
