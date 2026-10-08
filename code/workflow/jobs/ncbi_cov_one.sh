#!/usr/bin/env bash
#SBATCH -p cpu -c 40 --mem=64G -o logs/ncbi_cov_%j.out -e logs/ncbi_cov_%j.err
#
# Stage 17 worker — measured PacBio HiFi base coverage for the three NCBI genome
# submissions (wCchiA, wCchiB, wCchiA_unplaced).
#
# Maps the full CS_R HiFi read set (the metagenome that yielded both Wolbachia
# strains) against the exact FASTA files that will be deposited, and writes the
# per-contig and per-assembly mean depth.  The per-assembly mean is what goes
# into the mandatory "Genome coverage" field of the NCBI Genome Info table.
#
# Same pipeline run on ta (used 2026-09-23 while the cpu partition was saturated,
# reads streamed from hpc-c; identical minimap2/samtools flags):
#   ssh hpc-c "cat <hifi_cs_r>" | pigz -dc -p 6 \
#     | minimap2 -ax map-hifi -t 40 --sam-hit-only ncbi_submission_ref.fa - \
#     | samtools sort -@ 8 -m 1G -o hifi_vs_submission.bam -
set -euo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}"   # SLURM copies script to spool; locate repo via submit dir
source "$HERE/workflow/lib/common.sh"
source "$HOME/opt/miniforge3/etc/profile.d/conda.sh" && conda activate wolb

fq=$(cfg hifi_cs_r)
out="$ROOT/results/17_ncbi_cov"
ref="$out/ncbi_submission_ref.fa"
need "$fq" "$out/wCchiA.fsa" "$out/wCchiB.fsa" "$out/wCchiA_unplaced.fsa"
mkdir -p "$out"
cat "$out/wCchiA.fsa" "$out/wCchiB.fsa" "$out/wCchiA_unplaced.fsa" > "$ref"

# --sam-hit-only drops the ~99.9 % of reads that do not touch the 2.8 Mb reference,
# keeping the BAM (and the sort step) small; pigz decompresses the 22 GB fastq in
# parallel while minimap2 maps.
pigz -dc -p 6 "$fq" \
  | minimap2 -ax map-hifi -t 40 --sam-hit-only "$ref" - \
  | samtools sort -@ 8 -m 1G -o "$out/hifi_vs_submission.bam" -
samtools index -@ 4 "$out/hifi_vs_submission.bam"
mapped=$(samtools view -c -F 0x904 "$out/hifi_vs_submission.bam")

# Per-base depth over every position of every deposited contig (including zeros).
samtools depth -a -d 0 "$out/hifi_vs_submission.bam" > "$out/per_base_depth.tsv"

# Per-contig mean, then length-weighted per-assembly mean (= total aligned base /
# total bp).  Assembly labels follow the SeqID prefixes of the submission FASTAs.
awk -v OFS='\t' '
  { n[$1]++; s[$1]+=$3; if($3>0) c[$1]++ }
  END { print "seqid","length","mean_depth","frac_covered"
        for (k in n) printf "%s\t%d\t%.3f\t%.4f\n", k, n[k], s[k]/n[k], c[k]/n[k] }
' "$out/per_base_depth.tsv" | sort -k1,1 > "$out/per_contig_depth.tsv"

awk -v OFS='\t' -v mapped="$mapped" '
  NR==1 { next }
  { a = ($1 ~ /unplaced/) ? "wCchiA_unplaced" : ($1 ~ /^wCchiB/) ? "wCchiB" : "wCchiA"
    nc[a]++; bp[a]+=$2; bases[a]+=$2*$3 }
  END { print "assembly","n_contigs","total_bp","mean_depth","mapped_read_pairs"
        for (k in bp) printf "%s\t%d\t%d\t%.2f\t%s\n", k, nc[k], bp[k], bases[k]/bp[k], mapped }
' "$out/per_contig_depth.tsv" | sort -k1,1 > "$out/coverage_summary.tsv"

cat "$out/coverage_summary.tsv"
mark_done "$out/coverage"
