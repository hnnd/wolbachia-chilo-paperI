#!/usr/bin/env bash
# 从 fq 目录 + info.txt 生成 samples.tsv
set -euo pipefail
cd "$(dirname "$0")/.."
ROOT=$(pwd)
FQ_LAB=$(readlink -f ../fq/lab_seq) ; FQ_SRA=$(readlink -f ../fq/sra) ; INFO=$(readlink -f ../info.txt)
OUT=config/samples.tsv
printf "sample_id\tsource\tgroup\thost_plant\tgeo\tfq1\tfq2\n" > "$OUT"
# lab：<prefix>_<n>_1.fq.gz；prefix 末位 R/W = 水稻/茭草，group 取 prefix
for f1 in "$FQ_LAB"/*_1.fq.gz; do
  b=$(basename "$f1" _1.fq.gz)              # e.g. AQ_R_1
  f2="$FQ_LAB/${b}_2.fq.gz"
  [[ -f $f2 ]] || { echo "SKIP (missing pair): $f2" >&2; continue; }
  grp=$(echo "$b" | sed -E 's/_[0-9]+$//')  # AQ_R
  hp=W; [[ $grp == *_R ]] && hp=rice; [[ $grp == *_W ]] && hp=wateroat
  geo=$(echo "$grp" | sed -E 's/_[RW]$//')  # AQ
  printf "%s\tlab\t%s\t%s\t%s\t%s\t%s\n" "$b" "$grp" "$hp" "$geo" "$f1" "$f2" >> "$OUT"
done
# sra：<BioSample>_R1.fastq.gz；host/geo 从 info.txt 的 BioSample 行查（列号见 info.txt 表头）
# col 6 = BioSample, col 30 = geo_loc_name (verified from header)
for f1 in "$FQ_SRA"/*_R1.fastq.gz; do
  bs=$(basename "$f1" _R1.fastq.gz)
  f2="$FQ_SRA/${bs}_R2.fastq.gz"
  [[ -f $f2 ]] || { echo "missing pair: $f2" >&2; exit 1; }
  geo=$(awk -F'\t' -v b="$bs" '$6==b{print $30; exit}' "$INFO"); geo=${geo:-NA}
  printf "%s\tsra\tSRA\trice\t%s\t%s\t%s\n" "$bs" "$geo" "$f1" "$f2" >> "$OUT"
done
echo "samples: $(($(wc -l < "$OUT")-1))"
