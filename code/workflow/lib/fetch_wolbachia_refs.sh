#!/usr/bin/env bash
# Fetch representative Wolbachia reference genomes (verified accessions) and
# concatenate into $1; also build a wsp blast db next to it.
set -euo pipefail
OUT=$1
# Verified accessions (2026-06-28): A=wMel,wRi ; B=wPip ; parasitoid-source key control = wTpre
# wVitA/wVitB (Nasonia) can be appended later (one A + one B) once accessions confirmed.
ACC=(GCF_000008025.1 GCF_000022285.1 GCF_000073005.1 GCF_001439985.1)  # wMel,wRi,wPip,wTpre
: > "$OUT"
for a in "${ACC[@]}"; do
  datasets download genome accession "$a" --include genome --filename "/tmp/$a.zip"
  unzip -oq "/tmp/$a.zip" -d "/tmp/$a"
  cat /tmp/$a/ncbi_dataset/data/*/*.fna >> "$OUT"
  rm -rf "/tmp/$a" "/tmp/$a.zip"
done
seqkit stats "$OUT"
# wsp blast db (used in Stage 4 for wsp typing); built from the same refs for now.
makeblastdb -in "$OUT" -dbtype nucl -out "$(dirname "$OUT")/wsp_db" 2>/dev/null || true
