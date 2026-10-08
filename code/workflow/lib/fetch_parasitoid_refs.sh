#!/usr/bin/env bash
# Fetch parasitoid reference genome into $1 and build a blast db.
# Primary: Cotesia chilonis (WGS RJVT00000000; no standalone GCA -> resolve by taxon)
# Fallback: Cotesia congregata chromosome-level GCA_905319865.3
set -euo pipefail
OUT=$1; D=$(dirname "$OUT")
if datasets download genome taxon "Cotesia chilonis" --include genome --filename /tmp/cc.zip 2>/dev/null && \
   unzip -oq /tmp/cc.zip -d /tmp/cc 2>/dev/null && ls /tmp/cc/ncbi_dataset/data/*/*.fna >/dev/null 2>&1; then
  echo "using Cotesia chilonis assembly" >&2
else
  echo "fallback: Cotesia congregata GCA_905319865.3" >&2
  datasets download genome accession GCA_905319865.3 --include genome --filename /tmp/cc.zip
  rm -rf /tmp/cc; unzip -oq /tmp/cc.zip -d /tmp/cc
fi
cat /tmp/cc/ncbi_dataset/data/*/*.fna > "$OUT"
seqkit stats "$OUT"
rm -rf /tmp/cc /tmp/cc.zip
# COI barcode db (parasitoid species ID, Stage 6) — fetch separately when needed:
#   efetch -db nuccore -query '(Cotesia[Organism] OR Trichogramma[Organism]) AND COI' -format fasta > "$D/coi.fa"
#   makeblastdb -in "$D/coi.fa" -dbtype nucl -out "$D/coi_db"
makeblastdb -in "$OUT" -dbtype nucl -out "$D/parasitoid_refs" 2>/dev/null || true
