#!/usr/bin/env bash
# Build the combined dual-source mitochondrial reference for Stage 03c:
#   Chilo suppressalis mitogenome (NC_015612.1, RefSeq, verified via Entrez)
#   + the parasitoid mitogenomes already curated in refs/enemies/enemies_mito.fa
# Idempotent: skips work if the minibwa index already exists. Header convention
# is "Species|Accession" so Stage 03c can map idxstats contig -> species.
set -euo pipefail
source "$(dirname "$0")/common.sh"
source "$HOME/opt/miniforge3/etc/profile.d/conda.sh" && conda activate wolb

OUT="$ROOT/refs/mito_dualsource"; mkdir -p "$OUT"
CHILO_ACC="NC_015612.1"
chilo="$OUT/chilo_mito.fa"
combined="$OUT/mito_combined.fa"
prefix="$OUT/mito_combined"        # minibwa index prefix (.l2b/.mbw)

need "$ROOT/refs/enemies/enemies_mito.fa"

# 1. Chilo mitogenome (small RefSeq record; efetch is fine, no aria2c needed)
if [[ ! -s $chilo ]]; then
  log "fetching Chilo mito $CHILO_ACC via efetch"
  efetch -db nucleotide -id "$CHILO_ACC" -format fasta > "$chilo.raw" \
    || die "efetch failed for $CHILO_ACC (check network/proxy)"
  # Rename header to the Species|Accession convention used by enemies_mito.fa
  seqkit replace -p '^.*$' -r "Chilo_suppressalis|$CHILO_ACC" "$chilo.raw" > "$chilo"
  rm -f "$chilo.raw"
  need "$chilo"
  # Fail loud: edirect can return 0 despite a transient SSL/curl error, leaving a
  # truncated record. NC_015612.1 is 15,395 bp — guard against partial fetches.
  n=$(grep -c '^>' "$chilo"); L=$(seqkit fx2tab -nil "$chilo" | awk '{s+=$NF} END{print s}')
  [[ $n -eq 1 ]] || die "Chilo mito fetch malformed: $n seqs (expected 1) — rm $chilo and retry"
  [[ $L -ge 15000 && $L -le 16000 ]] || die "Chilo mito length ${L}bp outside expected ~15395 — truncated fetch; rm $chilo and retry"
fi

# 2. Combine Chilo + parasitoid mitos
if [[ ! -s $combined || $chilo -nt $combined || $ROOT/refs/enemies/enemies_mito.fa -nt $combined ]]; then
  log "building combined mito reference"
  cat "$chilo" "$ROOT/refs/enemies/enemies_mito.fa" > "$combined"
fi

# 3. minibwa index
if [[ ! -s $prefix.l2b || ! -s $prefix.mbw || $combined -nt $prefix.l2b ]]; then
  log "indexing combined mito reference with minibwa"
  minibwa index "$combined" "$prefix" || die "minibwa index failed"
fi
need "$prefix.l2b" "$prefix.mbw"
n=$(grep -c '^>' "$combined")
log "mito reference ready: $n contigs (1 Chilo + parasitoids) at $prefix"
