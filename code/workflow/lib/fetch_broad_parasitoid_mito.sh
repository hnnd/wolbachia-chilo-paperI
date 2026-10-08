#!/usr/bin/env bash
# Build a BROAD parasitoid mitogenome reference for Stage 03e per-sample census.
# The para_reads library is Cotesia chilonis genome only, so para_frac is Cotesia-
# centric and misses other parasitoids (Pachyneuron, tachinids, ichneumonids...).
# Here we fetch complete mitogenomes across the parasitoid families relevant to
# rice stem borers, so reads can be assigned to the right genus/family.
#
# Strategy: take all mitogenomes from the smaller families; cap the large ones
# (Tachinidae, Encyrtidae) to limit redundancy/cross-mapping. Headers normalised
# to "Genus_species|accession". Chilo mito is included so host reads have a home.
# Idempotent: skips fetch if the per-family fasta already exists; rebuild index
# only when inputs change.
set -euo pipefail
source "$(dirname "$0")/common.sh"
source "$HOME/opt/miniforge3/etc/profile.d/conda.sh" && conda activate wolb

OUT="$ROOT/refs/parasitoid_broad"; mkdir -p "$OUT/fam"
combined="$OUT/parasitoid_mito_broad.fa"
prefix="$OUT/parasitoid_mito_broad"
CHILO="$ROOT/refs/mito_dualsource/chilo_mito.fa"
need "$CHILO"

# family : cap (0 = all available)
declare -A CAP=(
  [Braconidae]=0 [Ichneumonidae]=0 [Pteromalidae]=0 [Eulophidae]=0
  [Chalcididae]=0 [Trichogrammatidae]=0 [Eupelmidae]=0 [Bethylidae]=0
  [Tachinidae]=25 [Encyrtidae]=15
)

for fam in "${!CAP[@]}"; do
  fa="$OUT/fam/$fam.fa"
  if [[ ! -s $fa ]]; then
    log "fetching $fam complete mitogenomes"
    timeout 300 bash -c "esearch -db nuccore -query '${fam}[Organism] AND mitochondrion[Title] AND complete genome[Title]' \
      | efetch -format fasta" > "$fa.raw" 2>/dev/null || { log "WARN: $fam fetch failed, skipping"; rm -f "$fa.raw"; continue; }
    cap=${CAP[$fam]}
    if [[ $cap -gt 0 ]]; then seqkit head -n "$cap" "$fa.raw" > "$fa"; else mv "$fa.raw" "$fa"; fi
    rm -f "$fa.raw"
  fi
  printf "  %-18s %s seqs\n" "$fam" "$(grep -c '^>' "$fa" 2>/dev/null || echo 0)"
done

# Combine + normalise headers to Genus_species|accession (parse "ACC Genus species ...").
# Strip qualifier prefixes (UNVERIFIED:, TPA_asm:, etc.) that sit before the genus.
# Drop sequences shorter than 12 kb (partial) and dedup identical sequences.
log "combining + normalising headers"
cat "$OUT"/fam/*.fa > "$OUT/all_raw.fa"
seqkit seq -m 12000 "$OUT/all_raw.fa" \
  | seqkit replace -p '^(\S+)\s+(UNVERIFIED:|TPA_asm:|TPA:|MAG:)\s+' -r '${1} ' \
  | seqkit replace -p '^(\S+)\s+(\S+)\s+(\S+).*$' -r '${2}_${3}|${1}' \
  > "$OUT/fam_renamed.fa"
# UNION with the curated enemies_mito.fa: it carries empirically-relevant taxa
# (Pachyneuron, Lypha, Trichomalopsis...) whose NCBI titles lack "complete genome"
# and are therefore missed by the family query. Already Species|Accession form.
need "$ROOT/refs/enemies/enemies_mito.fa"
cat "$ROOT/refs/enemies/enemies_mito.fa" "$OUT/fam_renamed.fa" \
  | seqkit rmdup -s 2>/dev/null > "$OUT/parasitoid_only.fa"
# Prepend Chilo mito (host control) — already in Species|Accession form.
cat "$CHILO" "$OUT/parasitoid_only.fa" > "$combined"
rm -f "$OUT/all_raw.fa" "$OUT/fam_renamed.fa"

n=$(grep -c '^>' "$combined")
[[ $n -ge 50 ]] || die "broad reference too small ($n seqs) — fetch likely failed"
log "broad reference: $n contigs (1 Chilo + $((n-1)) parasitoid mitos)"

# minibwa index
if [[ ! -s $prefix.l2b || ! -s $prefix.mbw || $combined -nt $prefix.l2b ]]; then
  log "indexing broad reference"
  minibwa index "$combined" "$prefix" || die "minibwa index failed"
fi
need "$prefix.l2b" "$prefix.mbw"
log "broad parasitoid mito reference ready at $prefix"
