#!/usr/bin/env bash
# Build the host-annotated Wolbachia MLST reference panel for evidence line B3
# (strain phylogeny). Downloads reference Wolbachia genomes spanning the competing
# host categories, then extracts the 5 pubMLST loci (gatB/coxA/hcpA/ftsZ/fbpA) from
# each by BLAST so our low-coverage Chilo strains can be placed among them.
#
# Discriminating design: PARASITOID (Hymenoptera) Wolbachia = H1 source; LEPIDOPTERA
# (esp. Ostrinia, the crambid moth closest to Chilo with native Wolbachia) = H0
# comparator. If the Chilo strain clusters with parasitoids and NOT with the
# closely-related corn-borer Wolbachia, that argues against native vertical origin.
set -euo pipefail
source "$(dirname "$0")/common.sh"
source "$HOME/opt/miniforge3/etc/profile.d/conda.sh" && conda activate wolb
cd "$ROOT"
out="refs/wolb_phylo"; gdir="$out/genomes"; ldir="$out/loci"
mkdir -p "$gdir" "$ldir"
mlstdb="$HOME/opt/miniforge3/envs/wolb/db/pubmlst/wolbachia"
loci=(gatB coxA hcpA ftsZ fbpA)

# Panel: accession  strain  host_species  host_type  supergroup(known else ?)
# host_type ∈ {PARASITOID, LEPIDOPTERA, FLY, MOSQUITO, NEMATODE}
panel=$(cat <<'EOF'
GCF_000008025.1	wMel	Drosophila_melanogaster	FLY	A
GCF_000022285.1	wRi	Drosophila_simulans	FLY	A
GCF_000376585.1	wNo	Drosophila_simulans	FLY	B
GCF_000073005.1	wPip	Culex_quinquefasciatus	MOSQUITO	B
GCA_001439985.1	wTpre	Trichogramma_pretiosum	PARASITOID	B
GCA_045765165.1	wTkay	Trichogramma_kaykai	PARASITOID	?
GCA_000204545.1	wVitB	Nasonia_vitripennis	PARASITOID	B
GCA_000174095.1	wUni	Muscidifurax_uniraptor	PARASITOID	A
GCA_006334525.1	wLcla	Leptopilina_clavipes	PARASITOID	?
GCA_039540065.1	wEnc	Encarsia_formosa	PARASITOID	?
GCA_023559125.1	wOfur	Ostrinia_furnacalis	LEPIDOPTERA	?
GCA_023559145.1	wOscap	Ostrinia_scapulalis	LEPIDOPTERA	?
GCA_040687725.1	wEel	Ephestia_elutella	LEPIDOPTERA	?
GCF_000008385.1	wBm	Brugia_malayi	NEMATODE	D
EOF
)

printf "strain\thost_species\thost_type\tsupergroup\taccession\n" > "$out/panel.tsv"
# Fresh per-locus reference files
for L in "${loci[@]}"; do : > "$ldir/$L.refs.fa"; done

while IFS=$'\t' read -r acc strain host htype sg; do
  [[ -z "${acc:-}" ]] && continue
  g="$gdir/$strain.fna"
  if [[ ! -s $g ]]; then
    log "fetch $strain ($acc, $host)"
    tmp="$gdir/.dl_$strain"; rm -rf "$tmp"; mkdir -p "$tmp"
    if datasets download genome accession "$acc" --include genome \
         --filename "$tmp/g.zip" >/dev/null 2>&1 && unzip -oq "$tmp/g.zip" -d "$tmp"; then
      cat "$tmp"/ncbi_dataset/data/"$acc"*/*.fna > "$g" 2>/dev/null || \
        find "$tmp" -name '*.fna' -exec cat {} + > "$g"
    fi
    rm -rf "$tmp"
  fi
  if [[ ! -s $g ]]; then log "WARN: $strain genome empty/failed — skipping"; continue; fi
  printf "%s\t%s\t%s\t%s\t%s\n" "$strain" "$host" "$htype" "$sg" "$acc" >> "$out/panel.tsv"

  # Extract each MLST locus: BLAST all pubMLST alleles vs this genome, take the
  # single best hit, pull that genomic region (revcomp if on minus strand).
  makeblastdb -in "$g" -dbtype nucl -out "$gdir/.db_$strain" >/dev/null 2>&1
  samtools faidx "$g"
  for L in "${loci[@]}"; do
    hit=$(blastn -query "$mlstdb/$L.tfa" -db "$gdir/.db_$strain" \
            -outfmt "6 sseqid sstart send bitscore length" -max_target_seqs 1 2>/dev/null \
          | sort -k4,4nr | head -1 || true)
    [[ -z "$hit" ]] && { log "  $strain $L: no hit"; continue; }
    read -r sid ss se bs ln <<<"$hit"
    if (( se >= ss )); then region="$sid:$ss-$se"; rc=0; else region="$sid:$se-$ss"; rc=1; fi
    seq=$(samtools faidx "$g" "$region" 2>/dev/null | seqkit seq -s -w0 2>/dev/null || true)
    [[ $rc -eq 1 ]] && seq=$(printf "%s\n" "$seq" | rev | tr 'ACGTacgtNn' 'TGCAtgcaNn')
    [[ -z "$seq" ]] && continue
    printf ">%s|%s\n%s\n" "$strain" "$htype" "$seq" >> "$ldir/$L.refs.fa"
  done
  rm -f "$gdir/.db_$strain".*
done <<<"$panel"

log "reference panel built: $(($(wc -l < "$out/panel.tsv")-1)) strains"
for L in "${loci[@]}"; do printf "  %s: %s refs\n" "$L" "$(grep -c '^>' "$ldir/$L.refs.fa" 2>/dev/null || echo 0)"; done
need "$out/panel.tsv" "$ldir/gatB.refs.fa"
