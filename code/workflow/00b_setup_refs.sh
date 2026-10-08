#!/usr/bin/env bash
# Stage 0b: download + index all references. Idempotent via .done sentinels.
# NOTE: cluster proxy throttles large downloads to ~15KB/s, but DIRECT internet
# (no proxy) reaches S3/NCBI at ~1GB/s. So big downloads use aria2c --no-proxy.
# Mapping index uses minibwa (minibwa index -> <prefix>.l2b/.mbw).
set -euo pipefail
source "$(dirname "$0")/lib/common.sh"
source "$HOME/opt/miniforge3/etc/profile.d/conda.sh" && conda activate wolb
cd "$ROOT"; mkdir -p refs/{host,wolbachia,parasitoid,kraken2_db}

# direct (no-proxy) multi-connection download via aria2c (in its own env)
aria2_dl(){ # url outdir outname
  env -u http_proxy -u https_proxy -u HTTP_PROXY -u HTTPS_PROXY -u all_proxy \
    conda run -n aria2 aria2c -x16 -s16 --no-proxy='*' -d "$2" -o "$3" "$1"
}

# 1) Host: Chilo suppressalis GCA_902850365.2
if [[ ! -f refs/host/chilo_suppressalis.done ]]; then
  datasets download genome accession GCA_902850365.2 --include genome --filename refs/host/chilo.zip
  unzip -o refs/host/chilo.zip -d refs/host/_dl
  cat refs/host/_dl/ncbi_dataset/data/*/*.fna > refs/host/chilo_suppressalis.fa
  minibwa index refs/host/chilo_suppressalis.fa refs/host/chilo_suppressalis
  mark_done refs/host/chilo_suppressalis
fi

# 2) Wolbachia representative strains (incl. parasitoid-source wTpre GCF_001439985.1)
if [[ ! -f refs/wolbachia/wolbachia_refs.done ]]; then
  bash workflow/lib/fetch_wolbachia_refs.sh refs/wolbachia/wolbachia_refs.fa
  minibwa index refs/wolbachia/wolbachia_refs.fa refs/wolbachia/wolbachia_refs
  mark_done refs/wolbachia/wolbachia_refs
fi

# 2b) Per-strain independent minibwa indexes. Merging all strains into one
#     reference distorts breadth/depth: conserved regions multimap, primary is
#     split arbitrarily across strains, undercounting the true strain's coverage
#     (measured wTpre 0.27->0.85; see memory/merged-ref-distorts-coverage).
#     Profiling maps each sample to every strain independently and keeps the
#     best-breadth strain. strain_list.txt drives the per-strain loop.
if [[ ! -f refs/wolbachia/strains.done ]]; then
  mkdir -p refs/wolbachia/strains
  samtools faidx refs/wolbachia/wolbachia_refs.fa
  cut -f1 refs/wolbachia/wolbachia_refs.fa.fai > refs/wolbachia/strains/strain_list.txt
  while read -r sn; do
    [[ -n "$sn" ]] || continue
    samtools faidx refs/wolbachia/wolbachia_refs.fa "$sn" > "refs/wolbachia/strains/$sn.fa"
    minibwa index "refs/wolbachia/strains/$sn.fa" "refs/wolbachia/strains/$sn"
  done < refs/wolbachia/strains/strain_list.txt
  mark_done refs/wolbachia/strains
fi

# 3) Parasitoid (Cotesia chilonis / fallback C. congregata) + COI
if [[ ! -f refs/parasitoid/parasitoid_refs.done ]]; then
  bash workflow/lib/fetch_parasitoid_refs.sh refs/parasitoid/parasitoid_refs.fa
  minibwa index refs/parasitoid/parasitoid_refs.fa refs/parasitoid/parasitoid_refs
  mark_done refs/parasitoid/parasitoid_refs
fi

# 4) Kraken2 PlusPF prebuilt DB (16 GB) — aria2c direct (proxy is too slow)
if [[ ! -f refs/kraken2_db.done ]]; then
  [[ -s refs/kraken2_db/db.tar.gz ]] || aria2_dl \
    https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_16gb_20240112.tar.gz refs/kraken2_db db.tar.gz
  tar -xzf refs/kraken2_db/db.tar.gz -C refs/kraken2_db
  mark_done refs/kraken2_db
fi

# 5) CheckM2 diamond DB (~3 GB; Stage 4). checkm2 fetches from its host; if the
#    proxy throttles it, set CHECKM2DB_URL to a direct URL and aria2_dl it instead.
if [[ ! -f refs/checkm2_db.done ]]; then
  CHECKM2DB=refs/checkm2_db; mkdir -p "$CHECKM2DB"
  env -u http_proxy -u https_proxy -u HTTP_PROXY -u HTTPS_PROXY -u all_proxy \
    conda run -n wolb-checkm2 checkm2 database --download --path "$CHECKM2DB"
  mark_done refs/checkm2_db
fi

# MD5 registry
{ echo "# Reference MD5 ($(date '+%F'))"; find refs -maxdepth 2 -name '*.fa' -exec md5sum {} \; ; } >> refs/REFERENCES.md
log "all references ready"
