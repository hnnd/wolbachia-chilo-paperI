#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/lib/common.sh"
source "$HOME/opt/miniforge3/etc/profile.d/conda.sh" && conda activate wolb
cd "$ROOT"; td=results/05_typing

# --- CheckM2 completeness (guard: only if DB is available) ---
# Detect DB via $CHECKM2DB env var pointing to an existing file, or the
# conventional local path refs/checkm2_db/uniref100.KO.1.dmnd.
_checkm2_db="${CHECKM2DB:-}"
_local_db="$ROOT/refs/checkm2_db/uniref100.KO.1.dmnd"
_checkm2_ok=0
_db_to_pass=""
if [[ -n "$_checkm2_db" && -f "$_checkm2_db" ]]; then
    _checkm2_ok=1
    _db_to_pass="${CHECKM2DB:-$_local_db}"
elif [[ -f "$_local_db" ]]; then
    _checkm2_ok=1
    _db_to_pass="$_local_db"
fi

# Guard: safe glob for assembly fastas (nullglob-safe via find)
shopt -s nullglob
_asms=( results/04_assembly/*/assembly.fasta )
shopt -u nullglob

if [[ ${#_asms[@]} -eq 0 ]]; then
    log "no assemblies yet — writing empty typing_summary.tsv and exiting"
    printf "sample\tST\tcompleteness\tcontamination\n" > "$td/typing_summary.tsv"
    exit 0
fi

# Build bins dir for CheckM2
mkdir -p "$td/bins"
for d in "${_asms[@]}"; do
    s=$(basename "$(dirname "$d")")
    ln -sf "$ROOT/$d" "$td/bins/$s.fasta"
done

if [[ $_checkm2_ok -eq 1 ]]; then
    log "running CheckM2 on ${#_asms[@]} assembly bins"
    conda run -n wolb-checkm2 checkm2 predict -i "$td/bins" -x fasta -o "$td/checkm2" -t "$(cfg threads)" --force --database-path "$_db_to_pass"
else
    log "CheckM2 DB not found (CHECKM2DB unset and $ROOT/refs/checkm2_db absent) — skipping; completeness/contamination will be NA"
fi

# Summarise MLST + completeness/contamination
{ printf "sample\tST\tcompleteness\tcontamination\n"
  for s in $(cat "$td/pos_samples.txt"); do
    st=$(awk '{print $3}' "$td/$s.mlst.tsv" 2>/dev/null | head -1)
    if [[ $_checkm2_ok -eq 1 && -f "$td/checkm2/quality_report.tsv" ]]; then
        cm=$(awk -F'\t' -v n="$s" '$1==n{print $2"\t"$3}' "$td/checkm2/quality_report.tsv")
    else
        cm="NA\tNA"
    fi
    printf "%s\t%s\t%s\n" "$s" "${st:-NA}" "${cm:-NA\tNA}"
  done; } > "$td/typing_summary.tsv"
need "$td/typing_summary.tsv"; column -t "$td/typing_summary.tsv" | head
