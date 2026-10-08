#!/usr/bin/env bash
# Stage 0a: install Miniforge (reuse existing) and create the analysis conda envs.
# Envs are split to avoid solver conflicts/slow solves:
#   wolb         core tools (Stage 1-2 critical path + assembly/typing/phylo)
#   wolb-checkm2 checkm2 only (Stage 4; bundles conflicting deps)
#   wolb-r       R + tidyverse (Stage 3 stats/plots)
set -euo pipefail
source "$(dirname "$0")/lib/common.sh"
MF="$HOME/opt/miniforge3"
# Reuse existing Miniforge at ~/miniforge3 (cluster pre-install) via symlink.
if [[ ! -x $MF/bin/conda && -x $HOME/miniforge3/bin/conda ]]; then
  log "found existing Miniforge at $HOME/miniforge3, symlinking to $MF"
  mkdir -p "$HOME/opt"
  ln -s "$HOME/miniforge3" "$MF"
fi
if [[ ! -x $MF/bin/conda ]]; then
  log "installing Miniforge to $MF"
  wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/mf.sh
  bash /tmp/mf.sh -b -p "$MF"
fi
source "$MF/etc/profile.d/conda.sh"
export CONDA_ALWAYS_YES=true   # non-interactive: never block on Confirm [Y/n]

# Create one env from a yaml; idempotent (skip if it already has the sentinel).
create_env(){
  local name=$1 yaml=$2
  if conda env list | awk '{print $1}' | grep -qx "$name"; then
    log "env $name exists, updating"
    mamba env update -y -n "$name" -f "$yaml"
  else
    log "creating env $name"
    mamba env create -y -f "$yaml"
  fi
  conda list -n "$name" > "$ROOT/envs/$name.lock.txt"
  mark_done "$ROOT/envs/$name"
  log "env $name ready"
}

create_env wolb         "$ROOT/envs/wolb.yaml"
create_env wolb-checkm2 "$ROOT/envs/wolb-checkm2.yaml"
create_env wolb-r       "$ROOT/envs/wolb-r.yaml"
log "all envs ready"
