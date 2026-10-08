#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/lib/common.sh"
source "$HOME/opt/miniforge3/etc/profile.d/conda.sh" && conda activate wolb
out="$ROOT/results/03_profile"
# 合并 metrics
csvtk concat -t "$out"/*.metrics.tsv > "$out/_metrics_all.tsv"
# 合并 hoststat
csvtk concat -t "$ROOT"/results/02_nonhost/*.hoststat.txt > "$out/_host_all.tsv"
# join 元数据 + 阳性判定
csvtk join -t -f sample "$out/_metrics_all.tsv" \
  <(csvtk rename -t -f 1 -n sample "$out/_host_all.tsv") \
  | csvtk join -t -f "sample;sample_id" - config/samples.tsv \
  | csvtk mutate2 -t -n wolb_pos -e "\$wolb_breadth>=$(cfg wolb_min_breadth) && \$wolb_rpm>=$(cfg wolb_min_rpm) ? 1:0" \
  | csvtk mutate2 -t -n para_pos -e "\$para_frac>=$(cfg para_min_frac) ? 1:0" \
  > "$out/per_sample_profile.tsv"
need "$out/per_sample_profile.tsv"

# --- optional de-novo parasitoid rescue (Stage 03f) --------------------------
# A confidently assembled parasitoid mitogenome (>=1 kb) marks a truly parasitized
# sample even when its Cotesia-referenced para_frac is below threshold — fly/
# ichneumonid parasitoids are systematically under-counted by the Cotesia reference
# (e.g. tachinid-parasitized HC_W_9). Best-effort: applied only when de novo exists.
dn="$out/../06_denovo/per_sample_denovo_id.tsv"
if [[ -s "$dn" ]]; then
  awk -F'\t' 'NR>1{ g=$7; ml=$10+0;
      p=(ml>=1000 && index("|Cotesia|Pexopsis|Scambus|Neotrichoporoides|Lypha|Brachymeria|Trichogramma|Trichomalopsis|Pachyneuron|Tetrastichus|","|"g"|"))?1:0;
      print $1"\t"p }' "$dn" > "$out/_denovo_para.tsv"
  awk -F'\t' -v OFS='\t' -v PT="$(cfg para_min_frac)" '
    NR==FNR { dp[$1]=$2; next }
    FNR==1 { for(i=1;i<=NF;i++) h[$i]=i; print; next }
    { s=$(h["sample"]); pf=$(h["para_frac"])+0; dv=(s in dp?dp[s]:0);
      $(h["para_pos"])=(pf>=PT || dv==1)?1:0; print }
  ' "$out/_denovo_para.tsv" "$out/per_sample_profile.tsv" > "$out/_recon.tsv" \
    && mv "$out/_recon.tsv" "$out/per_sample_profile.tsv"
  rm -f "$out/_denovo_para.tsv"
  echo "de-novo parasitoid rescue applied (Stage 03f present)"
fi
echo "wolb_pos / para_pos 计数："
csvtk -t freq -f wolb_pos "$out/per_sample_profile.tsv"
csvtk -t freq -f para_pos "$out/per_sample_profile.tsv"
