#!/usr/bin/env bash
#SBATCH -p cpu -c 16 --mem=32G -o logs/prof_%A_%a.out -e logs/prof_%A_%a.err
set -euo pipefail
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}"   # SLURM copies script to spool; locate repo via submit dir
source "$HERE/workflow/lib/common.sh"
source "$HOME/opt/miniforge3/etc/profile.d/conda.sh" && conda activate wolb
i=${SLURM_ARRAY_TASK_ID:?}; T=$(cfg threads)
sid=$(sample_field "$i" 1)
nh="$ROOT/results/02_nonhost"; out="$ROOT/results/03_profile"; mkdir -p "$out"
r1="$nh/${sid}_R1.fq.gz"; r2="$nh/${sid}_R2.fq.gz"; need "$r1" "$r2"
[[ -f $out/$sid.metrics.done ]] && exit 0
nonhost=$(($(zcat "$r1" | wc -l)/4*2))
# --- Wolbachia: 对每株各自独立比对，取 breadth 最高株 ---
# 把 4 株合并成单参考会因保守区多重比对，把 primary 随机分散到各株，系统性低估
# 真实感染株的 breadth/depth（实测 wTpre 0.27->0.85）。单株独立比对恢复真实信号。
# 见 memory/merged-ref-distorts-coverage。breadth=覆盖>=1x 比例，depth=全长 mean。
wolb_strain="NA"; wolb_breadth=0; wolb_depth=0; wolb_rpm=0
while read -r sn; do
  [[ -n "$sn" ]] || continue
  sb="$out/$sid.$sn.bam"
  minibwa map -x sr -t "$T" "refs/wolbachia/strains/$sn" "$r1" "$r2" \
    | samtools sort -@ "$T" -o "$sb" - && samtools index "$sb"
  bd=$(samtools depth -a "$sb" | awk '{t++; s+=$3; if($3>0)c++} END{printf "%.4f %.2f", t?c/t:0, t?s/t:0}') \
    || die "samtools depth failed on $sb"   # here-string would swallow this failure
  read b d <<<"$bd"
  m=$(samtools view -c -F 0x904 "$sb")
  r=$(awk -v m="$m" -v n="$nonhost" 'BEGIN{printf "%.2f", n?m/n*1e6:0}')
  if awk -v a="$b" -v c="$wolb_breadth" 'BEGIN{exit !(a>c)}'; then
    wolb_strain="$sn"; wolb_breadth="$b"; wolb_depth="$d"; wolb_rpm="$r"
  fi
  rm -f "$sb"*
done < refs/wolbachia/strains/strain_list.txt
# --- 寄生蜂 ---
minibwa map -x sr -t "$T" refs/parasitoid/parasitoid_refs "$r1" "$r2" \
  | samtools view -c -F 0x904 - > "$out/$sid.para.count"
para_reads=$(cat "$out/$sid.para.count")
para_frac=$(awk -v m=$para_reads -v n=$nonhost 'BEGIN{printf "%.6f", n?m/n:0}')
# Kraken2 taxonomic profile runs as a separate, high-memory pass (Stage 2c,
# workflow/jobs/kraken_one.sh) — the standard DB (~104GB) needs ~150G RAM,
# far more than this light 32G job; keeping it out keeps profiling fast.
printf "sample\twolb_strain\twolb_breadth\twolb_depth\twolb_rpm\tpara_reads\tpara_frac\n%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
  "$sid" "$wolb_strain" "$wolb_breadth" "$wolb_depth" "$wolb_rpm" "$para_reads" "$para_frac" > "$out/$sid.metrics.tsv"
rm -f "$out/$sid.para.count"
mark_done "$out/$sid.metrics"
