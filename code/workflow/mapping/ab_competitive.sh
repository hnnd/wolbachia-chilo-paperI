#!/bin/bash
# Q3: can A+B co-infection be established? Competitive mapping of the enriched
# Wolbachia read set to a COMBINED wMel(A)+wPip(B) reference, keeping only
# uniquely-placed reads (MAPQ>=20). Conserved-core reads map equally well to both
# references and get MAPQ 0, so they are excluded rather than double-counted --
# which is what the original independent two-reference mapping could not do.
set -u
export PATH="$HOME/miniforge3/envs/wolb/bin:$PATH"
cd /mnt/inaisfs/home/wangyunsheng/wangyunsheng/work/chilo/reseq/01.wolb
refA=refs/wolbachia/strains/NC_002978.6.fa
refB=refs/wolbachia/strains/NC_010981.1.fa
od=results/07_cophylo/ab_competitive; mkdir -p "$od"
comb="$od/AB_combined.fa"

# rename contigs so provenance is unambiguous after concatenation
awk '/^>/{print ">A_wMel"; next}{print}' "$refA"  > "$comb"
awk '/^>/{print ">B_wPip"; next}{print}' "$refB" >> "$comb"
lenA=$(awk '/^>A_/{f=1;next}/^>/{f=0}f{gsub(/[^ACGTNacgtn]/,"");n+=length($0)}END{print n}' "$comb")
lenB=$(awk '/^>B_/{f=1;next}/^>/{f=0}f{gsub(/[^ACGTNacgtn]/,"");n+=length($0)}END{print n}' "$comb")
echo "ref lengths: A_wMel=$lenA  B_wPip=$lenB"

out="$od/ab_competitive.tsv"
printf "sample\ttotal_reads\tmapq0_ambig\tuniqA\tuniqB\tBfrac_uniq\tbreadthA\tbreadthB\tdepthA_cov\tdepthB_cov\tsoloA_uniq\tsoloB_uniq\n" > "$out"

for s in AQ_R_12 AQ_R_8 AQ_R_9 CD_R_10 CD_R_3 CD_R_8 CD_W_1 NC2025_R_1 NC2025_R_6 NL_R_5; do
  w1="results/04_assembly/$s/w_R1.fq.gz"; w2="results/04_assembly/$s/w_R2.fq.gz"
  [[ -s "$w1" && -s "$w2" ]] || { echo "skip $s"; continue; }
  bam="$od/.$s.comb.bam"
  minimap2 -ax sr -t 8 "$comb" "$w1" "$w2" 2>/dev/null | samtools sort -@ 4 -o "$bam" - 2>/dev/null
  samtools index "$bam"
  tot=$(samtools view -c -F 0x900 "$bam")
  amb=$(samtools view -c -F 0x904 -q 0 "$bam" | cat)                      # placeholder, refined below
  amb=$(samtools view -F 0x904 "$bam" | awk '$5==0{n++}END{print n+0}')
  ua=$(samtools view -c -F 0x904 -q 20 "$bam" A_wMel)
  ub=$(samtools view -c -F 0x904 -q 20 "$bam" B_wPip)
  read -r brA dpA < <(samtools depth -a -Q 20 -r A_wMel "$bam" | awk -v L="$lenA" '{n++; if($3>0){c++; s+=$3}}END{printf "%.4f %.2f", (L? c/L:0), (c? s/c:0)}')
  read -r brB dpB < <(samtools depth -a -Q 20 -r B_wPip "$bam" | awk -v L="$lenB" '{n++; if($3>0){c++; s+=$3}}END{printf "%.4f %.2f", (L? c/L:0), (c? s/c:0)}')
  # for contrast: independent single-reference mapping, unique reads only
  sa=$(minimap2 -ax sr -t 8 "$refA" "$w1" "$w2" 2>/dev/null | samtools view -c -F 0x904 -q 20 -)
  sb=$(minimap2 -ax sr -t 8 "$refB" "$w1" "$w2" 2>/dev/null | samtools view -c -F 0x904 -q 20 -)
  bf=$(awk -v a="$ua" -v b="$ub" 'BEGIN{printf "%.4f", ((a+b)>0? b/(a+b):-1)}')
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$s" "$tot" "$amb" "$ua" "$ub" "$bf" "$brA" "$brB" "$dpA" "$dpB" "$sa" "$sb" >> "$out"
  echo "[done] $s tot=$tot ambig=$amb uniqA=$ua uniqB=$ub brA=$brA brB=$brB"
  rm -f "$bam" "$bam".bai
done
cp "$out" ./ab_competitive.tsv
cat "$out"
