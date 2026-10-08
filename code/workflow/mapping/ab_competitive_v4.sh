#SBATCH --partition=cpu
#SBATCH -c 32
#SBATCH --mem=64G
#SBATCH --time=05:00:00

# Paper I, open item 8: competitive mapping of the non-host read set of every
# Wolbachia-positive library against a COMBINED wCchiA + wCchiB reference,
# replacing the earlier wMel/wPip external-reference scheme. Reads from the
# conserved core map equally well to both strains and receive MAPQ 0; they are
# counted as ambiguous rather than assigned to either strain, which is exactly
# what independent single-reference mapping cannot do.

source ~/miniforge3/etc/profile.d/conda.sh
conda activate wolb
BASE=/mnt/inaisfs/home/wangyunsheng/wangyunsheng/work/chilo/reseq/01.wolb
THREADS=32

refA=$BASE/results/14_strainA_split/bins/strainA_main.fna
refB=$BASE/results/14_strainA_split/bins/strainB.fna
for f in "$refA" "$refB"; do
  n=$(grep -c '^>' "$f"); echo "$f : $n contig(s)"
  [ "$n" -eq 1 ] || { echo "FATAL: expected a single contig in $f"; exit 1; }
done

awk '/^>/{print ">A_wCchiA"; next}{print}' "$refA"  > comb.fa
awk '/^>/{print ">B_wCchiB"; next}{print}' "$refB" >> comb.fa
lenA=$(awk '/^>A_/{f=1;next}/^>/{f=0}f{gsub(/[^ACGTacgt]/,"");n+=length($0)}END{print n}' comb.fa)
lenB=$(awk '/^>B_/{f=1;next}/^>/{f=0}f{gsub(/[^ACGTacgt]/,"");n+=length($0)}END{print n}' comb.fa)
echo "non-N reference lengths: A_wCchiA=$lenA  B_wCchiB=$lenB"

out=ab_competitive_v4.tsv
printf "sample\tnonhost_reads\tmapped_primary\tambiguous_mapq0\tuniqA\tuniqB\tB_fraction\tbreadthA\tbreadthB\tdepthA_covered\tdepthB_covered\twinA_10kb\twinB_10kb\tpct_winA\tpct_winB\n" > "$out"

while read -r s; do
  [ -n "$s" ] || continue
  r1=$BASE/results/02_nonhost/${s}_R1.fq.gz
  r2=$BASE/results/02_nonhost/${s}_R2.fq.gz
  if [ ! -s "$r1" ] || [ ! -s "$r2" ]; then echo "[skip] $s : non-host reads absent"; continue; fi
  nh=$(( $(zcat "$r1" | wc -l) / 4 ))
  bam=.$s.bam
  minimap2 -ax sr -t $THREADS comb.fa "$r1" "$r2" 2>/dev/null | samtools sort -@ 6 -o "$bam" -
  samtools index "$bam"
  mapped=$(samtools view -c -F 0x904 "$bam")
  amb=$(samtools view -F 0x904 "$bam" | awk '$5==0{n++}END{print n+0}')
  ua=$(samtools view -c -F 0x904 -q 20 "$bam" A_wCchiA)
  ub=$(samtools view -c -F 0x904 -q 20 "$bam" B_wCchiB)
  samtools depth -a -Q 20 -r A_wCchiA "$bam" | awk -v L="$lenA" '{if($3>0){c++;s+=$3;w[int($2/10000)]=1}}END{printf "%.6f %.2f %d\n", (L?c/L:0), (c?s/c:0), length(w)}' > .statA
  samtools depth -a -Q 20 -r B_wCchiB "$bam" | awk -v L="$lenB" '{if($3>0){c++;s+=$3;w[int($2/10000)]=1}}END{printf "%.6f %.2f %d\n", (L?c/L:0), (c?s/c:0), length(w)}' > .statB
  read -r brA dpA wA < .statA
  read -r brB dpB wB < .statB
  bf=$(awk -v a="$ua" -v b="$ub" 'BEGIN{printf "%.4f", ((a+b)>0? b/(a+b) : -1)}')
  pwa=$(awk -v w="$wA" -v L="$lenA" 'BEGIN{printf "%.1f", 100*w/int(L/10000)}')
  pwb=$(awk -v w="$wB" -v L="$lenB" 'BEGIN{printf "%.1f", 100*w/int(L/10000)}')
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$s" "$nh" "$mapped" "$amb" "$ua" "$ub" "$bf" "$brA" "$brB" "$dpA" "$dpB" "$wA" "$wB" "$pwa" "$pwb" >> "$out"
  echo "[done] $s nonhost=$nh mapped=$mapped ambig=$amb uniqA=$ua uniqB=$ub brA=$brA brB=$brB Bfrac=$bf"
  rm -f "$bam" "$bam".bai
done < samples_competitive.txt

echo "=== result ==="
cat "$out"
ls -lh "$out" comb.fa || true
