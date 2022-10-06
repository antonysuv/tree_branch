#!/bin/bash
trees=$1
len=$2

split -d -l 1000 -a6 $trees chunk_ --additional-suffix=.tre --numeric-suffixes=1

sed -n 1p "$trees" > fixed_tree.tre

n=$(($(cat $trees | wc -l)/1000))

batch=alisim_iqtree_array.sh

echo '#!/bin/bash
#SBATCH --job-name=alisim
#SBATCH -N 1
#SBATCH --mem=1G
#SBATCH -n 1
#SBATCH -t 100:00:00
#SBATCH --array=1-'"$n"'
#SBATCH --requeue' >> $batch

echo 'date
start=$((${SLURM_ARRAY_TASK_ID} * 1000 - 1000 + 1))
stop=$((${SLURM_ARRAY_TASK_ID} * 1000))
slurm_id=$(printf "%06d" ${SLURM_ARRAY_TASK_ID})
split -d -l 1 -a6 chunk_$slurm_id.tre in_tree_ --additional-suffix=.tre --numeric-suffixes=$start
for tr in `seq -f '%06g' $start $stop`
do
myrand=$((`shuf -i 0-4294967295 -n1`))
file=in_tree_$tr.tre 
P1=$(seq 0 .001 0.979 | shuf | head -n1)
P2=$(seq 0 .001 0.979 | shuf | head -n1)
P3=$(seq 0 .001 0.979 | shuf | head -n1)
P4=$(seq 0 .001 0.979 | shuf | head -n1)
P5=$(seq 0 .001 0.979 | shuf | head -n1)
P6=$(seq 0 .001 0.979 | shuf | head -n1)
P7=$(seq 0 .001 0.979 | shuf | head -n1)
P8=$(seq 0 .001 0.979 | shuf | head -n1)
P9=$(seq 0 .001 0.979 | shuf | head -n1)
P10=$(seq 0 .001 0.979 | shuf | head -n1)
P11=$(seq 0 .001 0.979 | shuf | head -n1)
~/soft/iqtree-2.2.0-Linux/bin/iqtree2 --alisim alignment_$tr -m 12.12{$P1/$P2/$P3/$P4/$P5/$P6/$P7/$P8/$P9/$P10/$P11}+F{uniform/uniform/uniform/uniform}+G4{uniform} -t $file --length '$len' -af fasta -seed $myrand
~/soft/iqtree-2.2.0-Linux/bin/iqtree2 -s alignment_$tr.fa -t fixed_tree.tre -m 12.12+G4 -blmin 1e-300 -quiet -pre brls_$tr -safe -redo --tree-fix
rm brls_$tr.ckp.gz brls_$tr.iqtree brls_$tr.log $file
done

for tr in `seq -f '%06g' $start $stop`
do
cat in_tree_$tr.tre.log >> aln_${SLURM_ARRAY_TASK_ID}.log
rm in_tree_$tr.tre.log
cat alignment_$tr.fa >> aln_${SLURM_ARRAY_TASK_ID}.fa
rm alignment_$tr.fa
cat brls_$tr.treefile >> ml_${SLURM_ARRAY_TASK_ID}.tre
rm brls_$tr.treefile
done
date' >> $batch
sbatch $batch