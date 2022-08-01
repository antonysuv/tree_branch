#!/bin/bash
trees=$1
len=$2

n=$( cat $trees | wc -l )
sed -n 1p "$trees" > fixed_tree.tre
for i in `seq 1 $n`
do 
	sed -n "$i"p "$trees" > in_tree_"$i".tre
done

batch=alisim_iqtree_array.sh

echo '#!/bin/bash
#SBATCH --job-name=alisim
#SBATCH -N 1
#SBATCH --mem=1G
#SBATCH -n 1
#SBATCH -t 100:00:00
#SBATCH --array=1-'"$n"':1000
#SBATCH --requeue' >> $batch

echo 'date
range=$((${SLURM_ARRAY_TASK_ID}+999))
for tr in `seq ${SLURM_ARRAY_TASK_ID} $range`
do
myrand=$((`shuf -i 0-4294967295 -n1`))
file=in_tree_$tr.tre 
~/soft/iqtree-2.2.0-Linux/bin/iqtree2 --alisim alignment_$tr -m JC+G4{uniform} -t $file --length '$len' -af fasta -seed $myrand
~/soft/iqtree-2.2.0-Linux/bin/iqtree2 -s alignment_$tr.fa -t fixed_tree.tre -m JC+G4 -blmin 1e-300 -quiet -pre brls_$tr -safe -redo --tree-fix
rm brls_$tr.ckp.gz brls_$tr.iqtree brls_$tr.log $file
done

for tr in `seq ${SLURM_ARRAY_TASK_ID} $range`
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

