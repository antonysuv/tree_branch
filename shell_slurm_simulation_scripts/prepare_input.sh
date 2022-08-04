dir=$1
cd $dir
ls -v TRAIN/aln_*.fa
ls -v TRAIN/aln_*.fa | xargs cat > TRAIN/TRAIN.fa
ls -v TEST/aln_*.fa
ls -v TEST/aln_*.fa | xargs cat > TEST/TEST.fa
ls -v TEST/ml_*.tre
ls -v TEST/ml_*.tre | xargs cat > TEST/ML_all.tre
~/scripts/fasta2numeric.py --tr TRAIN/TRAIN.fa --te TEST/TEST.fa
cd ..