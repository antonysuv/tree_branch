dir=$1
cd $dir
cd TRAIN
~/scripts/produceTRAIN_alisim_iqtree_array_GTR.sh train_simulation.tre 1000
cd ..
cd TEST
~/scripts/produceTEST_alisim_iqtree_array_JC_GTR.sh test_simulation.tre 1000
cd ..
cd ..