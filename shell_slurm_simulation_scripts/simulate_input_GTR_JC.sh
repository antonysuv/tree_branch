dir=$1
cd $dir
cd TRAIN
~/scripts/produceTRAIN_alisim_iqtree_array.sh train_simulation.tre 1000
cd ..
cd TEST
~/scripts/produceTEST_alisim_iqtree_array_GTR_JC.sh test_simulation.tre 1000
cd ..
cd ..