dir=$1
cd $dir
~/scripts/produceKeras_MLP_BRANCH.py 
~/scripts/produceKeras_CNN_BRANCH.py
#~/scripts/produceKeras_CNNVI_BRANCH.py
cd ..