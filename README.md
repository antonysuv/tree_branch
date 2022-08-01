# Estimation of branch lengths using deep learning
Solarized dark             |  Solarized Ocean
:-------------------------:|:-------------------------:
![3D rotating gif](https://github.com/antonysuv/tree_branch/blob/main/img/mlp_intraslice_knn_5_dropout_0.0_epoch.gif) | ![3D rotating gif](https://github.com/antonysuv/tree_branch/blob/main/img/mlp_intraslice_knn_5_dropout_0.0_epoch.gif)
## Pipeline overview 
### 1) Generate newick trees with branch lengths  
```
alisim_input_generate.r -d exp,10 -l 1000 -t '((A,B),(C,D));' -n 50000 
```
```
Options:
	-t CHARACTER, --topo=CHARACTER
		Tree in newick e.g. ((A,B),C); or file

	-n NUMERIC, --nsim=NUMERIC
		Number of simulations

	-d CHARACTER, --distribution=CHARACTER
		Distribution

	-l NUMERIC, --len=NUMERIC
		Alignment length

	-f CHARACTER, --fname=CHARACTER
		File name
	
    -p NUMERIC, --prop=NUMERIC
		Proportion of TRAIN data to generate TEST data
	
    -h, --help
		Show this help message and exit
```
### 2) Generate multiple sequence alignment (MSA) for each tree with branch lengths 
```
iqtree2 --alisim alignment_$tr -m JC+G4{uniform} -t $file -af fasta
```
### 3) Perform training
#### 3a) Multi-Layer Perceptron (MLP) using site pattern distributios extracted from MSAs
```
python3.9 keras_MLP_BRANCH.py --tr TRAIN.npy --te TEST.npy --trl train_trees.unrooted.Y.txt  --tel test_trees.unrooted.Y.txt 
```
#### 3b) Convolutional Neural Network using MSAs
```
python3.9 keras_CNN_BRANCH.py --tr TRAIN.npy --te TEST.npy --trl train_trees.unrooted.Y.txt  --tel test_trees.unrooted.Y.txt 
```
#### 3c) Bayesian Convolutional Neural Network using MSAs
```
python3.9 keras_CNNVI_BRANCH.py --tr TRAIN.npy --te TEST.npy --trl train_trees.unrooted.Y.txt  --tel test_trees.unrooted.Y.txt 
```

# Summary table of performed experiments
| Experiment | MLP | CNN | M-PHATE | ML | Branch lengths | MSA model | N taxa | 
| --- | :---: | :---: | :---: | :---: | --- | --- | --- |
| Network generalization | ✅ | ✅ | ✅ | ❌ | ```exp,10``` | JC+G | 4 |
| Branch length heterogeneity (BL-space) | ✅ | ✅ | ❌ | ✅ | ```beta mixture``` | JC+G | 4 |
| Branch length heterogeneity (BL-space) | ✅ | ✅ | ❌ | ✅ | ```beta mixture``` | GTR+G | 4 |
| Exponential branch lengths | ✅ | ✅ | ❌ | ✅ | ```exp,10``` | JC+G | 4 (unrooted) |
| Exponential branch lengths | ✅ | ✅ | ❌ | ✅ | ```exp,10``` | GTR+G | 4 (unrooted) |
| Exponential branch lengths | ✅ | ✅ | ❌ | ✅ | ```exp,10``` | UNREST+G | 4 (rooted) |
| Birth-death branch lengths | ✅ | ✅ | ❌ | ✅ | ```bd,0.05,0.02,300,0.001``` | JC+G | 4 (unrooted) |
| Birth-death branch lengths | ✅ | ✅ | ❌ | ✅ | ```bd,0.05,0.02,300,0.001``` | GTR+G | 4 (unrooted) |
| Birth-death branch lengths | ✅ | ✅ | ❌ | ✅ | ```bd,0.05,0.02,300,0.001``` | UNREST+G | 4 (rooted) |
| Tree shape (balanced) | ✅ | ✅ | ❌ | ✅ | ```exp,10``` | GTR+G | 8 |
| Tree shape (pectinate) | ✅ | ✅ | ❌ | ✅ | ```exp,10``` | GTR+G | 8 |
| Model misspecification (train on GTR test on JC) | ✅ | ✅ | ❌ | ✅ | ```exp,10``` | JC+G | 4 (unrooted) |
| Model misspecification (train on JC test on GTR) | ✅ | ✅ | ❌ | ✅ | ```exp,10``` | JC+G | 4 (unrooted) |



# AliSim input simulation commands
### Network generalization in MPHATE  
```
alisim_input_generate.r -d exp,10 -l 1000 -t '(A,B,(C,D));' -n 50000
```

### BL space
#### 1000 sites
```
alisim_input_generate.r -d mixb -l 1000 -t '(A,B,(C,D));' -n 50000 -p 1
```

#### 5000 sites
```
alisim_input_generate.r -d mixb -l 5000 -t '(A,B,(C,D));' -n 50000 -p 1
```



### Exponential distribution 


#### average branch length = 0.01
```
alisim_input_generate.r -d exp,100 -l 1000 -t '(A,B,(C,D));' -n 50000
```

#### average branch length = 0.1
```
alisim_input_generate.r -d exp,10 -l 1000 -t '(A,B,(C,D));' -n 50000
```

#### average branch length = 1
```
alisim_input_generate.r -d exp,1 -l 1000 -t '(A,B,(C,D));' -n 50000
```

