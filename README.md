# Estimation of branch lengths using deep learning
## Pipeline overview 
### 1) Generate newick trees with branch lengths  
```
alisim_input_generate.r -d exp,10 -l 1000 -t '((A,B),(C,D));' -n 50000 -f trees
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
