# Estimation of branch lengths using deep learning
![Pipeline](https://github.com/antonysuv/tree_branch/blob/main/img/Fig1.png)
## Pipeline overview 
### 1) Generate newick trees with branch lengths  
```
alisim_input_generate.r -d bd,0.9,200,0.001 -t '((A,B),(C,D));' -m 33b 
```
```
Options:
	-t CHARACTER, --topo=CHARACTER
		Tree in newick e.g. ((A,B),C); or file

	-n NUMERIC, --nsim=NUMERIC
		Number of TRAIN simulations (default: 150000)

	-d CHARACTER, --distribution=CHARACTER
		Branch length distribution (unif,exp,mixb,bd,dir)

	-l NUMERIC, --len=NUMERIC
		Alignment length for IQTREE partiton file (default: 1000)

	-f CHARACTER, --fname=CHARACTER
		Name prefix for main file outputs (default: simulation)

	-p NUMERIC, --ntest=NUMERIC
		Number of TEST simulations (default: 10000)

	-m CHARACTER, --mdir=CHARACTER
		Name suffix for the output dir (default: MODEL)

	-h, --help
		Show this help message and exit

```
### 2) Generate multiple sequence alignment (MSA) for each tree with branch lengths with AliSim (example) 
```
iqtree2 --alisim alignment_$tr -m JC+G4{uniform} -t $file -af fasta
```
### 3) Perform training
#### 3a) Multi-Layer Perceptron (MLP) (+ ROE) using site pattern frequencies extracted from MSAs
```
python3.9 keras_MLP_BRANCH.py --tr TRAIN.npy --te TEST.npy --trl train_trees.unrooted.Y.txt  --tel test_trees.unrooted.Y.txt 
```
#### 3b) Convolutional Neural Network (+ ROE) using MSAs
```
python3.9 keras_CNN_BRANCH.py --tr TRAIN.npy --te TEST.npy --trl train_trees.unrooted.Y.txt  --tel test_trees.unrooted.Y.txt 
```

# Simulation of trees with branch lengths  

### Uniform model 
```
alisim_input_generate.r -d unif,a,b -t '(A,B,(C,D));'
```
| a | b | model |
| --- | :---: | :---: | 
| 0 | 0.001 | JC | 
| 0.001 | 0.01 | JC |
| 0.01 | 0.1 | JC |
| 0.1 | 1 | JC | 
| 1 | 10 | JC | 


### Exponential model
```
alisim_input_generate.r -d exp,theta -t '(A,B,(C,D));'
```
| theta | model | 
| --- | :---: |
| 1 | JC, GTR |   
| 10 | JC, GTR | 
| 100 | JC, GTR |  


### Branch length heterogeneity (BL-space) model
```
alisim_input_generate.r -d mixb -t '(A,B,(C,D));' 
```
| model |
| --- | 
| JC, GTR |



### Birth-death model 
```
alisim_input_generate.r -d bd,turnover,root_age,clock_rate -t '((A,B),(C,D));'
```
| root age | turnover(mu/lambda) | clock rate | model |
| --- | :---: | :---: | :---: | 
| 200 | 0, 0.5, 0.9 | 0.001 | 3.3b, UNREST |  
| 100 | 0, 0.5, 0.9 | 0.001 | 3.3b, UNREST |
| 50 | 0, 0.5, 0.9 | 0.001 | 3.3b, UNREST |
| 10 | 0, 0.5, 0.9 | 0.001 | 3.3b, UNREST |





