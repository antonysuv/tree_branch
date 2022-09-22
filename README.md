# Estimation of branch lengths using deep learning
Multilayer perceptron | Convolutional neural network
:-------------------------:|:-------------------------:
![3D rotating gif](https://github.com/antonysuv/tree_branch/blob/main/img/mlp_intraslice_knn_5_dropout_0.0_epoch.gif) | ![3D rotating gif](https://github.com/antonysuv/tree_branch/blob/main/img/cnn_intraslice_knn_2_dropout_0.0_epoch.gif)
![3D rotating gif](https://github.com/antonysuv/tree_branch/blob/main/img/mlp_intraslice_knn_5_dropout_0.0_valloss.gif) | ![3D rotating gif](https://github.com/antonysuv/tree_branch/blob/main/img/cnn_intraslice_knn_2_dropout_0.0_valloss.gif)
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

### Uniform model 
```
alisim_input_generate.r -d unif,a,b -t '(A,B,(C,D));'
```
| a | b | 
| --- | :---: | 
| 0 | 0.001 | 
| 0.001 | 0.1 |
| 0.1 | 1 | 
| 1 | 10| 

### BL-space
```
alisim_input_generate.r -d mixb -t '(A,B,(C,D));' 
```
### Exponential model
```
alisim_input_generate.r -d exp,theta -t '(A,B,(C,D));'
```
| theta |
| --- |
| 1 |  
| 10 |
| 100 | 


### Birth-death model 
```
alisim_input_generate.r -d bd,mu/lambda,root_age,clock_rate -t '((A,B),(C,D));'
```
| root age | mu/lambda | clock rate |
| --- | :---: | :---: | 
| 200 | 0, 0.5, 0.9 | 0.001 | 
| 100 | 0, 0.5, 0.9 | 0.001 |
| 50 | 0, 0.5, 0.9 | 0.001 |
| 10 | 0, 0.5, 0.9 | 0.001 |



### Network generalization in MPHATE (JC) 
```
alisim_input_generate.r -d exp,10 -l 1000 -t '(A,B,(C,D));' -n 150000
```

### BL space
#### 1000 sites (JC and GTR)
```
alisim_input_generate.r -d mixb -t '(A,B,(C,D));' 
```

#### 5000 sites (JC)
```
alisim_input_generate.r -d mixb -l 5000 -t '(A,B,(C,D));' -n 150000 -p 1
```



### Exponential distribution (JC and GTR)


#### average branch length = 0.01
```
alisim_input_generate.r -d exp,100 -l 1000 -t '(A,B,(C,D));' -n 150000
```

#### average branch length = 0.1
```
alisim_input_generate.r -d exp,10 -l 1000 -t '(A,B,(C,D));' -n 150000
```

#### average branch length = 1
```
alisim_input_generate.r -d exp,1 -l 1000 -t '(A,B,(C,D));' -n 150000
```

#### 25 taxa
```
set.seed(45)
phy=rtree(25, rooted = F, tip.label = LETTERS[1:25], br = NULL, equiprob = FALSE)
(G,P,(D,((W,C),(((X,I),((O,(B,E)),(N,K))),(((((J,M),((R,Q),T)),(A,S)),(L,(V,(U,F)))),(H,Y))))));
alisim_input_generate.r -d exp,10 -l 1000 -t '(G,P,(D,((W,C),(((X,I),((O,(B,E)),(N,K))),(((((J,M),((R,Q),T)),(A,S)),(L,(V,(U,F)))),(H,Y))))));' -n 50000
alisim_input_generate.r -d exp,10 -l 1000 -t '(G,P,(D,((W,C),(((X,I),((O,(B,E)),(N,K))),(((((J,M),((R,Q),T)),(A,S)),(L,(V,(U,F)))),(H,Y))))));' -n 150000 -p 0.0333333333 
```

#### 25 taxa
```
set.seed(45)
phy=rtree(25, rooted = F, tip.label = LETTERS[1:25], br = NULL, equiprob = FALSE)
(G,P,(D,((W,C),(((X,I),((O,(B,E)),(N,K))),(((((J,M),((R,Q),T)),(A,S)),(L,(V,(U,F)))),(H,Y))))));
alisim_input_generate.r -d exp,10 -l 1000 -t '(G,P,(D,((W,C),(((X,I),((O,(B,E)),(N,K))),(((((J,M),((R,Q),T)),(A,S)),(L,(V,(U,F)))),(H,Y))))));' -n 50000
```

#### 50 taxa
```
set.seed(132)
phy=rtree(50, rooted = F, br = rexp, equiprob = FALSE)
(t6,t45,((t25,t23),(((t15,t10),((t39,t50),t9)),(((t37,(t46,t2)),((((t30,(t28,(t26,(t49,t5)))),t8),((((t42,t44),t20),t36),((((t38,t18),t22),t41),t32))),((t4,(t21,t48)),(t17,t33)))),((t24,t47),((((t3,t35),((t1,(t29,t31)),(((t27,t11),(t7,t13)),t43))),((t19,(t40,t16)),t14)),(t12,t34)))))));
alisim_input_generate.r -d exp,10 -l 1000 -t '(t6,t45,((t25,t23),(((t15,t10),((t39,t50),t9)),(((t37,(t46,t2)),((((t30,(t28,(t26,(t49,t5)))),t8),((((t42,t44),t20),t36),((((t38,t18),t22),t41),t32))),((t4,(t21,t48)),(t17,t33)))),((t24,t47),((((t3,t35),((t1,(t29,t31)),(((t27,t11),(t7,t13)),t43))),((t19,(t40,t16)),t14)),(t12,t34)))))));' -n 50000
```