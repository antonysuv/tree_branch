# Estimation of branch lengths using deep learning

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
