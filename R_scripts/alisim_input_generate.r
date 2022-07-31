list_of_packages = c("TreeSim", "NELSI","MCMCpack","optparse","devtools","dplyr")
new_packages = list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
if (any("NELSI" == new_packages))
{
    install.packages(new_packages[new_packages!="NELSI"],repos = "http://cran.us.r-project.org")
    library("devtools")
    install_github("sebastianduchene/NELSI")
}else if (length(new_packages)>0)
{
    install.packages(new_packages,repos = "http://cran.us.r-project.org")
}    

library("dplyr")
library("TreeSim")
library("NELSI")
library("MCMCpack")
library("optparse")
options(scipen = 100000000)

option_list = list(
  make_option(c("-t", "--topo"), type="character", default="((A,B),C);",help="Tree in newick e.g. ((A,B),C); or file", metavar="character"),
  make_option(c("-n", "--nsim"), type="numeric", default=10, help="Number of simulations", metavar="numeric"),  
  make_option(c("-d", "--distribution"), type="character", default="unif,0,10", help="Branch length distribution", metavar="character"),
  make_option(c("-l", "--len"), type="numeric", default=100, help="Alignment length", metavar="numeric"),
  make_option(c("-f", "--fname"), type="character", default="simulation", help="File name", metavar="character"),
  make_option(c("-p", "--prop"), type="numeric", default=0.1, help="Proportion of TRAIN data to generate TEST data", metavar="numeric")  
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#####Example
#Rscript ~/Desktop/git_repos/tree_branch/R_scripts/alisim_input_generate.r -d exp,5 -l 1000 -t '((A,B),C);' -n 50000 -f test_trees
############

my_seed=runif(1,0,429496729)
set.seed(my_seed)

#Disable scientific notation in ape 
.write.tree2 <- function(phy, digits = 10, tree.prefix = "", check_tips)
{
    brl <- !is.null(phy$edge.length)
    nodelab <- !is.null(phy$node.label)
    #if (check_tips) phy$tip.label <- checkLabel(phy$tip.label)
    #if (nodelab) phy$node.label <- checkLabel(phy$node.label)
    f.d <- paste("%.", digits, "f", sep = "")
    cp <- function(x){
        STRING[k] <<- x
        k <<- k + 1
    }
    add.internal <- function(i) {
        cp("(")
        desc <- kids[[i]]
        for (j in desc) {
            if (j > n) add.internal(j)
            else add.terminal(ind[j])
            if (j != desc[length(desc)]) cp(",")
        }
        cp(")")
        if (nodelab && i > n) cp(phy$node.label[i - n]) # fixed by Naim Matasci (2010-12-07)
        if (brl) {
            cp(":")
            cp(sprintf(f.d, phy$edge.length[ind[i]]))
        }
    }
    add.terminal <- function(i) {
        cp(phy$tip.label[phy$edge[i, 2]])
        if (brl) {
            cp(":")
            cp(sprintf(f.d, phy$edge.length[i]))
        }
    }

    n <- length(phy$tip.label)

    ## borrowed from phangorn:
    parent <- phy$edge[, 1]
    children <- phy$edge[, 2]
    kids <- vector("list", n + phy$Nnode)
    for (i in 1:length(parent))
        kids[[parent[i]]] <- c(kids[[parent[i]]], children[i])

    ind <- match(1:max(phy$edge), phy$edge[, 2])

    LS <- 4*n + 5
    if (brl) LS <- LS + 4*n
    if (nodelab)  LS <- LS + n
    STRING <- character(LS)
    k <- 1
    cp(tree.prefix)
    cp("(")
    getRoot <- function(phy)
        phy$edge[, 1][!match(phy$edge[, 1], phy$edge[, 2], 0)][1]
    root <- getRoot(phy) # replaced n+1 with root - root has not be n+1
    desc <- kids[[root]]
    for (j in desc) {
        if (j > n) add.internal(j)
        else add.terminal(ind[j])
        if (j != desc[length(desc)]) cp(",")
    }

    if (is.null(phy$root.edge)) {
        cp(")")
        if (nodelab) cp(phy$node.label[1])
        cp(";")
    }
    else {
        cp(")")
        if (nodelab) cp(phy$node.label[1])
        cp(":")
        cp(sprintf(f.d, phy$root.edge))
        cp(";")
    }
    paste(STRING, collapse = "")
}

assignInNamespace(".write.tree2", .write.tree2, "ape")


#RevBayes template generate
get_revbayes_in = function(birth_rate,death_rate,root_age,nsim)
{    
    my_template= paste("
T <- readTrees(\"input_tree.tre\")[1]
taxa <- T.taxa()
moves = VectorMoves()
monitors = VectorMonitors()
root_age ~ dnUniform(",root_age,",",root_age+0.0001,")
timetree ~ dnBDP(lambda=",birth_rate,", mu=",death_rate,", rho=1, rootAge=root_age, samplingStrategy=\"uniform\", condition=\"nTaxa\", taxa=taxa)\n
#Set up a starting value
timetree.setValue(T)\n
#Node age moves
moves.append( mvNodeTimeSlideUniform(timetree, weight=30.0) )
moves.append( mvSlide(root_age, delta=2.0, tune=true, weight=10.0) )
moves.append( mvScale(root_age, lambda=2.0, tune=true, weight=10.0) )
moves.append( mvTreeScale(tree=timetree, rootAge=root_age, delta=1.0, tune=true, weight=3.0) )\n
monitors.append(mnFile(filename=\"revbayes_mcmc.trees\", printgen=1, timetree))
mymodel = model(timetree)
mymcmc = mcmcmc(mymodel, monitors, moves, nruns=1, nchains=1, tuneHeat=TRUE)
mymcmc.burnin(generations=5000,tuningInterval=100)
mymcmc.run(generations=",nsim,")
q()")
    write(my_template,"input_revbayes.Rev")    
}

#BL-space generator
mix_beta=function(nsim)
    {
      #Set up progress bar
      cat("\nWARNING: This distribution can only be used for unrooted quartet trees.\nThis command may fail if grid cell is empty, re-run the command\n")    
      tree = read.tree(text= "(A,B,(C,D));")
      tree_list=rep(tree,nsim)
      v_select=combn(1:5,2)
      boot=c()
      space_t=matrix(sample(rbeta(10000000,c(0.1,0.5,1),c(0.1,0.5,1))),ncol=5)
      #AS assymetry score (= pairwise distance PD) NB neigbour sum (= sum of neighboring branches NS) + L tree length   
      AS=apply(space_t[,v_select[1,]]-space_t[,v_select[2,]],1,sum)
      LB=2*apply(space_t,1,sum)+apply(space_t[,2:4],1,sum)
      pdf(file="raw_BL_space.pdf")
      plot(AS[1:200000],LB[1:200000],col=adjustcolor("black", alpha.f = 0.05),pch=16,xlab="PD",ylab = "LNS")
      dev.off()
      x=seq(-6,6,0.1)
      y=seq(0,13,0.1)
      m=matrix(1:(length(x)*length(y)),nrow=length(x),ncol=length(y))
      xint=findInterval(AS,x)
      yint=findInterval(LB,y)
      all_t=data.frame(space_t,AS=AS,LB=LB,XI=xint,YI=yint,Fact=m[cbind(xint,yint)])
      #Uniform sampling from tree space
      cat("\nStarting uniform sampling from BL-space\n")
      pb1 = txtProgressBar(min = 0,     
                             max = 10, 
                             style = 3,    
                             width = 50,   
                             char = "=")
    
      
      for (i in 1:10)
      {  
        sampletree=data.frame(all_t %>% group_by(Fact) %>% sample_n(size = 1,replace=F))
        boot=rbind(boot,sampletree)
        setTxtProgressBar(pb1, i)  
      }
      close(pb1)
      boot=boot[sample(nrow(boot),nsim),]
      pdf(file="uniform_sample_from_BL_space.pdf")
      plot(boot$AS,boot$LB,col=adjustcolor("black", alpha.f = 0.05),pch=16,xlab="PD",ylab = "LNS")
      dev.off()
      write.table(data.frame(PD=boot$AS,LNS=boot$LB),"bl_coordinates.txt",row.names=F,quote=F)
      #Assign branch lengths to trees
      cat("\nInitializing the trees\n")
      pb2 = txtProgressBar(min = 0,     
                             max = nsim, 
                             style = 3,    
                             width = 50,   
                             char = "=")
      for ( i in 1:length(tree_list))
      {
        tree_list[[i]]$edge.length = as.numeric(boot[i,c("X1","X2","X5","X3","X4")])
        setTxtProgressBar(pb2, i)
      } 
      close(pb2)
      return(tree_list)
    }  


#Simulate branch lengths function
sim.brls=function(tree,nsim,distr)
{
    
    if (!any(distr[1] == c("bd", "mixb")))
    {    
        cat("\nInitializing the trees")
        n_br = nrow(tree$edge)
        tree$edge.length=rep(0,n_br)
        tree_list=rep(tree,nsim)
        pb = txtProgressBar(min = 0,     
                         max = nsim, 
                         style = 3,    
                         width = 50,   
                         char = "=")   
        if (distr[1] == "exp") 
        {
            #Example: exp,0.5
            rate = as.numeric(distr[2])
            for ( i in 1:length(tree_list))
            {
                tree_list[[i]]$edge.length = rexp(n_br, rate = rate )
                setTxtProgressBar(pb, i)
            }    
        }    
        if (distr[1] == "unif")
        {
            #Example: unif,3,10
            lower = as.numeric(distr[2])
            upper = as.numeric(distr[3])
            for ( i in 1:length(tree_list))
            {
                tree_list[[i]]$edge.length = runif(n_br, min = lower, max = upper)
                setTxtProgressBar(pb, i)
            }    
        }
        if (distr[1] == "dir") 
        {
            #Example: dir,0.5
            #Concentration parameters
            alphas = rep(as.numeric(distr[2]),n_br)
            for ( i in 1:length(tree_list))
            {
                tree_list[[i]]$edge.length = rdirichlet(1, alpha = alphas)[1,]
                setTxtProgressBar(pb, i)
            }    
        }
        close(pb)
    }   
    
    if (distr[1] == "bd")
    {
        cat("\nStarting RevBayes to simulate BD trees and re-scale their branch lengths with NELSI")
        #Example: bd,5,1,200,0.01 lambda mu root_age clock_rate
        get_revbayes_in(birth_rate = as.numeric(distr[2]),death_rate = as.numeric(distr[3]),root_age=as.numeric(distr[4]),nsim = nsim)
        system("/Users/anton/Downloads/rb input_revbayes.Rev")
        rev_trees = read.tree("revbayes_mcmc.trees")
        rev_trees = rev_trees[2:length(rev_trees)]
        clock_scaled = lapply(rev_trees,simulate.clock,params = list(rate = as.numeric(distr[5]), noise = 0))
        #Dummy list 
        clock_scaled_unlisted = (unlist(lapply(clock_scaled,"[",1),recursive = F,use.names = F))
        tree_list=read.tree(text="(Dummy1,Dummy2);")
        tree_list=c(tree_list,tree_list)
        for (i in 1:length(clock_scaled))
        {
           tree_list[i]=clock_scaled_unlisted[i]
        }    
    }
    if (distr[1] == "mixb")
    {
       #Example: mixb
       tree_list=mix_beta(nsim) 
    }    
    

    return(tree_list)  
}

#Generate partition file for IQTREE simulatior 
get.nexus.part=function(nsim,aln_len,file)
{
    file = paste(file,".part",sep="")
    write("#nexus\nbegin sets;",file)
    write.table(data.frame('\tcharset',
                           paste("locus",1:nsim,sep="_"),
                           "=",
                           paste("DNA",",",sep=""),
                           paste(seq(1,aln_len*nsim,by=aln_len),
                               "-",
                                 paste(seq(aln_len,aln_len*nsim,by=aln_len),";",sep=""),
                           sep="")),
                           file,
                           append = T,
                           quote = F,
                           row.names = F,
                           col.names = F)
    write("end;",file,append=T)
}    


main_gen=function(nwk_tree,distr,nsim,aln_len,file)
{
   topo = nwk_tree
   trees = sim.brls(topo,nsim,distr)
   cat("\nGetting nexus partition file for IQTREE")
   get.nexus.part(nsim,aln_len,file)  
   cat("\nSaving trees and their branch lengths")
   write.table(unlist(lapply(trees,write.tree)),paste(file,".tre",sep=""),quote = F, row.names = F, col.names = F)
   write.table(t(data.frame(lapply(trees,"[[","edge.length"))),paste(file,".Y.txt",sep=""),quote = F, row.names = F, col.names = F)
   if (is.rooted(topo))
   {
      cat("\nRooted trees detected! Unrooting")
      trees_unroot = lapply(trees,unroot)
      write.table(unlist(lapply(trees_unroot,write.tree)),paste(file,".unrooted.tre",sep=""),quote = F, row.names = F, col.names = F)
      cat("\nSaving urooted trees and their branch lengths")
      write.table(t(data.frame(lapply(trees_unroot,"[[","edge.length"))),paste(file,".unrooted.Y.txt",sep=""),quote = F, row.names = F, col.names = F) 
      
   }    
   cat("\nDone!")    
}    



#Save input tree
if (strsplit(opt$topo,split="")[[1]][1] == "(")
{
    tree_string = read.tree(text = opt$topo)
}else{
    tree_string = read.tree(opt$topo)
}    

n_taxa = length(tree_string$tip.label)
#Create a directory structure
main_dir = paste("experiment_",paste(strsplit(opt$distribution,",")[[1]],collapse="_"),"_",n_taxa,"taxa",sep="")
dir.create(main_dir)
setwd(main_dir)
#Write tree in newick
write.tree(tree_string,"input_tree.tre")
write(my_seed,"seed.txt")
cat("\nGenerating TEST dataset")
dir.create("TEST")
setwd("TEST")

main_gen(nwk_tree = tree_string,
         distr = strsplit(opt$distribution,",")[[1]],
         nsim = round(opt$nsim*opt$prop),
         aln_len = opt$len,
         file = paste("test_",opt$fname,sep=""))
setwd("..")
cat("\nGenerating TRAIN dataset")
dir.create("TRAIN")
setwd("TRAIN")
  
main_gen(nwk_tree = tree_string,
         distr = strsplit(opt$distribution,",")[[1]],
         nsim = opt$nsim,
         aln_len = opt$len,
         file = paste("train_",opt$fname,sep=""))
setwd("../..")




