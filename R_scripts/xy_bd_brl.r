library('TreeSim')
library('NELSI')
library("MCMCpack")
library("optparse")
options(scipen = 100000000)

option_list = list(
  make_option(c("-n", "--n_tips"), type="numeric", default=NULL,help="number of tips on BD tree", metavar="numeric"),
  make_option(c("-b", "--birth_lambda"), type="numeric", default=NULL, help="birth rate lambda", metavar="numeric"),
  make_option(c("-d", "--death_mu"), type="numeric", default=NULL, help="death rate mu", metavar="numeric"),  
  make_option(c("-e", "--seed1"), type="numeric", default=NULL, help="seed for BD tree generation", metavar="numeric"),
  make_option(c("-c", "--n_calib_nodes"), type="numeric", default=NULL, help="number of random calibration nodes", metavar="numeric"),
  make_option(c("-f", "--ditance_fraction"), type="numeric", default=NULL, help="distance from true nodal age (upper and lower calibration bounds)", metavar="numeric"),
  make_option(c("-z", "--seed2"), type="numeric", default=NULL, help="seed for random calibration points", metavar="numeric"),  
  make_option(c("-r", "--clock_rate"), type="numeric", default=NULL, help="strict clock rate", metavar="numeric"),
  make_option(c("-m", "--n_mcmc_draws"), type="numeric", default=NULL, help="total number of MCMC draws", metavar="numeric"),
  make_option(c("-s", "--n_log_every"), type="numeric", default=NULL, help="sampling frequency from MCMC", metavar="numeric"),
  make_option(c("-a", "--aln_length"), type="numeric", default=NULL, help="MSA length", metavar="numeric"),  
  make_option(c("-x", "--xml_name"), type="character", default=NULL, help="input XML file name for BEAST2", metavar="character"),
  make_option(c("-t", "--tmplate_path"), type="character", default=NULL, help="path to BEAST2 XML template", metavar="character"),
  make_option(c("-p", "--beast_path"), type="character", default=NULL, help="path to BEAST2 executable", metavar="character")  
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


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


#Generate BD tree
bd_tree=function(n_tips,birth_lambda,death_mu,seed_n)
{
    set.seed(seed_n)
    phy = sim.bd.taxa(n = n_tips, numbsim = 1, lambda = birth_lambda, mu = death_mu, complete=F)[[1]]
    #phy$root.edge = NULL
    write.tree(phy,paste("main_phy_",birth_lambda,"_",death_mu,"_",seed_n,".tre",sep=""))
    return(phy)
} 

#MSA Block for BEAST2 
aln_block_beast2 = function(labls)
{
    dummy_data = paste('\t\t\t<sequence id="seq_',labls,'" spec=\"Sequence\" taxon=\"',labls,'\" totalcount="4" value="A"/>\n',sep="",collapse=" ")
    return(dummy_data)
}    

#Taxa Block for BEAST2 
taxonset_block_beast2 = function(labls)
{
    start_line = '<distribution id="all.prior" spec="beast.math.distributions.MRCAPrior" monophyletic="true" tree="@Tree.t:dummy">\n\t\t<taxonset id="all" spec="TaxonSet">\n'
    taxa = paste('\t\t\t<taxon id="',labls,'" spec="Taxon"/>\n',sep="",collapse=" ")
    end_line = '\t\t</taxonset>\n\t    </distribution>'
    taxa_block = paste(start_line,taxa,end_line,sep="",collapse=" ")
    return(taxa_block)
}    

#Fossil calibration Block for BEAST2 
prior_block_beast2 = function(labls,id_prior,lower_bound,upper_bound)
{
    start_line = paste('<distribution id="',id_prior,'.prior" spec="beast.math.distributions.MRCAPrior" monophyletic="true" tree="@Tree.t:dummy">\n',sep="",collapse=" ")
    start_line_2 = paste('\t\t<taxonset id="',id_prior,'" spec="TaxonSet">\n',sep="",collapse=" ")
    taxa_mrca = paste('\t\t\t<taxon idref="',labls,'"/>\n',sep="",collapse=" ")
    end_line = '\t\t</taxonset>\n'
    end_line_2 = paste('\t\t<Uniform id="Uniform.',id_prior,'" lower="',lower_bound,'" name="distr" upper="',upper_bound,'"/>\n\t    </distribution>',sep="",collapse=" ")
    prior_block = paste(start_line,start_line_2,taxa_mrca,end_line,end_line_2,sep="",collapse=" ")
    id_prior_logger = paste('<log idref="',id_prior,'.prior"/>\n',sep="",collapse=" ")
    prior_list=list(priors = prior_block,prior_ids = id_prior_logger)
    return(prior_list)           
}    

#Generate a random fossil calibration blocks plus a calibration block for root
draw_priors_beast2 = function(n,phy,ditance_fraction,seed_n)
{
    set.seed(seed_n)
    node_age = branching.times(phy)
    root_age = node_age[1]
    root_node = names(node_age[1])
    root_taxa = extract.clade(phy,node = as.numeric(root_node))$tip.label
    prior_block = prior_block_beast2(root_taxa,"root",root_age-root_age*ditance_fraction,root_age+root_age*ditance_fraction)
    
    if (n > 1)
    {
        node_list = sample(node_age[2:length(node_age)], size = n-1)
        for (i in 1:length(node_list))
        {
            node = node_list[i]
            node_taxa = extract.clade(phy,node = as.numeric(names(node)))$tip.label
            node_block = prior_block_beast2(node_taxa,paste("node_",names(node),sep=""),node-node*ditance_fraction,node+node*ditance_fraction)
            prior_block$priors = paste(prior_block$priors,node_block$priors,collapce="\n")
            prior_block$prior_ids = paste(prior_block$prior_ids,node_block$prior_ids,collapse="\n")       
        }
        write(c(root_age,node_list),"calibration_nodes.txt")
    }else{
        write(cbind(root_age,root_age-root_age*ditance_fraction,root_age+root_age*ditance_fraction),"calibration_nodes.txt")
    }
    
    return(prior_block)
}    

#Get XML input for BEAST2
get_xml_beast2 = function(tmplate_path,phy,birth_lambda,death_mu,list_age_priors,n_mcmc_draws,n_log_every,xml_name)
{
    my_template = readLines(tmplate_path)
    tip_labels = phy$tip.label
    taxalist=taxonset_block_beast2(tip_labels)
    phy_constrain_topology = phy
    phy_constrain_topology$edge.length = NULL
    phy_constrain_topology = write.tree(phy_constrain_topology)
    phy_starting_tree = write.tree(phy)
    dymmy_aln = aln_block_beast2(tip_labels)
    birth_rate = birth_lambda - death_mu
    death_rate = death_mu/birth_lambda
    
    my_template = gsub("DUMMYDATA",dymmy_aln,my_template)
    my_template = gsub("N_MCMC",n_mcmc_draws,my_template)
    my_template = gsub("N_STORE_MCMC",n_log_every,my_template)
    my_template = gsub("DEATHRATE",death_rate,my_template)
    my_template = gsub("BIRTHRATE",birth_rate,my_template)
    my_template = gsub("FIXEDTREETOPO",phy_constrain_topology,my_template)
    my_template = gsub("STARTTREEWITHBL",phy_starting_tree,my_template)
    my_template = gsub("TAXA_LIST",taxalist,my_template)
    my_template = gsub("PRIOR_LIST",list_age_priors$priors,my_template)
    my_template = gsub("ID_PRIOR",list_age_priors$prior_ids,my_template)
    write(my_template,xml_name)
}    
    
#Get nodal ages Y input for neural net
get_nodal_ages_table = function(phy_prior)
{
    age_table = t(data.frame(lapply(phy_prior,branching.times)))
    write.table(age_table,"Y_ages.txt",row.names=F,quote=F,col.names=F)
}    

#Get clock-scaled branch lengths and trees. Input for MSA simulatior INDELIBLE
get_clock_bls = function(phy_prior,clock_rate)
{
    clock_scaled = lapply(phy_prior,simulate.clock,params = list(rate = clock_rate, noise = 0))
    clock_scaled_newick = unlist(lapply(lapply(clock_scaled,"[[",1),write.tree))
    bl_table = t(data.frame(lapply(lapply(clock_scaled,"[[",2),function(x){col = x[,"length.subst"]; return(col)})))
    write.table(bl_table,"Y_brl.txt",row.names = F,quote = F,col.names = F)
    return(clock_scaled_newick)
}    

#Substitution model block for INDELIBLE control file
model_GTR_indelible = function(n_model_blocks,file)
{
  write(paste('[TYPE] NUCLEOTIDE 2\n[SETTINGS]\n [output] FASTA\n [randomseed] ',round(runif(1,1,100000))),file)
  model_id = 1
  for (block in 1:n_model_blocks)
  {
    #Gamma rate
    A = runif(1,0,5)
    #Nucl proportions DIRICHLET 
    options(digits = 5)
    Pi = format(rdirichlet(1, alpha=c(5,5,5,5)))
    
    write(paste('\n[MODEL] ','GTRModel',model_id,sep = ''),file,append=T)
    options(digits=2)
    model=paste(c('GTR ',format(runif(5,0,3))),sep = '')
    model_id=model_id+1
    write(paste(' [submodel] ',paste(model,collapse = ' '),'\n [rates] 0 ',A,' 0'),file,append=T)
    write(paste(' [statefreq]',paste(Pi,collapse = ' ')),file,append=T)  
    
  }
}

#TREE, PARTITIONS and EVOLVE blocks for INDELIBLE control file
tpe_indelible = function(newick_trees,aln_length,file)
{
    id_tree = paste("t",rep("_sim",times = length(newick_trees)),1:length(newick_trees),sep = "")
    write.table(data.frame('[TREE]',id_tree,newick_trees),file,append = T,quote = F,row.names = F,col.names = F)
    p_name = paste("p",1:length(newick_trees),sep = "")
    model_id = paste("GTRModel",1:length(newick_trees),sep = "")
    write.table(data.frame('[PARTITIONS]',p_name,"[",id_tree,model_id,aln_length,"]"),file,append = T,quote = F,row.names = F,col.names = F)
    #Set EVOLVE block
    write('[EVOLVE]',file,append = T)
    #Write INDELIBLE control file
    write.table(data.frame(p_name,1,apply(data.frame(id_tree,"_",model_id),1,paste,collapse="")),file,append = T,quote = F,row.names = F,col.names = F)
}    

main_gen = function(n_tips,birth_lambda,death_mu,seed_n1,n_calib_nodes,ditance_fraction,seed_n2,tmplate_path,n_mcmc_draws,n_log_every,xml_name,clock_rate,aln_length,beast_path)
{
    phy = bd_tree(n_tips,birth_lambda,death_mu,seed_n1)
    list_age_priors = draw_priors_beast2(n_calib_nodes,phy,ditance_fraction,seed_n2)
    get_xml_beast2(tmplate_path,phy,birth_lambda,death_mu,list_age_priors,n_mcmc_draws,n_log_every,xml_name)
    system(paste(beast_path,xml_name))
    phy_prior = read.nexus("ultrametric.trees")
    #Remove 0-length root edge 
    for (i in 1:length(phy_prior)){phy_prior[[i]]$root.edge = NULL}
    #Remove the first extra tree generated by BEAST2 
    phy_prior = phy_prior[2:length(phy_prior)]
    get_nodal_ages_table(phy_prior)
    trees_for_indelible=get_clock_bls(phy_prior,clock_rate)
    model_GTR_indelible(length(trees_for_indelible),"control.txt")
    tpe_indelible(trees_for_indelible,aln_length,"control.txt")
    
}

main_gen(n_tips = opt$n_tips,
         birth_lambda = opt$birth_lambda,
         death_mu = opt$death_mu,
         seed_n1 = opt$seed1,
         n_calib_nodes = opt$n_calib_nodes,
         ditance_fraction = opt$ditance_fraction,
         seed_n2 = opt$seed2,
         tmplate_path = opt$tmplate_path,
         n_mcmc_draws = opt$n_mcmc_draws,
         n_log_every = opt$n_log_every,
         xml_name = opt$xml_name,
         clock_rate = opt$clock_rate,
         aln_length = opt$aln_length,
         beast_path = opt$beast_path
        )
