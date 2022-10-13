library("ape")
library("ggpubr")
library("reshape")
library("optparse")


option_list = list(
  make_option(c("-m", "--ml"), type="character", default="TEST/ML_all.tre",help="path to ML branch lengths", metavar="character"),
  make_option(c("-c", "--cnn"), type="character", default="brls.predicted.cnn.sqrt.txt", help="path to CNN branch lengths", metavar="character"),  
  make_option(c("-p", "--mlp"), type="character", default="brls.predicted.mlp.reg.sqrt.txt", help="path to MLP branch lengths", metavar="character"),
  make_option(c("-t", "--truebl"), type="character", default="TEST/test_simulation.Y.txt", help="path to true branch lengths",metavar="character") 
)  

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


bls_extract = function(file_in)
{
    phy = read.tree(file_in)
    brls = data.frame(t(data.frame(lapply(phy,"[[","edge.length")))) 
    row.names(brls) = NULL
    
    return(brls)
}    


clean_outs=function(t1,t2)
{
    outs=c()
    for (i in 1:ncol(t1))
    {
        if (length(boxplot(t1[,i],plot=FALSE)$out) != 0)
        {
            outs = c(outs,-which(t1[,i] %in% boxplot(t1[,i],plot=FALSE)$out))
        }    
    }
    if (!is.null(outs))
    {
        t1=t1[outs,]
        t2=t2[outs,]
        return(list(t1,t2))
    }else{
        return(list(t1,t2))
    }       
}


# t1 = estimated, t2 = true 
#Tables with branch lengths 
data_prep=function(t1,t2,cl_out = FALSE,filename = "hist_bls.pdf" )
{
    if (cl_out == TRUE)
    {    
        t1t2cl=clean_outs(t1,t2)
        t1=t1t2cl[[1]]
        t2=t1t2cl[[2]]
    }    
    t1$trl = apply(t1,1,sum)
    t2$trl = apply(t2,1,sum)
    t1 = melt(t1)
    t2 = melt(t2)
    t1_t2 = data.frame(br_id = t1$variable, estimated_brl = t1$value, true_brl = t2$value)
    return(t1_t2)    
}    

# t1 = estimated, t2 = true 
get_hist = function(t1,t2,cl_out = FALSE,filename = "hist_bls.pdf")
{
    pl = data_prep(t1,t2,cl_out)
    pl_m = melt(pl)
    n_lab = length(unique((pl_m$br_id)))
    facet_labs = paste(c(rep("Branch",n_lab-1),"Tree length"),c(as.character(1:(n_lab-1)),""))
    p = gghistogram(pl_m, x = "value",y="..density..",fill = "variable",bins =50,palette = c("gold", "navy"),facet.by="br_id",scales = "free",add="median",add_density = TRUE, xlab = "branch length",panel.labs = list(br_id= facet_labs))
    ggsave(plot = p, width = 15, height = 15, dpi = 300, filename = filename)
}


# t1 = estimated, t2 = true 
get_cor = function(t1,t2,cl_out = FALSE,filename = "dot_bls.pdf") 
{

    pl = data_prep(t1,t2,cl_out)
    pl$er = (pl$estimated_brl-pl$true_brl)^2
    pl$aer = abs(pl$estimated_brl-pl$true_brl)
    mses = paste("MSE =",round(tapply(pl$er,pl$br_id,mean),6))
    maes = paste("\nMAE =",round(tapply(pl$aer,pl$br_id,mean),6))
    mses_labels = paste(c(rep("Branch",length(mses)-1),"Tree length"),c(as.character(1:(length(mses)-1)),""),":",mses,maes)
    pl_g = ggscatter(pl, x = "estimated_brl", y = "true_brl",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          size = 0.5,
          #xlim=c(0,0.2),
          #ylim=c(0,0.2),
          add.params = list(color = "blue",
                            fill = "lightgray")
                            )+
    #+stat_cor(method = "pearson")
    geom_abline(slope = 1,color = "red")+
    xlab("Estimated branch lengths")+ylab("True branch lengths")
  
  pl_g = facet(pl_g,facet.by = "br_id",scales = "free",panel.labs = list(br_id=mses_labels),ncol = 4)  
  ggsave(plot = pl_g, width = 15, height = 15, dpi = 300, filename = filename)
}

tt=read.table(opt$truebl)  
ml = bls_extract(opt$ml) 
cnn = read.table(opt$cnn)
mlp = read.table(opt$mlp)


get_cor(ml,tt,filename = "plot_ml_dot_raw.pdf")
get_cor(cnn,tt,filename = "plot_cnn_dot_raw.pdf")
get_cor(mlp,tt,filename = "plot_mlp_dot_raw.pdf")

get_hist(ml,tt,filename = "plot_ml_hist_raw.pdf")
get_hist(cnn,tt,filename = "plot_cnn_hist_raw.pdf" )
get_hist(mlp,tt,filename = "plot_mlp_hist_raw.pdf")

get_cor(ml,tt,cl_out = TRUE,filename = "plot_ml_dot_clean.pdf")
get_cor(cnn,tt,cl_out = TRUE,filename = "plot_cnn_dot_clean.pdf")
get_cor(mlp,tt,cl_out = TRUE,filename = "plot_mlp_dot_clean.pdf")

get_hist(ml,tt,cl_out = TRUE,filename = "plot_ml_hist_clean.pdf")
get_hist(cnn,tt,cl_out = TRUE,filename = "plot_cnn_hist_clean.pdf" )
get_hist(mlp,tt,cl_out = TRUE,filename = "plot_mlp_hist_clean.pdf")