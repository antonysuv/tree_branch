library("ape")
library("ggpubr")
library("reshape")


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
data_prep=function(t1,t2)
{
    t1t2cl=clean_outs(t1,t2)
    t1=t1t2cl[[1]]
    t2=t1t2cl[[2]]
    t1$trl = apply(t1,1,sum)
    t2$trl = apply(t2,1,sum)
    t1 = melt(t1)
    t2 = melt(t2)
    t1_t2 = data.frame(br_id = t1$variable, estimated_brl = t1$value, true_brl = t2$value)
    return(t1_t2)    
}    

# t1 = estimated, t2 = true 
get_hist = function(t1,t2)
{
    pl = data_prep(t1,t2)
    pl_m = melt(pl)
    n_lab = length(unique((pl_m$br_id)))
    facet_labs = paste(c(rep("Branch",n_lab-1),"Tree length"),c(as.character(1:(n_lab-1)),""))
    p = gghistogram(pl_m, x = "value",y="..density..",fill = "variable",bins =50,palette = c("gold", "navy"),facet.by="br_id",scales = "free",add="median",add_density = TRUE, xlab = "branch length",panel.labs = list(br_id= facet_labs))
    ggsave(plot = p, width = 15, height = 15, dpi = 300, filename = "hist_brs.pdf")
}
