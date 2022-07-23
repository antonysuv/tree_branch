library("ape")
library("ggpubr")
library("reshape")


#Extract branch lenghts from a file with newick trees
bls_extract = function(file_in)
{
    phy = read.tree(file_in)
    brls = data.frame(t(data.frame(lapply(phy,"[[","edge.length")))) 
    row.names(brls) = NULL
    return(brls)
}    


data_prep=function(t1,t2)
{
    t1$trl = apply(t1,1,sum)
    t2$trl = apply(t2,1,sum)
    t1 = melt(t1)
    t2 = melt(t2)
    t1_t2 = data.frame(br_id = t1$variable, estimated_brl = t1$value, true_brl = t2$value)
    return(t1_t2)    
}    



get_cor = function(t1,t2) 
{


    pl = data_prep(t1,t2)
    pl$er = (pl$estimated_brl-pl$true_brl)^2
    mses = paste("MSE =",round(tapply(pl$er,pl$br_id,mean),6))
    mses_labels = paste("Branch",1:length(mses),":",mses)
    pl_g = ggscatter(pl, x = "estimated_brl", y = "true_brl",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          size = 0.5,
          #xlim=c(0,1),
          #ylim=c(0,1),
          add.params = list(color = "blue",
                            fill = "lightgray")
                            )+
  stat_cor(method = "pearson")+geom_abline(slope = 1,color = "red")+xlab("Estimated branch lengths")+ylab("True branch lengths")
  
  pl_g = facet(pl_g,facet.by = "br_id",scales = "free",panel.labs = list(br_id=mses_labels),ncol = 4)  
  return(pl_g)
}

   + geom_density_2d_filled(contour_var = "ndensity")


get_cor = function(t1,t2) 
{


    pl = data_prep(t1,t2)
    pl$er = (pl$estimated_brl-pl$true_brl)^2
    mses = paste("MSE =",round(tapply(pl$er,pl$br_id,mean),6))
    mses_labels = paste("Branch",1:length(mses),":",mses)
    pl_g = ggscatter(pl, x = "estimated_brl", y = "true_brl",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                 # Add confidence interval
          size = 0.1,
          #xlim=c(0,1),
          #ylim=c(0,1),
          add.params = list(color = "blue",
                            fill = "lightgray")
                            )+
  stat_cor(method = "pearson")+geom_density_2d_filled(contour_var = "ndensity",alpha = 0.6)+geom_abline(slope = 1,color = "red",size = 0.1)+xlab("Estimated branch lengths")+ylab("True branch lengths")
  
  pl_g = facet(pl_g,facet.by = "br_id",scales = "free",panel.labs = list(br_id=mses_labels),ncol = 4)  
  return(pl_g)
}

   + geom_density_2d_filled(contour_var = "ndensity")

geom_density_2d_filled(contour_var = "ndensity",alpha = 0.6)









           
get_viol = function(t1,t2)
{    
    pl = data_prep(t1,t2)
    pl_m = melt(pl)
    pl_g = ggviolin(pl_m, x = "variable", y = "value", color = "variable",add = "boxplot",palette = c("black","orange"),add.params = list(fill = "white"))+
    stat_compare_means(aes(label = paste0("p = ", ..p.format..)))
    pl_g = facet(pl_g,facet.by = "br_id",scales = "free",ncol = 4) 
    return(pl_g)
}




get_main=function(t1,t2)
{
    pl1 = get_cor(t1,t2)
    pl2 = get_viol(t1,t2)
    figure = ggarrange(pl1,pl2,
                    labels = c("A", "B"),
                    ncol = 1, nrow = 2)
    return(figure)

}    


figure <- ggarrange(bxp, dp, lp,
                    labels = c("A", "B", "C"),
                    ncol = 2, nrow = 2)












#Error measures
mae = function(x,y)
{
    return(mean(abs(x-y)))
}

mse = function(x,y)
{
    return(mean((x-y)^2))
}

rmsd = function(x,y)
{
    return(sqrt(mse(x,y)))
}

get_corr_plot = function(a,b,tl)
{
    corr_reults = cor.test(a,b)
    plot(a,b,xlab = "True branch lengths", "Predicted branch lengths", main = tl)
    abline(c(0,1),col="red")
    
    
}    

