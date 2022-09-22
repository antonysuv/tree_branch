library("ape")
library("ggpubr")
library("reshape")


#Polygon t1 = estimated branch lengths, t2 = true branch lengths 
get_oct=function(t1,t2,bl_coord,param)
{
    ae = abs(t1-t2)
    names(ae)=paste("B",c(1,2,5,3,4),sep="")
    ae$mae=apply(ae,1,mean)
    ae$error=apply(ae[,1:5]/t2,1,mean)
    ae=cbind(ae,bl_coord)
    polyxy = data.frame(x=c(c(-6,-6,6,6),c(0,-4,-6,-6,-4,0,4,6,6,4,0)),y=c(c(0,13,13,0),c(0,2,5,8,11,13,11,8,5,2,0)),subid=c(rep(1,4),rep(2,11)))
    lower5 = sort(ae[,param])[round(nrow(ae)*0.05)]
    upper5 = sort(ae[,param],decreasing = TRUE)[round(nrow(ae)*0.05)]
    ae_lower = ae[ae[,param]<lower5,]
    ae_upper = ae[ae[,param]>upper5,]
    ae_lower$id = paste(param,"lower 5%")
    ae_upper$id = paste(param,"upper 5%")
    ae_combined = rbind(ae_lower,ae_upper)
    
    pl = ggplot(ae_combined, aes(x=PD, y=LNS))+
    stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE)+
    geom_polygon(data=polyxy, aes(x = x, y = y,subgroup = subid),color="white",fill="white")+scale_fill_viridis_c(option = "magma")
    pl = facet(pl,facet.by = "id",scales = "free",ncol = 2)  
    return(pl)
}    



polyxy = data.frame(x=c(c(-6,-6,6,6),c(0,-4,-6,-6,-4,0,4,6,6,4,0)),y=c(c(0,13,13,0),c(0,2,5,8,11,13,11,8,5,2,0)),subid=c(rep(1,4),rep(2,11)))
ggplot(polyxy, aes(x = x, y = y))+geom_polygon(aes(subgroup = subid),color="white",fill="white")
#geom_polygon(data=polyxy, aes(x = x, y = y,subgroup = subid),color="white",fill="white")

ggplot(gg[gg$gr<0.1,], aes(x=PD, y=LNS))+stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE)+geom_polygon(data=polyxy, aes(x = x, y = y,subgroup = subid),color="white",fill="white")+scale_fill_viridis_c(option = "magma")
ggplot(gg[gg$mlerr<0.1,], aes(x=PD, y=LNS))+stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE)+geom_polygon(data=polyxy, aes(x = x, y = y,subgroup = subid),color="white",fill="white")+scale_fill_viridis_c(option = "magma")


#Extract branch lenghts from a file with newick trees
bls_extract = function(file_in)
{
    phy = read.tree(file_in)
    brls = data.frame(t(data.frame(lapply(phy,"[[","edge.length")))) 
    row.names(brls) = NULL
    return(brls)
}    

#Tables with branch lengths 
data_prep=function(t1,t2)
{
    t1$trl = apply(t1,1,sum)
    t2$trl = apply(t2,1,sum)
    t1 = melt(t1)
    t2 = melt(t2)
    t1_t2 = data.frame(br_id = t1$variable, estimated_brl = t1$value, true_brl = t2$value)
    return(t1_t2)    
}    


# t1 = estimated, t2 = true 
get_cor = function(t1,t2) 
{


    pl = data_prep(t1,t2)
    pl$er = (pl$estimated_brl-pl$true_brl)^2
    pl$aer = abs(pl$estimated_brl-pl$true_brl)
    mses = paste("MSE =",round(tapply(pl$er,pl$br_id,mean),6))
    maes = paste("\nMAE =",round(tapply(pl$aer,pl$br_id,mean),6))
    mses_labels = paste(c(rep("Branch",length(mses)-1),"Tree length"),c(as.character(1:(length(mses)-1)),""),":",mses,maes)
    pl_g = ggscatter(pl, x = "estimated_brl", y = "true_brl",
          add = "reg.line",                                 # Add regression line
          conf.int = TRUE,                                  # Add confidence interval
          size = 0.5,
          xlim=c(0,0.2),
          ylim=c(0,0.2),
          add.params = list(color = "blue",
                            fill = "lightgray")
                            )+
    #+stat_cor(method = "pearson")
    geom_abline(slope = 1,color = "red")+
    xlab("Estimated branch lengths")+ylab("True branch lengths")
  
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
          #xlim=c(0,2),
          #ylim=c(0,2),
          add.params = list(color = "blue",
                            fill = "lightgray")
                            )+
    #+stat_cor(method = "pearson")
    +geom_density_2d_filled(contour_var = "ndensity",alpha = 0.6)+
    +geom_abline(slope = 1,color = "red",size = 0.1)+
    +xlab("Estimated branch lengths")+ylab("True branch lengths")
  
  pl_g = facet(pl_g,facet.by = "br_id",scales = "free",panel.labs = list(br_id=mses_labels),ncol = 4)  
  return(pl_g)
}

   + geom_density_2d_filled(contour_var = "ndensity")

geom_density_2d_filled(contour_var = "ndensity",alpha = 0.6)









           
get_viol = function(t1,t2)
{    
    pl = data_prep(t1,t2)
    pl_m = melt(pl)
    pl_g = ggviolin(pl_m, x = "variable", y = "value", color = "variable",add = "boxplot",palette = c("black","orange"),add.params = list(fill = "white"),ylim=c(0,10))+
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

