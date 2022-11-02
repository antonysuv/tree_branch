library("ape")
library("ggpubr")
library("reshape")
library("optparse")
library("plyr")
options(scipen = 100000000)
#Caution: warnings are disabled 
options(warn=-1)


option_list = list(
  make_option(c("-m", "--ml"), type="character", default="TEST/ML_all.tre",help="path to ML branch lengths", metavar="character"),
  make_option(c("-c", "--cnn"), type="character", default="brls.predicted.cnn.reg.sqrt.txt", help="path to CNN branch lengths", metavar="character"),  
  make_option(c("-p", "--mlp"), type="character", default="brls.predicted.mlp.reg.sqrt.txt", help="path to MLP branch lengths", metavar="character"),
  make_option(c("-s", "--cnn_noreg"), type="character", default="brls.predicted.cnn.sqrt.txt", help="path to CNN branch lengths", metavar="character"),  
  make_option(c("-r", "--mlp_noreg"), type="character", default="brls.predicted.mlp.sqrt.txt", help="path to MLP branch lengths", metavar="character"),
  make_option(c("-t", "--truebl"), type="character", default="TEST/test_simulation.Y.txt", help="path to true branch lengths",metavar="character"),
  make_option(c("-d", "--dir"), type="character", default="plots_blspace", help="create directory for outputs",metavar="character"),
  make_option(c("-b", "--blspace"), type="character", default="TEST/bl_coordinates.txt", help="path to BL space coordinates",metavar="character")  
)  

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)


########## Utility functions ##########
#Extract branch lengths from a set of newick trees
bls_extract = function(file_in)
{
    phy = read.tree(file_in)
    brls = data.frame(t(data.frame(lapply(phy,"[[","edge.length")))) 
    row.names(brls) = NULL
    
    return(brls)
}    



#Prepare tables with branch lengths 
#t1 = estimated, t2 = true bl_coor = branch length space coordinates output from alisim_input_generate.r mixb function 
data_prep=function(t1,t2,bl_coor, methodname = "Method")
{
  
    bl_coor$MSE = apply((t2-t1)^2,1,mean)
    bl_coor$MAE = apply(abs(t2-t1),1,mean)
    bl_coor$Method = methodname
    return(bl_coor)

}  

prep_all=function(cnn,cnn_roe,mlp,mlp_roe,ml,true_brl,bl_coor)
{
    d_cnn = data_prep(cnn,true_brl,bl_coor,methodname = "CNN") 
    d_cnn_roe = data_prep(cnn_roe,true_brl,bl_coor,methodname = "CNN-ROE")
    d_mlp= data_prep(mlp,true_brl,bl_coor,methodname = "MLP") 
    d_mlp_roe = data_prep(mlp_roe,true_brl,bl_coor,methodname = "MLP-ROE")
    d_ml = data_prep(ml,true_brl,bl_coor,methodname = "ML")
    d_master = rbind(d_cnn,d_cnn_roe,d_mlp,d_mlp_roe,d_ml)
    #Plotting order consistent CNN,CNN-ROE,MLP,MLP-ROE,ML
    d_master$Method=factor(d_master$Method,levels=c("CNN","CNN-ROE","MLP","MLP-ROE","ML"))
    return(d_master)
}    


get_blspace=function(t1,stat,minq,maxq,dir_name)
{
    
    
    my_q = quantile(t1[,stat],probs=seq(0,1,0.01))
    t1$Quantile_range = ifelse(t1[,stat]<my_q[minq],paste(stat,"<",minq),ifelse(t1[,stat]>my_q[maxq],paste(stat,">",maxq),paste(minq,"<",stat,"<",maxq)))
    polyxy = data.frame(x=c(c(-6,-6,6,6),c(0,-4,-6,-6,-4,0,4,6,6,4,0)),y=c(c(0,13,13,0),c(0,2,5,8,11,13,11,8,5,2,0)),subid=c(rep(1,4),rep(2,11)))
    
    for (i in unique(t1$Quantile_range))
    {    
        g = gsub("%","",gsub(">","greater",gsub(" ","",gsub("<","smaller",i))))
        p1 = ggplot(t1[t1$Quantile_range==i,],aes(x=PD,y=LNS))+
             stat_density_2d(geom = "raster",aes(fill = after_stat(density)),contour = FALSE,scales="free")+
             facet_wrap("Method",ncol=5,nrow=1)+
             scale_fill_viridis_c(option = "magma")+
             geom_polygon(data=polyxy, aes(x = x, y = y,subgroup = subid),color="white",fill="white")+
             ggtitle(paste("Quantile range:",i))
        ggsave(plot = p1, width = 10, height = 2.5, dpi = 300, filename = paste(dir_name,"/",g,".pdf",sep="",collapse=""))
        
    }
}    

                              
dir.create(opt$dir)
                     
cat("Read true\n")
true_brl = read.table(opt$truebl)  
cat("Read ml\n")
ml = bls_extract(opt$ml) 
cat("Read cnn ROE\n")
cnn_roe = read.table(opt$cnn)
cat("Read mlp ROE\n")
mlp_roe = read.table(opt$mlp)
cat("Read mlp\n")
mlp = read.table(opt$mlp_noreg)
cat("Read cnn\n")
cnn = read.table(opt$cnn_noreg)
cat("Read BL table\n")
bl_coord=read.table(opt$blspace,header=T)

d_master = prep_all(cnn,cnn_roe,mlp,mlp_roe,ml,true_brl,bl_coord) 


get_blspace(d_master,"MSE","25%","75%",opt$dir)
get_blspace(d_master,"MAE","25%","75%",opt$dir)
