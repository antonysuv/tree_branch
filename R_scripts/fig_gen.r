library("ape")
library("ggpubr")
library("reshape")
library("optparse")
library("plyr")
library("stringr")
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
  make_option(c("-d", "--dir"), type="character", default="plots_all", help="create directory for outputs",metavar="character")  
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

#Remove outliers 
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


#Prepare tables with branch lengths 
#t1 = estimated, t2 = true 
data_prep=function(t1,t2,cl_out = FALSE, methodname = "Method")
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
    t1_t2 = data.frame(Method = methodname, br_id = t1$variable,
                       estimated_brl = t1$value,
                       true_brl = t2$value,
                       residual = t2$value - t1$value,
                       se = (t1$value - t2$value)^2,
                       ae = abs(t1$value - t2$value))
    n_br = length(unique(t1_t2$br_id))
    br_ids = c(paste("Branch",1:(n_br-1),sep="_"),"Tree_length")
    t1_t2$br_id = rep(br_ids,table(t1_t2$br_id))
    t1_t2$br_id=factor(t1_t2$br_id,levels=str_sort(unique(t1_t2$br_id), numeric = TRUE))
    return(t1_t2)    
}


#Preparing all dataframes
prep_all=function(cnn,cnn_roe,mlp,mlp_roe,ml,true_brl,clean = FALSE)
{
    d_cnn = data_prep(cnn,true_brl,methodname = "CNN",cl_out = clean) 
    d_cnn_roe = data_prep(cnn_roe,true_brl,methodname = "CNN-ROE",cl_out = clean)
    d_mlp= data_prep(mlp,true_brl,methodname = "MLP",cl_out = clean) 
    d_mlp_roe = data_prep(mlp_roe,true_brl,methodname = "MLP-ROE",cl_out = clean)
    d_ml = data_prep(ml,true_brl,methodname = "ML",cl_out = clean)
    d_master = rbind(d_cnn,d_cnn_roe,d_mlp,d_mlp_roe,d_ml)
    #Plotting order consistent CNN,CNN-ROE,MLP,MLP-ROE,ML
    d_master$Method=factor(d_master$Method,levels=c("CNN","CNN-ROE","MLP","MLP-ROE","ML"))
    return(d_master)
}    

########## Plotting/Metadata storing functions ##########

# Over-under estimation violin plot 
get_violin=function(t1,filename = "violin.pdf")
{
    my_q = quantile(t1$residual,probs = seq(0, 1, 0.001))
    n_col = length(unique(t1$br_id))
    facet_labs = gsub("_"," ",unique(t1$br_id))
    p = ggviolin(t1, 
                 x = "Method",
                 y = "residual",
                 xlab=FALSE,
                 facet.by="br_id",
                 fill="Method",
                 palette =c("#999999", "#E69F00", "#56B4E9","#009e73","#f0e442"),
                 ncol=n_col,
                 panel.labs = list(br_id = facet_labs),
                 size=0.1,
                 draw_quantiles = 0.5)+
                 scale_y_continuous(limits = c(my_q["0.1%"], my_q["99.9%"]))+
                 theme(axis.text.x=element_blank())+
                 geom_hline(yintercept=0, linetype="solid", color = "red", size=0.3)+labs(fill = "Method:")
    
    fracs = tapply(t1$residual<0,list(t1$Method,t1$br_id),sum)/tapply(t1$residual>=0,list(t1$Method,t1$br_id),sum)
    n_comparisons = (tapply(t1$residual>=0,list(t1$Method,t1$br_id),sum)+tapply(t1$residual<0,list(t1$Method,t1$br_id),sum))[1,1]
    fracs = melt(fracs)
    names(fracs)=c("Method","br_id","frac")
    n_over = tapply(t1$residual<0,list(t1$Method,t1$br_id),sum)
    fracs$nover = melt(n_over)$value  
    fracs$pval = round(sapply(fracs$nover,function(f) binom.test(f,n=n_comparisons)$p.val),3)
    fracs$my_col = ifelse(fracs$frac<1 & fracs$pval < 0.05,"blue",ifelse(fracs$frac>1 & fracs$pval < 0.05,"red","black"))                       
    #target = c("CNN","CNN-ROE","MLP","MLP-ROE","ML")
    #coordinates for 5 methods per facet
    #Order of labels 
    fracs$x = c(1,2,3,4,5)
    fracs$y = my_q["99.9%"]
    p1 = p+geom_text(data = fracs,mapping = aes(x = x, y = y, label = round(frac,2)),color = fracs$my_col,size=2.5)
    ggsave(plot = p1, width = n_col+4, height = 5, dpi = 300, filename = filename) 
    return(p1)                          
}    

                            
# Corrleation plot across all branches and methods
get_cor=function(t1,filename = "corr.pdf")
{
    n_col = length(unique(t1$br_id))
    facet_labs = gsub("_"," ",unique(t1$br_id))
    p = ggscatter(t1, x = "estimated_brl", y = "true_brl",
                  facet.by=c("Method","br_id"),
                  ncol=n_col,
                  scales="free",
                  panel.labs = list(br_id = facet_labs),
                  size=0.1)+
                  geom_abline(slope = 1,color = "red",size=0.3)+
                  xlab("Estimated branch lengths")+
                  ylab("True branch lengths")+
                  theme(axis.text=element_text(size=6))
    stats = melt(tapply(t1$se,list(t1$Method,t1$br_id),mean))
    stats$mae = melt(tapply(t1$ae,list(t1$Method,t1$br_id),mean))$value 
    names(stats)=c("Method","br_id","MSE","MAE")
    stats$er=paste(paste("MSE:",round(stats$MSE,5),sep=""), paste("\nMAE:",round(stats$MAE,5),sep=""))
    p1 = p + geom_text(data = stats,mapping = aes(x = -Inf, y = Inf, label = er),hjust = -0.1, vjust = 1.1,color = "purple",size=2.5)
    ggsave(plot = p1, width = n_col, height = 6, dpi = 300, filename = filename) 
    return(p1)
}                            
                            


gghistogram_util = function(t1_m,my_facets=c("Method","br_id"),p_labs = list(br_id = ""),my_scales="free")
{
    p = gghistogram(t1_m, x = "value",y="..density..",
                   fill = "variable",
                   bins =50,
                   facet.by=my_facets,
                   palette = c("red", "navy"),
                   scales = my_scales,
                   add="median",
                   ncol=6,
                   color=NA,
                   add_density = TRUE,
                   panel.labs = p_labs,
                   add.params=list(linetype = "dashed",size=0.3))+
                   xlab("Branch length")+
                   theme(axis.text=element_text(size=6))
    return(p)
}    

                              
                              
                              
                              
                              
# Predicted vs. True branch length distribution plot                               
get_hist = function(t1,filename = "hist_bls")
{                            
   n_col = length(unique(t1$br_id))
   facet_labs = gsub("_"," ",unique(t1$br_id)) 
   t1_m = melt(t1[,c("Method","br_id","estimated_brl","true_brl")])
   t1_m$variable = ifelse(t1_m$variable=="true_brl","True","Estimated") 
   p = gghistogram_util(t1_m,p_labs = list(br_id = facet_labs))
   p_cl = ggpar(p, legend.title = "Branch length:")
   ggsave(plot = p_cl, width = n_col, height = 6, dpi = 300, filename = paste(filename,"hist_all.pdf",sep=""))
   p1 = gghistogram_util(t1_m[t1_m$br_id=="Tree_length",],"Method",p_labs = list(Method=c("CNN","CNN-ROE","MLP","MLP-ROE","ML")),my_scales="fixed")
   p1_cl = ggpar(p1, legend.title = "Branch length:")
   ggsave(plot = p1_cl, width = n_col, height = 3, dpi = 300, filename = paste(filename,"hist_treeL.pdf",sep="")) 
   return(p1_cl)   
}


#Get barplots and tables with stats 
get_stats=function(t1, filename = "bars_bls")
{
    my_stats= ddply(t1,.(Method,br_id),summarize,
                    MSE = round(mean(se),6),
                    MAE =  round(mean(ae),6),
                    Rho =  round(cor.test(estimated_brl,true_brl)$estimate,6),
                    P_rho =  round(cor.test(estimated_brl,true_brl)$p.value,6),
                    Bias =  round(sqrt(mean((quantile(true_brl,probs = seq(0, 1, 0.01))-quantile(estimated_brl,probs = seq(0, 1, 0.01)))^2)),6),
                    D =  round(ks.test(estimated_brl,true_brl)$statistic,6),
                    P_D =  round(ks.test(estimated_brl,true_brl)$p.value,6))

    n_col = length(unique(t1$br_id))
    facet_labs = gsub("_"," ",unique(t1$br_id))
    p1 = ggbarplot(my_stats,"Method","MSE",fill = "Method",color = "white",
                   facet.by="br_id",
                   scales="free",
                   label = F,
                   ylab="MSE",
                   xlab=F,
                   ncol=n_col,
                   panel.labs = list(br_id = facet_labs))+  
                   scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9","#009e73","#f0e442"))+
                   theme(axis.text.x=element_blank())+ 
                   labs(fill = "Method:")
    ggsave(plot = p1, width = n_col+4, height = 5, dpi = 300, filename = paste(filename,"_stats_mse.pdf",sep=""))
    
    my_stats_m = melt(my_stats)
    p2 = ggbarplot(my_stats_m,"Method","value",fill = "Method",color = "white",
                   facet.by=c("variable","br_id"),
                   scales="free",
                   label = TRUE,
                   xlab=F,
                   ylab=F,
                   lab.vjust=0.3,
                   lab.size = 2,
                   ncol=n_col,
                   panel.labs = list(br_id = facet_labs))+
                   theme(axis.text.x=element_blank(),axis.text=element_text(size=7))+
                   scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9","#009e73","#f0e442"))+
                   labs(fill = "Method:")  
    
   ggsave(plot = p2, width = 20, height = 20, dpi = 300, filename = paste(filename,"_stats_all.pdf",sep=""))  
    
   #Save table with stats
   write.csv(cast(my_stats_m,Method~variable~br_id),paste(filename,"_stats_all.csv",sep=""),row.names=T)
   #Save stats into separate tables  
   for (i in (unique(my_stats_m$variable)))
   { 
       write.csv(cast(my_stats_m[my_stats_m$variable==i,],Method~variable~br_id),paste(filename,"_stats_",i,".csv",sep=""),row.names=T)
   }
       return(p1)
}    


get_all_meta=function(t1,dir_name)
{
    p1 = get_violin(t1,filename = paste(dir_name,"/","violin_bls.pdf",sep=""))                              
    p2 = get_cor(t1,filename =  paste(dir_name,"/","corr_bls.pdf",sep=""))
    p3 = get_hist(t1,filename =  paste(dir_name,"/","hist_bls.pdf",sep=""))
    p4 = get_stats(t1, filename =  paste(dir_name,"/","bars_bls.pdf",sep=""))
}
                                                            
                              
                              
                              
                              
dir.create(opt$dir)
dir.create(paste(opt$dir,"/raw",sep=""))
dir.create(paste(opt$dir,"/clean",sep=""))                          

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
#true_brl = read.table("test_simulation.Y.txt");ml = bls_extract("ML_all.tre"); cnn_roe = read.table("brls.predicted.cnn.reg.sqrt.txt"); mlp_roe = read.table("brls.predicted.mlp.reg.sqrt.txt"); mlp = read.table("brls.predicted.mlp.sqrt.txt"); cnn = read.table("brls.predicted.cnn.sqrt.txt") 

d_master_raw = prep_all(cnn,cnn_roe,mlp,mlp_roe,ml,true_brl,clean = FALSE) 
d_master_clean = prep_all(cnn,cnn_roe,mlp,mlp_roe,ml,true_brl,clean = TRUE)                              
                              
get_all_meta(d_master_raw, dir_name = paste(opt$dir,"/raw",sep=""))
get_all_meta(d_master_clean, dir_name = paste(opt$dir,"/clean",sep=""))                              
                              
