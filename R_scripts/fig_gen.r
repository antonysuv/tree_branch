library("ape")
library("ggpubr")
library("reshape")
library("optparse")
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
data_prep=function(t1,t2,cl_out = FALSE)
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


# t1 = table after data-prep formatting 
cor_test = function(t1,methodname)
{
    my_stats=c()
    for (i in names(table(t1$br_id)))
    {
        sub_t1 = t1[t1$br_id==i,c("estimated_brl","true_brl")]
        d = sub_t1[,1]-sub_t1[,2]
        my_cdfbias = sqrt(mean((quantile(sub_t1[,2],probs = seq(0, 1, 0.01))-quantile(sub_t1[,1],probs = seq(0, 1, 0.01)))^2))
        my_ks = ks.test(sub_t1[,1],sub_t1[,2])
        my_mse = mean(d^2)
        my_mae = mean(abs(d))
        my_test = cor.test(sub_t1[,1],(sub_t1[,2]), method = "spearman",exact=FALSE)
        my_values = c(round(my_mse,6),
                          round(my_mae,6),
                          round(my_test$estimate,3),
                          round(my_test$p.value,6),
                          round(my_cdfbias,4),
                          my_ks$statistic,
                          round(my_ks$p.value,6)
                          )
        my_stats = c(my_stats,my_values)
    }
    return(c(methodname,my_stats))
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

# t1 = dd table  
get_bars=function(t1,filename = "bar.pdf")
{
    facet_labs = gsub("_"," ",unique(dd_t$branch))
    p1 = ggbarplot(t1,"variable","value",fill = "variable",color = "white",facet.by=c("test_stats_f","branch"),
             scales="free",
             label = TRUE,
             lab.vjust=0.3,
             lab.size = 2,
             sort.val = "none",
             ylab=F,
             panel.labs = list(branch = facet_labs),
             xlab=F)+theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9","#009e73","#f0e442"))+theme(axis.text.x=element_blank())+labs(fill = "Method")
   
    p2 = ggbarplot(t1[t1$test_stats=="MSE",],"variable","value",fill = "variable",color = "white",
                   facet.by=c("branch"),
                   scales="free",
                   label = F,
                   ylab=F,
                   xlab=F,
                   ncol=6,
                   panel.labs = list(branch = facet_labs))+ theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))+ scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9","#009e73","#f0e442"))+theme(axis.text.x=element_blank())+ labs(fill = "Method")
   
   ggsave(plot = p1, width = 20, height = 20, dpi = 300, filename = paste("all_",filename,sep="")) 
   ggsave(plot = p2, width = 10, height = 5, dpi = 300, filename = paste("mse_",filename,sep=""))
    
    

}    

ggbarplot(dd_t,"variable","value",fill = "variable",color = "white",facet.by=c("test_stats_f","branch"),scales="free",label = TRUE,lab.vjust=0.3,lab.size = 2,sort.val = "none",ylab=F,xlab=F)+ theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust=1))





dir.create(opt$dir)

cat("Read true\n")
tt = read.table(opt$truebl)  
cat("Read ml\n")
ml = bls_extract(opt$ml) 
cat("Read cnn ROE\n")
cnn = read.table(opt$cnn)
cat("Read mlp ROE\n")
mlp = read.table(opt$mlp)
cat("Read mlp\n")
mlp_noreg = read.table(opt$mlp_noreg)
cat("Read cnn\n")
cnn_noreg = read.table(opt$cnn_noreg)

#Write table with different stats
cat("Prepare data\n")
d_cnn_noreg = data_prep(cnn_noreg,tt); d_cnn = data_prep(cnn,tt); d_mlp_noreg = data_prep(mlp_noreg,tt); d_mlp = data_prep(mlp,tt); d_ml = data_prep(ml,tt)
dd = data.frame(rbind(cor_test(d_cnn_noreg,"CNN"),
                    cor_test(d_cnn,"CNN-ROE"),
                    cor_test(d_mlp_noreg,"MLP"),
                    cor_test(d_mlp,"MLP-ROE"),
                    cor_test(d_ml,"ML")
                   ))
#write.table(dd,"test_tab.txt",quote=F,row.names=F)
dd_t = data.frame(t(dd)) 
names(dd_t) = dd_t[1,]
dd_t=dd_t[-1,]
test_stats = rep(c("MSE","MAE","Rho","P_Rho","Bias","D","P_D"),(ncol(dd)-1)/7)
dd_t$test_stats = test_stats
branch_names = rep(c(paste("Branch",1:((ncol(dd)-2)/7),sep="_"),"Tree_length"),each=7)
dd_t$branch=branch_names
dd_t=melt(dd_t, measure.vars=c("CNN","CNN-ROE","MLP","MLP-ROE","ML"))
dd_t$value=as.numeric(dd_t$value)
dd_t$test_stats_f=factor(dd_t$test_stats,levels=c("MSE","MAE","Rho","P_Rho","Bias","D","P_D"))

ggbarplot(dd_t[dd_t$test_stats=="Bias",],"variable","value",fill = "variable",facet.by="branch",scales = "free",nrow = 1,label = TRUE,short.panel.labs=TRUE,lab.size = 3,xlab = FALSE,label.pos = "in",lab.nb.digits=4)





test_stats = c(" ",rep(c("MSE","MAE","Rho","P_Rho","Bias","D","P_D"),(ncol(dd)-1)/7))
test_stats = data.frame(t(data.frame(test_stats)))
branch_names = c("Method",rep(c(paste("Branch",1:((ncol(dd)-2)/7),sep="_"),"Tree_length"),each=7))
names(test_stats) = names(branch_names)
names(dd) = names(branch_names)
print(head(test_stats))
print(head(branch_names))
print(head(dd))
dd = rbind(branch_names,test_stats,dd) 
#names(dd) = c("Method",rep(c(paste("Branch",1:((ncol(dd)-2)/7),sep="_"),"Tree_length"),each=7))
write.csv(dd, paste(opt$dir,"/",opt$dir,".csv",sep=""),row.names=F,col.names=F)

get_cor(ml,tt,filename = paste(opt$dir,"/","plot_ml_dot_raw.pdf",sep=""))
get_cor(cnn,tt,filename =  paste(opt$dir,"/","plot_cnn_dot_raw.pdf",sep=""))
get_cor(mlp,tt,filename =  paste(opt$dir,"/","plot_mlp_dot_raw.pdf",sep=""))
get_cor(cnn_noreg,tt,filename =  paste(opt$dir,"/","plot_cnn_noreg_dot_raw.pdf",sep=""))
get_cor(mlp_noreg,tt,filename =  paste(opt$dir,"/","plot_mlp_noreg_dot_raw.pdf",sep=""))

get_hist(ml,tt,filename =  paste(opt$dir,"/","plot_ml_hist_raw.pdf",sep=""))
get_hist(cnn,tt,filename =  paste(opt$dir,"/","plot_cnn_hist_raw.pdf",sep=""))
get_hist(mlp,tt,filename =  paste(opt$dir,"/","plot_mlp_hist_raw.pdf",sep=""))
get_hist(cnn_noreg,tt,filename =  paste(opt$dir,"/","plot_cnn_hist_noreg_raw.pdf",sep=""))
get_hist(mlp_noreg,tt,filename =  paste(opt$dir,"/","plot_mlp_hist_noreg_raw.pdf",sep=""))

get_cor(ml,tt,cl_out = TRUE,filename =  paste(opt$dir,"/","plot_ml_dot_clean.pdf",sep=""))
get_cor(cnn,tt,cl_out = TRUE,filename =  paste(opt$dir,"/","plot_cnn_dot_clean.pdf",sep=""))
get_cor(mlp,tt,cl_out = TRUE,filename =  paste(opt$dir,"/","plot_mlp_dot_clean.pdf",sep=""))
get_cor(cnn_noreg,tt,cl_out = TRUE,filename =  paste(opt$dir,"/","plot_cnn_noreg_dot_clean.pdf",sep=""))
get_cor(mlp_noreg,tt,cl_out = TRUE,filename =  paste(opt$dir,"/","plot_mlp_noreg_dot_clean.pdf",sep=""))

get_hist(ml,tt,cl_out = TRUE,filename =  paste(opt$dir,"/","plot_ml_hist_clean.pdf",sep=""))
get_hist(cnn,tt,cl_out = TRUE,filename =  paste(opt$dir,"/","plot_cnn_hist_clean.pdf",sep=""))
get_hist(mlp,tt,cl_out = TRUE,filename =  paste(opt$dir,"/","plot_mlp_hist_clean.pdf",sep=""))
get_hist(cnn_noreg,tt,cl_out = TRUE,filename =  paste(opt$dir,"/","plot_cnn_hist_noreg_clean.pdf",sep=""))
get_hist(mlp_noreg,tt,cl_out = TRUE,filename =  paste(opt$dir,"/","plot_mlp_hist_noreg_clean.pdf",sep=""))