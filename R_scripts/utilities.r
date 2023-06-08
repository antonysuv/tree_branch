library("ape")
library("ggpubr")
library("reshape")

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

bls_extract = function(file_in)
{
    phy = read.tree(file_in)
    brls = data.frame(t(data.frame(lapply(phy,"[[","edge.length")))) 
    row.names(brls) = NULL
    
    return(brls)
}    

#Quantile plot
delta     = 0.001 
quantiles = 20
z.df     = data.frame(branch_length = seq(from=0, to=10, by=delta))
z.df$dencity = dexp(z.df$branch_length)
z.df$qt  = cut(pexp(z.df$branch_length),breaks=quantiles,labels=F)


gg = ggplot(z.df,aes(x=branch_length,y=dencity))+
  geom_area(aes(x=branch_length,y=dencity,group=qt,fill=qt),color="black")+
  scale_fill_gradient2(midpoint=median(unique(z.df$qt)), guide="none") +
  theme_bw()
plot(gg)



get_wrapped = function(true_bl_tab,cnn_tab,cnn_roe_tab,mlp_tab,mlp_roe_tab,ml_tab)
{
    true_bl_tab = read.table(true_bl_tab)
    cnn_tab = read.table(cnn_tab)
    cnn_roe_tab = read.table(cnn_roe_tab)
    mlp_tab = read.table(mlp_tab)
    mlp_roe_tab = read.table(mlp_roe_tab)
    ml_tab=bls_extract(ml_tab)
    
    d = data.frame(true_bl = as.vector(unlist(true_bl_tab)),
                  cnn = as.vector(unlist(cnn_tab)),
                  cnn_roe = as.vector(unlist(cnn_roe_tab)),
                  mlp = as.vector(unlist(mlp_tab)),
                  mlp_roe = as.vector(unlist(mlp_roe_tab)),
                  ml = as.vector(unlist(ml_tab)))
    return(d)
}

dd = get_wrapped("test_simulation.Y.txt",
    "brls.predicted.cnn.sqrt.txt",
    "brls.predicted.cnn.reg.sqrt.txt",
    "brls.predicted.mlp.sqrt.txt",
    "brls.predicted.mlp.reg.sqrt.txt",
    "ML_all.tre")

qs = qexp(seq(0,1,0.05))
#Replace Inf by big branch length
qs[!is.finite(qs)] = 1000
qs_tab = data.frame(left = qs[-length(qs)],right = qs[-1])

main_tab_mae = c()
main_tab_mse = c()
for (i in 1:nrow(qs_tab))
{
    q_result_mae = c()
    q_result_mse = c()
    qs_int = qs_tab[i,]
    dd_sub = dd[dd$true_bl >= qs_int$left & dd$true_bl < qs_int$right,]
    for (n in c("cnn","cnn_roe","mlp","mlp_roe","ml"))
    {
        mae = mean(abs(dd_sub$true_bl - dd_sub[,n]))
        mse = mean((dd_sub$true_bl - dd_sub[,n])^2)
        q_result_mae = c(q_result_mae,mae)
        q_result_mse = c(q_result_mse,mse)
    }
    main_tab_mae = rbind(main_tab_mae,q_result_mae)
    main_tab_mse = rbind(main_tab_mse,q_result_mse)

}

write.csv(main_tab_mae, "quantile_error_mae.csv")
write.csv(main_tab_mse, "quantile_error_mse.csv")