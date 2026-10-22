# Analyze results of simulations
# Focus on end-point results

library(resample) # To use colVars
library(ggplot2)

setwd("/Users/rexjiang/Desktop/novelty/ver_4_convergence/")

effect_all<-read.table("data_effect_all.txt",sep="\t");effect_all=as.matrix(effect_all)
n_regulator=ncol(effect_all);n_dim=nrow(effect_all)

X_ances=rep(1,n_regulator)
z_ances=effect_all%*%X_ances # Calculate ancestral character state

d_selection_all=c(5,10,15,20,25,30,35,40,45,50)

X_ances_new<-read.table(paste("./out_pdx/out_",sim,"/X_ances_diverged.txt",sep=""),sep="\t");X_ances_new=as.matrix(X_ances_new)

# Go through each set of repeated simulations
n_simulation=10
out_list=list() # List of matrices of summary statistics, each corresponding to a repeated set of simulations
for(sim in 1:n_simulation){
  out=matrix(0,nrow=length(d_selection_all),ncol=8) # Data matrix of summary statistics
  for(c in 1:length(d_selection_all)){
    d_selection=d_selection_all[c]
    fn0=paste("./out_pdx/out_",sim,"/selection_list_",d_selection,".txt",sep="")
    fn1=paste("./out_pdx/out_",sim,"/sim_out_X_",d_selection,".txt",sep="")
    fn2=paste("./out_pdx/out_",sim,"/sim_out_z_",d_selection,".txt",sep="")
    fn3=paste("./out_pdx/out_",sim,"/sim_out_fitness_",d_selection,".txt",sep="")
    sub_selection<-read.table(fn0,sep="\t");sub_selection=sub_selection[,1]
    X_out_all<-read.table(fn1,sep="\t");X_out_all=as.matrix(X_out_all)
    z_out_all<-read.table(fn2,sep="\t");z_out_all=as.matrix(z_out_all)
    fitness_out_all<-read.table(fn3,sep="\t");fitness_out_all=fitness_out_all[,1]
    
    # Variation among replicate lineages
    exp_var=colVars(log(X_out_all)) # Expression variance of each regulator
    out[c,1]=sum(exp_var) # Expression variance
    out[c,2]=out[c,1]*sqrt(2/(nrow(X_out_all)-1)) # SE of expression variance
    
    # Expression divergence from ancestor
    X_dist=rep(0,nrow(X_out_all))
    for(i in 1:nrow(X_out_all)){
      X_dist[i]=sqrt(sum((log(X_out_all[i,])-log(X_ances_new[i,]))^2)) # Euclidean distance to ancestral expression, computed for each lineage
    }
    out[c,3]=mean(X_dist) # Mean across lineages
    out[c,4]=sd(X_dist)/sqrt(nrow(X_out_all)) # SE
    out[c,5]=mean(as.vector(z_out_all)) # Mean trait value of traits under selection
    out[c,6]=sd(as.vector(z_out_all))/sqrt(nrow(z_out_all)*d_selection) # SE (calculated as SE of mean of all elements of a matrix)
    out[c,7]=mean(fitness_out_all) # Mean fitness
    out[c,8]=sd(fitness_out_all)/length(fitness_out_all) # SE of fitness
  }
  
  d=data.frame(d_selection_all,out)
  colnames(d)=c("d_selection","var","var_se","dist","dist_se","z","z_se","w","w_se")
  out_list[[sim]]=d
}

# Merge data file
d_merged=out_list[[1]]
for(sim in 2:n_simulation){
  d_merged=rbind(d_merged,out_list[[sim]])
}
write.table(d_merged,file="out_sum_merged_pdx.txt",sep="\t")

# Re-read merged data file
d<-read.table("out_sum_merged_pdx.txt",sep="\t",header=TRUE)
X_ances_new<-read.table(paste("./out_pdx/out_",sim,"/X_ances_diverged.txt",sep=""),sep="\t");X_ances_new=as.matrix(X_ances_new)
var_ref=sum(colVars(log(X_ances_new))) # Expression divergence at the end of stage 1

# Plot expression variance among lineages
g<-ggplot(d,aes(x=d_selection,y=var))
g=g+geom_point(size=3)+geom_smooth()
g=g+geom_hline(yintercept=var_ref,color="red",linetype="dashed",linewidth=1.5)
g=g+theme_classic()+theme(axis.text=element_text(size=12))+xlab("")+ylab("")
g=g+scale_x_continuous(limits=c(4,51),breaks=c(5,10,15,20,25,30,35,40,45,50))
g=g+scale_y_continuous(limits=c(0,60),breaks=c(0,10,20,30,40,50))

# Plot divergence from ancestral expression
h<-ggplot(d,aes(x=d_selection,y=dist))
h=h+geom_point(size=3)+geom_smooth()
h=h+theme_classic()+theme(axis.text=element_text(size=12))+xlab("")+ylab("")
h=h+scale_x_continuous(limits=c(4,51),breaks=c(5,10,15,20,25,30,35,40,45,50))
h=h+scale_y_continuous(limits=c(0,8.5),breaks=c(0,2,4,6,8))

# Plot fitness
f<-ggplot(d,aes(x=d_selection,y=w))
f=f+geom_point(size=3)+geom_smooth()
f=f+theme_classic()+theme(axis.text=element_text(size=12))+xlab("")+ylab("")
f=f+scale_x_continuous(limits=c(4,51),breaks=c(5,10,15,20,25,30,35,40,45,50))+scale_y_continuous(limits=c(0,1.05),breaks=c(0,0.2,0.4,0.6,0.8,1))

ggsave("./plot/plot_pdx_exp_var.pdf",plot=g,width=7,height=5)
ggsave("./plot/plot_pdx_exp_div.pdf",plot=h,width=7,height=5)
ggsave("./plot/plot_pdx_fitness.pdf",plot=f,width=7,height=5)
