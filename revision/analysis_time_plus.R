# Analyze results of simulations
# Examine trend over time
# For each simulation setting, calculate statistics over time

library(resample) # To use colVars
library(ggplot2)

setwd("/Users/rexjiang/Desktop/novelty/ver_4_convergence/")

# Function to calculate fitness (Gaussian function of a multivariate phenotype with w_max=1)
fitness_calc<-function(z,opt,omega) {
  Lambda<-solve(omega)
  D2<-t(z-opt)%*%Lambda%*%(z-opt)
  w<-exp(-D2/2)
  return(as.numeric(w))
}

# Initialize (set parameters used for all simulations)
# Read pre-generated matrices if already generated
effect_all<-read.table("data_effect_all.txt",sep="\t");effect_all=as.matrix(effect_all)
n_regulator=ncol(effect_all);n_dim=nrow(effect_all)

# Strength of selection
z_div_scale=2
omega_scale=100

n_step=200
n_rep=20

X_ances=rep(1,n_regulator)
z_ances=effect_all%*%X_ances # Calculate ancestral character state

d_selection_all=c(5,10,15,20,25,30,35,40,45,50)

n_simulation=10
# List of matrices of summary statistics, each element corresponding to a repeated set of simulations
var_time_all=list() 
for(sim in 1:n_simulation){
  var_time=matrix(0,nrow=length(d_selection_all),ncol=n_step) # Expression variance over time
  for(c in 1:length(d_selection_all)){
    d_selection=d_selection_all[c]
    fn1=paste("./out_plus/out_",sim,"/sim_out_X_time_",d_selection,".txt",sep="")
    sub_selection<-read.table(fn0,sep="\t");sub_selection=sub_selection[,1]
    X_time<-read.table(fn1,sep="\t");X_time=as.matrix(X_time)
    
    # Optimal phenotype and shape of adaptive landscape
    z_opt=z_ances[sub_selection]+z_div_scale
    omega=matrix(0,nrow=d_selection,ncol=d_selection);diag(omega)=omega_scale # Adaptive landscape (no correlational selection for now)
    
    for(t in 1:n_step){
      X_t=matrix(0,nrow=n_rep,ncol=n_regulator)
      for(n in 1:n_rep){
        X_t[n,]=X_time[which(X_time[,1]==n),t+1] # The n-th lineage at the end of the t-th time step
      }
      var_time[c,t]=sum(colVars(log(X_t)))
    }
  }
  var_time_all[[sim]]=var_time
}

# Convert to ggplot-ready format
d_merged=matrix(0,nrow=1,ncol=5)
for(sim in 1:n_simulation){
  for(c in 1:length(d_selection_all)){
    for(t in 1:n_step){
      row=c(sim,d_selection_all[c],t,var_time_all[[sim]][c,t])
      d_merged=rbind(d_merged,row)
    }
  }
}
d_merged=d_merged[-1,];rownames(d_merged)=NULL
colnames(d_merged)=c("sim","d_selection","t","var")
write.table(d_merged,file="out_time_merged_plus.txt",sep="\t")

d<-read.table("out_time_merged_plus.txt",sep="\t")
d$d_selection=factor(d$d_selection)

g<-ggplot(d,aes(x=t,y=var,colour=d_selection))
g=g+geom_point(alpha=0.5)
g=g+xlab("Time")+ylab("")

ggsave("./plot/plot_time_plus_var.pdf",plot=g,width=7,height=5)


