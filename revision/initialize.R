# Generate data files to be used across simulations

setwd("/Users/rexjiang/Desktop/novelty/ver_4_convergence/")
#library(mvtnorm)

# C-matrix
n_regulator=50 # The number of regulators available
n_dim=1000 # The total number of phenotypic dimensions (not all under selection)
effect_all=matrix(0,nrow=n_dim,ncol=n_regulator) # Data matrix containing all regulators' phenotypic effects per unit expression
# Effects of regulators on a random trait
effect_size_base_1=rexp(n_regulator/2,rate=2)
effect_size_base_2=-effect_size_base_1
effect_size_base=c(effect_size_base_1,effect_size_base_2)
for(n in 1:n_dim){
  effect_all[n,]=signif(sample(effect_size_base,n_regulator),5)
}
write.table(effect_all,file="data_effect_all.txt",sep="\t") # Save data for future use

# Regime of selection
d_selection_all=(1:10)*5 # Complexity levels to consider
n_simulation=5 # Number of regimes of selection to simulate for each d_selection
for(sim in 1:n_simulation){
  for(d_selection in d_selection_all){
    sub_selection=sort(sample(1:n_dim,d_selection)) # Randomly draw a set of traits subject to selection
    fn=paste("./out/out_",sim,"/selection_list_",d_selection,".txt",sep="")
    write.table(sub_selection,file=fn,sep="\t")
  }
}

# Random divergence of regulator expression
n_rep=20 # Number of lineages
sigma_x=0.1 # Evolutionary rate (assumed to be equal for all)
T_0=100 # Time
X_ances_new=matrix(0,nrow=n_rep,ncol=n_regulator)
for(i in 1:n_rep){
  delta=rnorm(n_regulator,mean=0,sd=sigma_x*sqrt(T_0))
  X_ances_new[i,]=X_ances*exp(delta)
}
fn=paste("./out_pdx/X_ances_diverged.txt",sep="")
write.table(X_ances_new,file=fn,sep="\t")

# Random divergence of C-matrix
effect_all<-read.table("data_effect_all.txt",sep="\t");effect_all=as.matrix(effect_all)
n_regulator=ncol(effect_all);n_dim=nrow(effect_all)
n_rep=20
T_0=100
# Little divergence
sigma_cis=0.001
for(n in 1:n_rep){
  effect_all_new=effect_all
  for(r in 1:n_regulator){
    delta=rnorm(n_dim,mean=0,sd=sigma_cis*sqrt(T_0))
    effect_all_new[,r]=effect_all_new[,r]+delta
  }
  fn=paste("./out_pdc_slow/effect_all_diverged_",n,".txt",sep="")
  write.table(effect_all_new,file=fn,sep="\t")
}
# Moderate divergence
sigma_cis=0.01
for(n in 1:n_rep){
  effect_all_new=effect_all
  for(r in 1:n_regulator){
    delta=rnorm(n_dim,mean=0,sd=sigma_cis*sqrt(T_0))
    effect_all_new[,r]=effect_all_new[,r]+delta
  }
  fn=paste("./out_pdc/effect_all_diverged_",n,".txt",sep="")
  write.table(effect_all_new,file=fn,sep="\t")
}
# Fast divergence
sigma_cis=0.1
for(n in 1:n_rep){
  effect_all_new=effect_all
  for(r in 1:n_regulator){
    delta=rnorm(n_dim,mean=0,sd=sigma_cis*sqrt(T_0))
    effect_all_new[,r]=effect_all_new[,r]+delta
  }
  fn=paste("./out_pdc_fast/effect_all_diverged_",n,".txt",sep="")
  write.table(effect_all_new,file=fn,sep="\t")
}

