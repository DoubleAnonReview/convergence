# Simulate the origin of a novel character identity under a sequential-fixation regime

#library(Matrix) # Used to check system solvability
#library(resample) # To use colVars
#library(ggplot2)

setwd("/Users/rexjiang/Desktop/novelty/ver_4_convergence/")

# Function to calculate fitness (Gaussian function of a multivariate phenotype with w_max=1)
fitness_calc<-function(z,opt,omega) {
  Lambda<-solve(omega)
  D2<-t(z-opt)%*%Lambda%*%(z-opt)
  w<-exp(-D2/2)
  return(as.numeric(w))
}

# Function to calculate fixation probability
# Parameters: ancestral fitness (wa), mutant fitness (wm), effective population size (Ne)
fix_prob_calc <- function(wa,wm,Ne){
  if(wa>0){
    s=(wm/wa)-1 # Coefficient of selection
    if(s==0){
      p=1/(2*Ne) # Neutral mutation
    }else{
      p=(1-exp(-2*s))/(1-exp(-4*Ne*s))
    }
  }else{ # The ancestral fitness is close to 0 (close enough to be recognized as zero by R)
    if(wm==0){ # Both ancestral and mutant fitness are close to 0
      p=1/(2*Ne) # The mutation is considered neutral
    }else{ # Ancestral fitness is close to 0 while mutant fitness isn't
      p=1 # The mutation is considered strongly beneficial
    }
  }
  return(p)
}

# Function to simulate evolution (step-by-step, with intermediate results saved)
# Each step should be interpreted as a relatively long time interval
# Arguments: ancestral expression levels of regulators (X_ances, a vector), per-unit-expression phenotypic effect of regulators (effect_all, a matrix), target state (z_opt, a vector), adaptive landscape (omega, a covariance matrix), 
# Arguments (continued): effective population size (Ne), simulation time (n_step), mutation rate per genome per time step (U), mutation effect size SD (sigma)
evo_sim <- function(X_ances,effect_all,z_opt,omega,Ne,n_step,U,sigma){
  X_out=matrix(0,nrow=length(X_ances),ncol=n_step) # Data matrix to store the output
  X=X_ances # 'Initialize' X
  for(t in 1:n_step){
    n_mut=rpois(1,lambda=2*Ne*U) # Number of mutations to occur over the time step
    if(n_mut>0){
      z=effect_all%*%log(X)
      w=fitness_calc(z,z_opt,omega)
      for(n in 1:n_mut){
        X_mut=X
        affect=sample(1:length(X),1) # Each mutation affects a single gene's expression
        effect=rnorm(1,mean=0,sd=sigma) # Effect of mutation (log scale)
        X_mut[affect]=X_mut[affect]*exp(effect) # Mutant expression level
        z_mut=effect_all%*%log(X_mut) # Mutant phenotype
        w_mut=fitness_calc(z_mut,z_opt,omega) # Mutant fitness
        fp=fix_prob_calc(w,w_mut,Ne) # Fixation probability
        fix=rbinom(n=1,size=1,prob=fp) # Whether the mutations fixes
        if(fix==1){
          X=X_mut;z=z_mut;w=w_mut # Pass the mutant's genotype & phenotype down to next step if the mutation is fixed
        }
      }
    }
    X_out[,t]=X # Write the expression vector at the end of the time step
  }
  return(X_out)
}

# Initialize (set parameters used for all simulations)
# Read pre-generated matrices if already generated
effect_all<-read.table("data_effect_all.txt",sep="\t");effect_all=as.matrix(effect_all)
n_regulator=ncol(effect_all);n_dim=nrow(effect_all)

# Set ancestral state
X_ances=rep(1,n_regulator) # An ancestral state where everything is lowly expressed
z_ances=effect_all%*%X_ances # Calculate ancestral character state

# Strength of selection
z_div_scale=2
omega_scale=100

# Number of replicate lineages for each setting 
n_rep=20

# Set parameters for simulation function (reduce steps and increase per-step mutation rate for new simulation algorithm)
Ne=1000;n_step=200;U=n_regulator*0.05;sigma=0.1

# Go through parameter combinations and simulate
d_selection_all=c(5,10,15,20,25,30,35,40,45,50) # Complexity levels to consider
n_simulation=10 # Number of regimes of selection to simulate for each d_selection
for(sim in 1:n_simulation){
  for(d_selection in d_selection_all){
    # List of traits subject to selection
    fn0=paste("./out/out_",sim,"/selection_list_",d_selection,".txt",sep="")
    sub_selection<-read.table(fn0,sep="\t")
    sub_selection=sub_selection[,1]
    #sub_selection
    effect_sub=effect_all[sub_selection,];z_ances_sub=z_ances[sub_selection] # Extract rows corresponding to traits under selection
    
    # Optimal phenotype and shape of adaptive landscape
    z_opt=z_ances[sub_selection]+z_div_scale
    omega=matrix(0,nrow=d_selection,ncol=d_selection);diag(omega)=omega_scale # Adaptive landscape (no correlational selection for now)
    
    X_time_list=list() # Data matrices for expression over time
    X_out_all=matrix(nrow=n_rep,ncol=n_regulator) # Data matrix for end-point expression
    z_out_all=matrix(nrow=n_rep,ncol=d_selection) # Data matrix for end-point phenotype
    fitness_out_all=rep(0,n_rep) # Data matrix for end-point fitness
    for(n in 1:n_rep){
      X_time=evo_sim(X_ances,effect_sub,z_opt,omega,Ne,n_step,U,sigma)
      X_time_list[[n]]=cbind(rep(n,nrow(X_time)),X_time)
      
      X_out=X_time[,n_step]
      z_out=effect_sub%*%log(X_out)
      fitness_out=fitness_calc(z_out,z_opt,omega)
      X_out_all[n,]=X_out
      z_out_all[n,]=z_out
      fitness_out_all[n]=fitness_out
      cat(paste(sim,"\t",d_selection,"\t",n,"\t",fitness_out,"\n",sep=""))
    }
    
    X_time_matrix=X_time_list[[1]]
    for(n in 2:n_rep){
      X_time_matrix=rbind(X_time_matrix,X_time_list[[n]])
    }
    
    fn1=paste("./out/out_",sim,"/sim_out_X_",d_selection,".txt",sep="")
    fn2=paste("./out/out_",sim,"/sim_out_z_",d_selection,".txt",sep="")
    fn3=paste("./out/out_",sim,"/sim_out_fitness_",d_selection,".txt",sep="")
    fn4=paste("./out/out_",sim,"/sim_out_X_time_",d_selection,".txt",sep="")
    write.table(X_out_all,file=fn1,sep="\t")
    write.table(z_out_all,file=fn2,sep="\t")
    write.table(fitness_out_all,file=fn3,sep="\t")
    write.table(X_time_matrix,file=fn4,sep="\t")
  }
}