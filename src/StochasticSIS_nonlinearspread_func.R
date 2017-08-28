# Stochastic SIS model
#Nonlinear spread (bets*S*I/(1+alpha*I^2)), so that spread decreases with increasing infectives (e.g. due to more resources spent on control),
#added to get interesting dynamics (outbreaks, no stable eq) -> See Liu, Levin, Iwasa 1986 for complete discussion, more general forms
#also has an "observation" component, determined by parameter rho, i.e. our probability of detecting infected individuals
# Thus can change disease dynamics (spread speed bet, recovery rate gamma, strength of intervention alpha) or observation intensity (rho)
#JLW 23/8/17

symbolizeOutput <- function(x,threshold,extremes=0,burnin=100){
  x <- x[burnin:length(x)] #Allow for brief burn-in
  
  if(extremes==0){
    x <- (x>=quantile(x,threshold))*1
  } else {
    if(threshold<0.5){stop("If using extremes=1, please use threshold >=0.5")}
    x <- (x>=quantile(x,threshold))+(x<=quantile(x,1-threshold))
  }
  x <- paste(x,collapse="")
  return(x)
}



stepSIS <- function(S,I,parameters){
  N <- parameters[["N"]]
  bet <- parameters[["beta"]]
  gam <- parameters[["gamma"]]
  rho <- parameters[["rho"]]
  alpha <- parameters[["alpha"]]
  
  #t1 <- rbinom(n=1,size=S,prob=1-exp(-bet*I/N))
  t1 <- rbinom(n=1,size=S,prob=1-exp(-bet*I/(1+alpha*I^2)))
  t2 <- rbinom(n=1,size=I,prob=1-exp(-gam))
  
  S  <- S - t1 + t2;
  I  <- I - t2 + t1;
  
  B <- rpois(n=1,lambda=rho*I+1e-6);
  
  return(c(S,I,B))
}


runSIS <- function(parameters,nsteps,full_out=F){
  N <- parameters[["N"]]
  I <- 1
  S <- N-1
  
  sim_list <- list()
  counter <- 1
  for(i in 1:nsteps){
    sim_list[[counter]] <- stepSIS(S,I,parameters)
    S <- sim_list[[counter]][1]
    I <- max(sim_list[[counter]][2],1)
    counter=counter+1
  }
  
  run_df <- as.data.frame(do.call(rbind.data.frame, sim_list),stringsAsFactors=FALSE)
  names(run_df) <- c("S","I","B")
  if(full_out==T){
    return(run_df)
  } else {
    return(run_df$B)
  }
    
}

varyParameter_SIS <- function(parameter_name,parameter_seq,parameters_list,sim_length=10000,outbreak_percentile=0.75){
  out_list <- list()
  counter <- 1
  for(par in parameter_seq){
    parameters_list[[parameter_name]] <- par
    run_raw <- runSIS(parameters_list,sim_length) 
    run_sym <- symbolizeOutput(run_raw,outbreak_percentile,extremes = 1)
    gEE_out <- getExcessEntropy(seq=run_sym,maxL=20)
    E <- gEE_out[1]
    hmu_est <- gEE_out[2]
    out_list[[counter]] <- c(par,hmu_est,E)
    counter <- counter+1
  }
  out_df <- as.data.frame(do.call(rbind.data.frame, out_list),stringsAsFactors=FALSE)
  names(out_df) <- c("par_value","hmu","E")
  return(out_df)
}










####SIRS model #####################################################################################

stepSIRS <- function(S,I,R,parameters){
  N <- parameters[["N"]]
  bet <- parameters[["beta"]]
  gam <- parameters[["gamma"]]
  rho <- parameters[["rho"]]
  alpha <- parameters[["alpha"]]
  mu <- parameters[["mu"]]
  
  #t1 <- rbinom(n=1,size=S,prob=1-exp(-bet*I/N))
  t1 <- rbinom(n=1,size=S,prob=1-exp(-bet*I/(1+alpha*I^2)))
  t2 <- rbinom(n=1,size=I,prob=1-exp(-gam))
  t3 <- rbinom(n=1,size=R,prob=1-exp(-mu))
  
  S  <- S - t1 + t3;
  I  <- I - t2 + t1;
  R  <- R - t3 + t2;
  
  B <- rpois(n=1,lambda=rho*I+1e-6);
  
  return(c(S,I,R,B))
}


runSIRS <- function(parameters,nsteps,full_out=F){
  N <- parameters[["N"]]
  I <- 1
  S <- N-1
  R <- 0
  
  sim_list <- list()
  counter <- 1
  for(i in 1:nsteps){
    sim_list[[counter]] <- stepSIRS(S,I,R,parameters)
    S <- sim_list[[counter]][1]
    I <- max(sim_list[[counter]][2],1)
    R <- sim_list[[counter]][3]
    counter <- counter+1
  }
  
  run_df <- as.data.frame(do.call(rbind.data.frame, sim_list),stringsAsFactors=FALSE)
  names(run_df) <- c("S","I","R","B")
  if(full_out==T){
    return(run_df)
  } else {
    return(run_df$B)
  }
  
}


varyParameter_SIRS <- function(parameter_name,parameter_seq,parameters_list,sim_length=10000,outbreak_percentile=0.75){
  out_list <- list()
  counter <- 1
  for(par in parameter_seq){
    parameters_list[[parameter_name]] <- par
    run_raw <- runSIRS(parameters_list,sim_length) 
    run_sym <- symbolizeOutput(run_raw,outbreak_percentile,extremes = 1)
    gEE_out <- getExcessEntropy(seq=run_sym,maxL=20)
    E <- gEE_out[1]
    hmu_est <- gEE_out[2]
    out_list[[counter]] <- c(par,hmu_est,E)
    counter <- counter+1
  }
  out_df <- as.data.frame(do.call(rbind.data.frame, out_list),stringsAsFactors=FALSE)
  names(out_df) <- c("par_value","hmu","E")
  return(out_df)
}
