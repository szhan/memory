# Stochastic SIS model
#Nonlinear spread (bets*S*I/(1+alpha*I^2)), so that spread decreases with increasing infectives (e.g. due to more resources spent on control),
#added to get interesting dynamics (outbreaks, no stable eq) -> See Liu, Levin, Iwasa 1986 for complete discussion, more general forms
#also has an "observation" component, determined by parameter rho, i.e. our probability of detecting infected individuals
# Thus can change disease dynamics (spread speed beta, recovery rate gamma, strength of intervention alpha) or observation intensity (rho)
#JLW 23/8/17

symbolizeOutput <- function(x,threshold){
  x <- x[100:length(x)] #Allow for brief burn-in
  x <- (x>=threshold)*1
  x <- paste(x,collapse="")
  return(x)
}

stepSIS <- function(S,I,parameters){
  N <- parameters[[1]]
  beta <- parameters[[2]]
  gam <- parameters[[3]]
  rho <- parameters[[4]]
  alpha <- parameters[[5]]
  
  #t1 <- rbinom(n=1,size=S,prob=1-exp(-beta*I/N))
  t1 <- rbinom(n=1,size=S,prob=1-exp(-beta*I/(1+alpha*I^2)))
  t2 <- rbinom(n=1,size=I,prob=1-exp(-gam))
  
  S  <- S - t1 + t2;
  I  <- I - t2 + t1;
  
  B <- rpois(n=1,lambda=rho*I+1e-6);
  
  return(c(S,I,B))
}


runSIS <- function(parameters,nsteps,full_out=F){
  N <- parameters[[1]]
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


varyRho <- function(len,pars,rho_seq){
  N <- pars[["N"]]
  beta <- pars[["beta"]]
  gam <- pars[["gamma"]]
  alpha <- pars[["alpha"]]
  out_list <- list()
  counter <- 1
  for(rho in rho_seq){
    parameters <- list(N,beta,gam,rho,alpha)
    run_raw <- runSIS(parameters,10000) 
    run_sym <- symbolizeOutput(run_raw,median(run_raw))
    gEE_out <- getExcessEntropy(seq=run_sym,maxL=20)
    E <- gEE_out[1]
    hmu_est <- gEE_out[2]
    out_list[[counter]] <- c(rho,hmu_est,E)
    counter <- counter+1
  }
  out_df <- as.data.frame(do.call(rbind.data.frame, out_list),stringsAsFactors=FALSE)
  names(out_df) <- c("rho","hmu","E")
  return(out_df)
}


varyBeta <- function(len,pars,beta_seq){
  N <- pars[["N"]]
  rho <- pars[["rho"]]
  gam <- pars[["gamma"]]
  alpha <- pars[["alpha"]]
  out_list <- list()
  counter <- 1
  for(beta in beta_seq){
    parameters <- list(N,beta,gam,rho,alpha)
    run_raw <- runSIS(parameters,10000) 
    run_sym <- symbolizeOutput(run_raw,median(run_raw))
    gEE_out <- getExcessEntropy(seq=run_sym,maxL=20)
    E <- gEE_out[1]
    hmu_est <- gEE_out[2]
    out_list[[counter]] <- c(beta,hmu_est,E)
    counter <- counter+1
  }
  out_df <- as.data.frame(do.call(rbind.data.frame, out_list),stringsAsFactors=FALSE)
  names(out_df) <- c("beta","hmu","E")
  return(out_df)
}

####SIRS model

stepSIRS <- function(S,I,R,parameters){
  N <- parameters[[1]]
  beta <- parameters[[2]]
  gam <- parameters[[3]]
  rho <- parameters[[4]]
  alpha <- parameters[[5]]
  mu <- parameters[[6]]
  
  #t1 <- rbinom(n=1,size=S,prob=1-exp(-beta*I/N))
  t1 <- rbinom(n=1,size=S,prob=1-exp(-beta*I/(1+alpha*I^2)))
  t2 <- rbinom(n=1,size=I,prob=1-exp(-gam))
  t3 <- rbinom(n=1,size=R,prob=1-exp(-mu))
  
  S  <- S - t1 + t3;
  I  <- I - t2 + t1;
  R  <- R - t3 + t2;
  
  B <- rpois(n=1,lambda=rho*I+1e-6);
  
  return(c(S,I,R,B))
}


runSIRS <- function(parameters,nsteps,full_out=F){
  N <- parameters[[1]]
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


varyRho_SIRS <- function(len,pars,rho_seq){
  N <- pars[["N"]]
  beta <- pars[["beta"]]
  gam <- pars[["gamma"]]
  alpha <- pars[["alpha"]]
  mu <- pars[["mu"]]
  out_list <- list()
  counter <- 1
  for(rho in rho_seq){
    parameters <- list(N,beta,gam,rho,alpha,mu)
    run_raw <- runSIS(parameters,10000) 
    run_sym <- symbolizeOutput(run_raw,median(run_raw))
    gEE_out <- getExcessEntropy(seq=run_sym,maxL=20)
    E <- gEE_out[1]
    hmu_est <- gEE_out[2]
    out_list[[counter]] <- c(rho,hmu_est,E)
    counter <- counter+1
  }
  out_df <- as.data.frame(do.call(rbind.data.frame, out_list),stringsAsFactors=FALSE)
  names(out_df) <- c("rho","hmu","E")
  return(out_df)
}


varyBeta_SIRS <- function(len,pars,beta_seq){
  N <- pars[["N"]]
  rho <- pars[["rho"]]
  gam <- pars[["gamma"]]
  alpha <- pars[["alpha"]]
  mu <- pars[["mu"]]
  out_list <- list()
  counter <- 1
  for(beta in beta_seq){
    parameters <- list(N,beta,gam,rho,alpha,mu)
    run_raw <- runSIS(parameters,10000) 
    run_sym <- symbolizeOutput(run_raw,median(run_raw))
    gEE_out <- getExcessEntropy(seq=run_sym,maxL=20)
    E <- gEE_out[1]
    hmu_est <- gEE_out[2]
    out_list[[counter]] <- c(beta,hmu_est,E)
    counter <- counter+1
  }
  out_df <- as.data.frame(do.call(rbind.data.frame, out_list),stringsAsFactors=FALSE)
  names(out_df) <- c("beta","hmu","E")
  return(out_df)
}