# Stochastic SIS model
#Nonlinear spread (bets*S*I/(1+alpha*I^2)), so that spread decreases with increasing infectives (e.g. due to more resources spent on control),
#added to get interesting dynamics (outbreaks, no stable eq) -> See Liu, Levin, Iwasa 1986 for complete discussion, more general forms
#also has an "observation" component, determined by parameter rho, i.e. our probability of detecting infected individuals
# Thus can change disease dynamics (spread speed beta, recovery rate gamma, strength of intervention alpha) or observation intensity (rho)
#JLW 23/8/17

stepSIS <- function(S,I,parameters){
  N <- parameters[[1]]
  beta <- parameters[[2]]
  gam <- parameters[[3]]
  rho <- parameters[[4]]
  alpha <- parameters[[5]]
  
  #t1 <- rbinom(n=1,size=S,prob=1-exp(-beta*I/N))
  t1 <- rbinom(n=1,size=S,prob=1-exp(-beta*I/(1+alpha*I^2)))
  t2 <- rbinom(n=1,size=I,prob=1-exp(-gam))
  
  S  = S - t1 + t2;
  I  = I - t2 + t1;
  
  B = rpois(n=1,lambda=rho*I+1e-6);
  
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

symbolizeOutput <- function(x,threshold){
  x <- (x>=threshold)*1
  x <- paste(x,collapse="")
  return(x)
}

