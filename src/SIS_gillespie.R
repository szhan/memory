# Stochastic SIS and SIRS models simulated w/ adaptive tau leaping algorithm from adaptivetau package
#
#Nonlinear spread (bets*S*I/(1+alpha*I^2)), so that spread decreases with increasing infectives (e.g. due to more resources spent on control),
#added to get interesting dynamics (outbreaks, no stable eq) -> See Liu, Levin, Iwasa 1986 for complete discussion, more general forms
#also has an "observation" component, determined by parameter rho, i.e. our probability of detecting infected individuals
# Thus can change disease dynamics (spread speed bet, recovery rate gamma, strength of intervention alpha) or observation intensity (rho)
#JLW 28/8/17


#install.packages("adaptivetau")
#library(adaptivetau)
#params = list(z=1e-6, beta=1e-5, mu=1e-2, gamma=1e-1)

#######
# SIS #
#######

#Linear spread#
#SISgRates <- function(x, p, t) {
#  return(c(x["S"]*(params$beta*x["I"] + params$z), #Infection
#           params$gamma*x["I"])) #Recovery
#}

SISgRates <- function(x, p, t) {
  return(c(x["S"]*(params$beta*x["I"]/(1+params$alpha*x["I"]^2) + params$z), #Infection
           params$gamma*x["I"])) #Recovery
}

SISg <- function(params,sim_length=10000,S0=1e5,I0=0){
  init.values <- c(S=S0,I=I0)
  transitions <-   list(c(S = -1, I = +1), c(I = -1, S = +1))
  sis_out <- ssa.adaptivetau(init.values, transitions, SISgRates, params, tf=sim_length)
  sis_out <- as.data.frame(sis_out,stringsAsFactors=FALSE)
  sis_out$B <- rpois(n=length(sis_out$I),lambda=(sis_out$I*params$rho+1e-6))
  return(sis_out)
}

SISg.exact <- function(params,sim_length=10000,S0=1e5,I0=0){
  init.values <- c(S=S0,I=I0)
  transitions <-   list(c(S = -1, I = +1), c(I = -1, S = +1))
  sis_out <- ssa.exact(init.values, transitions, SISgRates, params, tf=sim_length)
  sis_out <- as.data.frame(sis_out,stringsAsFactors=FALSE)
  sis_out$B <- rpois(n=length(sis_out$I),lambda=(sis_out$I*params$rho+1e-6))
  return(sis_out)
}


SISgRates_logistic <- function(x, p, t) {
  return(c(x["S"]*(params$beta*x["I"]/(1+params$alpha*x["I"]^2) + params$z), #Infection
           params$gamma*x["I"], #Recovery
           (params$r*(1-x["S"]/params$K))*x["S"]*(params$r*(1-x["S"]/params$K) > 0), #Density dependent growth
           -(params$r*(1-x["S"]/params$K))*x["S"]*(params$r*(1-x["S"]/params$K) < 0), #Density-dependent death
           params$m*x["I"])) #Mortality of diseased individuals
}

SISg_logistic <- function(params,sim_length=10000,S0=1e5,I0=0){
  init.values <- c(S=S0,I=I0)
  transitions <-   list(c(S = -1, I = +1), c(I = -1, S = +1), c(S = +1), c(S = -1), c(I =-1))
  sis_out <- ssa.adaptivetau(init.values, transitions, SISgRates_logistic, params, tf=sim_length)
  sis_out <- as.data.frame(sis_out,stringsAsFactors=FALSE)
  sis_out$B <- rpois(n=length(sis_out$I),lambda=(sis_out$I*params$rho+1e-6))
  return(sis_out)
}




########
# SIRS #
########

#Linear spread#
#SIRSgRates <- function(x, p, t) {
#  return(c(x["S"]*(params$beta*x["I"] + params$z), #Infection rate 
#           params$gamma*x["I"], # recovery rate
#           params$mu*x["R"])) #Immune loss rate
#}


SIRSgRates <- function(x, p, t) {
  return(c(x["S"]*(params$beta*x["I"]/(1+params$alpha*x["I"]^2) + params$z), #Infection rate 
           params$gamma*x["I"], # recovery rate
           params$mu*x["R"])) #Immune loss rate
}

SIRSg <- function(params,sim_length=10000,S0=1e5,I0=0,R0=0){
  init.values <- c(S=S0,I=I0,R=R0)
  transitions <-   list(c(S = -1, I = +1), c(I = -1, R = +1), c(R = -1, S = +1))
  sirs_out <- ssa.adaptivetau(init.values, transitions, SIRSgRates, params, tf=sim_length)
  sirs_out <- as.data.frame(sirs_out,stringsAsFactors=FALSE)
  sirs_out$B <- rpois(n=length(sirs_out$I),lambda=params$rho*sirs_out$I+1e-6)
  return(sirs_out)
}


########################################################################################
# SIE Model (environmental pool, see Sharp and Pastor 2011 in Ecological Applications) #
########################################################################################

SIEgRates <- function(x, p, t) {
  return(c(x["S"]*(params$beta*x["E"] + params$z), #Infection
           (params$r*(1-x["S"]/params$K))*x["S"]*(params$r*(1-x["S"]/params$K) > 0), #Density dependent growth
           -(params$r*(1-x["S"]/params$K))*x["S"]*(params$r*(1-x["S"]/params$K) < 0), #Density-dependent death
           params$m*x["I"], #Mortality of diseased individuals
           params$epsilon*x["I"],
           params$tau*x["E"]))
}

SIEg <- function(params,sim_length=10000,S0=1e5,I0=0,E0=0){
  init.values <- c(S=S0,I=I0,E=E0)
  transitions <-   list(c(S = -1, I = +1), c(S = +1), c(S = -1), c(I = -1), c(E = +1), c(E = -1))
  sis_out <- ssa.adaptivetau(init.values, transitions, SIEgRates, params, tf=sim_length)
  sis_out <- as.data.frame(sis_out,stringsAsFactors=FALSE)
  sis_out$B <- rpois(n=length(sis_out$I),lambda=(sis_out$I*params$rho+1e-6))
  return(sis_out)
}


SIEg.exact <- function(params,sim_length=10000,S0=1e5,I0=0,E0=0){
  init.values <- c(S=S0,I=I0,E=E0)
  transitions <-   list(c(S = -1, I = +1), c(S = +1), c(S = -1), c(I = -1), c(E = +1), c(E = -1))
  sis_out <- ssa.exact(init.values, transitions, SIEgRates, params, tf=sim_length)
  sis_out <- as.data.frame(sis_out,stringsAsFactors=FALSE)
  sis_out$B <- rpois(n=length(sis_out$I),lambda=(sis_out$I*params$rho+1e-6))
  return(sis_out)
}


###########################################################################################
# Seasonally forced stochastic SIR (Black and McKane 2010 Journal of Theoretical Biology) #
###########################################################################################

SIRSgRates_forced <- function(x, params, t) {
  return(c(x["S"]*(params$beta*x["I"]/(1+params$alpha*x["I"]^2) + params$z), #Infection rate 
           params$gamma*x["I"], # recovery rate
           params$mu*x["R"])) #Immune loss rate
}


SIRSg_forced <- function(params,sim_iter,on_interval=1,off_interval=1,S0=1e5,I0=0,R0=0){
  
  interval_is_on <- 0
  #sim_iter #Number of cycles
  #on_interval #length of school year
  #off_interval  #length of summer break
  
  init.values <- c(S=S0,I=I0,R=R0)
  transitions <-   list(c(S = -1, I = +1), c(I = -1, R = +1), c(R = -1, S = +1))
  
  time_counter <- 0
  sirs_out <- list()
  for(i in 1:sim_iter){
    params$beta <- params$beta0*(1+params$beta1*interval_is_on)
    interval_length = on_interval*interval_is_on + off_interval*(1-interval_is_on)
    interval_is_on = 1-interval_is_on
    
    sirs_out_i <- ssa.adaptivetau(init.values, transitions, SIRSgRates_forced, params, tf=interval_length)
    sirs_out_i <- as.data.frame(sirs_out_i,stringsAsFactors=FALSE)
    sirs_out_i$B <- rpois(n=length(sirs_out_i$I),lambda=params$rho*sirs_out_i$I+1e-6)
    sirs_out_i$time <- sirs_out_i$time + time_counter
    init.values <- c(S=sirs_out_i$S[length(sirs_out_i$S)],I=sirs_out_i$I[length(sirs_out_i$I)],R=sirs_out_i$R[length(sirs_out_i$R)])
    sirs_out <- rbind(sirs_out,sirs_out_i)
    time_counter <- sirs_out$time[length(sirs_out$time)]
  }
  return(sirs_out)
}








