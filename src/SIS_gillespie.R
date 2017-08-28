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

SIRSgRates <- function(x, p, t) {
  return(c(x["S"]*(params$beta*x["I"] + params$z), #Infection rate 
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

