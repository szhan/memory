library(CSSR)


#setwd("/media/jlweissman/OS/Storage/Grad/LabNotebook/CSSS/Memory/memory")
source("./src/StochasticSIS_nonlinearspread_func.R")
source("./src/entropy_calc.R")
source("./src/run_cssr_by_min_bic.R")


pars <- list(N=800,beta=2.7,gamma=5,rho=1,alpha=1e-3)
test <- runSIS(pars,100000) #full_out=TRUE for S,I in output
plot(1:100000,test[1:100000],type="l")
abline(h=median(test),col="red")
abline(h=quantile(test,0.6),col="red")
abline(h=quantile(test,0.4),col="red")
abline(h=quantile(test,0.75),col="red")
abline(h=quantile(test,0.25),col="red")

hist(test,breaks=100)
abline(v=median(test))
median(test)

test_sym <- symbolizeOutput(test,0.75,extremes = 1)
hmu <- getEntropyCurve(test_sym,20)
plot(abs(diff(hmu)))
plot(diff(diff(hmu)))
abline(h=0)
plot(hmu)

#getExcessEntropy(seq=test_sym,maxL=20) #From scratch, slower if hmu list already calculated
gEE_out <- getExcessEntropy(hmu=hmu)
E <- gEE_out[1]
hmu_est <- gEE_out[2]
cutoff_ind <- gEE_out[3]
plot(hmu)
abline(h=hmu_est)
abline(v=cutoff_ind)

pars <- list(N=800,beta=2,gamma=5,alpha=1)
len <- 100000
rho_out <- varyRho(len,pars,rho_seq=c(0.1,0.5,0.9,1))
plot(rho_out$hmu,rho_out$E,xlab=expression(h[mu]),ylab="E")
lines(rho_out$hmu,rho_out$E)

pars <- list(N=800,rho=1,gamma=5,alpha=1e-3)
len <- 100000
beta_out <- varyBeta(len,pars,beta_seq=seq(0.1,5,0.3))
plot(beta_out$hmu,beta_out$E,xlab=expression(h[mu]),ylab="E")
lines(beta_out$hmu,beta_out$E)

pars <- list(N=800,rho=1,gamma=5,alpha=1e-3)
beta_out <- varyParameter_SIS(parameter_name="beta",parameter_seq=seq(0.1,5,0.3),parameters_list=pars,sim_length=10000,outbreak_percentile=0.75)
plot(beta_out$hmu,beta_out$E,xlab=expression(h[mu]),ylab="E")
lines(beta_out$hmu,beta_out$E)

pars <- list(N=800,beta=2.7,gamma=5,alpha=1e-3)
rho_out <- varyParameter_SIS(parameter_name="rho",parameter_seq=seq(0.1,1,0.1),parameters_list=pars,sim_length=10000,outbreak_percentile=0.75)
plot(rho_out$hmu,rho_out$E,xlab=expression(h[mu]),ylab="E")
lines(rho_out$hmu,rho_out$E)

pars <- list(N=800,beta=2.7,rho=1,alpha=1e-3)
gamma_out <- varyParameter_SIS(parameter_name="gamma",parameter_seq=seq(0.1,6,0.3),parameters_list=pars,sim_length=10000,outbreak_percentile=0.75)
plot(gamma_out$hmu,gamma_out$E,xlab=expression(h[mu]),ylab="E")
lines(gamma_out$hmu,gamma_out$E)










##### SIRS

N <- 10000
beta <- 2
gam <- 5
rho <- 1
alpha <- 1
mu <- 1
parameters <- list(N,beta,gam,rho,alpha,mu)
test <- runSIRS(parameters,10000) #full_out=TRUE for S,I in output

plot(10:1000,test[10:1000],type="l")
abline(h=median(test))
abline(h=quantile(test,0.6))
abline(h=quantile(test,0.4))
abline(h=quantile(test,0.75))
abline(h=quantile(test,0.25))
test_sym <- symbolizeOutput(test,0.75,extremes = 1)

hist(test,breaks=100)
abline(v=median(test))
median(test)

#test_sym <- symbolizeOutput(test,0.5)
hmu <- getEntropyCurve(test_sym,20)
plot(abs(diff(hmu)))
plot(diff(hmu))
plot(diff(diff(hmu)))
abline(h=0)
plot(hmu)

gEE_out <- getExcessEntropy(hmu=hmu)
E <- gEE_out[1]
hmu_est <- gEE_out[2]
cutoff_ind <- gEE_out[3]
plot(hmu,xlab="L",ylab=expression(h[mu](L)))
abline(h=hmu_est)
abline(v=cutoff_ind)

pars <- list(N=800,beta=2,gamma=5,alpha=1,mu=0.1)
len <- 10000000
rho_out <- varyRho_SIRS(len,pars,rho_seq=seq(0.1,1,0.1))
plot(rho_out$hmu,rho_out$E,xlab=expression(h[mu]),ylab="E")
lines(rho_out$hmu,rho_out$E)


pars <- list(N=800,rho=1,gamma=5,alpha=1e-3,mu=10)
len <- 10000
beta_out <- varyBeta_SIRS(len,pars,beta_seq=seq(0.1,5,0.3))
plot(beta_out$hmu,beta_out$E,xlab=expression(h[mu]),ylab="E")
lines(beta_out$hmu,beta_out$E)


pars <- list(N=800,rho=1,gamma=5,alpha=1e-3,mu=0.1)
len <- 100000
beta_out <- varyParameter_SIRS(parameter_name="beta",parameter_seq=seq(0.1,3,0.3),parameters_list=pars,sim_length=10000,outbreak_percentile=0.75)
plot(beta_out$hmu,beta_out$E,xlab=expression(h[mu]),ylab="E")
lines(beta_out$hmu,beta_out$E)



f#library(devtools)
#install_github("dimalik/Hrate")
library(Hrate)
alpha <- unlist(strsplit(test_sym, split=""))
stabilization.rate <- stabilize.estimate(text = alpha, step.size = 100000, max.length = 1000000, every.word = 10, method="downsample",rate = 5)

print(get.stabilization(stabilization.rate))
print(get.criterion(stabilization.rate))

estimate <- get.estimate(text = alpha, every.word = 10, max.length = 50000)

plot(convergence.rate)


