library(CSSR)

setwd("/home/jake/memory")
source("./src/StochasticSIS_nonlinearspread_func.R")
source("./src/entropy_calc.R")
source("./src/run_cssr_by_min_bic.R")

N <- 800
beta <- 2
gam <- 5
rho <- 0.9
alpha <- 1
parameters <- list(N,beta,gam,rho,alpha)
test <- runSIS(parameters,10000) #full_out=TRUE for S,I in output

plot(1:1000,test[1:1000],type="l")
hist(test,breaks=100)
abline(v=median(test))
median(test)

test_sym <- symbolizeOutput(test,median(test))
test_sym
#dat <- test_sym
#out_file <- "example.out"
#alpha <- paste0(sort(unique(unlist(strsplit(dat, split="")))), collapse="")
#cssr_out <- run_cssr_by_min_bic(alpha, dat, out_file, threshold = 0.1)

hmu <- getEntropyCurve(test_sym,20)
plot(abs(diff(hmu)))
plot(diff(hmu))
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
text(0.4,0.25,expression(paste(rho,"=0.1")))
text(0.25,0.16,expression(paste(rho,"=0.5")))
text(0.23,0.14,expression(paste(rho,"=0.9")))
text(0.215,0.11,expression(paste(rho,"=1.0")))


pars <- list(N=800,rho=1,gamma=5,alpha=1)
len <- 100000
beta_out <- varyBeta(len,pars,beta_seq=seq(0.1,1,0.3))
plot(beta_out$hmu,beta_out$E,xlab=expression(h[mu]),ylab="E")
lines(beta_out$hmu,beta_out$E)



##### SIRS

N <- 10000
beta <- 2
gam <- 5
rho <- 1
alpha <- 1
mu <- 1
parameters <- list(N,beta,gam,rho,alpha,mu)
test <- runSIRS(parameters,1000000) #full_out=TRUE for S,I in output

plot(10:1000,test[10:1000],type="l")
hist(test,breaks=100)
abline(v=median(test))
median(test)

test_sym <- symbolizeOutput(test,median(test))
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
len <- 100000
rho_out <- varyRho_SIRS(len,pars,rho_seq=c(0.1,0.5,0.9,1))
plot(rho_out$hmu,rho_out$E,xlab=expression(h[mu]),ylab="E")
lines(rho_out$hmu,rho_out$E)


pars <- list(N=800,rho=1,gamma=5,alpha=1,mu=0.1)
len <- 100000
beta_out <- varyBeta_SIRS(len,pars,beta_seq=seq(0.1,3,0.3))
plot(beta_out$hmu,beta_out$E,xlab=expression(h[mu]),ylab="E")
lines(beta_out$hmu,beta_out$E)


library(Hrate)
alpha <- unlist(strsplit(test_sym, split=""))
stabilization.rate <- stabilize.estimate(text = alpha, step.size = 1000, max.length = 50000, every.word = 10, method="downsample",rate = 5)

print(get.stabilization(stabilization.rate))
print(get.criterion(stabilization.rate))

estimate <- get.estimate(text = alpha, every.word = 10, max.length = 50000)

plot(convergence.rate)