library(CSSR)
library(adaptivetau)

setwd("/media/jlweissman/OS/Storage/Grad/LabNotebook/CSSS/Memory/memory")
source("./src/SIS_gillespie.R")
source("./src/entropy_calc.R")
source("./src/run_cssr_by_min_bic.R")

#SIS

#SISg() for adaptive tau-leaping algorithm (fast)
#SISg.exact() for gillespie algorithm (slow)

params = list(z=1e-6, beta=1e-5, gamma=1e-3, rho=1, alpha=1e-6)
sis_out <- SISg(params,sim_length=10000,S0=800,I0=1)

plot(sis_out$time,sis_out$S,col="blue",ylim=c(1,1e3),type="l")
lines(sis_out$time,sis_out$I,col="red")


##################
# Stochastic SIR #
##################

#Gives stochastic amplification
#params <-  list(z=1e-5, beta=(1e-4)*(1.2), gamma=5e-1, rho=0.5, alpha=0,mu=1e-4)
#sirs_out <- SIRSg(params,sim_length=10000,S0=2e6,I0=1)

#Gives sharp outbreaks of ~1000, periods where cycles between outbreak and clearance as well as endemic periods
#params <-  list(z=1e-4, beta=(1e-3)*(1.2), gamma=5e0, rho=0.5, alpha=0,mu=1e-3)
#sirs_out <- SIRSg(params,sim_length=10000,S0=5e5,I0=1)

params <-  list(z=1e-4, beta=(1e-3)*(1.2), gamma=5e0, rho=0.5, alpha=0,mu=1e-3)
sirs_out <- SIRSg(params,sim_length=10001,S0=5e5,I0=1)

plot(sirs_out$time,sirs_out$S,col="blue",ylim=c(1,1e6),type="l",log="y",xlim=c(200,600),xlab="Time",ylab="Individuals",main="Stochastic SIR with Immigration")
lines(sirs_out$time,sirs_out$I,col="red")
lines(sirs_out$time,sirs_out$R,col="green")
lines(sirs_out$time,sirs_out$B,col="magenta")
legend(550,8e5,legend=c("S","I","R","B"),fill=c("blue","red","green","magenta"))

plot(sirs_out$time,sirs_out$I,col="red",ylim=c(0,2000),xlim=c(0,1000),type="l")
abline(h=500)

plot(sirs_out$time,sirs_out$I,col="red",ylim=c(0,2000),xlim=c(200,600),type="l",xlab="Time",ylab="Infected",main="Stochastic SIR with Immigration")
plot(sirs_out$time,sirs_out$I,col="red",ylim=c(0,500),xlim=c(200,600),type="l",xlab="Time",ylab="Infected",main="Stochastic SIR with Immigration")

so_y <- approx(sirs_out$time,sirs_out$I,xout=1:(floor(sirs_out$time[length(sirs_out$time)])-1))
plot(so_y$x,so_y$y,type="l",xlim=c(200,600))

par(mfrow=c(2,1))
plot(sirs_out$time,sirs_out$I,xlim=c(400,600),ylim=c(0,4000),type="l")
abline(h=500)
abline(h=5)
plot(so_y$x,so_y$y,type="l",xlim=c(400,600))
abline(h=500)
abline(h=5)
par(mfrow=c(1,1))

par(mfrow=c(2,1))
plot(sirs_out$time,sirs_out$I,xlim=c(400,600),ylim=c(0,10),type="l")
abline(h=5)
plot(so_y$x,so_y$y,type="l",xlim=c(400,600),ylim=c(0,10))
abline(h=5)
par(mfrow=c(1,1))

#Need lots of data
params <-  list(z=1e-4, beta=(1e-3)*(1.2), gamma=5e0, rho=0.5, alpha=0,mu=1e-3)
sirs_out <- SIRSg(params,sim_length=1e5,S0=5e5,I0=1)

x <- symbolizeSIR(sirs_out,extinction_thresh=5,epidemic_thresh=500)

hmu <- getEntropyCurve(x,50)
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




#############
# SIE Model #
#############

params = list(z=0, r=1.5, K=12, beta=0.8, m=0.5, epsilon=1e-1, tau=1e-2, rho=1)
sie_out <- SIEg(params,sim_length=100,S0=params$K,I0=1,E0=0)

par(mfrow=c(3,1))
plot(sie_out$time,sie_out$S,col="blue",type="l")
plot(sie_out$time,sie_out$I,col="red",type="l")
plot(sie_out$time,sie_out$E,col="green",type="l")
par(mfrow=c(1,1))

params = list(z=0, r=1.5, K=3000, beta=0.1, m=0.567, epsilon=0.111+0.4, tau=1, rho=1)
sie_out <- SIEg.exact(params,sim_length=100,S0=params$K,I0=5,E0=0)

par(mfrow=c(3,1))
plot(sie_out$time,sie_out$S,col="blue",type="l")
plot(sie_out$time,sie_out$I,col="red",type="l")
plot(sie_out$time,sie_out$E,col="green",type="l")
par(mfrow=c(1,1))


###########################################################################################
# Seasonally forced stochastic SIR (Black and McKane 2010 Journal of Theoretical Biology) #
###########################################################################################

params <-  list(z=1e-5, beta0=1e-4, beta1=0.2, gamma=5e-1, rho=0.5, alpha=0,mu=1e-4)
sirs_out <- SIRSg_forced(params,sim_iter=8000,on_interval=1,off_interval=0,S0=2e6,I0=1,R0=0)

plot(sirs_out$time,sirs_out$S,col="blue",type="l",log="y")
lines(sirs_out$time,sirs_out$I,col="red")
lines(sirs_out$time,sirs_out$R,col="green")
lines(sirs_out$time,sirs_out$B,col="magenta")

plot(sirs_out$time,sirs_out$I,col="red",ylim=c(0,1000),xlim=c(0,4000),type="l")





