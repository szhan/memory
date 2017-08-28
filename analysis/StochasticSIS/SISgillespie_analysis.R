library(CSSR)


#setwd("/media/jlweissman/OS/Storage/Grad/LabNotebook/CSSS/Memory/memory")
source("./src/SIS_gillespie.R")
source("./src/entropy_calc.R")
source("./src/run_cssr_by_min_bic.R")

params = list(z=1e-6, beta=1, gamma=1e-1, rho=1, alpha=1)
sis_out <- SISg(params,sim_length=10000,S0=1e5,I0=0)
plot(sis_out$time,sis_out$I)
plot(sis_out$time,sis_out$S)
