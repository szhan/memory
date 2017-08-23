library(CSSR)

setwd("/media/jlweissman/OS/Storage/Grad/LabNotebook/CSSS/Memory/memory")
source("./src/StochasticSIS_nonlinearspread_func.R")
source("./src/entropy_calc.R")
source("./src/run_cssr_by_min_bic.R")

N <- 800
beta <- 2
gam <- 5
rho <- 0.9
alpha <- 1
parameters <- list(N,beta,gam,rho,alpha)
test <- runSIS(parameters,100000) #full_out=TRUE for S,I in output

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




