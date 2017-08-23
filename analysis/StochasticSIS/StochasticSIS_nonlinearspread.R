
setwd("/media/jlweissman/OS/Storage/Grad/LabNotebook/CSSS/Memory/memory-master")
source("./src/StochasticSIS_nonlinearspread_func.R")
source("./src/entropy_calc.R")
source("./src/run_cssr_by_min_bic.R")

N <- 800
beta <- 2
gam <- 5
rho <- 0.9
alpha <- 1
parameters <- list(N,beta,gam,rho,alpha)
test <- runSIS(parameters,1000) #full_out=TRUE for S,I in output

plot(1:1000,test,type="l")
hist(test,breaks=100)
abline(v=median(test))
median(test)

test_sym <- symbolizeOutput(test,median(test))
test_sym
dat <- test_sym
out_file <- "example.out"
alpha <- paste0(sort(unique(unlist(strsplit(dat, split="")))), collapse="")
cssr_out <- run_cssr_by_min_bic(alpha, dat, out_file, threshold = 0.1)
