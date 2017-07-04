#library(devtools)
#install_github("szhan/CSSR")

require(CSSR)
require(DOT)
library(pryr)
library(foreach)
library(doParallel)
cl<-makeCluster(4)
registerDoParallel(cl)

print(mem_used())

allsubstr <- function(x, n){substring(x, 1:(nchar(x) - n + 1), n:nchar(x))}

shannonH <- function(p){-sum(p*log2(p))}

getHL <- function(L,seq){
  a <- table(allsubstr(seq,L))
  return(shannonH(a/sum(a)))
}

getEntropyCurve <- function(seq,maxL){
  HL <- lapply(1:maxL,getHL,seq=seq)
  hmu <- diff(unlist(HL))
  return(hmu)
}

get.hmu <- function(i,w,seq,plotit=FALSE){
  maxL <- 20
  seq <- substr(seq,i,i+w-1)
  hmu <- getEntropyCurve(seq,maxL)
  max_ind <- which(abs(diff(hmu)) %in% max(abs(diff(hmu))))
  min_ind <- which(abs(diff(hmu[1:max_ind])) %in% min(abs(diff(hmu[1:max_ind]))))
  
  if(plotit==TRUE){
    par(mfrow=c(2,1))
    plot(1:(min_ind),hmu[1:min_ind])
    plot(hmu)
    par(mfrow=c(1,1)) 
  }
  
  return(hmu[min_ind])
}

var.hmu <- function(sequence,w){
  print(nchar(sequence)-w+1)
  lapply(1:(nchar(sequence)-w+1),get.hmu,w=w,seq=sequence)
  hmu_vec <-unlist(lapply(1:(nchar(sequence)-w+1),get.hmu,w=w,seq=sequence))
  return(var(hmu_vec))
}

var.hmu.para <- function(sequence,w){
  hmu_vec = foreach(i=1:(nchar(sequence)-w+1),.export=ls(envir=globalenv()),.combine="c") %dopar% {  
    get.hmu(i,w,sequence)
  }
  return(var(hmu_vec))
}

var.hmu.para.samp <- function(sequence,w,n){
  window_sample <- sample(1:(nchar(sequence)-w+1),n,replace=FALSE)
  hmu_vec = foreach(i=1:length(window_sample),.export=ls(envir=globalenv()),.combine="c") %dopar% {  
    get.hmu(window_sample[i],w,sequence)
  }
  return(var(hmu_vec))
}

#######################################

#Non stationary#
pGen <- function(t,z){(0.5*t^3)/(t^3+z^3)}
t <- 1:2000
p_coin <- pGen(t,500)
plot(p_coin,xlab="time",ylab="p")
abline(v=1000)

p_coin <- (numeric(5000)+1)*.5
x <- rbinom(length(p_coin),1,p_coin)
x <- paste(x,collapse="")

hmu <- getEntropyCurve(x,20)
plot(abs(diff(hmu)))
plot(diff(hmu))
plot(hmu)


get.hmu(1,2000,x,plotit=TRUE)

#tnp <- system.time(var.hmu(x,1900))
#tp <- system.time(var.hmu.para(x,1900)) #Better
#####################################
n <- 100
intervals <- seq(100,1900,100)
hmu_var_vec <- numeric(length(intervals))
counter <- 1
for(i in intervals){
  hmu_var_vec[counter] <- var.hmu.para.samp(x,i,n)
  counter <- counter+1
}
plot(intervals,hmu_var_vec,xlab="w",ylab=expression(paste("Var(",h[mu],")")),main="100 Samples Per w")




























