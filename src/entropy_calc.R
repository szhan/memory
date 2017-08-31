# author: jlweissman
#Everything is in nats. Deal with it.

allsubstr <- function(x, n){substring(x, 1:(nchar(x) - n + 1), n:nchar(x))}

shannonH <- function(p){-sum(p*log(p))}

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

#Bias correction for estimated shannon entropies from Grassberger 1988 Phys. Lett. A 
#also see Grassberber 2003 arXiv for a better(?) correction and Bonachela et. al 2008 J. Phys. A: Math Theor. for another alternative
# or Nemenmann, Shafee, Bialek 2002 for a Bayesian approach
#n is a vector of empirical counts
shannonH.grassberger88 <- function(n){
  N <- sum(n)
  chi_n <- -(n/N)*log(n/N) + (1/(2*N)) - ((-1)^n)/(N*(n+1))
  return(sum(chi_n))
}

partialHarmonic <- function(start,end){return(sum(1/(start:end)))}

shannonH.bonachela <- function(n){
  N <- sum(n)
  return((1/(N+2))*sum((n+1)*sapply(n+2,partialHarmonic,end=N+2)))
}

getRecurrence <- function(n_max){
  em_gamma <- -digamma(1)
  g <- c(-em_gamma-log(2),2-em_gamma-log(2),2-em_gamma-log(2))
  counter <- 3
  n <- 1
  for(i in 1:(floor(n_max/2)-1)){
    g[(counter+1):(counter+2)] <- g[counter]+2/(2*n+1)
    counter <- counter+2
    n <- n+1
  }
  return(g)
}

g03_chi <- function(ni,N,g){return((ni/N)*(log(N)-g[ni]))}

shannonH.grassberger03 <- function(n,g){return(sum(sapply(n,g03_chi,N=sum(n),g=g)))}

getHL.corrected <- function(L,seq,method,g=0){
  n <- as.vector(table(allsubstr(seq,L)))
  
  if (method=="grassberger88"){
    return(shannonH.grassberger88(n))
    
  } else if (method=="grassberger03"){
    return(shannonH.grassberger03(n,g))
    
  } else if (method=="bonachela"){
    return(shannonH.bonachela(n))
    
  } else if (method %in% c("ML","MM","Jeffreys","Laplace","SG","minimax","CS","NSB","shrink")){
    return(entropy(n,method=method))
    
  } else {
    stop("Not a valid method")
  }

}

getEntropyCurve.corrected <- function(seq,maxL,method="grassberger88"){
  if(method=="grassberger03"){
    n_max <- max(as.vector(table(allsubstr(seq,1))))
    g <- getRecurrence(n_max)
    HL <- lapply(1:maxL,getHL.corrected,seq=seq,method=method,g=g)
  } else {
    HL <- lapply(1:maxL,getHL.corrected,seq=seq,method=method)
  }
  
  hmu <- diff(unlist(HL))
  return(hmu)
}

get.hmu.corrected <- function(i,w,seq,plotit=FALSE){
  maxL <- 20
  seq <- substr(seq,i,i+w-1)
  hmu <- getEntropyCurve.corrected(seq,maxL)
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

get.vw.w <- function(x,interval_jump,n){
  intervals <- seq(100,nchar(x)-interval_jump-1,interval_jump)
  hmu_var_vec <- numeric(length(intervals))
  counter <- 1
  for(i in intervals){
    hmu_var_vec[counter] <- var.hmu.para.samp(x,i,n)
    counter <- counter+1
  }
  plot(intervals,hmu_var_vec,xlab="w",ylab=expression(paste("Var(",h[mu],")")),main="100 Samples Per w",ylim=c(0,0.03))
  return(hmu_var_vec)
}

generateDataFromEM <- function(best_res,n_dat=10000){
  tp <- best_res$results$matrices$trans_prob
  te <- best_res$results$matrices$trans_stat + 1
  sp <- best_res$results$matrices$state_prob
  
  data <- numeric(n_dat)
  draw <- runif(1)
  state <- min(which(cumsum(sp)>=draw))
  for(i in 1:n_dat){
    draw <- runif(1)
    emission <- min(which(cumsum(tp[state,])>=draw))
    data[i] <- emission-1
    state <- te[state,emission]
  }
  data <- paste(data,collapse="")
  
  return(data)
}

getExcessEntropy <- function(seq=0,maxL=0,hmu=0){
  if(hmu==0){
    if(seq!=0 && maxL!=0){
      hmu <- getEntropyCurve(seq,maxL)
    } else {
      stop('Must provide sequence or pre-calculated hmu list')
    }
  } else {
    maxL <- length(hmu)
  }

  cutoff_ind <- which(diff(diff(hmu))<0)[1] #Where the concavity changes
  if(cutoff_ind==0){cutoff_ind <- maxL}
  hmu_est <- hmu[cutoff_ind]
  E <- sum(hmu[1:cutoff_ind]-hmu_est)
  return(c(E,hmu_est,cutoff_ind))
}


