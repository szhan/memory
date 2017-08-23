# author: jlweissman

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

  cutoff_ind <- which(diff(diff(hmu))<0)[1]-1 #Where the concavity changes
  if(cutoff_ind==0){cutoff_ind <- maxL}
  hmu_est <- hmu[cutoff_ind]
  E <- sum(hmu[1:cutoff_ind]-hmu_est)
  return(c(E,hmu_est,cutoff_ind))
}


