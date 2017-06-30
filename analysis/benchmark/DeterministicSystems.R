seasonF = function(t,A=1,m=1,O=0,k1=3,k2=1,Q=1){
  # NOTE:
  #  we use cosine rather than sine because cos(0) = 1
  #  such that the parameter O references the seasonal peak
  
  #  A = amplitude
  #  m = mean; if m<1 negative numbers are set to zero
  #  Q = number of peaks per year, can be < 1
  #  O = the peak (relative to 0)
  #  k1,k2  are shape parameters
  
  pmax(0, A*(m + cos(2*Q*pi*(t-O)/365))^k1, 0)^k2
}

seasonGenF = function(t, PAR){with(PAR,{
  sFi = function(i){seasonF(t,A[i],m[i],O[i],k1[i],k2[i],Q[i])} 
  apply(sapply(1:N,sFi),1,prod) 
})}


PAR = list(
  N=3,
  A=c(1,3,1),
  m=c(1,.5,1), 
  O=c(5,30,1), 
  k1=c(1,1,1), 
  k2=c(2,1,1), 
  Q =c(1,2,.5)
)

tt = c(0:7300)

plot(tt, seasonGenF(tt,PAR), type = "l")
