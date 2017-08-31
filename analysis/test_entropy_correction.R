#JLW 08/31/17

library(entropy)

setwd("/home/jake/memory")
source("./src/entropy_calc.R")


# General conclusions: 
# 1. the grassberger03 method appears to perform best. 
# 2. CS may beat is for very high block lengths but appears to behave inconsistently (drops then peaks when it should be flat)
# 3. the grassberger88 method (default) performs nearly as well as the grassberger03 method for these simple cases
# All methods (excluding perhaps bonachela) do better than the basic plug-in method for calculating the entropy at large block lengthd

# A note: I have left grassberger88 as the default rather than grassberger03 because I am not 100% confident in my implementation in the latter
# (someone should double check I got the calculation right)

#################
# Weighted Coin #
#################

p <- 0.6
n <- 1e4
x <- rbinom(n,size=1,prob=p)
x <- paste(x,collapse="")

hmu_curve <- getEntropyCurve(x,20)
hmu_curve_grassberger88 <- getEntropyCurve.corrected(x,20) #grassberger88 by default
hmu_curve_grassberger03 <- getEntropyCurve.corrected(x,20,method="grassberger03") #slight improvement over grassberger88
hmu_curve_bonachela <- getEntropyCurve.corrected(x,20,method="bonachela") #Bonachela doesn't work well w/ low probabilities... fails at med-large block size (overestimates instead of underestimating)

#Also various estimators from the "entropy" package
hmu_curve_ML <- getEntropyCurve.corrected(x,20,method="ML") #Just the base method (plug-in)
hmu_curve_MM <- getEntropyCurve.corrected(x,20,method="MM")
hmu_curve_Jeffreys <- getEntropyCurve.corrected(x,20,method="Jeffreys")
hmu_curve_Laplace <- getEntropyCurve.corrected(x,20,method="Laplace")
hmu_curve_SG <- getEntropyCurve.corrected(x,20,method="SG")
hmu_curve_minimax <- getEntropyCurve.corrected(x,20,method="minimax")
hmu_curve_CS <- getEntropyCurve.corrected(x,20,method="CS") #Maybe promising?
#hmu_curve_NSB <- getEntropyCurve.corrected(x,20,method="NSB") #Currently giving an error -> investigate
hmu_curve_shrink <- getEntropyCurve.corrected(x,20,method="shrink")

plot(hmu_curve,type="l",col="red",ylim=c(0,2),xlab="L",ylab=expression(h[mu](L)))
lines(hmu_curve_grassberger88,col="blue")
lines(hmu_curve_grassberger03,col="magenta")
lines(hmu_curve_bonachela,col="green") 
lines(hmu_curve_ML) 
lines(hmu_curve_MM) 
lines(hmu_curve_Jeffreys) 
lines(hmu_curve_Laplace) 
lines(hmu_curve_SG) 
lines(hmu_curve_minimax) 
lines(hmu_curve_CS) 
lines(hmu_curve_NSB) 
lines(hmu_curve_shrink) 

plot(hmu_curve,type="l",col="red",ylim=c(0,1),xlab="L",ylab=expression(h[mu](L)))
lines(hmu_curve_grassberger88,col="blue")
lines(hmu_curve_grassberger03,col="magenta")
lines(hmu_curve_ML) 
lines(hmu_curve_MM) 
lines(hmu_curve_Jeffreys) 
lines(hmu_curve_Laplace) 
lines(hmu_curve_SG) 
lines(hmu_curve_minimax) 
lines(hmu_curve_CS) 
#lines(hmu_curve_NSB) 
lines(hmu_curve_shrink) 

###########
# 2-cycle # #Should all be zero->good
###########

x <- rep(c(0,1),1e4)
x <- paste(x,collapse="")
hmu_curve <- getEntropyCurve(x,20)
hmu_curve_grassberger88 <- getEntropyCurve.corrected(x,20) #grassberger88 by default
hmu_curve_grassberger03 <- getEntropyCurve.corrected(x,20,method="grassberger03") 
hmu_curve_ML <- getEntropyCurve.corrected(x,20,method="ML")
hmu_curve_MM <- getEntropyCurve.corrected(x,20,method="MM")
hmu_curve_Jeffreys <- getEntropyCurve.corrected(x,20,method="Jeffreys")
hmu_curve_Laplace <- getEntropyCurve.corrected(x,20,method="Laplace")
hmu_curve_SG <- getEntropyCurve.corrected(x,20,method="SG")
hmu_curve_minimax <- getEntropyCurve.corrected(x,20,method="minimax")
hmu_curve_CS <- getEntropyCurve.corrected(x,20,method="CS")
hmu_curve_shrink <- getEntropyCurve.corrected(x,20,method="shrink")

plot(hmu_curve,type="l",col="red",ylim=c(0,1),xlab="L",ylab=expression(h[mu](L)))
lines(hmu_curve_grassberger88,col="blue")
lines(hmu_curve_grassberger03,col="magenta")
lines(hmu_curve_ML) 
lines(hmu_curve_MM) 
lines(hmu_curve_Jeffreys) 
lines(hmu_curve_Laplace) 
lines(hmu_curve_SG) 
lines(hmu_curve_minimax) 
lines(hmu_curve_CS) 
lines(hmu_curve_shrink)

###############
# Golden mean #
###############

x <- sample(c("1","01"),1e4,replace=T)
x <- paste(x,collapse="")
hmu_curve <- getEntropyCurve(x,20)
hmu_curve_grassberger88 <- getEntropyCurve.corrected(x,20) #grassberger88 by default
hmu_curve_grassberger03 <- getEntropyCurve.corrected(x,20,method="grassberger03") 
hmu_curve_ML <- getEntropyCurve.corrected(x,20,method="ML")
hmu_curve_MM <- getEntropyCurve.corrected(x,20,method="MM")
hmu_curve_Jeffreys <- getEntropyCurve.corrected(x,20,method="Jeffreys")
hmu_curve_Laplace <- getEntropyCurve.corrected(x,20,method="Laplace")
hmu_curve_SG <- getEntropyCurve.corrected(x,20,method="SG")
hmu_curve_minimax <- getEntropyCurve.corrected(x,20,method="minimax")
hmu_curve_CS <- getEntropyCurve.corrected(x,20,method="CS")
hmu_curve_shrink <- getEntropyCurve.corrected(x,20,method="shrink")

plot(hmu_curve,type="l",col="red",ylim=c(0,0.5),xlab="L",ylab=expression(h[mu](L)))
lines(hmu_curve_grassberger88,col="blue")
lines(hmu_curve_grassberger03,col="magenta")
lines(hmu_curve_ML) 
lines(hmu_curve_MM) 
lines(hmu_curve_Jeffreys) 
lines(hmu_curve_Laplace) 
lines(hmu_curve_SG) 
lines(hmu_curve_minimax) 
lines(hmu_curve_CS) 
lines(hmu_curve_shrink)


#Can then find E using this curve
gEE_out <- getExcessEntropy(hmu=hmu_curve_grassberger03)
E <- gEE_out[1]
hmu_est <- gEE_out[2]
cutoff_ind <- gEE_out[3]
plot(hmu_curve_grassberger03)
abline(h=hmu_est)
abline(v=cutoff_ind)

plot(hmu_curve_grassberger03*(1:20),xlab="L",ylab=expression(H(L)))




