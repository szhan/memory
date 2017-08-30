
setwd("/media/jlweissman/OS/Storage/Grad/LabNotebook/CSSS/Memory/memory")
source("./src/entropy_calc.R")

p <- 0.5
n <- 1e4
x <- rbinom(n,size=1,prob=p)
x <- paste(x,collapse="")

hmu_curve <- getEntropyCurve(x,20)
hmu_curve_grassberger88 <- getEntropyCurve.corrected(x,20) #grassberger88 by default
hmu_curve_grassberger03 <- getEntropyCurve.corrected(x,20,method="grassberger03") 
hmu_curve_bonachela <- getEntropyCurve.corrected(x,20,method="bonachela") #Bonachela doesn't seem to work.... need to investigate

plot(hmu_curve,type="l",col="red",ylim=c(0,2))
lines(hmu_curve_grassberger88,col="blue")
lines(hmu_curve_grassberger03,col="magenta")
lines(hmu_curve_bonachela,col="green") 

plot(hmu_curve,type="l",col="red",ylim=c(0,1))
lines(hmu_curve_grassberger88,col="blue")
lines(hmu_curve_grassberger03,col="magenta")

