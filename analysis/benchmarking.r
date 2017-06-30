# bench marking CSSR
rm(list = ls())

require('CSSR')

output <- runCSSR(alphabet="01", data="010101101010101", maxLength=15, 
                  isChi=FALSE, sigLevel=0.001, outputPrefix="../output/test")
