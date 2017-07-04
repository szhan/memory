# author: yliu11

# Partition data points as per one of two strategies:
# 1. Equiprobable bins
# 2. Extrema (lower 2.5%tile or upper 97.5%tile)
parti <- function(a, parti_type='maxent', n_bin=3) {
  switch(parti_type,
         maxent={
           thres <- quantile(a, probs=seq(0,1,length.out=n_bin + 1))
           b <- as.numeric(cut(a, breaks=thres, include.lowest=T)) - 1
         },
         extrema={
           thres <- quantile(a, probs=c(0, 0.025, 0.975, 1))
           b <- as.numeric(cut(a, breaks=thres, include.lowest=T)) - 1
         },
         stop("ERROR: Unrecognized parti_type")
  )
  return(b)
}

