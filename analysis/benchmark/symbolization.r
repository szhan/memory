# partitioning strategies 

# the function for symblization, but the name "parti" is more attractive
parti <- function(a, parti_type = 'maxent_bin', n_bin = 3) { 
  len_max <- length(a)
  switch(parti_type,
         maxent_bin={
           thresh <- quantile(a, probs = seq(0,1,length.out = n_bin + 1))
           b <- as.numeric(cut(a, breaks = thresh, include.lowest = T)) - 1
           },
         extrema={
           n_bin = 3
           thresh <- quantile(a, probs = c(0, 0.025, 0.975, 1))
           b <- as.numeric(cut(a, breaks = thresh, include.lowest = T)) - 1
         }, 
         stop("unsurported parti_type")
  )
  return(b)
}

# example use:
x <- rnorm(1000)
y <- parti(x, n_bin = 5) # the default is the "maxent_bin" option
plot(y)
y <- parti(x, parti_type = 'extrema')
plot(y)
