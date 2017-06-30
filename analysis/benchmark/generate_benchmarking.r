# Benchmarking
rm(list = ls())
## 1. generating idealized "memory" datasets ===========================

# Type I memory: AR(1) timeseries
# parameter value for p in AR(1) model
a_range <- c(0.05, 0.3, 0.5, 0.7, 0.9)
sd_range <- c(0.1, 0.5, 1)
len <- 4 # 10^len number of records

for(a_value in a_range){
  for(sd_value in sd_range){
    x <- arima.sim(list(order=c(1,0,0), sd = sd_value, ar=c(a_value)), n = 10^len)
    save(x, file = paste('./benchmark/TypeI_a', a_value, 'sd_value', sd_value,"benchmarking_data.rda", sep='_'))
  }
}

# Type II memory: 
N_range <- c(3, 5, 10, 15, 30)
sd_range <- c(0.1, 0.5, 1)
Iters <- 20
len <- 4 # 10^len number of records

rdirichlet <- function(a) {
  y <- rgamma(length(a), a, 1)
  return(y / sum(y))
}

for(N in N_range){
  for(sd_value in sd_range){
  wts <- rdirichlet(runif(n = N, min = 0, max = 10))
  x <- rnorm(n = N, mean = 0, sd = sd_value)
  for(i in 1:(10^len)){
    new <- x[1:3] * wts + rnorm(n = 1, mean = 0, sd = sd_value)
    x <- append(x, new, after = 0)
  }
  x <- rev(x)
  save(x, file = paste('./benchmark/TypeII_N', N_value, 'sd_value', sd_value, 
                       "benchmarking_data.rda", sep='_'))
  }
}  

## 2. pertubations =================
# + add iid white noise
# + add long (seasonal) cycles
# + add short (diurnal) cycles 
# + add long-term trend

# pert_type can be 'idd_noise', 'long_cycle', 'short_cycle', and 'trend'
add_pert <- function(a, pert_type = 'iid_noise') { 
  len_max <- length(a)
  switch(pert_type,
         iid_noise={b <- a + rnorm(n = len_max, mean = 0, sd = sd(a) + 0.2)},
         long_cycle={ # to mimic a long-term cycle, e.g. seasonal
           b <- a + (max(a) - min(a)) * sin((1:len_max)/50)
           }, 
         short_cycle={# to mimic a short-term cycle, e.g. diurnal
           b <- a + (max(a) - min(a)) * sin(1:len_max)
         },
         trend={# trend
           b <- a + seq(0, (max(a) - min(a)), length.out = len_max)
         },
         stop("unsurported pert_type")
  )
  return(b)
}

# example use:
x <- rnorm(1000)
y1 <- add_pert(x, pert_type = 'long_cycle')
plot(y1)
y2 <- add_pert(x, pert_type = 'short_cycle')
plot(y2)
y3 <- add_pert(x, pert_type = 'trend')
plot(y3)
