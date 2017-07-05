require(forecast)
require(pheatmap)
require(CSSR)

src_dir <- "../../src/"
source(paste0(src_dir, "symbolization.R", collapse=""))
source(paste0(src_dir, "entropy_calc.R", collapse=""))
source(paste0(src_dir, "run_cssr_by_min_bic.R", collapse=""))

setwd("/Users/szhan/Desktop/Learn/csss2017/memory/analysis/benchmark/")

## 1. Simulate stationary AR and MA time series
# General settings
steps <- 10^4 # Length of time series
ntimes <- 10 # Number of simulated time series per parameter set
sd_range <- c(0.1, 0.5, 1, 2) # Std dev for white noise
n_parts <- seq(from=2, to=6, by=2)
min_delta_bic <- function(n) { (n - 1) * log2(steps) }


# AR of order 1
a_range <- c(0.05, 0.3, 0.5, 0.7, 0.9)

ar_ts <- list()
for ( n in 1:ntimes ) {
  m <- matrix(nrow=length(a_range) * length(sd_range), ncol=steps)
  for ( i in 1:length(a_range) ) {
    for ( j in 1:length(sd_range) ) {
      k <- j + (i-1) * length(sd_range)
      m[k,] <- arima.sim(list(order=c(1,0,0), ar=a_range[i], sd=sd_range[j]), n=steps)
    }
  }
  ar_ts[[n]] <- m
}


# MA of order 5 to 30, by increments of 5
q_range <- seq(from=5, to=10, by=5)

rdirichlet <- function(a) {
  y <- rgamma(length(a), a, 1)
  return(y / sum(y))
}

ma_ts <- list()
for ( n in 1:ntimes ) {
  m <- matrix(nrow=length(q_range) * length(sd_range), ncol=steps)
  for ( i in 1:length(q_range) ) {
    for ( j in 1:length(sd_range) ) {
      ma_coef <- rdirichlet(runif(q_range[i], 0, 10))
      k <- j + (i-1) * length(sd_range)
      m[k,] <- arima.sim(list(order=c(0,0,q_range[i]), ma=ma_coef), sd=sd_range[j], n=steps)
    }
  }
  ma_ts[[n]] <- m
}


### Run CSSR (by min BIC) on time series data
run_cssr_ts <- function(sim_ts, par_range, part_strat) {
  nbr_sims <- length(sim_ts)
  nbr_psets <- nrow(sim_ts[[1]])
  nbr_strats <- length(part_strat)
  
  ent_rate <- matrix(rep(0, nbr_psets*nbr_strats), nrow=nbr_psets, ncol=nbr_strats)
  row_names <- c()
  for ( i in par_range ) {
    for ( j in sd_range ) {
      row_names <- append(row_names, paste0("par=", i, "_", "sd=", j, collapse=""))
    }
  }
  row.names(ent_rate) <- row_names
  colnames(ent_rate) <- paste0("b=", part_strat)
  
  rel_ent <- ent_rate
  rel_ent_rate <- ent_rate
  c_mu <- ent_rate
  variation <- ent_rate
  nbr_inf_states <- ent_rate
  
  for ( n in 1:nbr_sims ) {
    print(paste0("sim:", n))
    for ( i in 1:nbr_psets ) {
      for ( j in 1:nbr_strats ) {
        x <- sim_ts[[n]][i,]
        y <- paste0(parti(x, parti_type='maxent', n_bin=part_strat[j]), collapse="")
        s <- paste0(seq(0, part_strat[j]-1, 1), collapse="")
        r <- run_cssr_by_min_bic(alphabet=s, data=y, out_file="test", 
                                 threshold=min_delta_bic(part_strat[j]))
        
        ent_rate[i,j] <- ent_rate[i,j] + r$results$results$ent_rate
        rel_ent[i,j] <- rel_ent[i,j] + r$results$results$rel_ent
        rel_ent_rate[i,j] <- rel_ent_rate[i,j] + r$results$results$rel_ent_rate
        c_mu[i,j] <- c_mu[i,j] + r$results$results$c_mu
        variation[i,j] <- variation[i,j] + r$results$results$variation
        nbr_inf_states[i,j] <- nbr_inf_states[i,j] + r$results$results$nbr_inferred_states
      }
    }
  }
  
  # Take mean over 'steps' number of simulations
  ent_rate <- ent_rate / nbr_sims
  rel_ent <- rel_ent / nbr_sims
  rel_ent_rate <- rel_ent_rate / nbr_sims
  c_mu <- c_mu / nbr_sims
  variation <- variation / nbr_sims
  nbr_inf_states <- nbr_inf_states / nbr_sims
  
  res <- list(ent_rate=ent_rate,
              rel_ent=rel_ent,
              rel_ent_rate=rel_ent_rate,
              c_mu=c_mu,
              variation=variation,
              nbr_inf_states=nbr_inf_states)
  
  return(res)
}


cssr_ar <- run_cssr_ts(sim_ts=ar_ts, par_range=a_range, part_strat=n_parts)
cssr_ma <- run_cssr_ts(sim_ts=ma_ts, par_range=q_range, part_strat=n_parts)


### 2. Perturbations
add_pert <- function(a, pert_type='none') {
  len_max <- length(a)
  double_ampl <- max(a) - min(a)
  
  switch(pert_type,
         none={ b <- a },
         short_cycle={ b <- a + double_ampl * sin(1:len_max) },
         long_cycle={ b <- a + double_ampl * sin((1:len_max)/50) },
         up_trend={ b <- a + seq(0, double_ampl, length.out=len_max) },
         down_trend={ b <- a - seq(0, double_ampl, length.out=len_max) },
         stop("ERROR: Unrecognized pert_type")
  )
  
  return(b)
}

p_types <- c('none', 'long_cycle', 'short_cycle', 'up_trend')

