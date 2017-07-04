require(forecast)
require(pheatmap)
require(CSSR)

src_dir <- "../../src/"
source(paste0(src_dir, "symbolization.R", collapse=""))
source(paste0(src_dir, "entropy_calc.R", collapse=""))
source(paste0(src_dir, "run_cssr_by_min_bic.R", collapse=""))


## 1. Simulate stationary AR and MA time series
# AR of order 1
a_range <- c(0.05, 0.3, 0.5, 0.7, 0.9)
sd_range <- c(0.1, 0.5, 1, 5)
steps <- 10^4
ntimes <- 100

ar_ts <- list()
for ( n in 1:ntimes ) {
  m <- matrix(nrow=length(a_range) * length(sd_range), ncol=steps)
  for ( i in 1:length(a_range) ) {
    for ( j in 1:length(sd_range) ) {
      k <- j + (i - 1)*length(sd_range)
      m[k,] <- arima.sim(list(order=c(1,0,0), ar=a_range[i], sd=sd_range[j]), n=steps)
      }
  }
  ar_ts[[n]] <- m
}


# MA of order 5 to 30, by increments of 5
q_range <- seq(from=5, to=30, by=5)
sd_range <- c(0.1, 0.5, 1, 5)
steps <- 10^4

rdirichlet <- function(a) {
  y <- rgamma(length(a), a, 1)
  return(y / sum(y))
}


ma_ts <- matrix(nrow=length(q_range) * length(sd_range), ncol=steps)

m <- 1
for ( i in 1:length(q_range) ) {
  for ( j in 1:length(sd_range) ) {
    ma_coef <- rdirichlet(runif(q_range[i],0,10))
    x <- arima.sim(list(order=c(0,0,q_range[i]), ma=ma_coef), sd=sd_range[j], n=steps)
    ma_ts[m,] <- x
    m <- m + 1
  }
}


### Print AR(p=1) reults in heatmaps
ar_ent_rate <- matrix(nrow=nrow(ar_ts), ncol=length(n_parts))
row_names <- c()
for ( i in a_range ) {
  for ( j in sd_range ) {
    row_names <- append(row_names, paste0("a=", i, "_", "sd=", j, collapse=""))
  }
}
row.names(ar_ent_rate) <- row_names
colnames(ar_ent_rate) <- paste0("b=", n_parts)

ar_rel_ent <- ar_ent_rate
ar_rel_ent_rate <- ar_ent_rate
ar_c_mu <- ar_ent_rate
ar_variation <- ar_ent_rate
ar_nbr_inf_states <- ar_ent_rate



x <- ar_ts[i,]
y <- paste0(parti(x, parti_type='maxent', n_bin=n_parts[j]), collapse="")
s <- paste0(seq(0, n_parts[j] - 1, 1), collapse="")
r1 <- runCSSR(alphabet=s, data=y, maxLength=1, isChi=FALSE, sigLevel=0.001, outputPrefix="test")
# L <= log2(size of data) / h - 1
# h is bounded by log2(entropy rate)
L <- trunc(log2(length(x)) / r1$ent_rate - 1)
r2 <- runCSSR(alphabet=s, data=y, maxLength=L, isChi=FALSE, sigLevel=0.001, outputPrefix="test")


n_parts <- seq(from=2, to=10, by=2)
for ( i in 1:nrow(ar_ts) ) {
  for ( j in 1:length(n_parts) ) {
    
    ar_ent_rate[i,j] <- r2$ent_rate
    ar_rel_ent[i,j] <- r2$rel_ent
    ar_rel_ent_rate[i,j] <- r2$rel_ent_rate
    ar_c_mu[i,j] <- r2$c_mu
    ar_variation[i,j] <- r2$variation
    ar_nbr_inf_states[i,j] <- r2$nbr_inferred_states
  }
}


### Print MA results in heatmaps
ma_ent_rate <- matrix(nrow=nrow(ma_ts), ncol=length(n_parts))
row_names <- c()
for ( i in q_range ) {
  for ( j in sd_range ) {
    row_names <- append(row_names, paste0("q=", i, "_", "sd=", j, collapse=""))
  }
}
row.names(ma_ent_rate) <- row_names
colnames(ma_ent_rate) <- paste0("b=", n_parts)

ma_rel_ent <- ma_ent_rate
ma_rel_ent_rate <- ma_ent_rate
ma_c_mu <- ma_ent_rate
ma_variation <- ma_ent_rate
ma_nbr_inf_states <- ma_ent_rate

n_parts <- seq(from=2, to=10, by=2)
for ( i in 1:nrow(ma_ts) ) {
  for ( j in 1:length(n_parts) ) {
    x <- ma_ts[i,]
    y <- paste0(parti(x, parti_type='maxent', n_bin=n_parts[j]), collapse="")
    s <- paste0(seq(0, n_parts[j] - 1, 1), collapse="")
    r1 <- runCSSR(alphabet=s, data=y, maxLength=1, isChi=FALSE, sigLevel=0.001, outputPrefix="test")
    # L <= log2(size of data) / h - 1
    # h is bounded by log2(entropy rate)
    L <- trunc(log2(length(x)) / r1$ent_rate - 1)
    r2 <- runCSSR(alphabet=s, data=y, maxLength=L, isChi=FALSE, sigLevel=0.001, outputPrefix="test")
    
    ma_ent_rate[i,j] <- r2$ent_rate
    ma_rel_ent[i,j] <- r2$rel_ent
    ma_rel_ent_rate[i,j] <- r2$rel_ent_rate
    ma_c_mu[i,j] <- r2$c_mu
    ma_variation[i,j] <- r2$variation
    ma_nbr_inf_states[i,j] <- r2$nbr_inferred_states
  }
}


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


par(mfrow=c(2,2))
x <- ar_ts[1,]
x_range <- c(0,1000)

y1 <- add_pert(x, pert_type = 'long_cycle')
plot(y1, xlim=x_range, pch=16)
y2 <- add_pert(x, pert_type = 'short_cycle')
plot(y2, xlim=x_range, pch=16)
y3 <- add_pert(x, pert_type = 'up_trend')
plot(y3, xlim=x_range, pch=16)
y4 <- add_pert(x, pert_type = 'down_trend')
plot(y4, xlim=x_range, pch=16)


ts_oi <- ar_ts[1,]
p_types <- c('none', 'long_cycle', 'short_cycle', 'up_trend')
n_parts <- seq(from=2, to=10, by=2)

mat_ent_rate <- matrix(nrow=length(p_types), ncol=length(n_parts))
row.names(mat_ent_rate) <- p_types
colnames(mat_ent_rate) <- paste0("b=", n_parts)

mat_rel_ent <- mat_ent_rate
mat_rel_ent_rate <- mat_ent_rate
mat_c_mu <- mat_ent_rate
mat_variation <- mat_ent_rate
mat_nbr_inf_states <- mat_ent_rate

for ( i in 1:length(p_types) ) {
  for ( j in 1:length(n_parts) ) {
    x <- add_pert(ts_oi, p_types[i])
    y <- paste0(parti(x, parti_type='maxent', n_bin=n_parts[j]), collapse="")
    s <- paste0(seq(0, n_parts[j] - 1, 1), collapse="")
    r1 <- runCSSR(alphabet=s, data=y, maxLength=1, isChi=FALSE, sigLevel=0.001, outputPrefix="test")
    # L <= log2(size of data) / h - 1
    # h is bounded by log2(entropy rate)
    L <- trunc(log2(length(x)) / r1$ent_rate - 1)
    r2 <- runCSSR(alphabet=s, data=y, maxLength=L, isChi=FALSE, sigLevel=0.001, outputPrefix="test")
    
    mat_ent_rate[i,j] <- r2$ent_rate
    mat_rel_ent[i,j] <- r2$rel_ent
    mat_rel_ent_rate[i,j] <- r2$rel_ent_rate
    mat_c_mu[i,j] <- r2$c_mu
    mat_variation[i,j] <- r2$variation
    mat_nbr_inf_states[i,j] <- r2$nbr_inferred_states
  }
}

