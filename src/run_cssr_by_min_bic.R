# author: szhan

# Calculate the BIC score of an inferred minimum HMM,
#   given state probability distributions and transition matrices
compute_lnlik_cssr <- function(res, alphabet, data, exponent=5) {
  alpha_seq <- unlist(strsplit(alphabet, split=""))
  data_seq <- unlist(strsplit(data, split=""))

  state_prob <- res$matrices$state_prob
  trans_prob <- res$matrices$trans_prob
  trans_stat <- res$matrices$trans_stat
  nbr_observed_states <- res$results$alpha_size
  nbr_inferred_states <- res$results$nbr_inferred_states
  data_size <- res$info$data_size

  # Unpack matrices into observed state-specific transition probability matrices
  state_mat <- list()
  for ( i in 1:nbr_observed_states ) {
    tmp_mat <- matrix(data=rep(0, nbr_inferred_states^2),
                      nrow=nbr_inferred_states, 
                      ncol=nbr_inferred_states)
    for ( j in 1:nbr_inferred_states ) {
      tmp_mat[j, trans_stat[j,i] + 1] <- trans_prob[j,i]
    }
    state_mat[[i]] <- tmp_mat
  }

  # Compute log likelihood of data given inferred HMM
  base_mul_factor <- 10^exponent
  min_threshold <- 10^-exponent
  log_mul_factor <- 0
  
  lik <- state_prob
  for ( n in data_seq ) {
    idx <- grep(n, alpha_seq)
    lik <- lik %*% state_mat[[idx]]
    if ( sum(lik) < min_threshold ) {
      lik <- lik * base_mul_factor
      log_mul_factor <- log_mul_factor + log(base_mul_factor)
    }
  }
  
  lnlik <- log(sum(lik)) - log_mul_factor
  bic <- -2*lnlik + nbr_inferred_states * (nbr_observed_states-1) * log(data_size)
  
  return(bic)
}


# Run CSSR, sweeping over a range of max history length (1 to upper-bounded max history length),
#   retaining the results from the run with the minimum BIC
run_cssr_by_min_bic <- function(alphabet, data, out_file, threshold=0) {
  # Index corresponds to max history length
  bic_scores <- c()

  # Run CSSR on max history length of 1 to determine max history length
  ini_run <- runCSSR(alphabet=alpha, data=dat, maxLength=1, 
                  isChi=FALSE, sigLevel=0.001, outputPrefix=out_file)
  ini_bic <- compute_lnlik_cssr(res=ini_run, alphabet=alpha, data=dat)
  bic_scores <- ini_bic
  max_hist_len <- floor(log2(ini_run$info$data_size) / ini_run$results$ent_rate) - 1

  # Run CSSR over a range of max history length
  for ( i in 2:max_hist_len ) {
    run <- runCSSR(alphabet=alpha, data=dat, maxLength=i, 
                   isChi=FALSE, sigLevel=0.001, outputPrefix=out_file)
    bic <- compute_lnlik_cssr(res=run, alphabet=alpha, data=dat)
    bic_scores <- c(bic_scores, bic)
  }
  
  # Identify minimum max history length within threshold of minimum BIC across allowable max history lengths
  best_hist_len <- min(which(bic_scores <= min(bic_scores)+threshold))
  best_run <- runCSSR(alphabet=alpha, data=dat, maxLength=best_hist_len, 
                      isChi=FALSE, sigLevel=0.001, outputPrefix=out_file)
  

  return(list(results=best_run, scores=bic_scores, optim_hist_len=best_hist_len))
}

