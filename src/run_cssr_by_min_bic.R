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
  min_threshold <- 10^-exponent
  base_factor <- 10^exponent
  log2_base_factor <- log2(base_factor)
  log2_mul_factor <- log2(1)

  lik <- state_prob
  for ( n in data_seq ) {
    idx <- grep(n, alpha_seq)
    lik <- lik %*% state_mat[[idx]]
    if ( sum(lik) < min_threshold ) {
      lik <- lik * base_factor
      log2_mul_factor <- log2_mul_factor + log2_base_factor
    }
  }
  
  log2lik <- log2(sum(lik)) - log2_mul_factor
  bic <- -2*log2lik + nbr_inferred_states * (nbr_observed_states-1) * log2(data_size)
  
  return(bic)
}


# Run CSSR, sweeping over a range of max history length (1 to upper-bounded max history length),
#   retaining the results from the run with the minimum BIC
run_cssr_by_min_bic <- function(alphabet, data, out_file, threshold=0) {
  # Run CSSR on max history length of 1 to determine max history length
  ini_run <- runCSSR(alphabet=alphabet, data=data, maxLength=1, 
                     isChi=FALSE, sigLevel=0.001, outputPrefix=out_file)
  ini_bic <- compute_lnlik_cssr(res=ini_run, alphabet=alphabet, data=data)
  max_hist_len <- floor(log2(ini_run$info$data_size) / ini_run$results$ent_rate) - 1
  
  # Run CSSR over a range of max history length
  best_run <- ini_run
  best_bic <- ini_bic
  prev_bic <- ini_bic

  if ( max_hist_len >= 2 ) {
    for ( i in 2:max_hist_len ) {
      curr_run <- runCSSR(alphabet=alphabet, data=data, maxLength=i,
                          isChi=FALSE, sigLevel=0.001, outputPrefix=out_file)
      curr_bic <- compute_lnlik_cssr(res=curr_run, alphabet=alphabet, data=data)

      if ( (prev_bic - curr_bic) >= threshold ) {
        best_run <- curr_run
        best_bic <- curr_bic
        prev_bic <- curr_bic
      } else {
        break
      }
    }
  }

  return(list(results=best_run, bic=best_bic))
}

