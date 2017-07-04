require(CSSR)

data_file <- "../../data/demo_two-states.txt"
out_file <- paste0(data_file, ".out", collapse="")

dat <- readLines(data_file)
alpha <- paste0(sort(unique(unlist(strsplit(dat, split="")))), collapse="")

# Get results of CSSR run under optimal max history length determined by minimum BIC
best_res <- run_cssr_by_min_bic(alpha, dat, out_file)

# Show how BIC changes over max history length
plot(y=best_res$scores, x=1:length(best_res$scores), ylab="BIC", xlab="Max history length")

# Show structure of the inferred HMM
dot(best_res$results$results$dot_graph)

