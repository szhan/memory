require(plyr)


all_substr <- function(x, n){
	substring(x, 1:(nchar(x) - n + 1), n:nchar(x))
}


compute_shannon_entropy <- function(p) {
	-sum(p * log2(p))
}


get_HL <- function(L, seq) {
	a <- table(all_substr(seq, L))
	return(compute_shannon_entropy(a / sum(a)))
}


get_entropy_curve <- function(seq, maxL) {
	HL <- lapply(1:maxL, get_HL, seq=seq)
	hmu <- diff(unlist(HL))
	return(hmu)
}


get_hmu <- function(i, w, seq, plotit=FALSE){
	maxL <- 20
	seq <- substr(seq, i, i+w-1)
	
	hmu <- get_entropy_curve(seq, maxL)
	max_ind <- which(abs(diff(hmu)) %in% max(abs(diff(hmu))))
	min_ind <- which(abs(diff(hmu[1:max_ind])) %in% min(abs(diff(hmu[1:max_ind]))))

	if ( plotit == TRUE ) {
		par(mfrow=c(2,1))
		plot(1:(min_ind), hmu[1:min_ind])
		plot(hmu)
		par(mfrow=c(1,1))
	}

	return(hmu[min_ind])
}


var_hmu <- function(sequence, w){
	print(nchar(sequence) - w + 1)
	lapply(1:(nchar(sequence) - w + 1), get_hmu, w=w, seq=sequence)
  	hmu_vec <- unlist(lapply(1:(nchar(sequence)-w+1), get_hmu, w=w, seq=sequence))
	return(var(hmu_vec))
}


