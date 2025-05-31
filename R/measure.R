#---------------------------
# MEASUREMENT                     
#---------------------------
measure <- function(popsize, 
                    genomeLength, 
                    q_alphabeta, 
                    chromosome) {
  # Generate all random numbers at once
  p_alpha <- matrix(runif(popsize * genomeLength), nrow = popsize, ncol = genomeLength)
  # Compute the threshold matrix
  thresh <- 2 * (q_alphabeta[, 1, ])^2
  thresh <- t(thresh) # Adjust to [i, j] layout as in chromosome
  # Assign 0 or 1 based on comparison
  chromosome[p_alpha <= thresh] <- 0
  chromosome[p_alpha > thresh] <- 1
  return(chromosome)
}
