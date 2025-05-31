#---------------------------
# POPULATION INITIALIZATION                     
#---------------------------
generate_pop <- function(popsize, genomeLength, q_alphabeta, rot, theta, h, qubit_0) {
  # Pre-generate all random thetas in radians (uniformly between 0 and 2*pi)
  thetas <- runif(popsize * genomeLength, min = 0, max = 2 * pi)
  dim(thetas) <- c(genomeLength, popsize)
  
  # Precompute the qubit basis vector (length 2)
  qubit_0 <- as.numeric(qubit_0)
  
  # Precompute h
  h <- as.matrix(h)
  
  for (i in seq_len(popsize)) {
    for (j in seq_len(genomeLength)) {
      theta <- thetas[j, i]
      rot <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow = 2, byrow = TRUE)
      q_ab <- rot %*% h %*% qubit_0
      q_alphabeta[j, , i] <- as.numeric(q_ab)
    }
  }
  return(q_alphabeta)
}
