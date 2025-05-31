#---------------------------
# MEASUREMENT                     
#---------------------------
measure <- function(popsize,
                    genomeLength,
                    q_alphabeta,
                    chromosome) {
  for (i in (1:popsize)) {
    for (j in (1:genomeLength)) {
      p_alpha <- runif(1)
      if (p_alpha <= 2*q_alphabeta[j, 1, i]^2) chromosome[i, j] <- 0
      if (p_alpha > 2*q_alphabeta[j, 1, i]^2) chromosome[i, j] <- 1
    }
  }
  return(chromosome)
}
