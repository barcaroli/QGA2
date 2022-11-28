#----------
# MUTATION                   
#----------

mutation <- function(pop_mutation_rate, 
                     mutation_rate,
                     popsize,
                     chromosome,
                     solution_best,
                     q_alphabeta,
                     work_q_alphabeta,
                     genomeLength) {
  work_q_alphabeta <- q_alphabeta
  for (i in c(1:popsize)) {
    if (sum(chromosome[i, ] != solution_best) != 0) {
      rnd1 <- runif(1)
      if (rnd1 < pop_mutation_rate) {
        for (j in c(1:genomeLength)) {
          rnd2 <- runif(1)
          if (rnd2 < mutation_rate) {
            work_q_alphabeta[j, 1, i] <- q_alphabeta[j, 2, i]
            work_q_alphabeta[j, 2, i] <- q_alphabeta[j, 1, i]
          }
          if (rnd2 >= mutation_rate) {
            work_q_alphabeta[j, 1, i] <- q_alphabeta[j, 1, i]
            work_q_alphabeta[j, 2, i] <- q_alphabeta[j, 2, i]
          }
        }
      }
    }
  }
  q_alphabeta <- work_q_alphabeta
  return(q_alphabeta)
}
