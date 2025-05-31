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
  
  # Precompute which chromosomes differ from the best solution
  differing <- rowSums(chromosome != matrix(solution_best, nrow=popsize, ncol=length(solution_best), byrow=TRUE)) != 0
  
  # Random numbers for population-level mutation
  pop_rnd <- runif(popsize)
  mutate_pop <- (pop_rnd < pop_mutation_rate) & differing
  
  if (any(mutate_pop)) {
    # Only process chromosomes where mutation should occur
    idx <- which(mutate_pop)
    for (i in idx) {
      # Vector of random numbers for each gene
      gene_rnd <- runif(genomeLength)
      mutate_genes <- gene_rnd < mutation_rate
      keep_genes <- !mutate_genes
      
      # Swap where mutate_genes is TRUE, keep otherwise
      if (any(mutate_genes)) {
        work_q_alphabeta[mutate_genes, 1, i] <- q_alphabeta[mutate_genes, 2, i]
        work_q_alphabeta[mutate_genes, 2, i] <- q_alphabeta[mutate_genes, 1, i]
      }
      # These assignments are redundant unless q_alphabeta has changed elsewhere,so we can skip
      # work_q_alphabeta[keep_genes, 1, i] <- q_alphabeta[keep_genes, 1, i]
      # work_q_alphabeta[keep_genes, 2, i] <- q_alphabeta[keep_genes, 2, i]
    }
  }
  return(work_q_alphabeta)
}
