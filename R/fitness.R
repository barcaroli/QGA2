#---------------------------
# EVALUATION OF THE SOLUTION                   
#---------------------------
evaluate <- function(chromosome,
                     best_chromosome,
                     popsize,
                     Genome,
                     geneLength,
                     nvalues_sol,
                     generation,
                     eval_fitness,
                     eval_func_inputs) {
  # Pre-allocate fitness vector
  fitness <- numeric(popsize)
  fitness_total <- 0
  fitness_max <- -Inf
  the_best_chromosome <- 0
  solution_max <- NULL  # Ensure it's defined even if not used

  # Precompute binary weights for efficiency
  bin_weights <- 2^((geneLength - 1):0)
  
  # Vectorized evaluation of all solutions
  for (i in seq_len(popsize)) {
    # Reshape gene row into matrix
    solution1 <- matrix(chromosome[i, ], nrow = geneLength, ncol = Genome)
    # Vectorized binary to integer conversion and add 1
    solution <- colSums(solution1 * bin_weights) + 1
    # Evaluate fitness
    fitness[i] <- eval_fitness(solution, eval_func_inputs)
    fitness_total <- fitness_total + fitness[i]
    if (fitness[i] >= fitness_max) {
      fitness_max <- fitness[i]
      the_best_chromosome <- i
      solution_max <- chromosome[i, ]
    }
  }

  fitness_average <- fitness_total / popsize
  best_chromosome[generation] <- the_best_chromosome
  return(list(
    fitness = fitness,
    fitness_max = fitness_max,
    fitness_average = fitness_average,
    best_chromosome = best_chromosome,
    solution_max = solution_max
  ))
}
