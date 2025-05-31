#----------
# REPAIR                   
#----------
repair <- function(popsize,
                   chromosome,
                   geneLength,
                   genomeLength,
                   nvalues_sol,
                   Genome) {
  diff <- 2^geneLength - nvalues_sol
  acceptable_values <- seq_len(nvalues_sol)
  a <- matrix(seq_len(genomeLength), nrow = geneLength, ncol = Genome)
  
  for (i in seq_len(popsize)) {
    solution1 <- matrix(chromosome[i, ], nrow = geneLength, ncol = Genome)
    # Convert binary to integer (vectorized)
    powers <- 2^((geneLength - 1):0)
    solution <- as.vector(colSums(solution1 * powers)) + 1

    # Repair out-of-range solutions
    out_of_range <- which(!solution %in% acceptable_values)
    if (length(out_of_range) > 0) {
      solution[out_of_range] <- solution[out_of_range] - diff
    }

    # Write repaired solution back
    for (x in seq_len(Genome)) {
      y1 <- a[1, x]
      y2 <- a[geneLength, x]
      chromosome[i, y1:y2] <- as.binary(solution[x] - 1, n = geneLength)
    }
  }
  return(chromosome)
}
