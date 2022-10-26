#' Best stratification fitness evaluation function
#'
#' @description 
#' 
#' Fitness evaluation to be used for the optimization of a sampling frame stratification
#' 
#' @details
#' 
#' This function is the one that perform the evaluation of the fitness in the case of 
#' the optimization of a sampling frame stratification.
#' The function takes the set of solution considered at the k-th iteration of the Quantum
#' Genetic Algorithm, and for each of them determines the minimum sample size required to
#' be compliant with precision constraints on the target variables.
#' 
#' @param chromosome the current set of the solutions to be evaluated
#' @param frame the sampling frame
#' @param cv the set of precision constraints
#' 
#' @export
#'  
#' 
best_stratification <- function(chromosome,
                                frame,
                                cv) {
  fitness_total <- 0
  sum_sqr <- 0
  fitness_average <- -999999
  variance <- 0
  if (nstrat < 2^Genome_el) chromosome <- repair(chromosome)
  for (i in c(1:popsize)) {
    solution1 <- array(chromosome[i,],c(Genome_el,Genome))
    solution <- c(rep(0,Genome))
    for (x in c(1:Genome)) {
      for (y in c(1:Genome_el)) {
        solution[x] <- solution[x] + solution1[y,x]*2^(Genome_el - y) 
      }
    }
    solution <- solution + 1
    strata = aggrStrata2(dataset=frame,
                         model=NULL,
                         vett=solution,
                         dominio=1)
    fitness[i] <- -sum(bethel(strata, cv, realAllocation = TRUE))
    fitness_total <- fitness_total + fitness[i]
  }
  fitness_max <- -999999
  fitness_average <- fitness_total / popsize
  the_best_chromosome <- 0
  for (i in c(1:popsize)) {
    if (fitness[i] >= fitness_max) {
      fitness_max <- fitness[i]
      the_best_chromosome <- i
      solution_max <- chromosome[i, ]
    }
    best_chromosome[generation] <- the_best_chromosome
  }
  return(list(
    fitness = fitness,
    fitness_max = fitness_max,
    fitness_average = fitness_average,
    best_chromosome = best_chromosome,
    solution_max = solution_max
  ))
}