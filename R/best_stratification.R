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
#' Requires package SamplingStrata to execute functions aggrStrata and Bethel
#' 
#' @param chromosome the current set of the solutions to be evaluated
#' @param best_chromosome the vector of best chromosomes found so far
#' @param popsize the number of generated solutions (population) to be evaluated at each iteration
#' (default is 20)
#' @param Genome the length of the genome, i.e. the elements of the chromosome representing an 
#' individual in the population 
#' @param Genome_el the length of the single elements in the elements, where each position 
#' represents a qubit
#' @param nvalues_sol the number of possible values contained in each element of the solution 
#' (it must be nvalues_sol <= 2^Genome_el) 
#' @param generation index of the current generation
#' @param eval_func_inputs specific inputs for best stratification 
#' (list with sampling frame and precision constraints) 
#' 
#' @export
#'  
#' 
best_stratification <- function(chromosome,
                                best_chromosome,
                                popsize,
                                Genome,
                                Genome_el,
                                nvalues_sol,
                                generation,
                                eval_func_inputs) {
  require(SamplingStrata)
  nstrat <- nvalues_sol
  frame <- eval_func_inputs[[1]]
  cv <- eval_func_inputs[[2]]
  fitness <- array(0.0, c(1, popsize))
  fitness_total <- 0
  sum_sqr <- 0
  fitness_average <- -99999999
  variance <- 0
  # if (nstrat < 2^Genome_el) chromosome <- repair(chromosome)
  for (i in c(1:popsize)) {
    solution1 <- array(chromosome[i,],c(Genome_el,Genome))
    solution <- c(rep(0,Genome))
    for (x in c(1:Genome)) {
      for (y in c(1:Genome_el)) {
        solution[x] <- solution[x] + solution1[y,x]*2^(Genome_el - y) 
      }
    }
    solution <- solution + 1
    strata = SamplingStrata::aggrStrata2(dataset=frame,
                         model=NULL,
                         vett=solution,
                         dominio=1)
    fitness[i] <- -sum(SamplingStrata::bethel(strata, cv, realAllocation = TRUE))
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