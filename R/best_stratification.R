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
#' The function takes one of the solutions considered at the k-th iteration of the Quantum
#' Genetic Algorithm, and determines the minimum sample size required to
#' be compliant with precision constraints on the target variables.
#' 
#' Requires package SamplingStrata to execute functions aggrStrata and Bethel
#' 
#' @param solution the solution to be evaluated
#' @param eval_func_inputs specific inputs for best stratification 
#' (list with sampling frame and precision constraints) 
#' 
#' @export
#'  
#' 
best_stratification <- function(solution,
                                eval_func_inputs) {
  require(SamplingStrata)
  nstrat <- length(table(solution))
  frame <- eval_func_inputs[[1]]
  cv <- eval_func_inputs[[2]]
  strata = SamplingStrata::aggrStrata2(dataset=frame,
                         model=NULL,
                         vett=solution,
                         dominio=1)
  fitness <- -sum(SamplingStrata::bethel(strata, cv, realAllocation = TRUE))
  return(fitness)
}