#' Knapsack Problem fitness evaluation function
#'
#' @description 
#' 
#' Fitness evaluation to be used for the optimization of the knapsack problem
#' 
#' @details
#' 
#' This function is the one that performs the evaluation of the fitness in the case of 
#' the optimization of a Knapsack Problem.
#' The function takes one of the solutions considered at the k-th iteration of the Quantum
#' Genetic Algorithm, and 
#' 
#' 1. counts the number of items in knapsack;
#' 2. verifies if the total weights of the items in the knapsack is below the maximum allowed: 
#' if not, applies a penalty.
#' 
#' @param solution the solution to be evaluated
#' @param eval_func_inputs specific inputs for knapsack problem
#' (list with list of items with corresponding weights and allowed total weight) 
#' 
#' @export
#'  
#' 
KnapsackProblem <- function(solution,
                            eval_func_inputs) {
  solution <- solution - 1
  items <- eval_func_inputs[[1]]
  maxweight <- eval_func_inputs[[2]]
  tot_items <- sum(solution)
  # Penalization
  if (sum(items$weight[solution]) > maxweight) {
    tot_items <- tot_items - (sum(items$weight[solution]) - maxweight)  
  }
  return(tot_items)
}