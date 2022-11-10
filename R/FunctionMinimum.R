#' Function Minimum fitness evaluation function
#'
#' @description 
#' 
#' Fitness evaluation to be used for finding the minimum of a function
#' 
#' @details
#' 
#' This function is the one that performs the evaluation of the fitness in the case of 
#' searching the minimum of a function.
#' The function takes one of the solutions considered at the k-th iteration of the Quantum
#' Genetic Algorithm, and 
#' 
#' 1. convert the vector of binary values representing the current solution (value of the x);
#' 2. calculates 
#' 
#' @param solution the solution to be evaluated
#' @param eval_func_inputs specific inputs for the evaluation function
#' (n: the number of digits with required by the integer numbers of the x variable,
#'  y: the vector of values of the y variable) 
#'  
#' @value the value of the y variable
#' 
#' @export
#'  
#' 
FunctionMinimum <- function(solution,inputs) {
  n <- inputs[[1]] 
  y <- inputs[[2]]
  valuex <- 0
  for (w in c(1:n)) {
    valuex <- valuex + solution[w]*2^(n-w) 
  }
  valuex <- valuex + 1
  valuey <- y[valuex]
  return(-valuey)
}
