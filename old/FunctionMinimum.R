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
#' The function is trivial: takes one of the solutions considered at the k-th iteration 
#' (the value of the X), and reports the value of the y.
#'  
#' @param solution the solution to be evaluated
#' @param eval_func_inputs specific input for the evaluation function
#' (y: the vector of values of the y variable) 
#'  
#' @return the value of the y variable
#' 
#' @export
#'  
#' 
FunctionMinimum <- function(solution,y) {
  value <- y[solution]
  return(-value)
}
