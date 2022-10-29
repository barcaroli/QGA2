#' Traveller Salesman Problem fitness evaluation function
#'
#' @description 
#' 
#' Fitness evaluation to be used for the optimization of a Traveller Salesman Problem
#' 
#' @details
#' 
#' This function is the one that performs the evaluation of the fitness in the case of 
#' the optimization of a Traveller Salesman Problem.
#' The function takes one of the solutions considered at the k-th iteration of the Quantum
#' Genetic Algorithm, and determines the minimum sample size required to
#' be compliant with precision constraints on the target variables.
#' 
#' Input required is the dataframe containing, for each city, their geographical coordinates.
#' 
#' @param solution the solution to be evaluated
#' @param eval_func_inputs specific inputs for best stratification 
#' (list with sampling frame and precision constraints) 
#' 
#' @export
#'  
#' 
TravellerSalesman <- function(solution,distance) {
  l = 0.0  
  for (i in 2:length(solution)) {
    l = l+distance[solution[i-1], solution[i]]
  }
  l + distance[solution[1],solution[length(solution)]]
  penal <- ((nrow(distance)) - length(table(solution)))*sum(distance)/10
  cost <- -(penal+l)
  return(cost)
}

