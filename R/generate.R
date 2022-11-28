#' generate_pop
#'
#' @description 
#' 
#' Function to initialize first population
#' 
#' @param popsize the number of generated solutions (population) to be evaluated at each iteration
#' (default is 20)
#' @param genomeLength the length of the genome (or chromosome), representing a possible solution 
#' @param q_alphabeta the array containing the values of the amplitudes of the qubits
#' @param rot the rotation array
#' @param theta the rotation angle 
#' @param h the hadamard gate 
#' to each individual in the population (default is 1/(popsize+1))
#' @param qubit_0 the array of qubit 0
#' 
#' @export
#' 
#' @return A numeric vector giving the best solution obtained by the QGA

#---------------------------
# POPULATION INITIALIZATION                     
#---------------------------

generate_pop <- function(popsize,
                         genomeLength,
                         q_alphabeta,
                         rot,
                         theta,
                         h,
                         qubit_0) {
  for (i in c(1:popsize)) {
    for (j in c(1:genomeLength)) {
      theta <- runif(1) * 360
      theta <- pi*theta
      rot[1, 1] <- cos(theta)
      rot[1, 2] <- -sin(theta)
      rot[2, 1] <- sin(theta)
      rot[2, 2] <- cos(theta)
      q_alphabeta[j, 1, i] <- rot[1, 1] * h[1, 1] * qubit_0[1] + rot[1, 2] * h[1, 2] * qubit_0[2]
      q_alphabeta[j, 2, i] <- rot[2, 1] * h[2, 1] * qubit_0[1] + rot[2, 2] * h[2, 2] * qubit_0[2]
    }
  }
  return(q_alphabeta)
}