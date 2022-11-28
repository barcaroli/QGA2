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