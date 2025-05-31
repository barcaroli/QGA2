#--------------
# ROTATION                   
#--------------
rotation <- function(chromosome, best_chromosome, generation, genomeLength,
                     solution_best, q_alphabeta, work_q_alphabeta, popsize,
                     fitness, theta) {

  apply_rotation <- function(q, s, theta) {
    cs <- cos(s * theta)
    sn <- sin(s * theta)
    r1 <- cs * q[1] - sn * q[2]
    r2 <- sn * q[1] + cs * q[2]
    c(r1, r2)
  }

  for (i in seq_len(popsize)) {
    if (any(chromosome[i, ] != solution_best)) {
      for (j in seq_len(genomeLength)) {
        xb <- chromosome[i, j]
        bb <- chromosome[best_chromosome[generation], j]
        qa <- q_alphabeta[j, , i]

        # Determine if fitness[i] is better than best
        better <- fitness[i] < fitness[best_chromosome[generation]]

        # Select action
        s <- NULL
        if (better) {
          if (xb == 0 && bb == 1) {
            if (qa[1] * qa[2] >= 0) s <- 1
            else if (qa[1] * qa[2] < 0) s <- -1
            else if (qa[1] == 0) s <- 0
            else if (qa[2] == 0) s <- ifelse(runif(1) < 0.5, 1, -1)
          } else if (xb == 1 && bb == 0) {
            if (qa[1] * qa[2] >= 0) s <- -1
            else if (qa[1] * qa[2] < 0) s <- 1
            else if (qa[1] == 0) s <- ifelse(runif(1) < 0.5, 1, -1)
            else if (qa[2] == 0) s <- 0
          }
        } else { # fitness[i] >= best
          if (xb == 0 && bb == 1) {
            if (qa[1] * qa[2] >= 0) s <- -1
            else if (qa[1] * qa[2] < 0) s <- 1
            else if (qa[1] == 0) s <- ifelse(runif(1) < 0.5, 1, -1)
            else if (qa[2] == 0) s <- 0
          } else if (xb == 0 && bb == 1) {
            if (qa[1] * qa[2] >= 0) s <- 1
            else if (qa[1] * qa[2] < 0) s <- -1
            else if (qa[1] == 0) s <- 0
            else if (qa[2] == 0) s <- ifelse(runif(1) < 0.5, 1, -1)
          } else if (xb == 1 && bb == 1) {
            if (qa[1] * qa[2] >= 0) s <- 1
            else if (qa[1] * qa[2] < 0) s <- -1
            else if (qa[1] == 0) s <- 0
            else if (qa[2] == 0) s <- ifelse(runif(1) < 0.5, 1, -1)
          }
        }

        if (!is.null(s)) {
          rotated <- apply_rotation(qa, s, theta)
          q_alphabeta[j, 1, i] <- rotated[1]
          q_alphabeta[j, 2, i] <- rotated[2]
        }
      }
    }
  }
  q_alphabeta
}
