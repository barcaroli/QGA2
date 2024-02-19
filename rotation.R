#--------------
# ROTATION                   
#--------------

rotation <- function(chromosome,
                     best_chromosome,
                     generation,
                     genomeLength,
                     solution_best,
                     q_alphabeta,
                     work_q_alphabeta,
                     popsize,
                     fitness, 
                     theta) {
  rot <- array(0, c(2, 2))
  for (i in c(1:popsize)) {
    if (sum(chromosome[i, ] != solution_best) != 0) {
      for (j in c(1:genomeLength)) {
        # Han-Kim lookup table
        # f(x) > f(b) FALSE
        if (fitness[i] < fitness[best_chromosome[generation]]) {
          # x = 0 b = 1
          if (chromosome[i, j] == 0 & chromosome[best_chromosome[generation], j] == 1) {
            if (q_alphabeta[j, 1, i]*q_alphabeta[j, 2, i] >= 0) {
              rot[1, 1] <- cos(theta)
              rot[1, 2] <- -sin(theta)
              rot[2, 1] <- sin(theta)
              rot[2, 2] <- cos(theta)
              work_q_alphabeta[j, 1, i] <- (rot[1, 1] * q_alphabeta[j, 1, i]) + (rot[1, 2] * q_alphabeta[j, 2, i])
              work_q_alphabeta[j, 2, i] <- (rot[2, 1] * q_alphabeta[j, 1, i]) + (rot[2, 2] * q_alphabeta[j, 2, i])
              q_alphabeta[j, 1, i] <- work_q_alphabeta[j, 1, i]
              q_alphabeta[j, 2, i] <- work_q_alphabeta[j, 2, i]
            }
            if (q_alphabeta[j, 1, i]*q_alphabeta[j, 2, i] < 0) {
              rot[1, 1] <- cos(-theta)
              rot[1, 2] <- -sin(-theta)
              rot[2, 1] <- sin(-theta)
              rot[2, 2] <- cos(-theta)
              work_q_alphabeta[j, 1, i] <- (rot[1, 1] * q_alphabeta[j, 1, i]) + (rot[1, 2] * q_alphabeta[j, 2, i])
              work_q_alphabeta[j, 2, i] <- (rot[2, 1] * q_alphabeta[j, 1, i]) + (rot[2, 2] * q_alphabeta[j, 2, i])
              q_alphabeta[j, 1, i] <- work_q_alphabeta[j, 1, i]
              q_alphabeta[j, 2, i] <- work_q_alphabeta[j, 2, i]
            }
          }
          # x = 1 b = 0
          if (chromosome[i, j] == 1 & chromosome[best_chromosome[generation], j] == 0) {
            if (q_alphabeta[j, 1, i]*q_alphabeta[j, 2, i] >= 0) {
              rot[1, 1] <- cos(-theta)
              rot[1, 2] <- -sin(-theta)
              rot[2, 1] <- sin(-theta)
              rot[2, 2] <- cos(-theta)
              work_q_alphabeta[j, 1, i] <- (rot[1, 1] * q_alphabeta[j, 1, i]) + (rot[1, 2] * q_alphabeta[j, 2, i])
              work_q_alphabeta[j, 2, i] <- (rot[2, 1] * q_alphabeta[j, 1, i]) + (rot[2, 2] * q_alphabeta[j, 2, i])
              q_alphabeta[j, 1, i] <- work_q_alphabeta[j, 1, i]
              q_alphabeta[j, 2, i] <- work_q_alphabeta[j, 2, i]
            }
            if (q_alphabeta[j, 1, i]*q_alphabeta[j, 2, i] < 0) {
              rot[1, 1] <- cos(theta)
              rot[1, 2] <- -sin(theta)
              rot[2, 1] <- sin(theta)
              rot[2, 2] <- cos(theta)
              work_q_alphabeta[j, 1, i] <- (rot[1, 1] * q_alphabeta[j, 1, i]) + (rot[1, 2] * q_alphabeta[j, 2, i])
              work_q_alphabeta[j, 2, i] <- (rot[2, 1] * q_alphabeta[j, 1, i]) + (rot[2, 2] * q_alphabeta[j, 2, i])
              q_alphabeta[j, 1, i] <- work_q_alphabeta[j, 1, i]
              q_alphabeta[j, 2, i] <- work_q_alphabeta[j, 2, i]
            }
          }
        }
        # f(x) > f(b) TRUE
        if (fitness[i] >= fitness[best_chromosome[generation]]) {
          # x = 0 b = 1
          if (chromosome[i, j] == 0 & chromosome[best_chromosome[generation], j] == 1) {
            if (q_alphabeta[j, 1, i]*q_alphabeta[j, 2, i] >= 0) {
              rot[1, 1] <- cos(-theta)
              rot[1, 2] <- -sin(-theta)
              rot[2, 1] <- sin(-theta)
              rot[2, 2] <- cos(-theta)
              work_q_alphabeta[j, 1, i] <- (rot[1, 1] * q_alphabeta[j, 1, i]) + (rot[1, 2] * q_alphabeta[j, 2, i])
              work_q_alphabeta[j, 2, i] <- (rot[2, 1] * q_alphabeta[j, 1, i]) + (rot[2, 2] * q_alphabeta[j, 2, i])
              q_alphabeta[j, 1, i] <- work_q_alphabeta[j, 1, i]
              q_alphabeta[j, 2, i] <- work_q_alphabeta[j, 2, i]
            }
            if (q_alphabeta[j, 1, i]*q_alphabeta[j, 2, i] < 0) {
              rot[1, 1] <- cos(theta)
              rot[1, 2] <- -sin(theta)
              rot[2, 1] <- sin(theta)
              rot[2, 2] <- cos(theta)
              work_q_alphabeta[j, 1, i] <- (rot[1, 1] * q_alphabeta[j, 1, i]) + (rot[1, 2] * q_alphabeta[j, 2, i])
              work_q_alphabeta[j, 2, i] <- (rot[2, 1] * q_alphabeta[j, 1, i]) + (rot[2, 2] * q_alphabeta[j, 2, i])
              q_alphabeta[j, 1, i] <- work_q_alphabeta[j, 1, i]
              q_alphabeta[j, 2, i] <- work_q_alphabeta[j, 2, i]
            }
          }
          # x = 1 b = 0
          if (chromosome[i, j] == 1 & chromosome[best_chromosome[generation], j] == 0) {
            if (q_alphabeta[j, 1, i]*q_alphabeta[j, 2, i] >= 0) {
              rot[1, 1] <- cos(theta)
              rot[1, 2] <- -sin(theta)
              rot[2, 1] <- sin(theta)
              rot[2, 2] <- cos(theta)
              work_q_alphabeta[j, 1, i] <- (rot[1, 1] * q_alphabeta[j, 1, i]) + (rot[1, 2] * q_alphabeta[j, 2, i])
              work_q_alphabeta[j, 2, i] <- (rot[2, 1] * q_alphabeta[j, 1, i]) + (rot[2, 2] * q_alphabeta[j, 2, i])
              q_alphabeta[j, 1, i] <- work_q_alphabeta[j, 1, i]
              q_alphabeta[j, 2, i] <- work_q_alphabeta[j, 2, i]
            }
            if (q_alphabeta[j, 1, i]*q_alphabeta[j, 2, i] < 0) {
              rot[1, 1] <- cos(-theta)
              rot[1, 2] <- -sin(-theta)
              rot[2, 1] <- sin(-theta)
              rot[2, 2] <- cos(-theta)
              work_q_alphabeta[j, 1, i] <- (rot[1, 1] * q_alphabeta[j, 1, i]) + (rot[1, 2] * q_alphabeta[j, 2, i])
              work_q_alphabeta[j, 2, i] <- (rot[2, 1] * q_alphabeta[j, 1, i]) + (rot[2, 2] * q_alphabeta[j, 2, i])
              q_alphabeta[j, 1, i] <- work_q_alphabeta[j, 1, i]
              q_alphabeta[j, 2, i] <- work_q_alphabeta[j, 2, i]
            }
          }
          # if (chromosome[i, j] == 1 & chromosome[best_chromosome[generation], j] == 1) {
          #   if (q_alphabeta[j, 1, i]*q_alphabeta[j, 2, i] >= 0) {
          #     rot[1, 1] <- cos(theta)
          #     rot[1, 2] <- -sin(theta)
          #     rot[2, 1] <- sin(theta)
          #     rot[2, 2] <- cos(theta)
          #     work_q_alphabeta[j, 1, i] <- (rot[1, 1] * q_alphabeta[j, 1, i]) + (rot[1, 2] * q_alphabeta[j, 2, i])
          #     work_q_alphabeta[j, 2, i] <- (rot[2, 1] * q_alphabeta[j, 1, i]) + (rot[2, 2] * q_alphabeta[j, 2, i])
          #     q_alphabeta[j, 1, i] <- work_q_alphabeta[j, 1, i]
          #     q_alphabeta[j, 2, i] <- work_q_alphabeta[j, 2, i]
          #   }
          #   if (q_alphabeta[j, 1, i]*q_alphabeta[j, 2, i] < 0) {
          #     rot[1, 1] <- cos(-theta)
          #     rot[1, 2] <- -sin(-theta)
          #     rot[2, 1] <- sin(-theta)
          #     rot[2, 2] <- cos(-theta)
          #     work_q_alphabeta[j, 1, i] <- (rot[1, 1] * q_alphabeta[j, 1, i]) + (rot[1, 2] * q_alphabeta[j, 2, i])
          #     work_q_alphabeta[j, 2, i] <- (rot[2, 1] * q_alphabeta[j, 1, i]) + (rot[2, 2] * q_alphabeta[j, 2, i])
          #     q_alphabeta[j, 1, i] <- work_q_alphabeta[j, 1, i]
          #     q_alphabeta[j, 2, i] <- work_q_alphabeta[j, 2, i]
          #   }
          # }
        }
      }
    }
  }
  return(q_alphabeta)
}
