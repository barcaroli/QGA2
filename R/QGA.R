#' Quantum Genetic Algorithm
#'
#' @description 
#' 
#' Main function to execute a Quantum Genetic Algorithm
#' 
#' @details
#' 
#' This function is the 'engine', which performs the quantum genetic algorithm calling
#' two more functions that need to be developed or personalized:
#' 1. eval_fitness
#' 2. repair_function
#' 
#' @param popsize the number of generated solutions (population) to be evaluated at each iteration
#' (default is 20)
#' @param generation_max the number of iterations to be performed
#' (default is 200)
#' @param Genome the length of the genome, i.e. the elements of the chromosome representing an 
#' individual in the population 
#' @param Genome_el the length of the single elements in the elements, where each position 
#' represents a qubit
#' @param nvalues_sol the number of possible values contained in each element of the solution 
#' (it must be nvalues_sol <= 2^Genome_el) 
#' @param thetamax the maximum angle (expressed in radiants) to be used when applyin the rotation gate 
#' (default is pi * 0.05, where pi = 3.1415926535)
#' @param thetamin the minimum angle (expressed in radiants) to be used when applyin the rotation gate 
#' (default is pi * 0.025, where pi = 3.1415926535)
#' @param pop_mutation_rate_max maximum mutation rate to be used when applying the X-Pauli gate, applied 
#' to each individual in the population (default is 1/(popsize+1))
#' @param pop_mutation_rate_min minimum mutation rate to be used when applying the X-Pauli gate, applied 
#' to each individual in the population (default is 1/(popsize+1))
#' @param mutation_rate_max maximum mutation rate to be used when applying the X-Pauli gate, applied 
#' to each element of the chromosome  (default is 1/(Genome+1)))
#' @param mutation_rate_min minimum mutation rate to be used when applying the X-Pauli gate, applied 
#' to each element of the chromosome (default is 1/(Genome+1))
#' @param mutation_flag flag indicating if the mutation gate is to be applied or not (default is TRUE)
#' @param eval_fitness name of the function that will be used to evaluate the fitness of each solution
#' @param eval_func_inputs specific inputs required by the eval_fitness function
#' 
#' @export
#' 
#' @return A numeric vector giving the best solution obtained by the QGA
#' 
#' @examples
#' 


QGA <- function(popsize = 20,  
                generation_max = 200,
                nvalues_sol,
                Genome,
                Genome_el,
                thetamax = 3.1415926535 * 0.05,
                thetamin = 3.1415926535 * 0.025,
                pop_mutation_rate_max = 1/(popsize+1),
                pop_mutation_rate_min = 1/(popsize+1),
                mutation_rate_max = 1/(Genome+1),
                mutation_rate_min = 1/(Genome+1),
                mutation_flag = TRUE,
                eval_function,
                eval_func_inputs) {
  
  # Check
  if (nvalues_sol > 2^Genome_el) stop("Values of nvalues_sol and Genome_el not compatible")
  genomeLength <- Genome * Genome_el 
  
  
  #---------------------
  # WORK VARIABLES                                  
  #---------------------
  
  qubit_0 <- array(c(1, 0), c(2, 1))
  qubit_1 <- array(c(0, 1), c(2, 1))
  
  fitness <- array(0.0, c(1, popsize))
  
  q_alphabeta <- array(0.0, c(genomeLength, 2, popsize))
  work_q_alphabeta <- array(0.0, c(genomeLength, 2, popsize))
  
  chromosome <- array(0, c(popsize, genomeLength))
  best_chromosome <- array(0.0, generation_max)
  
  # Hadamard gate
  h <- array(c(1 / sqrt(2.0), 1 / sqrt(2.0), 1 / sqrt(2.0), -1 / sqrt(2.0)), c(2, 2))
  
  # Rotation Q-gate
  rot <- array(0.0, c(2, 2))
  
  #---------------------------
  # POPULATION INITIALIZATION                     
  #---------------------------
  
  generate_pop <- function() {
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
  
  #---------------------------
  # MEASUREMENT                     
  #---------------------------
  measure <- function() {
    for (i in (1:popsize)) {
      for (j in (1:genomeLength)) {
        p_alpha <- runif(1)
        if (p_alpha <= 2*q_alphabeta[j, 1, i]^2) chromosome[i, j] <- 0
        if (p_alpha > 2*q_alphabeta[j, 1, i]^2) chromosome[i, j] <- 1
      }
    }
    return(chromosome)
  }
  
  
  #--------------
  # ROTATION                   
  #--------------
  
  rotation <- function(fitness, theta) {
    rot <- array(0, c(2, 2))
    for (i in c(1:popsize)) {
      if (sum(chromosome[i, ] != solution_best) != 0) {
        for (j in c(1:genomeLength)) {
          if (fitness[i] < fitness[best_chromosome[generation]]) {
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
          if (fitness[i] >= fitness[best_chromosome[generation]]) {
            if (chromosome[i, j] == 0 & chromosome[best_chromosome[generation], j] == 1) {
              if (q_alphabeta[j, 1, i]*q_alphabeta[j, 2, i] >= 0) {
                rot[1, 1] <- cos(-0.05*pi)
                rot[1, 2] <- -sin(-0.05*pi)
                rot[2, 1] <- sin(-0.05*pi)
                rot[2, 2] <- cos(0.05*pi)
                work_q_alphabeta[j, 1, i] <- (rot[1, 1] * q_alphabeta[j, 1, i]) + (rot[1, 2] * q_alphabeta[j, 2, i])
                work_q_alphabeta[j, 2, i] <- (rot[2, 1] * q_alphabeta[j, 1, i]) + (rot[2, 2] * q_alphabeta[j, 2, i])
                q_alphabeta[j, 1, i] <- work_q_alphabeta[j, 1, i]
                q_alphabeta[j, 2, i] <- work_q_alphabeta[j, 2, i]
              }
              if (q_alphabeta[j, 1, i]*q_alphabeta[j, 2, i] < 0) {
                rot[1, 1] <- cos(0.05*pi)
                rot[1, 2] <- -sin(0.05*pi)
                rot[2, 1] <- sin(0.05*pi)
                rot[2, 2] <- cos(0.05*pi)
                work_q_alphabeta[j, 1, i] <- (rot[1, 1] * q_alphabeta[j, 1, i]) + (rot[1, 2] * q_alphabeta[j, 2, i])
                work_q_alphabeta[j, 2, i] <- (rot[2, 1] * q_alphabeta[j, 1, i]) + (rot[2, 2] * q_alphabeta[j, 2, i])
                q_alphabeta[j, 1, i] <- work_q_alphabeta[j, 1, i]
                q_alphabeta[j, 2, i] <- work_q_alphabeta[j, 2, i]
              }
            }
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
          }
        }
      }
    }
    return(q_alphabeta)
  }
  
  #----------
  # MUTATION                   
  #----------
  
  mutation <- function(pop_mutation_rate, mutation_rate) {
    work_q_alphabeta <- q_alphabeta
    for (i in c(1:popsize)) {
      if (sum(chromosome[i, ] != solution_best) != 0) {
        rnd1 <- runif(1)
        if (rnd1 < pop_mutation_rate) {
          for (j in c(1:genomeLength)) {
            rnd2 <- runif(1)
            if (rnd2 < mutation_rate) {
              work_q_alphabeta[j, 1, i] <- q_alphabeta[j, 2, i]
              work_q_alphabeta[j, 2, i] <- q_alphabeta[j, 1, i]
            }
            if (rnd2 >= mutation_rate) {
              work_q_alphabeta[j, 1, i] <- q_alphabeta[j, 1, i]
              work_q_alphabeta[j, 2, i] <- q_alphabeta[j, 2, i]
            }
          }
        }
      }
    }
    q_alphabeta <- work_q_alphabeta
    return(q_alphabeta)
  }
  
  #----------
  # PLOT                   
  #----------
  
  plot_Output <- function(res) {
    y1 = min(min(res$fitness_average),min(res$fitness_best))
    y2 = max(max(res$fitness_average),max(res$fitness_best))
  
    if (y1 > 0) {
      ymin = y1*0.8
      ymax = y2*1.2
    }
    
    if (y1 < 0 & (abs(y1) > abs(y2))) {
      ymin = y1*1.2
      ymax = y2*0.8
    }
  
    plot(fitness_average ~ generation,
      type = "l", data = res, col = "red",
      ylim = (c(ymin, ymax)), ylab = "Fitness"
    )
    title("QGA - Best stratification")
    points(fitness_best ~ generation, type = "l", data = res, col = "black")
    legend("bottomright",
      legend = c("Best fitness: BLACK", "Average fitness: RED"),
      ncol = 1, cex = 0.8, text.font = 1
    )
  }  
  

  #----------
  # EXECUTION                   
  #----------

  res <- NULL
  res$generation <- c(1:(generation_max + 1))
  res$fitness_average <- rep(0, (generation_max + 1))
  res$fitness_best <- rep(0, (generation_max + 1))
  res <- as.data.frame(res)
  
  fitness_best <- -999999
  solution_best <- rep(0, genomeLength)
  generation <- 1
  q_alphabeta <- generate_pop()
  # Show_population()
  chromosome <- measure()
  a <- eval_fitness(generation)
  fitness <- a$fitness
  fitness_max <- a$fitness_max
  fitness_average <- a$fitness_average
  best_chromosome <- a$best_chromosome
  solution_max <- a$solution_max
  if (fitness_max > fitness_best) {
    fitness_best <- fitness_max
    solution_best <- solution_max
  }
  res$fitness_average[generation] <- fitness_average
  res$fitness_best[generation] <- fitness_best
  plot_Output(res[c(1:generation), ])
  cat("\n***", generation, ",", fitness_average, ",", fitness_max,",",length(unique(solution_best)))
  
  while (generation <= generation_max) {
    theta <- thetamax - ((thetamax - thetamin) / generation_max) * generation
    if (theta < 0) theta <- 0
    q_alphabeta <- rotation(fitness, theta)
    # pop_mutation_rate <- pop_mutation_rate_min + ((pop_mutation_rate_max - pop_mutation_rate_min) / generation_max) * generation
    # mutation_rate <- mutation_rate_min + ((mutation_rate_max - mutation_rate_min) / generation_max) * generation
    pop_mutation_rate = pop_mutation_rate_max - ((pop_mutation_rate_max - pop_mutation_rate_min) / generation_max) * generation
    mutation_rate = mutation_rate_max - ((mutation_rate_max - mutation_rate_min) / generation_max) * generation
    if (mutation_flag == TRUE) q_alphabeta <- mutation(pop_mutation_rate, mutation_rate)
    generation <- generation + 1
    chromosome <- measure()
    chromosome <- repair(chromosome)
    a <- eval_fitness(generation)
    fitness <- a$fitness
    fitness_max <- a$fitness_max
    fitness_average <- a$fitness_average
    best_chromosome <- a$best_chromosome
    solution_max <- a$solution_max
    if (fitness_max > fitness_best) {
      fitness_best <- fitness_max
      solution_best <- solution_max
    }
    res$fitness_average[generation] <- fitness_average
    res$fitness_best[generation] <- fitness_best
    plot_Output(res[c(1:generation), ])
    cat("\n", generation, ",", fitness_average, ",", fitness_best)
  }
  
  solution1 <- array(solution_best,c(Genome_el,Genome))
  solution <- c(rep(0,Genome))
  for (x in c(1:Genome)) {
    for (y in c(1:Genome_el)) {
      solution[x] <- solution[x] + solution1[y,x]*2^(Genome_el - y) 
    }
  }
  solution <- solution + 1
  return(solution)
}
