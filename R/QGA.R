#' Quantum Genetic Algorithm
#'
#' @description 
#' 
#' Main function to execute a Quantum Genetic Algorithm
#' 
#' @details
#' 
#' This function is the 'engine', which performs the quantum genetic algorithm calling
#' the function for the evaluation of the fitness that is specific for the particulare
#' problem to be optmized.
#' 
#' @param popsize the number of generated solutions (population) to be evaluated at each iteration
#' (default is 20)
#' @param generation_max the number of iterations to be performed
#' (default is 200)
#' @param Genome the length of the genome (or chromosome), representing a possible solution 
#' @param nvalues_sol the number of possible integer values contained in each element (gene) of the solution 
#' @param thetainit the angle (expressed in radiants) to be used when applying the rotation gate
#' when starting the iterations 
#' (default is pi * 0.05, where pi = 3.1415926535)
#' @param thetaend the angle (expressed in radiants) to be used when applying the rotation gate 
#' at the end of the iterations
#' (default is pi * 0.025, where pi = 3.1415926535)
#' @param pop_mutation_rate_init initial mutation rate to be used when applying the X-Pauli gate, applied 
#' to each individual in the population (default is 1/(popsize+1))
#' @param pop_mutation_rate_end final mutation rate to be used when applying the X-Pauli gate, applied 
#' to each individual in the population (default is 1/(popsize+1))
#' @param mutation_rate_init initial mutation rate to be used when applying the X-Pauli gate, applied 
#' to each element of the chromosome  (default is 1/(Genome+1)))
#' @param mutation_rate_end final mutation rate to be used when applying the X-Pauli gate, applied 
#' to each element of the chromosome (default is 1/(Genome+1))
#' @param mutation_flag flag indicating if the mutation gate is to be applied or not (default is TRUE)
#' @param plotting flag indicating plotting during iterations
#' @param verbose flag indicating printing fitness during iterations
#' @param eval_fitness name of the function that will be used to evaluate the fitness of each solution
#' @param eval_func_inputs specific inputs required by the eval_fitness function
#' 
#' @export
#' 
#' @return A numeric vector giving the best solution obtained by the QGA
#' 
#' @examples 
#' #----------------------
#' library(QGA)
#' #----------------------
#' 
#' #----------------------------------------
#' # Fitness evaluation for Knapsack Problem
#----------------------------------------
#' KnapsackProblem <- function(solution,
#'                             eval_func_inputs) {
#'   solution <- solution - 1
#'   items <- eval_func_inputs[[1]]
#'   maxweight <- eval_func_inputs[[2]]
#'   tot_items <- sum(solution)
#'   # Penalization
#'   if (sum(items$weight[solution]) > maxweight) {
#'     tot_items <- tot_items - (sum(items$weight[solution]) - maxweight)  
#'   }
#'   return(tot_items)
#' }
#' 
#' #----------------------------------------
#' # Prepare data for fitness evaluation
#' items <- as.data.frame(list(Item = paste0("item",c(1:300)),
#'                             weight = rep(NA,300)))
#' set.seed(1234)
#' items$weight <- rnorm(300,mean=50,sd=20)
#' hist(items$weight)
#' sum(items$weight)
#' maxweight = sum(items$weight) / 2
#' maxweight
#' #----------------------
#' # Set parameters
#' popsize = 20
#' generation_max = 250
#' nvalues_sol = 2
#' Genome = nrow(items)
#' thetainit = 3.1415926535 * 0.05
#' thetaend = 3.1415926535 * 0.001
#' pop_mutation_rate_init = 1/(popsize + 1)
#' pop_mutation_rate_end = 1/(popsize + 1)
#' mutation_rate_init = 1/(Genome + 1)
#' mutation_rate_end = 1/(Genome + 1)
#' mutation_flag = TRUE
#' eval_fitness = KnapsackProblem
#' eval_func_inputs = list(items,
#'                         maxweight)
#' plotting = TRUE
#' verbose = FALSE
#' 
#' #----------------------
#' # Perform optimization
#' solution <- QGA(popsize,
#'                 generation_max,
#'                 nvalues_sol,
#'                 Genome,
#'                 thetainit,
#'                 thetaend,
#'                 pop_mutation_rate_init,
#'                 pop_mutation_rate_end,
#'                 mutation_rate_init,
#'                 mutation_rate_end,
#'                 mutation_flag,
#'                 plotting,
#'                 verbose,
#'                 eval_fitness,
#'                 eval_func_inputs)
#' 
#' #----------------------
#' # Analyze results
#' solution <- solution - 1
#' sum(solution)
#' sum(items$weight[solution])
#' maxweight
#' 
 
QGA <- function(popsize = 20,  
                generation_max = 200,
                nvalues_sol,
                Genome,
                thetainit = 3.1415926535 * 0.05,
                thetaend = 3.1415926535 * 0.025,
                pop_mutation_rate_init,
                pop_mutation_rate_end,
                mutation_rate_init,
                mutation_rate_end,
                mutation_flag = TRUE,
                plotting = TRUE,
                verbose = TRUE,
                eval_fitness,
                eval_func_inputs) {
  # check
  if (is.null(nvalues_sol)) stop("nvalues_sol parameter value missing!")
  if (is.null(Genome)) stop("Genome parameter value missing!")
  
  # default values
  if (is.null(pop_mutation_rate_init)) pop_mutation_rate_init = 1/(popsize+1)
  if (is.null(pop_mutation_rate_end)) pop_mutation_rate_end = 1/(popsize+1)
  if (is.null(mutation_rate_init)) mutation_rate_init = 1/(Genome+1)
  if (is.null(mutation_rate_end)) mutation_rate_end = 1/(Genome+1)
  
  # Calculate the number of (qu)bits necessary for each element in the genome/chromosome
  n = 0
  while (nvalues_sol > 2^n) {
    n = n+1
  }
  geneLength = n 
  genomeLength <- Genome * geneLength 
  
  #---------------------
  #  WORKING VARIABLES                                  
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
  q_alphabeta <- generate_pop(popsize,
                              genomeLength,
                              q_alphabeta,
                              rot,
                              theta,
                              h,
                              qubit_0)
  chromosome <- measure(popsize,
                        genomeLength,
                        q_alphabeta,
                        chromosome)
  chromosome <- repair(popsize,
                       chromosome,
                       geneLength,
                       genomeLength,
                       nvalues_sol,
                       Genome)
  a <- evaluate(chromosome,
                    best_chromosome,
                    popsize,
                    Genome,
                    geneLength,
                    nvalues_sol,
                    generation,
                    eval_fitness,
                    eval_func_inputs)
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
  if (plotting == TRUE) plot_Output(res[c(1:generation), ])
  if (verbose == TRUE) cat("\n", generation, ",", fitness_average, ",", fitness_max)
  
  if (verbose == FALSE) pb <- txtProgressBar(min = 0, max = generation_max, style = 3)
  while (generation <= generation_max) {
    if (verbose == FALSE) setTxtProgressBar(pb, generation)
    # cat("\n Iteration: ",generation)
    theta <- thetainit - ((thetainit - thetaend) / generation_max) * generation
    # switch_theta = generation_max * 0.25
    # if (generation < switch_theta) theta = thetainit
    # if (generation >= switch_theta) theta = thetaend
    
    if (theta < 0) theta <- 0
    q_alphabeta <- rotation(chromosome,
                            best_chromosome,
                            generation,
                            genomeLength,
                            solution_best,
                            q_alphabeta,
                            work_q_alphabeta,
                            popsize,
                            fitness, 
                            theta)
    generation <- generation + 1
    pop_mutation_rate = pop_mutation_rate_init - ((pop_mutation_rate_init - pop_mutation_rate_end) / generation_max) * generation
    mutation_rate = mutation_rate_init - ((mutation_rate_init - mutation_rate_end) / generation_max) * generation
    if (mutation_flag == TRUE) {
      q_alphabeta <- mutation(pop_mutation_rate, 
                              mutation_rate,
                              popsize,
                              chromosome,
                              solution_best,
                              q_alphabeta,
                              work_q_alphabeta,
                              genomeLength)      
    }
    chromosome <- measure(popsize,
                          genomeLength,
                          q_alphabeta,
                          chromosome)
    chromosome <- repair(popsize,
                         chromosome,
                         geneLength,
                         genomeLength,
                         nvalues_sol,
                         Genome)
    a <- evaluate(chromosome,
                      best_chromosome,
                      popsize,
                      Genome,
                      geneLength,
                      nvalues_sol,
                      generation,
                      eval_fitness,
                      eval_func_inputs)
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
    if (plotting == TRUE) plot_Output(res[c(1:generation), ])
    if (verbose == TRUE) cat("\n", generation, ",", fitness_average, ",", fitness_best)
  }
  if (verbose == FALSE) close(pb)
  cat("\n *** Best fitness: ",fitness_best)
  plot_Output(res)
  solution1 <- array(solution_best,c(geneLength,Genome))
  solution <- c(rep(0,Genome))
  for (x in c(1:Genome)) {
    for (y in c(1:geneLength)) {
      solution[x] <- solution[x] + solution1[y,x]*2^(geneLength - y) 
    }
  }
  solution <- solution + 1
  return(solution)
}
