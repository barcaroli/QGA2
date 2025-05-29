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
#' @param progress flag indicating progress bar during iterations
#' @param eval_fitness name of the function that will be used to evaluate the fitness of each solution
#' @param eval_func_inputs specific inputs required by the eval_fitness function
#' @param stop_limit value to stop the iterations if the fitness is higher
#' 
#' @export
#' 
#' @return A numeric vector (positive integers) giving the best solution obtained by the QGA
#' 
#' @examples 
#' #----------------------------------------
#' # Fitness evaluation for Knapsack Problem
#' #----------------------------------------
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
#' # Perform optimization
#' popsize = 20
#' Genome = nrow(items)
#' solutionQGA <- QGA(popsize = 20,
#'                 generation_max = 500,
#'                 nvalues_sol = 2,
#'                 Genome = nrow(items),
#'                 thetainit = 3.1415926535 * 0.05,
#'                 thetaend = 3.1415926535 * 0.025,
#'                 pop_mutation_rate_init = 1/(popsize + 1),
#'                 pop_mutation_rate_end = 1/(popsize + 1),
#'                 mutation_rate_init = 1,
#'                 mutation_rate_end = 1,
#'                 mutation_flag = TRUE,
#'                 plotting = FALSE,
#'                 verbose = FALSE,
#'                 progress = FALSE,
#'                 eval_fitness = KnapsackProblem,
#'                 eval_func_inputs = list(items,
#'                                         maxweight))
#' #----------------------
#' # Analyze results
#' solution <- solutionQGA[[1]]
#' solution <- solution - 1
#' sum(solution)
#' sum(items$weight[solution])
#' maxweight
#' 
 
QGA <- function(
  popsize = 20,
  generation_max = 200,
  nvalues_sol,
  Genome,
  thetainit = 3.1415926535 * 0.05,
  thetaend = 3.1415926535 * 0.025,
  pop_mutation_rate_init = NULL,
  pop_mutation_rate_end = NULL,
  mutation_rate_init = NULL,
  mutation_rate_end = NULL,
  mutation_flag = TRUE,
  plotting = TRUE,
  verbose = TRUE,
  progress = TRUE,
  eval_fitness,
  eval_func_inputs,
  stop_limit = NULL,
  plot_every = 10,      # Plot every n generations
  verbose_every = 10    # Print every n generations
) {
  # --- Parameter checks ---
  if (is.null(nvalues_sol)) stop("nvalues_sol parameter value missing!")
  if (is.null(Genome)) stop("Genome parameter value missing!")
  
  # --- Defaults ---
  if (is.null(pop_mutation_rate_init) & mutation_flag) pop_mutation_rate_init <- 1/(popsize+1)
  if (is.null(pop_mutation_rate_end) & mutation_flag) pop_mutation_rate_end <- 1/(popsize+1)
  if (is.null(mutation_rate_init) & mutation_flag) mutation_rate_init <- 1/(Genome+1)
  if (is.null(mutation_rate_end) & mutation_flag) mutation_rate_end <- 1/(Genome+1)
  if (is.null(stop_limit)) stop_limit <- Inf
  
  # --- Parallel setup ---
  numCores <- parallel::detectCores() - 1
  cl <- parallel::makeCluster(numCores)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))
  
  # --- Precompute constants ---
  theta_step <- (thetainit - thetaend) / generation_max
  pop_mut_step <- (pop_mutation_rate_init - pop_mutation_rate_end) / generation_max
  mut_step <- (mutation_rate_init - mutation_rate_end) / generation_max
  
  # --- Calculate required bits per gene ---
  geneLength <- ceiling(log2(nvalues_sol))
  genomeLength <- Genome * geneLength
  
  # --- Allocate working arrays --
  qubit_0 <- c(1, 0)
  qubit_1 <- c(0, 1)
  q_alphabeta <- array(0.0, c(genomeLength, 2, popsize))
  work_q_alphabeta <- array(0.0, c(genomeLength, 2, popsize))
  chromosome <- matrix(0, popsize, genomeLength)
  best_chromosome <- array(0.0, generation_max)
  h <- matrix(c(1, 1, 1, -1)/sqrt(2), 2, 2)
  rot <- matrix(0.0, 2, 2)
  
  # --- Results dataframe ---
  res <- data.frame(
    generation = seq_len(generation_max + 1),
    fitness_average = numeric(generation_max + 1),
    fitness_best = numeric(generation_max + 1)
  )
  
  fitness_best <- -Inf
  solution_best <- rep(0, genomeLength)
  generation <- 1
  
  # --- Initial Population ---
  theta <- thetainit
  q_alphabeta <- generate_pop(popsize, genomeLength, q_alphabeta, rot, theta, h, qubit_0)
  chromosome <- measure(popsize, genomeLength, q_alphabeta, chromosome)
  chromosome <- repair(popsize, chromosome, geneLength, genomeLength, nvalues_sol, Genome)
  
  # --- Parallel fitness evaluation ---
  fitness <- foreach(i = 1:popsize, .combine = c, .packages = character()) %dopar% {
    eval_fitness(chromosome[i, ], eval_func_inputs)
  }
  fitness_max <- max(fitness)
  fitness_average <- mean(fitness)
  best_idx <- which.max(fitness)
  best_chromosome <- chromosome[best_idx, ]
  solution_max <- chromosome[best_idx, ]
  
  if (fitness_max > fitness_best) {
    fitness_best <- fitness_max
    solution_best <- solution_max
  }
  res$fitness_average[generation] <- fitness_average
  res$fitness_best[generation] <- fitness_best
  if (plotting && (generation %% plot_every == 0)) plot_Output(res[1:generation, ])
  if (verbose && (generation %% verbose_every == 0)) cat("\n", generation, ",", fitness_average, ",", fitness_max)
  
  if (progress) pb <- txtProgressBar(min = 0, max = generation_max, style = 3)
  
  iter <- 1
  while (generation < generation_max && fitness_max < stop_limit) {
    if (progress) setTxtProgressBar(pb, generation)
    generation <- generation + 1
    iter <- iter + 1
    theta <- max(0, thetainit - theta_step * generation)
    pop_mutation_rate <- pop_mutation_rate_init - pop_mut_step * generation
    mutation_rate <- mutation_rate_init - mut_step * generation
    
    q_alphabeta <- rotation(
      chromosome, best_chromosome, generation, genomeLength, solution_best,
      q_alphabeta, work_q_alphabeta, popsize, fitness, theta
    )
    
    if (mutation_flag) {
      q_alphabeta <- mutation(
        pop_mutation_rate, mutation_rate, popsize, chromosome, solution_best,
        q_alphabeta, work_q_alphabeta, genomeLength
      )
    }
    
    chromosome <- measure(popsize, genomeLength, q_alphabeta, chromosome)
    chromosome <- repair(popsize, chromosome, geneLength, genomeLength, nvalues_sol, Genome)
    
    # --- Parallel fitness evaluation ---
    fitness <- foreach(i = 1:popsize, .combine = c, .packages = character()) %dopar% {
      eval_fitness(chromosome[i, ], eval_func_inputs)
    }
    fitness_max <- max(fitness)
    fitness_average <- mean(fitness)
    best_idx <- which.max(fitness)
    best_chromosome <- chromosome[best_idx, ]
    solution_max <- chromosome[best_idx, ]
    
    if (fitness_max > fitness_best) {
      fitness_best <- fitness_max
      solution_best <- solution_max
    }
    res$fitness_average[generation] <- fitness_average
    res$fitness_best[generation] <- fitness_best
    if (plotting && (generation %% plot_every == 0)) plot_Output(res[1:generation, ])
    if (verbose && (generation %% verbose_every == 0)) cat("\n", generation, ",", fitness_average, ",", fitness_best)
    if (fitness_max >= stop_limit) break
  }
  if (progress) close(pb)
  cat("\n *** Best fitness:", fitness_best)
  
  # --- Vectorized decoding of the solution ---
  solution1 <- matrix(solution_best, geneLength, Genome)
  powers <- 2^((geneLength - 1):0)
  solution <- as.integer(colSums(t(solution1) * powers)) + 1
  
  # --- Output ---
  out <- list(
    solution = solution,
    stats = res[1:generation, ]
  )
  return(out)
}
