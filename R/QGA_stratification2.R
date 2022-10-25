#########################################################
#                                                       #
#       QUANTUM GENETIC ALGORITHM - Version             #
#     to optimize sampling frame stratification         #
#                                                       #
#########################################################
pi = 3.1415926535
# devtools::install_github("d4ndo/binaryLogic")
library(binaryLogic)

#---------------------
# ALGORITHM PARAMETERS                                  
#---------------------
popsize <- 20 # Define here the population size 
nbitsol <- 3 # Define here the desired number of strata
Genome <- 150 # Define here the genome length
Genome_el <- 2 # Define here the length of the elementary unit in the genome (number of bits)
              # Note: it must be nbitsol <= 2^Genome_el
generation_max <- 2000 # Define here the maximum number of generations/iterations
thetamax <- pi * 0.05 # Define here the maximum rotation angle
thetamin <- pi * 0.025 # Define here the minimum rotation angle
pop_mutation_rate_max = 1/(popsize+1)  # Define here the mutation rate for each individual
pop_mutation_rate_min = 1/(popsize+1)  # Define here the mutation rate for each individual
mutation_rate_max = 1/(Genome+1)      # Define here the mutation rate for each qubit
mutation_rate_min = 1/(Genome+1)      # Define here the mutation rate for each qubit
mutation_flag = TRUE

# Check
if (nbitsol > 2^Genome_el) stop("Values of nbitsol and Genome_el not compatible")
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

#------------------
# REPAIR PROCEDURE                     
#------------------
repair <- function(chromosome) {
  diff = 2^Genome_el - nbitsol
  for (i in c(1:popsize)) {
    solution1 <- array(chromosome[i,],c(Genome_el,Genome))
    solution <- c(rep(0,Genome))
    for (x in c(1:Genome)) {
      for (y in c(1:Genome_el)) {
        solution[x] <- solution[x] + solution1[y,x]*2^(Genome_el - y) 
      }
    }
    solution <- solution + 1
    table(solution)
    sum(table(solution))
    t <- as.numeric(table(solution))
    length(t)
    if (length(t) > nbitsol) { 
      solution[solution %in% c(which(t==min(t))[1]:length(t)) & !(solution %in% c(1:diff))] <- solution[solution %in% c(which(t==min(t))[1]:length(t)) & !(solution %in% c(1:diff))] - diff
    }
    a = array(c(1:genomeLength),c(Genome_el,Genome))
    for (x in c(1:Genome)) {
      y1 = a[1,x]
      y2 = a[Genome_el,x]
      chromosome[i,c(y1:y2)] <- as.binary(solution[x]-1,n=Genome_el)
    }
  }  
  return(chromosome)
}


#-------------------
# FITNESS EVALUATION                      
#-------------------
eval_fitness <- function(generation) {
  fitness_total <- 0
  sum_sqr <- 0
  fitness_average <- -999999
  variance <- 0
  if (nbitsol < 2^Genome_el) chromosome <- repair(chromosome)
  for (i in c(1:popsize)) {
    solution1 <- array(chromosome[i,],c(Genome_el,Genome))
    solution <- c(rep(0,Genome))
    for (x in c(1:Genome)) {
      for (y in c(1:Genome_el)) {
        solution[x] <- solution[x] + solution1[y,x]*2^(Genome_el - y) 
      }
    }
    solution <- solution + 1
    table(solution)
    # strata <- aggrStrata(strata = atomic_strata, 
    #                      nvar = nvar, 
    #                      vett = solution, 
    #                      censiti = 0,
    #                      dominio = 1)
    strata = aggrStrata2(dataset=frame,
                          model=NULL,
                          vett=solution,
                          dominio=1)
    fitness[i] <- -sum(bethel(strata, cv, realAllocation = TRUE))
    fitness_total <- fitness_total + fitness[i]
  }
  fitness_max <- -999999
  fitness_average <- fitness_total / popsize
  the_best_chromosome <- 0
  for (i in c(1:popsize)) {
    if (fitness[i] >= fitness_max) {
      fitness_max <- fitness[i]
      the_best_chromosome <- i
      solution_max <- chromosome[i, ]
    }
    best_chromosome[generation] <- the_best_chromosome
  }
  return(list(
    fitness = fitness,
    fitness_max = fitness_max,
    fitness_average = fitness_average,
    best_chromosome = best_chromosome,
    solution_max = solution_max
  ))
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

# Prepare data for fitness evaluation
data(iris)
iris$id <- c(1:nrow(iris))
iris$dom <- 1
library(SamplingStrata)
frame <- buildFrameDF(
  df = iris,
  id = "id",
  domainvalue = "dom",
  X = c("id"),
  Y = c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
)
cv <- as.data.frame(list(
  DOM = "DOM1",
  CV1 = 0.03,
  CV2 = 0.03,
  CV3 = 0.03,
  CV4 = 0.03,
  domainvalue = 1
))

nvar = length(grep("Y",colnames(frame)))
source("aggrStrata.R")
source("aggrStrata2.R")
atomic_strata = aggrStrata2(dataset=frame,
                            model=NULL,
                            vett=c(1:150),
                            dominio=1)

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

strata <- aggrStrata2(dataset = frame, 
                      vett = solution, 
                      dominio = 1)
sum(bethel(strata, cv, realAllocation = TRUE))
# [1] 19.75488
iris$stratum <- solution
table(iris$Species, iris$stratum)
#             1  2  3
# setosa      0  0 50
# versicolor 48  2  0
# virginica   4 46  0



########################################################
#                                                      #
# Comparison with Genetic Algorithm (SamplingStrata)   #
#                                                      #
########################################################

sol <- optimStrata(method = "atomic",
            framesamp = frame,
            errors = cv,
            nbitsola = 3,
            pops = 20,
            iter = 1000)
sol$aggr_strata
sum(sol$aggr_strata$SOLUZ)
# # [1] 20.73226
# 
framenew <- sol$framenew
framenew <- merge(framenew,iris[,c("id","Species")],by.x=c("ID"),by.y=c("id"))
table(framenew$Species,framenew$LABEL)
#             1  2  3
# setosa     50  0  0
# versicolor  0 48  2
# virginica   0  5 45
