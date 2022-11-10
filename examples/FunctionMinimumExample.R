#---------------------------------------------
# Application of the Quantum Genetic Algorithm
# to find the minimum of a function
#---------------------------------------------

#----------------------
library(QGA)
#----------------------
# Prepare data for fitness evaluation
vals <- 64
x <- c(1:vals)
y <- rep(NA,vals)
for (i in c(1:vals)) {
  y[i] <- (x[i]-vals/2)^2
  y[i] <- y[i] + runif(1,min=-vals,max=vals)
}
plot(x,y)
which(y[x]==min(y))

#----------------------
# Fitness evaluation

evaluate <- function(solution,n) {
  valuex <- 0
  for (w in c(1:n)) {
    valuex <- valuex + solution[w]*2^(n-w) 
  }
  valuex <- valuex + 1
  valuey <- y[valuex]
  return(-valuey)
}
n = 0
while (vals > 2^n) {
  n = n+1
}
# solution <- c(0,1,1,1,1,1)
# evaluate(solution,n)


# Set parameters
popsize = 20
generation_max = 100
nvalues_sol = vals
Genome = 1
thetainit = 3.1415926535 * 0.05
thetaend = 3.1415926535 * 0.001
pop_mutation_rate_init = 1/(popsize + 1)
pop_mutation_rate_end = 1/(popsize + 1)
mutation_rate_init = 1/(Genome + 1)
mutation_rate_end = 1/(Genome + 1)
mutation_flag = TRUE
eval_fitness = evaluate
eval_func_inputs = n
plotting = TRUE
verbose = TRUE

#----------------------
# Perform optimization
solution <- QGA(popsize,
                generation_max,
                nvalues_sol,
                Genome,
                thetainit,
                thetaend,
                pop_mutation_rate_init,
                pop_mutation_rate_end,
                mutation_rate_init,
                mutation_rate_end,
                mutation_flag,
                plotting,
                verbose,
                eval_fitness,
                eval_func_inputs)

#----------------------
# Analyze results
solution <- solution - 1
sum(solution)
sum(items$weight[solution])
maxweight

#-----------------------------------------
# Compare with classical genetic algorithm
library(genalg)
evaluate <- function(solution) {
  tot_items <- sum(solution)
  # Penalization
  if (sum(items$weight[solution]) > maxweight) {
    tot_items <- tot_items - (sum(items$weight[solution]) - maxweight)  
  }
  return(-tot_items)
}
rbga.results <- rbga.bin(size=nrow(items),
         suggestions=NULL,
         popSize=20, 
         iters=1000, 
         elitism=NA, 
         zeroToOneRatio=round(maxweight*10/sum(items$weight)),
         evalFunc=evaluate)
plot(rbga.results)
filter = rbga.results$evaluations == min(rbga.results$evaluations)
bestObjectCount = sum(rep(1, rbga.results$popSize)[filter])
if (bestObjectCount > 1) {
  bestSolution = rbga.results$population[filter, ][1, 
  ]
} else {
  bestSolution = rbga.results$population[filter, ]
}
sum(bestSolution)
sum(items$weight[bestSolution])
maxweight
