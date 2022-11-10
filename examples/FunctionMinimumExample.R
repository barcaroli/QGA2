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
  y[i] <- y[i] + runif(1,min=-vals,max=vals) + 1000
}
plot(x,y)
which(y[x]==min(y))
y[which(y[x]==min(y))]


# Set parameters
popsize = 20
generation_max = 20
nvalues_sol = vals
Genome = 1
thetainit = 3.1415926535 * 0.05
thetaend = 3.1415926535 * 0.001
pop_mutation_rate_init = 1/(popsize + 1)
pop_mutation_rate_end = 1/(popsize + 1)
mutation_rate_init = 1/(Genome + 1)
mutation_rate_end = 1/(Genome + 1)
mutation_flag = TRUE
eval_fitness = FunctionMinimum
eval_func_inputs = y
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
solution
y[solution]
which(y[x]==min(y))
y[which(y[x]==min(y))]

#-----------------------------------------
# Compare with classical genetic algorithm
library(genalg)
evaluate <- function(solution) {
  solution <- round(solution)
  value <- y[solution]
  return(value)
}
solutionGA <- rbga(stringMin=c(1), 
                   stringMax=c(vals),
                   popSize=20, 
                   iters=20, 
                   elitism=NA, 
                   evalFunc=evaluate)
plot(solutionGA)
filter = solutionGA$evaluations == min(solutionGA$evaluations)
bestObjectCount = sum(rep(1, solutionGA$popSize)[filter])
if (bestObjectCount > 1) {
  bestSolution = solutionGA$population[filter, ][1]
} else {
  bestSolution = solutionGA$population[filter, ]
}
bestSolution = round(bestSolution)
bestSolution
y[bestSolution]

