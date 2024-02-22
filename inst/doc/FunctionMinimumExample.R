#---------------------------------------------
# Application of the Quantum Genetic Algorithm
# to find the minimum of a function
#---------------------------------------------

#----------------------
library(QGA)
#----------------------

#----------------------
# Fitness evaluation
evaluate1 <- function(solution,eval_func_inputs) {
  y <- eval_func_inputs
  value <- y[solution]
  return(-value)
}

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


#----------------------
# Set parameters
popsize = 20
Genome = 1
thetainit = 3.1415926535 * 0.05
thetaend = 3.1415926535 * 0.001
pop_mutation_rate_init = 1/(popsize + 1)
pop_mutation_rate_end = 1/(popsize + 1)
mutation_rate_init = 1/(Genome + 1)
mutation_rate_end = 1/(Genome + 1)

#----------------------
# Perform optimization
set.seed(1234)
solution <- QGA(popsize=20,
                generation_max=20,
                nvalues_sol=64,
                Genome=1,
                thetainit,
                thetaend,
                pop_mutation_rate_init,
                pop_mutation_rate_end,
                mutation_rate_init,
                mutation_rate_end,
                mutation_flag=FALSE,
                plotting=TRUE,
                verbose=FALSE,
                progress=FALSE,
                eval_fitness=evaluate1,
                eval_func_inputs=y)
QGA:::plot_Output(solution[[2]])

#----------------------
# Analyze results
solution <- solution[[1]]
solution
y[solution]
which(y[x]==min(y))
y[which(y[x]==min(y))]

#-----------------------------------------
# Compare with classical genetic algorithm
library(genalg)
evaluate2 <- function(solution) {
  solution <- round(solution)
  value <- y[solution]
  return(value)
}
solutionGA <- rbga(stringMin=c(1), 
                   stringMax=c(vals),
                   popSize=20, 
                   iters=20, 
                   elitism=NA, 
                   evalFunc=evaluate2)
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

