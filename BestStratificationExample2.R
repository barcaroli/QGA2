#---------------------------------------------
# Application of the Quantum Genetic Algorithm
# to the optimization of a sampling frame
#---------------------------------------------

#----------------------
require(SamplingStrata)
library(QGA)
#-------------------------------------------
# Fitness evaluation for best stratification
#-------------------------------------------
BestStratification <- function(solution,
                               eval_func_inputs) {
  frame <- eval_func_inputs[[1]]
  cv <- eval_func_inputs[[2]]
  strata = SamplingStrata::aggrStrata2(dataset=frame,
                                       model=NULL,
                                       vett=solution,
                                       dominio=1)
  fitness <- -sum(SamplingStrata::bethel(strata, cv, realAllocation = TRUE))
  return(fitness)
}

#-------------------------------------------
# Prepare data for fitness evaluation
data(iris)
iris$id <- c(1:nrow(iris))
iris$dom <- 1
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
nstrat = 3

#----------------------
# Set parameters
popsize = 20
generation_max = 1000
nvalues_sol = nstrat
Genome = nrow(iris)
thetainit = 3.1415926535 * 0.05
thetaend = 3.1415926535 * 0.025
pop_mutation_rate_init = 1/(popsize + 1)
pop_mutation_rate_end = 1/(popsize + 1)
mutation_rate_init = 1/(Genome + 1)
mutation_rate_end = 1/(Genome + 1)
mutation_flag = TRUE
eval_fitness = BestStratification
eval_func_inputs = list(frame, cv)

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
                plotting = TRUE,
                verbose = FALSE,
                eval_fitness,
                eval_func_inputs)
#----------------------
# Analyze results
table(solution)
strata <- aggrStrata2(dataset = frame, 
                      vett = solution, 
                      dominio = 1)
sum(bethel(strata, cv, realAllocation = TRUE))
iris$stratum <- solution
table(iris$Species, iris$stratum)

#-------------------------------
# Comparison with genalg
library(genalg)
evaluate <- function(solution) {
  solution <- round(solution)
  strata = SamplingStrata::aggrStrata2(dataset=frame,
                                       model=NULL,
                                       vett=solution,
                                       dominio=1)
  fitness <- sum(SamplingStrata::bethel(strata, cv, realAllocation = TRUE))
  return(fitness)
}
solutionGA <- rbga(stringMin=c(rep(1,nrow(iris))), 
                   stringMax=c(rep(nstrat,nrow(iris))),
                   popSize=20, 
                   iters=500, 
                   elitism=NA, 
                   evalFunc=evaluate)
plot(solutionGA)
filter = solutionGA$evaluations == min(solutionGA$evaluations)
bestObjectCount = sum(rep(1, solutionGA$popSize)[filter])
if (bestObjectCount > 1) {
  bestSolution = solutionGA$population[filter, ][1, 
  ]
} else {
  bestSolution = solutionGA$population[filter, ]
}
bestSolution <- round(bestSolution)
strata <- aggrStrata2(dataset = frame, 
                      vett = bestSolution, 
                      dominio = 1)
sum(bethel(strata, cv, realAllocation = TRUE))
iris$stratum <- bestSolution
table(iris$Species, iris$stratum)
