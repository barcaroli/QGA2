#---------------------------------------------
# Application of the Quantum Genetic Algorithm
# to the Traveler Salesman Problem
#---------------------------------------------

#----------------------
library(QGA)
#----------------------

#-------------------------------------------------
# Fitness evaluation for Traveler Salesman Problem
#-------------------------------------------------
TravellerSalesman <- function(solution,distance) {
  l = 0.0  
  for (i in 2:length(solution)) {
    l = l+distance[solution[i-1], solution[i]]
  }
  l = l + distance[solution[1],solution[length(solution)]]
  penal <- ((nrow(distance)) - length(table(solution)))*sum(distance)/10
  cost <- -(penal+l)
  return(cost)
}

#-------------------------------------------------
# Prepare data for fitness evaluation
cities <- read.csv(".\\examples\\cities.csv")
ncities <- 8
cities <- cities[c(1:ncities),]
distance <- as.matrix(dist(cities[,c(2:3)]))
#----------------------

# Set parameters
popsize = 20
generation_max = 1000
nvalues_sol = nrow(cities)
Genome = nrow(cities)
thetainit = 3.1415926535 * 0.01
thetaend = 3.1415926535 * 0.001
pop_mutation_rate_init = 1/(popsize + 1)
pop_mutation_rate_end = 1/(popsize + 1)
mutation_rate_init = 1/(Genome + 1)
mutation_rate_end = 1/(Genome + 1)
mutation_flag = FALSE
eval_fitness = TravellerSalesman
eval_func_inputs = distance

#----------------------
# Perform optimization
solutionQGA <- QGA(popsize,
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
                plotting = FALSE,
                verbose = FALSE,
                eval_fitness,
                eval_func_inputs)

#----------------------
# Analyze results
cities$city[solutionQGA]
cities_tsp <- cities[solutionQGA,]
plot(y~x,data=cities_tsp)
polygon(cities_tsp$x,cities_tsp$y,border="red")
text(x = cities_tsp$x, y = cities_tsp$y, labels = cities_tsp$city, cex=.75)
title("Best path")

#-----------------------------------------
# Compare with classical genetic algorithm

library(genalg)
evaluate <- function(solution) {
  solution <- round(solution)
  l = 0.0  
  for (i in 2:length(solution)) {
    l = l+distance[solution[i-1], solution[i]]
  }
  l = l + distance[solution[1],solution[length(solution)]]
  penal <- ((nrow(distance)) - length(table(solution)))*sum(distance)/10
  cost <- penal+l
  return(cost)
}
solutionGA <- rbga(stringMin=c(rep(1,nrow(cities))), 
                     stringMax=c(rep(nrow(cities),nrow(cities))),
                     popSize=20, 
                     iters=1000, 
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
cities$city[bestSolution]
l = 0.0  
for (i in 2:length(bestSolution)) {
  l = l+distance[bestSolution[i-1], bestSolution[i]]
}
l = l + distance[bestSolution[1],bestSolution[length(bestSolution)]]
l
cities_tsp <- cities[bestSolution,]
plot(y~x,data=cities_tsp)
polygon(cities_tsp$x,cities_tsp$y,border="red")
text(x = cities_tsp$x, y = cities_tsp$y, labels = cities_tsp$city, cex=.75)
title("Best path")

