#---------------------------------------------
# Application of the Quantum Genetic Algorithm
# to the Traveler Salesman Problem
#---------------------------------------------


#----------------------
library(QGA)
#----------------------
# Prepare data for fitness evaluation
cities <- read.csv("cities.csv")
ncities <- 8
cities <- cities[c(1:ncities),]
distance <- as.matrix(dist(cities[,c(2:3)]))
#----------------------
# Set parameters
popsize = 20
generation_max = 2000
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
                plotting = FALSE,
                verbose = FALSE,
                eval_fitness,
                eval_func_inputs)
#----------------------
# Analyze results
cities$city[solution]
cities_tsp <- cities[solution,]
plot(y~x,data=cities_tsp)
polygon(cities_tsp$x,cities_tsp$y,border="red")
text(x = cities_tsp$x, y = cities_tsp$y, labels = cities_tsp$city, cex=.75)
title("Best path")
