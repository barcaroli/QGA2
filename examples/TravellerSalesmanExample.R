#----------------------
library(QGA)
#----------------------
# Prepare data for fitness evaluation
cities <- read.csv("cities.csv")
distance <- as.matrix(dist(cities[,c(2:3)]))
#----------------------
# Set parameters
popsize = 20
generation_max = 3000
nvalues_sol = nrow(cities)
Genome = nrow(cities)
thetamax = 3.1415926535 * 0.01
thetamin = 3.1415926535 * 0.003
pop_mutation_rate_max = 1/(popsize + 1)
pop_mutation_rate_min = 1/(popsize + 1)
mutation_rate_max = 1/(Genome + 1)
mutation_rate_min = 1/(Genome + 1)
mutation_flag = FALSE
eval_fitness = TravellerSalesman
eval_func_inputs = distance
#----------------------
# Perform optimization
solution <- QGA(popsize,
                generation_max,
                nvalues_sol,
                Genome,
                thetamax,
                thetamin,
                pop_mutation_rate_max,
                pop_mutation_rate_min,
                mutation_rate_max,
                mutation_rate_min,
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
