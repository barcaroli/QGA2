#---------------------------------------------
# Application of the Quantum Genetic Algorithm
# to the Knapsack Problem
#---------------------------------------------

#----------------------
library(QGA)
#----------------------
# Prepare data for fitness evaluation
items <- as.data.frame(list(Item = paste0("item",c(1:300)),
                            weight = rep(NA,300)))
set.seed(1234)
items$weight <- rnorm(300,mean=50,sd=20)
hist(items$weight)
sum(items$weight)
maxweight = sum(items$weight) / 3
maxweight

#----------------------
# Set parameters
popsize = 20
generation_max = 30
nvalues_sol = 2
Genome = nrow(items)
thetainit = 3.1415926535 * 0.05
thetaend = 3.1415926535 * 0.001
pop_mutation_rate_init = 1/(popsize + 1)
pop_mutation_rate_end = 1/(popsize + 1)
mutation_rate_init = 1/(Genome + 1)
mutation_rate_end = 1/(Genome + 1)
mutation_flag = TRUE
eval_fitness = KnapsackProblem
eval_func_inputs = list(items,
                        maxweight)
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
