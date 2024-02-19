#---------------------------------------------
# Application of the Quantum Genetic Algorithm
# to find the minimum/maximum of a function
#---------------------------------------------

library(QGA)
#-------------------------------------------
# Fitness evaluation for best stratification
# Note: -value for the minimum
#        value for the maximum
#-------------------------------------------
FunctionEval <- function(solution,eval_func_inputs) {
  x_values = eval_func_inputs[[1]]
  y_values = eval_func_inputs[[2]]
  x = x_values[solution[1]]
  y = y_values[solution[2]]
  if (x <= length(x_values) & y <= length(y_values)) z = 0.5*(x*(1-x) + y*(1-y)) + 12 * cos(x*y)*sin(2*x+y)
  if (x > length(x_values) | y > length(y_values)) z = -100
  # cat("\nx: ",x,"  y: ",y,"  z: ",z,"\n")
  return(z)
}


#----------------------
# Set specific parameters
step = 0.00001
var1 = seq(-4,4,step)
var2 = seq(-4,4,step)
length(var1)
length(var1)^2
Genome = 2
nvalues_sol = length(var1)
eval_fitness = FunctionEval
eval_func_inputs = list(var1,var2)
n = 0
while (nvalues_sol > 2^n) {
  n = n+1
}
n
# Total number of qubits required:
Genome * n
# Visualize function in 3d 
x = seq(-4,4,0.1) # less values to permit visualization
y = seq(-4,4,0.1)
f = function(x,y) {r = 0.5*(x*(1-x) + y*(1-y)) + 12 * cos(x*y)*sin(2*x+y)}
z = outer(x,y,f)
persp(x, y, z, xlab='X Variable', ylab='Y Variable', zlab='Z Variable',
      main='3D Plot', col='orange', shade=.1, theta = 30, phi = 30, 
      ticktype='detailed')


# Find maximum and minimum of the function (only for a low value of step)
# x = var1 
# y = var2
# f = function(x,y) {r = 0.5*(x*(1-x) + y*(1-y)) + 12 * cos(x*y)*sin(2*x+y)}
# z = outer(x,y,f)
# max(z)
# min(z)


#----------------------
# Perform optimization
times <- 50
resQGA <- rep(NA,times)
popsize = 20
Genome = 2
for (i in c(1:times)) {
  cat("\n i: ",i)
  solution <- QGA(popsize,
                  generation_max = 100,
                  nvalues_sol = length(var1),
                  Genome,
                  thetainit = 3.1415926535 * 0.025,
                  thetaend = 3.1415926535 * 0.0025,
                  pop_mutation_rate_init = 1/(popsize + 1),
                  pop_mutation_rate_end = 1/(popsize + 1),
                  mutation_rate_init = 1/(Genome + 1),
                  mutation_rate_end = 1/(Genome + 1),
                  mutation_flag = TRUE,
                  plotting = FALSE,
                  verbose = FALSE,
                  eval_fitness = FunctionEval,
                  eval_func_inputs = list(var1,var2))
  solution
  x1 = var1[solution[1]]
  y1 = var2[solution[2]]
  z1 = 0.5*(x1*(1-x1) + y1*(1-y1)) + 12 * cos(x1*y1)*sin(2*x1+y1)
  resQGA[i] <- z1
}
resQGA 

#-------------------------------
# Comparison with genalg
library(genalg)
evaluate <- function(solution) {
  x = var1[round(solution[1])]
  y = var2[round(solution[2])]
  z = 0.5*(x*(1-x) + y*(1-y)) + 12 * cos(x*y)*sin(2*x+y)
  # cat("\nx: ",x,"  y: ",y,"  z: ",z,"\n")
  return(-z)
}

resgenalg <- rep(NA,times)
for (i in c(1:times)) {
  cat("\n i: ",i)
  solution_genalg <- rbga(stringMin=c(rep(1,2)), 
                          stringMax=c(rep(nvalues_sol,2)),
                          popSize=20, 
                          iters=100, 
                          elitism=NA, 
                          evalFunc=evaluate)
  plot(solution_genalg)
  filter = solution_genalg$evaluations == min(solution_genalg$evaluations)
  bestObjectCount = sum(rep(1, solution_genalg$popSize)[filter])
  if (bestObjectCount > 1) {
    bestSolution = solution_genalg$population[filter, ][1, 
    ]
  } else {
    bestSolution = solution_genalg$population[filter, ]
  }
  bestSolution <- round(bestSolution)
  bestSolution
  x2 = var1[round(bestSolution[1])]
  y2 = var2[round(bestSolution[2])]
  z2 = 0.5*(x2*(1-x2) + y2*(1-y2)) + 12 * cos(x2*y2)*sin(2*x2+y2)
  resgenalg[i] <- z2
}
resgenalg

results <- as.data.frame(list(algorithm=c(rep("QGA",times),rep("genalg",times)),
                              res=c(resQGA,resgenalg)))
boxplot(res~algorithm,data=results,col="orange")
max(z)
max(resQGA)
max(resgenalg)
summary(resQGA)
summary(resgenalg)
save.image("run_best_stratification.RData")
