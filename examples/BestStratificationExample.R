#---------------------------------------------
# Application of the Quantum Genetic Algorithm
# to the optimization of a sampling frame
#---------------------------------------------

#----------------------
library(SamplingStrata)
library(QGA)
#----------------------
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
generation_max = 500
nvalues_sol = nstrat
Genome = 150
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
# Comparison with SamplingStrata
# sol <-optimStrata(method = "atomic",
#                   framesamp = frame,
#                   nStrata = nstrat,
#                   errors = cv,
#                   pops = popsize,
#                   minnumstr = 1,
#                   iter = 500)
# sum(sol$aggr_strata$SOLUZ)
