library(SamplingStrata)
library(QGA)


# Prepare data for fitness evaluation
data(iris)
iris$id <- c(1:nrow(iris))
iris$dom <- 1
library(SamplingStrata)
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

nstrat = 2

# Execute optimization

solution <- QGA(
  popsize = 20,
  generation_max = 500,
  nvalues_sol = nstrat,
  Genome = nrow(frame),
  Genome_el = 1,
  thetamax = 3.1415926535 * 0.05,
  thetamin = 3.1415926535 * 0.025,
  pop_mutation_rate_max = 1/(popsize + 1),
  pop_mutation_rate_min = 1/(popsize + 1),
  mutation_rate_max = 1/(Genome + 1),
  mutation_rate_min = 1/(Genome + 1),
  mutation_flag = FALSE,
  eval_fitness = best_stratification,
  eval_func_inputs = list(frame, cv)
)

# Analyze results

strata <- aggrStrata2(dataset = frame, 
                      vett = solution, 
                      dominio = 1)
sum(bethel(strata, cv, realAllocation = TRUE))
# [1] 19.75488
iris$stratum <- solution
table(iris$Species, iris$stratum)

