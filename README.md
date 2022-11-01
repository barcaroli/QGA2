
# Quantum Genetic Algorithm R package

This package allows the application of the Quantum Genetic Algorithm
that was first proposed by Han and Kim in 2000.

In this package, each optimization problem is represented as a
maximization one, where each solution is a sequence of (qu)bits.
Following the quantum paradigm, these qubits are in a superposition
state: when measuring them, they collapse in a 0 or 1 state. After
measurement, the fitness of the solution is calculated as in usual
genetic algorithms.

The evolution at each iteration is oriented by the application of two
quantum gates to the amplitudes of the qubits:

1.  a rotation gate (always);
2.  an X-Pauli gate (optionally).

The rotation is based on the theta angle values: higher values allow a
quicker evolution, and lower values avoid local maxima.

The X-Pauli gate is equivalent to the classical mutation operator and
determines the swap between alfa and beta amplitudes of a given qubit.

The package has been developed in such a way as to permit a complete
separation between the “engine”, and the particular problem object of
the combinatorial optimization. So far, two specific problems have been
implemented:

1.  best stratification of a sampling frame;
2.  traveler salesman problem.

In the “examples” folder, the corresponding two applications are
contained. Once installed the package by executing

``` r
devtools::install_github("barcaroli/QGA")
```

they can be run, and their results analyzed.

For instance, in the case of the optimization of a sampling frame:

``` r
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
pop_mutation_rate_max = 1/(popsize + 1)
pop_mutation_rate_min = 1/(popsize + 1)
mutation_rate_max = 1/(Genome + 1)
mutation_rate_min = 1/(Genome + 1)
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
                pop_mutation_rate_max,
                pop_mutation_rate_min,
                mutation_rate_max,
                mutation_rate_min,
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

#-------------------------------------------------------------
# Comparison with SamplingStrata (classical genetic algorithm)
sol <-optimStrata(method = "atomic",
                  framesamp = frame,
                  nStrata = nstrat,
                  errors = cv,
                  pops = popsize,
                  minnumstr = 1,
                  iter = 500)
sum(sol$aggr_strata$SOLUZ)
```

## Reference

Kuk-Hyun Han and Jong-Hwan Kim, “Genetic quantum algorithm and its
application to combinatorial optimization problem,” Proceedings of the
2000 Congress on Evolutionary Computation. CEC00 (Cat. No.00TH8512),
2000, pp. 1354-1360 vol.2, doi: 10.1109/CEC.2000.870809.
