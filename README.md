
# QGA (Quantum Genetic Algorithm) R package

This package allows the application of the Quantum Genetic Algorithm
that was first proposed by Han and Kim in 2000.

This R development is a derivation of the Python implementation by
Lahoz-Beltra
(<https://github.com/ResearchCodesHub/QuantumGeneticAlgorithms>).

Each optimization problem is represented as a maximization one, where
each solution is a sequence of (qu)bits. Following the quantum paradigm,
these (qu)bits are in a superposition state: when measuring them, they
collapse in a 0 or 1 state. After measurement, the fitness of the
solution is calculated as in usual genetic algorithms.

The evolution at each iteration is oriented by the application of two
quantum gates to the amplitudes of the qubits:

1.  a rotation gate (always);
2.  a Pauli-X gate (optionally).

The rotation is based on the theta angle values: higher values allow a
quicker evolution, and lower values avoid local maxima.

The Pauli-X gate is equivalent to the classical mutation operator and
determines the swap between alfa and beta amplitudes of a given qubit.

The package has been developed in such a way as to permit a complete
separation between the “engine”, and the particular problem subject to
combinatorial optimization. So far, three specific problems have been
implemented, namely:

1.  best stratification of a sampling frame;
2.  traveler salesman problem;
3.  knapsack problem.

In the “inst/docs” folder, the corresponding examples of applications
are contained:

1.  BestStratificationExample.R
2.  TravelerSalesmanExample.R
3.  KnapsackProblemExample.R

Once installed the package, after executing

``` r
devtools::install_github("barcaroli/QGA")
```

they can be run, and their results analyzed.

In particular, QGA can be compared with the traditional genetic
algorithm.

In the case of the best stratification, QGA is compared with the GA
implemented in the package “SamplingStrata”. It can be verified that QGA
converges to a convenient solution more rapidly than SamplingStrata.

In the other cases, QGA is compared with the GA implemented in the
“genalg” package. QGA is more rapidly converging to a good solution in
the knapsack problem, while the opposite is in the case of the traveler
salesman.

## References

Kuk-Hyun Han and Jong-Hwan Kim, “Genetic quantum algorithm and its
application to combinatorial optimization problem,” Proceedings of the
2000 Congress on Evolutionary Computation. CEC00 (Cat. No.00TH8512),
2000, pp.1354-1360 vol.2, doi: 10.1109/CEC.2000.870809.

Lahoz-Beltra, Rafael. 2008. “QuantumGeneticAlgorithms.” GitHub
Repository.
<https://github.com/ResearchCodesHub/QuantumGeneticAlgorithms>.

Lahoz-Beltra, Rafael. 2016. “Quantum Genetic Algorithms for Computer
Scientists.” Computers 5(4). <https://doi.org/10.3390/computers5040024>.

Nowotniak, Robert. 2010. “Survey of Quantum-Inspired Evolutionary
Algorithms.” Technical University of Łódź, Computer Engineering
Department. <https://robert.nowotniak.com/papers/survey2010fimb.pdf>.

Zhang, Gexiang. 2011. “Quantum-Inspired Evolutionary Algorithms: A
Survey and Empirical Study.” Journal of Heuristics 17: 303–51.
<https://doi.org/10.1007/s10732-010-9136-0>
