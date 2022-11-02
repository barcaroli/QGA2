
# QGA (Quantum Genetic Algorithm) R package

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
separation between the “engine”, and the particular problem subject to
combinatorial optimization. So far, three specific problems have been
implemented:

1.  best stratification of a sampling frame;
2.  traveler salesman problem;
3.  knapsack problem.

In the “examples” folder, the corresponding three applications are
contained.

Once installed the package, by executing

``` r
devtools::install_github("d4ndo/binaryLogic")
devtools::install_github("barcaroli/QGA")
```

they can be run, and their results analyzed.

## Reference

Kuk-Hyun Han and Jong-Hwan Kim, “Genetic quantum algorithm and its
application to combinatorial optimization problem,” Proceedings of the
2000 Congress on Evolutionary Computation. CEC00 (Cat. No.00TH8512),
2000, pp. 1354-1360 vol.2, doi: 10.1109/CEC.2000.870809.
