# Inverse Ising

This repository contains Julia implementation of RISE, logRISE and RPLE algorithms for the inverse Ising problem.
Reference: Andrey Y. Lokhov, Marc Vuffray, Sidhant Misra, Michael Chertkov, "Optimal structure and parameter learning of Ising models" (2017)

### Prerequisites

For running the code, you need to have the latest version of Julia installed on your computer, as well as JuMP, Ipopt and StatsBase packages.

### Running

Specify desired parameters and file names in the arguments.csv file.

Then run `Inverse_Ising.jl` in the command line:

```
julia Inverse_Ising.jl
```

For small systems (e.g. N<=25), samples can be exaustively generated with `Gibbs_Sampler.jl`:

```
julia Gibbs_Sampler.jl input_adjacency.csv num_samples output_samples.csv
```

A small synthetic example of input and output files is provided in the folder "synthetic_example".

## D-Wave data set

The real data set generated on the D-Wave 2X quantum annealer "Ising" at Los Alamos National Laboratory and used in the paper for illustration is available in the folder "data_dwave".

## License

D-WISC is provided under a BSD-ish license with a "modifications must be indicated" clause. See the [LICENSE.md](LICENSE.md) file for the full text. This package is part of the Hybrid Quantum-Classical Computing suite, known internally as LA-CC-16-032.
