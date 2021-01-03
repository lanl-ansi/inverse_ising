# Inverse Ising

This repository contains Julia implementation of RISE, logRISE and RPLE algorithms for the inverse Ising problem.

## Prerequisites

For running the code, you need to have Julia installed on your computer, as well as the Julia packages JuMP, Ipopt and StatsBase packages. The code includes `Project.toml` and `Manifest.toml` files that can be used to configure a Julia environment that the code was tested with.

## Running

Specify desired parameters and file names in the arguments.csv file: reconstruction method (RISE, logRISE, RPLE), regularization coefficient c_lambda, symmetrization of reconstructed couplings (Y, N), name of the input sample file, name of the output parameters file.

Then run `Inverse_Ising.jl` in the command line:

```
julia Inverse_Ising.jl RISE 0.0 N synthetic_example/output_samples.csv output_model.csv
```

The input csv sample file should be in the histogram form, where each line is in the format "number of time a configuration has been sampled, configuration". The output file of reconstructed parameters has the form of a csv matrix that includes couplings (as off-diagonal entries) and magnetic fields (as diagonal entries).

For small systems (e.g. N<=25), samples can be exhaustively generated with `Gibbs_Sampler.jl`:

```
julia Gibbs_Sampler.jl synthetic_example/input_adjacency.csv 1000 output_samples.csv
```

See above for the formats of the adjacency matrix (containing parameters) and of the output file (containing samples in the histogram representation). A small synthetic example of input and output files is provided in the folder "synthetic_example".

## D-Wave data set

The real data set generated on the D-Wave 2X quantum annealer "Ising" at Los Alamos National Laboratory and used in the paper for illustration is available in the folder "data_dwave".

## Reference

If you find this code useful in your work, we kindly request that you cite the following [paper](https://arxiv.org/abs/1612.05024):
* A. Y. Lokhov, M. Vuffray, S. Misra, M. Chertkov (2016). Optimal structure and parameter learning of Ising models.
arXiv preprint arXiv:1612.05024.
```
@article{lokhov2016optimal,
  title={Optimal structure and parameter learning of Ising models},
  author={Lokhov, Andrey Y and Vuffray, Marc and Misra, Sidhant and Chertkov, Michael},
  journal={arXiv preprint arXiv:1612.05024},
  year={2016}
}
```

## License

D-WISC is provided under a BSD-ish license with a "modifications must be indicated" clause. See the [LICENSE.md](LICENSE.md) file for the full text. This package is part of the Hybrid Quantum-Classical Computing suite, known internally as LA-CC-16-032.
