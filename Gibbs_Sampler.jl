#!/usr/bin/env julia

###########
# Brute-force generation of samples (in the histogram representation) distributed according to Boltzmann Distribution for small systems
# Authors: Andrey Y. Lokhov, Marc Vuffray, Sidhant Misra, Michael Chertkov
# Paper: Optimal structure and parameter learning of Ising models, 2017
# Example of use: julia Gibbs_Sampler.jl input_adjacency.csv 1000000 output_samples.csv  
###########

using StatsBase
using CSV
using DataFrames
using LinearAlgebra

function int_to_spin(int_representation, spin_number)
  spin = 2 .* digits(int_representation, base = 2, pad = spin_number) .- 1
  return spin
end

function weigh_proba(int_representation, adj, prior)
  spin_number  = size(adj,1)
  current_spin = int_to_spin(int_representation, spin_number)
  return exp(((0.5) * current_spin' * adj * current_spin + prior * current_spin)[1])
end

function sample_generation(sample_number, adj, prior)
  spin_number   = size(adj,1)
  config_number = 2^spin_number

  items   = [i for i in 0:(config_number-1)]
  weights = [weigh_proba(i, adj, prior) for i in (0:config_number-1)]

  raw_sample = sample(items, StatsBase.Weights(weights), sample_number)
  raw_binning= countmap(raw_sample)

  spin_sample = [ vcat(raw_binning[i], int_to_spin(i, spin_number)) for i in keys(raw_binning)]
  spin_sample = hcat(spin_sample...)'
  return spin_sample
end

#Reading the input graph adjacency matrix weighted with couplings (off-diagional entries) and magnetic fields (diagonal entries).
file_adj          = ARGS[1]
adjacency_matrix  = convert(Matrix{Float64}, CSV.read(file_adj, DataFrame, datarow=1)) #couplings part
prior_vector      = transpose(diag(adjacency_matrix)) #priors, or magnetic fields part

#Reading number of samples to be generated
number_sample     = parse(Int64,ARGS[2])

#Reading the name of the ouput file
file_samples      = ARGS[3]

#Generation of samples
my_sample         = sample_generation(number_sample, adjacency_matrix, prior_vector)

#Writing generated samples in the histogram form: each line of the histogram of samples is in the format "number of time a configuration has been sampled, configuration"
CSV.write(file_samples, DataFrame(my_sample), writeheader=false)
