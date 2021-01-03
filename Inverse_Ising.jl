#!/usr/bin/env julia

###########
# Implementation of RISE, logRISE and RPLE algorithms
# Authors: Andrey Y. Lokhov, Marc Vuffray, Sidhant Misra, Michael Chertkov
# Paper: Optimal structure and parameter learning of Ising models, 2017
# Example of use: julia Inverse_Ising.jl    
###########

using CSV
using DataFrames

using JuMP
using Ipopt

#initialization of arguments (type of method, regularization coefficient, post symmetrization of couplings, input & output files).

adjacency = nothing
if length(ARGS) >= 5
    method              = strip(ARGS[1])
    regularizing_value  = parse(Float64, ARGS[2])
    symmetrization      = strip(ARGS[3])
    file_samples_histo  = strip(ARGS[4])
    file_reconstruction = "$(strip(ARGS[5]))"
    if length(ARGS) >= 6
        adjacency = convert(Matrix{Float64}, CSV.read(ARGS[6], DataFrame, datarow=1))
    end
else
    #reading arguments in the argument file
    args = CSV.read("arguments.csv", DataFrame, datarow=1)
    method              = strip(args[1])
    regularizing_value  = args[2]
    symmetrization      = strip(args[3])
    file_samples_histo  = strip(args[4])
    file_reconstruction = "$(strip(args[5]))"
    if length(args) >= 6
        adjacency = convert(Matrix{Float64}, CSV.read(ARGS[6], DataFrame, datarow=1))
    end
end

#Initialization of the histogram of samples and extraction of number of spins and configuarions.
#Each line of the histogram of samples is in the format "number of time a configuration has been sampled, configuration".
samples_histo       = convert(Matrix{Float64}, CSV.read(file_samples_histo, DataFrame, datarow=1))
(num_conf, num_row) = size(samples_histo)
num_spins           = num_row - 1

#Declaration of the matrix of reconstructed couplings (off-diagional entries) and magnetic fields (diagonal entries).
reconstruction = Array{Float64}(undef, num_spins, num_spins)

#Extraction of the total number of samples contained in the histogram.
num_samples = 0
for k=1:num_conf
    global num_samples += samples_histo[k,1]
end

#Initialization of the regularizing parameter from regularization coefficient. Declaration of the RPLE and RISE objective funcions.
lambda           = regularizing_value*sqrt(log((num_spins^2)/0.05)/num_samples)
RISEobjective(h) = exp(-h)
RPLEobjective(h) = log(1 + exp(-2h))

#Loop for the reconstruction of couplings around "current_spin"
for current_spin = 1:num_spins
    #Printing out progress
    println("Reconstructing the parameters adjacent to node ", current_spin);

    #Construction of the statistics around node "current_spin" in the histogram format.
    #An element (k,i) reads "configuration[k,current_spin] * configuration[k,i]" if i != current_spin and "configuration[k,current_spin]" otherwise.
    nodal_stat  = [ samples_histo[k, 1 + current_spin] * (i == current_spin ? 1 : samples_histo[k, 1 + i]) for k=1:num_conf , i=1:num_spins]

    #Declaration in JuMP of the optimization model "m" and the convex solver. Here the convex solver is choosen to be Ipopt.
    #The tolerance tol is chosen based on experiments in the paper, remove for using a default tolerance of the solver, e.g. solver = IpoptSolver()
    #The option print_level=0 disables the output of the Ipopt solver, remove for reading a default detailed information on the convergence, e.g. solver = IpoptSolver()
    #m = Model(with_optimizer(Ipopt.Optimizer, tol=1e-12, print_level=0))
    m = Model(with_optimizer(Ipopt.Optimizer, tol=1e-12))

    #Initialization of the loss function for JuMP. RISE is choosen by default unless the user specifies RPLE in the arguments.
    JuMP.register(m, :IIPobjective, 1, (method == "RPLE" ? RPLEobjective : RISEobjective), autodiff=true)

    #Declaration in JuMP of "x", the array of variables for couplings and magnetic fields and "z", the array of slack variables for the l1 norm.
    #The magnetic field variable is x[current_spin] and the coupling variables are x[i] for i!= current_spin. The slack variable z[current_spin] is uneccessary.
    @variable(m, x[1:num_spins])
    @variable(m, z[1:num_spins])

    #Declaration in JuMP of the objective function: (log(RISE) + l1, RISE + l1, RPLE + l1). RISE is choosen by default unless the user specifies RPLE or logRISE in the arguments.
    if method == "logRISE"
        @NLobjective(m, Min,
            log(sum((samples_histo[k,1]/num_samples)* IIPobjective(sum(x[i]*nodal_stat[k,i] for i=1:num_spins)) for k=1:num_conf)) +
            lambda*sum(z[j] for j=1:num_spins if current_spin!=j)
            )
    else
        @NLobjective(m, Min,
            sum((samples_histo[k,1]/num_samples)* IIPobjective(sum(x[i]*nodal_stat[k,i] for i=1:num_spins)) for k=1:num_conf) +
            lambda*sum(z[j] for j=1:num_spins if current_spin!=j)
            )
    end

    #Declaration in JuMP of slack constraints for the l1 penalty.
    for j in 1:num_spins
        @constraint(m, z[j] >=  x[j]) #z_plus
        @constraint(m, z[j] >= -x[j]) #z_minus
        if adjacency != nothing && adjacency[current_spin,j] == 0
            @constraint(m, x[j] == 0.0)
            @constraint(m, z[j] == 0.0)
        end
    end

    #Lauching convex optimization, printing results and updating the matrix of reconstructed parameters accordingly.
    optimize!(m)
    spin_values = [JuMP.value(x[i]) for i in 1:num_spins]
    println(current_spin, " = ", spin_values)
    for i in 1:num_spins
        reconstruction[current_spin,i] = spin_values[i]
    end

    # release memory for next iteration
    nodal_stat = 0
    m = 0
    GC.gc()
end

#symmetrization of the couplings. No symmetrization is choosen by defaut unless the user specifies "Y" in the arguments.
if symmetrization == "Y"
    reconstruction = 0.5*(reconstruction + transpose(reconstruction))
end

println(typeof(reconstruction))
println(reconstruction)

println(typeof(file_reconstruction))
println(file_reconstruction)

#ouputing reconstruction in the CSV file
CSV.write(file_reconstruction, DataFrame(reconstruction), writeheader=false)
